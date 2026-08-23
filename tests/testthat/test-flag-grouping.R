# Tests for flag_group: scoping the CROSS-ROW quality checks by a grouping
# variable (multivariate / multi-outcome meta-analysis). Default flag_group = NULL
# must reproduce the whole-dataset behaviour exactly; when set, cross-row checks
# (D1/D2/D3, G, H, E4/E6/E7) run WITHIN each group.

# A trial contributing several outcomes repeats its study_id across outcome groups
# (legitimate); within one outcome, study S1 appears twice (a true duplicate). Row 4
# is a large within-'dep' outlier.
grp_dat <- function() data.frame(
  study_id   = c("S1", "S2", "S3", "S1", "S1", "S2", "S3"),
  outcome    = c("dep", "dep", "dep", "dep", "anx", "anx", "anx"),
  mean_exp   = c(10, 11, 10.5, 55, 20, 21, 20.5),
  mean_nexp  = c(8, 9, 8.5, 9, 18, 19, 18.5),
  mean_sd_exp  = rep(3, 7), mean_sd_nexp = rep(3, 7),
  n_exp = rep(30, 7), n_nexp = rep(30, 7)
)

fl <- function(dat, fo = list()) {
  s <- summary(convert_df(dat, measure = "g", verbose = FALSE),
               flags = TRUE, flag_options = fo)
  as.character(s$flags_crude)
}

test_that(".build_group_key: NULL / absent column -> one pool; else the column value", {
  d <- data.frame(o = c("a", "a", "b"), s = c("x", "y", "x"))
  expect_equal(metaConvert:::.build_group_key(d, NULL, n = 3), rep("__all__", 3))
  expect_equal(metaConvert:::.build_group_key(d, "missing_col", n = 3), rep("__all__", 3))
  expect_equal(metaConvert:::.build_group_key(d, "o"), c("a", "a", "b"))
  # multi-column interaction. The INTERNAL separator is "\r", not the " / " the
  # user is shown -- see the collision test below for why.
  expect_equal(metaConvert:::.build_group_key(d, c("o", "s")),
               c("a\rx", "a\ry", "b\rx"))
  # ...and it still READS as " / " wherever a message prints it
  expect_equal(metaConvert:::.group_label(
                 metaConvert:::.build_group_key(d, c("o", "s"))),
               c("a / x", "a / y", "b / x"))
  # NA -> shared "<NA>" bucket
  d2 <- data.frame(o = c("a", NA, NA))
  expect_equal(metaConvert:::.build_group_key(d2, "o"), c("a", "<NA>", "<NA>"))
})

test_that(".build_group_key: a value containing ' / ' cannot collide two distinct rows", {
  # Audit 2.0.1 finding #40. Joining multi-column keys with " / " -- a sequence a
  # grouping value can itself contain -- merged distinct rows into one pool, so
  # D1/D2/D3, G, H, E4 and E6/E7/E8 silently compared rows that are not comparable.
  #   (g1 = "x / y", g2 = "z")  and  (g1 = "x", g2 = "y / z")  both gave "x / y / z".
  d <- data.frame(g1 = c("x / y", "x"), g2 = c("z", "y / z"))
  k <- metaConvert:::.build_group_key(d, c("g1", "g2"))
  expect_false(identical(k[1], k[2]))            # was TRUE: the whole bug
  expect_equal(length(unique(k)), 2L)            # two pools, not one
})

test_that(".group_label: neutralises both the internal separator and the flag-merge separator", {
  # "; " is the sequence the merged flag string is split on. A study_id or group
  # value carrying it truncates the message and, for the re-split Tier-1 checks,
  # mis-routes the leading fragment into BOTH the crude and adjusted scopes
  # (it carries no quoted column name, so it matches every scope). Audit #41.
  expect_equal(metaConvert:::.group_label("HAM-D; total score"), "HAM-D, total score")
  expect_equal(metaConvert:::.group_label("a\rb"), "a / b")
  expect_equal(metaConvert:::.group_label("a\rb; c"), "a / b, c")
})

test_that(".row_ref: sanitises a study_id containing '; ' (covers V23, V36 and H)", {
  expect_equal(metaConvert:::.row_ref(1, "Huang; 2017"), "Huang, 2017 (row 1)")
  expect_false(grepl("; ", metaConvert:::.row_ref(1, "Huang; 2017"), fixed = TRUE))
  # unaffected when there is nothing to sanitise
  expect_equal(metaConvert:::.row_ref(2, c("A", "Xu 2017")), "Xu 2017 (row 2)")
  # falls back to the bare row number with no study_id
  expect_equal(metaConvert:::.row_ref(3, NULL), "row 3")
})

test_that(".by_group: NULL-equivalent single group leaves messages untagged; multi tags them", {
  gk1 <- rep("__all__", 3)
  out1 <- metaConvert:::.by_group(gk1, 3, function(idx) lapply(idx, function(i) "msg"))
  expect_equal(unlist(out1), rep("msg", 3))                 # no group tag with 1 group
  gk2 <- c("a", "a", "b")
  out2 <- metaConvert:::.by_group(gk2, 3, function(idx) lapply(seq_along(idx), function(j) "msg"))
  expect_true(all(grepl("within group '", unlist(out2))))   # tagged with >1 group
  expect_true(any(grepl("within group 'a'", out2[[1]])))
  expect_true(any(grepl("within group 'b'", out2[[3]])))
})

test_that("H duplication: whole dataset flags every repeated study_id (multivariate false positive)", {
  f <- fl(grp_dat())                                   # flag_group = NULL
  # S1 (rows 1,4,5), S2 (2,6), S3 (3,7) all repeat across outcomes -> all 7 flagged
  expect_equal(sum(grepl("Duplicate study_id", f)), 7)
})

test_that("H duplication: grouped by outcome flags only the WITHIN-group repeat", {
  f <- fl(grp_dat(), list(flag_group = "outcome"))
  dup <- grepl("Duplicate study_id", f)
  # Only S1 twice within 'dep' (rows 1 and 4) qualifies
  expect_equal(which(dup), c(1L, 4L))
  expect_true(all(grepl("within group 'dep'", f[c(1, 4)])))
  # S2/S3 (once per outcome) and S1-in-anx are NOT duplicates
  expect_false(any(dup[c(2, 3, 5, 6, 7)]))
})

test_that("D1 ES outlier is judged within the group", {
  f <- fl(grp_dat(), list(flag_group = "outcome"))
  # row 4 (g huge) is the outlier within the 4-row 'dep' group
  expect_true(grepl("ES outlier", f[4]))
  expect_true(grepl("within group 'dep'", f[4]))
  # the 3 anx rows are homogeneous -> no ES-outlier flag among them
  expect_false(any(grepl("ES outlier", f[5:7])))
})

test_that("multi-column flag_group works end to end", {
  d <- grp_dat(); d$subgroup <- "a"
  f <- fl(d, list(flag_group = c("outcome", "subgroup")))
  expect_equal(which(grepl("Duplicate study_id", f)), c(1L, 4L))
})

test_that("default flag_group = NULL is unchanged vs not naming the option at all", {
  d <- grp_dat()
  expect_identical(fl(d, list()), fl(d, list(flag_group = NULL)))
})
