# Regression tests for aggregate_df() data-plumbing bugs (independent of
# metafor). Each targets a specific defect fixed in main_aggregate_df.R that the
# original test-agg.R could not catch (it used only pre-sorted single-column
# integer ids, single-column col_* args, no NA, no factor, no 'y' column).

## ---- weighted-mean alignment (unsorted cluster ids) ------------------------

test_that("col_weighted_mean is attached to the correct study when ids are unsorted", {
  # appearance order (Z, A) differs from sorted order (A, Z): the old
  # cbind(unique(agg), ...) mislabelled the values.
  df <- data.frame(
    agg = c("Z", "A", "Z", "A"),
    es  = c(0.1, 0.2, 0.3, 0.4),
    se  = c(0.5, 0.5, 0.5, 0.5),
    n   = c(100, 1, 100, 1),
    W   = c(3, 1, 5, 2)
  )
  truth <- c(A = weighted.mean(c(1, 1), c(1, 2)),
             Z = weighted.mean(c(100, 100), c(3, 5)))  # = A:1, Z:100

  for (dep in c("subgroups", "outcomes")) {
    r <- aggregate_df(df, dependence = dep, agg_fact = "agg", es = "es", se = "se",
                      cor_unit = 0.5, col_weighted_mean = "n", weights = "W")
    got <- setNames(r$n, as.character(r$agg))
    expect_equal(got[names(truth)], truth, info = dep)
    expect_type(r$n, "double")   # cbind used to coerce to character
  }
})

test_that("col_weighted_mean accepts MULTIPLE columns", {
  df <- data.frame(
    agg   = rep(c(1, 2, 3), each = 3),
    es    = rnorm(9), se = abs(rnorm(9)) + 0.3,
    n     = c(10, 20, 30, 40, 50, 60, 70, 80, 90),
    score = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    W     = c(1, 2, 3, 1, 2, 3, 1, 2, 3)
  )
  r <- aggregate_df(df, dependence = "subgroups", agg_fact = "agg",
                    es = "es", se = "se",
                    col_weighted_mean = c("n", "score"), weights = "W")
  expect_true(all(c("n", "score") %in% colnames(r)))
  split_df <- split(df, df$agg)
  exp_n     <- sapply(split_df, function(d) weighted.mean(d$n, d$W))
  exp_score <- sapply(split_df, function(d) weighted.mean(d$score, d$W))
  got <- setNames(seq_len(nrow(r)), as.character(r$agg))
  expect_equal(setNames(r$n, as.character(r$agg))[names(exp_n)], exp_n)
  expect_equal(setNames(r$score, as.character(r$agg))[names(exp_score)], exp_score)
})

## ---- 'y' column collision --------------------------------------------------

test_that("a user column named 'y' survives aggregation (no placeholder collision)", {
  df <- data.frame(
    agg = c(1, 1, 2, 2),
    es  = c(0.1, 0.2, 0.3, 0.4),
    se  = c(0.5, 0.5, 0.5, 0.5),
    y   = c(11, 12, 13, 14)
  )
  for (dep in c("subgroups", "outcomes")) {
    r <- aggregate_df(df, dependence = dep, agg_fact = "agg", es = "es", se = "se",
                      cor_unit = 0.5, col_mean = "y")
    expect_true("y" %in% colnames(r), info = dep)
    got <- setNames(r$y, as.character(r$agg))
    expect_equal(got[c("1", "2")], c("1" = 11.5, "2" = 13.5), info = dep)
  }
})

## ---- factor cluster ids do not inject phantom rows -------------------------

test_that("factor agg_fact with unused levels produces no phantom NA rows", {
  ids_chr <- rep(c("a", "b", "c"), c(2, 3, 2))
  df_chr  <- data.frame(agg = ids_chr, es = rnorm(7), se = abs(rnorm(7)) + 0.3)
  df_fac  <- df_chr
  df_fac$agg <- factor(df_chr$agg, levels = c("a", "b", "c", "d", "e"))

  for (dep in c("subgroups", "outcomes")) {
    r_fac <- aggregate_df(df_fac, dependence = dep, agg_fact = "agg",
                          es = "es", se = "se", cor_unit = 0.5)
    r_chr <- aggregate_df(df_chr, dependence = dep, agg_fact = "agg",
                          es = "es", se = "se", cor_unit = 0.5)
    expect_equal(nrow(r_fac), 3, info = dep)
    expect_false(any(is.na(r_fac$agg)), info = dep)
    # factor and character keys yield identical effect sizes
    expect_equal(setNames(r_fac$es, as.character(r_fac$agg))[c("a", "b", "c")],
                 setNames(r_chr$es, as.character(r_chr$agg))[c("a", "b", "c")],
                 info = dep)
  }
})

## ---- all-NA groups: NA (not 0 / NaN) across every dependence type ----------

test_that("an all-NA aggregated column returns NA consistently for every dependence", {
  # cluster 1 has all-NA 'v'; cluster 2 has values.
  base <- data.frame(
    agg      = c(1, 1, 2, 2),
    time_agg = c(1, 2, 1, 2),
    es       = c(0.2, 0.3, 0.4, 0.5),
    se       = c(0.3, 0.3, 0.3, 0.3),
    v        = c(NA, NA, 10, 20)
  )
  for (dep in c("subgroups", "outcomes", "times")) {
    r_sum  <- aggregate_df(base, dependence = dep, agg_fact = "agg", es = "es",
                           se = "se", cor_unit = 0.5, col_sum = "v")
    r_mean <- aggregate_df(base, dependence = dep, agg_fact = "agg", es = "es",
                           se = "se", cor_unit = 0.5, col_mean = "v")
    s <- setNames(r_sum$v,  as.character(r_sum$agg))
    m <- setNames(r_mean$v, as.character(r_mean$agg))
    expect_true(is.na(s[["1"]]), info = paste(dep, "sum"))   # not 0
    expect_true(is.na(m[["1"]]), info = paste(dep, "mean"))  # not NaN
    expect_equal(s[["2"]], 30, info = paste(dep, "sum"))
    expect_equal(m[["2"]], 15, info = paste(dep, "mean"))
  }
})

test_that("all-NA col_min / col_max return NA without a base-R min/max warning", {
  base <- data.frame(
    agg = c(1, 1, 2, 2),
    es  = c(0.2, 0.3, 0.4, 0.5),
    se  = c(0.3, 0.3, 0.3, 0.3),
    v   = c(NA, NA, 10, 20)
  )
  for (dep in c("outcomes", "times")) {
    dat <- base
    if (dep == "times") dat$time_agg <- c(1, 2, 1, 2)
    expect_silent(
      r <- aggregate_df(dat, dependence = dep, agg_fact = "agg", es = "es",
                        se = "se", cor_unit = 0.5, col_min = "v")
    )
    mn <- setNames(r$v, as.character(r$agg))
    expect_true(is.na(mn[["1"]]), info = dep)
    expect_equal(mn[["2"]], 10, info = dep)
  }
})

## ---- NA es/se coherence ----------------------------------------------------

test_that("a cluster reduced to a single valid ES returns that ES (not NA)", {
  df <- data.frame(
    agg = c(1, 1, 2, 2),
    es  = c(0.3, NA, 0.4, 0.6),
    se  = c(0.2, 0.25, 0.3, 0.3)
  )
  r <- aggregate_df(df, dependence = "subgroups", agg_fact = "agg",
                    es = "es", se = "se")
  got <- setNames(r$es, as.character(r$agg))
  expect_equal(got[["1"]], 0.3)              # the lone non-NA estimate
  # and se is coherent with a single-study estimate, not a 2-study SE
  expect_equal(setNames(r$se, as.character(r$agg))[["1"]], 0.2)
})
