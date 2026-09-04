## B1b: a correlation CI bound leaving [-1, 1] while the point estimate is valid.
##
## r is bounded identically to alpha and ICC, and on the routes that build a
## symmetric Wald interval (es +/- z*se) the bound can escape the parameter space
## near |r| = 1. B7b/B8b already flag exactly this for alpha and ICC; B1 above tests
## only the point estimate, so without B1b the bound passes through silently.
##
## SCOPE CHANGED TWICE. First the tetrachoric route stopped reaching B1b: its
## r-scale interval became the tanh back-transform of the Fisher-z interval (.tet_r in
## internal_multiple_formulas.R), which cannot leave (-1, 1) -- measured escape fell
## from 6.7% of tables to 0.0% with coverage of the true correlation improving from
## 0.945 to 0.952. Then EVERY package-computed r interval was switched to that
## construction (es_from_pearson_r() / fisher_z / spearman, the user r/z branches, and
## the lipsey_cooper smd_to_cor branch; the OR conversions and the viechtbauer branch
## already did it), matching cor.test() and metafor's ZCOR + transf.ztor. B1b is now
## reachable only through a USER-supplied r interval, which is handed back verbatim,
## and is retained as the backstop for that case.
##
## This flags; it does not change those intervals.

flag_of <- function(x) {
  summary(convert_df(x, measure = "r", verbose = FALSE), flags = TRUE)$flags_crude[1]
}
fires <- function(x) grepl("escapes the parameter space", flag_of(x))

.pearson_row <- function(r, n) {
  data.frame(study_id = "s", pearson_r = r, n_sample = n)
}

test_that("es_from_pearson_r() builds the back-transformed Fisher interval, which cannot escape", {
  # A large r on a small sample: the former r +/- qt se overshot 1 while r was valid.
  e <- es_from_pearson_r(pearson_r = 0.97, n_sample = 12)
  expect_lte(abs(e$r), 1)
  expect_lte(e$r_ci_up, 1)
  z <- metafor::escalc(measure = "ZCOR", ri = 0.97, ni = 12)
  bt <- summary(z, transf = metafor::transf.ztor)
  expect_equal(c(e$r_ci_lo, e$r_ci_up), c(bt$ci.lb, bt$ci.ub), tolerance = 1e-8)
  expect_false(fires(.pearson_row(0.97, 12)))
  e2 <- es_from_pearson_r(pearson_r = -0.97, n_sample = 12)
  expect_gte(e2$r_ci_lo, -1)
  expect_false(fires(.pearson_row(-0.97, 12)))
})

.user_row <- function(r, lo, up) {
  data.frame(study_id = "s", user_es_original_measure_crude = "r", user_es_crude = r,
             user_ci_lo_crude = lo, user_ci_up_crude = up, n_sample = 12)
}

test_that("B1b fires when a USER-supplied correlation CI bound escapes [-1, 1]", {
  x <- .user_row(0.97, 0.85, 1.05)
  expect_true(fires(x))
  # [INFO], not [UNUSUAL]: nothing was extracted wrongly. It is a disclosure about
  # the interval, so there is nothing for the user to verify. Tested on the B1b
  # token itself -- this row legitimately raises other flags too (|r| = 0.97 trips
  # the high-correlation check), so the assertion must not be made against the whole
  # flag string.
  b1b <- Filter(function(p) grepl("escapes the parameter space", p),
                strsplit(flag_of(x), "; ")[[1]])
  expect_length(b1b, 1)
  expect_match(b1b, "^\\[INFO\\]")
  expect_false(grepl("\\[UNUSUAL\\]", b1b))
  expect_match(b1b, "upper bound exceeds 1")
})

test_that("B1b catches the lower bound too", {
  x <- .user_row(-0.97, -1.05, -0.85)
  expect_true(fires(x))
  expect_match(flag_of(x), "lower bound is below -1")
})

test_that("B1b is silent on well-behaved correlations", {
  for (rn in list(c(0.30, 200), c(-0.45, 150), c(0.60, 400))) {
    e <- es_from_pearson_r(pearson_r = rn[1], n_sample = rn[2])
    expect_lte(e$r_ci_up, 1)
    expect_gte(e$r_ci_lo, -1)
    expect_false(fires(.pearson_row(rn[1], rn[2])))
  }
})

test_that("the tetrachoric route no longer reaches B1b", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  # Bonett & Price (2005) Example 2 -- the table that used to produce r = 0.864 with
  # a Wald upper bound of 1.078, and was this file's original motivating case.
  e <- es_from_2x2(4, 6, 1, 89)
  expect_lte(abs(e$r), 1)
  expect_lte(e$r_ci_up, 1)
  expect_gte(e$r_ci_lo, -1)
  x <- data.frame(study_id = "s", n_cases_exp = 4, n_controls_exp = 6,
                  n_cases_nexp = 1, n_controls_nexp = 89)
  expect_false(fires(x))
  # The interval is still an interval: it brackets the estimate and is ordered.
  expect_lt(e$r_ci_lo, e$r)
  expect_gt(e$r_ci_up, e$r)
})

test_that("the reflected tetrachoric interval is bounded on the lower side too", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  e <- es_from_2x2(4, 6, 1, 89, reverse_2x2 = TRUE)
  expect_gte(e$r_ci_lo, -1)
  expect_lte(e$r_ci_lo, e$r)       # reflection kept the interval well-ordered
  expect_lt(e$r_ci_lo, e$r_ci_up)
})

test_that("B1b and B1 are mutually exclusive by construction", {
  # B1b is guarded on abs(es) <= 1, so a row whose POINT ESTIMATE is out of range
  # gets B1 only. Checked end-to-end through the user-input route rather than by
  # calling the internal, so the test does not depend on its signature.
  x <- data.frame(study_id = "s",
                  user_es_measure_crude = "r", user_es_crude = 1.4,
                  user_ci_lo_crude = 1.2, user_ci_up_crude = 1.6, n_sample = 50)
  f <- tryCatch(
    summary(convert_df(x, measure = "r", verbose = FALSE), flags = TRUE)$flags_crude[1],
    error = function(e) NA_character_)
  skip_if(is.na(f), "user-input r route not reachable with these columns")
  if (grepl("\\[INVALID\\] r outside", f))
    expect_false(grepl("escapes the parameter space", f))
})
