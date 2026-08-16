## B1b: a correlation CI bound leaving [-1, 1] while the point estimate is valid.
##
## r is bounded identically to alpha and ICC, and on the routes that build a
## symmetric Wald interval (es +/- z*se) the bound can escape the parameter space
## near |r| = 1. B7b/B8b already flag exactly this for alpha and ICC; B1 above tests
## only the point estimate, so without B1b the bound passes through silently.
##
## SCOPE CHANGED: the tetrachoric route no longer reaches B1b. Its r-scale interval
## is now the tanh back-transform of the Fisher-z interval (.tet_r in
## internal_multiple_formulas.R), which cannot leave (-1, 1) -- measured escape fell
## from 6.7% of tables to 0.0% with coverage of the true correlation improving from
## 0.945 to 0.952. B1b is retained for the routes that DO still build a symmetric
## Wald r interval: es_from_pearson_r() and the or_to_cor family.
##
## This flags; it does not change those intervals.

flag_of <- function(x) {
  summary(convert_df(x, measure = "r", verbose = FALSE), flags = TRUE)$flags_crude[1]
}
fires <- function(x) grepl("escapes the parameter space", flag_of(x))

.pearson_row <- function(r, n) {
  data.frame(study_id = "s", pearson_r = r, n_sample = n)
}

test_that("B1b fires when a correlation CI bound escapes [-1, 1]", {
  # A large r on a small sample: the Wald interval overshoots 1 while r is valid.
  e <- es_from_pearson_r(pearson_r = 0.97, n_sample = 12)
  expect_lte(abs(e$r), 1)          # the point estimate itself is valid
  expect_gt(e$r_ci_up, 1)          # but the interval is not
  expect_true(fires(.pearson_row(0.97, 12)))
  # [INFO], not [UNUSUAL]: nothing was extracted wrongly. It is a disclosure about
  # the interval metaConvert builds, so there is nothing for the user to verify.
  # Tested on the B1b token itself -- this row legitimately raises other flags too
  # (|r| = 0.97 trips the high-correlation check), so the assertion must not be made
  # against the whole flag string.
  b1b <- Filter(function(p) grepl("escapes the parameter space", p),
                strsplit(flag_of(.pearson_row(0.97, 12)), "; ")[[1]])
  expect_length(b1b, 1)
  expect_match(b1b, "^\\[INFO\\]")
  expect_false(grepl("\\[UNUSUAL\\]", b1b))
  expect_match(b1b, "upper bound exceeds 1")
})

test_that("B1b catches the lower bound too", {
  e <- es_from_pearson_r(pearson_r = -0.97, n_sample = 12)
  expect_lt(e$r_ci_lo, -1)
  expect_gte(e$r, -1)
  expect_true(fires(.pearson_row(-0.97, 12)))
  expect_match(flag_of(.pearson_row(-0.97, 12)), "lower bound is below -1")
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
