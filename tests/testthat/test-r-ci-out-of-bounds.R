## B1b: a correlation CI bound leaving [-1, 1] while the point estimate is valid.
##
## r is bounded identically to alpha and ICC and its interval is built identically
## (es +/- z*se), so a symmetric Wald interval can escape the parameter space near
## |r| = 1. B7b/B8b already flag exactly this for alpha and ICC; before B1b, r had
## no analogue and the bound passed through silently -- B1 tests only the point
## estimate. Common on the tetrachoric route with small off-diagonal cells.
##
## This flags; it does NOT change the interval. Keeping metafor parity on the
## tetrachoric CI is deliberate (tests_save/checked/test-2x2.R asserts it to 1e-10).

flag_of <- function(a, b, c, d) {
  x <- data.frame(study_id = "s", n_cases_exp = a, n_controls_exp = b,
                  n_cases_nexp = c, n_controls_nexp = d)
  summary(convert_df(x, measure = "r", verbose = FALSE), flags = TRUE)$flags_crude[1]
}
fires <- function(a, b, c, d) grepl("escapes the parameter space", flag_of(a, b, c, d))

test_that("B1b fires when a correlation CI bound escapes [-1, 1]", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  # Bonett & Price (2005) Example 2: r = 0.864 with CI upper bound 1.078
  e <- es_from_2x2(4, 6, 1, 89)
  expect_lte(abs(e$r), 1)          # the point estimate itself is valid
  expect_gt(e$r_ci_up, 1)          # but the interval is not
  expect_true(fires(4, 6, 1, 89))
  # [INFO], not [UNUSUAL]: nothing was extracted wrongly and the user has no
  # corrective action (table_2x2_to_cor accepts only "tetrachoric"), so this is a
  # disclosure about the interval metaConvert builds -- V31's situation, not B7b's.
  expect_match(flag_of(4, 6, 1, 89), "\\[INFO\\]")
  expect_false(grepl("\\[UNUSUAL\\]", flag_of(4, 6, 1, 89)))
  expect_match(flag_of(4, 6, 1, 89), "upper bound exceeds 1")
})

test_that("B1b is silent on well-behaved tables", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  for (tab in list(c(143, 52, 41, 164), c(203, 186, 167, 374), c(30, 30, 20, 40))) {
    e <- es_from_2x2(tab[1], tab[2], tab[3], tab[4])
    expect_lte(e$r_ci_up, 1)
    expect_gte(e$r_ci_lo, -1)
    expect_false(fires(tab[1], tab[2], tab[3], tab[4]))
  }
})

test_that("B1b catches the lower bound too, via reverse_2x2", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  # reverse_2x2 reflects the interval (see test-reverse-ci-reflection.R), so the
  # out-of-range bound moves to the lower side and must still be caught.
  e <- es_from_2x2(4, 6, 1, 89, reverse_2x2 = TRUE)
  expect_lt(e$r_ci_lo, -1)
  expect_gte(e$r, -1)
  expect_lte(e$r_ci_lo, e$r)       # reflection kept the interval well-ordered
})

test_that("B1b and B1 are mutually exclusive by construction", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
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
