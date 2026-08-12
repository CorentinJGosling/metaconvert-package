# Extracted from test-r-ci-out-of-bounds.R:48

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "metaConvert", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
flag_of <- function(a, b, c, d) {
  x <- data.frame(study_id = "s", n_cases_exp = a, n_controls_exp = b,
                  n_cases_nexp = c, n_controls_nexp = d)
  summary(convert_df(x, measure = "r", verbose = FALSE), flags = TRUE)$flags_crude[1]
}
fires <- function(a, b, c, d) grepl("escapes the parameter space", flag_of(a, b, c, d))

# test -------------------------------------------------------------------------
e <- es_from_2x2(4, 6, 1, 89, reverse_2x2 = TRUE)
expect_lt(e$r_ci_lo, -1)
expect_gte(e$r, -1)
expect_lte(e$r_ci_lo, e$r)
