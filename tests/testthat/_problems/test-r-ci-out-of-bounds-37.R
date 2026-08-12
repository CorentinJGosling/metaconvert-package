# Extracted from test-r-ci-out-of-bounds.R:37

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
for (tab in list(c(143, 52, 41, 164), c(203, 186, 167, 374), c(30, 30, 20, 40))) {
    e <- es_from_2x2(tab[1], tab[2], tab[3], tab[4])
    expect_lte(e$r_ci_up, 1)
    expect_gte(e$r_ci_lo, -1)
    expect_false(fires(tab[1], tab[2], tab[3], tab[4]))
  }
