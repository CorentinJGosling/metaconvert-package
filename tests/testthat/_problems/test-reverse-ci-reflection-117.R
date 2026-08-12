# Extracted from test-reverse-ci-reflection.R:117

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "metaConvert", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
tetra_grid <- function() {
  g <- expand.grid(
    a = c(3, 12, 40, 143, 220),
    b = c(2, 9, 52, 130),
    c = c(1, 15, 41, 160),
    d = c(4, 20, 164, 300)
  )
  g[, c("a", "b", "c", "d")]
}

# test -------------------------------------------------------------------------
fwd <- es_from_2x2_sum(n_cases_exp = 143, n_exp = 195,
                         n_cases_nexp = 41, n_nexp = 205,
                         table_2x2_to_cor = "tetrachoric", reverse_2x2 = FALSE)
rev <- es_from_2x2_sum(n_cases_exp = 143, n_exp = 195,
                         n_cases_nexp = 41, n_nexp = 205,
                         table_2x2_to_cor = "tetrachoric", reverse_2x2 = TRUE)
expect_equal(rev$r_ci_lo, -fwd$r_ci_up, tolerance = 0)
expect_equal(rev$r_ci_up, -fwd$r_ci_lo, tolerance = 0)
expect_lt(rev$r_ci_lo, rev$r_ci_up)
