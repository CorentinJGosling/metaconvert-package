# Extracted from test-reverse-ci-reflection.R:143

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
f <- metaConvert:::.phi_to_cor(
    phi = .5, n_sample = 400,
    n_cases_exp = 143, n_controls_exp = 52,
    n_cases_nexp = 41, n_controls_nexp = 164,
    phi_to_cor = "tetrachoric", reverse_phi = FALSE)
r <- metaConvert:::.phi_to_cor(
    phi = .5, n_sample = 400,
    n_cases_exp = 143, n_controls_exp = 52,
    n_cases_nexp = 41, n_controls_nexp = 164,
    phi_to_cor = "tetrachoric", reverse_phi = TRUE)
expect_equal(as.numeric(r[1]), -as.numeric(f[1]), tolerance = 0)
expect_equal(as.numeric(r[3]), -as.numeric(f[4]), tolerance = 0)
expect_equal(as.numeric(r[4]), -as.numeric(f[3]), tolerance = 0)
expect_equal(as.numeric(r[7]), -as.numeric(f[8]), tolerance = 0)
expect_equal(as.numeric(r[8]), -as.numeric(f[7]), tolerance = 0)
expect_lt(as.numeric(r[3]), as.numeric(r[4]))
expect_lt(as.numeric(r[7]), as.numeric(r[8]))
