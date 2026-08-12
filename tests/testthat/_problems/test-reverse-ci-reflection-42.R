# Extracted from test-reverse-ci-reflection.R:42

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
fwd <- es_from_2x2(143, 52, 41, 164, table_2x2_to_cor = "tetrachoric",
                     reverse_2x2 = FALSE)
rev <- es_from_2x2(143, 52, 41, 164, table_2x2_to_cor = "tetrachoric",
                     reverse_2x2 = TRUE)
expect_equal(fwd$r, 0.74587459, tolerance = 1e-7)
expect_equal(fwd$r_ci_lo, 0.65912073, tolerance = 1e-7)
expect_equal(fwd$r_ci_up, 0.83262844, tolerance = 1e-7)
expect_equal(rev$r, -0.74587459, tolerance = 1e-7)
expect_equal(rev$r_ci_lo, -0.83262844, tolerance = 1e-7)
expect_equal(rev$r_ci_up, -0.65912073, tolerance = 1e-7)
