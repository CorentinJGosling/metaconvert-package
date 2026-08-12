# Extracted from test-mapply-memo.R:58

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "metaConvert", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
res <- es_from_2x2(n_cases_exp     = c(15, 15),
                     n_controls_exp  = c(25, 25),
                     n_cases_nexp    = c(8, 8),
                     n_controls_nexp = c(32, 32),
                     reverse_2x2     = c(FALSE, TRUE))
expect_equal(res$r[1], -res$r[2], tolerance = 1e-12)
expect_false(isTRUE(all.equal(res$r[1], res$r[2])))
