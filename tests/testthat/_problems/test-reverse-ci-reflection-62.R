# Extracted from test-reverse-ci-reflection.R:62

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
g <- tetra_grid()
fwd <- suppressWarnings(es_from_2x2(g$a, g$b, g$c, g$d,
                                      table_2x2_to_cor = "tetrachoric",
                                      reverse_2x2 = rep(FALSE, nrow(g))))
rev <- suppressWarnings(es_from_2x2(g$a, g$b, g$c, g$d,
                                      table_2x2_to_cor = "tetrachoric",
                                      reverse_2x2 = rep(TRUE, nrow(g))))
ok <- !is.na(fwd$r) & !is.na(fwd$r_ci_lo) & !is.na(fwd$r_ci_up) &
        !is.na(rev$r) & !is.na(rev$r_ci_lo) & !is.na(rev$r_ci_up) &
        !is.na(fwd$z_ci_lo) & !is.na(rev$z_ci_lo)
expect_gt(sum(ok), 200)
