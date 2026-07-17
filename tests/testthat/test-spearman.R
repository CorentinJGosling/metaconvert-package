# =============================================================================
# Spearman rho input (es_from_spearman_rho) — promoted from tests_save/checked/
# (P26, audit digest [25]).
#
# The archived file's SE / propagation / reversal cases duplicated what
# tests/testthat/test-spearman-fixes.R now anchors against hand-derived
# published formulas and a metafor::conv.delta oracle; those duplicates are
# dropped here. This file keeps the archived file's UNIQUE cases:
#   - the Rupinski & Dunlap (1996) point conversion (incl. the r_s = 0 and
#     r_s = 1 fixed points), with a frozen literal regression pin,
#   - the Fisher z of the converted r,
#   - vectorization,
#   - convert_df pipeline integration for measure = "r" and measure = "g".
# =============================================================================

test_that("es_from_spearman_rho converts to Pearson r (Rupinski & Dunlap 1996)", {
  # published transform: r_p = 2 * sin(pi/6 * r_s)
  res <- es_from_spearman_rho(spearman_r = 0.55, n_sample = 100)
  expect_equal(res$r, 2 * sin(pi / 6 * 0.55), tolerance = 1e-10)
  # frozen literal pin (computed externally) so a transform regression cannot
  # hide behind a re-typed formula:
  expect_equal(res$r, 0.5680306894078453, tolerance = 1e-12)

  # fixed points of the transform: r_s = 0 -> 0, r_s = 1 -> 2*sin(pi/6) = 1
  res0 <- es_from_spearman_rho(spearman_r = 0, n_sample = 50)
  expect_equal(res0$r, 0, tolerance = 1e-10)
  res1 <- es_from_spearman_rho(spearman_r = 1, n_sample = 50)
  expect_equal(res1$r, 1, tolerance = 1e-10)
})

test_that("es_from_spearman_rho produces Fisher's z of the converted r", {
  res <- es_from_spearman_rho(spearman_r = 0.40, n_sample = 80)
  # z = atanh(r_p) with r_p = 2*sin(pi/6 * 0.40); frozen external values:
  expect_equal(res$r, 0.4158233816355186, tolerance = 1e-12)
  expect_equal(res$z, 0.4426315850892963, tolerance = 1e-12)
  expect_equal(res$z, atanh(res$r), tolerance = 1e-10)
})

test_that("es_from_spearman_rho handles vectorized input", {
  res <- es_from_spearman_rho(
    spearman_r = c(0.3, 0.5, 0.7),
    n_sample = c(50, 100, 200)
  )
  expect_equal(nrow(res), 3)
  expect_equal(res$info_used, rep("spearman_r", 3))
  # vectorization consistency: each row equals its own scalar call
  row2 <- es_from_spearman_rho(spearman_r = 0.5, n_sample = 100)
  expect_equal(res$r[2], row2$r, tolerance = 1e-12)
  expect_equal(res$r_se[2], row2$r_se, tolerance = 1e-12)
})

test_that("Spearman integrates into convert_df hierarchy for measure='r'", {
  dat <- data.frame(spearman_r = 0.55, n_sample = 100)
  mc <- convert_df(dat, measure = "r", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)
  expect_equal(s$info_used, "spearman_r")
  expect_equal(s$es, 0.5680306894078453, tolerance = 1e-10)
})

test_that("Spearman feeds into d/g hierarchy", {
  dat <- data.frame(spearman_r = 0.55, n_sample = 100)
  mc <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)
  expect_equal(s$info_used, "spearman_r")
  expect_false(is.na(s$es))
})
