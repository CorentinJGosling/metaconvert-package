# =============================================================================
# es_disattenuate — promoted from tests_save/checked/ (P26, audit digest [31]).
#
# Adapted to the phase-1 fixes in R/es_disattenuate.R:
#   - P13/P20a: rows with |r_corrected| > 0.999 now warn and NA ALL
#     derived-scale outputs (Fisher z, its SE/CI, and the r CI); only the
#     point estimate and its first-order SE are preserved. The archived
#     expectation that z stays finite on a clamped row is inverted.
# External anchors: psychmeta::correct_r and psych::correct.cor.
# =============================================================================

test_that("es_disattenuate corrects correlation correctly", {
  res <- es_disattenuate(r = 0.50, r_se = 0.05,
                          reliability_x = 0.85, reliability_y = 0.80,
                          n_sample = 100)
  A <- sqrt(0.85 * 0.80)
  expect_equal(res$r_corrected, 0.50 / A, tolerance = 1e-10)
  expect_equal(res$r_corrected_se, 0.05 / A, tolerance = 1e-10)
  expect_equal(res$attenuation_factor, A, tolerance = 1e-10)
})

test_that("es_disattenuate with perfect reliability returns original", {
  res <- es_disattenuate(r = 0.50, r_se = 0.05,
                          reliability_x = 1.0, reliability_y = 1.0,
                          n_sample = 100)
  expect_equal(res$r_corrected, 0.50, tolerance = 1e-10)
  expect_equal(res$r_corrected_se, 0.05, tolerance = 1e-10)
})

test_that("es_disattenuate provides Fisher's z transformation", {
  res <- es_disattenuate(r = 0.50, r_se = 0.05,
                          reliability_x = 0.85, reliability_y = 0.80,
                          n_sample = 100)
  expect_false(is.na(res$z_corrected))
  expect_false(is.na(res$z_corrected_se))
  # Fisher's z should be larger (in magnitude) than corrected r for |r| < 1
  expect_true(abs(res$z_corrected) >= abs(res$r_corrected))
})

test_that("es_disattenuate handles vectorized input", {
  res <- es_disattenuate(
    r = c(0.40, 0.60),
    r_se = c(0.06, 0.04),
    reliability_x = c(0.85, 0.90),
    reliability_y = c(0.80, 0.85),
    n_sample = c(80, 150)
  )
  expect_equal(nrow(res), 2)
  expect_equal(res$r_corrected[1], 0.40 / sqrt(0.85 * 0.80), tolerance = 1e-10)
  expect_equal(res$r_corrected[2], 0.60 / sqrt(0.90 * 0.85), tolerance = 1e-10)
})

test_that("es_disattenuate corrects for one measure only", {
  # Correct only for X (set Y reliability to 1)
  res_x <- es_disattenuate(r = 0.50, r_se = 0.05,
                            reliability_x = 0.85, reliability_y = 1.0,
                            n_sample = 100)
  expect_equal(res_x$r_corrected, 0.50 / sqrt(0.85), tolerance = 1e-10)
})

test_that("Fisher's z SE uses delta method from r_corrected_se", {
  res <- es_disattenuate(r = 0.50, r_se = 0.05,
                          reliability_x = 0.85, reliability_y = 0.80,
                          n_sample = 100)
  A <- sqrt(0.85 * 0.80)
  r_c <- 0.50 / A          # non-extreme (0.606), no bounding involved
  r_c_se <- 0.05 / A
  expected_z_se <- r_c_se / (1 - r_c^2)
  expect_equal(res$z_corrected_se, expected_z_se, tolerance = 1e-10)
})

test_that("Fisher's z SE differs with different attenuation factors", {
  # Same n_sample, different reliabilities -> different z SE
  res1 <- es_disattenuate(r = 0.30, r_se = 0.08,
                           reliability_x = 1.0, reliability_y = 1.0,
                           n_sample = 100)
  res2 <- es_disattenuate(r = 0.30, r_se = 0.08,
                           reliability_x = 0.50, reliability_y = 0.50,
                           n_sample = 100)
  # With the old buggy code, these would be equal (both 1/sqrt(97))
  expect_false(abs(res1$z_corrected_se - res2$z_corrected_se) < 1e-10)
  # Study with lower reliability has larger z SE (more uncertainty)
  expect_true(res2$z_corrected_se > res1$z_corrected_se)
})

test_that("Fisher's z SE with perfect reliability equals delta method on observed r", {
  r <- 0.50; r_se <- 0.05
  res <- es_disattenuate(r = r, r_se = r_se,
                          reliability_x = 1.0, reliability_y = 1.0,
                          n_sample = 100)
  # With A=1, r_corrected = r, r_corrected_se = r_se
  # delta method: z_se = r_se / (1 - r^2)
  expected <- r_se / (1 - r^2)
  expect_equal(res$z_corrected_se, expected, tolerance = 1e-10)
})

test_that("Fisher's z SE does not require n_sample", {
  # n_sample is missing -> z SE should still be computed (from r_corrected_se)
  res <- es_disattenuate(r = 0.50, r_se = 0.05,
                          reliability_x = 0.85, reliability_y = 0.80)
  expect_false(is.na(res$z_corrected_se))
  expect_true(is.finite(res$z_corrected_se))
})

test_that("es_disattenuate handles corrected r > 1: point kept, derived scales NA'd", {
  # High r with low reliabilities -> corrected r exceeds 1. Post-P13/P20a
  # behaviour: a warning fires, the point estimate and first-order SE are
  # preserved for inspection, and EVERY derived-scale output (Fisher z, its
  # SE/CI, and the corrected-r CI) is NA — a bounded atanh(0.999) would be a
  # pure artifact of the internal constant.
  expect_warning(
    res <- es_disattenuate(r = 0.80, r_se = 0.03,
                           reliability_x = 0.60, reliability_y = 0.60,
                           n_sample = 100),
    "Fisher's z outputs are set to NA"
  )
  expect_true(res$r_corrected > 1)
  expect_equal(res$r_corrected, 0.80 / sqrt(0.60 * 0.60), tolerance = 1e-10)
  expect_equal(res$r_corrected_se, 0.03 / sqrt(0.60 * 0.60), tolerance = 1e-10)
  expect_true(is.na(res$z_corrected))
  expect_true(is.na(res$z_corrected_se))
  expect_true(is.na(res$z_corrected_ci_lo))
  expect_true(is.na(res$z_corrected_ci_up))
  expect_true(is.na(res$r_corrected_ci_lo))
  expect_true(is.na(res$r_corrected_ci_up))
})

test_that("es_disattenuate derives r_se from n_sample when r_se is omitted", {
  r <- 0.50; n <- 100; rel_x <- 0.85; rel_y <- 0.80
  res <- es_disattenuate(r = r, reliability_x = rel_x, reliability_y = rel_y,
                         n_sample = n)
  A <- sqrt(rel_x * rel_y)
  r_se_derived <- sqrt((1 - r^2)^2 / (n - 1))       # large-sample Pearson SE
  expect_equal(res$r_corrected_se, r_se_derived / A, tolerance = 1e-10)
  expect_true(is.finite(res$z_corrected_se))
  # identical to passing that r_se explicitly
  res2 <- es_disattenuate(r = r, r_se = r_se_derived,
                          reliability_x = rel_x, reliability_y = rel_y, n_sample = n)
  expect_equal(res$r_corrected_se, res2$r_corrected_se, tolerance = 1e-10)
})


# ==============================================================================
# Cross-validation against psychmeta::correct_r()
# ==============================================================================

test_that("es_disattenuate matches psychmeta::correct_r (both measures)", {
  skip_if_not_installed("psychmeta")
  r_obs <- 0.50; rel_x <- 0.85; rel_y <- 0.80

  mc <- es_disattenuate(r = r_obs, r_se = 0.05,
                         reliability_x = rel_x, reliability_y = rel_y,
                         n_sample = 100)
  pm <- psychmeta::correct_r(rxyi = r_obs, rxx = rel_x, ryy = rel_y,
                              correction = "meas")
  expect_equal(mc$r_corrected, pm$correlations$rtp, tolerance = 1e-10)
})

test_that("es_disattenuate matches psychmeta::correct_r (one-sided X)", {
  skip_if_not_installed("psychmeta")
  r_obs <- 0.50; rel_x <- 0.85

  mc <- es_disattenuate(r = r_obs, r_se = 0.05,
                         reliability_x = rel_x, reliability_y = 1.0,
                         n_sample = 100)
  pm <- psychmeta::correct_r(rxyi = r_obs, rxx = rel_x, ryy = 1.0,
                              correction = "meas")
  expect_equal(mc$r_corrected, pm$correlations$rtp, tolerance = 1e-10)
})

test_that("es_disattenuate matches psychmeta::correct_r (vectorized)", {
  skip_if_not_installed("psychmeta")
  r_obs <- c(0.40, 0.60, 0.25)
  rel_x <- c(0.85, 0.90, 0.75)
  rel_y <- c(0.80, 0.85, 0.70)

  mc <- es_disattenuate(r = r_obs, r_se = rep(0.05, 3),
                         reliability_x = rel_x, reliability_y = rel_y,
                         n_sample = rep(100, 3))
  pm <- psychmeta::correct_r(rxyi = r_obs, rxx = rel_x, ryy = rel_y,
                              correction = "meas")
  expect_equal(mc$r_corrected, pm$correlations$rtp, tolerance = 1e-10)
})

test_that("es_disattenuate matches psychmeta::correct_r (small r, high reliability)", {
  skip_if_not_installed("psychmeta")
  r_obs <- 0.10; rel_x <- 0.95; rel_y <- 0.95

  mc <- es_disattenuate(r = r_obs, r_se = 0.08,
                         reliability_x = rel_x, reliability_y = rel_y,
                         n_sample = 50)
  pm <- psychmeta::correct_r(rxyi = r_obs, rxx = rel_x, ryy = rel_y,
                              correction = "meas")
  expect_equal(mc$r_corrected, pm$correlations$rtp, tolerance = 1e-10)
})

test_that("es_disattenuate matches psychmeta::correct_r (negative r)", {
  skip_if_not_installed("psychmeta")
  r_obs <- -0.35; rel_x <- 0.80; rel_y <- 0.85

  mc <- es_disattenuate(r = r_obs, r_se = 0.06,
                         reliability_x = rel_x, reliability_y = rel_y,
                         n_sample = 80)
  pm <- psychmeta::correct_r(rxyi = r_obs, rxx = rel_x, ryy = rel_y,
                              correction = "meas")
  expect_equal(mc$r_corrected, pm$correlations$rtp, tolerance = 1e-10)
})


# ==============================================================================
# Cross-validation against psych::correct.cor()
# ==============================================================================

test_that("es_disattenuate matches psych::correct.cor", {
  skip_if_not_installed("psych")
  r_obs <- 0.50; rel_x <- 0.85; rel_y <- 0.80

  mc <- es_disattenuate(r = r_obs, r_se = 0.05,
                         reliability_x = rel_x, reliability_y = rel_y,
                         n_sample = 100)
  cor_mat <- matrix(c(1, r_obs, r_obs, 1), nrow = 2)
  corrected_mat <- psych::correct.cor(x = cor_mat, y = c(rel_x, rel_y))
  # correct.cor puts disattenuated r in upper triangle
  expect_equal(mc$r_corrected, corrected_mat[1, 2], tolerance = 1e-10)
})
