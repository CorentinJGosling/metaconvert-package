# =============================================================================
# Coverage for the es_disattenuate() audit fixes (AUDIT-FIX-TRACKER P12, P13,
# P20; digest [28], [29], [34]).
#
# Design rule for this file: every expected value is either a hand-derived
# formula written independently below (classical Spearman disattenuation
# r_c = r / sqrt(rel_x * rel_y); large-sample Pearson SE
# sqrt((1 - r^2)^2 / (n - 1)); delta-method Fisher z), or an external oracle
# (psychmeta::correct_r) -- never a re-typed copy of the package source.
# =============================================================================

# -----------------------------------------------------------------------------
# P12: scalar n_sample with a vector r must recycle -- pre-fix, rows 2..n of
# r_corrected_se were silently NA (n_sample[need_se] NA-padded the index).
# -----------------------------------------------------------------------------
test_that("P12: scalar n_sample recycles across a vector r (SE present for ALL rows)", {
  r <- c(0.3, 0.4, 0.5)
  res <- es_disattenuate(
    r = r, reliability_x = 0.8, reliability_y = 0.8, n_sample = 100
  )
  # Hand oracle: SE(r) = sqrt((1 - r^2)^2 / (n - 1)) (Cooper et al. 2019),
  # then r_c_se = SE(r) / sqrt(rel_x * rel_y) = SE(r) / 0.8.
  expected_se <- sqrt((1 - r^2)^2 / (100 - 1)) / 0.8
  expect_false(any(is.na(res$r_corrected_se)))
  expect_equal(res$r_corrected_se, expected_se, tolerance = 1e-10)
})

test_that("P12: scalar reliabilities recycle alongside a vector r", {
  r <- c(0.20, 0.45)
  res <- es_disattenuate(
    r = r, r_se = 0.05, reliability_x = 0.9, reliability_y = 0.7
  )
  A <- sqrt(0.9 * 0.7)
  expect_equal(res$r_corrected, r / A, tolerance = 1e-10)
  expect_equal(res$r_corrected_se, rep(0.05 / A, 2), tolerance = 1e-10)
  expect_equal(res$attenuation_factor, rep(A, 2), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# P13: rows with an impossible corrected r (|r_c| beyond the extreme threshold)
# must NA the z-scale outputs too, not just the r-scale CI. Pre-fix they
# emitted z_corrected = atanh(0.9999) = 4.9517 (a pure clamp artifact),
# z_corrected_se ~ 250, and finite z CIs.
# -----------------------------------------------------------------------------
test_that("P13: extreme corrected r yields NA for ALL z outputs (and the r CI)", {
  # r_c = 0.80 / sqrt(0.6 * 0.6) = 4/3 > 1: no valid Fisher z exists.
  expect_warning(
    res <- es_disattenuate(
      r = 0.80, r_se = 0.03, reliability_x = 0.6, reliability_y = 0.6
    )
  )
  # point estimate and its first-order SE are preserved (faithful arithmetic)
  expect_equal(res$r_corrected, 0.80 / 0.6, tolerance = 1e-10)
  expect_equal(res$r_corrected_se, 0.03 / 0.6, tolerance = 1e-10)
  # every derived-scale output is NA
  expect_true(is.na(res$z_corrected))
  expect_true(is.na(res$z_corrected_se))
  expect_true(is.na(res$z_corrected_ci_lo))
  expect_true(is.na(res$z_corrected_ci_up))
  expect_true(is.na(res$r_corrected_ci_lo))
  expect_true(is.na(res$r_corrected_ci_up))
})

test_that("P13: in a mixed vector only the extreme row is NA'd; clean row keeps its z", {
  expect_warning(
    res <- es_disattenuate(
      r = c(0.50, 0.80), r_se = c(0.05, 0.03),
      reliability_x = c(0.9, 0.6), reliability_y = c(0.9, 0.6)
    )
  )
  # clean row, hand-derived: r_c = 0.5/0.9, z = atanh(r_c), z_se = r_c_se/(1-r_c^2)
  rc <- 0.50 / 0.9
  rc_se <- 0.05 / 0.9
  expect_equal(res$z_corrected[1], atanh(rc), tolerance = 1e-10)
  expect_equal(res$z_corrected_se[1], rc_se / (1 - rc^2), tolerance = 1e-10)
  expect_false(is.na(res$r_corrected_ci_lo[1]))
  # extreme row fully NA on derived scales
  expect_true(is.na(res$z_corrected[2]))
  expect_true(is.na(res$z_corrected_se[2]))
  expect_true(is.na(res$z_corrected_ci_lo[2]))
  expect_true(is.na(res$z_corrected_ci_up[2]))
  expect_true(is.na(res$r_corrected_ci_lo[2]))
  expect_true(is.na(res$r_corrected_ci_up[2]))
})

test_that("P13: Inf corrected r (reliability 0) is also fully suppressed on derived scales", {
  w <- testthat::capture_warnings(
    res <- es_disattenuate(
      r = 0.50, r_se = 0.05, reliability_x = 0, reliability_y = 0.8
    )
  )
  expect_true(length(w) >= 1)
  expect_true(is.na(res$z_corrected))
  expect_true(is.na(res$z_corrected_se))
  expect_true(is.na(res$z_corrected_ci_lo))
  expect_true(is.na(res$z_corrected_ci_up))
})

# -----------------------------------------------------------------------------
# P20a: rows with 0.999 < |r_c| < 0.9999 used to get a warning claiming
# "values clamped to +/-0.9999" when nothing was clamped. The warning must now
# describe what actually happens (CI / Fisher z suppressed as unreliable), and
# the suppression must apply.
# -----------------------------------------------------------------------------
test_that("P20a: r_c = 0.9995 warns without claiming clamping; CI and z suppressed", {
  # A = sqrt(0.81 * 1) = 0.9; r_c = 0.89955 / 0.9 = 0.9995
  w <- testthat::capture_warnings(
    res <- es_disattenuate(
      r = 0.89955, r_se = 0.02, reliability_x = 0.81, reliability_y = 1
    )
  )
  expect_true(length(w) >= 1)
  expect_false(any(grepl("clamp", w, ignore.case = TRUE)))
  # point estimate preserved, derived scales suppressed
  expect_equal(res$r_corrected, 0.89955 / 0.9, tolerance = 1e-10)
  expect_true(is.na(res$r_corrected_ci_lo))
  expect_true(is.na(res$r_corrected_ci_up))
  expect_true(is.na(res$z_corrected))
  expect_true(is.na(res$z_corrected_se))
})

test_that("P20a: r_c just below the extreme threshold is untouched (no warning)", {
  # A = 0.9; r_c = 0.89 / 0.9 = 0.98888... < 0.999
  w <- testthat::capture_warnings(
    res <- es_disattenuate(
      r = 0.89, r_se = 0.02, reliability_x = 0.81, reliability_y = 1
    )
  )
  expect_length(w, 0)
  rc <- 0.89 / 0.9
  expect_equal(res$z_corrected, atanh(rc), tolerance = 1e-10)
  expect_false(is.na(res$r_corrected_ci_lo))
  expect_false(is.na(res$r_corrected_ci_up))
})

# -----------------------------------------------------------------------------
# P20b: an observed |r| > 1 must trigger an input-validity warning (mirroring
# the reliability checks). Pre-fix it passed through silently, with only the
# reliability-blaming clamp warning firing.
# -----------------------------------------------------------------------------
test_that("P20b: observed r outside [-1, 1] triggers an input-validity warning", {
  w <- testthat::capture_warnings(
    res <- es_disattenuate(
      r = 1.5, r_se = 0.05, reliability_x = 0.81, reliability_y = 1
    )
  )
  expect_true(any(grepl("observed r", w, ignore.case = TRUE)))
  # arithmetic is still faithful on the point estimate
  expect_equal(res$r_corrected, 1.5 / 0.9, tolerance = 1e-10)
})

test_that("P20b: a valid negative r does not trigger the input warning", {
  w <- testthat::capture_warnings(
    es_disattenuate(r = -0.35, r_se = 0.06,
                    reliability_x = 0.80, reliability_y = 0.85)
  )
  expect_length(w, 0)
})

# -----------------------------------------------------------------------------
# Regression: a clean case must keep all its (verified-correct) outputs.
# Hand-derived oracle written independently below; cross-checked against
# psychmeta::correct_r for the point estimate.
# -----------------------------------------------------------------------------
test_that("regression: clean case outputs match hand-derived formulas exactly", {
  res <- es_disattenuate(
    r = 0.5, r_se = 0.05, reliability_x = 0.8, reliability_y = 0.9
  )
  A <- sqrt(0.8 * 0.9)                 # attenuation factor
  rc <- 0.5 / A                        # Spearman (1904) disattenuation
  rc_se <- 0.05 / A                    # first-order approximation
  zc <- atanh(rc)                      # Fisher z of corrected r
  zc_se <- rc_se / (1 - rc^2)          # delta method, dz/dr = 1/(1 - r^2)
  zlo <- zc - qnorm(0.975) * zc_se
  zup <- zc + qnorm(0.975) * zc_se
  expect_equal(res$r_corrected, rc, tolerance = 1e-10)
  expect_equal(res$r_corrected_se, rc_se, tolerance = 1e-10)
  expect_equal(res$z_corrected, zc, tolerance = 1e-10)
  expect_equal(res$z_corrected_se, zc_se, tolerance = 1e-10)
  expect_equal(res$z_corrected_ci_lo, zlo, tolerance = 1e-10)
  expect_equal(res$z_corrected_ci_up, zup, tolerance = 1e-10)
  expect_equal(res$r_corrected_ci_lo, tanh(zlo), tolerance = 1e-10)
  expect_equal(res$r_corrected_ci_up, tanh(zup), tolerance = 1e-10)
  expect_equal(res$attenuation_factor, A, tolerance = 1e-10)
})

test_that("regression: clean case point estimate matches psychmeta::correct_r", {
  skip_if_not_installed("psychmeta")
  res <- es_disattenuate(
    r = 0.5, r_se = 0.05, reliability_x = 0.8, reliability_y = 0.9
  )
  pm <- psychmeta::correct_r(rxyi = 0.5, rxx = 0.8, ryy = 0.9,
                             correction = "meas")
  expect_equal(res$r_corrected, pm$correlations$rtp, tolerance = 1e-10)
})
