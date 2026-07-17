# Externally anchored tests for es_from_cronbach_alpha() and es_from_icc().
# These replace the circular archived tests (tests_save/checked/test-alpha.R /
# test-icc.R, whose expectations recomputed the source formulas): the oracles
# here are metafor::escalc (ABT/ARAW) and an independent hand derivation of the
# ICC(3,1) variance from the F-statistic route -- NOT the package's own code.
# P16 edge-guard tests failed on the pre-fix tree (Inf/NaN emissions).

skip_if_not_installed("metafor")
library(metafor)

TOL <- 1e-12

# ---------------------------------------------------------------------------
# ALPHA -- Bonett transform vs metafor ABT (sign conventions differ: metaConvert
# stores +log(1 - alpha), metafor stores -log(1 - alpha); SEs must be identical)
# ---------------------------------------------------------------------------
test_that("alpha bonett: ES/SE match metafor ABT (opposite sign convention)", {
  grid <- list(c(0.80, 10, 100), c(0.60, 5, 40), c(0.95, 20, 250))
  for (g in grid) {
    a <- g[1]; k <- g[2]; n <- g[3]
    pkg <- es_from_cronbach_alpha(cronbach_alpha = a, n_sample = n, n_items = k)
    ora <- escalc(measure = "ABT", ai = a, mi = k, ni = n)
    expect_equal(pkg$alpha, -as.numeric(ora$yi), tolerance = TOL)
    expect_equal(pkg$alpha_se, sqrt(as.numeric(ora$vi)), tolerance = TOL)
  }
})

test_that("alpha raw: ES/SE match metafor ARAW", {
  grid <- list(c(0.80, 10, 100), c(0.60, 5, 40), c(0.95, 20, 250))
  for (g in grid) {
    a <- g[1]; k <- g[2]; n <- g[3]
    pkg <- es_from_cronbach_alpha(cronbach_alpha = a, n_sample = n, n_items = k,
                                  alpha_to_es = "raw")
    ora <- escalc(measure = "ARAW", ai = a, mi = k, ni = n)
    expect_equal(pkg$alpha, as.numeric(ora$yi), tolerance = TOL)
    expect_equal(pkg$alpha_se, sqrt(as.numeric(ora$vi)), tolerance = TOL)
  }
})

test_that("alpha bonett CI back-transform (with the bounds swap) brackets alpha", {
  pkg <- es_from_cronbach_alpha(cronbach_alpha = 0.8, n_sample = 100, n_items = 10)
  # ln(1 - a) is DECREASING in a: the transformed-scale UPPER bound maps to the
  # alpha-scale LOWER bound (the vignette recipe)
  a_lo <- 1 - exp(pkg$alpha_ci_up)
  a_up <- 1 - exp(pkg$alpha_ci_lo)
  expect_lt(a_lo, 0.8)
  expect_gt(a_up, 0.8)
  expect_lt(a_lo, a_up)
})

# ---------------------------------------------------------------------------
# ICC -- transformed SE vs an INDEPENDENT derivation from the F-statistic route:
# ICC3 = (F0 - 1)/(F0 + k - 1) with F0 on df (n - 1), (n - 1)(k - 1);
# Var(ln F0) = 2/(n-1) + 2/((n-1)(k-1)); ln(1 - ICC3) = ln k - ln(F0 + k - 1);
# delta method => Var(ln(1 - ICC)) = [(1 + (k-1) rho)/k]^2 * (2/(n-1) + 2/((n-1)(k-1)))
# ---------------------------------------------------------------------------
test_that("icc bonett SE matches the hand-derived F-route variance (consistency)", {
  grid <- list(c(0.8, 2, 50), c(0.6, 3, 30), c(0.9, 4, 100)) # rho, k, n
  for (g in grid) {
    rho <- g[1]; k <- g[2]; n <- g[3]
    pkg <- es_from_icc(icc = rho, n_sample = n, n_measurements = k,
                       icc_type = "consistency")
    var_lnF <- 2 / (n - 1) + 2 / ((n - 1) * (k - 1))
    se_hand <- sqrt(((1 + (k - 1) * rho) / k)^2 * var_lnF)
    expect_equal(pkg$icc, log(1 - rho), tolerance = TOL)
    expect_equal(pkg$icc_se, se_hand, tolerance = TOL)
  }
})

test_that("icc raw SE = (1 - rho) * transformed SE (delta-method identity)", {
  rho <- 0.8; k <- 2; n <- 50
  bon <- es_from_icc(icc = rho, n_sample = n, n_measurements = k,
                     icc_type = "consistency")
  raw <- es_from_icc(icc = rho, n_sample = n, n_measurements = k,
                     icc_type = "consistency", icc_to_es = "raw")
  expect_equal(raw$icc, rho, tolerance = TOL)
  expect_equal(raw$icc_se, (1 - rho) * bon$icc_se, tolerance = TOL)
})

# ---------------------------------------------------------------------------
# P16 -- edge cases degrade to NA (pre-fix: Inf/NaN/-Inf)
# ---------------------------------------------------------------------------
test_that("P16: alpha edge cases return NA, valid rows unaffected", {
  res <- es_from_cronbach_alpha(
    cronbach_alpha = c(1.0, 1.2, -0.5, 0.8, 0.8, 0.8),
    n_sample = c(100, 100, 100, 2, 1, 100),
    n_items = c(10, 10, 10, 10, 10, 1)
  )
  # alpha = 1 (bonett: log(0)), alpha > 1, n <= 2, k = 1 -> NA; alpha < 0 valid
  expect_true(is.na(res$alpha[1]) && is.na(res$alpha_se[1]))
  expect_true(is.na(res$alpha[2]) && is.na(res$alpha_se[2]))
  expect_equal(res$alpha[3], log(1.5))          # negative alpha is valid
  expect_true(all(is.na(res$alpha_se[4:6])))
  expect_false(any(is.infinite(res$alpha_se), na.rm = TRUE))
  expect_false(any(is.nan(res$alpha)))
  # raw scale: alpha = 1 is a representable boundary (es = 1, se = 0)
  raw1 <- es_from_cronbach_alpha(cronbach_alpha = 1, n_sample = 100, n_items = 10,
                                 alpha_to_es = "raw")
  expect_equal(raw1$alpha, 1)
  expect_equal(raw1$alpha_se, 0)
})

test_that("P16: icc edge cases return NA, valid rows unaffected", {
  res <- es_from_icc(
    icc = c(1, 1.3, -1.3, 0.8, 0.8, 0.8),
    n_sample = c(50, 50, 50, 1, 50, 50),
    n_measurements = c(2, 2, 2, 2, 1, 2)
  )
  expect_true(all(is.na(res$icc[1:3])))          # icc = 1, |icc| > 1
  expect_true(all(is.na(res$icc_se[4:5])))       # n <= 1, k < 2
  expect_false(any(is.infinite(res$icc), na.rm = TRUE))
  expect_false(any(is.infinite(res$icc_se), na.rm = TRUE))
  ok <- es_from_icc(icc = 0.8, n_sample = 50, n_measurements = 2)
  expect_equal(ok$icc, log(0.2), tolerance = TOL)
})

# ---------------------------------------------------------------------------
# Regression pins (verified-correct current values, 1e-12)
# ---------------------------------------------------------------------------
test_that("regression pins: one full output row per function per method", {
  a_bon <- es_from_cronbach_alpha(cronbach_alpha = 0.85, n_sample = 200, n_items = 10)
  expect_equal(a_bon$alpha, log(0.15), tolerance = TOL)
  expect_equal(a_bon$alpha_se, sqrt(20 / (9 * 198)), tolerance = TOL)
  expect_equal(a_bon$alpha_ci_lo, a_bon$alpha - qnorm(.975) * a_bon$alpha_se, tolerance = TOL)

  i_bon <- es_from_icc(icc = 0.80, n_sample = 50, n_measurements = 2)
  expect_equal(i_bon$icc, log(0.2), tolerance = TOL)
  expect_equal(i_bon$icc_se, sqrt(2 * 1.8^2 / (2 * 1 * 49)), tolerance = TOL)

  i_raw <- es_from_icc(icc = 0.80, n_sample = 50, n_measurements = 2, icc_to_es = "raw")
  expect_equal(i_raw$icc, 0.80, tolerance = TOL)
  expect_equal(i_raw$icc_se, 0.2 * sqrt(2 * 1.8^2 / (2 * 1 * 49)), tolerance = TOL)
})
