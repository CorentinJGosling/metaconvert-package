test_that("reliability_change_score computes correctly", {
  res <- reliability_change_score(reliability = 0.85, r_pre_post = 0.60)
  expected <- (0.85 - 0.60) / (1 - 0.60)
  expect_equal(res$rel_change, expected, tolerance = 1e-10)
})

test_that("reliability_change_score handles edge cases", {
  # Perfect reliability
  res <- reliability_change_score(reliability = 1.0, r_pre_post = 0.50)
  expect_equal(res$rel_change, 1.0, tolerance = 1e-10)

  # Reliability equals pre-post correlation (change score unreliable)
  res_zero <- reliability_change_score(reliability = 0.70, r_pre_post = 0.70)
  expect_equal(res_zero$rel_change, 0.0, tolerance = 1e-10)
})

test_that("reliability_change_score handles vectorized input", {
  res <- reliability_change_score(
    reliability = c(0.85, 0.90),
    r_pre_post = c(0.60, 0.70)
  )
  expect_equal(nrow(res), 2)
  expect_equal(res$rel_change[1], (0.85 - 0.60) / (1 - 0.60), tolerance = 1e-10)
  expect_equal(res$rel_change[2], (0.90 - 0.70) / (1 - 0.70), tolerance = 1e-10)
})

test_that("reliability_change_score: low reliability and high r_pre_post yields negative", {
  # When rel < r_pre_post, change score reliability is negative (unreliable)
  res <- reliability_change_score(reliability = 0.50, r_pre_post = 0.70)
  expect_true(res$rel_change < 0)
  expect_equal(res$rel_change, (0.50 - 0.70) / (1 - 0.70), tolerance = 1e-10)
})


test_that("compute_sem computes correctly", {
  res <- compute_sem(sd = 10, icc = 0.85, n_sample = 100)
  expected_sem <- 10 * sqrt(1 - 0.85)
  expect_equal(res$sem, expected_sem, tolerance = 1e-10)
  expect_false(is.na(res$sem_se))
  expect_true(res$sem_ci_lo < res$sem)
  expect_true(res$sem_ci_up > res$sem)
})

test_that("compute_sem: zero ICC gives SEM = SD", {
  res <- compute_sem(sd = 10, icc = 0, n_sample = 100)
  expect_equal(res$sem, 10, tolerance = 1e-10)
})

test_that("compute_sem: perfect ICC gives SEM = 0", {
  res <- compute_sem(sd = 10, icc = 1, n_sample = 100)
  expect_equal(res$sem, 0, tolerance = 1e-10)
})

test_that("compute_sem handles vectorized input", {
  res <- compute_sem(
    sd = c(10, 15),
    icc = c(0.85, 0.90),
    n_sample = c(100, 200)
  )
  expect_equal(nrow(res), 2)
  expect_equal(res$sem[1], 10 * sqrt(0.15), tolerance = 1e-10)
  expect_equal(res$sem[2], 15 * sqrt(0.10), tolerance = 1e-10)
})

test_that("compute_sem accepts custom icc_se", {
  res_default <- compute_sem(sd = 10, icc = 0.85, n_sample = 100)
  res_custom <- compute_sem(sd = 10, icc = 0.85, n_sample = 100, icc_se = 0.05)
  # SEM should be the same
  expect_equal(res_default$sem, res_custom$sem, tolerance = 1e-10)
  # But SEM SE should differ
  expect_false(res_default$sem_se == res_custom$sem_se)
})

test_that("compute_sem: delta-method SE validated against manual calculation", {
  sd_val <- 10; icc_val <- 0.85; n <- 100; icc_se_val <- 0.05
  res <- compute_sem(sd = sd_val, icc = icc_val, n_sample = n, icc_se = icc_se_val)

  # Manual delta-method variance: Var(SEM) = SD^2/(4(1-ICC)) * Var(ICC) + (1-ICC) * Var(SD)
  var_sd <- sd_val^2 / (2 * (n - 1))
  var_icc <- icc_se_val^2
  var_sem <- (sd_val^2 / (4 * (1 - icc_val))) * var_icc + (1 - icc_val) * var_sd
  expect_equal(res$sem_se, sqrt(var_sem), tolerance = 1e-10)
})

test_that("compute_sem: larger n gives smaller SE", {
  res_small <- compute_sem(sd = 10, icc = 0.85, n_sample = 30)
  res_large <- compute_sem(sd = 10, icc = 0.85, n_sample = 500)
  # Point estimate is the same regardless of n
  expect_equal(res_small$sem, res_large$sem, tolerance = 1e-10)
  # SE decreases with larger sample
  expect_true(res_large$sem_se < res_small$sem_se)
})

test_that("compute_sem default (no icc_se) uses the exact same-sample variance", {
  # SEM = SD*sqrt(1-ICC) = sqrt(MSE); Var(SEM) = SEM^2 / (2 (n-1)(k-1)), default k=2
  sd_val <- 10; icc_val <- 0.85; n <- 100
  res <- compute_sem(sd = sd_val, icc = icc_val, n_sample = n)
  sem <- sd_val * sqrt(1 - icc_val)
  expected_se <- sqrt(sem^2 / (2 * (n - 1) * (2 - 1)))
  expect_equal(res$sem_se, expected_se, tolerance = 1e-10)
})

test_that("compute_sem: n_measurements (k) affects the same-sample SE", {
  n <- 100; sd_val <- 10; icc_val <- 0.85
  sem <- sd_val * sqrt(1 - icc_val)
  res_k2 <- compute_sem(sd = sd_val, icc = icc_val, n_sample = n, n_measurements = 2)
  res_k3 <- compute_sem(sd = sd_val, icc = icc_val, n_sample = n, n_measurements = 3)
  expect_equal(res_k3$sem_se, sqrt(sem^2 / (2 * (n - 1) * (3 - 1))), tolerance = 1e-10)
  # more measurements -> more error df -> smaller SE
  expect_true(res_k3$sem_se < res_k2$sem_se)
})

test_that("compute_sem: exact same-sample SE is far below the old independence-delta", {
  # Regression guard for the Cov(SD,ICC) fix: the same-sample exact variance is
  # 2-6x smaller than the former independence-delta approximation.
  sd_val <- 10; icc_val <- 0.85; n <- 100
  res <- compute_sem(sd = sd_val, icc = icc_val, n_sample = n)   # exact chi-square
  var_sd <- sd_val^2 / (2 * (n - 1))
  var_icc_old <- (1 - icc_val^2)^2 / (n - 1)
  old_se <- sqrt((sd_val^2 / (4 * (1 - icc_val))) * var_icc_old + (1 - icc_val) * var_sd)
  expect_true(res$sem_se < 0.75 * old_se)
})


test_that("compute_sdc computes correctly", {
  sem_val <- 3.87
  res <- compute_sdc(sem = sem_val, sem_se = 0.5)
  multiplier <- qnorm(0.975) * sqrt(2)
  expected_sdc <- multiplier * sem_val
  expect_equal(res$sdc, expected_sdc, tolerance = 1e-10)
  expect_equal(res$sdc_se, multiplier * 0.5, tolerance = 1e-10)
})

test_that("compute_sdc chains from compute_sem", {
  sem_res <- compute_sem(sd = 10, icc = 0.85, n_sample = 100)
  sdc_res <- compute_sdc(sem = sem_res$sem, sem_se = sem_res$sem_se)
  multiplier <- qnorm(0.975) * sqrt(2)
  expect_equal(sdc_res$sdc, multiplier * sem_res$sem, tolerance = 1e-10)
  expect_true(sdc_res$sdc > sem_res$sem)  # SDC should be larger than SEM
})

test_that("compute_sdc: CI is correctly constructed", {
  res <- compute_sdc(sem = 5, sem_se = 0.3)
  expect_equal(res$sdc_ci_lo, max(0, res$sdc - qnorm(0.975) * res$sdc_se), tolerance = 1e-10)
  expect_equal(res$sdc_ci_up, res$sdc + qnorm(0.975) * res$sdc_se, tolerance = 1e-10)
  expect_true(res$sdc_ci_lo < res$sdc)
  expect_true(res$sdc_ci_up > res$sdc)
})

test_that("compute_sdc: SDC without sem_se gives NA for uncertainty", {
  res <- compute_sdc(sem = 5)
  multiplier <- qnorm(0.975) * sqrt(2)
  expect_equal(res$sdc, multiplier * 5, tolerance = 1e-10)
  expect_true(is.na(res$sdc_se))
})


# ==============================================================================
# Full SEM -> SDC pipeline cross-check
# ==============================================================================

test_that("full SEM-SDC pipeline: known values from de Vet et al. (2011)", {
  # de Vet et al. (2011) Table 6.2 example-like scenario:
  # SD = 10, ICC = 0.70 -> SEM = 10 * sqrt(0.30) = 5.477
  # SDC = 1.96 * sqrt(2) * 5.477 = 15.18 (rounded)
  sem_res <- compute_sem(sd = 10, icc = 0.70, n_sample = 50)
  sdc_res <- compute_sdc(sem = sem_res$sem, sem_se = sem_res$sem_se)

  expect_equal(sem_res$sem, 10 * sqrt(0.30), tolerance = 1e-10)
  expect_equal(sdc_res$sdc, qnorm(0.975) * sqrt(2) * 10 * sqrt(0.30), tolerance = 1e-10)
  # SDC should be approximately 2.77 * SEM
  expect_equal(sdc_res$sdc / sem_res$sem, qnorm(0.975) * sqrt(2), tolerance = 1e-10)
})
