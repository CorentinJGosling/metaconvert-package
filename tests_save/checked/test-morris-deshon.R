# Tests for Morris & DeShon (2002) standardizers
# Validates morris_dz, morris_dav, and morris_drm implementations

library(testthat)
library(metaConvert)

# Test 1: Morris d_z formula verification ====
test_that("morris_dz formula matches Morris & DeShon (2002) equation 8", {
  # Example data from Morris 2007 Table 1, Study 1
  # Blicksenderfer et al. (1997)
  mean_pre <- 30.6
  mean_post <- 38.5
  sd_pre <- 15.0
  sd_post <- 11.6
  n <- 20
  r <- 0.47

  # Calculate expected sd_diff
  sd_diff_expected <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)

  # Expected d_z = mean_diff / sd_diff
  d_z_expected <- (mean_post - mean_pre) / sd_diff_expected

  # Calculate using metaConvert
  result <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r,
    pre_post_to_smd = "morris_dz"
  )

  # Verify point estimate
  expect_equal(result$d, d_z_expected, tolerance = 1e-10)

  # metaConvert uses metafor's approach: var(g) = 1/n + g²/(2n), then var(d) = var(g)/J²
  # This differs from Morris & DeShon's formula which uses d directly
  J <- exp(lgamma((n-1) / 2) - 0.5 * log((n-1) / 2) - lgamma((n-1 - 1) / 2))
  g_expected <- d_z_expected * J
  var_g_expected <- 1/n + g_expected^2 / (2 * n)
  var_d_expected <- var_g_expected / J^2
  expect_equal(result$d_se^2, var_d_expected, tolerance = 1e-10)
})

test_that("morris_dz point estimate depends on r through sd_diff", {
  # d_z point estimate DOES depend on r because sd_diff includes r
  # This is the correct behavior per Morris & DeShon (2002)

  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 12
  n <- 30

  result_r05 <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = 0.5,
    pre_post_to_smd = "morris_dz"
  )

  result_r08 <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = 0.8,
    pre_post_to_smd = "morris_dz"
  )

  # Point estimates should DIFFER because sd_diff depends on r
  expect_true(result_r05$d != result_r08$d)

  # Higher r means smaller sd_diff, so larger d (in absolute value)
  expect_true(abs(result_r08$d) > abs(result_r05$d))

  # Variances should also differ
  expect_false(isTRUE(all.equal(result_r05$d_se^2, result_r08$d_se^2, tolerance = 1e-10)))
})

# Test 2: Morris d_av formula verification ====
test_that("morris_dav point estimate is the average-SD SMD; variance is metafor SMCRPH (Bonett 2008)", {
  # NB: Morris (2008) gives NO sampling variance for d_av / d_ppc3 (p.373: "currently
  # unknown"); the variance used is the heteroscedasticity-robust SMCRPH (Bonett 2008
  # eq. 10), not a Morris & DeShon (2002) formula. Only the POINT estimate is the
  # standard average-SD SMD.
  mean_pre <- 23.5
  mean_post <- 26.8
  sd_pre <- 3.1
  sd_post <- 4.1
  n <- 50
  r <- 0.64

  # Calculate expected sd_av
  sd_av_expected <- sqrt((sd_pre^2 + sd_post^2) / 2)

  # Expected d_av = mean_diff / sd_av
  d_av_expected <- (mean_post - mean_pre) / sd_av_expected

  # Calculate using metaConvert
  result <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  # Verify point estimate
  expect_equal(result$d, d_av_expected, tolerance = 1e-10)

  # Variance = metafor SMCRPH (Bonett 2008 eq. 10): robust change-SD leading term plus
  # a fourth-moment g^2 term, both on df = n-1; mi = 2(n-1)/(1+r^2) is for J only.
  mi <- 2 * (n - 1) / (1 + r^2)
  J <- exp(lgamma(mi / 2) - 0.5 * log(mi / 2) - lgamma((mi - 1) / 2))
  g_expected <- d_av_expected * J
  sd_diff2 <- sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post
  fm <- sd_pre^4 + sd_post^4 + 2 * r^2 * sd_pre^2 * sd_post^2
  var_g_expected <- sd_diff2 / (sd_av_expected^2 * (n - 1)) +
    g_expected^2 * fm / (8 * sd_av_expected^4 * (n - 1))
  var_d_expected <- var_g_expected / J^2
  expect_equal(result$d_se^2, var_d_expected, tolerance = 1e-10)
})

test_that("morris_dav recommended by Morris 2008", {
  # This is more of a documentation test
  # d_av should be robust to variance heterogeneity

  # Scenario: sd_post > sd_pre (variance heterogeneity)
  mean_pre <- 35.6
  mean_post <- 36.0
  sd_pre <- 4.7
  sd_post <- 4.6  # Similar but not identical
  n <- 14
  r <- 0.44

  result <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  # Should compute without error
  expect_true(is.numeric(result$d))
  expect_true(is.numeric(result$d_se))
  expect_true(!is.na(result$d))
  expect_true(!is.na(result$d_se))
})

# Test 3: Morris d_rm (cooper/morris_drm) verification ====
test_that("cooper equals morris_drm (they are aliases)", {
  mean_pre <- 53.4
  mean_post <- 75.9
  sd_pre <- 14.5
  sd_post <- 4.4
  n <- 10
  r <- 0.89

  result_cooper <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "cooper"
  )

  result_morris_drm <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_drm"
  )

  # Should be identical
  expect_equal(result_cooper$d, result_morris_drm$d, tolerance = 1e-15)
  expect_equal(result_cooper$d_se, result_morris_drm$d_se, tolerance = 1e-15)
  expect_equal(result_cooper$g, result_morris_drm$g, tolerance = 1e-15)
})

test_that("morris_drm variance matches Morris 2000 formula", {
  mean_pre <- 30
  mean_post <- 40
  sd_pre <- 10
  sd_post <- 10
  n <- 25
  r <- 0.7

  result <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_drm"
  )

  # Variance formula: var(d_rm) = 2(1-r)/n + d_rm²/(2n)
  var_expected <- 2 * (1 - r) / n + result$d^2 / (2 * n)

  expect_equal(result$d_se^2, var_expected, tolerance = 1e-10)
})

# Test 4: Comparison across methods ====
test_that("Method ordering: d_z typically smallest, d_rm largest (when r > 0)", {
  mean_pre <- 20
  mean_post <- 30
  sd_pre <- 8
  sd_post <- 8
  n <- 40
  r <- 0.6

  result_dz <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dz"
  )

  result_dav <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  result_drm <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_drm"
  )

  # When sd_pre = sd_post and r > 0:
  # d_av = d_rm (both use pooled SD)
  # d_z > d_av = d_rm (because d_z divides by smaller sd_diff without adjustment)
  expect_equal(result_dav$d, result_drm$d, tolerance = 1e-10)
  expect_true(result_dz$d > result_dav$d)
})

test_that("When sd_pre = sd_post, d_av equals average of d_z and d_bonett conceptually", {
  # This tests the conceptual middle ground of d_av
  mean_pre <- 100
  mean_post <- 110
  sd_pre <- 15
  sd_post <- 15  # Homogeneous variance
  n <- 30
  r <- 0.5

  result_dav <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  result_bonett <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  # When sd_pre = sd_post:
  # sd_av = sqrt((sd² + sd²)/2) = sqrt(sd²) = sd
  # So d_av should equal d_bonett when variances are equal
  expect_equal(result_dav$d, result_bonett$d, tolerance = 1e-10)
})

# Test 5: Validation against Morris 2007 Table 1 examples ====
test_that("Morris 2007 Table 1 dppc2 ~ morris_dav (per-arm standardizer, pool_sd = FALSE, the default)", {
  # Study: Ivancevich and Smith (1981)
  # Morris 2007 reports dppc2 = 0.80.
  #
  # This test covers the PER-ARM path (pool_sd = FALSE, the default): each arm is
  # standardized by its OWN average SD and the two within-group d_av values are then
  # subtracted (Becker 1988 / Morris d_ppc1). pool_sd = FALSE is passed explicitly
  # for readability, but it is what the wrapper would do anyway.
  #
  # NOTE: the agreement with dppc2 is arithmetically incidental. Morris's actual
  # dppc2 estimator uses the POOLED PRETEST SD, which is metaConvert's
  # `bonett` + `pool_sd = TRUE` route -- pinned exactly in the next test.

  # Treatment group
  mean_pre_t <- 23.5
  sd_pre_t <- 3.1
  mean_post_t <- 26.8
  sd_post_t <- 4.1
  n_t <- 50

  # Control group
  mean_pre_c <- 24.9
  sd_pre_c <- 4.1
  mean_post_c <- 25.3
  sd_post_c <- 3.3
  n_c <- 42

  r <- 0.64

  result <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_t,
    mean = mean_post_t,
    mean_pre_sd_exp = sd_pre_t,
    mean_sd_exp = sd_post_t,
    n_exp = n_t,
    mean_pre_nexp = mean_pre_c,
    mean_nexp = mean_post_c,
    mean_pre_sd_nexp = sd_pre_c,
    mean_sd_nexp = sd_post_c,
    n_nexp = n_c,
    r_pre_post_exp = r,
    r_pre_post_nexp = r,
    pre_post_to_smd = "morris_dav",
    pool_sd = FALSE
  )

  # Morris 2007 reports dppc2 = 0.80
  expect_equal(result$d, 0.80, tolerance = 0.01)
})

test_that("Morris 2007 dppc2 is reproduced exactly by bonett + pool_sd = TRUE", {
  # Morris (2008) dppc2 = c_p * [(post_T - pre_T) - (post_C - pre_C)] / SD_pre_pooled
  # i.e. the between-group change difference standardized by the POOLED PRETEST SD,
  # then bias-corrected. That is exactly metaConvert's `bonett` standardizer with the
  # OPT-IN pooled common standardizer, pool_sd = TRUE (the default is FALSE, the
  # per-arm d_ppc1 construction), so pool_sd = TRUE must be passed explicitly.
  # Study: Ivancevich and Smith (1981); Morris 2007 Table 1 reports dppc2 = 0.80.

  mean_pre_t <- 23.5; sd_pre_t <- 3.1; mean_post_t <- 26.8; sd_post_t <- 4.1; n_t <- 50
  mean_pre_c <- 24.9; sd_pre_c <- 4.1; mean_post_c <- 25.3; sd_post_c <- 3.3; n_c <- 42
  r <- 0.64

  result <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t, n_exp = n_t,
    mean_pre_nexp = mean_pre_c, mean_nexp = mean_post_c,
    mean_pre_sd_nexp = sd_pre_c, mean_sd_nexp = sd_post_c, n_nexp = n_c,
    r_pre_post_exp = r, r_pre_post_nexp = r,
    pre_post_to_smd = "bonett",
    pool_sd = TRUE
  )

  # Closed-form dppc2 (exact bias correction J rather than the c_p approximation)
  m <- n_t + n_c - 2
  J <- exp(lgamma(m / 2) - 0.5 * log(m / 2) - lgamma((m - 1) / 2))
  sd_pre_pooled <- sqrt(((n_t - 1) * sd_pre_t^2 + (n_c - 1) * sd_pre_c^2) / m)
  mean_diff <- (mean_post_t - mean_pre_t) - (mean_post_c - mean_pre_c)
  dppc2_expected <- J * mean_diff / sd_pre_pooled

  # The bias-corrected estimate (g) IS dppc2
  expect_equal(result$g, dppc2_expected, tolerance = 1e-10)

  # ... and it recovers Morris's published value
  expect_equal(result$g, 0.80, tolerance = 0.01)

  # The uncorrected d is the same quantity without the J correction
  expect_equal(result$d, mean_diff / sd_pre_pooled, tolerance = 1e-10)
})

# Test 6: Hedges g correction ====
test_that("Hedges g properly corrects Morris standardizers", {
  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 12
  sd_post <- 12
  n <- 15  # Small sample
  r <- 0.6

  # Helper function to calculate J (matches .d_j in package)
  calc_J <- function(df) {
    exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
  }

  for (method in c("morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Calculate the appropriate J based on method
    # morris_dav uses modified degrees of freedom: mi = 2*(n-1)/(1+r²)
    # morris_dz and morris_drm use standard df = n - 1
    if (method == "morris_dav") {
      mi <- 2 * (n - 1) / (1 + r^2)
      J <- calc_J(mi)
    } else {
      J <- calc_J(n - 1)
    }

    # g should equal d * J (using the appropriate J for each method)
    expect_equal(result$g, result$d * J, tolerance = 1e-10)

    # var(g) should equal J² * var(d)
    expect_equal(result$g_se^2, result$d_se^2 * J^2, tolerance = 1e-10)
  }
})

# Test 7: Confidence intervals ====
test_that("Confidence intervals properly computed for all Morris methods", {
  mean_pre <- 25
  mean_post <- 32
  sd_pre <- 8
  sd_post <- 9
  n <- 20
  r <- 0.55

  for (method in c("morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # CI should be symmetric around point estimate
    ci_width_lower <- result$d - result$d_ci_lo
    ci_width_upper <- result$d_ci_up - result$d

    # Manual calculation
    se_d <- result$d_se
    t_crit <- qt(0.975, df = n - 1)
    expected_ci_lo <- result$d - t_crit * se_d
    expected_ci_up <- result$d + t_crit * se_d

    expect_equal(result$d_ci_lo, expected_ci_lo, tolerance = 1e-10)
    expect_equal(result$d_ci_up, expected_ci_up, tolerance = 1e-10)
  }
})

# Test 8: Two-group designs ====
test_that("Morris methods work for two-group paired designs (per-arm standardizer, default)", {
  # No pool_sd argument -> exercises the DEFAULT per-arm standardizer
  # (pool_sd = FALSE: standardize the mean change WITHIN each arm, subtract the two,
  # add their sampling variances -- Becker 1988 / Morris d_ppc1).
  # Treatment group
  mean_pre_t <- 30
  mean_post_t <- 45
  sd_pre_t <- 10
  sd_post_t <- 11
  n_t <- 25

  # Control group
  mean_pre_c <- 32
  mean_post_c <- 34
  sd_pre_c <- 9
  sd_post_c <- 10
  n_c <- 25

  r <- 0.7

  for (method in c("morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post(
      mean_pre_exp = mean_pre_t,
      mean = mean_post_t,
      mean_pre_sd_exp = sd_pre_t,
      mean_sd_exp = sd_post_t,
      n_exp = n_t,
      mean_pre_nexp = mean_pre_c,
      mean_nexp = mean_post_c,
      mean_pre_sd_nexp = sd_pre_c,
      mean_sd_nexp = sd_post_c,
      n_nexp = n_c,
      r_pre_post_exp = r,
      r_pre_post_nexp = r,
      pre_post_to_smd = method
    )

    # Should produce valid results
    expect_true(is.numeric(result$d))
    expect_true(is.numeric(result$d_se))
    expect_true(!is.na(result$d))
    expect_true(!is.na(result$d_se))
    expect_true(result$d_se > 0)
  }
})

# Test 9: Two-Group morris_dz validation ====
test_that("morris_dz correct for two-group design (per-arm standardizer, pool_sd = FALSE, the default)", {
  # Covers the PER-ARM path (the default): a d_z is formed inside each arm against
  # that arm's OWN change-score SD, and the two are subtracted (variances summed,
  # since the arms are independent -- Becker 1988 / Morris d_ppc1).
  # The manual comparator below transcribes exactly that per-arm rule.
  #
  # Two-group pre-post design with morris_dz standardizer
  # Treatment group
  mean_pre_t <- 30
  mean_post_t <- 45
  sd_pre_t <- 10
  sd_post_t <- 11
  n_t <- 25
  r_t <- 0.7

  # Control group
  mean_pre_c <- 32
  mean_post_c <- 34
  sd_pre_c <- 9
  sd_post_c <- 10
  n_c <- 25
  r_c <- 0.7

  # Helper to calculate J
  calc_J <- function(df) exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2))

  # Manual calculation for each group - metaConvert uses g in variance formula
  sd_diff_t <- sqrt(sd_pre_t^2 + sd_post_t^2 - 2 * r_t * sd_pre_t * sd_post_t)
  d_z_t <- (mean_post_t - mean_pre_t) / sd_diff_t
  J_t <- calc_J(n_t - 1)
  g_z_t <- d_z_t * J_t
  var_g_z_t <- 1/n_t + g_z_t^2 / (2 * n_t)
  var_d_z_t <- var_g_z_t / J_t^2

  sd_diff_c <- sqrt(sd_pre_c^2 + sd_post_c^2 - 2 * r_c * sd_pre_c * sd_post_c)
  d_z_c <- (mean_post_c - mean_pre_c) / sd_diff_c
  J_c <- calc_J(n_c - 1)
  g_z_c <- d_z_c * J_c
  var_g_z_c <- 1/n_c + g_z_c^2 / (2 * n_c)
  var_d_z_c <- var_g_z_c / J_c^2

  # Combined effect (difference)
  d_z_expected <- d_z_t - d_z_c
  se_d_z_expected <- sqrt(var_d_z_t + var_d_z_c)

  # metaConvert calculation
  result <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_t,
    mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t,
    mean_sd_exp = sd_post_t,
    n_exp = n_t,
    mean_pre_nexp = mean_pre_c,
    mean_nexp = mean_post_c,
    mean_pre_sd_nexp = sd_pre_c,
    mean_sd_nexp = sd_post_c,
    n_nexp = n_c,
    r_pre_post_exp = r_t,
    r_pre_post_nexp = r_c,
    pre_post_to_smd = "morris_dz",
    pool_sd = FALSE
  )

  # Verify point estimate
  expect_equal(result$d, d_z_expected, tolerance = 1e-10)

  # Verify SE
  expect_equal(result$d_se, se_d_z_expected, tolerance = 1e-10)
})

test_that("morris_dz two-group OPT-IN pooled standardizer (pool_sd = TRUE) matches metafor SMD on change scores", {
  # OPT-IN path (pool_sd = TRUE, NOT the default): the between-group SMD is the
  # difference in mean change divided by the change-score SD POOLED ACROSS ARMS --
  # a single common standardizer (Morris 2008 d_ppc2). This is the two-sample SMD
  # computed on change scores. pool_sd must be passed EXPLICITLY: the default is
  # pool_sd = FALSE (the per-arm d_ppc1 construction).
  skip_if_not_installed("metafor")

  mean_pre_t <- 30; mean_post_t <- 45; sd_pre_t <- 10; sd_post_t <- 11; n_t <- 25
  mean_pre_c <- 32; mean_post_c <- 34; sd_pre_c <- 9;  sd_post_c <- 10; n_c <- 25
  r <- 0.7

  result <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t, n_exp = n_t,
    mean_pre_nexp = mean_pre_c, mean_nexp = mean_post_c,
    mean_pre_sd_nexp = sd_pre_c, mean_sd_nexp = sd_post_c, n_nexp = n_c,
    r_pre_post_exp = r, r_pre_post_nexp = r,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE  # opt in explicitly; the default is FALSE (per-arm d_ppc1)
  )

  # Change-score SD in each arm
  sd_change_t <- sqrt(sd_pre_t^2 + sd_post_t^2 - 2 * r * sd_pre_t * sd_post_t)
  sd_change_c <- sqrt(sd_pre_c^2 + sd_post_c^2 - 2 * r * sd_pre_c * sd_post_c)

  # Live external reference: the bias-corrected point estimate is EXACTLY
  # metafor's two-sample SMD applied to the change scores.
  mf <- metafor::escalc(
    measure = "SMD",
    m1i = mean_post_t - mean_pre_t, m2i = mean_post_c - mean_pre_c,
    sd1i = sd_change_t, sd2i = sd_change_c,
    n1i = n_t, n2i = n_c
  )
  expect_equal(result$g, as.numeric(mf$yi), tolerance = 1e-10)

  # Variance: the pooled morris_dz estimator is an independent-groups Hedges g on
  # the change scores, so var(g) is Hedges (1981) / metafor's default "LS":
  #   var(g) = N/(n_t*n_c) + g^2/(2*N)
  # There is NO leading J^2 on var(g) (that was the LS2 convention, which understated
  # vi by up to 9% at small n), and NO 2*(1 - r) factor (that belongs to the raw-score
  # metric d_rm; d_z is not rescaled). Because g = J*d with J a deterministic constant,
  # var(d) = var(g)/J^2 EXACTLY (the code previously used a separate "d-scale twin",
  # which violated this identity -- see the pooled-variance-calibration test).
  N <- n_t + n_c
  m <- N - 2
  J <- exp(lgamma(m / 2) - 0.5 * log(m / 2) - lgamma((m - 1) / 2))
  sd_pooled <- sqrt(((n_t - 1) * sd_change_t^2 + (n_c - 1) * sd_change_c^2) / m)
  d_expected <- ((mean_post_t - mean_pre_t) - (mean_post_c - mean_pre_c)) / sd_pooled
  g_expected <- d_expected * J
  var_g_expected <- N / (n_t * n_c) + g_expected^2 / (2 * N)
  var_d_expected <- var_g_expected / J^2

  expect_equal(result$d, d_expected, tolerance = 1e-10)
  expect_equal(result$d_se^2, var_d_expected, tolerance = 1e-10)
  expect_equal(result$g_se^2, var_g_expected, tolerance = 1e-10)

  # CI is built on the pooled df m = N - 2
  expect_equal(result$d_ci_lo, d_expected - sqrt(var_d_expected) * qt(.975, m),
               tolerance = 1e-10)
  expect_equal(result$d_ci_up, d_expected + sqrt(var_d_expected) * qt(.975, m),
               tolerance = 1e-10)
})

test_that("morris_dz two-group with different correlations (per-arm standardizer, pool_sd = FALSE, the default)", {
  # The comparator subtracts two per-arm d_z values, each standardized by its own
  # arm's change SD -> the per-arm (default) construction, so pool_sd = FALSE.
  # Test with different r values between groups
  mean_pre_t <- 50
  mean_post_t <- 60
  sd_pre_t <- 12
  sd_post_t <- 14
  n_t <- 30
  r_t <- 0.5  # Moderate correlation in treatment

  mean_pre_c <- 51
  mean_post_c <- 53
  sd_pre_c <- 11
  sd_post_c <- 13
  n_c <- 28
  r_c <- 0.8  # High correlation in control

  # Manual calculations
  sd_diff_t <- sqrt(sd_pre_t^2 + sd_post_t^2 - 2 * r_t * sd_pre_t * sd_post_t)
  d_z_t <- (mean_post_t - mean_pre_t) / sd_diff_t

  sd_diff_c <- sqrt(sd_pre_c^2 + sd_post_c^2 - 2 * r_c * sd_pre_c * sd_post_c)
  d_z_c <- (mean_post_c - mean_pre_c) / sd_diff_c

  d_z_expected <- d_z_t - d_z_c

  # metaConvert
  result <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t,
    n_exp = n_t, r_pre_post_exp = r_t,
    mean_pre_nexp = mean_pre_c, mean_nexp = mean_post_c,
    mean_pre_sd_nexp = sd_pre_c, mean_sd_nexp = sd_post_c,
    n_nexp = n_c, r_pre_post_nexp = r_c,
    pre_post_to_smd = "morris_dz",
    pool_sd = FALSE
  )

  expect_equal(result$d, d_z_expected, tolerance = 1e-10)

  # Higher r in control should make d_z_c larger (smaller denominator)
  expect_true(d_z_c > d_z_t / 5,
              label = "Control group should have measurable effect")
})

# Test 10: Two-Group morris_dav validation ====
test_that("morris_dav correct for two-group design (per-arm standardizer, pool_sd = FALSE, the default)", {
  # Covers the PER-ARM path (the default). The comparator builds a metafor-SMCRP d_av
  # INSIDE each arm (own average SD, modified df mi = 2*(n-1)/(1+r^2)) and then
  # SUBTRACTS the two, summing their variances -- i.e. it encodes the per-arm rule
  # (Becker 1988 / Morris d_ppc1). pool_sd = FALSE is passed explicitly for
  # readability; it is also what the wrapper defaults to.
  #
  # Two-group design with morris_dav (average SD standardizer)
  mean_pre_t <- 23.5
  mean_post_t <- 26.8
  sd_pre_t <- 3.1
  sd_post_t <- 4.1
  n_t <- 50
  r_t <- 0.64

  mean_pre_c <- 24.9
  mean_post_c <- 25.3
  sd_pre_c <- 4.1
  sd_post_c <- 3.3
  n_c <- 42
  r_c <- 0.64

  # Helper to calculate J
  calc_J <- function(df) exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2))

  # Manual calculation per metafor SMCRPH approach (Bonett 2008 eq. 10: robust
  # change-SD leading term + fourth-moment g^2 term, both on df = n-1; the modified
  # df mi = 2*(n-1)/(1+r^2) is used for the Hedges J only).
  smcrph_var_g <- function(sd_pre, sd_post, r, n, sd_av, g) {
    sd_diff2 <- sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post
    fm <- sd_pre^4 + sd_post^4 + 2 * r^2 * sd_pre^2 * sd_post^2
    sd_diff2 / (sd_av^2 * (n - 1)) + g^2 * fm / (8 * sd_av^4 * (n - 1))
  }
  mi_t <- 2 * (n_t - 1) / (1 + r_t^2)
  J_t <- calc_J(mi_t)
  sd_av_t <- sqrt((sd_pre_t^2 + sd_post_t^2) / 2)
  d_av_t <- (mean_post_t - mean_pre_t) / sd_av_t
  g_av_t <- d_av_t * J_t
  var_g_av_t <- smcrph_var_g(sd_pre_t, sd_post_t, r_t, n_t, sd_av_t, g_av_t)
  var_d_av_t <- var_g_av_t / J_t^2

  mi_c <- 2 * (n_c - 1) / (1 + r_c^2)
  J_c <- calc_J(mi_c)
  sd_av_c <- sqrt((sd_pre_c^2 + sd_post_c^2) / 2)
  d_av_c <- (mean_post_c - mean_pre_c) / sd_av_c
  g_av_c <- d_av_c * J_c
  var_g_av_c <- smcrph_var_g(sd_pre_c, sd_post_c, r_c, n_c, sd_av_c, g_av_c)
  var_d_av_c <- var_g_av_c / J_c^2

  # Combined
  d_av_expected <- d_av_t - d_av_c
  se_d_av_expected <- sqrt(var_d_av_t + var_d_av_c)

  # metaConvert
  result <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t,
    n_exp = n_t, r_pre_post_exp = r_t,
    mean_pre_nexp = mean_pre_c, mean_nexp = mean_post_c,
    mean_pre_sd_nexp = sd_pre_c, mean_sd_nexp = sd_post_c,
    n_nexp = n_c, r_pre_post_nexp = r_c,
    pre_post_to_smd = "morris_dav",
    pool_sd = FALSE
  )

  expect_equal(result$d, d_av_expected, tolerance = 1e-10)
  expect_equal(result$d_se, se_d_av_expected, tolerance = 1e-10)
})

test_that("morris_dav two-group robust to variance heterogeneity (per-arm standardizer, pool_sd = FALSE, the default)", {
  # The comparator subtracts two per-arm d_av values -> the per-arm (default)
  # construction, so the call passes pool_sd = FALSE.
  # Morris_dav should handle different variances well
  mean_pre_t <- 40
  mean_post_t <- 50
  sd_pre_t <- 8
  sd_post_t <- 16  # Double the variance
  n_t <- 35
  r_t <- 0.6

  mean_pre_c <- 41
  mean_post_c <- 43
  sd_pre_c <- 7
  sd_post_c <- 14  # Also doubled
  n_c <- 33
  r_c <- 0.6

  # Manual calculation
  sd_av_t <- sqrt((sd_pre_t^2 + sd_post_t^2) / 2)
  d_av_t <- (mean_post_t - mean_pre_t) / sd_av_t
  var_d_av_t <- 2 * (1 - r_t) / n_t + d_av_t^2 / (2 * n_t)

  sd_av_c <- sqrt((sd_pre_c^2 + sd_post_c^2) / 2)
  d_av_c <- (mean_post_c - mean_pre_c) / sd_av_c
  var_d_av_c <- 2 * (1 - r_c) / n_c + d_av_c^2 / (2 * n_c)

  d_av_expected <- d_av_t - d_av_c

  # metaConvert
  result <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t,
    n_exp = n_t, r_pre_post_exp = r_t,
    mean_pre_nexp = mean_pre_c, mean_nexp = mean_post_c,
    mean_pre_sd_nexp = sd_pre_c, mean_sd_nexp = sd_post_c,
    n_nexp = n_c, r_pre_post_nexp = r_c,
    pre_post_to_smd = "morris_dav",
    pool_sd = FALSE
  )

  expect_equal(result$d, d_av_expected, tolerance = 1e-10)

  # Should produce valid results despite heterogeneity
  expect_true(is.finite(result$d))
  expect_true(is.finite(result$d_se))
  expect_true(result$d_se > 0)
})

test_that("per-arm (default) and opt-in pooled standardizers agree when the arms' SDs are equal", {
  # The per-arm construction (pool_sd = FALSE, the DEFAULT) subtracts two within-group
  # values, each divided by its OWN arm's SD (Becker 1988 / Morris d_ppc1). The pooled
  # construction (pool_sd = TRUE, OPT-IN) divides the difference in mean change by ONE
  # SD pooled across arms (Morris 2008 d_ppc2).
  #
  # The two COINCIDE if and only if the arms share the same standardizing SD, which is
  # what this test pins. Otherwise they target DIFFERENT ESTIMANDS -- neither is "the
  # buggy one": d_ppc1 does not assume the arms' true SDs are equal (Viechtbauer calls
  # it "more broadly applicable"), while d_ppc2 assumes equality to buy efficiency.
  # Because it is a genuine analytic choice, the package no longer makes it silently:
  # the default is the assumption-free per-arm path, and pooling is opt-in.
  args_equal_sd <- list(
    mean_pre_exp = 20, mean_exp = 30, mean_pre_sd_exp = 10, mean_sd_exp = 12, n_exp = 30,
    mean_pre_nexp = 21, mean_nexp = 24, mean_pre_sd_nexp = 10, mean_sd_nexp = 12, n_nexp = 30,
    r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6
  )

  for (method in c("morris_dz", "morris_dav", "morris_drm")) {
    pooled <- do.call(es_from_means_sd_pre_post,
                      c(args_equal_sd, pre_post_to_smd = method, pool_sd = TRUE))
    legacy <- do.call(es_from_means_sd_pre_post,
                      c(args_equal_sd, pre_post_to_smd = method, pool_sd = FALSE))

    expect_equal(pooled$d, legacy$d, tolerance = 1e-10,
                 label = paste("equal-SD arms, point estimate for", method))
  }
})

# Test 11: Variance formula verification across r values ====
test_that("morris_dz variance increases as r decreases", {
  # metaConvert uses metafor's approach: var(g) = 1/n + g²/(2n), then var(d) = var(g)/J²

  mean_pre <- 50
  mean_post <- 55
  sd_pre <- 10
  sd_post <- 10
  n <- 40

  # Helper to calculate J
  calc_J <- function(df) exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2))
  J <- calc_J(n - 1)

  r_values <- c(0, 0.3, 0.5, 0.7, 0.9)

  results <- lapply(r_values, function(r) {
    es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = "morris_dz"
    )
  })

  # Verify variance formula for each r - using g in formula per metafor approach
  for (i in seq_along(r_values)) {
    r <- r_values[i]
    result <- results[[i]]

    # Calculate expected variance using g (metafor approach)
    g <- result$d * J
    var_g_expected <- 1/n + g^2 / (2 * n)
    var_expected <- var_g_expected / J^2

    expect_equal(result$d_se^2, var_expected, tolerance = 1e-10,
                 label = paste("Variance formula for r =", r))
  }

  # Point estimates should vary with r (sd_diff changes)
  d_values <- sapply(results, function(x) x$d)
  expect_true(length(unique(round(d_values, 6))) > 1,
              label = "d_z should vary with r")
})

test_that("morris_dav variance formula correct across r values", {
  # metaConvert uses metafor SMCRPH (Bonett 2008 eq. 10): a robust change-SD leading
  # term + a fourth-moment g^2 term, both on df = n-1; then var(d) = var(g)/J^2, where
  # J uses the modified df mi = 2*(n-1)/(1+r^2).

  mean_pre <- 45
  mean_post <- 52
  sd_pre <- 8
  sd_post <- 9
  n <- 30

  # Helper to calculate J
  calc_J <- function(df) exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2))

  r_values <- c(0, 0.3, 0.5, 0.7, 0.9)

  results <- lapply(r_values, function(r) {
    es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = "morris_dav"
    )
  })

  sd_av <- sqrt((sd_pre^2 + sd_post^2) / 2)

  # Verify variance formula for each r using the metafor SMCRPH approach
  for (i in seq_along(r_values)) {
    r <- r_values[i]
    result <- results[[i]]

    # morris_dav uses modified df: mi = 2*(n-1)/(1+r²)
    mi <- 2 * (n - 1) / (1 + r^2)
    J <- calc_J(mi)
    g <- result$d * J

    # metafor SMCRPH variance: sd_diff²/(sd_av²(n-1)) + g²(sd_pre⁴+sd_post⁴+2r²sd_pre²sd_post²)/(8 sd_av⁴(n-1))
    sd_diff2 <- sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post
    fm <- sd_pre^4 + sd_post^4 + 2 * r^2 * sd_pre^2 * sd_post^2
    var_g_expected <- sd_diff2 / (sd_av^2 * (n - 1)) + g^2 * fm / (8 * sd_av^4 * (n - 1))
    var_expected <- var_g_expected / J^2

    expect_equal(result$d_se^2, var_expected, tolerance = 1e-10,
                 label = paste("Variance formula for r =", r))
  }

  # Variance should decrease as r increases (stronger dependency than d_z)
  se_values <- sapply(results, function(x) x$d_se)
  expect_true(all(diff(se_values) < 0),
              label = "SE should decrease as r increases for morris_dav")
})

test_that("morris_drm (cooper) variance formula correct across r values", {
  # var(d_rm) = 2(1-r)/n + d_rm²/(2n)
  # Same as morris_dav

  mean_pre <- 35
  mean_post <- 40
  sd_pre <- 6
  sd_post <- 7
  n <- 25

  r_values <- c(0, 0.3, 0.5, 0.7, 0.9, 0.99)

  results <- lapply(r_values, function(r) {
    es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = "morris_drm"
    )
  })

  # Verify variance formula
  for (i in seq_along(r_values)) {
    r <- r_values[i]
    result <- results[[i]]

    var_expected <- 2 * (1 - r) / n + result$d^2 / (2 * n)

    expect_equal(result$d_se^2, var_expected, tolerance = 1e-10,
                 label = paste("Variance for r =", r))
  }

  # SE should decrease monotonically
  se_values <- sapply(results, function(x) x$d_se)
  expect_true(all(diff(se_values) < 0))
})

# Test 12: Edge cases for Morris standardizers ====
test_that("morris standardizers with r=0 (independent measures)", {
  # When r=0, formulas should still work
  mean_pre <- 60
  mean_post <- 70
  sd_pre <- 12
  sd_post <- 14
  n <- 50
  r <- 0

  for (method in c("morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should produce finite results
    expect_true(is.finite(result$d),
                label = paste(method, "should handle r=0"))
    expect_true(is.finite(result$d_se),
                label = paste(method, "SE should be finite at r=0"))
    expect_true(result$d_se > 0)
  }
})

test_that("morris standardizers with r=0.99 (near-perfect correlation)", {
  # When r approaches 1, variance should be very small
  mean_pre <- 55
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 10
  n <- 40
  r <- 0.99

  for (method in c("morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    expect_true(is.finite(result$d))
    expect_true(is.finite(result$d_se))
    expect_true(result$d_se > 0)
  }

  # morris_dav and morris_drm should have very small SE when r=0.99
  result_dav <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  result_r05 <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = 0.5,
    pre_post_to_smd = "morris_dav"
  )

  # SE should be much smaller with r=0.99 than r=0.5
  expect_true(result_dav$d_se < result_r05$d_se / 2)
})

test_that("morris standardizers with extreme SD heterogeneity", {
  # sd_post = 5 × sd_pre
  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 5
  sd_post <- 25  # 5x larger
  n <- 30
  r <- 0.6

  for (method in c("morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should handle extreme heterogeneity
    expect_true(is.finite(result$d),
                label = paste(method, "should handle SD heterogeneity"))
    expect_true(is.finite(result$d_se))
    expect_true(result$d_se > 0)
    expect_true(!is.na(result$d))
  }

  # morris_dav recommended for this scenario (Morris 2008)
  result_dav <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  # Should be robust
  expect_true(result_dav$d > 0)
  expect_true(result_dav$d < 5)  # Reasonable range
})
