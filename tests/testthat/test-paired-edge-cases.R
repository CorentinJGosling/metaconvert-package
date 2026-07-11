# Edge case tests for paired statistics formulas
# Tests boundary conditions and extreme values

library(testthat)
library(metaConvert)

# Test 1: r = 0 (independent measures) ====
test_that("r = 0 reduces to independent groups formula - morris_drm", {
  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 10
  n <- 30
  r <- 0  # Independent

  result <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_drm"
  )

  # When r = 0, variance should be: 2(1-0)/n + d²/(2n) = 2/n + d²/(2n)
  var_expected <- 2/n + result$d^2 / (2 * n)

  expect_equal(result$d_se^2, var_expected, tolerance = 1e-10)

  # d_rm when r=0 should be larger than when r>0 (less precision)
  result_r05 <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = 0.5,
    pre_post_to_smd = "morris_drm"
  )

  # Variance with r=0 should be larger
  expect_true(result$d_se > result_r05$d_se)
})

test_that("r = 0 makes d_z equal to independent groups effect size", {
  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 12
  sd_post <- 12
  n <- 25
  r <- 0

  result_dz <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dz"
  )

  # When r = 0 and sd_pre = sd_post:
  # sd_diff = sqrt(sd² + sd² - 0) = sqrt(2)*sd
  sd_diff_expected <- sqrt(sd_pre^2 + sd_post^2)
  d_expected <- (mean_post - mean_pre) / sd_diff_expected

  expect_equal(result_dz$d, d_expected, tolerance = 1e-10)
})

# Test 2: r = 0.99 (near-perfect correlation) ====
test_that("r = 0.99 gives smallest variance - morris_drm", {
  mean_pre <- 30
  mean_post <- 38
  sd_pre <- 8
  sd_post <- 8
  n <- 40

  result_r099 <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = 0.99,
    pre_post_to_smd = "morris_drm"
  )

  result_r05 <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = 0.5,
    pre_post_to_smd = "morris_drm"
  )

  # Higher correlation = smaller variance (var = 2(1-r)/n + d²/(2n))
  expect_true(result_r099$d_se < result_r05$d_se)

  # Variance should be very small when r is very high
  # With r=0.99, n_exp=40: var ≈ 2*0.01/40 = 0.0005, se ≈ 0.02 (ignoring d² term)
  # But d is also small (multiplied by sqrt(0.02)=0.14), so se remains small
  expect_true(result_r099$d_se < 0.15)
})

test_that("r = 0.99 makes d_rm much smaller than d_z", {
  # When r is very high, morris_drm multiplies by sqrt(2*(1-r)) which is very small
  # This makes d_rm < d_z (opposite of the old test expectation)
  mean_pre <- 20
  mean_post <- 25
  sd_pre <- 6
  sd_post <- 6
  n <- 35
  r <- 0.99

  result_drm <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_drm"
  )

  result_dz <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dz"
  )

  # d_rm should be much smaller: d_rm = d_z * sqrt(2*(1-0.99)) = d_z * 0.1414
  # So d_z should be about 7x larger than d_rm
  expect_true(result_dz$d > result_drm$d * 5)

  # Verify the relationship
  expected_ratio <- sqrt(2 * (1 - r))
  expect_equal(result_drm$d / result_dz$d, expected_ratio, tolerance = 1e-10)
})

# Test 3: Negative correlation ====
test_that("r = -0.5 handled correctly - all methods", {
  mean_pre <- 45
  mean_post <- 52
  sd_pre <- 9
  sd_post <- 10
  n <- 30
  r <- -0.5

  for (method in c("bonett", "morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should not error
    expect_true(is.numeric(result$d))
    expect_true(!is.na(result$d))

    # Variance should be positive
    expect_true(result$d_se > 0)
  }
})

test_that("Negative r increases variance compared to positive r", {
  mean_pre <- 35
  mean_post <- 42
  sd_pre <- 8
  sd_post <- 8
  n <- 25

  result_neg <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = -0.3,
    pre_post_to_smd = "morris_drm"
  )

  result_pos <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = 0.3,
    pre_post_to_smd = "morris_drm"
  )

  # Negative correlation increases variance
  # var = 2(1-r)/n, so 2(1-(-0.3))/n > 2(1-0.3)/n
  expect_true(result_neg$d_se > result_pos$d_se)
})

# Test 4: Homogeneous variance (sd_pre = sd_post) ====
test_that("sd_pre = sd_post creates algebraic simplifications", {
  mean_pre <- 50
  mean_post <- 60
  sd <- 10  # Same for both
  n <- 30
  r <- 0.6

  result_dav <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd, mean_sd_exp = sd,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  result_bonett <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd, mean_sd_exp = sd,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  # When sd_pre = sd_post:
  # sd_av = sqrt((sd² + sd²)/2) = sqrt(sd²) = sd
  # So d_av should equal d_bonett (both use sd as denominator)
  expect_equal(result_dav$d, result_bonett$d, tolerance = 1e-10)
})

test_that("sd_pre = sd_post: d_z and d_rm have known relationship", {
  mean_pre <- 40
  mean_post <- 52
  sd <- 12
  n <- 35
  r <- 0.7

  result_dz <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd, mean_sd_exp = sd,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dz"
  )

  result_drm <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd, mean_sd_exp = sd,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_drm"
  )

  # When sd_pre = sd_post:
  # sd_diff = sqrt(2*sd²*(1-r)) = sd*sqrt(2*(1-r))
  # d_z = mean_diff / (sd*sqrt(2*(1-r)))
  # d_rm = mean_diff / (sd*sqrt(2*(1-r))) * sqrt(2*(1-r)) = mean_diff / sd
  # Therefore: d_rm = d_z * sqrt(2*(1-r))

  expected_ratio <- sqrt(2 * (1 - r))
  actual_ratio <- result_drm$d / result_dz$d

  expect_equal(actual_ratio, expected_ratio, tolerance = 1e-10)
})

# Test 5: Extreme effect sizes ====
test_that("Very large effect sizes (d > 5) handled correctly", {
  mean_pre <- 10
  mean_post <- 100
  sd_pre <- 15
  sd_post <- 15
  n <- 50
  r <- 0.5

  for (method in c("bonett", "morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should produce valid results
    expect_true(result$d > 5)
    expect_true(is.finite(result$d))
    expect_true(is.finite(result$d_se))
    expect_true(result$d_se > 0)
  }
})

test_that("Negative effect sizes (d < -5) handled correctly", {
  mean_pre <- 100
  mean_post <- 10
  sd_pre <- 15
  sd_post <- 15
  n <- 50
  r <- 0.5

  for (method in c("bonett", "morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should produce valid negative results
    expect_true(result$d < -5)
    expect_true(is.finite(result$d))
    expect_true(is.finite(result$d_se))
    expect_true(result$d_se > 0)
  }
})

test_that("Very small effect sizes (|d| < 0.01) handled correctly", {
  mean_pre <- 50.00
  mean_post <- 50.05
  sd_pre <- 10
  sd_post <- 10
  n <- 100
  r <- 0.6

  for (method in c("bonett", "morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should produce valid small results
    expect_true(abs(result$d) < 0.1)
    expect_true(is.finite(result$d))
    expect_true(result$d_se > 0)
  }
})

# Test 6: Small sample sizes ====
test_that("Small samples (n = 5) handled correctly", {
  mean_pre <- 30
  mean_post <- 40
  sd_pre <- 8
  sd_post <- 9
  n <- 5
  r <- 0.6

  for (method in c("bonett", "morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should produce results with high variance
    expect_true(is.finite(result$d))
    expect_true(result$d_se > 0)

    # Variance should be relatively large for small n
    expect_true(result$d_se > 0.2)
  }
})

test_that("Hedges correction more important for small samples", {
  mean_pre <- 35
  mean_post <- 45
  sd_pre <- 10
  sd_post <- 10
  r <- 0.5

  # Small sample
  result_small <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = 8, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  # Large sample
  result_large <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = 100, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  # Calculate J values using EXACT formula (matches .d_j function)
  # For morris_dav, use modified df: mi = 2*(n-1)/(1+r²)
  calc_J <- function(df) {
    exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
  }

  mi_small <- 2 * (8 - 1) / (1 + r^2)
  mi_large <- 2 * (100 - 1) / (1 + r^2)
  J_small <- calc_J(mi_small)
  J_large <- calc_J(mi_large)

  # Ratio g/d should equal J
  ratio_small <- result_small$g / result_small$d
  ratio_large <- result_large$g / result_large$d

  expect_equal(ratio_small, J_small, tolerance = 1e-10)
  expect_equal(ratio_large, J_large, tolerance = 1e-10)

  # J should be smaller (more correction) for small n
  expect_true(J_small < J_large)
})

# Test 7: Extreme variance heterogeneity ====
test_that("sd_post >> sd_pre handled correctly", {
  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 5
  sd_post <- 25  # 5x larger
  n <- 30
  r <- 0.4

  for (method in c("bonett", "morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should handle extreme heterogeneity
    expect_true(is.finite(result$d))
    expect_true(is.finite(result$d_se))
    expect_true(result$d_se > 0)
  }
})

test_that("sd_pre >> sd_post handled correctly", {
  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 25
  sd_post <- 5  # 5x smaller
  n <- 30
  r <- 0.4

  for (method in c("bonett", "morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Should handle extreme heterogeneity
    expect_true(is.finite(result$d))
    expect_true(is.finite(result$d_se))
    expect_true(result$d_se > 0)
  }
})

# Test 8: Zero effect (mean_pre = mean_post) ====
test_that("Zero effect (mean_pre = mean_post) produces d = 0", {
  mean <- 50
  sd_pre <- 10
  sd_post <- 11
  n <- 40
  r <- 0.6

  for (method in c("bonett", "morris_dz", "morris_dav", "morris_drm")) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean,
      mean_exp = mean,
      mean_pre_sd_exp = sd_pre,
      mean_sd_exp = sd_post,
      n_exp = n,
      r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Effect size should be exactly zero
    expect_equal(result$d, 0, tolerance = 1e-15)
    expect_equal(result$g, 0, tolerance = 1e-15)

    # Variance should still be positive
    expect_true(result$d_se > 0)

    # CI should be symmetric around zero
    expect_equal(result$d_ci_lo, -result$d_ci_up, tolerance = 1e-10)
  }
})
