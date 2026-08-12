# Input Format Equivalence Tests
# Validates that SD, SE, and CI inputs produce identical results
# This ensures internal consistency across different input formats

library(testthat)
library(metaConvert)

# Helper function to convert SD to SE and CI
.convert_sd_to_se_ci <- function(mean, sd, n) {
  se <- sd / sqrt(n)
  t_crit <- qt(0.975, n - 1)
  ci_lo <- mean - t_crit * se
  ci_up <- mean + t_crit * se

  list(se = se, ci_lo = ci_lo, ci_up = ci_up)
}

# Test 1: Two-Group SD ↔ SE ↔ CI Equivalence (bonett) ====
test_that("SD/SE/CI inputs equivalent for bonett standardizer (two-group)", {
  # Base data
  mean_pre_exp <- 50
  mean_exp <- 60
  sd_pre_exp <- 10
  sd_post_exp <- 12
  n_exp <- 30

  mean_pre_nexp <- 51
  mean_post_nexp <- 53
  sd_pre_nexp <- 9
  sd_post_nexp <- 11
  n_nexp <- 28

  r_exp <- 0.6
  r_nexp <- 0.65

  # Convert SD to SE and CI for experimental group
  conv_pre_exp <- .convert_sd_to_se_ci(mean_pre_exp, sd_pre_exp, n_exp)
  conv_post_exp <- .convert_sd_to_se_ci(mean_exp, sd_post_exp, n_exp)

  # Convert SD to SE and CI for control group
  conv_pre_nexp <- .convert_sd_to_se_ci(mean_pre_nexp, sd_pre_nexp, n_nexp)
  conv_post_nexp <- .convert_sd_to_se_ci(mean_post_nexp, sd_post_nexp, n_nexp)

  # Calculate from SD
  result_sd <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_sd_exp = sd_pre_exp, mean_sd_exp = sd_post_exp,
    mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_post_nexp,
    mean_pre_sd_nexp = sd_pre_nexp, mean_sd_nexp = sd_post_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_exp, r_pre_post_nexp = r_nexp,
    pre_post_to_smd = "bonett"
  )

  # Calculate from SE
  result_se <- es_from_means_se_pre_post(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_se_exp = conv_pre_exp$se, mean_se_exp = conv_post_exp$se,
    mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_post_nexp,
    mean_pre_se_nexp = conv_pre_nexp$se, mean_se_nexp = conv_post_nexp$se,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_exp, r_pre_post_nexp = r_nexp,
    pre_post_to_smd = "bonett"
  )

  # Calculate from CI
  result_ci <- es_from_means_ci_pre_post(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_ci_lo_exp = conv_pre_exp$ci_lo, mean_pre_ci_up_exp = conv_pre_exp$ci_up,
    mean_ci_lo_exp = conv_post_exp$ci_lo, mean_ci_up_exp = conv_post_exp$ci_up,
    mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_post_nexp,
    mean_pre_ci_lo_nexp = conv_pre_nexp$ci_lo, mean_pre_ci_up_nexp = conv_pre_nexp$ci_up,
    mean_ci_lo_nexp = conv_post_nexp$ci_lo, mean_ci_up_nexp = conv_post_nexp$ci_up,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_exp, r_pre_post_nexp = r_nexp,
    pre_post_to_smd = "bonett"
  )

  # All should be identical
  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10,
               label = "SD and SE should give same d")
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8,
               label = "SD and CI should give same d")

  expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10,
               label = "SD and SE should give same SE")
  expect_equal(result_sd$d_se, result_ci$d_se, tolerance = 1e-8,
               label = "SD and CI should give same SE")

  # Test g as well
  expect_equal(result_sd$g, result_se$g, tolerance = 1e-10)
  expect_equal(result_sd$g, result_ci$g, tolerance = 1e-8)

  # Test MD
  expect_equal(result_sd$md, result_se$md, tolerance = 1e-10)
  expect_equal(result_sd$md, result_ci$md, tolerance = 1e-8)
})

# Test 2: Two-Group SD ↔ SE ↔ CI Equivalence (cooper) ====
test_that("SD/SE/CI inputs equivalent for cooper standardizer (two-group)", {
  mean_pre_exp <- 45
  mean_exp <- 55
  sd_pre_exp <- 8
  sd_post_exp <- 9
  n_exp <- 35

  mean_pre_nexp <- 46
  mean_post_nexp <- 47
  sd_pre_nexp <- 7
  sd_post_nexp <- 8
  n_nexp <- 33

  r_exp <- 0.7
  r_nexp <- 0.7

  # Convert
  conv_pre_exp <- .convert_sd_to_se_ci(mean_pre_exp, sd_pre_exp, n_exp)
  conv_post_exp <- .convert_sd_to_se_ci(mean_exp, sd_post_exp, n_exp)
  conv_pre_nexp <- .convert_sd_to_se_ci(mean_pre_nexp, sd_pre_nexp, n_nexp)
  conv_post_nexp <- .convert_sd_to_se_ci(mean_post_nexp, sd_post_nexp, n_nexp)

  # Calculate from each format
  result_sd <- es_from_means_sd_pre_post(
    mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r_exp, r_nexp,
    pre_post_to_smd = "cooper"
  )

  result_se <- es_from_means_se_pre_post(
    mean_pre_exp, mean_exp, conv_pre_exp$se, conv_post_exp$se,
    mean_pre_nexp, mean_post_nexp, conv_pre_nexp$se, conv_post_nexp$se,
    n_exp, n_nexp, r_exp, r_nexp,
    pre_post_to_smd = "cooper"
  )

  result_ci <- es_from_means_ci_pre_post(
    mean_pre_exp, mean_exp,
    conv_pre_exp$ci_lo, conv_pre_exp$ci_up,
    conv_post_exp$ci_lo, conv_post_exp$ci_up,
    mean_pre_nexp, mean_post_nexp,
    conv_pre_nexp$ci_lo, conv_pre_nexp$ci_up,
    conv_post_nexp$ci_lo, conv_post_nexp$ci_up,
    n_exp, n_nexp, r_exp, r_nexp,
    pre_post_to_smd = "cooper"
  )

  # Verify equivalence
  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)
  expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10)
  expect_equal(result_sd$g, result_se$g, tolerance = 1e-10)
  expect_equal(result_sd$md, result_se$md, tolerance = 1e-10)
})

# Test 3: Two-Group SD ↔ SE ↔ CI Equivalence (morris_dz) ====
test_that("SD/SE/CI inputs equivalent for morris_dz standardizer (two-group)", {
  mean_pre_exp <- 40
  mean_exp <- 50
  sd_pre_exp <- 10
  sd_post_exp <- 10
  n_exp <- 40

  mean_pre_nexp <- 41
  mean_post_nexp <- 42
  sd_pre_nexp <- 9
  sd_post_nexp <- 9
  n_nexp <- 38

  r <- 0.5

  conv_pre_exp <- .convert_sd_to_se_ci(mean_pre_exp, sd_pre_exp, n_exp)
  conv_post_exp <- .convert_sd_to_se_ci(mean_exp, sd_post_exp, n_exp)
  conv_pre_nexp <- .convert_sd_to_se_ci(mean_pre_nexp, sd_pre_nexp, n_nexp)
  conv_post_nexp <- .convert_sd_to_se_ci(mean_post_nexp, sd_post_nexp, n_nexp)

  result_sd <- es_from_means_sd_pre_post(
    mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_dz"
  )

  result_se <- es_from_means_se_pre_post(
    mean_pre_exp, mean_exp, conv_pre_exp$se, conv_post_exp$se,
    mean_pre_nexp, mean_post_nexp, conv_pre_nexp$se, conv_post_nexp$se,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_dz"
  )

  result_ci <- es_from_means_ci_pre_post(
    mean_pre_exp, mean_exp,
    conv_pre_exp$ci_lo, conv_pre_exp$ci_up,
    conv_post_exp$ci_lo, conv_post_exp$ci_up,
    mean_pre_nexp, mean_post_nexp,
    conv_pre_nexp$ci_lo, conv_pre_nexp$ci_up,
    conv_post_nexp$ci_lo, conv_post_nexp$ci_up,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_dz"
  )

  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)
  expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10)
})

# Test 4: Two-Group SD ↔ SE ↔ CI Equivalence (morris_dav) ====
test_that("SD/SE/CI inputs equivalent for morris_dav standardizer (two-group)", {
  mean_pre_exp <- 35
  mean_exp <- 42
  sd_pre_exp <- 6
  sd_post_exp <- 7
  n_exp <- 25

  mean_pre_nexp <- 36
  mean_post_nexp <- 37
  sd_pre_nexp <- 6
  sd_post_nexp <- 7
  n_nexp <- 25

  r <- 0.6

  conv_pre_exp <- .convert_sd_to_se_ci(mean_pre_exp, sd_pre_exp, n_exp)
  conv_post_exp <- .convert_sd_to_se_ci(mean_exp, sd_post_exp, n_exp)
  conv_pre_nexp <- .convert_sd_to_se_ci(mean_pre_nexp, sd_pre_nexp, n_nexp)
  conv_post_nexp <- .convert_sd_to_se_ci(mean_post_nexp, sd_post_nexp, n_nexp)

  result_sd <- es_from_means_sd_pre_post(
    mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_dav"
  )

  result_se <- es_from_means_se_pre_post(
    mean_pre_exp, mean_exp, conv_pre_exp$se, conv_post_exp$se,
    mean_pre_nexp, mean_post_nexp, conv_pre_nexp$se, conv_post_nexp$se,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_dav"
  )

  result_ci <- es_from_means_ci_pre_post(
    mean_pre_exp, mean_exp,
    conv_pre_exp$ci_lo, conv_pre_exp$ci_up,
    conv_post_exp$ci_lo, conv_post_exp$ci_up,
    mean_pre_nexp, mean_post_nexp,
    conv_pre_nexp$ci_lo, conv_pre_nexp$ci_up,
    conv_post_nexp$ci_lo, conv_post_nexp$ci_up,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_dav"
  )

  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)
  expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10)
})

# Test 5: Two-Group SD ↔ SE ↔ CI Equivalence (morris_drm) ====
test_that("SD/SE/CI inputs equivalent for morris_drm standardizer (two-group)", {
  mean_pre_exp <- 55
  mean_exp <- 65
  sd_pre_exp <- 12
  sd_post_exp <- 14
  n_exp <- 50

  mean_pre_nexp <- 56
  mean_post_nexp <- 58
  sd_pre_nexp <- 11
  sd_post_nexp <- 13
  n_nexp <- 48

  r <- 0.8

  conv_pre_exp <- .convert_sd_to_se_ci(mean_pre_exp, sd_pre_exp, n_exp)
  conv_post_exp <- .convert_sd_to_se_ci(mean_exp, sd_post_exp, n_exp)
  conv_pre_nexp <- .convert_sd_to_se_ci(mean_pre_nexp, sd_pre_nexp, n_nexp)
  conv_post_nexp <- .convert_sd_to_se_ci(mean_post_nexp, sd_post_nexp, n_nexp)

  result_sd <- es_from_means_sd_pre_post(
    mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_drm"
  )

  result_se <- es_from_means_se_pre_post(
    mean_pre_exp, mean_exp, conv_pre_exp$se, conv_post_exp$se,
    mean_pre_nexp, mean_post_nexp, conv_pre_nexp$se, conv_post_nexp$se,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_drm"
  )

  result_ci <- es_from_means_ci_pre_post(
    mean_pre_exp, mean_exp,
    conv_pre_exp$ci_lo, conv_pre_exp$ci_up,
    conv_post_exp$ci_lo, conv_post_exp$ci_up,
    mean_pre_nexp, mean_post_nexp,
    conv_pre_nexp$ci_lo, conv_pre_nexp$ci_up,
    conv_post_nexp$ci_lo, conv_post_nexp$ci_up,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "morris_drm"
  )

  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)
  expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10)

  # Verify morris_drm = cooper
  result_cooper <- es_from_means_sd_pre_post(
    mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r, r,
    pre_post_to_smd = "cooper"
  )

  expect_equal(result_sd$d, result_cooper$d, tolerance = 1e-15,
               label = "morris_drm should equal cooper")
})

# Test 6: Single-Group SD ↔ SE ↔ CI Equivalence (bonett) ====
test_that("SD/SE/CI inputs equivalent for bonett standardizer (single-group)", {
  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 12
  n <- 30
  r <- 0.6

  conv_pre <- .convert_sd_to_se_ci(mean_pre, sd_pre, n)
  conv_post <- .convert_sd_to_se_ci(mean_post, sd_post, n)

  result_sd <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  result_se <- es_from_means_se_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_se_exp = conv_pre$se, mean_se_exp = conv_post$se,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  result_ci <- es_from_means_ci_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_ci_lo_exp = conv_pre$ci_lo, mean_pre_ci_up_exp = conv_pre$ci_up,
    mean_ci_lo_exp = conv_post$ci_lo, mean_ci_up_exp = conv_post$ci_up,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)
  expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10)
  expect_equal(result_sd$d_se, result_ci$d_se, tolerance = 1e-8)
})

# Test 7: Single-Group SD ↔ SE ↔ CI for all standardizers ====
test_that("SD/SE/CI inputs equivalent for bonett and cooper (single-group)", {
  # Note: SE and CI functions only support bonett and cooper standardizers
  # Morris standardizers (morris_dz, morris_dav, morris_drm) only work with SD input
  mean_pre <- 45
  mean_post <- 52
  sd_pre <- 8
  sd_post <- 9
  n <- 25
  r <- 0.7

  conv_pre <- .convert_sd_to_se_ci(mean_pre, sd_pre, n)
  conv_post <- .convert_sd_to_se_ci(mean_post, sd_post, n)

  # Test bonett and cooper with all three input formats
  standardizers <- c("bonett", "cooper")

  for (std in standardizers) {
    result_sd <- es_from_means_sd_pre_post_single_group(
      mean_pre, mean_post, sd_pre, sd_post, n, r,
      pre_post_to_smd = std
    )

    result_se <- es_from_means_se_pre_post_single_group(
      mean_pre, mean_post, conv_pre$se, conv_post$se, n, r,
      pre_post_to_smd = std
    )

    result_ci <- es_from_means_ci_pre_post_single_group(
      mean_pre, mean_post,
      conv_pre$ci_lo, conv_pre$ci_up,
      conv_post$ci_lo, conv_post$ci_up,
      n, r, pre_post_to_smd = std
    )

    expect_equal(result_sd$d, result_se$d, tolerance = 1e-10,
                 label = paste("SD=SE for", std))
    expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8,
                 label = paste("SD=CI for", std))
    expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10,
                 label = paste("SE match for", std))
  }

  # Test Morris standardizers only with SD input (SE/CI not supported)
  morris_standardizers <- c("morris_dz", "morris_dav", "morris_drm")

  for (std in morris_standardizers) {
    result_sd <- es_from_means_sd_pre_post_single_group(
      mean_pre, mean_post, sd_pre, sd_post, n, r,
      pre_post_to_smd = std
    )

    # Should produce valid results
    expect_true(is.finite(result_sd$d),
                label = paste(std, "should work with SD input"))
    expect_true(is.finite(result_sd$d_se))
  }
})

# Test 8: Mean Change SD ↔ SE ↔ CI Equivalence (Two-Group) ====
test_that("Mean change SD/SE/CI inputs equivalent (two-group)", {
  # Construct mean change data
  mean_change_exp <- 10
  sd_change_exp <- 5
  n_exp <- 30
  r_exp <- 0.6

  mean_change_nexp <- 2
  sd_change_nexp <- 4
  n_nexp <- 28
  r_nexp <- 0.65

  # Convert to SE and CI
  conv_exp <- .convert_sd_to_se_ci(mean_change_exp, sd_change_exp, n_exp)
  conv_nexp <- .convert_sd_to_se_ci(mean_change_nexp, sd_change_nexp, n_nexp)

  # Calculate from SD
  result_sd <- es_from_mean_change_sd(
    mean_change_exp = mean_change_exp, mean_change_sd_exp = sd_change_exp,
    mean_change_nexp = mean_change_nexp, mean_change_sd_nexp = sd_change_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_exp, r_pre_post_nexp = r_nexp
  )

  # Calculate from SE
  result_se <- es_from_mean_change_se(
    mean_change_exp = mean_change_exp, mean_change_se_exp = conv_exp$se,
    mean_change_nexp = mean_change_nexp, mean_change_se_nexp = conv_nexp$se,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_exp, r_pre_post_nexp = r_nexp
  )

  # Calculate from CI
  result_ci <- es_from_mean_change_ci(
    mean_change_exp = mean_change_exp,
    mean_change_ci_lo_exp = conv_exp$ci_lo, mean_change_ci_up_exp = conv_exp$ci_up,
    mean_change_nexp = mean_change_nexp,
    mean_change_ci_lo_nexp = conv_nexp$ci_lo, mean_change_ci_up_nexp = conv_nexp$ci_up,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_exp, r_pre_post_nexp = r_nexp
  )

  # All should match
  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)
  expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10)
  expect_equal(result_sd$d_se, result_ci$d_se, tolerance = 1e-8)
})

# Test 9: Mean Change SD ↔ SE ↔ CI Equivalence (Single-Group) ====
test_that("Mean change SD/SE/CI inputs equivalent (single-group)", {
  mean_change <- 8
  sd_change <- 4
  n <- 25
  r <- 0.7

  conv <- .convert_sd_to_se_ci(mean_change, sd_change, n)

  result_sd <- es_from_mean_change_sd_single_group(
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r
  )

  result_se <- es_from_mean_change_se_single_group(
    mean_change_exp = mean_change,
    mean_change_se_exp = conv$se,
    n_exp = n,
    r_pre_post_exp = r
  )

  result_ci <- es_from_mean_change_ci_single_group(
    mean_change_exp = mean_change,
    mean_change_ci_lo_exp = conv$ci_lo,
    mean_change_ci_up_exp = conv$ci_up,
    n_exp = n,
    r_pre_post_exp = r
  )

  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)
  expect_equal(result_sd$d_se, result_se$d_se, tolerance = 1e-10)
})

# Test 10: Cross-validation with different r values ====
test_that("Input format equivalence holds across different r values", {
  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 10
  sd_post <- 10
  n <- 30

  r_values <- c(0, 0.3, 0.5, 0.7, 0.9)

  for (r in r_values) {
    conv_pre <- .convert_sd_to_se_ci(mean_pre, sd_pre, n)
    conv_post <- .convert_sd_to_se_ci(mean_post, sd_post, n)

    result_sd <- es_from_means_sd_pre_post_single_group(
      mean_pre, mean_post, sd_pre, sd_post, n, r,
      pre_post_to_smd = "cooper"
    )

    result_se <- es_from_means_se_pre_post_single_group(
      mean_pre, mean_post, conv_pre$se, conv_post$se, n, r,
      pre_post_to_smd = "cooper"
    )

    result_ci <- es_from_means_ci_pre_post_single_group(
      mean_pre, mean_post,
      conv_pre$ci_lo, conv_pre$ci_up,
      conv_post$ci_lo, conv_post$ci_up,
      n, r, pre_post_to_smd = "cooper"
    )

    expect_equal(result_sd$d, result_se$d, tolerance = 1e-10,
                 label = paste("SD=SE for r =", r))
    expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8,
                 label = paste("SD=CI for r =", r))
  }
})

# Test 11: Extreme sample sizes ====
test_that("Input format equivalence with small sample size", {
  mean_pre <- 30
  mean_post <- 35
  sd_pre <- 5
  sd_post <- 6
  n <- 10  # Small sample
  r <- 0.6

  conv_pre <- .convert_sd_to_se_ci(mean_pre, sd_pre, n)
  conv_post <- .convert_sd_to_se_ci(mean_post, sd_post, n)

  # Use bonett (supported by SE/CI functions)
  result_sd <- es_from_means_sd_pre_post_single_group(
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  result_se <- es_from_means_se_pre_post_single_group(
    mean_pre, mean_post, conv_pre$se, conv_post$se, n, r,
    pre_post_to_smd = "bonett"
  )

  result_ci <- es_from_means_ci_pre_post_single_group(
    mean_pre, mean_post,
    conv_pre$ci_lo, conv_pre$ci_up,
    conv_post$ci_lo, conv_post$ci_up,
    n, r, pre_post_to_smd = "bonett"
  )

  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)

  # CI conversion should handle t-distribution correctly
  expect_true(is.finite(result_ci$d))
  expect_true(is.finite(result_ci$d_se))
})

test_that("Input format equivalence with large sample size", {
  mean_pre <- 60
  mean_post <- 65
  sd_pre <- 12
  sd_post <- 14
  n <- 200  # Large sample
  r <- 0.75

  conv_pre <- .convert_sd_to_se_ci(mean_pre, sd_pre, n)
  conv_post <- .convert_sd_to_se_ci(mean_post, sd_post, n)

  result_sd <- es_from_means_sd_pre_post_single_group(
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  result_se <- es_from_means_se_pre_post_single_group(
    mean_pre, mean_post, conv_pre$se, conv_post$se, n, r,
    pre_post_to_smd = "bonett"
  )

  result_ci <- es_from_means_ci_pre_post_single_group(
    mean_pre, mean_post,
    conv_pre$ci_lo, conv_pre$ci_up,
    conv_post$ci_lo, conv_post$ci_up,
    n, r, pre_post_to_smd = "bonett"
  )

  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_ci$d, tolerance = 1e-8)
})
