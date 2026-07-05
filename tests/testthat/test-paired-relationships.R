# Algebraic relationship tests for paired statistics
# Verifies mathematical relationships between different standardizers

library(testthat)
library(metaConvert)

# Test 1: Bonett vs Cooper algebraic relationship ====
test_that("Bonett vs Cooper when sd_pre = sd_post", {
  mean_pre <- 40
  mean_post <- 52
  sd <- 11  # Same for both
  n <- 35
  r <- 0.65

  result_bonett <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd, mean_sd_exp = sd,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  result_cooper <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd, mean_sd_exp = sd,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "cooper"
  )

  # When sd_pre = sd_post:
  # d_bonett = mean_diff / sd
  # d_cooper = mean_diff / (sd*sqrt(2*(1-r))) * sqrt(2*(1-r)) = mean_diff / sd
  # So they should be equal!
  expect_equal(result_bonett$d, result_cooper$d, tolerance = 1e-10)
})

test_that("Bonett > Cooper when sd_post > sd_pre", {
  mean_pre <- 30
  mean_post <- 40
  sd_pre <- 8
  sd_post <- 12
  n <- 30
  r <- 0.5

  result_bonett <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  result_cooper <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "cooper"
  )

  # Bonett uses smaller sd_pre in denominator, so d should be larger
  expect_true(result_bonett$d > result_cooper$d)
})

# Test 2: d_z vs d_av relationship ====
test_that("d_z vs d_av when sd_pre = sd_post", {
  mean_pre <- 45
  mean_post <- 58
  sd <- 10
  n <- 40
  r <- 0.6

  result_dz <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd, mean_sd_exp = sd,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dz"
  )

  result_dav <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd, mean_sd_exp = sd,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  # When sd_pre = sd_post = sd:
  # sd_diff = sd*sqrt(2*(1-r)) < sd (when r > 0)
  # d_z = mean_diff / sd_diff > mean_diff / sd = d_av
  # Therefore: d_z > d_av when r > 0

  expect_true(result_dz$d > result_dav$d)

  # Ratio: d_z / d_av = sd / sd_diff = 1 / sqrt(2*(1-r))
  ratio <- result_dz$d / result_dav$d
  expected_ratio <- 1 / sqrt(2 * (1 - r))

  expect_equal(ratio, expected_ratio, tolerance = 1e-10)
})

test_that("When sd_pre = sd_post, d_z largest, d_av = d_rm", {
  # When sd_pre = sd_post and r > 0:
  # d_av = d_rm (both use pooled SD)
  # d_z > d_av = d_rm (because d_z divides by smaller sd_diff without adjustment)
  mean_pre <- 35
  mean_post <- 47
  sd_pre <- 9
  sd_post <- 9
  n <- 30
  r <- 0.8  # High correlation

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
  # d_z > d_av = d_rm
  expect_equal(result_dav$d, result_drm$d, tolerance = 1e-10)
  expect_true(result_dz$d > result_dav$d)
})

# Test 3: All methods yield same sign ====
test_that("All methods agree on direction (positive effect)", {
  mean_pre <- 30
  mean_post <- 45
  sd_pre <- 10
  sd_post <- 11
  n <- 35
  r <- 0.55

  methods <- c("bonett", "morris_dz", "morris_dav", "morris_drm")
  results <- lapply(methods, function(method) {
    es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )
  })

  # All should be positive
  for (result in results) {
    expect_true(result$d > 0)
  }
})

test_that("All methods agree on direction (negative effect)", {
  mean_pre <- 60
  mean_post <- 48
  sd_pre <- 12
  sd_post <- 13
  n <- 40
  r <- 0.6

  methods <- c("bonett", "morris_dz", "morris_dav", "morris_drm")
  results <- lapply(methods, function(method) {
    es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )
  })

  # All should be negative
  for (result in results) {
    expect_true(result$d < 0)
  }
})

# Test 4: Variance formulas for d vs g ====
test_that("var(g) = J² × var(d) for all methods", {
  mean_pre <- 42
  mean_post <- 54
  sd_pre <- 11
  sd_post <- 12
  n <- 28
  r <- 0.62

  # Helper function to calculate J (matches .d_j in package)
  calc_J <- function(df) {
    exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
  }

  methods <- c("bonett", "morris_dz", "morris_dav", "morris_drm")
  for (method in methods) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Calculate the appropriate J based on method
    # morris_dav uses modified degrees of freedom: mi = 2*(n-1)/(1+r²)
    # Other methods use standard df = n - 1
    if (method == "morris_dav") {
      mi <- 2 * (n - 1) / (1 + r^2)
      J <- calc_J(mi)
    } else {
      J <- calc_J(n - 1)
    }

    # Check relationship: var(g) = J² × var(d)
    expect_equal(result$g_se^2, result$d_se^2 * J^2, tolerance = 1e-10,
                 label = paste("Method:", method))

    # Also verify g = d * J
    expect_equal(result$g, result$d * J, tolerance = 1e-10,
                 label = paste("Method:", method))
  }
})

# Test 5: Method ordering when r varies ====
test_that("As r increases, variance decreases for morris_drm", {
  mean_pre <- 38
  mean_post <- 47
  sd_pre <- 9
  sd_post <- 10
  n <- 35

  r_values <- c(0.2, 0.5, 0.8)
  variances <- sapply(r_values, function(r) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = "morris_drm"
    )
    result$d_se^2
  })

  # Variance should decrease as r increases
  expect_true(variances[1] > variances[2])
  expect_true(variances[2] > variances[3])
})

test_that("d_rm constant when sd_pre = sd_post (r cancels out)", {
  # When sd_pre = sd_post = sd:
  # sd_diff = sd * sqrt(2*(1-r))
  # d_rm = mean_diff / sd_diff * sqrt(2*(1-r)) = mean_diff / sd (constant!)
  mean_pre <- 50
  mean_post <- 60
  sd <- 12
  n <- 40

  r_values <- c(0.3, 0.6, 0.9)
  effect_sizes <- sapply(r_values, function(r) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd, mean_sd_exp = sd,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = "morris_drm"
    )
    result$d
  })

  # d_rm should be constant (r factors cancel)
  expect_equal(effect_sizes[1], effect_sizes[2], tolerance = 1e-10)
  expect_equal(effect_sizes[2], effect_sizes[3], tolerance = 1e-10)
})

test_that("d_z depends on r through sd_diff (Morris & DeShon 2002)", {
  # d_z = mean_diff / sd_diff where sd_diff = sqrt(sd_pre² + sd_post² - 2*r*sd_pre*sd_post)
  # Higher r means smaller sd_diff, thus LARGER d_z
  mean_pre <- 45
  mean_post <- 53
  sd_pre <- 10
  sd_post <- 11
  n <- 30

  r_values <- c(0.2, 0.5, 0.8)
  effect_sizes <- sapply(r_values, function(r) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = "morris_dz"
    )
    result$d
  })

  # Effect sizes should INCREASE with r (smaller sd_diff denominator)
  expect_true(effect_sizes[1] < effect_sizes[2])
  expect_true(effect_sizes[2] < effect_sizes[3])
})

# Test 6: Confidence interval relationships ====
test_that("Wider CIs for smaller n across all methods", {
  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 10
  sd_post <- 10
  r <- 0.6

  n_small <- 15
  n_large <- 100

  methods <- c("bonett", "morris_dz", "morris_dav", "morris_drm")
  for (method in methods) {
    result_small <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n_small, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    result_large <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n_large, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # CI width for small n
    ci_width_small <- result_small$d_ci_up - result_small$d_ci_lo

    # CI width for large n
    ci_width_large <- result_large$d_ci_up - result_large$d_ci_lo

    # Small n should have wider CIs
    expect_true(ci_width_small > ci_width_large,
                label = paste("Method:", method))
  }
})

# Test 7: Transformation consistency ====
test_that("d to g transformation consistent across methods", {
  mean_pre <- 32
  mean_post <- 41
  sd_pre <- 8
  sd_post <- 9
  n <- 25
  r <- 0.58

  # Helper function to calculate J
  calc_J <- function(df) {
    exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
  }

  methods <- c("bonett", "morris_dz", "morris_dav", "morris_drm")
  for (method in methods) {
    result <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = mean_pre, mean_exp = mean_post,
      mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
      n_exp = n, r_pre_post_exp = r,
      pre_post_to_smd = method
    )

    # Calculate the appropriate J based on method
    # morris_dav uses modified degrees of freedom: mi = 2*(n-1)/(1+r²)
    # Other methods use standard df = n - 1
    if (method == "morris_dav") {
      mi <- 2 * (n - 1) / (1 + r^2)
      J <- calc_J(mi)
    } else {
      J <- calc_J(n - 1)
    }

    # Verify transformation
    expect_equal(result$g, result$d * J, tolerance = 1e-10,
                 label = paste("Method:", method))
  }
})

# Test 8: Special case: r = 0.5 ====
test_that("When r = 0.5, sqrt(2*(1-r)) = 1.0 for morris_drm", {
  mean_pre <- 35
  mean_post <- 45
  sd_pre <- 10
  sd_post <- 10
  n <- 30
  r <- 0.5

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

  # When r = 0.5 and sd_pre = sd_post:
  # sqrt(2*(1-0.5)) = 1.0
  # So d_rm = d_z * 1.0 = d_z
  expect_equal(result_drm$d, result_dz$d, tolerance = 1e-10)
})
