# smd_var = "hedges_olkin" on every pre/post call below (2.1.0). These tests pin
# agreement with metafor's SMCC/SMCR/SMCRH/SMCRPH variances and the Bonett (2008) /
# Morris & DeShon (2002) formulas, which are the Hedges-Olkin ("LS") form. Since 2.1.0
# the pre/post routes honour smd_var like the two-group routes, and the package
# default "borenstein" applies J^2 to the d-scale variance instead; the metafor form is
# reached with smd_var = "hedges_olkin", which is what these pins now request. The
# exception is the cooper / morris_drm blocks, whose Morris & DeShon (2002) oracle is
# the d-scale variance 2(1-r)/n + d^2/(2n): that is the J^2-outside ("borenstein")
# convention, so those calls request the default explicitly.

# Comprehensive External Validation Tests
# Validates metaConvert pre-post functions against metafor and esc packages
# Tests both effect sizes and standard errors across all standardizers

library(testthat)
library(metaConvert)

# Check for required packages
has_metafor <- requireNamespace("metafor", quietly = TRUE)
has_esc <- requireNamespace("esc", quietly = TRUE)

# =============================================================================
# METAFOR VALIDATIONS
# =============================================================================

# Test 1: Bonett standardizer vs metafor SMCRH (two-group) ====
# NOTE: the metafor comparator below is built by computing each arm's SMCRH
# separately (each arm standardized by its OWN baseline SD) and SUBTRACTING the
# two, and by SUMMING the two per-arm variances. That construction IS the legacy
# per-arm standardizer, i.e. pool_sd = FALSE. The current default (pool_sd = TRUE)
# standardizes the mean-change difference by a SD pooled ACROSS arms (Morris 2008),
# which no per-arm-subtraction comparator can reproduce when the arms' SDs differ.
# So this test pins the legacy per-arm path and must pass pool_sd = FALSE.
test_that("bonett two-group matches metafor SMCRH for effect size and SE (legacy per-arm standardizer, pool_sd = FALSE)", {
  skip_if_not(has_metafor, "metafor package not available")

  # Two-group pre-post design
  mean_pre_exp <- 50
  mean_post_exp <- 60
  sd_pre_exp <- 10
  sd_post_exp <- 12
  n_exp <- 30
  r_exp <- 0.6

  mean_pre_nexp <- 51
  mean_post_nexp <- 53
  sd_pre_nexp <- 9
  sd_post_nexp <- 11
  n_nexp <- 28
  r_nexp <- 0.65

  # metaConvert calculation
  # pool_sd = FALSE: match the per-arm rule the metafor comparator replicates
  # (arm-specific baseline SD as standardizer, then subtract the two arms).
  result_mc <- es_from_means_sd_pre_post(smd_var = "hedges_olkin", 
    mean_pre_exp, mean_post_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r_exp, r_nexp,
    pre_post_to_smd = "bonett",
    pool_sd = FALSE
  )

  # metafor SMCRH computes (m1 - m2) / sd1, so to get (post - pre) / sd_pre:
  # We use m1=pre, m2=post, sd1=sd_pre, then NEGATE the result
  result_metafor_exp <- metafor::escalc(
    measure = "SMCRH",
    m1i = mean_pre_exp, m2i = mean_post_exp,
    sd1i = sd_pre_exp, sd2i = sd_post_exp,
    ni = n_exp, ri = r_exp,
    data = data.frame(
      m1i = mean_pre_exp, m2i = mean_post_exp,
      sd1i = sd_pre_exp, sd2i = sd_post_exp,
      ni = n_exp, ri = r_exp
    )
  )

  result_metafor_nexp <- metafor::escalc(
    measure = "SMCRH",
    m1i = mean_pre_nexp, m2i = mean_post_nexp,
    sd1i = sd_pre_nexp, sd2i = sd_post_nexp,
    ni = n_nexp, ri = r_nexp,
    data = data.frame(
      m1i = mean_pre_nexp, m2i = mean_post_nexp,
      sd1i = sd_pre_nexp, sd2i = sd_post_nexp,
      ni = n_nexp, ri = r_nexp
    )
  )

  # Negate metafor results to match (post - pre) convention, then compute difference
  # metaConvert: g = g_exp - g_nexp = (post_exp - pre_exp)/sd_pre - (post_nexp - pre_nexp)/sd_pre
  # metafor: yi = (pre - post)/sd_pre, so -yi = (post - pre)/sd_pre
  g_metafor <- -as.numeric(result_metafor_exp$yi) - (-as.numeric(result_metafor_nexp$yi))
  se_g_metafor <- sqrt(as.numeric(result_metafor_exp$vi) + as.numeric(result_metafor_nexp$vi))

  expect_equal(result_mc$g, g_metafor, tolerance = 1e-6,
               label = "Bonett g matches metafor SMCRH")
  expect_equal(result_mc$g_se, se_g_metafor, tolerance = 1e-6,
               label = "Bonett g SE matches metafor SMCRH")
})

test_that("bonett single-group matches metafor SMCRH for effect size and SE", {
  skip_if_not(has_metafor, "metafor package not available")

  mean_pre <- 45
  mean_post <- 52
  sd_pre <- 8
  sd_post <- 9
  n <- 25
  r <- 0.7

  # metaConvert
  result_mc <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  # metafor SMCRH computes (m1 - m2) / sd1
  # To get (post - pre) / sd_pre, use m1=pre, m2=post, sd1=sd_pre, then NEGATE
  result_metafor <- metafor::escalc(
    measure = "SMCRH",
    m1i = mean_pre, m2i = mean_post,
    sd1i = sd_pre, sd2i = sd_post,
    ni = n, ri = r,
    data = data.frame(
      m1i = mean_pre, m2i = mean_post,
      sd1i = sd_pre, sd2i = sd_post,
      ni = n, ri = r
    )
  )

  # metafor SMCRH returns Hedges' g (bias-corrected), negate to match post-pre
  expect_equal(result_mc$g, -as.numeric(result_metafor$yi), tolerance = 1e-6,
               label = "Single-group bonett g matches metafor SMCRH")
  expect_equal(result_mc$g_se, sqrt(as.numeric(result_metafor$vi)), tolerance = 1e-6,
               label = "Single-group bonett g SE matches metafor SMCRH")
})

# Test 2: Mean change consistency (manual calculation) ====
test_that("mean change produces correct MD and SE", {
  # Construct mean change data
  mean_change_exp <- 10
  sd_change_exp <- 5
  n_exp <- 30

  mean_change_nexp <- 2
  sd_change_nexp <- 4.5
  n_nexp <- 28

  # metaConvert calculation
  result_mc <- es_from_mean_change_sd(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change_exp,
    mean_change_sd_exp = sd_change_exp,
    n_exp = n_exp,
    mean_change_nexp = mean_change_nexp,
    mean_change_sd_nexp = sd_change_nexp,
    n_nexp = n_nexp
  )

  # Manual calculation
  md_expected <- mean_change_exp - mean_change_nexp
  se_expected <- sqrt(sd_change_exp^2 / n_exp + sd_change_nexp^2 / n_nexp)

  expect_equal(result_mc$md, md_expected, tolerance = 1e-6,
               label = "Mean change MD matches manual calculation")
  expect_equal(result_mc$md_se, se_expected, tolerance = 1e-6,
               label = "Mean change SE matches manual calculation")
})

# Test 3: Cooper/morris_drm validation (conceptual - metafor doesn't have direct equivalent) ====
test_that("cooper standardizer produces valid results matching manual calculation", {
  # Since metafor doesn't have a direct equivalent to cooper/morris_drm,
  # we validate against the published Morris & DeShon (2002) formula

  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 12
  n <- 30
  r <- 0.6

  # metaConvert
  result_mc <- es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "cooper"
  )

  # Manual calculation per Morris & DeShon (2002)
  sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
  d_rm_manual <- (mean_post - mean_pre) / sd_change * sqrt(2 * (1 - r))
  var_rm_manual <- 2 * (1 - r) / n + d_rm_manual^2 / (2 * n)
  se_rm_manual <- sqrt(var_rm_manual)

  expect_equal(result_mc$d, d_rm_manual, tolerance = 1e-10,
               label = "Cooper d matches Morris & DeShon (2002) formula")
  expect_equal(result_mc$d_se, se_rm_manual, tolerance = 1e-10,
               label = "Cooper SE matches Morris & DeShon (2002) formula")
})

# =============================================================================
# VARIANCE RELATIONSHIP TESTS
# =============================================================================

# Test 4: Hedges' g variance relates to Cohen's d variance by J² ====
test_that("Hedges g variance equals Cohen d variance times J²", {
  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 12
  n <- 30
  r <- 0.6

  result <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  # Calculate J
  df <- n - 1
  J <- exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))

  # Variance relationship
  var_g_expected <- result$d_se^2 * J^2

  expect_equal(result$g_se^2, var_g_expected, tolerance = 1e-10,
               label = "var(g) = J² × var(d)")
})

# Test 5: SE decreases as correlation increases (Cooper standardizer) ====
test_that("Cooper SE decreases monotonically as r increases", {
  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 10
  sd_post <- 10
  n <- 30

  r_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)

  results <- lapply(r_values, function(r) {
    es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
      mean_pre, mean_post, sd_pre, sd_post, n, r,
      pre_post_to_smd = "cooper"
    )
  })

  se_values <- sapply(results, function(x) x$d_se)

  # SE should decrease monotonically
  expect_true(all(diff(se_values) < 0),
              label = "Cooper SE decreases as r increases")
})

# Test 6: Two-group vs single-group consistency ====
test_that("Two-group reduces to single-group when control has zero change", {
  mean_pre_exp <- 50
  mean_post_exp <- 60
  sd_pre_exp <- 10
  sd_post_exp <- 12
  n_exp <- 30
  r_exp <- 0.6

  # Control with zero change
  mean_pre_nexp <- 55
  mean_post_nexp <- 55  # No change
  sd_pre_nexp <- 10
  sd_post_nexp <- 10
  n_nexp <- 30
  r_nexp <- 0.6

  # Two-group
  result_two <- es_from_means_sd_pre_post(smd_var = "hedges_olkin", 
    mean_pre_exp, mean_post_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r_exp, r_nexp,
    pre_post_to_smd = "bonett"
  )

  # Single-group (experimental only)
  result_single <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp, mean_post_exp, sd_pre_exp, sd_post_exp,
    n_exp, r_exp,
    pre_post_to_smd = "bonett"
  )

  # Should be approximately equal (accounting for control group variance)
  expect_equal(abs(result_two$d), abs(result_single$d), tolerance = 0.2,
               label = "Two-group with zero control change approximates single-group")
})

# Test 7: Mean change consistency ====
test_that("mean change from means/SD produces consistent results", {
  mean_pre <- 45
  mean_post <- 52
  sd_pre <- 8
  sd_post <- 9
  n <- 25
  r <- 0.7

  # Calculate mean change
  mean_change <- mean_post - mean_pre
  sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)

  # From means/SD (any standardizer should give same MD)
  result_means <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  # From mean change
  result_mc <- es_from_mean_change_sd_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n
  )

  expect_equal(result_means$mdw, result_mc$mdw, tolerance = 1e-6,
               label = "MDw from means/SD matches MDw from mean change")
  expect_equal(result_means$mdw_se, result_mc$mdw_se, tolerance = 1e-6,
               label = "MDw SE from means/SD matches SE from mean change")
})
