# Cross-validation tests for paired statistics
# Verifies that different input formats yield identical results

library(testthat)
library(metaConvert)

# Test 1: Paired t-test matches means/SD calculations ====
test_that("Paired t-test matches means/SD - Cooper/morris_drm", {
  # Setup: Known paired data
  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 12
  sd_post <- 14
  n <- 30
  r <- 0.65

  # Calculate paired t statistic manually
  mean_diff <- mean_post - mean_pre
  sd_diff <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
  se_diff <- sd_diff / sqrt(n)
  paired_t <- mean_diff / se_diff

  # Method 1: From paired t-test
  result_t <- es_from_paired_t_single_group(
    paired_t_exp = paired_t,
    n_exp = n,
    r_pre_post_exp = r
  )

  # Method 2: From means/SD with cooper
  result_means <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r,
    pre_post_to_smd = "cooper"
  )

  # Should be identical (within numerical precision)
  expect_equal(result_t$d, result_means$d, tolerance = 1e-10)
  expect_equal(result_t$d_se^2, result_means$d_se^2, tolerance = 1e-10)
  expect_equal(result_t$g, result_means$g, tolerance = 1e-10)
})

test_that("Paired t-test formula derivation is correct", {
  # Verify: paired_t = mean_diff / (sd_diff / sqrt(n))
  # Therefore: d_rm = paired_t * sqrt(2*(1-r)/n)

  paired_t <- 2.5
  n <- 20
  r <- 0.7

  result <- es_from_paired_t_single_group(
    paired_t_exp = paired_t,
    n_exp = n,
    r_pre_post_exp = r
  )

  # Manual calculation
  d_expected <- paired_t * sqrt(2 * (1 - r) / n)

  expect_equal(result$d, d_expected, tolerance = 1e-10)
})

# Test 2: Mean change matches pre-post for all methods ====
# Note: bonett method is skipped because it requires baseline SD
# which is conceptually different from change SD

test_that("Mean change matches pre-post - cooper/morris_drm", {
  mean_change <- 8.5
  sd_change <- 12.3
  n <- 40
  r <- 0.6

  # Method 1: Mean change function
  result_change <- es_from_mean_change_sd_single_group(
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r
  )

  # Method 2: Pre-post with trick (mean_pre = 0, mean_exp = mean_change)
  result_prepost <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = 0,
    mean_exp = mean_change,
    mean_pre_sd_exp = 0,
    mean_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r,
    pre_post_to_smd = "cooper"
  )

  # Should be identical
  expect_equal(result_change$d, result_prepost$d, tolerance = 1e-10)
  expect_equal(result_change$d_se^2, result_prepost$d_se^2, tolerance = 1e-10)
})

test_that("Mean change equals direct d_z when r = 0.5", {
  # When r = 0.5, sqrt(2*(1-0.5)) = 1.0
  # So d_rm = mean_change/sd_change * 1.0 = d_z

  mean_change <- 7.2
  sd_change <- 10.5
  n <- 35
  r <- 0.5

  result <- es_from_mean_change_sd_single_group(
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r
  )

  # This should equal d_z (no correlation adjustment)
  d_z_expected <- mean_change / sd_change

  expect_equal(result$d, d_z_expected, tolerance = 1e-10)
})

# Test 3: Two-group paired matches difference of single-group ====
# These two tests pin the per-arm standardizer (pool_sd = FALSE, the DEFAULT):
# each arm is standardized by its OWN SD and the two within-group values are then
# subtracted, their independent sampling variances added. This is Morris (2008)
# d_ppc1 (from Becker 1988) -- the construction metafor users get from
# escalc(measure = "SMCR") per arm followed by `yi = yT - yC; vi = vT + vC`.
# The comparator below (single-group exp minus single-group nexp) IS that
# construction by definition, so the identity is EXACT under pool_sd = FALSE.
# It does NOT hold under the opt-in pooled standardizer (pool_sd = TRUE, Morris
# 2008 d_ppc2: one SD pooled across arms), which is covered by its own test below.
# pool_sd is passed explicitly here to document the path under test.
test_that("Two-group paired equals difference of two single-group - bonett (per-arm standardizer, pool_sd = FALSE)", {
  # Setup two groups
  mean_pre_t <- 40
  mean_post_t <- 52
  sd_pre_t <- 11
  sd_post_t <- 13
  n_t <- 25

  mean_pre_c <- 42
  mean_post_c <- 44
  sd_pre_c <- 10
  sd_post_c <- 11
  n_c <- 25

  r <- 0.68

  # Method 1: Two-group function on the per-arm path (the default; passed
  # explicitly so the test documents which standardizer it pins)
  result_two_group <- es_from_means_sd_pre_post(
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
    r_pre_post_exp = r,
    r_pre_post_nexp = r,
    pre_post_to_smd = "bonett",
    pool_sd = FALSE
  )

  # Method 2: Difference of two single-group calls
  result_t <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre_t,
    mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t,
    mean_sd_exp = sd_post_t,
    n_exp = n_t,
    r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  result_c <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre_c,
    mean_exp = mean_post_c,
    mean_pre_sd_exp = sd_pre_c,
    mean_sd_exp = sd_post_c,
    n_exp = n_c,
    r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  # On the per-arm path the two-group d IS d_treatment - d_control by
  # construction, so this is an exact identity (not an approximation).
  expected_d_diff <- result_t$d - result_c$d

  expect_equal(result_two_group$d, expected_d_diff, tolerance = 1e-10)
})

test_that("Two-group paired equals difference of single-group - morris_dav (per-arm standardizer, pool_sd = FALSE)", {
  mean_pre_t <- 35
  mean_post_t <- 48
  sd_pre_t <- 9
  sd_post_t <- 10
  n_t <- 30

  mean_pre_c <- 37
  mean_post_c <- 39
  sd_pre_c <- 8
  sd_post_c <- 9
  n_c <- 30

  r <- 0.72

  # Two-group function on the per-arm path -- the default (see note above)
  result_two_group <- es_from_means_sd_pre_post(
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
    r_pre_post_exp = r,
    r_pre_post_nexp = r,
    pre_post_to_smd = "morris_dav",
    pool_sd = FALSE
  )

  # Single-group calls
  result_t <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t,
    n_exp = n_t, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  result_c <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre_c, mean_exp = mean_post_c,
    mean_pre_sd_exp = sd_pre_c, mean_sd_exp = sd_post_c,
    n_exp = n_c, r_pre_post_exp = r,
    pre_post_to_smd = "morris_dav"
  )

  # Exact identity on the per-arm path (see note above)
  expected_d_diff <- result_t$d - result_c$d

  expect_equal(result_two_group$d, expected_d_diff, tolerance = 1e-10)
})

# Test 3b: pool_sd is FALSE by default; TRUE is an opt-in ====
# The two-group pre/post standardizer is a genuine analytic CHOICE, so the
# package does not make it silently:
#   pool_sd = FALSE (DEFAULT) -- standardize the change WITHIN each arm by that
#     arm's own SD, subtract, and ADD the two independent sampling variances.
#     Morris (2008) d_ppc1, from Becker (1988). This is what metafor users build
#     with escalc(measure = "SMCR") per arm + `yi = yT - yC; vi = vT + vC`
#     (metafor-project.org/doku.php/analyses:morris2008). It does NOT assume the
#     arms' true standardizing SDs are equal, which is why Viechtbauer calls it
#     "more broadly applicable".
#   pool_sd = TRUE (opt-in) -- divide the between-group mean change by a SINGLE
#     SD pooled across arms. Morris (2008) d_ppc2 (eq. 8-9): more efficient, but
#     it assumes equal true arm SDs.
# The two coincide under homoscedasticity and target different estimands
# otherwise. The first test pins the default; the second validates the pooled
# path (which is exercised only with pool_sd = TRUE passed EXPLICITLY).
test_that("pool_sd defaults to FALSE (per-arm d_ppc1, Becker 1988 / Morris d_ppc1)", {
  mean_pre_t <- 40
  mean_post_t <- 52
  sd_pre_t <- 11
  sd_post_t <- 13
  n_t <- 25

  mean_pre_c <- 42
  mean_post_c <- 44
  sd_pre_c <- 10
  sd_post_c <- 11
  n_c <- 25

  r <- 0.68

  call_pp <- function(...) {
    es_from_means_sd_pre_post(
      mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
      mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t, n_exp = n_t,
      mean_pre_nexp = mean_pre_c, mean_nexp = mean_post_c,
      mean_pre_sd_nexp = sd_pre_c, mean_sd_nexp = sd_post_c, n_nexp = n_c,
      r_pre_post_exp = r, r_pre_post_nexp = r, ...
    )
  }

  res_default <- call_pp(pre_post_to_smd = "bonett")
  res_perarm <- call_pp(pre_post_to_smd = "bonett", pool_sd = FALSE)
  res_pooled <- call_pp(pre_post_to_smd = "bonett", pool_sd = TRUE)

  # (1) Omitting pool_sd must give the PER-ARM result: the default is FALSE.
  expect_equal(res_default$d, res_perarm$d, tolerance = 1e-10)
  expect_equal(res_default$d_se, res_perarm$d_se, tolerance = 1e-10)

  # (2) The default is the difference of the two within-arm standardized mean
  #     changes, with their independent sampling variances ADDED (d_ppc1).
  res_t <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
    mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t,
    n_exp = n_t, r_pre_post_exp = r, pre_post_to_smd = "bonett"
  )
  res_c <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre_c, mean_exp = mean_post_c,
    mean_pre_sd_exp = sd_pre_c, mean_sd_exp = sd_post_c,
    n_exp = n_c, r_pre_post_exp = r, pre_post_to_smd = "bonett"
  )

  expect_equal(res_default$d, res_t$d - res_c$d, tolerance = 1e-10)
  expect_equal(res_default$d_se^2, res_t$d_se^2 + res_c$d_se^2, tolerance = 1e-10)

  # (3) The two constructions are genuinely different, so (1) really does
  #     discriminate between them (guards against the default silently
  #     switching to the pooled path).
  expect_false(isTRUE(all.equal(res_pooled$d, res_perarm$d, tolerance = 1e-6)))
})

test_that("pool_sd = TRUE (opt-in) pools the standardizer across arms - Morris d_ppc2", {
  mean_pre_t <- 40
  mean_post_t <- 52
  sd_pre_t <- 11
  sd_post_t <- 13
  n_t <- 25

  mean_pre_c <- 42
  mean_post_c <- 44
  sd_pre_c <- 10
  sd_post_c <- 11
  n_c <- 25

  r <- 0.68

  call_pp <- function(...) {
    es_from_means_sd_pre_post(
      mean_pre_exp = mean_pre_t, mean_exp = mean_post_t,
      mean_pre_sd_exp = sd_pre_t, mean_sd_exp = sd_post_t, n_exp = n_t,
      mean_pre_nexp = mean_pre_c, mean_nexp = mean_post_c,
      mean_pre_sd_nexp = sd_pre_c, mean_sd_nexp = sd_post_c, n_nexp = n_c,
      r_pre_post_exp = r, r_pre_post_nexp = r, ...
    )
  }

  # (1) bonett + pool_sd = TRUE is, by definition of the estimand, the
  #     between-group mean change over the pooled BASELINE SD.
  res_pooled <- call_pp(pre_post_to_smd = "bonett", pool_sd = TRUE)

  sd_pooled_baseline <- sqrt(
    ((n_t - 1) * sd_pre_t^2 + (n_c - 1) * sd_pre_c^2) / (n_t + n_c - 2)
  )
  d_estimand <- ((mean_post_t - mean_pre_t) - (mean_post_c - mean_pre_c)) /
    sd_pooled_baseline

  expect_equal(res_pooled$d, d_estimand, tolerance = 1e-10)

  # (2) External check: morris_dz + pool_sd = TRUE is a two-sample SMD computed
  #     on the change scores, so its point estimate must equal metafor's
  #     escalc(measure = "SMD") fed the change means/SDs. (Only the point
  #     estimate is pinned: metaConvert's variance follows the Hedges (1981)
  #     J^2 form, which differs from every escalc vtype by a small-sample
  #     df/J^2 convention.)
  skip_if_not_installed("metafor")

  res_dz <- call_pp(pre_post_to_smd = "morris_dz", pool_sd = TRUE)

  sd_change_t <- sqrt(sd_pre_t^2 + sd_post_t^2 - 2 * r * sd_pre_t * sd_post_t)
  sd_change_c <- sqrt(sd_pre_c^2 + sd_post_c^2 - 2 * r * sd_pre_c * sd_post_c)

  ref <- metafor::escalc(
    measure = "SMD",
    m1i = mean_post_t - mean_pre_t, m2i = mean_post_c - mean_pre_c,
    sd1i = sd_change_t, sd2i = sd_change_c,
    n1i = n_t, n2i = n_c
  )

  expect_equal(res_dz$g, as.numeric(ref$yi), tolerance = 1e-10)
})

# Test 4: Different input formats for same data ====
test_that("SD, SE, and CI inputs yield identical results - cooper", {
  mean_pre <- 45
  mean_post <- 55
  sd_pre <- 11
  sd_post <- 12
  n <- 50
  r <- 0.63

  # Method 1: From SD
  result_sd <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "cooper"
  )

  # Method 2: From SE
  se_pre <- sd_pre / sqrt(n)
  se_post <- sd_post / sqrt(n)

  result_se <- es_from_means_se_pre_post_single_group(
    mean_pre_exp = mean_pre, mean_exp = mean_post,
    mean_pre_se_exp = se_pre, mean_se_exp = se_post,
    n_exp = n, r_pre_post_exp = r,
    pre_post_to_smd = "cooper"
  )

  # Should be identical
  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d_se^2, result_se$d_se^2, tolerance = 1e-10)
})

test_that("Mean change from SD vs SE vs pval - all match", {
  mean_change <- 6.8
  sd_change <- 9.5
  n <- 45
  r <- 0.58

  # Method 1: From SD
  result_sd <- es_from_mean_change_sd_single_group(
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r
  )

  # Method 2: From SE
  se_change <- sd_change / sqrt(n)
  result_se <- es_from_mean_change_se_single_group(
    mean_change_exp = mean_change,
    mean_change_se_exp = se_change,
    n_exp = n,
    r_pre_post_exp = r
  )

  # Method 3: From p-value
  t_stat <- mean_change / se_change
  pval <- 2 * pt(abs(t_stat), df = n - 1, lower.tail = FALSE)

  result_pval <- es_from_mean_change_pval_single_group(
    mean_change_exp = mean_change,
    mean_change_pval_exp = pval,
    n_exp = n,
    r_pre_post_exp = r
  )

  # All should match
  expect_equal(result_sd$d, result_se$d, tolerance = 1e-10)
  expect_equal(result_sd$d, result_pval$d, tolerance = 1e-6)
  expect_equal(result_sd$d_se^2, result_se$d_se^2, tolerance = 1e-10)
})

# Test 5: Consistency check with external packages ====
# NOTE: This validation has been moved to test-external-validation-packages.R
# and test-external-validation-comprehensive.R which provide systematic
# comparisons against metafor, TOSTER, and other packages.
