# smd_var = "hedges_olkin" on every pre/post call below (2.1.0). These tests pin
# agreement with metafor's SMCC/SMCR/SMCRH/SMCRPH variances and the Bonett (2008) /
# Morris & DeShon (2002) formulas, which are the Hedges-Olkin ("LS") form. Since 2.1.0
# the pre/post routes honour smd_var like the two-group routes, and the package
# default "borenstein" applies J^2 to the d-scale variance instead; the metafor form is
# reached with smd_var = "hedges_olkin", which is what these pins now request.

### Tests for single-group within-group effect sizes (dw, gw, mdw) -----
# These tests validate the new single-group functions against metafor and TOSTER
#
# SCOPE NOTE (pool_sd / Morris 2008 default flip):
# Every function exercised here is a *_single_group() route. A single-group
# pre/post design has one arm, so there is no second arm's SD to pool with and
# nothing for a common standardizer to average over: the pool_sd argument does
# not exist on these functions, and the flip of its default (FALSE -> TRUE) in
# convert_df() and the seven two-group wrappers cannot reach this file. That is
# why -- unlike the two-group pre/post test files -- no call below needs an
# explicit pool_sd = FALSE. The standardizer here is always the single arm's own
# SD (pre SD for "bonett", change SD for "cooper"), which is unchanged.
#
# Likewise untouched by the paired-t variance harmonization: this file does not
# exercise es_from_paired_t_single_group() (the one single-group route whose
# morris_dz variance moved to the metafor SMCC convention). That function is
# covered by test-EXTERNAL-PAIRED-SINGLE-GROUP.R, test-paired-cross-validation.R
# and test-prepost-fixes.R.
#
# metafor is an Import; TOSTER is Suggests-only, so the external comparators are
# guarded per-test with skip_if_not_installed() rather than attached with
# library() at the top of the file (all comparator calls below are :: qualified).

library(testthat)
library(metaConvert)


### DW/GW - Pre-Post Means SD - Bonett method -----
test_that("DW - Means/SD - bonett - single group", {
  skip_if_not_installed("metafor")

  # Create test data
  set.seed(123)
  n <- 50
  mean_pre <- 100
  mean_post <- 95
  sd_pre <- 15
  sd_post <- 14
  r_pre_post <- 0.6

  # metafor calculation for within-group Hedges' g (SMCRH standardizes by pre SD)
  metafor_res <- metafor::escalc(
    m1i = mean_post,
    sd1i = sd_pre,
    m2i = mean_pre,
    sd2i = sd_post,
    ni = n,
    ri = r_pre_post,
    measure = "SMCRH"
  )

  # Convert metafor's g back to d (Bonett uses d = (mean_post - mean_pre) / sd_pre)
  J <- .d_j(n - 1)
  metafor_d <- as.numeric(metafor_res$yi) / J
  metafor_var_d <- metafor_res$vi / (J^2)
  metafor_se_d <- sqrt(metafor_var_d)

  # metaConvert calculation
  es <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Test d
  expect_equal(es$d, metafor_d, tolerance = 1e-10)
  expect_equal(es$d_se, metafor_se_d, tolerance = 1e-10)

  # Test g
  expect_equal(es$g, as.numeric(metafor_res$yi), tolerance = 1e-10)
  expect_equal(es$g_se, sqrt(metafor_res$vi), tolerance = 1e-10)

  # Test mdw
  expected_mdw <- mean_post - mean_pre
  expected_mdw_var <- (sd_pre^2 + sd_post^2 - 2 * r_pre_post * sd_pre * sd_post) / n
  expected_mdw_se <- sqrt(expected_mdw_var)

  expect_equal(es$mdw, expected_mdw, tolerance = 1e-10)
  expect_equal(es$mdw_se, expected_mdw_se, tolerance = 1e-10)

  # Test info_used
  expect_equal(es$info_used, "means_sd_pre_post_single_group")
})

test_that("GW - Means/SD - bonett - single group", {
  skip_if_not_installed("metafor")

  # Create test data
  set.seed(456)
  n <- 75
  mean_pre <- 50
  mean_post <- 48
  sd_pre <- 10
  sd_post <- 11
  r_pre_post <- 0.7

  # metafor calculation
  metafor_res <- metafor::escalc(
    m1i = mean_post,
    sd1i = sd_pre,
    m2i = mean_pre,
    sd2i = sd_post,
    ni = n,
    ri = r_pre_post,
    measure = "SMCRH"
  )

  # metaConvert calculation
  es <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Test g
  expect_equal(es$g, as.numeric(metafor_res$yi), tolerance = 1e-10)
  expect_equal(es$g_se, sqrt(metafor_res$vi), tolerance = 1e-10)
})

### DW/GW - Pre-Post Means SD - Cooper method -----
test_that("DW - Means/SD - cooper - single group", {
  skip_if_not_installed("TOSTER")

  # Create test data with simulated paired scores
  set.seed(789)
  n <- 60
  d_true <- 0.4

  scores_pre <- rnorm(n, 0, 1)
  scores_post <- scores_pre - rnorm(n, d_true, 1)

  mean_pre <- mean(scores_pre)
  mean_post <- mean(scores_post)
  sd_pre <- sd(scores_pre)
  sd_post <- sd(scores_post)
  r_pre_post <- cor(scores_pre, scores_post)

  # TOSTER calculation (Cooper method with bias_correction = FALSE)
  toster_res <- TOSTER::smd_calc(
    formula = c(scores_post, scores_pre) ~ rep(c("A", "B"), each = n),
    bias_correction = FALSE,
    rm_correction = TRUE,
    paired = TRUE,
    smd_ci = "t"
  )

  # metaConvert calculation
  es <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "cooper",
    reverse_means_pre_post = FALSE
  )

  # Test d
  # The point estimate is the load-bearing external check here: metaConvert and
  # TOSTER agree on d_rm to machine precision (~5e-17 at this seed).
  expect_equal(es$d, toster_res$estimate, tolerance = 1e-6)
  # The SE is NOT an exact cross-package check. metaConvert's d_rm variance and
  # TOSTER's rm_correction SE use different conventions and genuinely disagree
  # (~10.8% at this seed: 0.1131 vs 0.1268), so the 5e-1 tolerance only asserts
  # the SEs are of the same order. Verified identical before and after the
  # pool_sd / paired-t variance work -- this gap is a long-standing convention
  # difference, not a regression, and the tolerance is left as found. Tightening
  # it would require first deciding which d_rm SE convention the package wants.
  expect_equal(es$d_se, toster_res$SE, tolerance = 5e-1)

  # Test info_used
  expect_equal(es$info_used, "means_sd_pre_post_single_group")
})

test_that("GW - Means/SD - cooper - single group", {
  skip_if_not_installed("TOSTER")

  # Create test data with simulated paired scores
  set.seed(321)
  n <- 80
  d_true <- 0.5

  scores_pre <- rnorm(n, 0, 1)
  scores_post <- scores_pre - rnorm(n, d_true, 1)

  mean_pre <- mean(scores_pre)
  mean_post <- mean(scores_post)
  sd_pre <- sd(scores_pre)
  sd_post <- sd(scores_post)
  r_pre_post <- cor(scores_pre, scores_post)

  # TOSTER calculation (Cooper method with bias_correction = TRUE for g)
  toster_res <- TOSTER::smd_calc(
    formula = c(scores_post, scores_pre) ~ rep(c("A", "B"), each = n),
    bias_correction = TRUE,
    rm_correction = TRUE,
    paired = TRUE,
    smd_ci = "t"
  )

  # metaConvert calculation
  es <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "cooper",
    reverse_means_pre_post = FALSE
  )

  # Test g
  # Point estimate agrees with TOSTER to machine precision (~1e-16 at this seed).
  expect_equal(es$g, toster_res$estimate, tolerance = 1e-6)
  # See the d_se note in the test above: the SE conventions differ (~2.3% here,
  # 0.0857 vs 0.0877), so this is an order-of-magnitude check only. Left as found.
  expect_equal(es$g_se, toster_res$SE, tolerance = 5e-1)
})

### Pre-Post Means SE -----
test_that("DW - Means/SE - bonett - single group", {
  # Create test data
  set.seed(111)
  n <- 50
  mean_pre <- 100
  mean_post <- 95
  sd_pre <- 15
  sd_post <- 14
  se_pre <- sd_pre / sqrt(n)
  se_post <- sd_post / sqrt(n)
  r_pre_post <- 0.6

  # Reference calculation using SD function
  es_sd <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Test SE function
  es_se <- es_from_means_se_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_se_exp = se_pre,
    mean_se_exp = se_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Should be identical
  expect_equal(es_se$d, es_sd$d, tolerance = 1e-10)
  expect_equal(es_se$g, es_sd$g, tolerance = 1e-10)
  expect_equal(es_se$mdw, es_sd$mdw, tolerance = 1e-10)
  expect_equal(es_se$info_used, "means_se_pre_post_single_group")
})

### Pre-Post Means CI -----
test_that("DW - Means/CI - bonett - single group", {
  # Create test data
  set.seed(222)
  n <- 50
  mean_pre <- 100
  mean_post <- 95
  sd_pre <- 15
  sd_post <- 14
  se_pre <- sd_pre / sqrt(n)
  se_post <- sd_post / sqrt(n)
  r_pre_post <- 0.6

  # Create CIs
  ci_lo_pre <- mean_pre - qt(0.975, n - 1) * se_pre
  ci_up_pre <- mean_pre + qt(0.975, n - 1) * se_pre
  ci_lo_post <- mean_post - qt(0.975, n - 1) * se_post
  ci_up_post <- mean_post + qt(0.975, n - 1) * se_post

  # Reference calculation using SD function
  es_sd <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Test CI function
  es_ci <- es_from_means_ci_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_ci_lo_exp = ci_lo_pre,
    mean_pre_ci_up_exp = ci_up_pre,
    mean_ci_lo_exp = ci_lo_post,
    mean_ci_up_exp = ci_up_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Should be very close (minor differences due to CI->SE conversion)
  expect_equal(es_ci$d, es_sd$d, tolerance = 1e-8)
  expect_equal(es_ci$g, es_sd$g, tolerance = 1e-8)
  expect_equal(es_ci$mdw, es_sd$mdw, tolerance = 1e-8)
  expect_equal(es_ci$info_used, "means_ci_pre_post_single_group")
})

### Mean Change SD -----
test_that("DW - Mean Change SD - single group", {
  # Create test data
  set.seed(333)
  n <- 50
  mean_pre <- 100
  mean_post <- 95
  sd_pre <- 15
  sd_post <- 14
  r_pre_post <- 0.6

  # Calculate mean change
  mean_change <- mean_post - mean_pre
  sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r_pre_post * sd_pre * sd_post)

  # Pre-post calculation (Cooper method, which is appropriate for change scores)
  es_prepost <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "cooper",
    reverse_means_pre_post = FALSE
  )

  # Mean change calculation
  es_change <- es_from_mean_change_sd_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = FALSE
  )

  # Should be identical
  expect_equal(es_change$d, es_prepost$d, tolerance = 1e-10)
  expect_equal(es_change$g, es_prepost$g, tolerance = 1e-10)
  expect_equal(es_change$mdw, es_prepost$mdw, tolerance = 1e-10)
  expect_equal(es_change$info_used, "mean_change_sd_single_group")
})

### Mean Change SE -----
test_that("DW - Mean Change SE - single group", {
  # Create test data
  set.seed(444)
  n <- 50
  mean_change <- 5
  sd_change <- 8
  se_change <- sd_change / sqrt(n)
  r_pre_post <- 0.6

  # Reference calculation using SD function
  es_sd <- es_from_mean_change_sd_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = FALSE
  )

  # Test SE function
  es_se <- es_from_mean_change_se_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_se_exp = se_change,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = FALSE
  )

  # Should be identical
  expect_equal(es_se$d, es_sd$d, tolerance = 1e-10)
  expect_equal(es_se$g, es_sd$g, tolerance = 1e-10)
  expect_equal(es_se$mdw, es_sd$mdw, tolerance = 1e-10)
  expect_equal(es_se$info_used, "mean_change_se_single_group")
})

### Mean Change CI -----
test_that("DW - Mean Change CI - single group", {
  # Create test data
  set.seed(555)
  n <- 50
  mean_change <- 5
  sd_change <- 8
  se_change <- sd_change / sqrt(n)
  r_pre_post <- 0.6

  # Create CI
  ci_lo <- mean_change - qt(0.975, n - 1) * se_change
  ci_up <- mean_change + qt(0.975, n - 1) * se_change

  # Reference calculation using SD function
  es_sd <- es_from_mean_change_sd_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = FALSE
  )

  # Test CI function
  es_ci <- es_from_mean_change_ci_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_ci_lo_exp = ci_lo,
    mean_change_ci_up_exp = ci_up,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = FALSE
  )

  # Should be very close
  expect_equal(es_ci$d, es_sd$d, tolerance = 1e-8)
  expect_equal(es_ci$g, es_sd$g, tolerance = 1e-8)
  expect_equal(es_ci$mdw, es_sd$mdw, tolerance = 1e-8)
  expect_equal(es_ci$info_used, "mean_change_ci_single_group")
})

### Mean Change p-value -----
test_that("DW - Mean Change p-value - single group", {
  # Create test data
  set.seed(666)
  n <- 50
  mean_change <- 5
  sd_change <- 8
  r_pre_post <- 0.6

  # Calculate p-value from t-test
  se_change <- sd_change / sqrt(n)
  t_stat <- mean_change / se_change
  pval <- 2 * pt(abs(t_stat), df = n - 1, lower.tail = FALSE)

  # Reference calculation using SD function
  es_sd <- es_from_mean_change_sd_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = FALSE
  )

  # Test p-value function
  es_pval <- es_from_mean_change_pval_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_pval_exp = pval,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = FALSE
  )

  # Should be identical
  expect_equal(es_pval$d, es_sd$d, tolerance = 1e-10)
  expect_equal(es_pval$g, es_sd$g, tolerance = 1e-10)
  expect_equal(es_pval$mdw, es_sd$mdw, tolerance = 1e-10)
  expect_equal(es_pval$info_used, "mean_change_pval_single_group")
})

### REVERSE - Pre-Post -----
test_that("REVERSE - Pre-Post - bonett - single group", {
  # Create test data
  set.seed(777)
  n <- 50
  mean_pre <- 100
  mean_post <- 95
  sd_pre <- 15
  sd_post <- 14
  r_pre_post <- 0.6

  # Normal direction
  es_normal <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Reversed direction
  es_reverse <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = TRUE
  )

  # Test d reversal
  expect_equal(es_normal$d, -es_reverse$d, tolerance = 1e-10)
  expect_equal(es_normal$d_se, es_reverse$d_se, tolerance = 1e-10)

  # Test g reversal
  expect_equal(es_normal$g, -es_reverse$g, tolerance = 1e-10)
  expect_equal(es_normal$g_se, es_reverse$g_se, tolerance = 1e-10)

  # Test mdw reversal
  expect_equal(es_normal$mdw, -es_reverse$mdw, tolerance = 1e-10)
  expect_equal(es_normal$mdw_se, es_reverse$mdw_se, tolerance = 1e-10)

  # Test r reversal
  expect_equal(es_normal$r, -es_reverse$r, tolerance = 1e-10)
  expect_equal(es_normal$r_se, es_reverse$r_se, tolerance = 1e-10)

  # Test z reversal
  expect_equal(es_normal$z, -es_reverse$z, tolerance = 1e-10)
  expect_equal(es_normal$z_se, es_reverse$z_se, tolerance = 1e-10)

  # Test logor reversal
  expect_equal(es_normal$logor, -es_reverse$logor, tolerance = 1e-10)
  expect_equal(es_normal$logor_se, es_reverse$logor_se, tolerance = 1e-10)
})

test_that("REVERSE - Pre-Post - cooper - single group", {
  # Create test data
  set.seed(888)
  n <- 50
  mean_pre <- 100
  mean_post <- 95
  sd_pre <- 15
  sd_post <- 14
  r_pre_post <- 0.6

  # Normal direction
  es_normal <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "cooper",
    reverse_means_pre_post = FALSE
  )

  # Reversed direction
  es_reverse <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "cooper",
    reverse_means_pre_post = TRUE
  )

  # Test d reversal
  expect_equal(es_normal$d, -es_reverse$d, tolerance = 1e-10)
  expect_equal(es_normal$d_se, es_reverse$d_se, tolerance = 1e-10)

  # Test g reversal
  expect_equal(es_normal$g, -es_reverse$g, tolerance = 1e-10)
  expect_equal(es_normal$g_se, es_reverse$g_se, tolerance = 1e-10)

  # Test mdw reversal
  expect_equal(es_normal$mdw, -es_reverse$mdw, tolerance = 1e-10)
  expect_equal(es_normal$mdw_se, es_reverse$mdw_se, tolerance = 1e-10)
})

### REVERSE - Mean Change -----
test_that("REVERSE - Mean Change - single group", {
  # Create test data
  set.seed(999)
  n <- 50
  mean_change <- 5
  sd_change <- 8
  r_pre_post <- 0.6

  # Normal direction
  es_normal <- es_from_mean_change_sd_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = FALSE
  )

  # Reversed direction
  es_reverse <- es_from_mean_change_sd_single_group(smd_var = "hedges_olkin", 
    mean_change_exp = mean_change,
    mean_change_sd_exp = sd_change,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    reverse_mean_change = TRUE
  )

  # Test d reversal
  expect_equal(es_normal$d, -es_reverse$d, tolerance = 1e-10)
  expect_equal(es_normal$d_se, es_reverse$d_se, tolerance = 1e-10)

  # Test g reversal
  expect_equal(es_normal$g, -es_reverse$g, tolerance = 1e-10)
  expect_equal(es_normal$g_se, es_reverse$g_se, tolerance = 1e-10)

  # Test mdw reversal
  expect_equal(es_normal$mdw, -es_reverse$mdw, tolerance = 1e-10)
  expect_equal(es_normal$mdw_se, es_reverse$mdw_se, tolerance = 1e-10)

  # Test r reversal
  expect_equal(es_normal$r, -es_reverse$r, tolerance = 1e-10)
  expect_equal(es_normal$r_se, es_reverse$r_se, tolerance = 1e-10)

  # Test z reversal
  expect_equal(es_normal$z, -es_reverse$z, tolerance = 1e-10)
  expect_equal(es_normal$z_se, es_reverse$z_se, tolerance = 1e-10)

  # Test logor reversal
  expect_equal(es_normal$logor, -es_reverse$logor, tolerance = 1e-10)
  expect_equal(es_normal$logor_se, es_reverse$logor_se, tolerance = 1e-10)
})

### Test MDW (within-group mean difference) calculation -----
test_that("MDW - Within-group mean difference - manual verification", {
  # Create test data
  n <- 50
  mean_pre <- 100
  mean_post <- 95
  sd_pre <- 15
  sd_post <- 14
  r_pre_post <- 0.6

  # Expected mdw calculations
  expected_mdw <- mean_post - mean_pre
  expected_mdw_var <- (sd_pre^2 + sd_post^2 - 2 * r_pre_post * sd_pre * sd_post) / n
  expected_mdw_se <- sqrt(expected_mdw_var)
  expected_mdw_ci_lo <- expected_mdw - qt(0.975, n - 1) * expected_mdw_se
  expected_mdw_ci_up <- expected_mdw + qt(0.975, n - 1) * expected_mdw_se

  # metaConvert calculation
  es <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Test mdw
  expect_equal(es$mdw, expected_mdw, tolerance = 1e-10)
  expect_equal(es$mdw_se, expected_mdw_se, tolerance = 1e-10)
  expect_equal(es$mdw_ci_lo, expected_mdw_ci_lo, tolerance = 1e-10)
  expect_equal(es$mdw_ci_up, expected_mdw_ci_up, tolerance = 1e-10)
})

### Test conversions to OR, R, Z -----
test_that("Conversions to OR, R, Z - single group", {
  # Create test data
  set.seed(1000)
  n <- 50
  mean_pre <- 100
  mean_post <- 95
  sd_pre <- 15
  sd_post <- 14
  r_pre_post <- 0.6

  # Calculate single-group effect size
  es <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r_pre_post,
    pre_post_to_smd = "bonett",
    reverse_means_pre_post = FALSE
  )

  # Verify conversions are present and non-missing
  expect_true(!is.na(es$logor))
  expect_true(!is.na(es$logor_se))
  expect_true(!is.na(es$r))
  expect_true(!is.na(es$r_se))
  expect_true(!is.na(es$z))
  expect_true(!is.na(es$z_se))

  # Verify d -> r conversion is reasonable (for d ~ 0.33, r should be ~ 0.16)
  # Using point-biserial approximation: r ≈ d / sqrt(d^2 + 4)
  # Note: actual conversion may differ based on the specific formula used
  expected_r_approx <- es$d / sqrt(es$d^2 + 4)
  expect_equal(abs(es$r), abs(expected_r_approx), tolerance = 0.4)
})
