# =============================================================================
# External Validation Tests: PAIRED SINGLE-GROUP Functions
# =============================================================================
#
# Purpose: Comprehensive external validation of all functions in
#          R/es_from_PAIRED_SINGLE_GROUP.R against external packages
#
# Functions tested:
#   1. es_from_means_sd_pre_post_single_group()
#   2. es_from_means_se_pre_post_single_group()
#   3. es_from_means_ci_pre_post_single_group()
#   4. es_from_mean_change_sd_single_group()
#   5. es_from_mean_change_se_single_group()
#   6. es_from_mean_change_ci_single_group()
#   7. es_from_mean_change_pval_single_group()
#   8. es_from_paired_t_single_group()
#
# External validation sources:
#   - metafor::escalc (SMCRH, SMCRP, SMCC, MN)
#   - TOSTER::smd_calc (with raw data generation)
#   - Morris & DeShon (2002) manual formulas
#
# Dataset: metaumbrella::df.SMC (real meta-analysis data)
#
# Tolerance standards:
#   - 1e-9: bonett vs metafor SMCRH (exact match)
#   - 1e-6: morris_dav vs metafor SMCRP (modified df)
#   - 1e-10: Manual formula implementations
#   - 1e-9: TOSTER g (effect size only, NOT SE - see note below)
#
# NOTE on TOSTER:
#   - TOSTER uses different variance formulas for Hedges' g SE
#   - Effect sizes (g) match exactly (tolerance 1e-9)
#   - Standard errors differ by ~2-5% (both formulas are valid)
#   - We test g but NOT g_se against TOSTER
#
# =============================================================================

library(testthat)
library(metaConvert)

# Check for required packages
skip_if_not_installed("metaumbrella")
skip_if_not_installed("metafor")
skip_if_not_installed("TOSTER")

# Load required packages
suppressWarnings(library(metaumbrella))
suppressWarnings(library(metafor))
suppressWarnings(library(TOSTER))

# Prepare base dataset from metaumbrella::df.SMC
# Use first 5 studies for comprehensive tests
data("df.SMC", package = "metaumbrella")
df_raw <- df.SMC[1:5, ]

# Create standardized column names for single-group design (experimental group only)
df_base <- data.frame(
  study_id = seq_len(nrow(df_raw)),
  n_exp = df_raw$n_cases,
  mean_pre_exp = df_raw$mean_pre_cases,
  mean_exp = df_raw$mean_cases,  # Post-intervention mean
  mean_pre_sd_exp = df_raw$sd_pre_cases,
  mean_sd_exp = df_raw$sd_cases,  # Post-intervention SD
  r_pre_post_exp = df_raw$pre_post_cor
)

# =============================================================================
# SECTION 1: es_from_means_sd_pre_post_single_group() - SINGLE GROUP
# =============================================================================

# -----------------------------------------------------------------------------
# Test 1.1: bonett matches metafor SMCRH
# -----------------------------------------------------------------------------
test_that("MEANS-SD-SG: bonett single-group matches metafor SMCRH (g + g_se)", {
  # Prepare data (use only experimental group)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    mean_pre_exp = df_base$mean_pre_exp,
    mean_exp = df_base$mean_exp,
    mean_sd_exp = df_base$mean_sd_exp,
    mean_pre_sd_exp = df_base$mean_pre_sd_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post_single_group",
    measure = "gw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "bonett"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post_single_group")

  # metafor SMCRH calculation
  yi_mf <- vi_mf <- rep(NA, nrow(res))
  for (i in seq_len(nrow(res))) {
    esc <- escalc(
      measure = "SMCRH",
      m1i = res$mean_exp[i],
      m2i = res$mean_pre_exp[i],
      sd1i = res$mean_pre_sd_exp[i],
      sd2i = res$mean_sd_exp[i],
      ni = res$n_exp[i],
      ri = res$r_pre_post_exp[i]
    )
    yi_mf[i] <- esc$yi
    vi_mf[i] <- esc$vi
  }

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, yi_mf, tolerance = 1e-9)
  expect_equal(es_mc$se_crude, sqrt(vi_mf), tolerance = 1e-9)
})

# -----------------------------------------------------------------------------
# Test 1.2: cooper single-group matches Morris & DeShon 2002 formula (d + d_se)
# -----------------------------------------------------------------------------
test_that("MEANS-SD-SG: cooper single-group matches Morris & DeShon 2002 (d + d_se)", {
  # Prepare data (use only experimental group)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    mean_pre_exp = df_base$mean_pre_exp,
    mean_exp = df_base$mean_exp,
    mean_sd_exp = df_base$mean_sd_exp,
    mean_pre_sd_exp = df_base$mean_pre_sd_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post_single_group",
    measure = "dw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "cooper"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post_single_group")

  # Manual calculation using Morris & DeShon (2002) formula
  # d_rm = (mean_post - mean_pre) / sd_change * sqrt(2 * (1 - r))
  # where sd_change = sqrt(sd_pre^2 + sd_post^2 - 2*r*sd_pre*sd_post)
  mean_diff <- res$mean_exp - res$mean_pre_exp
  sd_change <- sqrt(
    res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
    2 * res$r_pre_post_exp * res$mean_pre_sd_exp * res$mean_sd_exp
  )
  d_rm_manual <- (mean_diff / sd_change) * sqrt(2 * (1 - res$r_pre_post_exp))

  # Variance: var(d_rm) = 2*(1-r)/n + d_rm^2/(2*n)
  d_rm_var_manual <- 2 * (1 - res$r_pre_post_exp) / res$n_exp +
                     d_rm_manual^2 / (2 * res$n_exp)
  d_rm_se_manual <- sqrt(d_rm_var_manual)

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, d_rm_manual, tolerance = 1e-10)
  expect_equal(es_mc$se_crude, d_rm_se_manual, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 1.3: morris_dav single-group matches metafor SMCRP (g + g_se)
# -----------------------------------------------------------------------------
test_that("MEANS-SD-SG: morris_dav single-group matches metafor SMCRP (g + g_se)", {
  # Prepare data (use only experimental group)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    mean_pre_exp = df_base$mean_pre_exp,
    mean_exp = df_base$mean_exp,
    mean_sd_exp = df_base$mean_sd_exp,
    mean_pre_sd_exp = df_base$mean_pre_sd_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post_single_group",
    measure = "gw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dav"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post_single_group")

  # metafor SMCRP calculation
  yi_mf <- vi_mf <- rep(NA, nrow(res))
  for (i in seq_len(nrow(res))) {
    esc <- escalc(
      measure = "SMCRP",
      m1i = res$mean_exp[i],
      m2i = res$mean_pre_exp[i],
      sd1i = res$mean_sd_exp[i],
      sd2i = res$mean_pre_sd_exp[i],
      ni = res$n_exp[i],
      ri = res$r_pre_post_exp[i]
    )
    yi_mf[i] <- esc$yi
    vi_mf[i] <- esc$vi
  }

  # Test effect sizes and SEs
  # Use 1e-6 tolerance due to modified df in SMCRP
  expect_equal(es_mc$es_crude, yi_mf, tolerance = 1e-6)
  expect_equal(es_mc$se_crude, sqrt(vi_mf), tolerance = 1e-6)
})

# -----------------------------------------------------------------------------
# Test 1.4: morris_dav single-group matches metafor SMCRP manual formula
# -----------------------------------------------------------------------------
test_that("MEANS-SD-SG: morris_dav single-group matches metafor SMCRP manual (d + d_se)", {
  # Prepare data (use only experimental group)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    mean_pre_exp = df_base$mean_pre_exp,
    mean_exp = df_base$mean_exp,
    mean_sd_exp = df_base$mean_sd_exp,
    mean_pre_sd_exp = df_base$mean_pre_sd_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post_single_group",
    measure = "gw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dav"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post_single_group")

  # Manual calculation following metafor source code
  # d_av = (mean_post - mean_pre) / sqrt((sd_pre^2 + sd_post^2) / 2)
  mean_diff <- res$mean_exp - res$mean_pre_exp
  sd_av <- sqrt((res$mean_pre_sd_exp^2 + res$mean_sd_exp^2) / 2)
  d_av_manual <- mean_diff / sd_av

  # Modified df for bias correction (Bonett 2008)
  mi <- 2 * (res$n_exp - 1) / (1 + res$r_pre_post_exp^2)
  J <- exp(lgamma(mi / 2) - log(sqrt(mi / 2)) - lgamma((mi - 1) / 2))
  g_av_manual <- d_av_manual * J

  # Variance (metaConvert formula matching metafor SMCRP)
  # var(g) = 2*(1-r)/n + g^2*(1+r^2)/(4*n)
  g_av_var_manual <- 2 * (1 - res$r_pre_post_exp) / res$n_exp +
    g_av_manual^2 * (1 + res$r_pre_post_exp^2) / (4 * res$n_exp)
  g_av_se_manual <- sqrt(g_av_var_manual)

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, g_av_manual, tolerance = 1e-10)
  expect_equal(es_mc$se_crude, g_av_se_manual, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 1.5: morris_dz single-group matches metafor SMCC manual formula (d + d_se)
# -----------------------------------------------------------------------------
test_that("MEANS-SD-SG: morris_dz single-group matches metafor SMCC manual (d + d_se)", {
  # Prepare data (use only experimental group)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    mean_pre_exp = df_base$mean_pre_exp,
    mean_exp = df_base$mean_exp,
    mean_sd_exp = df_base$mean_sd_exp,
    mean_pre_sd_exp = df_base$mean_pre_sd_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post_single_group",
    measure = "gw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dz"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post_single_group")

  # Manual calculation using SMCC formula
  # d_z = (mean_post - mean_pre) / sd_change
  # where sd_change = sqrt(sd_pre^2 + sd_post^2 - 2*r*sd_pre*sd_post)
  mean_diff <- res$mean_exp - res$mean_pre_exp
  sd_change <- sqrt(
    res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
    2 * res$r_pre_post_exp * res$mean_pre_sd_exp * res$mean_sd_exp
  )
  d_z_manual <- mean_diff / sd_change

  # Bias correction
  J <- exp(lgamma((res$n_exp - 1) / 2) - log(sqrt((res$n_exp - 1) / 2)) -
           lgamma((res$n_exp - 2) / 2))
  g_z_manual <- d_z_manual * J

  # Variance: var(g_z) = 1/n + g_z^2/(2*n)
  g_z_var_manual <- 1 / res$n_exp + g_z_manual^2 / (2 * res$n_exp)
  g_z_se_manual <- sqrt(g_z_var_manual)

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, g_z_manual, tolerance = 1e-10)
  expect_equal(es_mc$se_crude, g_z_se_manual, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 1.6: cooper single-group matches TOSTER (g only, not SE)
# -----------------------------------------------------------------------------
test_that("MEANS-SD-SG: cooper single-group matches TOSTER (g only)", {
  skip_if_not_installed("MASS")

  # Generate raw paired data with known correlation structure
  scenarios <- list(
    list(n = 30, mu1 = 50, mu2 = 60, sd1 = 10, sd2 = 12, r = 0.7, seed = 42),
    list(n = 50, mu1 = 20, mu2 = 28, sd1 = 8, sd2 = 10, r = 0.8, seed = 7),
    list(n = 100, mu1 = 0, mu2 = 0.5, sd1 = 1, sd2 = 1.2, r = 0.3, seed = 99)
  )

  for (sc in scenarios) {
    set.seed(sc$seed)
    Sigma <- matrix(c(sc$sd1^2, sc$r * sc$sd1 * sc$sd2,
                       sc$r * sc$sd1 * sc$sd2, sc$sd2^2), 2, 2)
    dat <- MASS::mvrnorm(sc$n, c(sc$mu1, sc$mu2), Sigma)
    pre <- dat[, 1]
    post <- dat[, 2]

    # TOSTER g_rm: rm_correction=TRUE applies sqrt(2(1-r)) factor (= cooper/morris_drm)
    # Sign convention: TOSTER computes x-y (pre-post), metaConvert computes post-pre
    toster_g <- TOSTER::smd_calc(
      x = pre, y = post,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = TRUE, smd_ci = "nct"
    )

    # metaConvert from summary stats via convert_df pipeline
    res <- data.frame(
      n_exp = length(pre),
      mean_pre_exp = mean(pre),
      mean_exp = mean(post),
      mean_pre_sd_exp = sd(pre),
      mean_sd_exp = sd(post),
      r_pre_post_exp = cor(pre, post)
    )

    es_mc <- summary(convert_df(res,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_sd_pre_post_single_group",
      measure = "gw",
      pre_post_to_smd = "cooper"
    ), digits = 11)

    expect_equal(es_mc$es_crude, -as.numeric(toster_g$estimate), tolerance = 1e-9,
                 label = paste("cooper g, seed =", sc$seed))
  }
})

# -----------------------------------------------------------------------------
# Test 1.7: morris_dz single-group matches TOSTER (g only, not SE)
# -----------------------------------------------------------------------------
test_that("MEANS-SD-SG: morris_dz single-group matches TOSTER (g only)", {
  skip_if_not_installed("MASS")

  # Generate raw paired data with known correlation structure
  scenarios <- list(
    list(n = 30, mu1 = 50, mu2 = 60, sd1 = 10, sd2 = 12, r = 0.7, seed = 42),
    list(n = 50, mu1 = 20, mu2 = 28, sd1 = 8, sd2 = 10, r = 0.8, seed = 7),
    list(n = 100, mu1 = 0, mu2 = 0.5, sd1 = 1, sd2 = 1.2, r = 0.3, seed = 99)
  )

  for (sc in scenarios) {
    set.seed(sc$seed)
    Sigma <- matrix(c(sc$sd1^2, sc$r * sc$sd1 * sc$sd2,
                       sc$r * sc$sd1 * sc$sd2, sc$sd2^2), 2, 2)
    dat <- MASS::mvrnorm(sc$n, c(sc$mu1, sc$mu2), Sigma)
    pre <- dat[, 1]
    post <- dat[, 2]

    # TOSTER g_z: rm_correction=FALSE gives d_z (standardized by SD of differences)
    # Sign convention: TOSTER computes x-y (pre-post), metaConvert computes post-pre
    toster_g <- TOSTER::smd_calc(
      x = pre, y = post,
      paired = TRUE, rm_correction = FALSE,
      bias_correction = TRUE, smd_ci = "nct"
    )

    # metaConvert from summary stats via convert_df pipeline
    res <- data.frame(
      n_exp = length(pre),
      mean_pre_exp = mean(pre),
      mean_exp = mean(post),
      mean_pre_sd_exp = sd(pre),
      mean_sd_exp = sd(post),
      r_pre_post_exp = cor(pre, post)
    )

    es_mc <- summary(convert_df(res,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_sd_pre_post_single_group",
      measure = "gw",
      pre_post_to_smd = "morris_dz"
    ), digits = 11)

    expect_equal(es_mc$es_crude, -as.numeric(toster_g$estimate), tolerance = 1e-9,
                 label = paste("morris_dz g, seed =", sc$seed))
  }
})

# =============================================================================
# SECTION 2: es_from_means_se_pre_post_single_group() - SINGLE GROUP
# =============================================================================

test_that("MEANS-SE-SG: all methods match SD version after SE→SD conversion", {
  methods <- c("bonett", "cooper", "morris_dav", "morris_dz")

  for (method in methods) {
    # Prepare data with SD (use only experimental group)
    res_sd <- data.frame(
      study_id = df_base$study_id[1:3],
      n_exp = df_base$n_exp[1:3],
      mean_pre_exp = df_base$mean_pre_exp[1:3],
      mean_exp = df_base$mean_exp[1:3],
      mean_sd_exp = df_base$mean_sd_exp[1:3],
      mean_pre_sd_exp = df_base$mean_pre_sd_exp[1:3],
      r_pre_post_exp = df_base$r_pre_post[1:3]
    )

    # From SD
    es_sd <- summary(convert_df(res_sd,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_sd_pre_post_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Prepare data with SE (remove SD columns to force SE usage)
    res_se <- res_sd
    res_se$mean_pre_se_exp <- res_se$mean_pre_sd_exp / sqrt(res_se$n_exp)
    res_se$mean_se_exp <- res_se$mean_sd_exp / sqrt(res_se$n_exp)
    res_se$mean_pre_sd_exp <- NULL
    res_se$mean_sd_exp <- NULL

    # From SE
    es_se <- summary(convert_df(res_se,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_se_pre_post_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Verify correct hierarchies were used
    expect_equal(unique(es_sd$info_used_crude), "means_sd_pre_post_single_group")
    expect_equal(unique(es_se$info_used_crude), "means_se_pre_post_single_group")

    # Test equivalence
    expect_equal(es_sd$es_crude, es_se$es_crude, tolerance = 1e-10,
                 label = paste0("SE→SD equivalence for method: ", method))
    expect_equal(es_sd$se_crude, es_se$se_crude, tolerance = 1e-10,
                 label = paste0("SE→SD SE equivalence for method: ", method))
  }
})

# =============================================================================
# SECTION 3: es_from_means_ci_pre_post_single_group() - SINGLE GROUP
# =============================================================================

test_that("MEANS-CI-SG: all methods match SD version after CI→SD conversion", {
  methods <- c("bonett", "cooper", "morris_dav", "morris_dz")

  for (method in methods) {
    # Prepare data with SD (use only experimental group)
    res_sd <- data.frame(
      study_id = df_base$study_id[1:3],
      n_exp = df_base$n_exp[1:3],
      mean_pre_exp = df_base$mean_pre_exp[1:3],
      mean_exp = df_base$mean_exp[1:3],
      mean_sd_exp = df_base$mean_sd_exp[1:3],
      mean_pre_sd_exp = df_base$mean_pre_sd_exp[1:3],
      r_pre_post_exp = df_base$r_pre_post[1:3]
    )

    # From SD
    es_sd <- summary(convert_df(res_sd,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_sd_pre_post_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Prepare data with CI (remove SD and SE columns to force CI usage)
    res_ci <- res_sd
    se_pre <- res_ci$mean_pre_sd_exp / sqrt(res_ci$n_exp)
    se_post <- res_ci$mean_sd_exp / sqrt(res_ci$n_exp)
    res_ci$mean_pre_ci_lo_exp <- res_ci$mean_pre_exp - qt(0.975, res_ci$n_exp - 1) * se_pre
    res_ci$mean_pre_ci_up_exp <- res_ci$mean_pre_exp + qt(0.975, res_ci$n_exp - 1) * se_pre
    res_ci$mean_ci_lo_exp <- res_ci$mean_exp - qt(0.975, res_ci$n_exp - 1) * se_post
    res_ci$mean_ci_up_exp <- res_ci$mean_exp + qt(0.975, res_ci$n_exp - 1) * se_post
    res_ci$mean_pre_sd_exp <- NULL
    res_ci$mean_sd_exp <- NULL

    # From CI
    es_ci <- summary(convert_df(res_ci,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_ci_pre_post_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Verify correct hierarchies were used
    expect_equal(unique(es_sd$info_used_crude), "means_sd_pre_post_single_group")
    expect_equal(unique(es_ci$info_used_crude), "means_ci_pre_post_single_group")

    # Test equivalence (use 1e-9 tolerance for CI due to t-distribution)
    expect_equal(es_sd$es_crude, es_ci$es_crude, tolerance = 1e-9,
                 label = paste0("CI→SD equivalence for method: ", method))
    expect_equal(es_sd$se_crude, es_ci$se_crude, tolerance = 1e-9,
                 label = paste0("CI→SD SE equivalence for method: ", method))
  }
})

# =============================================================================
# SECTION 4: es_from_mean_change_sd_single_group()
# =============================================================================

# -----------------------------------------------------------------------------
# Test 4.1: cooper matches Morris & DeShon 2002 formula
# -----------------------------------------------------------------------------
test_that("MC-SD-SG: cooper matches Morris & DeShon 2002 (d + d_se)", {
  # Prepare data (compute mean change from pre/post)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # Compute mean change and SD
  res$mean_change_exp <- df_base$mean_exp - df_base$mean_pre_exp
  res$mean_change_sd_exp <- sqrt(
    df_base$mean_pre_sd_exp^2 + df_base$mean_sd_exp^2 -
    2 * res$r_pre_post_exp * df_base$mean_pre_sd_exp * df_base$mean_sd_exp
  )

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "mean_change_sd_single_group",
    measure = "dw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "cooper"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "mean_change_sd_single_group")

  # Manual calculation using Morris & DeShon (2002) formula
  # d_rm = (mean_change / sd_change) * sqrt(2 * (1 - r))
  d_rm_manual <- (res$mean_change_exp / res$mean_change_sd_exp) *
                 sqrt(2 * (1 - res$r_pre_post_exp))

  # Variance: var(d_rm) = 2*(1-r)/n + d_rm^2/(2*n)
  d_rm_var_manual <- 2 * (1 - res$r_pre_post_exp) / res$n_exp +
                     d_rm_manual^2 / (2 * res$n_exp)
  d_rm_se_manual <- sqrt(d_rm_var_manual)

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, d_rm_manual, tolerance = 1e-10)
  expect_equal(es_mc$se_crude, d_rm_se_manual, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 4.2: morris_dz matches metafor SMCC manual formula
# -----------------------------------------------------------------------------
test_that("MC-SD-SG: morris_dz matches metafor SMCC manual (d + d_se)", {
  # Prepare data (compute mean change from pre/post)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # Compute mean change and SD
  res$mean_change_exp <- df_base$mean_exp - df_base$mean_pre_exp
  res$mean_change_sd_exp <- sqrt(
    df_base$mean_pre_sd_exp^2 + df_base$mean_sd_exp^2 -
    2 * res$r_pre_post_exp * df_base$mean_pre_sd_exp * df_base$mean_sd_exp
  )

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "mean_change_sd_single_group",
    measure = "gw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dz"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "mean_change_sd_single_group")

  # Manual calculation using SMCC formula
  # d_z = mean_change / sd_change
  d_z_manual <- res$mean_change_exp / res$mean_change_sd_exp

  # Bias correction
  J <- exp(lgamma((res$n_exp - 1) / 2) - log(sqrt((res$n_exp - 1) / 2)) -
           lgamma((res$n_exp - 2) / 2))
  g_z_manual <- d_z_manual * J

  # Variance: var(g_z) = 1/n + g_z^2/(2*n)
  g_z_var_manual <- 1 / res$n_exp + g_z_manual^2 / (2 * res$n_exp)
  g_z_se_manual <- sqrt(g_z_var_manual)

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, g_z_manual, tolerance = 1e-10)
  expect_equal(es_mc$se_crude, g_z_se_manual, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 4.3: MDw matches metafor MN
# -----------------------------------------------------------------------------
test_that("MC-SD-SG: MDw matches metafor MN (mdw + mdw_se)", {
  # Prepare data (compute mean change from pre/post)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # Compute mean change and SD
  res$mean_change_exp <- df_base$mean_exp - df_base$mean_pre_exp
  res$mean_change_sd_exp <- sqrt(
    df_base$mean_pre_sd_exp^2 + df_base$mean_sd_exp^2 -
    2 * res$r_pre_post_exp * df_base$mean_pre_sd_exp * df_base$mean_sd_exp
  )

  # metaConvert using convert_df - get MDw
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "mean_change_sd_single_group",
    measure = "mdw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "cooper"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "mean_change_sd_single_group")

  # metafor MN calculation
  yi_mf <- vi_mf <- rep(NA, nrow(res))
  for (i in seq_len(nrow(res))) {
    esc <- escalc(
      measure = "MN",
      mi = res$mean_change_exp[i],
      sdi = res$mean_change_sd_exp[i],
      ni = res$n_exp[i]
    )
    yi_mf[i] <- esc$yi
    vi_mf[i] <- esc$vi
  }

  # Test mean difference and SE
  # Use 1e-6 tolerance for metafor MN
  expect_equal(es_mc$es_crude, yi_mf, tolerance = 1e-6)
  expect_equal(es_mc$se_crude, sqrt(vi_mf), tolerance = 1e-6)
})

# -----------------------------------------------------------------------------
# Test 4.4: Equivalence with means_sd_pre_post_single_group (mean_pre=0)
# -----------------------------------------------------------------------------
test_that("MC-SD-SG: equivalence with means_sd_pre_post_single_group (mean_pre=0)", {
  methods <- c("cooper", "morris_dz")

  for (method in methods) {
    # Prepare mean change data
    res_mc <- data.frame(
      study_id = df_base$study_id[1:3],
      n_exp = df_base$n_exp[1:3],
      r_pre_post_exp = df_base$r_pre_post[1:3]
    )

    res_mc$mean_change_exp <- df_base$mean_exp[1:3] - df_base$mean_pre_exp[1:3]
    res_mc$mean_change_sd_exp <- sqrt(
      df_base$mean_pre_sd_exp[1:3]^2 + df_base$mean_sd_exp[1:3]^2 -
      2 * res_mc$r_pre_post_exp * df_base$mean_pre_sd_exp[1:3] * df_base$mean_sd_exp[1:3]
    )

    # From mean change
    es_mc <- summary(convert_df(res_mc,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "mean_change_sd_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Prepare pre/post data with mean_pre=0
    res_pp <- data.frame(
      study_id = df_base$study_id[1:3],
      n_exp = df_base$n_exp[1:3],
      mean_pre_exp = 0,
      mean_exp = res_mc$mean_change_exp,
      mean_pre_sd_exp = 0,
      mean_sd_exp = res_mc$mean_change_sd_exp,
      r_pre_post_exp = res_mc$r_pre_post_exp
    )

    # From pre/post with mean_pre=0. The zero baseline SD is a deliberate
    # degenerate construction for this equivalence check; the centralized
    # input validation would (rightly) zap SD = 0 as a data-entry error, so
    # disable input correction for this synthetic dataset.
    es_pp <- summary(convert_df(res_pp,
      verbose = FALSE,
      correct_inputs = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_sd_pre_post_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Verify correct hierarchies were used
    expect_equal(unique(es_mc$info_used_crude), "mean_change_sd_single_group")
    expect_equal(unique(es_pp$info_used_crude), "means_sd_pre_post_single_group")

    # Test equivalence (machine precision)
    expect_equal(es_mc$es_crude, es_pp$es_crude, tolerance = 1e-15,
                 label = paste0("MC ≡ pre/post equivalence for method: ", method))
    expect_equal(es_mc$se_crude, es_pp$se_crude, tolerance = 1e-15,
                 label = paste0("MC ≡ pre/post SE equivalence for method: ", method))
  }
})

# =============================================================================
# SECTION 5: es_from_mean_change_se_single_group()
# =============================================================================

test_that("MC-SE-SG: all methods match SD version after SE→SD conversion", {
  methods <- c("cooper", "morris_dz")

  for (method in methods) {
    # Prepare mean change data with SD
    res_sd <- data.frame(
      study_id = df_base$study_id[1:3],
      n_exp = df_base$n_exp[1:3],
      r_pre_post_exp = df_base$r_pre_post[1:3]
    )

    res_sd$mean_change_exp <- df_base$mean_exp[1:3] - df_base$mean_pre_exp[1:3]
    res_sd$mean_change_sd_exp <- sqrt(
      df_base$mean_pre_sd_exp[1:3]^2 + df_base$mean_sd_exp[1:3]^2 -
      2 * res_sd$r_pre_post_exp * df_base$mean_pre_sd_exp[1:3] * df_base$mean_sd_exp[1:3]
    )

    # From SD
    es_sd <- summary(convert_df(res_sd,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "mean_change_sd_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Prepare data with SE (remove SD columns to force SE usage)
    res_se <- res_sd
    res_se$mean_change_se_exp <- res_se$mean_change_sd_exp / sqrt(res_se$n_exp)
    res_se$mean_change_sd_exp <- NULL

    # From SE
    es_se <- summary(convert_df(res_se,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "mean_change_se_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Verify correct hierarchies were used
    expect_equal(unique(es_sd$info_used_crude), "mean_change_sd_single_group")
    expect_equal(unique(es_se$info_used_crude), "mean_change_se_single_group")

    # Test equivalence
    expect_equal(es_sd$es_crude, es_se$es_crude, tolerance = 1e-10,
                 label = paste0("SE→SD equivalence for method: ", method))
    expect_equal(es_sd$se_crude, es_se$se_crude, tolerance = 1e-10,
                 label = paste0("SE→SD SE equivalence for method: ", method))
  }
})

# =============================================================================
# SECTION 6: es_from_mean_change_ci_single_group()
# =============================================================================

test_that("MC-CI-SG: all methods match SD version after CI→SD conversion", {
  methods <- c("cooper", "morris_dz")

  for (method in methods) {
    # Prepare mean change data with SD
    res_sd <- data.frame(
      study_id = df_base$study_id[1:3],
      n_exp = df_base$n_exp[1:3],
      r_pre_post_exp = df_base$r_pre_post[1:3]
    )

    res_sd$mean_change_exp <- df_base$mean_exp[1:3] - df_base$mean_pre_exp[1:3]
    res_sd$mean_change_sd_exp <- sqrt(
      df_base$mean_pre_sd_exp[1:3]^2 + df_base$mean_sd_exp[1:3]^2 -
      2 * res_sd$r_pre_post_exp * df_base$mean_pre_sd_exp[1:3] * df_base$mean_sd_exp[1:3]
    )

    # From SD
    es_sd <- summary(convert_df(res_sd,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "mean_change_sd_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Prepare data with CI (remove SD and SE columns to force CI usage)
    res_ci <- res_sd
    se_change <- res_ci$mean_change_sd_exp / sqrt(res_ci$n_exp)
    res_ci$mean_change_ci_lo_exp <- res_ci$mean_change_exp - qt(0.975, res_ci$n_exp - 1) * se_change
    res_ci$mean_change_ci_up_exp <- res_ci$mean_change_exp + qt(0.975, res_ci$n_exp - 1) * se_change
    res_ci$mean_change_sd_exp <- NULL

    # From CI
    es_ci <- summary(convert_df(res_ci,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "mean_change_ci_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Verify correct hierarchies were used
    expect_equal(unique(es_sd$info_used_crude), "mean_change_sd_single_group")
    expect_equal(unique(es_ci$info_used_crude), "mean_change_ci_single_group")

    # Test equivalence (use 1e-9 tolerance for CI due to t-distribution)
    expect_equal(es_sd$es_crude, es_ci$es_crude, tolerance = 1e-9,
                 label = paste0("CI→SD equivalence for method: ", method))
    expect_equal(es_sd$se_crude, es_ci$se_crude, tolerance = 1e-9,
                 label = paste0("CI→SD SE equivalence for method: ", method))
  }
})

# =============================================================================
# SECTION 7: es_from_mean_change_pval_single_group()
# =============================================================================

test_that("MC-PVAL-SG: all methods match SD version after pval→t→SE→SD conversion", {
  methods <- c("cooper", "morris_dz")

  for (method in methods) {
    # Prepare mean change data with SD
    res_sd <- data.frame(
      study_id = df_base$study_id[1:3],
      n_exp = df_base$n_exp[1:3],
      r_pre_post_exp = df_base$r_pre_post[1:3]
    )

    res_sd$mean_change_exp <- df_base$mean_exp[1:3] - df_base$mean_pre_exp[1:3]
    res_sd$mean_change_sd_exp <- sqrt(
      df_base$mean_pre_sd_exp[1:3]^2 + df_base$mean_sd_exp[1:3]^2 -
      2 * res_sd$r_pre_post_exp * df_base$mean_pre_sd_exp[1:3] * df_base$mean_sd_exp[1:3]
    )

    # From SD
    es_sd <- summary(convert_df(res_sd,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "mean_change_sd_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Compute p-value from mean change and SD
    t_stat <- res_sd$mean_change_exp / (res_sd$mean_change_sd_exp / sqrt(res_sd$n_exp))
    pval <- 2 * pt(abs(t_stat), res_sd$n_exp - 1, lower.tail = FALSE)

    # Prepare data with pval (remove SD columns to force pval usage)
    res_pval <- data.frame(
      study_id = res_sd$study_id,
      n_exp = res_sd$n_exp,
      mean_change_exp = res_sd$mean_change_exp,
      mean_change_pval_exp = pval,
      r_pre_post_exp = res_sd$r_pre_post_exp
    )

    # From pval
    es_pval <- summary(convert_df(res_pval,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "mean_change_pval_single_group",
      measure = "gw",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Verify correct hierarchies were used
    expect_equal(unique(es_sd$info_used_crude), "mean_change_sd_single_group")
    expect_equal(unique(es_pval$info_used_crude), "mean_change_pval_single_group")

    # Test equivalence (use 1e-9 tolerance for pval conversion chain)
    expect_equal(es_sd$es_crude, es_pval$es_crude, tolerance = 1e-9,
                 label = paste0("pval→t→SE→SD equivalence for method: ", method))
    expect_equal(es_sd$se_crude, es_pval$se_crude, tolerance = 1e-9,
                 label = paste0("pval→t→SE→SD SE equivalence for method: ", method))
  }
})

# =============================================================================
# SECTION 8: es_from_paired_t_single_group()
# =============================================================================

# -----------------------------------------------------------------------------
# Test 8.1: cooper matches Morris & DeShon 2002 formula
# -----------------------------------------------------------------------------
test_that("PAIRED-T-SG: cooper matches Morris & DeShon 2002 (d + d_se)", {
  # Prepare data (compute paired t from pre/post data)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # Compute paired t-statistic
  mean_diff <- df_base$mean_exp - df_base$mean_pre_exp
  sd_change <- sqrt(
    df_base$mean_pre_sd_exp^2 + df_base$mean_sd_exp^2 -
    2 * res$r_pre_post_exp * df_base$mean_pre_sd_exp * df_base$mean_sd_exp
  )
  res$paired_t_exp <- mean_diff / (sd_change / sqrt(res$n_exp))

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "paired_t_single_group",
    measure = "dw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "cooper"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "paired_t_single_group")

  # Manual calculation using Morris & DeShon (2002) formula
  # d_rm = t * sqrt(2*(1-r)/n)
  d_rm_manual <- res$paired_t_exp * sqrt(2 * (1 - res$r_pre_post_exp) / res$n_exp)

  # Variance: var(d_rm) = 2*(1-r)/n + d_rm^2/(2*n)
  d_rm_var_manual <- 2 * (1 - res$r_pre_post_exp) / res$n_exp +
                     d_rm_manual^2 / (2 * res$n_exp)
  d_rm_se_manual <- sqrt(d_rm_var_manual)

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, d_rm_manual, tolerance = 1e-10)
  expect_equal(es_mc$se_crude, d_rm_se_manual, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 8.2: morris_dz matches metafor SMCC manual formula
# -----------------------------------------------------------------------------
test_that("PAIRED-T-SG: morris_dz matches metafor SMCC manual (d + d_se)", {
  # Prepare data (compute paired t from pre/post data)
  res <- data.frame(
    study_id = df_base$study_id,
    n_exp = df_base$n_exp,
    r_pre_post_exp = df_base$r_pre_post
  )

  # Compute paired t-statistic
  mean_diff <- df_base$mean_exp - df_base$mean_pre_exp
  sd_change <- sqrt(
    df_base$mean_pre_sd_exp^2 + df_base$mean_sd_exp^2 -
    2 * res$r_pre_post_exp * df_base$mean_pre_sd_exp * df_base$mean_sd_exp
  )
  res$paired_t_exp <- mean_diff / (sd_change / sqrt(res$n_exp))

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "paired_t_single_group",
    measure = "gw",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dz"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "paired_t_single_group")

  # Manual calculation using SMCC formula
  # d_z = t / sqrt(n)
  d_z_manual <- res$paired_t_exp / sqrt(res$n_exp)

  # Bias correction
  J <- exp(lgamma((res$n_exp - 1) / 2) - log(sqrt((res$n_exp - 1) / 2)) -
           lgamma((res$n_exp - 2) / 2))
  g_z_manual <- d_z_manual * J

  # Variance: var(g_z) = J^2 * var(d_z) = J^2/n + g_z^2/(2*n)
  # Note: paired_t_single_group uses g_se = d_se * J, so J^2 applies to full variance
  g_z_var_manual <- J^2 / res$n_exp + g_z_manual^2 / (2 * res$n_exp)
  g_z_se_manual <- sqrt(g_z_var_manual)

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, g_z_manual, tolerance = 1e-10)
  expect_equal(es_mc$se_crude, g_z_se_manual, tolerance = 1e-10)
})
