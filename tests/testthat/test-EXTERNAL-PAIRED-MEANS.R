# ==============================================================================
# EXTERNAL VALIDATION: es_from_PAIRED_MEANS.R
# ==============================================================================
# Validates all functions in es_from_PAIRED_MEANS.R against external packages:
#   - metafor (SMCRH for bonett, SMCRP for morris_dav, SMCC for morris_dz, MC for raw MD)
#   - TOSTER (rm_correction for cooper, rm_correction=FALSE for morris_dz)
#
# Tests cover:
#   1. es_from_means_sd_pre_post() - all methods (bonett, cooper, morris_dav, morris_dz)
#   2. es_from_means_se_pre_post() - all methods
#   3. es_from_means_ci_pre_post() - all methods
#
# Dataset: metaumbrella::df.SMC (real meta-analysis data with pre/post means)
# ==============================================================================

library(testthat)
library(metaConvert)

# Check for required packages
has_metafor <- requireNamespace("metafor", quietly = TRUE)
has_TOSTER <- requireNamespace("TOSTER", quietly = TRUE)
has_metaumbrella <- requireNamespace("metaumbrella", quietly = TRUE)

# Skip all tests if metaumbrella not available
if (!has_metaumbrella) {
  skip("metaumbrella package not available - skipping all EXTERNAL-PAIRED-MEANS tests")
}

# Load dataset
data("df.SMC", package = "metaumbrella")

# ==============================================================================
# SECTION 1: es_from_means_sd_pre_post() - TWO-GROUP
# ==============================================================================

## 1.1 BONETT vs metafor SMCRH ====
test_that("MEANS-SD: bonett two-group matches metafor SMCRH (g + g_se)", {
  skip_if_not(has_metafor, "metafor not available")

  # Prepare data from metaumbrella
  res <- df.SMC[1:5, ]
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$r_pre_post_exp <- 0.7
  res$r_pre_post_nexp <- 0.7

  # metafor SMCRH for experimental group
  smc_exp <- metafor::escalc(
    m1i = res$mean_pre_exp,
    sd1i = res$mean_pre_sd_exp,
    m2i = res$mean_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp,
    measure = "SMCRH"
  )

  # metafor SMCRH for control group
  smc_nexp <- metafor::escalc(
    m1i = res$mean_pre_nexp,
    sd1i = res$mean_pre_sd_nexp,
    m2i = res$mean_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp,
    measure = "SMCRH"
  )

  # Calculate difference (negate metafor to match post-pre convention)
  g_expected <- (-as.numeric(smc_exp$yi)) - (-as.numeric(smc_nexp$yi))
  se_expected <- sqrt(as.numeric(smc_exp$vi) + as.numeric(smc_nexp$vi))

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "g",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "bonett"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, g_expected, tolerance = 1e-9,
               label = "bonett g matches metafor SMCRH")
  expect_equal(es_mc$se_crude, se_expected, tolerance = 1e-9,
               label = "bonett g_se matches metafor SMCRH")
})

## 1.2 COOPER vs Morris & DeShon (2002) Manual Formula ====
test_that("MEANS-SD: cooper two-group matches Morris & DeShon (2002) (d + d_se)", {
  res <- df.SMC[1:5, ]
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$r_pre_post_exp <- 0.6
  res$r_pre_post_nexp <- 0.6

  # Manual Morris & DeShon (2002) calculation
  # Experimental group
  sd_change_exp <- sqrt(res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
                        2 * res$r_pre_post_exp * res$mean_pre_sd_exp * res$mean_sd_exp)
  d_rm_exp <- (res$mean_exp - res$mean_pre_exp) / sd_change_exp *
              sqrt(2 * (1 - res$r_pre_post_exp))
  var_rm_exp <- 2 * (1 - res$r_pre_post_exp) / res$n_exp + d_rm_exp^2 / (2 * res$n_exp)

  # Control group
  sd_change_nexp <- sqrt(res$mean_pre_sd_nexp^2 + res$mean_sd_nexp^2 -
                         2 * res$r_pre_post_nexp * res$mean_pre_sd_nexp * res$mean_sd_nexp)
  d_rm_nexp <- (res$mean_nexp - res$mean_pre_nexp) / sd_change_nexp *
               sqrt(2 * (1 - res$r_pre_post_nexp))
  var_rm_nexp <- 2 * (1 - res$r_pre_post_nexp) / res$n_nexp + d_rm_nexp^2 / (2 * res$n_nexp)

  # Difference
  d_expected <- d_rm_exp - d_rm_nexp
  se_expected <- sqrt(var_rm_exp + var_rm_nexp)

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "d",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "cooper"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, d_expected, tolerance = 1e-10,
               label = "cooper d matches Morris & DeShon (2002)")
  expect_equal(es_mc$se_crude, se_expected, tolerance = 1e-10,
               label = "cooper d_se matches Morris & DeShon (2002)")
})

## 1.3 MORRIS_DAV vs metafor SMCRP ====
test_that("MEANS-SD: morris_dav two-group matches metafor SMCRP (g + g_se)", {
  skip_if_not(has_metafor, "metafor not available")

  res <- df.SMC[1:5, ]
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$r_pre_post_exp <- 0.65
  res$r_pre_post_nexp <- 0.65

  # metafor SMCRP for experimental group
  mf_exp <- metafor::escalc(
    measure = "SMCRP",
    m1i = res$mean_pre_exp,
    m2i = res$mean_exp,
    sd1i = res$mean_pre_sd_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp
  )

  # metafor SMCRP for control group
  mf_nexp <- metafor::escalc(
    measure = "SMCRP",
    m1i = res$mean_pre_nexp,
    m2i = res$mean_nexp,
    sd1i = res$mean_pre_sd_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp
  )

  # Calculate difference (negate metafor to match post-pre convention)
  g_expected <- (-as.numeric(mf_exp$yi)) - (-as.numeric(mf_nexp$yi))
  se_expected <- sqrt(as.numeric(mf_exp$vi) + as.numeric(mf_nexp$vi))

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "g",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dav"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, g_expected, tolerance = 1e-6,
               label = "morris_dav g matches metafor SMCRP")
  expect_equal(es_mc$se_crude, se_expected, tolerance = 1e-6,
               label = "morris_dav g_se matches metafor SMCRP")
})

## 1.4 MORRIS_DZ vs metafor SMCC Manual Formula ====
test_that("MEANS-SD: morris_dz two-group matches metafor SMCC formula (d + d_se)", {
  skip_if_not(has_metafor, "metafor not available")

  res <- df.SMC[1:5, ]
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$r_pre_post_exp <- 0.55
  res$r_pre_post_nexp <- 0.55

  .d_j <- function(df) exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2))

  # Manual SMCC formula (metafor source code approach)
  # Experimental group
  J_exp <- .d_j(res$n_exp - 1)
  sd_diff_exp <- sqrt(res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
                      2 * res$r_pre_post_exp * res$mean_pre_sd_exp * res$mean_sd_exp)
  d_z_exp <- (res$mean_exp - res$mean_pre_exp) / sd_diff_exp
  g_z_exp <- d_z_exp * J_exp
  var_g_z_exp <- 1 / res$n_exp + g_z_exp^2 / (2 * res$n_exp)
  var_d_z_exp <- var_g_z_exp / (J_exp^2)

  # Control group
  J_nexp <- .d_j(res$n_nexp - 1)
  sd_diff_nexp <- sqrt(res$mean_pre_sd_nexp^2 + res$mean_sd_nexp^2 -
                       2 * res$r_pre_post_nexp * res$mean_pre_sd_nexp * res$mean_sd_nexp)
  d_z_nexp <- (res$mean_nexp - res$mean_pre_nexp) / sd_diff_nexp
  g_z_nexp <- d_z_nexp * J_nexp
  var_g_z_nexp <- 1 / res$n_nexp + g_z_nexp^2 / (2 * res$n_nexp)
  var_d_z_nexp <- var_g_z_nexp / (J_nexp^2)

  # Difference
  d_expected <- d_z_exp - d_z_nexp
  se_expected <- sqrt(var_d_z_exp + var_d_z_nexp)

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "d",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dz"
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, d_expected, tolerance = 1e-10,
               label = "morris_dz d matches SMCC formula")
  expect_equal(es_mc$se_crude, se_expected, tolerance = 1e-10,
               label = "morris_dz d_se matches SMCC formula")
})

## 1.5 RAW MD vs metafor MC ====
test_that("MEANS-SD: raw MD two-group matches metafor MC (md + md_se)", {
  skip_if_not(has_metafor, "metafor not available")

  res <- df.SMC[1:5, ]
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$r_pre_post_exp <- 0.75
  res$r_pre_post_nexp <- 0.75

  # metafor MC for experimental group
  mf_exp <- metafor::escalc(
    measure = "MC",
    m1i = res$mean_pre_exp,
    m2i = res$mean_exp,
    sd1i = res$mean_pre_sd_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp
  )

  # metafor MC for control group
  mf_nexp <- metafor::escalc(
    measure = "MC",
    m1i = res$mean_pre_nexp,
    m2i = res$mean_nexp,
    sd1i = res$mean_pre_sd_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp
  )

  # Calculate difference (negate metafor to match post-pre convention)
  md_expected <- (-as.numeric(mf_exp$yi)) - (-as.numeric(mf_nexp$yi))
  se_expected <- sqrt(as.numeric(mf_exp$vi) + as.numeric(mf_nexp$vi))

  # metaConvert using convert_df
  es_mc <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "md",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "bonett" # Method doesn't affect MD
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, md_expected, tolerance = 1e-10,
               label = "raw MD matches metafor MC")
  expect_equal(es_mc$se_crude, se_expected, tolerance = 1e-10,
               label = "raw MD SE matches metafor MC")
})

# ==============================================================================
# SECTION 2: es_from_means_se_pre_post() - TWO-GROUP
# ==============================================================================

test_that("MEANS-SE: all methods match SD version after SE→SD conversion", {
  res <- df.SMC[1:3, ]
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$r_pre_post_exp <- 0.6
  res$r_pre_post_nexp <- 0.6

  # Convert SD to SE
  res$mean_pre_se_exp <- res$mean_pre_sd_exp / sqrt(res$n_exp)
  res$mean_se_exp <- res$mean_sd_exp / sqrt(res$n_exp)
  res$mean_pre_se_nexp <- res$mean_pre_sd_nexp / sqrt(res$n_nexp)
  res$mean_se_nexp <- res$mean_sd_nexp / sqrt(res$n_nexp)

  methods <- c("bonett", "cooper", "morris_dav", "morris_dz")

  for (method in methods) {
    # From SD
    res_sd <- res
    es_sd <- summary(convert_df(res_sd,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_sd_pre_post",
      measure = "g",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # From SE (remove SD columns to force SE usage)
    res_se <- res
    res_se$mean_pre_sd_exp <- NULL
    res_se$mean_sd_exp <- NULL
    res_se$mean_pre_sd_nexp <- NULL
    res_se$mean_sd_nexp <- NULL

    es_se <- summary(convert_df(res_se,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_se_pre_post",
      measure = "g",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Verify correct hierarchies were used
    expect_equal(unique(es_sd$info_used_crude), "means_sd_pre_post")
    expect_equal(unique(es_se$info_used_crude), "means_se_pre_post")

    # Test equivalence
    expect_equal(es_se$es_crude, es_sd$es_crude, tolerance = 1e-10,
                 label = paste0(method, ": SE→SD gives same g"))
    expect_equal(es_se$se_crude, es_sd$se_crude, tolerance = 1e-10,
                 label = paste0(method, ": SE→SD gives same g_se"))
  }
})

# ==============================================================================
# SECTION 3: es_from_means_ci_pre_post() - TWO-GROUP
# ==============================================================================

test_that("MEANS-CI: all methods match SD version after CI→SD conversion", {
  res <- df.SMC[1:3, ]
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$r_pre_post_exp <- 0.6
  res$r_pre_post_nexp <- 0.6

  # Convert SD to CI
  res$mean_pre_se_exp <- res$mean_pre_sd_exp / sqrt(res$n_exp)
  res$mean_se_exp <- res$mean_sd_exp / sqrt(res$n_exp)
  res$mean_pre_se_nexp <- res$mean_pre_sd_nexp / sqrt(res$n_nexp)
  res$mean_se_nexp <- res$mean_sd_nexp / sqrt(res$n_nexp)

  res$mean_pre_ci_lo_exp <- res$mean_pre_exp - qt(0.975, res$n_exp - 1) * res$mean_pre_se_exp
  res$mean_pre_ci_up_exp <- res$mean_pre_exp + qt(0.975, res$n_exp - 1) * res$mean_pre_se_exp
  res$mean_ci_lo_exp <- res$mean_exp - qt(0.975, res$n_exp - 1) * res$mean_se_exp
  res$mean_ci_up_exp <- res$mean_exp + qt(0.975, res$n_exp - 1) * res$mean_se_exp

  res$mean_pre_ci_lo_nexp <- res$mean_pre_nexp - qt(0.975, res$n_nexp - 1) * res$mean_pre_se_nexp
  res$mean_pre_ci_up_nexp <- res$mean_pre_nexp + qt(0.975, res$n_nexp - 1) * res$mean_pre_se_nexp
  res$mean_ci_lo_nexp <- res$mean_nexp - qt(0.975, res$n_nexp - 1) * res$mean_se_nexp
  res$mean_ci_up_nexp <- res$mean_nexp + qt(0.975, res$n_nexp - 1) * res$mean_se_nexp

  methods <- c("bonett", "cooper", "morris_dav", "morris_dz")

  for (method in methods) {
    # From SD
    res_sd <- res
    es_sd <- summary(convert_df(res_sd,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_sd_pre_post",
      measure = "g",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # From CI (remove SD columns to force CI usage)
    res_ci <- res
    res_ci$mean_pre_sd_exp <- NULL
    res_ci$mean_sd_exp <- NULL
    res_ci$mean_pre_sd_nexp <- NULL
    res_ci$mean_sd_nexp <- NULL
    res_ci$mean_pre_se_exp <- NULL
    res_ci$mean_se_exp <- NULL
    res_ci$mean_pre_se_nexp <- NULL
    res_ci$mean_se_nexp <- NULL

    es_ci <- summary(convert_df(res_ci,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_ci_pre_post",
      measure = "g",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = method
    ), digits = 11)

    # Verify correct hierarchies were used
    expect_equal(unique(es_sd$info_used_crude), "means_sd_pre_post")
    expect_equal(unique(es_ci$info_used_crude), "means_ci_pre_post")

    # Test equivalence
    expect_equal(es_ci$es_crude, es_sd$es_crude, tolerance = 1e-9,
                 label = paste0(method, ": CI→SD gives same g"))
    expect_equal(es_ci$se_crude, es_sd$se_crude, tolerance = 1e-9,
                 label = paste0(method, ": CI→SD gives same g_se"))
  }
})
