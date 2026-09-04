# smd_var = "hedges_olkin" on every pre/post call below (2.1.0). These tests pin
# agreement with metafor's SMCC/SMCR/SMCRH/SMCRPH variances and the Bonett (2008) /
# Morris & DeShon (2002) formulas, which are the Hedges-Olkin ("LS") form. Since 2.1.0
# the pre/post routes honour smd_var like the two-group routes, and the package
# default "borenstein" applies J^2 to the d-scale variance instead; the metafor form is
# reached with smd_var = "hedges_olkin", which is what these pins now request. The
# exception is the cooper / morris_drm blocks, whose Morris & DeShon (2002) oracle is
# the d-scale variance 2(1-r)/n + d^2/(2n): that is the J^2-outside ("borenstein")
# convention, so those calls request the default explicitly.

# Long-running file (> 25 s): executed locally and in CI (devtools::test and
# testthat set NOT_CRAN=true); skipped wholesale on CRAN to respect check-time
# limits. The fast pre/post files still run on CRAN.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
# ==============================================================================
# EXTERNAL VALIDATION: es_from_PAIRED_MEANS.R
# ==============================================================================
# Validates all functions in es_from_PAIRED_MEANS.R against external packages:
#   - metafor (SMCRH for bonett, SMCRPH for morris_dav, SMCC for morris_dz, MC for raw MD)
#   - TOSTER (rm_correction for cooper, rm_correction=FALSE for morris_dz)
#
# Tests cover:
#   1. es_from_means_sd_pre_post() - all methods (bonett, cooper, morris_dav, morris_dz)
#   2. es_from_means_se_pre_post() - all methods
#   3. es_from_means_ci_pre_post() - all methods
#
# Dataset: metaumbrella::df.SMC (real meta-analysis data with pre/post means)
#
# ------------------------------------------------------------------------------
# NOTE ON pool_sd (Sections 1 and 1b)
# ------------------------------------------------------------------------------
# Every STANDARDIZED comparator in Section 1 (tests 1.1-1.4) is built the same
# way: a within-group effect size is computed for EACH arm separately (each arm
# standardized by its OWN SD) and the two are then SUBTRACTED, adding the two
# arms' sampling variances (the arms are independent). metafor's SMCRH / SMCRPH /
# SMCC are single-group (pre-post) estimators, so an arm-wise call plus a
# subtraction is the only way to build a between-group value out of them.
#
# That construction IS the `pool_sd = FALSE` path -- which is the DEFAULT, and
# the long-standing behaviour of the package. It is Morris (2008) d_ppc1 (from
# Becker 1988), and it is exactly what metafor users do for this design; see
# https://www.metafor-project.org/doku.php/analyses:morris2008 where Viechtbauer
# computes escalc(measure = "SMCR") per arm and then `yi = yT - yC; vi = vT + vC`.
# It is "more broadly applicable" (Viechtbauer) precisely because it does NOT
# assume the two arms' true standardizing SDs are equal.
#
# `pool_sd = TRUE` is the OPT-IN alternative: the difference in mean change is
# divided by ONE SD pooled across arms = Morris (2008) d_ppc2 (eq. 8-9). Morris
# recommends it (it is more efficient) but it does assume equal true arm SDs. The
# two coincide when that assumption holds and target different estimands
# otherwise, so the package does not make the choice silently.
#
# Section 1 therefore passes `pool_sd = FALSE` explicitly (redundant with the
# default, but it documents which estimand the external comparator encodes).
# Section 1b passes `pool_sd = TRUE` explicitly to give the opt-in pooled path
# its own external coverage.
# (Test 1.5, raw MD vs metafor MC, is UNSTANDARDIZED and hence pool_sd-invariant;
# it deliberately keeps running on the default.)
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
test_that("MEANS-SD: bonett two-group matches metafor SMCRH (g + g_se) (per-arm d_ppc1 standardizer, pool_sd = FALSE, the default)", {
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

  # Calculate difference (negate metafor to match post-pre convention).
  # SMCRH standardizes each arm by ITS OWN baseline SD, so subtracting the two
  # arm-wise values (and adding their variances) replicates the per-arm d_ppc1
  # rule -> pool_sd = FALSE below (the default, passed explicitly for clarity).
  g_expected <- (-as.numeric(smc_exp$yi)) - (-as.numeric(smc_nexp$yi))
  se_expected <- sqrt(as.numeric(smc_exp$vi) + as.numeric(smc_nexp$vi))

  # metaConvert using convert_df
  es_mc <- summary(convert_df(smd_var = "hedges_olkin", res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "g",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "bonett",
    pool_sd = FALSE # the DEFAULT (per-arm d_ppc1) -- what the comparator encodes
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
test_that("MEANS-SD: cooper two-group matches Morris & DeShon (2002) (d + d_se) (per-arm d_ppc1 standardizer, pool_sd = FALSE, the default)", {
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

  # Manual Morris & DeShon (2002) calculation.
  # This is the PUBLISHED d_rm formula applied ARM-WISE (each arm standardized by
  # its own change SD) and then subtracted -- i.e. the per-arm d_ppc1 rule, which
  # is why the metaConvert call below passes pool_sd = FALSE (the default).
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
  es_mc <- summary(convert_df(smd_var = "borenstein", res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "d",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "cooper",
    pool_sd = FALSE # the DEFAULT (per-arm d_ppc1) -- what the comparator encodes
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, d_expected, tolerance = 1e-10,
               label = "cooper d matches Morris & DeShon (2002)")
  expect_equal(es_mc$se_crude, se_expected, tolerance = 1e-10,
               label = "cooper d_se matches Morris & DeShon (2002)")
})

## 1.3 MORRIS_DAV vs metafor SMCRPH ====
test_that("MEANS-SD: morris_dav two-group matches metafor SMCRPH (g + g_se) (per-arm d_ppc1 standardizer, pool_sd = FALSE, the default)", {
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

  # metafor SMCRPH (heteroscedasticity-robust average-SD standardizer = Bonett 2008
  # eq. 10) for the experimental arm. d_av now uses the robust SMCRPH variance, not
  # the homoscedastic SMCRP form.
  mf_exp <- metafor::escalc(
    measure = "SMCRPH",
    m1i = res$mean_pre_exp,
    m2i = res$mean_exp,
    sd1i = res$mean_pre_sd_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp
  )

  # metafor SMCRPH for control group
  mf_nexp <- metafor::escalc(
    measure = "SMCRPH",
    m1i = res$mean_pre_nexp,
    m2i = res$mean_nexp,
    sd1i = res$mean_pre_sd_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp
  )

  # Calculate difference (negate metafor to match post-pre convention).
  # SMCRPH standardizes each arm by ITS OWN average pre/post SD, so subtracting the
  # two arm-wise values (and adding their variances) replicates the per-arm d_ppc1
  # rule -> pool_sd = FALSE below (the default, passed explicitly for clarity).
  g_expected <- (-as.numeric(mf_exp$yi)) - (-as.numeric(mf_nexp$yi))
  se_expected <- sqrt(as.numeric(mf_exp$vi) + as.numeric(mf_nexp$vi))

  # metaConvert using convert_df
  es_mc <- summary(convert_df(smd_var = "hedges_olkin", res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "g",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dav",
    pool_sd = FALSE # the DEFAULT (per-arm d_ppc1) -- what the comparator encodes
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, g_expected, tolerance = 1e-6,
               label = "morris_dav g matches metafor SMCRPH")
  expect_equal(es_mc$se_crude, se_expected, tolerance = 1e-6,
               label = "morris_dav g_se matches metafor SMCRPH")
})

## 1.4 MORRIS_DZ vs metafor SMCC ====
test_that("MEANS-SD: morris_dz two-group matches metafor SMCC (d + d_se) (per-arm d_ppc1 standardizer, pool_sd = FALSE, the default)", {
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

  # This test previously hand-transcribed the SMCC formula, which made it a
  # restatement of metaConvert's own source rather than an external check. It now
  # calls metafor::escalc(measure = "SMCC") LIVE, for each arm, and subtracts:
  #   - correct = FALSE  -> yi is the UNCORRECTED d_z  (the point estimate)
  #   - correct = TRUE   -> yi is g_z, vi is var(g_z)  (the variance)
  # The Hedges factor J is recovered purely from metafor (the ratio of the two
  # yi's), so var(d) = var(g)/J^2 needs no formula transcribed from the package.
  # SMCC standardizes each arm by ITS OWN change SD, so the arm-wise subtraction
  # below is the per-arm d_ppc1 rule -> pool_sd = FALSE in the convert_df call
  # (the default, passed explicitly for clarity).
  smcc <- function(arm, correct) {
    metafor::escalc(
      measure = "SMCC",
      m1i  = res[[paste0("mean_pre_", arm)]],
      m2i  = res[[paste0("mean_", arm)]],
      sd1i = res[[paste0("mean_pre_sd_", arm)]],
      sd2i = res[[paste0("mean_sd_", arm)]],
      ni   = res[[paste0("n_", arm)]],
      ri   = res[[paste0("r_pre_post_", arm)]],
      correct = correct
    )
  }

  g_exp  <- smcc("exp", TRUE)
  d_exp  <- smcc("exp", FALSE)
  g_nexp <- smcc("nexp", TRUE)
  d_nexp <- smcc("nexp", FALSE)

  # Hedges' J, derived from metafor itself (g / d), not transcribed
  J_exp  <- as.numeric(g_exp$yi) / as.numeric(d_exp$yi)
  J_nexp <- as.numeric(g_nexp$yi) / as.numeric(d_nexp$yi)

  # Difference (negate metafor to match the post-pre convention)
  d_expected <- (-as.numeric(d_exp$yi)) - (-as.numeric(d_nexp$yi))
  se_expected <- sqrt(as.numeric(g_exp$vi) / J_exp^2 +
                      as.numeric(g_nexp$vi) / J_nexp^2)

  # metaConvert using convert_df
  es_mc <- summary(convert_df(smd_var = "hedges_olkin", res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "d",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dz",
    pool_sd = FALSE # the DEFAULT (per-arm d_ppc1) -- what the comparator encodes
  ), digits = 11)

  # Verify correct hierarchy was used
  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # Test effect sizes and SEs
  expect_equal(es_mc$es_crude, d_expected, tolerance = 1e-10,
               label = "morris_dz d matches metafor SMCC")
  expect_equal(es_mc$se_crude, se_expected, tolerance = 1e-10,
               label = "morris_dz d_se matches metafor SMCC")
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
  # NB: no pool_sd here -- the raw MD is unstandardized, so it is invariant to
  # pool_sd. This test therefore runs on the DEFAULT (pool_sd = FALSE) and pins
  # that the pooling option never leaks into the MD path.
  es_mc <- summary(convert_df(smd_var = "hedges_olkin", res,
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
# SECTION 1b: POOLED STANDARDIZER (the OPT-IN pool_sd = TRUE path)
# ==============================================================================
# Section 1 validates the DEFAULT per-arm path (d_ppc1). This section validates
# the OPT-IN pooled path (Morris 2008 d_ppc2), which callers must now request
# explicitly with pool_sd = TRUE.
#
# With pool_sd = TRUE + pre_post_to_smd = "morris_dz", the between-group SMD is
# the difference in mean CHANGE divided by a single change-SD pooled across arms.
# That is exactly a two-sample Hedges' g computed ON THE CHANGE SCORES, i.e.
# metafor::escalc(measure = "SMD") fed the per-arm change means and change SDs --
# giving a real external comparator for the pooled path (no formula from the
# package is transcribed here).
# ==============================================================================

test_that("MEANS-SD: pooled morris_dz (opt-in pool_sd = TRUE) matches metafor SMD on change scores", {
  skip_if_not_installed("metafor")

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
  r <- 0.55
  res$r_pre_post_exp <- r
  res$r_pre_post_nexp <- r

  # Per-arm change means and change SDs (definitional, not a package formula)
  change_exp <- res$mean_exp - res$mean_pre_exp
  change_nexp <- res$mean_nexp - res$mean_pre_nexp
  sd_change_exp <- sqrt(res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
                        2 * r * res$mean_pre_sd_exp * res$mean_sd_exp)
  sd_change_nexp <- sqrt(res$mean_pre_sd_nexp^2 + res$mean_sd_nexp^2 -
                         2 * r * res$mean_pre_sd_nexp * res$mean_sd_nexp)

  # External comparator: ordinary two-sample SMD on the change scores
  mf <- metafor::escalc(
    measure = "SMD",
    m1i = change_exp, m2i = change_nexp,
    sd1i = sd_change_exp, sd2i = sd_change_nexp,
    n1i = res$n_exp, n2i = res$n_nexp,
    vtype = "LS2"
  )

  # metaConvert on the OPT-IN pooled path (pool_sd = TRUE passed explicitly --
  # it is no longer the default, so this must be requested)
  es_mc <- summary(convert_df(smd_var = "hedges_olkin", res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "g",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE
  ), digits = 11)

  expect_equal(unique(es_mc$info_used_crude), "means_sd_pre_post")

  # POINT ESTIMATE: exact agreement with metafor
  expect_equal(es_mc$es_crude, as.numeric(mf$yi), tolerance = 1e-9,
               label = "pooled morris_dz g matches metafor SMD on change scores")

  # SE: metaConvert uses var(g) = J^2 * (N/(n1*n2) + g^2/(2*m)) with m = N - 2;
  # metafor's LS2 uses J^2*(1/n1 + 1/n2) + g^2/(2*N). The two differ ONLY in the
  # small-sample g^2 term (J^2 and m = N-2 vs N), which is a sub-1% effect -- so a
  # 2% relative check is the tight bound here, NOT a loosened tolerance. It is a
  # real regression guard: the pre-fix variance carried a spurious 2*(1-r_avg)
  # factor on the sampling term, which at r = 0.55 scales the SE by sqrt(0.9)
  # (~5% off) and at small r by up to ~41% -- both far outside this band.
  expect_equal(es_mc$se_crude, sqrt(as.numeric(mf$vi)), tolerance = 0.02,
               label = "pooled morris_dz g_se agrees with metafor SMD (LS2) within 2%")
})

test_that("MEANS-SD: pool_sd defaults to FALSE (per-arm d_ppc1, Becker 1988 / Morris d_ppc1) and is not equivalent to the pooled path", {
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

  run <- function(...) {
    summary(convert_df(smd_var = "hedges_olkin", res,
      verbose = FALSE,
      es_selected = "hierarchy",
      hierarchy = "means_sd_pre_post",
      measure = "g",
      smd_to_cor = "viechtbauer",
      pre_post_to_smd = "bonett",
      ...
    ), digits = 11)
  }

  es_default <- run()
  es_pooled <- run(pool_sd = TRUE)
  es_per_arm <- run(pool_sd = FALSE)

  # The default IS the per-arm standardizer (d_ppc1): each arm is standardized by
  # its OWN SD, the two within-arm SMDs are subtracted, and their sampling
  # variances are ADDED (the arms are independent). This is the long-standing
  # behaviour and the construction metafor users build by hand from SMCR/SMCRH.
  expect_equal(es_default$es_crude, es_per_arm$es_crude, tolerance = 1e-12)
  expect_equal(es_default$se_crude, es_per_arm$se_crude, tolerance = 1e-12)

  # ... and the pooled standardizer (d_ppc2) is a genuinely DIFFERENT estimand on
  # this fixture (the arms' SDs are unequal), so it must be opted into explicitly
  # and the pool_sd = TRUE in Section 1b is load-bearing, not cosmetic.
  expect_false(isTRUE(all.equal(es_default$es_crude, es_pooled$es_crude,
                                tolerance = 1e-6)))
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
    es_sd <- summary(convert_df(smd_var = "hedges_olkin", res_sd,
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

    es_se <- summary(convert_df(smd_var = "hedges_olkin", res_se,
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
    es_sd <- summary(convert_df(smd_var = "hedges_olkin", res_sd,
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

    es_ci <- summary(convert_df(smd_var = "hedges_olkin", res_ci,
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

} # end NOT_CRAN gate
