# =============================================================================
# External Validation Tests: PROPORTION SINGLE-GROUP Functions
# =============================================================================
#
# Purpose: Comprehensive external validation of proportion functions in
#          R/es_from_PROP.R against metafor::escalc
#          (promoted from tests_save/checked/ — P26/P27, audit digest [25]).
#
# Functions tested:
#   1. es_from_prop_single_group()
#   2. es_from_prop_single_group_counts()
#
# External validation sources:
#   - metafor::escalc (PLO = logit, PFT = Freeman-Tukey, PR = raw proportion)
#
# Tolerance standards:
#   - 1e-10: Logit ES + SE against metafor PLO (exact formula match)
#   - 1e-10: Raw proportion ES + SE against metafor PR (exact formula match)
#   - 1e-10: Freeman-Tukey ES + SE against metafor PFT (exact match, both use
#            1/sqrt(4n+2))
#   - 1e-10: BOUNDARY (p = 0 or 1) SEs against metafor PR/PLO (P27: the
#            archived file asserted only boundary point estimates, which is
#            exactly where the pre-P5 SE deviated by sqrt((n+1)/n); the SEs
#            bit-match metafor now that the continuity-corrected denominator
#            n + 1 flows into the SE)
#
# =============================================================================

library(testthat)
library(metaConvert)

skip_if_not_installed("metafor")
library(metafor)

# =============================================================================
# SECTION 1: es_from_prop_single_group() — Direct Function vs metafor
# =============================================================================

# Test data: 5 studies with varied proportions and sample sizes
props <- c(0.30, 0.15, 0.05, 0.50, 0.75)
ns    <- c(100, 200, 300, 150, 80)
cases <- round(props * ns)

# -----------------------------------------------------------------------------
# Test 1.1: Logit (PLO) — ES + SE against metafor
# -----------------------------------------------------------------------------
test_that("PROP-SG: logit matches metafor PLO (ES + SE)", {
  mfr <- escalc(xi = cases, ni = ns, measure = "PLO")

  mc <- es_from_prop_single_group(prop = props, n_sample = ns, prop_to_es = "logit")

  expect_equal(mc$prop, as.numeric(mfr$yi), tolerance = 1e-10)
  expect_equal(mc$prop_se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 1.2: Freeman-Tukey (PFT) — ES + SE against metafor
# -----------------------------------------------------------------------------
test_that("PROP-SG: Freeman-Tukey matches metafor PFT (ES + SE)", {
  mfr <- escalc(xi = cases, ni = ns, measure = "PFT")

  mc <- es_from_prop_single_group(prop = props, n_sample = ns, prop_to_es = "freeman_tukey")

  expect_equal(mc$prop, as.numeric(mfr$yi), tolerance = 1e-10)
  expect_equal(mc$prop_se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 1.3: Raw proportion (PR) — ES + SE against metafor
# -----------------------------------------------------------------------------
test_that("PROP-SG: raw proportion matches metafor PR (ES + SE)", {
  mfr <- escalc(xi = cases, ni = ns, measure = "PR")

  mc <- es_from_prop_single_group(prop = props, n_sample = ns, prop_to_es = "raw")

  expect_equal(mc$prop, as.numeric(mfr$yi), tolerance = 1e-10)
  expect_equal(mc$prop_se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
})

# =============================================================================
# SECTION 2: es_from_prop_single_group_counts() — Counts Input
# =============================================================================

# -----------------------------------------------------------------------------
# Test 2.1: Counts produce identical results to proportion input
# -----------------------------------------------------------------------------
test_that("PROP-COUNTS: counts input matches proportion input exactly", {
  mc_prop <- es_from_prop_single_group(prop = cases / ns, n_sample = ns)
  mc_counts <- es_from_prop_single_group_counts(n_cases = cases, n_sample = ns)

  # All numeric columns must match (except info_used)
  numeric_cols <- setdiff(names(mc_prop), "info_used")
  for (col in numeric_cols) {
    expect_equal(mc_counts[[col]], mc_prop[[col]], tolerance = 1e-10,
                 label = paste("counts vs prop for column:", col))
  }

  # info_used differs as expected
  expect_true(all(mc_prop$info_used == "prop_single_group"))
  expect_true(all(mc_counts$info_used == "prop_single_group_counts"))
})

# -----------------------------------------------------------------------------
# Test 2.2: Counts logit + FT match metafor exactly
# -----------------------------------------------------------------------------
test_that("PROP-COUNTS: logit from counts matches metafor", {
  mfr_plo <- escalc(xi = cases, ni = ns, measure = "PLO")

  mc <- es_from_prop_single_group_counts(n_cases = cases, n_sample = ns, prop_to_es = "logit")

  expect_equal(mc$prop, as.numeric(mfr_plo$yi), tolerance = 1e-10)
  expect_equal(mc$prop_se, sqrt(as.numeric(mfr_plo$vi)), tolerance = 1e-10)
})

test_that("PROP-COUNTS: Freeman-Tukey from counts matches metafor", {
  mfr_pft <- escalc(xi = cases, ni = ns, measure = "PFT")

  mc <- es_from_prop_single_group_counts(n_cases = cases, n_sample = ns, prop_to_es = "freeman_tukey")

  expect_equal(mc$prop, as.numeric(mfr_pft$yi), tolerance = 1e-10)
})

# =============================================================================
# SECTION 3: Edge Cases — metafor Validation
# =============================================================================

# -----------------------------------------------------------------------------
# Test 3.1: Zero prevalence (0/n) — continuity correction
# -----------------------------------------------------------------------------
test_that("PROP-EDGE: zero prevalence matches metafor PLO + PFT + PR (ES AND SE)", {
  zero_cases <- c(0, 0, 0, 0)
  zero_ns <- c(20, 50, 100, 200)

  mfr_plo <- escalc(xi = zero_cases, ni = zero_ns, measure = "PLO")
  mfr_pft <- escalc(xi = zero_cases, ni = zero_ns, measure = "PFT")
  mfr_pr  <- escalc(xi = zero_cases, ni = zero_ns, measure = "PR")

  mc_logit <- es_from_prop_single_group_counts(n_cases = zero_cases, n_sample = zero_ns,
                                                prop_to_es = "logit")
  mc_ft <- es_from_prop_single_group_counts(n_cases = zero_cases, n_sample = zero_ns,
                                             prop_to_es = "freeman_tukey")
  mc_raw <- es_from_prop_single_group_counts(n_cases = zero_cases, n_sample = zero_ns,
                                              prop_to_es = "raw")

  # Logit point estimates match (both apply continuity correction)
  expect_equal(mc_logit$prop, as.numeric(mfr_plo$yi), tolerance = 1e-10)

  # Freeman-Tukey handles zeros naturally
  expect_equal(mc_ft$prop, as.numeric(mfr_pft$yi), tolerance = 1e-10)

  # Raw proportion: (0 + 0.5) / (n + 1) continuity correction
  expected_corrected <- 0.5 / (zero_ns + 1)
  expect_equal(mc_raw$prop, expected_corrected, tolerance = 1e-10)

  # P27: boundary SEs must ALSO bit-match the metafor oracle. Pre-P5, the SE
  # used the raw n with the corrected p and exceeded sqrt(vi) by exactly
  # sqrt((n+1)/n) — the one comparison the archived file skipped.
  expect_equal(mc_raw$prop_se, sqrt(as.numeric(mfr_pr$vi)), tolerance = 1e-10)
  expect_equal(mc_logit$prop_se, sqrt(as.numeric(mfr_plo$vi)), tolerance = 1e-10)
  expect_equal(mc_ft$prop_se, sqrt(as.numeric(mfr_pft$vi)), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 3.2: Full prevalence (n/n) — continuity correction
# -----------------------------------------------------------------------------
test_that("PROP-EDGE: full prevalence matches metafor PLO + PFT + PR (ES AND SE)", {
  full_ns <- c(20, 50, 100, 200)
  full_cases <- full_ns  # 100% prevalence (includes x = 20, n = 20)

  mfr_plo <- escalc(xi = full_cases, ni = full_ns, measure = "PLO")
  mfr_pft <- escalc(xi = full_cases, ni = full_ns, measure = "PFT")
  mfr_pr  <- escalc(xi = full_cases, ni = full_ns, measure = "PR")

  mc_logit <- es_from_prop_single_group_counts(n_cases = full_cases, n_sample = full_ns,
                                                prop_to_es = "logit")
  mc_ft <- es_from_prop_single_group_counts(n_cases = full_cases, n_sample = full_ns,
                                             prop_to_es = "freeman_tukey")
  mc_raw <- es_from_prop_single_group_counts(n_cases = full_cases, n_sample = full_ns,
                                              prop_to_es = "raw")

  expect_equal(mc_logit$prop, as.numeric(mfr_plo$yi), tolerance = 1e-10)
  expect_equal(mc_ft$prop, as.numeric(mfr_pft$yi), tolerance = 1e-10)

  # Raw proportion: (n + 0.5) / (n + 1) continuity correction
  expected_corrected <- (full_cases + 0.5) / (full_ns + 1)
  expect_equal(mc_raw$prop, expected_corrected, tolerance = 1e-10)

  # P27: boundary SEs at p = 1 bit-match metafor (see Test 3.1)
  expect_equal(mc_raw$prop_se, sqrt(as.numeric(mfr_pr$vi)), tolerance = 1e-10)
  expect_equal(mc_logit$prop_se, sqrt(as.numeric(mfr_plo$vi)), tolerance = 1e-10)
  expect_equal(mc_ft$prop_se, sqrt(as.numeric(mfr_pft$vi)), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 3.3: Near-boundary proportions (small and large)
# -----------------------------------------------------------------------------
test_that("PROP-EDGE: near-boundary proportions match metafor", {
  near_props <- c(0.001, 0.01, 0.99, 0.999)
  near_ns <- c(1000, 500, 500, 1000)
  near_cases <- round(near_props * near_ns)

  mfr_plo <- escalc(xi = near_cases, ni = near_ns, measure = "PLO")
  mfr_pft <- escalc(xi = near_cases, ni = near_ns, measure = "PFT")

  mc_logit <- es_from_prop_single_group(prop = near_cases / near_ns, n_sample = near_ns,
                                         prop_to_es = "logit")
  mc_ft <- es_from_prop_single_group(prop = near_cases / near_ns, n_sample = near_ns,
                                      prop_to_es = "freeman_tukey")

  expect_equal(mc_logit$prop, as.numeric(mfr_plo$yi), tolerance = 1e-10)
  expect_equal(mc_ft$prop, as.numeric(mfr_pft$yi), tolerance = 1e-10)
})

# =============================================================================
# SECTION 4: Pipeline Integration (convert_df + summary)
# =============================================================================

# -----------------------------------------------------------------------------
# Test 4.1: Full pipeline end-to-end with default raw method
# -----------------------------------------------------------------------------
test_that("PROP-PIPELINE: convert_df + summary end-to-end (raw)", {
  dat <- data.frame(
    study_id = 1:5,
    prop = props,
    n_sample = ns
  )

  es <- summary(convert_df(dat, verbose = FALSE, measure = "prop"), digits = 11)

  expect_equal(unique(es$info_used_crude), "prop_single_group")
  expect_equal(es$es_crude, props, tolerance = 1e-10)

  expected_se <- sqrt(props * (1 - props) / ns)
  expect_equal(es$se_crude, expected_se, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 4.2: Full pipeline with logit method
# -----------------------------------------------------------------------------
test_that("PROP-PIPELINE: convert_df + summary with logit method", {
  dat <- data.frame(
    study_id = 1:5,
    prop = props,
    n_sample = ns
  )

  mfr <- escalc(xi = cases, ni = ns, measure = "PLO")

  es <- summary(convert_df(dat, verbose = FALSE, measure = "prop", prop_to_es = "logit"), digits = 11)

  expect_equal(es$es_crude, as.numeric(mfr$yi), tolerance = 1e-10)
  expect_equal(es$se_crude, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 4.3: Full pipeline with freeman_tukey method
# -----------------------------------------------------------------------------
test_that("PROP-PIPELINE: convert_df + summary with freeman_tukey method", {
  dat <- data.frame(
    study_id = 1:5,
    prop = props,
    n_sample = ns
  )

  mfr <- escalc(xi = cases, ni = ns, measure = "PFT")

  es <- summary(convert_df(dat, verbose = FALSE, measure = "prop", prop_to_es = "freeman_tukey"), digits = 11)

  expect_equal(es$es_crude, as.numeric(mfr$yi), tolerance = 1e-10)
  expect_equal(es$se_crude, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 4.4: Mixed input (prop + n_cases) — hierarchy attribution
# -----------------------------------------------------------------------------
test_that("PROP-PIPELINE: mixed input hierarchy attribution", {
  dat <- data.frame(
    study_id = 1:4,
    prop = c(0.30, NA, 0.50, NA),
    n_cases = c(NA, 40, NA, 15),
    n_sample = c(100, 200, 150, 300)
  )

  es <- summary(convert_df(dat, verbose = FALSE, measure = "prop"), digits = 11)

  expect_equal(es$info_used_crude,
               c("prop_single_group", "prop_single_group_counts",
                 "prop_single_group", "prop_single_group_counts"))
  expect_equal(es$es_crude, c(0.30, 40/200, 0.50, 15/300), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 4.5: reverse_prop through full pipeline
# -----------------------------------------------------------------------------
test_that("PROP-PIPELINE: reverse_prop flips proportion", {
  dat_fwd <- data.frame(
    study_id = 1:3,
    prop = c(0.30, 0.15, 0.75),
    n_sample = c(100, 200, 80)
  )
  dat_rev <- dat_fwd
  dat_rev$reverse_prop <- TRUE

  es_fwd <- summary(convert_df(dat_fwd, verbose = FALSE, measure = "prop"), digits = 11)
  es_rev <- summary(convert_df(dat_rev, verbose = FALSE, measure = "prop"), digits = 11)

  # Reversed proportion = 1 - original
  expect_equal(es_rev$es_crude, 1 - es_fwd$es_crude, tolerance = 1e-10)
})

# =============================================================================
# SECTION 5: SE Variance Validation against metafor
# =============================================================================

# -----------------------------------------------------------------------------
# Test 5.1: Logit SE matches metafor sqrt(vi) for PLO
# -----------------------------------------------------------------------------
test_that("PROP-SE: logit SE matches metafor PLO variance", {
  test_props <- c(0.10, 0.25, 0.50, 0.75, 0.90)
  test_ns    <- c(50, 100, 200, 300, 500)
  test_cases <- round(test_props * test_ns)

  mfr <- escalc(xi = test_cases, ni = test_ns, measure = "PLO")
  mc  <- es_from_prop_single_group(prop = test_cases / test_ns, n_sample = test_ns,
                                    prop_to_es = "logit")

  expect_equal(mc$prop_se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 5.2: Freeman-Tukey SE matches metafor PFT variance
# -----------------------------------------------------------------------------
test_that("PROP-SE: Freeman-Tukey SE matches metafor PFT variance", {
  test_ns <- c(10, 50, 100, 500, 5000)
  test_props <- rep(0.30, length(test_ns))
  test_cases <- round(test_props * test_ns)

  mfr <- escalc(xi = test_cases, ni = test_ns, measure = "PFT")
  mc <- es_from_prop_single_group(prop = test_cases / test_ns, n_sample = test_ns,
                                   prop_to_es = "freeman_tukey")

  expect_equal(mc$prop_se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)

  # Verify the formula directly
  expected_se <- 1 / sqrt(4 * test_ns + 2)
  expect_equal(mc$prop_se, expected_se, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Test 5.3: Raw proportion SE matches metafor sqrt(vi) for PR
# -----------------------------------------------------------------------------
test_that("PROP-SE: raw proportion SE matches metafor PR variance", {
  test_props <- c(0.10, 0.25, 0.50, 0.75, 0.90)
  test_ns    <- c(50, 100, 200, 300, 500)
  test_cases <- round(test_props * test_ns)

  mfr <- escalc(xi = test_cases, ni = test_ns, measure = "PR")
  mc  <- es_from_prop_single_group(prop = test_cases / test_ns, n_sample = test_ns,
                                    prop_to_es = "raw")

  expect_equal(mc$prop_se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
})

# =============================================================================
# SECTION 6: Extreme Sample Sizes with metafor Validation
# =============================================================================

test_that("PROP-EXTREME: very small and very large n match metafor", {
  extreme_cases <- c(3, 5, 50, 500, 5000)
  extreme_ns    <- c(10, 20, 200, 2000, 20000)

  mfr_plo <- escalc(xi = extreme_cases, ni = extreme_ns, measure = "PLO")
  mfr_pft <- escalc(xi = extreme_cases, ni = extreme_ns, measure = "PFT")
  mfr_pr  <- escalc(xi = extreme_cases, ni = extreme_ns, measure = "PR")

  mc_logit <- es_from_prop_single_group(prop = extreme_cases / extreme_ns,
                                         n_sample = extreme_ns, prop_to_es = "logit")
  mc_ft <- es_from_prop_single_group(prop = extreme_cases / extreme_ns,
                                      n_sample = extreme_ns, prop_to_es = "freeman_tukey")
  mc_raw <- es_from_prop_single_group(prop = extreme_cases / extreme_ns,
                                       n_sample = extreme_ns, prop_to_es = "raw")

  expect_equal(mc_logit$prop, as.numeric(mfr_plo$yi), tolerance = 1e-10)
  expect_equal(mc_ft$prop, as.numeric(mfr_pft$yi), tolerance = 1e-10)
  expect_equal(mc_raw$prop, as.numeric(mfr_pr$yi), tolerance = 1e-10)

  expect_equal(mc_logit$prop_se, sqrt(as.numeric(mfr_plo$vi)), tolerance = 1e-10)
  expect_equal(mc_raw$prop_se, sqrt(as.numeric(mfr_pr$vi)), tolerance = 1e-10)
})
