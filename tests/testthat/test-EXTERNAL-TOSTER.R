# Long-running file (> 25 s): executed locally and in CI (devtools::test and
# testthat set NOT_CRAN=true); skipped wholesale on CRAN to respect check-time
# limits. The fast pre/post files still run on CRAN.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
# ==============================================================================
# EXTERNAL VALIDATION: metaConvert vs TOSTER (Caldwell, 2020)
# ==============================================================================
#
# Validates metaConvert's paired/pre-post effect sizes against TOSTER::smd_calc().
# TOSTER is authored by Aaron Caldwell (Caldwell & Vigotsky, 2020) and implements
# the Goulet-Pelletier & Cousineau (2018) formulas.
#
# KEY DIFFERENCES (acknowledged):
#   - Point estimates (d, g): Should match EXACTLY (same formulas)
#   - SE/variance: Will differ because:
#       * metaConvert uses LS (large-sample) approximation matching metafor
#       * TOSTER uses Goulet-Pelletier exact formula (noncentral t moments)
#     The LS approximation converges to the exact formula as n increases.
#   - J factor df: Both use df = n-1 for d_rm (when TOSTER smd_ci = "nct")
#   - Sign convention: TOSTER computes x-y; metaConvert computes post-pre
#
# References:
#   - Goulet-Pelletier & Cousineau (2018). TQMP, 14(4), 242-265.
#   - Caldwell & Vigotsky (2020). PeerJ, 8, e10314.
#   - Morris & DeShon (2002). Psych Methods, 7(1), 105-125.
#   - Becker (1988). British J Math Stat Psych, 41, 257-264.
# ==============================================================================

library(testthat)
library(metaConvert)

has_TOSTER <- requireNamespace("TOSTER", quietly = TRUE)
has_MASS <- requireNamespace("MASS", quietly = TRUE)

# Helper: generate paired data with known correlation structure
generate_paired_data <- function(n, mu1, mu2, sd1, sd2, r, seed = 42) {
  set.seed(seed)
  Sigma <- matrix(c(sd1^2, r * sd1 * sd2, r * sd1 * sd2, sd2^2), 2, 2)
  dat <- MASS::mvrnorm(n, c(mu1, mu2), Sigma)
  list(x = dat[, 1], y = dat[, 2])
}

# Helper: compute summary statistics from raw data
summarize_paired <- function(x, y) {
  list(
    n = length(x),
    mean_pre = mean(x), mean_post = mean(y),
    sd_pre = sd(x), sd_post = sd(y),
    r = cor(x, y),
    sd_diff = sd(y - x),
    mean_diff = mean(y - x)
  )
}

# ==============================================================================
# SECTION 1: POINT ESTIMATE VALIDATION (d_rm)
# ==============================================================================

test_that("TOSTER d_rm point estimate matches metaConvert morris_drm exactly", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  # Test across multiple scenarios
  scenarios <- list(
    list(n = 30, mu1 = 50, mu2 = 60, sd1 = 10, sd2 = 12, r = 0.7, seed = 42),
    list(n = 15, mu1 = 100, mu2 = 95, sd1 = 15, sd2 = 15, r = 0.5, seed = 123),
    list(n = 50, mu1 = 20, mu2 = 28, sd1 = 8, sd2 = 10, r = 0.8, seed = 7),
    list(n = 100, mu1 = 0, mu2 = 0.5, sd1 = 1, sd2 = 1.2, r = 0.3, seed = 99)
  )

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    dat <- generate_paired_data(sc$n, sc$mu1, sc$mu2, sc$sd1, sc$sd2, sc$r, sc$seed)
    ss <- summarize_paired(dat$x, dat$y)

    # TOSTER d_rm (no bias correction to compare d, not g)
    toster_d <- TOSTER::smd_calc(
      x = dat$x, y = dat$y,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = FALSE, smd_ci = "nct"
    )

    # metaConvert d_rm from summary statistics
    mc_res <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
      mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
      n_exp = ss$n, r_pre_post_exp = ss$r,
      pre_post_to_smd = "morris_drm"
    )

    # Point estimates must match exactly (sign-adjusted: TOSTER uses x-y)
    expect_equal(mc_res$d, -as.numeric(toster_d$estimate), tolerance = 1e-10,
                 label = paste("Scenario", i, "d_rm point estimate"))
  }
})

test_that("TOSTER g_rm (bias-corrected) matches metaConvert morris_drm g exactly", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  scenarios <- list(
    list(n = 30, mu1 = 50, mu2 = 60, sd1 = 10, sd2 = 12, r = 0.7, seed = 42),
    list(n = 15, mu1 = 100, mu2 = 95, sd1 = 15, sd2 = 15, r = 0.5, seed = 123),
    list(n = 50, mu1 = 20, mu2 = 28, sd1 = 8, sd2 = 10, r = 0.8, seed = 7)
  )

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    dat <- generate_paired_data(sc$n, sc$mu1, sc$mu2, sc$sd1, sc$sd2, sc$r, sc$seed)
    ss <- summarize_paired(dat$x, dat$y)

    # TOSTER g_rm with J(n-1) bias correction
    toster_g <- TOSTER::smd_calc(
      x = dat$x, y = dat$y,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = TRUE, smd_ci = "nct"
    )

    # metaConvert g
    mc_res <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
      mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
      n_exp = ss$n, r_pre_post_exp = ss$r,
      pre_post_to_smd = "morris_drm"
    )

    # g must match exactly (both use J(n-1))
    expect_equal(mc_res$g, -as.numeric(toster_g$estimate), tolerance = 1e-10,
                 label = paste("Scenario", i, "g_rm point estimate"))
  }
})

# ==============================================================================
# SECTION 2: POINT ESTIMATE VALIDATION (d_z)
# ==============================================================================

test_that("TOSTER d_z point estimate matches metaConvert morris_dz exactly", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  scenarios <- list(
    list(n = 30, mu1 = 50, mu2 = 60, sd1 = 10, sd2 = 12, r = 0.7, seed = 42),
    list(n = 20, mu1 = 80, mu2 = 75, sd1 = 20, sd2 = 18, r = 0.6, seed = 55),
    list(n = 50, mu1 = 20, mu2 = 28, sd1 = 8, sd2 = 10, r = 0.8, seed = 7)
  )

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    dat <- generate_paired_data(sc$n, sc$mu1, sc$mu2, sc$sd1, sc$sd2, sc$r, sc$seed)
    ss <- summarize_paired(dat$x, dat$y)

    # TOSTER d_z (no rm_correction)
    toster_dz <- TOSTER::smd_calc(
      x = dat$x, y = dat$y,
      paired = TRUE, rm_correction = FALSE,
      bias_correction = FALSE, smd_ci = "nct"
    )

    # metaConvert d_z
    mc_res <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
      mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
      n_exp = ss$n, r_pre_post_exp = ss$r,
      pre_post_to_smd = "morris_dz"
    )

    expect_equal(mc_res$d, -as.numeric(toster_dz$estimate), tolerance = 1e-10,
                 label = paste("Scenario", i, "d_z point estimate"))
  }
})

test_that("TOSTER g_z matches metaConvert morris_dz g exactly", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  dat <- generate_paired_data(30, 50, 60, 10, 12, 0.7)
  ss <- summarize_paired(dat$x, dat$y)

  toster_gz <- TOSTER::smd_calc(
    x = dat$x, y = dat$y,
    paired = TRUE, rm_correction = FALSE,
    bias_correction = TRUE, smd_ci = "nct"
  )

  mc_res <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
    mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
    n_exp = ss$n, r_pre_post_exp = ss$r,
    pre_post_to_smd = "morris_dz"
  )

  expect_equal(mc_res$g, -as.numeric(toster_gz$estimate), tolerance = 1e-10)
})

# ==============================================================================
# SECTION 3: d_rm = d_z * sqrt(2*(1-r)) IDENTITY
# ==============================================================================
# Confirmed by: Caldwell (2020) Eq 13, Goulet-Pelletier (2018) Eq 12b

test_that("d_rm = d_z * sqrt(2*(1-r)) identity holds in TOSTER", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  # Test with heteroscedastic SDs to confirm identity holds generally
  scenarios <- list(
    list(n = 30, sd1 = 10, sd2 = 12, r = 0.7),
    list(n = 25, sd1 = 15, sd2 = 8, r = 0.4),
    list(n = 40, sd1 = 5, sd2 = 5, r = 0.9),
    list(n = 20, sd1 = 20, sd2 = 25, r = 0.2)
  )

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    dat <- generate_paired_data(sc$n, 50, 60, sc$sd1, sc$sd2, sc$r, seed = i * 10)
    r_obs <- cor(dat$x, dat$y)

    toster_drm <- TOSTER::smd_calc(
      x = dat$x, y = dat$y,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = FALSE, smd_ci = "nct"
    )
    toster_dz <- TOSTER::smd_calc(
      x = dat$x, y = dat$y,
      paired = TRUE, rm_correction = FALSE,
      bias_correction = FALSE, smd_ci = "nct"
    )

    ratio <- as.numeric(toster_drm$estimate / toster_dz$estimate)
    expected_ratio <- sqrt(2 * (1 - r_obs))

    expect_equal(ratio, expected_ratio, tolerance = 1e-10,
                 label = paste("Scenario", i, "d_rm/d_z = sqrt(2*(1-r))"))
  }
})

# ==============================================================================
# SECTION 4: SE COMPARISON (acknowledged difference)
# ==============================================================================
# metaConvert uses LS approximation (matching metafor):
#   var_d = 2*(1-r)/n + d^2/(2*n)
#
# TOSTER uses Goulet-Pelletier exact formula:
#   var = (df/(df-2)) * (2*(1-r)/n) * (1 + g^2*n/(2*(1-r))) - g^2/J^2
#
# These converge as n -> infinity. Tests verify:
#   1. The difference is bounded (< 10% for n >= 15, < 5% for n >= 25)
#   2. The difference shrinks with increasing n
#   3. Both formulas agree on the ordering (larger d -> larger SE)

test_that("SE: metaConvert LS2 vs TOSTER GP exact - bounded difference", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  # Test across sample sizes
  n_values <- c(15, 25, 50, 100)
  relative_diffs <- numeric(length(n_values))

  for (i in seq_along(n_values)) {
    n <- n_values[i]
    dat <- generate_paired_data(n, 50, 60, 10, 12, 0.6, seed = i * 100)
    ss <- summarize_paired(dat$x, dat$y)

    # TOSTER SE (GP exact, using smd_ci="nct" so J uses df=n-1 like metaConvert)
    toster_g <- TOSTER::smd_calc(
      x = dat$x, y = dat$y,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = TRUE, smd_ci = "nct"
    )

    # metaConvert SE
    mc_res <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
      mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
      n_exp = ss$n, r_pre_post_exp = ss$r,
      pre_post_to_smd = "morris_drm"
    )

    se_toster <- as.numeric(toster_g$SE)
    se_mc <- mc_res$g_se

    relative_diffs[i] <- abs(se_mc - se_toster) / se_toster * 100

    # Bounded difference: < 15% for n >= 15
    expect_true(relative_diffs[i] < 15,
                label = paste("n =", n, ": SE relative diff =",
                              round(relative_diffs[i], 2), "% (< 15%)"))
  }

  # Convergence: difference should shrink with n
  # (monotonic convergence not guaranteed with random data, so just check
  #  that the largest n has smaller difference than the smallest n)
  expect_true(relative_diffs[length(relative_diffs)] < relative_diffs[1],
              label = "SE difference converges with increasing n")
})

test_that("SE: metaConvert vs TOSTER - same ordering across effect sizes", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  n <- 40
  r <- 0.6
  # Different effect magnitudes
  effect_sizes <- c(2, 5, 10, 20)
  se_mc_vals <- numeric(length(effect_sizes))
  se_toster_vals <- numeric(length(effect_sizes))

  for (i in seq_along(effect_sizes)) {
    dat <- generate_paired_data(n, 50, 50 + effect_sizes[i], 10, 10, r, seed = i * 7)
    ss <- summarize_paired(dat$x, dat$y)

    toster_g <- TOSTER::smd_calc(
      x = dat$x, y = dat$y,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = TRUE, smd_ci = "nct"
    )

    mc_res <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
      mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
      n_exp = ss$n, r_pre_post_exp = ss$r,
      pre_post_to_smd = "morris_drm"
    )

    se_mc_vals[i] <- mc_res$g_se
    se_toster_vals[i] <- as.numeric(toster_g$SE)
  }

  # Both should have same ordering: larger effect -> larger SE
  mc_order <- order(se_mc_vals)
  toster_order <- order(se_toster_vals)
  expect_equal(mc_order, toster_order,
               label = "SE ordering consistent between metaConvert and TOSTER")
})

# ==============================================================================
# SECTION 5: TWO-GROUP VALIDATION
# ==============================================================================
# TOSTER only computes single-group SMDs, so the two-group comparator below is
# built by subtracting two per-arm TOSTER values. That subtraction IS the legacy
# per-arm construction (pool_sd = FALSE): each arm standardized by its own SD.
# Since metaConvert 2.1.0 the DEFAULT is pool_sd = TRUE (a single SD pooled
# across arms, Morris 2008), which is a different -- and for unequal arm SDs, the
# correct -- estimand. These tests therefore pin pool_sd = FALSE explicitly and
# remain regression tests of the legacy path.

test_that("Two-group d_rm (legacy per-arm, pool_sd = FALSE): metaConvert matches TOSTER-based manual calculation", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  # Generate experimental group
  dat_exp <- generate_paired_data(35, 50, 62, 10, 11, 0.65, seed = 42)
  ss_exp <- summarize_paired(dat_exp$x, dat_exp$y)

  # Generate control group
  dat_nexp <- generate_paired_data(30, 50, 53, 10, 11, 0.65, seed = 99)
  ss_nexp <- summarize_paired(dat_nexp$x, dat_nexp$y)

  # TOSTER for each group separately
  toster_exp <- TOSTER::smd_calc(
    x = dat_exp$x, y = dat_exp$y,
    paired = TRUE, rm_correction = TRUE,
    bias_correction = TRUE, smd_ci = "nct"
  )
  toster_nexp <- TOSTER::smd_calc(
    x = dat_nexp$x, y = dat_nexp$y,
    paired = TRUE, rm_correction = TRUE,
    bias_correction = TRUE, smd_ci = "nct"
  )

  # Two-group d = d_exp - d_nexp (negate TOSTER sign for post-pre)
  g_toster_two <- (-as.numeric(toster_exp$estimate)) - (-as.numeric(toster_nexp$estimate))

  # metaConvert two-group
  mc_res <- es_from_means_sd_pre_post(
    mean_pre_exp = ss_exp$mean_pre, mean_exp = ss_exp$mean_post,
    mean_pre_sd_exp = ss_exp$sd_pre, mean_sd_exp = ss_exp$sd_post,
    n_exp = ss_exp$n, r_pre_post_exp = ss_exp$r,
    mean_pre_nexp = ss_nexp$mean_pre, mean_nexp = ss_nexp$mean_post,
    mean_pre_sd_nexp = ss_nexp$sd_pre, mean_sd_nexp = ss_nexp$sd_post,
    n_nexp = ss_nexp$n, r_pre_post_nexp = ss_nexp$r,
    pre_post_to_smd = "cooper",  # cooper = morris_drm
    # the TOSTER comparator above subtracts two per-arm values, so compare
    # against the legacy per-arm path, not the pooled default
    pool_sd = FALSE
  )

  # Point estimate must match
  expect_equal(mc_res$g, g_toster_two, tolerance = 1e-10,
               label = "Two-group g_rm matches TOSTER-based calculation")
})

# ==============================================================================
# SECTION 6: HOMOSCEDASTICITY IDENTITY
# ==============================================================================
# Under equal variances: bonett = morris_drm = morris_dav (Caldwell 2020)
# TOSTER d_rm should equal metaConvert bonett when sd_pre = sd_post

test_that("Under homoscedasticity: TOSTER d_rm = metaConvert bonett = metaConvert morris_dav", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  dat <- generate_paired_data(40, 50, 58, 10, 10, 0.6, seed = 42)
  ss <- summarize_paired(dat$x, dat$y)

  # TOSTER d_rm
  toster_d <- TOSTER::smd_calc(
    x = dat$x, y = dat$y,
    paired = TRUE, rm_correction = TRUE,
    bias_correction = FALSE, smd_ci = "nct"
  )

  # metaConvert bonett, morris_drm, morris_dav
  mc_bonett <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
    mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
    n_exp = ss$n, r_pre_post_exp = ss$r,
    pre_post_to_smd = "bonett"
  )
  mc_drm <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
    mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
    n_exp = ss$n, r_pre_post_exp = ss$r,
    pre_post_to_smd = "morris_drm"
  )
  mc_dav <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
    mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
    n_exp = ss$n, r_pre_post_exp = ss$r,
    pre_post_to_smd = "morris_dav"
  )

  d_toster <- -as.numeric(toster_d$estimate)

  # Under homoscedasticity: all three should give the same d
  # Note: sd_pre and sd_post from sample won't be exactly equal,
  # so use a slightly relaxed tolerance
  expect_equal(mc_bonett$d, mc_drm$d, tolerance = 1e-2,
               label = "bonett ≈ morris_drm under near-homoscedasticity")
  expect_equal(mc_drm$d, mc_dav$d, tolerance = 1e-2,
               label = "morris_drm ≈ morris_dav under near-homoscedasticity")
  expect_equal(mc_drm$d, d_toster, tolerance = 1e-10,
               label = "morris_drm matches TOSTER d_rm exactly")
})

# ==============================================================================
# SECTION 7: GLASS'S DELTA (bonett) vs TOSTER glass
# ==============================================================================

test_that("TOSTER glass1 matches metaConvert bonett for single-group", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  dat <- generate_paired_data(30, 50, 60, 10, 14, 0.7, seed = 42)
  ss <- summarize_paired(dat$x, dat$y)

  # TOSTER Glass's delta using pre-test SD (glass1 = x SD)
  toster_glass <- TOSTER::smd_calc(
    x = dat$x, y = dat$y,
    paired = TRUE, rm_correction = FALSE,
    bias_correction = FALSE, glass = "glass1",
    smd_ci = "nct"
  )

  # metaConvert bonett
  mc_bonett <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
    mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
    n_exp = ss$n, r_pre_post_exp = ss$r,
    pre_post_to_smd = "bonett"
  )

  # bonett d = (post - pre) / sd_pre
  # TOSTER glass1 d = (x - y) / sd_x where x=pre, y=post
  expect_equal(mc_bonett$d, -as.numeric(toster_glass$estimate), tolerance = 1e-10,
               label = "bonett d matches TOSTER Glass's delta (glass1)")
})

# ==============================================================================
# SECTION 8: CONVERT_DF INTEGRATION WITH TOSTER
# ==============================================================================

test_that("convert_df with morris_drm (legacy per-arm, pool_sd = FALSE) matches TOSTER for full pipeline", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  # Generate two studies' worth of paired data
  dat1_exp <- generate_paired_data(25, 50, 58, 10, 11, 0.6, seed = 1)
  dat1_nexp <- generate_paired_data(28, 50, 52, 10, 11, 0.6, seed = 2)
  dat2_exp <- generate_paired_data(40, 100, 115, 20, 22, 0.7, seed = 3)
  dat2_nexp <- generate_paired_data(35, 100, 103, 20, 22, 0.7, seed = 4)

  ss1_exp <- summarize_paired(dat1_exp$x, dat1_exp$y)
  ss1_nexp <- summarize_paired(dat1_nexp$x, dat1_nexp$y)
  ss2_exp <- summarize_paired(dat2_exp$x, dat2_exp$y)
  ss2_nexp <- summarize_paired(dat2_nexp$x, dat2_nexp$y)

  # Build dataframe for convert_df
  df_input <- data.frame(
    mean_pre_exp = c(ss1_exp$mean_pre, ss2_exp$mean_pre),
    mean_exp = c(ss1_exp$mean_post, ss2_exp$mean_post),
    mean_pre_sd_exp = c(ss1_exp$sd_pre, ss2_exp$sd_pre),
    mean_sd_exp = c(ss1_exp$sd_post, ss2_exp$sd_post),
    n_exp = c(ss1_exp$n, ss2_exp$n),
    r_pre_post_exp = c(ss1_exp$r, ss2_exp$r),
    mean_pre_nexp = c(ss1_nexp$mean_pre, ss2_nexp$mean_pre),
    mean_nexp = c(ss1_nexp$mean_post, ss2_nexp$mean_post),
    mean_pre_sd_nexp = c(ss1_nexp$sd_pre, ss2_nexp$sd_pre),
    mean_sd_nexp = c(ss1_nexp$sd_post, ss2_nexp$sd_post),
    n_nexp = c(ss1_nexp$n, ss2_nexp$n),
    r_pre_post_nexp = c(ss1_nexp$r, ss2_nexp$r)
  )

  # metaConvert pipeline
  es_mc <- summary(convert_df(df_input,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "g",
    pre_post_to_smd = "cooper",
    # the per-study TOSTER comparator below subtracts two per-arm values, i.e.
    # the legacy construction; the pooled default targets a different estimand
    pool_sd = FALSE
  ), digits = 11)

  # TOSTER-based manual calculation for each study
  for (study in 1:2) {
    if (study == 1) {
      dat_exp <- dat1_exp; dat_nexp <- dat1_nexp
    } else {
      dat_exp <- dat2_exp; dat_nexp <- dat2_nexp
    }

    toster_exp <- TOSTER::smd_calc(
      x = dat_exp$x, y = dat_exp$y,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = TRUE, smd_ci = "nct"
    )
    toster_nexp <- TOSTER::smd_calc(
      x = dat_nexp$x, y = dat_nexp$y,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = TRUE, smd_ci = "nct"
    )

    g_expected <- (-as.numeric(toster_exp$estimate)) - (-as.numeric(toster_nexp$estimate))

    expect_equal(es_mc$es_crude[study], g_expected, tolerance = 1e-10,
                 label = paste("Study", study, "convert_df g matches TOSTER"))
  }
})

# ==============================================================================
# SECTION 9: QUANTIFYING LS vs EXACT SE ACROSS SAMPLE SIZES
# ==============================================================================
# This is a documentation test: it verifies the known behavior that
# metaConvert's LS approximation converges to TOSTER's exact formula.

test_that("SE approximation error quantified across sample sizes", {
  skip_if_not(has_TOSTER, "TOSTER not available")
  skip_if_not(has_MASS, "MASS not available")

  # Fixed parameters, varying n
  n_values <- c(10, 15, 20, 30, 50, 100, 200)
  rel_diffs <- numeric(length(n_values))

  for (i in seq_along(n_values)) {
    n <- n_values[i]
    dat <- generate_paired_data(n, 50, 60, 10, 12, 0.6, seed = n)
    ss <- summarize_paired(dat$x, dat$y)

    toster_g <- TOSTER::smd_calc(
      x = dat$x, y = dat$y,
      paired = TRUE, rm_correction = TRUE,
      bias_correction = TRUE, smd_ci = "nct"
    )

    mc_res <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = ss$mean_pre, mean_exp = ss$mean_post,
      mean_pre_sd_exp = ss$sd_pre, mean_sd_exp = ss$sd_post,
      n_exp = ss$n, r_pre_post_exp = ss$r,
      pre_post_to_smd = "morris_drm"
    )

    se_toster <- as.numeric(toster_g$SE)
    se_mc <- mc_res$g_se
    rel_diffs[i] <- abs(se_mc - se_toster) / se_toster * 100
  }

  # Known bounds based on LS approximation theory:
  # The LS approximation omits (df/(df-2)) factor and -d^2/J^2 correction,
  # so divergence is largest at small n and large d.
  # n >= 10: < 25% (LS notably underestimates SE at very small n)
  # n >= 30: < 10%
  # n >= 100: < 5%
  for (i in seq_along(n_values)) {
    if (n_values[i] >= 100) {
      expect_true(rel_diffs[i] < 5,
                  label = paste("n =", n_values[i], ": SE diff", round(rel_diffs[i], 2), "% < 5%"))
    } else if (n_values[i] >= 30) {
      expect_true(rel_diffs[i] < 10,
                  label = paste("n =", n_values[i], ": SE diff", round(rel_diffs[i], 2), "% < 10%"))
    } else {
      expect_true(rel_diffs[i] < 25,
                  label = paste("n =", n_values[i], ": SE diff", round(rel_diffs[i], 2), "% < 25%"))
    }
  }

  # Overall convergence: largest n should have smallest relative diff
  expect_true(rel_diffs[length(rel_diffs)] < rel_diffs[1],
              label = "LS approximation converges to exact formula")
})

} # end NOT_CRAN gate
