library(testthat)

# ==============================================================================
# Tests for the quality/plausibility flag system
# ==============================================================================

test_that(".default_flag_options returns expected structure", {
  opts <- metaConvert:::.default_flag_options()
  expect_type(opts, "list")
  expect_true("smd_max" %in% names(opts))
  expect_true("r_max" %in% names(opts))
  expect_true("log_or_max" %in% names(opts))
  expect_true("n_min" %in% names(opts))
  expect_true("iqr_mult" %in% names(opts))
  expect_true("dispersion_max_smd" %in% names(opts))
  expect_true("dispersion_max_r" %in% names(opts))
  expect_true("dispersion_max_logor" %in% names(opts))
  expect_true("overlap_min" %in% names(opts))
  expect_true("enable_cross_row" %in% names(opts))
  expect_true("sd_ratio_max" %in% names(opts))
  expect_equal(opts$smd_max, 3)
  expect_equal(opts$sd_ratio_max, 10)
  expect_equal(opts$overlap_min, 0.85)
  expect_equal(opts$dispersion_max_smd, 0.5)
  expect_equal(opts$dispersion_max_r, 0.15)
  expect_equal(opts$dispersion_max_logor, 1.0)
})


# ==============================================================================
# Category A: Numeric Integrity
# ==============================================================================

test_that("[A1] flags Inf/NaN ES values", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(Inf, -Inf, NaN, 0.5, NA),
    se = c(0.1, 0.1, 0.1, 0.1, NA),
    ci_lo = c(-0.1, -0.3, -0.1, 0.3, NA),
    ci_up = c(0.3, 0.1, 0.3, 0.7, NA)
  )
  expect_true(any(grepl("ES is non-finite or undefined", flags[[1]])))
  expect_true(any(grepl("ES is non-finite or undefined", flags[[2]])))
  expect_true(any(grepl("ES is non-finite or undefined", flags[[3]])))
  expect_false(any(grepl("ES is non-finite or undefined", flags[[4]])))
  expect_equal(length(flags[[5]]), 0)
})

test_that("[A2] flags Inf/NaN SE values", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(0.5, 0.5),
    se = c(Inf, NaN),
    ci_lo = c(0.3, 0.3),
    ci_up = c(0.7, 0.7)
  )
  expect_true(any(grepl("SE is non-finite or undefined", flags[[1]])))
  expect_true(any(grepl("SE is non-finite or undefined", flags[[2]])))
})

test_that("[A3] flags negative SE", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(0.5),
    se = c(-0.1),
    ci_lo = c(0.3),
    ci_up = c(0.7)
  )
  expect_true(any(grepl("Negative SE", flags[[1]])))
})

test_that("[A4] flags inverted CI", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(0.5),
    se = c(0.1),
    ci_lo = c(0.8),
    ci_up = c(0.2)
  )
  expect_true(any(grepl("Inverted CI: lower", flags[[1]])))
})

test_that("[A5] flags ES outside CI", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(2.0),
    se = c(0.1),
    ci_lo = c(0.3),
    ci_up = c(0.7)
  )
  expect_true(any(grepl("ES outside its CI", flags[[1]])))
})

test_that("No A-flags for normal values", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(0.5, NA),
    se = c(0.1, NA),
    ci_lo = c(0.3, NA),
    ci_up = c(0.7, NA)
  )
  expect_equal(length(flags[[1]]), 0)
  expect_equal(length(flags[[2]]), 0)
})


# ==============================================================================
# Category B: Bounds Violations
# ==============================================================================

test_that("[B1] flags correlation outside [-1, 1]", {
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(1.05, -1.2, 0.8),
    se = c(0.1, 0.1, 0.1),
    ci_lo = c(0.9, -1.4, 0.6),
    ci_up = c(1.2, -1.0, 1.0),
    measure = "r", exp = FALSE
  )
  expect_true(any(grepl("r outside", flags[[1]])))
  expect_true(any(grepl("r outside", flags[[2]])))
  expect_false(any(grepl("r outside", flags[[3]])))
})

test_that("[B3] flags NNT magnitude less than 1", {
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(0.5, -0.3, 5, 0),
    se = c(NA, NA, NA, NA),
    ci_lo = c(NA, NA, NA, NA),
    ci_up = c(NA, NA, NA, NA),
    measure = "nnt", exp = FALSE
  )
  expect_true(any(grepl("NNT magnitude < 1", flags[[1]])))
  expect_true(any(grepl("NNT magnitude < 1", flags[[2]])))
  expect_false(any(grepl("NNT magnitude < 1", flags[[3]])))
  expect_false(any(grepl("NNT magnitude < 1", flags[[4]])))  # NNT=0 is excluded
})

test_that("[B4] flags proportion outside [0, 1]", {
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(-0.1, 1.2, 0.5),
    se = c(0.1, 0.1, 0.1),
    ci_lo = c(-0.3, 1.0, 0.3),
    ci_up = c(0.1, 1.4, 0.7),
    measure = "prop", exp = FALSE
  )
  expect_true(any(grepl("Proportion outside", flags[[1]])))
  expect_true(any(grepl("Proportion outside", flags[[2]])))
  expect_false(any(grepl("Proportion outside", flags[[3]])))
})

test_that("[B4] does NOT flag logit-scale proportions outside [0, 1]", {
  # logit(0.3) = -0.847, logit(0.8) = 1.386 — valid on logit scale
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(-0.847, 1.386, -6.9),
    se = c(0.1, 0.1, 0.5),
    ci_lo = c(-1.0, 1.2, -7.9),
    ci_up = c(-0.6, 1.6, -5.9),
    measure = "prop", exp = FALSE,
    prop_to_es = "logit"
  )
  expect_false(any(grepl("Proportion outside", flags[[1]])))
  expect_false(any(grepl("Proportion outside", flags[[2]])))
  expect_false(any(grepl("Proportion outside", flags[[3]])))
})

test_that("[B4] does NOT flag Freeman-Tukey-scale proportions outside [0, 1]", {
  # FT(0.8, n=100) ~ 1.107 — valid on FT scale (range [0, pi/2])
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(1.107, 1.539),
    se = c(0.05, 0.016),
    ci_lo = c(1.0, 1.5),
    ci_up = c(1.2, 1.6),
    measure = "prop", exp = FALSE,
    prop_to_es = "freeman_tukey"
  )
  expect_false(any(grepl("Proportion outside", flags[[1]])))
  expect_false(any(grepl("Proportion outside", flags[[2]])))
})


# ==============================================================================
# Category C: Plausibility
# ==============================================================================

test_that("[C1] flags large SMD", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(4.0, -5.0, 0.5, 2.9),
    se = c(0.5, 0.5, 0.1, 0.3),
    measure = "d", opts = opts,
    n_sample = rep(NA, 4)
  )
  expect_true(any(grepl("Large SMD", flags[[1]])))
  expect_true(any(grepl("Large SMD", flags[[2]])))
  expect_false(any(grepl("Large SMD", flags[[3]])))
  expect_false(any(grepl("Large SMD", flags[[4]])))  # 2.9 < 3
})

test_that("[C1] respects custom smd_max threshold", {
  opts <- metaConvert:::.default_flag_options()
  opts$smd_max <- 5
  flags <- metaConvert:::.flag_plausibility(
    es = c(4.0, 6.0),
    se = c(0.5, 0.5),
    measure = "g", opts = opts,
    n_sample = rep(NA, 2)
  )
  expect_false(any(grepl("Large SMD", flags[[1]])))  # 4.0 < 5
  expect_true(any(grepl("Large SMD", flags[[2]])))   # 6.0 > 5
})

test_that("[C2] flags zero SE", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.5, 0.5),
    se = c(0.0, 0.1),
    measure = "d", opts = opts,
    n_sample = rep(NA, 2)
  )
  expect_true(any(grepl("SE is zero", flags[[1]])))
  expect_false(any(grepl("SE is zero", flags[[2]])))
})

test_that("[C3] removed — SE/ES ratio check no longer exists", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.1, 0.5),
    se = c(5.0, 0.3),
    measure = "d", opts = opts,
    n_sample = rep(NA, 2)
  )
  # C3 removed: genuine SE problems caught by C2 (SE=0) and F1/F2 (SE vs N)
  expect_false(any(grepl("SE very large", flags[[1]])))
  expect_false(any(grepl("SE very large", flags[[2]])))
})

test_that("[C4] flags high correlation", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.98, 0.80),
    se = c(0.01, 0.05),
    measure = "r", opts = opts,
    n_sample = rep(NA, 2)
  )
  expect_true(any(grepl("High correlation", flags[[1]])))
  expect_false(any(grepl("High correlation", flags[[2]])))
})

test_that("[C6] small sample is suppressed by default (informational flag)", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.5, 0.5, 0.5),
    se = c(0.1, 0.1, 0.1),
    measure = "d", opts = opts,
    n_sample = c(3, 10, 100)
  )
  expect_false(any(grepl("Small sample", flags[[1]])))
  expect_false(any(grepl("Small sample", flags[[2]])))
  expect_false(any(grepl("Small sample", flags[[3]])))
})

test_that("[C6] small sample fires when enable_informational = TRUE", {
  opts <- metaConvert:::.default_flag_options()
  opts$enable_informational <- TRUE
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.5, 0.5, 0.5),
    se = c(0.1, 0.1, 0.1),
    measure = "d", opts = opts,
    n_sample = c(3, 10, 100)
  )
  expect_true(any(grepl("Small sample", flags[[1]])))
  expect_false(any(grepl("Small sample", flags[[2]])))  # n=10, not < 10
  expect_false(any(grepl("Small sample", flags[[3]])))
})

test_that("[C9] flags extreme transformed proportion - logit scale", {
  opts <- metaConvert:::.default_flag_options()
  # logit(0.001) ~ -6.9, logit(0.999) ~ 6.9, logit(0.3) ~ -0.847
  flags <- metaConvert:::.flag_plausibility(
    es = c(-6.9, 6.9, -0.847),
    se = c(0.5, 0.5, 0.1),
    measure = "prop", opts = opts,
    n_sample = c(100, 100, 100),
    prop_to_es = "logit"
  )
  expect_true(any(grepl("Extreme transformed proportion", flags[[1]])))
  expect_true(any(grepl("Extreme transformed proportion", flags[[2]])))
  expect_false(any(grepl("Extreme transformed proportion", flags[[3]])))
})

test_that("[C9] flags extreme transformed proportion - Freeman-Tukey scale", {
  opts <- metaConvert:::.default_flag_options()
  # FT(0.001, n=1000) ~ 0.032, FT(0.999, n=1000) ~ 1.539, FT(0.3, n=100) ~ 0.580
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.032, 1.539, 0.580),
    se = c(0.016, 0.016, 0.050),
    measure = "prop", opts = opts,
    n_sample = c(1000, 1000, 100),
    prop_to_es = "freeman_tukey"
  )
  expect_true(any(grepl("Extreme transformed proportion", flags[[1]])))
  expect_true(any(grepl("Extreme transformed proportion", flags[[2]])))
  expect_false(any(grepl("Extreme transformed proportion", flags[[3]])))
})

test_that("[C9] does not flag raw proportions (B4 handles those)", {
  opts <- metaConvert:::.default_flag_options()
  # Raw proportion 0.01 is in [0, 1] — not on transformed scale
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.01, 0.99),
    se = c(0.01, 0.01),
    measure = "prop", opts = opts,
    n_sample = c(100, 100),
    prop_to_es = "raw"
  )
  expect_false(any(grepl("Extreme transformed proportion", flags[[1]])))
  expect_false(any(grepl("Extreme transformed proportion", flags[[2]])))
})


# ==============================================================================
# Category D: Cross-Row Outliers
# ==============================================================================

test_that("[D1] flags ES outlier", {
  opts <- metaConvert:::.default_flag_options()
  # All values near 0.5, one extreme outlier at 20
  es <- c(0.5, 0.4, 0.6, 0.5, 0.45, 20.0)
  se <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  flags <- metaConvert:::.flag_cross_row_outliers(es, se, opts)
  expect_true(any(grepl("ES outlier", flags[[6]])))
  expect_false(any(grepl("ES outlier", flags[[1]])))
})

test_that("[D2] flags SE outlier", {
  opts <- metaConvert:::.default_flag_options()
  # SE values with natural variance so IQR > 0, plus one extreme outlier
  es <- c(0.5, 0.4, 0.6, 0.5, 0.45, 0.5)
  se <- c(0.10, 0.15, 0.12, 0.20, 0.11, 50.0)
  flags <- metaConvert:::.flag_cross_row_outliers(es, se, opts)
  expect_true(any(grepl("SE outlier", flags[[6]])))
  expect_false(any(grepl("SE outlier", flags[[1]])))
})

test_that("Cross-row checks need >= 4 values", {
  opts <- metaConvert:::.default_flag_options()
  # Only 3 values: IQR method should not fire
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = c(0.5, 0.5, 100),
    se = c(0.1, 0.1, 0.1),
    opts
  )
  expect_equal(length(flags[[3]]), 0)
})


# ==============================================================================
# Category E: Internal Consistency
# ==============================================================================

test_that("[E1] flags high dispersion", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3, 3),
    dispersion = c(0.8, 0.2),
    overlap = c(0.5, 0.5),
    diff_mm = c(0.5, 0.5),
    opts = opts, measure = "d"
  )
  expect_true(any(grepl("High dispersion", flags[[1]])))
  expect_false(any(grepl("High dispersion", flags[[2]])))
})

test_that("[E2] flags zero CI overlap", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3, 3),
    dispersion = c(0.1, 0.1),
    overlap = c(0.0, 0.5),
    diff_mm = c(0.5, 0.5),
    opts = opts, measure = "d"
  )
  expect_true(any(grepl("Zero CI overlap", flags[[1]])))
  expect_false(any(grepl("Zero CI overlap", flags[[2]])))
})

test_that("[E2b] flags low CI overlap", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3, 3, 3),
    dispersion = c(0.1, 0.1, 0.1),
    overlap = c(0.5, 0.9, 0.83),
    diff_mm = c(0.5, 0.5, 0.5),
    opts = opts, measure = "d"
  )
  expect_true(any(grepl("Low CI overlap", flags[[1]])))   # 0.50 < 0.85
  expect_false(any(grepl("Low CI overlap", flags[[2]])))  # 0.90 >= 0.85
  expect_true(any(grepl("Low CI overlap", flags[[3]])))   # 0.83 < 0.85 (new threshold)
})

test_that("overlap_min_max formula returns intersection / union for nested CIs", {
  # Regression: the prior |max_lo - min_up| / |max_up - min_lo| formula
  # was correct for shifted CIs but wrong for nested CIs (could return > 1
  # or hide nesting at 100%). Verify the new formula reports proper I/U.
  #
  # Build a dataset where one method produces a wider CI than another at
  # the same point estimate, so the narrower CI (mostly) nests inside the
  # wider. Row 1: means_sd path (recomputed SE = 0.341) vs user_input path.
  # Since Step G the user CI is preserved verbatim ([-2.33, -1.22], centred
  # at -1.775 rather than rebuilt symmetric around -1.70), so the expected
  # intersection/union is 0.822 (was 0.831 with the old symmetric rebuild).
  df <- data.frame(
    n_exp = 82, n_nexp = 41,
    mean_exp = -9.4, mean_sd_exp = 2.1,
    mean_nexp = -7.7, mean_sd_nexp = 1.6,
    user_es_original_measure_crude = "md",
    user_es_target_measure_crude = "md",
    user_es_crude = -1.70,
    user_ci_lo_crude = -2.33,
    user_ci_up_crude = -1.22
  )
  res <- as.data.frame(summary(convert_df(df, measure = "md"), flags = TRUE))
  overlap <- as.numeric(res$overlap_min_max_crude[1])
  expect_gte(overlap, 0)
  expect_lte(overlap, 1)
  expect_equal(overlap, 0.822, tolerance = 0.01)
})

test_that("[E3] flags large min-max difference", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3, 3),
    dispersion = c(0.1, 0.1),
    overlap = c(0.9, 0.9),
    diff_mm = c(2.0, 0.5),
    opts = opts, measure = "d"
  )
  expect_true(any(grepl("Large min-max", flags[[1]])))   # 2.0 > 1.0
  expect_false(any(grepl("Large min-max", flags[[2]])))  # 0.5 < 1.0
})

test_that("E-checks skip rows with < 2 estimations", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(1, 0),
    dispersion = c(0.8, 0.8),
    overlap = c(0.0, 0.0),
    diff_mm = c(2.0, 2.0),
    opts = opts, measure = "d"
  )
  expect_equal(length(flags[[1]]), 0)
  expect_equal(length(flags[[2]]), 0)
})


# ==============================================================================
# Integration: summary.metaConvert with flags
# ==============================================================================

test_that("summary with split_adjusted=FALSE includes flags column", {
  res <- summary(convert_df(df.haza, measure = "g", verbose = FALSE,
                            split_adjusted = FALSE))
  expect_true("flags" %in% colnames(res))
  expect_type(res$flags, "character")
})

test_that("summary with split_adjusted=TRUE (default) includes flags_crude and flags_adjusted", {
  res <- summary(convert_df(df.haza, measure = "g", verbose = FALSE))
  expect_true("flags_crude" %in% colnames(res))
  expect_true("flags_adjusted" %in% colnames(res))
  expect_type(res$flags_crude, "character")
  expect_type(res$flags_adjusted, "character")
})

test_that("flags = FALSE omits flags columns", {
  res <- summary(convert_df(df.haza, measure = "g", verbose = FALSE,
                            split_adjusted = FALSE),
                 flags = FALSE)
  expect_false("flags" %in% colnames(res))

  res2 <- summary(convert_df(df.haza, measure = "g", verbose = FALSE),
                  flags = FALSE)
  expect_false("flags_crude" %in% colnames(res2))
  expect_false("flags_adjusted" %in% colnames(res2))
})

test_that("flag_options can be overridden in summary()", {
  obj <- convert_df(df.haza, measure = "g", verbose = FALSE, split_adjusted = FALSE)
  res_strict <- summary(obj, flag_options = list(smd_max = 0.01))
  res_loose <- summary(obj, flag_options = list(smd_max = 100))
  n_strict <- sum(grepl("Large SMD", res_strict$flags))
  n_loose <- sum(grepl("Large SMD", res_loose$flags))
  expect_true(n_strict >= n_loose)
})

test_that("flag_options stored in convert_df are used by summary", {
  obj <- convert_df(df.haza, measure = "g", verbose = FALSE,
                    split_adjusted = FALSE,
                    flag_options = list(smd_max = 0.01))
  res <- summary(obj)
  n_flagged <- sum(grepl("Large SMD", res$flags))
  expect_true(n_flagged > 0)
})

test_that("flags column contains empty strings for NA rows", {
  res <- summary(convert_df(df.haza, measure = "g", verbose = FALSE,
                            split_adjusted = FALSE))
  na_rows <- which(is.na(res$es))
  if (length(na_rows) > 0) {
    expect_true(all(res$flags[na_rows] == ""))
  }
})


# ==============================================================================
# New tests: Previously untested flags and edge cases
# ==============================================================================

test_that("[B2] flags non-positive OR/RR on natural scale (exp=TRUE)", {
  # When exp=TRUE, ES column contains exp(logor) which must be > 0.
  # es=0 corresponds to exp(-Inf)=0, es=-1 is impossible for exp().
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(0, -1, 2.5, NA),
    se = c(0.1, 0.1, 0.1, NA),
    ci_lo = c(NA, NA, NA, NA),
    ci_up = c(NA, NA, NA, NA),
    measure = "logor", exp = TRUE
  )
  expect_true(any(grepl("non-positive", flags[[1]])))   # es=0 flagged
  expect_true(any(grepl("non-positive", flags[[2]])))   # es=-1 flagged
  expect_false(any(grepl("non-positive", flags[[3]])))  # es=2.5 valid
  expect_equal(length(flags[[4]]), 0)                # NA skipped
})

test_that("[B2] does not flag when exp=FALSE (log scale)", {
  # On log scale, negative values are legitimate (OR < 1)
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(-0.5, 0.5),
    se = c(0.1, 0.1),
    ci_lo = c(NA, NA),
    ci_up = c(NA, NA),
    measure = "logor", exp = FALSE
  )
  expect_false(any(grepl("non-positive", flags[[1]])))
  expect_false(any(grepl("non-positive", flags[[2]])))
})

test_that("[B2b] flags non-positive OR/RR/IRR CI bounds on natural scale", {
  # ES is positive but the CI lower bound is <= 0 on the natural scale.
  # Row 1: ci_lo = 0  → flag
  # Row 2: ci_lo = -0.5 → flag (lower)
  # Row 3: ci_up = -0.1 → flag (upper)
  # Row 4: valid CI, no flag
  flags <- metaConvert:::.flag_bounds_violations(
    es    = c(2.0, 1.5, 1.2, 2.0),
    se    = c(0.1, 0.1, 0.1, 0.1),
    ci_lo = c(0,   -0.5, 0.5, 0.5),
    ci_up = c(4.0,  3.0, -0.1, 4.0),
    measure = "logor", exp = TRUE
  )
  expect_true(any(grepl("CI bound non-positive", flags[[1]])))
  expect_true(any(grepl("CI lower", flags[[2]])))
  expect_true(any(grepl("CI upper", flags[[3]])))
  expect_false(any(grepl("CI bound non-positive", flags[[4]])))
})

test_that("[B2b] does not flag on log scale (exp=FALSE)", {
  # On log scale, negative CI bounds are legitimate
  flags <- metaConvert:::.flag_bounds_violations(
    es    = c(0.5),
    se    = c(0.1),
    ci_lo = c(-0.5),
    ci_up = c(1.0),
    measure = "logor", exp = FALSE
  )
  expect_false(any(grepl("CI bound non-positive", flags[[1]])))
})

test_that("[B6] flags RD point estimate outside [-1, 1]", {
  flags <- metaConvert:::.flag_bounds_violations(
    es    = c(1.5, -1.2, 0.4, NA),
    se    = c(0.1, 0.1, 0.1, NA),
    ci_lo = c(1.4, -1.3, 0.3, NA),
    ci_up = c(1.6, -1.1, 0.5, NA),
    measure = "rd", exp = FALSE
  )
  expect_true(any(grepl("RD outside", flags[[1]])))
  expect_true(any(grepl("RD outside", flags[[2]])))
  expect_false(any(grepl("RD outside", flags[[3]])))
})

test_that("[B6b] flags RD CI bound outside [-1, 1] even when point estimate is valid", {
  # Row 1: rd = 0.5 (valid), ci_lo = -1.2 (invalid) → flag
  # Row 2: rd = 0.7 (valid), ci_up = 1.4 (invalid) → flag
  # Row 3: rd = 0.4 with valid CI → no flag
  flags <- metaConvert:::.flag_bounds_violations(
    es    = c(0.5, 0.7, 0.4),
    se    = c(0.1, 0.1, 0.1),
    ci_lo = c(-1.2, 0.5, 0.3),
    ci_up = c( 0.9, 1.4, 0.5),
    measure = "rd", exp = FALSE
  )
  expect_true(any(grepl("CI bound outside \\[-1, 1\\]", flags[[1]])))
  expect_true(any(grepl("CI lower", flags[[1]])))
  expect_true(any(grepl("CI bound outside \\[-1, 1\\]", flags[[2]])))
  expect_true(any(grepl("CI upper", flags[[2]])))
  expect_false(any(grepl("CI bound outside", flags[[3]])))
})

test_that("[C5] flags large log OR/RR magnitude", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(6.0, -6.0, 4.0),
    se = c(0.5, 0.5, 0.5),
    measure = "logor", opts = opts,
    n_sample = rep(NA, 3)
  )
  expect_true(any(grepl("Large log", flags[[1]])))   # |6.0| > 5
  expect_true(any(grepl("Large log", flags[[2]])))   # |-6.0| > 5
  expect_false(any(grepl("Large log", flags[[3]])))  # |4.0| < 5
})

test_that("[C1] exact boundary: ES=3.0 with smd_max=3 is NOT flagged (strict >)", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(3.0, 3.001),
    se = c(0.1, 0.1),
    measure = "d", opts = opts,
    n_sample = rep(NA, 2)
  )
  expect_false(any(grepl("Large SMD", flags[[1]])))  # 3.0 is NOT > 3
  expect_true(any(grepl("Large SMD", flags[[2]])))   # 3.001 is > 3
})

test_that("[C4] exact boundary: r=0.95 with r_max=0.95 is NOT flagged (strict >)", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.95, 0.951),
    se = c(0.01, 0.01),
    measure = "r", opts = opts,
    n_sample = rep(NA, 2)
  )
  expect_false(any(grepl("High correlation", flags[[1]])))  # 0.95 is NOT > 0.95
  expect_true(any(grepl("High correlation", flags[[2]])))   # 0.951 is > 0.95
})

test_that("[E2b] exact boundary: overlap=overlap_min is NOT flagged (strict <)", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3, 3),
    dispersion = c(0.1, 0.1),
    overlap = c(opts$overlap_min, opts$overlap_min - 0.01),
    diff_mm = c(0.5, 0.5),
    opts = opts, measure = "d"
  )
  expect_false(any(grepl("Low CI overlap", flags[[1]])))  # exactly at threshold: not flagged
  expect_true(any(grepl("Low CI overlap", flags[[2]])))   # just below threshold: flagged
})

test_that("enable_cross_row=FALSE disables D1/D2", {
  opts <- metaConvert:::.default_flag_options()
  opts$enable_cross_row <- FALSE
  # Even with an extreme outlier, no D flags should fire
  es <- c(0.5, 0.4, 0.6, 0.5, 0.45, 20.0)
  se <- c(0.10, 0.15, 0.12, 0.20, 0.11, 50.0)
  flags <- metaConvert:::.flag_cross_row_outliers(es, se, opts)
  for (i in seq_along(flags)) {
    expect_equal(length(flags[[i]]), 0)
  }
})

test_that("[C3-removed] ES=0 with large SE does NOT trigger", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_plausibility(
    es = c(0.0),
    se = c(999.0),
    measure = "d", opts = opts,
    n_sample = rep(NA, 1)
  )
  # C3 removed entirely
  expect_false(any(grepl("SE very large", flags[[1]])))
})

test_that("[D1] n=4 with IQR > 0 activates cross-row detection", {
  opts <- metaConvert:::.default_flag_options()
  # Exactly 4 values: IQR method should activate
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = c(0.5, 0.5, 0.6, 100.0),
    se = c(0.1, 0.1, 0.1, 0.1),
    opts
  )
  expect_true(any(grepl("ES outlier", flags[[4]])))  # 100.0 is extreme
})

test_that("[D1] n=4 with IQR=0 does NOT activate cross-row detection", {
  opts <- metaConvert:::.default_flag_options()
  # 4 identical values + 1 outlier: but first 4 are identical so IQR=0
  # Wait: with 5 values total the IQR is computed across all 5, not just the first 4
  # Use 4 identical values plus the outlier contributes to IQR
  # Actually, let's test pure IQR=0 case: all values identical
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = c(0.5, 0.5, 0.5, 0.5),
    se = c(0.1, 0.1, 0.1, 0.1),
    opts
  )
  for (i in 1:4) {
    expect_false(any(grepl("ES outlier", flags[[i]])))
  }
})

test_that("Multiple flag categories can trigger on single row", {
  # Row with Inf ES should trigger A1 but not crash other checks
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(Inf),
    se = c(-0.1),
    ci_lo = c(0.8),
    ci_up = c(0.2)
  )
  # Should get both A1 (Inf ES) and A3 (negative SE) and A4 (inverted CI)
  expect_true(any(grepl("ES is non-finite or undefined", flags[[1]])))
  expect_true(any(grepl("Negative SE", flags[[1]])))
  expect_true(any(grepl("Inverted CI: lower", flags[[1]])))
})

test_that("[E1] uses measure-specific dispersion thresholds", {
  opts <- metaConvert:::.default_flag_options()
  # For measure="r", dispersion_max_r=0.15
  flags_r <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3, 3),
    dispersion = c(0.20, 0.10),
    overlap = c(0.9, 0.9),
    diff_mm = c(0.1, 0.1),
    opts = opts, measure = "r"
  )
  expect_true(any(grepl("High dispersion", flags_r[[1]])))   # 0.20 > 0.15
  expect_false(any(grepl("High dispersion", flags_r[[2]])))  # 0.10 < 0.15

  # For measure="d", dispersion_max_smd=0.5
  flags_d <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3, 3),
    dispersion = c(0.20, 0.60),
    overlap = c(0.9, 0.9),
    diff_mm = c(0.1, 0.1),
    opts = opts, measure = "d"
  )
  expect_false(any(grepl("High dispersion", flags_d[[1]])))  # 0.20 < 0.5
  expect_true(any(grepl("High dispersion", flags_d[[2]])))   # 0.60 > 0.5

  # For measure="logor", dispersion_max_logor=1.0
  flags_or <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3, 3),
    dispersion = c(1.2, 0.8),
    overlap = c(0.9, 0.9),
    diff_mm = c(0.1, 0.1),
    opts = opts, measure = "logor"
  )
  expect_true(any(grepl("High dispersion", flags_or[[1]])))   # 1.2 > 1.0
  expect_false(any(grepl("High dispersion", flags_or[[2]])))  # 0.8 < 1.0
})


# ==============================================================================
# Sign reversal: systematic vs partial
# ==============================================================================

test_that("systematic sign reversal (all rows) suppresses all E-flags", {
  opts <- metaConvert:::.default_flag_options()
  opts$ignore_sign_reversal <- TRUE
  # All 3 rows have min ≈ -max (systematic reversal)
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(2, 2, 2),
    dispersion = c(0.8, 0.8, 0.8),
    overlap = c(0, 0, 0),
    diff_mm = c(1.0, 1.0, 1.0),
    opts = opts, measure = "d",
    min_es = c(-0.5, -0.3, -0.8),
    max_es = c(0.5, 0.3, 0.8)
  )
  # All E-flags should be suppressed
  for (i in 1:3) {
    expect_equal(length(flags[[i]]), 0,
      info = paste("Row", i, "should have no E-flags"))
  }
})

test_that("partial sign reversal (some rows) keeps all E-flags active", {
  opts <- metaConvert:::.default_flag_options()
  opts$ignore_sign_reversal <- TRUE
  # Row 1 and 2 have reversal, row 3 does NOT
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(2, 2, 2),
    dispersion = c(0.8, 0.8, 0.8),
    overlap = c(0, 0, 0),
    diff_mm = c(1.0, 1.0, 1.0),
    opts = opts, measure = "d",
    min_es = c(-0.5, -0.3, 0.1),
    max_es = c(0.5, 0.3, 0.8)
  )
  # E-flags should fire on ALL rows (including reversed ones)
  for (i in 1:3) {
    expect_true(length(flags[[i]]) > 0,
      info = paste("Row", i, "should have E-flags"))
  }
})

test_that("ignore_sign_reversal = FALSE keeps E-flags even when all reversed", {
  opts <- metaConvert:::.default_flag_options()
  opts$ignore_sign_reversal <- FALSE
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(2, 2),
    dispersion = c(0.8, 0.8),
    overlap = c(0, 0),
    diff_mm = c(1.0, 1.0),
    opts = opts, measure = "d",
    min_es = c(-0.5, -0.3),
    max_es = c(0.5, 0.3)
  )
  for (i in 1:2) {
    expect_true(length(flags[[i]]) > 0,
      info = paste("Row", i, "should have E-flags when ignore_sign_reversal=FALSE"))
  }
})

test_that("single-method rows are ignored when checking systematic reversal", {
  opts <- metaConvert:::.default_flag_options()
  opts$ignore_sign_reversal <- TRUE
  # Row 1: multi-method, reversed. Row 2: single method. Row 3: multi-method, reversed.
  # All multi-method rows are reversed → systematic
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(2, 1, 2),
    dispersion = c(0.8, NA, 0.8),
    overlap = c(0, NA, 0),
    diff_mm = c(1.0, NA, 1.0),
    opts = opts, measure = "d",
    min_es = c(-0.5, NA, -0.3),
    max_es = c(0.5, NA, 0.3)
  )
  # All suppressed (single-method row has nothing to flag anyway)
  for (i in 1:3) {
    expect_equal(length(flags[[i]]), 0)
  }
})

# ==============================================================================
# Category V: Input Data Validation (.validate_input_data)
# ==============================================================================

test_that("[V1] flags and sets negative n_exp to NA", {
  dat <- data.frame(n_exp = c(-5, 10, NA), n_nexp = c(10, 10, 10))
  result <- metaConvert:::.validate_input_data(dat)
  expect_true(is.na(result$data$n_exp[1]))
  expect_equal(result$data$n_exp[2], 10)
  expect_true(is.na(result$data$n_exp[3]))  # was already NA
  expect_true(grepl("Negative input.*n_exp", result$issues[1]))
  expect_equal(result$issues[2], "")
  expect_equal(result$issues[3], "")
})

test_that("[V1] flags negative SD values", {
  dat <- data.frame(mean_sd_exp = c(-1.5, 2.0), mean_sd_nexp = c(1.0, -0.5))
  result <- metaConvert:::.validate_input_data(dat)
  expect_true(is.na(result$data$mean_sd_exp[1]))
  expect_equal(result$data$mean_sd_exp[2], 2.0)
  expect_equal(result$data$mean_sd_nexp[1], 1.0)
  expect_true(is.na(result$data$mean_sd_nexp[2]))
  expect_true(grepl("Negative input.*mean_sd_exp", result$issues[1]))
  expect_true(grepl("Negative input.*mean_sd_nexp", result$issues[2]))
})

test_that("[V1] flags negative SE values", {
  dat <- data.frame(user_se_crude = c(-0.1, 0.2))
  result <- metaConvert:::.validate_input_data(dat)
  expect_true(is.na(result$data$user_se_crude[1]))
  expect_equal(result$data$user_se_crude[2], 0.2)
  expect_true(grepl("Negative input.*user_se_crude", result$issues[1]))
})

test_that("[V2] flags inverted CI (lo > up) and sets triplet to NA", {
  dat <- data.frame(
    mean_exp = c(5.0, 5.0),
    mean_ci_lo_exp = c(6.0, 4.0),  # Row 1: lo=6 > up=4 → inverted
    mean_ci_up_exp = c(4.0, 6.0)
  )
  result <- metaConvert:::.validate_input_data(dat)
  # Row 1: all three columns set to NA
  expect_true(is.na(result$data$mean_exp[1]))
  expect_true(is.na(result$data$mean_ci_lo_exp[1]))
  expect_true(is.na(result$data$mean_ci_up_exp[1]))
  # Row 2: untouched
  expect_equal(result$data$mean_exp[2], 5.0)
  expect_equal(result$data$mean_ci_lo_exp[2], 4.0)
  expect_equal(result$data$mean_ci_up_exp[2], 6.0)
  expect_true(grepl("Inverted CI for.*mean_exp", result$issues[1]))
  expect_equal(result$issues[2], "")
})

test_that("[V3] flags value outside CI and sets triplet to NA", {
  dat <- data.frame(
    md = c(10.0, 5.0),
    md_ci_lo = c(1.0, 4.0),
    md_ci_up = c(3.0, 6.0)  # Row 1: val=10 > up=3 → outside
  )
  result <- metaConvert:::.validate_input_data(dat)
  expect_true(is.na(result$data$md[1]))
  expect_true(is.na(result$data$md_ci_lo[1]))
  expect_true(is.na(result$data$md_ci_up[1]))
  expect_equal(result$data$md[2], 5.0)
  expect_true(grepl("Value outside CI.*md", result$issues[1]))
  expect_equal(result$issues[2], "")
})

test_that("[V4] flags asymmetric CI but preserves data (warn only)", {
  # md=5, CI=[2, 9]: lower_dist=3, upper_dist=4, avg=3.5, asym=28.6% > 10%
  dat <- data.frame(
    md = c(5.0, 5.0),
    md_ci_lo = c(2.0, 4.0),
    md_ci_up = c(9.0, 6.0)  # Row 1 asymmetric, Row 2 symmetric
  )
  result <- metaConvert:::.validate_input_data(dat, max_asymmetry = 10)
  # Data is preserved (warn only)
  expect_equal(result$data$md[1], 5.0)
  expect_equal(result$data$md_ci_lo[1], 2.0)
  expect_equal(result$data$md_ci_up[1], 9.0)
  expect_equal(result$data$md[2], 5.0)
  # Flag is recorded for row 1 only
  expect_true(grepl("Asymmetric CI.*md", result$issues[1]))
  expect_equal(result$issues[2], "")
})

test_that("[V4] respects custom max_asymmetry threshold", {
  # md=5, CI=[4, 7]: lower_dist=1, upper_dist=2, avg=1.5, asym=66.7%
  dat <- data.frame(md = 5.0, md_ci_lo = 4.0, md_ci_up = 7.0)
  # With very high threshold, no flag
  result_high <- metaConvert:::.validate_input_data(dat, max_asymmetry = 100)
  expect_equal(result_high$issues[1], "")
  expect_equal(result_high$data$md[1], 5.0)
  # With low threshold, flags it but preserves data
  result_low <- metaConvert:::.validate_input_data(dat, max_asymmetry = 10)
  expect_true(grepl("Asymmetric CI", result_low$issues[1]))
  expect_equal(result_low$data$md[1], 5.0)
})

test_that("[V4] skips zero-width CI (no division by zero)", {
  # md=5, CI=[5, 5]: avg_dist=0, should be skipped
  dat <- data.frame(md = 5.0, md_ci_lo = 5.0, md_ci_up = 5.0)
  result <- metaConvert:::.validate_input_data(dat, max_asymmetry = 10)
  expect_equal(result$issues[1], "")
  expect_equal(result$data$md[1], 5.0)
})

test_that("[V4] log-transforms exp-scale CIs (OR) before asymmetry check", {
  # OR=2, CI=[1.2, 4.0]
  # On natural scale: asymmetric (0.8 below, 2.0 above)
  # On log scale: log(2)=0.693, log(1.2)=0.182, log(4)=1.386
  #   lower_dist=0.511, upper_dist=0.693, avg=0.602, asym=30.2% > 10%
  dat <- data.frame(or = 2.0, or_ci_lo = 1.2, or_ci_up = 4.0)
  result <- metaConvert:::.validate_input_data(dat, max_asymmetry = 10)
  expect_true(grepl("Asymmetric CI.*or", result$issues[1]))
  # Data preserved
  expect_equal(result$data$or[1], 2.0)
})

test_that("[V4] no flag for symmetric exp-scale CI (OR)", {
  # OR=2, CI symmetrically constructed on log scale
  # log(2)=0.693, se=0.3 => CI = exp(0.693 +/- 1.96*0.3) = [exp(0.105), exp(1.281)]
  lo <- exp(log(2) - 1.96 * 0.3)
  up <- exp(log(2) + 1.96 * 0.3)
  dat <- data.frame(or = 2.0, or_ci_lo = lo, or_ci_up = up)
  result <- metaConvert:::.validate_input_data(dat, max_asymmetry = 10)
  expect_equal(result$issues[1], "")
})

test_that("Multiple V-flags on one row are semicolon-separated", {
  dat <- data.frame(
    n_exp = c(-5),
    mean_sd_exp = c(-1.0),
    mean_sd_nexp = c(1.0)
  )
  result <- metaConvert:::.validate_input_data(dat)
  expect_true(grepl("Negative input.*n_exp", result$issues[1]))
  expect_true(grepl("Negative input.*mean_sd_exp", result$issues[1]))
  expect_true(grepl("; ", result$issues[1]))
})

test_that(".validate_input_data skips columns not present in data", {
  dat <- data.frame(some_column = c(-5, 10))
  result <- metaConvert:::.validate_input_data(dat)
  # No known columns, so no issues
  expect_equal(result$issues[1], "")
  expect_equal(result$issues[2], "")
})

test_that(".validate_input_data emits console message for bad rows", {
  dat <- data.frame(n_exp = c(-5, 10))
  expect_message(
    metaConvert:::.validate_input_data(dat),
    "Input data validation"
  )
})

test_that(".validate_input_data is silent when all data is valid", {
  dat <- data.frame(n_exp = c(10, 20), mean_sd_exp = c(1.0, 2.0))
  expect_silent(metaConvert:::.validate_input_data(dat))
})

test_that(".validate_input_data respects verbose=FALSE", {
  dat <- data.frame(n_exp = c(-5, 10))
  expect_silent(metaConvert:::.validate_input_data(dat, verbose = FALSE))
  # Issues still recorded even with verbose=FALSE
  result <- suppressMessages(metaConvert:::.validate_input_data(dat, verbose = FALSE))
  expect_true(grepl("Negative input", result$issues[1]))
})

test_that(".positive_columns returns a non-empty character vector", {
  cols <- metaConvert:::.positive_columns()
  expect_type(cols, "character")
  expect_true(length(cols) > 20)
  expect_true("n_exp" %in% cols)
  expect_true("user_se_crude" %in% cols)
})

test_that(".ci_triplets returns a list of triplets with val/lo/up/scale", {
  triplets <- metaConvert:::.ci_triplets()
  expect_type(triplets, "list")
  expect_true(length(triplets) > 10)
  for (tri in triplets) {
    expect_true(all(c("val", "lo", "up", "scale") %in% names(tri)))
    expect_true(tri$scale %in% c("additive", "exp", "user"))
  }
  # Check that exp-scale triplets include OR and RR
  exp_vals <- vapply(triplets[vapply(triplets, function(t) t$scale == "exp", logical(1))],
                     function(t) t$val, character(1))
  expect_true("or" %in% exp_vals)
  expect_true("rr" %in% exp_vals)
})


# ==============================================================================
# Integration: centralized validation flows through to summary() flags
# ==============================================================================

test_that("convert_df does not crash when user provides negative SE", {
  dat <- data.frame(
    measure = c("g", "g"),
    user_es_measure_crude = c("my_es", "my_es"),
    user_es_crude = c(0.5, 0.8),
    user_se_crude = c(-0.1, 0.2)
  )
  obj <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = FALSE)
  # Validation issues stored on the object
  input_val <- attr(obj, "input_validation")
  expect_true(grepl("Negative input", input_val[1]))
  expect_equal(input_val[2], "")

  res <- summary(obj)
  # Row 2 should still produce a valid ES
  expect_false(is.na(res$es[2]))
})

test_that("V-flags appear in flags column of summary()", {
  dat <- data.frame(
    measure = c("g", "g"),
    user_es_measure_crude = c("my_es", "my_es"),
    user_es_crude = c(0.5, 0.8),
    user_se_crude = c(-0.1, 0.2)
  )
  suppressWarnings(
    obj <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = FALSE)
  )
  res <- summary(obj, flags = TRUE)
  expect_true("flags" %in% colnames(res))
  # Row 1 should have V1 flag
  expect_true(grepl("Negative input", res$flags[1]))
  # Row 2 should NOT have V1 flag
  expect_false(grepl("Negative input", res$flags[2]))
})

test_that("V-flags appear first in the flags string (before A-E flags)", {
  dat <- data.frame(
    measure = c("g", "g"),
    user_es_measure_crude = c("my_es", "my_es"),
    user_es_crude = c(0.5, 0.8),
    user_se_crude = c(-0.1, 0.2)
  )
  suppressWarnings(
    obj <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = FALSE)
  )
  res <- summary(obj, flags = TRUE)
  # Row 1 should have at least the V1 flag
  v_flag_row <- res$flags[1]
  expect_true(nchar(v_flag_row) > 0)
  # The first flag should be a V-flag (validation), carrying one of the four
  # family prefixes: [INVALID] (data replaced with NA), [UNUSUAL] (data
  # preserved, plausibility warning), [DISCORDANT] (cross-method or cross-row
  # inconsistency), or [INFO] (soft informational). A negative user_se_crude
  # is V1 (hard error — value set to NA), so it must be [INVALID].
  first_flag <- strsplit(v_flag_row, "; ")[[1]][1]
  expect_true(grepl("^\\[(INVALID|UNUSUAL|DISCORDANT|INFO)\\] (Negative input|Inverted CI for|Value outside CI|Asymmetric CI|Out-of-range)", first_flag))
})

test_that("flags=FALSE still works with centralized validation", {
  dat <- data.frame(
    measure = c("g", "g"),
    user_es_measure_crude = c("my_es", "my_es"),
    user_es_crude = c(0.5, 0.8),
    user_se_crude = c(-0.1, 0.2)
  )
  suppressWarnings(
    obj <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = FALSE)
  )
  res <- summary(obj, flags = FALSE)
  expect_false("flags" %in% colnames(res))
})

test_that("V-flags work with split_adjusted=TRUE", {
  dat <- data.frame(
    measure = c("g", "g"),
    user_es_measure_crude = c("my_es", "my_es"),
    user_es_crude = c(0.5, 0.8),
    user_se_crude = c(-0.1, 0.2)
  )
  suppressWarnings(
    obj <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = TRUE)
  )
  res <- summary(obj, flags = TRUE)
  expect_true("flags_crude" %in% colnames(res))
  expect_true("flags_adjusted" %in% colnames(res))
  # Row 1 crude should have V1 flag (user_se_crude is crude-scoped)
  expect_true(grepl("Negative input", res$flags_crude[1]))
})

test_that("V-flags are scoped: crude-only flags don't appear in flags_adjusted", {
  dat <- data.frame(
    md = 2, md_ci_lo = 1, md_ci_up = 4,
    n_exp = 20, n_nexp = 20
  )
  obj <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = TRUE)
  res <- summary(obj, flags = TRUE)
  # md is a crude column — V4 should only be in flags_crude
  expect_true(grepl("Asymmetric CI", res$flags_crude[1]))
  expect_false(grepl("Asymmetric CI", res$flags_adjusted[1]))
})

test_that("V-flags for shared columns (n_exp) appear in both crude and adjusted", {
  dat <- data.frame(
    measure = c("g"),
    user_es_measure_crude = c("my_es"),
    user_es_crude = c(0.5),
    user_se_crude = c(0.1),
    n_exp = c(-5),
    n_nexp = c(20)
  )
  obj <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = TRUE)
  res <- summary(obj, flags = TRUE)
  # n_exp is shared — V1 should appear in both
  expect_true(grepl("Negative input.*n_exp", res$flags_crude[1]))
  expect_true(grepl("Negative input.*n_exp", res$flags_adjusted[1]))
})

test_that("V4 asymmetry flags flow through to summary()", {
  dat <- data.frame(
    measure = c("g", "g"),
    user_es_measure_crude = c("my_es", "my_es"),
    user_es_crude = c(0.5, 0.8),
    user_ci_lo_crude = c(0.1, 0.6),
    user_ci_up_crude = c(1.5, 1.0)  # Row 1: asym=66.7% (lo_dist=0.4, up_dist=1.0)
  )
  suppressWarnings(
    obj <- convert_df(dat, measure = "g", verbose = FALSE,
                      split_adjusted = FALSE, max_asymmetry = 10)
  )
  res <- summary(obj, flags = TRUE)
  expect_true(grepl("Asymmetric CI", res$flags[1]))
  expect_false(grepl("Asymmetric CI", res$flags[2]))
  # V4 is warn-only: ES should still be computed for row 1
  expect_false(is.na(res$es[1]))
})

test_that("V4 asymmetry skips log-transform for ratio-scale user_es (OR/RR)", {
  # OR = 2.0, CI [1.2, 3.5] — naturally asymmetric on additive scale
  # but symmetric on log scale: log(2.0)=0.693, log(1.2)=0.182, log(3.5)=1.253
  # log-scale distances: 0.693-0.182=0.511, 1.253-0.693=0.560 → ~9.2% asymmetry
  dat <- data.frame(
    user_es_original_measure_crude = c("or", "or"),
    user_es_crude = c(2.0, 1.5),
    user_ci_lo_crude = c(1.2, 1.1),
    user_ci_up_crude = c(3.5, 2.1)
  )
  suppressWarnings(
    obj <- convert_df(dat, measure = "logor", verbose = FALSE,
                      split_adjusted = FALSE, max_asymmetry = 10)
  )
  res <- summary(obj, flags = TRUE)
  # Row 1: ~9.2% asymmetry on log scale → should NOT fire at 10% threshold
  expect_false(grepl("Asymmetric CI", res$flags[1]))
})

test_that("V4 asymmetry still fires for genuinely asymmetric ratio-scale user_es", {
  # OR = 5.0, CI [0.5, 100] — highly asymmetric even on log scale
  # log(5)=1.609, log(0.5)=-0.693, log(100)=4.605
  # log-scale: 1.609-(-0.693)=2.302, 4.605-1.609=2.996 → ~26.2% asymmetry
  dat <- data.frame(
    user_es_original_measure_crude = c("or"),
    user_es_crude = c(5.0),
    user_ci_lo_crude = c(0.5),
    user_ci_up_crude = c(100)
  )
  suppressWarnings(
    obj <- convert_df(dat, measure = "logor", verbose = FALSE,
                      split_adjusted = FALSE, max_asymmetry = 10)
  )
  res <- summary(obj, flags = TRUE)
  expect_true(grepl("Asymmetric CI", res$flags[1]))
})

test_that("flags=TRUE and guidance=TRUE produce independent columns", {
  res <- summary(convert_df(df.haza, measure = "g", verbose = FALSE,
                            split_adjusted = FALSE),
                 flags = TRUE, guidance = TRUE)
  expect_true("flags" %in% colnames(res))
  expect_true("es_guidance" %in% colnames(res))
})


# ==============================================================================
# Fix: C5 with exp=TRUE checks on the log scale (not skipped, not natural-scale)
# ==============================================================================

test_that("[C5-fix] exp=TRUE checks the log of the natural-scale value", {
  # OR=10 -> log 2.3 < 5: not flagged. OR=200 -> log 5.3 > 5: flagged.
  # (Previously C5 was skipped entirely under exp=TRUE, silently disabling
  # the check for measure = "or"/"rr"/"irr"/"hr" users.)
  flags <- metaConvert:::.flag_plausibility(
    es = c(10, 200),
    se = c(0.3, 0.5),
    measure = "logor",
    opts = metaConvert:::.default_flag_options(),
    n_sample = c(100, 100),
    exp = TRUE
  )
  expect_false(any(grepl("Large log", flags[[1]])))
  expect_true(any(grepl("Large log OR", flags[[2]])))
})

test_that("[C5-fix] exp=TRUE ignores non-positive natural-scale values (B2's job)", {
  flags <- metaConvert:::.flag_plausibility(
    es = c(0, -2),
    se = c(0.3, 0.5),
    measure = "logor",
    opts = metaConvert:::.default_flag_options(),
    n_sample = c(100, 100),
    exp = TRUE
  )
  expect_false(any(grepl("Large log", flags[[1]])))
  expect_false(any(grepl("Large log", flags[[2]])))
})

test_that("[C5-fix] loghr is included in C5", {
  flags <- metaConvert:::.flag_plausibility(
    es = c(6),
    se = c(0.5),
    measure = "loghr",
    opts = metaConvert:::.default_flag_options(),
    n_sample = c(100),
    exp = FALSE
  )
  expect_true(any(grepl("Large log HR", flags[[1]])))
})

test_that("[C5-fix] exp=FALSE still flags large logOR", {
  flags <- metaConvert:::.flag_plausibility(
    es = c(6),
    se = c(0.5),
    measure = "logor",
    opts = metaConvert:::.default_flag_options(),
    n_sample = c(100),
    exp = FALSE
  )
  expect_true(any(grepl("Large log OR", flags[[1]])))
})

test_that("[C5-fix] logirr is now included in C5", {
  flags <- metaConvert:::.flag_plausibility(
    es = c(6),
    se = c(0.5),
    measure = "logirr",
    opts = metaConvert:::.default_flag_options(),
    n_sample = c(100),
    exp = FALSE
  )
  expect_true(any(grepl("Large log IRR", flags[[1]])))
})

test_that("[C3-removed] exp=TRUE no longer relevant (C3 removed)", {
  # C3 removed entirely — this test confirms no SE/ES ratio flag exists
  flags <- metaConvert:::.flag_plausibility(
    es = c(100),
    se = c(0.5),
    measure = "logor",
    opts = metaConvert:::.default_flag_options(),
    n_sample = c(50),
    exp = TRUE
  )
  expect_false(any(grepl("SE very large", flags[[1]])))
})


# ==============================================================================
# Fix: B7/B8 Bonett scale bounds
# ==============================================================================

test_that("[B7-fix] Bonett alpha: positive value implies negative alpha (UNUSUAL)", {
  # On the Bonett scale, ln(1-alpha) > 0 corresponds to alpha < 0 — possible
  # (negative inter-item covariance) but a serious data-quality signal, so it
  # is flagged [UNUSUAL] rather than [INVALID].
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(0.5, -0.5, 1.5),
    se = c(0.1, 0.1, 0.1),
    ci_lo = c(0.3, -0.7, 1.3),
    ci_up = c(0.7, -0.3, 1.7),
    measure = "alpha",
    exp = FALSE,
    alpha_to_es = "bonett"
  )
  expect_true(any(grepl("implies a negative Cronbach's alpha", flags[[1]])))  # 0.5 > 0
  expect_true(any(grepl("\\[UNUSUAL\\]", flags[[1]])))
  expect_false(any(grepl("\\[INVALID\\]", flags[[1]])))
  expect_false(any(grepl("alpha", flags[[2]], ignore.case = TRUE)))  # -0.5 < 0, valid
  expect_true(any(grepl("implies a negative Cronbach's alpha", flags[[3]])))  # 1.5 > 0
})

test_that("[B7-fix] raw alpha: only es > 1 is impossible", {
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(0.5, 1.2),
    se = c(0.1, 0.1),
    ci_lo = c(0.3, 1.0),
    ci_up = c(0.7, 1.4),
    measure = "alpha",
    exp = FALSE,
    alpha_to_es = "raw"
  )
  expect_false(any(grepl("Alpha", flags[[1]])))  # 0.5 < 1, valid on raw
  expect_true(any(grepl("Alpha exceeds 1", flags[[2]])))  # 1.2 > 1
})

test_that("[B8-fix] Bonett ICC: positive value implies negative ICC (UNUSUAL)", {
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(0.3, -1.0),
    se = c(0.1, 0.1),
    ci_lo = c(0.1, -1.2),
    ci_up = c(0.5, -0.8),
    measure = "icc",
    exp = FALSE,
    icc_to_es = "bonett"
  )
  expect_true(any(grepl("implies a negative ICC", flags[[1]])))  # 0.3 > 0
  expect_true(any(grepl("\\[UNUSUAL\\]", flags[[1]])))
  expect_false(any(grepl("\\[INVALID\\]", flags[[1]])))
  expect_false(any(grepl("ICC", flags[[2]])))  # -1.0 < 0, valid
})


# ==============================================================================
# V7: 2x2 additive consistency
# ==============================================================================

test_that("[V7] flags 2x2 additive inconsistency and sets the three cells to NA", {
  dat <- data.frame(
    n_cases_exp = c(10, 15),
    n_controls_exp = c(20, 25),
    n_exp = c(35, 40),  # Row 1: 10+20=30 != 35; Row 2: 15+25=40 ✓
    n_cases_nexp = c(5, 8),
    n_controls_nexp = c(25, 22),
    n_nexp = c(30, 30)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  # Row 1 should have the inconsistency flag
  expect_true(grepl("2x2 inconsistency", result$issues[1]))
  expect_true(grepl("n_cases_exp \\+ n_controls_exp != n_exp", result$issues[1]))
  # Row 1's three cells should now be NA (cannot determine which is wrong)
  expect_true(is.na(result$data$n_cases_exp[1]))
  expect_true(is.na(result$data$n_controls_exp[1]))
  expect_true(is.na(result$data$n_exp[1]))
  # Row 2 should NOT have the flag and its data should be preserved
  expect_false(grepl("2x2 inconsistency", result$issues[2]))
  expect_equal(result$data$n_cases_exp[2], 15)
  expect_equal(result$data$n_exp[2], 40)
})

test_that("[V7] without correct_inputs the data is preserved and the message shows actual sum/total", {
  dat <- data.frame(
    n_cases_exp = c(10),
    n_controls_exp = c(20),
    n_exp = c(35),  # 10 + 20 = 30 != 35
    n_cases_nexp = c(5),
    n_controls_nexp = c(25),
    n_nexp = c(30)
  )
  result <- metaConvert:::.validate_input_data(
    dat, verbose = FALSE, correct_inputs = FALSE
  )
  expect_true(grepl("2x2 inconsistency", result$issues[1]))
  expect_true(grepl("sum=30, total=35", result$issues[1]))
  # Data preserved when correct_inputs = FALSE
  expect_equal(result$data$n_cases_exp[1], 10)
  expect_equal(result$data$n_exp[1], 35)
})


# ==============================================================================
# V8: Zero cell in 2x2 table
# ==============================================================================

test_that("[V8] flags zero cell in 2x2 table when enable_informational = TRUE", {
  dat <- data.frame(
    n_cases_exp = c(0, 10),
    n_controls_exp = c(20, 20),
    n_cases_nexp = c(5, 5),
    n_controls_nexp = c(25, 25)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE,
                                                enable_informational = TRUE)
  expect_true(grepl("Zero cell in 2x2 table", result$issues[1]))
  expect_true(grepl("n_cases_exp.*= 0", result$issues[1]))
  expect_false(grepl("Zero cell", result$issues[2]))
  # Data should be preserved (warn only)
  expect_equal(result$data$n_cases_exp[1], 0)
})

test_that("[V8] zero cell flags suppressed by default (enable_informational = FALSE)", {
  dat <- data.frame(
    n_cases_exp = c(0, 10),
    n_controls_exp = c(20, 20),
    n_cases_nexp = c(5, 5),
    n_controls_nexp = c(25, 25)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("Zero cell", result$issues[1]))
})

# ==============================================================================
# V9: Possible change score as endpoint
# ==============================================================================

test_that("[V9] flags negative mean with small SD as possible change score", {
  dat <- data.frame(
    mean_exp = c(-3.6, 12.5),
    mean_sd_exp = c(0.25, 5.0),
    mean_nexp = c(-2.95, 15.0),
    mean_sd_nexp = c(0.20, 6.0)
  )
  # Only fires with enable_informational = TRUE
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE,
                                                enable_informational = TRUE)
  # Row 1: negative means + small SDs (SD/|mean| < 0.3) — flagged
  expect_true(grepl("Possible change score.*mean_exp", result$issues[1]))
  expect_true(grepl("Possible change score.*mean_nexp", result$issues[1]))
  # Row 2: positive means — no flag
  expect_false(grepl("change score", result$issues[2]))
})

test_that("[V9] suppressed by default (enable_informational = FALSE)", {
  dat <- data.frame(
    mean_exp = c(-3.6),
    mean_sd_exp = c(0.25),
    mean_nexp = c(-2.95),
    mean_sd_nexp = c(0.20)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("change score", result$issues[1]))
})

test_that("[V9] does NOT fire when SD is large relative to mean", {
  dat <- data.frame(
    mean_exp = c(-5.0),
    mean_sd_exp = c(10.0),
    mean_nexp = c(-3.0),
    mean_sd_nexp = c(8.0)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE,
                                                enable_informational = TRUE)
  # SD/|mean| = 2.0 and 2.67 — both > 0.3, no flag
  expect_false(grepl("change score", result$issues[1]))
})

test_that("[V9] preserves data (warn only, no NA)", {
  dat <- data.frame(
    mean_exp = c(-3.6),
    mean_sd_exp = c(0.25),
    mean_nexp = c(3.0),
    mean_sd_nexp = c(1.0)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE,
                                                enable_informational = TRUE)
  # Data preserved
  expect_equal(result$data$mean_exp[1], -3.6)
  expect_equal(result$data$mean_sd_exp[1], 0.25)
  # Only mean_exp flagged (negative + small SD), not mean_nexp (positive)
  expect_true(grepl("Possible change score.*mean_exp", result$issues[1]))
  expect_false(grepl("mean_nexp", result$issues[1]))
})


# ==============================================================================
# V10: Zero SD/SE detection
# ==============================================================================

test_that("[V10] flags zero SD and sets to NA", {
  dat <- data.frame(
    mean_exp = c(10.0, 10.0),
    mean_sd_exp = c(0, 5.0),
    mean_nexp = c(12.0, 12.0),
    mean_sd_nexp = c(3.0, 0)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  # Row 1: sd_exp=0 → flagged and set to NA
  expect_true(grepl("Zero SD.*mean_sd_exp", result$issues[1]))
  expect_true(is.na(result$data$mean_sd_exp[1]))
  # Row 1: sd_nexp=3.0 preserved
  expect_equal(result$data$mean_sd_nexp[1], 3.0)
  # Row 2: sd_nexp=0 → flagged and set to NA
  expect_true(grepl("Zero SD.*mean_sd_nexp", result$issues[2]))
  expect_true(is.na(result$data$mean_sd_nexp[2]))
  # Row 2: sd_exp=5.0 preserved
  expect_equal(result$data$mean_sd_exp[2], 5.0)
})

test_that("[V10] flags zero input SE and sets to NA", {
  dat <- data.frame(
    mean_exp = c(10.0),
    mean_se_exp = c(0),
    mean_nexp = c(12.0),
    mean_se_nexp = c(1.5)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("Zero SE.*mean_se_exp", result$issues[1]))
  expect_true(is.na(result$data$mean_se_exp[1]))
  expect_equal(result$data$mean_se_nexp[1], 1.5)
})

test_that("[V10] non-zero SD/SE are not affected", {
  dat <- data.frame(
    mean_exp = c(10.0),
    mean_sd_exp = c(5.0),
    mean_nexp = c(12.0),
    mean_sd_nexp = c(3.0)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("Zero SD", result$issues[1]))
  expect_equal(result$data$mean_sd_exp[1], 5.0)
  expect_equal(result$data$mean_sd_nexp[1], 3.0)
})

test_that("[V11] flags and sets p-value > 1 to NA", {
  dat <- data.frame(
    student_t_pval = c(0.05, 1.5, -0.1),
    anova_f_pval = c(0.8, 1.01, 0.5)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  # p = 1.5 should be set to NA (> 1)
  expect_true(is.na(result$data$student_t_pval[2]))
  expect_true(grepl("Out-of-range p-value.*student_t_pval", result$issues[2]))
  # p = -0.1 should be set to NA (< 0)
  expect_true(is.na(result$data$student_t_pval[3]))
  expect_true(grepl("Out-of-range p-value.*student_t_pval", result$issues[3]))
  # p = 0.05 should be untouched
  expect_equal(result$data$student_t_pval[1], 0.05)
  # anova_f_pval = 1.01 should be set to NA
  expect_true(is.na(result$data$anova_f_pval[2]))
  expect_true(grepl("Out-of-range p-value.*anova_f_pval", result$issues[2]))
})

# ==============================================================================
# V12: Extreme SD ratio between arms
# ==============================================================================

test_that("[V12] flags extreme SD ratio between arms", {
  dat <- data.frame(
    mean_sd_exp = c(55, 5, 8),
    mean_sd_nexp = c(5, 5, 10)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE, sd_ratio_max = 10)
  # Row 1: ratio = 11, should fire
  expect_true(grepl("Extreme SD ratio", result$issues[1]))
  # Row 2: ratio = 1, should not fire
  expect_false(grepl("Extreme SD ratio", result$issues[2]))
  # Row 3: ratio = 1.25, should not fire
  expect_false(grepl("Extreme SD ratio", result$issues[3]))
})

test_that("[V12] preserves data (warn only, no NA)", {
  dat <- data.frame(
    mean_sd_exp = c(100),
    mean_sd_nexp = c(5)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE, sd_ratio_max = 10)
  expect_equal(result$data$mean_sd_exp[1], 100)
  expect_equal(result$data$mean_sd_nexp[1], 5)
  expect_true(grepl("Extreme SD ratio", result$issues[1]))
})

test_that("[V12] respects configurable threshold", {
  dat <- data.frame(
    mean_sd_exp = c(30),
    mean_sd_nexp = c(5)
  )
  # ratio = 6: does not fire at default 10
  result_default <- metaConvert:::.validate_input_data(dat, verbose = FALSE, sd_ratio_max = 10)
  expect_false(grepl("Extreme SD ratio", result_default$issues[1]))
  # ratio = 6: fires at threshold = 5
  result_strict <- metaConvert:::.validate_input_data(dat, verbose = FALSE, sd_ratio_max = 5)
  expect_true(grepl("Extreme SD ratio", result_strict$issues[1]))
})

test_that("[V12] checks all SD pair types", {
  dat <- data.frame(
    mean_pre_sd_exp = c(50),
    mean_pre_sd_nexp = c(3)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE, sd_ratio_max = 10)
  expect_true(grepl("Extreme SD ratio.*mean_pre_sd", result$issues[1]))
})

test_that("[V12] skips when either SD is NA or zero", {
  dat <- data.frame(
    mean_sd_exp = c(50, 50),
    mean_sd_nexp = c(NA, 0)
  )
  # Zero SD will be caught by V10 first and set to NA, so V12 won't fire
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE, sd_ratio_max = 10)
  expect_false(grepl("Extreme SD ratio", result$issues[1]))
  expect_false(grepl("Extreme SD ratio", result$issues[2]))
})

# ==============================================================================
# V14: Identical group SDs
# ==============================================================================

test_that("[V14] flags identical SDs across arms when informational enabled", {
  dat <- data.frame(
    mean_sd_exp = c(3.5, 5.0),
    mean_sd_nexp = c(3.5, 6.0)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE,
                                                enable_informational = TRUE)
  # Row 1: identical SDs, should fire
  expect_true(grepl("Identical SDs", result$issues[1]))
  # Row 2: different SDs, should not fire
  expect_false(grepl("Identical SDs", result$issues[2]))
})

test_that("[V14] suppressed by default (enable_informational = FALSE)", {
  dat <- data.frame(
    mean_sd_exp = c(3.5),
    mean_sd_nexp = c(3.5)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE,
                                                enable_informational = FALSE)
  expect_false(grepl("Identical SDs", result$issues[1]))
})

test_that("[V14] preserves data (info only, no NA)", {
  dat <- data.frame(
    mean_sd_exp = c(3.5),
    mean_sd_nexp = c(3.5)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE,
                                                enable_informational = TRUE)
  expect_equal(result$data$mean_sd_exp[1], 3.5)
  expect_equal(result$data$mean_sd_nexp[1], 3.5)
})

# ==============================================================================
# D2: SE outlier normalization by sqrt(N)
# ==============================================================================

test_that("[D2] SE variation from sample size alone does not flag (normalized)", {
  opts <- metaConvert:::.default_flag_options()
  # 6 studies: identical true effect (d=0.5), varying N (equal groups)
  # SE = sqrt(1/n_exp + 1/n_nexp + d^2/(2*N)) matches expected formula exactly
  n_exps  <- c(15, 25, 40, 60, 100, 250)
  n_nexps <- c(15, 25, 40, 60, 100, 250)
  ns <- n_exps + n_nexps
  ses <- sapply(seq_along(ns), function(j)
    sqrt(1/n_exps[j] + 1/n_nexps[j] + 0.5^2/(2*ns[j])))
  # Raw SEs range widely but SE/expected_SE = 1 for all — no outlier
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = rep(0.5, 6), se = ses, opts = opts,
    n_total = ns, measure = "d",
    n_exp = n_exps, n_nexp = n_nexps
  )
  for (i in seq_along(flags)) {
    expect_false(any(grepl("SE outlier", flags[[i]])),
      info = paste("Study with N =", ns[i], "should not be SE outlier"))
  }
})

test_that("[D2] large effect size does not cause false SE outlier for d/g", {
  opts <- metaConvert:::.default_flag_options()
  # Husby 1993 scenario: d=-2.36 with N=36, SE=0.436 is mathematically correct
  # Other studies have smaller |d|, so raw SE is lower, but SE/expected_SE ≈ 1
  d_vals <- c(-0.74, -1.11, -2.36, -0.31, 0.47, -0.63, 0.00)
  n_exps  <- c(50, 46, 20, 42, 9, 35, 28)
  n_nexps <- c(30, 40, 16, 40, 8, 35, 27)
  ses <- sapply(seq_along(d_vals), function(j)
    sqrt(1/n_exps[j] + 1/n_nexps[j] + d_vals[j]^2/(2*(n_exps[j]+n_nexps[j]))))
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = d_vals, se = ses, opts = opts,
    n_total = n_exps + n_nexps, measure = "d",
    n_exp = n_exps, n_nexp = n_nexps
  )
  # No study should be flagged — all SEs match their expected values
  for (i in seq_along(flags)) {
    expect_false(any(grepl("SE outlier", flags[[i]])),
      info = paste("Study", i, "d =", d_vals[i], "should not be SE outlier"))
  }
})

test_that("[D2] genuine SE outlier still detected after normalization", {
  opts <- metaConvert:::.default_flag_options()
  # 5 normal studies + 1 with implausibly large SE for its sample size
  n_exps  <- c(25, 25, 25, 25, 25, 25)
  n_nexps <- c(25, 25, 25, 25, 25, 25)
  ns <- n_exps + n_nexps
  ses <- c(0.20, 0.21, 0.19, 0.22, 0.20, 2.5)  # last one is 12x too large
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = rep(0.5, 6), se = ses, opts = opts,
    n_total = ns, measure = "d",
    n_exp = n_exps, n_nexp = n_nexps
  )
  expect_true(any(grepl("SE outlier", flags[[6]])))
  expect_false(any(grepl("SE outlier", flags[[1]])))
})

test_that("[D2] falls back to raw SE when n_total is NULL", {
  opts <- metaConvert:::.default_flag_options()
  # Need IQR > 0, so use some variation
  ses <- c(0.15, 0.20, 0.25, 0.18, 0.22, 5.0)
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = rep(0.5, 6), se = ses, opts = opts,
    n_total = NULL
  )
  # SE=5.0 should be flagged as outlier (raw comparison, no normalization)
  expect_true(any(grepl("SE outlier", flags[[6]])))
})

test_that("[D2] logOR: rare-event studies not falsely flagged as SE outliers", {
  opts <- metaConvert:::.default_flag_options()
  # Scenario: 6 studies, all with n_exp=100, n_nexp=100.
  # 4 common-event studies (event rate ~50%, SE_logOR ~ 0.28)
  # 1 rare-event study (event rate ~5%, SE_logOR ~ 0.65)
  # 1 genuine outlier (impossible SE for its N)
  # With old sqrt(N) normalization, the rare-event study would be flagged.
  # With 2x2 reconstruction, it should NOT be flagged because its SE matches
  # the expected SE given its cell counts.

  # Common events: logOR ~ 0.5, SE from sqrt(1/50+1/50+1/50+1/50) ~ 0.283
  # Rare events: logOR ~ 0.7, approx 5 cases per group
  #   SE ~ sqrt(1/5 + 1/95 + 1/5 + 1/95) ~ 0.64
  # Genuine outlier: SE=10.0 is wildly impossible for N=200 (data entry error)
  es_vals  <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.7, 0.5)
  se_vals  <- c(0.28, 0.28, 0.29, 0.27, 0.28, 0.28, 0.28, 0.64, 10.0)
  n_exps   <- rep(100, 9)
  n_nexps  <- rep(100, 9)
  ns       <- n_exps + n_nexps

  flags <- metaConvert:::.flag_cross_row_outliers(
    es = es_vals, se = se_vals, opts = opts,
    n_total = ns, measure = "logor",
    n_exp = n_exps, n_nexp = n_nexps
  )
  # Rare-event study (8): NOT an outlier — SE matches its event rate
  expect_false(any(grepl("SE outlier", flags[[8]])),
    info = "Rare-event study should not be falsely flagged with 2x2 reconstruction")
  # Genuine outlier (9): SE=10.0 is impossible for N=200
  expect_true(any(grepl("SE outlier", flags[[9]])),
    info = "Genuine SE outlier should still be detected")
})

test_that("[D2] extreme group imbalance does not crash d/g normalization", {
  opts <- metaConvert:::.default_flag_options()
  # Extreme imbalance: n_exp=500, n_nexp=5
  # SE = sqrt(1/500 + 1/5 + 0.5^2/(2*505)) ≈ 0.456
  n_exps  <- c(50, 50, 50, 50, 500)
  n_nexps <- c(50, 50, 50, 50, 5)
  ns <- n_exps + n_nexps
  ses <- sapply(seq_along(ns), function(j)
    sqrt(1/n_exps[j] + 1/n_nexps[j] + 0.5^2/(2*ns[j])))
  # Should not crash and should not flag (SEs match expected)
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = rep(0.5, 5), se = ses, opts = opts,
    n_total = ns, measure = "d",
    n_exp = n_exps, n_nexp = n_nexps
  )
  # No crash is already a pass; additionally verify no false positive
  for (i in seq_along(flags)) {
    expect_false(any(grepl("SE outlier", flags[[i]])),
      info = paste("Study", i, "should not be SE outlier"))
  }
})

test_that("[F1/F2] logIRR is included in standardized measures", {
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(0.5), se = c(0.0001), measure = "logirr",
    n_total = c(100)
  )
  # SE=0.0001 with N=100: floor = 0.1/sqrt(100) = 0.01
  # 0.0001 < 0.01 → should flag
  expect_true(any(grepl("Implausibly small SE", flags[[1]])))
})

test_that("[F3] flags SE below Cox/Wald theoretical floor for loghr", {
  # Feng-style per-10%-PDC HR pooled with binary HRs:
  # log(HR)=log(0.93)=-0.0726, SE_log = (log(0.94) - log(0.92))/3.92 ≈ 0.0055
  # n_exp=52514, n_nexp=227 → floor = sqrt(1/52514 + 1/227) ≈ 0.0665
  # 0.0055 < 0.0665 → should flag F3
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(log(0.93)),
    se = c((log(0.94) - log(0.92)) / 3.92),
    measure = "loghr",
    n_total = c(52741),
    n_exp = c(52514),
    n_nexp = c(227)
  )
  expect_true(any(grepl("SE below theoretical minimum given group sizes", flags[[1]])))
  expect_true(any(grepl("per-unit-slope model, wrong N, or wrong SE", flags[[1]])))
})

test_that("[F3] does not flag legitimate large RCT with tight CI (common outcome)", {
  # 500k-arm RCT with common outcome: HR=0.85, SE_log=0.0153
  # Floor = sqrt(2/500000) ≈ 0.002. 0.0153 >> 0.002 → no F3 flag.
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(log(0.85)), se = c(0.0153),
    measure = "loghr",
    n_total = c(1000000), n_exp = c(500000), n_nexp = c(500000)
  )
  expect_false(any(grepl("SE below theoretical minimum", flags[[1]])))
})

test_that("[F3] does not flag rare-outcome large cohort (Nielsen-style)", {
  # Nielsen 2012 Pass A: HR=0.81, SE_log=0.038, n_exp=18721, n_nexp=277204
  # Floor = sqrt(1/18721 + 1/277204) ≈ 0.0076. 0.038 >> 0.0076 → no flag.
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(log(0.81)), se = c(0.038),
    measure = "loghr",
    n_total = c(295925), n_exp = c(18721), n_nexp = c(277204)
  )
  expect_false(any(grepl("SE below theoretical minimum", flags[[1]])))
})

test_that("[F3] applies to loghr / logirr but not to logor / logrr / d / g", {
  # The Cox/Wald floor sqrt(1/n_exp + 1/n_nexp) is the minimum of
  # SE_logHR = sqrt(1/d_exp + 1/d_nexp) when d = n. It applies to loghr and
  # logirr (Cox / Poisson rate models). For binary 2x2 OR and RR the SE
  # formula depends on non-event counts too and can approach 0 as events -> N,
  # so the floor does not apply.

  # loghr with SE between F1 floor (0.1/sqrt(300)=0.0058) and F3 floor
  # (sqrt(2/150)=0.1155): F1 silent, F3 should fire.
  flags_loghr <- metaConvert:::.flag_se_sample_size(
    es = c(log(0.5)), se = c(0.05),
    measure = "loghr",
    n_total = c(300), n_exp = c(150), n_nexp = c(150)
  )
  expect_true(any(grepl("SE below theoretical minimum given group sizes", flags_loghr[[1]])))

  # logor with same numerics: F3 should NOT fire (binary 2x2 SE_logOR has
  # no useful N-only floor).
  flags_logor <- metaConvert:::.flag_se_sample_size(
    es = c(log(0.5)), se = c(0.05),
    measure = "logor",
    n_total = c(300), n_exp = c(150), n_nexp = c(150)
  )
  expect_false(any(grepl("theoretical minimum given group sizes", flags_logor[[1]])))

  # logrr with same numerics: F3 should NOT fire (binary 2x2 SE_logRR has
  # no useful N-only floor).
  flags_logrr <- metaConvert:::.flag_se_sample_size(
    es = c(log(0.5)), se = c(0.05),
    measure = "logrr",
    n_total = c(300), n_exp = c(150), n_nexp = c(150)
  )
  expect_false(any(grepl("theoretical minimum given group sizes", flags_logrr[[1]])))

  # d-measure with same numerics: F3 should NOT fire (Cox/Wald floor doesn't
  # apply to standardized mean differences).
  flags_d <- metaConvert:::.flag_se_sample_size(
    es = c(0.5), se = c(0.05),
    measure = "d",
    n_total = c(300), n_exp = c(150), n_nexp = c(150)
  )
  expect_false(any(grepl("theoretical minimum given group sizes", flags_d[[1]])))
})

test_that("[F3] does not duplicate F1's 'Implausibly small SE' message", {
  # When SE is below BOTH F1's 0.1/sqrt(N) floor AND F3's Cox/Wald floor,
  # only F1 should fire (F3 explicitly skips to avoid duplicate messages).
  # Tested on loghr (where F3 is active); F3 is gated out for logor/logrr
  # so this dedup logic only matters for loghr/logirr.
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(log(0.5)), se = c(0.005),
    measure = "loghr",
    n_total = c(300), n_exp = c(150), n_nexp = c(150)
  )
  f1_count <- sum(grepl("Implausibly small SE", flags[[1]]))
  f3_count <- sum(grepl("SE below theoretical minimum", flags[[1]]))
  expect_equal(f1_count, 1)
  expect_equal(f3_count, 0)
})

test_that("[F3] silent when n_exp / n_nexp missing (falls back to F1)", {
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(log(0.93)), se = c(0.0055),
    measure = "loghr",
    n_total = c(52741),
    n_exp = c(NA_real_), n_nexp = c(NA_real_)
  )
  expect_false(any(grepl("SE below theoretical minimum", flags[[1]])))
})

# ==============================================================================
# E5: Direction disagreement
# ==============================================================================

test_that("[E5] flags direction disagreement across methods", {
  # Create a result dataframe that mimics what .generate_df() produces
  res <- data.frame(
    row_id = c(1, 2),
    es = c(0.3, -0.2),
    se = c(0.1, 0.1),
    es_ci_lo = c(0.1, -0.4),
    es_ci_up = c(0.5, 0.0),
    n_estimations = c(3, 3),
    dispersion_es = c(0.2, 0.1),
    overlap_min_max = c(0.5, 0.9),
    diff_min_max = c(0.6, 0.3),
    info_used = c("means_sd", "student_t"),
    min_info = c("student_t", "means_sd"),
    max_info = c("means_sd", "student_t"),
    min_es_value = c(-0.1, -0.3),  # Row 1: min negative, max positive → disagreement
    max_es_value = c(0.5, -0.1),   # Row 2: both negative → no disagreement
    measure = c("g", "g"),
    stringsAsFactors = FALSE
  )
  raw_data <- data.frame(row_id = c(1, 2), n_sample = c(50, 50))
  opts <- metaConvert:::.default_flag_options()

  flagged <- metaConvert:::.flag_es_quality(res, "g", FALSE, "", raw_data, opts)
  expect_true(grepl("Direction disagreement", flagged$flags[1]))
  expect_true(grepl("check group labeling", flagged$flags[1]))
  expect_false(grepl("Direction disagreement", flagged$flags[2]))
})

test_that("[E5] skips direction check for unsigned measures", {
  # prop is not a signed measure — should not check direction
  res <- data.frame(
    row_id = c(1),
    es = c(0.5),
    se = c(0.05),
    es_ci_lo = c(0.4),
    es_ci_up = c(0.6),
    n_estimations = c(2),
    dispersion_es = c(0.1),
    overlap_min_max = c(0.8),
    diff_min_max = c(0.2),
    info_used = c("prop_single_group"),
    min_info = c("prop_single_group"),
    max_info = c("prop_single_group_counts"),
    min_es_value = c(-0.1),
    max_es_value = c(0.5),
    measure = c("prop"),
    stringsAsFactors = FALSE
  )
  raw_data <- data.frame(row_id = c(1), n_sample = c(100))
  opts <- metaConvert:::.default_flag_options()

  flagged <- metaConvert:::.flag_es_quality(res, "prop", FALSE, "", raw_data, opts)
  expect_false(grepl("Direction disagreement", flagged$flags[1]))
})

test_that("[E5] suppressed by ignore_sign_reversal when all rows are reversed", {
  # All rows have min ≈ -max (systematic reversal)
  res <- data.frame(
    row_id = c(1, 2, 3),
    es = c(0.3, 0.2, 0.5),
    se = c(0.1, 0.1, 0.1),
    es_ci_lo = c(0.1, 0.0, 0.3),
    es_ci_up = c(0.5, 0.4, 0.7),
    n_estimations = c(2, 2, 2),
    dispersion_es = c(0.1, 0.1, 0.1),
    overlap_min_max = c(0.5, 0.5, 0.5),
    diff_min_max = c(0.2, 0.2, 0.2),
    info_used = c("means_sd", "means_sd", "means_sd"),
    min_info = c("student_t", "student_t", "student_t"),
    max_info = c("means_sd", "means_sd", "means_sd"),
    min_es_value = c(-0.5, -0.3, -0.8),  # min ≈ -max for all rows
    max_es_value = c(0.5, 0.3, 0.8),
    stringsAsFactors = FALSE
  )
  raw_data <- data.frame(row_id = c(1, 2, 3), n_sample = c(50, 50, 50))
  opts <- metaConvert:::.default_flag_options()
  opts$ignore_sign_reversal <- TRUE

  flagged <- metaConvert:::.flag_es_quality(res, "g", FALSE, "", raw_data, opts)
  for (i in 1:3) {
    expect_false(grepl("Direction disagreement", flagged$flags[i]),
      info = paste("Row", i, "E5 should be suppressed when all reversed + ignore_sign_reversal"))
  }

  # But with ignore_sign_reversal = FALSE, E5 should fire
  opts$ignore_sign_reversal <- FALSE
  flagged2 <- metaConvert:::.flag_es_quality(res, "g", FALSE, "", raw_data, opts)
  expect_true(grepl("Direction disagreement", flagged2$flags[1]),
    info = "E5 should fire when ignore_sign_reversal=FALSE")
})


# ==============================================================================
# A6: SE/CI width consistency
# ==============================================================================

test_that("[A6] flags CI width inconsistent with SE", {
  # Expected width = 2 * 1.96 * 0.1 = 0.392; actual = 1.0 → 155% deviation
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(0.5),
    se = c(0.1),
    ci_lo = c(0.0),
    ci_up = c(1.0),
    measure = "g"
  )
  expect_true(any(grepl("CI width inconsistent with SE", flags[[1]])))
})

test_that("[A6] does not flag consistent CI/SE", {
  se_val <- 0.2
  expected_half <- stats::qnorm(0.975) * se_val
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(0.5),
    se = c(se_val),
    ci_lo = c(0.5 - expected_half),
    ci_up = c(0.5 + expected_half),
    measure = "g"
  )
  expect_false(any(grepl("CI width inconsistent", flags[[1]])))
})

test_that("[A6] skips for NNT (discontinuous CIs)", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(5.0),
    se = c(2.0),
    ci_lo = c(-20),
    ci_up = c(3.0),
    measure = "nnt"
  )
  expect_false(any(grepl("CI width inconsistent", flags[[1]])))
})

test_that("[A6] skips when exp=TRUE (mixed scales)", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(2.5),   # natural-scale OR
    se = c(0.3),   # log-scale SE
    ci_lo = c(1.2),
    ci_up = c(5.1),
    measure = "logor",
    exp = TRUE
  )
  expect_false(any(grepl("CI width inconsistent", flags[[1]])))
})

test_that("[A6] tighter tolerance for large N catches 7% deviation", {
  # With N=1000, tolerance = max(0.05, 0.30/sqrt(1000)) ≈ 5%
  # A 7% deviation should be flagged at N=1000 but NOT at N=9 (tolerance ~10%)
  se_val <- 0.1
  z <- stats::qnorm(0.975)
  expected_width <- 2 * z * se_val  # 0.392
  # 7% wider than expected
  actual_width <- expected_width * 1.07
  es_val <- 0.5
  ci_lo_val <- es_val - actual_width / 2
  ci_up_val <- es_val + actual_width / 2

  # Large N: should flag (tol ~5%, deviation 7%)
  flags_large <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(ci_lo_val), ci_up = c(ci_up_val),
    measure = "logor", n_total = c(1000)
  )
  expect_true(any(grepl("CI width inconsistent", flags_large[[1]])),
    info = "7% deviation should be flagged at N=1000 (tol ~5%)")

  # Small N: should NOT flag (tol ~10%, deviation 7%)
  flags_small <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(ci_lo_val), ci_up = c(ci_up_val),
    measure = "logor", n_total = c(9)
  )
  expect_false(any(grepl("CI width inconsistent", flags_small[[1]])),
    info = "7% deviation should NOT be flagged at N=9 (tol ~10%)")
})

test_that("[A6] uses qt() for d/g with small N", {
  # With N=8 (df=6), qt(.975,6) = 2.447 vs qnorm(.975) = 1.96
  # Build CI using qt() — should match A6 expectations
  se_val <- 0.4
  n_val <- 8
  z_crit <- stats::qt(0.975, n_val - 2)  # 2.447
  expected_half <- z_crit * se_val  # 0.979
  es_val <- 0.5
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - expected_half),
    ci_up = c(es_val + expected_half),
    measure = "g", n_total = c(n_val)
  )
  expect_false(any(grepl("CI width inconsistent", flags[[1]])),
    info = "CI built with qt() should not be flagged for d/g measure")

  # Same CI but using qnorm() — should be flagged because it's too narrow
  z_norm <- stats::qnorm(0.975)
  narrow_half <- z_norm * se_val  # 0.784
  flags_narrow <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - narrow_half),
    ci_up = c(es_val + narrow_half),
    measure = "g", n_total = c(n_val)
  )
  # Expected width with qt: 1.958; actual with qnorm: 1.568; deviation ~20%
  expect_true(any(grepl("CI width inconsistent", flags_narrow[[1]])),
    info = "CI built with qnorm should be flagged at small N for d/g (qt expected)")
})

# ==============================================================================
# n_min default changed to 10
# ==============================================================================

test_that("n_min default is now 10", {
  opts <- metaConvert:::.default_flag_options()
  expect_equal(opts$n_min, 10)
})

test_that("[C6] n=7 is suppressed by default but fires with enable_informational", {
  # Default: suppressed
  flags_default <- metaConvert:::.flag_plausibility(
    es = c(0.5),
    se = c(0.3),
    measure = "g",
    opts = metaConvert:::.default_flag_options(),
    n_sample = c(7)
  )
  expect_false(any(grepl("Small sample size", flags_default[[1]])))

  # With enable_informational = TRUE: fires
  opts <- metaConvert:::.default_flag_options()
  opts$enable_informational <- TRUE
  flags_info <- metaConvert:::.flag_plausibility(
    es = c(0.5),
    se = c(0.3),
    measure = "g",
    opts = opts,
    n_sample = c(7)
  )
  expect_true(any(grepl("Small sample size", flags_info[[1]])))
})


# ==============================================================================
# Category G: Cross-Row Direction Conflict
# ==============================================================================

test_that("[G1] flags all rows on both sides of direction conflict", {
  opts <- metaConvert:::.default_flag_options()
  # 3 sig-positive (ci_lo > 0), 2 sig-negative (ci_up < 0), 1 non-significant
  es     <- c(0.48, 0.38, 0.77, -1.38, -0.63, 0.00)
  ci_lo  <- c(0.05, 0.10, 0.27, -2.04, -1.16, -0.41)
  ci_up  <- c(0.91, 0.66, 1.27, -0.73, -0.10,  0.42)
  info   <- rep("means_sd", 6)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, info_used = info, measure = "g"
  )

  # Sig-positive rows (1, 2, 3) should be flagged
  expect_true(any(grepl("significantly positive", flags[[1]])))
  expect_true(any(grepl("significantly positive", flags[[2]])))
  expect_true(any(grepl("significantly positive", flags[[3]])))
  # Sig-negative rows (4, 5) should be flagged
  expect_true(any(grepl("significantly negative", flags[[4]])))
  expect_true(any(grepl("significantly negative", flags[[5]])))
  # Non-significant row (6) should NOT be flagged
  expect_equal(length(flags[[6]]), 0)
})

test_that("[G1] silent when a lone outlier is below BOTH count and lone_pct", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_min <- 2
  # 5 sig-positive, 1 sig-negative (1/6 = 16.7% < lone_pct 20%): the single
  # negative clears neither the count floor (2) nor the lone-share gate (20%),
  # so the pool stays silent.
  es     <- c(0.5, 0.4, 0.6, 0.55, 0.45, -0.8)
  ci_lo  <- c(0.1, 0.1, 0.2, 0.15, 0.05, -1.5)
  ci_up  <- c(0.9, 0.7, 1.0, 0.95, 0.85, -0.1)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  for (i in seq_along(flags)) {
    expect_equal(length(flags[[i]]), 0)
  }
})

test_that("[G1] lone significant outlier fires when it is >= lone_pct of a small pool", {
  opts <- metaConvert:::.default_flag_options()  # lone_pct default = 20
  # Dauw scenario: 1 sig-positive (1/5 = 20%), 2 sig-negative (40%), 2 non-sig.
  es     <- c(0.90, -1.20, -1.18, -0.40, -1.00)
  ci_lo  <- c(0.46, -2.22, -1.85, -1.62, -2.68)
  ci_up  <- c(1.34, -0.18, -0.51,  0.82,  0.68)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "md"
  )
  # The lone positive (row 1) and both significant negatives (rows 2, 3) fire.
  expect_true(any(grepl("significantly positive", flags[[1]])))
  expect_true(any(grepl("significantly negative", flags[[2]])))
  expect_true(any(grepl("significantly negative", flags[[3]])))
  # Non-significant rows stay silent.
  expect_equal(length(flags[[4]]), 0)
  expect_equal(length(flags[[5]]), 0)
})

test_that("[G1] lone outlier below lone_pct in a large pool stays silent", {
  opts <- metaConvert:::.default_flag_options()
  # 1 sig-positive (1/10 = 10% < 20%), 2 sig-negative, 7 non-sig: the lone
  # positive clears neither the count floor nor the 20% lone gate.
  es     <- c(0.5, -0.8, -0.9, rep(0.1, 7))
  ci_lo  <- c(0.1, -1.5, -1.6, rep(-0.5, 7))
  ci_up  <- c(0.9, -0.1, -0.2, rep(0.7, 7))

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  for (i in seq_along(flags)) {
    expect_equal(length(flags[[i]]), 0)
  }
})

test_that("[G1] lone-outlier branch is disabled with direction_conflict_lone_pct = Inf", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_lone_pct <- Inf  # revert to pure count-floor behaviour
  # Same Dauw scenario as above, but the lone positive no longer qualifies.
  es     <- c(0.90, -1.20, -1.18, -0.40, -1.00)
  ci_lo  <- c(0.46, -2.22, -1.85, -1.62, -2.68)
  ci_up  <- c(1.34, -0.18, -0.51,  0.82,  0.68)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "md"
  )
  for (i in seq_along(flags)) {
    expect_equal(length(flags[[i]]), 0)
  }
})

test_that("[G1] minority side gets [UNUSUAL], majority side gets [INFO]", {
  opts <- metaConvert:::.default_flag_options()
  # 2 sig-positive (majority), 1 sig-negative (minority), 1 non-sig
  es     <- c(0.5, 0.6, -0.8, 0.1)
  ci_lo  <- c(0.1, 0.2, -1.5, -0.4)
  ci_up  <- c(0.9, 1.0, -0.1,  0.6)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  # majority = positive (2 > 1): positive rows -> INFO, the lone negative -> UNUSUAL
  expect_true(any(grepl("[INFO]",    flags[[1]], fixed = TRUE)))
  expect_true(any(grepl("[INFO]",    flags[[2]], fixed = TRUE)))
  expect_true(any(grepl("[UNUSUAL]", flags[[3]], fixed = TRUE)))
  expect_equal(length(flags[[4]]), 0)
})

test_that("[G1] Dauw-style lone outlier is [UNUSUAL], its opposite cluster is [INFO]", {
  opts <- metaConvert:::.default_flag_options()
  es     <- c(0.90, -1.20, -1.18, -0.40, -1.00)  # 1 sig-pos (minority), 2 sig-neg (majority)
  ci_lo  <- c(0.46, -2.22, -1.85, -1.62, -2.68)
  ci_up  <- c(1.34, -0.18, -0.51,  0.82,  0.68)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "md"
  )
  expect_true(any(grepl("[UNUSUAL]", flags[[1]], fixed = TRUE)))  # Dauw outlier
  expect_true(any(grepl("[INFO]",    flags[[2]], fixed = TRUE)))  # majority cluster
  expect_true(any(grepl("[INFO]",    flags[[3]], fixed = TRUE)))
})

test_that("[G1] tie (no majority) makes every significant row [UNUSUAL]", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_min <- 1
  es     <- c(0.5, -0.8)
  ci_lo  <- c(0.1, -1.5)
  ci_up  <- c(0.9, -0.1)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  expect_true(any(grepl("[UNUSUAL]", flags[[1]], fixed = TRUE)))
  expect_true(any(grepl("[UNUSUAL]", flags[[2]], fixed = TRUE)))
})

test_that("[G1] direction_conflict_info = FALSE suppresses the majority [INFO] flags", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_info <- FALSE
  es     <- c(0.5, 0.6, -0.8)  # majority positive
  ci_lo  <- c(0.1, 0.2, -1.5)
  ci_up  <- c(0.9, 1.0, -0.1)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  expect_equal(length(flags[[1]]), 0)  # majority rows suppressed
  expect_equal(length(flags[[2]]), 0)
  expect_true(any(grepl("[UNUSUAL]", flags[[3]], fixed = TRUE)))  # outlier still flagged
})

test_that("[G1] silent when below direction_conflict_pct threshold", {
  opts <- metaConvert:::.default_flag_options()
  # 100 studies: 1 sig-positive (1%), 1 sig-negative (1%), rest non-significant
  # Below 10% pct threshold, so should NOT fire
  es     <- c(0.5, -0.8, rep(0.1, 98))
  ci_lo  <- c(0.1, -1.5, rep(-0.5, 98))
  ci_up  <- c(0.9, -0.1, rep(0.7, 98))

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  for (i in seq_along(flags)) {
    expect_equal(length(flags[[i]]), 0)
  }
})

test_that("[G1] fires with direction_conflict_min = 1", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_min <- 1
  # 1 sig-positive, 1 sig-negative
  es     <- c(0.5, -0.8)
  ci_lo  <- c(0.1, -1.5)
  ci_up  <- c(0.9, -0.1)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  expect_true(any(grepl("significantly positive", flags[[1]])))
  expect_true(any(grepl("significantly negative", flags[[2]])))
})

test_that("[G1] silent for non-signed measures", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_min <- 1
  es     <- c(0.5, -0.8)
  ci_lo  <- c(0.1, -1.5)
  ci_up  <- c(0.9, -0.1)

  for (m in c("nnt", "rd", "prop", "alpha", "icc")) {
    flags <- metaConvert:::.flag_cross_row_direction_conflict(
      es, ci_lo, ci_up, opts, measure = m
    )
    for (i in seq_along(flags)) {
      expect_equal(length(flags[[i]]), 0, info = paste("measure:", m))
    }
  }
})

test_that("[G1] respects enable_cross_row = FALSE", {
  opts <- metaConvert:::.default_flag_options()
  opts$enable_cross_row <- FALSE
  opts$direction_conflict_min <- 1
  es     <- c(0.5, -0.8)
  ci_lo  <- c(0.1, -1.5)
  ci_up  <- c(0.9, -0.1)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  for (i in seq_along(flags)) {
    expect_equal(length(flags[[i]]), 0)
  }
})

test_that("[G1] handles NA CI values gracefully", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_min <- 1
  es     <- c(0.5, -0.8, NA)
  ci_lo  <- c(0.1, -1.5, NA)
  ci_up  <- c(0.9, -0.1, NA)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "g"
  )
  expect_true(any(grepl("significantly positive", flags[[1]])))
  expect_true(any(grepl("significantly negative", flags[[2]])))
  expect_equal(length(flags[[3]]), 0)
})

test_that("[G1] handles exp=TRUE with null at 1 for logor", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_min <- 1
  # On exponentiated scale: OR > 1 is "positive", OR < 1 is "negative"
  es     <- c(2.5, 0.3)
  ci_lo  <- c(1.2, 0.1)  # ci_lo > 1 means sig positive
  ci_up  <- c(5.0, 0.8)  # ci_up < 1 means sig negative

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "logor", exp = TRUE
  )
  expect_true(any(grepl("significantly positive", flags[[1]])))
  expect_true(any(grepl("significantly negative", flags[[2]])))
})

test_that("[G1] includes method traceability in message", {
  opts <- metaConvert:::.default_flag_options()
  opts$direction_conflict_min <- 1
  es     <- c(0.5, -0.8)
  ci_lo  <- c(0.1, -1.5)
  ci_up  <- c(0.9, -0.1)
  info   <- c("means_sd", "student_t")

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, info_used = info, measure = "g"
  )
  expect_true(any(grepl("from means_sd", flags[[1]])))
  expect_true(any(grepl("from student_t", flags[[2]])))
})

test_that("[G1] message includes opposing count", {
  opts <- metaConvert:::.default_flag_options()
  # 3 sig-positive (majority), 2 sig-negative (minority)
  es     <- c(0.5, 0.4, 0.6, -0.8, -0.7)
  ci_lo  <- c(0.1, 0.1, 0.2, -1.5, -1.3)
  ci_up  <- c(0.9, 0.7, 1.0, -0.1, -0.1)

  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es, ci_lo, ci_up, opts, measure = "d"
  )
  # Majority (positive) row -> [INFO], mentions the 2 opposing rows
  expect_true(any(grepl("2 significant row(s) point the opposite way", flags[[1]], fixed = TRUE)))
  # Minority (negative) row -> [UNUSUAL], mentions the 3 majority rows
  expect_true(any(grepl("3 significant row(s) favour the other direction", flags[[4]], fixed = TRUE)))
})


# ==============================================================================
# Severity prefixes
# ==============================================================================

test_that("Category A/B flags use [INVALID] prefix", {
  # A1: ES is Inf
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(Inf), se = c(0.1), ci_lo = c(0), ci_up = c(1)
  )
  expect_true(grepl("^\\[INVALID\\]", flags[[1]][1]))

  # B1: r > 1
  flags <- metaConvert:::.flag_bounds_violations(
    es = c(1.5), se = c(0.1), ci_lo = c(1.0), ci_up = c(2.0),
    measure = "r", exp = FALSE
  )
  expect_true(grepl("^\\[INVALID\\]", flags[[1]][1]))
})

test_that("Category C/D/F flags use [UNUSUAL] prefix", {
  opts <- metaConvert:::.default_flag_options()
  # C1: large SMD
  flags <- metaConvert:::.flag_plausibility(
    es = c(5.0), se = c(0.5), measure = "d",
    opts = opts, n_sample = c(50)
  )
  expect_true(grepl("^\\[UNUSUAL\\]", flags[[1]][1]))

  # D1: ES outlier (needs 4+ studies, IQR > 0, deviation > min_dev=1.0 for d)
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = c(0.1, 0.2, 0.3, 0.15, 0.25, 50.0),
    se = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    opts = opts, measure = "d"
  )
  expect_true(grepl("^\\[UNUSUAL\\]", flags[[6]][1]))
})

test_that("Category E flags use [DISCORDANT] prefix", {
  # E flags: cross-method consistency checks fire when multiple ES estimates
  # for the same row disagree. These are now classified as the Discordant
  # family in the framework taxonomy.
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(3), dispersion = c(1.0),
    overlap = c(0.9), diff_mm = c(0.1),
    opts = opts, measure = "d"
  )
  expect_true(grepl("^\\[DISCORDANT\\]", flags[[1]][1]))
})

test_that("Category G flags use [UNUSUAL] (minority/outlier) and [INFO] (majority) prefixes", {
  # Cross-row direction conflict is split by leverage: the minority/outlier
  # side escalates to [UNUSUAL]; the majority side is surfaced as [INFO] context.
  opts <- metaConvert:::.default_flag_options()
  # Clear majority: 2 sig-positive (majority), 1 sig-negative (minority)
  flags <- metaConvert:::.flag_cross_row_direction_conflict(
    es = c(0.5, 0.6, -0.8), ci_lo = c(0.1, 0.2, -1.5), ci_up = c(0.9, 1.0, -0.1),
    opts = opts, measure = "g"
  )
  expect_true(grepl("^\\[INFO\\]",    flags[[1]][1]))  # majority side
  expect_true(grepl("^\\[UNUSUAL\\]", flags[[3]][1]))  # minority/outlier side

  # Tie (no majority): every significant row is an outlier -> [UNUSUAL]
  opts$direction_conflict_min <- 1
  flags2 <- metaConvert:::.flag_cross_row_direction_conflict(
    es = c(0.5, -0.8), ci_lo = c(0.1, -1.5), ci_up = c(0.9, -0.1),
    opts = opts, measure = "g"
  )
  expect_true(grepl("^\\[UNUSUAL\\]", flags2[[1]][1]))
})

test_that("V-flags use [INVALID], [UNUSUAL], or [INFO] prefix by family", {
  # Validation flags carry one of four family prefixes:
  #   [INVALID]    — data is invalid and will be replaced with NA
  #   [UNUSUAL]    — data preserved, plausibility warning (e.g. extreme SD
  #                  ratio between arms)
  #   [DISCORDANT] — within-row or cross-row inconsistency
  #   [INFO]       — soft informational flag (e.g. asymmetric CI)

  # V1: negative input → [INVALID] (value is set to NA, cannot be used)
  dat <- data.frame(n_exp = c(-5))
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("^\\[INVALID\\]", result$issues[[1]][1]))

  # V4: asymmetric CI → [INFO] (data preserved; often just rounding)
  dat <- data.frame(
    mean_exp = c(10), mean_ci_lo_exp = c(5), mean_ci_up_exp = c(20)
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  if (length(result$issues[[1]]) > 0 && nchar(result$issues[[1]][1]) > 0) {
    expect_true(grepl("^\\[INFO\\]", result$issues[[1]][1]))
  }
})


# ==============================================================================
# Category H: Cross-row study duplication
# ==============================================================================

test_that("[H1] duplicate study_id flags every row in the duplicate group", {
  flags <- metaConvert:::.flag_cross_row_duplicates(
    study_id = c("ELOQUENT-2", "ELOQUENT-1", "ELOQUENT-2"),
    opts     = list(enable_cross_row = TRUE)
  )
  expect_true(any(grepl("Duplicate study_id 'ELOQUENT-2'", flags[[1]], fixed = TRUE)))
  expect_true(any(grepl("ELOQUENT-2 (row 3)", flags[[1]], fixed = TRUE)))
  expect_equal(length(flags[[2]]), 0)
  expect_true(any(grepl("Duplicate study_id 'ELOQUENT-2'", flags[[3]], fixed = TRUE)))
  expect_true(any(grepl("ELOQUENT-2 (row 1)", flags[[3]], fixed = TRUE)))
  # Severity prefix
  expect_true(grepl("^\\[INFO\\]", flags[[1]][1]))
})

test_that("[H1] unique study_ids produce no duplicate flag", {
  flags <- metaConvert:::.flag_cross_row_duplicates(
    study_id = c("A", "B", "C"),
    opts     = list(enable_cross_row = TRUE)
  )
  expect_equal(length(flags[[1]]), 0)
  expect_equal(length(flags[[2]]), 0)
  expect_equal(length(flags[[3]]), 0)
})

test_that("[H1] NA / empty study_id values are ignored (not treated as duplicates)", {
  flags <- metaConvert:::.flag_cross_row_duplicates(
    study_id = c(NA, "", NA, "A"),
    opts     = list(enable_cross_row = TRUE)
  )
  expect_true(all(vapply(flags, length, integer(1)) == 0))
})

test_that("[H1] enable_cross_row = FALSE disables the flag", {
  flags <- metaConvert:::.flag_cross_row_duplicates(
    study_id = c("X", "X"),
    opts     = list(enable_cross_row = FALSE)
  )
  expect_equal(length(flags[[1]]), 0)
  expect_equal(length(flags[[2]]), 0)
})

test_that("[H1] NULL study_id (column absent) returns no flags", {
  flags <- metaConvert:::.flag_cross_row_duplicates(
    study_id = NULL,
    opts     = list(enable_cross_row = TRUE)
  )
  expect_equal(length(flags), 0)
})

test_that("[H1] three-way duplicate names both other rows in each message", {
  flags <- metaConvert:::.flag_cross_row_duplicates(
    study_id = c("T", "T", "T"),
    opts     = list(enable_cross_row = TRUE)
  )
  expect_true(grepl("T (row 2), T (row 3)", flags[[1]], fixed = TRUE))
  expect_true(grepl("T (row 1), T (row 3)", flags[[2]], fixed = TRUE))
  expect_true(grepl("T (row 1), T (row 2)", flags[[3]], fixed = TRUE))
})

test_that("[H1] end-to-end: summary(convert_df()) emits H1 flag when study_id is duplicated", {
  df <- data.frame(
    study_id        = c("ELOQUENT-2", "ELOQUENT-1", "ELOQUENT-2"),
    n_cases_exp     = c(129, 119, 13),
    n_controls_exp  = c(189, 252, 18),
    n_exp           = c(318, 371, 31),
    n_cases_nexp    = c(81, 82, 8),
    n_controls_nexp = c(236, 289, 21),
    n_nexp          = c(317, 371, 29)
  )
  res <- summary(convert_df(df, measure = "rr"), flags = TRUE)
  expect_true(grepl("Duplicate study_id 'ELOQUENT-2'", res$flags_crude[1], fixed = TRUE))
  expect_true(grepl("Duplicate study_id 'ELOQUENT-2'", res$flags_crude[3], fixed = TRUE))
  expect_false(grepl("Duplicate study_id", res$flags_crude[2], fixed = TRUE))
})


# ==============================================================================
# Fixes from the 2026-06 code review (exp-scale handling, F1 r floor,
# negative reliability, V22 negative correlation)
# ==============================================================================

test_that("[F1-fix] high |r| with a legitimately small SE is not flagged", {
  # r = 0.9, n = 100: true SE ~ (1 - r^2)/sqrt(n - 1) = 0.019. The old fixed
  # floor 0.3/sqrt(N) = 0.03 would have flagged it as implausibly small.
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(0.9), se = c(0.02), measure = "r", n_total = c(100)
  )
  expect_false(any(grepl("Implausibly small SE", flags[[1]])))
})

test_that("[F1-fix] variance-entered-as-SE still caught for moderate r", {
  # r = 0.3, n = 100: variance = ((1-r^2)^2)/(n-1) ~ 0.0084, below the scaled
  # floor 0.3*(1-r^2)/sqrt(N) ~ 0.027.
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(0.3), se = c(0.0084), measure = "r", n_total = c(100)
  )
  expect_true(any(grepl("Implausibly small SE", flags[[1]])))
})

test_that("[V25] negative Cronbach's alpha is UNUSUAL (preserved), not INVALID", {
  dat <- data.frame(cronbach_alpha = c(-0.2, 0.8), n_sample = c(50, 50),
                    n_items = c(10, 10))
  out <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("\\[UNUSUAL\\] Negative alpha", out$issues[1]))
  expect_false(grepl("Out-of-range alpha", out$issues[1]))
  expect_equal(out$data$cronbach_alpha[1], -0.2)  # preserved
  expect_equal(out$issues[2], "")
})

test_that("[V25] ICC in [-1, 0) is UNUSUAL; ICC < -1 is INVALID", {
  dat <- data.frame(icc = c(-0.3, -1.5, 0.9), n_sample = rep(50, 3),
                    n_measurements = rep(2, 3))
  out <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("\\[UNUSUAL\\] Negative ICC", out$issues[1]))
  expect_equal(out$data$icc[1], -0.3)                       # preserved
  expect_true(grepl("Out-of-range ICC", out$issues[2]))     # impossible
  expect_true(is.na(out$data$icc[2]))                       # zapped
  # row 3 (icc = 0.9) is valid -> no V25/V11 flag, but it carries the always-on V31
  # agreement-type-SE informational note (icc_type defaults to "agreement")
  expect_false(grepl("Out-of-range|Negative ICC", out$issues[3]))
  expect_true(grepl("\\[INFO\\] ICC agreement-type SE", out$issues[3]))
})

test_that("[V22-fix] a negative cov_outcome_r counts as provided", {
  dat <- data.frame(ancova_mean_sd_pooled = c(2.5, 2.5),
                    cov_outcome_r = c(-0.5, NA),
                    n_exp = c(50, 50), n_nexp = c(50, 50))
  out <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("cov_outcome_r missing", out$issues[1]))
  expect_true(grepl("cov_outcome_r missing", out$issues[2]))
})

test_that("[D1-fix] exp=TRUE runs the IQR outlier check on the log scale", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = c(0.9, 1.0, 1.1, 1.05, 50), se = rep(0.2, 5),
    opts = opts, measure = "logor", exp = TRUE
  )
  expect_true(any(grepl("ES outlier", flags[[5]])))
  # message shows the natural-scale value with the natural-scale label
  expect_true(any(grepl("or = 50", flags[[5]])))
  expect_false(any(grepl("ES outlier", flags[[1]])))
})

test_that("[E5-fix] direction disagreement fires around null = 1 when exp=TRUE", {
  res <- data.frame(row_id = 1, es = 1.0, se = 0.2,
                    es_ci_lo = 0.7, es_ci_up = 1.4,
                    n_estimations = 2, dispersion_es = 0.1,
                    overlap_min_max = 0.9, diff_min_max = 1.5,
                    min_es_value = 0.5, max_es_value = 2.0,
                    min_info = "2x2", max_info = "user_input_crude",
                    info_used = "2x2")
  raw <- data.frame(row_id = 1)
  out <- metaConvert:::.flag_es_quality(res, measure = "logor", exp = TRUE,
                                        suffix = "", raw_data = raw,
                                        opts = metaConvert:::.default_flag_options())
  expect_true(grepl("Direction disagreement", out$flags))
})

test_that("[E3-fix] min-max difference is evaluated on the log scale when exp=TRUE", {
  # natural diff = 2.1 - 0.9 = 1.2 < diff_max_logor = 2, and log diff = 0.85:
  # neither scale should flag (the OLD code compared 1.2 on the wrong scale).
  # By contrast min=0.2, max=2.5 has natural diff 2.3 (old: flag) and log
  # diff 2.53 (new: flag on the correct scale).
  opts <- metaConvert:::.default_flag_options()
  f <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(2, 2), dispersion = c(NA, NA), overlap = c(1, 1),
    diff_mm = c(1.2, 2.3), opts = opts, measure = "logor",
    min_es = c(0.9, 0.2), max_es = c(2.1, 2.5),
    exp = TRUE, es = c(1.4, 0.7)
  )
  expect_false(any(grepl("Large min-max", f[[1]])))
  expect_true(any(grepl("Large min-max", f[[2]])))
  expect_true(any(grepl("log scale", f[[2]])))
})

test_that("[reversal-fix] systematic 1/x reversal detected on the natural scale", {
  expect_true(metaConvert:::.rows_all_reversed(
    n_estimations = c(2, 2), min_es = c(0.5, 0.25), max_es = c(2, 4),
    on_exp_scale = TRUE
  ))
  expect_false(metaConvert:::.rows_all_reversed(
    n_estimations = c(2, 2), min_es = c(0.5, 1.2), max_es = c(2, 3),
    on_exp_scale = TRUE
  ))
})

test_that("[A6-fix] user rows accept both z- and t-based CIs (construction rarely reported)", {
  # A user row whose CI matches the t-based width is NOT flagged [DISCORDANT]
  # by default: the source's CI construction is rarely reported and a t CI is
  # legitimate. (Deliberate trade-off: a paired-vs-independent SE error of the
  # same magnitude hides in the t-z gap — opt in via enable_informational.)
  se_val <- 46.24; n_val <- 16; es_val <- 100
  half_t <- stats::qt(0.975, n_val - 2) * se_val   # ~ z-width * 1.09
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - half_t), ci_up = c(es_val + half_t),
    info_used = c("user_input_crude"),
    measure = "mdw", n_total = c(n_val)
  )
  expect_false(any(grepl("\\[DISCORDANT\\] CI width inconsistent", flags[[1]])))

  # With enable_informational = TRUE the t-consistent width emits an [INFO] note.
  flags_info <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - half_t), ci_up = c(es_val + half_t),
    info_used = c("user_input_crude"),
    measure = "mdw", n_total = c(n_val),
    enable_informational = TRUE
  )
  expect_true(any(grepl("\\[INFO\\] CI width matches a t-based CI", flags_info[[1]])))

  # The same user row with a clean Wald-z CI passes silently.
  half_z <- stats::qnorm(0.975) * se_val
  flags_z <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - half_z), ci_up = c(es_val + half_z),
    info_used = c("user_input_crude"),
    measure = "mdw", n_total = c(n_val)
  )
  expect_false(any(grepl("CI width", flags_z[[1]])))

  # A width matching NEITHER convention (e.g. SD entered as SE) is still flagged.
  flags_bad <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - 3 * half_z), ci_up = c(es_val + 3 * half_z),
    info_used = c("user_input_crude"),
    measure = "mdw", n_total = c(n_val)
  )
  expect_true(any(grepl("\\[DISCORDANT\\] CI width inconsistent", flags_bad[[1]])))
})

test_that("[A6-fix] package-computed mdw rows expect t-based CIs", {
  se_val <- 0.4; n_val <- 8; es_val <- 0.5
  half_t <- stats::qt(0.975, n_val - 2) * se_val
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - half_t), ci_up = c(es_val + half_t),
    info_used = c("means_sd_pre_post_single_group"),
    measure = "mdw", n_total = c(n_val)
  )
  expect_false(any(grepl("CI width inconsistent", flags[[1]])))
})

test_that("[E5-guard] rounding-level straddles of the null are suppressed", {
  base <- data.frame(row_id = 1, se = 0.05,
                     n_estimations = 2, dispersion_es = 0.02,
                     overlap_min_max = 0.95, diff_min_max = 0.1,
                     min_info = "2x2_sum", max_info = "chisq_pval",
                     info_used = "2x2")
  raw <- data.frame(row_id = 1)
  opts <- metaConvert:::.default_flag_options()
  # RR 0.97 vs 1.05: straddles 1 but |log| < 0.05 on the min side -> suppressed
  res1 <- cbind(base, es = 1.0, es_ci_lo = 0.9, es_ci_up = 1.11,
                min_es_value = 0.971, max_es_value = 1.045)
  out1 <- metaConvert:::.flag_es_quality(res1, measure = "logrr", exp = TRUE,
                                         suffix = "", raw_data = raw, opts = opts)
  expect_false(grepl("Direction disagreement", out1$flags))
  # RR 0.5 vs 2.0: genuine disagreement -> fires
  res2 <- cbind(base, es = 1.0, es_ci_lo = 0.7, es_ci_up = 1.4,
                min_es_value = 0.5, max_es_value = 2.0)
  out2 <- metaConvert:::.flag_es_quality(res2, measure = "logrr", exp = TRUE,
                                         suffix = "", raw_data = raw, opts = opts)
  expect_true(grepl("Direction disagreement", out2$flags))
  # Additive scale: -0.02 vs +0.03 suppressed, -0.41 vs +0.39 fires
  res3 <- cbind(base, es = 0.01, es_ci_lo = -0.1, es_ci_up = 0.12,
                min_es_value = -0.02, max_es_value = 0.03)
  out3 <- metaConvert:::.flag_es_quality(res3, measure = "g", exp = FALSE,
                                         suffix = "", raw_data = raw, opts = opts)
  expect_false(grepl("Direction disagreement", out3$flags))
  res4 <- cbind(base, es = 0.1, es_ci_lo = -0.3, es_ci_up = 0.5,
                min_es_value = -0.41, max_es_value = 0.39)
  out4 <- metaConvert:::.flag_es_quality(res4, measure = "g", exp = FALSE,
                                         suffix = "", raw_data = raw, opts = opts)
  expect_true(grepl("Direction disagreement", out4$flags))
})

test_that("[A6-df] within-group measures use the paired df (n_pairs - 1)", {
  # 7 paired subjects: the package builds the mdw CI with qt(.975, 6). With
  # n_total = 14 (n_exp + n_nexp double-counting the same people), the old
  # n_total - 2 = 12 expectation produced a ~12% width gap and a false flag.
  se_val <- 32.7; n_pairs <- 7; es_val <- 86.6
  half <- stats::qt(0.975, n_pairs - 1) * se_val
  flags <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - half), ci_up = c(es_val + half),
    info_used = c("means_sd_pre_post_single_group"),
    measure = "mdw", n_total = c(14), n_exp = c(n_pairs)
  )
  expect_false(any(grepl("CI width inconsistent", flags[[1]])))

  # A genuinely wrong CI (e.g. built from an SD instead of an SE) still flags.
  flags_bad <- metaConvert:::.flag_numeric_integrity(
    es = c(es_val), se = c(se_val),
    ci_lo = c(es_val - 3 * half), ci_up = c(es_val + 3 * half),
    info_used = c("means_sd_pre_post_single_group"),
    measure = "mdw", n_total = c(14), n_exp = c(n_pairs)
  )
  expect_true(any(grepl("CI width inconsistent", flags_bad[[1]])))
})
