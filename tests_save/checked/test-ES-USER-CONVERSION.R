# =============================================================================
# PASSTHROUGH TESTS
# =============================================================================

test_that("USER passthrough: unknown type produces passthrough with warning", {
  dat <- data.frame(
    user_es_original_measure_crude = "custom_metric",
    user_es_crude = 0.5,
    user_se_crude = 0.1
  )
  expect_warning(
    res <- convert_df(dat, measure = "g"),
    "Unknown user_es_original_measure_crude"
  )
  s <- summary(res, digits = 11)
  expect_equal(s$es_crude, 0.5, tolerance = 1e-10)
  expect_equal(s$se_crude, 0.1, tolerance = 1e-10)
  expect_equal(s$info_used_crude, "user_input_crude")
})

test_that("USER passthrough: NA type (no original measure) produces passthrough", {
  dat <- data.frame(
    user_es_original_measure_crude = NA_character_,
    user_es_crude = 0.5,
    user_se_crude = 0.1
  )
  res <- convert_df(dat, measure = "g", verbose = FALSE)
  s <- summary(res, digits = 11)
  expect_equal(s$es_crude, 0.5, tolerance = 1e-10)
  expect_equal(s$info_used_crude, "user_input_crude")
})


# =============================================================================
# PIPELINE-VS-PIPELINE CONVERSION TESTS (CRUDE)
# Each test compares: direct input through convert_df vs USER input through
# convert_df. Both paths must produce identical ES/SE and correct info_used.
# =============================================================================

test_that("USER vs direct pipeline: logor -> g", {
  logor_val <- log(2.5)
  logor_se <- 0.3

  # Direct path: or + logor_se columns
  dat_direct <- data.frame(
    or = exp(logor_val), logor_se = logor_se,
    n_exp = 50, n_nexp = 50
  )
  res_direct <- convert_df(dat_direct, measure = "g", verbose = FALSE)
  s_direct <- summary(res_direct, digits = 11)

  # USER path
  dat_user <- data.frame(
    user_es_original_measure_crude = "logor",
    user_es_crude = logor_val,
    user_se_crude = logor_se,
    n_exp = 50, n_nexp = 50
  )
  res_user <- convert_df(dat_user, measure = "g", verbose = FALSE)
  s_user <- summary(res_user, digits = 11)

  # Values must match exactly
  expect_equal(s_user$es_crude, s_direct$es_crude, tolerance = 1e-10)
  expect_equal(s_user$se_crude, s_direct$se_crude, tolerance = 1e-10)

  # Info used: direct = "or_se", USER = "user_input_crude"
  expect_equal(s_direct$info_used_crude, "or_se")
  expect_equal(s_user$info_used_crude, "user_input_crude")
})

test_that("USER vs direct pipeline: natural OR + log-scale SE -> g", {
  or_val <- 2.5
  logor_se <- 0.3  # SE on log scale (SE of logOR)

  # Direct path: or + logor_se (log-scale SE)
  dat_direct <- data.frame(
    or = or_val, logor_se = logor_se,
    n_exp = 50, n_nexp = 50
  )
  res_direct <- convert_df(dat_direct, measure = "g", verbose = FALSE)
  s_direct <- summary(res_direct, digits = 11)

  # USER path: natural-scale OR + log-scale SE
  dat_user <- data.frame(
    user_es_original_measure_crude = "or",
    user_es_crude = or_val,
    user_se_crude = logor_se,
    n_exp = 50, n_nexp = 50
  )
  res_user <- convert_df(dat_user, measure = "g", verbose = FALSE)
  s_user <- summary(res_user, digits = 11)

  expect_equal(s_user$es_crude, s_direct$es_crude, tolerance = 1e-10)
  expect_equal(s_user$se_crude, s_direct$se_crude, tolerance = 1e-10)
  expect_equal(s_direct$info_used_crude, "or_se")
  expect_equal(s_user$info_used_crude, "user_input_crude")
})

test_that("USER vs direct pipeline: d -> g", {
  d_val <- 0.5
  n_e <- 50
  n_ne <- 50
  # Use the same d_se formula as .es_from_d() computes from sample sizes
  d_se <- sqrt((n_e + n_ne) / (n_e * n_ne) + d_val^2 / (2 * (n_e + n_ne)))

  # Direct path: cohen_d (SE computed from sample sizes by es_from_cohen_d)
  dat_direct <- data.frame(
    cohen_d = d_val,
    n_exp = n_e, n_nexp = n_ne
  )
  res_direct <- convert_df(dat_direct, measure = "g", verbose = FALSE)
  s_direct <- summary(res_direct, digits = 11)

  # USER path: provide the formula-matching d_se
  dat_user <- data.frame(
    user_es_original_measure_crude = "d",
    user_es_crude = d_val,
    user_se_crude = d_se,
    n_exp = n_e, n_nexp = n_ne
  )
  res_user <- convert_df(dat_user, measure = "g", verbose = FALSE)
  s_user <- summary(res_user, digits = 11)

  expect_equal(s_user$es_crude, s_direct$es_crude, tolerance = 1e-10)
  expect_equal(s_user$se_crude, s_direct$se_crude, tolerance = 1e-10)
  expect_equal(s_direct$info_used_crude, "cohen_d")
  expect_equal(s_user$info_used_crude, "user_input_crude")
})

test_that("USER vs direct pipeline: r -> g", {
  r_val <- 0.3

  # Direct path: pearson_r column
  dat_direct <- data.frame(
    pearson_r = r_val,
    n_sample = 100
  )
  res_direct <- convert_df(dat_direct, measure = "g", verbose = FALSE)
  s_direct <- summary(res_direct, digits = 11)

  # USER path: r with user SE (preserved, not recomputed from sample size)
  dat_user <- data.frame(
    user_es_original_measure_crude = "r",
    user_es_crude = r_val,
    user_se_crude = sqrt((1 - r_val^2)^2 / 99),
    n_sample = 100
  )
  res_user <- convert_df(dat_user, measure = "g", verbose = FALSE)
  s_user <- summary(res_user, digits = 11)

  expect_equal(s_user$es_crude, s_direct$es_crude, tolerance = 1e-10)
  expect_equal(s_user$se_crude, s_direct$se_crude, tolerance = 1e-10)
  expect_equal(s_direct$info_used_crude, "pearson_r")
  expect_equal(s_user$info_used_crude, "user_input_crude")
})

test_that("USER vs direct pipeline: z -> g (ES matches, SE preserved from user)", {
  z_val <- atanh(0.3)
  z_se_user <- sqrt(1 / 97)

  # Direct path: fisher_z column (recomputes SE from n)
  dat_direct <- data.frame(
    fisher_z = z_val,
    n_sample = 100
  )
  res_direct <- convert_df(dat_direct, measure = "g", verbose = FALSE)
  s_direct <- summary(res_direct, digits = 11)

  # USER path (preserves user-provided SE)
  dat_user <- data.frame(
    user_es_original_measure_crude = "z",
    user_es_crude = z_val,
    user_se_crude = z_se_user,
    n_sample = 100
  )
  res_user <- convert_df(dat_user, measure = "g", verbose = FALSE)
  s_user <- summary(res_user, digits = 11)

  # ES (point estimate) matches because both use the same r->d conversion
  expect_equal(s_user$es_crude, s_direct$es_crude, tolerance = 1e-10)
  # SE differs slightly because USER path preserves user's z_se while
  # direct pipeline recomputes r_se from sample size formula
  expect_equal(s_direct$info_used_crude, "fisher_z")
  expect_equal(s_user$info_used_crude, "user_input_crude")
})

test_that("USER vs direct pipeline: rd is not converted to g in either path", {
  rd_val <- 0.15
  rd_se <- 0.05

  # Direct path: rd + rd_se columns
  dat_direct <- data.frame(
    rd = rd_val, rd_se = rd_se,
    baseline_risk = 0.3,
    n_exp = 50, n_nexp = 50
  )
  s_direct <- summary(convert_df(dat_direct, measure = "g", verbose = FALSE), digits = 11)

  # USER path
  dat_user <- data.frame(
    user_es_original_measure_crude = "rd",
    user_es_crude = rd_val,
    user_se_crude = rd_se,
    baseline_risk = 0.3,
    n_exp = 50, n_nexp = 50
  )
  s_user <- summary(convert_df(dat_user, measure = "g", verbose = FALSE), digits = 11)

  # A risk difference is a ratio-family measure: it is not converted to an SMD
  # via EITHER entry route (the raw-column and user-input paths stay consistent).
  expect_true(is.na(s_direct$es_crude))
  expect_true(is.na(s_user$es_crude))
})

test_that("USER vs direct pipeline: logrr -> logor", {
  rr_val <- 1.5
  logrr_se <- 0.2

  # Direct path: rr + logrr_se columns
  dat_direct <- data.frame(
    rr = rr_val, logrr_se = logrr_se,
    baseline_risk = 0.2,
    n_exp = 50, n_nexp = 50,
    n_cases = 30, n_controls = 70
  )
  res_direct <- convert_df(dat_direct, measure = "logor", verbose = FALSE)
  s_direct <- summary(res_direct, digits = 11)

  # USER path
  dat_user <- data.frame(
    user_es_original_measure_crude = "logrr",
    user_es_crude = log(rr_val),
    user_se_crude = logrr_se,
    baseline_risk = 0.2,
    n_exp = 50, n_nexp = 50,
    n_cases = 30, n_controls = 70
  )
  res_user <- convert_df(dat_user, measure = "logor", verbose = FALSE)
  s_user <- summary(res_user, digits = 11)

  expect_equal(s_user$es_crude, s_direct$es_crude, tolerance = 1e-10)
  expect_equal(s_user$se_crude, s_direct$se_crude, tolerance = 1e-10)
  expect_equal(s_direct$info_used_crude, "rr_se")
  expect_equal(s_user$info_used_crude, "user_input_crude")
})

test_that("USER vs direct pipeline: natural RR + log-scale SE -> logrr", {
  rr_val <- 1.8
  logrr_se <- 0.2  # SE on log scale (SE of logRR)

  # Direct path: rr + logrr_se (log-scale SE)
  dat_direct <- data.frame(
    rr = rr_val, logrr_se = logrr_se,
    baseline_risk = 0.2,
    n_cases = 30, n_controls = 70
  )
  res_direct <- convert_df(dat_direct, measure = "logrr", verbose = FALSE)
  s_direct <- summary(res_direct, digits = 11)

  # USER path: natural-scale RR + log-scale SE
  dat_user <- data.frame(
    user_es_original_measure_crude = "rr",
    user_es_crude = rr_val,
    user_se_crude = logrr_se,
    baseline_risk = 0.2,
    n_cases = 30, n_controls = 70
  )
  res_user <- convert_df(dat_user, measure = "logrr", verbose = FALSE)
  s_user <- summary(res_user, digits = 11)

  expect_equal(s_user$es_crude, s_direct$es_crude, tolerance = 1e-10)
  expect_equal(s_user$se_crude, s_direct$se_crude, tolerance = 1e-10)
  expect_equal(s_direct$info_used_crude, "rr_se")
  expect_equal(s_user$info_used_crude, "user_input_crude")
})

test_that("USER natural-scale RR with CI (no SE): SE derived on log scale", {
  # Roorda 1998 case: RR=2.7, CI=[0.125, 58.239] on natural scale
  # The correct SE is derived from log-transformed CI:
  #   log(58.239) - log(0.125) / 3.92 = 1.567
  # NOT from natural-scale CI: (58.239 - 0.125) / 3.92 = 14.82
  rr_val <- 2.7
  ci_lo <- 0.125173636511331
  ci_up <- 58.2391005260929
  expected_logrr <- log(rr_val)
  expected_se <- (log(ci_up) - log(ci_lo)) / (2 * qnorm(.975))

  dat <- data.frame(
    user_es_original_measure_crude = "rr",
    user_es_crude = rr_val,
    user_ci_lo_crude = ci_lo,
    user_ci_up_crude = ci_up,
    n_exp = 9, n_nexp = 8
  )
  res <- convert_df(dat, measure = "logrr", verbose = FALSE)
  s <- summary(res, digits = 11)

  expect_equal(s$es_crude, expected_logrr, tolerance = 1e-4)
  expect_equal(s$se_crude, expected_se, tolerance = 1e-2)
  # SE should match 2x2 result (~1.567), NOT the wrong value (~5.491)
  expect_true(s$se_crude < 2.0)
})

test_that("USER rr/logrr are ratio-family: not converted to SMD or correlation", {
  # A user-entered risk ratio must behave like the risk-ratio pipeline measure
  # (es_from_rr_se): it converts to OR / NNT / RD, never to a d/g/r/z. This
  # mirrors the risk-difference user path.
  not_produced <- function(x, col) !(col %in% names(x)) || all(is.na(x[[col]]))
  for (om in c("rr", "logrr")) {
    val <- if (om == "rr") 1.5 else log(1.5)
    res <- es_from_user_crude(
      user_es_original_measure_crude = om,
      user_es_crude = val,
      user_se_crude = 0.2,
      baseline_risk = 0.2,
      n_exp = 50, n_nexp = 50,
      n_cases = 30, n_controls = 70
    )
    # The risk-ratio conversion itself is produced ...
    expect_false(is.na(res$logrr), info = om)
    # ... but no SMD / correlation value.
    expect_true(not_produced(res, "d"), info = om)
    expect_true(not_produced(res, "g"), info = om)
    expect_true(not_produced(res, "r"), info = om)
    expect_true(not_produced(res, "z"), info = om)
  }
})

test_that("USER rr through convert_df: g is NA but OR is still produced", {
  dat <- data.frame(
    user_es_original_measure_crude = "rr",
    user_es_crude = 1.5,
    user_se_crude = 0.2,
    baseline_risk = 0.2,
    n_exp = 50, n_nexp = 50,
    n_cases = 30, n_controls = 70
  )
  s_g <- summary(convert_df(dat, measure = "g", verbose = FALSE), digits = 11)
  expect_true(is.na(s_g$es_crude))

  s_or <- summary(convert_df(dat, measure = "or", verbose = FALSE), digits = 11)
  expect_false(is.na(s_or$es_crude))
  expect_equal(s_or$info_used_crude, "user_input_crude")
})

test_that("USER natural-scale OR with CI (no SE): SE derived on log scale", {
  or_val <- 3.0
  ci_lo <- 0.5
  ci_up <- 18.0
  expected_logor <- log(or_val)
  expected_se <- (log(ci_up) - log(ci_lo)) / (2 * qnorm(.975))

  dat <- data.frame(
    user_es_original_measure_crude = "or",
    user_es_crude = or_val,
    user_ci_lo_crude = ci_lo,
    user_ci_up_crude = ci_up,
    n_exp = 20, n_nexp = 20
  )
  res <- convert_df(dat, measure = "logor", verbose = FALSE)
  s <- summary(res, digits = 11)

  expect_equal(s$es_crude, expected_logor, tolerance = 1e-4)
  expect_equal(s$se_crude, expected_se, tolerance = 1e-2)
})


# =============================================================================
# ADJUSTED VARIANT TESTS
# The conversion formulas are scale transformations — they don't depend on
# whether the estimate is crude or adjusted. We verify this by comparing
# adjusted USER output against crude USER output with identical values.
# =============================================================================

test_that("USER adjusted vs crude: logor -> g produces identical values", {
  dat_crude <- data.frame(
    user_es_original_measure_crude = "logor",
    user_es_crude = log(2.5),
    user_se_crude = 0.3,
    n_exp = 50, n_nexp = 50
  )
  res_crude <- convert_df(dat_crude, measure = "g", verbose = FALSE)
  s_crude <- summary(res_crude, digits = 11)

  dat_adj <- data.frame(
    user_es_original_measure_adj = "logor",
    user_es_adj = log(2.5),
    user_se_adj = 0.3,
    n_exp = 50, n_nexp = 50
  )
  res_adj <- convert_df(dat_adj, measure = "g", verbose = FALSE, split_adjusted = TRUE)
  s_adj <- summary(res_adj, digits = 11)

  expect_equal(s_adj$es_adjusted, s_crude$es_crude, tolerance = 1e-10)
  expect_equal(s_adj$se_adjusted, s_crude$se_crude, tolerance = 1e-10)
  expect_equal(s_adj$info_used_adjusted, "user_input_adj")
})

test_that("USER adjusted vs crude: natural OR + log-scale SE -> g produces identical values", {
  dat_crude <- data.frame(
    user_es_original_measure_crude = "or",
    user_es_crude = 2.5,
    user_se_crude = 0.3,  # SE on log scale
    n_exp = 50, n_nexp = 50
  )
  res_crude <- convert_df(dat_crude, measure = "g", verbose = FALSE)
  s_crude <- summary(res_crude, digits = 11)

  dat_adj <- data.frame(
    user_es_original_measure_adj = "or",
    user_es_adj = 2.5,
    user_se_adj = 0.3,  # SE on log scale
    n_exp = 50, n_nexp = 50
  )
  res_adj <- convert_df(dat_adj, measure = "g", verbose = FALSE, split_adjusted = TRUE)
  s_adj <- summary(res_adj, digits = 11)

  expect_equal(s_adj$es_adjusted, s_crude$es_crude, tolerance = 1e-10)
  expect_equal(s_adj$se_adjusted, s_crude$se_crude, tolerance = 1e-10)
  expect_equal(s_adj$info_used_adjusted, "user_input_adj")
})

test_that("USER adjusted vs crude: d -> g produces identical values", {
  d_se <- sqrt(100 / (50 * 50) + 0.5^2 / 200)  # formula-matching SE

  dat_crude <- data.frame(
    user_es_original_measure_crude = "d",
    user_es_crude = 0.5,
    user_se_crude = d_se,
    n_exp = 50, n_nexp = 50
  )
  res_crude <- convert_df(dat_crude, measure = "g", verbose = FALSE)
  s_crude <- summary(res_crude, digits = 11)

  dat_adj <- data.frame(
    user_es_original_measure_adj = "d",
    user_es_adj = 0.5,
    user_se_adj = d_se,
    n_exp = 50, n_nexp = 50
  )
  res_adj <- convert_df(dat_adj, measure = "g", verbose = FALSE, split_adjusted = TRUE)
  s_adj <- summary(res_adj, digits = 11)

  expect_equal(s_adj$es_adjusted, s_crude$es_crude, tolerance = 1e-10)
  expect_equal(s_adj$se_adjusted, s_crude$se_crude, tolerance = 1e-10)
  expect_equal(s_adj$info_used_adjusted, "user_input_adj")
})

test_that("USER adjusted vs crude: r -> g produces identical values", {
  dat_crude <- data.frame(
    user_es_original_measure_crude = "r",
    user_es_crude = 0.3,
    user_se_crude = sqrt((1 - 0.3^2)^2 / 99),
    n_sample = 100
  )
  res_crude <- convert_df(dat_crude, measure = "g", verbose = FALSE)
  s_crude <- summary(res_crude, digits = 11)

  dat_adj <- data.frame(
    user_es_original_measure_adj = "r",
    user_es_adj = 0.3,
    user_se_adj = sqrt((1 - 0.3^2)^2 / 99),
    n_sample = 100
  )
  res_adj <- convert_df(dat_adj, measure = "g", verbose = FALSE, split_adjusted = TRUE)
  s_adj <- summary(res_adj, digits = 11)

  expect_equal(s_adj$es_adjusted, s_crude$es_crude, tolerance = 1e-10)
  expect_equal(s_adj$se_adjusted, s_crude$se_crude, tolerance = 1e-10)
  expect_equal(s_adj$info_used_adjusted, "user_input_adj")
})

test_that("USER adjusted vs crude: rd -> g is NA in both (removal is symmetric)", {
  dat_crude <- data.frame(
    user_es_original_measure_crude = "rd",
    user_es_crude = 0.15,
    user_se_crude = 0.05,
    baseline_risk = 0.3,
    n_exp = 50, n_nexp = 50
  )
  s_crude <- summary(convert_df(dat_crude, measure = "g", verbose = FALSE), digits = 11)

  dat_adj <- data.frame(
    user_es_original_measure_adj = "rd",
    user_es_adj = 0.15,
    user_se_adj = 0.05,
    baseline_risk = 0.3,
    n_exp = 50, n_nexp = 50
  )
  s_adj <- summary(convert_df(dat_adj, measure = "g", verbose = FALSE, split_adjusted = TRUE), digits = 11)

  # RD -> SMD is removed for both the crude and adjusted user-input scopes.
  expect_true(is.na(s_crude$es_crude))
  expect_true(is.na(s_adj$es_adjusted))
})

test_that("USER adjusted vs crude: logrr -> logor produces identical values", {
  dat_crude <- data.frame(
    user_es_original_measure_crude = "logrr",
    user_es_crude = log(1.5),
    user_se_crude = 0.2,
    baseline_risk = 0.2,
    n_exp = 50, n_nexp = 50,
    n_cases = 30, n_controls = 70
  )
  res_crude <- convert_df(dat_crude, measure = "logor", verbose = FALSE)
  s_crude <- summary(res_crude, digits = 11)

  dat_adj <- data.frame(
    user_es_original_measure_adj = "logrr",
    user_es_adj = log(1.5),
    user_se_adj = 0.2,
    baseline_risk = 0.2,
    n_exp = 50, n_nexp = 50,
    n_cases = 30, n_controls = 70
  )
  res_adj <- convert_df(dat_adj, measure = "logor", verbose = FALSE, split_adjusted = TRUE)
  s_adj <- summary(res_adj, digits = 11)

  expect_equal(s_adj$es_adjusted, s_crude$es_crude, tolerance = 1e-10)
  expect_equal(s_adj$se_adjusted, s_crude$se_crude, tolerance = 1e-10)
  expect_equal(s_adj$info_used_adjusted, "user_input_adj")
})


# =============================================================================
# FEATURE TESTS
# =============================================================================

test_that("USER conversion: logrr is ratio-family and is NOT converted to g", {
  dat <- data.frame(
    user_es_original_measure_crude = "logrr",
    user_es_crude = log(1.5),
    user_se_crude = 0.2,
    baseline_risk = 0.2,
    n_exp = 50, n_nexp = 50,
    n_cases = 30, n_controls = 70
  )
  # RR is a ratio-family measure: no standardized mean difference is produced.
  s_g <- summary(convert_df(dat, measure = "g", verbose = FALSE), digits = 11)
  expect_true(is.na(s_g$es_crude))

  # The ratio conversion RR -> OR still works.
  s_or <- summary(convert_df(dat, measure = "logor", verbose = FALSE), digits = 11)
  expect_false(is.na(s_or$es_crude))
  expect_equal(s_or$info_used_crude, "user_input_crude")
})

test_that("USER conversion: CI-only input derives ES and SE", {
  ci_lo <- log(1.5)
  ci_up <- log(4.0)
  # Derived: es = (ci_lo + ci_up)/2, se = (ci_up - ci_lo)/(2*qnorm(.975))
  derived_es <- (ci_up + ci_lo) / 2
  derived_se <- (ci_up - ci_lo) / (2 * qnorm(.975))

  # Reference: provide the derived values directly
  dat_ref <- data.frame(
    user_es_original_measure_crude = "logor",
    user_es_crude = derived_es,
    user_se_crude = derived_se,
    n_exp = 50, n_nexp = 50
  )
  res_ref <- convert_df(dat_ref, measure = "g", verbose = FALSE)
  s_ref <- summary(res_ref, digits = 11)

  # CI-only input
  dat_ci <- data.frame(
    user_es_original_measure_crude = "logor",
    user_ci_lo_crude = ci_lo,
    user_ci_up_crude = ci_up,
    n_exp = 50, n_nexp = 50
  )
  res_ci <- convert_df(dat_ci, measure = "g", verbose = FALSE)
  s_ci <- summary(res_ci, digits = 11)

  expect_equal(s_ci$es_crude, s_ref$es_crude, tolerance = 1e-10)
  expect_equal(s_ci$se_crude, s_ref$se_crude, tolerance = 1e-10)
})

test_that("USER conversion: missing companion columns produce NA gracefully", {
  dat <- data.frame(
    user_es_original_measure_crude = "d",
    user_es_crude = 0.5,
    user_se_crude = 0.15
    # No n_exp, n_nexp -> .es_from_d needs sample sizes for J correction
  )
  res <- convert_df(dat, measure = "g", verbose = FALSE)
  s <- summary(res, digits = 11)

  expect_true(is.na(s$es_crude))
})

test_that("USER conversion: mixed types per row", {
  dat <- data.frame(
    user_es_original_measure_crude = c("logor", "d", "custom_metric"),
    user_es_crude = c(log(2.5), 0.5, 0.8),
    user_se_crude = c(0.3, 0.15, 0.1),
    n_exp = c(50, 50, NA),
    n_nexp = c(50, 50, NA)
  )
  expect_warning(
    res <- convert_df(dat, measure = "g"),
    "Unknown user_es_original_measure_crude"
  )
  s <- summary(res, digits = 11)

  # Row 1: logor converted to g
  expect_false(is.na(s$es_crude[1]))
  # Row 2: d converted to g
  expect_false(is.na(s$es_crude[2]))
  # Row 3: passthrough, raw value as g
  expect_equal(s$es_crude[3], 0.8, tolerance = 1e-10)

  expect_equal(s$info_used_crude[1], "user_input_crude")
  expect_equal(s$info_used_crude[2], "user_input_crude")
  expect_equal(s$info_used_crude[3], "user_input_crude")
})

test_that("USER conversion: standalone function defaults target to g", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "logor",
    user_es_crude = log(2.5),
    user_se_crude = 0.3,
    n_exp = 50, n_nexp = 50
  )

  expect_true("g" %in% colnames(res))
  expect_false(is.na(res$g))
  expect_equal(res$info_used, "user_input_crude")
})

test_that("USER conversion: md passthrough (no pooled SD for conversion)", {
  dat <- data.frame(
    user_es_original_measure_crude = "md",
    user_es_crude = 5.0,
    user_se_crude = 1.2
  )
  res <- convert_df(dat, measure = "md", verbose = FALSE)
  s <- summary(res, digits = 11)

  expect_equal(s$es_crude, 5.0, tolerance = 1e-10)
  expect_equal(s$se_crude, 1.2, tolerance = 1e-10)
  expect_equal(s$info_used_crude, "user_input_crude")
})

test_that("USER conversion: case-insensitive type matching", {
  dat_lower <- data.frame(
    user_es_original_measure_crude = "logor",
    user_es_crude = log(2.5),
    user_se_crude = 0.3,
    n_exp = 50, n_nexp = 50
  )
  res_lower <- convert_df(dat_lower, measure = "g", verbose = FALSE)
  s_lower <- summary(res_lower, digits = 11)

  dat_mixed <- data.frame(
    user_es_original_measure_crude = "LogOR",
    user_es_crude = log(2.5),
    user_se_crude = 0.3,
    n_exp = 50, n_nexp = 50
  )
  res_mixed <- convert_df(dat_mixed, measure = "g", verbose = FALSE)
  s_mixed <- summary(res_mixed, digits = 11)

  expect_equal(s_mixed$es_crude, s_lower$es_crude, tolerance = 1e-10)
  expect_equal(s_mixed$se_crude, s_lower$se_crude, tolerance = 1e-10)
})

# =============================================================================
# SE PRESERVATION TESTS (Issue A fix)
# =============================================================================

test_that("USER r: user-provided SE is preserved, not recomputed from n", {
  r_val <- 0.5
  user_se <- 0.12  # deliberately different from formula-based SE
  formula_se <- sqrt((1 - r_val^2)^2 / (30 - 1))  # ~0.075

  res <- es_from_user_crude(
    user_es_original_measure_crude = "r",
    user_es_crude = r_val,
    user_se_crude = user_se,
    n_sample = 30
  )

  # r_se must be the user's value, not the formula-based value
  expect_equal(res$r, r_val)
  expect_equal(res$r_se, user_se)
  expect_true(abs(res$r_se - formula_se) > 0.01)  # they differ

  # z_se should be derived from user's r_se via delta method
  expected_z_se <- user_se / (1 - r_val^2)
  expect_equal(res$z_se, expected_z_se, tolerance = 1e-10)
})

test_that("USER z: user-provided SE is preserved, not recomputed from n", {
  z_val <- atanh(0.5)
  user_z_se <- 0.2  # deliberately different from 1/sqrt(n-3)
  formula_z_se <- sqrt(1 / (30 - 3))  # ~0.192

  res <- es_from_user_crude(
    user_es_original_measure_crude = "z",
    user_es_crude = z_val,
    user_se_crude = user_z_se,
    n_sample = 30
  )

  # z_se must be the user's value
  expect_equal(res$z, z_val)
  expect_equal(res$z_se, user_z_se)

  # r_se should be derived from user's z_se via inverse delta method
  r_val <- tanh(z_val)
  expected_r_se <- user_z_se * (1 - r_val^2)
  expect_equal(res$r_se, expected_r_se, tolerance = 1e-10)
})

# =============================================================================
# NNT CHAINING TEST (Issue E fix)
# =============================================================================

test_that("USER nnt: chains to es_from_rd_se for ratio conversions (OR/RR/RD), not SMD", {
  nnt_val <- 10       # RD = 0.1
  nnt_se <- 2
  baseline_risk <- 0.3  # treatment_risk = 0.3 - 0.1 = 0.2 (valid: 0 < pt < 1)

  res <- es_from_user_crude(
    user_es_original_measure_crude = "nnt",
    user_es_crude = nnt_val,
    user_se_crude = nnt_se,
    baseline_risk = baseline_risk,
    n_exp = 50, n_nexp = 50
  )

  # NNT columns preserve user values
  expect_equal(res$nnt, nnt_val)
  expect_equal(res$nnt_se, nnt_se)

  # RD should be derived
  expect_equal(res$rd, 1 / nnt_val, tolerance = 1e-10)
  expect_false(is.na(res$rd_se))

  # With baseline_risk the ratio-family conversion is produced (nnt -> rd -> OR/RR)...
  expect_false(is.na(res$logor))
  expect_false(is.na(res$logrr))
  # ... but a risk difference is not converted to an SMD or correlation. (The
  # user wrapper seeds the target-measure column as NA, so check for no VALUE
  # rather than column absence.)
  not_produced <- function(x, col) !(col %in% names(x)) || all(is.na(x[[col]]))
  expect_true(not_produced(res, "d"))
  expect_true(not_produced(res, "g"))
  expect_true(not_produced(res, "r"))
  expect_true(not_produced(res, "z"))
})

test_that("USER nnt: without baseline_risk, only rd/nnt columns produced", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "nnt",
    user_es_crude = 5,
    user_se_crude = 1
  )

  expect_equal(res$nnt, 5)
  expect_equal(res$rd, 0.2, tolerance = 1e-10)
  # Without baseline_risk, logor is NA
  expect_true(is.na(res$logor) || !("logor" %in% names(res)))
})

# =============================================================================
# BACKWARDS COMPAT
# =============================================================================

test_that("USER backwards compat: old user_es_measure_crude column works", {
  # New column name
  dat_new <- data.frame(
    user_es_original_measure_crude = "logor",
    user_es_crude = log(2.5),
    user_se_crude = 0.3,
    n_exp = 50, n_nexp = 50
  )
  res_new <- convert_df(dat_new, measure = "g", verbose = FALSE)
  s_new <- summary(res_new, digits = 11)

  # Old column name (backwards compat)
  dat_old <- data.frame(
    user_es_measure_crude = "logor",
    user_es_crude = log(2.5),
    user_se_crude = 0.3,
    n_exp = 50, n_nexp = 50
  )
  res_old <- convert_df(dat_old, measure = "g", verbose = FALSE)
  s_old <- summary(res_old, digits = 11)

  expect_equal(s_old$es_crude, s_new$es_crude, tolerance = 1e-10)
  expect_equal(s_old$se_crude, s_new$se_crude, tolerance = 1e-10)
})

# =============================================================================
# BUG FIX: NNT CI should be discontinuous, not symmetric (Altman 1998)
# =============================================================================

test_that("USER NNT: CI computed via RD inversion, not symmetric normal", {
  # NNT = 10, SE = 20 → RD = 0.1, RD_SE ≈ 0.2
  # RD CI crosses zero: [0.1 - 1.96*0.2, 0.1 + 1.96*0.2] = [-0.292, 0.492]
  # → NNT CI should be NA (discontinuous per Altman 1998)
  res <- es_from_user_crude(
    user_es_original_measure_crude = "nnt",
    user_es_crude = 10, user_se_crude = 20,
    baseline_risk = 0.3
  )
  # Point estimate and SE preserved
  expect_equal(res$nnt, 10)
  expect_equal(res$nnt_se, 20)
  # CI should be NA because underlying RD CI crosses zero
  expect_true(is.na(res$nnt_ci_lo))
  expect_true(is.na(res$nnt_ci_up))
})

test_that("USER NNT: CI is proper when RD CI does NOT cross zero", {
  # NNT = 5, SE = 1 → RD = 0.2, RD_SE = 0.04
  # RD CI: [0.2 - 1.96*0.04, 0.2 + 1.96*0.04] = [0.1216, 0.2784]
  # → NNT CI: [1/0.2784, 1/0.1216] = [3.59, 8.22] (inverted bounds)
  res <- es_from_user_crude(
    user_es_original_measure_crude = "nnt",
    user_es_crude = 5, user_se_crude = 1,
    baseline_risk = 0.3
  )
  expect_equal(res$nnt, 5)
  expect_equal(res$nnt_se, 1)
  # CI should exist (not NA) and be computed via RD inversion
  expect_false(is.na(res$nnt_ci_lo))
  expect_false(is.na(res$nnt_ci_up))
  # CI bounds from RD inversion: 1 / rd_ci_up and 1 / rd_ci_lo
  rd <- 1 / 5
  rd_se <- 1 / 5^2
  rd_ci_lo <- rd - qnorm(.975) * rd_se
  rd_ci_up <- rd + qnorm(.975) * rd_se
  expect_equal(res$nnt_ci_lo, 1 / rd_ci_up, tolerance = 1e-10)
  expect_equal(res$nnt_ci_up, 1 / rd_ci_lo, tolerance = 1e-10)
})

# =============================================================================
# BUG FIX: md mixed rows preserve SE for rows without sample sizes
# =============================================================================

test_that("USER md: SE preserved for rows without sample sizes in mixed batch", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = c("md", "md"),
    user_es_crude = c(4.0, 3.0),
    user_se_crude = c(1.2, 0.8),
    n_exp = c(30, NA),
    n_nexp = c(30, NA),
    user_es_target_measure_crude = "g"
  )
  # Row 1 (with sample sizes): full conversion
  expect_false(is.na(res$md[1]))
  expect_false(is.na(res$md_se[1]))
  expect_false(is.na(res$g[1]))

  # Row 2 (without sample sizes): md + SE should be preserved (not NA)
  expect_equal(res$md[2], 3.0)
  expect_equal(res$md_se[2], 0.8)
  expect_false(is.na(res$md_ci_lo[2]))
  expect_false(is.na(res$md_ci_up[2]))
  # d/g should be NA (can't convert without sample sizes)
  expect_true(is.na(res$g[2]))
})

# =============================================================================
# BUG FIX: g without sample sizes preserves output
# =============================================================================

test_that("USER g: output preserved when sample sizes are missing", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "g",
    user_es_crude = 0.5, user_se_crude = 0.15,
    user_es_target_measure_crude = "g"
  )
  # g and SE should NOT be NA — user provided them directly
  expect_equal(res$g, 0.5)
  expect_equal(res$g_se, 0.15)
  expect_false(is.na(res$g_ci_lo))
  expect_false(is.na(res$g_ci_up))
  # d should also be populated (via J ≈ 1 fallback)
  expect_false(is.na(res$d))
  # logOR should be populated (d → logOR doesn't need sample sizes)
  expect_false(is.na(res$logor))
})

test_that("USER g: with sample sizes, exact J factor is used (not fallback)", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "g",
    user_es_crude = 0.5, user_se_crude = 0.15,
    n_exp = 30, n_nexp = 30,
    user_es_target_measure_crude = "g"
  )
  # With proper sample sizes, J should be the exact correction factor
  df <- 30 + 30 - 2
  J <- exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
  expect_equal(res$d, 0.5 / J, tolerance = 1e-10)
  expect_equal(res$g, 0.5, tolerance = 1e-10)
})

# =============================================================================
# IMPROVEMENT: irr natural-scale support
# =============================================================================

test_that("USER irr: natural-scale IRR with CI produces logirr output", {
  # IRR = 0.7, CI = [0.5, 0.98]
  res <- es_from_user_crude(
    user_es_original_measure_crude = "irr",
    user_es_crude = 0.7,
    user_ci_lo_crude = 0.5, user_ci_up_crude = 0.98,
    user_es_target_measure_crude = "logirr"
  )
  # Should produce logirr output (log-transformed)
  expect_equal(res$logirr, log(0.7), tolerance = 1e-10)
  expected_se <- (log(0.98) - log(0.5)) / (2 * qnorm(.975))
  expect_equal(res$logirr_se, expected_se, tolerance = 1e-10)
})

test_that("USER irr: natural-scale IRR with log-scale SE produces logirr output", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "irr",
    user_es_crude = 0.7, user_se_crude = 0.2,
    user_es_target_measure_crude = "logirr"
  )
  expect_equal(res$logirr, log(0.7), tolerance = 1e-10)
  expect_equal(res$logirr_se, 0.2, tolerance = 1e-10)
})

# =============================================================================
# IMPROVEMENT: SE scale warning for ratio measures
# =============================================================================

test_that("USER or: warns when SE appears to be on natural scale", {
  expect_warning(
    es_from_user_crude(
      user_es_original_measure_crude = "or",
      user_es_crude = 2.5, user_se_crude = 3.0  # SE > ES → likely natural scale
    ),
    "LOG scale"
  )
})

test_that("USER or: no warning when SE is plausibly on log scale", {
  expect_no_warning(
    es_from_user_crude(
      user_es_original_measure_crude = "or",
      user_es_crude = 2.5, user_se_crude = 0.3  # SE < ES → likely log scale
    )
  )
})


# =============================================================================
# WITHIN-GROUP TARGET TESTS (md->mdw, d->dw, g->gw + direct mdw/dw/gw)
# Verify that user_input_crude populates the within-group output columns when
# the pipeline targets a within-group measure. Regression test for the bug
# where convert_df(..., measure="mdw") with user_es_original="md" silently
# dropped user_input_crude as a method.
# =============================================================================

test_that("USER md -> mdw: target columns populated, point preserved", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "md",
    user_es_target_measure_crude = "mdw",
    user_es_crude = 86.62,
    user_ci_lo_crude = -4.01, user_ci_up_crude = 177.25
  )
  expect_equal(res$mdw, 86.62, tolerance = 1e-6)
  expected_se <- (177.25 - (-4.01)) / (2 * qnorm(.975))
  expect_equal(res$mdw_se, expected_se, tolerance = 1e-6)
  expect_equal(res$mdw_ci_lo, -4.01, tolerance = 1e-6)
  expect_equal(res$mdw_ci_up, 177.25, tolerance = 1e-6)
  # md slot also populated (point-preserving mirror)
  expect_equal(res$md, 86.62, tolerance = 1e-6)
})

test_that("USER mdw -> mdw: direct passthrough", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "mdw",
    user_es_target_measure_crude = "mdw",
    user_es_crude = 50, user_se_crude = 10
  )
  expect_equal(res$mdw, 50, tolerance = 1e-6)
  expect_equal(res$mdw_se, 10, tolerance = 1e-6)
  expect_equal(res$md, 50, tolerance = 1e-6)
  expect_equal(res$md_se, 10, tolerance = 1e-6)
})

test_that("USER d -> dw: target columns populated", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "d",
    user_es_target_measure_crude = "dw",
    user_es_crude = 0.5, user_se_crude = 0.15,
    n_exp = 30, n_nexp = 30
  )
  expect_equal(res$dw, 0.5, tolerance = 1e-6)
  expect_equal(res$dw_se, 0.15, tolerance = 1e-6)
  expect_false(is.na(res$d))  # original d slot also kept
})

test_that("USER g -> gw: target columns populated", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "g",
    user_es_target_measure_crude = "gw",
    user_es_crude = 0.6, user_se_crude = 0.2,
    n_exp = 25, n_nexp = 25
  )
  expect_equal(res$gw, 0.6, tolerance = 1e-6)
  expect_equal(res$gw_se, 0.2, tolerance = 1e-6)
  expect_false(is.na(res$g))
})

test_that("convert_df(measure='mdw') surfaces user_input_crude as a method", {
  # Regression test: build a tiny one-row sheet with user_es_crude populated.
  # The summary should report 2 methods available (user_input_crude +
  # whichever method the team-extracted pre/post supports), not 1.
  dat <- data.frame(
    user_es_original_measure_crude = "md",
    user_es_target_measure_crude   = "mdw",
    user_es_crude    = 86.62,
    user_ci_lo_crude = -4.01,
    user_ci_up_crude = 177.25,
    mean_pre_exp = 8.10,  mean_pre_sd_exp  = 2.74,
    mean_exp     = 94.62, mean_sd_exp      = 91.91,
    n_exp        = 27
  )
  res <- convert_df(dat, measure = "mdw", verbose = FALSE)
  s <- as.data.frame(summary(res, flags = TRUE))
  # min_info/max_info should be populated (no "< 2 types..." sentinel)
  expect_false(grepl("< 2 types", s$min_info_crude[1], fixed = TRUE))
  expect_false(grepl("< 2 types", s$max_info_crude[1], fixed = TRUE))
  # One side should be user_input_crude, the other the means-pre-post path
  methods <- c(s$min_info_crude[1], s$max_info_crude[1])
  expect_true("user_input_crude" %in% methods)
  expect_true(any(grepl("means.*pre_post", methods)))
})


# =============================================================================
# HAZARD RATIO (HR) TESTS
# HR is a recognised measure but only via user_input (no raw-data methods,
# because HR estimation requires time-to-event survival data).
# =============================================================================

test_that("USER hr (CI-only) -> loghr: SE recovered on log scale", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "hr",
    user_es_target_measure_crude = "hr",
    user_es_crude = 0.75,
    user_ci_lo_crude = 0.58, user_ci_up_crude = 0.97
  )
  expect_equal(res$loghr, log(0.75), tolerance = 1e-6)
  expected_se <- (log(0.97) - log(0.58)) / (2 * qnorm(.975))
  expect_equal(res$loghr_se, expected_se, tolerance = 1e-6)
  # Step G preserves the user's CI verbatim on the entered-measure slot
  # (log-transformed for natural-scale ratio input): the source's reported
  # bounds are kept, not a symmetric rebuild around log(0.75).
  expect_equal(res$loghr_ci_lo, log(0.58), tolerance = 1e-6)
  expect_equal(res$loghr_ci_up, log(0.97), tolerance = 1e-6)
  expect_equal(res$info_used, "user_input_crude")
})

test_that("USER hr (with SE on log scale) -> loghr: ES log-transformed", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "hr",
    user_es_target_measure_crude = "hr",
    user_es_crude = 2.0,
    user_se_crude = 0.15
  )
  expect_equal(res$loghr, log(2.0), tolerance = 1e-6)
  expect_equal(res$loghr_se, 0.15, tolerance = 1e-6)
})

test_that("USER loghr -> loghr: direct passthrough on log scale", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "loghr",
    user_es_target_measure_crude = "loghr",
    user_es_crude = -0.3,
    user_se_crude = 0.1
  )
  expect_equal(res$loghr, -0.3, tolerance = 1e-6)
  expect_equal(res$loghr_se, 0.1, tolerance = 1e-6)
  expect_equal(res$loghr_ci_lo, -0.3 - qnorm(.975) * 0.1, tolerance = 1e-6)
  expect_equal(res$loghr_ci_up, -0.3 + qnorm(.975) * 0.1, tolerance = 1e-6)
})

test_that("convert_df(measure='hr') returns user_input_crude as sole method", {
  dat <- data.frame(
    user_es_original_measure_crude = "hr",
    user_es_target_measure_crude   = "hr",
    user_es_crude    = 0.75,
    user_ci_lo_crude = 0.58,
    user_ci_up_crude = 0.97
  )
  res <- convert_df(dat, measure = "hr", verbose = FALSE)
  s <- as.data.frame(summary(res, flags = TRUE))
  expect_equal(s$info_used_crude, "user_input_crude")
  # exp=TRUE for measure='hr' -> back-transformed display
  expect_equal(s$es_crude, 0.75, tolerance = 1e-4)
  # CI should be the input CI back-transformed from log scale
  expect_equal(s$es_ci_lo_crude, 0.58, tolerance = 1e-3)
  expect_equal(s$es_ci_up_crude, 0.97, tolerance = 1e-3)
})


# ==============================================================================
# Step G: the user's CI is preserved verbatim on the entered-measure slot
# (2026-06 fix: converters' t-rebuilt CIs silently replaced the source CI and
# manufactured a z-recover/t-rebuild SE-CI inconsistency that A6 then flagged)
# ==============================================================================

test_that("[StepG] user md + CI: output CI is the user CI, not a t-rebuild", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "md",
    user_es_crude = 86.62,
    user_ci_lo_crude = -4.01, user_ci_up_crude = 177.25,
    user_es_target_measure_crude = "mdw",
    n_exp = 7, n_nexp = 7
  )
  expect_equal(res$md_ci_lo[1], -4.01)
  expect_equal(res$md_ci_up[1], 177.25)
  # mirrored within-group twin carries the same user CI
  expect_equal(res$mdw_ci_lo[1], -4.01)
  expect_equal(res$mdw_ci_up[1], 177.25)
  # SE recovered with Wald-z from the same CI -> slot is internally coherent
  expect_equal(res$md_se[1], (177.25 - (-4.01)) / (2 * qnorm(.975)), tolerance = 1e-10)
})

test_that("[StepG] user d + SE only: CI is z-based, not the converter's t-rebuild", {
  d <- 0.5; se <- 0.2
  res <- es_from_user_crude(
    user_es_original_measure_crude = "d",
    user_es_crude = d, user_se_crude = se,
    user_es_target_measure_crude = "g",
    n_exp = 10, n_nexp = 10
  )
  expect_equal(res$d_ci_lo[1], d - qnorm(.975) * se, tolerance = 1e-10)
  expect_equal(res$d_ci_up[1], d + qnorm(.975) * se, tolerance = 1e-10)
  expect_equal(res$dw_ci_lo[1], d - qnorm(.975) * se, tolerance = 1e-10)
})

test_that("[StepG] natural OR with SE and CI: logor slot keeps log(user CI)", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "or",
    user_es_crude = 2.5, user_se_crude = 0.3,
    user_ci_lo_crude = 1.4, user_ci_up_crude = 4.5,
    user_es_target_measure_crude = "or",
    n_exp = 50, n_nexp = 50
  )
  expect_equal(res$logor_ci_lo[1], log(1.4), tolerance = 1e-10)
  expect_equal(res$logor_ci_up[1], log(4.5), tolerance = 1e-10)
  expect_equal(res$logor_se[1], 0.3)
})

test_that("[StepG] natural OR with SE only: converter CI untouched (no scale mixing)", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "or",
    user_es_crude = 2.5, user_se_crude = 0.3,
    user_es_target_measure_crude = "or",
    n_exp = 50, n_nexp = 50
  )
  expect_equal(res$logor_ci_lo[1], log(2.5) - qnorm(.975) * 0.3, tolerance = 1e-6)
  expect_equal(res$logor_ci_up[1], log(2.5) + qnorm(.975) * 0.3, tolerance = 1e-6)
})

test_that("[StepG] NNT + SE only keeps the RD-inverted CI (no symmetric overwrite)", {
  res <- es_from_user_crude(
    user_es_original_measure_crude = "nnt",
    user_es_crude = 8, user_se_crude = 2,
    user_es_target_measure_crude = "nnt",
    baseline_risk = 0.3
  )
  symmetric_lo <- 8 - qnorm(.975) * 2
  expect_false(isTRUE(all.equal(res$nnt_ci_lo[1], symmetric_lo)))
})

test_that("[StepG] end-to-end: user md CI survives convert_df + summary unchanged", {
  df <- data.frame(
    study_id = c("Elzanaty 2017", "Park 2024"),
    user_es_original_measure_crude = "md",
    user_es_crude = c(86.62, 207.00),
    user_ci_lo_crude = c(-4.01, 121.91),
    user_ci_up_crude = c(177.25, 292.09),
    n_exp = c(7, 5), n_nexp = c(7, 5)
  )
  res <- summary(convert_df(df, measure = "mdw", verbose = FALSE),
                 flags = TRUE)
  expect_equal(as.numeric(res$es_ci_lo_crude), c(-4.01, 121.91), tolerance = 1e-3)
  expect_equal(as.numeric(res$es_ci_up_crude), c(177.25, 292.09), tolerance = 1e-3)
  # the z-coherent slot no longer trips the A6 round-trip artifact
  expect_false(any(grepl("CI width inconsistent", res$flags_crude)))
})
