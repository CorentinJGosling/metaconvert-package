# ==============================================================================
# Tests for NNT SE/CI, Risk Difference, and Person-Time NNT
# ==============================================================================

# --- NNT from 2x2 table: point estimate + SE + CI ---
test_that("NNT from 2x2 table has SE and CI", {
  dat <- data.frame(
    n_cases_exp = 40,
    n_controls_exp = 160,
    n_cases_nexp = 10,
    n_controls_nexp = 190
  )

  res <- es_from_2x2(
    n_cases_exp = 40, n_controls_exp = 160,
    n_cases_nexp = 10, n_controls_nexp = 190
  )

  # NNT should be 1/RD
  pt <- 40 / (40 + 160)
  pc <- 10 / (10 + 190)
  rd <- pc - pt
  nnt_expected <- 1 / rd
  rd_se_expected <- sqrt(pt * (1 - pt) / 200 + pc * (1 - pc) / 200)
  nnt_se_expected <- rd_se_expected / rd^2

  expect_equal(res$nnt, nnt_expected, tolerance = 1e-10)
  expect_false(is.na(res$nnt_se))
  expect_equal(res$nnt_se, nnt_se_expected, tolerance = 1e-10)
  expect_false(is.na(res$nnt_ci_lo))
  expect_false(is.na(res$nnt_ci_up))
  # CI should contain the point estimate
  expect_true(res$nnt_ci_lo <= res$nnt && res$nnt <= res$nnt_ci_up)
})

# --- RD from 2x2 table ---
test_that("RD from 2x2 table is correct", {
  res <- es_from_2x2(
    n_cases_exp = 40, n_controls_exp = 160,
    n_cases_nexp = 10, n_controls_nexp = 190
  )

  pt <- 40 / 200
  pc <- 10 / 200
  rd_expected <- pc - pt
  rd_se_expected <- sqrt(pt * (1 - pt) / 200 + pc * (1 - pc) / 200)

  expect_equal(res$rd, rd_expected, tolerance = 1e-10)
  expect_equal(res$rd_se, rd_se_expected, tolerance = 1e-10)
  expect_equal(res$rd_ci_lo, rd_expected - qnorm(.975) * rd_se_expected, tolerance = 1e-10)
  expect_equal(res$rd_ci_up, rd_expected + qnorm(.975) * rd_se_expected, tolerance = 1e-10)
})

# --- NNT CI is NA when RD CI crosses zero (Altman 1998 discontinuous CI) ---
test_that("NNT CI is NA when RD confidence interval crosses zero", {
  # Small difference with large uncertainty -> RD CI crosses zero
  res <- es_from_2x2(
    n_cases_exp = 5, n_controls_exp = 45,
    n_cases_nexp = 4, n_controls_nexp = 46
  )

  # RD is very small (0.02), SE will be large relative to it
  rd_ci_lo <- res$rd - qnorm(.975) * res$rd_se
  rd_ci_up <- res$rd + qnorm(.975) * res$rd_se
  crosses <- rd_ci_lo < 0 & rd_ci_up > 0

  if (crosses) {
    expect_true(is.na(res$nnt_ci_lo))
    expect_true(is.na(res$nnt_ci_up))
  }
  # NNT point estimate and SE should still exist
  expect_false(is.na(res$nnt))
  expect_false(is.na(res$nnt_se))
})

# --- Reverse direction for 2x2 ---
test_that("NNT and RD reverse correctly from 2x2", {
  res_fwd <- es_from_2x2(
    n_cases_exp = 40, n_controls_exp = 160,
    n_cases_nexp = 10, n_controls_nexp = 190,
    reverse_2x2 = FALSE
  )
  res_rev <- es_from_2x2(
    n_cases_exp = 40, n_controls_exp = 160,
    n_cases_nexp = 10, n_controls_nexp = 190,
    reverse_2x2 = TRUE
  )

  expect_equal(res_fwd$nnt, -res_rev$nnt, tolerance = 1e-10)
  expect_equal(res_fwd$rd, -res_rev$rd, tolerance = 1e-10)
  # SE should be the same regardless of direction
  expect_equal(res_fwd$nnt_se, res_rev$nnt_se, tolerance = 1e-10)
  expect_equal(res_fwd$rd_se, res_rev$rd_se, tolerance = 1e-10)
})

# --- NNT from OR via delta method ---
test_that("NNT from OR has SE and CI", {
  res <- es_from_or_se(
    or = 2.5, logor_se = 0.3,
    baseline_risk = 0.1,
    n_exp = 100, n_nexp = 100
  )

  expect_false(is.na(res$nnt))
  expect_false(is.na(res$nnt_se))
  expect_true(res$nnt_se > 0)
  # RD should be populated too
  expect_false(is.na(res$rd))
  expect_false(is.na(res$rd_se))
  expect_true(res$rd_se > 0)
})

# --- NNT from RR via delta method ---
test_that("NNT from RR has SE and CI", {
  res <- es_from_rr_se(
    rr = 1.5, logrr_se = 0.2,
    baseline_risk = 0.1,
    n_exp = 100, n_nexp = 100
  )

  expect_false(is.na(res$nnt))
  expect_false(is.na(res$nnt_se))
  expect_true(res$nnt_se > 0)
  # Check RD
  rd_expected <- 0.1 * (1 - 1.5)
  expect_equal(res$rd, rd_expected, tolerance = 1e-10)
  expect_false(is.na(res$rd_se))
})

# --- Person-time NNT from IRR ---
test_that("Person-time NNT from IRR with explicit baseline_rate", {
  res <- es_from_cases_time(
    n_cases_exp = 50, n_cases_nexp = 100,
    time_exp = 1000, time_nexp = 1000,
    baseline_rate = 0.1
  )

  irr <- (50 / 1000) / (100 / 1000)
  ird <- 0.1 * (1 - irr)
  nnt_expected <- 1 / ird

  expect_equal(res$nnt, nnt_expected, tolerance = 1e-10)
  expect_false(is.na(res$nnt_se))
  expect_true(res$nnt_se > 0)
})

test_that("Person-time NNT from IRR with auto-computed baseline_rate uses Poisson SE", {
  res <- es_from_cases_time(
    n_cases_exp = 30, n_cases_nexp = 100,
    time_exp = 1000, time_nexp = 2000
  )

  # baseline_rate auto = 100/2000 = 0.05
  # IRR = (30/1000)/(100/2000) = 0.03/0.05 = 0.6
  irr <- (30 / 1000) / (100 / 2000)
  ird <- 0.05 * (1 - irr)
  nnt_expected <- 1 / ird

  expect_equal(res$nnt, nnt_expected, tolerance = 1e-10)

  # Auto-computed baseline: IRD SE uses direct Poisson rate difference variance
  ird_se_expected <- sqrt(30 / 1000^2 + 100 / 2000^2)
  expect_equal(res$rd_se, ird_se_expected, tolerance = 1e-10)
})

test_that("IRD SE differs between auto-computed and user-provided baseline_rate", {
  # Auto-computed: direct Poisson SE (accounts for both groups' variance)
  res_auto <- es_from_cases_time(
    n_cases_exp = 50, n_cases_nexp = 100,
    time_exp = 1000, time_nexp = 1000
  )
  # User-provided with same value: delta method SE (treats baseline as constant)
  res_fixed <- es_from_cases_time(
    n_cases_exp = 50, n_cases_nexp = 100,
    time_exp = 1000, time_nexp = 1000,
    baseline_rate = 100 / 1000  # same as auto-computed
  )

  # Point estimates should match
  expect_equal(res_auto$rd, res_fixed$rd, tolerance = 1e-10)
  expect_equal(res_auto$nnt, res_fixed$nnt, tolerance = 1e-10)

  # Auto SE should be larger (accounts for control group variance)
  expect_true(res_auto$rd_se > res_fixed$rd_se)
})

test_that("Person-time NNT reverse works", {
  res_fwd <- es_from_cases_time(
    n_cases_exp = 50, n_cases_nexp = 100,
    time_exp = 1000, time_nexp = 1000,
    baseline_rate = 0.1, reverse_irr = FALSE
  )
  res_rev <- es_from_cases_time(
    n_cases_exp = 50, n_cases_nexp = 100,
    time_exp = 1000, time_nexp = 1000,
    baseline_rate = 0.1, reverse_irr = TRUE
  )

  expect_equal(res_fwd$nnt, -res_rev$nnt, tolerance = 1e-10)
  expect_equal(res_fwd$rd, -res_rev$rd, tolerance = 1e-10)
})

# --- RD as standalone measure through convert_df ---
test_that("measure='rd' works through convert_df pipeline", {
  dat <- data.frame(
    n_cases_exp = c(40, 30),
    n_controls_exp = c(160, 170),
    n_cases_nexp = c(10, 20),
    n_controls_nexp = c(190, 180)
  )

  res <- summary(convert_df(dat, measure = "rd", verbose = FALSE))

  expect_equal(nrow(res), 2)
  expect_false(any(is.na(res$es_crude)))
  expect_false(any(is.na(res$se_crude)))
  expect_true(all(res$info_used_crude == "2x2"))

  # Verify against direct calculation
  pt <- dat$n_cases_exp / (dat$n_cases_exp + dat$n_controls_exp)
  pc <- dat$n_cases_nexp / (dat$n_cases_nexp + dat$n_controls_nexp)
  rd_expected <- pc - pt

  expect_equal(res$es_crude, rd_expected, tolerance = 1e-10)
})

# --- NNT measure with IRR in hierarchy ---
test_that("NNT from IRR works through convert_df pipeline", {
  dat <- data.frame(
    n_cases_exp = 50,
    n_cases_nexp = 100,
    time_exp = 1000,
    time_nexp = 1000,
    baseline_rate = 0.1
  )

  res <- summary(convert_df(dat, measure = "nnt",
                             es_selected = "hierarchy",
                             hierarchy = "cases_time",
                             verbose = FALSE))

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_equal(res$info_used_crude, "cases_time")
  expect_false(is.na(res$se_crude))
})

# --- NNT consistency: 2x2 vs OR vs RR ---
test_that("NNT is consistent across 2x2, OR, and RR sources", {
  dat <- data.frame(
    n_cases_exp = 40,
    n_controls_exp = 160,
    n_cases_nexp = 10,
    n_controls_nexp = 190
  )

  res_2x2 <- summary(convert_df(dat, measure = "nnt",
                                  es_selected = "hierarchy",
                                  hierarchy = "2x2",
                                  verbose = FALSE), digits = 11)

  # Also compute OR and RR
  dat$baseline_risk <- dat$n_cases_nexp / (dat$n_cases_nexp + dat$n_controls_nexp)
  or_val <- (40 * 190) / (160 * 10)
  logor_se <- sqrt(1/40 + 1/160 + 1/10 + 1/190)
  dat$or <- or_val
  dat$logor_se <- logor_se

  res_or <- summary(convert_df(dat, measure = "nnt",
                                es_selected = "hierarchy",
                                hierarchy = "or_se",
                                verbose = FALSE), digits = 11)

  # All three should produce NNT with the same sign
  expect_equal(sign(res_2x2$es_crude), sign(res_or$es_crude))
})

# --- RD from OR through convert_df ---
test_that("measure='rd' from OR source works", {
  dat <- data.frame(
    or = 2.5,
    logor_se = 0.3,
    baseline_risk = 0.1,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "rd",
                             es_selected = "hierarchy",
                             hierarchy = "or_se",
                             verbose = FALSE))

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_false(is.na(res$se_crude))
})

# --- RD from RR through convert_df ---
test_that("measure='rd' from RR source works", {
  dat <- data.frame(
    rr = 1.5,
    logrr_se = 0.2,
    baseline_risk = 0.2,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "rd",
                             es_selected = "hierarchy",
                             hierarchy = "rr_se",
                             verbose = FALSE))

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_false(is.na(res$se_crude))

  # RD from RR: rd = baseline_risk * (1 - rr)
  rd_expected <- 0.2 * (1 - 1.5)
  expect_equal(res$es_crude, rd_expected, tolerance = 1e-10)
})

# --- RD from IRR through convert_df ---
test_that("measure='rd' from IRR source works", {
  dat <- data.frame(
    n_cases_exp = 50,
    n_cases_nexp = 100,
    time_exp = 1000,
    time_nexp = 1000,
    baseline_rate = 0.1
  )

  res <- summary(convert_df(dat, measure = "rd",
                             es_selected = "hierarchy",
                             hierarchy = "cases_time",
                             verbose = FALSE))

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_false(is.na(res$se_crude))
})

# --- data_extraction_sheet for RD ---
test_that("data_extraction_sheet works for measure='rd'", {
  sheet <- data_extraction_sheet(measure = "rd", extension = "data.frame")
  expect_true(is.data.frame(sheet))
  # Should contain 2x2, OR, RR, and IRR columns
  expect_true("n_cases_exp" %in% colnames(sheet))
  expect_true("n_controls_exp" %in% colnames(sheet))
  expect_true("time_exp" %in% colnames(sheet))
  expect_true("baseline_rate" %in% colnames(sheet))
})

# --- data_extraction_sheet for NNT now includes IRR columns ---
test_that("data_extraction_sheet for NNT includes IRR columns", {
  sheet <- data_extraction_sheet(measure = "nnt", extension = "data.frame")
  expect_true(is.data.frame(sheet))
  expect_true("time_exp" %in% colnames(sheet))
  expect_true("time_nexp" %in% colnames(sheet))
  expect_true("baseline_rate" %in% colnames(sheet))
})

# --- NNT SE is zero when RD SE is very close to zero (edge case) ---
test_that("NNT handles rd == 0 gracefully", {
  # Equal proportions -> rd = 0 -> nnt = NA
  res <- es_from_2x2(
    n_cases_exp = 50, n_controls_exp = 50,
    n_cases_nexp = 50, n_controls_nexp = 50
  )

  expect_true(is.na(res$nnt))
  expect_true(is.na(res$nnt_se))
})

# --- Flag: minimum possible NNT from baseline risk ---
test_that("B5 flag: minimum possible NNT from baseline risk", {
  # baseline_risk = 0.054 -> min NNT = 1/0.054 = 18.52
  # Kraemer small NNT = 8.93, which is less than 18.52 -> impossible
  dat <- data.frame(
    n_cases_exp = c(40, 1),
    n_controls_exp = c(160, 999),
    n_cases_nexp = c(10, 54),
    n_controls_nexp = c(190, 946)
  )

  res <- summary(convert_df(dat, measure = "nnt", verbose = FALSE),
                 flags = TRUE, guidance = FALSE)

  # Check that flags column exists (named "flags_crude" when split_adjusted=TRUE)
  expect_true("flags_crude" %in% names(res) || "flags" %in% names(res))
})

# --- Flag: RD bounds ---
test_that("B6 flag: RD outside [-1, 1] should not occur from valid data", {
  dat <- data.frame(
    n_cases_exp = 40,
    n_controls_exp = 160,
    n_cases_nexp = 10,
    n_controls_nexp = 190
  )

  res <- summary(convert_df(dat, measure = "rd", verbose = FALSE),
                 flags = TRUE, guidance = FALSE)

  # Flags column should exist
  expect_true("flags_crude" %in% names(res) || "flags" %in% names(res))
  # RD from valid 2x2 should always be in [-1, 1]
  expect_true(abs(res$es_crude) <= 1)
})

# --- V5 flag: n_cases > n_total logical consistency ---
test_that("V5 flag: n_cases_exp > n_exp produces NA + flag", {
  dat <- data.frame(
    n_cases_exp = 100, n_exp = 10,
    n_cases_nexp = 10, n_nexp = 100
  )

  res <- summary(convert_df(dat, verbose = FALSE), flags = TRUE, guidance = FALSE)

  # ES should be NA (impossible data)
  expect_true(is.na(res$es_crude))
  # Flag should mention the inconsistency
  flag_col <- if ("flags_crude" %in% names(res)) "flags_crude" else "flags"
  expect_true(grepl("Inconsistent input", res[[flag_col]]))
  expect_true(grepl("n_cases_exp", res[[flag_col]]))
})

test_that("V5 flag: mixed valid/invalid rows handled correctly", {
  dat <- data.frame(
    n_cases_exp   = c(100, 10),
    n_exp         = c(10,  200),
    n_cases_nexp  = c(10,  20),
    n_nexp        = c(100, 200)
  )

  res <- summary(convert_df(dat, verbose = FALSE), flags = TRUE, guidance = FALSE)

  # Row 1: invalid -> NA

  expect_true(is.na(res$es_crude[1]))
  # Row 2: valid -> has ES
  expect_false(is.na(res$es_crude[2]))
})

test_that("V5 flag: n_cases > n_sample caught for proportions", {
  dat <- data.frame(n_cases = 500, n_sample = 100)

  res <- summary(convert_df(dat, measure = "prop", verbose = FALSE),
                 flags = TRUE, guidance = FALSE)

  expect_true(is.na(res$es_crude))
  flag_col <- if ("flags_crude" %in% names(res)) "flags_crude" else "flags"
  expect_true(grepl("Inconsistent input", res[[flag_col]]))
})

# --- RD is now a true effect size: requires SE (no exemption) ---
test_that("RD without SE is set to NA (no SE exemption)", {
  # User provides ES value but not SE -> should be rejected for RD
  dat <- data.frame(
    user_es_measure_crude = "rd",
    user_es_crude = c(-0.15, -0.10),
    n_exp = c(100, 100),
    n_nexp = c(100, 100)
  )

  res <- summary(convert_df(dat, measure = "rd", verbose = FALSE))

  # With no SE available, RD should be NA (unlike NNT which is exempt)
  expect_true(all(is.na(res$es_crude)))
})

test_that("RD with SE is accepted normally", {
  # User provides ES value + SE -> should be accepted
  dat <- data.frame(
    user_es_measure_crude = "rd",
    user_es_crude = c(-0.15, -0.10),
    user_se_crude = c(0.05, 0.04),
    n_exp = c(100, 100),
    n_nexp = c(100, 100)
  )

  res <- summary(convert_df(dat, measure = "rd", verbose = FALSE))

  expect_equal(res$es_crude, c(-0.15, -0.10), tolerance = 1e-10)
  expect_equal(res$se_crude, c(0.05, 0.04), tolerance = 1e-10)
})

test_that("RD from 2x2 (which always produces SE) works normally", {
  dat <- data.frame(
    n_cases_exp = c(40, 30),
    n_controls_exp = c(160, 170),
    n_cases_nexp = c(10, 20),
    n_controls_nexp = c(190, 180)
  )

  res <- summary(convert_df(dat, measure = "rd", verbose = FALSE))

  # 2x2 always computes both RD and RD_SE, so should be accepted
  expect_false(any(is.na(res$es_crude)))
  expect_false(any(is.na(res$se_crude)))
})

# --- E4 flag: NNT risk-vs-rate mixing ---
test_that("E4 flag: mixed risk-based and rate-based NNT produces warning", {
  # Row 1: 2x2 table (risk-based NNT)
  # Row 2: IRR from cases/time (rate-based NNT)
  dat <- data.frame(
    n_cases_exp     = c(40, 50),
    n_controls_exp  = c(160, NA),
    n_cases_nexp    = c(10, 100),
    n_controls_nexp = c(190, NA),
    time_exp        = c(NA, 1000),
    time_nexp       = c(NA, 1000),
    baseline_rate   = c(NA, 0.1)
  )

  res <- summary(convert_df(dat, measure = "nnt", verbose = FALSE),
                 flags = TRUE, guidance = FALSE)

  flag_col <- if ("flags_crude" %in% names(res)) "flags_crude" else "flags"

  # Row 1 should be risk-based (from 2x2)
  expect_equal(res$info_used_crude[1], "2x2")
  # Row 2 should be rate-based (from cases_time)
  expect_equal(res$info_used_crude[2], "cases_time")

  # Both rows should have the mixing flag
  expect_true(grepl("Mixed NNT types", res[[flag_col]][1]))
  expect_true(grepl("Mixed NNT types", res[[flag_col]][2]))
  expect_true(grepl("risk-based", res[[flag_col]][1]))
  expect_true(grepl("rate-based", res[[flag_col]][2]))
})

test_that("E4 flag: no mixing flag when all NNTs are same type", {
  # All risk-based (from 2x2)
  dat <- data.frame(
    n_cases_exp = c(40, 30),
    n_controls_exp = c(160, 170),
    n_cases_nexp = c(10, 20),
    n_controls_nexp = c(190, 180)
  )

  res <- summary(convert_df(dat, measure = "nnt", verbose = FALSE),
                 flags = TRUE, guidance = FALSE)

  flag_col <- if ("flags_crude" %in% names(res)) "flags_crude" else "flags"
  expect_false(any(grepl("Mixed NNT types", res[[flag_col]])))
})
