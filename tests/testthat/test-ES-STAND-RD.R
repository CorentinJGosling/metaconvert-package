# ==============================================================================
# Tests for es_from_stand_RD.R: es_from_rd_se, es_from_rd_ci, es_from_rd_pval
# ==============================================================================

# --- es_from_rd_se: basic RD + NNT (always computable) ---
test_that("es_from_rd_se: RD output matches input", {
  res <- es_from_rd_se(
    rd = 0.15, rd_se = 0.05,
    n_exp = 100, n_nexp = 100
  )

  expect_equal(res$rd, 0.15, tolerance = 1e-10)
  expect_equal(res$rd_se, 0.05, tolerance = 1e-10)
  expect_equal(res$rd_ci_lo, 0.15 - qnorm(.975) * 0.05, tolerance = 1e-10)
  expect_equal(res$rd_ci_up, 0.15 + qnorm(.975) * 0.05, tolerance = 1e-10)
})

test_that("es_from_rd_se: NNT is 1/RD with correct SE", {
  rd <- 0.15
  rd_se <- 0.05
  res <- es_from_rd_se(rd = rd, rd_se = rd_se, n_exp = 100, n_nexp = 100)

  nnt_expected <- 1 / rd
  nnt_se_expected <- rd_se / rd^2

  expect_equal(res$nnt, nnt_expected, tolerance = 1e-10)
  expect_equal(res$nnt_se, nnt_se_expected, tolerance = 1e-10)
})

test_that("es_from_rd_se: without baseline_risk, OR/RR/d/g/r/z are NA", {
  res <- es_from_rd_se(
    rd = 0.15, rd_se = 0.05,
    n_exp = 100, n_nexp = 100
  )

  # NNT and RD should be available
  expect_false(is.na(res$nnt))
  expect_false(is.na(res$rd))

  # OR, RR, d, g, r, z should be NA (no baseline_risk)
  expect_true(is.na(res$logor))
  expect_true(is.na(res$logrr))
  expect_true(is.na(res$d))
  expect_true(is.na(res$g))
  expect_true(is.na(res$r))
  expect_true(is.na(res$z))
})

# --- es_from_rd_se: full conversion with baseline_risk ---
test_that("es_from_rd_se: OR conversion with baseline_risk", {
  rd <- 0.15
  rd_se <- 0.05
  br <- 0.30
  pt <- br - rd  # 0.15

  res <- es_from_rd_se(
    rd = rd, rd_se = rd_se,
    baseline_risk = br, n_exp = 100, n_nexp = 100
  )

  or_expected <- (pt * (1 - br)) / (br * (1 - pt))
  logor_expected <- log(or_expected)
  logor_se_expected <- rd_se / (pt * (1 - pt))

  expect_equal(res$logor, logor_expected, tolerance = 1e-10)
  expect_equal(res$logor_se, logor_se_expected, tolerance = 1e-10)
  expect_false(is.na(res$logor_ci_lo))
  expect_false(is.na(res$logor_ci_up))
})

test_that("es_from_rd_se: RR conversion with baseline_risk", {
  rd <- 0.15
  rd_se <- 0.05
  br <- 0.30
  pt <- br - rd  # 0.15

  res <- es_from_rd_se(
    rd = rd, rd_se = rd_se,
    baseline_risk = br, n_exp = 100, n_nexp = 100
  )

  rr_expected <- 1 - rd / br  # 0.5
  logrr_expected <- log(rr_expected)
  logrr_se_expected <- rd_se / abs(pt)

  expect_equal(res$logrr, logrr_expected, tolerance = 1e-10)
  expect_equal(res$logrr_se, logrr_se_expected, tolerance = 1e-10)
})

test_that("es_from_rd_se: d/g/r/z conversion via logOR", {
  rd <- 0.15
  rd_se <- 0.05
  br <- 0.30
  pt <- br - rd

  res <- es_from_rd_se(
    rd = rd, rd_se = rd_se,
    baseline_risk = br, n_exp = 100, n_nexp = 100
  )

  or <- (pt * (1 - br)) / (br * (1 - pt))
  logor <- log(or)
  logor_se <- rd_se / (pt * (1 - pt))

  d_expected <- logor * sqrt(3) / pi
  d_se_expected <- sqrt(logor_se^2 * 3 / pi^2)

  expect_equal(res$d, d_expected, tolerance = 1e-10)
  expect_equal(res$d_se, d_se_expected, tolerance = 1e-10)
  expect_false(is.na(res$g))
  expect_false(is.na(res$g_se))
  expect_false(is.na(res$r))
  expect_false(is.na(res$r_se))
  expect_false(is.na(res$z))
  expect_false(is.na(res$z_se))
})

# --- Edge cases ---
test_that("es_from_rd_se: RD = 0 gives NNT = NA", {
  res <- es_from_rd_se(
    rd = 0, rd_se = 0.05,
    n_exp = 100, n_nexp = 100
  )

  expect_true(is.na(res$nnt))
  expect_true(is.na(res$nnt_se))
  # RD itself should still be valid
  expect_equal(res$rd, 0, tolerance = 1e-10)
})

test_that("es_from_rd_se: treatment_risk out of (0,1) gives NA for OR/RR", {
  # rd = 0.5, baseline_risk = 0.3 -> treatment_risk = -0.2 (invalid)
  res <- es_from_rd_se(
    rd = 0.5, rd_se = 0.05,
    baseline_risk = 0.3, n_exp = 100, n_nexp = 100
  )

  # OR/RR should be NA (invalid treatment risk)
  expect_true(is.na(res$logor))
  expect_true(is.na(res$logrr))
  expect_true(is.na(res$d))

  # NNT should still work
  expect_false(is.na(res$nnt))
  expect_equal(res$nnt, 1 / 0.5, tolerance = 1e-10)
})

test_that("es_from_rd_se: baseline_risk = 0 gives RR = NA", {
  res <- es_from_rd_se(
    rd = 0.15, rd_se = 0.05,
    baseline_risk = 0, n_exp = 100, n_nexp = 100
  )

  expect_true(is.na(res$logrr))
  expect_true(is.na(res$logor))
})

test_that("es_from_rd_se: negative RD gives negative NNT", {
  res <- es_from_rd_se(
    rd = -0.10, rd_se = 0.03,
    n_exp = 100, n_nexp = 100
  )

  expect_true(res$nnt < 0)
  expect_equal(res$nnt, 1 / (-0.10), tolerance = 1e-10)
})

# --- Reversal ---
test_that("es_from_rd_se: reverse flips RD, NNT, OR, and d signs", {
  rd <- 0.15
  rd_se <- 0.05
  br <- 0.30

  res_fwd <- es_from_rd_se(rd = rd, rd_se = rd_se,
                            baseline_risk = br, n_exp = 100, n_nexp = 100,
                            reverse_rd = FALSE)
  res_rev <- es_from_rd_se(rd = rd, rd_se = rd_se,
                            baseline_risk = br, n_exp = 100, n_nexp = 100,
                            reverse_rd = TRUE)

  # Signs flip

  expect_equal(res_fwd$rd, -res_rev$rd, tolerance = 1e-10)
  expect_equal(res_fwd$nnt, -res_rev$nnt, tolerance = 1e-10)
  expect_equal(res_fwd$logor, -res_rev$logor, tolerance = 1e-10)
  expect_equal(res_fwd$d, -res_rev$d, tolerance = 1e-10)

  # SEs remain the same
  expect_equal(res_fwd$rd_se, res_rev$rd_se, tolerance = 1e-10)
  expect_equal(res_fwd$nnt_se, res_rev$nnt_se, tolerance = 1e-10)
  expect_equal(res_fwd$logor_se, res_rev$logor_se, tolerance = 1e-10)
  expect_equal(res_fwd$d_se, res_rev$d_se, tolerance = 1e-10)
})

# --- NNT CI: discontinuous when crossing zero ---
test_that("es_from_rd_se: NNT CI is NA when RD CI crosses zero", {
  # Large SE relative to RD -> CI crosses zero
  res <- es_from_rd_se(
    rd = 0.02, rd_se = 0.10,
    n_exp = 50, n_nexp = 50
  )

  rd_ci_lo <- 0.02 - qnorm(.975) * 0.10
  rd_ci_up <- 0.02 + qnorm(.975) * 0.10
  expect_true(rd_ci_lo < 0 & rd_ci_up > 0)

  expect_true(is.na(res$nnt_ci_lo))
  expect_true(is.na(res$nnt_ci_up))
  # Point estimate and SE still available
  expect_false(is.na(res$nnt))
  expect_false(is.na(res$nnt_se))
})

test_that("es_from_rd_se: NNT CI is present when RD CI does not cross zero", {
  res <- es_from_rd_se(
    rd = 0.20, rd_se = 0.03,
    n_exp = 200, n_nexp = 200
  )

  rd_ci_lo <- 0.20 - qnorm(.975) * 0.03
  rd_ci_up <- 0.20 + qnorm(.975) * 0.03
  expect_true(rd_ci_lo > 0)

  expect_false(is.na(res$nnt_ci_lo))
  expect_false(is.na(res$nnt_ci_up))
  # CI should contain the point estimate
  expect_true(res$nnt_ci_lo <= res$nnt & res$nnt <= res$nnt_ci_up)
})

# --- es_from_rd_ci ---
test_that("es_from_rd_ci: derives SE from CI correctly", {
  rd <- 0.15
  rd_se <- 0.05
  ci_lo <- rd - qnorm(.975) * rd_se
  ci_up <- rd + qnorm(.975) * rd_se

  res_ci <- es_from_rd_ci(
    rd = rd, rd_ci_lo = ci_lo, rd_ci_up = ci_up,
    baseline_risk = 0.30, n_exp = 100, n_nexp = 100
  )
  res_se <- es_from_rd_se(
    rd = rd, rd_se = rd_se,
    baseline_risk = 0.30, n_exp = 100, n_nexp = 100
  )

  # All ES values should match exactly
  expect_equal(res_ci$rd, res_se$rd, tolerance = 1e-10)
  expect_equal(res_ci$rd_se, res_se$rd_se, tolerance = 1e-10)
  expect_equal(res_ci$nnt, res_se$nnt, tolerance = 1e-10)
  expect_equal(res_ci$logor, res_se$logor, tolerance = 1e-10)
  expect_equal(res_ci$d, res_se$d, tolerance = 1e-10)
  expect_equal(res_ci$g, res_se$g, tolerance = 1e-10)
  expect_equal(res_ci$r, res_se$r, tolerance = 1e-10)

  expect_equal(res_ci$info_used, "rd_ci")
})

# --- es_from_rd_pval ---
test_that("es_from_rd_pval: derives SE from p-value correctly", {
  rd <- 0.15
  pval <- 0.01
  z_rd <- qnorm(pval / 2, lower.tail = FALSE)
  rd_se_derived <- abs(rd / z_rd)

  res_pval <- es_from_rd_pval(
    rd = rd, rd_pval = pval,
    baseline_risk = 0.30, n_exp = 100, n_nexp = 100
  )
  res_se <- es_from_rd_se(
    rd = rd, rd_se = rd_se_derived,
    baseline_risk = 0.30, n_exp = 100, n_nexp = 100
  )

  # All ES values should match
  expect_equal(res_pval$rd, res_se$rd, tolerance = 1e-10)
  expect_equal(res_pval$rd_se, res_se$rd_se, tolerance = 1e-10)
  expect_equal(res_pval$nnt, res_se$nnt, tolerance = 1e-10)
  expect_equal(res_pval$logor, res_se$logor, tolerance = 1e-10)
  expect_equal(res_pval$d, res_se$d, tolerance = 1e-10)

  expect_equal(res_pval$info_used, "rd_pval")
})

# --- Vectorized input ---
test_that("es_from_rd_se: vectorized input works", {
  res <- es_from_rd_se(
    rd = c(0.10, 0.20, -0.05),
    rd_se = c(0.03, 0.04, 0.02),
    baseline_risk = c(0.30, 0.40, 0.20),
    n_exp = c(100, 150, 80),
    n_nexp = c(100, 150, 80)
  )

  expect_equal(nrow(res), 3)
  expect_equal(res$rd, c(0.10, 0.20, -0.05), tolerance = 1e-10)
  expect_false(any(is.na(res$logor)))
  expect_false(any(is.na(res$nnt)))
})

# --- Pipeline: convert_df with measure = "rd" ---
test_that("convert_df pipeline: RD from rd+rd_se", {
  dat <- data.frame(
    rd = 0.15,
    rd_se = 0.05,
    baseline_risk = 0.30,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "rd", verbose = FALSE))

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_false(is.na(res$se_crude))
  expect_equal(res$es_crude, 0.15, tolerance = 1e-10)
  expect_equal(res$se_crude, 0.05, tolerance = 1e-10)
  expect_equal(res$info_used_crude, "rd_se")
})

test_that("convert_df pipeline: RD from rd+CI", {
  rd_se <- 0.05
  dat <- data.frame(
    rd = 0.15,
    rd_ci_lo = 0.15 - qnorm(.975) * rd_se,
    rd_ci_up = 0.15 + qnorm(.975) * rd_se,
    baseline_risk = 0.30,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "rd",
                             es_selected = "hierarchy",
                             hierarchy = "rd_ci",
                             verbose = FALSE))

  expect_equal(nrow(res), 1)
  expect_equal(res$es_crude, 0.15, tolerance = 1e-10)
  expect_equal(res$se_crude, rd_se, tolerance = 1e-4)
  expect_equal(res$info_used_crude, "rd_ci")
})

test_that("convert_df pipeline: RD from rd+pval", {
  dat <- data.frame(
    rd = 0.15,
    rd_pval = 0.01,
    baseline_risk = 0.30,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "rd",
                             es_selected = "hierarchy",
                             hierarchy = "rd_pval",
                             verbose = FALSE))

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_false(is.na(res$se_crude))
  expect_equal(res$info_used_crude, "rd_pval")
})

# --- Pipeline: NNT from rd ---
test_that("convert_df pipeline: NNT from rd+SE", {
  dat <- data.frame(
    rd = 0.15,
    rd_se = 0.05,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "nnt",
                             es_selected = "hierarchy",
                             hierarchy = "rd_se",
                             verbose = FALSE), digits = 11)

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  nnt_expected <- 1 / 0.15
  expect_equal(res$es_crude, nnt_expected, tolerance = 1e-10)
  expect_equal(res$info_used_crude, "rd_se")
})

# --- Pipeline: OR from rd ---
test_that("convert_df pipeline: OR from rd+SE (with baseline_risk)", {
  rd <- 0.15
  br <- 0.30
  pt <- br - rd

  dat <- data.frame(
    rd = rd,
    rd_se = 0.05,
    baseline_risk = br,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "or",
                             es_selected = "hierarchy",
                             hierarchy = "rd_se",
                             verbose = FALSE), digits = 11)

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_equal(res$info_used_crude, "rd_se")

  # measure="or" sets exp=TRUE, so es_crude is exp(logor) = OR
  standalone <- es_from_rd_se(rd = rd, rd_se = 0.05,
                              baseline_risk = br, n_exp = 100, n_nexp = 100)
  expect_equal(res$es_crude, exp(standalone$logor), tolerance = 1e-10)
})

# --- Pipeline: RR from rd ---
test_that("convert_df pipeline: RR from rd+SE (with baseline_risk)", {
  rd <- 0.15
  br <- 0.30

  dat <- data.frame(
    rd = rd,
    rd_se = 0.05,
    baseline_risk = br,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "rr",
                             es_selected = "hierarchy",
                             hierarchy = "rd_se",
                             verbose = FALSE), digits = 11)

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_equal(res$info_used_crude, "rd_se")

  # measure="rr" sets exp=TRUE, so es_crude is exp(logrr) = RR
  standalone <- es_from_rd_se(rd = rd, rd_se = 0.05,
                              baseline_risk = br, n_exp = 100, n_nexp = 100)
  expect_equal(res$es_crude, exp(standalone$logrr), tolerance = 1e-10)
})

# --- Pipeline: d/g from rd ---
test_that("convert_df pipeline: d from rd+SE (with baseline_risk)", {
  rd <- 0.15
  br <- 0.30
  pt <- br - rd

  dat <- data.frame(
    rd = rd,
    rd_se = 0.05,
    baseline_risk = br,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "d",
                             es_selected = "hierarchy",
                             hierarchy = "rd_se",
                             verbose = FALSE), digits = 11)

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_equal(res$info_used_crude, "rd_se")

  # Verify against standalone function
  standalone <- es_from_rd_se(rd = rd, rd_se = 0.05,
                              baseline_risk = br, n_exp = 100, n_nexp = 100)
  expect_equal(res$es_crude, standalone$d, tolerance = 1e-10)
})

test_that("convert_df pipeline: g from rd+SE (with baseline_risk)", {
  dat <- data.frame(
    rd = 0.15,
    rd_se = 0.05,
    baseline_risk = 0.30,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "g",
                             es_selected = "hierarchy",
                             hierarchy = "rd_se",
                             verbose = FALSE), digits = 11)

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_equal(res$info_used_crude, "rd_se")
})

# --- Pipeline: r/z from rd ---
test_that("convert_df pipeline: r from rd+SE (with baseline_risk)", {
  dat <- data.frame(
    rd = 0.15,
    rd_se = 0.05,
    baseline_risk = 0.30,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "r",
                             es_selected = "hierarchy",
                             hierarchy = "rd_se",
                             verbose = FALSE), digits = 11)

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$es_crude))
  expect_equal(res$info_used_crude, "rd_se")
})

# --- Round-trip consistency: 2x2 → RD matches direct RD input ---
test_that("Round-trip: 2x2 RD matches direct RD input", {
  # Compute RD from a 2x2 table
  n_ce <- 40; n_coe <- 160; n_cn <- 10; n_con <- 190
  pt <- n_ce / (n_ce + n_coe)
  pc <- n_cn / (n_cn + n_con)
  rd <- pc - pt
  rd_se <- sqrt(pt * (1 - pt) / (n_ce + n_coe) + pc * (1 - pc) / (n_cn + n_con))

  # Direct from 2x2
  res_2x2 <- es_from_2x2(
    n_cases_exp = n_ce, n_controls_exp = n_coe,
    n_cases_nexp = n_cn, n_controls_nexp = n_con
  )

  # From standalone RD input
  res_rd <- es_from_rd_se(
    rd = rd, rd_se = rd_se,
    baseline_risk = pc,
    n_exp = n_ce + n_coe, n_nexp = n_cn + n_con
  )

  # RD should match
  expect_equal(res_rd$rd, res_2x2$rd, tolerance = 1e-10)
  expect_equal(res_rd$rd_se, res_2x2$rd_se, tolerance = 1e-10)

  # NNT should match
  expect_equal(res_rd$nnt, res_2x2$nnt, tolerance = 1e-10)
  expect_equal(res_rd$nnt_se, res_2x2$nnt_se, tolerance = 1e-10)
})

# --- Hierarchy position: rd_se preferred over 2x2 when first in hierarchy ---
test_that("Hierarchy: rd_se can be selected over 2x2", {
  dat <- data.frame(
    rd = 0.15,
    rd_se = 0.05,
    n_cases_exp = 40,
    n_controls_exp = 160,
    n_cases_nexp = 10,
    n_controls_nexp = 190,
    n_exp = 200,
    n_nexp = 200
  )

  res <- summary(convert_df(dat, measure = "rd",
                             es_selected = "hierarchy",
                             hierarchy = "rd_se",
                             verbose = FALSE))

  expect_equal(res$info_used_crude, "rd_se")
})

# --- Data extraction sheet ---
test_that("data_extraction_sheet includes rd columns for measure='rd'", {
  sheet <- data_extraction_sheet(measure = "rd", extension = "data.frame")
  expect_true("rd" %in% colnames(sheet))
  expect_true("rd_se" %in% colnames(sheet))
  expect_true("rd_ci_lo" %in% colnames(sheet))
  expect_true("rd_ci_up" %in% colnames(sheet))
  expect_true("rd_pval" %in% colnames(sheet))
  expect_true("reverse_rd" %in% colnames(sheet))
})

test_that("data_extraction_sheet includes rd columns for measure='nnt'", {
  sheet <- data_extraction_sheet(measure = "nnt", extension = "data.frame")
  expect_true("rd" %in% colnames(sheet))
  expect_true("rd_se" %in% colnames(sheet))
})

test_that("data_extraction_sheet includes rd columns for measure='or'", {
  sheet <- data_extraction_sheet(measure = "or", extension = "data.frame")
  expect_true("rd" %in% colnames(sheet))
  expect_true("rd_se" %in% colnames(sheet))
})

# --- info_used values ---
test_that("info_used is correctly set for each function", {
  res_se <- es_from_rd_se(rd = 0.10, rd_se = 0.03,
                           n_exp = 100, n_nexp = 100)
  res_ci <- es_from_rd_ci(rd = 0.10,
                           rd_ci_lo = 0.04, rd_ci_up = 0.16,
                           n_exp = 100, n_nexp = 100)
  res_pval <- es_from_rd_pval(rd = 0.10, rd_pval = 0.001,
                               n_exp = 100, n_nexp = 100)

  expect_equal(res_se$info_used, "rd_se")
  expect_equal(res_ci$info_used, "rd_ci")
  expect_equal(res_pval$info_used, "rd_pval")
})

# --- Missing input → all NA ---
test_that("es_from_rd_se: NA input gives NA output", {
  res <- es_from_rd_se(rd = NA_real_, rd_se = 0.05,
                        n_exp = 100, n_nexp = 100)

  expect_true(is.na(res$rd))
  expect_true(is.na(res$nnt))
  expect_true(is.na(res$d))
})

test_that("es_from_rd_se: NA SE gives NA for everything except RD point estimate", {
  res <- es_from_rd_se(rd = 0.15, rd_se = NA_real_,
                        baseline_risk = 0.30,
                        n_exp = 100, n_nexp = 100)

  # RD point estimate preserved but SE-dependent fields are NA
  expect_true(is.na(res$rd_se))
})

# --- Guidance: near-miss detection ---
test_that("Guidance suggests rd methods when partially filled", {
  dat <- data.frame(
    rd = 0.15,
    n_exp = 100,
    n_nexp = 100
  )

  res <- summary(convert_df(dat, measure = "rd", verbose = FALSE),
                 guidance = TRUE, flags = FALSE)

  guidance_col <- if ("es_guidance_crude" %in% names(res)) "es_guidance_crude" else "es_guidance"
  # Should suggest adding rd_se (or CI or pval)
  expect_true(grepl("rd_se|rd_ci|rd_pval", res[[guidance_col]]))
})
