# Red-half regression tests for the AUDIT-2026-08-28 "unit_type / correlation-family"
# cluster. Every test below asserts the CORRECT behaviour, so each one FAILS on the
# current tree and passes once the corresponding defect is repaired.
#
# Findings covered (AUDIT-2026-08-28-findings.md):
#   632-641 + 1113-1124  unit_type default "raw_scale" overwritten with NA by missing()
#   412-421              small_margin_prop range-checked nowhere
#   599-608              es_from_beta_std sign-flips d/g on a negative sd_dv
#   1022-1033            smd_var = "hedges_olkin" leaks into the biserial r/z variance


# --- helpers: independent closed forms -------------------------------------------
# Mathur & VanderWeele (2020) r -> d, raw-units increase:
#     d = r * unit_increase_iv / (sd_iv * sqrt(1 - r^2))
.mathur_d_raw <- function(r, sd_iv, unit_increase_iv) {
  r * unit_increase_iv / (sd_iv * sqrt(1 - r^2))
}

# Bonett & Price (2005) OR -> r, as printed in the paper:
#     c = (1 - |p1. - p.1| / 5 - (0.5 - pmin)^2) / 2 ;  r = cos(pi / (1 + OR^c))
.bonett_r <- function(or, logor_se, n_exp, n_cases, n_sample, pmin) {
  cc <- (1 - abs(n_exp / n_sample - n_cases / n_sample) / 5 - (0.5 - pmin)^2) / 2
  r <- cos(pi / (1 + or^cc))
  r_se <- logor_se * (pi * cc * or^cc) * sin(pi / (1 + or^cc)) / (1 + or^cc)^2
  list(r = r, r_se = r_se)
}


# ---------------------------------------------------------------------------------
# AUDIT lines 632-641 (duplicated at 1113-1124).
# DEFECT: seven exported routes declare `unit_type = "raw_scale"` and then immediately
# run `if (missing(unit_type)) unit_type <- rep(NA, n)`; the mathur eligibility gate
# requires !is.na(unit_type), so relying on the documented default silently returns
# NA for d/d_se/g/g_se/logor instead of the raw-units conversion.
test_that("AUDIT-unit-type: the documented default reaches all seven cor_to_smd='mathur' routes", {
  sd_iv <- 2
  incr <- 1

  # --- es_from_pearson_r ---------------------------------------------------------
  p_def <- es_from_pearson_r(pearson_r = 0.4, n_sample = 100, cor_to_smd = "mathur",
                             sd_iv = sd_iv, unit_increase_iv = incr)
  p_raw <- es_from_pearson_r(pearson_r = 0.4, n_sample = 100, cor_to_smd = "mathur",
                             sd_iv = sd_iv, unit_increase_iv = incr, unit_type = "raw_scale")
  expect_false(is.na(p_def$d))
  expect_equal(p_def$d, .mathur_d_raw(0.4, sd_iv, incr))   # 0.2182178902
  expect_equal(p_def[, c("d", "d_se", "g", "g_se")], p_raw[, c("d", "d_se", "g", "g_se")])

  # --- es_from_fisher_z ----------------------------------------------------------
  f_def <- es_from_fisher_z(fisher_z = atanh(0.4), n_sample = 100, cor_to_smd = "mathur",
                            sd_iv = sd_iv, unit_increase_iv = incr)
  f_raw <- es_from_fisher_z(fisher_z = atanh(0.4), n_sample = 100, cor_to_smd = "mathur",
                            sd_iv = sd_iv, unit_increase_iv = incr, unit_type = "raw_scale")
  expect_false(is.na(f_def$d))
  expect_equal(f_def$d, .mathur_d_raw(0.4, sd_iv, incr))
  expect_equal(f_def[, c("d", "d_se", "g", "g_se")], f_raw[, c("d", "d_se", "g", "g_se")])

  # --- es_from_spearman_rho ------------------------------------------------------
  # Rupinski & Dunlap (1996): r_p = 2 * sin(pi/6 * r_s), then the Mathur raw-units step.
  rp <- 2 * sin(pi / 6 * 0.4)
  s_def <- es_from_spearman_rho(spearman_r = 0.4, n_sample = 100, cor_to_smd = "mathur",
                                sd_iv = sd_iv, unit_increase_iv = incr)
  s_raw <- es_from_spearman_rho(spearman_r = 0.4, n_sample = 100, cor_to_smd = "mathur",
                                sd_iv = sd_iv, unit_increase_iv = incr, unit_type = "raw_scale")
  expect_false(is.na(s_def$d))
  expect_equal(s_def$d, .mathur_d_raw(rp, sd_iv, incr))    # 0.2286137171
  expect_equal(s_def[, c("d", "d_se", "g", "g_se")], s_raw[, c("d", "d_se", "g", "g_se")])

  # --- es_from_linreg_t ----------------------------------------------------------
  # r = t / sqrt(t^2 + df), df = n_sample - n_covariates - 2, then Mathur raw units.
  rt <- 3 / sqrt(3^2 + (100 - 2 - 2))
  t_def <- es_from_linreg_t(linreg_t = 3, n_sample = 100, n_covariates = 2,
                            cor_to_smd = "mathur", sd_iv = sd_iv, unit_increase_iv = incr)
  t_raw <- es_from_linreg_t(linreg_t = 3, n_sample = 100, n_covariates = 2,
                            cor_to_smd = "mathur", sd_iv = sd_iv, unit_increase_iv = incr,
                            unit_type = "raw_scale")
  expect_false(is.na(t_def$d))
  expect_equal(t_def$d, .mathur_d_raw(rt, sd_iv, incr))    # 0.1530931089
  expect_equal(t_def[, c("d", "d_se", "g", "g_se")], t_raw[, c("d", "d_se", "g", "g_se")])

  # --- es_from_linreg_b_se -------------------------------------------------------
  bse_def <- es_from_linreg_b_se(linreg_b = 0.5, linreg_b_se = 0.15, n_sample = 100,
                                 n_covariates = 2, cor_to_smd = "mathur",
                                 sd_iv = sd_iv, unit_increase_iv = incr)
  bse_raw <- es_from_linreg_b_se(linreg_b = 0.5, linreg_b_se = 0.15, n_sample = 100,
                                 n_covariates = 2, cor_to_smd = "mathur",
                                 sd_iv = sd_iv, unit_increase_iv = incr,
                                 unit_type = "raw_scale")
  expect_false(is.na(bse_def$d))
  expect_equal(bse_def[, c("d", "d_se", "g", "g_se")], bse_raw[, c("d", "d_se", "g", "g_se")])

  # --- es_from_linreg_b_ci -------------------------------------------------------
  bci_def <- es_from_linreg_b_ci(linreg_b = 0.5, linreg_b_ci_lo = 0.2, linreg_b_ci_up = 0.8,
                                 n_sample = 100, n_covariates = 2, cor_to_smd = "mathur",
                                 sd_iv = sd_iv, unit_increase_iv = incr)
  bci_raw <- es_from_linreg_b_ci(linreg_b = 0.5, linreg_b_ci_lo = 0.2, linreg_b_ci_up = 0.8,
                                 n_sample = 100, n_covariates = 2, cor_to_smd = "mathur",
                                 sd_iv = sd_iv, unit_increase_iv = incr,
                                 unit_type = "raw_scale")
  expect_false(is.na(bci_def$d))
  expect_equal(bci_def[, c("d", "d_se", "g", "g_se")], bci_raw[, c("d", "d_se", "g", "g_se")])

  # --- es_from_linreg_b_pval -----------------------------------------------------
  bpv_def <- es_from_linreg_b_pval(linreg_b = 0.5, linreg_b_pval = 0.01, n_sample = 100,
                                   n_covariates = 2, cor_to_smd = "mathur",
                                   sd_iv = sd_iv, unit_increase_iv = incr)
  bpv_raw <- es_from_linreg_b_pval(linreg_b = 0.5, linreg_b_pval = 0.01, n_sample = 100,
                                   n_covariates = 2, cor_to_smd = "mathur",
                                   sd_iv = sd_iv, unit_increase_iv = incr,
                                   unit_type = "raw_scale")
  expect_false(is.na(bpv_def$d))
  expect_equal(bpv_def[, c("d", "d_se", "g", "g_se")], bpv_raw[, c("d", "d_se", "g", "g_se")])
})


# ---------------------------------------------------------------------------------
# AUDIT lines 1113-1124, **Verifier correction** (the convert_df() path the original
# finding missed).
# DEFECT: .validate_unit_type(column = TRUE) neutralises an unrecognised cell such as
# "SD" to NA while warning "Those rows are treated as raw units"; the mathur gate then
# DROPS the row, so the warning states the opposite of what happens. A cell that is
# already NA in the column IS treated as raw units (the merge loop fills it from the
# argument default), so the neutralised cell must behave the same way.
test_that("AUDIT-unit-type-convert_df: a neutralised unit_type cell is estimated as raw units, as its own warning promises", {
  d <- data.frame(
    study_id = paste0("s", 1:3),
    pearson_r = 0.4, n_sample = 100, sd_iv = 2, unit_increase_iv = 1,
    unit_type = c("raw_scale", "SD", NA),
    stringsAsFactors = FALSE
  )

  # The neutralisation warning is emitted; that part is already correct.
  expect_warning(
    convert_df(d, measure = "g", cor_to_smd = "mathur", verbose = FALSE),
    "treated as raw units"
  )

  res <- suppressWarnings(convert_df(d, measure = "g", cor_to_smd = "mathur", verbose = FALSE))
  s <- as.data.frame(suppressWarnings(summary(res)))

  ref <- s$es_crude[s$study_id == "s1"]   # explicit "raw_scale"
  neu <- s$es_crude[s$study_id == "s2"]   # "SD" -> neutralised to NA
  blk <- s$es_crude[s$study_id == "s3"]   # already NA in the column

  # Independent value: Hedges' g of the Mathur raw-units d at r = .4, sd_iv = 2, incr = 1.
  expect_false(is.na(ref))
  expect_false(is.na(blk))
  expect_equal(neu, ref)
  expect_equal(neu, blk)
})


# ---------------------------------------------------------------------------------
# AUDIT lines 412-421.
# DEFECT: small_margin_prop (documented as a proportion in (0, 0.5]) is validated
# nowhere -- absent from .positive_columns() and .bounded_columns(), and unguarded in
# .or_to_cor()'s bonett branch. A percentage (30) drives Bonett's c negative, which
# reverses r to -1 and returns a NEGATIVE standard error (-1.4e-221, i.e. a ~5e442
# inverse-variance weight); 0.9 or -0.3 give a silently wrong r with no flag at all.
test_that("AUDIT-small-margin-prop: an out-of-range small_margin_prop is rejected, never sign-reversing r or the SE", {
  # ---- direct route call: a percentage must not produce r = -1 with a negative SE ----
  bad <- es_from_or_se(or = 3, logor_se = 0.2, small_margin_prop = 30,
                       n_exp = 40, n_cases = 30, n_sample = 100, or_to_cor = "bonett")
  # A standard error is never negative. (Today: -1.4148309e-221.)
  expect_true(is.na(bad$r_se) || bad$r_se >= 0)
  # OR = 3 > 1 is a positive association; the Bonett r cannot be negative, let alone -1.
  expect_true(is.na(bad$r) || bad$r > 0)
  expect_true(is.na(bad$z_se) || is.finite(bad$z_se))

  # ---- through convert_df(): correct_inputs = TRUE must neutralise the bad cell, so
  # the margin-derived pmin is used and every row returns the SAME correct value ----
  x <- data.frame(
    study_id = paste0("s", 1:4),
    or = 3, logor_se = 0.2, n_exp = 40, n_nexp = 60,
    n_cases = 30, n_controls = 70, n_sample = 100,
    small_margin_prop = c(NA, 0.9, 30, -0.3),
    stringsAsFactors = FALSE
  )
  s <- as.data.frame(suppressWarnings(
    summary(convert_df(x, measure = "r", verbose = FALSE), flags = TRUE)
  ))

  # Independent expectation: the smallest margin is min(40, 60, 30, 70) / 100 = 0.30,
  # so c = (1 - |0.40 - 0.30|/5 - (0.5 - 0.30)^2) / 2 = 0.47.
  ok <- .bonett_r(or = 3, logor_se = 0.2, n_exp = 40, n_cases = 30,
                  n_sample = 100, pmin = 0.30)
  expect_equal(ok$r, 0.386434174938, tolerance = 1e-9)      # sanity on the closed form
  expect_equal(ok$r_se, 0.0637479990183, tolerance = 1e-9)

  expect_equal(s$es_crude, rep(ok$r, 4), tolerance = 1e-8)
  expect_equal(s$se_crude, rep(ok$r_se, 4), tolerance = 1e-8)

  # And the three out-of-range cells must be named by a Tier-1 input flag.
  fl <- paste(s$flags_crude, if ("flags_adjusted" %in% names(s)) s$flags_adjusted else "")
  expect_false(grepl("small_margin_prop", fl[1]))
  expect_true(all(grepl("small_margin_prop", fl[2:4])))
})


# ---------------------------------------------------------------------------------
# AUDIT lines 599-608.
# DEFECT: es_from_beta_std() computes unstd_beta <- beta_std * (sd_dv / sd_dummy) with
# no .positive_or_na() guard. A negative sd_dv negates unstd_beta while the pooled SD
# (which uses sd_dv^2) stays positive, so d/g come back exactly sign-reversed with a
# healthy positive SE and a correctly ordered CI -- verbatim the failure class
# R/internal_guards.R exists to close. sd_dv is also absent from .positive_columns().
test_that("AUDIT-beta-std-sd-dv: a negative sd_dv is neutralised, not silently sign-flipped", {
  n_exp <- 20; n_nexp <- 22; beta_std <- 0.35; sd_dv <- 0.98

  # Independent closed form, from the documented formula in ?es_from_beta_unstd:
  #   sd_dummy = sqrt((n_exp - n_exp^2/N) / (N - 1))
  #   unstd_beta = beta_std * sd_dv / sd_dummy
  #   sd_pooled  = sqrt((sd_dv^2 (N-1) - unstd_beta^2 n_exp n_nexp / N) / (N - 2))
  #   d = unstd_beta / sd_pooled
  N <- n_exp + n_nexp
  sd_dummy <- sqrt((n_exp - n_exp^2 / N) / (N - 1))
  ub <- beta_std * sd_dv / sd_dummy
  sd_pooled <- sqrt((sd_dv^2 * (N - 1) - ub^2 * n_exp * n_nexp / N) / (N - 2))
  d_ok <- ub / sd_pooled

  pos <- es_from_beta_std(beta_std = beta_std, sd_dv = sd_dv, n_exp = n_exp, n_nexp = n_nexp)
  expect_equal(pos$d, d_ok)                                 # control: valid input untouched
  expect_equal(pos$d, 0.730083921904, tolerance = 1e-9)

  # A negative SD is arithmetically impossible: the route must emit NA, not -d.
  neg <- es_from_beta_std(beta_std = beta_std, sd_dv = -sd_dv, n_exp = n_exp, n_nexp = n_nexp)
  expect_true(is.na(neg$d))
  expect_true(is.na(neg$d_se))
  expect_true(is.na(neg$g))
  expect_true(is.na(neg$g_se))

  # Through convert_df(): sd_dv belongs in .positive_columns(), so the cell is set to NA
  # under correct_inputs = TRUE and a Tier-1 flag names it. Nothing may be pooled.
  dd <- data.frame(study_id = "a", beta_std = beta_std, sd_dv = -sd_dv,
                   n_exp = n_exp, n_nexp = n_nexp, stringsAsFactors = FALSE)
  r <- as.data.frame(suppressWarnings(
    summary(convert_df(dd, measure = "g", verbose = FALSE), flags = TRUE)
  ))
  expect_true(is.na(r$es_crude[1]))
  expect_true(grepl("sd_dv", r$flags_crude[1]))
})


# ---------------------------------------------------------------------------------
# AUDIT lines 1022-1033 (CONFIRMED by the verifier).
# DEFECT: .smd_to_cor() normalises Soper's closed-form biserial variance by a
# hard-coded Borenstein/LS2 vd_crude while .es_from_d() builds the default d_se under
# whichever smd_var convention is in force. Under the opt-in smd_var = "hedges_olkin"
# the crude prec_ratio is 1/J^2 != 1, so r_se^2 drifts off metafor's exact RBIS
# variance and z_se^2 stops being the stabilised 1/(N-1).
test_that("AUDIT-prec-ratio: smd_var='hedges_olkin' must not leak into the biserial r/z variance", {
  n1 <- 30; n2 <- 30; dval <- 0.5
  N <- n1 + n2

  bor <- es_from_cohen_d(cohen_d = dval, n_exp = n1, n_nexp = n2, smd_var = "borenstein")
  ho  <- es_from_cohen_d(cohen_d = dval, n_exp = n1, n_nexp = n2, smd_var = "hedges_olkin")
  mf  <- metafor::escalc(measure = "RBIS", di = dval, n1i = n1, n2i = n2)

  # Point estimate is untouched by the convention (control: passes today).
  expect_equal(ho$r, as.numeric(mf$yi))
  expect_equal(bor$r, as.numeric(mf$yi))

  # Soper's biserial variance is exact and carries no smd_var convention: both arms
  # must reproduce metafor's RBIS vi (0.0227358445786558 at n = 30/30, d = 0.5).
  expect_equal(bor$r_se^2, as.numeric(mf$vi))
  expect_equal(ho$r_se^2, as.numeric(mf$vi))

  # The viechtbauer transform is variance-stabilising: z_se^2 must be 1/(N - 1).
  expect_equal(bor$z_se^2, 1 / (N - 1))
  expect_equal(ho$z_se^2, 1 / (N - 1))

  # The SMD variance itself does differ by the convention -- that is the documented
  # behaviour and must stay (control: passes today).
  expect_gt(ho$d_se^2, bor$d_se^2)
})
