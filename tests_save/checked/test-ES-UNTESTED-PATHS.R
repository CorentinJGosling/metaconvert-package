# ==============================================================================
# Smoke tests for calculation paths that previously lacked a dedicated test file.
#
# Each test drives convert_df() + summary() with the minimal columns required to
# exercise one orphan es_from_*() function, and asserts that the resulting ES and
# SE are finite. These are *not* formula-correctness tests — the intent is to
# catch silent regressions (NA output, errors, stop() from the method counter).
# Formula correctness for most of these is already covered indirectly by
# test-INTERNAL-FULL-LIFECYCLE.R and test-formula-verification.R.
#
# Also includes:
#   - correct_inputs = FALSE preservation + Tier-1 flag merging (previously
#     untested behaviour for the validation toggle added in the flag work).
#   - stand_RR delta chain: rr + logrr_se must produce NNT/RD with non-NA SE.
#   - stand_IRR person-time NNT: cases_time inputs must yield rate-based NNT.
# ==============================================================================

library(testthat)
library(metaConvert)

# Helper: pull the single computed value for one method out of the metaConvert
# object by name (avoids depending on hierarchy selection).
.method_row <- function(mc_obj, method_name, col) {
  stopifnot(inherits(mc_obj, "metaConvert"))
  df <- mc_obj[[method_name]]
  if (is.null(df) || !col %in% colnames(df)) return(NA_real_)
  df[[col]]
}

# ==============================================================================
# Orphan path 1: es_from_chisq / es_from_phi (phi_chisq_list_L8)
# ==============================================================================
test_that("chi-square and phi inputs produce finite SMD", {
  dat <- data.frame(
    study_id = "chi1",
    chisq = 10.0,
    phi = 0.3,
    n_sample = 100,
    n_cases = 50,
    n_exp = 50,
    n_nexp = 50
  )
  res <- convert_df(dat, measure = "g", verbose = FALSE)
  g_chisq <- .method_row(res, "es_chisq", "g")
  g_phi   <- .method_row(res, "es_phi", "g")
  expect_true(is.finite(g_chisq), info = "chi-square -> g")
  expect_true(is.finite(g_phi),   info = "phi -> g")
  expect_true(g_chisq > 0 && g_phi > 0)
})

# ==============================================================================
# Orphan path 2: es_from_etasq (anova_list_L11)
# ==============================================================================
test_that("eta-squared input produces finite SMD", {
  dat <- data.frame(
    study_id = "eta1",
    etasq = 0.10,
    n_exp = 40,
    n_nexp = 40
  )
  res <- convert_df(dat, measure = "g", verbose = FALSE)
  g_eta <- .method_row(res, "es_etasq", "g")
  expect_true(is.finite(g_eta))
  expect_true(g_eta > 0)
})

# ==============================================================================
# Orphan path 3: es_from_pt_bis_r (point-biserial correlation)
# ==============================================================================
test_that("point-biserial correlation input produces finite r and g", {
  dat <- data.frame(
    study_id = "pb1",
    pt_bis_r = 0.35,
    n_exp = 30,
    n_nexp = 30
  )
  res <- convert_df(dat, measure = "r", verbose = FALSE)
  r_pb <- .method_row(res, "es_r_point_bis", "r")
  expect_true(is.finite(r_pb))
  # es_from_pt_bis_r() adjusts the point-biserial to a Pearson r by accounting
  # for the dichotomization of one variable (Lipsey & Wilson, 2001) — the
  # returned r is larger in magnitude than the raw pt-biserial. We only assert
  # sign and plausible magnitude, not equality.
  expect_true(r_pb > 0 && r_pb < 1)
  expect_true(r_pb >= 0.35)
})

# ==============================================================================
# Orphan path 4: es_from_ancova_means_sd (ANCOVA adjusted means + SDs)
# ==============================================================================
test_that("ANCOVA adjusted means + SDs produce finite SMD", {
  dat <- data.frame(
    study_id = "anc1",
    ancova_mean_exp = 10,
    ancova_mean_nexp = 12,
    ancova_mean_sd_exp = 3,
    ancova_mean_sd_nexp = 3,
    n_exp = 50,
    n_nexp = 50,
    cov_outcome_r = 0.5,
    n_cov_ancova = 1
  )
  res <- convert_df(dat, measure = "g", verbose = FALSE)
  g_ancova <- .method_row(res, "es_ancova_means_sd", "g")
  expect_true(is.finite(g_ancova))
  expect_true(g_ancova < 0)  # exp mean (10) < nexp mean (12)
})

# ==============================================================================
# Orphan path 5: es_from_ancova_md_sd (ANCOVA adjusted mean-difference)
# ==============================================================================
test_that("ANCOVA mean difference input produces finite SMD", {
  dat <- data.frame(
    study_id = "ancmd1",
    ancova_md = 2.0,
    ancova_md_sd = 4.0,
    n_exp = 50,
    n_nexp = 50,
    cov_outcome_r = 0.5,
    n_cov_ancova = 1
  )
  res <- convert_df(dat, measure = "g", verbose = FALSE)
  g_ancmd <- .method_row(res, "es_ancova_md_sd", "g")
  expect_true(is.finite(g_ancmd))
  expect_true(g_ancmd > 0)
})

# ==============================================================================
# Orphan path 6: es_from_beta_std (standardized regression beta)
# ==============================================================================
test_that("standardized beta input produces finite SMD", {
  dat <- data.frame(
    study_id = "beta1",
    beta_std = 0.25,
    sd_dv = 10,
    n_exp = 60,
    n_nexp = 60
  )
  res <- convert_df(dat, measure = "g", verbose = FALSE)
  g_beta <- .method_row(res, "es_std_beta", "g")
  expect_true(is.finite(g_beta))
  expect_true(g_beta > 0)
})

# ==============================================================================
# stand_RR delta chain: rr + logrr_se must yield RD, NNT with finite SEs
# ==============================================================================
test_that("es_from_rr_se delta chain produces finite RD and NNT with SE", {
  dat <- data.frame(
    study_id = "rr1",
    rr = 0.7,
    logrr_se = 0.15,
    baseline_risk = 0.3,
    n_exp = 200,
    n_nexp = 200,
    n_cases = 100,
    n_controls = 100
  )
  res <- convert_df(dat, measure = "rd", verbose = FALSE)
  rd      <- .method_row(res, "es_rr_se", "rd")
  rd_se   <- .method_row(res, "es_rr_se", "rd_se")
  nnt     <- .method_row(res, "es_rr_se", "nnt")
  nnt_se  <- .method_row(res, "es_rr_se", "nnt_se")
  expect_true(all(is.finite(c(rd, rd_se, nnt, nnt_se))))
  # RD = baseline_risk * (1 - rr) = 0.3 * 0.3 = 0.09
  expect_equal(rd, 0.09, tolerance = 1e-6)
  # NNT = 1 / 0.09 ~ 11.11
  expect_equal(nnt, 1 / 0.09, tolerance = 1e-6)
})

# ==============================================================================
# stand_IRR (es_from_cases_time) person-time NNT from cases + time
# ==============================================================================
test_that("es_from_cases_time produces finite person-time NNT via IRR", {
  dat <- data.frame(
    study_id = "irr1",
    n_cases_exp = 10,
    n_cases_nexp = 20,
    time_exp = 1000,
    time_nexp = 1000
  )
  # baseline_rate auto-computed from n_cases_nexp / time_nexp = 0.02
  res <- convert_df(dat, measure = "nnt", verbose = FALSE)
  nnt    <- .method_row(res, "es_cases_time", "nnt")
  nnt_se <- .method_row(res, "es_cases_time", "nnt_se")
  logirr <- .method_row(res, "es_cases_time", "logirr")
  expect_true(is.finite(logirr))
  expect_true(is.finite(nnt))
  expect_true(is.finite(nnt_se))
  # IRR = (10/1000) / (20/1000) = 0.5, so logIRR = log(0.5)
  expect_equal(logirr, log(0.5), tolerance = 1e-6)
  # IRD = baseline_rate * (1 - IRR) = 0.02 * 0.5 = 0.01, so NNT = 100
  expect_equal(nnt, 100, tolerance = 1e-4)
})

# ==============================================================================
# correct_inputs = FALSE: invalid values preserved AND surfaced in es_flags
# ==============================================================================
test_that("correct_inputs = FALSE preserves raw values and still flags them", {
  # A row with a negative SD — invalid, must be flagged either way.
  dat <- data.frame(
    study_id = "bad1",
    mean_exp = 10, mean_nexp = 12,
    mean_sd_exp = -3,              # invalid: SDs must be >= 0
    mean_sd_nexp = 3,
    n_exp = 40, n_nexp = 40
  )

  # correct_inputs = TRUE (default): value replaced with NA
  res_corr <- convert_df(dat, measure = "g", verbose = FALSE,
                         correct_inputs = TRUE)
  raw_corr <- attr(res_corr, "raw_data")
  expect_true(is.na(raw_corr$mean_sd_exp),
              info = "correct_inputs=TRUE must set negative SD to NA")

  # correct_inputs = FALSE: value preserved, flag still present
  res_keep <- convert_df(dat, measure = "g", verbose = FALSE,
                         correct_inputs = FALSE)
  raw_keep <- attr(res_keep, "raw_data")
  expect_equal(raw_keep$mean_sd_exp, -3,
               info = "correct_inputs=FALSE must preserve raw invalid value")

  # Tier-1 validation issues are stored as an attribute and must mention the
  # negative SD regardless of correct_inputs
  issues_corr <- attr(res_corr, "input_validation")
  issues_keep <- attr(res_keep, "input_validation")
  expect_true(length(issues_corr) >= 1 &&
              any(grepl("Negative input", unlist(issues_corr))))
  expect_true(length(issues_keep) >= 1 &&
              any(grepl("Negative input", unlist(issues_keep))))

  # After summary(), the flags column must carry the Tier-1 message.
  # With split_adjusted = TRUE (default) the column is flags_crude; otherwise
  # it's flags. Accept either shape so the test is robust to the default.
  s_keep <- summary(res_keep, flags = TRUE)
  flag_col <- intersect(c("flags_crude", "flags"), colnames(s_keep))[1]
  expect_false(is.na(flag_col),
               info = "summary() must expose a flags column when flags=TRUE")
  expect_true(any(grepl("Negative input", s_keep[[flag_col]])),
              info = "summary() must propagate Tier-1 flags to the flags column")
})
