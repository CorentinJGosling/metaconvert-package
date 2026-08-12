# Fixes from the 2026-07-16 formula audit (AUDIT-FIX-TRACKER.md):
# the package must not flag its own correct output (P6/P7/P8), must keep
# flagging genuinely wrong CIs, and several consistency fixes (P17, P22-P25,
# V31). Each test failed on the pre-fix tree (TDD protocol).

fc_flag_col <- function(s) {
  col <- grep("^flags", names(s), value = TRUE)[1]
  as.character(s[[col]])
}

# ---------------------------------------------------------------------------
# P6 - A6 must accept the package's own Welch-df md CI (band: pooled..min-1 df)
# ---------------------------------------------------------------------------
test_that("P6: Welch md CI on clean unequal-variance data is not flagged", {
  df <- data.frame(study_id = "s1", n_exp = 5, n_nexp = 45,
                   mean_exp = 10, mean_sd_exp = 12, mean_nexp = 4, mean_sd_nexp = 1.5)
  s <- suppressMessages(summary(convert_df(df, measure = "md", verbose = FALSE), flags = TRUE))
  expect_false(grepl("CI width inconsistent", fc_flag_col(s)[1]))

  # moderate imbalance case from the audit
  df2 <- data.frame(study_id = "s1", n_exp = 15, n_nexp = 30,
                    mean_exp = 8, mean_sd_exp = 6, mean_nexp = 6, mean_sd_nexp = 1.5)
  s2 <- suppressMessages(summary(convert_df(df2, measure = "md", verbose = FALSE), flags = TRUE))
  expect_false(grepl("CI width inconsistent", fc_flag_col(s2)[1]))
})

test_that("P6: a genuinely wrong md CI still fires A6", {
  # user-entered md whose CI is twice the width its SE implies
  df <- data.frame(study_id = "s1", n_exp = 20, n_nexp = 20,
                   user_es_original_measure_crude = "md",
                   user_es_crude = 2, user_se_crude = 1,
                   user_ci_lo_crude = 2 - 2 * qnorm(.975) * 1,
                   user_ci_up_crude = 2 + 2 * qnorm(.975) * 1)
  s <- suppressMessages(summary(convert_df(df, measure = "md", verbose = FALSE), flags = TRUE))
  expect_true(grepl("CI width inconsistent", fc_flag_col(s)[1]))
})

# ---------------------------------------------------------------------------
# P7 - A6 must expect qt(.975, n_nexp - 1) for Glass rows
# ---------------------------------------------------------------------------
test_that("P7: clean Glass rows are not flagged; pooled rows unchanged", {
  df <- data.frame(study_id = "s1", n_exp = 10, n_nexp = 10,
                   mean_exp = 2.3, mean_sd_exp = 1.2, mean_nexp = 1.9, mean_sd_nexp = 0.9)
  s_glass <- suppressMessages(summary(
    convert_df(df, measure = "g", smd_denom = "glass", verbose = FALSE), flags = TRUE))
  expect_false(grepl("CI width inconsistent", fc_flag_col(s_glass)[1]))

  # unbalanced case from the audit (n = 50/8: qt(7) vs qt(56) differ by 13%)
  df2 <- data.frame(study_id = "s1", n_exp = 50, n_nexp = 8,
                    mean_exp = 2.3, mean_sd_exp = 1.2, mean_nexp = 1.9, mean_sd_nexp = 0.9)
  s_glass2 <- suppressMessages(summary(
    convert_df(df2, measure = "g", smd_denom = "glass", verbose = FALSE), flags = TRUE))
  expect_false(grepl("CI width inconsistent", fc_flag_col(s_glass2)[1]))

  # pooled rows keep the pooled expectation (clean row, still silent)
  s_pooled <- suppressMessages(summary(
    convert_df(df, measure = "g", verbose = FALSE), flags = TRUE))
  expect_false(grepl("CI width inconsistent", fc_flag_col(s_pooled)[1]))
})

# ---------------------------------------------------------------------------
# P8 - A6 must accept the deliberate [0,1] clamp on raw proportion CIs
# ---------------------------------------------------------------------------
test_that("P8: clamped raw-prop CIs are not flagged; interior rows unchanged", {
  df <- data.frame(study_id = 1:2, prop = c(0.05, 0.30), n_sample = c(20, 100))
  s <- suppressMessages(summary(convert_df(df, measure = "prop", verbose = FALSE), flags = TRUE))
  fl <- fc_flag_col(s)
  expect_false(grepl("CI width inconsistent", fl[1]))  # clamped at 0
  expect_false(grepl("CI width inconsistent", fl[2]))  # interior

  # boundary row via counts (continuity corrected, CI clamped at 0)
  df2 <- data.frame(study_id = 1, n_cases = 0, n_sample = 20)
  s2 <- suppressMessages(summary(convert_df(df2, measure = "prop", verbose = FALSE), flags = TRUE))
  expect_false(grepl("CI width inconsistent", fc_flag_col(s2)[1]))
})

# ---------------------------------------------------------------------------
# P17 - raw-scale Wald CI escaping [0/-1, 1] gets an [UNUSUAL] bound flag
# ---------------------------------------------------------------------------
test_that("P17: raw alpha/ICC CI > 1 flagged [UNUSUAL]; bonett scale untouched", {
  df <- data.frame(study_id = 1, cronbach_alpha = 0.9, n_sample = 10, n_items = 3)
  s_raw <- suppressMessages(summary(
    convert_df(df, measure = "alpha", alpha_to_es = "raw", verbose = FALSE), flags = TRUE))
  expect_true(grepl("CI upper bound exceeds 1", fc_flag_col(s_raw)[1]))

  s_bon <- suppressMessages(summary(
    convert_df(df, measure = "alpha", verbose = FALSE), flags = TRUE))
  expect_false(grepl("CI upper bound exceeds 1", fc_flag_col(s_bon)[1]))

  df_icc <- data.frame(study_id = 1, icc = 0.9, n_sample = 10, n_measurements = 2,
                       icc_type = "consistency")
  s_icc <- suppressMessages(summary(
    convert_df(df_icc, measure = "icc", icc_to_es = "raw", verbose = FALSE), flags = TRUE))
  expect_true(grepl("ICC CI bound outside", fc_flag_col(s_icc)[1]))
})

# ---------------------------------------------------------------------------
# V31 - agreement-type ICC rows carry the anti-conservative-SE info note
# ---------------------------------------------------------------------------
test_that("V31: agreement rows get the [INFO] note, consistency rows do not", {
  df <- data.frame(study_id = 1:2, icc = c(0.8, 0.8), n_sample = 50,
                   n_measurements = 2, icc_type = c("agreement", "consistency"))
  s <- suppressMessages(summary(convert_df(df, measure = "icc", verbose = FALSE), flags = TRUE))
  fl <- fc_flag_col(s)
  expect_true(grepl("between-rater variance", fl[1]))
  expect_false(grepl("between-rater variance", fl[2]))
})

# ---------------------------------------------------------------------------
# P25 - E7 fires when a pool contains ONLY paired-t rows under pool_sd = TRUE
# ---------------------------------------------------------------------------
test_that("P25: all-paired-t pool under pool_sd = TRUE is disclosed", {
  df <- data.frame(study_id = c("a", "b", "c"),
                   paired_t_exp = c(2.5, 3.1, 2.2), paired_t_nexp = c(1.1, 0.9, 1.4),
                   n_exp = 20, n_nexp = 20, r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6)
  s_pool <- suppressMessages(summary(
    convert_df(df, measure = "g", pool_sd = TRUE, verbose = FALSE), flags = TRUE))
  expect_true(all(grepl("Per-arm standardizer", fc_flag_col(s_pool))))

  s_def <- suppressMessages(summary(
    convert_df(df, measure = "g", verbose = FALSE), flags = TRUE))
  expect_false(any(grepl("Per-arm standardizer", fc_flag_col(s_def))))
})

# ---------------------------------------------------------------------------
# P22 - deprecated 'measure' alias warns and does not override new-style arg
# ---------------------------------------------------------------------------
test_that("P22: alias warns; explicit user_es_target_measure_crude wins", {
  expect_warning(
    old_style <- es_from_user_crude(measure = "g", user_es_measure_crude = "d",
                                    user_es_crude = 0.5, user_se_crude = 0.1,
                                    n_exp = 20, n_nexp = 20),
    "deprecated"
  )
  new_style <- es_from_user_crude(user_es_target_measure_crude = "g",
                                  user_es_original_measure_crude = "d",
                                  user_es_crude = 0.5, user_se_crude = 0.1,
                                  n_exp = 20, n_nexp = 20)
  expect_equal(old_style$g, new_style$g)
  # both supplied: the new-style argument wins (alias ignored, still warns? no -
  # the guard skips the alias silently when the new-style arg is present)
  both <- es_from_user_crude(measure = "d", user_es_target_measure_crude = "g",
                             user_es_original_measure_crude = "d",
                             user_es_crude = 0.5, user_se_crude = 0.1,
                             n_exp = 20, n_nexp = 20)
  expect_equal(both$g, new_style$g)
})

# ---------------------------------------------------------------------------
# P23 - es_from_pt_bis_r honours smd_to_cor
# ---------------------------------------------------------------------------
test_that("P23: smd_to_cor changes the r returned by es_from_pt_bis_r", {
  a <- es_from_pt_bis_r(pt_bis_r = 0.35, n_exp = 13, n_nexp = 19,
                        smd_to_cor = "viechtbauer")
  b <- es_from_pt_bis_r(pt_bis_r = 0.35, n_exp = 13, n_nexp = 19,
                        smd_to_cor = "lipsey_cooper")
  expect_false(isTRUE(all.equal(a$r, b$r)))
})

# ---------------------------------------------------------------------------
# P24 - per-row smd_var NA falls back to the default instead of erroring
# ---------------------------------------------------------------------------
test_that("P24: NA smd_var defaults to borenstein (mirrors smd_denom)", {
  res <- es_from_cohen_d(cohen_d = c(0.5, 0.6), n_exp = c(10, 10), n_nexp = c(12, 12),
                         smd_var = c("borenstein", NA))
  ref <- es_from_cohen_d(cohen_d = c(0.5, 0.6), n_exp = c(10, 10), n_nexp = c(12, 12))
  expect_equal(res$g_se, ref$g_se)
})
