# Regression test for the D2 SE-outlier normalization on arm-size-less d/g pools (A8).
#
# The d/g balanced-arms 2/sqrt(N) fallback was nested inside if(any(has_both)), so a
# Cohen's d pool built entirely from pearson_r/fisher_z inputs (which carry only
# n_sample; n_exp/n_nexp stay NA) skipped normalization and entered RAW, N-dependent SEs
# into the IQR -- falsely flagging a legitimately small-N study as an SE outlier. Fix:
# run the fallback at branch level (like the logOR/logRR branch).

test_that("small-N study in a pure pearson_r d pool is NOT flagged as an SE outlier", {
  df <- data.frame(pearson_r = rep(0.30, 6), n_sample = c(300, 350, 400, 450, 500, 25))
  s <- suppressWarnings(suppressMessages(
    summary(convert_df(df, measure = "d"), flags = TRUE)
  ))
  fl_col <- grep("flag", colnames(s), value = TRUE)[1]
  # the N=25 row used to be flagged; after the fix no row is an SE outlier
  expect_false(any(grepl("SE outlier", s[[fl_col]])))
})

test_that("a genuinely pathological SE in an arm-size-less d pool is still flagged", {
  opts <- metaConvert:::.default_flag_options()
  n_total <- rep(400, 6)                       # uniform N; 2/sqrt(N) normalization uniform
  norm <- c(0.85, 0.92, 1.0, 1.08, 1.15, 10)   # target normalized SEs: realistic spread + 1 outlier
  se <- norm * (2 / sqrt(n_total))             # so se_for_iqr == norm after the 2/sqrt(N) fix
  fl <- metaConvert:::.flag_cross_row_outliers(
    es = rep(0.3, 6), se = se, opts = opts, measure = "d",
    n_total = n_total, n_exp = rep(NA_real_, 6), n_nexp = rep(NA_real_, 6)
  )
  expect_true(any(grepl("SE outlier", fl[[6]])))            # pathological row still fires
  expect_false(any(grepl("SE outlier", unlist(fl[1:5]))))   # normal rows do not
})
