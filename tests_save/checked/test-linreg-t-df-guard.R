# Regression test for es_from_linreg_t non-positive-df handling (A7).
#
# A row with n_sample - n_covariates - 2 <= 0 makes the partial correlation undefined.
# es_from_linreg_t used stop(), which aborted the ENTIRE convert_df() run (it is called
# on the whole column) even when the offending row was a non-regression row that merely
# carried a stray n_covariates. Now: NA the bad rows + warn, valid rows continue.

test_that("a stray n_covariates on a means+SD row no longer aborts convert_df", {
  df <- data.frame(mean_exp = 5, mean_nexp = 4, mean_sd_exp = 1, mean_sd_nexp = 1,
                   n_exp = 20, n_nexp = 20, n_sample = 40, n_covariates = 39)  # df = -1
  expect_no_error(
    s <- suppressWarnings(suppressMessages(summary(convert_df(df, measure = "d"))))
  )
  es_col <- grep("^es_crude$|^es$", colnames(s), value = TRUE)[1]
  expect_false(is.na(s[[es_col]][1]))  # the means_sd Cohen's d is still computed
})

test_that("one bad regression row does not abort the run", {
  df <- data.frame(linreg_t = c(2.5, 3.0), n_sample = c(100, 40), n_covariates = c(2, 39))
  expect_no_error(
    suppressWarnings(suppressMessages(convert_df(df, measure = "r")))
  )
})

test_that("direct call warns on bad rows, NAs them, and keeps valid rows exact", {
  expect_warning(
    e <- es_from_linreg_t(linreg_t = c(6.19, 2.0), n_sample = c(232, 5),
                          n_covariates = c(6, 6)),
    "non-positive residual df"
  )
  # partial correlation rp = t / sqrt(t^2 + df); row 1 df = 232 - 6 - 2 = 224
  expect_equal(e$rp[1], 6.19 / sqrt(6.19^2 + 224), tolerance = 1e-4)
  expect_false(is.na(e$rp[1]))  # valid row computed
  expect_true(is.na(e$rp[2]))   # bad row (df = 5 - 6 - 2 = -3) NA, not an abort
})
