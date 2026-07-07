# ==============================================================================
# Regression: for the ratio measures, the non-auto selection modes
# (hierarchy / minimum / maximum) use the SAME method ordering as
# es_selected = "auto".
#
# convert_df() previously defined the method ordering twice -- once per measure
# in the es_selected == "auto" block, and once as a single generic fallback for
# every other mode. The two drifted: e.g. for measure = "rr" the non-auto
# fallback ranked a reconstructed 2x2 above the directly reported risk ratio,
# whereas auto correctly preferred rr_se. The non-auto block now mirrors auto
# per measure. This test locks that parity in.
#
# NB: this parity applies to the FALLBACK ordering, i.e. when the hierarchy
# string does not list an applicable method. For the ratio measures the default
# hierarchy string ("means_sd > means_se > means_ci") never applies, so the
# fallback ordering governs. For d/g/r/z the default string DOES apply and
# hierarchy mode intentionally diverges from auto -- that divergence is the
# point of hierarchy mode and is NOT tested here.
# ==============================================================================

library(testthat)
library(metaConvert)

test_that("hierarchy prefers the reported RR over a reconstructed 2x2 (mirrors auto)", {
  # A row carrying BOTH a directly reported risk ratio AND a 2x2 table.
  dat <- data.frame(
    rr = 1.5, logrr_se = 0.2, baseline_risk = 0.2,
    n_cases_exp = 30, n_controls_exp = 70,
    n_cases_nexp = 20, n_controls_nexp = 80,
    n_exp = 100, n_nexp = 100
  )
  a <- summary(convert_df(dat, measure = "rr", es_selected = "auto",
                          verbose = FALSE), digits = 11)
  h <- summary(convert_df(dat, measure = "rr", es_selected = "hierarchy",
                          verbose = FALSE), digits = 11)

  expect_equal(a$info_used_crude, "rr_se")          # reported RR, not the 2x2
  expect_equal(h$info_used_crude, a$info_used_crude) # hierarchy fallback == auto
})

test_that("hierarchy prefers the reported OR over a reconstructed 2x2 (mirrors auto)", {
  dat <- data.frame(
    or = 2.0, logor_se = 0.25,
    n_cases_exp = 30, n_controls_exp = 70,
    n_cases_nexp = 20, n_controls_nexp = 80,
    n_exp = 100, n_nexp = 100
  )
  a <- summary(convert_df(dat, measure = "or", es_selected = "auto",
                          verbose = FALSE), digits = 11)
  h <- summary(convert_df(dat, measure = "or", es_selected = "hierarchy",
                          verbose = FALSE), digits = 11)

  # A reported-OR method wins (info_used starts with "or"), not the 2x2.
  expect_true(startsWith(a$info_used_crude, "or"))
  expect_false(identical(a$info_used_crude, "2x2"))
  expect_equal(h$info_used_crude, a$info_used_crude) # hierarchy fallback == auto
})

test_that("ratio measures select the same method in auto and hierarchy (fallback parity)", {
  # No continuous data, so the default hierarchy string never applies and the
  # non-auto fallback ordering (now mirroring auto) governs every ratio measure.
  dat <- data.frame(
    rr = 1.5, logrr_se = 0.2,
    or = 2.0, logor_se = 0.25,
    rd = 0.1, rd_se = 0.04,
    baseline_risk = 0.2,
    n_cases_exp = 30, n_controls_exp = 70,
    n_cases_nexp = 20, n_controls_nexp = 80,
    n_exp = 100, n_nexp = 100
  )
  for (m in c("or", "logor", "rr", "logrr", "nnt", "rd")) {
    a <- summary(convert_df(dat, measure = m, es_selected = "auto",
                            verbose = FALSE))
    h <- summary(convert_df(dat, measure = m, es_selected = "hierarchy",
                            verbose = FALSE))
    expect_equal(h$info_used_crude, a$info_used_crude,
                 info = sprintf("measure = '%s'", m))
  }
})
