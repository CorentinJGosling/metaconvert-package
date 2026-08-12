# Regression test for the reverse-flag CI on ratio conversions (A6).
#
# On reverse the point estimate is negated, so a symmetric log-scale CI must be NEGATED
# AND SWAPPED: new_lo = -old_up, new_up = -old_lo. The converted-CI blocks had swapped
# WITHOUT negating, producing an inverted (lo > up), wrong-signed interval that did not
# bracket the reversed estimate. Sites: es_from_stand_RR.R (RR->OR logor CI) and
# es_from_stand_OR.R (OR->RR logrr CI, OR->r CI, OR->z CI).

test_that("RR -> OR: reversed logor CI is ordered, brackets the estimate, and is negate+swap of forward", {
  f <- es_from_rr_se(rr = 2.0, logrr_se = 0.2, n_cases = 80, n_controls = 120,
                     n_exp = 100, n_nexp = 100, rr_to_or = "metaumbrella")
  r <- es_from_rr_se(rr = 2.0, logrr_se = 0.2, n_cases = 80, n_controls = 120,
                     n_exp = 100, n_nexp = 100, rr_to_or = "metaumbrella", reverse_rr = TRUE)
  expect_lte(r$logor_ci_lo[1], r$logor_ci_up[1])
  expect_gte(r$logor[1], r$logor_ci_lo[1]); expect_lte(r$logor[1], r$logor_ci_up[1])
  expect_equal(r$logor[1], -f$logor[1], tolerance = 1e-9)
  expect_equal(r$logor_ci_lo[1], -f$logor_ci_up[1], tolerance = 1e-9)
  expect_equal(r$logor_ci_up[1], -f$logor_ci_lo[1], tolerance = 1e-9)
})

test_that("OR -> RR/r/z: reversed converted CIs are ordered, bracket the estimate, negate+swap", {
  f <- es_from_or_se(or = 2.5, logor_se = 0.25, n_cases = 90, n_controls = 110,
                     n_exp = 100, n_nexp = 100, n_sample = 200,
                     or_to_rr = "metaumbrella_cases", or_to_cor = "pearson")
  r <- es_from_or_se(or = 2.5, logor_se = 0.25, n_cases = 90, n_controls = 110,
                     n_exp = 100, n_nexp = 100, n_sample = 200,
                     or_to_rr = "metaumbrella_cases", or_to_cor = "pearson", reverse_or = TRUE)
  # logrr
  expect_lte(r$logrr_ci_lo[1], r$logrr_ci_up[1])
  expect_equal(r$logrr_ci_lo[1], -f$logrr_ci_up[1], tolerance = 1e-9)
  expect_equal(r$logrr_ci_up[1], -f$logrr_ci_lo[1], tolerance = 1e-9)
  # r (tanh is odd, so negate+swap is correct even for the asymmetric r interval)
  expect_lte(r$r_ci_lo[1], r$r_ci_up[1])
  expect_gte(r$r[1], r$r_ci_lo[1]); expect_lte(r$r[1], r$r_ci_up[1])
  expect_equal(r$r_ci_lo[1], -f$r_ci_up[1], tolerance = 1e-9)
  expect_equal(r$r_ci_up[1], -f$r_ci_lo[1], tolerance = 1e-9)
  # z
  expect_lte(r$z_ci_lo[1], r$z_ci_up[1])
  expect_gte(r$z[1], r$z_ci_lo[1]); expect_lte(r$z[1], r$z_ci_up[1])
  expect_equal(r$z_ci_lo[1], -f$z_ci_up[1], tolerance = 1e-9)
  expect_equal(r$z_ci_up[1], -f$z_ci_lo[1], tolerance = 1e-9)
})

test_that("non-reversed converted CIs remain ordered and unchanged", {
  f <- es_from_rr_se(rr = 2.0, logrr_se = 0.2, n_cases = 80, n_controls = 120,
                     n_exp = 100, n_nexp = 100, rr_to_or = "metaumbrella")
  expect_lte(f$logor_ci_lo[1], f$logor_ci_up[1])
  expect_gte(f$logor[1], f$logor_ci_lo[1]); expect_lte(f$logor[1], f$logor_ci_up[1])
})

test_that("convert_df(measure='or') summary OR CI is ordered under reverse_rr (end-to-end)", {
  df <- data.frame(rr = 2.0, logrr_se = 0.2, n_cases = 80, n_controls = 120,
                   n_exp = 100, n_nexp = 100, reverse_rr = TRUE)
  s <- suppressWarnings(suppressMessages(
    summary(convert_df(df, measure = "or", rr_to_or = "metaumbrella"))
  ))
  lo <- grep("ci_lo", colnames(s), value = TRUE)[1]
  up <- grep("ci_up", colnames(s), value = TRUE)[1]
  expect_lte(s[[lo]][1], s[[up]][1])   # natural-scale OR CI ordered (was inverted)
})
