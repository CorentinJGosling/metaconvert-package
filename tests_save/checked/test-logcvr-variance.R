# Regression test for the lnCVR sampling-variance fix (A3).
#
# The Var(ln s) term in the lnCVR variance (Senior et al. 2020 eq. 16) carries a 1/2
# factor (Var(ln s) ~= 0.5*n/(n-1)^2 for a normal sample). It had been written without
# the 1/2, making the lnCVR variance ~2x too large (SE ~40% too wide). The sibling
# logVR SE (same Var(ln s) quantity) was already correct, so this pins both.

test_that("logCVR SE matches metafor escalc CVR after the 1/2-factor fix", {
  skip_if_not_installed("metafor")
  m1 <- 2.3; s1 <- 1.2; n1 <- 55; m2 <- 1.9; s2 <- 0.9; n2 <- 55
  es <- es_variab_from_means_sd(mean_exp = m1, mean_sd_exp = s1, n_exp = n1,
                                mean_nexp = m2, mean_sd_nexp = s2, n_nexp = n2)
  mf <- metafor::escalc(measure = "CVR", m1i = m1, sd1i = s1, n1i = n1,
                        m2i = m2, sd2i = s2, n2i = n2)
  se_mf <- sqrt(as.numeric(mf$vi[1]))
  # metaConvert's IND2 form adds a small 4th-moment term metafor omits, so allow ~3%.
  expect_equal(es$logcvr_se[1], se_mf, tolerance = 0.03)
  expect_equal(es$logcvr[1], as.numeric(mf$yi[1]), tolerance = 0.02)
  # guard against a regression back to the ~1.3x-too-large (no-1/2) value
  expect_lt(es$logcvr_se[1], 1.15 * se_mf)
})

test_that("logVR SE is unaffected by the logCVR fix", {
  s1 <- 1.2; n1 <- 55; s2 <- 0.9; n2 <- 55
  es <- es_variab_from_means_sd(mean_exp = 2.3, mean_sd_exp = s1, n_exp = n1,
                                mean_nexp = 1.9, mean_sd_nexp = s2, n_nexp = n2)
  # metaConvert uses the n/(n-1)^2 small-sample form (with its 1/2 factor); this was
  # already correct before the logCVR fix and must stay exact.
  expect_equal(es$logvr_se[1],
               sqrt(0.5 * (n2 / ((n2 - 1)^2) + n1 / ((n1 - 1)^2))),
               tolerance = 1e-9)
  # ~within 1% of metafor's 1/(n-1) form (a documented small-sample convention diff)
  if (requireNamespace("metafor", quietly = TRUE)) {
    mf <- metafor::escalc(measure = "VR", sd1i = s1, n1i = n1, sd2i = s2, n2i = n2)
    expect_equal(es$logvr_se[1], sqrt(as.numeric(mf$vi[1])), tolerance = 0.02)
  }
})
