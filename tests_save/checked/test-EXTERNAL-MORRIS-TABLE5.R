# =============================================================================
# PUBLISHED external anchor: Morris (2008), Table 5.
#
# Morris SB (2008), Organizational Research Methods 11(2):364-386, analyses five
# studies from Carlson & Schmidt (1999) and prints d_ppc1 and d_ppc2 in Table 5.
# Viechtbauer reproduces both on the metafor-project site
# (https://www.metafor-project.org/doku.php/analyses:morris2008), which also
# pins down what "the metafor way" actually is for this design:
#
#   d_ppc1 (Becker 1988): compute a standardized mean change WITHIN each arm
#     (escalc measure = "SMCR"), then SUBTRACT the two and ADD their sampling
#     variances (the arms are independent). This is metaConvert's pool_sd = FALSE.
#
#   d_ppc2 (Morris 2008 eq. 8-9): divide the difference in mean change by a
#     single pretest SD pooled across arms. This is metaConvert's pool_sd = TRUE.
#
# Viechtbauer notes d_ppc2 "assumes that the true pretest SDs are equal for the
# two groups", while d_ppc1 "does not make that assumption and therefore is more
# broadly applicable (but may be slightly less efficient)". Morris recommends
# d_ppc2. Both are supported; pool_sd = FALSE is the default.
#
# This file pins the POINT ESTIMATES of both columns against the published values.
# (Sampling variances are deliberately NOT pinned to Viechtbauer's numbers: he
# uses the homoscedastic SMCR/eq.-25-style forms, whereas metaConvert uses the
# heteroscedasticity-robust SMCRH-style forms. The two agree exactly when
# SD_pre = SD_post and diverge otherwise -- which is asserted below.)
# =============================================================================

datT <- data.frame(
  m_pre   = c(30.6, 23.5, 0.5, 53.4, 35.6),
  m_post  = c(38.5, 26.8, 0.7, 75.9, 36.0),
  sd_pre  = c(15.0, 3.1, 0.1, 14.5, 4.7),
  sd_post = c(11.6, 4.1, 0.1, 4.4, 4.6),
  ni      = c(20, 50, 9, 10, 14),
  ri      = c(0.47, 0.64, 0.77, 0.89, 0.44)
)
datC <- data.frame(
  m_pre   = c(23.1, 24.9, 0.6, 55.7, 34.8),
  m_post  = c(19.7, 25.3, 0.6, 60.7, 33.4),
  sd_pre  = c(13.8, 4.1, 0.2, 17.3, 3.1),
  sd_post = c(14.8, 3.3, 0.2, 17.9, 6.9),
  ni      = c(20, 42, 9, 11, 14),
  ri      = c(0.47, 0.64, 0.77, 0.89, 0.44)
)

mc_fit <- function(pool_sd) {
  es_from_means_sd_pre_post(
    n_exp = datT$ni, n_nexp = datC$ni,
    mean_pre_exp = datT$m_pre, mean_exp = datT$m_post,
    mean_pre_sd_exp = datT$sd_pre, mean_sd_exp = datT$sd_post,
    mean_pre_nexp = datC$m_pre, mean_nexp = datC$m_post,
    mean_pre_sd_nexp = datC$sd_pre, mean_sd_nexp = datC$sd_post,
    r_pre_post_exp = datT$ri, r_pre_post_nexp = datC$ri,
    pre_post_to_smd = "bonett", pool_sd = pool_sd
  )
}

test_that("pool_sd = FALSE reproduces Morris (2008) Table 5, column d_ppc1", {
  skip_if_not_installed("metafor")
  # Viechtbauer's own computation, verbatim from the metafor-project page.
  eT <- metafor::escalc(measure = "SMCR", m1i = m_post, m2i = m_pre,
                        sd1i = sd_pre, ni = ni, ri = ri, data = datT)
  eC <- metafor::escalc(measure = "SMCR", m1i = m_post, m2i = m_pre,
                        sd1i = sd_pre, ni = ni, ri = ri, data = datC)
  d_ppc1 <- as.numeric(eT$yi) - as.numeric(eC$yi)

  expect_equal(mc_fit(FALSE)$g, d_ppc1, tolerance = 1e-10)
  # and against the values printed on the page / in Morris Table 5
  expect_equal(round(mc_fit(FALSE)$g, 2), c(0.74, 0.95, 1.81, 1.15, 0.51),
               tolerance = 1e-8)
})

test_that("pool_sd = TRUE reproduces Morris (2008) Table 5, column d_ppc2", {
  skip_if_not_installed("metafor")
  # Viechtbauer's own computation, verbatim from the metafor-project page.
  sd_pool <- sqrt((with(datT, (ni - 1) * sd_pre^2) +
                   with(datC, (ni - 1) * sd_pre^2)) / (datT$ni + datC$ni - 2))
  d_ppc2 <- metafor:::.cmicalc(datT$ni + datC$ni - 2) *
    (with(datT, m_post - m_pre) - with(datC, m_post - m_pre)) / sd_pool

  expect_equal(mc_fit(TRUE)$g, as.numeric(d_ppc2), tolerance = 1e-10)
  expect_equal(round(mc_fit(TRUE)$g, 2), c(0.77, 0.80, 1.20, 1.05, 0.44),
               tolerance = 1e-8)
})

test_that("the pooled variance reduces EXACTLY to Viechtbauer's when SD_pre = SD_post", {
  skip_if_not_installed("metafor")
  # Viechtbauer's d_ppc2 sampling variance (metafor-project page):
  #   vi = 2*(1 - r)*(1/nT + 1/nC) + yi^2/(2*N)
  # It assumes SD_pre = SD_post WITHIN each arm. metaConvert instead takes the
  # numerator variance from the empirical pooled change SD, so the two coincide
  # exactly under that assumption and diverge when it fails.
  res <- mc_fit(TRUE)
  vi_v <- 2 * (1 - datT$ri) * (1 / datT$ni + 1 / datC$ni) +
    res$g^2 / (2 * (datT$ni + datC$ni))

  # Study 3 is the one study with sd_pre == sd_post in BOTH arms (0.1/0.1, 0.2/0.2):
  # Viechtbauer's assumption holds exactly, so the variances must agree exactly.
  expect_equal(res$g_se[3]^2, vi_v[3], tolerance = 1e-8)

  # Study 5's control arm goes 3.1 -> 6.9: the assumption fails badly, and
  # metaConvert must report a LARGER variance (the homoscedastic form understates it).
  expect_gt(res$g_se[5]^2, vi_v[5] * 1.5)
})

test_that("pool_sd = FALSE is the default", {
  expect_equal(mc_fit(FALSE)$g, mc_fit(pool_sd = formals(es_from_means_sd_pre_post)$pool_sd)$g,
               tolerance = 1e-12)
  expect_false(isTRUE(formals(es_from_means_sd_pre_post)$pool_sd))
  expect_false(isTRUE(formals(convert_df)$pool_sd))
})
