## Agreement with metafor::escalc().
##
## metaConvert and metafor compute an overlapping set of effect sizes from the same
## inputs, and where they overlap the numbers must match to machine precision -- users
## routinely feed a metaConvert summary() straight into rma(). This file pins that
## agreement measure by measure so a change to a conversion formula (or to a sampling
## variance) cannot silently break it. metafor is an Imports, so it is always available.
##
## Where the two packages legitimately differ, the difference is asserted explicitly
## rather than omitted -- see the mean-difference block.

tol <- 1e-10

test_that("continuous two-group measures match metafor", {
  m1 <- 10; m2 <- 8; s1 <- 3; s2 <- 3.2; n1 <- 40; n2 <- 60
  mc <- es_from_means_sd(mean_exp = m1, mean_nexp = m2, mean_sd_exp = s1,
                         mean_sd_nexp = s2, n_exp = n1, n_nexp = n2)

  mf <- metafor::escalc("SMD", m1i = m1, m2i = m2, sd1i = s1, sd2i = s2, n1i = n1, n2i = n2)
  expect_equal(mc$g, mf$yi[1], tolerance = tol, ignore_attr = TRUE)
  # metafor's SMD variance is the Hedges-Olkin one; metaConvert reaches it via smd_var
  expect_equal(
    es_from_means_sd(mean_exp = m1, mean_nexp = m2, mean_sd_exp = s1, mean_sd_nexp = s2,
                     n_exp = n1, n_nexp = n2, smd_var = "hedges_olkin")$g_se,
    sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE)

  mf <- metafor::escalc("RBIS", m1i = m1, m2i = m2, sd1i = s1, sd2i = s2, n1i = n1, n2i = n2)
  expect_equal(mc$r,    mf$yi[1],       tolerance = tol, ignore_attr = TRUE)
  expect_equal(mc$r_se, sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE)

  mf <- metafor::escalc("SMD1", m1i = m1, m2i = m2, sd2i = s2, n1i = n1, n2i = n2)
  expect_equal(
    es_from_means_sd(mean_exp = m1, mean_nexp = m2, mean_sd_exp = s1, mean_sd_nexp = s2,
                     n_exp = n1, n_nexp = n2, smd_denom = "control")$g,
    mf$yi[1], tolerance = tol, ignore_attr = TRUE)

  expect_equal(es_from_cohen_d(cohen_d = 0.6, n_exp = n1, n_nexp = n2)$logor,
               metafor::escalc("D2ORL", di = 0.6, n1i = n1, n2i = n2)$yi[1],
               tolerance = tol, ignore_attr = TRUE)
})

test_that("the mean difference matches the metafor variance appropriate to each route", {
  m1 <- 10; m2 <- 8; s1 <- 3; s2 <- 3.2; n1 <- 40; n2 <- 60
  mc <- es_from_means_sd(mean_exp = m1, mean_nexp = m2, mean_sd_exp = s1,
                         mean_sd_nexp = s2, n_exp = n1, n_nexp = n2)
  # per-arm SDs -> Welch (unpooled), i.e. metafor's DEFAULT vtype = "LS"
  mf <- metafor::escalc("MD", m1i = m1, m2i = m2, sd1i = s1, sd2i = s2, n1i = n1, n2i = n2)
  expect_equal(mc$md,    mf$yi[1],       tolerance = tol, ignore_attr = TRUE)
  expect_equal(mc$md_se, sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE)
  # a single pooled SD -> the homoscedastic variance, metafor's vtype = "HO"
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  mfho <- metafor::escalc("MD", m1i = m1, m2i = m2, sd1i = s1, sd2i = s2,
                          n1i = n1, n2i = n2, vtype = "HO")
  expect_equal(es_from_md_sd(md = m1 - m2, md_sd = sp, n_exp = n1, n_nexp = n2)$md_se,
               sqrt(mfho$vi[1]), tolerance = tol, ignore_attr = TRUE)
})

test_that("binary measures from a 2x2 table match metafor", {
  a <- 12; b <- 28; c_ <- 10; d_ <- 50
  mc <- es_from_2x2(n_cases_exp = a, n_controls_exp = b,
                    n_cases_nexp = c_, n_controls_nexp = d_)
  for (nm in c("OR", "RR")) {
    mf <- metafor::escalc(nm, ai = a, bi = b, ci = c_, di = d_)
    lo <- if (nm == "OR") "logor" else "logrr"
    expect_equal(mc[[lo]],               mf$yi[1],       tolerance = tol, ignore_attr = TRUE)
    expect_equal(mc[[paste0(lo, "_se")]], sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE)
  }
  mf <- metafor::escalc("RD", ai = a, bi = b, ci = c_, di = d_)
  # metaConvert reports the risk difference control-minus-exposed (so that a
  # protective effect gives a positive NNT); metafor reports exposed-minus-control
  expect_equal(-mc$rd,   mf$yi[1],       tolerance = tol, ignore_attr = TRUE)
  expect_equal(mc$rd_se, sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE)

  skip_if_not_installed("mvtnorm")   # metafor Suggests it; RTET needs it
  mf <- metafor::escalc("RTET", ai = a, bi = b, ci = c_, di = d_)
  expect_equal(mc$r,    mf$yi[1],       tolerance = tol, ignore_attr = TRUE)
  expect_equal(mc$r_se, sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE)
})

test_that("correlations and incidence rates match metafor", {
  r0 <- 0.42; ns <- 100
  mc <- es_from_pearson_r(pearson_r = r0, n_sample = ns)
  mf <- metafor::escalc("ZCOR", ri = r0, ni = ns)
  expect_equal(mc$z,    mf$yi[1],       tolerance = tol, ignore_attr = TRUE)
  expect_equal(mc$z_se, sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE)
  expect_equal(mc$r_se, sqrt(metafor::escalc("COR", ri = r0, ni = ns)$vi[1]),
               tolerance = tol, ignore_attr = TRUE)

  mci <- es_from_cases_time(n_cases_exp = 15, time_exp = 100,
                            n_cases_nexp = 25, time_nexp = 120)
  mf <- metafor::escalc("IRR", x1i = 15, t1i = 100, x2i = 25, t2i = 120)
  expect_equal(mci$logirr,    mf$yi[1],       tolerance = tol, ignore_attr = TRUE)
  expect_equal(mci$logirr_se, sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE)
})

test_that("single-group proportions match metafor, including the zero-cell correction", {
  k <- 22; ns <- 100
  for (cfg in list(list("raw", "PR"), list("logit", "PLO"), list("freeman_tukey", "PFT"))) {
    mc <- es_from_prop_single_group(prop = k / ns, n_sample = ns, prop_to_es = cfg[[1]])
    mf <- metafor::escalc(cfg[[2]], xi = k, ni = ns)
    expect_equal(mc$prop,    mf$yi[1],       tolerance = tol, ignore_attr = TRUE,
                 label = paste("prop", cfg[[1]]))
    expect_equal(mc$prop_se, sqrt(mf$vi[1]), tolerance = tol, ignore_attr = TRUE,
                 label = paste("prop_se", cfg[[1]]))
  }
  # prop = 0 must use metafor's (x + 1/2) / (n + 1) convention bit-for-bit
  expect_equal(es_from_prop_single_group(prop = 0, n_sample = 50, prop_to_es = "raw")$prop,
               metafor::escalc("PR", xi = 0, ni = 50, add = 1/2, to = "only0")$yi[1],
               tolerance = tol, ignore_attr = TRUE)
  expect_equal(es_from_prop_single_group(prop = 0, n_sample = 50, prop_to_es = "logit")$prop,
               metafor::escalc("PLO", xi = 0, ni = 50, add = 1/2, to = "only0")$yi[1],
               tolerance = tol, ignore_attr = TRUE)
})
