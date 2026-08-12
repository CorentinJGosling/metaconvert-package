# RR => NNT ---------------
test_that("RR + SE => meta-NNT", {
  dat <- metaumbrella::df.OR
  rr <- with(
    dat,
    metaumbrella:::.estimate_rr_from_n(
      n_cases_exp = n_cases_exp,
      n_cases_nexp = n_cases_nexp,
      n_exp = n_exp, n_nexp = n_nexp
    )
  )
  dat$rr <- rr$value
  dat$logrr_se <- rr$se
  dat$rr_ci_lo <- exp(log(dat$rr) - qnorm(.975) * dat$logrr_se)
  dat$rr_ci_up <- exp(log(dat$rr) + qnorm(.975) * dat$logrr_se)
  z <- log(dat$rr) / dat$logrr_se
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp

  nnt_meta = meta::nnt(x = dat$rr, p.c = dat$baseline_risk, sm = "RR")
  es_mcv_nnt = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_se",
                                  verbose = FALSE), digits=11)

  expect_true(all(es_mcv_nnt$info_used_crude=="rr_se"))
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
})
test_that("RR + CI => meta-NNT", {
  dat <- metaumbrella::df.OR
  rr <- with(
    dat,
    metaumbrella:::.estimate_rr_from_n(
      n_cases_exp = n_cases_exp,
      n_cases_nexp = n_cases_nexp,
      n_exp = n_exp, n_nexp = n_nexp
    )
  )
  dat$rr <- rr$value
  dat$logrr_se <- rr$se
  dat$rr_ci_lo <- exp(log(dat$rr) - qnorm(.975) * dat$logrr_se)
  dat$rr_ci_up <- exp(log(dat$rr) + qnorm(.975) * dat$logrr_se)
  z <- log(dat$rr) / dat$logrr_se
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp

  nnt_meta = meta::nnt(x = dat$rr, p.c = dat$baseline_risk, sm = "RR")
  es_mcv_nnt = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_ci",
                                  verbose = FALSE), digits=11)

  expect_true(all(es_mcv_nnt$info_used_crude=="rr_ci"))
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
})
test_that("RR + pval => meta-NNT", {
  dat <- metaumbrella::df.OR
  rr <- with(
    dat,
    metaumbrella:::.estimate_rr_from_n(
      n_cases_exp = n_cases_exp,
      n_cases_nexp = n_cases_nexp,
      n_exp = n_exp, n_nexp = n_nexp
    )
  )
  dat$rr <- rr$value
  dat$logrr_se <- rr$se
  dat$rr_ci_lo <- exp(log(dat$rr) - qnorm(.975) * dat$logrr_se)
  dat$rr_ci_up <- exp(log(dat$rr) + qnorm(.975) * dat$logrr_se)
  z <- log(dat$rr) / dat$logrr_se
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp

  nnt_meta = meta::nnt(x = dat$rr, p.c = dat$baseline_risk, sm = "RR")
  es_mcv_nnt = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_pval",
                                  verbose = FALSE), digits=11)

  expect_true(all(es_mcv_nnt$info_used_crude=="rr_pval"))
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
})

# OR => NNT ---------------
test_that("OR => meta-NNT", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  nnt_meta = meta::nnt(x = dat$or, p.c = dat$baseline_risk, sm = "OR")
  es_mcv_nnt = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or",
                                  verbose = FALSE), digits=11)

  expect_true(all(es_mcv_nnt$info_used_crude=="or"))
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
})

test_that("OR + SE => meta-NNT", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  nnt_meta = meta::nnt(x = dat$or, p.c = dat$baseline_risk, sm = "OR")
  es_mcv_nnt = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_se",
                                  verbose = FALSE), digits=11)

  expect_true(all(es_mcv_nnt$info_used_crude=="or_se"))
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
})
test_that("RR + CI => meta-NNT", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  nnt_meta = meta::nnt(x = dat$or, p.c = dat$baseline_risk, sm = "OR")
  es_mcv_nnt = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_ci",
                                  verbose = FALSE), digits=11)

  expect_true(all(es_mcv_nnt$info_used_crude=="or_ci"))
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
})
test_that("RR + pval => meta-NNT", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  nnt_meta = meta::nnt(x = dat$or, p.c = dat$baseline_risk, sm = "OR")
  es_mcv_nnt = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_pval",
                                  verbose = FALSE), digits=11)

  expect_true(all(es_mcv_nnt$info_used_crude=="or_pval"))
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
  expect_equal(es_mcv_nnt$es_crude, nnt_meta$NNT, tolerance = 1e-10)
})

# 2x2 to NNT -----
test_that("NNT from 2x2", {
  dat = data.frame(
    n_cases_exp = 40,
    n_controls_exp = 160,
    n_cases_nexp = 10,
    n_controls_nexp = 190
  )

  epinnt = epiR::epi.2by2(dat = matrix(c(40, 10, 160, 190), 2),
                          method = "cross.sectional",
                          conf.level = 0.95, units = 100,
                          interpret = FALSE)
  epinnt$massoc.detail$NNT.strata.score$est
  res = summary(convert_df(dat, measure = "nnt", verbose = FALSE),
                digits = 11)
  expect_equal(res$es_crude,
               -epinnt$massoc.detail$NNT.strata.score$est,
               tolerance = 1e-10)
})

# Internal NNT from 2x2 v OR v RR ------------
test_that("2x2 v OR v RR => NNT", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  rr <- with(
    dat,
    metaumbrella:::.estimate_rr_from_n(
      n_cases_exp = n_cases_exp,
      n_cases_nexp = n_cases_nexp,
      n_exp = n_exp, n_nexp = n_nexp
    )
  )
  dat$rr <- rr$value
  dat$logrr_se <- rr$se
  dat$rr_ci_lo <- exp(log(dat$rr) - qnorm(.975) * dat$logrr_se)
  dat$rr_ci_up <- exp(log(dat$rr) + qnorm(.975) * dat$logrr_se)
  z <- log(dat$rr) / dat$logrr_se
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$rr_pval != 1, ]

  es_mcv_nnt_or = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_se",
                                     verbose = FALSE), digits=11)
  es_mcv_nnt_rr = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_se",
                                     verbose = FALSE), digits=11)
  es_mcv_nnt_2x2 = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="2x2",
                                     verbose = FALSE), digits=11)

  expect_true(all(es_mcv_nnt_or$info_used_crude=="or_se"))
  expect_true(all(es_mcv_nnt_rr$info_used_crude=="rr_se"))
  expect_true(all(es_mcv_nnt_2x2$info_used_crude=="2x2"))

  expect_equal(es_mcv_nnt_or$es_crude, es_mcv_nnt_rr$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_or$es_crude, es_mcv_nnt_2x2$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_2x2$es_crude, es_mcv_nnt_rr$es_crude, tolerance = 1e-10)
})
# REVERSE -----
test_that("NNT - Reverse", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  rr <- with(
    dat,
    metaumbrella:::.estimate_rr_from_n(
      n_cases_exp = n_cases_exp,
      n_cases_nexp = n_cases_nexp,
      n_exp = n_exp, n_nexp = n_nexp
    )
  )
  dat$rr <- rr$value
  dat$logrr_se <- rr$se
  dat$rr_ci_lo <- exp(log(dat$rr) - qnorm(.975) * dat$logrr_se)
  dat$rr_ci_up <- exp(log(dat$rr) + qnorm(.975) * dat$logrr_se)
  z <- log(dat$rr) / dat$logrr_se
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5); dat = dat[dat$rr_pval != 1, ]

  dat$reverse_2x2 = FALSE
  dat$reverse_rr = FALSE
  dat$reverse_rr_pval = FALSE
  dat$reverse_or = FALSE
  dat$reverse_or_pval = FALSE

  es_mcv_nnt_or = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or",
                                        verbose = FALSE), digits=11)
  es_mcv_nnt_or_se = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_se",
                                        verbose = FALSE), digits=11)
  es_mcv_nnt_or_ci = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_ci",
                                        verbose = FALSE), digits=11)
  es_mcv_nnt_or_pval = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_pval",
                                        verbose = FALSE), digits=11)
  es_mcv_nnt_rr_se = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_se",
                                     verbose = FALSE), digits=11)
  es_mcv_nnt_rr_ci = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_ci",
                                     verbose = FALSE), digits=11)
  es_mcv_nnt_rr_pval = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_pval",
                                     verbose = FALSE), digits=11)
  es_mcv_nnt_2x2 = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="2x2",
                                      verbose = FALSE), digits=11)

  dat$reverse_2x2 = TRUE
  dat$reverse_rr = TRUE
  dat$reverse_rr_pval = TRUE
  dat$reverse_or = TRUE
  dat$reverse_or_pval = TRUE

  es_mcv_nnt_or_rv = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or",
                                     verbose = FALSE), digits=11)
  es_mcv_nnt_or_se_rv = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_se",
                                        verbose = FALSE), digits=11)
  es_mcv_nnt_or_ci_rv = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_ci",
                                        verbose = FALSE), digits=11)
  es_mcv_nnt_or_pval_rv = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="or_pval",
                                          verbose = FALSE), digits=11)
  es_mcv_nnt_rr_se_rv = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_se",
                                        verbose = FALSE), digits=11)
  es_mcv_nnt_rr_ci_rv = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_ci",
                                        verbose = FALSE), digits=11)
  es_mcv_nnt_rr_pval_rv = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="rr_pval",
                                          verbose = FALSE), digits=11)
  es_mcv_nnt_2x2_rv = summary(convert_df(dat, measure="nnt", es_selected = "hierarchy", hierarchy="2x2",
                                      verbose = FALSE), digits=11)


  expect_equal(unique(es_mcv_nnt_or$info_used_crude), "or")
  expect_equal(unique(es_mcv_nnt_or_se$info_used_crude), "or_se")
  expect_equal(unique(es_mcv_nnt_or_ci$info_used_crude), "or_ci")
  expect_equal(unique(es_mcv_nnt_or_pval$info_used_crude), "or_pval")
  expect_equal(unique(es_mcv_nnt_rr_se$info_used_crude), "rr_se")
  expect_equal(unique(es_mcv_nnt_rr_ci$info_used_crude), "rr_ci")
  expect_equal(unique(es_mcv_nnt_rr_pval$info_used_crude), "rr_pval")
  expect_equal(unique(es_mcv_nnt_2x2$info_used_crude), "2x2")

  expect_equal(unique(es_mcv_nnt_or_rv$info_used_crude), "or")
  expect_equal(unique(es_mcv_nnt_or_se_rv$info_used_crude), "or_se")
  expect_equal(unique(es_mcv_nnt_or_ci_rv$info_used_crude), "or_ci")
  expect_equal(unique(es_mcv_nnt_or_pval_rv$info_used_crude), "or_pval")
  expect_equal(unique(es_mcv_nnt_rr_se_rv$info_used_crude), "rr_se")
  expect_equal(unique(es_mcv_nnt_rr_ci_rv$info_used_crude), "rr_ci")
  expect_equal(unique(es_mcv_nnt_rr_pval_rv$info_used_crude), "rr_pval")
  expect_equal(unique(es_mcv_nnt_2x2_rv$info_used_crude), "2x2")

  expect_equal(es_mcv_nnt_or$es_crude, -es_mcv_nnt_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_or_se$es_crude, -es_mcv_nnt_or_se_rv$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_or_ci$es_crude, -es_mcv_nnt_or_ci_rv$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_or_pval$es_crude, -es_mcv_nnt_or_pval_rv$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_rr_se$es_crude, -es_mcv_nnt_rr_se_rv$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_rr_ci$es_crude, -es_mcv_nnt_rr_ci_rv$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_rr_pval$es_crude, -es_mcv_nnt_rr_pval_rv$es_crude, tolerance = 1e-10)
  expect_equal(es_mcv_nnt_2x2$es_crude, -es_mcv_nnt_2x2_rv$es_crude, tolerance = 1e-10)
})

