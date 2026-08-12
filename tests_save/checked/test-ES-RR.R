## RR => OR ============
test_that("effectsize SE/CI/p - Grant CI", {
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
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp

  dat$reverse_rr <- FALSE


  or_es <- suppressWarnings(
    log(effectsize::riskratio_to_oddsratio(RR = dat$rr, p0 = dat$baseline_risk))
  )

  es.mcv_or_se <- summary(convert_df(dat,
                                     verbose = FALSE,
                                     es_selected = "hierarchy", hierarchy = "rr_se", measure = "logor",
                                     rr_to_or = "grant"
  ), digits = 11)
  es.mcv_or_ci <- summary(convert_df(dat,
                                     verbose = FALSE,
                                     es_selected = "hierarchy", hierarchy = "rr_ci", measure = "logor",
                                     rr_to_or = "grant"
  ), digits = 11)
  es.mcv_or_p <- summary(convert_df(dat,
                                    verbose = FALSE,
                                    es_selected = "hierarchy", hierarchy = "rr_pval", measure = "logor",
                                    rr_to_or = "grant"
  ), digits = 11)

  row_se = which(!is.nan(or_es) & es.mcv_or_se$info_used_crude=="rr_se")
  row_ci =  which(!is.nan(or_es) & es.mcv_or_ci$info_used_crude=="rr_ci")
  row_p = which(!is.nan(or_es) & es.mcv_or_p$info_used_crude=="rr_pval")

  expect_equal(unique(es.mcv_or_se$info_used_crude[row_se]), "rr_se")
  expect_equal(unique(es.mcv_or_ci$info_used_crude[row_ci]), "rr_ci")
  expect_equal(unique(es.mcv_or_p$info_used_crude[row_p]), "rr_pval")

  expect_equal(es.mcv_or_se$es_crude[row_se],
               or_es[row_se],
               tolerance = 1e-10
  )
  expect_equal(es.mcv_or_ci$es_crude[row_ci],
               or_es[row_ci],
               tolerance = 1e-10
  )
  expect_equal(es.mcv_or_p$es_crude[row_p],
               or_es[row_p],
               tolerance = 1e-10
  )
})


test_that("RR+SE/CI/p to OR - TRANSPOSE 2x2", {
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
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp

  dat$reverse_rr <- FALSE


  es.mcv_rr_se <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se",
                                     measure = "logor",
                                     rr_to_or = "transpose"
  ), digits = 11)
  es.mcv_rr_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci",
                                     measure = "logor",
                                     rr_to_or = "transpose"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval",
                                    measure = "logor",
                                    rr_to_or = "transpose"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_se$info_used_crude), "rr_se")
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "rr_ci")
  expect_equal(unique(es.mcv_rr_p$info_used_crude), "rr_pval")

  expect_equal(es.mcv_rr_se$es_crude, log(dat$rr), tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci$es_crude, log(dat$rr), tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$es_crude, log(dat$rr), tolerance = 1e-10)

  expect_equal(es.mcv_rr_se$se_crude, dat$logrr_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci$se_crude, dat$logrr_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$se_crude, dat$logrr_se, tolerance = 1e-10)
})
test_that("RR+SE/CI/p to RR - DIPIE 2x2 - OR = 1", {
  dat <- metaumbrella::df.OR[1, ]
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
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp

  dat$reverse_rr <- FALSE

  res = estimraw::estim_raw(
    es = dat$rr, lb = dat$rr_ci_lo, ub = dat$rr_ci_up,
    m1 = dat$n_exp, m2 = dat$n_nexp, dec = 11, measure = "rr"
  )
  res
  expect_true(with(res, c/(c+d))[1] == with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp)))

  resor1 = with(res[1,],
                es_from_2x2(n_cases_exp = a, n_controls_exp=b,
                            n_cases_nexp = c, n_controls_nexp=d))

  es.mcv_rr_se <- summary(convert_df(dat,
                                     verbose = FALSE,
                                     es_selected = "hierarchy", hierarchy = "rr_se",
                                     measure = "logor",
                                     rr_to_or = "dipietrantonj"
  ), digits = 11)

  es.mcv_rr_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci",
                                     measure = "logor",
                                     rr_to_or = "dipietrantonj"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval",
                                    measure = "logor",
                                    rr_to_or = "dipietrantonj"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_se$info_used_crude), "rr_se")
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "rr_ci")
  expect_equal(unique(es.mcv_rr_p$info_used_crude), "rr_pval")

  expect_equal(es.mcv_rr_ci$es_crude, resor1$logor, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$es_crude, resor1$logor, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$es_crude, resor1$logor, tolerance = 1e-10)


  expect_equal(es.mcv_rr_ci$se_crude, resor1$logor_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$se_crude, resor1$logor_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$se_crude, resor1$logor_se, tolerance = 1e-10)
})

test_that("OR+SE/CI/p to RR - DIPIE 2x2 - noBR = 1", {
  dat <- metaumbrella::df.OR[1, ]
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
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat$n_cases_nexp = dat$n_cases_exp = NA
  dat$baseline_risk <- NA

  dat$reverse_rr <- FALSE

  res = estimraw::estim_raw(
    es = dat$rr, lb = dat$rr_ci_lo, ub = dat$rr_ci_up,
    m1 = dat$n_exp, m2 = dat$n_nexp, dec = 11, measure = "rr"
  )

  resor1 = with(res[1,],
                es_from_2x2(n_cases_exp = a, n_controls_exp=b,
                            n_cases_nexp = c, n_controls_nexp=d))

  es.mcv_rr_se <- summary(convert_df(dat,
                                     verbose = FALSE,
                                     es_selected = "hierarchy", hierarchy = "rr_se",
                                     measure = "logor",
                                     rr_to_or = "dipietrantonj"
  ), digits = 11)

  es.mcv_rr_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci",
                                     measure = "logor",
                                     rr_to_or = "dipietrantonj"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval",
                                    measure = "logor",
                                    rr_to_or = "dipietrantonj"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_se$info_used_crude), "rr_se")
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "rr_ci")
  expect_equal(unique(es.mcv_rr_p$info_used_crude), "rr_pval")

  expect_equal(es.mcv_rr_ci$es_crude, resor1$logor, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$es_crude, resor1$logor, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$es_crude, resor1$logor, tolerance = 1e-10)


  expect_equal(es.mcv_rr_ci$se_crude, resor1$logor_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$se_crude, resor1$logor_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$se_crude, resor1$logor_se, tolerance = 1e-10)
})


test_that("Internal SE/CI/p - Grant CI", {
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
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp

  dat$reverse_rr <- FALSE

  dat$logrr <- log(dat$rr)
  dat$rr <- NA
  es.mcv_lor_se <- summary(convert_df(dat,
                                      verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", measure = "logor",
                                      rr_to_or = "grant"
  ), digits = 11)
  es.mcv_lor_ci <- summary(convert_df(dat,
                                      verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", measure = "logor",
                                      rr_to_or = "grant"
  ), digits = 11)
  es.mcv_lor_p <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", measure = "logor",
                                     rr_to_or = "grant"
  ), digits = 11)
  expect_equal(es.mcv_lor_se$es_crude, es.mcv_lor_ci$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_lor_se$se_crude, es.mcv_lor_ci$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_lor_se$es_crude, es.mcv_lor_p$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_lor_se$se_crude, es.mcv_lor_p$se_crude, tolerance = 1e-10)

  es.mcv_or_se <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", measure = "or",
                                     rr_to_or = "grant"
  ), digits = 11)
  es.mcv_or_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", measure = "or",
                                     rr_to_or = "grant"
  ), digits = 11)
  es.mcv_or_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", measure = "or",
                                    rr_to_or = "grant"
  ), digits = 11)
  expect_equal(es.mcv_or_ci$es_crude, exp(es.mcv_lor_ci$es_crude), tolerance = 1e-10)
  expect_equal(es.mcv_or_ci$se_crude, es.mcv_lor_ci$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se$es_crude, exp(es.mcv_lor_se$es_crude), tolerance = 1e-10)
  expect_equal(es.mcv_or_se$se_crude, es.mcv_lor_se$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_p$es_crude, exp(es.mcv_lor_p$es_crude), tolerance = 1e-10)
  expect_equal(es.mcv_or_p$se_crude, es.mcv_lor_p$se_crude, tolerance = 1e-10)
})

## Internal ==========
test_that("ES from rr + CI", {
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
  dat$logrr <- log(rr$value)
  dat$logrr_se <- rr$se
  dat$rr_ci_lo <- exp(log(dat$rr) - qnorm(.975) * dat$logrr_se)
  dat$rr_ci_up <- exp(log(dat$rr) + qnorm(.975) * dat$logrr_se)
  z <- log(dat$rr) / dat$logrr_se
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- runif(nrow(dat), 0.1, 0.5)

  dat$reverse_rr <- FALSE
  dat$reverse_rr_pval <- FALSE

  es.mcv_rr_se <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "rr_se", measure = "logrr"
  ), digits = 11)
  es.mcv_rr_ci <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "rr_ci", measure = "logrr"
  ), digits = 11)

  expect_equal(unique(es.mcv_rr_se$info_used_crude), "rr_se")
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "rr_ci")

  expect_equal(es.mcv_rr_se$es_crude, es.mcv_rr_ci$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$se_crude, es.mcv_rr_ci$se_crude, tolerance = 1e-10)

  dat$logrr <- log(dat$rr)
  dat$rr <- NA
  dat$logrr_ci_lo <- log(dat$rr_ci_lo)
  dat$rr_ci_lo <- NA
  dat$logrr_ci_up <- log(dat$rr_ci_up)
  dat$rr_ci_up <- NA


  es.mcv_log <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", measure = "logrr"), digits = 11)
  expect_equal(es.mcv_rr_ci$es_crude, es.mcv_log$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci$se_crude, es.mcv_log$se_crude, tolerance = 1e-10)
})

test_that("ES from rr+p-val", {
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
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- runif(nrow(dat), 0.1, 0.5)

  dat$reverse_rr <- FALSE
  dat$reverse_rr_pval <- FALSE

  es.mcv_rr_se <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "rr_se", measure = "logrr"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "rr_pval", measure = "logrr"
  ), digits = 11)

  row <- which(es.mcv_rr_se$info_used_crude == "rr_se" & es.mcv_rr_p$info_used_crude == "rr_pval")

  expect_equal(unique(es.mcv_rr_se$info_used_crude[row]), "rr_se")
  expect_equal(unique(es.mcv_rr_p$info_used_crude[row]), "rr_pval")

  expect_equal(abs(es.mcv_rr_se$es_crude[row]), abs(es.mcv_rr_p$es_crude[row]), tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$se_crude[row], es.mcv_rr_p$se_crude[row], tolerance = 1e-10)

  dat$logrr <- log(dat$rr)
  dat$rr <- NA
  es.mcv_d_log <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", measure = "logrr"), digits = 11)
  expect_equal(es.mcv_rr_p$es_crude[row], es.mcv_d_log$es_crude[row], tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$se_crude[row], es.mcv_d_log$se_crude[row], tolerance = 1e-10)
})

test_that("RR - Reverse", {
  dat <- metaumbrella::df.OR[1:5, ]
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
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp

  dat$reverse_rr <- FALSE

  es.mcv_rr_se_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", measure = "logrr"), digits = 11)
  es.mcv_rr_ci_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", measure = "logrr"), digits = 11)

  es.mcv_rr_se_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", rr_to_or = "grant", measure = "logor"), digits = 11)
  es.mcv_rr_ci_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", rr_to_or = "grant", measure = "logor"), digits = 11)

  es.mcv_rr_se_rr4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", rr_to_or = "metaumbrella", measure = "logor"), digits = 11)
  es.mcv_rr_ci_rr4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", rr_to_or = "metaumbrella", measure = "logor"), digits = 11)

  es.mcv_rr_se_rr6 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", rr_to_or = "transpose", measure = "logor"), digits = 11)
  es.mcv_rr_ci_rr6 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", rr_to_or = "transpose", measure = "logor"), digits = 11)

  es.mcv_rr_se_rr7 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", rr_to_or = "dipietrantonj", measure = "logor"), digits = 11)
  es.mcv_rr_ci_rr7 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", rr_to_or = "dipietrantonj", measure = "logor"), digits = 11)

  es.mcv_rr_se_nnt <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", measure = "nnt"), digits = 11)
  es.mcv_rr_ci_nnt <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", measure = "nnt"), digits = 11)

  dat$reverse_rr <- TRUE

  es.mcv_rr_se_rr2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", rr_to_or = "grant", measure = "logor"), digits = 11)
  es.mcv_rr_ci_rr2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", rr_to_or = "grant", measure = "logor"), digits = 11)

  es.mcv_rr_se_rr4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", rr_to_or = "metaumbrella", measure = "logor"), digits = 11)
  es.mcv_rr_ci_rr4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", rr_to_or = "metaumbrella", measure = "logor"), digits = 11)

  es.mcv_rr_se_rr6_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", rr_to_or = "transpose", measure = "logor"), digits = 11)
  es.mcv_rr_ci_rr6_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", rr_to_or = "transpose", measure = "logor"), digits = 11)

  es.mcv_rr_se_rr7_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", rr_to_or = "dipietrantonj", measure = "logor"), digits = 11)
  es.mcv_rr_ci_rr7_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", rr_to_or = "dipietrantonj", measure = "logor"), digits = 11)

  es.mcv_rr_se_nnt_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", measure = "nnt"), digits = 11)
  es.mcv_rr_ci_nnt_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", measure = "nnt"), digits = 11)

  es.mcv_rr_se_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_se", measure = "logrr"), digits = 11)
  es.mcv_rr_ci_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_ci", measure = "logrr"), digits = 11)


  expect_true(all(
    c(es.mcv_rr_se_or$info_used_crude, es.mcv_rr_se_or_rv$info_used_crude,
      es.mcv_rr_se_rr2$info_used_crude, es.mcv_rr_se_rr2_rv$info_used_crude,
      es.mcv_rr_se_rr4$info_used_crude, es.mcv_rr_se_rr4_rv$info_used_crude,
      es.mcv_rr_se_rr6$info_used_crude, es.mcv_rr_se_rr6_rv$info_used_crude,
      es.mcv_rr_se_rr7$info_used_crude, es.mcv_rr_se_rr7_rv$info_used_crude,
      es.mcv_rr_se_nnt$info_used_crude, es.mcv_rr_se_nnt_rv$info_used_crude) == "rr_se"))

  expect_true(all(
    c(es.mcv_rr_ci_or$info_used_crude, es.mcv_rr_ci_or_rv$info_used_crude,
      es.mcv_rr_ci_rr2$info_used_crude, es.mcv_rr_ci_rr2_rv$info_used_crude,
      es.mcv_rr_ci_rr4$info_used_crude, es.mcv_rr_ci_rr4_rv$info_used_crude,
      es.mcv_rr_ci_rr6$info_used_crude, es.mcv_rr_ci_rr6_rv$info_used_crude,
      es.mcv_rr_ci_rr7$info_used_crude, es.mcv_rr_ci_rr7_rv$info_used_crude,
      es.mcv_rr_ci_nnt$info_used_crude, es.mcv_rr_ci_nnt_rv$info_used_crude) == "rr_ci"))


  expect_equal(es.mcv_rr_se_rr2$es_crude, -es.mcv_rr_se_rr2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_rr2$es_crude, -es.mcv_rr_ci_rr2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se_rr2$se_crude, es.mcv_rr_se_rr2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_rr2$se_crude, es.mcv_rr_ci_rr2_rv$se_crude, tolerance = 1e-10)


  expect_equal(es.mcv_rr_se_rr4$es_crude, -es.mcv_rr_se_rr4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_rr4$es_crude, -es.mcv_rr_ci_rr4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se_rr4$se_crude, es.mcv_rr_se_rr4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_rr4$se_crude, es.mcv_rr_ci_rr4_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_rr_se_rr6$es_crude, -es.mcv_rr_se_rr6_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_rr6$es_crude, -es.mcv_rr_ci_rr6_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se_rr6$se_crude, es.mcv_rr_se_rr6_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_rr6$se_crude, es.mcv_rr_ci_rr6_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_rr_se_rr7$es_crude, -es.mcv_rr_se_rr7_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_rr7$es_crude, -es.mcv_rr_ci_rr7_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se_rr7$se_crude, es.mcv_rr_se_rr7_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_rr7$se_crude, es.mcv_rr_ci_rr7_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_rr_se_nnt$es_crude, -es.mcv_rr_se_nnt_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_nnt$es_crude, -es.mcv_rr_ci_nnt_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se_nnt$se_crude, es.mcv_rr_se_nnt_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci_nnt$se_crude, es.mcv_rr_ci_nnt_rv$se_crude, tolerance = 1e-10)

})


test_that("RR - pval Reverse", {
  dat <- metaumbrella::df.OR[1:5, ]
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
  dat$rr_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$rr_pval != 1, ]
  dat$baseline_risk <- dat$n_cases_nexp / dat$n_nexp
  dat$reverse_rr <- FALSE
  dat$reverse_rr_pval <- FALSE

  es.mcv_or_pval_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", rr_to_or = "grant", measure = "logor"), digits = 11)
  es.mcv_or_pval_rr4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", rr_to_or = "metaumbrella", measure = "logor"), digits = 11)
  es.mcv_or_pval_rr6 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", rr_to_or = "transpose", measure = "logor"), digits = 11)
  es.mcv_or_pval_rr7 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", rr_to_or = "dipietrantonj", measure = "logor"), digits = 11)
  es.mcv_or_pval_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", measure = "logrr"), digits = 11)
  es.mcv_or_pval_nnt <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", measure = "nnt"), digits = 11)

  dat$reverse_rr_pval <- TRUE

  es.mcv_or_pval_rr2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", rr_to_or = "grant", measure = "logor"), digits = 11)
  es.mcv_or_pval_rr4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", rr_to_or = "metaumbrella", measure = "logor"), digits = 11)
  es.mcv_or_pval_rr6_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", rr_to_or = "transpose", measure = "logor"), digits = 11)
  es.mcv_or_pval_rr7_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", rr_to_or = "dipietrantonj", measure = "logor"), digits = 11)
  es.mcv_or_pval_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", measure = "logrr"), digits = 11)
  es.mcv_or_pval_nnt_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "rr_pval", measure = "nnt"), digits = 11)

  expect_true(all(
    c(es.mcv_or_pval_or$info_used_crude, es.mcv_or_pval_or_rv$info_used_crude,
      es.mcv_or_pval_rr2$info_used_crude, es.mcv_or_pval_rr2_rv$info_used_crude,
      es.mcv_or_pval_rr4$info_used_crude, es.mcv_or_pval_rr4_rv$info_used_crude,
      es.mcv_or_pval_rr6$info_used_crude, es.mcv_or_pval_rr6_rv$info_used_crude,
      es.mcv_or_pval_rr7$info_used_crude, es.mcv_or_pval_rr7_rv$info_used_crude,
      es.mcv_or_pval_nnt$info_used_crude, es.mcv_or_pval_nnt_rv$info_used_crude) == "rr_pval"))

  expect_equal(es.mcv_or_pval_or$es_crude, -es.mcv_or_pval_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_or$se_crude, es.mcv_or_pval_or_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_pval_rr2$es_crude, -es.mcv_or_pval_rr2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr4$es_crude, -es.mcv_or_pval_rr4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr6$es_crude, -es.mcv_or_pval_rr6_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr7$es_crude, -es.mcv_or_pval_rr7_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_nnt$es_crude, -es.mcv_or_pval_nnt_rv$es_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_pval_rr2$se_crude, es.mcv_or_pval_rr2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr4$se_crude, es.mcv_or_pval_rr4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr6$se_crude, es.mcv_or_pval_rr6_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr7$se_crude, es.mcv_or_pval_rr7_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_nnt$se_crude, es.mcv_or_pval_nnt_rv$se_crude, tolerance = 1e-10)

})
