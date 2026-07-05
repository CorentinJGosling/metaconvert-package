# OR to SMD ----
test_that("1. OR+SE to SMD", {
  skip_if_not_installed("esc")
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se

  es.mcv_d <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "or_se",
    measure = "d"
  ), digits = 11)
  expect_equal(unique(es.mcv_d$info_used_crude), "or_se")

  # test OR esc (scalar-only API, loop row-by-row) ---
  esc_or2d <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r <- esc::convert_or2d(or = dat$or[i], se = dat$logor_se[i],
                           totaln = dat$n_cases[i] + dat$n_controls[i])
    data.frame(es = r$es, se = r$se, stringsAsFactors = FALSE)
  }))
  expect_equal(es.mcv_d$es_crude, esc_or2d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, esc_or2d$se, tolerance = 1e-10)

  # test logOR esc ---
  dat$logor <- log(dat$or)
  dat$or <- NA
  es.mcv_d_log <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "d"), digits = 11)
  expect_equal(es.mcv_d$es_crude, es.mcv_d_log$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d_log$se_crude, tolerance = 1e-10)

  # test logOR mfr ---
  mfr <- metafor::escalc(
    ai = n_cases_exp, bi = n_controls_exp,
    ci = n_cases_nexp, di = n_controls_nexp,
    measure = "OR2DL", data = dat
  )
  expect_equal(es.mcv_d$es_crude, as.numeric(as.character(mfr$yi)), tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, sqrt(mfr$vi), tolerance = 1e-10)

  # test compute.es ---
  comp_d <- compute.es::lores(
    lor = dat$logor,
    var.lor = dat$logor_se^2,
    verbose = FALSE,
    n.1 = dat$n_exp, n.2 = dat$n_nexp,
    dig = 11
  )
  expect_equal(es.mcv_d$es_crude, comp_d$d, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude^2, comp_d$var.d, tolerance = 1e-10)
})

# OR to COR ----
test_that("OR+SE to R (lipsey_cooper 1)", {
  skip_if_not_installed("compute.es")
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
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))
  dat$n_sample <- dat$n_cases + dat$n_controls
  dat$n_exp <- dat$n_nexp <- dat$n_cases_exp <- dat$n_controls_exp <-
    dat$n_cases_nexp <- dat$n_controls_nexp <- dat$n_cases <- dat$n_controls <- NA
  dat$n_exp <- dat$n_nexp <- dat$n_sample / 2
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
    measure = "r",
    or_to_cor = "lipsey_cooper"
  ), digits = 12)

  expect_equal(unique(es.mcv_r$info_used_crude), "or_se")

  comp_r <- compute.es::lores(
    lor = log(dat$or),
    var.lor = dat$logor_se^2,
    verbose = FALSE,
    n.1 = dat$n_exp, n.2 = dat$n_nexp,
    dig = 12
  )
  expect_equal(es.mcv_r$es_crude, comp_r$r, tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude^2, comp_r$var.r, tolerance = 1e-10)
})
test_that("OR+SE to Z (lipsey_cooper 1)", {
  skip_if_not_installed("compute.es")
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
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))
  dat$n_sample <- dat$n_cases + dat$n_controls
  dat$n_exp <- dat$n_nexp <- dat$n_cases_exp <- dat$n_controls_exp <-
    dat$n_cases_nexp <- dat$n_controls_nexp <- dat$n_cases <- dat$n_controls <- NA
  dat$n_exp <- dat$n_nexp <- dat$n_sample / 2
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
    measure = "z",
    or_to_cor = "lipsey_cooper"
  ), digits = 12)

  expect_equal(unique(es.mcv_r$info_used_crude), "or_se")

  comp_r <- compute.es::lores(
    lor = log(dat$or),
    var.lor = dat$logor_se^2, verbose = FALSE,
    n.1 = dat$n_exp, n.2 = dat$n_nexp,
    dig = 12
  )
  expect_equal(es.mcv_r$es_crude, comp_r$fisher.z, tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude^2, comp_r$var.z, tolerance = 5e-2)
})

test_that("OR+SE to R (lipsey_cooper 2)", {
  skip_if_not_installed("esc")
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
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))
  dat$n_sample <- dat$n_cases + dat$n_controls
  dat$n_exp <- dat$n_nexp <- dat$n_cases_exp <- dat$n_controls_exp <-
    dat$n_cases_nexp <- dat$n_controls_nexp <- dat$n_cases <- dat$n_controls <- NA
  dat$n_exp <- dat$n_nexp <- dat$n_sample / 2
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
    measure = "r",
    or_to_cor = "lipsey_cooper"
  ), digits = 11)

  expect_equal(unique(es.mcv_r$info_used_crude), "or_se")

  # esc scalar-only API, loop row-by-row
  esc_d2r <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    or2d <- esc::convert_or2d(or = dat$or[i], se = dat$logor_se[i],
                              totaln = dat$n_exp[i] + dat$n_nexp[i], es.type = "d")
    d2r <- esc::convert_d2r(d = or2d$es, se = or2d$se,
                            grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i])
    data.frame(es = d2r$es, se = d2r$se, stringsAsFactors = FALSE)
  }))
  expect_equal(es.mcv_r$es_crude, esc_d2r$es, tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, esc_d2r$se, tolerance = 5e-2)
})

test_that("OR+SE to Z (lipsey_cooper 2)", {
  skip_if_not_installed("esc")
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
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))
  dat$n_sample <- dat$n_cases + dat$n_controls
  dat$n_exp <- dat$n_nexp <- dat$n_cases_exp <- dat$n_controls_exp <-
    dat$n_cases_nexp <- dat$n_controls_nexp <- dat$n_cases <- dat$n_controls <- NA
  dat$n_exp <- dat$n_nexp <- dat$n_sample / 2
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
    measure = "z",
    or_to_cor = "lipsey_cooper"
  ), digits = 11)

  expect_equal(unique(es.mcv_r$info_used_crude), "or_se")

  # esc scalar-only API, loop row-by-row
  esc_d2r <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    or2d <- esc::convert_or2d(or = dat$or[i], se = dat$logor_se[i],
                              totaln = dat$n_exp[i] + dat$n_nexp[i], es.type = "d")
    d2r <- esc::convert_d2r(d = or2d$es, se = or2d$se,
                            grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i])
    data.frame(es = d2r$es, se = d2r$se, stringsAsFactors = FALSE)
  }))
  expect_equal(es.mcv_r$es_crude, atanh(esc_d2r$es), tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, esc_d2r$se, tolerance = 1e-10)
})

test_that("OR+SE to R (digby)", {
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
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))
  dat$n_sample <- dat$n_cases + dat$n_controls
  dat$n_exp <- dat$n_nexp <- dat$n_cases_exp <- dat$n_controls_exp <-
    dat$n_cases_nexp <- dat$n_controls_nexp <- dat$n_cases <- dat$n_controls <- NA
  dat$n_exp <- dat$n_nexp <- dat$n_sample / 2

  res <- metafor::conv.delta(log(dat$or), dat$logor_se^2,
    transf = metafor::transf.lnortortet.digby,
    var.names = c("yi", "vi"), digits=11
  )

  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
    measure = "r",
    or_to_cor = "digby"
  ), digits = 11)

  expect_equal(unique(es.mcv_r$info_used_crude), "or_se")
  expect_equal(es.mcv_r$es_crude, as.numeric(res$yi), tolerance = 1e-8)
  expect_equal(es.mcv_r$se_crude^2, res$vi, tolerance = 1e-8)
})
test_that("OR+SE to R (pearson)", {
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
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))
  dat$n_sample <- dat$n_cases + dat$n_controls
  dat$n_exp <- dat$n_nexp <- dat$n_cases_exp <- dat$n_controls_exp <-
    dat$n_cases_nexp <- dat$n_controls_nexp <- dat$n_cases <- dat$n_controls <- NA
  dat$n_exp <- dat$n_nexp <- dat$n_sample / 2

  res <- metafor::conv.delta(log(dat$or), dat$logor_se^2,
                             transf = metafor::transf.lnortortet.pearson,
                             var.names = c("yi", "vi")
  )

  es.mcv_r <- summary(convert_df(dat,
                                 verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
                                 measure = "r",
                                 or_to_cor = "pearson"
  ), digits = 11)

  expect_equal(unique(es.mcv_r$info_used_crude), "or_se")
  expect_equal(es.mcv_r$es_crude, as.numeric(res$yi), tolerance = 1e-8)
  expect_equal(es.mcv_r$se_crude, sqrt(res$vi), tolerance = 1e-8)
})


test_that("OR+SE to R (bonett)", {
  dat = data.frame(study=1)
  dat$n_cases_exp = 143+0.5
  dat$n_controls_exp = 52+0.5
  dat$n_cases_nexp = 41+0.5
  dat$n_controls_nexp = 164+0.5
  dat$n_exp = dat$n_cases_exp + dat$n_controls_exp
  dat$n_nexp = dat$n_cases_nexp + dat$n_controls_nexp
  dat$n_cases = dat$n_cases_exp + dat$n_cases_nexp
  dat$n_controls = dat$n_controls_exp + dat$n_controls_nexp
  dat$n_sample = dat$n_exp + dat$n_cases_nexp + dat$n_controls_nexp
    dat$small_margin_prop = with(dat, min(
    n_cases/n_sample, n_controls/n_sample,
    n_exp/n_sample, n_nexp/n_sample
  ))
  dat$or = with(dat,
        suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp)))
  dat$logor_se = with(dat, suppressWarnings(
    sqrt(1/n_cases_exp +  1/n_controls_exp +
           1/n_cases_nexp + 1/n_controls_nexp)))

  r = .741
  r_se = .0446

  es.mcv_r <- summary(convert_df(dat,
                                 verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
                                 measure = "r",
                                 or_to_cor = "bonett"
  ), digits = 11)



  expect_equal(unique(es.mcv_r$info_used_crude), "or_se")
  expect_equal(es.mcv_r$es_crude, r, tolerance = 1e-3)
  expect_equal(es.mcv_r$se_crude, r_se, tolerance = 1e-4)
})


# OR to RR ----
test_that("OR+SE/CI/p to RR  - GRANT CI", {
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
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  es.mcv_rr_orG2 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or",
    measure = "logrr",
    or_to_rr = "grant"
  ), digits = 11)
  es.mcv_rr_seG2 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
    measure = "logrr",
    or_to_rr = "grant"
  ), digits = 11)
  es.mcv_rr_ciG2 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci",
    measure = "logrr",
    or_to_rr = "grant"
  ), digits = 11)
  es.mcv_rr_pG2 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval",
    measure = "logrr",
    or_to_rr = "grant"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_orG2$info_used_crude), "or")
  expect_equal(unique(es.mcv_rr_seG2$info_used_crude), "or_se")
  expect_equal(unique(es.mcv_rr_ciG2$info_used_crude), "or_ci")
  expect_equal(unique(es.mcv_rr_pG2$info_used_crude), "or_pval")

  expect_equal(es.mcv_rr_seG2$es_crude, es.mcv_rr_orG2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_seG2$es_crude, es.mcv_rr_ciG2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr_seG2$es_crude, es.mcv_rr_pG2$es_crude, tolerance = 1e-10)
})

test_that("OR+SE/CI/p to RR - TRANSPOSE 2x2", {
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
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  es.mcv_rr_or <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "or",
                                     measure = "logrr",
                                     or_to_rr = "transpose"
  ), digits = 11)
  es.mcv_rr_se <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se",
                                     measure = "logrr",
                                     or_to_rr = "transpose"
  ), digits = 11)
  es.mcv_rr_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci",
                                     measure = "logrr",
                                     or_to_rr = "transpose"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval",
                                    measure = "logrr",
                                    or_to_rr = "transpose"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_or$info_used_crude), "or")
  expect_equal(unique(es.mcv_rr_se$info_used_crude), "or_se")
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "or_ci")
  expect_equal(unique(es.mcv_rr_p$info_used_crude), "or_pval")

  expect_equal(es.mcv_rr_or$es_crude, log(dat$or), tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$es_crude, log(dat$or), tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci$es_crude, log(dat$or), tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$es_crude, log(dat$or), tolerance = 1e-10)

  # expect_equal(es.mcv_rr_or$se_crude, dat$logor_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$se_crude, dat$logor_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_ci$se_crude, dat$logor_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$se_crude, dat$logor_se, tolerance = 1e-10)
})
test_that("OR+SE/CI/p to RR - DIPIE 2x2 - BR = 1", {
  dat <- metaumbrella::df.OR[2, ]
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  res = estimraw::estim_raw(
    es = dat$or, lb = dat$or_ci_lo, ub = dat$or_ci_up,
    m1 = dat$n_exp, m2 = dat$n_nexp, dec = 11, measure = "or"
  )
  expect_true(with(res$solution1[1,],
                   c/(c+d)) == dat$n_cases_nexp/(dat$n_cases_nexp+dat$n_controls_nexp))

  resor1 = with(res$solution1[1,],
                es_from_2x2(n_cases_exp = a, n_controls_exp=b,
                            n_cases_nexp = c, n_controls_nexp=d))

  es.mcv_rr_se <- summary(convert_df(dat,
                                     verbose = FALSE,
                                     es_selected = "hierarchy", hierarchy = "or_se",
                                     measure = "logrr",
                                     or_to_rr = "dipietrantonj"
  ), digits = 11)

  es.mcv_rr_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci",
                                     measure = "logrr",
                                     or_to_rr = "dipietrantonj"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval",
                                    measure = "logrr",
                                    or_to_rr = "dipietrantonj"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "or_ci")
  expect_equal(unique(es.mcv_rr_p$info_used_crude), "or_pval")
  expect_equal(unique(es.mcv_rr_se$info_used_crude), "or_se")

  expect_equal(es.mcv_rr_ci$es_crude, resor1$logrr, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$es_crude, resor1$logrr, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$es_crude, resor1$logrr, tolerance = 1e-10)


  expect_equal(es.mcv_rr_ci$se_crude, resor1$logrr_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$se_crude, resor1$logrr_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$se_crude, resor1$logrr_se, tolerance = 1e-10)
})
test_that("OR+SE/CI/p to RR - DIPIE 2x2 - BR = 2", {
  dat <- metaumbrella::df.OR[1, ]
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp))

  res = estimraw::estim_raw(
    es = dat$or, lb = dat$or_ci_lo, ub = dat$or_ci_up,
    m1 = dat$n_exp, m2 = dat$n_nexp, dec = 11, measure = "or"
  )
  expect_true(with(res$solution2[1,], c/(c+d)) == with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp)))

  resor2 = with(res$solution2[1,],
                es_from_2x2(n_cases_exp = a, n_controls_exp=b,
                            n_cases_nexp = c, n_controls_nexp=d))

  es.mcv_rr_se <- summary(convert_df(dat,
                                     verbose = FALSE,
                                     es_selected = "hierarchy", hierarchy = "or_se",
                                     measure = "logrr",
                                     or_to_rr = "dipietrantonj"
  ), digits = 11)

  es.mcv_rr_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci",
                                     measure = "logrr",
                                     or_to_rr = "dipietrantonj"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval",
                                    measure = "logrr",
                                    or_to_rr = "dipietrantonj"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "or_ci")
  expect_equal(unique(es.mcv_rr_p$info_used_crude), "or_pval")
  expect_equal(unique(es.mcv_rr_se$info_used_crude), "or_se")

  expect_equal(es.mcv_rr_ci$es_crude, resor2$logrr, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$es_crude, resor2$logrr, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$es_crude, resor2$logrr, tolerance = 1e-10)


  expect_equal(es.mcv_rr_ci$se_crude, resor2$logrr_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$se_crude, resor2$logrr_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$se_crude, resor2$logrr_se, tolerance = 1e-10)
})

test_that("OR+SE/CI/p to RR - DIPIE 2x2 - noBR = 1", {
  dat <- metaumbrella::df.OR[2, ]
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- NA
  dat$n_cases_exp = dat$n_cases_nexp = dat$n_controls_exp = dat$n_controls_nexp = NA

  res = estimraw::estim_raw(
    es = dat$or, lb = dat$or_ci_lo, ub = dat$or_ci_up,
    m1 = dat$n_exp, m2 = dat$n_nexp, dec = 11, measure = "or"
  )
  # expect_true(with(res$solution1[1,], c/(c+d)) == with(dat, n_cases_nexp/(n_cases_nexp+n_controls_nexp)))

  resor1 = with(res$solution1[1,],
                es_from_2x2(n_cases_exp = a, n_controls_exp=b,
                            n_cases_nexp = c, n_controls_nexp=d))

  es.mcv_rr_se <- summary(convert_df(dat,
                                     verbose = FALSE,
                                     es_selected = "hierarchy", hierarchy = "or_se",
                                     measure = "logrr",
                                     or_to_rr = "dipietrantonj"
  ), digits = 11)

  es.mcv_rr_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci",
                                     measure = "logrr",
                                     or_to_rr = "dipietrantonj"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval",
                                    measure = "logrr",
                                    or_to_rr = "dipietrantonj"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "or_ci")
  expect_equal(unique(es.mcv_rr_p$info_used_crude), "or_pval")
  expect_equal(unique(es.mcv_rr_se$info_used_crude), "or_se")

  expect_equal(es.mcv_rr_ci$es_crude, resor1$logrr, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$es_crude, resor1$logrr, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$es_crude, resor1$logrr, tolerance = 1e-10)


  expect_equal(es.mcv_rr_ci$se_crude, resor1$logrr_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_p$se_crude, resor1$logrr_se, tolerance = 1e-10)
  expect_equal(es.mcv_rr_se$se_crude, resor1$logrr_se, tolerance = 1e-10)
})
test_that("OR+SE/CI/p to RR - DIPIE 2x2 - no BR = 2", {
  dat <- metaumbrella::df.OR[1, ]
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$baseline_risk <- NA
  dat$n_cases_exp = dat$n_cases_nexp = dat$n_controls_exp = dat$n_controls_nexp = NA

  res = estimraw::estim_raw(
    es = dat$or, lb = dat$or_ci_lo, ub = dat$or_ci_up,
    m1 = dat$n_exp, m2 = dat$n_nexp, dec = 11, measure = "or"
  )

  resor2 = with(res$solution2[1,],
                es_from_2x2(n_cases_exp = a, n_controls_exp=b,
                            n_cases_nexp = c, n_controls_nexp=d))

  es.mcv_rr_se <- summary(convert_df(dat,
                                     verbose = FALSE,
                                     es_selected = "hierarchy", hierarchy = "or_se",
                                     measure = "logrr",
                                     or_to_rr = "dipietrantonj"
  ), digits = 11)

  es.mcv_rr_ci <- summary(convert_df(dat,
                                     verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci",
                                     measure = "logrr",
                                     or_to_rr = "dipietrantonj"
  ), digits = 11)
  es.mcv_rr_p <- summary(convert_df(dat,
                                    verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval",
                                    measure = "logrr",
                                    or_to_rr = "dipietrantonj"
  ), digits = 11)
  expect_equal(unique(es.mcv_rr_ci$info_used_crude), "or_ci")
  expect_equal(unique(es.mcv_rr_p$info_used_crude), "or_pval")
  expect_equal(unique(es.mcv_rr_se$info_used_crude), "or_se")

  expect_false(es.mcv_rr_ci$es_crude == resor2$logrr)
  expect_false(es.mcv_rr_p$es_crude==resor2$logrr)
  expect_false(es.mcv_rr_se$es_crude==resor2$logrr)

  expect_false(es.mcv_rr_ci$se_crude==resor2$logrr_se)
  expect_false(es.mcv_rr_p$se_crude==resor2$logrr_se)
  expect_false(es.mcv_rr_se$se_crude==resor2$logrr_se)
})

# INTERNAL VERIFICATIONS ----
test_that("OR+SE = OR+CI", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)

  # compare with SE from metafor
  mfr_se <- metafor::conv.wald(
    out = log(dat$or),
    ci.lb = log(dat$or_ci_lo),
    ci.ub = log(dat$or_ci_up)
  )

  expect_equal(dat$logor_se, sqrt(mfr_se$vi), tolerance = 1e-10)

  es.mcv_or_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "logor"), digits = 11)
  es.mcv_or_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", measure = "logor"), digits = 11)
  expect_equal(unique(es.mcv_or_se$info_used_crude), "or_se")
  expect_equal(unique(es.mcv_or_ci$info_used_crude), "or_ci")

  expect_equal(es.mcv_or_se$es_crude, es.mcv_or_ci$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se$se_crude, es.mcv_or_ci$se_crude, tolerance = 1e-10)

  dat$logor <- log(dat$or)
  dat$or <- NA
  dat$logor_ci_lo <- log(dat$or_ci_lo)
  dat$or_ci_lo <- NA
  dat$logor_ci_up <- log(dat$or_ci_up)
  dat$or_ci_up <- NA

  es.mcv_exp <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", measure = "or"), digits = 11)
  expect_equal(es.mcv_or_ci$es_crude, log(es.mcv_exp$es_crude), tolerance = 1e-10)
  expect_equal(es.mcv_or_ci$se_crude, es.mcv_exp$se_crude, tolerance = 1e-10)
})

test_that("OR+SE = OR+p-val", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]

  # compare with pval from metafor
  mfr_se <- metafor::conv.wald(
    out = log(dat$or),
    pval = dat$or_pval
  )

  expect_equal(dat$logor_se, sqrt(mfr_se$vi), tolerance = 1e-10)

  es.mcv_or_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "logor"), digits = 11)
  es.mcv_or_p <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", measure = "logor"), digits = 11)
  expect_equal(unique(es.mcv_or_se$info_used_crude), "or_se")
  expect_equal(unique(es.mcv_or_p$info_used_crude), "or_pval")

  expect_equal(abs(es.mcv_or_se$es_crude), abs(es.mcv_or_p$es_crude), tolerance = 1e-10)
  expect_equal(es.mcv_or_se$se_crude, es.mcv_or_p$se_crude, tolerance = 1e-10)

  dat$logor <- log(dat$or)
  dat$or <- NA
  es.mcv_d_exp <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", measure = "or"), digits = 11)
  expect_equal(es.mcv_or_p$es_crude, log(es.mcv_d_exp$es_crude), tolerance = 1e-10)
  expect_equal(es.mcv_or_p$se_crude, es.mcv_d_exp$se_crude, tolerance = 1e-10)
})


# REVERSE ----

test_that("OR - Reverse", {
  dat <- metaumbrella::df.OR[1:5, ]
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$small_margin_prop = 0.10
  dat$baseline_risk = 0.10
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$reverse_or <- FALSE
  dat$reverse_or_pval <- FALSE

  es.mcv_or_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", measure = "d"), digits = 11)
  es.mcv_or_se_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "d"), digits = 11)
  es.mcv_or_ci_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", measure = "d"), digits = 11)

  es.mcv_or_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", measure = "g"), digits = 11)
  es.mcv_or_se_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "g"), digits = 11)
  es.mcv_or_ci_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", measure = "g"), digits = 11)

  es.mcv_or_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "lipsey_cooper", measure = "r"), digits = 11)
  es.mcv_or_se_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "lipsey_cooper", measure = "r"), digits = 11)
  es.mcv_or_ci_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "lipsey_cooper", measure = "r"), digits = 11)

  es.mcv_or_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "lipsey_cooper", measure = "z"), digits = 11)
  es.mcv_or_se_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "lipsey_cooper", measure = "z"), digits = 11)
  es.mcv_or_ci_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "lipsey_cooper", measure = "z"), digits = 11)

  es.mcv_or_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "digby", measure = "r"), digits = 11)
  es.mcv_or_se_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "digby", measure = "r"), digits = 11)
  es.mcv_or_ci_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "digby", measure = "r"), digits = 11)

  es.mcv_or_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "digby", measure = "z"), digits = 11)
  es.mcv_or_se_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "digby", measure = "z"), digits = 11)
  es.mcv_or_ci_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "digby", measure = "z"), digits = 11)

  es.mcv_or_r3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "pearson", measure = "r"), digits = 11)
  es.mcv_or_se_r3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "pearson", measure = "r"), digits = 11)
  es.mcv_or_ci_r3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "pearson", measure = "r"), digits = 11)

  es.mcv_or_z3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "pearson", measure = "z"), digits = 11)
  es.mcv_or_se_z3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "pearson", measure = "z"), digits = 11)
  es.mcv_or_ci_z3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "pearson", measure = "z"), digits = 11)

  es.mcv_or_r4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "bonett", measure = "r"), digits = 11)
  es.mcv_or_se_r4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "bonett", measure = "r"), digits = 11)
  es.mcv_or_ci_r4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "bonett", measure = "r"), digits = 11)

  es.mcv_or_z4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "bonett", measure = "z"), digits = 11)
  es.mcv_or_se_z4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "bonett", measure = "z"), digits = 11)
  es.mcv_or_ci_z4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "bonett", measure = "z"), digits = 11)

  es.mcv_or_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", measure = "logor"), digits = 11)
  es.mcv_or_se_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "logor"), digits = 11)
  es.mcv_or_ci_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", measure = "logor"), digits = 11)



  es.mcv_or_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "grant", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "grant", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "grant", measure = "logrr"), digits = 11)


  es.mcv_or_rr4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "metaumbrella_cases", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "metaumbrella_cases", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "metaumbrella_cases", measure = "logrr"), digits = 11)

  es.mcv_or_rr5 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "metaumbrella_exp", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr5 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "metaumbrella_exp", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr5 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "metaumbrella_exp", measure = "logrr"), digits = 11)

  es.mcv_or_rr6 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "transpose", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr6 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "transpose", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr6 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "transpose", measure = "logrr"), digits = 11)

  es.mcv_or_rr7 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "dipietrantonj", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr7 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "dipietrantonj", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr7 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "dipietrantonj", measure = "logrr"), digits = 11)

  dat$reverse_or <- TRUE


  es.mcv_or_rr2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or",   or_to_rr = "grant",    measure = "logrr"), digits = 11)
  es.mcv_or_se_rr2_rv <- summary(convert_df(dat, verbose = FALSE,   es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "grant", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "grant", measure = "logrr"), digits = 11)

  es.mcv_or_rr4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "metaumbrella_cases", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "metaumbrella_cases", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "metaumbrella_cases", measure = "logrr"), digits = 11)

  es.mcv_or_rr5_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "metaumbrella_exp", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr5_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "metaumbrella_exp", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr5_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "metaumbrella_exp", measure = "logrr"), digits = 11)

  es.mcv_or_rr6_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "transpose", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr6_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "transpose", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr6_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "transpose", measure = "logrr"), digits = 11)

  es.mcv_or_rr7_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_rr = "dipietrantonj", measure = "logrr"), digits = 11)
  es.mcv_or_se_rr7_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_rr = "dipietrantonj", measure = "logrr"), digits = 11)
  es.mcv_or_ci_rr7_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_rr = "dipietrantonj", measure = "logrr"), digits = 11)

  es.mcv_or_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", measure = "d"), digits = 11)
  es.mcv_or_se_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "d"), digits = 11)
  es.mcv_or_ci_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", measure = "d"), digits = 11)

  es.mcv_or_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", measure = "g"), digits = 11)
  es.mcv_or_se_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "g"), digits = 11)
  es.mcv_or_ci_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", measure = "g"), digits = 11)

  es.mcv_or_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "lipsey_cooper", measure = "r"), digits = 11)
  es.mcv_or_se_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "lipsey_cooper", measure = "r"), digits = 11)
  es.mcv_or_ci_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "lipsey_cooper", measure = "r"), digits = 11)

  es.mcv_or_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "lipsey_cooper", measure = "z"), digits = 11)
  es.mcv_or_se_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "lipsey_cooper", measure = "z"), digits = 11)
  es.mcv_or_ci_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "lipsey_cooper", measure = "z"), digits = 11)

  es.mcv_or_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "digby", measure = "r"), digits = 11)
  es.mcv_or_se_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "digby", measure = "r"), digits = 11)
  es.mcv_or_ci_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "digby", measure = "r"), digits = 11)

  es.mcv_or_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "digby", measure = "z"), digits = 11)
  es.mcv_or_se_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "digby", measure = "z"), digits = 11)
  es.mcv_or_ci_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "digby", measure = "z"), digits = 11)

  es.mcv_or_r3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "pearson", measure = "r"), digits = 11)
  es.mcv_or_se_r3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "pearson", measure = "r"), digits = 11)
  es.mcv_or_ci_r3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "pearson", measure = "r"), digits = 11)

  es.mcv_or_z3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "pearson", measure = "z"), digits = 11)
  es.mcv_or_se_z3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "pearson", measure = "z"), digits = 11)
  es.mcv_or_ci_z3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "pearson", measure = "z"), digits = 11)

  es.mcv_or_r4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "bonett", measure = "r"), digits = 11)
  es.mcv_or_se_r4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "bonett", measure = "r"), digits = 11)
  es.mcv_or_ci_r4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "bonett", measure = "r"), digits = 11)

  es.mcv_or_z4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", or_to_cor = "bonett", measure = "z"), digits = 11)
  es.mcv_or_se_z4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", or_to_cor = "bonett", measure = "z"), digits = 11)
  es.mcv_or_ci_z4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", or_to_cor = "bonett", measure = "z"), digits = 11)

  es.mcv_or_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or", measure = "logor"), digits = 11)
  es.mcv_or_se_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_se", measure = "logor"), digits = 11)
  es.mcv_or_ci_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_ci", measure = "logor"), digits = 11)


  expect_true(all(
    c(es.mcv_or_d$info_used_crude, es.mcv_or_d_rv$info_used_crude,
      es.mcv_or_g$info_used_crude, es.mcv_or_g_rv$info_used_crude,
      es.mcv_or_or$info_used_crude, es.mcv_or_or_rv$info_used_crude,

      es.mcv_or_r1$info_used_crude, es.mcv_or_r1_rv$info_used_crude,
      es.mcv_or_r2$info_used_crude, es.mcv_or_r2_rv$info_used_crude,
      es.mcv_or_r3$info_used_crude, es.mcv_or_r3_rv$info_used_crude,
      es.mcv_or_r4$info_used_crude, es.mcv_or_r4_rv$info_used_crude,
      es.mcv_or_z1$info_used_crude, es.mcv_or_z1_rv$info_used_crude,
      es.mcv_or_z2$info_used_crude, es.mcv_or_z2_rv$info_used_crude,
      es.mcv_or_z3$info_used_crude, es.mcv_or_z3_rv$info_used_crude,
      es.mcv_or_z4$info_used_crude, es.mcv_or_z4_rv$info_used_crude,

      es.mcv_or_rr2$info_used_crude, es.mcv_or_rr2_rv$info_used_crude,
      es.mcv_or_rr4$info_used_crude, es.mcv_or_rr4_rv$info_used_crude,
      es.mcv_or_rr5$info_used_crude, es.mcv_or_rr5_rv$info_used_crude,
      es.mcv_or_rr6$info_used_crude, es.mcv_or_rr6_rv$info_used_crude,
      es.mcv_or_rr7$info_used_crude, es.mcv_or_rr7_rv$info_used_crude) == "or"))
  expect_true(all(
    c(es.mcv_or_se_d$info_used_crude, es.mcv_or_se_d_rv$info_used_crude,
      es.mcv_or_se_g$info_used_crude, es.mcv_or_se_g_rv$info_used_crude,
      es.mcv_or_se_or$info_used_crude, es.mcv_or_se_or_rv$info_used_crude,

      es.mcv_or_se_r1$info_used_crude, es.mcv_or_se_r1_rv$info_used_crude,
      es.mcv_or_se_r2$info_used_crude, es.mcv_or_se_r2_rv$info_used_crude,
      es.mcv_or_se_r3$info_used_crude, es.mcv_or_se_r3_rv$info_used_crude,
      es.mcv_or_se_r4$info_used_crude, es.mcv_or_se_r4_rv$info_used_crude,
      es.mcv_or_se_z1$info_used_crude, es.mcv_or_se_z1_rv$info_used_crude,
      es.mcv_or_se_z2$info_used_crude, es.mcv_or_se_z2_rv$info_used_crude,
      es.mcv_or_se_z3$info_used_crude, es.mcv_or_se_z3_rv$info_used_crude,
      es.mcv_or_se_z4$info_used_crude, es.mcv_or_se_z4_rv$info_used_crude,

      es.mcv_or_se_rr2$info_used_crude, es.mcv_or_se_rr2_rv$info_used_crude,
      es.mcv_or_se_rr4$info_used_crude, es.mcv_or_se_rr4_rv$info_used_crude,
      es.mcv_or_se_rr5$info_used_crude, es.mcv_or_se_rr5_rv$info_used_crude,
      es.mcv_or_se_rr6$info_used_crude, es.mcv_or_se_rr6_rv$info_used_crude,
      es.mcv_or_se_rr7$info_used_crude, es.mcv_or_se_rr7_rv$info_used_crude) == "or_se"))
  expect_true(all(
    c(es.mcv_or_ci_d$info_used_crude, es.mcv_or_ci_d_rv$info_used_crude,
      es.mcv_or_ci_g$info_used_crude, es.mcv_or_ci_g_rv$info_used_crude,
      es.mcv_or_ci_or$info_used_crude, es.mcv_or_ci_or_rv$info_used_crude,

      es.mcv_or_ci_r1$info_used_crude, es.mcv_or_ci_r1_rv$info_used_crude,
      es.mcv_or_ci_r2$info_used_crude, es.mcv_or_ci_r2_rv$info_used_crude,
      es.mcv_or_ci_r3$info_used_crude, es.mcv_or_ci_r3_rv$info_used_crude,
      es.mcv_or_ci_r4$info_used_crude, es.mcv_or_ci_r4_rv$info_used_crude,
      es.mcv_or_ci_z1$info_used_crude, es.mcv_or_ci_z1_rv$info_used_crude,
      es.mcv_or_ci_z2$info_used_crude, es.mcv_or_ci_z2_rv$info_used_crude,
      es.mcv_or_ci_z3$info_used_crude, es.mcv_or_ci_z3_rv$info_used_crude,
      es.mcv_or_ci_z4$info_used_crude, es.mcv_or_ci_z4_rv$info_used_crude,

      es.mcv_or_ci_rr2$info_used_crude, es.mcv_or_ci_rr2_rv$info_used_crude,
      es.mcv_or_ci_rr4$info_used_crude, es.mcv_or_ci_rr4_rv$info_used_crude,
      es.mcv_or_ci_rr5$info_used_crude, es.mcv_or_ci_rr5_rv$info_used_crude,
      es.mcv_or_ci_rr6$info_used_crude, es.mcv_or_ci_rr6_rv$info_used_crude,
      es.mcv_or_ci_rr7$info_used_crude, es.mcv_or_ci_rr7_rv$info_used_crude) == "or_ci"))

  expect_equal(es.mcv_or_d$es_crude, -es.mcv_or_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_d$es_crude, -es.mcv_or_se_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_d$es_crude, -es.mcv_or_ci_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_d$se_crude, es.mcv_or_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_d$se_crude, es.mcv_or_se_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_d$se_crude, es.mcv_or_ci_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_g$es_crude, -es.mcv_or_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_g$es_crude, -es.mcv_or_se_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_g$es_crude, -es.mcv_or_ci_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_g$se_crude, es.mcv_or_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_g$se_crude, es.mcv_or_se_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_g$se_crude, es.mcv_or_ci_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_r1$es_crude, -es.mcv_or_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_r1$es_crude, -es.mcv_or_se_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_r1$es_crude, -es.mcv_or_ci_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_r1$se_crude, es.mcv_or_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_r1$se_crude, es.mcv_or_se_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_r1$se_crude, es.mcv_or_ci_r1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_r2$es_crude, -es.mcv_or_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_r2$es_crude, -es.mcv_or_se_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_r2$es_crude, -es.mcv_or_ci_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_r2$se_crude, es.mcv_or_r2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_r2$se_crude, es.mcv_or_se_r2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_r2$se_crude, es.mcv_or_ci_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_r3$es_crude, -es.mcv_or_r3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_r3$es_crude, -es.mcv_or_se_r3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_r3$es_crude, -es.mcv_or_ci_r3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_r3$se_crude, es.mcv_or_r3_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_r3$se_crude, es.mcv_or_se_r3_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_r3$se_crude, es.mcv_or_ci_r3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_r4$es_crude, -es.mcv_or_r4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_r4$es_crude, -es.mcv_or_se_r4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_r4$es_crude, -es.mcv_or_ci_r4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_r4$se_crude, es.mcv_or_r4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_r4$se_crude, es.mcv_or_se_r4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_r4$se_crude, es.mcv_or_ci_r4_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_z1$es_crude, -es.mcv_or_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_z1$es_crude, -es.mcv_or_se_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_z1$es_crude, -es.mcv_or_ci_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_z1$se_crude, es.mcv_or_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_z1$se_crude, es.mcv_or_se_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_z1$se_crude, es.mcv_or_ci_z1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_z2$es_crude, -es.mcv_or_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_z2$es_crude, -es.mcv_or_se_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_z2$es_crude, -es.mcv_or_ci_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_z2$se_crude, es.mcv_or_z2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_z2$se_crude, es.mcv_or_se_z2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_z2$se_crude, es.mcv_or_ci_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_z3$es_crude, -es.mcv_or_z3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_z3$es_crude, -es.mcv_or_se_z3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_z3$es_crude, -es.mcv_or_ci_z3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_z3$se_crude, es.mcv_or_z3_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_z3$se_crude, es.mcv_or_se_z3_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_z3$se_crude, es.mcv_or_ci_z3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_se_z4$es_crude, -es.mcv_or_se_z4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_z4$es_crude, -es.mcv_or_ci_z4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_z4$se_crude, es.mcv_or_z4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_z4$se_crude, es.mcv_or_se_z4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_z4$se_crude, es.mcv_or_ci_z4_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_se_or$es_crude, -es.mcv_or_se_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_or$es_crude, -es.mcv_or_ci_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_or$se_crude, es.mcv_or_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_or$se_crude, es.mcv_or_se_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_or$se_crude, es.mcv_or_ci_or_rv$se_crude, tolerance = 1e-10)


  expect_equal(es.mcv_or_se_rr2$es_crude, -es.mcv_or_se_rr2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr2$es_crude, -es.mcv_or_rr2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr2$es_crude, -es.mcv_or_ci_rr2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr2$se_crude, es.mcv_or_rr2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_rr2$se_crude, es.mcv_or_se_rr2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr2$se_crude, es.mcv_or_ci_rr2_rv$se_crude, tolerance = 1e-10)


  expect_equal(es.mcv_or_se_rr4$es_crude, -es.mcv_or_se_rr4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr4$es_crude, -es.mcv_or_rr4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr4$es_crude, -es.mcv_or_ci_rr4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr4$se_crude, es.mcv_or_rr4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_rr4$se_crude, es.mcv_or_se_rr4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr4$se_crude, es.mcv_or_ci_rr4_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_se_rr5$es_crude, -es.mcv_or_se_rr5_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr5$es_crude, -es.mcv_or_rr5_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr5$es_crude, -es.mcv_or_ci_rr5_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr5$se_crude, es.mcv_or_rr5_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_rr5$se_crude, es.mcv_or_se_rr5_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr5$se_crude, es.mcv_or_ci_rr5_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_se_rr6$es_crude, -es.mcv_or_se_rr6_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr6$es_crude, -es.mcv_or_rr6_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr6$es_crude, -es.mcv_or_ci_rr6_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr6$se_crude, es.mcv_or_rr6_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_rr6$se_crude, es.mcv_or_se_rr6_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr6$se_crude, es.mcv_or_ci_rr6_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_se_rr7$es_crude, -es.mcv_or_se_rr7_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr7$es_crude, -es.mcv_or_rr7_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr7$es_crude, -es.mcv_or_ci_rr7_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_rr7$se_crude, es.mcv_or_rr7_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_se_rr7$se_crude, es.mcv_or_se_rr7_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_ci_rr7$se_crude, es.mcv_or_ci_rr7_rv$se_crude, tolerance = 1e-10)


})

test_that("OR - pval Reverse", {
  dat <- metaumbrella::df.OR[1:5, ]
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$baseline_risk = dat$n_cases_nexp / dat$n_nexp
  dat$or <- or$value
  dat$logor_se <- or$se
  dat$or_ci_lo <- exp(log(dat$value) - qnorm(.975) * dat$logor_se)
  dat$or_ci_up <- exp(log(dat$value) + qnorm(.975) * dat$logor_se)
  z <- log(dat$or) / dat$logor_se
  dat$or_pval <- 1 - 2 * abs(pnorm(z) - 0.5)
  dat = dat[dat$or_pval != 1, ]
  dat$reverse_or <- FALSE
  dat$reverse_or_pval <- FALSE

  es.mcv_or_pval_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "grant", measure = "logrr"), digits = 11)
  es.mcv_or_pval_rr4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "metaumbrella_cases", measure = "logrr"), digits = 11)
  es.mcv_or_pval_rr5 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "metaumbrella_exp", measure = "logrr"), digits = 11)
  es.mcv_or_pval_rr6 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "transpose", measure = "logrr"), digits = 11)
  es.mcv_or_pval_rr7 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "dipietrantonj", measure = "logrr"), digits = 11)
  es.mcv_or_pval_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", measure = "d"), digits = 11)
  es.mcv_or_pval_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", measure = "g"), digits = 11)
  es.mcv_or_pval_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "lipsey_cooper", measure = "r"), digits = 11)
  es.mcv_or_pval_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "lipsey_cooper", measure = "z"), digits = 11)
  es.mcv_or_pval_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "digby", measure = "r"), digits = 11)
  es.mcv_or_pval_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "digby", measure = "z"), digits = 11)
  es.mcv_or_pval_r3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "pearson", measure = "r"), digits = 11)
  es.mcv_or_pval_z3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "pearson", measure = "z"), digits = 11)
  es.mcv_or_pval_r4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "bonett", measure = "r"), digits = 11)
  es.mcv_or_pval_z4 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "bonett", measure = "z"), digits = 11)
  es.mcv_or_pval_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", measure = "logor"), digits = 11)

  dat$reverse_or_pval <- TRUE

  es.mcv_or_pval_rr2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "grant", measure = "logrr"), digits = 11)
  es.mcv_or_pval_rr4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "metaumbrella_cases", measure = "logrr"), digits = 11)
  es.mcv_or_pval_rr5_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "metaumbrella_exp", measure = "logrr"), digits = 11)
  es.mcv_or_pval_rr6_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "transpose", measure = "logrr"), digits = 11)
  es.mcv_or_pval_rr7_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_rr = "dipietrantonj", measure = "logrr"), digits = 11)
  es.mcv_or_pval_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", measure = "d"), digits = 11)
  es.mcv_or_pval_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", measure = "g"), digits = 11)
  es.mcv_or_pval_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "lipsey_cooper", measure = "r"), digits = 11)
  es.mcv_or_pval_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "lipsey_cooper", measure = "z"), digits = 11)
  es.mcv_or_pval_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "digby", measure = "r"), digits = 11)
  es.mcv_or_pval_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "digby", measure = "z"), digits = 11)
  es.mcv_or_pval_r3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "pearson", measure = "r"), digits = 11)
  es.mcv_or_pval_z3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "pearson", measure = "z"), digits = 11)
  es.mcv_or_pval_r4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "bonett", measure = "r"), digits = 11)
  es.mcv_or_pval_z4_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", or_to_cor = "bonett", measure = "z"), digits = 11)
  es.mcv_or_pval_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "or_pval", measure = "logor"), digits = 11)

  expect_true(all(
    c(es.mcv_or_pval_d$info_used_crude, es.mcv_or_pval_d_rv$info_used_crude,
      es.mcv_or_pval_g$info_used_crude, es.mcv_or_pval_g_rv$info_used_crude,
      es.mcv_or_pval_or$info_used_crude, es.mcv_or_pval_or_rv$info_used_crude,
      es.mcv_or_pval_r1$info_used_crude, es.mcv_or_pval_r1_rv$info_used_crude,
      es.mcv_or_pval_r2$info_used_crude, es.mcv_or_pval_r2_rv$info_used_crude,
      es.mcv_or_pval_r3$info_used_crude, es.mcv_or_pval_r3_rv$info_used_crude,
      es.mcv_or_pval_r4$info_used_crude, es.mcv_or_pval_r4_rv$info_used_crude,
      es.mcv_or_pval_z1$info_used_crude, es.mcv_or_pval_z1_rv$info_used_crude,
      es.mcv_or_pval_z2$info_used_crude, es.mcv_or_pval_z2_rv$info_used_crude,
      es.mcv_or_pval_z3$info_used_crude, es.mcv_or_pval_z3_rv$info_used_crude,
      es.mcv_or_pval_z4$info_used_crude, es.mcv_or_pval_z4_rv$info_used_crude,
      es.mcv_or_pval_rr2$info_used_crude, es.mcv_or_pval_rr2_rv$info_used_crude,
      es.mcv_or_pval_rr4$info_used_crude, es.mcv_or_pval_rr4_rv$info_used_crude,
      es.mcv_or_pval_rr5$info_used_crude, es.mcv_or_pval_rr5_rv$info_used_crude,
      es.mcv_or_pval_rr6$info_used_crude, es.mcv_or_pval_rr6_rv$info_used_crude,
      es.mcv_or_pval_rr7$info_used_crude, es.mcv_or_pval_rr7_rv$info_used_crude) == "or_pval"))

  expect_equal(es.mcv_or_pval_d$es_crude, -es.mcv_or_pval_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_g$es_crude, -es.mcv_or_pval_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_r1$es_crude, -es.mcv_or_pval_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_r2$es_crude, -es.mcv_or_pval_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_r3$es_crude, -es.mcv_or_pval_r3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_r4$es_crude, -es.mcv_or_pval_r4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_z1$es_crude, -es.mcv_or_pval_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_z2$es_crude, -es.mcv_or_pval_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_z3$es_crude, -es.mcv_or_pval_z3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_z4$es_crude, -es.mcv_or_pval_z4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_or$es_crude, -es.mcv_or_pval_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_d$se_crude, es.mcv_or_pval_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_g$se_crude, es.mcv_or_pval_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_or$se_crude, es.mcv_or_pval_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_r1$se_crude, es.mcv_or_pval_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_r2$se_crude, es.mcv_or_pval_r2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_r3$se_crude, es.mcv_or_pval_r3_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_r4$se_crude, es.mcv_or_pval_r4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_z1$se_crude, es.mcv_or_pval_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_z2$se_crude, es.mcv_or_pval_z2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_z3$se_crude, es.mcv_or_pval_z3_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_z4$se_crude, es.mcv_or_pval_z4_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_pval_rr2$es_crude, -es.mcv_or_pval_rr2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr4$es_crude, -es.mcv_or_pval_rr4_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr5$es_crude, -es.mcv_or_pval_rr5_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr6$es_crude, -es.mcv_or_pval_rr6_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr7$es_crude, -es.mcv_or_pval_rr7_rv$es_crude, tolerance = 1e-10)

  expect_equal(es.mcv_or_pval_rr2$se_crude, es.mcv_or_pval_rr2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr4$se_crude, es.mcv_or_pval_rr4_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr5$se_crude, es.mcv_or_pval_rr5_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr6$se_crude, es.mcv_or_pval_rr6_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_pval_rr7$se_crude, es.mcv_or_pval_rr7_rv$se_crude, tolerance = 1e-10)

  })
