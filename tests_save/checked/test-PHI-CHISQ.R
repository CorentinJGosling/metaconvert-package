# test_that("Internal X2-pval RR", {
#   tab <- matrix(c(50,20,50,80), nrow=2)
#   phi <- sqrt(chisq.test(tab, correct=FALSE)$statistic[[1]] / sum(tab))
#
#   res_na = es_from_phi(phi=phi, n_sample=sum(tab))
#
#   res_mcv_phi = es_from_phi(phi=phi, n_sample=sum(tab),
#                             n_cases=sum(tab[,1]), n_exp=sum(tab[1,]))
#   res_mcv_or = es_from_2x2(n_cases_exp = tab[1,1], n_cases_nexp = tab[2,1],
#                            n_controls_exp = tab[1,2], n_controls_nexp = tab[2,2])
#
#   res_mfr_phi = summary(escalc(measure="RTET", ai=tab[1,1], bi=tab[1,2],
#                            ci=tab[2,1], di=tab[2,2],
#                            digits=12))
#
#
#   if (is.na(res_na$r)) message("1. pass correct NA when missing cases/exp")
#   if (res_mcv_phi$r == res_mcv_or$r && res_mcv_phi$r_se == res_mcv_or$r_se) {
#     message("2. pass internal values from phi and 2x2 tables")
#   }
#   if (res_mcv_phi$r == res_mfr_phi$yi && res_mcv_phi$r_se^2 == res_mfr_phi$vi) {
#     message("3. pass chisq metafor comparison")
#   }
#
#
#   chi <- chisq.test(tab, correct=FALSE)$statistic
#
#   res_na = es_from_chisq(chisq=chi, n_sample=sum(tab))
#   if (is.na(res_na$r)) message("4. pass correct NA when missing cases/exp")
#
#   res_mcv_chisq = es_from_chisq(chisq=chi, n_sample=sum(tab),
#                               n_cases=sum(tab[,1]), n_exp=sum(tab[1,]))
#
#
#   if (res_mcv_chisq$r == res_mcv_or$r && res_mcv_chisq$r_se == res_mcv_or$r_se) {
#     message("5. pass internal values from chisq and 2x2 tables")
#   }
#   if (res_mcv_chisq$r == res_mfr_phi$yi && res_mcv_chisq$r_se^2 == res_mfr_phi$vi) {
#     message("6. pass phi metafor comparison")
#   }
# })

# PHI TO OR/RR -----
test_that("phi=>OR", {
  skip_if_not_installed("mvtnorm")
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:150, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))

  cont <- suppressWarnings(metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  ))
  dat <- dat[which(!is.na(cont$ai)), ]
  cont <- metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  )
  or_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi, ci = cont$ci, di = cont$di,
    measure = "OR",
    digits = 12
  )

  es.mcv_or1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "logor"),
                        digits = 12
  )

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "d"),
                      digits = 12
  )
  es.mcv_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "r"),
                      digits = 12
  )
  # esc >= 0.5.1 on R >= 4.3 rejects vector inputs, so call row by row
  esc_phi_row <- function(es.type) {
    do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
      data.frame(esc::esc_phi(phi = dat$phi[i], totaln = dat$n_sample[i], es.type = es.type))
    }))
  }
  comp_res_d <- esc_phi_row("d")
  comp_res_r <- esc_phi_row("r")


  ## test ES
  expect_equal(unique(es.mcv_d$info_used_crude), "phi")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 5e-1)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 5e-1)

  expect_equal(unique(es.mcv_r$info_used_crude), "phi")
  expect_equal(es.mcv_r$es_crude, comp_res_r$es, tolerance = 4e-1)
  expect_equal(es.mcv_r$se_crude, comp_res_r$se, tolerance = 3e-1)

  expect_equal(unique(es.mcv_or1$info_used_crude), "phi")
  expect_equal(es.mcv_or1$es_crude, as.numeric(or_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or1$se_crude^2, or_mfr$vi, tolerance = 1e-10)
})
test_that("phi=>RR", {
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:150, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))

  cont <- suppressWarnings(metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  ))
  dat <- dat[which(!is.na(cont$ai)), ]
  cont <- metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  )
  or_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi, ci = cont$ci, di = cont$di,
    measure = "RR",
    digits = 12
  )

  es.mcv_or1 <- summary(convert_df(
    dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi",
    measure = "logrr"),
    digits = 12
  )

  expect_equal(unique(es.mcv_or1$info_used_crude), "phi")
  expect_equal(es.mcv_or1$es_crude, as.numeric(or_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or1$se_crude^2, or_mfr$vi, tolerance = 1e-10)
})
# PHI TO SMD/R-Z -----
test_that("phi=>SMD/R-Z", {
  skip_if_not_installed("mvtnorm")
  set.seed(1234)

  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:150, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))

  cont <- suppressWarnings(metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  ))
  dat <- dat[which(!is.na(cont$ai)), ]
  cont <- metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  )
  or_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi, ci = cont$ci, di = cont$di,
    measure = "OR",
    digits = 12
  )
  es.mcv_or <- summary(convert_df(
    dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi",
    measure = "logor"),
    digits = 12
  )
  expect_equal(unique(es.mcv_or$info_used_crude), "phi")
  expect_equal(es.mcv_or$es_crude, as.numeric(or_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude^2, or_mfr$vi, tolerance = 1e-10)

  d_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi, ci = cont$ci, di = cont$di,
    measure = "OR2DL",
    digits = 12
  )
  es.mcv_d <- summary(convert_df(
    dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi",
    measure = "d"),
    digits = 12
  )

  expect_equal(unique(es.mcv_d$info_used_crude), "phi")
  expect_equal(es.mcv_d$es_crude, as.numeric(d_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude^2, d_mfr$vi, tolerance = 1e-10)

  r_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi, ci = cont$ci, di = cont$di,
    measure = "RTET",
    digits = 12
  )
  es.mcv_r <- summary(convert_df(
    dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi",
    measure = "r"),
    digits = 12
  )

  expect_equal(unique(es.mcv_r$info_used_crude), "phi")
  expect_equal(es.mcv_r$es_crude, as.numeric(r_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude^2, r_mfr$vi, tolerance = 1e-10)
})

# CHISQ TO SMD/OR/RR -----

test_that("X2 => OR", {
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:150, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))
  dat$reverse_chisq <- FALSE
  cont <- suppressWarnings(metafor::conv.2x2(
    x2i = dat$chisq, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  ))
  dat <- dat[which(!is.na(cont$ai)), ]

  cont <- metafor::conv.2x2(
    x2i = dat$chisq, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases,
    correct = FALSE
  )
  or_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi,
    ci = cont$ci, di = cont$di,
    measure = "OR",
    digits = 12
  )
  es.mcv_or <- summary(
    convert_df(dat,
      verbose = FALSE,
      es_selected = "hierarchy", hierarchy = "chisq", measure = "logor"
    ),
    digits = 12
  )


  expect_equal(unique(es.mcv_or$info_used_crude), "chisq")
  expect_equal(es.mcv_or$es_crude, as.numeric(or_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude^2, or_mfr$vi, tolerance = 1e-10)
})
test_that("X2 => RR", {
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:150, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))
  dat$reverse_chisq <- FALSE
  cont <- suppressWarnings(metafor::conv.2x2(
    x2i = dat$chisq, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  ))
  dat <- dat[which(!is.na(cont$ai)), ]

  cont <- metafor::conv.2x2(
    x2i = dat$chisq, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases,
    correct = FALSE
  )
  or_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi,
    ci = cont$ci, di = cont$di,
    measure = "RR",
    digits = 12
  )
  es.mcv_or <- summary(
    convert_df(dat,
               verbose = FALSE,
               es_selected = "hierarchy", hierarchy = "chisq", measure = "logrr"
    ),
    digits = 12
  )

  expect_equal(unique(es.mcv_or$info_used_crude), "chisq")
  expect_equal(es.mcv_or$es_crude, as.numeric(or_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude^2, or_mfr$vi, tolerance = 1e-10)
})

test_that("Internal X2-pval OR", {
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:150, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))
  dat$chisq_pval <- pchisq(dat$chisq, df = 1, lower.tail = FALSE)
  dat$reverse_chisq <- FALSE

  es.mcv_or1 <- summary(
    convert_df(
      dat,
      verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "or"
    ),
    digits = 11
  )
  es.mcv_or2 <- summary(
    convert_df(
      dat,
      verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "or"
    ),
    digits = 11
  )
  row1 <- which(es.mcv_or1$info_used_crude == "chisq_pval")
  row2 <- which(es.mcv_or2$info_used_crude == "chisq")
  expect_equal(unique(es.mcv_or1$info_used_crude[row1]), "chisq_pval")
  expect_equal(unique(es.mcv_or2$info_used_crude[row1]), "chisq")
  expect_equal(es.mcv_or1$es_crude[row1], es.mcv_or2$es_crude[row2], tolerance = 1e-10)
  expect_equal(es.mcv_or1$se_crude[row1], es.mcv_or2$se_crude[row2], tolerance = 1e-10)

})
# chi TO SMD/R-Z -----
test_that("chi=>SMD/R-Z", {
  skip_if_not_installed("mvtnorm")
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:150, ]
  dat$chisq <- runif(nrow(dat), 0.5, 3)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))
  dat$reverse_chisq=FALSE
  cont <- suppressWarnings(metafor::conv.2x2(
    x2i = dat$chisq, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases,
    correct=FALSE
  ))
  dat <- dat[which(!is.na(cont$ai)), ]
  or_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi, ci = cont$ci, di = cont$di,
    measure = "OR",
    digits = 12
  )
  es.mcv_or <- summary(convert_df(
    dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq",
    measure = "logor"),
    digits = 12
  )
  expect_equal(unique(es.mcv_or$info_used_crude), "chisq")
  expect_equal(es.mcv_or$es_crude, as.numeric(or_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude^2, or_mfr$vi, tolerance = 1e-10)

  d_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi, ci = cont$ci, di = cont$di,
    measure = "OR2DL",
    digits = 12
  )
  es.mcv_d <- summary(convert_df(
    dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq",
    measure = "d"),
    digits = 12
  )

  expect_equal(unique(es.mcv_d$info_used_crude), "chisq")
  expect_equal(es.mcv_d$es_crude, as.numeric(d_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude^2, d_mfr$vi, tolerance = 1e-10)


  r_mfr <- metafor::escalc(
    ai = cont$ai, bi = cont$bi, ci = cont$ci, di = cont$di,
    measure = "RTET",
    digits = 12
  )
  es.mcv_r <- summary(convert_df(
    dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq",
    measure = "r"),
    digits = 12
  )

  expect_equal(unique(es.mcv_r$info_used_crude), "chisq")
  expect_equal(es.mcv_r$es_crude, as.numeric(r_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude^2, r_mfr$vi, tolerance = 1e-10)
})
test_that("Internal X2-pval RR", {
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:150, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))
  dat$chisq_pval <- pchisq(dat$chisq, df = 1, lower.tail = FALSE)
  dat$reverse_chisq <- FALSE

  es.mcv_or1 <- summary(
    convert_df(
      dat,
      verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "rr"
    ),
    digits = 11
  )
  es.mcv_or2 <- summary(
    convert_df(
      dat,
      verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "rr"
    ),
    digits = 11
  )
  row1 <- which(es.mcv_or1$info_used_crude == "chisq_pval")
  row2 <- which(es.mcv_or2$info_used_crude == "chisq")
  expect_equal(unique(es.mcv_or1$info_used_crude[row1]), "chisq_pval")
  expect_equal(unique(es.mcv_or2$info_used_crude[row1]), "chisq")
  expect_equal(es.mcv_or1$es_crude[row1], es.mcv_or2$es_crude[row2], tolerance = 1e-10)
  expect_equal(es.mcv_or1$se_crude[row1], es.mcv_or2$se_crude[row2], tolerance = 1e-10)

  # View(cbind(es.mcv_or1$es_crude[row1], es.mcv_or2$es_crude[row2], dat$chisq_pval[row1]))
})

# REVERSE -----

test_that("reverse-phi", {
  skip_if_not_installed("mvtnorm")
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:50, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))
  dat$chisq_pval <- pchisq(dat$chisq, df = 1, lower.tail = FALSE)

  cont <- suppressWarnings(metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  ))
  dat <- dat[which(!is.na(cont$ai)), ]


  dat$reverse_phi <- FALSE

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "logor"), digits = 11)
  es.mcv_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "z"), digits = 11)
  es.mcv_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "r"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "g"), digits = 11)

  dat$reverse_phi <- TRUE

  es.mcv_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "d"), digits = 11)
  es.mcv_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "logor"), digits = 11)
  es.mcv_z_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "z"), digits = 11)
  es.mcv_r_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "r"), digits = 11)
  es.mcv_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "phi", measure = "g"), digits = 11)

  expect_true(all(c(es.mcv_d$info_used_crude, es.mcv_or$info_used_crude,
                    es.mcv_g$info_used_crude,
                    es.mcv_z$info_used_crude, es.mcv_r$info_used_crude,
                    es.mcv_d_r$info_used_crude, es.mcv_or_r$info_used_crude,
                    es.mcv_g_r$info_used_crude,
                    es.mcv_z_r$info_used_crude, es.mcv_r_r$info_used_crude) == "phi"))

  expect_equal(es.mcv_d$es_crude, -es.mcv_d_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, -es.mcv_or_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z$es_crude, -es.mcv_z_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z$se_crude, es.mcv_z_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r$es_crude, -es.mcv_r_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, es.mcv_r_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$es_crude, -es.mcv_g_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_r$se_crude, tolerance = 1e-10)
})


test_that("reverse-chisq", {
  skip_if_not_installed("mvtnorm")
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:50, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))
  dat$chisq_pval <- pchisq(dat$chisq, df = 1, lower.tail = FALSE)
  cont <- suppressWarnings(metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  ))
  dat <- dat[which(!is.na(cont$ai)), ]

  dat$reverse_chisq <- FALSE
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "logor"), digits = 11)
  es.mcv_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "z"), digits = 11)
  es.mcv_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "r"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "g"), digits = 11)

  dat$reverse_chisq <- TRUE

  es.mcv_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "d"), digits = 11)
  es.mcv_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "logor"), digits = 11)
  es.mcv_z1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "z"), digits = 11)
  es.mcv_r1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "r"), digits = 11)
  es.mcv_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq", measure = "g"), digits = 11)

  expect_true(all(c(es.mcv_d$info_used_crude, es.mcv_or$info_used_crude,
                    es.mcv_g$info_used_crude,
                    es.mcv_z1$info_used_crude, es.mcv_r1$info_used_crude,
                    es.mcv_d_r$info_used_crude, es.mcv_or_r$info_used_crude,
                    es.mcv_g_r$info_used_crude,
                    es.mcv_z1_r$info_used_crude, es.mcv_r1_r$info_used_crude) == "chisq"))

  expect_equal(es.mcv_d$es_crude, -es.mcv_d_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, -es.mcv_or_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$es_crude, -es.mcv_z1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$se_crude, es.mcv_z1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$es_crude, -es.mcv_r1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$se_crude, es.mcv_r1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$es_crude, -es.mcv_g_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_r$se_crude, tolerance = 1e-10)
})


test_that("reverse-chisq chisq_pval", {
  skip_if_not_installed("mvtnorm")
  set.seed(1234)
  # drop the original df.haza counts: the synthetic phi/chisq + n_cases/n_sample
  # built below would contradict them and trip the 2x2 consistency checks
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp),
                select = -c(n_cases_exp, n_cases_nexp, n_controls_exp,
                            n_controls_nexp, n_cases, n_controls, n_sample))
  dat <- dat[1:50, ]
  dat$phi <- abs(rnorm(nrow(dat), 0, 0.35))
  dat$phi[abs(dat$phi) > 0.5] <- 0.5
  dat$chisq <- runif(nrow(dat), 0.1, 5)
  dat$n_exp = dat$n_exp + 50
  dat$n_nexp = dat$n_nexp + 50
  dat$n_sample <- dat$n_exp + dat$n_nexp
  dat$n_cases <- dat$n_exp + round(runif(nrow(dat), 15,30))
  dat$chisq_pval <- pchisq(dat$chisq, df = 1, lower.tail = FALSE)
  cont <- suppressWarnings(metafor::conv.2x2(
    ri = dat$phi, ni = dat$n_sample,
    n1i = dat$n_exp, n2i = dat$n_cases
  ))
  dat <- dat[which(!is.na(cont$ai)), ]

  dat$reverse_chisq_pval <- FALSE
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "logor"), digits = 11)
  es.mcv_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "z"), digits = 11)
  es.mcv_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "r"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "g"), digits = 11)

  dat$reverse_chisq_pval <- TRUE

  es.mcv_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "d"), digits = 11)
  es.mcv_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "logor"), digits = 11)
  es.mcv_z1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "z"), digits = 11)
  es.mcv_r1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "r"), digits = 11)
  es.mcv_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "chisq_pval", measure = "g"), digits = 11)

  expect_true(all(c(es.mcv_d$info_used_crude, es.mcv_or$info_used_crude,
                    es.mcv_g$info_used_crude,
                    es.mcv_z1$info_used_crude, es.mcv_r1$info_used_crude,
                    es.mcv_d_r$info_used_crude, es.mcv_or_r$info_used_crude, es.mcv_g_r$info_used_crude,
                    es.mcv_z1_r$info_used_crude, es.mcv_r1_r$info_used_crude) == "chisq_pval"))

  expect_equal(es.mcv_d$es_crude, -es.mcv_d_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, -es.mcv_or_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$es_crude, -es.mcv_z1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$se_crude, es.mcv_z1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$es_crude, -es.mcv_r1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$se_crude, es.mcv_r1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$es_crude, -es.mcv_g_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_r$se_crude, tolerance = 1e-10)
})
