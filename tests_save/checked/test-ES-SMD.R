### SMD to G/OR/R ------------
test_that("d => G", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))

  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$cohen_d <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$d)
  dat$reverse_d <- FALSE
  es.mfr <- metafor::escalc(
    m1i = dat$mean_exp, sd1i = dat$mean_sd_exp, n1i = dat$n_exp,
    m2i = dat$mean_nexp, sd2i = dat$mean_sd_nexp, n2i = dat$n_nexp,
    measure = "SMD", vtype = "LS2"
  )

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "d"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "g"), digits = 11)

  # SE of D
  comp_es <- with(dat, compute.es::des(d = cohen_d, n.1 = n_exp, n.2 = n_nexp, dig = 12, verbose = FALSE))
  expect_equal(es.mcv_d$es_crude, comp_es$d, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude^2, comp_es$var.d, tolerance = 1e-10)

  # G + var against METAFOR
  expect_equal(unique(es.mcv_g$info_used_crude), "cohen_d")
  expect_equal(es.mcv_g$es_crude, as.numeric(as.character(es.mfr$yi)), tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude^2, as.numeric(as.character(es.mfr$vi)), tolerance = 1e-10)
})
test_that("d => OR", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))

  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$cohen_d <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$d)
  dat$reverse_d <- FALSE
  es.mfr <- metafor::escalc(
    m1i = dat$mean_exp, sd1i = dat$mean_sd_exp, n1i = dat$n_exp,
    m2i = dat$mean_nexp, sd2i = dat$mean_sd_nexp, n2i = dat$n_nexp,
    measure = "SMD", vtype = "LS2"
  )

  comp_es <- with(dat, compute.es::des(d = cohen_d, n.1 = n_exp, n.2 = n_nexp, dig = 12, verbose = FALSE))

  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "logor"), digits = 11)
  mfr_or <- metafor::escalc(
    di = cohen_d, n1i = n_exp, n2i = n_nexp,
    measure = "D2ORL", data = dat
  )


  # D to OR
  expect_equal(unique(es.mcv_or$info_used_crude), "cohen_d")
  expect_equal(es.mcv_or$es_crude, comp_es$lOR, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_ci_lo_crude, comp_es$l.lor, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_ci_up_crude, comp_es$u.lor, tolerance = 1e-10)

  expect_equal(es.mcv_or$es_crude, as.numeric(mfr_or$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude^2, as.numeric(mfr_or$vi), tolerance = 1e-10)
})

test_that("d => R/Z - LIPSEY-COOPER", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))

  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$cohen_d <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$d)
  dat$reverse_d <- FALSE
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "d"), digits = 11)


  es.mcv_r_co <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "r",
    smd_to_cor = "lipsey_cooper"
  ), digits = 12)
  es.mcv_z_co <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "z",
    smd_to_cor = "lipsey_cooper"
  ), digits = 12)

  # esc >= 0.5.1 on R >= 4.3 rejects vector inputs, so convert row by row
  esc_r <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    data.frame(esc::convert_d2r(
      d = dat$cohen_d[i], se = es.mcv_d$se_crude[i],
      grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i]
    ))
  }))
  comp_es <- with(dat, compute.es::des(
    d = cohen_d, n.1 = n_exp, n.2 = n_nexp,
    dig = 12, verbose = FALSE
  ))

  # D to R
  expect_equal(unique(es.mcv_r_co$info_used_crude), "cohen_d")
  expect_equal(es.mcv_r_co$es_crude, comp_es$r, tolerance = 1e-10)
  expect_equal(es.mcv_r_co$se_crude^2, comp_es$var.r, tolerance = 1e-10)

  expect_equal(unique(es.mcv_z_co$info_used_crude), "cohen_d")
  expect_equal(es.mcv_z_co$es_crude, esc_r$fishers.z, tolerance = 1e-10)
  # esc::convert_d2r()'s z-SE (esc_r$se) is NOT the Fisher transform of its own r
  # variance; metaConvert intentionally uses the Fisher-consistent delta value
  # se(z) = se(r)/(1 - r^2). The z POINT estimate above still matches esc.
  expect_equal(es.mcv_z_co$se_crude,
               es.mcv_r_co$se_crude / (1 - es.mcv_r_co$es_crude^2),
               tolerance = 1e-9)
})

test_that("d => R/Z - VIECHTBAUER", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))

  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$cohen_d <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$d)
  dat$reverse_d <- FALSE

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "d"), digits = 11)
  es.mcv_rb <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "r",
    smd_to_cor = "viechtbauer"
  ), digits = 11)
  es.mcv_zb <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "z",
    smd_to_cor = "viechtbauer"
  ), digits = 11)

  mfr_rb <- metafor::escalc(
    di = cohen_d, n1i = n_exp, n2i = n_nexp,
    measure = "RBIS", data = dat
  )
  mfr_zb <- suppressWarnings(
    metafor::escalc(
      di = cohen_d, n1i = n_exp, n2i = n_nexp,
      measure = "ZBIS", data = dat
    )
  )

  zbis <- function(rb, n1, n2) {
    n <- n1 + n2
    p <- n1 / n
    fzp <- dnorm(qnorm(p))
    a <- sqrt(fzp) / (p*(1-p))^(1/4)
    rb = ifelse(rb > 1, 1,
                 ifelse(rb < -1, -1, rb))
    zrb <- (a/2) * log((1+a*rb)/(1-a*rb))


    ci_1 <- zrb - qnorm(.975) * sqrt(1/(n-1))
    ci_2 <- zrb + qnorm(.975) * sqrt(1/(n-1))
    cbind(zrb,
          ci_1,
          ci_2)
  }
  res_z = zbis(as.numeric(mfr_rb$yi), dat$n_exp, dat$n_nexp)


  expect_equal(unique(es.mcv_rb$info_used_crude), "cohen_d")
  expect_equal(es.mcv_rb$es_crude, as.numeric(mfr_rb$yi), tolerance = 1e-10)
  expect_equal(es.mcv_rb$se_crude^2, as.numeric(mfr_rb$vi), tolerance = 1e-10)

  expect_equal(unique(es.mcv_zb$info_used_crude), "cohen_d")
  expect_equal(es.mcv_zb$es_crude, res_z[,1], tolerance = 1e-10)
  expect_equal(es.mcv_zb$es_ci_lo_crude, res_z[,2], tolerance = 1e-10)
  expect_equal(es.mcv_zb$es_ci_up_crude, res_z[,3], tolerance = 1e-10)
})

test_that("g => g + SE", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$hedges_g <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$g)
  dat$cohen_d <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$d)
  dat$reverse_g <- FALSE
  dat$reverse_d <- FALSE
  es.mfr <- metafor::escalc(
    m1i = dat$mean_exp, sd1i = dat$mean_sd_exp, n1i = dat$n_exp,
    m2i = dat$mean_nexp, sd2i = dat$mean_sd_nexp, n2i = dat$n_nexp,
    measure = "SMD", vtype = "LS2"
  )

  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "g"), digits = 11)
  es.mcv_or_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "logor"), digits = 11)
  es.mcv_r_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "r"), digits = 11)
  es.mcv_z_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "z"), digits = 11)
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "g"), digits = 11)
  es.mcv_or_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "logor"), digits = 11)
  es.mcv_r_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "r"), digits = 11)
  es.mcv_z_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "z"), digits = 11)

  expect_equal(unique(es.mcv_g$info_used_crude), "hedges_g")
  expect_equal(unique(es.mcv_d$info_used_crude), "cohen_d")
  expect_equal(es.mcv_g$es_crude, es.mcv_d$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_d$se_crude, tolerance = 1e-10)

  expect_equal(as.numeric(es.mfr$yi), es.mcv_d$es_crude, tolerance = 1e-10)
  expect_equal(as.numeric(es.mfr$vi), es.mcv_d$se_crude^2, tolerance = 1e-10)

  # G to OR
  expect_equal(unique(es.mcv_or_d$info_used_crude), "cohen_d")
  expect_equal(unique(es.mcv_or_g$info_used_crude), "hedges_g")
  expect_equal(es.mcv_or_d$es_crude, es.mcv_or_g$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or_d$se_crude, es.mcv_or_g$se_crude, tolerance = 1e-10)

  expect_equal(unique(es.mcv_r_d$info_used_crude), "cohen_d")
  expect_equal(unique(es.mcv_r_g$info_used_crude), "hedges_g")
  expect_equal(es.mcv_r_d$es_crude, es.mcv_r_g$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r_d$se_crude, es.mcv_r_g$se_crude, tolerance = 1e-10)

  expect_equal(unique(es.mcv_z_d$info_used_crude), "cohen_d")
  expect_equal(unique(es.mcv_z_g$info_used_crude), "hedges_g")
  expect_equal(es.mcv_z_d$es_crude, es.mcv_z_g$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z_d$se_crude, es.mcv_z_g$se_crude, tolerance = 1e-10)
})

test_that("es_from_cohen_d - REVERSE", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$hedges_g <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$g)
  dat$cohen_d <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$d)
  dat$reverse_g <- FALSE
  dat$reverse_d <- FALSE

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "d"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "g"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "logor"), digits = 11)
  es.mcv_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "r",
                                 smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "r",
                                 smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "z",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "z",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)

  dat$reverse_d <- TRUE
  es.mcv_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "d"), digits = 11)
  es.mcv_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "g"), digits = 11)
  es.mcv_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "logor"), digits = 11)
  es.mcv_r1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "r",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "r",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_z1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "z",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d", measure = "z",
                                    smd_to_cor = "lipsey_cooper"), digits = 11)

  expect_equal(es.mcv_d$es_crude, -es.mcv_d_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$es_crude, -es.mcv_g_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, -es.mcv_or_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$es_crude, -es.mcv_z1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$se_crude, es.mcv_z1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$es_crude, -es.mcv_z2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$se_crude, es.mcv_z2_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$es_crude, -es.mcv_r1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$se_crude, es.mcv_r1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$es_crude, -es.mcv_r2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$se_crude, es.mcv_r2_r$se_crude, tolerance = 1e-10)
})


test_that("es_from_hedges_g - REVERSE", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$hedges_g <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$g)
  dat$cohen_d <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$d)
  dat$reverse_g <- FALSE
  dat$reverse_d <- FALSE

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "d"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "g"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "logor"), digits = 11)
  es.mcv_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "r",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "r",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "z",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "z",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)

  dat$reverse_g <- TRUE
  es.mcv_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "d"), digits = 11)
  es.mcv_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "g"), digits = 11)
  es.mcv_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "logor"), digits = 11)
  es.mcv_r1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "r",
                                    smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "r",
                                    smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_z1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "z",
                                    smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "hedges_g", measure = "z",
                                    smd_to_cor = "lipsey_cooper"), digits = 11)


  expect_equal(es.mcv_d$es_crude, -es.mcv_d_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$es_crude, -es.mcv_g_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, -es.mcv_or_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$es_crude, -es.mcv_z1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$se_crude, es.mcv_z1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$es_crude, -es.mcv_z2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$se_crude, es.mcv_z2_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$es_crude, -es.mcv_r1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$se_crude, es.mcv_r1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$es_crude, -es.mcv_r2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$se_crude, es.mcv_r2_r$se_crude, tolerance = 1e-10)
})
#
test_that("es_from_cohen_d_adj - REVERSE", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$hedges_g <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$g)
  dat$cohen_d <- with(dat, es_from_means_sd(
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp, n_exp = n_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp, n_nexp = n_nexp
  )$d)
  dat$cohen_d_adj <- dat$cohen_d
  dat$cov_outcome_r <- 0.3
  dat$n_cov_ancova <- 4
  dat$reverse_d <- FALSE
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "d"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "g"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "logor"), digits = 11)
  es.mcv_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "r",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "r",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "z",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "z",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)

  dat$reverse_d <- TRUE
  es.mcv_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "d"), digits = 11)
  es.mcv_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "g"), digits = 11)
  es.mcv_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "logor"), digits = 11)
  es.mcv_r1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "r",
                                    smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "r",
                                    smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_z1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "z",
                                    smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "cohen_d_adj", measure = "z",
                                    smd_to_cor = "lipsey_cooper"), digits = 11)


  expect_equal(unique(es.mcv_d$info_used_adjusted), "cohen_d_adj")
  expect_equal(es.mcv_d$es_adjusted, -es.mcv_d_r$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_adjusted, es.mcv_d_r$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_g$es_adjusted, -es.mcv_g_r$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_adjusted, es.mcv_g_r$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_adjusted, -es.mcv_or_r$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_adjusted, es.mcv_or_r$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_z1$es_crude, -es.mcv_z1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$se_crude, es.mcv_z1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$es_crude, -es.mcv_z2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$se_crude, es.mcv_z2_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$es_crude, -es.mcv_r1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$se_crude, es.mcv_r1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$es_crude, -es.mcv_r2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$se_crude, es.mcv_r2_r$se_crude, tolerance = 1e-10)
})
