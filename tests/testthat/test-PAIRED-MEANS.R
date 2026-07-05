library(testthat)
library(metaConvert)

### SMD from PRE POST MEANS/SD-----
test_that("D - Means/SD - bonett", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE

  smc_exp <- metafor::escalc(
    m1i = res$mean_pre_exp,
    sd1i = res$mean_pre_sd_exp,
    m2i = res$mean_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp,
    measure = "SMCRH"
  )
  smc_nexp <- metafor::escalc(
    m1i = res$mean_pre_nexp,
    sd1i = res$mean_pre_sd_nexp,
    m2i = res$mean_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp,
    measure = "SMCRH"
  )

  J_exp = metaConvert:::.d_j(res$n_exp - 1)
  J_nexp = metaConvert:::.d_j(res$n_nexp - 1)


  smcc_vd <- (-as.numeric(as.character(smc_exp$yi)) / J_exp) -
             (-as.numeric(as.character(smc_nexp$yi)) / J_nexp)
  se_vd <- sqrt(smc_exp$vi / (J_exp^2) + smc_nexp$vi / J_nexp^2)
  ## metaconvert
  es.mcv_d <- summary(convert_df(res,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = "means_sd_pre_post",
    measure = "d",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "bonett"
  ), digits = 11)
  ## d
  expect_equal(unique(es.mcv_d$info_used_crude), "means_sd_pre_post")
  expect_equal(es.mcv_d$es_crude, as.numeric(as.character(smcc_vd)),
    tolerance = 1e-10
  )
  expect_equal(es.mcv_d$se_crude, se_vd, tolerance = 1e-10)
})
test_that("D - Means/SD - cooper", {
  p = 0.4
  n = 98
  d = 0.31323
  dat = data.frame(n_exp=NA)
  dat$n_exp = round(p * n)
  dat$n_nexp = n - dat$n_exp
  scores_grp1 <- rnorm(dat$n_exp, d, 1)
  scores_grp2 <- rnorm(dat$n_nexp, 0, 1)
  dat$mean_exp <- mean(scores_grp1)
  dat$mean_nexp <- mean(scores_grp2)
  dat$mean_sd_exp <- sd(scores_grp1)
  dat$mean_sd_nexp <- sd(scores_grp2)

  scores_grp1_pre <- scores_grp1 - rnorm(dat$n_exp, 0, 1)
  scores_grp2_pre <- scores_grp2 - rnorm(dat$n_nexp, 0, 1)
  dat$mean_pre_exp <- mean(scores_grp1_pre)
  dat$mean_pre_nexp <- mean(scores_grp2_pre)
  dat$mean_pre_sd_exp <- sd(scores_grp1_pre)
  dat$mean_pre_sd_nexp <- sd(scores_grp2_pre)

  dat$r_pre_post_exp <- cor.test(~scores_grp1_pre+scores_grp1)$estimate
  dat$r_pre_post_nexp <- cor.test(~scores_grp2_pre+scores_grp2)$estimate
  dat$reverse_means_pre_post <- FALSE

  res_g1 = TOSTER::smd_calc(
    formula = c(scores_grp1_pre, scores_grp1) ~
      rep(c("A", "B"), each = length(scores_grp1_pre)),
    bias_correction = FALSE, rm_correction = TRUE,
    paired = TRUE, smd_ci = "t")
  res_g2 = TOSTER::smd_calc(
    formula = c(scores_grp2_pre, scores_grp2) ~
      rep(c("A", "B"), each = length(scores_grp2_pre)),
    bias_correction = FALSE, rm_correction = TRUE,
    paired = TRUE, smd_ci = "t")


  es.mcv_d <- summary(convert_df(dat,
                                 verbose = FALSE,
                                 es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
                                 measure = "d",
                                 smd_to_cor = "viechtbauer",
                                 pre_post_to_smd = "cooper"
  ), digits = 11)
  ## d

  es.mcv_d$es_crude
  expect_equal(unique(es.mcv_d$info_used_crude), "means_sd_pre_post")
  expect_equal(es.mcv_d$es_crude, (-res_g1$estimate) - (-res_g2$estimate),
               tolerance = 1e-6)
  expect_equal(es.mcv_d$se_crude, sqrt(res_g1$SE^2 + res_g2$SE^2),
               tolerance = 5e-1)
})


test_that("G - Means/SD - bonett", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE

  smc_exp <- metafor::escalc(
    m1i = res$mean_pre_exp,
    sd1i = res$mean_pre_sd_exp,
    m2i = res$mean_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp,
    measure = "SMCRH"
  )
  smc_nexp <- metafor::escalc(
    m1i = res$mean_pre_nexp,
    sd1i = res$mean_pre_sd_nexp,
    m2i = res$mean_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp,
    measure = "SMCRH"
  )

  smcc_vg <- (-as.numeric(as.character(smc_exp$yi))) - (-as.numeric(as.character(smc_nexp$yi)))
  se_vg <- sqrt(smc_exp$vi + smc_nexp$vi)
  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "g",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "bonett"
  ), digits = 11)
  ## d
  expect_equal(unique(es.mcv_g$info_used_crude), "means_sd_pre_post")
  expect_equal(es.mcv_g$es_crude, as.numeric(as.character(smcc_vg)),
    tolerance = 1e-10
  )
  expect_equal(es.mcv_g$se_crude, se_vg, tolerance = 1e-10)
})
test_that("G - Means/SD - cooper", {
  p = 0.4
  n = 98
  d = 0.5
  dat = data.frame(n_exp=NA)
  dat$n_exp = round(p * n)
  dat$n_nexp = n - dat$n_exp
  scores_grp1 <- rnorm(dat$n_exp, d, 1)
  scores_grp2 <- rnorm(dat$n_nexp, 0, 1)
  dat$mean_exp <- mean(scores_grp1)
  dat$mean_nexp <- mean(scores_grp2)
  dat$mean_sd_exp <- sd(scores_grp1)
  dat$mean_sd_nexp <- sd(scores_grp2)

  scores_grp1_pre <- scores_grp1 - rnorm(dat$n_exp, 0, 1)
  scores_grp2_pre <- scores_grp2 - rnorm(dat$n_nexp, 0, 1)
  dat$mean_pre_exp <- mean(scores_grp1_pre)
  dat$mean_pre_nexp <- mean(scores_grp2_pre)
  dat$mean_pre_sd_exp <- sd(scores_grp1_pre)
  dat$mean_pre_sd_nexp <- sd(scores_grp2_pre)

  dat$r_pre_post_exp <- cor.test(~scores_grp1_pre+scores_grp1)$estimate
  dat$r_pre_post_nexp <- cor.test(~scores_grp2_pre+scores_grp2)$estimate
  dat$reverse_means_pre_post <- FALSE

  res_g1 = TOSTER::smd_calc(
    formula = c(scores_grp1_pre, scores_grp1) ~
      rep(c("A", "B"), each = length(scores_grp1_pre)),
    bias_correction = TRUE, rm_correction = TRUE,
    paired = TRUE, smd_ci = "t")
  res_g2 = TOSTER::smd_calc(
    formula = c(scores_grp2_pre, scores_grp2) ~
      rep(c("A", "B"), each = length(scores_grp2_pre)),
    bias_correction = TRUE, rm_correction = TRUE,
    paired = TRUE, smd_ci = "t")


  es.mcv_d <- summary(convert_df(dat,
                                 verbose = FALSE,
                                 es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
                                 measure = "g",
                                 smd_to_cor = "viechtbauer",
                                 pre_post_to_smd = "cooper"
  ), digits = 11)
  ## d

  expect_equal(unique(es.mcv_d$info_used_crude), "means_sd_pre_post")
  expect_equal(es.mcv_d$es_crude, (-res_g1$estimate) - (-res_g2$estimate),
               tolerance = 1e-6)
  expect_equal(es.mcv_d$se_crude, sqrt(res_g1$SE^2 + res_g2$SE^2),
               tolerance = 5e-1)
})

test_that("G - SE - bonett", {
  res <- metaumbrella::df.SMC
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls

  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls

  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_pre_se_nexp <- res$mean_pre_sd_nexp / sqrt(res$n_nexp)

  res$mean_sd_nexp <- res$sd_controls
  res$mean_se_nexp <- res$mean_sd_nexp / sqrt(res$n_nexp)

  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases

  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_pre_se_exp <- res$mean_pre_sd_exp / sqrt(res$n_exp)

  res$mean_sd_exp <- res$sd_cases
  res$mean_se_exp <- res$mean_sd_exp / sqrt(res$n_exp)

  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE

  smc_exp <- metafor::escalc(
    m1i = res$mean_pre_exp,
    sd1i = res$mean_pre_sd_exp,
    m2i = res$mean_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp,
    measure = "SMCRH"
  )
  smc_nexp <- metafor::escalc(
    m1i = res$mean_pre_nexp,
    sd1i = res$mean_pre_sd_nexp,
    m2i = res$mean_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp,
    measure = "SMCRH"
  )

  smcc_vg <- (-as.numeric(as.character(smc_exp$yi))) - (-as.numeric(as.character(smc_nexp$yi)))
  se_vg <- sqrt(smc_exp$vi + smc_nexp$vi)

  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "g",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "bonett"
  ), digits = 11)
  ## d
  expect_equal(unique(es.mcv_g$info_used_crude), "means_se_pre_post")
  expect_equal(es.mcv_g$es_crude, as.numeric(as.character(smcc_vg)),
    tolerance = 1e-10
  )
  expect_equal(es.mcv_g$se_crude, se_vg, tolerance = 1e-10)
})

test_that("G - CI - bonett", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls

  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls

  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_pre_se_nexp <- res$mean_pre_sd_nexp / sqrt(res$n_nexp)
  res$mean_pre_ci_lo_nexp <- res$mean_pre_nexp - qt(.975, res$n_nexp - 1) * res$mean_pre_se_nexp
  res$mean_pre_ci_up_nexp <- res$mean_pre_nexp + qt(.975, res$n_nexp - 1) * res$mean_pre_se_nexp

  res$mean_sd_nexp <- res$sd_controls
  res$mean_se_nexp <- res$mean_sd_nexp / sqrt(res$n_nexp)
  res$mean_ci_lo_nexp <- res$mean_nexp - qt(.975, res$n_nexp - 1) * res$mean_se_nexp
  res$mean_ci_up_nexp <- res$mean_nexp + qt(.975, res$n_nexp - 1) * res$mean_se_nexp

  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases

  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_pre_se_exp <- res$mean_pre_sd_exp / sqrt(res$n_exp)
  res$mean_pre_ci_lo_exp <- res$mean_pre_exp - qt(.975, res$n_exp - 1) * res$mean_pre_se_exp
  res$mean_pre_ci_up_exp <- res$mean_pre_exp + qt(.975, res$n_exp - 1) * res$mean_pre_se_exp

  res$mean_sd_exp <- res$sd_cases
  res$mean_se_exp <- res$mean_sd_exp / sqrt(res$n_exp)
  res$mean_ci_lo_exp <- res$mean_exp - qt(.975, res$n_exp - 1) * res$mean_se_exp
  res$mean_ci_up_exp <- res$mean_exp + qt(.975, res$n_exp - 1) * res$mean_se_exp

  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE

  smc_exp <- metafor::escalc(
    m1i = res$mean_pre_exp,
    sd1i = res$mean_pre_sd_exp,
    m2i = res$mean_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp,
    measure = "SMCRH"
  )
  smc_nexp <- metafor::escalc(
    m1i = res$mean_pre_nexp,
    sd1i = res$mean_pre_sd_nexp,
    m2i = res$mean_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp,
    measure = "SMCRH"
  )

  smcc_vg <- (-as.numeric(as.character(smc_exp$yi))) - (-as.numeric(as.character(smc_nexp$yi)))
  se_vg <- sqrt(smc_exp$vi + smc_nexp$vi)
  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "g",
    smd_to_cor = "viechtbauer",
    pre_post_to_smd = "bonett"
  ), digits = 11)
  ## d
  expect_equal(unique(es.mcv_g$info_used_crude), "means_ci_pre_post")
  expect_equal(es.mcv_g$es_crude, as.numeric(as.character(smcc_vg)),
    tolerance = 1e-10
  )
  expect_equal(es.mcv_g$se_crude, se_vg, tolerance = 1e-10)
})


### MD from PRE POST -----
test_that("MD - SD", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE

  smc_exp_md <- metafor::escalc(
    m1i = res$mean_pre_exp,
    sd1i = res$mean_pre_sd_exp,
    m2i = res$mean_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp,
    measure = "MC"
  )
  smc_nexp_md <- metafor::escalc(
    m1i = res$mean_pre_nexp,
    sd1i = res$mean_pre_sd_nexp,
    m2i = res$mean_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp,
    measure = "MC"
  )
  smcc_md <- (-smc_exp_md$yi) - (-smc_nexp_md$yi)
  se_md <- sqrt(smc_exp_md$vi + smc_nexp_md$vi)

  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
    measure = "md"
  ), digits = 11)
  expect_equal(unique(es.mcv_g$info_used_crude), "means_sd_pre_post")
  expect_equal(es.mcv_g$es_crude, as.numeric(as.character(smcc_md)),
    tolerance = 1e-10
  )
  expect_equal(es.mcv_g$se_crude, se_md, tolerance = 1e-10)
})
test_that("MD - SE", {
  res <- metaumbrella::df.SMC
  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls

  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls

  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_pre_se_nexp <- res$mean_pre_sd_nexp / sqrt(res$n_nexp)

  res$mean_sd_nexp <- res$sd_controls
  res$mean_se_nexp <- res$mean_sd_nexp / sqrt(res$n_nexp)

  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases

  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_pre_se_exp <- res$mean_pre_sd_exp / sqrt(res$n_exp)

  res$mean_sd_exp <- res$sd_cases
  res$mean_se_exp <- res$mean_sd_exp / sqrt(res$n_exp)

  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE

  smc_exp_md <- metafor::escalc(
    m1i = res$mean_pre_exp,
    sd1i = res$mean_pre_sd_exp,
    m2i = res$mean_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp,
    measure = "MC"
  )
  smc_nexp_md <- metafor::escalc(
    m1i = res$mean_pre_nexp,
    sd1i = res$mean_pre_sd_nexp,
    m2i = res$mean_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp,
    measure = "MC"
  )
  smcc_md <- (-smc_exp_md$yi) - (-smc_nexp_md$yi)
  se_md <- sqrt(smc_exp_md$vi + smc_nexp_md$vi)

  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se_pre_post",
    measure = "md"
  ), digits = 11)
  expect_equal(unique(es.mcv_g$info_used_crude), "means_se_pre_post")
  expect_equal(es.mcv_g$es_crude, as.numeric(as.character(smcc_md)),
    tolerance = 1e-10
  )
  expect_equal(es.mcv_g$se_crude, se_md, tolerance = 1e-10)
})
test_that("MD - CI", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls

  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls

  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_pre_se_nexp <- res$mean_pre_sd_nexp / sqrt(res$n_nexp)
  res$mean_pre_ci_lo_nexp <- res$mean_pre_nexp - qt(.975, res$n_nexp - 1) * res$mean_pre_se_nexp
  res$mean_pre_ci_up_nexp <- res$mean_pre_nexp + qt(.975, res$n_nexp - 1) * res$mean_pre_se_nexp

  res$mean_sd_nexp <- res$sd_controls
  res$mean_se_nexp <- res$mean_sd_nexp / sqrt(res$n_nexp)
  res$mean_ci_lo_nexp <- res$mean_nexp - qt(.975, res$n_nexp - 1) * res$mean_se_nexp
  res$mean_ci_up_nexp <- res$mean_nexp + qt(.975, res$n_nexp - 1) * res$mean_se_nexp

  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases

  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_pre_se_exp <- res$mean_pre_sd_exp / sqrt(res$n_exp)
  res$mean_pre_ci_lo_exp <- res$mean_pre_exp - qt(.975, res$n_exp - 1) * res$mean_pre_se_exp
  res$mean_pre_ci_up_exp <- res$mean_pre_exp + qt(.975, res$n_exp - 1) * res$mean_pre_se_exp

  res$mean_sd_exp <- res$sd_cases
  res$mean_se_exp <- res$mean_sd_exp / sqrt(res$n_exp)
  res$mean_ci_lo_exp <- res$mean_exp - qt(.975, res$n_exp - 1) * res$mean_se_exp
  res$mean_ci_up_exp <- res$mean_exp + qt(.975, res$n_exp - 1) * res$mean_se_exp

  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE

  smc_exp_md <- metafor::escalc(
    m1i = res$mean_pre_exp,
    sd1i = res$mean_pre_sd_exp,
    m2i = res$mean_exp,
    sd2i = res$mean_sd_exp,
    ni = res$n_exp,
    ri = res$r_pre_post_exp,
    measure = "MC"
  )
  smc_nexp_md <- metafor::escalc(
    m1i = res$mean_pre_nexp,
    sd1i = res$mean_pre_sd_nexp,
    m2i = res$mean_nexp,
    sd2i = res$mean_sd_nexp,
    ni = res$n_nexp,
    ri = res$r_pre_post_nexp,
    measure = "MC"
  )
  smcc_md <- (-smc_exp_md$yi) - (-smc_nexp_md$yi)
  se_md <- sqrt(smc_exp_md$vi + smc_nexp_md$vi)

  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci_pre_post",
    measure = "md"
  ), digits = 11)
  expect_equal(unique(es.mcv_g$info_used_crude), "means_ci_pre_post")
  expect_equal(es.mcv_g$es_crude, as.numeric(as.character(smcc_md)),
    tolerance = 1e-10
  )
  expect_equal(es.mcv_g$se_crude, se_md, tolerance = 1e-10)
})


### SMD from CHANGE -------
test_that("D - pre/post v change", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.5
  res$reverse_means_pre_post <- FALSE

  res$mean_change_exp <- res$mean_exp - res$mean_pre_exp
  res$mean_change_nexp <- res$mean_nexp - res$mean_pre_nexp
  res$mean_change_sd_exp <-
    sqrt(res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
      2 * res$r_pre_post_exp * res$mean_pre_sd_exp * res$mean_sd_exp)
  res$mean_change_sd_nexp <-
    sqrt(res$mean_pre_sd_nexp^2 + res$mean_sd_nexp^2 -
      2 * res$r_pre_post_nexp * res$mean_pre_sd_nexp * res$mean_sd_nexp)

  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
    pre_post_to_smd = "cooper", measure = "d"
  ), digits = 11)
  es.mcv_g_c <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd",
    measure = "d"
  ), digits = 11)
  expect_equal(unique(es.mcv_g$info_used_crude), "means_sd_pre_post")
  expect_equal(unique(es.mcv_g_c$info_used_crude), "mean_change_sd")
  expect_equal(es.mcv_g$es_crude, es.mcv_g_c$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_c$se_crude, tolerance = 1e-10)
})

test_that("G - pre/post v change", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.5
  res$reverse_means_pre_post <- FALSE

  res$mean_change_exp <- res$mean_exp - res$mean_pre_exp
  res$mean_change_nexp <- res$mean_nexp - res$mean_pre_nexp
  res$mean_change_sd_exp <-
    sqrt(res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
      2 * res$r_pre_post_exp * res$mean_pre_sd_exp * res$mean_sd_exp)
  res$mean_change_sd_nexp <-
    sqrt(res$mean_pre_sd_nexp^2 + res$mean_sd_nexp^2 -
      2 * res$r_pre_post_nexp * res$mean_pre_sd_nexp * res$mean_sd_nexp)

  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
    pre_post_to_smd = "cooper", measure = "g"
  ), digits = 11)
  es.mcv_g_c <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd",
    measure = "g"
  ), digits = 11)
  expect_equal(unique(es.mcv_g$info_used_crude), "means_sd_pre_post")
  expect_equal(unique(es.mcv_g_c$info_used_crude), "mean_change_sd")
  expect_equal(es.mcv_g$es_crude, es.mcv_g_c$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_c$se_crude, tolerance = 1e-10)
})

test_that("MD - pre/post v change", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE

  res$mean_change_exp <- res$mean_exp - res$mean_pre_exp
  res$mean_change_nexp <- res$mean_nexp - res$mean_pre_nexp
  res$mean_change_sd_exp <-
    sqrt(res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
      2 * res$r_pre_post_exp * res$mean_pre_sd_exp * res$mean_sd_exp)
  res$mean_change_sd_nexp <-
    sqrt(res$mean_pre_sd_nexp^2 + res$mean_sd_nexp^2 -
      2 * res$r_pre_post_nexp * res$mean_pre_sd_nexp * res$mean_sd_nexp)

  ## metaconvert
  es.mcv_g <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
    pre_post_to_smd = "cooper", measure = "md"
  ), digits = 11)
  es.mcv_g_c <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd",
    measure = "md"
  ), digits = 11)
  expect_equal(unique(es.mcv_g$info_used_crude), "means_sd_pre_post")
  expect_equal(unique(es.mcv_g_c$info_used_crude), "mean_change_sd")
  expect_equal(es.mcv_g$es_crude, es.mcv_g_c$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_c$se_crude, tolerance = 1e-10)
})


### REVERSE ------

test_that("REVERSE - pre/post - BONNETT", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls

  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls

  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_pre_se_nexp <- res$mean_pre_sd_nexp / sqrt(res$n_nexp)
  res$mean_pre_ci_lo_nexp <- res$mean_pre_nexp - qt(.975, res$n_nexp - 1) * res$mean_pre_se_nexp
  res$mean_pre_ci_up_nexp <- res$mean_pre_nexp + qt(.975, res$n_nexp - 1) * res$mean_pre_se_nexp

  res$mean_sd_nexp <- res$sd_controls
  res$mean_se_nexp <- res$mean_sd_nexp / sqrt(res$n_nexp)
  res$mean_ci_lo_nexp <- res$mean_nexp - qt(.975, res$n_nexp - 1) * res$mean_se_nexp
  res$mean_ci_up_nexp <- res$mean_nexp + qt(.975, res$n_nexp - 1) * res$mean_se_nexp

  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases

  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_pre_se_exp <- res$mean_pre_sd_exp / sqrt(res$n_exp)
  res$mean_pre_ci_lo_exp <- res$mean_pre_exp - qt(.975, res$n_exp - 1) * res$mean_pre_se_exp
  res$mean_pre_ci_up_exp <- res$mean_pre_exp + qt(.975, res$n_exp - 1) * res$mean_pre_se_exp

  res$mean_sd_exp <- res$sd_cases
  res$mean_se_exp <- res$mean_sd_exp / sqrt(res$n_exp)
  res$mean_ci_lo_exp <- res$mean_exp - qt(.975, res$n_exp - 1) * res$mean_se_exp
  res$mean_ci_up_exp <- res$mean_exp + qt(.975, res$n_exp - 1) * res$mean_se_exp

  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE
  dat <- res

  es.mcv_m_sd_md <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "md"), digits = 11)
  es.mcv_m_se_md <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "md"), digits = 11)
  es.mcv_m_ci_md <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "md"), digits = 11)

  es.mcv_m_sd_d <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "d"), digits = 11)
  es.mcv_m_se_d <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "d"), digits = 11)
  es.mcv_m_ci_d <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "d"), digits = 11)

  es.mcv_m_sd_g <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "g"), digits = 11)
  es.mcv_m_se_g <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "g"), digits = 11)
  es.mcv_m_ci_g <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "g"), digits = 11)

  es.mcv_m_sd_r <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "r"), digits = 11)
  es.mcv_m_se_r <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "r"), digits = 11)
  es.mcv_m_ci_r <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "r"), digits = 11)

  es.mcv_m_sd_z <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "z"), digits = 11)
  es.mcv_m_se_z <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "z"), digits = 11)
  es.mcv_m_ci_z <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "z"), digits = 11)

  es.mcv_m_sd_or <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "logor"), digits = 11)
  es.mcv_m_se_or <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "logor"), digits = 11)
  es.mcv_m_ci_or <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "logor"), digits = 11)

  dat$reverse_means_pre_post <- TRUE

  es.mcv_m_sd_md_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "md"), digits = 11)
  es.mcv_m_se_md_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "md"), digits = 11)
  es.mcv_m_ci_md_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "md"), digits = 11)

  es.mcv_m_sd_d_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "d"), digits = 11)
  es.mcv_m_se_d_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "d"), digits = 11)
  es.mcv_m_ci_d_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "d"), digits = 11)

  es.mcv_m_sd_g_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "g"), digits = 11)
  es.mcv_m_se_g_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "g"), digits = 11)
  es.mcv_m_ci_g_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "g"), digits = 11)

  es.mcv_m_sd_r_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "r"), digits = 11)
  es.mcv_m_se_r_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "r"), digits = 11)
  es.mcv_m_ci_r_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "r"), digits = 11)

  es.mcv_m_sd_z_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "z"), digits = 11)
  es.mcv_m_se_z_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "z"), digits = 11)
  es.mcv_m_ci_z_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "z"), digits = 11)

  es.mcv_m_sd_or_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "logor"), digits = 11)
  es.mcv_m_se_or_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "logor"), digits = 11)
  es.mcv_m_ci_or_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "bonett",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "logor"), digits = 11)

  expect_equal(es.mcv_m_sd_d$info_used_crude, es.mcv_m_sd_d_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_d$es_crude, -es.mcv_m_sd_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$es_crude, -es.mcv_m_se_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$es_crude, -es.mcv_m_ci_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_d$se_crude, es.mcv_m_sd_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$se_crude, es.mcv_m_se_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$se_crude, es.mcv_m_ci_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_md$info_used_crude, es.mcv_m_sd_md_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_md$es_crude, -es.mcv_m_sd_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_md$es_crude, -es.mcv_m_se_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_md$es_crude, -es.mcv_m_ci_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_md$se_crude, es.mcv_m_sd_md_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_md$se_crude, es.mcv_m_se_md_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_md$se_crude, es.mcv_m_ci_md_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_g$info_used_crude, es.mcv_m_sd_g_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_g$es_crude, -es.mcv_m_sd_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_g$es_crude, -es.mcv_m_se_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_g$es_crude, -es.mcv_m_ci_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_g$se_crude, es.mcv_m_sd_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_g$se_crude, es.mcv_m_se_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_g$se_crude, es.mcv_m_ci_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_r$info_used_crude, es.mcv_m_sd_r_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_r$es_crude, -es.mcv_m_sd_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_r$es_crude, -es.mcv_m_se_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_r$es_crude, -es.mcv_m_ci_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_r$se_crude, es.mcv_m_sd_r_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_r$se_crude, es.mcv_m_se_r_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_r$se_crude, es.mcv_m_ci_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_z$info_used_crude, es.mcv_m_sd_z_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_z$es_crude, -es.mcv_m_sd_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_z$es_crude, -es.mcv_m_se_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_z$es_crude, -es.mcv_m_ci_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_z$se_crude, es.mcv_m_sd_z_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_z$se_crude, es.mcv_m_se_z_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_z$se_crude, es.mcv_m_ci_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_or$info_used_crude, es.mcv_m_sd_or_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_or$es_crude, -es.mcv_m_sd_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_or$es_crude, -es.mcv_m_se_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_or$es_crude, -es.mcv_m_ci_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_or$se_crude, es.mcv_m_sd_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_or$se_crude, es.mcv_m_se_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_or$se_crude, es.mcv_m_ci_or_rv$se_crude, tolerance = 1e-10)
})

test_that("REVERSE - pre/post - COOPER", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls

  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls

  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_pre_se_nexp <- res$mean_pre_sd_nexp / sqrt(res$n_nexp)
  res$mean_pre_ci_lo_nexp <- res$mean_pre_nexp - qt(.975, res$n_nexp - 1) * res$mean_pre_se_nexp
  res$mean_pre_ci_up_nexp <- res$mean_pre_nexp + qt(.975, res$n_nexp - 1) * res$mean_pre_se_nexp

  res$mean_sd_nexp <- res$sd_controls
  res$mean_se_nexp <- res$mean_sd_nexp / sqrt(res$n_nexp)
  res$mean_ci_lo_nexp <- res$mean_nexp - qt(.975, res$n_nexp - 1) * res$mean_se_nexp
  res$mean_ci_up_nexp <- res$mean_nexp + qt(.975, res$n_nexp - 1) * res$mean_se_nexp

  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases

  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_pre_se_exp <- res$mean_pre_sd_exp / sqrt(res$n_exp)
  res$mean_pre_ci_lo_exp <- res$mean_pre_exp - qt(.975, res$n_exp - 1) * res$mean_pre_se_exp
  res$mean_pre_ci_up_exp <- res$mean_pre_exp + qt(.975, res$n_exp - 1) * res$mean_pre_se_exp

  res$mean_sd_exp <- res$sd_cases
  res$mean_se_exp <- res$mean_sd_exp / sqrt(res$n_exp)
  res$mean_ci_lo_exp <- res$mean_exp - qt(.975, res$n_exp - 1) * res$mean_se_exp
  res$mean_ci_up_exp <- res$mean_exp + qt(.975, res$n_exp - 1) * res$mean_se_exp

  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_means_pre_post <- FALSE
  dat <- res

  es.mcv_m_sd_md <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "md"), digits = 11)
  es.mcv_m_se_md <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "md"), digits = 11)
  es.mcv_m_ci_md <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "md"), digits = 11)

  es.mcv_m_sd_d <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "d"), digits = 11)
  es.mcv_m_se_d <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "d"), digits = 11)
  es.mcv_m_ci_d <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "d"), digits = 11)

  es.mcv_m_sd_g <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "g"), digits = 11)
  es.mcv_m_se_g <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "g"), digits = 11)
  es.mcv_m_ci_g <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "g"), digits = 11)

  es.mcv_m_sd_r <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "r"), digits = 11)
  es.mcv_m_se_r <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "r"), digits = 11)
  es.mcv_m_ci_r <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "r"), digits = 11)

  es.mcv_m_sd_z <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "z"), digits = 11)
  es.mcv_m_se_z <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "z"), digits = 11)
  es.mcv_m_ci_z <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "z"), digits = 11)

  es.mcv_m_sd_or <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "logor"), digits = 11)
  es.mcv_m_se_or <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "logor"), digits = 11)
  es.mcv_m_ci_or <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "logor"), digits = 11)

  dat$reverse_means_pre_post <- TRUE

  es.mcv_m_sd_md_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "md"), digits = 11)
  es.mcv_m_se_md_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "md"), digits = 11)
  es.mcv_m_ci_md_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "md"), digits = 11)

  es.mcv_m_sd_d_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "d"), digits = 11)
  es.mcv_m_se_d_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "d"), digits = 11)
  es.mcv_m_ci_d_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "d"), digits = 11)

  es.mcv_m_sd_g_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "g"), digits = 11)
  es.mcv_m_se_g_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "g"), digits = 11)
  es.mcv_m_ci_g_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "g"), digits = 11)

  es.mcv_m_sd_r_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "r"), digits = 11)
  es.mcv_m_se_r_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "r"), digits = 11)
  es.mcv_m_ci_r_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "r"), digits = 11)

  es.mcv_m_sd_z_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "z"), digits = 11)
  es.mcv_m_se_z_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "z"), digits = 11)
  es.mcv_m_ci_z_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "z"), digits = 11)

  es.mcv_m_sd_or_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_sd_pre_post", measure = "logor"), digits = 11)
  es.mcv_m_se_or_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_se_pre_post", measure = "logor"), digits = 11)
  es.mcv_m_ci_or_rv <- summary(convert_df(dat, verbose = FALSE, pre_post_to_smd = "cooper",es_selected = "hierarchy", hierarchy = "means_ci_pre_post", measure = "logor"), digits = 11)

  expect_equal(es.mcv_m_sd_d$info_used_crude, es.mcv_m_sd_d_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_d$es_crude, -es.mcv_m_sd_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$es_crude, -es.mcv_m_se_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$es_crude, -es.mcv_m_ci_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_d$se_crude, es.mcv_m_sd_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$se_crude, es.mcv_m_se_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$se_crude, es.mcv_m_ci_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_md$info_used_crude, es.mcv_m_sd_md_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_md$es_crude, -es.mcv_m_sd_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_md$es_crude, -es.mcv_m_se_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_md$es_crude, -es.mcv_m_ci_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_md$se_crude, es.mcv_m_sd_md_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_md$se_crude, es.mcv_m_se_md_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_md$se_crude, es.mcv_m_ci_md_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_g$info_used_crude, es.mcv_m_sd_g_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_g$es_crude, -es.mcv_m_sd_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_g$es_crude, -es.mcv_m_se_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_g$es_crude, -es.mcv_m_ci_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_g$se_crude, es.mcv_m_sd_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_g$se_crude, es.mcv_m_se_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_g$se_crude, es.mcv_m_ci_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_r$info_used_crude, es.mcv_m_sd_r_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_r$es_crude, -es.mcv_m_sd_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_r$es_crude, -es.mcv_m_se_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_r$es_crude, -es.mcv_m_ci_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_r$se_crude, es.mcv_m_sd_r_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_r$se_crude, es.mcv_m_se_r_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_r$se_crude, es.mcv_m_ci_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_z$info_used_crude, es.mcv_m_sd_z_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_z$es_crude, -es.mcv_m_sd_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_z$es_crude, -es.mcv_m_se_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_z$es_crude, -es.mcv_m_ci_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_z$se_crude, es.mcv_m_sd_z_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_z$se_crude, es.mcv_m_se_z_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_z$se_crude, es.mcv_m_ci_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_or$info_used_crude, es.mcv_m_sd_or_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_or$es_crude, -es.mcv_m_sd_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_or$es_crude, -es.mcv_m_se_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_or$es_crude, -es.mcv_m_ci_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_or$se_crude, es.mcv_m_sd_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_or$se_crude, es.mcv_m_se_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_or$se_crude, es.mcv_m_ci_or_rv$se_crude, tolerance = 1e-10)
})
test_that("REVERSE - change - bonett", {
  res <- metaumbrella::df.SMC

  res$n_exp <- res$n_cases
  res$n_nexp <- res$n_controls
  res$mean_pre_nexp <- res$mean_pre_controls
  res$mean_nexp <- res$mean_controls
  res$mean_pre_sd_nexp <- res$sd_pre_controls
  res$mean_sd_nexp <- res$sd_controls
  res$mean_pre_exp <- res$mean_pre_cases
  res$mean_exp <- res$mean_cases
  res$mean_pre_sd_exp <- res$sd_pre_cases
  res$mean_sd_exp <- res$sd_cases
  res$r_pre_post_nexp <- res$r_pre_post_exp <- 0.8
  res$reverse_mean_change <- FALSE

  res$mean_change_exp <- res$mean_exp - res$mean_pre_exp
  res$mean_change_nexp <- res$mean_nexp - res$mean_pre_nexp
  res$mean_change_sd_exp <-
    sqrt(res$mean_pre_sd_exp^2 + res$mean_sd_exp^2 -
      2 * res$r_pre_post_exp * res$mean_pre_sd_exp * res$mean_sd_exp)
  res$mean_change_sd_nexp <-
    sqrt(res$mean_pre_sd_nexp^2 + res$mean_sd_nexp^2 -
      2 * res$r_pre_post_nexp * res$mean_pre_sd_nexp * res$mean_sd_nexp)
  dat <- res

  es.mcv_m_sd_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "md"), digits = 11)
  es.mcv_m_sd_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "d"), digits = 11)
  es.mcv_m_sd_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "g"), digits = 11)
  es.mcv_m_sd_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "r"), digits = 11)
  es.mcv_m_sd_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "z"), digits = 11)
  es.mcv_m_sd_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "logor"), digits = 11)

  dat$reverse_mean_change <- TRUE

  es.mcv_m_sd_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "md"), digits = 11)
  es.mcv_m_sd_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "d"), digits = 11)
  es.mcv_m_sd_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "g"), digits = 11)
  es.mcv_m_sd_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "r"), digits = 11)
  es.mcv_m_sd_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "z"), digits = 11)
  es.mcv_m_sd_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "mean_change_sd", measure = "logor"), digits = 11)

  expect_equal(es.mcv_m_sd_d$info_used_crude, es.mcv_m_sd_d_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_d$es_crude, -es.mcv_m_sd_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_d$se_crude, es.mcv_m_sd_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_md$info_used_crude, es.mcv_m_sd_md_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_md$es_crude, -es.mcv_m_sd_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_md$se_crude, es.mcv_m_sd_md_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_g$info_used_crude, es.mcv_m_sd_g_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_g$es_crude, -es.mcv_m_sd_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_g$se_crude, es.mcv_m_sd_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_r$info_used_crude, es.mcv_m_sd_r_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_r$es_crude, -es.mcv_m_sd_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_r$se_crude, es.mcv_m_sd_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_z$info_used_crude, es.mcv_m_sd_z_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_z$es_crude, -es.mcv_m_sd_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_z$se_crude, es.mcv_m_sd_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_or$info_used_crude, es.mcv_m_sd_or_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_or$es_crude, -es.mcv_m_sd_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_or$se_crude, es.mcv_m_sd_or_rv$se_crude, tolerance = 1e-10)
})
