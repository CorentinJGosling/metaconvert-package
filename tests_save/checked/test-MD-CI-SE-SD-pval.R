# MD => SMD ===============
test_that("MD => d", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp) & mean_sd_nexp != 0)

  dat_MD <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "MD", vtype="HO", data = dat, digits = 11
  ))

  dat$md <- dat_MD$yi
  dat$md_sd <- sqrt(dat_MD$vi) / sqrt(1 / dat$n_exp + 1 / dat$n_nexp)

  dat$reverse_md <- FALSE
  dat$reverse_means <- FALSE

  es.mcv_md_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "d"), digits = 11)
  es.mcv_means_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_md_sd$info_used_crude), "md_sd")
  expect_equal(unique(es.mcv_means_sd$info_used_crude), "means_sd")
  expect_equal(es.mcv_md_sd$es_crude, es.mcv_means_sd$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_means_sd$se_crude, tolerance = 1e-10)
})

# Means => MD ===============
test_that("MD from means", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_md <- FALSE
  dat$reverse_means <- FALSE

  dat_MD <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "MD", data = dat, digits = 11
  ))

  es.mcv_md_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "md"), digits = 11)

  expect_equal(unique(es.mcv_md_sd$info_used_crude), "means_sd")
  expect_equal(es.mcv_md_sd$es_crude, as.numeric(dat_MD$yi), tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude^2, dat_MD$vi, tolerance = 1e-10)
})

# MD + SD => MD + ? ===========
test_that("D - MD+SD vs MD+SE/CI/p", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp))

  dat_MD <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "MD", data = dat, digits = 11
  ))

  dat$md <- dat_MD$yi
  dat$md_sd <- sqrt(dat_MD$vi) / sqrt(1 / dat$n_exp + 1 / dat$n_nexp)
  dat$md_se <- sqrt(dat_MD$vi)
  dat$md_ci_lo <- with(dat, md - qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_ci_up <- with(dat, md + qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_t <- abs(dat$md / dat$md_se)
  dat$md_pval <- 2 * pt(dat$md_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)

  dat <- subset(dat, !is.na(md_pval) & md_pval != 1)
  dat$reverse_md <- FALSE
  dat$reverse_md_pval <- FALSE

  es.mcv_md_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "d"), digits = 11)
  es.mcv_md_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "d"), digits = 11)
  es.mcv_md_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "d"), digits = 11)
  es.mcv_md_pval <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "d"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_md_sd$info_used_crude), "md_sd")
  expect_equal(unique(es.mcv_md_se$info_used_crude), "md_se")
  expect_equal(unique(es.mcv_md_ci$info_used_crude), "md_ci")
  expect_equal(unique(es.mcv_md_pval$info_used_crude), "md_pval")

  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_se$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_se$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_ci$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_ci$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_pval$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_pval$se_crude, tolerance = 1e-10)
})
test_that("OR - MD+SD vs MD+SE/CI/p", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp))

  dat_MD <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "MD", data = dat, digits = 11
  ))

  dat$md <- dat_MD$yi
  dat$md_sd <- sqrt(dat_MD$vi) / sqrt(1 / dat$n_exp + 1 / dat$n_nexp)
  dat$md_se <- sqrt(dat_MD$vi)
  dat$md_ci_lo <- with(dat, md - qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_ci_up <- with(dat, md + qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_t <- abs(dat$md / dat$md_se)
  dat$md_pval <- 2 * pt(dat$md_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)

  dat <- subset(dat, !is.na(md_pval) & md_pval != 1)
  dat$reverse_md <- FALSE
  dat$reverse_md_pval <- FALSE

  es.mcv_md_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "or"), digits = 11)
  es.mcv_md_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "or"), digits = 11)
  es.mcv_md_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "or"), digits = 11)
  es.mcv_md_pval <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "or"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_md_sd$info_used_crude), "md_sd")
  expect_equal(unique(es.mcv_md_se$info_used_crude), "md_se")
  expect_equal(unique(es.mcv_md_ci$info_used_crude), "md_ci")
  expect_equal(unique(es.mcv_md_pval$info_used_crude), "md_pval")

  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_se$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_se$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_ci$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_ci$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_pval$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_pval$se_crude, tolerance = 1e-10)
})
test_that("Z - MD+SD vs MD+SE/CI/p", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp))

  dat_MD <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "MD", data = dat, digits = 11
  ))

  dat$md <- dat_MD$yi
  dat$md_sd <- sqrt(dat_MD$vi) / sqrt(1 / dat$n_exp + 1 / dat$n_nexp)
  dat$md_se <- sqrt(dat_MD$vi)
  dat$md_ci_lo <- with(dat, md - qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_ci_up <- with(dat, md + qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_t <- abs(dat$md / dat$md_se)
  dat$md_pval <- 2 * pt(dat$md_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)

  dat <- subset(dat, !is.na(md_pval) & md_pval != 1)
  dat$reverse_md <- FALSE
  dat$reverse_md_pval <- FALSE

  es.mcv_md_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "z"), digits = 11)
  es.mcv_md_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "z"), digits = 11)
  es.mcv_md_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "z"), digits = 11)
  es.mcv_md_pval <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "z"), digits = 11)

  ## test ES
  # z declined on an out-of-range biserial r -- see helper-z-declined.R
  expect_route(es.mcv_md_sd$info_used_crude, "md_sd")
  # z declined on an out-of-range biserial r -- see helper-z-declined.R
  expect_route(es.mcv_md_se$info_used_crude, "md_se")
  # z declined on an out-of-range biserial r -- see helper-z-declined.R
  expect_route(es.mcv_md_ci$info_used_crude, "md_ci")
  # z declined on an out-of-range biserial r -- see helper-z-declined.R
  expect_route(es.mcv_md_pval$info_used_crude, "md_pval")

  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_se$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_se$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_ci$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_ci$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$es_crude, es.mcv_md_pval$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd$se_crude, es.mcv_md_pval$se_crude, tolerance = 1e-10)
})

test_that("MD - Reverse standard", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp))

  dat_MD <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "MD", data = dat, digits = 11
  ))

  dat$md <- dat_MD$yi
  dat$md_sd <- sqrt(dat_MD$vi) / sqrt(1 / dat$n_exp + 1 / dat$n_nexp)
  dat$md_se <- sqrt(dat_MD$vi)
  dat$md_ci_lo <- with(dat, md - qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_ci_up <- with(dat, md + qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_t <- abs(dat$md / dat$md_se)
  dat$md_pval <- 2 * pt(dat$md_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)

  dat <- subset(dat, !is.na(md_pval) & md_pval != 1)
  dat$reverse_md <- FALSE
  dat$reverse_md_pval <- FALSE

  es.mcv_md_sd_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "md"), digits = 11)
  es.mcv_md_se_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "md"), digits = 11)
  es.mcv_md_ci_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "md"), digits = 11)

  es.mcv_md_sd_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "d"), digits = 11)
  es.mcv_md_se_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "d"), digits = 11)
  es.mcv_md_ci_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "d"), digits = 11)

  es.mcv_md_sd_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "g"), digits = 11)
  es.mcv_md_se_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "g"), digits = 11)
  es.mcv_md_ci_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "g"), digits = 11)

  es.mcv_md_sd_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "r",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_md_se_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "r",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_md_ci_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "r",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)

  es.mcv_md_sd_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "z",
                                       smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_md_se_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "z",
                                       smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_md_ci_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "z",
                                       smd_to_cor = "lipsey_cooper"), digits = 11)

  es.mcv_md_sd_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "r",
                                        smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_md_se_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "r",
                                        smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_md_ci_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "r",
                                        smd_to_cor = "viechtbauer"), digits = 11)

  es.mcv_md_sd_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "z",
                                        smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_md_se_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "z",
                                        smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_md_ci_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "z",
                                        smd_to_cor = "viechtbauer"), digits = 11)

  es.mcv_md_sd_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "logor"), digits = 11)
  es.mcv_md_se_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "logor"), digits = 11)
  es.mcv_md_ci_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "logor"), digits = 11)

  dat$reverse_md <- TRUE

  es.mcv_md_sd_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "md"), digits = 11)
  es.mcv_md_se_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "md"), digits = 11)
  es.mcv_md_ci_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "md"), digits = 11)

  es.mcv_md_sd_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "d"), digits = 11)
  es.mcv_md_se_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "d"), digits = 11)
  es.mcv_md_ci_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "d"), digits = 11)

  es.mcv_md_sd_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "g"), digits = 11)
  es.mcv_md_se_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "g"), digits = 11)
  es.mcv_md_ci_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "g"), digits = 11)

  es.mcv_md_sd_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "r",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_md_se_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "r",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_md_ci_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "r",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)

  es.mcv_md_sd_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "z",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_md_se_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "z",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_md_ci_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "z",
                                        smd_to_cor = "lipsey_cooper"), digits = 11)

  es.mcv_md_sd_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "r",
                                        smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_md_se_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "r",
                                        smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_md_ci_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "r",
                                        smd_to_cor = "viechtbauer"), digits = 11)

  es.mcv_md_sd_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "z",
                                        smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_md_se_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "z",
                                        smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_md_ci_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "z",
                                        smd_to_cor = "viechtbauer"), digits = 11)

  es.mcv_md_sd_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_sd", measure = "logor"), digits = 11)
  es.mcv_md_se_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_se", measure = "logor"), digits = 11)
  es.mcv_md_ci_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_ci", measure = "logor"), digits = 11)


  expect_route(c(es.mcv_md_se_d$info_used_crude,
    es.mcv_md_se_md$info_used_crude,
    es.mcv_md_se_g$info_used_crude,
    es.mcv_md_se_r1$info_used_crude,
    es.mcv_md_se_r2$info_used_crude,
    es.mcv_md_se_z1$info_used_crude,
    es.mcv_md_se_z2$info_used_crude,
    es.mcv_md_se_d_rv$info_used_crude,
    es.mcv_md_se_md_rv$info_used_crude,
    es.mcv_md_se_g_rv$info_used_crude,
    es.mcv_md_se_r1_rv$info_used_crude,
    es.mcv_md_se_r2_rv$info_used_crude,
    es.mcv_md_se_z1_rv$info_used_crude,
    es.mcv_md_se_z2_rv$info_used_crude), "md_se")

  expect_route(c(es.mcv_md_sd_d$info_used_crude,
    es.mcv_md_sd_md$info_used_crude,
    es.mcv_md_sd_g$info_used_crude,
    es.mcv_md_sd_r1$info_used_crude,
    es.mcv_md_sd_r2$info_used_crude,
    es.mcv_md_sd_z1$info_used_crude,
    es.mcv_md_sd_z2$info_used_crude,
    es.mcv_md_sd_d_rv$info_used_crude,
    es.mcv_md_sd_md_rv$info_used_crude,
    es.mcv_md_sd_g_rv$info_used_crude,
    es.mcv_md_sd_r1_rv$info_used_crude,
    es.mcv_md_sd_r2_rv$info_used_crude,
    es.mcv_md_sd_z1_rv$info_used_crude,
    es.mcv_md_sd_z2_rv$info_used_crude), "md_sd")

  expect_route(c(es.mcv_md_ci_d$info_used_crude,
    es.mcv_md_ci_md$info_used_crude,
    es.mcv_md_ci_g$info_used_crude,
    es.mcv_md_ci_r1$info_used_crude,
    es.mcv_md_ci_r2$info_used_crude,
    es.mcv_md_ci_z1$info_used_crude,
    es.mcv_md_ci_z2$info_used_crude,
    es.mcv_md_ci_d_rv$info_used_crude,
    es.mcv_md_ci_md_rv$info_used_crude,
    es.mcv_md_ci_g_rv$info_used_crude,
    es.mcv_md_ci_r1_rv$info_used_crude,
    es.mcv_md_ci_r2_rv$info_used_crude,
    es.mcv_md_ci_z1_rv$info_used_crude,
    es.mcv_md_ci_z2$info_used_crude), "md_ci")

  expect_equal(es.mcv_md_sd_d$es_crude, -es.mcv_md_sd_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_d$es_crude, -es.mcv_md_se_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_d$es_crude, -es.mcv_md_ci_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd_d$se_crude, es.mcv_md_sd_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_d$se_crude, es.mcv_md_se_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_d$se_crude, es.mcv_md_ci_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_sd_md$info_used_crude, es.mcv_md_sd_md_rv$info_used_crude)
  expect_equal(es.mcv_md_sd_md$es_crude, -es.mcv_md_sd_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_md$es_crude, -es.mcv_md_se_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_md$es_crude, -es.mcv_md_ci_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd_md$se_crude, es.mcv_md_sd_md_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_md$se_crude, es.mcv_md_se_md_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_md$se_crude, es.mcv_md_ci_md_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_sd_g$es_crude, -es.mcv_md_sd_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_g$es_crude, -es.mcv_md_se_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_g$es_crude, -es.mcv_md_ci_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd_g$se_crude, es.mcv_md_sd_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_g$se_crude, es.mcv_md_se_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_g$se_crude, es.mcv_md_ci_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_sd_r1$es_crude, -es.mcv_md_sd_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_r1$es_crude, -es.mcv_md_se_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_r1$es_crude, -es.mcv_md_ci_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd_r1$se_crude, es.mcv_md_sd_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_r1$se_crude, es.mcv_md_se_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_r1$se_crude, es.mcv_md_ci_r1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_sd_r2$es_crude, -es.mcv_md_sd_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_r2$es_crude, -es.mcv_md_se_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_r2$es_crude, -es.mcv_md_ci_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd_r2$se_crude, es.mcv_md_sd_r2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_r2$se_crude, es.mcv_md_se_r2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_r2$se_crude, es.mcv_md_ci_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_sd_z1$es_crude, -es.mcv_md_sd_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_z1$es_crude, -es.mcv_md_se_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_z1$es_crude, -es.mcv_md_ci_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd_z1$se_crude, es.mcv_md_sd_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_z1$se_crude, es.mcv_md_se_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_z1$se_crude, es.mcv_md_ci_z1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_sd_z2$es_crude, -es.mcv_md_sd_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_z2$es_crude, -es.mcv_md_se_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_z2$es_crude, -es.mcv_md_ci_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd_z2$se_crude, es.mcv_md_sd_z2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_z2$se_crude, es.mcv_md_se_z2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_z2$se_crude, es.mcv_md_ci_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_sd_or$es_crude, -es.mcv_md_sd_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_or$es_crude, -es.mcv_md_se_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_or$es_crude, -es.mcv_md_ci_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_sd_or$se_crude, es.mcv_md_sd_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_se_or$se_crude, es.mcv_md_se_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_ci_or$se_crude, es.mcv_md_ci_or_rv$se_crude, tolerance = 1e-10)
})

test_that("MD - Reverse pval", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp))

  dat_MD <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "MD", data = dat, digits = 11
  ))

  dat$md <- dat_MD$yi
  dat$md_sd <- sqrt(dat_MD$vi) / sqrt(1 / dat$n_exp + 1 / dat$n_nexp)
  dat$md_se <- sqrt(dat_MD$vi)
  dat$md_ci_lo <- with(dat, md - qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_ci_up <- with(dat, md + qt(.975, n_exp + n_nexp - 2) * md_se)
  dat$md_t <- abs(dat$md / dat$md_se)
  dat$md_pval <- 2 * pt(dat$md_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)

  dat <- subset(dat, !is.na(md_pval) & md_pval != 1)
  dat$reverse_md <- FALSE
  dat$reverse_md <- FALSE

  es.mcv_md_pval_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "logor"), digits = 11)
  es.mcv_md_pval_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval",
                                          smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_md_pval_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval",
                                          smd_to_cor = "viechtbauer",  measure = "r"), digits = 11)
  es.mcv_md_pval_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval",
                                          smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)
  es.mcv_md_pval_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval",
                                          smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_md_pval_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "g"), digits = 11)
  es.mcv_md_pval_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "d"), digits = 11)
  es.mcv_md_pval_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "md"), digits = 11)

  dat$reverse_md <- TRUE
  es.mcv_md_pval_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "md"), digits = 11)
  es.mcv_md_pval_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "g"), digits = 11)
  es.mcv_md_pval_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval",
                                          smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_md_pval_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval",
                                          smd_to_cor = "viechtbauer",  measure = "r"), digits = 11)
  es.mcv_md_pval_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval",
                                          smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)
  es.mcv_md_pval_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval",
                                          smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_md_pval_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "logor"), digits = 11)
  es.mcv_md_pval_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "md_pval", measure = "d"), digits = 11)

  expect_route(c(es.mcv_md_pval_d$info_used_crude,
                    es.mcv_md_pval_md$info_used_crude,
                    es.mcv_md_pval_g$info_used_crude,
                    es.mcv_md_pval_r1$info_used_crude,
                    es.mcv_md_pval_r2$info_used_crude,
                    es.mcv_md_pval_z1$info_used_crude,
                    es.mcv_md_pval_z2$info_used_crude,
                    es.mcv_md_pval_d_rv$info_used_crude,
                    es.mcv_md_pval_md_rv$info_used_crude,
                    es.mcv_md_pval_g_rv$info_used_crude,
                    es.mcv_md_pval_r1_rv$info_used_crude,
                    es.mcv_md_pval_r2_rv$info_used_crude,
                    es.mcv_md_pval_z1_rv$info_used_crude,
                    es.mcv_md_pval_z2$info_used_crude), "md_pval")

  expect_equal(es.mcv_md_pval_d$es_crude, -es.mcv_md_pval_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_d$se_crude, es.mcv_md_pval_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_pval_md$es_crude, -es.mcv_md_pval_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_md$se_crude, es.mcv_md_pval_md_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_d$es_crude, -es.mcv_md_pval_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_d$se_crude, es.mcv_md_pval_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_pval_g$se_crude, es.mcv_md_pval_g_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_g$es_crude, -es.mcv_md_pval_g_rv$es_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_pval_r1$es_crude, -es.mcv_md_pval_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_r1$se_crude, es.mcv_md_pval_r1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_pval_z1$es_crude, -es.mcv_md_pval_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_z1$se_crude, es.mcv_md_pval_z1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_pval_r2$es_crude, -es.mcv_md_pval_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_r2$se_crude, es.mcv_md_pval_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_pval_z2$es_crude, -es.mcv_md_pval_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_z2$se_crude, es.mcv_md_pval_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_md_pval_or$se_crude, es.mcv_md_pval_or_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_md_pval_or$es_crude, -es.mcv_md_pval_or_rv$es_crude, tolerance = 1e-10)
})
