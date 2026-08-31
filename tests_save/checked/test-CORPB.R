# RPB to COR (VIECHTBAUER) --------
test_that("rpb to R bis (viecht)", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp))
  dat$reverse_means <- FALSE

  dat$pt_bis_r <- as.numeric(metafor::escalc(
    vtype = "CS", measure = "RPB",
    m1i = mean_exp, m2i = mean_nexp, sd1i = mean_sd_exp, sd2i = mean_sd_nexp,
    n1i = n_exp, n2i = n_nexp, data = dat
  )$yi)

  dat$pt_bis_r_pval <- with(dat, 2 * pt(pt_bis_r * sqrt((n_exp + n_nexp - 2) / (1 - pt_bis_r^2)),
    n_exp + n_nexp - 2,
    lower.tail = FALSE
  ))
  dat <- subset(dat, pt_bis_r_pval < 1)

  es.mcv_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "r"), digits = 11)
  es.mcv_mean <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "r",
    smd_to_cor = "viechtbauer"
  ), digits = 11)

  expect_equal(unique(es.mcv_r$info_used_crude), "pt_bis_r")
  expect_equal(unique(es.mcv_mean$info_used_crude), "means_sd")
  expect_equal(es.mcv_r$es_crude, es.mcv_mean$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, es.mcv_mean$se_crude, tolerance = 1e-10)
})
test_that("rpb to Z (viecht)", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp))
  dat$reverse_means <- FALSE
  dat$pt_bis_r <- as.numeric(metafor::escalc(
    vtype = "CS", measure = "RPB",
    m1i = mean_exp, m2i = mean_nexp, sd1i = mean_sd_exp, sd2i = mean_sd_nexp,
    n1i = n_exp, n2i = n_nexp, data = dat
  )$yi)
  dat$pt_bis_r_pval <- with(dat, 2 * pt(pt_bis_r * sqrt((n_exp + n_nexp - 2) / (1 - pt_bis_r^2)), n_exp + n_nexp - 2, lower.tail = FALSE))
  dat <- subset(dat, pt_bis_r_pval < 1)


  es.mcv_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "z"), digits = 11)
  es.mcv_mean <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "z",
    smd_to_cor = "viechtbauer"
  ), digits = 11)

  # z declined on an out-of-range biserial r -- see helper-z-declined.R
  expect_route(es.mcv_r$info_used_crude, "pt_bis_r")
  # z declined on an out-of-range biserial r -- see helper-z-declined.R
  expect_route(es.mcv_mean$info_used_crude, "means_sd")
  expect_equal(es.mcv_r$es_crude, es.mcv_mean$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, es.mcv_mean$se_crude, tolerance = 1e-10)
})


# RPB to SMD (VIECHTBAUER) --------
test_that("rpb to SMD (viecht)", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp))
  dat$reverse_means <- FALSE
  dat$pt_bis_r <- as.numeric(metafor::escalc(
    vtype = "CS", measure = "RPB",
    m1i = mean_exp, m2i = mean_nexp, sd1i = mean_sd_exp, sd2i = mean_sd_nexp,
    n1i = n_exp, n2i = n_nexp, data = dat
  )$yi)
  dat$pt_bis_r_pval <- with(dat, 2 * pt(pt_bis_r * sqrt((n_exp + n_nexp - 2) / (1 - pt_bis_r^2)), n_exp + n_nexp - 2, lower.tail = FALSE))
  dat <- subset(dat, pt_bis_r_pval < 1)


  es.mcv_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "d"), digits = 11)
  es.mcv_mean <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_r$info_used_crude), "pt_bis_r")
  expect_equal(unique(es.mcv_mean$info_used_crude), "means_sd")
  expect_equal(es.mcv_r$es_crude, es.mcv_mean$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, es.mcv_mean$se_crude, tolerance = 1e-10)
})

test_that("rpb pval to SMD (viecht)", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp))
  dat$pt_bis_r <- as.numeric(metafor::escalc(
    vtype = "CS", measure = "RPB",
    m1i = mean_exp, m2i = mean_nexp, sd1i = mean_sd_exp, sd2i = mean_sd_nexp,
    n1i = n_exp, n2i = n_nexp, data = dat
  )$yi)
  dat$pt_bis_r_pval <- with(dat, 2 * pt(pt_bis_r * sqrt((n_exp + n_nexp - 2) / (1 - pt_bis_r^2)), n_exp + n_nexp - 2, lower.tail = FALSE))
  dat <- subset(dat, pt_bis_r_pval < 1)
  dat$reverse_means <- FALSE

  es.mcv_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval", measure = "d"), digits = 11)
  es.mcv_mean <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_r$info_used_crude), "pt_bis_r_pval")
  expect_equal(unique(es.mcv_mean$info_used_crude), "means_sd")
  expect_equal(es.mcv_r$es_crude, es.mcv_mean$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, es.mcv_mean$se_crude, tolerance = 1e-10)
})

### REVERSE ------
test_that("point-bis pval - Reverse", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp))
  dat$pt_bis_r <- as.numeric(metafor::escalc(
    vtype = "CS", measure = "RPB",
    m1i = mean_exp, m2i = mean_nexp, sd1i = mean_sd_exp, sd2i = mean_sd_nexp,
    n1i = n_exp, n2i = n_nexp, data = dat
  )$yi)
  dat$pt_bis_r_pval <- with(dat, 2 * pt(pt_bis_r * sqrt((n_exp + n_nexp - 2) / (1 - pt_bis_r^2)), n_exp + n_nexp - 2, lower.tail = FALSE))
  dat <- subset(dat, pt_bis_r_pval < 1)
  dat$reverse_means <- FALSE

  dat$reverse_pt_bis_r <- FALSE
  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval", measure = "logor"), digits = 11)
  es.mcv_t_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval",
                                    smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval",
                                    smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval",
                                    smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval",
                                    smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  dat$reverse_pt_bis_r_pval <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval", measure = "logor"), digits = 11)
  es.mcv_t_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)
  expect_route(c(es.mcv_t_d$info_used_crude, es.mcv_t_g$info_used_crude,
                    es.mcv_t_or$info_used_crude, es.mcv_t_r1$info_used_crude,
                    es.mcv_t_r2$info_used_crude, es.mcv_t_z1$info_used_crude,
                    es.mcv_t_z2$info_used_crude), "pt_bis_r_pval")

  expect_equal(es.mcv_t_d$es_crude, -es.mcv_t_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_crude, es.mcv_t_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_crude, -es.mcv_t_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_crude, es.mcv_t_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r1$es_crude, -es.mcv_t_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r1$se_crude, es.mcv_t_r1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z1$es_crude, -es.mcv_t_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z1$se_crude, es.mcv_t_z1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r2$es_crude, -es.mcv_t_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$se_crude, es.mcv_t_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z2$es_crude, -es.mcv_t_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$se_crude, es.mcv_t_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_crude, -es.mcv_t_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_crude, es.mcv_t_or_rv$se_crude, tolerance = 1e-10)
})

test_that("point-bis - Reverse", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp))
  dat$pt_bis_r <- as.numeric(metafor::escalc(
    vtype = "CS", measure = "RPB",
    m1i = mean_exp, m2i = mean_nexp, sd1i = mean_sd_exp, sd2i = mean_sd_nexp,
    n1i = n_exp, n2i = n_nexp, data = dat
  )$yi)
  dat$pt_bis_r_pval <- with(dat, 2 * pt(pt_bis_r * sqrt((n_exp + n_nexp - 2) / (1 - pt_bis_r^2)), n_exp + n_nexp - 2, lower.tail = FALSE))
  dat <- subset(dat, pt_bis_r_pval < 1)
  dat$reverse_pt_bis_r <- FALSE

  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "logor"), digits = 11)
  es.mcv_t_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r",
                                    smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r",
                                    smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r",
                                    smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r",
                                    smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  dat$reverse_pt_bis_r <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "logor"), digits = 11)
  es.mcv_t_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)
  expect_route(c(es.mcv_t_d$info_used_crude, es.mcv_t_g$info_used_crude,
                    es.mcv_t_or$info_used_crude, es.mcv_t_r1$info_used_crude,
                    es.mcv_t_r2$info_used_crude, es.mcv_t_z1$info_used_crude,
                    es.mcv_t_z2$info_used_crude), "pt_bis_r")

  expect_equal(es.mcv_t_d$es_crude, -es.mcv_t_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_crude, es.mcv_t_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_crude, -es.mcv_t_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_crude, es.mcv_t_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r1$es_crude, -es.mcv_t_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r1$se_crude, es.mcv_t_r1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z1$es_crude, -es.mcv_t_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z1$se_crude, es.mcv_t_z1_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r2$es_crude, -es.mcv_t_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$se_crude, es.mcv_t_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z2$es_crude, -es.mcv_t_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$se_crude, es.mcv_t_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_crude, -es.mcv_t_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_crude, es.mcv_t_or_rv$se_crude, tolerance = 1e-10)
})
