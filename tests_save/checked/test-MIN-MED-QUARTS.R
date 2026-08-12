### Cohen's d ==============
test_that("min-med-max", {
  df.haza$min_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$min_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$max_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$max_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE

  res_exp <- with(dat, meta:::mean_sd_range(
    n = n_exp, median = med_exp, min = min_exp,
    max = max_exp, method.mean = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_range(
    n = n_nexp, median = med_nexp, min = min_nexp,
    max = max_nexp, method.mean = "Wan"
  ))

  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max",
    measure = "d"
  ), digits = 11)

  expect_equal(unique(es.mcv_d$info_used_crude), "med_min_max")
  expect_equal(es.mcv_d$es_crude, es$d, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$d_se, tolerance = 1e-10)
})

test_that("q1-med-q3", {
  df.haza$q1_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$q1_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$q3_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$q3_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE

  res_exp <- with(dat, meta:::mean_sd_iqr(
    n = n_exp, median = med_exp,
    q1 = q1_exp, q3 = q3_exp, method.mean = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_iqr(
    n = n_nexp, median = med_nexp,
    q1 = q1_nexp, q3 = q3_nexp, method.mean = "Wan"
  ))

  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "d"), digits = 11)
  expect_equal(unique(es.mcv_d$info_used_crude), "med_quarts")
  expect_equal(es.mcv_d$es_crude, es$d, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$d_se, tolerance = 1e-10)
})
test_that("min-q1-med-q3-max ", {
  df.haza$min_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$min_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$max_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$max_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$q1_exp <- round(runif(nrow(df.haza), 11, 19))
  df.haza$q1_nexp <- round(runif(nrow(df.haza), 11, 19))
  df.haza$q3_exp <- round(runif(nrow(df.haza), 26, 29))
  df.haza$q3_nexp <- round(runif(nrow(df.haza), 26, 29))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE
  res_exp <- with(dat, meta:::mean_sd_iqr_range(
    n = n_exp, min = min_exp, max = max_exp,
    median = med_exp,
    q1 = q1_exp, q3 = q3_exp,
    method.mean = "Wan",
    method.sd = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_iqr_range(
    n = n_nexp, min = min_nexp, max = max_nexp,
    median = med_nexp,
    q1 = q1_nexp, q3 = q3_nexp, method.mean = "Wan",
    method.sd = "Wan"
  ))
  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "d"), digits = 11)
  expect_equal(unique(es.mcv_d$info_used_crude), "med_min_max_quarts")
  expect_equal(es.mcv_d$es_crude, es$d, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$d_se, tolerance = 1e-10)
})

### Hedges' g ==============
test_that("min-med-max", {
  df.haza$min_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$min_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$max_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$max_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE

  res_exp <- with(dat, meta:::mean_sd_range(
    n = n_exp, median = med_exp, min = min_exp,
    max = max_exp, method.mean = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_range(
    n = n_nexp, median = med_nexp, min = min_nexp,
    max = max_nexp, method.mean = "Wan"
  ))

  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat,
                                 verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max",
                                 measure = "g"
  ), digits = 11)

  expect_equal(unique(es.mcv_d$info_used_crude), "med_min_max")
  expect_equal(es.mcv_d$es_crude, es$g, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$g_se, tolerance = 1e-10)
})

test_that("q1-med-q3", {
  df.haza$q1_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$q1_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$q3_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$q3_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE

  res_exp <- with(dat, meta:::mean_sd_iqr(
    n = n_exp, median = med_exp,
    q1 = q1_exp, q3 = q3_exp, method.mean = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_iqr(
    n = n_nexp, median = med_nexp,
    q1 = q1_nexp, q3 = q3_nexp, method.mean = "Wan"
  ))

  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "g"), digits = 11)
  expect_equal(unique(es.mcv_d$info_used_crude), "med_quarts")
  expect_equal(es.mcv_d$es_crude, es$g, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$g_se, tolerance = 1e-10)
})
test_that("min-q1-med-q3-max ", {
  df.haza$min_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$min_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$max_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$max_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$q1_exp <- round(runif(nrow(df.haza), 11, 19))
  df.haza$q1_nexp <- round(runif(nrow(df.haza), 11, 19))
  df.haza$q3_exp <- round(runif(nrow(df.haza), 26, 29))
  df.haza$q3_nexp <- round(runif(nrow(df.haza), 26, 29))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE
  res_exp <- with(dat, meta:::mean_sd_iqr_range(
    n = n_exp, min = min_exp, max = max_exp,
    median = med_exp,
    q1 = q1_exp, q3 = q3_exp,
    method.mean = "Wan",
    method.sd = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_iqr_range(
    n = n_nexp, min = min_nexp, max = max_nexp,
    median = med_nexp,
    q1 = q1_nexp, q3 = q3_nexp, method.mean = "Wan",
    method.sd = "Wan"
  ))
  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "g"), digits = 11)
  expect_equal(unique(es.mcv_d$info_used_crude), "med_min_max_quarts")
  expect_equal(es.mcv_d$es_crude, es$g, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$g_se, tolerance = 1e-10)
})

### MD ==============
test_that("min-med-max", {
  df.haza$min_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$min_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$max_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$max_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE

  res_exp <- with(dat, meta:::mean_sd_range(
    n = n_exp, median = med_exp, min = min_exp,
    max = max_exp, method.mean = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_range(
    n = n_nexp, median = med_nexp, min = min_nexp,
    max = max_nexp, method.mean = "Wan"
  ))

  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat,
                                 verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max",
                                 measure = "md"
  ), digits = 11)

  expect_equal(unique(es.mcv_d$info_used_crude), "med_min_max")
  expect_equal(es.mcv_d$es_crude, es$md, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$md_se, tolerance = 1e-10)
})

test_that("q1-med-q3", {
  df.haza$q1_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$q1_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$q3_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$q3_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE

  res_exp <- with(dat, meta:::mean_sd_iqr(
    n = n_exp, median = med_exp,
    q1 = q1_exp, q3 = q3_exp, method.mean = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_iqr(
    n = n_nexp, median = med_nexp,
    q1 = q1_nexp, q3 = q3_nexp, method.mean = "Wan"
  ))

  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "md"), digits = 11)
  expect_equal(unique(es.mcv_d$info_used_crude), "med_quarts")
  expect_equal(es.mcv_d$es_crude, es$md, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$md_se, tolerance = 1e-10)
})

test_that("min-q1-med-q3-max ", {
  df.haza$min_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$min_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$max_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$max_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$q1_exp <- round(runif(nrow(df.haza), 11, 19))
  df.haza$q1_nexp <- round(runif(nrow(df.haza), 11, 19))
  df.haza$q3_exp <- round(runif(nrow(df.haza), 26, 29))
  df.haza$q3_nexp <- round(runif(nrow(df.haza), 26, 29))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE
  res_exp <- with(dat, meta:::mean_sd_iqr_range(
    n = n_exp, min = min_exp, max = max_exp,
    median = med_exp,
    q1 = q1_exp, q3 = q3_exp,
    method.mean = "Wan",
    method.sd = "Wan"
  ))
  res_nexp <- with(dat, meta:::mean_sd_iqr_range(
    n = n_nexp, min = min_nexp, max = max_nexp,
    median = med_nexp,
    q1 = q1_nexp, q3 = q3_nexp, method.mean = "Wan",
    method.sd = "Wan"
  ))
  es <- es_from_means_sd(
    mean_exp = res_exp$mean, mean_nexp = res_nexp$mean,
    mean_sd_exp = res_exp$sd, mean_sd_nexp = res_nexp$sd,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp
  )

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "md"), digits = 11)
  expect_equal(unique(es.mcv_d$info_used_crude), "med_min_max_quarts")
  expect_equal(es.mcv_d$es_crude, es$md, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es$md_se, tolerance = 1e-10)
})

test_that("Reverse", {
  df.haza$min_exp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$min_nexp <- round(runif(nrow(df.haza), 0, 10))
  df.haza$max_exp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$max_nexp <- round(runif(nrow(df.haza), 30, 40))
  df.haza$q1_exp <- round(runif(nrow(df.haza), 11, 19))
  df.haza$q1_nexp <- round(runif(nrow(df.haza), 11, 19))
  df.haza$q3_exp <- round(runif(nrow(df.haza), 26, 29))
  df.haza$q3_nexp <- round(runif(nrow(df.haza), 26, 29))
  df.haza$med_exp <- round(runif(nrow(df.haza), 20, 25))
  df.haza$med_nexp <- round(runif(nrow(df.haza), 20, 25))
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_med <- FALSE
  es.mcv_m_sd_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "md"), digits = 11)
  es.mcv_m_se_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "md"), digits = 11)
  es.mcv_m_ci_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "md"), digits = 11)

  es.mcv_m_sd_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "d"), digits = 11)
  es.mcv_m_se_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "d"), digits = 11)
  es.mcv_m_ci_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "d"), digits = 11)

  es.mcv_m_sd_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "g"), digits = 11)
  es.mcv_m_se_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "g"), digits = 11)
  es.mcv_m_ci_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "g"), digits = 11)

  es.mcv_m_sd_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "r"), digits = 11)
  es.mcv_m_se_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "r"), digits = 11)
  es.mcv_m_ci_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "r"), digits = 11)

  es.mcv_m_sd_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "z"), digits = 11)
  es.mcv_m_se_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "z"), digits = 11)
  es.mcv_m_ci_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "z"), digits = 11)

  es.mcv_m_sd_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "logor"), digits = 11)
  es.mcv_m_se_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "logor"), digits = 11)
  es.mcv_m_ci_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "logor"), digits = 11)

  dat$reverse_med <- TRUE

  es.mcv_m_sd_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "md"), digits = 11)
  es.mcv_m_se_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "md"), digits = 11)
  es.mcv_m_ci_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "md"), digits = 11)

  es.mcv_m_sd_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "d"), digits = 11)
  es.mcv_m_se_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "d"), digits = 11)
  es.mcv_m_ci_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "d"), digits = 11)

  es.mcv_m_sd_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "g"), digits = 11)
  es.mcv_m_se_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "g"), digits = 11)
  es.mcv_m_ci_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "g"), digits = 11)

  es.mcv_m_sd_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "r"), digits = 11)
  es.mcv_m_se_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "r"), digits = 11)
  es.mcv_m_ci_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "r"), digits = 11)

  es.mcv_m_sd_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "z"), digits = 11)
  es.mcv_m_se_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "z"), digits = 11)
  es.mcv_m_ci_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "z"), digits = 11)

  es.mcv_m_sd_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max", measure = "logor"), digits = 11)
  es.mcv_m_se_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_quarts", measure = "logor"), digits = 11)
  es.mcv_m_ci_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "med_min_max_quarts", measure = "logor"), digits = 11)

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

  expect_equal(es.mcv_m_sd_d$info_used_crude, es.mcv_m_sd_d_rv$info_used_crude)
  expect_equal(es.mcv_m_sd_d$es_crude, -es.mcv_m_sd_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$es_crude, -es.mcv_m_se_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$es_crude, -es.mcv_m_ci_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_d$se_crude, es.mcv_m_sd_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$se_crude, es.mcv_m_se_d_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$se_crude, es.mcv_m_ci_d_rv$se_crude, tolerance = 1e-10)

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
