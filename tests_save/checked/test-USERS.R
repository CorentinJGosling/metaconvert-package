test_that("es user - crude", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  dat$user_es_crude <- rnorm(nrow(dat), 0, 2)
  dat$user_se_crude <- abs(rnorm(nrow(dat), 0, 1))
  dat$user_ci_lo_crude <- dat$user_es_crude - qnorm(.975) * dat$user_se_crude
  dat$user_ci_up_crude <- dat$user_es_crude + qnorm(.975) * dat$user_se_crude
  dat$user_es_measure_crude <- "ABC"

  es.mcv_d <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "user_input_crude",
    measure = "d"
  ), digits = 11)

  expect_equal(unique(es.mcv_d$info_used_crude), "user_input_crude")
  expect_equal(es.mcv_d$measure_crude, rep("ABC", nrow(dat)), tolerance = 1e-10)
  expect_equal(es.mcv_d$es_crude, dat$user_es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, dat$user_se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$es_ci_lo_crude, dat$user_ci_lo_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$es_ci_up_crude, dat$user_ci_up_crude, tolerance = 1e-10)
})
test_that("es user - crude 2", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  dat$user_es_crude <- rnorm(nrow(dat), 0, 2)
  dat$user_se_crude <- abs(rnorm(nrow(dat), 0, 1))
  dat$user_ci_lo_crude <- dat$user_es_crude - qnorm(.975) * dat$user_se_crude
  dat$user_ci_up_crude <- dat$user_es_crude + qnorm(.975) * dat$user_se_crude
  dat$user_es_measure_crude <- "ABC"

  es.mcv_d1 <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "user_input_crude",
    measure = "d"
  ), digits = 11)

  dat$user_se_crude <- NA
  es.mcv_d2 <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "user_input_crude",
    measure = "d"
  ), digits = 11)

  expect_equal(unique(es.mcv_d1$info_used_crude), "user_input_crude")
  expect_equal(unique(es.mcv_d2$info_used_crude), "user_input_crude")
  expect_equal(es.mcv_d1$measure_crude, rep("ABC", nrow(dat)), tolerance = 1e-10)
  expect_equal(es.mcv_d2$measure_crude, rep("ABC", nrow(dat)), tolerance = 1e-10)
  expect_equal(es.mcv_d1$es_crude, es.mcv_d2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d1$se_crude, es.mcv_d2$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d1$es_ci_lo_crude, es.mcv_d2$es_ci_lo_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d1$es_ci_up_crude, es.mcv_d2$es_ci_up_crude, tolerance = 1e-10)
})

test_that("es user - adjusted", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp))

  dat$user_es_adj <- rnorm(nrow(dat), 0, 2)
  dat$user_se_adj <- abs(rnorm(nrow(dat), 0, 1))
  dat$user_ci_lo_adj <- dat$user_es_adj - qnorm(.975) * dat$user_se_adj
  dat$user_ci_up_adj <- dat$user_es_adj + qnorm(.975) * dat$user_se_adj
  dat$user_es_measure_adj <- "ABC"

  es.mcv_d <- summary(convert_df(dat,
                                 verbose = FALSE,
                                 es_selected = "hierarchy", hierarchy = "user_input_adj",
                                 measure = "d"
  ), digits = 11)

  expect_equal(unique(es.mcv_d$info_used_adj), "user_input_adj")
  expect_equal(es.mcv_d$measure_adj, rep("ABC", nrow(dat)), tolerance = 1e-10)
  expect_equal(es.mcv_d$es_adjusted, dat$user_es_adj, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_adjusted, dat$user_se_adj, tolerance = 1e-10)
  expect_equal(es.mcv_d$es_ci_lo_adjusted, dat$user_ci_lo_adj, tolerance = 1e-10)
  expect_equal(es.mcv_d$es_ci_up_adjusted, dat$user_ci_up_adj, tolerance = 1e-10)
})
test_that("es user - adjusted2", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
                  !is.na(n_exp) & !is.na(n_nexp))

  dat$user_es_adj <- rnorm(nrow(dat), 0, 2)
  dat$user_se_adj <- abs(rnorm(nrow(dat), 0, 1))
  dat$user_ci_lo_adj <- dat$user_es_adj - qnorm(.975) * dat$user_se_adj
  dat$user_ci_up_adj <- dat$user_es_adj + qnorm(.975) * dat$user_se_adj
  dat$user_es_measure_adj <- "ABC"

  es.mcv_d1 <- summary(convert_df(dat,
                                 verbose = FALSE,
                                 es_selected = "hierarchy", hierarchy = "user_input_adj",
                                 measure = "d"
  ), digits = 11)

  dat$user_se_adj <- NA
  es.mcv_d2 <- summary(convert_df(dat,
                                  verbose = FALSE,
                                  es_selected = "hierarchy", hierarchy = "user_input_adj",
                                  measure = "d"
  ), digits = 11)

  expect_equal(unique(es.mcv_d1$info_used_adjusted), "user_input_adj")
  expect_equal(unique(es.mcv_d2$info_used_adjusted), "user_input_adj")
  expect_equal(es.mcv_d1$measure_adjusted, rep("ABC", nrow(dat)), tolerance = 1e-10)
  expect_equal(es.mcv_d2$measure_adjusted, rep("ABC", nrow(dat)), tolerance = 1e-10)
  expect_equal(es.mcv_d1$es_adjusted, es.mcv_d2$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d1$se_adjusted, es.mcv_d2$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d1$es_ci_lo_adjusted, es.mcv_d2$es_ci_lo_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d1$es_ci_up_adjusted, es.mcv_d2$es_ci_up_adjusted, tolerance = 1e-10)
})
