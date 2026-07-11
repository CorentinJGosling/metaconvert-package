# Long-running file (> 25 s): executed locally and in CI (devtools::test and
# testthat set NOT_CRAN=true); skipped wholesale on CRAN to respect check-time
# limits. The fast pre/post files still run on CRAN.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
library(testthat)
library(metaConvert)

# T-test to SMD/OR/COR
test_that("1. paired t-test => D", {
  res <- data.frame(1)

  mean_pre_nexp <- rnorm(40, 80, 15)
  res$mean_pre_sd_nexp <- sd(mean_pre_nexp)
  res$mean_pre_nexp <- mean(mean_pre_nexp)

  mean_nexp <- rnorm(40, 85, 15)
  res$mean_sd_nexp <- sd(mean_nexp)
  res$mean_nexp <- mean(mean_nexp)

  mean_pre_exp <- rnorm(33, 80, 15)
  res$mean_pre_sd_exp <- sd(mean_pre_exp)
  res$mean_pre_exp <- mean(mean_pre_exp)

  mean_exp <- rnorm(33, 95, 15)
  res$mean_sd_exp <- sd(mean_exp)
  res$mean_exp <- mean(mean_exp)

  res$n_exp <- length(mean_pre_exp)
  res$n_nexp <- length(mean_pre_nexp)

  res$r_pre_post_nexp <- cor.test(mean_nexp, mean_pre_nexp)$estimate
  res$r_pre_post_exp <- cor.test(mean_exp, mean_pre_exp)$estimate

  res$paired_t_exp <- t.test(mean_exp, mean_pre_exp, paired = TRUE, alternative = "two.sided")$statistic
  res$paired_t_nexp <- t.test(mean_nexp, mean_pre_nexp, paired = TRUE, alternative = "two.sided")$statistic

  res1 <- summary(convert_df(res, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "d"), digits = 11)
  res2 <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
    pre_post_to_smd = "cooper", measure = "d"
  ), digits = 11)
  expect_equal(unique(res1$info_used_crude), "paired_t")
  expect_equal(unique(res2$info_used_crude), "means_sd_pre_post")

  expect_equal(res1$es_crude, res2$es_crude, tolerance = 1e-10)
  expect_equal(res1$se_crude, res2$se_crude, tolerance = 1e-10)
})


test_that("REVERSE - paired t-test ", {
  dat <- df.haza
  dat$paired_t_exp <- runif(nrow(dat), -3, 3)
  dat$paired_t_nexp <- runif(nrow(dat), -3, 3)
  dat$r_pre_post_exp <- dat$r_pre_post_nexp <- 0.8
  dat <- subset(dat, !is.na(n_exp) & !is.na(n_nexp))

  dat$reverse_paired_t <- FALSE
  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "z"), digits = 11)

  dat$reverse_paired_t <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "z"), digits = 11)

  expect_equal(es.mcv_t_d$es_crude, -es.mcv_t_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_crude, es.mcv_t_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_crude, -es.mcv_t_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_crude, es.mcv_t_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_crude, -es.mcv_t_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_crude, es.mcv_t_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_crude, -es.mcv_t_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_crude, es.mcv_t_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_crude, -es.mcv_t_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_crude, es.mcv_t_or_rv$se_crude, tolerance = 1e-10)
})

test_that("2. paired t-test p-value", {
  res <- data.frame(1)

  mean_pre_nexp <- rnorm(40, 80, 5)
  res$mean_pre_sd_nexp <- sd(mean_pre_nexp)
  res$mean_pre_nexp <- mean(mean_pre_nexp)

  mean_nexp <- rnorm(40, 100, 5)
  res$mean_sd_nexp <- sd(mean_nexp)
  res$mean_nexp <- mean(mean_nexp)

  mean_pre_exp <- rnorm(33, 80, 5)
  res$mean_pre_sd_exp <- sd(mean_pre_exp)
  res$mean_pre_exp <- mean(mean_pre_exp)

  mean_exp <- rnorm(33, 105, 5)
  res$mean_sd_exp <- sd(mean_exp)
  res$mean_exp <- mean(mean_exp)

  res$n_exp <- length(mean_pre_exp)
  res$n_nexp <- length(mean_pre_nexp)

  res$r_pre_post_nexp <- cor.test(mean_nexp, mean_pre_nexp)$estimate
  res$r_pre_post_exp <- cor.test(mean_exp, mean_pre_exp)$estimate


  res$paired_t_pval_exp <- t.test(mean_pre_exp, mean_exp,
    paired = TRUE, alternative = "two.sided"
  )$p.value
  res$paired_t_pval_nexp <- t.test(mean_pre_nexp, mean_nexp,
    paired = TRUE, alternative = "two.sided"
  )$p.value

  res1 <- summary(convert_df(res, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "d"), digits = 11)
  res2 <- summary(convert_df(res,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
    pre_post_to_smd = "cooper", measure = "d"
  ), digits = 11)
  expect_equal(unique(res1$info_used_crude), "paired_t_pval")
  expect_equal(unique(res2$info_used_crude), "means_sd_pre_post")

  expect_equal(abs(res1$es_crude), abs(res2$es_crude), tolerance = 1e-10)
  expect_equal(res1$se_crude, res2$se_crude, tolerance = 1e-10)
})

test_that("2-bis. paired t-test p-value 2", {
  dat <- df.haza
  dat$paired_t_exp <- abs(runif(nrow(dat), -5, 5))
  dat$paired_t_nexp <- abs(runif(nrow(dat), -5, 5))

  dat$paired_t_pval_exp <- 2 * pt(dat$paired_t_exp, dat$n_exp - 1, lower.tail = FALSE)
  dat$paired_t_pval_nexp <- 2 * pt(dat$paired_t_nexp, dat$n_nexp - 1, lower.tail = FALSE)
  dat$r_pre_post_exp <- dat$r_pre_post_nexp <- 0.8
  dat <- subset(dat, !is.na(n_exp) & !is.na(n_nexp))

  res1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "d"), digits = 11)
  res2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "d"), digits = 11)
  expect_equal(unique(res1$info_used_crude), "paired_t_pval")
  expect_equal(unique(res2$info_used_crude), "paired_t")

  expect_equal(res1$es_crude, res2$es_crude, tolerance = 1e-10)
  expect_equal(res1$se_crude, res2$se_crude, tolerance = 1e-10)
})

test_that("REVERSE - paired t-test p-value", {
  dat <- df.haza
  # p-values must be valid two-sided probabilities in (0, 1); direction is
  # carried by reverse_paired_t_pval. (Negative t encoded as a one-tailed
  # complement p > 1 is an invalid p-value: the centralized input validation
  # flags it and sets it to NA, making the row fall back to other methods
  # that do not respond to the reverse flag.)
  dat$paired_t_exp <- abs(runif(nrow(dat), -3, 3))
  dat$paired_t_nexp <- abs(runif(nrow(dat), -3, 3))
  dat$paired_t_pval_exp <- 2 * pt(dat$paired_t_exp, dat$n_exp - 1, lower.tail = FALSE)
  dat$paired_t_pval_nexp <- 2 * pt(dat$paired_t_nexp, dat$n_nexp - 1, lower.tail = FALSE)
  dat <- subset(dat, paired_t_pval_exp != 1 & paired_t_pval_nexp != 1 & !is.na(n_exp) & !is.na(n_nexp))

  dat$reverse_paired_t_pval <- FALSE
  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "z"), digits = 11)

  dat$reverse_paired_t_pval <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t_pval", measure = "z"), digits = 11)

  expect_equal(es.mcv_t_d$es_crude, -es.mcv_t_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_crude, es.mcv_t_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_crude, -es.mcv_t_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_crude, es.mcv_t_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_crude, -es.mcv_t_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_crude, es.mcv_t_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_crude, -es.mcv_t_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_crude, es.mcv_t_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_crude, -es.mcv_t_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_crude, es.mcv_t_or_rv$se_crude, tolerance = 1e-10)
})


test_that("3. paired f", {
  dat <- df.haza
  dat$paired_t_exp <- abs(runif(nrow(dat), -5, 5))
  dat$paired_t_nexp <- abs(runif(nrow(dat), -5, 5))

  dat$paired_f_exp <- dat$paired_t_exp^2
  dat$paired_f_nexp <- dat$paired_t_nexp^2
  dat$r_pre_post_exp <- dat$r_pre_post_nexp <- 0.8
  dat <- subset(dat, !is.na(n_exp) & !is.na(n_nexp))

  res1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "d"), digits = 11)
  res2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "d"), digits = 11)
  expect_equal(unique(res1$info_used_crude), "paired_f")
  expect_equal(unique(res2$info_used_crude), "paired_t")

  expect_equal(res1$es_crude, res2$es_crude, tolerance = 1e-10)
  expect_equal(res1$se_crude, res2$se_crude, tolerance = 1e-10)
})


test_that("REVERSE - paired f value", {
  dat <- df.haza
  dat$paired_t_exp <- abs(runif(nrow(dat), -5, 5))
  dat$paired_t_nexp <- abs(runif(nrow(dat), -5, 5))

  dat$paired_f_exp <- dat$paired_t_exp^2
  dat$paired_f_nexp <- dat$paired_t_nexp^2
  dat$r_pre_post_exp <- dat$r_pre_post_nexp <- 0.8
  dat <- subset(dat, !is.na(n_exp) & !is.na(n_nexp))

  dat$reverse_paired_f <- FALSE
  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "z"), digits = 11)

  dat$reverse_paired_f <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f", measure = "z"), digits = 11)

  expect_equal(es.mcv_t_d$es_crude, -es.mcv_t_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_crude, es.mcv_t_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_crude, -es.mcv_t_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_crude, es.mcv_t_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_crude, -es.mcv_t_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_crude, es.mcv_t_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_crude, -es.mcv_t_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_crude, es.mcv_t_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_crude, -es.mcv_t_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_crude, es.mcv_t_or_rv$se_crude, tolerance = 1e-10)
})

test_that("4. paired f p-value", {
  dat <- df.haza
  dat$paired_t_exp <- abs(runif(nrow(dat), -5, 5))
  dat$paired_t_nexp <- abs(runif(nrow(dat), -5, 5))

  dat$paired_f_pval_exp <- 2 * pt(dat$paired_t_exp, dat$n_exp - 1, lower.tail = FALSE)
  dat$paired_f_pval_nexp <- 2 * pt(dat$paired_t_nexp, dat$n_nexp - 1, lower.tail = FALSE)
  dat$r_pre_post_exp <- dat$r_pre_post_nexp <- 0.8
  dat <- subset(dat, !is.na(n_exp) & !is.na(n_nexp))

  res1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "d"), digits = 11)
  res2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_t", measure = "d"), digits = 11)
  expect_equal(unique(res1$info_used_crude), "paired_f_pval")
  expect_equal(unique(res2$info_used_crude), "paired_t")

  expect_equal(res1$es_crude, res2$es_crude, tolerance = 1e-10)
  expect_equal(res1$se_crude, res2$se_crude, tolerance = 1e-10)
})


test_that("REVERSE - anova pval", {
  dat <- df.haza
  dat$paired_t_exp <- abs(runif(nrow(dat), -5, 5))
  dat$paired_t_nexp <- abs(runif(nrow(dat), -5, 5))

  dat$paired_f_pval_exp <- 2 * pt(dat$paired_t_exp, dat$n_exp - 1, lower.tail = FALSE)
  dat$paired_f_pval_nexp <- 2 * pt(dat$paired_t_nexp, dat$n_nexp - 1, lower.tail = FALSE)
  dat$r_pre_post_exp <- dat$r_pre_post_nexp <- 0.8
  dat <- subset(dat, !is.na(n_exp) & !is.na(n_nexp))

  dat$reverse_paired_f_pval <- FALSE
  es.mcv_f_p_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "d"), digits = 11)
  es.mcv_f_p_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "g"), digits = 11)
  es.mcv_f_p_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "logor"), digits = 11)
  es.mcv_f_p_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "r"), digits = 11)
  es.mcv_f_p_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "z"), digits = 11)

  dat$reverse_paired_f_pval <- TRUE
  es.mcv_f_p_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "d"), digits = 11)
  es.mcv_f_p_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "g"), digits = 11)
  es.mcv_f_p_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "logor"), digits = 11)
  es.mcv_f_p_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "r"), digits = 11)
  es.mcv_f_p_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "paired_f_pval", measure = "z"), digits = 11)

  expect_equal(es.mcv_f_p_d$es_crude, -es.mcv_f_p_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_p_d$se_crude, es.mcv_f_p_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_p_g$es_crude, -es.mcv_f_p_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_p_g$se_crude, es.mcv_f_p_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_p_r$es_crude, -es.mcv_f_p_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_p_r$se_crude, es.mcv_f_p_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_p_z$es_crude, -es.mcv_f_p_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_p_z$se_crude, es.mcv_f_p_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_p_or$es_crude, -es.mcv_f_p_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_p_or$se_crude, es.mcv_f_p_or_rv$se_crude, tolerance = 1e-10)
})

} # end NOT_CRAN gate
