# esc >= 0.5.1 on R >= 4.3 no longer accepts vector inputs in esc_t/esc_f
# (their internal `missing(x) || is.na(x)` guards error on length > 1),
# so the comparators are called row by row.
esc_by_row <- function(.fn, ..., es.type) {
  args <- data.frame(...)
  do.call(rbind, lapply(seq_len(nrow(args)), function(i) {
    data.frame(do.call(.fn, c(as.list(args[i, , drop = FALSE]), list(es.type = es.type))))
  }))
}

test_that("pre-test. t-test and d coincides from raw data", {
  scg1 <- rnorm(250, 2, 1)
  scg2 <- rnorm(250, 0, 1)

  dat <- data.frame(
    mean_exp = mean(scg1),
    mean_nexp = mean(scg2),
    mean_sd_exp = sd(scg1),
    mean_sd_nexp = sd(scg2),
    n_exp = length(scg1),
    n_nexp = length(scg2)
  )

  t <- t.test(c(scg1, scg2) ~ rep(c("A", "B"), each = 250), var.equal = TRUE)

  dat$student_t <- t$statistic

  es.mcv_means <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_t <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "d"), digits = 11)
  expect_equal(unique(es.mcv_means$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_t$info_used_crude), "student_t")

  expect_equal(es.mcv_means$es_crude, es.mcv_t$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_means$se_crude, es.mcv_t$se_crude, tolerance = 1e-10)
})
test_that("pre-test. t-test pval and d coincides from raw data", {
  scg1 <- rnorm(250, 2, 1)
  scg2 <- rnorm(250, 0, 1)

  dat <- data.frame(
    mean_exp = mean(scg1),
    mean_nexp = mean(scg2),
    mean_sd_exp = sd(scg1),
    mean_sd_nexp = sd(scg2),
    n_exp = length(scg1),
    n_nexp = length(scg2)
  )

  t <- t.test(c(scg1, scg2) ~ rep(c("A", "B"), each = 250), var.equal = TRUE)

  dat$student_t_pval <- t$p.value

  es.mcv_means <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_t <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "d"), digits = 11)
  expect_equal(unique(es.mcv_means$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_t$info_used_crude), "student_t_pval")
  expect_equal(es.mcv_means$es_crude, es.mcv_t$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_means$se_crude, es.mcv_t$se_crude, tolerance = 1e-10)
})
test_that("pre-test. f-test and d coincides from raw data", {
  scg1 <- rnorm(250, 2, 1)
  scg2 <- rnorm(250, 0, 1)

  dat <- data.frame(
    mean_exp = mean(scg1),
    mean_nexp = mean(scg2),
    mean_sd_exp = sd(scg1),
    mean_sd_nexp = sd(scg2),
    n_exp = length(scg1),
    n_nexp = length(scg2)
  )

  f <- summary(aov(c(scg1, scg2) ~ rep(c("A", "B"), each = 250)))

  dat$anova_f <- f[[1]][1, 4]

  es.mcv_means <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_t <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "d"), digits = 11)
  expect_equal(unique(es.mcv_means$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_t$info_used_crude), "anova_f")

  expect_equal(es.mcv_means$es_crude, es.mcv_t$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_means$se_crude, es.mcv_t$se_crude, tolerance = 1e-10)
})
test_that("pre-test. f-test pval and d coincides from raw data", {
  scg1 <- rnorm(250, 2, 1)
  scg2 <- rnorm(250, 0, 1)

  dat <- data.frame(
    mean_exp = mean(scg1),
    mean_nexp = mean(scg2),
    mean_sd_exp = sd(scg1),
    mean_sd_nexp = sd(scg2),
    n_exp = length(scg1),
    n_nexp = length(scg2)
  )

  f <- summary(aov(c(scg1, scg2) ~ rep(c("A", "B"), each = 250)))

  dat$anova_f_pval <- f[[1]][1, 5]

  es.mcv_means <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_t <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "d"), digits = 11)
  expect_equal(unique(es.mcv_means$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_t$info_used_crude), "anova_f_pval")
  expect_equal(es.mcv_means$es_crude, es.mcv_t$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_means$se_crude, es.mcv_t$se_crude, tolerance = 1e-10)
})

#
test_that("pre-test. rpb and t coincides from raw data", {

  scg1=rnorm(50); scg2=rnorm(50)

  dat <- data.frame(
    mean_exp=mean(scg1), mean_sd_exp=sd(scg1),
    mean_nexp=mean(scg2), mean_sd_nexp=sd(scg2),
    n_exp = length(scg1),
    n_nexp = length(scg2)
  )

  t = t.test(c(scg1, scg2)~rep(c("A", "B"), each=50), var.equal=TRUE)$statistic
  dat$student_t = t

  corr = cor.test(~c(scg1, scg2)+rep(c(1, 0), each = 50))$estimate
  dat$pt_bis_r = corr

  es.mcv_means <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_t <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "d"), digits = 11)
  es.mcv_rpb <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pt_bis_r", measure = "d"), digits = 11)
  expect_equal(unique(es.mcv_means$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_t$info_used_crude), "student_t")
  expect_equal(unique(es.mcv_rpb$info_used_crude), "pt_bis_r")
  expect_equal(es.mcv_means$es_crude, es.mcv_t$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_means$se_crude, es.mcv_t$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_means$es_crude, es.mcv_rpb$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_means$se_crude, es.mcv_rpb$se_crude, tolerance = 1e-10)
})

# T-test to SMD/OR/COR
test_that("1. t-test value", {
  dat <- df.haza
  dat$student_t <- runif(nrow(dat), -3, 3)
  dat <- subset(dat, !is.na(student_t) & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_student_t <- FALSE


  comp_res_d <- esc_by_row(esc::esc_t, t = dat$student_t, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "d")
  comp_res_or <- esc_by_row(esc::esc_t, t = dat$student_t, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "or")
  t_mfr <- metafor::escalc(
    ti = student_t, n1i = n_exp, n2i = n_nexp,
    measure = "SMD", data = dat,
    vtype = "LS2"
  )
  r_mfr <- metafor::escalc(
    ti = student_t, n1i = n_exp, n2i=  n_nexp,
    measure = "RBIS", data = dat
  )
  zbis <- function(rb, n1, n2) {
    n <- n1 + n2
    p <- n1 / n
    fzp <- dnorm(qnorm(p))
    a <- sqrt(fzp) / (p*(1-p))^(1/4)
    zrb <- (a/2) * log((1+a*rb)/(1-a*rb))

    ci_1 <- zrb - qnorm(.975) * sqrt(1/(n-1))
    ci_2 <- zrb + qnorm(.975) * sqrt(1/(n-1))
    cbind(zrb,
          ci_1,
          ci_2)
  }
  res_z = zbis(as.numeric(r_mfr$yi), dat$n_exp, dat$n_nexp)
  res_z[,1][res_z[,1]>1] <- 1

  r_mfr$yi[r_mfr$yi>1] <- 1
  # z_mfr <- metafor::escalc(
  #   ti = student_t, n1i = n_exp, n2i=  n_nexp,
  #   measure = "ZBIS", data = dat
  # )
  # z_mfr$yi[z_mfr$yi>1] <- 1
  ## metaconvert
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "d"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "g"), digits = 11)
  es.mcv_logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "logor"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "or"), digits = 11)
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "r",
    smd_to_cor = "viechtbauer"
  ), digits = 11)
  es.mcv_z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "z",
    smd_to_cor = "viechtbauer"
  ), digits = 11)

  ## t => SMD d
  expect_equal(unique(es.mcv_d$info_used_crude), "student_t")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 1e-10)

  ## t => SMD g
  expect_equal(es.mcv_g$es_crude, as.numeric(t_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude^2, as.numeric(t_mfr$vi), tolerance = 1e-10)

  ## t => OR
  expect_equal(unique(es.mcv_or$info_used_crude), "student_t")
  expect_equal(es.mcv_or$es_crude, comp_res_or$es, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, comp_res_or$se, tolerance = 1e-10)

  ## t => logOR
  expect_equal(unique(es.mcv_logor$info_used_crude), "student_t")
  expect_equal(es.mcv_logor$es_crude, log(comp_res_or$es), tolerance = 1e-10)
  expect_equal(es.mcv_logor$se_crude, comp_res_or$se, tolerance = 1e-10)

  ## t => R
  expect_equal(unique(es.mcv_r$info_used_crude), "student_t")
  expect_equal(es.mcv_r$es_crude, as.numeric(r_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude^2, as.numeric(r_mfr$vi), tolerance = 1e-10)

  ## t => Z : FAILURE, escalc uses different calculations
  # expect_equal(unique(es.mcv_r$info_used_crude), "student_t")
  # expect_equal(es.mcv_z$es_crude, as.numeric(z_mfr$yi), tolerance = Inf)
  # expect_equal(es.mcv_z$se_crude^2, as.numeric(z_mfr$vi), tolerance = Inf)

  ## t => Z : Original formula
  expect_equal(es.mcv_z$es_crude, res_z[,1], tolerance = 1e-10)
  expect_equal(es.mcv_z$es_ci_lo_crude, res_z[,2], tolerance = 1e-10)
  expect_equal(es.mcv_z$es_ci_up_crude, res_z[,3], tolerance = 1e-10)
})
test_that("REVERSE - t-test value", {
  dat <- df.haza
  dat$student_t <- runif(nrow(dat), -3, 3)
  dat <- subset(dat, !is.na(student_t) & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_student_t <- FALSE
  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "logor"), digits = 11)
  es.mcv_t_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t",
                                    smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t",
                                    smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t",
                                    smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t",
                                    smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  dat$reverse_student_t <- TRUE

  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "logor"), digits = 11)
  es.mcv_t_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t",
                                    smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t",
                                    smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t",
                                    smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  expect_true(all(c(es.mcv_t_d$info_used_crude, es.mcv_t_g$info_used_crude,
  es.mcv_t_or$info_used_crude, es.mcv_t_r1$info_used_crude,
  es.mcv_t_r2$info_used_crude, es.mcv_t_z1$info_used_crude,
  es.mcv_t_z2$info_used_crude) == "student_t"))

  expect_equal(es.mcv_t_r1$es_crude, -es.mcv_t_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r1$se_crude, es.mcv_t_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$es_crude, -es.mcv_t_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$se_crude, es.mcv_t_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z1$es_crude, -es.mcv_t_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z1$se_crude, es.mcv_t_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$es_crude, -es.mcv_t_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$se_crude, es.mcv_t_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_d$es_crude, -es.mcv_t_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_crude, es.mcv_t_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_crude, -es.mcv_t_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_crude, es.mcv_t_g_rv$se_crude, tolerance = 1e-10)


  expect_equal(es.mcv_t_or$es_crude, -es.mcv_t_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_crude, es.mcv_t_or_rv$se_crude, tolerance = 1e-10)
})

test_that("2. t-test p-value", {
  dat <- df.haza
  dat$student_t <- runif(nrow(dat), -3, 3)
  dat$student_t_pval <- 2 * pt(dat$student_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)
  dat <- subset(dat, dat$student_t_pval < 1 & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_student_t_pval <- FALSE

  t_mfr <- metafor::escalc(
    pi = student_t_pval, n1i = n_exp, n2i = n_nexp,
    measure = "SMD", data = dat,
    vtype = "LS2"
  )
  r_mfr <- metafor::escalc(
    ti = student_t, n1i = n_exp, n2i=  n_nexp,
    measure = "RBIS", data = dat
  )
  z_mfr <- metafor::escalc(
    ti = student_t, n1i = n_exp, n2i=  n_nexp,
    measure = "ZBIS", data = dat
  )

  comp_res_p <- esc_by_row(esc::esc_t, t = dat$student_t, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "d")
  comp_res_d <- esc_by_row(esc::esc_t, p = dat$student_t_pval, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "d")
  comp_res_or <- esc_by_row(esc::esc_t, p = dat$student_t_pval, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "or")

  ## metaconvert
  es.mcv_t <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t", measure = "d"), digits = 11)
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "d"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "g"), digits = 11)
  es.mcv_logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "logor"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "or"), digits = 11)
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "r",
    smd_to_cor = "viechtbauer"
  ), digits = 11)
  es.mcv_z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "z",
    smd_to_cor = "viechtbauer"
  ), digits = 11)

  # t => SMD d
  expect_equal(unique(es.mcv_d$info_used_crude), "student_t_pval")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 1e-10)

  ## t pval => SMD g
  expect_equal(es.mcv_g$es_crude, as.numeric(t_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude^2, as.numeric(t_mfr$vi), tolerance = 1e-10)

  ## t pval => OR
  expect_equal(unique(es.mcv_or$info_used_crude), "student_t_pval")
  expect_equal(es.mcv_or$es_crude, comp_res_or$es, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, comp_res_or$se, tolerance = 1e-10)

  expect_equal(unique(es.mcv_logor$info_used_crude), "student_t_pval")
  expect_equal(es.mcv_logor$es_crude, log(comp_res_or$es), tolerance = 1e-10)
  expect_equal(es.mcv_logor$se_crude, comp_res_or$se, tolerance = 1e-10)

  ## t => R
  expect_equal(unique(es.mcv_r$info_used_crude), "student_t_pval")
  expect_equal(es.mcv_r$es_crude, as.numeric(r_mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude^2, as.numeric(r_mfr$vi), tolerance = 1e-10)


  ## t => Z : Original formula
  zbis <- function(rb, n1, n2) {
    n <- n1 + n2
    p <- n1 / n
    fzp <- dnorm(qnorm(p))
    a <- sqrt(fzp) / (p*(1-p))^(1/4)
    zrb <- (a/2) * log((1+a*rb)/(1-a*rb))

    ci_1 <- zrb - qnorm(.975) * sqrt(1/(n-1))
    ci_2 <- zrb + qnorm(.975) * sqrt(1/(n-1))
    cbind(zrb,
          ci_1,
          ci_2)
  }
  res_z = zbis(as.numeric(r_mfr$yi), dat$n_exp, dat$n_nexp)
  expect_equal(es.mcv_z$es_crude, res_z[,1], tolerance = 1e-10)
  expect_equal(es.mcv_z$es_ci_lo_crude, res_z[,2], tolerance = 1e-10)
  expect_equal(es.mcv_z$es_ci_up_crude, res_z[,3], tolerance = 1e-10)
})
test_that("REVESE - t-test p-value", {
  dat <- df.haza
  dat$student_t <- runif(nrow(dat), -3, 3)
  dat$student_t_pval <- 2 * pt(dat$student_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)
  dat <- subset(dat, dat$student_t_pval < 1 & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_student_t_pval <- FALSE

  es.mcv_t_p_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "d"), digits = 11)
  es.mcv_t_p_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "g"), digits = 11)
  es.mcv_t_p_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "logor"), digits = 11)
  es.mcv_t_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  dat$reverse_student_t_pval <- TRUE
  es.mcv_t_p_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "d"), digits = 11)
  es.mcv_t_p_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "g"), digits = 11)
  es.mcv_t_p_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval", measure = "logor"), digits = 11)
  es.mcv_t_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "student_t_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  expect_true(all(c(es.mcv_t_p_d$info_used_crude, es.mcv_t_p_g$info_used_crude,
                    es.mcv_t_p_or$info_used_crude, es.mcv_t_r1$info_used_crude,
                    es.mcv_t_r2$info_used_crude, es.mcv_t_z1$info_used_crude,
                    es.mcv_t_z2$info_used_crude) == "student_t_pval"))

  expect_equal(es.mcv_t_r1$es_crude, -es.mcv_t_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r1$se_crude, es.mcv_t_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$es_crude, -es.mcv_t_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$se_crude, es.mcv_t_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z1$es_crude, -es.mcv_t_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z1$se_crude, es.mcv_t_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$es_crude, -es.mcv_t_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$se_crude, es.mcv_t_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_p_d$es_crude, -es.mcv_t_p_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_p_d$se_crude, es.mcv_t_p_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_p_g$es_crude, -es.mcv_t_p_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_p_g$se_crude, es.mcv_t_p_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_p_or$es_crude, -es.mcv_t_p_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_p_or$se_crude, es.mcv_t_p_or_rv$se_crude, tolerance = 1e-10)
})

test_that("3. anova f", {
  dat <- df.haza
  dat$student_t <- runif(nrow(dat), -3, 3)
  dat$anova_f <- dat$student_t^2
  dat <- subset(dat, !is.na(anova_f) & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_anova_f <- FALSE

  comp_res_d <- esc_by_row(esc::esc_f, f = dat$anova_f, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "d")
  comp_res_or <- esc_by_row(esc::esc_f, f = dat$anova_f, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "or")
  comp_res_r <- esc_by_row(esc::esc_f, f = dat$anova_f, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "r")

  ## metaconvert
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "or"), digits = 11)
  es.mcv_logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "logor"), digits = 11)
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "r",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)
  es.mcv_z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "z",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_d$info_used_crude), "anova_f")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 1e-10)

  expect_equal(unique(es.mcv_or$info_used_crude), "anova_f")
  expect_equal(exp(es.mcv_logor$es_crude), comp_res_or$es, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, comp_res_or$se, tolerance = 1e-10)

  expect_equal(unique(es.mcv_logor$info_used_crude), "anova_f")
  expect_equal(es.mcv_logor$es_crude, log(comp_res_or$es), tolerance = 1e-10)
  expect_equal(es.mcv_logor$se_crude, comp_res_or$se, tolerance = 1e-10)

  expect_equal(unique(es.mcv_r$info_used_crude), "anova_f")
  expect_equal(es.mcv_r$es_crude, comp_res_r$es, tolerance = 1e-1)
  expect_equal(es.mcv_r$se_crude, comp_res_r$se, tolerance = 1.5e-1)

  expect_equal(unique(es.mcv_z$info_used_crude), "anova_f")
  expect_equal(es.mcv_z$es_crude, atanh(comp_res_r$es), tolerance = 1e-1)
  expect_equal(es.mcv_z$se_crude, comp_res_r$se, tolerance = 1e-1)
})
test_that("REVERSE - anova f", {
  dat <- df.haza
  dat$student_t <- runif(nrow(dat), -3, 3)
  dat$anova_f <- dat$student_t^2
  dat <- subset(dat, !is.na(anova_f) & !is.na(n_exp) & !is.na(n_nexp))

  dat$reverse_anova_f <- FALSE
  es.mcv_f_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "d"), digits = 11)
  es.mcv_f_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "g"), digits = 11)
  es.mcv_f_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "logor"), digits = 11)
  es.mcv_f_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "r"), digits = 11)
  es.mcv_f_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "z"), digits = 11)
  es.mcv_t_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  dat$reverse_anova_f <- TRUE
  es.mcv_f_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "d"), digits = 11)
  es.mcv_f_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "g"), digits = 11)
  es.mcv_f_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f", measure = "logor"), digits = 11)
  es.mcv_t_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  expect_true(all(c(es.mcv_f_d$info_used_crude, es.mcv_f_g$info_used_crude,
                    es.mcv_f_or$info_used_crude, es.mcv_t_r1$info_used_crude,
                    es.mcv_t_r2$info_used_crude, es.mcv_t_z1$info_used_crude,
                    es.mcv_t_z2$info_used_crude) == "anova_f"))

  expect_equal(es.mcv_t_r1$es_crude, -es.mcv_t_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r1$se_crude, es.mcv_t_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$es_crude, -es.mcv_t_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$se_crude, es.mcv_t_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z1$es_crude, -es.mcv_t_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z1$se_crude, es.mcv_t_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$es_crude, -es.mcv_t_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$se_crude, es.mcv_t_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_d$es_crude, -es.mcv_f_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_d$se_crude, es.mcv_f_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_g$es_crude, -es.mcv_f_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_g$se_crude, es.mcv_f_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_or$es_crude, -es.mcv_f_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_or$se_crude, es.mcv_f_or_rv$se_crude, tolerance = 1e-10)
})

test_that("4. anova p-value", {
  dat <- df.haza
  dat$student_t <- runif(nrow(dat), -3, 3)
  dat$anova_f_pval <- 2 * pt(dat$student_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)
  dat <- subset(dat, anova_f_pval < 1 & !is.na(anova_f_pval) & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_anova_f_pval <- FALSE

  comp_res_d <- esc_by_row(esc::esc_t, p = dat$anova_f_pval, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "d")
  comp_res_or <- esc_by_row(esc::esc_t, p = dat$anova_f_pval, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "or")
  comp_res_r <- esc_by_row(esc::esc_t, p = dat$anova_f_pval, grp1n = dat$n_exp, grp2n = dat$n_nexp, es.type = "r")
  ## metaconvert
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "or"), digits = 11)
  es.mcv_logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "logor"), digits = 11)
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "r",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)
  es.mcv_z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "z",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)

  expect_equal(unique(es.mcv_d$info_used_crude), "anova_f_pval")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 1e-10)

  expect_equal(unique(es.mcv_or$info_used_crude), "anova_f_pval")
  expect_equal(exp(es.mcv_logor$es_crude), comp_res_or$es, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, comp_res_or$se, tolerance = 1e-10)

  expect_equal(unique(es.mcv_logor$info_used_crude), "anova_f_pval")
  expect_equal(es.mcv_logor$es_crude, log(comp_res_or$es), tolerance = 1e-10)
  expect_equal(es.mcv_logor$se_crude, comp_res_or$se, tolerance = 1e-10)

  expect_equal(unique(es.mcv_r$info_used_crude), "anova_f_pval")
  expect_equal(es.mcv_r$es_crude, comp_res_r$es, tolerance = 1e-1)
  expect_equal(es.mcv_r$se_crude, comp_res_r$se, tolerance = 1.5e-1)

  expect_equal(unique(es.mcv_z$info_used_crude), "anova_f_pval")
  expect_equal(es.mcv_z$es_crude, atanh(comp_res_r$es), tolerance = 1e-1)
  expect_equal(es.mcv_z$se_crude, comp_res_r$se, tolerance = 1e-1)
})


test_that("REVERSE - anova pval", {
  dat <- df.haza
  dat$student_t <- runif(nrow(dat), -3, 3)
  dat$anova_f_pval <- 2 * pt(dat$student_t, dat$n_exp + dat$n_nexp - 2, lower.tail = FALSE)
  dat <- subset(dat, anova_f_pval < 1 & !is.na(anova_f_pval) & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_anova_f_pval <- FALSE
  es.mcv_f_p_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "d"), digits = 11)
  es.mcv_f_p_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "g"), digits = 11)
  es.mcv_f_p_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "logor"), digits = 11)
  es.mcv_t_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  dat$reverse_anova_f_pval <- TRUE
  es.mcv_f_p_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "d"), digits = 11)
  es.mcv_f_p_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "g"), digits = 11)
  es.mcv_f_p_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval", measure = "logor"), digits = 11)
  es.mcv_t_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "anova_f_pval",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  expect_true(all(c(es.mcv_f_p_d$info_used_crude, es.mcv_f_p_g$info_used_crude,
                    es.mcv_f_p_or$info_used_crude, es.mcv_t_r1$info_used_crude,
                    es.mcv_t_r2$info_used_crude, es.mcv_t_z1$info_used_crude,
                    es.mcv_t_z2$info_used_crude) == "anova_f_pval"))

  expect_equal(es.mcv_t_r1$es_crude, -es.mcv_t_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r1$se_crude, es.mcv_t_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$es_crude, -es.mcv_t_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$se_crude, es.mcv_t_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z1$es_crude, -es.mcv_t_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z1$se_crude, es.mcv_t_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$es_crude, -es.mcv_t_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$se_crude, es.mcv_t_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_p_d$es_crude, -es.mcv_f_p_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_p_d$se_crude, es.mcv_f_p_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_p_g$es_crude, -es.mcv_f_p_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_p_g$se_crude, es.mcv_f_p_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_f_p_or$es_crude, -es.mcv_f_p_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_f_p_or$se_crude, es.mcv_f_p_or_rv$se_crude, tolerance = 1e-10)
})

test_that("5. etasq", {
  dat <- subset(df.haza, select = -c(n_cases_exp, n_cases_nexp, n_controls_exp, n_controls_nexp))
  dat$etasq <- runif(nrow(dat), 0.01, 0.3)
  dat <- subset(dat, !is.na(etasq) & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_etasq <- FALSE
  dat$n_exp <- dat$n_nexp <- 30
  dat$n_sample <- dat$n_exp + dat$n_nexp

  comp_res_d <- esc::cohens_d(eta = dat$etasq)
  comp_res_or <- esc::log_odds(eta = dat$etasq)
  comp_res_r <- esc::pearsons_r(eta = dat$etasq)

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "logor"), digits = 11)
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "r",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)
  es.mcv_z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "z",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_d$info_used_crude), "etasq")
  expect_equal(unique(es.mcv_or$info_used_crude), "etasq")
  expect_equal(unique(es.mcv_r$info_used_crude), "etasq")
  expect_equal(unique(es.mcv_z$info_used_crude), "etasq")
  # NEWS.md 2.0.1: es_from_etasq() now inverts eta^2 = F / (F + N - 2) exactly rather
  # than using the large-sample, equal-n limit 2*sqrt(eta/(1-eta)) that esc still uses.
  # The two differ by exactly sqrt((N-2)/N) at equal arm sizes (1.7% at n = 30/30) and
  # by up to 39% at 10/90. The exact form is the one that reproduces Cohen's d computed
  # from the raw data, and it agrees with es_from_anova_f() and es_from_pt_bis_r().
  # This test therefore now pins the DELIBERATE divergence from esc.
  N <- 60L
  df_err <- N - 2L
  exact_d <- sqrt(dat$etasq / (1 - dat$etasq) * df_err * (1 / 30 + 1 / 30))
  expect_equal(exact_d, comp_res_d * sqrt(df_err / N), tolerance = 1e-10)

  expect_equal(es.mcv_d$es_crude, exact_d, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, exact_d * pi / sqrt(3), tolerance = 1e-10)
  # r / z on the lipsey_cooper scale, derived from the same exact d
  a_lc <- (30 + 30)^2 / (30 * 30)
  exact_r <- exact_d / sqrt(exact_d^2 + a_lc)
  expect_equal(es.mcv_r$es_crude, exact_r, tolerance = 1e-10)
  expect_equal(es.mcv_z$es_crude, atanh(exact_r), tolerance = 1e-10)
})

test_that("REVERSE - etasq", {
  dat <- df.haza
  dat$etasq <- runif(nrow(dat), 0.01, 0.3)
  dat <- subset(dat, !is.na(etasq) & !is.na(n_exp) & !is.na(n_nexp))
  dat$reverse_etasq <- FALSE

  es.mcv_eta_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "d"), digits = 11)
  es.mcv_eta_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "g"), digits = 11)
  es.mcv_eta_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "logor"), digits = 11)
  es.mcv_t_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  dat$reverse_etasq <- TRUE
  es.mcv_eta_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "d"), digits = 11)
  es.mcv_eta_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "g"), digits = 11)
  es.mcv_eta_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq", measure = "logor"), digits = 11)
  es.mcv_t_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  expect_true(all(c(es.mcv_eta_d$info_used_crude, es.mcv_eta_g$info_used_crude,
                    es.mcv_eta_or$info_used_crude, es.mcv_t_r1$info_used_crude,
                    es.mcv_t_r2$info_used_crude, es.mcv_t_z1$info_used_crude,
                    es.mcv_t_z2$info_used_crude) == "etasq"))

  expect_equal(es.mcv_t_r1$es_crude, -es.mcv_t_r1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r1$se_crude, es.mcv_t_r1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$es_crude, -es.mcv_t_r2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$se_crude, es.mcv_t_r2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z1$es_crude, -es.mcv_t_z1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z1$se_crude, es.mcv_t_z1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$es_crude, -es.mcv_t_z2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$se_crude, es.mcv_t_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_eta_d$es_crude, -es.mcv_eta_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_eta_d$se_crude, es.mcv_eta_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_eta_g$es_crude, -es.mcv_eta_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_eta_g$se_crude, es.mcv_eta_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_eta_or$es_crude, -es.mcv_eta_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_eta_or$se_crude, es.mcv_eta_or_rv$se_crude, tolerance = 1e-10)
})

test_that("REVERSE - etasq adjusted", {
  dat <- df.haza
  dat$etasq_adj <- runif(nrow(dat), 0.01, 0.3)
  dat <- subset(dat, !is.na(etasq_adj) & !is.na(n_exp) & !is.na(n_nexp))
  dat$cov_outcome_r <- 0.3
  dat$n_cov_ancova <- 0.3
  # 2.0.1: the ADJUSTED eta-squared route is flipped by its own reverse_etasq_adj
  # column, not by reverse_etasq (which now governs only the crude route). A dataset
  # that relied on reverse_etasq flipping both must set both -- see NEWS 2.0.1.
  dat$reverse_etasq <- FALSE
  dat$reverse_etasq_adj <- FALSE
  es.mcv_eta_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj", measure = "d"), digits = 11)
  es.mcv_eta_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj", measure = "g"), digits = 11)
  es.mcv_eta_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj", measure = "logor"), digits = 11)
  es.mcv_t_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  dat$reverse_etasq <- TRUE
  dat$reverse_etasq_adj <- TRUE
  es.mcv_eta_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj", measure = "d"), digits = 11)
  es.mcv_eta_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj", measure = "g"), digits = 11)
  es.mcv_eta_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj", measure = "logor"), digits = 11)
  es.mcv_t_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj",
                                       smd_to_cor = "viechtbauer", measure = "r"), digits = 11)
  es.mcv_t_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj",
                                       smd_to_cor = "lipsey_cooper",  measure = "r"), digits = 11)
  es.mcv_t_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj",
                                       smd_to_cor = "viechtbauer", measure = "z"), digits = 11)
  es.mcv_t_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "etasq_adj",
                                       smd_to_cor = "lipsey_cooper",  measure = "z"), digits = 11)

  expect_equal(unique(es.mcv_eta_d$info_used_adjusted), "etasq_adj")
  expect_equal(unique(es.mcv_eta_g$info_used_adjusted), "etasq_adj")
  expect_equal(unique(es.mcv_eta_or$info_used_adjusted), "etasq_adj")

  expect_true(all(c(es.mcv_eta_d$info_used_adjusted, es.mcv_eta_g$info_used_adjusted,
                    es.mcv_eta_or$info_used_adjusted, es.mcv_t_r1$info_used_adjusted,
                    es.mcv_t_r2$info_used_adjusted, es.mcv_t_z1$info_used_adjusted,
                    es.mcv_t_z2$info_used_adjusted) == "etasq_adj"))

  expect_equal(es.mcv_t_r1$es_adjusted, -es.mcv_t_r1_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_r1$se_adjusted, es.mcv_t_r1_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$es_adjusted, -es.mcv_t_r2_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_r2$se_adjusted, es.mcv_t_r2_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_z1$es_adjusted, -es.mcv_t_z1_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_z1$se_adjusted, es.mcv_t_z1_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$es_adjusted, -es.mcv_t_z2_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_z2$se_adjusted, es.mcv_t_z2_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_eta_d$es_adjusted, -es.mcv_eta_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_eta_d$se_adjusted, es.mcv_eta_d_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_eta_g$es_adjusted, -es.mcv_eta_g_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_eta_g$se_adjusted, es.mcv_eta_g_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_eta_or$es_adjusted, -es.mcv_eta_or_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_eta_or$se_adjusted, es.mcv_eta_or_rv$se_adjusted, tolerance = 1e-10)
})
