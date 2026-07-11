# Long-running file (> 25 s): executed locally and in CI (devtools::test and
# testthat set NOT_CRAN=true); skipped wholesale on CRAN to respect check-time
# limits. The fast pre/post files still run on CRAN.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
# esc >= 0.5.1 on R >= 4.3 no longer accepts vector inputs
# (its internal `missing(x) || is.na(x)` guards error on length > 1),
# so the esc comparators are called row by row.
esc_by_row <- function(.fn, ..., es.type) {
  args <- data.frame(...)
  do.call(rbind, lapply(seq_len(nrow(args)), function(i) {
    data.frame(do.call(.fn, c(as.list(args[i, , drop = FALSE]), list(es.type = es.type))))
  }))
}

# MEANS TO SMD ----------------------------------------------------
test_that("raw_means + SD => d", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))

  dat$reverse_means <- FALSE
  es.esc <- esc_by_row(esc::esc_mean_sd,
    grp1m = dat$mean_exp, grp2m = dat$mean_nexp,
    grp1sd = dat$mean_sd_exp, grp2sd = dat$mean_sd_nexp,
    grp1n = dat$n_exp, grp2n = dat$n_nexp,
    es.type = "d"
  )
  es.mcv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  expect_equal(unique(es.mcv$info_used_crude), "means_sd")
  expect_equal(es.mcv$es_crude, es.esc$es, tolerance = 1e-10)
  expect_equal(es.mcv$se_crude, es.esc$se, tolerance = 1e-10)
})

test_that("raw_means + SE => d", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))

  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$reverse_means <- FALSE
  es.esc <- with(dat, esc_by_row(esc::esc_mean_se,
    grp1m = mean_exp, grp2m = mean_nexp,
    grp1se = mean_se_exp, grp2se = mean_se_nexp,
    grp1n = n_exp, grp2n = n_nexp,
    es.type = "d"
  ))
  es.esc <- data.frame(es.esc)
  es.mcv_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_sd$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_se$info_used_crude), "means_se")
  expect_equal(es.mcv_se$es_crude, es.esc$es, tolerance = 1e-1)
  expect_equal(es.mcv_se$se_crude, es.esc$se, tolerance = 1e-2)
  expect_equal(es.mcv_se$es_crude, es.mcv_sd$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_se$se_crude, es.mcv_sd$se_crude, tolerance = 1e-10)
})

test_that("raw_means + SD => g", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE

  es.mfr <- metafor::escalc(
    m1i = dat$mean_exp, sd1i = dat$mean_sd_exp, n1i = dat$n_exp,
    m2i = dat$mean_nexp, sd2i = dat$mean_sd_nexp, n2i = dat$n_nexp,
    measure = "SMD", vtype = "LS2"
  )
  es.mcv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "g"), digits = 11)
  expect_equal(unique(es.mcv$info_used_crude), "means_sd")
  expect_equal(es.mcv$es_crude, as.numeric(as.character(es.mfr$yi)), tolerance = 1e-10)
  expect_equal(es.mcv$se_crude^2, as.numeric(as.character(es.mfr$vi)), tolerance = 1e-10)
})
# MEANS TO MD ----------------------------------------------------
test_that("raw_means + SD => md", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE
  es.mfr <- metafor::escalc(
    m1i = dat$mean_exp, sd1i = dat$mean_sd_exp, n1i = dat$n_exp,
    m2i = dat$mean_nexp, sd2i = dat$mean_sd_nexp, n2i = dat$n_nexp,
    measure = "MD"
  )
  es.mcv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "md"), digits = 11)
  expect_equal(unique(es.mcv$info_used_crude), "means_sd")
  expect_equal(es.mcv$es_crude, as.numeric(as.character(es.mfr$yi)), tolerance = 1e-10)
  expect_equal(es.mcv$se_crude^2, as.numeric(as.character(es.mfr$vi)), tolerance = 1e-10)
})

# MEANS TO OR ----------------------------------------------------
test_that("raw_means + SD => d + OR", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE

  comp_res_d <- esc_by_row(esc::esc_mean_sd,
    grp1m = dat$mean_exp, grp1sd = dat$mean_sd_exp, grp1n = dat$n_exp,
    grp2m = dat$mean_nexp, grp2sd = dat$mean_sd_nexp, grp2n = dat$n_nexp,
    es.type = "d"
  )
  comp_res_or <- esc_by_row(esc::esc_mean_sd,
    grp1m = dat$mean_exp,
    grp1sd = dat$mean_sd_exp,
    grp1n = dat$n_exp,
    grp2m = dat$mean_nexp,
    grp2sd = dat$mean_sd_nexp,
    grp2n = dat$n_nexp,
    es.type = "or"
  )
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "logor"), digits = 11)

  es.mfr <- metafor::escalc(
    m1i = dat$mean_exp, sd1i = dat$mean_sd_exp, n1i = dat$n_exp,
    m2i = dat$mean_nexp, sd2i = dat$mean_sd_nexp, n2i = dat$n_nexp,
    measure = "D2ORL"
  )
  ## test ES
  expect_equal(unique(es.mcv_d$info_used_crude), "means_sd")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 1e-10)

  expect_equal(unique(es.mcv_or$info_used_crude), "means_sd")
  expect_equal(es.mcv_or$es_crude, log(comp_res_or$es), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, comp_res_or$se, tolerance = 1e-10)

  expect_equal(es.mcv_or$es_crude, as.numeric(as.character(es.mfr$yi)), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude^2, as.numeric(as.character(es.mfr$vi)), tolerance = 1e-10)
})

# MEANS TO COR ----------------------------------------------------

test_that("raw_means + SD => R (LIPSEY-COOPER)", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE

  es.esc <- with(dat, esc_by_row(esc::esc_mean_sd,
    grp1m = mean_exp, grp2m = mean_nexp,
    grp1sd = mean_sd_exp, grp2sd = mean_sd_nexp,
    grp1n = n_exp, grp2n = n_nexp,
    es.type = "r"
  ))
  es.esc <- data.frame(es.esc)

  comp_es <- with(dat, compute.es::mes(
    m.1 = mean_exp, m.2 = mean_nexp, sd.1 = mean_sd_exp, sd.2 = mean_sd_nexp,
    n.1 = n_exp, n.2 = n_nexp,
    dig = 12, verbose = FALSE
  ))
  es.mcv <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "r",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)

  expect_equal(unique(es.mcv$info_used_crude), "means_sd")
  expect_equal(es.mcv$es_crude, es.esc$es, tolerance = 1e-10)
  expect_equal(es.mcv$se_crude^2, comp_es$var.r, tolerance = 1e-10)
})
test_that("raw_means + SD => Z (LIPSEY-COOPER)", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE
  es.esc <- with(dat, esc_by_row(esc::esc_mean_sd,
    grp1m = mean_exp, grp2m = mean_nexp,
    grp1sd = mean_sd_exp, grp2sd = mean_sd_nexp,
    grp1n = n_exp, grp2n = n_nexp,
    es.type = "r"
  ))

  es.esc <- data.frame(es.esc)
  comp_es <- with(dat, compute.es::mes(
    m.1 = mean_exp, m.2 = mean_nexp, sd.1 = mean_sd_exp, sd.2 = mean_sd_nexp, n.1 = n_exp, n.2 = n_nexp,
    dig = 12, verbose = FALSE
  ))
  es.mcv_z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "z",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)

  expect_equal(unique(es.mcv_z$info_used_crude), "means_sd")

  expect_equal(unique(es.mcv_z$info_used_crude), "means_sd")
  expect_equal(es.mcv_z$es_crude, es.esc$fishers.z, tolerance = 1e-10)
  expect_equal(es.mcv_z$es_ci_lo_crude, es.esc$ci.lo.z, tolerance = 1e-10)
  expect_equal(es.mcv_z$es_ci_up_crude, es.esc$ci.hi.z, tolerance = 1e-10)
})


test_that("raw_means + SD => R (VIECHTBAUER)", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE
  mfr <- metafor::escalc(
    measure = "RBIS",
    m1i = mean_exp, m2i = mean_nexp,
    sd1i = mean_sd_exp, sd2i = mean_sd_nexp,
    n1i = n_exp, n2i = n_nexp,
    data = dat, digits = 16
  )
  # mfr_r <- ifelse(as.numeric(mfr$yi) > 1, 1, ifelse(as.numeric(mfr$yi) < -1, -1, as.numeric(mfr$yi)))

  es.mcv_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", smd_to_cor = "viechtbauer", measure = "r"), digits = 16)
  expect_equal(unique(es.mcv_r$info_used_crude), "means_sd")
  expect_equal(es.mcv_r$es_crude, as.numeric(mfr$yi), tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude^2, mfr$vi, tolerance = 1e-10)
})

test_that("raw_means + SD => Z (VIECHTBAUER)", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE
  mfr_rb <- metafor::escalc(
    measure = "RBIS",
    m1i = mean_exp, m2i = mean_nexp,
    sd1i = mean_sd_exp, sd2i = mean_sd_nexp,
    n1i = n_exp, n2i = n_nexp,
    data = dat, digits = 13
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

  es.mcv_z <- summary(
    convert_df(dat,
      verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd",
      smd_to_cor = "viechtbauer", measure = "z"
    ),
    digits = 13
  )
  expect_equal(unique(es.mcv_z$info_used_crude), "means_sd")
  expect_equal(es.mcv_z$es_crude, res_z[,1], tolerance = 1e-10)
  expect_equal(es.mcv_z$es_ci_lo_crude, res_z[,2], tolerance = 1e-10)
  expect_equal(es.mcv_z$es_ci_up_crude, res_z[,3], tolerance = 1e-10)
})


# INTERNAL CALCULATIONS -------
test_that("means_sd = means_sd_pooled", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE
  dat$mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * mean_sd_exp^2 + (n_nexp - 1) * mean_sd_nexp^2) / (n_exp + n_nexp - 2)))

  es.mcv_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_sd_pooled <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_sd$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_sd_pooled$info_used_crude), "means_sd_pooled")
  expect_equal(es.mcv_sd$es_crude, es.mcv_sd_pooled$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_sd$se_crude, es.mcv_sd_pooled$se_crude, tolerance = 1e-10)
})

test_that("means_sd = means_se", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE
  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)


  es.mcv_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_sd$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_se$info_used_crude), "means_se")
  expect_equal(es.mcv_sd$es_crude, es.mcv_se$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_sd$se_crude, es.mcv_se$se_crude, tolerance = 1e-10)
})


test_that("means_sd = means_ci", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE
  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$mean_ci_lo_exp <- dat$mean_exp - dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$mean_ci_up_exp <- dat$mean_exp + dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_up_nexp <- dat$mean_nexp + dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  es.mcv_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_sd$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_ci$info_used_crude), "means_ci")
  expect_equal(es.mcv_sd$es_crude, es.mcv_ci$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_sd$se_crude, es.mcv_ci$se_crude, tolerance = 1e-10)
})

test_that("means_sd = means_plot", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means_plot <- FALSE
  dat$reverse_means <- FALSE

  dat$plot_mean_exp <- dat$mean_exp
  dat$plot_mean_nexp <- dat$mean_nexp
  dat$plot_mean_sd_lo_exp <- dat$mean_exp - dat$mean_sd_exp
  dat$plot_mean_sd_up_exp <- dat$mean_exp + dat$mean_sd_exp
  dat$plot_mean_sd_lo_nexp <- dat$mean_nexp - dat$mean_sd_nexp
  dat$plot_mean_sd_up_nexp <- dat$mean_nexp + dat$mean_sd_nexp

  # View(dat[14, c("plot_mean_exp", "plot_mean_nexp", "plot_mean_se_lo_exp",
  #                "plot_mean_se_up_exp", "plot_mean_se_lo_nexp", "plot_mean_se_up_nexp")])
  #
  # x= dat[14, ]
  es.mcv_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_plot <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_sd$info_used_crude), "means_sd")
  expect_equal(unique(es.mcv_plot$info_used_crude), "means_plot")
  expect_equal(es.mcv_sd$es_crude, es.mcv_plot$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_sd$se_crude, es.mcv_plot$se_crude, tolerance = 1e-10)
})
test_that("means_se = means_plot", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means_plot <- FALSE
  dat$reverse_means <- FALSE
  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$plot_mean_sd_lo_exp <- NA
  dat$plot_mean_sd_up_exp <- NA
  dat$plot_mean_sd_lo_nexp <- NA
  dat$plot_mean_sd_up_nexp <- NA
  dat$plot_mean_exp <- dat$mean_exp
  dat$plot_mean_nexp <- dat$mean_nexp
  dat$plot_mean_se_lo_exp <- dat$mean_exp - dat$mean_se_exp
  dat$plot_mean_se_up_exp <- dat$mean_exp + dat$mean_se_exp
  dat$plot_mean_se_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp
  dat$plot_mean_se_up_nexp <- dat$mean_nexp + dat$mean_se_nexp

  es.mcv_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "d"), digits = 12)
  es.mcv_plot <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "d"), digits = 12)

  expect_equal(unique(es.mcv_se$info_used_crude), "means_se")
  expect_equal(unique(es.mcv_plot$info_used_crude), "means_plot")
  expect_equal(es.mcv_se$es_crude, es.mcv_plot$es_crude, tolerance = 1e-8)
  expect_equal(es.mcv_se$se_crude, es.mcv_plot$se_crude, tolerance = 1e-8)
})
test_that("means_ci = means_plot", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_means <- FALSE
  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$mean_ci_lo_exp <- dat$mean_exp - dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$mean_ci_up_exp <- dat$mean_exp + dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_up_nexp <- dat$mean_nexp + dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$plot_mean_sd_lo_exp <- NA
  dat$plot_mean_sd_up_exp <- NA
  dat$plot_mean_sd_lo_nexp <- NA
  dat$plot_mean_sd_up_nexp <- NA
  dat$plot_mean_se_lo_exp <- NA
  dat$plot_mean_se_up_exp <- NA
  dat$plot_mean_se_lo_nexp <- NA
  dat$plot_mean_se_up_nexp <- NA
  dat$plot_mean_exp <- dat$mean_exp
  dat$plot_mean_nexp <- dat$mean_nexp
  dat$plot_mean_ci_lo_exp <- dat$mean_ci_lo_exp
  dat$plot_mean_ci_up_exp <- dat$mean_ci_up_exp
  dat$plot_mean_ci_lo_nexp <- dat$mean_ci_lo_nexp
  dat$plot_mean_ci_up_nexp <- dat$mean_ci_up_nexp

  es.mcv_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "d"), digits = 11)
  es.mcv_plot <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_ci$info_used_crude), "means_ci")
  expect_equal(unique(es.mcv_plot$info_used_crude), "means_plot")
  expect_equal(es.mcv_ci$es_crude, es.mcv_plot$es_crude, tolerance = 1e-8)
  expect_equal(es.mcv_ci$se_crude, es.mcv_plot$se_crude, tolerance = 1e-8)
})

# REVERSE ----------------------
test_that("means - Reverse standard", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  # dat = dat[1:50, ]
  dat$reverse_means <- FALSE
  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$mean_ci_lo_exp <- dat$mean_exp - dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$mean_ci_up_exp <- dat$mean_exp + dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_up_nexp <- dat$mean_nexp + dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  es.mcv_m_sd_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "md"), digits = 11)
  es.mcv_m_se_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "md"), digits = 11)
  es.mcv_m_ci_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "md"), digits = 11)
  es.mcv_m_pooled_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "md"), digits = 11)

  es.mcv_m_sd_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_m_se_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "d"), digits = 11)
  es.mcv_m_ci_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "d"), digits = 11)
  es.mcv_m_pooled_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "d"), digits = 11)

  es.mcv_m_sd_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "g"), digits = 11)
  es.mcv_m_se_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "g"), digits = 11)
  es.mcv_m_ci_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "g"), digits = 11)
  es.mcv_m_pooled_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "g"), digits = 11)

  es.mcv_m_sd_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "r"), digits = 11)
  es.mcv_m_se_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "r"), digits = 11)
  es.mcv_m_ci_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "r"), digits = 11)
  es.mcv_m_pooled_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "r"), digits = 11)

  es.mcv_m_sd_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "z"), digits = 11)
  es.mcv_m_se_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "z"), digits = 11)
  es.mcv_m_ci_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "z"), digits = 11)
  es.mcv_m_pooled_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "z"), digits = 11)

  es.mcv_m_sd_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "logor"), digits = 11)
  es.mcv_m_se_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "logor"), digits = 11)
  es.mcv_m_ci_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "logor"), digits = 11)
  es.mcv_m_pooled_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "logor"), digits = 11)

  dat$reverse_means <- TRUE

  es.mcv_m_sd_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "md"), digits = 11)
  es.mcv_m_se_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "md"), digits = 11)
  es.mcv_m_ci_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "md"), digits = 11)
  es.mcv_m_pooled_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "md"), digits = 11)

  es.mcv_m_sd_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "d"), digits = 11)
  es.mcv_m_se_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "d"), digits = 11)
  es.mcv_m_ci_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "d"), digits = 11)
  es.mcv_m_pooled_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "d"), digits = 11)

  es.mcv_m_sd_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "g"), digits = 11)
  es.mcv_m_se_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "g"), digits = 11)
  es.mcv_m_ci_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "g"), digits = 11)
  es.mcv_m_pooled_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "g"), digits = 11)

  es.mcv_m_sd_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "r"), digits = 11)
  es.mcv_m_se_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "r"), digits = 11)
  es.mcv_m_ci_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "r"), digits = 11)
  es.mcv_m_pooled_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "r"), digits = 11)

  es.mcv_m_sd_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "z"), digits = 11)
  es.mcv_m_se_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "z"), digits = 11)
  es.mcv_m_ci_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "z"), digits = 11)
  es.mcv_m_pooled_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "z"), digits = 11)

  es.mcv_m_sd_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd", measure = "logor"), digits = 11)
  es.mcv_m_se_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_se", measure = "logor"), digits = 11)
  es.mcv_m_ci_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_ci", measure = "logor"), digits = 11)
  es.mcv_m_pooled_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_sd_pooled", measure = "logor"), digits = 11)

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

test_that("MD - Reverse pval", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$reverse_plot_means <- FALSE
  dat$plot_mean_exp <- dat$mean_exp
  dat$plot_mean_nexp <- dat$mean_nexp
  dat$plot_mean_sd_lo_exp <- dat$mean_exp - dat$mean_sd_exp
  dat$plot_mean_sd_up_exp <- dat$mean_exp + dat$mean_sd_exp
  dat$plot_mean_sd_lo_nexp <- dat$mean_nexp - dat$mean_sd_nexp
  dat$plot_mean_sd_up_nexp <- dat$mean_nexp + dat$mean_sd_nexp
  dat <- dat[1:20, ]

  es.mcv_m_plot_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "logor"), digits = 11)
  es.mcv_m_plot_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "z"), digits = 11)
  es.mcv_m_plot_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "r"), digits = 11)
  es.mcv_m_plot_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "g"), digits = 11)
  es.mcv_m_plot_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "d"), digits = 11)
  es.mcv_m_plot_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "md"), digits = 11)

  dat$reverse_plot_means <- TRUE
  es.mcv_m_plot_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "md"), digits = 11)
  es.mcv_m_plot_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "g"), digits = 11)
  es.mcv_m_plot_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "r"), digits = 11)
  es.mcv_m_plot_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "z"), digits = 11)
  es.mcv_m_plot_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "logor"), digits = 11)
  es.mcv_m_plot_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "means_plot", measure = "d"), digits = 11)

  expect_equal(es.mcv_m_plot_d$info_used_crude, es.mcv_m_plot_d_rv$info_used_crude)
  expect_equal(es.mcv_m_plot_d$es_crude, -es.mcv_m_plot_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_d$se_crude, es.mcv_m_plot_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_md$info_used_crude, es.mcv_m_plot_md_rv$info_used_crude)
  expect_equal(es.mcv_m_plot_md$es_crude, -es.mcv_m_plot_md_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_md$se_crude, es.mcv_m_plot_md_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_g$info_used_crude, es.mcv_m_plot_g_rv$info_used_crude)
  expect_equal(es.mcv_m_plot_g$es_crude, -es.mcv_m_plot_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_g$se_crude, es.mcv_m_plot_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_r$info_used_crude, es.mcv_m_plot_r_rv$info_used_crude)
  expect_equal(es.mcv_m_plot_r$es_crude, -es.mcv_m_plot_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_r$se_crude, es.mcv_m_plot_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_z$info_used_crude, es.mcv_m_plot_z_rv$info_used_crude)
  expect_equal(es.mcv_m_plot_z$es_crude, -es.mcv_m_plot_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_z$se_crude, es.mcv_m_plot_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_or$info_used_crude, es.mcv_m_plot_or_rv$info_used_crude)
  expect_equal(es.mcv_m_plot_or$es_crude, -es.mcv_m_plot_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_or$se_crude, es.mcv_m_plot_or_rv$se_crude, tolerance = 1e-10)
})

} # end NOT_CRAN gate
