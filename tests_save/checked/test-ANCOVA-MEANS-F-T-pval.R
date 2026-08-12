## ANCOVA MEANS ====
test_that("adjusted means + pooled sd adjusted", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_sd_pooled <- with(
    dat,
    sqrt(((n_exp - 1) * ancova_mean_sd_exp^2 +
      (n_nexp - 1) * ancova_mean_sd_nexp^2) /
      (n_exp + n_nexp - 2))
  )
  dat$reverse_ancova_means <- FALSE

  comp_es <- compute.es::a.mes(
    m.1.adj = dat$ancova_mean_exp, m.2.adj = dat$ancova_mean_nexp,
    sd.adj = dat$ancova_mean_sd_pooled,
    n.1 = dat$n_exp, n.2 = dat$n_nexp,
    R = dat$cov_outcome_r, q = dat$n_cov_ancova,
    level = 95, dig = 12, verbose = FALSE
  )

  es.mcv.d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "d"), digits = 11)
  es.mcv.or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "or"), digits = 11)
  es.mcv.logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "logor"), digits = 11)
  es.mcv.r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj",
    smd_to_cor = "lipsey_cooper", measure = "r"
  ), digits = 11)
  es.mcv.z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj",
    smd_to_cor = "lipsey_cooper", measure = "z"
  ), digits = 11)

  ## test ES
  ## SMD
  expect_equal(unique(es.mcv.d$info_used_adjusted), "ancova_means_sd_pooled_adj")
  expect_equal(es.mcv.d$es_adjusted, comp_es$d, tolerance = 1e-10)
  expect_equal(es.mcv.d$se_adjusted, sqrt(comp_es$var.d), tolerance = 1e-10)
  # logOR
  expect_equal(unique(es.mcv.logor$info_used_adjusted), "ancova_means_sd_pooled_adj")
  expect_equal(es.mcv.logor$es_adjusted, comp_es$lOR, tolerance = 1e-10)
  expect_equal(es.mcv.logor$se_adjusted, (comp_es$u.lor - comp_es$l.lor) / (2 * qnorm(0.975)), tolerance = 1e-10)
  # or
  expect_equal(unique(es.mcv.or$info_used_adjusted), "ancova_means_sd_pooled_adj")
  expect_equal(exp(comp_es$lOR), es.mcv.or$es_adjusted, tolerance = 1e-10)
  expect_equal(exp(comp_es$l.lor), es.mcv.or$es_ci_lo_adjusted, tolerance = 1e-10)
  expect_equal(exp(comp_es$u.lor), es.mcv.or$es_ci_up_adjusted, tolerance = 1e-10)
  # R
  expect_equal(unique(es.mcv.r$info_used_adjusted), "ancova_means_sd_pooled_adj")
  expect_equal(es.mcv.r$es_adjusted, comp_es$r, tolerance = 1e-10)
  expect_equal(es.mcv.r$se_adjusted^2, comp_es$var.r, tolerance = 1e-10)
  # Z
  expect_equal(unique(es.mcv.z$info_used_adjusted), "ancova_means_sd_pooled_adj")
  expect_equal(es.mcv.z$es_adjusted, comp_es$fisher.z, tolerance = 1e-10)
  expect_equal(es.mcv.z$se_adjusted^2, comp_es$var.z, tolerance = 1e-1)
})

test_that("adjusted means + sd pooled crude", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$mean_sd_pooled <- with(
    dat,
    sqrt(((n_exp - 1) * mean_sd_exp^2 +
      (n_nexp - 1) * mean_sd_nexp^2) /
      (n_exp + n_nexp - 2))
  )
  dat$reverse_ancova_means <- FALSE

  comp_es <- compute.es::a.mes2(
    m.1.adj = dat$ancova_mean_exp, m.2.adj = dat$ancova_mean_nexp,
    s.pooled = dat$mean_sd_pooled,
    n.1 = dat$n_exp, n.2 = dat$n_nexp,
    R = dat$cov_outcome_r, q = dat$n_cov_ancova,
    level = 95, cer = 0.2, dig = 12, verbose = FALSE
  )

  ## metaconvert
  es.mcv.d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "d"), digits = 12)
  es.mcv.or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "or"), digits = 12)
  es.mcv.logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "logor"), digits = 12)
  es.mcv.r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled",
    smd_to_cor = "lipsey_cooper", measure = "r"
  ), digits = 12)
  es.mcv.z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled",
    smd_to_cor = "lipsey_cooper", measure = "z"
  ), digits = 12)

  ## test ES
  ## SMD
  expect_equal(unique(es.mcv.d$info_used_adjusted), "ancova_means_sd_pooled")
  expect_equal(es.mcv.d$es_adjusted, comp_es$d, tolerance = 1e-10)
  expect_equal(es.mcv.d$se_adjusted, sqrt(comp_es$var.d), tolerance = 1e-10)
  # logOR
  expect_equal(unique(es.mcv.logor$info_used_adjusted), "ancova_means_sd_pooled")
  expect_equal(es.mcv.logor$es_adjusted, comp_es$lOR, tolerance = 1e-10)
  expect_equal(es.mcv.logor$se_adjusted, (comp_es$u.lor - comp_es$l.lor) / (2 * qnorm(0.975)), tolerance = 1e-10)
  # or
  expect_equal(unique(es.mcv.or$info_used_adjusted), "ancova_means_sd_pooled")
  expect_equal(exp(comp_es$lOR), es.mcv.or$es_adjusted, tolerance = 1e-10)
  expect_equal(exp(comp_es$l.lor), es.mcv.or$es_ci_lo_adjusted, tolerance = 1e-10)
  expect_equal(exp(comp_es$u.lor), es.mcv.or$es_ci_up_adjusted, tolerance = 1e-10)
  # R
  expect_equal(unique(es.mcv.r$info_used_adjusted), "ancova_means_sd_pooled")
  expect_equal(es.mcv.r$es_adjusted, comp_es$r, tolerance = 1e-10)
  expect_equal(es.mcv.r$se_adjusted^2, comp_es$var.r, tolerance = 1e-10)
  # Z
  expect_equal(unique(es.mcv.z$info_used_adjusted), "ancova_means_sd_pooled")
  expect_equal(es.mcv.z$es_adjusted, comp_es$fisher.z, tolerance = 1e-10)
  expect_equal(es.mcv.z$se_adjusted^2, comp_es$var.z, tolerance = 1e-1)
})

test_that("adjusted means + sd/se/ci => SMD", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$ancova_mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$ancova_mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  df_model <- dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova
  dat$ancova_mean_ci_lo_exp <- dat$mean_exp - dat$ancova_mean_se_exp * qt(0.975, df_model)
  dat$ancova_mean_ci_up_exp <- dat$mean_exp + dat$ancova_mean_se_exp * qt(0.975, df_model)
  dat$ancova_mean_ci_lo_nexp <- dat$mean_nexp - dat$ancova_mean_se_nexp * qt(0.975, df_model)
  dat$ancova_mean_ci_up_nexp <- dat$mean_nexp + dat$ancova_mean_se_nexp * qt(0.975, df_model)
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * ancova_mean_sd_exp^2 + (n_nexp - 1) * ancova_mean_sd_nexp^2) / (n_exp + n_nexp - 2)))

  dat$reverse_ancova_means <- FALSE
  ## metaconvert
  es.mcv_d_sd_comb <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "d"), digits = 11)
  es.mcv_d_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "d"), digits = 11)
  es.mcv_d_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "d"), digits = 11)
  es.mcv_d_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "d"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_d_sd_comb$info_used_adjusted), "ancova_means_sd_pooled_adj")
  expect_equal(unique(es.mcv_d_sd$info_used_adjusted), "ancova_means_sd")
  expect_equal(unique(es.mcv_d_se$info_used_adjusted), "ancova_means_se")
  expect_equal(unique(es.mcv_d_ci$info_used_adjusted), "ancova_means_ci")

  expect_equal(es.mcv_d_sd_comb$es_adjusted, es.mcv_d_sd$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d_sd_comb$se_adjusted, es.mcv_d_sd$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_d_sd_comb$es_adjusted, es.mcv_d_se$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d_sd_comb$se_adjusted, es.mcv_d_se$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_d_sd_comb$es_adjusted, es.mcv_d_ci$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d_sd_comb$se_adjusted, es.mcv_d_ci$se_adjusted, tolerance = 1e-10)
})
test_that("adjusted means + sd/se/ci => OR", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$ancova_mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$ancova_mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  df_model <- dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova
  dat$ancova_mean_ci_lo_exp <- dat$mean_exp - dat$ancova_mean_se_exp * qt(0.975, df_model)
  dat$ancova_mean_ci_up_exp <- dat$mean_exp + dat$ancova_mean_se_exp * qt(0.975, df_model)
  dat$ancova_mean_ci_lo_nexp <- dat$mean_nexp - dat$ancova_mean_se_nexp * qt(0.975, df_model)
  dat$ancova_mean_ci_up_nexp <- dat$mean_nexp + dat$ancova_mean_se_nexp * qt(0.975, df_model)
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * ancova_mean_sd_exp^2 + (n_nexp - 1) * ancova_mean_sd_nexp^2) / (n_exp + n_nexp - 2)))

  dat$reverse_ancova_means <- FALSE
  ## metaconvert
  es.mcv_d_sd_comb <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "logor"), digits = 11)
  es.mcv_d_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "logor"), digits = 11)
  es.mcv_d_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "logor"), digits = 11)
  es.mcv_d_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "logor"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_d_sd_comb$info_used_adjusted), "ancova_means_sd_pooled_adj")
  expect_equal(unique(es.mcv_d_sd$info_used_adjusted), "ancova_means_sd")
  expect_equal(unique(es.mcv_d_se$info_used_adjusted), "ancova_means_se")
  expect_equal(unique(es.mcv_d_ci$info_used_adjusted), "ancova_means_ci")

  expect_equal(es.mcv_d_sd_comb$es_adjusted, es.mcv_d_sd$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d_sd_comb$se_adjusted, es.mcv_d_sd$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_d_sd_comb$es_adjusted, es.mcv_d_se$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d_sd_comb$se_adjusted, es.mcv_d_se$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_d_sd_comb$es_adjusted, es.mcv_d_ci$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_d_sd_comb$se_adjusted, es.mcv_d_ci$se_adjusted, tolerance = 1e-10)
})

test_that("means + indiv/pooled SD => MD", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$mean_sd_pooled <- with(
    dat,
    sqrt(((n_exp - 1) * mean_sd_exp^2 +
            (n_nexp - 1) * mean_sd_nexp^2) /
           (n_exp + n_nexp - 2))
  )
  dat$reverse_ancova_means <- FALSE


  ## metaconvert
  es.mcv.pool <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled",
                                    measure = "md"), digits = 12)
  es.mcv.indiv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd",
                                     measure = "md"), digits = 12)

  expect_equal(unique(es.mcv.indiv$info_used_adjusted), "ancova_means_sd")
  expect_equal(unique(es.mcv.pool$info_used_adjusted), "ancova_means_sd_pooled")
  expect_equal(es.mcv.indiv$es_adjusted, es.mcv.pool$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv.indiv$se_adjusted, es.mcv.pool$se_adjusted, tolerance = 5e-1)
})
test_that("means + indiv/pooled SD => MD", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) & !is.na(n_exp) &
                  !is.na(mean_nexp) & !is.na(mean_sd_nexp) & !is.na(n_nexp))
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_sd_pooled <- with(
    dat,
    sqrt(((n_exp - 1) * mean_sd_exp^2 +
            (n_nexp - 1) * mean_sd_nexp^2) /
           (n_exp + n_nexp - 2))
  )
  dat$reverse_ancova_means <- FALSE


  ## metaconvert
  es.mcv.pool <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "md"), digits = 12)
  es.mcv.indiv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "md"), digits = 12)

  expect_equal(unique(es.mcv.indiv$info_used_adjusted), "ancova_means_sd")
  expect_equal(unique(es.mcv.pool$info_used_adjusted), "ancova_means_sd_pooled_adj")
  expect_equal(es.mcv.indiv$es_adjusted, es.mcv.pool$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv.indiv$se_adjusted, es.mcv.pool$se_adjusted, tolerance = 1e-1)
})

## ANCOVA statistics ====
test_that("ANCOVA t", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$reverse_ancova_t <- FALSE
  comp_es <- compute.es::a.tes(
    t = dat$ancova_t, n.1 = dat$n_exp, n.2 = dat$n_nexp,
    R = dat$cov_outcome_r,
    q = dat$n_cov_ancova,
    level = 95, cer = 0.2, dig = 12, verbose = FALSE
  )

  ## metaconvert
  es.mcv.d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "d"), digits = 12)
  es.mcv.or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "or"), digits = 12)
  es.mcv.logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "logor"), digits = 12)
  es.mcv.r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t",
    smd_to_cor = "lipsey_cooper", measure = "r"
  ), digits = 12)
  es.mcv.z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t",
    smd_to_cor = "lipsey_cooper", measure = "z"
  ), digits = 12)

  ## test ES
  ## SMD
  expect_equal(unique(es.mcv.d$info_used_adjusted), "ancova_t")
  expect_equal(es.mcv.d$es_adjusted, comp_es$d, tolerance = 1e-10)
  expect_equal(es.mcv.d$se_adjusted, sqrt(comp_es$var.d), tolerance = 1e-10)
  # logOR
  expect_equal(unique(es.mcv.logor$info_used_adjusted), "ancova_t")
  expect_equal(es.mcv.logor$es_adjusted, comp_es$lOR, tolerance = 1e-10)
  expect_equal(es.mcv.logor$se_adjusted, (comp_es$u.lor - comp_es$l.lor) / (2 * qnorm(0.975)), tolerance = 1e-10)
  # or
  expect_equal(unique(es.mcv.or$info_used_adjusted), "ancova_t")
  expect_equal(exp(comp_es$lOR), es.mcv.or$es_adjusted, tolerance = 1e-10)
  expect_equal(exp(comp_es$l.lor), es.mcv.or$es_ci_lo_adjusted, tolerance = 1e-10)
  expect_equal(exp(comp_es$u.lor), es.mcv.or$es_ci_up_adjusted, tolerance = 1e-10)
  # R
  expect_equal(unique(es.mcv.r$info_used_adjusted), "ancova_t")
  expect_equal(es.mcv.r$es_adjusted, comp_es$r, tolerance = 1e-10)
  expect_equal(es.mcv.r$se_adjusted^2, comp_es$var.r, tolerance = 1e-10)
  # Z
  expect_equal(unique(es.mcv.z$info_used_adjusted), "ancova_t")
  expect_equal(es.mcv.z$es_adjusted, comp_es$fisher.z, tolerance = 1e-10)
  expect_equal(es.mcv.z$se_adjusted^2, comp_es$var.z, tolerance = 1e-1)
})
test_that("ANCOVA f", {
  dat <- subset(df.haza, !is.na(mean_sd_exp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$reverse_ancova_f <- FALSE

  comp_es <- compute.es::a.fes(
    f = dat$ancova_f, n.1 = dat$n_exp, n.2 = dat$n_nexp,
    R = dat$cov_outcome_r,
    q = dat$n_cov_ancova,
    level = 95, cer = 0.2, dig = 12, verbose = FALSE
  )

  ## metaconvert
  es.mcv.d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "d"), digits = 12)
  es.mcv.or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "or"), digits = 12)
  es.mcv.logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "logor"), digits = 12)
  es.mcv.r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f",
    smd_to_cor = "lipsey_cooper", measure = "r"
  ), digits = 12)
  es.mcv.z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f",
    smd_to_cor = "lipsey_cooper", measure = "z"
  ), digits = 12)

  ## test ES
  ## SMD
  expect_equal(unique(es.mcv.d$info_used_adjusted), "ancova_f")
  expect_equal(es.mcv.d$es_adjusted, comp_es$d, tolerance = 1e-10)
  expect_equal(es.mcv.d$se_adjusted, sqrt(comp_es$var.d), tolerance = 1e-10)
  # logOR
  expect_equal(unique(es.mcv.logor$info_used_adjusted), "ancova_f")
  expect_equal(es.mcv.logor$es_adjusted, comp_es$lOR, tolerance = 1e-10)
  expect_equal(es.mcv.logor$se_adjusted, (comp_es$u.lor - comp_es$l.lor) / (2 * qnorm(0.975)), tolerance = 1e-10)
  # or
  expect_equal(unique(es.mcv.or$info_used_adjusted), "ancova_f")
  expect_equal(exp(comp_es$lOR), es.mcv.or$es_adjusted, tolerance = 1e-10)
  expect_equal(exp(comp_es$l.lor), es.mcv.or$es_ci_lo_adjusted, tolerance = 1e-10)
  expect_equal(exp(comp_es$u.lor), es.mcv.or$es_ci_up_adjusted, tolerance = 1e-10)
  # R
  expect_equal(unique(es.mcv.r$info_used_adjusted), "ancova_f")
  expect_equal(es.mcv.r$es_adjusted, comp_es$r, tolerance = 1e-10)
  expect_equal(es.mcv.r$se_adjusted^2, comp_es$var.r, tolerance = 1e-10)
  # Z
  expect_equal(unique(es.mcv.z$info_used_adjusted), "ancova_f")
  expect_equal(es.mcv.z$es_adjusted, comp_es$fisher.z, tolerance = 1e-10)
  expect_equal(es.mcv.z$se_adjusted^2, comp_es$var.z, tolerance = 1e-1)
})

test_that("ANCOVA t_pval", {
  df.haza$ancova_t <- runif(nrow(df.haza), -5, 5)
  df.haza$ancova_f <- df.haza$ancova_t^2
  df.haza$ancova_f_pval <- 2 * pt(abs(df.haza$ancova_t), df.haza$n_exp + df.haza$n_nexp - 2 - df.haza$n_cov_ancova, lower.tail = FALSE)
  df.haza$ancova_t_pval <- df.haza$ancova_f_pval
  df.haza$n_cov_ancova <- round(runif(nrow(df.haza), 1, 5))
  df.haza$cov_outcome_r <- 0.3
  df.haza$reverse_ancova_f_pval <- FALSE
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp) & !is.na(ancova_f_pval))

  comp_es <- compute.es::a.pes(
    p = dat$ancova_t_pval, n.1 = dat$n_exp, n.2 = dat$n_nexp,
    R = dat$cov_outcome_r,
    q = dat$n_cov_ancova,
    level = 95, cer = 0.2, dig = 12, verbose = FALSE
  )

  ## metaconvert
  es.mcv.d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "d"), digits = 12)
  es.mcv.or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "or"), digits = 12)
  es.mcv.logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "logor"), digits = 12)
  es.mcv.r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval",
    smd_to_cor = "lipsey_cooper", measure = "r"
  ), digits = 12)
  es.mcv.z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval",
    smd_to_cor = "lipsey_cooper", measure = "z"
  ), digits = 12)

  ## test ES
  ## SMD
  expect_equal(unique(es.mcv.d$info_used_adjusted), "ancova_t_pval")
  expect_equal(es.mcv.d$es_adjusted, abs(comp_es$d), tolerance = 1e-10)
  expect_equal(es.mcv.d$se_adjusted, sqrt(comp_es$var.d), tolerance = 1e-10)
  # logOR
  expect_equal(unique(es.mcv.logor$info_used_adjusted), "ancova_t_pval")
  expect_equal(es.mcv.logor$es_adjusted, abs(comp_es$lOR), tolerance = 1e-10)
  expect_equal(es.mcv.logor$se_adjusted, (comp_es$u.lor - comp_es$l.lor) / (2 * qnorm(0.975)), tolerance = 1e-10)
  # or
  expect_equal(unique(es.mcv.or$info_used_adjusted), "ancova_t_pval")
  expect_equal(log(es.mcv.or$es_adjusted), abs(comp_es$lOR), tolerance = 1e-10)
  expect_equal(es.mcv.or$se_adjusted, (comp_es$u.lor - comp_es$l.lor) / (2 * qnorm(0.975)), tolerance = 1e-10)
  # R
  expect_equal(unique(es.mcv.r$info_used_adjusted), "ancova_t_pval")
  expect_equal(es.mcv.r$es_adjusted, abs(comp_es$r), tolerance = 1e-10)
  expect_equal(es.mcv.r$se_adjusted^2, comp_es$var.r, tolerance = 1e-10)
  # Z
  expect_equal(unique(es.mcv.z$info_used_adjusted), "ancova_t_pval")
  expect_equal(es.mcv.z$es_adjusted, abs(comp_es$fisher.z), tolerance = 1e-10)
  expect_equal(es.mcv.z$se_adjusted^2, comp_es$var.z, tolerance = 1e-1)
})
test_that("ANCOVA f_pval", {
  df.haza$ancova_t <- runif(nrow(df.haza), -5, 5)
  df.haza$ancova_f <- df.haza$ancova_t^2
  df.haza$ancova_f_pval <- 2 * pt(abs(df.haza$ancova_t), df.haza$n_exp + df.haza$n_nexp - 2 - df.haza$n_cov_ancova, lower.tail = FALSE)
  df.haza$ancova_t_pval <- df.haza$ancova_f_pval
  df.haza$n_cov_ancova <- round(runif(nrow(df.haza), 1, 5))
  df.haza$cov_outcome_r <- 0.3
  df.haza$reverse_ancova_f_pval <- FALSE
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp) & !is.na(ancova_f_pval))

  comp_es <- compute.es::a.pes(
    p = dat$ancova_f_pval, n.1 = dat$n_exp, n.2 = dat$n_nexp,
    R = dat$cov_outcome_r,
    q = dat$n_cov_ancova,
    level = 95, cer = 0.2, dig = 12, verbose = FALSE
  )

  ## metaconvert
  es.mcv.d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "d"), digits = 12)
  es.mcv.or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "or"), digits = 12)
  es.mcv.logor <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "logor"), digits = 12)
  es.mcv.r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval",
    smd_to_cor = "lipsey_cooper", measure = "r"
  ), digits = 12)
  es.mcv.z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval",
    smd_to_cor = "lipsey_cooper", measure = "z"
  ), digits = 12)

  ## test ES
  ## SMD
  expect_equal(unique(es.mcv.d$info_used_adjusted), "ancova_f_pval")
  expect_equal(es.mcv.d$es_adjusted, abs(comp_es$d), tolerance = 1e-10)
  expect_equal(es.mcv.d$se_adjusted, sqrt(comp_es$var.d), tolerance = 1e-10)
  # logOR
  expect_equal(unique(es.mcv.logor$info_used_adjusted), "ancova_f_pval")
  expect_equal(es.mcv.logor$es_adjusted, abs(comp_es$lOR), tolerance = 1e-10)
  expect_equal(es.mcv.logor$se_adjusted, (comp_es$u.lor - comp_es$l.lor) / (2 * qnorm(0.975)), tolerance = 1e-10)
  # or
  expect_equal(unique(es.mcv.or$info_used_adjusted), "ancova_f_pval")
  expect_equal(log(es.mcv.or$es_adjusted), abs(comp_es$lOR), tolerance = 1e-10)
  expect_equal(es.mcv.or$se_adjusted, (comp_es$u.lor - comp_es$l.lor) / (2 * qnorm(0.975)), tolerance = 1e-10)
  # R
  expect_equal(unique(es.mcv.r$info_used_adjusted), "ancova_f_pval")
  expect_equal(es.mcv.r$es_adjusted, abs(comp_es$r), tolerance = 1e-10)
  expect_equal(es.mcv.r$se_adjusted^2, comp_es$var.r, tolerance = 1e-10)
  # Z
  expect_equal(unique(es.mcv.z$info_used_adjusted), "ancova_f_pval")
  expect_equal(es.mcv.z$es_adjusted, abs(comp_es$fisher.z), tolerance = 1e-10)
  expect_equal(es.mcv.z$se_adjusted^2, comp_es$var.z, tolerance = 1e-1)
})

## PLOT ====

test_that("adjusted plot means + plot SD => d", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_sd_exp))
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$ancova_mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp - 1)
  dat$ancova_mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp - 1)
  dat$ancova_mean_ci_lo_exp <- dat$mean_exp - dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_up_exp <- dat$mean_exp + dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_lo_nexp <- dat$mean_nexp - dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$ancova_mean_ci_up_nexp <- dat$mean_nexp + dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$plot_ancova_mean_exp <- dat$ancova_mean_exp
  dat$plot_ancova_mean_nexp <- dat$ancova_mean_nexp
  dat$plot_ancova_mean_sd_lo_exp <- dat$ancova_mean_exp - dat$ancova_mean_sd_exp
  dat$plot_ancova_mean_sd_up_exp <- dat$ancova_mean_exp + dat$ancova_mean_sd_exp
  dat$plot_ancova_mean_sd_lo_nexp <- dat$ancova_mean_nexp - dat$ancova_mean_sd_nexp
  dat$plot_ancova_mean_sd_up_nexp <- dat$ancova_mean_nexp + dat$ancova_mean_sd_nexp
  dat <- subset(dat, !is.na(plot_ancova_mean_sd_lo_exp) & !is.na(plot_ancova_mean_sd_up_exp) &
    !is.na(plot_ancova_mean_sd_lo_nexp) & !is.na(plot_ancova_mean_sd_up_nexp))
  dat$reverse_ancova_means <- FALSE
  dat$reverse_plot_ancova_means <- FALSE

  es.mcv_sd <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "d"), digits = 11)
  es.mcv_plot <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_sd$info_used_adjusted), "ancova_means_sd")
  expect_equal(unique(es.mcv_plot$info_used_adjusted), "ancova_means_plot")
  expect_equal(es.mcv_sd$es_adjusted, es.mcv_plot$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_sd$se_adjusted, es.mcv_plot$se_adjusted, tolerance = 1e-10)
})
test_that("adjusted plot means + plot SE => d", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_sd_exp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$ancova_mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$ancova_mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$ancova_mean_ci_lo_exp <- dat$mean_exp - dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_up_exp <- dat$mean_exp + dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_lo_nexp <- dat$mean_nexp - dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$ancova_mean_ci_up_nexp <- dat$mean_nexp + dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$ancova_mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * ancova_mean_sd_exp^2 + (n_nexp - 1) * ancova_mean_sd_nexp^2) / (n_exp + n_nexp - 2)))
  dat$mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * mean_sd_exp^2 + (n_nexp - 1) * mean_sd_nexp^2) / (n_exp + n_nexp - 2)))
  dat <- subset(dat, !is.na(ancova_mean_sd_pooled) & ancova_pval != 1)
  dat$plot_ancova_mean_sd_lo_exp <- NA
  dat$plot_ancova_mean_sd_up_exp <- NA
  dat$plot_ancova_mean_sd_lo_nexp <- NA
  dat$plot_ancova_mean_sd_up_nexp <- NA
  dat$plot_ancova_mean_exp <- dat$ancova_mean_exp
  dat$plot_ancova_mean_nexp <- dat$ancova_mean_nexp
  dat$plot_ancova_mean_se_lo_exp <- dat$ancova_mean_exp - dat$ancova_mean_se_exp
  dat$plot_ancova_mean_se_up_exp <- dat$ancova_mean_exp + dat$ancova_mean_se_exp
  dat$plot_ancova_mean_se_lo_nexp <- dat$ancova_mean_nexp - dat$ancova_mean_se_nexp
  dat$plot_ancova_mean_se_up_nexp <- dat$ancova_mean_nexp + dat$ancova_mean_se_nexp
  dat$reverse_ancova_means <- FALSE
  dat$reverse_plot_ancova_means <- FALSE

  es.mcv_se <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "d"), digits = 11)
  es.mcv_plot <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_se$info_used_adjusted), "ancova_means_se")
  expect_equal(unique(es.mcv_plot$info_used_adjusted), "ancova_means_plot")
  expect_equal(es.mcv_se$es_adjusted, es.mcv_plot$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_se$se_adjusted, es.mcv_plot$se_adjusted, tolerance = 1e-10)
})

test_that("adjusted plot means + plot CI => d", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_sd_exp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$ancova_mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$ancova_mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$ancova_mean_ci_lo_exp <- dat$mean_exp - dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_up_exp <- dat$mean_exp + dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_lo_nexp <- dat$mean_nexp - dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$ancova_mean_ci_up_nexp <- dat$mean_nexp + dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat <- subset(dat, !is.na(ancova_mean_ci_lo_exp) & !is.na(ancova_mean_ci_up_exp) & !is.na(ancova_mean_ci_lo_nexp) & !is.na(ancova_mean_ci_up_nexp))
  dat$plot_ancova_mean_sd_lo_exp <- NA
  dat$plot_ancova_mean_sd_up_exp <- NA
  dat$plot_ancova_mean_sd_lo_nexp <- NA
  dat$plot_ancova_mean_sd_up_nexp <- NA
  dat$plot_ancova_mean_se_lo_exp <- NA
  dat$plot_ancova_mean_se_up_exp <- NA
  dat$plot_ancova_mean_se_lo_nexp <- NA
  dat$plot_ancova_mean_se_up_nexp <- NA
  dat$plot_ancova_mean_exp <- dat$ancova_mean_exp
  dat$plot_ancova_mean_nexp <- dat$ancova_mean_nexp
  dat$plot_ancova_mean_ci_lo_exp <- dat$ancova_mean_ci_lo_exp
  dat$plot_ancova_mean_ci_up_exp <- dat$ancova_mean_ci_up_exp
  dat$plot_ancova_mean_ci_lo_nexp <- dat$ancova_mean_ci_lo_nexp
  dat$plot_ancova_mean_ci_up_nexp <- dat$ancova_mean_ci_up_nexp
  dat$reverse_ancova_means <- FALSE
  dat$reverse_plot_ancova_means <- FALSE

  es.mcv_ci <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "d"), digits = 11)
  es.mcv_plot <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "d"), digits = 11)

  expect_equal(unique(es.mcv_ci$info_used_adjusted), "ancova_means_ci")
  expect_equal(unique(es.mcv_plot$info_used_adjusted), "ancova_means_plot")
  expect_equal(es.mcv_ci$es_adjusted, es.mcv_plot$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_ci$se_adjusted, es.mcv_plot$se_adjusted, tolerance = 1e-10)
})

## REVERSE ====
test_that("ancova f - Reverse", {
  dat <- subset(df.haza, !is.na(mean_sd_exp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$reverse_ancova_f <- FALSE

  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "z"), digits = 11)

  dat$reverse_ancova_f <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f", measure = "z"), digits = 11)

  expect_equal(unique(es.mcv_t_d$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_g$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_or$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_r$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_z$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_d_rv$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_g_rv$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_or_rv$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_r_rv$info_used_adjusted), "ancova_f")
  expect_equal(unique(es.mcv_t_z_rv$info_used_adjusted), "ancova_f")

  expect_equal(es.mcv_t_d$es_adjusted, -es.mcv_t_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_adjusted, es.mcv_t_d_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_adjusted, -es.mcv_t_g_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_adjusted, es.mcv_t_g_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_adjusted, -es.mcv_t_r_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_adjusted, es.mcv_t_r_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_adjusted, -es.mcv_t_z_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_adjusted, es.mcv_t_z_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_adjusted, -es.mcv_t_or_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_adjusted, es.mcv_t_or_rv$se_adjusted, tolerance = 1e-10)
})
test_that("ancova t - Reverse", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$reverse_ancova_t <- FALSE

  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "z"), digits = 11)

  dat$reverse_ancova_t <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t", measure = "z"), digits = 11)

  expect_equal(unique(es.mcv_t_d$info_used_adjusted), unique(es.mcv_t_d_rv$info_used_adjusted))
  expect_equal(unique(es.mcv_t_g$info_used_adjusted), unique(es.mcv_t_g_rv$info_used_adjusted))
  expect_equal(unique(es.mcv_t_or$info_used_adjusted), unique(es.mcv_t_or_rv$info_used_adjusted))
  expect_equal(unique(es.mcv_t_r$info_used_adjusted), unique(es.mcv_t_r_rv$info_used_adjusted))
  expect_equal(unique(es.mcv_t_z$info_used_adjusted), unique(es.mcv_t_z_rv$info_used_adjusted))

  expect_equal(es.mcv_t_d$es_adjusted, -es.mcv_t_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_adjusted, es.mcv_t_d_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_adjusted, -es.mcv_t_g_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_adjusted, es.mcv_t_g_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_adjusted, -es.mcv_t_r_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_adjusted, es.mcv_t_r_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_adjusted, -es.mcv_t_z_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_adjusted, es.mcv_t_z_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_adjusted, -es.mcv_t_or_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_adjusted, es.mcv_t_or_rv$se_adjusted, tolerance = 1e-10)
})


test_that("ancova f pval - Reverse", {
  df.haza$ancova_t <- runif(nrow(df.haza), -5, 5)
  df.haza$ancova_f <- df.haza$ancova_t^2
  df.haza$ancova_f_pval <- 2 * pt(abs(df.haza$ancova_t), df.haza$n_exp + df.haza$n_nexp - 2 - df.haza$n_cov_ancova, lower.tail = FALSE)
  df.haza$ancova_t_pval <- df.haza$ancova_f_pval
  df.haza$n_cov_ancova <- round(runif(nrow(df.haza), 1, 5))
  df.haza$cov_outcome_r <- 0.3
  df.haza$reverse_ancova_f_pval <- FALSE
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp) & !is.na(ancova_f_pval))

  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "z"), digits = 11)

  dat$reverse_ancova_f_pval <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_f_pval", measure = "z"), digits = 11)

  expect_equal(es.mcv_t_d$es_adjusted, -es.mcv_t_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_adjusted, es.mcv_t_d_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_adjusted, -es.mcv_t_g_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_adjusted, es.mcv_t_g_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_adjusted, -es.mcv_t_r_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_adjusted, es.mcv_t_r_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_adjusted, -es.mcv_t_z_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_adjusted, es.mcv_t_z_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_adjusted, -es.mcv_t_or_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_adjusted, es.mcv_t_or_rv$se_adjusted, tolerance = 1e-10)
})

test_that("ancova t pval - Reverse", {
  df.haza$ancova_t <- runif(nrow(df.haza), -5, 5)
  df.haza$ancova_f <- df.haza$ancova_t^2
  df.haza$ancova_f_pval <- 2 * pt(abs(df.haza$ancova_t), df.haza$n_exp + df.haza$n_nexp - 2 - df.haza$n_cov_ancova, lower.tail = FALSE)
  df.haza$ancova_t_pval <- df.haza$ancova_f_pval
  df.haza$n_cov_ancova <- round(runif(nrow(df.haza), 1, 5))
  df.haza$cov_outcome_r <- 0.3
  df.haza$reverse_ancova_t_pval <- FALSE
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp) & !is.na(ancova_t_pval))

  es.mcv_t_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "z"), digits = 11)

  dat$reverse_ancova_t_pval <- TRUE
  es.mcv_t_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_t_pval", measure = "z"), digits = 11)

  expect_equal(es.mcv_t_d$es_adjusted, -es.mcv_t_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_d$se_adjusted, es.mcv_t_d_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_adjusted, -es.mcv_t_g_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_adjusted, es.mcv_t_g_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_adjusted, -es.mcv_t_r_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_adjusted, es.mcv_t_r_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_adjusted, -es.mcv_t_z_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_adjusted, es.mcv_t_z_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_adjusted, -es.mcv_t_or_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_adjusted, es.mcv_t_or_rv$se_adjusted, tolerance = 1e-10)
})



test_that("ancova means - Reverse", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_sd_exp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$ancova_mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp - 1)
  dat$ancova_mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp - 1)
  dat$ancova_mean_ci_lo_exp <- dat$mean_exp - dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_up_exp <- dat$mean_exp + dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_lo_nexp <- dat$mean_nexp - dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$ancova_mean_ci_up_nexp <- dat$mean_nexp + dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$ancova_mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * ancova_mean_sd_exp^2 + (n_nexp - 1) * ancova_mean_sd_nexp^2) / (n_exp + n_nexp - 2)))
  dat$mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * mean_sd_exp^2 + (n_nexp - 1) * mean_sd_nexp^2) / (n_exp + n_nexp - 2)))
  dat <- subset(dat, !is.na(ancova_mean_sd_pooled) & ancova_pval != 1)
  dat$reverse_ancova_means <- FALSE
  dat$reverse_plot_ancova_means <- FALSE

  es.mcv_m_sd_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "md"), digits = 11)
  es.mcv_m_se_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "md"), digits = 11)
  es.mcv_m_ci_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "md"), digits = 11)
  es.mcv_m_pooled_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "md"), digits = 11)
  es.mcv_m_pooled_md_adj <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "md"), digits = 11)

  es.mcv_m_sd_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "d"), digits = 11)
  es.mcv_m_se_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "d"), digits = 11)
  es.mcv_m_ci_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "d"), digits = 11)
  es.mcv_m_pooled_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "d"), digits = 11)
  es.mcv_m_pooled_d_adj <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "d"), digits = 11)

  es.mcv_m_sd_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "g"), digits = 11)
  es.mcv_m_se_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "g"), digits = 11)
  es.mcv_m_ci_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "g"), digits = 11)
  es.mcv_m_pooled_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "g"), digits = 11)
  es.mcv_m_pooled_g_adj <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "g"), digits = 11)

  es.mcv_m_sd_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "r"), digits = 11)
  es.mcv_m_se_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "r"), digits = 11)
  es.mcv_m_ci_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "r"), digits = 11)
  es.mcv_m_pooled_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "r"), digits = 11)
  es.mcv_m_pooled_r_adj <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "r"), digits = 11)

  es.mcv_m_sd_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "z"), digits = 11)
  es.mcv_m_se_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "z"), digits = 11)
  es.mcv_m_ci_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "z"), digits = 11)
  es.mcv_m_pooled_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "z"), digits = 11)
  es.mcv_m_pooled_z_adj <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "z"), digits = 11)

  es.mcv_m_sd_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "logor"), digits = 11)
  es.mcv_m_se_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "logor"), digits = 11)
  es.mcv_m_ci_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "logor"), digits = 11)
  es.mcv_m_pooled_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "logor"), digits = 11)
  es.mcv_m_pooled_or_adj <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "logor"), digits = 11)

  dat$reverse_ancova_means <- TRUE

  es.mcv_m_sd_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "md"), digits = 11)
  es.mcv_m_se_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "md"), digits = 11)
  es.mcv_m_ci_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "md"), digits = 11)
  es.mcv_m_pooled_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "md"), digits = 11)
  es.mcv_m_pooled_md_adj_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "md"), digits = 11)

  es.mcv_m_sd_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "d"), digits = 11)
  es.mcv_m_se_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "d"), digits = 11)
  es.mcv_m_ci_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "d"), digits = 11)
  es.mcv_m_pooled_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "d"), digits = 11)
  es.mcv_m_pooled_d_adj_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "d"), digits = 11)

  es.mcv_m_sd_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "g"), digits = 11)
  es.mcv_m_se_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "g"), digits = 11)
  es.mcv_m_ci_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "g"), digits = 11)
  es.mcv_m_pooled_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "g"), digits = 11)
  es.mcv_m_pooled_g_adj_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "g"), digits = 11)

  es.mcv_m_sd_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "r"), digits = 11)
  es.mcv_m_se_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "r"), digits = 11)
  es.mcv_m_ci_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "r"), digits = 11)
  es.mcv_m_pooled_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "r"), digits = 11)
  es.mcv_m_pooled_r_adj_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "r"), digits = 11)

  es.mcv_m_sd_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "z"), digits = 11)
  es.mcv_m_se_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "z"), digits = 11)
  es.mcv_m_ci_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "z"), digits = 11)
  es.mcv_m_pooled_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "z"), digits = 11)
  es.mcv_m_pooled_z_adj_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "z"), digits = 11)

  es.mcv_m_sd_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd", measure = "logor"), digits = 11)
  es.mcv_m_se_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_se", measure = "logor"), digits = 11)
  es.mcv_m_ci_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_ci", measure = "logor"), digits = 11)
  es.mcv_m_pooled_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled", measure = "logor"), digits = 11)
  es.mcv_m_pooled_or_adj_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_sd_pooled_adj", measure = "logor"), digits = 11)

  expect_equal(es.mcv_m_sd_d$info_used_adjusted, es.mcv_m_sd_d_rv$info_used_adjusted)
  expect_equal(es.mcv_m_sd_d$es_adjusted, -es.mcv_m_sd_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$es_adjusted, -es.mcv_m_se_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$es_adjusted, -es.mcv_m_ci_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_d$se_adjusted, es.mcv_m_sd_d_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$se_adjusted, es.mcv_m_se_d_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$se_adjusted, es.mcv_m_ci_d_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_md$info_used_adjusted, es.mcv_m_sd_md_rv$info_used_adjusted)
  expect_equal(es.mcv_m_sd_md$es_adjusted, -es.mcv_m_sd_md_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_md$es_adjusted, -es.mcv_m_se_md_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_md$es_adjusted, -es.mcv_m_ci_md_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_md$se_adjusted, es.mcv_m_sd_md_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_md$se_adjusted, es.mcv_m_se_md_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_md$se_adjusted, es.mcv_m_ci_md_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_d$info_used_adjusted, es.mcv_m_sd_d_rv$info_used_adjusted)
  expect_equal(es.mcv_m_sd_d$es_adjusted, -es.mcv_m_sd_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$es_adjusted, -es.mcv_m_se_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$es_adjusted, -es.mcv_m_ci_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_d$se_adjusted, es.mcv_m_sd_d_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_d$se_adjusted, es.mcv_m_se_d_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_d$se_adjusted, es.mcv_m_ci_d_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_g$info_used_adjusted, es.mcv_m_sd_g_rv$info_used_adjusted)
  expect_equal(es.mcv_m_sd_g$es_adjusted, -es.mcv_m_sd_g_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_g$es_adjusted, -es.mcv_m_se_g_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_g$es_adjusted, -es.mcv_m_ci_g_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_g$se_adjusted, es.mcv_m_sd_g_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_g$se_adjusted, es.mcv_m_se_g_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_g$se_adjusted, es.mcv_m_ci_g_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_r$info_used_adjusted, es.mcv_m_sd_r_rv$info_used_adjusted)
  expect_equal(es.mcv_m_sd_r$es_adjusted, -es.mcv_m_sd_r_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_r$es_adjusted, -es.mcv_m_se_r_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_r$es_adjusted, -es.mcv_m_ci_r_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_r$se_adjusted, es.mcv_m_sd_r_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_r$se_adjusted, es.mcv_m_se_r_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_r$se_adjusted, es.mcv_m_ci_r_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_z$info_used_adjusted, es.mcv_m_sd_z_rv$info_used_adjusted)
  expect_equal(es.mcv_m_sd_z$es_adjusted, -es.mcv_m_sd_z_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_z$es_adjusted, -es.mcv_m_se_z_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_z$es_adjusted, -es.mcv_m_ci_z_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_z$se_adjusted, es.mcv_m_sd_z_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_z$se_adjusted, es.mcv_m_se_z_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_z$se_adjusted, es.mcv_m_ci_z_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_sd_or$info_used_adjusted, es.mcv_m_sd_or_rv$info_used_adjusted)
  expect_equal(es.mcv_m_sd_or$es_adjusted, -es.mcv_m_sd_or_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_or$es_adjusted, -es.mcv_m_se_or_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_or$es_adjusted, -es.mcv_m_ci_or_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_sd_or$se_adjusted, es.mcv_m_sd_or_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_se_or$se_adjusted, es.mcv_m_se_or_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_ci_or$se_adjusted, es.mcv_m_ci_or_rv$se_adjusted, tolerance = 1e-10)
})


test_that("ancova means plot - Reverse", {
  df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")] <-
    suppressWarnings(apply(
      df.haza[, c("mean_exp", "mean_sd_exp", "n_exp", "mean_nexp", "mean_sd_nexp", "n_nexp")],
      2, function(x) as.numeric(as.character(x))
    ))
  dat <- subset(df.haza, !is.na(mean_sd_exp))
  dat$ancova_t <- runif(nrow(dat), -5, 5)
  dat$ancova_f <- dat$ancova_t^2
  dat$ancova_pval <- 2 * pt(abs(dat$ancova_t), dat$n_exp + dat$n_nexp - 2 - dat$n_cov_ancova, lower.tail = FALSE)
  dat$n_cov_ancova <- round(runif(nrow(dat), 1, 5))
  dat$cov_outcome_r <- 0.3
  dat$ancova_mean_exp <- dat$mean_exp
  dat$ancova_mean_nexp <- dat$mean_nexp
  dat$ancova_mean_sd_exp <- dat$mean_sd_exp
  dat$ancova_mean_sd_nexp <- dat$mean_sd_nexp
  dat$ancova_mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp - 1)
  dat$ancova_mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp - 1)
  dat$ancova_mean_ci_lo_exp <- dat$mean_exp - dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_up_exp <- dat$mean_exp + dat$ancova_mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$ancova_mean_ci_lo_nexp <- dat$mean_nexp - dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$ancova_mean_ci_up_nexp <- dat$mean_nexp + dat$ancova_mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$ancova_mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * ancova_mean_sd_exp^2 + (n_nexp - 1) * ancova_mean_sd_nexp^2) / (n_exp + n_nexp - 2)))
  dat$mean_sd_pooled <- with(dat, sqrt(((n_exp - 1) * mean_sd_exp^2 + (n_nexp - 1) * mean_sd_nexp^2) / (n_exp + n_nexp - 2)))
  dat <- subset(dat, !is.na(ancova_mean_sd_pooled) & ancova_pval != 1)
  dat$reverse_ancova_means <- FALSE
  dat$reverse_plot_ancova_means <- FALSE

  dat$plot_ancova_mean_exp <- dat$ancova_mean_exp
  dat$plot_ancova_mean_nexp <- dat$ancova_mean_nexp
  dat$plot_ancova_mean_sd_lo_exp <- dat$ancova_mean_exp - dat$ancova_mean_sd_exp
  dat$plot_ancova_mean_sd_up_exp <- dat$ancova_mean_exp + dat$ancova_mean_sd_exp
  dat$plot_ancova_mean_sd_lo_nexp <- dat$ancova_mean_nexp - dat$ancova_mean_sd_nexp
  dat$plot_ancova_mean_sd_up_nexp <- dat$ancova_mean_nexp + dat$ancova_mean_sd_nexp

  es.mcv_m_plot_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "logor"), digits = 11)
  es.mcv_m_plot_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "z"), digits = 11)
  es.mcv_m_plot_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "r"), digits = 11)
  es.mcv_m_plot_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "g"), digits = 11)
  es.mcv_m_plot_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "d"), digits = 11)
  es.mcv_m_plot_md <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "md"), digits = 11)

  dat$reverse_plot_ancova_means <- TRUE
  es.mcv_m_plot_md_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "md"), digits = 11)
  es.mcv_m_plot_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "g"), digits = 11)
  es.mcv_m_plot_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "r"), digits = 11)
  es.mcv_m_plot_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "z"), digits = 11)
  es.mcv_m_plot_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "logor"), digits = 11)
  es.mcv_m_plot_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "ancova_means_plot", measure = "d"), digits = 11)

  expect_equal(es.mcv_m_plot_d$es_adjusted, -es.mcv_m_plot_d_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_d$se_adjusted, es.mcv_m_plot_d_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_md$es_adjusted, -es.mcv_m_plot_md_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_md$se_adjusted, es.mcv_m_plot_md_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_g$se_adjusted, es.mcv_m_plot_g_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_g$es_adjusted, -es.mcv_m_plot_g_rv$es_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_r$es_adjusted, -es.mcv_m_plot_r_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_r$se_adjusted, es.mcv_m_plot_r_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_z$es_adjusted, -es.mcv_m_plot_z_rv$es_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_z$se_adjusted, es.mcv_m_plot_z_rv$se_adjusted, tolerance = 1e-10)

  expect_equal(es.mcv_m_plot_or$se_adjusted, es.mcv_m_plot_or_rv$se_adjusted, tolerance = 1e-10)
  expect_equal(es.mcv_m_plot_or$es_adjusted, -es.mcv_m_plot_or_rv$es_adjusted, tolerance = 1e-10)
})
