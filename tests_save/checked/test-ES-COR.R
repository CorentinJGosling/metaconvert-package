test_that("R var", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  df.Z.mfr <- metafor::escalc(measure = "COR",
                              ri = pearson_r, ni = n_sample,
                              data = dat, digits = 12)

  es.mcv_d1 <- summary(convert_df(dat,
                                  verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
                                  cor_to_smd = "cooper", measure = "r"
  ), digits = 11)

  expect_equal(unique(es.mcv_d1$info_used_crude), "pearson_r")
  expect_equal(es.mcv_d1$es_crude, as.numeric(df.Z.mfr$yi))
  expect_equal(es.mcv_d1$se_crude, sqrt(as.numeric(df.Z.mfr$vi)))
})

test_that("Z var", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  df.Z.mfr <- metafor::escalc(measure = "ZCOR",
                              ri = pearson_r, ni = n_sample,
                              data = dat, digits = 12)

  es.mcv_d1 <- summary(convert_df(dat,
                                  verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
                                  cor_to_smd = "cooper", measure = "z"
  ), digits = 11)

  expect_equal(unique(es.mcv_d1$info_used_crude), "pearson_r")
  expect_equal(es.mcv_d1$es_crude, as.numeric(df.Z.mfr$yi))
  expect_equal(es.mcv_d1$se_crude, sqrt(as.numeric(df.Z.mfr$vi)))
})
test_that("R var", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  df.Z.mfr <- metafor::escalc(measure = "COR",
                              ri = pearson_r, ni = n_sample,
                              data = dat, digits = 12)

  es.mcv_d1 <- summary(convert_df(dat,
                                  verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
                                  cor_to_smd = "cooper", measure = "r"
  ), digits = 11)

  expect_equal(unique(es.mcv_d1$info_used_crude), "fisher_z")
  expect_equal(es.mcv_d1$es_crude, as.numeric(df.Z.mfr$yi))
  expect_equal(es.mcv_d1$se_crude, sqrt(as.numeric(df.Z.mfr$vi)))
})

test_that("Z var", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  df.Z.mfr <- metafor::escalc(measure = "ZCOR",
                              ri = pearson_r, ni = n_sample,
                              data = dat, digits = 12)

  es.mcv_d1 <- summary(convert_df(dat,
                                  verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
                                  cor_to_smd = "cooper", measure = "z"
  ), digits = 11)

  expect_equal(unique(es.mcv_d1$info_used_crude), "fisher_z")
  expect_equal(es.mcv_d1$es_crude, as.numeric(df.Z.mfr$yi))
  expect_equal(es.mcv_d1$se_crude, sqrt(as.numeric(df.Z.mfr$vi)))
})


test_that("R to Z", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  es.mcv_d1 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
    cor_to_smd = "cooper", measure = "z"
  ), digits = 11)

  expect_equal(es.mcv_d1$es_crude, dat$fisher_z)
})
# Z to R--------
test_that("Z to R", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  es.mcv_d1 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
    cor_to_smd = "cooper", measure = "r"
  ), digits = 11)

  expect_equal(es.mcv_d1$es_crude, dat$pearson_r)
})

# R/Z to SMD (COOPER) --------
test_that("R to SMD", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)

  es.mcv_d <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
    cor_to_smd = "cooper", measure = "d"
  ), digits = 11)
  es.mcv_or <- summary(convert_df(dat,
    verbose = FALSE,
    cor_to_smd = "cooper", es_selected = "hierarchy", hierarchy = "pearson_r", measure = "logor"
  ), digits = 11)


  expect_equal(unique(es.mcv_d$info_used_crude), "pearson_r")
  expect_equal(es.mcv_d$es_crude, effectsize::r_to_d(r = dat$pearson_r))
  expect_equal(es.mcv_or$es_crude, effectsize::r_to_oddsratio(r = dat$pearson_r, log = TRUE))
})
test_that("Z to SMD", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  es.mcv_d1 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
    cor_to_smd = "cooper", measure = "d"
  ), digits = 11)
  es.mcv_d2 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
    cor_to_smd = "cooper", measure = "d"
  ), digits = 11)

  expect_equal(unique(es.mcv_d1$info_used_crude), "pearson_r")
  expect_equal(unique(es.mcv_d2$info_used_crude), "fisher_z")
  expect_equal(es.mcv_d1$es_crude, es.mcv_d2$es_crude)
  expect_equal(es.mcv_d1$se_crude, es.mcv_d2$se_crude)
})

# R/Z to SMD (VIECHT) --------
test_that("R to SMD", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  es.mcv_d <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
    cor_to_smd = "viechtbauer", measure = "g"
  ), digits = 11)

  dat <- metafor::escalc(
    measure = "COR", ri = pearson_r, ni = n_sample, vtype = "LS", data = dat,
    digits = 12, var.names = c("ri", "vri")
  )
  # These rows SUPPLY n_exp / n_nexp, so the biserial r is inverted at the study's own
  # split (transf.rtod's n1i / n2i), the inverse of .smd_to_cor()'s forward map; the
  # balanced-arms constant (no n1i / n2i) is the reference only for n_sample-only rows.
  dat <- metafor::conv.delta(yi = ri, vi = vri, data = dat, transf = metafor::transf.rtod,
                             n1i = n_exp, n2i = n_nexp, var.names = c("yi", "vi"))

  expect_equal(unique(es.mcv_d$info_used_crude), "pearson_r")
  expect_equal(es.mcv_d$es_crude, as.numeric(dat$yi), tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude^2, dat$vi, tolerance = 1e-10)
})


test_that("Z to SMD", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)

  es.mcv_d <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
    cor_to_smd = "viechtbauer", measure = "g"
  ), digits = 11)

  dat <- metafor::escalc(
    measure = "COR", ri = pearson_r, ni = n_sample, vtype = "LS", data = dat,
    digits = 12, var.names = c("ri", "vri")
  )
  # These rows SUPPLY n_exp / n_nexp, so the biserial r is inverted at the study's own
  # split (transf.rtod's n1i / n2i), the inverse of .smd_to_cor()'s forward map; the
  # balanced-arms constant (no n1i / n2i) is the reference only for n_sample-only rows.
  dat <- metafor::conv.delta(yi = ri, vi = vri, data = dat, transf = metafor::transf.rtod,
                             n1i = n_exp, n2i = n_nexp, var.names = c("yi", "vi"))

  expect_equal(unique(es.mcv_d$info_used_crude), "fisher_z")
  expect_equal(es.mcv_d$es_crude, as.numeric(dat$yi), tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude^2, dat$vi, tolerance = 1e-10)
})

# R/Z to SMD (Mathur) --------
test_that("R to SMD", {
  dat <- data.frame(1)
  dat$n_sample <- 50;
  dat$pearson_r <- -0.50;
  dat$sd_iv = 0.24;
  dat$unit_increase_iv = 1;
  dat$unit_type = "raw_data"

  es.mcv_d <- summary(convert_df(dat,
                                 verbose = FALSE,
                                 es_selected = "hierarchy", hierarchy = "pearson_r",
                                 unit_type = "raw_data",
                                 cor_to_smd = "mathur",
                                 measure = "d"
  ), digits = 11)

  res = MetaUtility::r_to_d(r=dat$pearson_r,
                            sx = dat$sd_iv,
                            delta = dat$unit_increase_iv,
                            N = dat$n_sample)
  dat$pearson_r * dat$unit_increase_iv / (dat$sd_iv * sqrt((1 - dat$pearson_r^2)))
  expect_equal(unique(es.mcv_d$info_used_crude), "pearson_r")
  expect_equal(es.mcv_d$es_crude, res$d, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, res$se, tolerance = 1e-10)
})


### REVERSE R/Z------
test_that("R - Reverse", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)
  dat$reverse_pearson_r <- FALSE
  dat$sd_iv = abs(rnorm(nrow(dat)));
  dat$unit_increase_iv = 1;
  dat$unit_type = "raw_data"

  es.mcv_t_d1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
                                    cor_to_smd = "cooper",measure = "d"), digits = 11)
  es.mcv_t_d2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
                                    cor_to_smd = "viechtbauer",measure = "d"), digits = 11)
  es.mcv_t_d3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
                                    cor_to_smd = "mathur",measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r", measure = "z"), digits = 11)

  dat$reverse_pearson_r <- TRUE
  es.mcv_t_d1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
                                       cor_to_smd = "cooper",measure = "d"), digits = 11)
  es.mcv_t_d2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
                                       cor_to_smd = "viechtbauer",measure = "d"), digits = 11)
  es.mcv_t_d3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r",
                                       cor_to_smd = "mathur", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "pearson_r", measure = "z"), digits = 11)

  expect_equal(es.mcv_t_d1$es_crude, -es.mcv_t_d1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d1$se_crude, es.mcv_t_d1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d2$es_crude, -es.mcv_t_d2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d2$se_crude, es.mcv_t_d2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d3$es_crude, -es.mcv_t_d3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d3$se_crude, es.mcv_t_d3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_crude, -es.mcv_t_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_crude, es.mcv_t_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_crude, -es.mcv_t_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_crude, es.mcv_t_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_crude, -es.mcv_t_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_crude, es.mcv_t_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_crude, -es.mcv_t_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_crude, es.mcv_t_or_rv$se_crude, tolerance = 1e-10)
})
test_that("Z - Reverse", {
  dat <- subset(df.haza, !is.na(n_exp) & !is.na(n_nexp))
  dat$n_sample <- dat$n_exp + dat$n_nexp

  dat$pearson_r <- runif(nrow(dat), -0.5, 0.5)
  dat$fisher_z <- atanh(dat$pearson_r)
  dat$reverse_fisher_z <- FALSE
  dat$sd_iv = abs(rnorm(nrow(dat)));
  dat$unit_increase_iv = 1;
  dat$unit_type = "raw_data"

  es.mcv_t_d1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
                                    cor_to_smd = "cooper",measure = "d"), digits = 11)
  es.mcv_t_d2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
                                    cor_to_smd = "viechtbauer",measure = "d"), digits = 11)
  es.mcv_t_d3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
                                    cor_to_smd = "mathur",measure = "d"), digits = 11)
  es.mcv_t_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z", measure = "g"), digits = 11)
  es.mcv_t_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z", measure = "logor"), digits = 11)
  es.mcv_t_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z", measure = "r"), digits = 11)
  es.mcv_t_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z", measure = "z"), digits = 11)

  dat$reverse_fisher_z <- TRUE
  es.mcv_t_d1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
                                       cor_to_smd = "cooper",measure = "d"), digits = 11)
  es.mcv_t_d2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
                                       cor_to_smd = "viechtbauer",measure = "d"), digits = 11)
  es.mcv_t_d3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z",
                                       cor_to_smd = "mathur", measure = "d"), digits = 11)
  es.mcv_t_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z", measure = "g"), digits = 11)
  es.mcv_t_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z", measure = "logor"), digits = 11)
  es.mcv_t_r_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z", measure = "r"), digits = 11)
  es.mcv_t_z_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "fisher_z", measure = "z"), digits = 11)

  row = which(!is.na(es.mcv_t_d3$es_crude))
  expect_equal(es.mcv_t_d1$es_crude, -es.mcv_t_d1_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d1$se_crude, es.mcv_t_d1_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d2$es_crude, -es.mcv_t_d2_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d2$se_crude, es.mcv_t_d2_rv$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_d3$es_crude[row], -es.mcv_t_d3_rv$es_crude[row], tolerance = 1e-10)
  expect_equal(es.mcv_t_d3$se_crude, es.mcv_t_d3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_g$es_crude, -es.mcv_t_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_g$se_crude, es.mcv_t_g_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_r$es_crude, -es.mcv_t_r_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_r$se_crude, es.mcv_t_r_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_z$es_crude, -es.mcv_t_z_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_z$se_crude, es.mcv_t_z_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_t_or$es_crude, -es.mcv_t_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_t_or$se_crude, es.mcv_t_or_rv$se_crude, tolerance = 1e-10)
})
