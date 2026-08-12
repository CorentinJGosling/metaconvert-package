### VR ====================
test_that("VR from means SD", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    mean_sd_exp != 0 & mean_sd_nexp != 0 &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  dat_VR <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "VR", data = dat, digits = 11
  ))
  dat_CVR <- suppressWarnings(metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "CVR", data = dat, digits = 11
  ))

  es.mcv_vr_sd <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_sd", measure = "logvr"
  ), digits = 11)

  expect_equal(unique(es.mcv_vr_sd$info_used_crude), "variability_means_sd")
  expect_equal(es.mcv_vr_sd$es_crude, as.numeric(dat_VR$yi), tolerance = 1e-10)
  expect_equal(es.mcv_vr_sd$se_crude^2, dat_VR$vi, tolerance = 1e-1)
})
test_that("VR from means SE", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    mean_sd_exp != 0 & mean_sd_nexp != 0 &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)

  dat_VR <- metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "VR", data = dat, digits = 11
  )
  dat_CVR <- metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "CVR", data = dat, digits = 11
  )

  es.mcv_vr_sd <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_se", measure = "logvr"
  ), digits = 11)

  expect_equal(unique(es.mcv_vr_sd$info_used_crude), "variability_means_se")
  expect_equal(es.mcv_vr_sd$es_crude, as.numeric(dat_VR$yi), tolerance = 1e-10)
  expect_equal(es.mcv_vr_sd$se_crude^2, dat_VR$vi, tolerance = 1e-1)
})
test_that("VR from means CI", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    mean_sd_exp != 0 & mean_sd_nexp != 0 &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$mean_ci_lo_exp <- dat$mean_exp - dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$mean_ci_up_exp <- dat$mean_exp + dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_up_nexp <- dat$mean_nexp + dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)

  dat_VR <- metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "VR", data = dat, digits = 11
  )
  dat_CVR <- metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "CVR", data = dat, digits = 11
  )
  es.mcv_vr_sd <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_ci", measure = "logvr"
  ), digits = 11)

  expect_equal(unique(es.mcv_vr_sd$info_used_crude), "variability_means_ci")
  expect_equal(es.mcv_vr_sd$es_crude, as.numeric(dat_VR$yi), tolerance = 1e-10)
  expect_equal(es.mcv_vr_sd$se_crude^2, dat_VR$vi, tolerance = 1e-1)
})
### CVR ====================
test_that("CVR from means SD", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    mean_sd_exp != 0 & mean_sd_nexp != 0 &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$mean_ci_lo_exp <- dat$mean_exp - dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$mean_ci_up_exp <- dat$mean_exp + dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_up_nexp <- dat$mean_nexp + dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)

  dat_CVR <- metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "CVR", data = dat, digits = 11
  )

  es.mcv_vr_sd <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_sd", measure = "logcvr"
  ), digits = 11)

  expect_equal(unique(es.mcv_vr_sd$info_used_crude), "variability_means_sd")
  expect_equal(es.mcv_vr_sd$es_crude, as.numeric(dat_CVR$yi), tolerance = 1e-1)
  expect_equal(es.mcv_vr_sd$se_crude^2, dat_CVR$vi, tolerance = 1e-1)
})
test_that("CVR from means SE", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    mean_sd_exp != 0 & mean_sd_nexp != 0 &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$mean_ci_lo_exp <- dat$mean_exp - dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$mean_ci_up_exp <- dat$mean_exp + dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_up_nexp <- dat$mean_nexp + dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)

  dat_CVR <- metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "CVR", data = dat, digits = 11
  )

  es.mcv_vr_sd <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_se", measure = "logcvr"
  ), digits = 11)

  expect_equal(unique(es.mcv_vr_sd$info_used_crude), "variability_means_se")
  expect_equal(es.mcv_vr_sd$es_crude, as.numeric(dat_CVR$yi), tolerance = 1e-1)
  expect_equal(es.mcv_vr_sd$se_crude^2, dat_CVR$vi, tolerance = 1e-1)
})
test_that("CVR from means CI", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    mean_sd_exp != 0 & mean_sd_nexp != 0 &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  # dat[which(dat$mean_sd_exp == 0), ]$mean_sd_exp <- 0.1
  # dat[which(dat$mean_sd_nexp == 0), ]$mean_sd_nexp <- 0.1
  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$mean_ci_lo_exp <- dat$mean_exp - dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$mean_ci_up_exp <- dat$mean_exp + dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_up_nexp <- dat$mean_nexp + dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)

  dat_CVR <- metafor::escalc(
    m1i = as.numeric(as.character(mean_exp)),
    sd1i = as.numeric(as.character(mean_sd_exp)),
    n1i = as.numeric(as.character(n_exp)),
    m2i = as.numeric(as.character(mean_nexp)),
    sd2i = as.numeric(as.character(mean_sd_nexp)),
    n2i = as.numeric(as.character(n_nexp)),
    measure = "CVR", data = dat, digits = 11
  )
  es.mcv_vr_sd <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_ci", measure = "logcvr"
  ), digits = 11)

  expect_equal(unique(es.mcv_vr_sd$info_used_crude), "variability_means_ci")
  expect_equal(es.mcv_vr_sd$es_crude, as.numeric(dat_CVR$yi), tolerance = 1e-1)
  expect_equal(es.mcv_vr_sd$se_crude^2, dat_CVR$vi, tolerance = 1e-1)
})

### REVERSE =========================

test_that("var - reverse", {
  dat <- subset(df.haza, !is.na(mean_exp) & !is.na(mean_sd_exp) &
    mean_sd_exp != 0 & mean_sd_nexp != 0 &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  # dat[which(dat$mean_sd_exp == 0), ]$mean_sd_exp <- 0.1
  # dat[which(dat$mean_sd_nexp == 0), ]$mean_sd_nexp <- 0.1
  dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
  dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)
  dat$mean_ci_lo_exp <- dat$mean_exp - dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_lo_nexp <- dat$mean_nexp - dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)
  dat$mean_ci_up_exp <- dat$mean_exp + dat$mean_se_exp * qt(0.975, dat$n_exp - 1)
  dat$mean_ci_up_nexp <- dat$mean_nexp + dat$mean_se_nexp * qt(0.975, dat$n_nexp - 1)

  dat$reverse_means_variability <- FALSE
  es.mcv_vr_ci <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_ci", measure = "logvr"
  ), digits = 11)
  es.mcv_vr_se <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_se", measure = "logvr"
  ), digits = 11)
  es.mcv_vr_sd <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_sd", measure = "logvr"
  ), digits = 11)
  es.mcv_cvr_ci <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_ci", measure = "logcvr"
  ), digits = 11)
  es.mcv_cvr_se <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_se", measure = "logcvr"
  ), digits = 11)
  es.mcv_cvr_sd <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_sd", measure = "logcvr"
  ), digits = 11)

  dat$reverse_means_variability <- TRUE
  es.mcv_vr_ci_rv <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_ci", measure = "logvr"
  ), digits = 11)
  es.mcv_vr_se_rv <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_se", measure = "logvr"
  ), digits = 11)
  es.mcv_vr_sd_rv <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_sd", measure = "logvr"
  ), digits = 11)
  es.mcv_cvr_ci_rv <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_ci", measure = "logcvr"
  ), digits = 11)
  es.mcv_cvr_se_rv <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_se", measure = "logcvr"
  ), digits = 11)
  es.mcv_cvr_sd_rv <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "variability_means_sd", measure = "logcvr"
  ), digits = 11)

  expect_equal(es.mcv_vr_sd$info_used_crude, es.mcv_vr_sd_rv$info_used_crude)
  expect_equal(es.mcv_vr_sd$es_crude, -es.mcv_vr_sd_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_vr_sd$se_crude, es.mcv_vr_sd_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_vr_se$info_used_crude, es.mcv_vr_se_rv$info_used_crude)
  expect_equal(es.mcv_vr_se$es_crude, -es.mcv_vr_se_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_vr_se$se_crude, es.mcv_vr_se_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_vr_ci$info_used_crude, es.mcv_vr_ci_rv$info_used_crude)
  expect_equal(es.mcv_vr_ci$es_crude, -es.mcv_vr_ci_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_vr_ci$se_crude, es.mcv_vr_ci_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_cvr_sd$info_used_crude, es.mcv_cvr_sd_rv$info_used_crude)
  expect_equal(es.mcv_cvr_sd$es_crude, -es.mcv_cvr_sd_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_cvr_sd$se_crude, es.mcv_cvr_sd_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_cvr_se$info_used_crude, es.mcv_cvr_se_rv$info_used_crude)
  expect_equal(es.mcv_cvr_se$es_crude, -es.mcv_cvr_se_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_cvr_se$se_crude, es.mcv_cvr_se_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_cvr_ci$info_used_crude, es.mcv_cvr_ci_rv$info_used_crude)
  expect_equal(es.mcv_cvr_ci$es_crude, -es.mcv_cvr_ci_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_cvr_ci$se_crude, es.mcv_cvr_ci_rv$se_crude, tolerance = 1e-10)
})
