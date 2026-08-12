# cases/times to IRR ----
test_that("1. cases/times to IRR", {
  dat <- metaumbrella::df.IRR

  es.mcv_irr <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "cases_time",
    measure = "irr"
  ), digits = 11)

  df.IRR.mfr <- metafor::escalc(
    measure = "IRR",
    x1i = n_cases_exp, t1i = time_exp,
    x2i = n_cases_nexp, t2i = time_nexp,
    data = dat
  )

  expect_equal(es.mcv_irr$es_crude, exp(as.numeric(as.character(df.IRR.mfr$yi))), tolerance = 1e-10)
  expect_equal(es.mcv_irr$se_crude, sqrt(df.IRR.mfr$vi), tolerance = 1e-10)

  es.mcv_logirr <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "cases_time",
    measure = "logirr"
  ), digits = 11)
  expect_equal(es.mcv_logirr$es_crude, as.numeric(as.character(df.IRR.mfr$yi)), tolerance = 1e-10)
  expect_equal(es.mcv_logirr$se_crude, sqrt(df.IRR.mfr$vi), tolerance = 1e-10)
})

# REVERSE ----

test_that("IRR - Reverse", {
  dat <- metaumbrella::df.IRR
  dat$reverse_irr <- FALSE
  es.mcv_irr <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "cases_time",
    measure = "irr"
  ), digits = 11)
  dat$reverse_irr <- TRUE

  es.mcv_irr_rv <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "cases_time",
    measure = "irr"
  ), digits = 11)

  expect_equal(es.mcv_irr$info_used_crude, es.mcv_irr_rv$info_used_crude)

  expect_equal(es.mcv_irr$es_crude, 1 / es.mcv_irr_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_irr$se_crude, es.mcv_irr_rv$se_crude, tolerance = 1e-10)
})
test_that("logIRR - Reverse", {
  dat <- metaumbrella::df.IRR
  dat$reverse_irr <- FALSE
  es.mcv_irr <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "cases_time",
    measure = "logirr"
  ), digits = 11)
  dat$reverse_irr <- TRUE

  es.mcv_irr_rv <- summary(convert_df(dat,
    verbose = FALSE,
    es_selected = "hierarchy", hierarchy = "cases_time",
    measure = "logirr"
  ), digits = 11)

  expect_equal(es.mcv_irr$info_used_crude, es.mcv_irr_rv$info_used_crude)

  expect_equal(es.mcv_irr$es_crude, -es.mcv_irr_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_irr$se_crude, es.mcv_irr_rv$se_crude, tolerance = 1e-10)
})
