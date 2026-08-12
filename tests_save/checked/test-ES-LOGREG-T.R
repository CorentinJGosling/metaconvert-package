# logreg_t standalone ----
test_that("logreg_t - standalone function", {
  # OR path
  res_or <- es_from_logreg_t(
    or = 2.12, logreg_t = 3.21,
    n_cases = 50, n_controls = 150
  )
  expect_equal(res_or$info_used, "logreg_t")
  expect_equal(res_or$logor, log(2.12), tolerance = 1e-10)
  expect_equal(res_or$logor_se, abs(log(2.12)) / 3.21, tolerance = 1e-10)

  # RR path
  res_rr <- es_from_logreg_t(
    rr = 1.5, logreg_t = 2.8,
    n_exp = 100, n_nexp = 100
  )
  expect_equal(res_rr$info_used, "logreg_t")
  expect_equal(res_rr$logrr, log(1.5), tolerance = 1e-10)
  expect_equal(res_rr$logrr_se, abs(log(1.5)) / 2.8, tolerance = 1e-10)

  # logor input should match OR input
  res_logor <- es_from_logreg_t(
    logor = log(2.12), logreg_t = 3.21,
    n_cases = 50, n_controls = 150
  )
  expect_equal(res_or$logor, res_logor$logor, tolerance = 1e-10)
  expect_equal(res_or$logor_se, res_logor$logor_se, tolerance = 1e-10)
  expect_equal(res_or$d, res_logor$d, tolerance = 1e-10)
  expect_equal(res_or$d_se, res_logor$d_se, tolerance = 1e-10)
})

# logreg_t OR path equivalence ----
test_that("logreg_t - OR path: equivalence with or_se", {
  dat <- metaumbrella::df.OR
  or <- with(dat, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))
  dat$or <- or$value
  dat$logor_se <- or$se

  # Derive Wald t-statistic from known OR and SE
  dat$logreg_t <- log(dat$or) / dat$logor_se

  # Remove rows where OR=1 (log(OR)=0 makes Wald t=0, SE unrecoverable)
  dat <- dat[dat$or != 1, ]

  # Clear 2x2 table columns so only or+logor_se and logreg_t paths are available
  dat$n_cases_exp <- dat$n_cases_nexp <- dat$n_controls_exp <- dat$n_controls_nexp <- NA

  # --- logor measure ---
  es_logreg <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", measure = "logor"
  ), digits = 11)
  es_or_se <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "or_se", measure = "logor"
  ), digits = 11)

  expect_equal(unique(es_logreg$info_used_crude), "logreg_t")
  expect_equal(unique(es_or_se$info_used_crude), "or_se")
  expect_equal(es_logreg$es_crude, es_or_se$es_crude, tolerance = 1e-10)
  expect_equal(es_logreg$se_crude, es_or_se$se_crude, tolerance = 1e-10)

  # --- d measure: compare with esc ---
  if (requireNamespace("esc", quietly = TRUE)) {
    es_logreg_d <- summary(convert_df(dat,
      verbose = FALSE, es_selected = "hierarchy",
      hierarchy = "logreg_t", measure = "d"
    ), digits = 11)

    # esc::convert_or2d() doesn't accept vectors (is.na(se) on length>1 fails),
    # so call it row-by-row
    esc_or2d <- do.call(rbind, mapply(function(or_i, se_i, n_i) {
      data.frame(esc::convert_or2d(or = or_i, se = se_i, totaln = n_i))
    }, dat$or, dat$logor_se, dat$n_cases + dat$n_controls, SIMPLIFY = FALSE))
    expect_equal(es_logreg_d$es_crude, esc_or2d$es, tolerance = 1e-10)
    expect_equal(es_logreg_d$se_crude, esc_or2d$se, tolerance = 1e-10)
  }
})

# logreg_t RR path equivalence ----
test_that("logreg_t - RR path: equivalence with rr_se", {
  src <- metaumbrella::df.OR

  # Compute RR and its SE from the 2x2 table
  rr_val <- with(src,
    (n_cases_exp / (n_cases_exp + n_controls_exp)) /
    (n_cases_nexp / (n_cases_nexp + n_controls_nexp))
  )
  logrr_se_val <- with(src,
    sqrt(1 / n_cases_exp - 1 / (n_cases_exp + n_controls_exp) +
         1 / n_cases_nexp - 1 / (n_cases_nexp + n_controls_nexp))
  )

  # Remove rows where RR=1 (log(RR)=0 makes Wald t=0, SE unrecoverable)
  keep <- rr_val != 1
  rr_val <- rr_val[keep]
  logrr_se_val <- logrr_se_val[keep]
  src <- src[keep, ]

  # Build a clean dataframe with only RR-related columns
  dat <- data.frame(
    rr = rr_val,
    logrr_se = logrr_se_val,
    logreg_t = log(rr_val) / logrr_se_val,
    n_exp = src$n_exp,
    n_nexp = src$n_nexp,
    n_cases = src$n_cases,
    n_controls = src$n_controls,
    baseline_risk = with(src, n_cases_nexp / (n_cases_nexp + n_controls_nexp))
  )

  es_logreg <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", measure = "logrr"
  ), digits = 11)
  es_rr_se <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "rr_se", measure = "logrr"
  ), digits = 11)

  expect_equal(unique(es_logreg$info_used_crude), "logreg_t")
  expect_equal(unique(es_rr_se$info_used_crude), "rr_se")
  expect_equal(es_logreg$es_crude, es_rr_se$es_crude, tolerance = 1e-10)
  expect_equal(es_logreg$se_crude, es_rr_se$se_crude, tolerance = 1e-10)
})

# logreg_t reverse ----
test_that("logreg_t - reverse flag", {
  src <- metaumbrella::df.OR[1:5, ]
  or <- with(src, metaumbrella:::.estimate_or_from_n(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp
  ))

  # Build clean dataframe with only or + logreg_t
  dat <- data.frame(
    or = or$value,
    logreg_t = log(or$value) / or$se,
    n_cases = src$n_cases,
    n_controls = src$n_controls,
    n_exp = src$n_exp,
    n_nexp = src$n_nexp,
    baseline_risk = with(src, n_cases_nexp / (n_cases_nexp + n_controls_nexp))
  )

  # Without reverse
  es_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", measure = "d"), digits = 11)
  es_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", measure = "g"), digits = 11)
  es_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", measure = "logor"), digits = 11)
  es_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", or_to_cor = "bonett", measure = "r"), digits = 11)
  es_z <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", or_to_cor = "bonett", measure = "z"), digits = 11)

  # With reverse
  dat$reverse_logreg_t <- TRUE

  es_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", measure = "d"), digits = 11)
  es_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", measure = "g"), digits = 11)
  es_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", measure = "logor"), digits = 11)
  es_r_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", or_to_cor = "bonett", measure = "r"), digits = 11)
  es_z_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "logreg_t", or_to_cor = "bonett", measure = "z"), digits = 11)

  # info_used check
  expect_true(all(c(
    es_d$info_used_crude, es_g$info_used_crude, es_or$info_used_crude,
    es_r$info_used_crude, es_z$info_used_crude,
    es_d_r$info_used_crude, es_g_r$info_used_crude, es_or_r$info_used_crude,
    es_r_r$info_used_crude, es_z_r$info_used_crude
  ) == "logreg_t"))

  # ES negates, SE unchanged
  expect_equal(es_d$es_crude, -es_d_r$es_crude, tolerance = 1e-10)
  expect_equal(es_d$se_crude, es_d_r$se_crude, tolerance = 1e-10)
  expect_equal(es_g$es_crude, -es_g_r$es_crude, tolerance = 1e-10)
  expect_equal(es_g$se_crude, es_g_r$se_crude, tolerance = 1e-10)
  expect_equal(es_or$es_crude, -es_or_r$es_crude, tolerance = 1e-10)
  expect_equal(es_or$se_crude, es_or_r$se_crude, tolerance = 1e-10)
  expect_equal(es_r$es_crude, -es_r_r$es_crude, tolerance = 1e-10)
  expect_equal(es_r$se_crude, es_r_r$se_crude, tolerance = 1e-10)
  expect_equal(es_z$es_crude, -es_z_r$es_crude, tolerance = 1e-10)
  expect_equal(es_z$se_crude, es_z_r$se_crude, tolerance = 1e-10)
})
