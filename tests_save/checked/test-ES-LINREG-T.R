# linreg_t standalone ----
test_that("linreg_t - standalone function", {
  # Aloe & Thompson (2013) Table 1, Cole et al. (2004): t=6.19, n=232, k=6
  res <- es_from_linreg_t(linreg_t = 6.19, n_sample = 232, n_covariates = 6)

  df <- 232 - 6 - 2  # = 224 (residual df: n - intercept - focal - covariates)
  rp_expected <- 6.19 / sqrt(6.19^2 + df)
  rp_se_expected <- sqrt((1 - rp_expected^2)^2 / df)
  zp_expected <- atanh(rp_expected)
  zp_se_expected <- sqrt(1 / (df - 1))  # = 1 / (n - n_covariates - 3)

  expect_equal(res$info_used, "linreg_t")
  expect_equal(res$rp, rp_expected, tolerance = 1e-10)
  expect_equal(res$rp_se, rp_se_expected, tolerance = 1e-10)
  expect_equal(res$rp_ci_lo, rp_expected - qt(.975, df) * rp_se_expected, tolerance = 1e-10)
  expect_equal(res$rp_ci_up, rp_expected + qt(.975, df) * rp_se_expected, tolerance = 1e-10)
  expect_equal(res$zp, zp_expected, tolerance = 1e-10)
  expect_equal(res$zp_se, zp_se_expected, tolerance = 1e-10)
  expect_equal(res$zp_ci_lo, zp_expected - qnorm(.975) * zp_se_expected, tolerance = 1e-10)
  expect_equal(res$zp_ci_up, zp_expected + qnorm(.975) * zp_se_expected, tolerance = 1e-10)

  # Simple case: t=2.5, n=100, k=2
  res2 <- es_from_linreg_t(linreg_t = 2.5, n_sample = 100, n_covariates = 2)

  df2 <- 100 - 2 - 2  # = 96
  rp2 <- 2.5 / sqrt(2.5^2 + df2)
  rp_se2 <- sqrt((1 - rp2^2)^2 / df2)
  zp2 <- atanh(rp2)
  zp_se2 <- sqrt(1 / (df2 - 1))

  expect_equal(res2$info_used, "linreg_t")
  expect_equal(res2$rp, rp2, tolerance = 1e-10)
  expect_equal(res2$rp_se, rp_se2, tolerance = 1e-10)
  expect_equal(res2$zp, zp2, tolerance = 1e-10)
  expect_equal(res2$zp_se, zp_se2, tolerance = 1e-10)

  # Negative t should give negative rp
  res_neg <- es_from_linreg_t(linreg_t = -3.0, n_sample = 50, n_covariates = 1)
  expect_true(res_neg$rp < 0)
  expect_equal(res_neg$rp, -3.0 / sqrt(9 + 47), tolerance = 1e-10)  # df = 50 - 1 - 2 = 47
})

# linreg_t ground-truth anchor (guards against df regression) ----
test_that("linreg_t - rp equals the true partial correlation (lm) and metafor PCOR", {
  skip_if_not_installed("metafor")

  set.seed(2024)
  N <- 40; k <- 3
  Z <- matrix(rnorm(N * k), N)
  X <- rnorm(N)
  Y <- 0.5 * X + Z %*% rnorm(k) + rnorm(N)
  fit <- lm(Y ~ X + Z)
  tval <- summary(fit)$coefficients["X", "t value"]

  # Ground truth 1: partial correlation via residualization
  rp_true <- cor(resid(lm(X ~ Z)), resid(lm(Y ~ Z)))
  # Ground truth 2: metafor's PCOR (mi = total predictors including the focal one)
  mf <- metafor::escalc(measure = "PCOR", ti = tval, ni = N, mi = k + 1)

  res <- es_from_linreg_t(linreg_t = tval, n_sample = N, n_covariates = k)

  expect_equal(res$rp, rp_true, tolerance = 1e-8)
  expect_equal(res$rp, as.numeric(mf$yi), tolerance = 1e-8)
  # Fisher-z variance of a partial correlation controlling for k covariates = 1 / (n - k - 3)
  expect_equal(res$zp_se, sqrt(1 / (N - k - 3)), tolerance = 1e-10)
})

# linreg_t formula verification ----
test_that("linreg_t - formula verification with multiple studies", {
  # Multiple studies from Aloe & Thompson (2013) Table 1
  t_vals <- c(6.19, 2.5, -1.8, 4.32, 0.75)
  n_vals <- c(232, 100, 50, 500, 30)
  k_vals <- c(6, 2, 1, 10, 3)

  res <- es_from_linreg_t(
    linreg_t = t_vals, n_sample = n_vals, n_covariates = k_vals
  )

  for (i in seq_along(t_vals)) {
    df_i <- n_vals[i] - k_vals[i] - 2
    rp_i <- t_vals[i] / sqrt(t_vals[i]^2 + df_i)
    rp_se_i <- sqrt((1 - rp_i^2)^2 / df_i)
    zp_i <- atanh(rp_i)
    zp_se_i <- sqrt(1 / (df_i - 1))

    expect_equal(res$rp[i], rp_i, tolerance = 1e-10)
    expect_equal(res$rp_se[i], rp_se_i, tolerance = 1e-10)
    expect_equal(res$zp[i], zp_i, tolerance = 1e-10)
    expect_equal(res$zp_se[i], zp_se_i, tolerance = 1e-10)
  }

  # All should have same info_used
  expect_equal(res$info_used, rep("linreg_t", 5))
})

# linreg_t convert_df integration ----
test_that("linreg_t - convert_df integration for rp/zp measures", {
  # Build a dataframe with linreg_t data
  dat <- data.frame(
    linreg_t = c(6.19, 2.5, -1.8, 4.32),
    n_sample = c(232, 100, 50, 500),
    n_covariates = c(6, 2, 1, 10)
  )

  # --- rp measure ---
  es_rp <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_t", measure = "rp"
  ), digits = 11)

  expect_equal(unique(es_rp$info_used_crude), "linreg_t")
  expect_equal(nrow(es_rp), 4)

  # Verify values match standalone function
  standalone <- es_from_linreg_t(
    linreg_t = dat$linreg_t,
    n_sample = dat$n_sample,
    n_covariates = dat$n_covariates
  )
  expect_equal(es_rp$es_crude, standalone$rp, tolerance = 1e-10)
  expect_equal(es_rp$se_crude, standalone$rp_se, tolerance = 1e-10)

  # --- zp measure ---
  es_zp <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_t", measure = "zp"
  ), digits = 11)

  expect_equal(unique(es_zp$info_used_crude), "linreg_t")
  expect_equal(es_zp$es_crude, standalone$zp, tolerance = 1e-10)
  expect_equal(es_zp$se_crude, standalone$zp_se, tolerance = 1e-10)
})

# linreg_t reverse ----
test_that("linreg_t - reverse flag", {
  dat <- data.frame(
    linreg_t = c(6.19, 2.5, -1.8, 4.32, 0.75),
    n_sample = c(232, 100, 50, 500, 30),
    n_covariates = c(6, 2, 1, 10, 3)
  )

  # Without reverse - rp measure
  es_rp <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_t", measure = "rp"), digits = 11)
  es_zp <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_t", measure = "zp"), digits = 11)

  # With reverse
  dat$reverse_linreg_t <- TRUE

  es_rp_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_t", measure = "rp"), digits = 11)
  es_zp_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_t", measure = "zp"), digits = 11)

  # info_used check
  expect_true(all(c(
    es_rp$info_used_crude, es_zp$info_used_crude,
    es_rp_r$info_used_crude, es_zp_r$info_used_crude
  ) == "linreg_t"))

  # ES negates, SE unchanged
  expect_equal(es_rp$es_crude, -es_rp_r$es_crude, tolerance = 1e-10)
  expect_equal(es_rp$se_crude, es_rp_r$se_crude, tolerance = 1e-10)
  expect_equal(es_zp$es_crude, -es_zp_r$es_crude, tolerance = 1e-10)
  expect_equal(es_zp$se_crude, es_zp_r$se_crude, tolerance = 1e-10)
})
