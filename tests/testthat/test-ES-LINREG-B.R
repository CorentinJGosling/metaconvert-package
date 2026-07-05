# linreg_b_se standalone ----
test_that("linreg_b_se - standalone function", {
  # b=1.5, se=0.6 => t = 2.5, n=100, k=2
  res <- es_from_linreg_b_se(linreg_b = 1.5, linreg_b_se = 0.6,
                              n_sample = 100, n_covariates = 2)

  t_val <- 1.5 / 0.6  # = 2.5
  df <- 100 - 2 - 1   # = 97
  rp_expected <- t_val / sqrt(t_val^2 + df)
  rp_se_expected <- sqrt((1 - rp_expected^2)^2 / df)
  zp_expected <- atanh(rp_expected)
  zp_se_expected <- sqrt(1 / (df - 2))

  expect_equal(res$info_used, "linreg_b_se")
  expect_equal(res$rp, rp_expected, tolerance = 1e-10)
  expect_equal(res$rp_se, rp_se_expected, tolerance = 1e-10)
  expect_equal(res$rp_ci_lo, rp_expected - qt(.975, df) * rp_se_expected, tolerance = 1e-10)
  expect_equal(res$rp_ci_up, rp_expected + qt(.975, df) * rp_se_expected, tolerance = 1e-10)
  expect_equal(res$zp, zp_expected, tolerance = 1e-10)
  expect_equal(res$zp_se, zp_se_expected, tolerance = 1e-10)

  # Negative b should give negative rp
  res_neg <- es_from_linreg_b_se(linreg_b = -2.0, linreg_b_se = 0.8,
                                  n_sample = 50, n_covariates = 1)
  expect_true(res_neg$rp < 0)
  t_neg <- -2.0 / 0.8
  df_neg <- 50 - 1 - 1
  expect_equal(res_neg$rp, t_neg / sqrt(t_neg^2 + df_neg), tolerance = 1e-10)
})

# linreg_b_se matches linreg_t ----
test_that("linreg_b_se - matches linreg_t with equivalent input", {
  b <- 1.5; se <- 0.6; n <- 100; k <- 2
  t_val <- b / se

  res_b <- es_from_linreg_b_se(linreg_b = b, linreg_b_se = se,
                                n_sample = n, n_covariates = k)
  res_t <- es_from_linreg_t(linreg_t = t_val, n_sample = n, n_covariates = k)

  expect_equal(res_b$rp, res_t$rp, tolerance = 1e-10)
  expect_equal(res_b$rp_se, res_t$rp_se, tolerance = 1e-10)
  expect_equal(res_b$zp, res_t$zp, tolerance = 1e-10)
  expect_equal(res_b$zp_se, res_t$zp_se, tolerance = 1e-10)
  expect_equal(res_b$d, res_t$d, tolerance = 1e-10)
  expect_equal(res_b$info_used, "linreg_b_se")
  expect_equal(res_t$info_used, "linreg_t")
})

# linreg_b_ci standalone ----
test_that("linreg_b_ci - standalone function", {
  b <- 1.5; se <- 0.6; n <- 100; k <- 2
  df <- n - k - 1
  ci_lo <- b - qt(.975, df) * se
  ci_up <- b + qt(.975, df) * se

  res_ci <- es_from_linreg_b_ci(linreg_b = b, linreg_b_ci_lo = ci_lo,
                                 linreg_b_ci_up = ci_up,
                                 n_sample = n, n_covariates = k)
  res_se <- es_from_linreg_b_se(linreg_b = b, linreg_b_se = se,
                                 n_sample = n, n_covariates = k)

  expect_equal(res_ci$info_used, "linreg_b_ci")
  expect_equal(res_ci$rp, res_se$rp, tolerance = 1e-8)
  expect_equal(res_ci$rp_se, res_se$rp_se, tolerance = 1e-8)
  expect_equal(res_ci$zp, res_se$zp, tolerance = 1e-8)
  expect_equal(res_ci$zp_se, res_se$zp_se, tolerance = 1e-8)
})

# linreg_b_pval standalone ----
test_that("linreg_b_pval - standalone function", {
  b <- 1.5; se <- 0.6; n <- 100; k <- 2
  df <- n - k - 1
  t_val <- b / se
  pval <- 2 * pt(-abs(t_val), df)

  res_pval <- es_from_linreg_b_pval(linreg_b = b, linreg_b_pval = pval,
                                     n_sample = n, n_covariates = k)
  res_se <- es_from_linreg_b_se(linreg_b = b, linreg_b_se = se,
                                 n_sample = n, n_covariates = k)

  expect_equal(res_pval$info_used, "linreg_b_pval")
  expect_equal(res_pval$rp, res_se$rp, tolerance = 1e-8)
  expect_equal(res_pval$rp_se, res_se$rp_se, tolerance = 1e-8)
  expect_equal(res_pval$zp, res_se$zp, tolerance = 1e-8)
  expect_equal(res_pval$zp_se, res_se$zp_se, tolerance = 1e-8)
})

# linreg_b_pval sign handling ----
test_that("linreg_b_pval - negative b produces negative rp", {
  b_pos <- 1.5; b_neg <- -1.5; pval <- 0.01; n <- 100; k <- 2

  res_pos <- es_from_linreg_b_pval(linreg_b = b_pos, linreg_b_pval = pval,
                                    n_sample = n, n_covariates = k)
  res_neg <- es_from_linreg_b_pval(linreg_b = b_neg, linreg_b_pval = pval,
                                    n_sample = n, n_covariates = k)

  expect_true(res_pos$rp > 0)
  expect_true(res_neg$rp < 0)
  expect_equal(abs(res_pos$rp), abs(res_neg$rp), tolerance = 1e-10)
  expect_equal(res_pos$rp_se, res_neg$rp_se, tolerance = 1e-10)
})

# linreg_b multiple studies ----
test_that("linreg_b - multiple studies", {
  b_vals <- c(1.5, -0.8, 2.3, 0.1)
  se_vals <- c(0.6, 0.3, 0.9, 0.5)
  n_vals <- c(100, 50, 200, 80)
  k_vals <- c(2, 1, 5, 3)

  res <- es_from_linreg_b_se(linreg_b = b_vals, linreg_b_se = se_vals,
                              n_sample = n_vals, n_covariates = k_vals)

  expect_equal(length(res$rp), 4)
  expect_equal(res$info_used, rep("linreg_b_se", 4))

  # Verify each manually
  for (i in 1:4) {
    t_i <- b_vals[i] / se_vals[i]
    df_i <- n_vals[i] - k_vals[i] - 1
    rp_i <- t_i / sqrt(t_i^2 + df_i)
    expect_equal(res$rp[i], rp_i, tolerance = 1e-10)
  }
})

# linreg_b_se convert_df integration ----
test_that("linreg_b_se - convert_df integration for rp/zp measures", {
  dat <- data.frame(
    linreg_b = c(1.5, -0.8, 2.3),
    linreg_b_se = c(0.6, 0.3, 0.9),
    n_sample = c(100, 50, 200),
    n_covariates = c(2, 1, 5)
  )

  # --- rp measure ---
  es_rp <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_b_se", measure = "rp"
  ), digits = 11)

  expect_equal(unique(es_rp$info_used_crude), "linreg_b_se")
  expect_equal(nrow(es_rp), 3)

  standalone <- es_from_linreg_b_se(
    linreg_b = dat$linreg_b, linreg_b_se = dat$linreg_b_se,
    n_sample = dat$n_sample, n_covariates = dat$n_covariates
  )
  expect_equal(es_rp$es_crude, standalone$rp, tolerance = 1e-10)
  expect_equal(es_rp$se_crude, standalone$rp_se, tolerance = 1e-10)

  # --- zp measure ---
  es_zp <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_b_se", measure = "zp"
  ), digits = 11)

  expect_equal(unique(es_zp$info_used_crude), "linreg_b_se")
  expect_equal(es_zp$es_crude, standalone$zp, tolerance = 1e-10)
  expect_equal(es_zp$se_crude, standalone$zp_se, tolerance = 1e-10)
})

# linreg_b_ci convert_df integration ----
test_that("linreg_b_ci - convert_df integration", {
  b <- c(1.5, -0.8); se <- c(0.6, 0.3)
  n <- c(100, 50); k <- c(2, 1)
  df <- n - k - 1
  ci_lo <- b - qt(.975, df) * se
  ci_up <- b + qt(.975, df) * se

  dat <- data.frame(
    linreg_b = b,
    linreg_b_ci_lo = ci_lo,
    linreg_b_ci_up = ci_up,
    n_sample = n,
    n_covariates = k
  )

  es_rp <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_b_ci", measure = "rp"
  ), digits = 11)

  expect_equal(unique(es_rp$info_used_crude), "linreg_b_ci")
  expect_equal(nrow(es_rp), 2)

  # Should match the SE-based results
  standalone_se <- es_from_linreg_b_se(
    linreg_b = b, linreg_b_se = se,
    n_sample = n, n_covariates = k
  )
  expect_equal(es_rp$es_crude, standalone_se$rp, tolerance = 1e-8)
})

# linreg_b_pval convert_df integration ----
test_that("linreg_b_pval - convert_df integration", {
  b <- c(1.5, -0.8); se <- c(0.6, 0.3)
  n <- c(100, 50); k <- c(2, 1)
  df <- n - k - 1
  t_vals <- b / se
  pvals <- 2 * pt(-abs(t_vals), df)

  dat <- data.frame(
    linreg_b = b,
    linreg_b_pval = pvals,
    n_sample = n,
    n_covariates = k
  )

  es_rp <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_b_pval", measure = "rp"
  ), digits = 11)

  expect_equal(unique(es_rp$info_used_crude), "linreg_b_pval")

  # Should match the SE-based results
  standalone_se <- es_from_linreg_b_se(
    linreg_b = b, linreg_b_se = se,
    n_sample = n, n_covariates = k
  )
  expect_equal(es_rp$es_crude, standalone_se$rp, tolerance = 1e-8)
})

# linreg_b reverse flag ----
test_that("linreg_b - reverse flag", {
  dat <- data.frame(
    linreg_b = c(1.5, -0.8),
    linreg_b_se = c(0.6, 0.3),
    n_sample = c(100, 50),
    n_covariates = c(2, 1)
  )

  es_rp <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_b_se", measure = "rp"), digits = 11)

  dat$reverse_linreg_b <- TRUE
  es_rp_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_b_se", measure = "rp"), digits = 11)

  expect_equal(es_rp$es_crude, -es_rp_r$es_crude, tolerance = 1e-10)
  expect_equal(es_rp$se_crude, es_rp_r$se_crude, tolerance = 1e-10)
})

# linreg_b_pval reverse flag ----
test_that("linreg_b_pval - reverse flag", {
  b <- c(1.5, -0.8); se <- c(0.6, 0.3)
  n <- c(100, 50); k <- c(2, 1)
  df <- n - k - 1
  pvals <- 2 * pt(-abs(b / se), df)

  dat <- data.frame(
    linreg_b = b,
    linreg_b_pval = pvals,
    n_sample = n,
    n_covariates = k
  )

  es_rp <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_b_pval", measure = "rp"), digits = 11)

  dat$reverse_linreg_b_pval <- TRUE
  es_rp_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy",
    hierarchy = "linreg_b_pval", measure = "rp"), digits = 11)

  expect_equal(es_rp$es_crude, -es_rp_r$es_crude, tolerance = 1e-10)
  expect_equal(es_rp$se_crude, es_rp_r$se_crude, tolerance = 1e-10)
})

# linreg_b hierarchy with linreg_t ----
test_that("linreg_b - hierarchy with linreg_t", {
  # Row 1: only linreg_b_se available; Row 2: only linreg_t available
  dat <- data.frame(
    linreg_b = c(1.5, NA),
    linreg_b_se = c(0.6, NA),
    linreg_t = c(NA, 2.5),
    n_sample = c(100, 100),
    n_covariates = c(2, 2)
  )

  es_rp <- summary(convert_df(dat, verbose = FALSE, measure = "rp"), digits = 11)
  expect_equal(es_rp$info_used_crude, c("linreg_b_se", "linreg_t"))

  # When linreg_t is also available in row 1, auto hierarchy should prefer linreg_t
  dat$linreg_t[1] <- 1.5 / 0.6  # same t as b/se
  es_rp2 <- summary(convert_df(dat, verbose = FALSE, measure = "rp"), digits = 11)
  expect_equal(nrow(es_rp2), 2)
  expect_true(all(!is.na(es_rp2$es_crude)))
})

# linreg_b mixed input types ----
test_that("linreg_b - mixed input types across rows", {
  b <- c(1.5, -0.8, 2.3)
  se <- c(0.6, 0.3, 0.9)
  n <- c(100, 50, 200)
  k <- c(2, 1, 5)
  df <- n - k - 1

  dat <- data.frame(
    # Row 1: b + se
    linreg_b = b,
    linreg_b_se = c(se[1], NA, NA),
    # Row 2: b + CI
    linreg_b_ci_lo = c(NA, b[2] - qt(.975, df[2]) * se[2], NA),
    linreg_b_ci_up = c(NA, b[2] + qt(.975, df[2]) * se[2], NA),
    # Row 3: b + pval
    linreg_b_pval = c(NA, NA, 2 * pt(-abs(b[3] / se[3]), df[3])),
    n_sample = n,
    n_covariates = k
  )

  es_rp <- summary(convert_df(dat, verbose = FALSE, measure = "rp"), digits = 11)
  expect_equal(nrow(es_rp), 3)
  expect_true(all(!is.na(es_rp$es_crude)))
  expect_equal(es_rp$info_used_crude, c("linreg_b_se", "linreg_b_ci", "linreg_b_pval"))
})

# linreg_b auto hierarchy ----
test_that("linreg_b - auto hierarchy (es_selected = 'auto')", {
  dat <- data.frame(
    linreg_b = c(1.5, -0.8),
    linreg_b_se = c(0.6, 0.3),
    n_sample = c(100, 50),
    n_covariates = c(2, 1)
  )

  es_rp <- summary(convert_df(dat, verbose = FALSE, measure = "rp"), digits = 11)
  expect_equal(nrow(es_rp), 2)
  expect_true(all(!is.na(es_rp$es_crude)))
})
