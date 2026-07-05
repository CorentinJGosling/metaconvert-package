# Validate ANCOVA functions against lm() output
# This test simulates ANCOVA data, fits lm(), and checks that
# es_from_ancova_* functions reproduce the model's MD, SE, t-statistic,
# Cohen's d, and Cohen's d SE.

test_that("ANCOVA md/md_se/d/d_se match lm() output", {
  set.seed(42)
  n_exp <- 60
  n_nexp <- 50
  n <- n_exp + n_nexp

  # Simulate data: two groups, one covariate
  # Orthogonalize covariate to group so that lm SE = sigma * sqrt(1/n1+1/n2) exactly
  # This eliminates the design-matrix source of discrepancy
  group <- c(rep(1, n_exp), rep(0, n_nexp))
  covariate_raw <- rnorm(n, mean = 50, sd = 10)
  covariate <- covariate_raw
  covariate[group == 1] <- covariate_raw[group == 1] - mean(covariate_raw[group == 1]) + mean(covariate_raw)
  covariate[group == 0] <- covariate_raw[group == 0] - mean(covariate_raw[group == 0]) + mean(covariate_raw)

  # True effect = 5, covariate effect = 0.8
  y <- 5 * group + 0.8 * covariate + rnorm(n, sd = 8)

  # Fit ANCOVA
  fit <- lm(y ~ group + covariate)
  smry <- summary(fit)
  coefs <- smry$coefficients

  # Extract all quantities from the lm summary
  lm_md       <- coefs["group", "Estimate"]
  lm_md_se    <- coefs["group", "Std. Error"]
  lm_t        <- coefs["group", "t value"]
  lm_pval     <- coefs["group", "Pr(>|t|)"]
  lm_F        <- lm_t^2
  lm_resid_sd <- smry$sigma  # = sqrt(RSS / (n - 3))

  # With orthogonal covariate, lm SE = sigma * sqrt(1/n1 + 1/n2) exactly
  expect_equal(lm_md_se, lm_resid_sd * sqrt(1/n_exp + 1/n_nexp), tolerance = 1e-10)

  # Covariate-outcome correlation (multiple R from covariate-only model)
  fit_cov <- lm(y ~ covariate)
  cov_outcome_r <- sqrt(summary(fit_cov)$r.squared)

  # Crude (unadjusted) pooled SD
  sd_exp  <- sd(y[group == 1])
  sd_nexp <- sd(y[group == 0])
  crude_sd_pooled <- sqrt(((n_exp - 1) * sd_exp^2 + (n_nexp - 1) * sd_nexp^2) / (n - 2))

  # Adjusted group means (from the ANCOVA model at mean covariate)
  grand_mean_cov <- mean(covariate)
  adj_mean_exp  <- coef(fit)["(Intercept)"] + coef(fit)["group"] * 1 + coef(fit)["covariate"] * grand_mean_cov
  adj_mean_nexp <- coef(fit)["(Intercept)"] + coef(fit)["group"] * 0 + coef(fit)["covariate"] * grand_mean_cov

  # Adjusted MD from means (sanity check)
  expect_equal(unname(adj_mean_exp - adj_mean_nexp), lm_md)

  # Adjusted per-group residual SDs (from model residuals)
  resid_exp  <- residuals(fit)[group == 1]
  resid_nexp <- residuals(fit)[group == 0]
  adj_sd_exp  <- sqrt(sum(resid_exp^2) / (n_exp - 1))
  adj_sd_nexp <- sqrt(sum(resid_nexp^2) / (n_nexp - 1))

  # CI from lm
  df_ancova <- n - 2 - 1  # n - 2 groups - 1 covariate
  lm_ci_lo <- lm_md - qt(0.975, df_ancova) * lm_md_se
  lm_ci_up <- lm_md + qt(0.975, df_ancova) * lm_md_se

  # ---- Test es_from_ancova_md_sd ----
  res_md_sd <- es_from_ancova_md_sd(
    ancova_md = lm_md,
    ancova_md_sd = lm_resid_sd,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  # MD and MD SE
  expect_equal(res_md_sd$md, lm_md, tolerance = 1e-10)
  expect_equal(res_md_sd$md_se, lm_md_se, tolerance = 1e-10,
               label = "md_se = ancova_md_sd * sqrt(1/n1 + 1/n2) = lm SE")

  # Cohen's d and SE
  expected_d <- lm_md / (lm_resid_sd / sqrt(1 - cov_outcome_r^2))
  expect_equal(res_md_sd$d, expected_d, tolerance = 1e-10,
               label = "Cohen's d from ancova_md_sd")
  expect_false(is.na(res_md_sd$d_se),
               label = "d_se should not be NA")

  # ---- Test es_from_ancova_md_se ----
  res_md_se_fn <- es_from_ancova_md_se(
    ancova_md = lm_md,
    ancova_md_se = lm_md_se,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  expect_equal(res_md_se_fn$md, lm_md, tolerance = 1e-10)
  # Round-trip: SE -> SD -> SE should preserve the value
  expect_equal(res_md_se_fn$md_se, lm_md_se, tolerance = 1e-10,
               label = "md_se from ancova_md_se should round-trip exactly")
  expect_equal(res_md_se_fn$md / res_md_se_fn$md_se, lm_t, tolerance = 1e-10,
               label = "md/md_se reproduces lm t-stat (SE round-trip)")
  expect_false(is.na(res_md_se_fn$d))
  expect_false(is.na(res_md_se_fn$d_se))

  # ---- Test es_from_ancova_md_ci ----
  res_md_ci <- es_from_ancova_md_ci(
    ancova_md = lm_md,
    ancova_md_ci_lo = lm_ci_lo,
    ancova_md_ci_up = lm_ci_up,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  expect_equal(res_md_ci$md, lm_md, tolerance = 1e-10)
  expect_equal(res_md_ci$md_se, lm_md_se, tolerance = 1e-10,
               label = "md_se from CI should recover lm SE exactly")
  expect_equal(res_md_ci$d, res_md_se_fn$d, tolerance = 1e-10,
               label = "d from CI should match d from SE")
  expect_equal(res_md_ci$d_se, res_md_se_fn$d_se, tolerance = 1e-10,
               label = "d_se from CI should match d_se from SE")

  # ---- Test es_from_ancova_md_pval ----
  res_md_pval <- es_from_ancova_md_pval(
    ancova_md = lm_md,
    ancova_md_pval = lm_pval,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  expect_equal(res_md_pval$md, lm_md, tolerance = 1e-10)
  expect_equal(res_md_pval$md_se, lm_md_se, tolerance = 1e-10,
               label = "md_se from pval should recover lm SE exactly")
  expect_equal(res_md_pval$d, res_md_se_fn$d, tolerance = 1e-10,
               label = "d from pval should match d from SE")
  expect_equal(res_md_pval$d_se, res_md_se_fn$d_se, tolerance = 1e-10,
               label = "d_se from pval should match d_se from SE")

  # ---- Test es_from_ancova_means_sd_pooled_adj ----
  res_pooled_adj <- es_from_ancova_means_sd_pooled_adj(
    ancova_mean_exp = unname(adj_mean_exp),
    ancova_mean_nexp = unname(adj_mean_nexp),
    ancova_mean_sd_pooled = lm_resid_sd,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  expect_equal(res_pooled_adj$md, lm_md, tolerance = 1e-10)
  expect_equal(res_pooled_adj$md_se, lm_md_se, tolerance = 1e-10,
               label = "md_se from pooled_adj should use adjusted SD")
  expect_equal(res_pooled_adj$d, res_md_sd$d, tolerance = 1e-10,
               label = "d from pooled_adj should match d from md_sd")
  expect_equal(res_pooled_adj$d_se, res_md_sd$d_se, tolerance = 1e-10,
               label = "d_se from pooled_adj should match d_se from md_sd")

  # ---- Test es_from_ancova_means_sd_pooled_crude ----
  res_pooled_crude <- es_from_ancova_means_sd_pooled_crude(
    ancova_mean_exp = unname(adj_mean_exp),
    ancova_mean_nexp = unname(adj_mean_nexp),
    mean_sd_pooled = crude_sd_pooled,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  expect_equal(res_pooled_crude$md, lm_md, tolerance = 1e-10)
  expect_false(is.na(res_pooled_crude$d))
  expect_false(is.na(res_pooled_crude$d_se))
  expect_false(is.na(res_pooled_crude$md_se))

  # ---- Test es_from_ancova_means_sd ----
  res_means_sd <- es_from_ancova_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    ancova_mean_exp = unname(adj_mean_exp),
    ancova_mean_nexp = unname(adj_mean_nexp),
    ancova_mean_sd_exp = adj_sd_exp,
    ancova_mean_sd_nexp = adj_sd_nexp,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1
  )

  expect_equal(res_means_sd$md, lm_md, tolerance = 1e-10)
  # md_se should use the adjusted SDs directly (not crude SDs)
  expected_means_sd_se <- sqrt(adj_sd_exp^2/n_exp + adj_sd_nexp^2/n_nexp)
  expect_equal(res_means_sd$md_se, expected_means_sd_se, tolerance = 1e-10,
               label = "md_se should use adjusted per-group SDs")
  expect_false(is.na(res_means_sd$d))
  expect_false(is.na(res_means_sd$d_se))

  # ---- Test es_from_ancova_t vs es_from_ancova_means_sd ----
  # With orthogonal covariate, lm_t = lm_md / (sigma * sqrt(1/n1+1/n2)),
  # so the only remaining discrepancy between ancova_t and ancova_means_sd
  # is the pooling denominator:
  #   - sigma = sqrt(RSS / (n-3))       [lm uses n - p degrees of freedom]
  #   - ancova_means_sd pools with n-2   [(n1-1) + (n2-1), ignoring covariate df]
  # Exact ratio: d_from_means_sd / d_from_t = sqrt((n-2) / (n-3))
  res_t <- es_from_ancova_t(
    ancova_t = lm_t,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  # d from ancova_t: d = t * sqrt(1/n1 + 1/n2) * sqrt(1 - r^2)
  expected_d_from_t <- lm_t * sqrt(1/n_exp + 1/n_nexp) * sqrt(1 - cov_outcome_r^2)
  expect_equal(res_t$d, expected_d_from_t, tolerance = 1e-10,
               label = "Cohen's d from ancova_t matches formula")
  expect_false(is.na(res_t$d_se))

  # ancova_t uses sigma (df = n-3), ancova_means_sd pools per-group SDs (df = n-2)
  # The ratio is exactly sqrt((n-2)/(n-3))
  n_cov <- 1
  pooling_ratio <- sqrt((n - 2) / (n - 2 - n_cov))
  expect_equal(res_means_sd$d, res_t$d * pooling_ratio, tolerance = 1e-2,
               label = "d_means_sd / d_t = sqrt((n-2)/(n-3)) [pooling denominator]")

  # ---- Test es_from_ancova_f (using lm F-stat) ----
  res_f <- es_from_ancova_f(
    ancova_f = lm_F,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  # F = t^2, so ancova_f and ancova_t must give identical d and d_se
  expect_equal(res_f$d, res_t$d, tolerance = 1e-10,
               label = "Cohen's d from ancova_f should match d from ancova_t")
  expect_equal(res_f$d_se, res_t$d_se, tolerance = 1e-10,
               label = "d_se from ancova_f should match d_se from ancova_t")

  # ---- Cross-function d consistency ----
  # Functions sharing the same SD input should give identical d
  expect_equal(res_md_sd$d, res_pooled_adj$d, tolerance = 1e-10)
  expect_equal(res_md_sd$d_se, res_pooled_adj$d_se, tolerance = 1e-10)

  # With orthogonal covariate: ancova_t(lm_t) = ancova_md_sd(lm_md, sigma) exactly
  expect_equal(res_t$d, res_md_sd$d, tolerance = 1e-10,
               label = "d from ancova_t = d from md_sd (orthogonal covariate)")
  expect_equal(res_t$d_se, res_md_sd$d_se, tolerance = 1e-10,
               label = "d_se from ancova_t = d_se from md_sd (orthogonal covariate)")

  # F -> t -> d chain should be exact
  expect_equal(res_f$d, res_t$d, tolerance = 1e-10)
  expect_equal(res_f$d_se, res_t$d_se, tolerance = 1e-10)

  # Delegation chain (md_se -> md_sd, md_ci -> md_se -> md_sd, md_pval -> md_se -> md_sd)
  expect_equal(res_md_se_fn$d, res_md_ci$d, tolerance = 1e-10)
  expect_equal(res_md_se_fn$d, res_md_pval$d, tolerance = 1e-10)
  expect_equal(res_md_se_fn$d_se, res_md_ci$d_se, tolerance = 1e-10)
  expect_equal(res_md_se_fn$d_se, res_md_pval$d_se, tolerance = 1e-10)
})

test_that("adjusted md_se is smaller than crude md_se (ANCOVA precision gain)", {
  set.seed(123)
  n_exp <- 80
  n_nexp <- 80
  n <- n_exp + n_nexp

  # Strong covariate: high r should give big precision gain
  group <- c(rep(1, n_exp), rep(0, n_nexp))
  covariate <- rnorm(n, mean = 0, sd = 1)
  y <- 3 * group + 2.5 * covariate + rnorm(n, sd = 3)

  fit <- lm(y ~ group + covariate)
  smry <- summary(fit)
  coefs <- smry$coefficients
  lm_md       <- coefs["group", "Estimate"]
  lm_md_se    <- coefs["group", "Std. Error"]
  lm_resid_sd <- smry$sigma

  fit_cov <- lm(y ~ covariate)
  cov_outcome_r <- sqrt(summary(fit_cov)$r.squared)

  # Crude SE (no covariate adjustment)
  sd_exp  <- sd(y[group == 1])
  sd_nexp <- sd(y[group == 0])
  crude_sd_pooled <- sqrt(((n_exp - 1) * sd_exp^2 + (n_nexp - 1) * sd_nexp^2) / (n - 2))
  crude_md_se <- crude_sd_pooled * sqrt(1/n_exp + 1/n_nexp)

  # Adjusted SE from es_from_ancova_md_sd (using the fix)
  res <- es_from_ancova_md_sd(
    ancova_md = lm_md,
    ancova_md_sd = lm_resid_sd,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = 1,
    n_exp = n_exp, n_nexp = n_nexp
  )

  # The adjusted SE should be SMALLER than the crude SE
  expect_true(res$md_se < crude_md_se,
              label = "Adjusted md_se should be smaller than crude md_se (ANCOVA precision gain)")

  # md_se should equal ancova_md_sd * sqrt(1/n1 + 1/n2) exactly
  expect_equal(res$md_se, lm_resid_sd * sqrt(1/n_exp + 1/n_nexp), tolerance = 1e-10)

  # d and d_se should be present
  expect_false(is.na(res$d))
  expect_false(is.na(res$d_se))
  expect_true(res$d_se > 0)
})

test_that("es_from_ancova_means_ci uses model df (not per-group df) to recover SE", {
  set.seed(99)
  n_exp <- 60
  n_nexp <- 50
  n <- n_exp + n_nexp
  n_cov <- 1

  # Simulate ANCOVA data with orthogonal covariate
  group <- c(rep(1, n_exp), rep(0, n_nexp))
  covariate_raw <- rnorm(n, mean = 50, sd = 10)
  covariate <- covariate_raw
  covariate[group == 1] <- covariate_raw[group == 1] - mean(covariate_raw[group == 1]) + mean(covariate_raw)
  covariate[group == 0] <- covariate_raw[group == 0] - mean(covariate_raw[group == 0]) + mean(covariate_raw)
  y <- 5 * group + 0.8 * covariate + rnorm(n, sd = 8)

  # Fit ANCOVA
  fit <- lm(y ~ group + covariate)
  smry <- summary(fit)
  sigma <- smry$sigma

  # Covariate-outcome correlation
  fit_cov <- lm(y ~ covariate)
  cov_outcome_r <- sqrt(summary(fit_cov)$r.squared)

  # Adjusted group means
  grand_mean_cov <- mean(covariate)
  adj_mean_exp  <- coef(fit)["(Intercept)"] + coef(fit)["group"] * 1 + coef(fit)["covariate"] * grand_mean_cov
  adj_mean_nexp <- coef(fit)["(Intercept)"] + coef(fit)["group"] * 0 + coef(fit)["covariate"] * grand_mean_cov

  # Adjusted per-group residual SDs
  resid_exp  <- residuals(fit)[group == 1]
  resid_nexp <- residuals(fit)[group == 0]
  adj_sd_exp  <- sqrt(sum(resid_exp^2) / (n_exp - 1))
  adj_sd_nexp <- sqrt(sum(resid_nexp^2) / (n_nexp - 1))

  # Compute adjusted per-group SEs
  adj_se_exp  <- adj_sd_exp / sqrt(n_exp)
  adj_se_nexp <- adj_sd_nexp / sqrt(n_nexp)

  # Build CIs using model df (the correct df)
  df_model <- n - 2 - n_cov
  ci_lo_exp  <- unname(adj_mean_exp)  - qt(0.975, df_model) * adj_se_exp
  ci_up_exp  <- unname(adj_mean_exp)  + qt(0.975, df_model) * adj_se_exp
  ci_lo_nexp <- unname(adj_mean_nexp) - qt(0.975, df_model) * adj_se_nexp
  ci_up_nexp <- unname(adj_mean_nexp) + qt(0.975, df_model) * adj_se_nexp

  # ---- Test es_from_ancova_means_ci recovers the SE exactly ----
  res_ci <- es_from_ancova_means_ci(
    n_exp = n_exp, n_nexp = n_nexp,
    ancova_mean_exp = unname(adj_mean_exp),
    ancova_mean_nexp = unname(adj_mean_nexp),
    ancova_mean_ci_lo_exp = ci_lo_exp,
    ancova_mean_ci_up_exp = ci_up_exp,
    ancova_mean_ci_lo_nexp = ci_lo_nexp,
    ancova_mean_ci_up_nexp = ci_up_nexp,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = n_cov
  )

  # Compare to es_from_ancova_means_sd (the reference)
  res_sd <- es_from_ancova_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    ancova_mean_exp = unname(adj_mean_exp),
    ancova_mean_nexp = unname(adj_mean_nexp),
    ancova_mean_sd_exp = adj_sd_exp,
    ancova_mean_sd_nexp = adj_sd_nexp,
    cov_outcome_r = cov_outcome_r,
    n_cov_ancova = n_cov
  )

  # With model df, CI -> SE -> SD round-trip should be exact
  expect_equal(res_ci$md, res_sd$md, tolerance = 1e-10,
               label = "md from CI should match md from SD")
  expect_equal(res_ci$md_se, res_sd$md_se, tolerance = 1e-10,
               label = "md_se from CI (model df) should recover the original SE exactly")
  expect_equal(res_ci$d, res_sd$d, tolerance = 1e-10,
               label = "d from CI should match d from SD")
  expect_equal(res_ci$d_se, res_sd$d_se, tolerance = 1e-10,
               label = "d_se from CI should match d_se from SD")

  # ---- Verify that using per-group df (the old bug) would give a DIFFERENT result ----
  # Old buggy SE recovery: divide by t with n_i - 1 df
  old_se_exp  <- (ci_up_exp - ci_lo_exp) / (2 * qt(0.975, n_exp - 1))
  old_se_nexp <- (ci_up_nexp - ci_lo_nexp) / (2 * qt(0.975, n_nexp - 1))
  old_md_se <- sqrt(old_se_exp^2 * n_exp / n_exp + old_se_nexp^2 * n_nexp / n_nexp)

  # Per-group df is smaller -> larger t-quantile -> smaller recovered SE
  expect_true(old_se_exp < adj_se_exp,
              label = "Old buggy SE (per-group df) underestimates the true SE")
  expect_true(old_se_nexp < adj_se_nexp,
              label = "Old buggy SE (per-group df) underestimates the true SE")
})
