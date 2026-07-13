# =============================================================================
# Calibration of the POOLED two-group pre/post variances.
#
# Published references exist for pooled d_z (= escalc(SMD,"LS") on change scores)
# and pooled d_av (Bonett 2008 eq. 19); d_rm follows from d_z by the exact rescaling
# d_rm = d_z*sqrt(2(1-r)). The robust-bonett departure from Viechtbauer's homoscedastic
# vi under heteroscedasticity has no closed-form two-group reference, so Monte-Carlo
# calibration is part of its specification: those tests ARE the evidence.
#
# Each branch follows metafor's LS *pattern* (an empirical leading term plus a g^2 term
# over 2N), with J = J(nu) at the standardizer df and Var(d) = Var(g)/J^2 (an exact
# identity). The leading term is heteroscedasticity-robust (built from the EMPIRICAL
# pooled change SD, not the 2*sigma^2*(1-r) identity that assumes SD_pre = SD_post).
# =============================================================================

# Slow (Monte Carlo): only run when explicitly asked.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {

# -----------------------------------------------------------------------------
# 1. Bit-exact agreement with metafor for the one branch that HAS a reference
# -----------------------------------------------------------------------------
test_that("pooled morris_dz is BIT-EXACT with escalc(measure = 'SMD') on change scores", {
  skip_if_not_installed("metafor")
  # Once the change SD is pooled across arms, the estimator IS an
  # independent-groups Hedges g on the change scores. metaConvert's .d_j is
  # identical to metafor's cm(), so agreement should be exact, not approximate.
  grid <- expand.grid(n1 = c(10, 29, 90), n2 = c(10, 34, 90))
  for (i in seq_len(nrow(grid))) {
    n1 <- grid$n1[i]; n2 <- grid$n2[i]
    res <- es_from_mean_change_sd(
      n_exp = n1, n_nexp = n2,
      mean_change_exp = -12.4, mean_change_sd_exp = 6.52,
      mean_change_nexp = -3.5, mean_change_sd_nexp = 10.67,
      r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6,
      pre_post_to_smd = "morris_dz", pool_sd = TRUE
    )
    ref <- metafor::escalc(
      measure = "SMD", vtype = "LS",
      m1i = -12.4, sd1i = 6.52, n1i = n1,
      m2i = -3.5, sd2i = 10.67, n2i = n2
    )
    expect_equal(res$g, as.numeric(ref$yi), tolerance = 1e-12,
                 label = sprintf("g, n = %d/%d", n1, n2))
    expect_equal(res$g_se^2, as.numeric(ref$vi), tolerance = 1e-12,
                 label = sprintf("vi, n = %d/%d", n1, n2))
  }
})

# -----------------------------------------------------------------------------
# 2. Monte-Carlo calibration of all four pooled branches
# -----------------------------------------------------------------------------
# Simulates bivariate-normal (pre, post) data, computes each estimator and its
# package variance, and checks that E[vi] tracks the TRUE sampling variance and
# that the 95% CI covers. Tolerances are deliberately loose (bias within 10%,
# coverage within [0.93, 0.97]) so the test guards against a real regression --
# e.g. reinstating the leading J^2 (which costs ~9% at n = 10/arm) or the
# homoscedastic numerator (which costs ~43% on bonett at q = 0.64) -- without
# being flaky.
mc_calibrate <- function(n1, n2, r, delta_raw, sd_pre, sd_post, method, nsim = 4000) {
  S <- matrix(c(sd_pre^2, r * sd_pre * sd_post,
                r * sd_pre * sd_post, sd_post^2), 2)
  g <- numeric(nsim); vi <- numeric(nsim)
  for (i in seq_len(nsim)) {
    A <- MASS::mvrnorm(n1, c(0, delta_raw), S)
    B <- MASS::mvrnorm(n2, c(0, 0), S)
    res <- es_from_means_sd_pre_post(
      n_exp = n1, n_nexp = n2,
      mean_pre_exp = mean(A[, 1]), mean_exp = mean(A[, 2]),
      mean_pre_sd_exp = sd(A[, 1]), mean_sd_exp = sd(A[, 2]),
      mean_pre_nexp = mean(B[, 1]), mean_nexp = mean(B[, 2]),
      mean_pre_sd_nexp = sd(B[, 1]), mean_sd_nexp = sd(B[, 2]),
      r_pre_post_exp = r, r_pre_post_nexp = r,
      pre_post_to_smd = method, pool_sd = TRUE
    )
    g[i] <- res$g; vi[i] <- res$g_se^2
  }
  ok <- is.finite(g) & is.finite(vi) & vi > 0
  g <- g[ok]; vi <- vi[ok]
  list(
    bias = mean(vi) / var(g) - 1,
    coverage = mean(abs(g - mean(g)) <= qnorm(.975) * sqrt(vi))
  )
}

test_that("all four pooled variances are calibrated under HOMOSCEDASTICITY", {
  skip_if_not_installed("MASS")
  set.seed(20260712)
  for (method in c("bonett", "morris_dz", "morris_drm", "morris_dav")) {
    out <- mc_calibrate(30, 30, r = 0.5, delta_raw = 0.8,
                        sd_pre = 1, sd_post = 1, method = method)
    expect_lt(abs(out$bias), 0.10,
              label = sprintf("|E[vi]/Var(g) - 1| for %s (%.1f%%)", method, 100 * out$bias))
    expect_gt(out$coverage, 0.93, label = sprintf("coverage for %s", method))
    expect_lt(out$coverage, 0.97, label = sprintf("coverage for %s", method))
  }
})

test_that("all four pooled variances stay calibrated under SD HETEROSCEDASTICITY", {
  skip_if_not_installed("MASS")
  # SD_pre / SD_post = 0.64 -- the median in the authors' own PETRA data, and the
  # regime metaConvert's V18 flag warns about. The OLD homoscedastic numerator
  # understated Var(g_bonett) by ~43% here and dropped coverage to 0.86; this is
  # the test that pins the fix.
  set.seed(20260713)
  for (method in c("bonett", "morris_dz", "morris_drm", "morris_dav")) {
    out <- mc_calibrate(30, 30, r = 0.5, delta_raw = 0.8,
                        sd_pre = 0.64, sd_post = 1.0, method = method)
    expect_lt(abs(out$bias), 0.10,
              label = sprintf("|E[vi]/Var(g) - 1| for %s under q = 0.64 (%.1f%%)",
                              method, 100 * out$bias))
    expect_gt(out$coverage, 0.93,
              label = sprintf("coverage for %s under q = 0.64", method))
  }
})

# -----------------------------------------------------------------------------
# 3. Structural properties the formulas must have
# -----------------------------------------------------------------------------
test_that("Var(g) = J^2 * Var(d) holds exactly on every pooled branch", {
  # No leading J^2 on var_g (the LS convention: var_g = leading_term + g^2/(2N)), but
  # because g = J*d with J a deterministic constant, var_d MUST equal var_g/J^2 exactly.
  # The earlier code violated this (it set var_d to a separate "d-scale twin"); this
  # test guards the restored identity. If someone re-introduces a leading J^2 on var_g
  # (the LS2 pattern var_g = J^2 * (leading + d^2/2N)) the g/d relation still holds but
  # the calibration tests above break -- so both are needed.
  for (method in c("bonett", "morris_dz", "morris_drm", "morris_dav")) {
    r <- es_from_means_sd_pre_post(
      n_exp = 12, n_nexp = 11,
      mean_pre_exp = 10, mean_exp = 16, mean_pre_sd_exp = 4, mean_sd_exp = 5,
      mean_pre_nexp = 10, mean_nexp = 12, mean_pre_sd_nexp = 4, mean_sd_nexp = 5,
      r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
      pre_post_to_smd = method, pool_sd = TRUE
    )
    J <- r$g / r$d
    expect_equal(r$g_se^2, J^2 * r$d_se^2, tolerance = 1e-10,
                 label = sprintf("var_g = J^2 * var_d for %s", method))
    expect_equal(r$g, J * r$d, tolerance = 1e-10,
                 label = sprintf("g = J*d for %s", method))
  }
})

test_that("morris_dav uses the Cousineau (2020) effective df, not N - 2", {
  # nu = 2m/(1+r^2) -- the same modified df metafor SMCRP uses per arm, and the
  # same one metaConvert's own single-group dav branch already uses.
  n1 <- 12; n2 <- 11; r <- 0.5
  m <- n1 + n2 - 2
  nu <- 2 * m / (1 + r^2)
  res <- es_from_means_sd_pre_post(
    n_exp = n1, n_nexp = n2,
    mean_pre_exp = 10, mean_exp = 16, mean_pre_sd_exp = 4, mean_sd_exp = 5,
    mean_pre_nexp = 10, mean_nexp = 12, mean_pre_sd_nexp = 4, mean_sd_nexp = 5,
    r_pre_post_exp = r, r_pre_post_nexp = r,
    pre_post_to_smd = "morris_dav", pool_sd = TRUE
  )
  expect_equal(res$g / res$d, metaConvert:::.d_j(nu), tolerance = 1e-10)
  # and NOT the unmodified df
  expect_false(isTRUE(all.equal(res$g / res$d, metaConvert:::.d_j(m), tolerance = 1e-6)))
})

test_that("pooled drm = pooled dz * sqrt(2(1-r)) in both g and SE", {
  r <- 0.6
  args <- list(
    n_exp = 29, n_nexp = 34,
    mean_change_exp = -12.4, mean_change_sd_exp = 6.52,
    mean_change_nexp = -3.5, mean_change_sd_nexp = 10.67,
    r_pre_post_exp = r, r_pre_post_nexp = r, pool_sd = TRUE
  )
  dz <- do.call(es_from_mean_change_sd, c(args, pre_post_to_smd = "morris_dz"))
  drm <- do.call(es_from_mean_change_sd, c(args, pre_post_to_smd = "morris_drm"))
  k <- sqrt(2 * (1 - r))
  expect_equal(drm$g, dz$g * k, tolerance = 1e-10)
  expect_equal(drm$g_se, dz$g_se * k, tolerance = 1e-10)
})

} # end NOT_CRAN gate
