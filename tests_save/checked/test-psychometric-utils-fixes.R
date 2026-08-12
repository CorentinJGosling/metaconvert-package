# =============================================================================
# Coverage for the compute_sem() / reliability_change_score() / compute_sdc()
# audit fixes (AUDIT-FIX-TRACKER P14, P18, P19, P21; digest [30], [32], [33],
# [35]).
#
# Design rule for this file: expected values are hand-derived formulas written
# independently below (SEM = SD*sqrt(1-ICC); same-sample Var(SEM) =
# SEM^2/(2 df) with df = (n-1)(k-1) from the scaled chi-square distribution of
# MSE; exact chi-square CI [SEM*sqrt(df/qchisq(.975,df)),
# SEM*sqrt(df/qchisq(.025,df))]; general Lord change-score reliability) --
# never re-typed copies of the package source.
# =============================================================================

# -----------------------------------------------------------------------------
# P18: n_measurements = 1 makes the same-sample df = (n-1)(k-1) = 0. Pre-fix
# this silently emitted sem_se = Inf, CI [0, Inf], no warning.
# -----------------------------------------------------------------------------
test_that("P18: n_measurements = 1 yields NA SE/CI with a warning (was silent Inf)", {
  w <- testthat::capture_warnings(
    res <- compute_sem(sd = 10, icc = 0.85, n_sample = 50, n_measurements = 1)
  )
  expect_true(any(grepl("n_measurements", w)))
  # SEM point estimate does not depend on k -- must remain intact
  expect_equal(res$sem, 10 * sqrt(1 - 0.85), tolerance = 1e-10)
  expect_true(is.na(res$sem_se))
  expect_true(is.na(res$sem_ci_lo))
  expect_true(is.na(res$sem_ci_up))
})

# -----------------------------------------------------------------------------
# P18: icc > 1 pre-fix rode the icc >= 1 branch and emitted the incoherent
# pair (sem = NaN, sem_se = 0) -- a maximally confident SE around an
# impossible point estimate.
# -----------------------------------------------------------------------------
test_that("P18: icc > 1 yields NA sem and NA sem_se with a warning (was NaN with SE = 0)", {
  w <- testthat::capture_warnings(
    res <- compute_sem(sd = 10, icc = 1.05, n_sample = 50)
  )
  expect_true(any(grepl("ICC", w)))
  expect_true(is.na(res$sem))
  expect_false(is.nan(res$sem))        # a deliberate NA_real_, not an arithmetic NaN
  expect_true(is.na(res$sem_se))
  expect_false(isTRUE(res$sem_se == 0))
})

test_that("P18 regression: icc = 1 exactly stays the coherent degenerate case (0, 0)", {
  res <- compute_sem(sd = 10, icc = 1, n_sample = 50)
  expect_equal(res$sem, 0, tolerance = 1e-10)
  expect_equal(res$sem_se, 0, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# P19: the same-sample branch (icc_se missing) knows df = (n-1)(k-1); its CI
# must be the exact chi-square interval, not the symmetric Wald interval
# (which undercovers: 89.7% at n = 10, k = 2, analytic).
# -----------------------------------------------------------------------------
test_that("P19: same-sample CI is the exact chi-square interval (n = 10, k = 2)", {
  res <- compute_sem(sd = 10, icc = 0.85, n_sample = 10, n_measurements = 2)
  sem <- 10 * sqrt(1 - 0.85)
  df <- (10 - 1) * (2 - 1)
  # Hand oracle from the primary definition: MSE/sigma_e^2 ~ chi-square_df/df,
  # so a 95% interval for sigma_e = SEM_true is
  # [SEM*sqrt(df/qchisq(.975, df)), SEM*sqrt(df/qchisq(.025, df))].
  expect_equal(res$sem_ci_lo, sem * sqrt(df / qchisq(0.975, df)), tolerance = 1e-10)
  expect_equal(res$sem_ci_up, sem * sqrt(df / qchisq(0.025, df)), tolerance = 1e-10)
  # the delta-method SE itself is unchanged: Var(SEM) = SEM^2 / (2 df)
  expect_equal(res$sem_se, sem / sqrt(2 * df), tolerance = 1e-10)
})

test_that("P19: same-sample chi-square CI also used at larger df (n = 100, k = 2)", {
  res <- compute_sem(sd = 10, icc = 0.85, n_sample = 100)
  sem <- 10 * sqrt(1 - 0.85)
  df <- (100 - 1) * (2 - 1)
  expect_equal(res$sem_ci_lo, sem * sqrt(df / qchisq(0.975, df)), tolerance = 1e-10)
  expect_equal(res$sem_ci_up, sem * sqrt(df / qchisq(0.025, df)), tolerance = 1e-10)
})

test_that("P19: user-supplied icc_se keeps the truncated Wald CI (df unknown)", {
  res <- compute_sem(sd = 10, icc = 0.85, n_sample = 50, icc_se = 0.05)
  sem <- 10 * sqrt(1 - 0.85)
  # Hand oracle: bivariate delta method with Cov = 0,
  # Var(SEM) = (SD^2/(4(1-ICC))) * icc_se^2 + (1-ICC) * SD^2/(2(n-1))
  var_sem <- (10^2 / (4 * (1 - 0.85))) * 0.05^2 + (1 - 0.85) * 10^2 / (2 * (50 - 1))
  se <- sqrt(var_sem)
  expect_equal(res$sem_se, se, tolerance = 1e-10)
  expect_equal(res$sem_ci_lo, max(0, sem - qnorm(0.975) * se), tolerance = 1e-10)
  expect_equal(res$sem_ci_up, sem + qnorm(0.975) * se, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# P14: compute_sem's icc_se must be a RAW-scale ICC standard error, but
# es_from_icc()'s default (bonett) returns a ln(1-ICC)-scale SE under the same
# column name. A heuristic warning must fire when icc_se > (1 - icc) -- a
# raw-scale SE cannot plausibly exceed that -- and never alter the value.
# -----------------------------------------------------------------------------
test_that("P14: Bonett-scale icc_se (0.264 at icc = 0.85) triggers the scale warning", {
  w <- testthat::capture_warnings(
    res <- compute_sem(sd = 10, icc = 0.85, n_sample = 50, icc_se = 0.264)
  )
  expect_true(any(grepl("raw", w, ignore.case = TRUE)))
  # warning only -- the value must NOT be altered: SE still computed from 0.264
  var_sem <- (10^2 / (4 * (1 - 0.85))) * 0.264^2 + (1 - 0.85) * 10^2 / (2 * (50 - 1))
  expect_equal(res$sem_se, sqrt(var_sem), tolerance = 1e-10)
})

test_that("P14: a plausible raw-scale icc_se (0.0396) does not trigger the scale warning", {
  w <- testthat::capture_warnings(
    compute_sem(sd = 10, icc = 0.85, n_sample = 50, icc_se = 0.0396)
  )
  expect_length(w, 0)
})

# -----------------------------------------------------------------------------
# P21: reliability_change_score() returned a negative "reliability" silently
# when r_pre_post > reliability (population-impossible under CTT; feeds
# sqrt(negative) = NaN downstream in es_disattenuate).
# -----------------------------------------------------------------------------
test_that("P21: negative change-score reliability warns; value preserved", {
  w <- testthat::capture_warnings(
    res <- reliability_change_score(reliability = 0.5, r_pre_post = 0.7)
  )
  expect_true(any(grepl("negative", w, ignore.case = TRUE)))
  # Lord formula by hand: (0.5 - 0.7) / (1 - 0.7) = -2/3; value is preserved
  expect_equal(res$rel_change, (0.5 - 0.7) / (1 - 0.7), tolerance = 1e-10)
})

test_that("P21 regression: Lord equal-variance special case is unchanged", {
  # General Lord (1963): (s1^2 rel1 + s2^2 rel2 - 2 s1 s2 r12) /
  # (s1^2 + s2^2 - 2 s1 s2 r12); with s1 = s2, rel1 = rel2 = 0.85, r12 = 0.60
  # this reduces to (0.85 - 0.60) / (1 - 0.60) = 0.625.
  w <- testthat::capture_warnings(
    res <- reliability_change_score(reliability = 0.85, r_pre_post = 0.60)
  )
  expect_length(w, 0)
  expect_equal(res$rel_change, 0.625, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# compute_sdc() chaining regression: the documented chain from compute_sem()
# must be unaffected by the fixes (it consumes only sem and sem_se).
# -----------------------------------------------------------------------------
test_that("compute_sdc chaining from compute_sem is unchanged on a clean case", {
  sem_res <- compute_sem(sd = 10, icc = 0.85, n_sample = 100)
  res <- compute_sdc(sem = sem_res$sem, sem_se = sem_res$sem_se)
  m <- qnorm(0.975) * sqrt(2)
  sem <- 10 * sqrt(1 - 0.85)
  sem_se <- sem / sqrt(2 * (100 - 1) * (2 - 1))   # same-sample exact variance
  expect_equal(res$sdc, m * sem, tolerance = 1e-10)
  expect_equal(res$sdc_se, m * sem_se, tolerance = 1e-10)
  expect_equal(res$sdc_ci_lo, max(0, m * sem - qnorm(0.975) * m * sem_se),
               tolerance = 1e-10)
  expect_equal(res$sdc_ci_up, m * sem + qnorm(0.975) * m * sem_se,
               tolerance = 1e-10)
})
