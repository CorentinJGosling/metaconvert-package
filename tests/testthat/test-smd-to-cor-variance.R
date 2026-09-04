## .smd_to_cor(): the d -> r/z conversion must (a) preserve the point estimate a crude
## two-group study would report, and (b) carry the ACTUAL precision of the supplied d.
##
## Background. r is a deterministic monotone function of d given the arm sizes, so
## Var(r) must inherit Var(d). The default "viechtbauer" branch uses Soper's closed
## form, which is a function of n, p and r only; to leading order that closed form IS
## the delta propagation of the CRUDE d variance, so applying it verbatim to a d whose
## precision differs from the crude design (ANCOVA, pre-post, Glass, user-supplied SE)
## silently reports the crude study's precision instead. The implementation therefore
## rescales it by vd / vd_crude, which is exactly 1 on crude rows.

## ---------------------------------------------------------------------------
## 1. Crude rows must be untouched -- and still bit-identical to metafor RBIS
## ---------------------------------------------------------------------------
test_that("crude d -> r matches metafor measure = 'RBIS' bit-for-bit", {
  skip_if_not_installed("metafor")
  for (nn in list(c(60, 60), c(20, 60), c(15, 45))) {
    for (d in c(0.2, 0.7, 1.3)) {
      e  <- es_from_cohen_d(cohen_d = d, n_exp = nn[1], n_nexp = nn[2])
      mf <- metafor::escalc(measure = "RBIS", di = d, n1i = nn[1], n2i = nn[2])
      lab <- sprintf("n = %d/%d, d = %.1f", nn[1], nn[2], d)
      expect_equal(e$r,      as.numeric(mf$yi), tolerance = 1e-12, info = lab)
      expect_equal(e$r_se^2, as.numeric(mf$vi), tolerance = 1e-12, info = lab)
    }
  }
})

test_that("the crude precision ratio is exactly 1 (no drift on ordinary rows)", {
  # vd_crude inside .smd_to_cor must be the same expression .es_from_d() uses
  for (nn in list(c(30, 30), c(12, 48))) {
    for (d in c(0, 0.5, 1.1)) {
      e <- es_from_cohen_d(cohen_d = d, n_exp = nn[1], n_nexp = nn[2])
      expect_equal(e$d_se^2,
                   (nn[1] + nn[2]) / (nn[1] * nn[2]) + d^2 / (2 * (nn[1] + nn[2])),
                   tolerance = 1e-14)
    }
  }
})

## ---------------------------------------------------------------------------
## 2. The r point estimate must not depend on the source study's covariate count
## ---------------------------------------------------------------------------
test_that("r and z do not depend on n_cov_ancova", {
  crude <- es_from_cohen_d(cohen_d = 0.7, n_exp = 20, n_nexp = 20)
  for (q in c(0, 1, 3, 5)) {
    adj <- es_from_cohen_d_adj(cohen_d_adj = 0.7, cov_outcome_r = 0.6,
                               n_cov_ancova = q, n_exp = 20, n_nexp = 20)
    # identical marginal d -> identical marginal correlation, whatever q
    expect_equal(adj$d, crude$d, tolerance = 1e-12, info = paste("q =", q))
    expect_equal(adj$r, crude$r, tolerance = 1e-12, info = paste("q =", q))
    expect_equal(adj$z, crude$z, tolerance = 1e-12, info = paste("q =", q))
  }
})

test_that("r reproduces the finite-sample point-biserial identity at df = N - 2", {
  # r_pb = d / sqrt(d^2 + h), h = (N-2)(1/n1 + 1/n2); r_bis = sqrt(pq)/f * r_pb
  n1 <- 25; n2 <- 35; d <- 0.6
  N <- n1 + n2; p <- n1 / N
  h <- (N - 2) / n1 + (N - 2) / n2
  expected <- sqrt(p * (1 - p)) / dnorm(qnorm(p, lower.tail = FALSE)) *
    d / sqrt(d^2 + h)
  expect_equal(es_from_cohen_d(cohen_d = d, n_exp = n1, n_nexp = n2)$r,
               expected, tolerance = 1e-12)
})

## ---------------------------------------------------------------------------
## 3. r_se / z_se must track d_se exactly
## ---------------------------------------------------------------------------
test_that("ANCOVA r_se and z_se carry the Cooper eq. 12.26 (1 - R^2) shrink", {
  # hold d fixed at 0.5 while varying R, so the point estimate (and hence Soper's
  # closed form) is constant and the whole change is the precision ratio
  ref <- NULL
  for (rr in c(0, 0.5, 0.7, 0.9)) {
    e <- es_from_ancova_md_sd(ancova_md = 0.5, ancova_md_sd = 1 * sqrt(1 - rr^2),
                              cov_outcome_r = rr, n_cov_ancova = 1,
                              n_exp = 60, n_nexp = 60)
    expect_equal(e$d, 0.5, tolerance = 1e-10)
    if (is.null(ref)) { ref <- e; next }
    lab <- paste("cov_outcome_r =", rr)
    expect_equal(e$r_se / ref$r_se, e$d_se / ref$d_se, tolerance = 1e-10, info = lab)
    expect_equal(e$z_se / ref$z_se, e$d_se / ref$d_se, tolerance = 1e-10, info = lab)
    # and it must actually shrink, not stay frozen
    expect_lt(e$r_se, ref$r_se)
  }
})

test_that("pre-post rows carry their own precision into r and z", {
  pp <- es_from_means_sd_pre_post(
    n_exp = 60, n_nexp = 60,
    mean_pre_exp = 0, mean_exp = 0.4, mean_pre_sd_exp = 1, mean_sd_exp = 1,
    mean_pre_nexp = 0, mean_nexp = 0, mean_pre_sd_nexp = 1, mean_sd_nexp = 1,
    r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8)
  indep <- es_from_cohen_d(cohen_d = pp$d, n_exp = 60, n_nexp = 60)

  expect_equal(pp$r, indep$r, tolerance = 1e-10)          # same d -> same r
  expect_lt(pp$d_se, indep$d_se)                          # but far more precise
  expect_equal(pp$r_se / indep$r_se, pp$d_se / indep$d_se, tolerance = 1e-10)
  expect_equal(pp$z_se / indep$z_se, pp$d_se / indep$d_se, tolerance = 1e-10)
})

test_that("scaling preserves the variance-stabilising relation between r and z", {
  # z is built so that (dz/dr)^2 * vr is constant; scaling both by the same factor
  # must leave vz/vr unchanged relative to a crude row with the same d
  crude <- es_from_cohen_d(cohen_d = 0.5, n_exp = 60, n_nexp = 60)
  adj   <- es_from_ancova_md_sd(ancova_md = 0.5, ancova_md_sd = sqrt(1 - 0.7^2),
                                cov_outcome_r = 0.7, n_cov_ancova = 1,
                                n_exp = 60, n_nexp = 60)
  expect_equal(adj$z_se^2 / adj$r_se^2, crude$z_se^2 / crude$r_se^2, tolerance = 1e-10)
})

## ---------------------------------------------------------------------------
## 4. lipsey_cooper branch
## ---------------------------------------------------------------------------
test_that("lipsey_cooper still propagates vd and stays Fisher-consistent", {
  e <- es_from_ancova_md_sd(ancova_md = 0.5, ancova_md_sd = sqrt(1 - 0.7^2),
                            cov_outcome_r = 0.7, n_cov_ancova = 1,
                            n_exp = 60, n_nexp = 60, smd_to_cor = "lipsey_cooper")
  p <- 0.5; a <- 1 / (p * (1 - p))
  expect_equal(e$z_se^2, e$d_se^2 / (e$d^2 + a), tolerance = 1e-12)
  expect_equal(e$z_se^2, e$r_se^2 / (1 - e$r^2)^2, tolerance = 1e-8)
})

test_that("lipsey_cooper r CI is the back-transformed z interval", {
  # Formerly r +/- qt(.975, N - 2 - q) r_se on the r scale -- the one package-computed r
  # interval still built that way, and one that could leave [-1, 1] at small n. It is
  # now tanh() of the z interval like every other correlation route.
  n1 <- 15; n2 <- 15; q <- 5
  e <- es_from_ancova_md_sd(ancova_md = 1, ancova_md_sd = 1, cov_outcome_r = 0.3,
                            n_cov_ancova = q, n_exp = n1, n_nexp = n2,
                            smd_to_cor = "lipsey_cooper")
  expect_equal(c(e$r_ci_lo, e$r_ci_up), tanh(c(e$z_ci_lo, e$z_ci_up)), tolerance = 1e-12)
  expect_equal(c(e$z_ci_lo, e$z_ci_up), e$z + c(-1, 1) * qnorm(.975) * e$z_se, tolerance = 1e-12)
  expect_true(abs(e$r_ci_lo) < 1 && abs(e$r_ci_up) < 1)
})

## ---------------------------------------------------------------------------
## 5. Calibration: the r interval must actually cover at its nominal rate
## ---------------------------------------------------------------------------
test_that("ANCOVA r and z 95% CIs are calibrated", {
  skip_on_cran()
  reps <- 1500; n1 <- n2 <- 60
  set.seed(20240)
  M <- matrix(NA_real_, reps, 6)
  for (i in seq_len(reps)) {
    x <- rnorm(n1 + n2)                       # random covariate: no construction bias
    grp <- rep(c(1, 0), c(n1, n2))
    y <- 0.8 * grp + x + rnorm(n1 + n2)
    fit <- lm(y ~ grp + x); s <- summary(fit)$sigma
    r_w <- sqrt(max(0, 1 - s^2 / summary(lm(y ~ grp))$sigma^2))
    a <- es_from_ancova_md_sd(ancova_md = unname(coef(fit)[2]), ancova_md_sd = s,
                              cov_outcome_r = r_w, n_cov_ancova = 1,
                              n_exp = n1, n_nexp = n2)
    M[i, ] <- c(a$r, a$r_ci_lo, a$r_ci_up, a$z, a$z_ci_lo, a$z_ci_up)
  }
  # population marginal biserial r implied by the design
  d_true <- 0.8 / sqrt(1 + 1)
  p <- 0.5; h <- (n1 + n2 - 2) * (1 / n1 + 1 / n2)
  r_true <- sqrt(p * (1 - p)) / dnorm(qnorm(p, lower.tail = FALSE)) *
    d_true / sqrt(d_true^2 + h)

  cov_r <- mean(M[, 2] <= r_true & r_true <= M[, 3])
  cov_z <- mean(M[, 5] <= mean(M[, 4]) & mean(M[, 4]) <= M[, 6])
  # before the precision-ratio fix these ran at ~0.99 (SE ~1.38x too large)
  expect_gt(cov_r, 0.925); expect_lt(cov_r, 0.972)
  expect_gt(cov_z, 0.925); expect_lt(cov_z, 0.972)
})
