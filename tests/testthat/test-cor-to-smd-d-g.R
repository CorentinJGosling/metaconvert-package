## `cor_to_smd` must report BOTH d and g, and report both accurately.
##
## The invariant that holds across every route is g = J * d, so a user asking for
## measure = "d" gets the uncorrected SMD and measure = "g" the bias-corrected one.
##
## Which quantity anchors the pair differs by route, and for "viechtbauer" that is
## a question about estimator bias rather than about the algebra of the transform.
## metafor::transf.rtod() is a population map with no df term, which invites the
## conclusion that its output is an uncorrected d -- but simulation says otherwise
## (bivariate normal, rho = 0.5, median split, true delta = 0.870126, nrep = 2e5;
## bias of the reported g against delta):
##
##   n      Hedges g from raw data   g = transf.rtod   g = J * transf.rtod
##   10           -0.004                  +0.033             -0.055
##   25           -0.000                  +0.012             -0.017
##  200           -0.000                  +0.002             -0.002
##
## transf.rtod(r_hat) belongs in g: E[r_hat] is biased down, which partly cancels
## the upward bias the pooled-SD denominator gives a directly computed d, so this
## route does not carry the full Hedges bias and applying J on top over-corrects.

test_that("g = J * d holds for every cor_to_smd route", {
  for (n in c(25, 50, 100, 300)) {
    J <- .d_j(n - 2)
    for (m in c("viechtbauer", "cooper")) {
      res <- es_from_pearson_r(pearson_r = 0.5, n_sample = n, cor_to_smd = m)
      expect_equal(res$g / res$d, J, tolerance = 1e-10, info = paste(m, "n =", n))
      expect_equal(res$g_se / res$d_se, J, tolerance = 1e-10, info = paste(m, "n =", n))
    }
  }
})

test_that("viechtbauer anchors g on transf.rtod; cooper anchors d on the point-biserial map", {
  for (r in c(-0.6, -0.2, 0, 0.3, 0.75)) {
    v <- es_from_pearson_r(pearson_r = r, n_sample = 200, cor_to_smd = "viechtbauer")
    expect_equal(v$g, metafor::transf.rtod(r), tolerance = 1e-12)
    expect_equal(v$d, metafor::transf.rtod(r) / .d_j(198), tolerance = 1e-12)

    cp <- es_from_pearson_r(pearson_r = r, n_sample = 200, cor_to_smd = "cooper")
    expect_equal(cp$d, 2 * r / sqrt(1 - r^2), tolerance = 1e-12)
  }
})

test_that("transf.rtod at a median split equals the analytic median-split SMD", {
  # Anchors the viechtbauer estimand to a quantity derived independently of the
  # package: dichotomising a standard bivariate normal at the median gives
  # 2*L*rho / sqrt(1 - L^2*rho^2) with L = dnorm(0)/0.5.
  L <- stats::dnorm(0) / 0.5
  for (rho in c(0.1, 0.3, 0.5, 0.7)) {
    analytic <- 2 * L * rho / sqrt(1 - L^2 * rho^2)
    expect_equal(metafor::transf.rtod(rho), analytic, tolerance = 1e-10,
                 info = paste("rho =", rho))
  }
})

test_that(".cor_to_smd_vec matches the per-row .cor_to_smd it replaced", {
  set.seed(99)
  N <- 300
  r <- runif(N, -0.9, 0.9); rse <- runif(N, 0.01, 0.15)
  n <- sample(20:500, N, TRUE); sdiv <- runif(N, 0.5, 2)
  uinc <- rep(1, N); utype <- sample(c("sd", "value"), N, TRUE)
  meth <- sample(c("viechtbauer", "cooper", "mathur"), N, TRUE)

  old <- t(mapply(.cor_to_smd, r = r, r_se = rse, n_sample = n, sd_iv = sdiv,
                  unit_increase_iv = uinc, unit_type = utype, cor_to_smd = meth))
  old <- cbind(unlist(old[, 1]), unlist(old[, 2]))
  new <- .cor_to_smd_vec(r, rse, uinc, sdiv, utype, n, meth)

  expect_identical(is.na(old), is.na(new))
  # The "viechtbauer" branch (conv.delta) is bit-identical, but the closed-form
  # "cooper" / "mathur" branches can land 1-2 ULP apart because R evaluates
  # length-1 and length-n arithmetic on different code paths -- and which rows
  # drift is compiler/ISA dependent, so tolerance = 0 here would be a platform
  # flake. Measured worst case 4.3e-16 relative.
  expect_equal(old[, 1], new[, 1], tolerance = 1e-12)
  expect_equal(old[, 2], new[, 2], tolerance = 1e-12)
})
