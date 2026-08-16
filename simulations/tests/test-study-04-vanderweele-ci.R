## =============================================================================
## Roadmap item 1.9 -- study 04's vanderweele_sqrt_or() built an interval its own
## source paper forbids.
##
## VanderWeele (2020) Biometrics 76(3):746-752, p.748:
##
##   "because the square-root transformation ... is an approximation, it cannot be
##    applied directly to the confidence interval of the odds ratio to obtain 95%
##    coverage over repeated samples of the true risk ratio. However ... if the
##    square-root transformation of the lower limit of the odds ratio confidence
##    interval is divided by 1.25 and the square-root transformation of the upper
##    limit of the odds ratio confidence interval is multiplied by 1.25, then this
##    resulting transformed confidence interval will have at least 95% coverage of
##    the true risk ratio provided the outcome probabilities do indeed fall between
##    w = 20% and u = 80% ... The coverage of the transformed interval will in
##    general be conservative."
##
## The candidate used to build `logrr +/- qnorm(.975) * logor_se/2`, i.e. exactly
## the direct transformation the first sentence rules out. Its coverage column
## therefore described an operation the paper rejects. The POINT estimate was and
## remains sqrt(OR), which Corollary 1 (p.747) proves optimal -- so only the
## interval was ever at issue, and this file pins that distinction.
## =============================================================================

## qnorm/pnorm etc. come from stats; the runner has already sourced the study.
z975 <- stats::qnorm(.975)

## --- the constant is a formula, not a magic number -------------------------

test_that("the bias ratio is 1.25 exactly at the paper's [0.2, 0.8] band", {
  expect_equal(vw_bias_ratio(0.2, 0.8), 1.25)
  # 1/sqrt(1 - 4v^2) for [0.5 - v, 0.5 + v]; v = 0.3 -> 1/sqrt(0.64) = 1.25
  expect_equal(vw_bias_ratio(0.2, 0.8), 1 / sqrt(1 - 4 * 0.3^2))
  # a narrower band is a weaker widening; the limit at v -> 0 is no widening
  expect_lt(vw_bias_ratio(0.4, 0.6), 1.25)
  expect_gt(vw_bias_ratio(0.1, 0.9), 1.25)
  expect_equal(vw_bias_ratio(0.5 - 1e-9, 0.5 + 1e-9), 1, tolerance = 1e-12)
  # and it is only derived for a band symmetric about 0.5
  expect_error(vw_bias_ratio(0.2, 0.7))
  expect_error(vw_bias_ratio(0.3, 0.2))
  expect_error(vw_bias_ratio(0, 1))
})

## --- the point estimate is untouched ---------------------------------------

test_that("the point estimate is still sqrt(OR) (Corollary 1)", {
  or <- c(0.25, 0.5, 1, 2, 4, 16)
  got <- vanderweele_sqrt_or(or, rep(0.3, length(or)))
  expect_equal(exp(got$logrr), sqrt(or))
  expect_equal(got$logrr, log(or) / 2)
  # the paper's own worked example, p.749: OR 2.3 -> RR 1.5
  expect_equal(round(exp(vanderweele_sqrt_or(2.3, 0.2)$logrr), 1), 1.5)
})

test_that("the reported SE is still the delta-method SE of the estimate", {
  # The paper gives no variance. The SE column documents the point estimate, and
  # is deliberately NOT the basis of the interval -- performance() reads coverage
  # from ci_lo/ci_up and se_ratio from se, independently.
  got <- vanderweele_sqrt_or(or = c(2, 4), logor_se = c(0.3, 0.5))
  expect_equal(got$logrr_se, c(0.15, 0.25))
})

## --- the interval now follows the paper ------------------------------------

test_that("the bounds are the transformed OR limits, divided/multiplied by 1.25", {
  or <- 2.3; se <- 0.2
  or_lo <- exp(log(or) - z975 * se)
  or_up <- exp(log(or) + z975 * se)

  got <- vanderweele_sqrt_or(or, se)
  # natural scale, stated exactly as the paper states it
  expect_equal(exp(got$logrr_ci_lo), sqrt(or_lo) / 1.25)
  expect_equal(exp(got$logrr_ci_up), sqrt(or_up) * 1.25)
  # log scale
  expect_equal(got$logrr_ci_lo, log(or_lo) / 2 - log(1.25))
  expect_equal(got$logrr_ci_up, log(or_up) / 2 + log(1.25))
})

test_that("the paper's worked example reproduces its published interval", {
  # p.749: OR 2.3 (95% CI 1.5 to 3.4) -> RR 1.5, "a conservative 95% confidence
  # interval for the risk ratio estimate of 1.5 is (1.0 to 2.3)".
  lo <- sqrt(1.5) / 1.25
  up <- sqrt(3.4) * 1.25
  expect_equal(round(lo, 1), 1.0)
  expect_equal(round(up, 1), 2.3)
})

test_that("the forbidden symmetric interval is what changed, by exactly 2*log(1.25)", {
  or <- c(0.3, 1, 2.5, 9); se <- c(0.2, 0.35, 0.15, 0.6)
  got <- vanderweele_sqrt_or(or, se)

  naive_lo <- log(or) / 2 - z975 * se / 2      # the pre-fix bounds
  naive_up <- log(or) / 2 + z975 * se / 2

  expect_equal(naive_lo - got$logrr_ci_lo, rep(log(1.25), length(or)))
  expect_equal(got$logrr_ci_up - naive_up, rep(log(1.25), length(or)))
  # so the interval is strictly wider, by the same amount at every OR and SE
  expect_equal((got$logrr_ci_up - got$logrr_ci_lo) - (naive_up - naive_lo),
               rep(2 * log(1.25), length(or)))
  expect_true(all(got$logrr_ci_lo < got$logrr_ci_up))
})

test_that("the interval is NOT reconstructible from the reported SE", {
  # Guards the documented inconsistency: if someone 'tidies' the function by
  # rebuilding the CI from logrr_se, this fails.
  got <- vanderweele_sqrt_or(2.5, 0.4)
  expect_false(isTRUE(all.equal(got$logrr_ci_lo, got$logrr - z975 * got$logrr_se)))
  expect_false(isTRUE(all.equal(got$logrr_ci_up, got$logrr + z975 * got$logrr_se)))
})

test_that("the interval always brackets the point estimate and is vectorised", {
  set.seed(19)
  or <- exp(stats::runif(200, log(0.1), log(10)))
  se <- stats::runif(200, 0.05, 1.2)
  got <- vanderweele_sqrt_or(or, se)
  expect_equal(nrow(got), 200L)
  expect_true(all(got$logrr_ci_lo < got$logrr))
  expect_true(all(got$logrr > got$logrr_ci_lo))
  expect_true(all(got$logrr_ci_up > got$logrr))
  expect_true(all(is.finite(unlist(got))))
})

## --- does the fix do what the paper claims? --------------------------------

test_that("inside the [0.2, 0.8] band the fixed interval covers and the old one did not", {
  # This is the claim the coverage column exists to test, so measure it rather
  # than assert the arithmetic twice. Cell inside the guaranteed region:
  # p0 = br = 0.30, p1 = rr * br = 0.60 -- both within [0.2, 0.8].
  set.seed(1909)
  br <- 0.30; rr <- 2; n_arm <- 300; nrep <- 2000
  a <- stats::rbinom(nrep, n_arm, rr * br); c_ <- stats::rbinom(nrep, n_arm, br)
  b <- n_arm - a; d <- n_arm - c_
  keep <- pmin(a, b, c_, d) > 0
  a <- a[keep]; b <- b[keep]; c_ <- c_[keep]; d <- d[keep]

  or <- (a * d) / (b * c_)
  se <- sqrt(1/a + 1/b + 1/c_ + 1/d)
  target <- log(rr)                       # true log RR

  fixed <- vanderweele_sqrt_or(or, se)
  cov_fixed <- mean(fixed$logrr_ci_lo <= target & fixed$logrr_ci_up >= target)

  naive_lo <- log(or) / 2 - z975 * se / 2
  naive_up <- log(or) / 2 + z975 * se / 2
  cov_naive <- mean(naive_lo <= target & naive_up >= target)

  # The paper promises AT LEAST 95%, "in general conservative", inside the band.
  expect_gte(cov_fixed, 0.95)
  # And the operation it forbids does not deliver it -- which is precisely why
  # the pre-fix coverage column was measuring a straw man.
  expect_lt(cov_naive, 0.95)
  expect_gt(cov_fixed, cov_naive)
})

test_that("outside the band the guarantee is not claimed and need not hold", {
  # p0 = br = 0.05 is below w = 0.2, so Theorem 1's premise fails. Recorded as a
  # fact about the grid, not a defect: this is why the coverage column for this
  # arm must be read separately inside and outside the guaranteed region.
  set.seed(1910)
  br <- 0.05; rr <- 2; n_arm <- 300; nrep <- 2000
  a <- stats::rbinom(nrep, n_arm, rr * br); c_ <- stats::rbinom(nrep, n_arm, br)
  b <- n_arm - a; d <- n_arm - c_
  keep <- pmin(a, b, c_, d) > 0
  a <- a[keep]; b <- b[keep]; c_ <- c_[keep]; d <- d[keep]

  or <- (a * d) / (b * c_)
  se <- sqrt(1/a + 1/b + 1/c_ + 1/d)
  fixed <- vanderweele_sqrt_or(or, se)
  cov_out <- mean(fixed$logrr_ci_lo <= log(rr) & fixed$logrr_ci_up >= log(rr))

  # No assertion on the level -- only that the cell is genuinely outside the
  # premise, so that a reader cannot mistake it for a failure of the method.
  expect_lt(br, 0.2)
  expect_true(is.finite(cov_out))
  message(sprintf("    [outside band: br = %.2f, rr = %.1f] coverage = %.3f", br, rr, cov_out))
})

## --- the grid split the README has to respect ------------------------------

test_that("the guaranteed region of study 04's grid is the documented handful of cells", {
  g <- build_grid_04()
  in_band <- g$br >= 0.2 & g$br <= 0.8 & g$rr * g$br >= 0.2 & g$rr * g$br <= 0.8
  cells <- unique(g[in_band, c("br", "rr")])
  cells <- cells[order(cells$br, cells$rr), ]

  expect_equal(cells$br, c(0.30, 0.30, 0.30, 0.50, 0.50, 0.50))
  expect_equal(cells$rr, c(0.75, 1, 2, 0.5, 0.75, 1))
  # every br <= 0.15 cell fails the premise outright (p0 < 0.2)
  expect_true(all(!in_band[g$br <= 0.15]))
  # and the guaranteed region is a minority of the grid
  expect_lt(mean(in_band), 0.5)
})

## --- stopifnot() fallback for a bare R (no testthat) ------------------------
if (!requireNamespace("testthat", quietly = TRUE)) {
  .g <- vanderweele_sqrt_or(2.3, 0.2)
  stopifnot(
    isTRUE(all.equal(vw_bias_ratio(0.2, 0.8), 1.25)),
    isTRUE(all.equal(exp(.g$logrr), sqrt(2.3))),
    isTRUE(all.equal(exp(.g$logrr_ci_lo),
                     sqrt(exp(log(2.3) - z975 * 0.2)) / 1.25)),
    isTRUE(all.equal(exp(.g$logrr_ci_up),
                     sqrt(exp(log(2.3) + z975 * 0.2)) * 1.25))
  )
}
