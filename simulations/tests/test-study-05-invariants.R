## =============================================================================
## Roadmap item 3.1 -- study 05's vanderweele_rr_squared computed 2*log(OR)
## instead of 2*log(RR), and (item 1.9's mirror) built the symmetric interval the
## source paper rules out.
##
## THE COLLISION. gen_2x2_rr() (study 04) returns theta_sample on the log RR
## scale. gen_2x2_or() overwrites it with the log OR, correctly -- that is this
## study's target. The candidate then read theta_sample and doubled it, so it
## squared the ODDS ratio. Intent is not ambiguous: the comment says "OR ~= RR^2"
## and the SE on the next line is 2 * logrr_se, which is only right for the RR
## reading.
##
## WHY IT IS WORSE THAN ONE BAD ROW. Scored against theta_sample = log OR, the
## bias is 2*logOR - logOR = log OR EXACTLY, so the cell's mean bias is just the
## mean |log OR| of the design -- an artifact that tracks the grid rather than the
## method. Correcting it does not merely shrink the number, it inverts the
## baseline-risk gradient the write-up interpreted.
##
## The fix gives the log RR its own column, logrr_sample, which gen_2x2_or() sets
## before the overwrite and never reuses.
## =============================================================================

z975 <- stats::qnorm(.975)

.cond05 <- function(br = 0.30, rr = 2, n = 300, br_guess = br)
  list(n = n, rr = rr, br = br, p_exp = 0.5, br_guess = br_guess)

## --- the two scales are kept apart -----------------------------------------

test_that("gen_2x2_or() keeps the log RR in its own column and targets the log OR", {
  set.seed(3101)
  d <- gen_2x2_or(.cond05(), 500)

  expect_true("logrr_sample" %in% names(d))
  expect_equal(d$theta_sample, log(d$or))          # target: log OR
  expect_equal(d$logrr_sample, log(d$rr))          # kept apart: log RR
  # and they really are different scales on this data
  expect_false(isTRUE(all.equal(d$theta_sample, d$logrr_sample)))
  # the RR column is the one gen_2x2_rr() produced, not something recomputed
  expect_equal(exp(d$logrr_sample), (d$a / (d$a + d$b)) / (d$c / (d$c + d$d)))
})

test_that("the population target is the log OR implied by (rr, br)", {
  d <- gen_2x2_or(.cond05(br = 0.30, rr = 2), 5)
  r1 <- 0.60; r0 <- 0.30
  expect_equal(unique(d$theta_pop), log((r1 / (1 - r1)) / (r0 / (1 - r0))))
})

test_that("gen_2x2_or()'s invariant fires if the two scales are ever merged again", {
  # Simulate the regression: hand the checker a frame whose logrr_sample has been
  # overwritten with the log OR, and assert the assertion would catch it.
  set.seed(3102)
  d <- gen_2x2_or(.cond05(), 50)
  d_bad <- d; d_bad$logrr_sample <- log(d_bad$or)
  expect_false(isTRUE(all.equal(d_bad$logrr_sample, log(d_bad$rr))))
  expect_true(isTRUE(all.equal(d$logrr_sample, log(d$rr))))
})

## --- the candidate squares the RIGHT quantity -------------------------------

test_that("the candidate is 2*log(RR), not 2*log(OR)", {
  set.seed(3103)
  d <- gen_2x2_or(.cond05(), 400)
  got <- vanderweele_rr_squared(d$logrr_sample, d$logrr_se)

  expect_equal(got$logor, 2 * log(d$rr))
  expect_equal(exp(got$logor), d$rr^2)
  expect_false(isTRUE(all.equal(got$logor, 2 * log(d$or))))   # the old behaviour
  expect_equal(got$logor_se, 2 * d$logrr_se)                  # SE matches the RR reading
})

test_that("the estimator wired into the study reads logrr_sample", {
  # End-to-end through estimate_rr_to_or(), which is what run_study() calls -- so
  # this one needs the package's es_from_rr_*() routes, not just the study file.
  skip_if_not(isTRUE(get0("METACONVERT_AVAILABLE", ifnotfound = FALSE)),
              "metaConvert not loadable")
  set.seed(3104)
  cond <- .cond05()
  d <- gen_2x2_or(cond, 200)
  # rr_to_or = "grant" warns per row here (rr_ci_up * br_guess >= 1 -- the item-1.4
  # domain guard doing its job, and one of the things this study measures). Quiet
  # THIS call only; never wrap the suite itself, which hides real assertions.
  out <- suppressWarnings(estimate_rr_to_or(d, cond))
  vw <- out[out$method == "vanderweele_rr_squared", ]

  expect_equal(nrow(vw), 200L)
  expect_equal(vw$est, 2 * log(d$rr))
  expect_false(isTRUE(all.equal(vw$est, 2 * log(d$or))))
})

test_that("the buggy reading reproduces its own signature bias exactly", {
  # bias = 2*logOR - logOR = logOR. Pinning it makes the artifact identifiable in
  # any pre-fix aggregate still in circulation.
  set.seed(3105)
  d <- gen_2x2_or(.cond05(br = 0.30, rr = 2, n = 300), 5000)
  bias_buggy <- mean(2 * d$theta_sample - d$theta_sample)
  expect_equal(bias_buggy, mean(d$theta_sample))
  expect_equal(bias_buggy, mean(abs(d$theta_sample)), tolerance = 1e-8)  # all logOR > 0 here

  bias_fixed <- mean(vanderweele_rr_squared(d$logrr_sample, d$logrr_se)$logor -
                       d$theta_sample)
  expect_lt(abs(bias_fixed), abs(bias_buggy) / 5)
})

## --- the mirrored interval (item 1.9's counterpart) -------------------------

test_that("the RR->OR bias-ratio bound is the SQUARE of the OR->RR one", {
  expect_equal(vw_bias_ratio_sq(0.2, 0.8), 1.5625)
  expect_equal(vw_bias_ratio_sq(0.2, 0.8), vw_bias_ratio(0.2, 0.8)^2)
  for (w in c(0.1, 0.2, 0.3, 0.4)) {
    v <- (1 - 2 * w) / 2
    expect_equal(vw_bias_ratio_sq(w, 1 - w), 1 / (1 - 4 * v^2))
  }
})

test_that("that bound is exact, by closed form and by numerical maximisation", {
  # RR^2/OR = p1(1-p1) / (p0(1-p0)). Derived here, not quoted from the paper, so
  # it is checked rather than asserted.
  g <- expand.grid(p0 = seq(0.2, 0.8, by = 0.002), p1 = seq(0.2, 0.8, by = 0.002))
  rr <- g$p1 / g$p0
  or <- (g$p1 * (1 - g$p0)) / (g$p0 * (1 - g$p1))
  expect_equal(rr^2 / or, (g$p1 * (1 - g$p1)) / (g$p0 * (1 - g$p0)))
  expect_lte(max(rr^2 / or), 1.5625 + 1e-9)
  expect_gt(max(rr^2 / or), 1.5624)              # the bound is attained, not loose
  expect_gte(min(rr^2 / or), 1 / 1.5625 - 1e-9)
})

test_that("the interval is the transformed RR limits, widened -- not symmetric", {
  logrr <- log(1.8); se <- 0.25
  got <- vanderweele_rr_squared(logrr, se)

  rr_lo <- exp(logrr - z975 * se); rr_up <- exp(logrr + z975 * se)
  expect_equal(exp(got$logor_ci_lo), rr_lo^2 / 1.5625)
  expect_equal(exp(got$logor_ci_up), rr_up^2 * 1.5625)

  # not reconstructible from the reported SE (guards a future "tidy-up")
  expect_false(isTRUE(all.equal(got$logor_ci_lo, got$logor - z975 * got$logor_se)))
  expect_false(isTRUE(all.equal(got$logor_ci_up, got$logor + z975 * got$logor_se)))
  # wider than the symmetric interval by exactly 2*log(1.5625)
  naive_w <- 2 * z975 * got$logor_se
  expect_equal((got$logor_ci_up - got$logor_ci_lo) - naive_w, 2 * log(1.5625))
})

test_that("inside the band the widened interval covers and the symmetric one did not", {
  # Measured, not asserted: this is the claim the coverage column reports.
  set.seed(3106)
  cond <- .cond05(br = 0.50, rr = 0.5, n = 300)
  d <- gen_2x2_or(cond, 8000)
  tgt <- d$theta_pop[1]

  got <- vanderweele_rr_squared(d$logrr_sample, d$logrr_se)
  cov_wide <- mean(got$logor_ci_lo <= tgt & got$logor_ci_up >= tgt)

  naive_lo <- got$logor - z975 * got$logor_se
  naive_up <- got$logor + z975 * got$logor_se
  cov_naive <- mean(naive_lo <= tgt & naive_up >= tgt)

  expect_gte(cov_wide, 0.95)
  expect_lt(cov_naive, 0.95)
  expect_gt(cov_wide, cov_naive)
})

test_that("the candidate is vectorised and never returns a broken interval", {
  set.seed(3107)
  logrr <- stats::runif(300, log(0.2), log(5))
  se <- stats::runif(300, 0.05, 1)
  got <- vanderweele_rr_squared(logrr, se)
  expect_equal(nrow(got), 300L)
  expect_true(all(got$logor_ci_lo < got$logor))
  expect_true(all(got$logor_ci_up > got$logor))
  expect_true(all(is.finite(unlist(got))))
})

## --- stopifnot() fallback for a bare R (no testthat) ------------------------
if (!requireNamespace("testthat", quietly = TRUE)) {
  set.seed(3199)
  .d <- gen_2x2_or(.cond05(), 100)
  .g <- vanderweele_rr_squared(.d$logrr_sample, .d$logrr_se)
  stopifnot(
    isTRUE(all.equal(.d$theta_sample, log(.d$or))),
    isTRUE(all.equal(.d$logrr_sample, log(.d$rr))),
    isTRUE(all.equal(.g$logor, 2 * log(.d$rr))),
    !isTRUE(all.equal(.g$logor, 2 * log(.d$or))),
    isTRUE(all.equal(vw_bias_ratio_sq(0.2, 0.8), 1.5625))
  )
}
