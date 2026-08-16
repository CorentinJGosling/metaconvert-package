## =============================================================================
## Roadmap item 3.3 -- performance() published cells built from a handful of
## replications while advertising n_valid in the hundreds.
##
## One row summarises one condition x method cell, but its columns rest on THREE
## different subsets of the replications:
##
##   bias / emp_se / rmse        replications whose POINT ESTIMATE is finite
##   mod_se / se_ratio           ... whose STANDARD ERROR is finite
##   coverage / ci_width         ... whose INTERVAL is finite
##
## Only the first count was reported. sum(ok_se) and sum(ok_ci) were computed,
## used as denominators, and thrown away -- so a row could say n_valid = 996 while
## its coverage came from 5 replications, and look in every published table exactly
## like a cell built from all 1000.
##
## Measured over the shipped aggregates before the fix: 54 of 19,129 cells
## published an SE ratio or coverage from < 100 replications, 30 from < 30, 4 from
## < 10, the worst from 3 (study 05 `grant`, n_valid 440). 124 cells had a CI count
## below half their n_valid. Concentrated in study 05's `grant` and study 03a's
## `phi (r)` / `phi (z)`.
##
## Both counts are now returned, and each family is withheld below min_valid
## (default 30) keyed on ITS OWN count. The point-estimate family has no floor by
## design: n_valid has always been printed beside bias, and no shipped cell reports
## a bias from fewer than 100 valid estimates (0 of 57,026).
## =============================================================================

## A cell of `n` replications in which the point estimate is always fine, the SE
## is finite for `n_se` of them and the interval for `n_ci`.
.cell <- function(n = 200, n_se = n, n_ci = n, target = 0, seed = 33) {
  set.seed(seed)
  est <- stats::rnorm(n, target, 1)
  se  <- rep(1, n);  se[seq_len(n - n_se)]  <- NA_real_
  lo  <- est - 1.96; up <- est + 1.96
  if (n_ci < n) { lo[seq_len(n - n_ci)] <- NA_real_; up[seq_len(n - n_ci)] <- NA_real_ }
  list(est = est, se = se, ci_lo = lo, ci_up = up, target = target)
}
.perf <- function(x, ...) performance(x$est, x$se, x$ci_lo, x$ci_up, x$target, ...)

## --- the counts themselves --------------------------------------------------

test_that("performance() reports a count for each of the three families", {
  p <- .perf(.cell(n = 200, n_se = 200, n_ci = 200), min_valid = 0)
  expect_true(all(c("n_valid", "n_valid_se", "n_valid_ci") %in% names(p)))
  expect_equal(p$n_sim, 200)
  expect_equal(p$n_valid, 200)
  expect_equal(p$n_valid_se, 200)
  expect_equal(p$n_valid_ci, 200)
})

test_that("the SE and CI counts track their own subsets, not n_valid", {
  p <- .perf(.cell(n = 1000, n_se = 5, n_ci = 5), min_valid = 0)
  expect_equal(p$n_valid, 1000)      # every point estimate is fine ...
  expect_equal(p$n_valid_se, 5)      # ... and this is what se_ratio rests on
  expect_equal(p$n_valid_ci, 5)
  # the exact shape of the shipped defect: n_valid in the hundreds, coverage on 5
  expect_gt(p$n_valid / p$n_valid_ci, 100)
})

test_that("the counts move independently of each other", {
  p <- .perf(.cell(n = 300, n_se = 300, n_ci = 40), min_valid = 0)
  expect_equal(p$n_valid_se, 300)
  expect_equal(p$n_valid_ci, 40)

  p2 <- .perf(.cell(n = 300, n_se = 40, n_ci = 300), min_valid = 0)
  expect_equal(p2$n_valid_se, 40)
  expect_equal(p2$n_valid_ci, 300)
})

test_that("a non-finite point estimate removes the replication from ALL THREE", {
  # ok_se and ok_ci are conjunctions with ok, so a dropped estimate cannot leave a
  # usable SE behind it.
  x <- .cell(n = 100)
  x$est[1:10] <- NA_real_
  p <- .perf(x, min_valid = 0)
  expect_equal(p$n_valid, 90)
  expect_equal(p$n_valid_se, 90)
  expect_equal(p$n_valid_ci, 90)
})

test_that("a negative or infinite SE is not counted as valid", {
  x <- .cell(n = 100)
  x$se[1:3] <- -1; x$se[4:5] <- Inf; x$se[6] <- NaN
  p <- .perf(x, min_valid = 0)
  expect_equal(p$n_valid, 100)
  expect_equal(p$n_valid_se, 94)
  expect_equal(p$n_valid_ci, 100)   # the interval is unaffected
})

test_that("a half-missing interval does not count", {
  x <- .cell(n = 100)
  x$ci_up[1:7] <- NA_real_          # lower bound still finite
  p <- .perf(x, min_valid = 0)
  expect_equal(p$n_valid_ci, 93)
  expect_equal(p$n_valid_se, 100)
})

test_that("the counts are reported even on the give-up path", {
  # n_valid < 2 returns early; before the fix that branch had no counts at all, so
  # rbind-ing it beside a normal row would have failed once the columns were added.
  x <- .cell(n = 50)
  x$est[3:50] <- NA_real_           # 2 estimates -> still early-return territory
  x$est[2] <- NA_real_              # ... now 1
  p <- .perf(x, min_valid = 0)
  expect_equal(p$n_valid, 1)
  expect_equal(p$n_valid_se, 1)
  expect_equal(p$n_valid_ci, 1)
  expect_true(is.na(p$bias))
  expect_equal(ncol(p), ncol(.perf(.cell(n = 50), min_valid = 0)))
})

## --- the suppression floor ---------------------------------------------------

test_that("below min_valid a family is withheld, and its count still explains why", {
  p <- .perf(.cell(n = 1000, n_se = 5, n_ci = 5))    # default min_valid = 30
  expect_true(is.na(p$se_ratio))
  expect_true(is.na(p$mod_se))
  expect_true(is.na(p$se_ratio_mcse))
  expect_true(is.na(p$coverage))
  expect_true(is.na(p$coverage_mcse))
  expect_true(is.na(p$ci_width))
  # the counts are NOT withheld -- the suppression must be explainable
  expect_equal(p$n_valid_se, 5)
  expect_equal(p$n_valid_ci, 5)
  # and the point-estimate family is untouched
  expect_false(is.na(p$bias))
  expect_false(is.na(p$emp_se))
  expect_false(is.na(p$rmse))
  expect_equal(p$n_valid, 1000)
})

test_that("each family is keyed on its OWN count", {
  # 1000 usable SEs, 5 usable intervals: keep se_ratio, drop coverage.
  p <- .perf(.cell(n = 1000, n_se = 1000, n_ci = 5))
  expect_false(is.na(p$se_ratio))
  expect_false(is.na(p$mod_se))
  expect_true(is.na(p$coverage))
  expect_true(is.na(p$ci_width))

  # and the mirror
  p2 <- .perf(.cell(n = 1000, n_se = 5, n_ci = 1000))
  expect_true(is.na(p2$se_ratio))
  expect_false(is.na(p2$coverage))
})

test_that("the floor fires exactly at the boundary, not near it", {
  expect_true(is.na(.perf(.cell(n = 500, n_ci = 29))$coverage))
  expect_false(is.na(.perf(.cell(n = 500, n_ci = 30))$coverage))
  expect_false(is.na(.perf(.cell(n = 500, n_ci = 31))$coverage))

  expect_true(is.na(.perf(.cell(n = 500, n_se = 29))$se_ratio))
  expect_false(is.na(.perf(.cell(n = 500, n_se = 30))$se_ratio))
})

test_that("min_valid is honoured as given, including 0 and Inf", {
  x <- .cell(n = 500, n_se = 10, n_ci = 10)
  expect_false(is.na(.perf(x, min_valid = 0)$coverage))     # disabled
  expect_false(is.na(.perf(x, min_valid = 10)$coverage))    # at the boundary
  expect_true(is.na(.perf(x, min_valid = 11)$coverage))
  expect_true(is.na(.perf(x, min_valid = Inf)$coverage))    # withhold everything
  expect_false(is.na(.perf(x, min_valid = Inf)$bias))       # ... except the point family
})

test_that("the default floor is SIM_DEFAULTS$min_valid", {
  expect_equal(SIM_DEFAULTS$min_valid, 30L)
  expect_equal(eval(formals(performance)$min_valid), SIM_DEFAULTS$min_valid)
})

## --- regression: a complete cell is unchanged --------------------------------

test_that("a fully populated cell is numerically identical to before the fix", {
  # The suppression must be inert wherever it is not needed, or every shipped
  # number moves. Recomputed here from first principles rather than from the
  # function under test.
  x <- .cell(n = 400)
  p <- .perf(x)
  ok <- rep(TRUE, 400)

  expect_equal(p$bias, mean(x$est - x$target))
  expect_equal(p$emp_se, stats::sd(x$est))
  expect_equal(p$bias_mcse, stats::sd(x$est) / sqrt(400))
  expect_equal(p$emp_se_mcse, stats::sd(x$est) / sqrt(2 * 399))
  expect_equal(p$mod_se, sqrt(mean(x$se^2)))
  expect_equal(p$se_ratio, sqrt(mean(x$se^2)) / stats::sd(x$est))
  expect_equal(p$rmse, sqrt(mean((x$est - x$target)^2)))
  cov_i <- x$ci_lo <= x$target & x$ci_up >= x$target
  expect_equal(p$coverage, mean(cov_i))
  expect_equal(p$coverage_mcse, sqrt(mean(cov_i) * (1 - mean(cov_i)) / 400))
  expect_equal(p$ci_width, mean(x$ci_up - x$ci_lo))
  expect_equal(p$nonest_rate, 0)
})

test_that("coverage_mcse uses the CI count, not n_valid", {
  # It divided by sum(ok_ci) before and must still do so -- if it silently switched
  # to n_valid the MCSE of a partial cell would be understated.
  p <- .perf(.cell(n = 1000, n_ci = 100), min_valid = 0)
  expect_equal(p$coverage_mcse, sqrt(p$coverage * (1 - p$coverage) / 100))
  expect_false(isTRUE(all.equal(p$coverage_mcse,
                                sqrt(p$coverage * (1 - p$coverage) / 1000))))
})

## --- through summarise_raw() -------------------------------------------------

test_that("summarise_raw() carries the counts and forwards min_valid", {
  raw <- data.frame(
    cond = rep(c("a", "b"), each = 100),
    method = "m",
    est = stats::rnorm(200),
    se = 1,
    ci_lo = -2, ci_up = 2,
    theta_pop = 0
  )
  raw$se[1:95] <- NA_real_                      # cell "a" keeps 5 usable SEs
  raw$ci_lo[1:95] <- NA_real_

  out <- summarise_raw(raw, "cond", c(population = "theta_pop"))
  expect_true(all(c("n_valid_se", "n_valid_ci") %in% names(out)))
  a <- out[out$cond == "a", ]; b <- out[out$cond == "b", ]
  expect_equal(a$n_valid, 100); expect_equal(a$n_valid_se, 5)
  expect_equal(b$n_valid_se, 100)
  expect_true(is.na(a$se_ratio))                # withheld at the default floor
  expect_false(is.na(b$se_ratio))

  relaxed <- summarise_raw(raw, "cond", c(population = "theta_pop"), min_valid = 0)
  expect_false(is.na(relaxed$se_ratio[relaxed$cond == "a"]))
})

## --- stopifnot() fallback for a bare R (no testthat) ------------------------
if (!requireNamespace("testthat", quietly = TRUE)) {
  .p <- performance(rnorm(100), c(rep(NA, 95), rep(1, 5)),
                    rep(-2, 100), rep(2, 100), 0)
  stopifnot(
    all(c("n_valid_se", "n_valid_ci") %in% names(.p)),
    .p$n_valid == 100, .p$n_valid_se == 5,
    is.na(.p$se_ratio), !is.na(.p$bias)
  )
}
