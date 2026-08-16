# =============================================================================
# Roadmap item 1.4 -- rr_to_or = "grant" must not return a finite effect size with
# a missing variance.
#
# Grant's conversion is
#     OR = RR (1 - BR) / (1 - RR * BR)
# which is defined only while RR * BR < 1. metaConvert applies it to THREE values --
# the point estimate and both confidence limits -- and derives the standard error
# from the transformed interval width. The upper limit is the largest of the three,
# so it leaves the domain first.
#
# The result was an asymmetric, silent failure: a row could return a plausible log
# odds ratio beside a NaN standard error and a half-open interval, with every log()
# wrapped in suppressWarnings(). A finite estimate with no usable variance is worse
# than no estimate, because it looks poolable and is not -- and in an
# inverse-variance meta-analysis a missing SE quietly drops the study instead of
# reporting that it could not be converted.
#
# The row is now returned as NA throughout, with a warning naming the condition.
#
# The MIRROR direction is deliberately untouched: .or_to_rr()'s grant branch computes
# or / (1 - BR + BR*or), whose denominator is positive for every or > 0 and
# BR in (0, 1), so it has no domain boundary. The last test pins that asymmetry so a
# future "fix for symmetry" does not add a guard that can never fire.
# =============================================================================

.grant_rr_to_or <- function(rr, lo, up, br) {
  es_from_rr_ci(rr = rr, rr_ci_lo = lo, rr_ci_up = up,
                baseline_risk = br, rr_to_or = "grant")
}

.quartet <- function(x) c(x$logor, x$logor_se, x$logor_ci_lo, x$logor_ci_up)

# --- the defect ------------------------------------------------------------

test_that("an out-of-domain upper limit returns NA throughout, not a mixed row", {
  # rr_ci_up * br = 2.5 * 0.5 = 1.25 >= 1, while rr * br = 0.75 < 1: the point
  # estimate is in the domain and the upper limit is not. This is the exact shape
  # that used to yield logor = 1.0986 with logor_se = NaN.
  res <- suppressWarnings(.grant_rr_to_or(1.5, 0.9, 2.5, 0.5))
  expect_true(all(is.na(.quartet(res))))
})

test_that("it warns, and the warning names the offending value and the condition", {
  expect_warning(.grant_rr_to_or(1.5, 0.9, 2.5, 0.5), "rr_to_or = 'grant'")
  w <- tryCatch(.grant_rr_to_or(1.5, 0.9, 2.5, 0.5),
                warning = function(e) conditionMessage(e))
  expect_match(w, "rr_ci_up = 2.5")
  expect_match(w, "baseline_risk = 0.5")
  expect_match(w, "RR \\* baseline_risk < 1")
  # It must point at a way out.
  expect_match(w, "metaumbrella")
})

test_that("an out-of-domain POINT ESTIMATE also returns NA throughout", {
  # rr * br = 3 * 0.5 = 1.5: everything is out of the domain.
  res <- suppressWarnings(.grant_rr_to_or(3, 2.2, 4.1, 0.5))
  expect_true(all(is.na(.quartet(res))))
})

test_that("no input can produce a partly-finite quartet", {
  # The property that matters, swept rather than spot-checked: whatever the inputs,
  # the four outputs are either all finite or all NA.
  set.seed(9)
  mixed <- 0L; n <- 0L
  for (i in 1:400) {
    rr <- exp(rnorm(1, 0, 0.5)); se <- runif(1, .1, .5); br <- runif(1, .05, .9)
    res <- suppressWarnings(
      .grant_rr_to_or(rr, rr * exp(-1.96 * se), rr * exp(1.96 * se), br))
    v <- .quartet(res); n <- n + 1L
    if (any(is.finite(v)) && !all(is.finite(v))) mixed <- mixed + 1L
  }
  expect_gt(n, 300)
  expect_equal(mixed, 0L)
})

# --- no regression in the domain -------------------------------------------

test_that("an in-domain conversion is unchanged and fully finite", {
  # NB es_from_rr_ci() derives logrr_se from the CI width and delegates to
  # es_from_rr_se(), which rebuilds a LOG-SYMMETRIC interval. An asymmetric input
  # interval therefore does not round-trip, so the fixture below is built symmetric
  # on the log scale and the expectation is exact.
  g <- function(x, br) log(x * (1 - br) / (1 - x * br))
  rr <- 1.2; se <- 0.1; br <- 0.2; z <- qnorm(.975)
  lo <- rr * exp(-z * se); up <- rr * exp(z * se)   # up * br = 0.29, in domain

  res <- .grant_rr_to_or(rr, lo, up, br)
  expect_true(all(is.finite(.quartet(res))))
  expect_equal(res$logor, g(rr, br), tolerance = 1e-10)
  expect_equal(res$logor_ci_lo, g(lo, br), tolerance = 1e-10)
  expect_equal(res$logor_ci_up, g(up, br), tolerance = 1e-10)
  expect_equal(res$logor_se, (g(up, br) - g(lo, br)) / (2 * z), tolerance = 1e-10)
})

test_that("an in-domain conversion emits no warning", {
  expect_no_warning(.grant_rr_to_or(1.2, 1.0, 1.45, 0.2))
  expect_no_warning(.grant_rr_to_or(1.05, 0.90, 1.22, 0.05))
})

test_that("the other rr_to_or methods are untouched by this guard", {
  # They have no domain restriction, so the same inputs must still convert.
  for (m in c("metaumbrella", "transpose")) {
    res <- suppressWarnings(es_from_rr_se(
      rr = 1.5, logrr_se = 0.25, baseline_risk = 0.5,
      n_cases = 60, n_controls = 140, rr_to_or = m))
    expect_true(is.finite(res$logor),
                info = paste0("rr_to_or = '", m, "' must still produce an estimate"))
  }
})

# --- the asymmetry is deliberate -------------------------------------------

test_that("the or_to_rr mirror has no domain boundary and is left alone", {
  # or / (1 - BR + BR*or): the denominator is 1 - BR + BR*or, and with BR in (0, 1)
  # and or > 0 both terms are positive, so it can never reach zero.
  for (or in c(0.01, 0.1, 1, 5, 50, 500)) {
    for (br in c(0.01, 0.05, 0.5, 0.95, 0.99)) {
      res <- suppressWarnings(es_from_or_ci(
        or = or, or_ci_lo = or * 0.6, or_ci_up = or * 1.8,
        baseline_risk = br, or_to_rr = "grant"))
      expect_true(is.finite(res$logrr),
                  info = paste0("or = ", or, ", baseline_risk = ", br))
      expect_true(is.finite(res$logrr_se),
                  info = paste0("or = ", or, ", baseline_risk = ", br))
    }
  }
})
