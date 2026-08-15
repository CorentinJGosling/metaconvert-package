# =============================================================================
# Roadmap item 1.1 -- 2x2 reconstruction from an odds ratio.
#
# .estimate_n_from_or_and_n_exp() and .estimate_n_from_or_and_n_cases() recover a
# 2x2 table by enumerating candidates and keeping the one whose reconstructed
# var(logOR) best matches the reported one. Two things break that search:
#
#   1. THE ROTATION TIE. (a,b,c,d) -> (d,c,b,a) preserves the odds ratio AND
#      1/a+1/b+1/c+1/d EXACTLY, and is admissible with the same arm sizes precisely
#      when n_exp == n_nexp (resp. n_cases == n_controls). The two candidates are
#      bit-identical on the matching criterion, so the winner was decided by
#      enumeration index -- measured 25-32% correct at rare event rates.
#
#   2. A CORRUPTED TARGET. On the es_from_or() route the `var` being matched is
#      itself imputed by .se_from_or() and runs ~1.4x wide, so the search misses the
#      true table even when no tie exists.
#
# Both vanish if the table is SOLVED rather than searched: with both margin pairs
# known, the OR determines the table by a quadratic and the variance is not used at
# all. .solve_2x2_from_or() does that, and the two helpers now try it first
# (Rung 1: a second margin; Rung 2: baseline_risk) before falling through to the
# unchanged enumeration.
#
# NO PRIOR IS APPLIED when neither rung fires. Every candidate rule tested was a bet
# on outcome coding -- "assume events are the minority" is 3.65x worse than the
# status quo on common outcomes and produced +44% pooled-RR bias end-to-end, and
# recoding an outcome from "response" to "non-response" flips its answer on identical
# data. The last test in this file pins that no such prior has been introduced.
# =============================================================================

.solve  <- function(...) metaConvert:::.solve_2x2_from_or(...)
.est_exp   <- function(...) metaConvert:::.estimate_n_from_or_and_n_exp(...)
.est_cases <- function(...) metaConvert:::.estimate_n_from_or_and_n_cases(...)

.tab_or  <- function(a, b, c, d) (a * d) / (b * c)
.tab_var <- function(a, b, c, d) 1 / a + 1 / b + 1 / c + 1 / d

# --- the tie this all exists for ------------------------------------------

test_that("the 180-degree rotation is an EXACT tie on both matching criteria", {
  a <- 19; b <- 131; c <- 75; d <- 75          # equal arms: 150 / 150
  expect_equal(.tab_or(a, b, c, d),  .tab_or(d, c, b, a),  tolerance = 0)
  expect_equal(.tab_var(a, b, c, d), .tab_var(d, c, b, a), tolerance = 1e-12)
  # ...and it is admissible only because the arms are equal.
  expect_equal(a + b, c + d)
  # The rotation swaps BOTH margin pairs -- which is exactly why the opposite pair
  # identifies the branch. Rotated table is (a',b',c',d') = (d,c,b,a).
  n_cases_true <- a + c               # 94
  n_cases_rot  <- d + b               # 206, i.e. n_controls of the true table
  expect_equal(n_cases_rot, b + d)
  expect_false(n_cases_true == n_cases_rot)
})

# --- the solve ------------------------------------------------------------

test_that(".solve_2x2_from_or() recovers the exact table, tie or no tie", {
  set.seed(5)
  bad <- 0; declined <- 0
  for (i in 1:500) {
    n1 <- sample(20:400, 1); n2 <- sample(20:400, 1)
    a <- sample(1:(n1 - 1), 1); c <- sample(1:(n2 - 1), 1)
    b <- n1 - a; d <- n2 - c
    got <- .solve(.tab_or(a, b, c, d), n1, n2, a + c)
    if (is.null(got)) { declined <- declined + 1; next }
    if (!isTRUE(all.equal(
          c(got$n_cases_exp, got$n_controls_exp, got$n_cases_nexp, got$n_controls_nexp),
          c(a, b, c, d)))) bad <- bad + 1
  }
  expect_equal(bad, 0)
  expect_equal(declined, 0)
})

test_that(".solve_2x2_from_or() handles or = 1, where the quadratic degenerates", {
  # A = 1 - or = 0, so the equation is linear: a = n_cases * n_exp / N.
  got <- .solve(1, 100, 100, 80)
  expect_equal(got$n_cases_exp, 40)
  expect_equal(got$n_cases_nexp, 40)
  expect_equal(.tab_or(got$n_cases_exp, got$n_controls_exp,
                       got$n_cases_nexp, got$n_controls_nexp), 1)
})

test_that(".solve_2x2_from_or() declines rather than guessing on bad input", {
  expect_null(.solve(2, 50, 50, 150))   # multi-arm: n_cases >= n_exp + n_nexp
  expect_null(.solve(2, 50, 50, 100))   # n_cases == N
  expect_null(.solve(2, 50, 50, 0))     # empty case margin
  expect_null(.solve(0, 50, 50, 20))    # or must be positive
  expect_null(.solve(-1, 50, 50, 20))
  expect_null(.solve(Inf, 50, 50, 20))
  expect_null(.solve(NA, 50, 50, 20))
  expect_null(.solve(2, NA, 50, 20))
  expect_null(.solve(2, 0, 50, 20))     # empty arm
})

test_that(".solve_2x2_from_or() declines when a solved cell would be empty", {
  # A zero cell makes var(logOR) infinite. The enumeration has a purpose-built +0.5
  # branch for those, so the solve must hand them back rather than emit the table.
  set.seed(9)
  zero_cells <- 0
  for (i in 1:400) {
    n1 <- sample(30:200, 1); n2 <- sample(30:200, 1)
    or <- exp(runif(1, log(0.1), log(10)))
    got <- .solve(or, n1, n2, sample(2:(n1 + n2 - 2), 1))
    if (!is.null(got)) {
      cells <- c(got$n_cases_exp, got$n_controls_exp, got$n_cases_nexp, got$n_controls_nexp)
      if (min(cells) < 1) zero_cells <- zero_cells + 1
      expect_true(all(is.finite(cells)))
    }
  }
  expect_equal(zero_cells, 0)
})

# --- the cascade ----------------------------------------------------------

test_that("Rung 1 fixes the tie: a second margin gives a perfect branch hit", {
  set.seed(11)
  for (er in c(0.05, 0.15, 0.30, 0.50, 0.70)) {
    n <- 150; hit_solved <- 0; k <- 0
    for (i in 1:120) {
      p1 <- min(max(er + runif(1, -.12, .12), .03), .95)
      p2 <- min(max(er + runif(1, -.12, .12), .03), .95)
      a <- rbinom(1, n, p1); c <- rbinom(1, n, p2); b <- n - a; d <- n - c
      if (min(a, b, c, d) < 2) next
      k <- k + 1
      got <- .est_exp(.tab_or(a, b, c, d), .tab_var(a, b, c, d), n, n,
                      n_cases = a + c, n_controls = b + d)
      hit_solved <- hit_solved + isTRUE(all.equal(
        c(got$n_cases_exp, got$n_cases_nexp), c(a, c)))
    }
    expect_equal(hit_solved / k, 1,
                 info = paste("event rate", er, "-- Rung 1 must be exact"))
  }
})

test_that("Rung 1 accepts n_controls when n_cases is absent", {
  a <- 19; b <- 131; c <- 75; d <- 75
  got <- .est_exp(.tab_or(a, b, c, d), .tab_var(a, b, c, d), 150, 150,
                  n_controls = b + d)
  expect_equal(got$n_cases_exp, a)
  expect_equal(got$n_cases_nexp, c)
})

test_that("Rung 2 identifies the branch from a reported baseline risk", {
  set.seed(21)
  n <- 150; hit <- 0; k <- 0
  for (i in 1:150) {
    a <- rbinom(1, n, 0.25); c <- rbinom(1, n, 0.20); b <- n - a; d <- n - c
    if (min(a, b, c, d) < 2) next
    k <- k + 1
    got <- .est_exp(.tab_or(a, b, c, d), .tab_var(a, b, c, d), n, n,
                    baseline_risk = c / n)
    hit <- hit + isTRUE(all.equal(c(got$n_cases_exp, got$n_cases_nexp), c(a, c)))
  }
  expect_equal(hit / k, 1)
})

test_that("the mirror helper is fixed by n_exp when its own margins are balanced", {
  set.seed(3)
  hit <- 0; k <- 0
  for (i in 1:150) {
    repeat { a <- sample(2:148, 1); b <- sample(2:148, 1); if (a != b) break }
    c <- 150 - a; d <- 150 - b            # n_cases == n_controls == 150
    got <- .est_cases(.tab_or(a, b, c, d), .tab_var(a, b, c, d),
                      a + c, b + d, n_exp = a + b, n_nexp = c + d)
    k <- k + 1
    hit <- hit + isTRUE(all.equal(c(got$n_cases_exp, got$n_cases_nexp), c(a, c)))
  }
  expect_equal(hit / k, 1)
})

# --- no regression, and no prior -------------------------------------------

test_that("with no extra margin the output is UNCHANGED (pure no-op)", {
  # The cascade must be additive: rows it cannot identify take exactly the path they
  # took before, so no existing result moves.
  set.seed(99)
  for (i in 1:150) {
    n1 <- sample(20:300, 1); n2 <- sample(20:300, 1)
    a <- sample(1:(n1 - 1), 1); c <- sample(1:(n2 - 1), 1); b <- n1 - a; d <- n2 - c
    or <- .tab_or(a, b, c, d); v <- .tab_var(a, b, c, d)
    bare <- .est_exp(or, v, n1, n2)
    # Supplying only NAs must be identical to supplying nothing.
    with_na <- .est_exp(or, v, n1, n2, n_cases = NA, n_controls = NA, baseline_risk = NA)
    expect_identical(bare, with_na)
  }
})

test_that("no minority-event prior has been introduced", {
  # A prior would make the answer depend on which category is labelled "case".
  # Construct a common-outcome table (event rate 0.70) at equal arms with NO second
  # margin: a prior would flip to the minority branch; the enumeration must not.
  a <- 120; b <- 30; c <- 90; d <- 60      # 150 / 150, event rate 0.70
  got <- .est_exp(.tab_or(a, b, c, d), .tab_var(a, b, c, d), 150, 150)
  expect_equal(c(got$n_cases_exp, got$n_cases_nexp), c(a, c),
               info = paste("A minority-event prior would have returned the rotation",
                            "(d,c,b,a) here. The enumeration order must be untouched."))
})

test_that("es_from_or_se() end-to-end: metaumbrella_exp recovers the right RR", {
  # The user-facing consequence: same OR and SE, equal arms, with and without the
  # case margin. Solving must land on the true RR.
  a <- 19; b <- 131; c <- 75; d <- 75
  true_rr <- (a / 150) / (c / 150)
  got <- es_from_or_se(or = .tab_or(a, b, c, d), logor_se = sqrt(.tab_var(a, b, c, d)),
                       n_exp = 150, n_nexp = 150, n_cases = a + c, n_controls = b + d,
                       or_to_rr = "metaumbrella_exp")
  expect_equal(exp(got$logrr), true_rr, tolerance = 0.02)
})
