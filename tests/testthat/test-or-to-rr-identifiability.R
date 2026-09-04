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

test_that(".solve_2x2_from_or() declines when a solved cell would be (near-)empty", {
  # A cell below 0.5 is smaller than any real table carries, raw (>= 1) or
  # continuity-corrected (>= 0.5), and its 1/cell term would dominate the Woolf
  # variance. The enumeration has a purpose-built +0.5 branch for those, so the solve
  # must hand them back rather than emit the table. Cells are otherwise the EXACT
  # root and may be fractional: the solved table's odds ratio must be the reported one.
  set.seed(9)
  small_cells <- 0; checked <- 0
  for (i in 1:400) {
    n1 <- sample(30:200, 1); n2 <- sample(30:200, 1)
    or <- exp(runif(1, log(0.1), log(10)))
    got <- .solve(or, n1, n2, sample(2:(n1 + n2 - 2), 1))
    if (!is.null(got)) {
      cells <- c(got$n_cases_exp, got$n_controls_exp, got$n_cases_nexp, got$n_controls_nexp)
      if (min(cells) < 0.5) small_cells <- small_cells + 1
      expect_true(all(is.finite(cells)))
      expect_equal(.tab_or(cells[1], cells[2], cells[3], cells[4]), or, tolerance = 1e-10)
      checked <- checked + 1
    }
  }
  expect_equal(small_cells, 0)
  expect_gt(checked, 300)
})

test_that("the solved table reproduces a REPORTED odds ratio exactly, not a rounded neighbour", {
  # OR = 2.00 on margins 40/60, 30/70 has no integer table: rounding the root to
  # 16/24/14/46 gave a table with OR = 2.19 and log RR 0.539; the exact root
  # a = 15.538 keeps OR = 2 and gives log RR 0.477.
  got <- .solve(2, 40, 60, 30)
  expect_equal(.tab_or(got$n_cases_exp, got$n_controls_exp, got$n_cases_nexp, got$n_controls_nexp), 2,
               tolerance = 1e-12)
  expect_equal(got$n_cases_exp, (170 - sqrt(170^2 - 4 * 2400)) / 2, tolerance = 1e-12)
  lrr <- log((got$n_cases_exp / 40) / (got$n_cases_nexp / 60))
  expect_equal(lrr, 0.4771998608, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(lrr, 0.5389965007)))   # the rounded-table answer
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

test_that("each reconstruction helper is defined exactly once across R/", {
  # R/estimate_n_from_es.R used to hold a SECOND definition of both OR helpers. With no
  # Collate field in DESCRIPTION, files are sourced alphabetically, so
  # internal_multiple_formulas.R ("i") overwrote estimate_n_from_es.R ("e") and the
  # copies there were dead -- but editable, and by the end they had drifted from the
  # live versions. They are deleted; this stops a duplicate reappearing unnoticed.
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = FALSE)
  r_files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(length(r_files) == 0, "not a source checkout")

  for (fn in c(".estimate_n_from_or_and_n_cases", ".estimate_n_from_or_and_n_exp",
               ".estimate_n_from_irr", ".estimate_n_from_rr", ".solve_2x2_from_or")) {
    n <- sum(vapply(r_files, function(f) {
      sum(grepl(paste0("^", gsub("\\.", "\\\\.", fn), "\\s*<-\\s*function"),
                readLines(f, warn = FALSE)))
    }, integer(1)))
    expect_equal(n, 1L, info = paste0(fn, " is defined ", n, " time(s) across R/"))
  }
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

# =============================================================================
# Roadmap item 1.8 -- the gate on the non-identified cell.
#
# Item 1.1 (above) SOLVES the reconstruction whenever a margin from the other pair
# is available. It does nothing for the documented minimal input -- or + logor_se +
# n_exp + n_nexp -- which contains no such margin. At EQUAL margins that input is
# genuinely non-identified: the 180-degree rotation reproduces the odds ratio and its
# standard error exactly while implying a different risk ratio, so the enumeration is
# tied and its winner is an artifact of loop order.
#
# Measured on the _exp parameterisation, 782 usable draws at n_exp == n_nexp == 50
# with no second margin: the true table comes back 37.0% of the time, its rotation
# 66.6% (4.4% of tables are their own rotation), genuinely wrong 62.3%. The cliff is
# sharp, not gradual -- one participant of imbalance makes the rotation inadmissible:
#
#   |n_exp - n_nexp|    0       1       2       5      10
#   hit rate         0.384   0.995   0.987   0.995   0.982     (exact SE)
#                    0.379   0.955   0.967   0.950   0.951     (OR rounded to 2dp)
#   _cases mirror    0.745   0.995   0.995   0.993   0.993
#
# So the gate tests EXACT equality and nothing looser, and fires only when the exact
# solve did not. Rows it fires on return NA for the RR quartet with a warning naming
# the column that would break the tie; every other row is untouched.
#
# Only the RR outputs are withheld: OR, D, G, R, Z, RD and NNT never touch the
# reconstruction.
# =============================================================================

.tied <- function(...) metaConvert:::.rotation_tied(...)

# Count warnings without stopping at the first (mapply raises one per gated row).
.catch_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr,
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  list(value = val, warnings = w)
}

# A concrete non-identified case, used throughout: true table 1/49/25/25.
# true OR = 0.0204, true RR = 0.040. Ungated, the enumeration returns the rotation
# 25/49/25/1, whose RR is 25 -- a 625-fold error.
.case <- list(a = 1, b = 49, c = 25, d = 25, n = 50)
.case_or <- with(.case, .tab_or(a, b, c, d))
.case_se <- with(.case, sqrt(.tab_var(a, b, c, d)))
.case_rr <- with(.case, (a / n) / (c / n))

# --- the "solved" marker ---------------------------------------------------

test_that("the helpers report whether they SOLVED or searched", {
  or <- .case_or; v <- .case_se^2

  # Rung 1, via n_cases / via n_controls
  expect_true(attr(.est_exp(or, v, 50, 50, n_cases = 26), "solved"))
  expect_true(attr(.est_exp(or, v, 50, 50, n_controls = 74), "solved"))
  # Rung 2, via baseline_risk
  expect_true(attr(.est_exp(or, v, 50, 50, baseline_risk = 25 / 50), "solved"))
  # Fall-through to the enumeration
  expect_false(attr(.est_exp(or, v, 50, 50), "solved"))
  expect_false(attr(.est_exp(or, v, 50, 50, n_cases = NA, baseline_risk = NA), "solved"))

  # The mirror helper, whose Rung 1 takes n_exp
  a <- 19; b <- 131; c <- 75; d <- 75
  or2 <- .tab_or(a, b, c, d); v2 <- .tab_var(a, b, c, d)
  expect_true(attr(.est_cases(or2, v2, a + c, b + d, n_exp = a + b, n_nexp = c + d), "solved"))
  expect_false(attr(.est_cases(or2, v2, a + c, b + d), "solved"))
})

test_that("the marker survives the enumeration branch's column assignments", {
  # The enumeration writes each cell into `res` with $<- AFTER the attribute is set.
  # If that dropped the attribute the marker would read NULL, and !isTRUE(NULL) is
  # TRUE, so the gate would still fire -- but for the wrong reason, and a solved row
  # would be indistinguishable from a searched one. Pin FALSE, not missing.
  got <- .est_exp(.case_or, .case_se^2, 50, 50)
  expect_identical(attr(got, "solved"), FALSE)
  expect_false(is.null(attr(got, "solved")))
  expect_false(anyNA(unlist(got)))          # the enumeration did return a table
})

test_that("a second margin that does NOT describe this table falls through as unsolved", {
  # .solve_2x2_from_or() declines a case margin >= N (a multi-arm trial reporting a
  # margin over arms not in n_exp + n_nexp). Keying the gate on what the helper
  # actually did, rather than on "was a second margin supplied", catches this.
  got <- .est_exp(.case_or, .case_se^2, 50, 50, n_cases = 500)
  expect_false(attr(got, "solved"))
})

# --- .rotation_tied() ------------------------------------------------------

test_that(".rotation_tied() is exact equality and nothing looser", {
  expect_true(.tied(50, 50))
  expect_true(.tied(50L, 50))
  expect_false(.tied(50, 51))
  expect_false(.tied(50, 49))
  expect_false(.tied(500, 501))       # a 0.2% imbalance is still not a tie
  expect_false(.tied(NA, 50))
  expect_false(.tied(50, NA))
  expect_false(.tied(NA, NA))
  expect_false(.tied(NaN, NaN))
  expect_false(.tied(Inf, Inf))
  expect_false(.tied(numeric(0), 50))
  expect_false(.tied(c(50, 50), c(50, 50)))   # never silently vectorises
})

# --- the gate, metaumbrella_exp -------------------------------------------

test_that("equal arms + minimal input: the RR quartet is NA and a warning names the fix", {
  res <- .catch_warnings(
    es_from_or_se(or = .case_or, logor_se = .case_se, n_exp = 50, n_nexp = 50,
                  or_to_rr = "metaumbrella_exp"))

  expect_true(is.na(res$value$logrr))
  expect_true(is.na(res$value$logrr_se))
  expect_true(is.na(res$value$logrr_ci_lo))
  expect_true(is.na(res$value$logrr_ci_up))

  expect_length(res$warnings, 1)
  expect_match(res$warnings[1], "metaumbrella_exp", fixed = TRUE)
  expect_match(res$warnings[1], "not identified when n_exp == n_nexp", fixed = TRUE)
  expect_match(res$warnings[1], "here both are 50", fixed = TRUE)
  expect_match(res$warnings[1], "n_cases", fixed = TRUE)
  expect_match(res$warnings[1], "baseline_risk", fixed = TRUE)
  # No hit-rate percentage is quoted: the two parameterisations measured differently.
  expect_false(grepl("%", res$warnings[1], fixed = TRUE))
})

test_that("the withheld value is the one that used to be wrong", {
  # Pin what the gate prevents, so a future "restore the old behaviour" change has to
  # confront the number. The enumeration returns the rotated table -- cells
  # (n_cases_exp, n_cases_nexp, n_controls_exp, n_controls_nexp) = (25, 49, 25, 1),
  # i.e. the true 1/49/25/25 turned 180 degrees.
  got <- .est_exp(.case_or, .case_se^2, 50, 50)
  expect_equal(c(got$n_cases_exp, got$n_cases_nexp,
                 got$n_controls_exp, got$n_controls_nexp), c(25, 49, 25, 1))
  wrong_rr <- (got$n_cases_exp / 50) / (got$n_cases_nexp / 50)
  expect_equal(wrong_rr, 25 / 49)             # 0.510, against a true RR of 0.04
  expect_equal(wrong_rr / .case_rr, 12.755, tolerance = 1e-3)   # 12.8-fold, unflagged
  # Both tables reproduce the reported inputs exactly -- which is the whole problem.
  expect_equal(.tab_or(got$n_cases_exp, got$n_controls_exp,
                       got$n_cases_nexp, got$n_controls_nexp), .case_or)
  expect_equal(.tab_var(got$n_cases_exp, got$n_controls_exp,
                        got$n_cases_nexp, got$n_controls_nexp), .case_se^2,
               tolerance = 1e-12)
})

test_that("only the RR columns are withheld", {
  res <- suppressWarnings(
    es_from_or_se(or = .case_or, logor_se = .case_se, n_exp = 50, n_nexp = 50,
                  or_to_rr = "metaumbrella_exp"))
  # These come from the log OR directly and never touch the reconstruction.
  for (col in c("logor", "logor_se", "logor_ci_lo", "logor_ci_up",
                "d", "d_se", "g", "g_se", "r", "r_se", "z", "z_se")) {
    expect_false(is.na(res[[col]]), info = paste(col, "must survive the RR gate"))
  }
  expect_equal(res$logor, log(.case_or))
  expect_equal(res$logor_se, .case_se)
})

test_that("RD and NNT survive the gate (measured on the _cases mirror)", {
  # For metaumbrella_exp a supplied baseline_risk identifies the table, so the gate
  # cannot fire alongside an RD. The _cases parameterisation does not consume
  # baseline_risk (Rung 2 is deliberately not wired there -- it is 0.9965, not exact),
  # so this is the configuration where a gated RR coexists with a computed RD/NNT.
  res <- suppressWarnings(
    es_from_or_se(or = 2, logor_se = 0.3, n_cases = 60, n_controls = 60,
                  baseline_risk = 0.2, or_to_rr = "metaumbrella_cases"))
  expect_true(is.na(res$logrr))
  expect_false(is.na(res$rd))
  expect_false(is.na(res$nnt))
  expect_false(is.na(res$rd_se))
  # and the RD is the OR-based one, untouched by the reconstruction
  treatment_risk <- (2 * 0.2) / (1 - 0.2 + 2 * 0.2)
  expect_equal(res$rd, 0.2 - treatment_risk)
})

test_that("each recovery path returns the TRUE risk ratio, silently", {
  for (extra in list(list(n_cases = 26), list(n_controls = 74),
                     list(baseline_risk = 25 / 50),
                     list(n_cases = 26, n_controls = 74))) {
    res <- .catch_warnings(do.call(es_from_or_se, c(
      list(or = .case_or, logor_se = .case_se, n_exp = 50, n_nexp = 50,
           or_to_rr = "metaumbrella_exp"), extra)))
    expect_equal(exp(res$value$logrr), .case_rr, tolerance = 1e-8,
                 info = paste("recovery via", paste(names(extra), collapse = "+")))
    expect_length(res$warnings, 0)
  }
})

test_that("unequal arms are untouched and silent", {
  # One participant of imbalance is enough: the rotation is inadmissible, so the
  # search is not tied and today's behaviour stands.
  res <- .catch_warnings(
    es_from_or_se(or = .case_or, logor_se = .case_se, n_exp = 50, n_nexp = 51,
                  or_to_rr = "metaumbrella_exp"))
  expect_false(is.na(res$value$logrr))
  expect_length(res$warnings, 0)
})

test_that("solved rows at EQUAL arms are not gated", {
  # The gate keys on "was it solved", not on "are the margins equal". A row with equal
  # arms AND a case margin is fully identified and must produce a value.
  set.seed(1808)
  k <- 0
  for (i in 1:40) {
    n <- sample(30:200, 1)
    a <- rbinom(1, n, runif(1, .05, .60)); c <- rbinom(1, n, runif(1, .05, .60))
    b <- n - a; d <- n - c
    if (min(a, b, c, d) < 2) next
    k <- k + 1
    res <- .catch_warnings(
      es_from_or_se(or = .tab_or(a, b, c, d), logor_se = sqrt(.tab_var(a, b, c, d)),
                    n_exp = n, n_nexp = n, n_cases = a + c, n_controls = b + d,
                    or_to_rr = "metaumbrella_exp"))
    expect_false(is.na(res$value$logrr))
    expect_length(res$warnings, 0)
    expect_equal(exp(res$value$logrr), (a / n) / (c / n), tolerance = 1e-8)
  }
  expect_gt(k, 20)
})

# --- the gate, metaumbrella_cases (the mirror) -----------------------------

test_that("the mirror gates at n_cases == n_controls and names n_exp", {
  res <- .catch_warnings(
    es_from_or_se(or = 2, logor_se = 0.3, n_cases = 60, n_controls = 60,
                  or_to_rr = "metaumbrella_cases"))
  expect_true(all(is.na(c(res$value$logrr, res$value$logrr_se,
                          res$value$logrr_ci_lo, res$value$logrr_ci_up))))
  expect_length(res$warnings, 1)
  expect_match(res$warnings[1], "metaumbrella_cases", fixed = TRUE)
  expect_match(res$warnings[1], "not identified when n_cases == n_controls", fixed = TRUE)
  expect_match(res$warnings[1], "n_exp", fixed = TRUE)
  # baseline_risk is NOT offered as a fix here: Rung 2 is wired for _exp only.
  expect_false(grepl("baseline_risk", res$warnings[1], fixed = TRUE))
})

test_that("the mirror recovers from n_exp and leaves unbalanced case margins alone", {
  a <- 19; b <- 131; c <- 75; d <- 75          # n_cases 94, n_controls 206
  or <- .tab_or(a, b, c, d); se <- sqrt(.tab_var(a, b, c, d))

  # unbalanced case margins: not tied, no gate
  res1 <- .catch_warnings(
    es_from_or_se(or = or, logor_se = se, n_cases = a + c, n_controls = b + d,
                  or_to_rr = "metaumbrella_cases"))
  expect_false(is.na(res1$value$logrr))
  expect_length(res1$warnings, 0)

  # balanced case margins + n_exp: solved, no gate, exact
  a2 <- 40; b2 <- 60; c2 <- 60; d2 <- 40       # n_cases = n_controls = 100
  res2 <- .catch_warnings(
    es_from_or_se(or = .tab_or(a2, b2, c2, d2), logor_se = sqrt(.tab_var(a2, b2, c2, d2)),
                  n_cases = a2 + c2, n_controls = b2 + d2,
                  n_exp = a2 + b2, n_nexp = c2 + d2, or_to_rr = "metaumbrella_cases"))
  expect_false(is.na(res2$value$logrr))
  expect_length(res2$warnings, 0)
  expect_equal(exp(res2$value$logrr), (a2 / (a2 + b2)) / (c2 / (c2 + d2)), tolerance = 1e-8)
})

# --- vectorised behaviour --------------------------------------------------

test_that("in a mixed vector the NA lands on the gated row only", {
  res <- .catch_warnings(
    es_from_or_se(or = rep(.case_or, 4), logor_se = rep(.case_se, 4),
                  n_exp = c(50, 50, 50, 50), n_nexp = c(50, 51, 50, 50),
                  n_cases = c(NA, NA, 26, NA), n_controls = c(NA, NA, 74, NA),
                  or_to_rr = "metaumbrella_exp"))
  v <- res$value
  expect_true(is.na(v$logrr[1]))            # gated
  expect_false(is.na(v$logrr[2]))           # unequal arms
  expect_false(is.na(v$logrr[3]))           # solved
  expect_true(is.na(v$logrr[4]))            # gated
  expect_equal(exp(v$logrr[3]), .case_rr, tolerance = 1e-8)
  expect_length(res$warnings, 2)
  # the whole quartet moves together, row by row
  expect_true(all(is.na(c(v$logrr_se[1], v$logrr_ci_lo[1], v$logrr_ci_up[1]))))
  expect_true(all(is.na(c(v$logrr_se[4], v$logrr_ci_lo[4], v$logrr_ci_up[4]))))
  expect_false(anyNA(c(v$logrr_se[2], v$logrr_ci_lo[2], v$logrr_ci_up[2])))
})

test_that("reverse_or on a gated row cannot resurrect a value", {
  res <- suppressWarnings(
    es_from_or_se(or = .case_or, logor_se = .case_se, n_exp = 50, n_nexp = 50,
                  reverse_or = TRUE, or_to_rr = "metaumbrella_exp"))
  expect_true(all(is.na(c(res$logrr, res$logrr_se, res$logrr_ci_lo, res$logrr_ci_up))))
  expect_equal(res$logor, -log(.case_or))   # the OR still reverses normally
})

# --- other entry points ----------------------------------------------------

test_that("the gate reaches es_from_or_ci() and es_from_or_pval()", {
  ci_lo <- exp(log(.case_or) - qnorm(.975) * .case_se)
  ci_up <- exp(log(.case_or) + qnorm(.975) * .case_se)

  r1 <- .catch_warnings(es_from_or_ci(or = .case_or, or_ci_lo = ci_lo, or_ci_up = ci_up,
                                      n_exp = 50, n_nexp = 50,
                                      or_to_rr = "metaumbrella_exp"))
  expect_true(is.na(r1$value$logrr))
  expect_gte(length(r1$warnings), 1)

  pv <- 2 * pnorm(abs(log(.case_or)) / .case_se, lower.tail = FALSE)
  r2 <- .catch_warnings(es_from_or_pval(or = .case_or, or_pval = pv,
                                        n_exp = 50, n_nexp = 50,
                                        or_to_rr = "metaumbrella_exp"))
  expect_true(is.na(r2$value$logrr))
  expect_gte(length(r2$warnings), 1)
})

test_that("es_from_or() is unaffected: its own inputs always identify the table", {
  # es_from_or() requires n_cases + n_controls, so with n_exp/n_nexp present the
  # solve fires (Rung 1) and the gate cannot. This also pins that the item-1.1 solve
  # still bypasses the imputed .se_from_or() variance.
  res <- .catch_warnings(es_from_or(or = .case_or, n_cases = 26, n_controls = 74,
                                    n_exp = 50, n_nexp = 50,
                                    or_to_rr = "metaumbrella_exp"))
  expect_false(is.na(res$value$logrr))
  expect_equal(exp(res$value$logrr), .case_rr, tolerance = 1e-8)
  expect_length(res$warnings, 0)
})

test_that("the other or_to_rr methods are untouched at equal margins", {
  for (m in c("transpose", "grant", "dipietrantonj")) {
    res <- .catch_warnings(
      es_from_or_se(or = .case_or, logor_se = .case_se, n_exp = 50, n_nexp = 50,
                    baseline_risk = 0.5, or_to_rr = m))
    expect_false(any(grepl("not identified when", res$warnings, fixed = TRUE)),
                 info = paste(m, "must not be gated"))
  }
})

# --- convert_df() ----------------------------------------------------------

test_that("convert_df() propagates the gate to the summary output", {
  df <- data.frame(or = rep(.case_or, 2), logor_se = rep(.case_se, 2),
                   n_exp = c(50, 50), n_nexp = c(50, 51))
  res <- suppressWarnings(summary(suppressWarnings(
    convert_df(df, measure = "rr", or_to_rr = "metaumbrella_exp"))))
  # The estimate is es_crude (on the log scale for measure = "rr"); the `rr` column in
  # the summary is the echoed INPUT column, all-NA here, and checking it would pass
  # for the wrong reason.
  expect_true("es_crude" %in% names(res))
  expect_true(is.na(res$es_crude[1]))       # equal arms -> withheld
  expect_false(is.na(res$es_crude[2]))      # unequal arms -> unchanged
  expect_true(is.na(res$info_used_crude[1]))
  expect_equal(res$info_used_crude[2], "or_se")
})

test_that("convert_df() at the default or_to_rr gates on balanced case margins", {
  # The package default is metaumbrella_cases, so the tied configuration there is
  # n_cases == n_controls.
  df <- data.frame(or = c(2, 2), logor_se = c(0.3, 0.3),
                   n_cases = c(60, 60), n_controls = c(60, 90))
  res <- suppressWarnings(summary(suppressWarnings(convert_df(df, measure = "rr"))))
  expect_true(is.na(res$es_crude[1]))
  expect_false(is.na(res$es_crude[2]))
})

# --- properties over a sweep ----------------------------------------------

test_that("the gate is deterministic: every tied row is withheld, none other is", {
  set.seed(1809)
  tied_na <- tied_n <- untied_na <- untied_n <- 0L
  tied_warn <- untied_warn <- 0L
  for (i in 1:120) {
    n <- sample(30:250, 1)
    a <- rbinom(1, n, runif(1, .03, .60)); c <- rbinom(1, n, runif(1, .03, .60))
    b <- n - a; d <- n - c
    if (min(a, b, c, d) < 2) next
    or <- .tab_or(a, b, c, d); se <- sqrt(.tab_var(a, b, c, d))

    r_tied <- .catch_warnings(
      es_from_or_se(or = or, logor_se = se, n_exp = n, n_nexp = n,
                    or_to_rr = "metaumbrella_exp"))
    tied_n <- tied_n + 1L
    tied_na <- tied_na + is.na(r_tied$value$logrr)
    tied_warn <- tied_warn + (length(r_tied$warnings) == 1L)

    # same table, one extra participant in the control arm -> no longer tied
    r_untied <- .catch_warnings(
      es_from_or_se(or = or, logor_se = se, n_exp = n, n_nexp = n + 1,
                    or_to_rr = "metaumbrella_exp"))
    untied_n <- untied_n + 1L
    untied_na <- untied_na + is.na(r_untied$value$logrr)
    untied_warn <- untied_warn + length(r_untied$warnings)
  }
  expect_gt(tied_n, 50)
  expect_equal(tied_na, tied_n)          # every tied row withheld
  expect_equal(tied_warn, tied_n)        # exactly one warning each
  expect_equal(untied_na, 0)             # no untied row withheld
  expect_equal(untied_warn, 0)           # and none warned
})

test_that("no row ever returns a finite RR beside a non-finite SE", {
  set.seed(1810)
  for (i in 1:80) {
    n1 <- sample(20:200, 1); n2 <- sample(20:200, 1)
    a <- sample(1:(n1 - 1), 1); c <- sample(1:(n2 - 1), 1); b <- n1 - a; d <- n2 - c
    for (m in c("metaumbrella_exp", "metaumbrella_cases")) {
      res <- suppressWarnings(
        es_from_or_se(or = .tab_or(a, b, c, d), logor_se = sqrt(.tab_var(a, b, c, d)),
                      n_exp = n1, n_nexp = n2, n_cases = NA, n_controls = NA,
                      or_to_rr = m))
      quartet <- c(res$logrr, res$logrr_se, res$logrr_ci_lo, res$logrr_ci_up)
      expect_true(all(is.finite(quartet)) || all(is.na(quartet)),
                  info = paste(m, "quartet must be all-finite or all-NA"))
    }
  }
})

# =============================================================================
# Roadmap item 1.1 follow-up -- the exact solve must not destroy a CORRECTED table.
#
# Surfaced by the 3.5 regeneration, not by any test: study 04's metaumbrella_cases
# mean |bias| went 0.006 -> 0.031, concentrated entirely at br = 0.01. Cause: the
# solve ended with a <- round(a), but a table that has received a +0.5 continuity
# correction -- which es_from_2x2() itself emits for a zero cell, and which study
# 04 supplies -- has HALF-INTEGER cells. Nothing in the inputs reveals it, because
# adding 0.5 to all four cells adds exactly 1 to every margin.
#
# Worked case: the corrected table 0.5/50.5/2.5/48.5 was reconstructed as
# 1/50/2/49, giving log RR -0.693 against the true -1.609 -- an error of 0.916.
# At br = 0.01, 75.8% of study 04's replications carry the correction.
#
# The rule is now consistency-gated: a whole-count (or, on the half-integer
# signature, a corrected) table is taken only if its OR rounds back to the reported
# OR at the reported precision -- i.e. it could be the source of the printed number.
# Otherwise the exact root is kept, so a reported OR that no integer table reproduces
# is no longer replaced by a neighbour with a different OR (see the "reported odds
# ratio" test above). Measured: integer tables with a 2 dp OR 97.3% exact (residual
# errors are the ambiguous small-OR cases where two tables round to the same 2 dp
# value), 3 dp OR 100%, corrected tables 100%.
# =============================================================================

test_that("a continuity-corrected table is reconstructed exactly, not rounded away", {
  ac <- 0.5; bc <- 50.5; cc <- 2.5; dc <- 48.5      # the worked case
  or <- (ac * dc) / (bc * cc)
  got <- .solve(or, ac + bc, cc + dc, ac + cc)
  expect_false(is.null(got))
  expect_equal(c(got$n_cases_exp, got$n_controls_exp, got$n_cases_nexp, got$n_controls_nexp),
               c(ac, bc, cc, dc))
  true_lrr <- log((ac / (ac + bc)) / (cc / (cc + dc)))
  got_lrr <- log((got$n_cases_exp / (got$n_cases_exp + got$n_controls_exp)) /
                 (got$n_cases_nexp / (got$n_cases_nexp + got$n_controls_nexp)))
  expect_equal(got_lrr, true_lrr)
  # the pre-fix answer, pinned so a revert is unmistakable
  expect_false(isTRUE(all.equal(got$n_cases_exp, 1)))
})

test_that("the half-integer guard admits a corrected 0.5 cell but still refuses a zero", {
  # 0.5 is the CORRECTED value of a legitimate zero: its variance is finite, so the
  # old ">= 1" rule would have discarded exactly the tables this branch exists for.
  ac <- 0.5; bc <- 30.5; cc <- 1.5; dc <- 29.5
  expect_false(is.null(.solve((ac * dc) / (bc * cc), ac + bc, cc + dc, ac + cc)))
  # a genuine zero margin is still handed back to the enumeration
  expect_null(.solve(2, 50, 50, 0))
})

test_that("ordinary integer tables are unaffected by the new rule", {
  # The reason round() is there at all: a realistically rounded OR.
  for (tab in list(c(12, 38, 6, 44), c(19, 131, 75, 75), c(40, 60, 60, 40),
                   c(3, 97, 1, 99), c(55, 45, 30, 70))) {
    a <- tab[1]; b <- tab[2]; c <- tab[3]; d <- tab[4]
    got <- .solve(round(.tab_or(a, b, c, d), 2), a + b, c + d, a + c)
    expect_false(is.null(got), info = paste(tab, collapse = "/"))
    expect_equal(got$n_cases_exp, a, info = paste(tab, collapse = "/"))
    expect_equal(got$n_cases_nexp, c, info = paste(tab, collapse = "/"))
  }
})

test_that("recovery rates hold on both populations, measured not assumed", {
  set.seed(1101)
  # corrected tables: must be exact
  err <- solved <- 0; n_b <- 0
  for (i in 1:200) {
    ne <- sample(30:150, 1); nn <- sample(30:150, 1); c0 <- sample(1:6, 1)
    ac <- 0.5; bc <- ne + 0.5; cc <- c0 + 0.5; dc <- nn - c0 + 0.5
    or <- (ac * dc) / (bc * cc)
    got <- .solve(or, ac + bc, cc + dc, ac + cc)
    n_b <- n_b + 1
    if (is.null(got)) next
    solved <- solved + 1
    tru <- log((ac / (ac + bc)) / (cc / (cc + dc)))
    g <- log((got$n_cases_exp / (got$n_cases_exp + got$n_controls_exp)) /
             (got$n_cases_nexp / (got$n_cases_nexp + got$n_controls_nexp)))
    err <- max(err, abs(g - tru))
  }
  expect_equal(solved, n_b)          # every corrected table is now solvable
  expect_lt(err, 1e-8)               # and exactly so

  # integer tables with a 2 dp OR: recovery must stay high
  hit <- k <- 0
  for (i in 1:300) {
    ne <- sample(20:200, 1); nn <- sample(20:200, 1)
    a <- sample(2:(ne - 2), 1); c <- sample(2:(nn - 2), 1)
    got <- .solve(round(.tab_or(a, ne - a, c, nn - c), 2), ne, nn, a + c)
    if (is.null(got)) next
    k <- k + 1
    hit <- hit + (got$n_cases_exp == a)
  }
  expect_gt(k, 250)
  expect_gt(hit / k, 0.90)           # measured 0.973; 0.90 leaves Monte Carlo room

  # integer tables with a 3 dp OR: the printed number now pins the table
  hit3 <- k3 <- 0
  for (i in 1:300) {
    ne <- sample(20:200, 1); nn <- sample(20:200, 1)
    a <- sample(2:(ne - 2), 1); c <- sample(2:(nn - 2), 1)
    got <- .solve(round(.tab_or(a, ne - a, c, nn - c), 3), ne, nn, a + c)
    if (is.null(got)) next
    k3 <- k3 + 1
    hit3 <- hit3 + (got$n_cases_exp == a)
  }
  expect_equal(hit3, k3)
})

test_that("es_from_or_se() end-to-end on a corrected table gives the right RR", {
  ac <- 0.5; bc <- 50.5; cc <- 2.5; dc <- 48.5
  or <- (ac * dc) / (bc * cc)
  se <- sqrt(1/ac + 1/bc + 1/cc + 1/dc)
  true_rr <- (ac / (ac + bc)) / (cc / (cc + dc))
  got <- es_from_or_se(or = or, logor_se = se,
                       n_exp = ac + bc, n_nexp = cc + dc,
                       n_cases = ac + cc, n_controls = bc + dc,
                       or_to_rr = "metaumbrella_cases")
  expect_equal(exp(got$logrr), true_rr, tolerance = 1e-8)
})
