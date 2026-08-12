# Regression tests for the RR -> OR reconstruction (.metaumbrella_rr_se_to_or).
#
# Guards against the arm-total sign error in var(log RR): the reconstruction picks
# the candidate 2x2 table whose implied variance matches the reported logrr_se^2.
# var(log RR) = 1/a - 1/n1 + 1/c - 1/n2 (arm totals SUBTRACTED) -- NOT the log-OR
# pattern 1/a + 1/b + 1/c + 1/d. With the wrong (added) signs the wrong table is
# selected and the recovered OR is biased (e.g. +9.1% for a=30,b=70,c=10,d=90).
#
# Strategy: take KNOWN 2x2 tables, compute their exact (RR, logrr_se, n_cases,
# n_controls), reconstruct, and assert the recovered OR == the true OR = a*d/(b*c).

test_that("RR->OR (metaumbrella) recovers the true OR of a known 2x2 table", {
  round_trip_or <- function(a, b, c, d) {
    n1 <- a + b; n2 <- c + d
    rr <- (a / n1) / (c / n2)
    logrr_se <- sqrt(1 / a - 1 / n1 + 1 / c - 1 / n2)
    es <- es_from_rr_se(
      rr = rr, logrr_se = logrr_se,
      n_cases = a + c, n_controls = b + d,
      n_exp = n1, n_nexp = n2,
      rr_to_or = "metaumbrella"
    )
    exp(es$logor[1])
  }

  # Canonical motivating case: true OR = (30*90)/(70*10) = 3.857142...
  expect_equal(round_trip_or(30, 70, 10, 90), (30 * 90) / (70 * 10), tolerance = 0.01)

  # A spread of well-conditioned tables (all cells >= 10), both RR>1 and RR<1.
  tables <- list(
    c(40, 60, 20, 80),   # OR = 2.667
    c(25, 75, 50, 50),   # OR < 1 (protective)
    c(60, 40, 30, 70),   # OR = 3.5
    c(15, 85, 45, 55),   # OR < 1
    c(50, 150, 20, 180), # unbalanced margins
    c(90, 110, 40, 160)
  )
  for (tb in tables) {
    a <- tb[1]; b <- tb[2]; c <- tb[3]; d <- tb[4]
    true_or <- (a * d) / (b * c)
    expect_equal(round_trip_or(a, b, c, d), true_or, tolerance = 0.02,
                 info = sprintf("table a=%d b=%d c=%d d=%d (true OR=%.4f)", a, b, c, d, true_or))
  }
})

test_that("RR->OR reconstruction is far more accurate than the pre-fix (added-sign) formula", {
  # Recompute what the OLD buggy var_sim (arm totals ADDED) would have selected, and
  # assert the shipping reconstruction is materially closer to the truth. This pins
  # the direction of the fix so a future regression to the added form is caught.
  buggy_recon <- function(rr, logrr_se, n_cases, n_controls) {
    ncn <- 1:(n_cases - 1)
    ncne <- round(ncn * ((rr * (n_cases + n_controls)) / (n_cases + (rr - 1) * ncn) - 1))
    a <- n_cases - ncn; ce <- n_controls - ncne
    ok <- which(ncn > 0 & ncne > 0 & a > 0 & ce > 0)
    ncn <- ncn[ok]; ncne <- ncne[ok]; a <- a[ok]; ce <- ce[ok]
    n1 <- a + ce; n2 <- ncn + ncne
    var_added <- 1 / a + 1 / n1 + 1 / ncn + 1 / n2          # the OLD (buggy) form
    j <- order((var_added - logrr_se^2)^2)[1]
    (a[j] * ncne[j]) / (ce[j] * ncn[j])                     # OR of the buggy pick
  }

  a <- 60; b <- 40; c <- 30; d <- 70
  n1 <- a + b; n2 <- c + d
  rr <- (a / n1) / (c / n2)
  logrr_se <- sqrt(1 / a - 1 / n1 + 1 / c - 1 / n2)
  true_or <- (a * d) / (b * c)

  fixed_or <- exp(es_from_rr_se(rr = rr, logrr_se = logrr_se, n_cases = a + c,
                                n_controls = b + d, n_exp = n1, n_nexp = n2,
                                rr_to_or = "metaumbrella")$logor[1])
  buggy_or <- buggy_recon(rr, logrr_se, a + c, b + d)

  expect_lt(abs(fixed_or - true_or), abs(buggy_or - true_or))  # fix is strictly closer
  expect_equal(fixed_or, true_or, tolerance = 0.02)
})

test_that("OR->RR (metaumbrella) direction is unaffected and still recovers the true RR", {
  # The sibling OR->2x2 sweep uses var(log OR) = 1/a+1/b+1/c+1/d (all added), which is
  # CORRECT and must not be touched by the RR->OR fix. Guard it here.
  a <- 40; b <- 60; c <- 20; d <- 80
  n1 <- a + b; n2 <- c + d
  or <- (a * d) / (b * c)
  logor_se <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  true_rr <- (a / n1) / (c / n2)
  es <- es_from_or_se(or = or, logor_se = logor_se,
                      n_cases = a + c, n_controls = b + d,
                      n_exp = n1, n_nexp = n2, or_to_rr = "metaumbrella_cases")
  expect_equal(exp(es$logrr[1]), true_rr, tolerance = 0.02)
})
