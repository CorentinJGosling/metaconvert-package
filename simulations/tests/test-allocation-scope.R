## =============================================================================
## Roadmap item 4.2 -- studies 04 and 05 are measured at 1:1 allocation ONLY, and
## the README must say so.
##
## WHAT THE ITEM ORIGINALLY CLAIMED, AND WHY THAT PART IS DEAD. The roadmap said
## `metaumbrella_exp`'s coverage of 0.497 was "conditional on p_exp = 0.5", because
## equal arms are the non-identified case for that reconstruction (item 1.1). After
## 1.1 / 1.8 / 1.10 the route SOLVES the table rather than searching for it, and the
## solve is exact whatever the arms look like -- so 0.497 no longer exists (the
## shipped aggregate has 0.976, bit-identical to `metaumbrella_cases`) and there is
## no reporting correction left to make. That half of the item is retired, not fixed.
##
## WHAT SURVIVES IS A DESIGN LIMITATION, AND IT POINTS AT A DIFFERENT FINDING.
## `p_exp = 0.5` is not a neutral choice: it is the allocation most FAVOURABLE to the
## sparsity diagnostic of item 4.1. P(any zero cell) depends on the two arms
## separately and is driven by the smaller one, so balance maximises the dense region
## -- the region in which that diagnostic finds every route calibrated. Everything
## downstream of sparsity (se_ratio, coverage, and the non-estimability of the two
## routes that carry the selection caveat) is therefore conditional on 1:1, while the
## point-estimate no-go map is not.
##
## These tests pin BOTH halves: the premise of the caveat (the grid really does hold
## one allocation), the fact that makes it worth stating (the dense region moves, and
## 0.5 is its maximum), the mechanism (the corrected-table share is minimised at
## balance, and that is what `metafor::conv.2x2` fails on), and the part that is NOT
## conditional (the reconstruction is exact at every allocation).
##
## IF SOMEONE LATER ADDS `p_exp` TO `build_grid_04()`: the first test below fails on
## purpose. Delete it, delete the README caveat it guards, and re-point the rest of
## this file at the regenerated aggregates.
## =============================================================================

test_that("the caveat's premise holds: studies 04 and 05 run at one allocation", {
  g <- suppressMessages(build_grid_04())
  expect_equal(sort(unique(g$p_exp)), 0.5)

  # and the SHIPPED aggregates agree, so a stale file cannot quietly contradict the
  # README paragraph
  for (f in c("04_or_to_rr_nrep1000.csv", "05_rr_to_or_nrep1000.csv")) {
    agg <- utils::read.csv(file.path(.sim_root, "data", "aggregated", f))
    expect_equal(sort(unique(agg$p_exp)), 0.5, info = f)
  }
})

test_that("the caveat is not vacuous: the dense region moves with allocation", {
  # This is the whole reason the paragraph exists. If the dense/sparse split were
  # allocation-invariant, holding p_exp fixed would cost nothing and the README
  # would be adding a hedge with no content behind it.
  g <- suppressMessages(build_grid_04())
  r <- dense_region_by_allocation(g)

  # 0.5 is the MAXIMUM of the dense region, not a middling value
  expect_equal(r$p_exp[which.max(r$dense)], 0.5)
  # the shipped split, reproduced through the new helper
  expect_equal(r$dense[r$p_exp == 0.5], 105L)
  expect_equal(r$sparse[r$p_exp == 0.5], 255L)
  # and it falls off in BOTH directions -- the two are not interchangeable, which is
  # why a single unbalanced level would not have covered this either
  expect_equal(r$dense[r$p_exp == 0.25], 65L)
  expect_equal(r$dense[r$p_exp == 0.75], 85L)
  expect_lt(r$dense[r$p_exp == 0.10], 40L)
  expect_lt(r$dense[r$p_exp == 0.90], 60L)
  # a third of the dense region is lost by moving to 3:1 in the thinner direction
  expect_lt(r$dense[r$p_exp == 0.25] / r$dense[r$p_exp == 0.5], 0.7)
})

test_that("the mechanism is the corrected-table share, and it is minimised at 1:1", {
  # Closed form: with arms n1 and n0 there are (n1+1)(n0+1) tables and (n1-1)(n0-1)
  # of them are strictly interior. The interior product is maximised at n1 = n0 for a
  # fixed total, so the share reaching a route ALREADY continuity-corrected is
  # smallest at balance.
  for (n in c(50, 300)) {
    s <- corrected_table_share(n, c(0.10, 0.25, 0.50, 0.75, 0.90))
    expect_equal(which.min(s), 3L, info = paste("n =", n))
    expect_true(all(diff(s[1:3]) < 0) && all(diff(s[3:5]) > 0), info = paste("n =", n))
  }
  # the numbers the README quotes
  expect_equal(round(corrected_table_share(50, 0.50), 4), 0.1479)
  expect_equal(round(corrected_table_share(50, 0.10), 4), 0.3623)
})

test_that("conv.2x2 fails only on corrected tables, which is what ties it to allocation", {
  skip_if_not_installed("metafor")
  # Deterministic: enumerate EVERY 2x2 with total n at an allocation, apply the
  # Haldane correction where a cell is 0, and ask conv.2x2 to reconstruct it. No
  # seed, no sampling -- the whole table space.
  enum <- function(n, p_exp) {
    n1 <- round(n * p_exp); n0 <- n - n1
    g <- expand.grid(a = 0:n1, c = 0:n0)
    g$b <- n1 - g$a; g$d <- n0 - g$c
    zero <- g$a == 0 | g$b == 0 | g$c == 0 | g$d == 0
    a <- g$a + .5 * zero; b <- g$b + .5 * zero
    cc <- g$c + .5 * zero; d <- g$d + .5 * zero
    or <- (a * d) / (b * cc)
    keep <- is.finite(or) & or > 0
    a <- a[keep]; b <- b[keep]; cc <- cc[keep]; d <- d[keep]
    or <- or[keep]; zero <- zero[keep]
    rec <- suppressWarnings(metafor::conv.2x2(ori = or, ni = a + b + cc + d,
                                              n1i = a + b, n2i = a + cc))
    bad <- is.na(rec$ai) | rec$ai <= 0 | rec$bi <= 0 | rec$ci <= 0 | rec$di <= 0
    c(all = mean(bad), corrected = mean(bad[zero]), integer = mean(bad[!zero]))
  }

  for (p in c(0.10, 0.50, 0.90)) {
    r <- enum(50, p)
    # it reconstructs EVERY integer table ...
    expect_equal(unname(r["integer"]), 0, info = paste("p_exp =", p))
    # ... and fails on roughly half the corrected ones, at every allocation
    expect_gt(unname(r["corrected"]), 0.3)
    expect_lt(unname(r["corrected"]), 0.8)
  }
  # so its overall failure rate is the corrected share times a constant, and is
  # therefore lowest exactly where this study measured it
  f <- vapply(c(0.10, 0.50, 0.90), function(p) enum(50, p)[["all"]], numeric(1))
  expect_equal(which.min(f), 2L)
  expect_gt(f[1] / f[2], 2)
  expect_gt(f[3] / f[2], 2)
})

test_that("what is NOT conditional: the 2x2 solve is exact at every allocation", {
  skip_if_not(isTRUE(get0("METACONVERT_AVAILABLE", ifnotfound = FALSE)),
              "metaConvert not loadable")
  # The original 4.2 worry was that equal arms are the non-identified case. With both
  # margins supplied -- which study 04 always does -- the OR determines the table by a
  # quadratic and the tie never arises, so the point estimates are allocation-free.
  # This is the sentence in the README that says the no-go MAP is not conditional.
  set.seed(4020)
  for (p_exp in c(0.10, 0.25, 0.50, 0.75, 0.90)) {
    n <- 300; n1 <- round(n * p_exp); n0 <- n - n1
    err <- numeric(0)
    for (i in 1:40) {
      br <- stats::runif(1, 0.05, 0.6)
      rr <- exp(stats::runif(1, log(0.3), log(2)))
      if (rr * br >= 1) next
      a <- stats::rbinom(1, n1, rr * br); c_ <- stats::rbinom(1, n0, br)
      b <- n1 - a; d <- n0 - c_
      zero <- a == 0 || b == 0 || c_ == 0 || d == 0
      a <- a + .5 * zero; b <- b + .5 * zero
      c_ <- c_ + .5 * zero; d <- d + .5 * zero
      truth <- log((a / (a + b)) / (c_ / (c_ + d)))
      for (m in c("metaumbrella_cases", "metaumbrella_exp")) {
        got <- suppressWarnings(es_from_or_se(
          or = (a * d) / (b * c_), logor_se = sqrt(1/a + 1/b + 1/c_ + 1/d),
          n_exp = a + b, n_nexp = c_ + d, n_cases = a + c_, n_controls = b + d,
          n_sample = n, baseline_risk = c_ / (c_ + d), or_to_rr = m))
        err <- c(err, abs(got$logrr - truth))
      }
    }
    expect_true(all(is.finite(err)), info = paste("p_exp =", p_exp))
    expect_lt(max(err), 1e-8, label = paste("max |log RR error| at p_exp =", p_exp))
  }
})

## --- stopifnot() fallback for a bare R (no testthat) ------------------------
if (!requireNamespace("testthat", quietly = TRUE)) {
  .g <- suppressMessages(build_grid_04())
  .r <- dense_region_by_allocation(.g)
  stopifnot(
    identical(sort(unique(.g$p_exp)), 0.5),
    .r$p_exp[which.max(.r$dense)] == 0.5,
    .r$dense[.r$p_exp == 0.5] == 105L,
    which.min(corrected_table_share(50, c(0.10, 0.25, 0.50, 0.75, 0.90))) == 3L
  )
}
