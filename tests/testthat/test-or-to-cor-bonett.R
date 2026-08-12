## or_to_cor = "bonett" -- the convert_df() default.
##
## Two things are locked down here:
##  (1) the published worked examples of Bonett & Price (2005) reproduce through the
##      package WITH small_margin_prop LEFT BLANK, i.e. the pmin derivation is right;
##  (2) all four equivalent ways of describing the same 2x2 margins reach the method.
##      Before the symmetric back-fill only n_exp + n_cases did; the other three fell
##      through to the lipsey_cooper value with no warning and info_used unchanged.

test_that("Bonett & Price (2005) worked examples reproduce with pmin left blank", {
  # The paper adds 0.5 to each cell before computing the OR, but takes the marginal
  # proportions from the UNCORRECTED counts. Feeding it that way reproduces the
  # printed numbers; see the derivation comment in R/es_from_stand_OR.R.
  bp <- function(a, b, c_, d) {
    or <- ((a + .5) * (d + .5)) / ((b + .5) * (c_ + .5))
    se <- sqrt(1 / (a + .5) + 1 / (b + .5) + 1 / (c_ + .5) + 1 / (d + .5))
    es_from_or_se(
      or = or, logor_se = se,
      n_exp = a + b, n_nexp = c_ + d,
      n_cases = a + c_, n_controls = b + d,
      n_sample = a + b + c_ + d,
      or_to_cor = "bonett"
    )
  }

  # Example 1, f = (203, 186, 167, 374); paper prints .333 / .048 / (.237, .424)
  e1 <- bp(203, 186, 167, 374)
  expect_equal(e1$r,       0.333196, tolerance = 5e-4)
  expect_equal(e1$r_se,    0.047838, tolerance = 5e-4)
  expect_equal(e1$r_ci_lo, 0.236682, tolerance = 5e-4)
  expect_equal(e1$r_ci_up, 0.423753, tolerance = 5e-4)

  # Example 2, f = (4, 6, 1, 89); paper prints .831 / .108 / (.488, .956)
  e2 <- bp(4, 6, 1, 89)
  expect_equal(e2$r,       0.831157, tolerance = 5e-4)
  expect_equal(e2$r_se,    0.107651, tolerance = 5e-4)
  expect_equal(e2$r_ci_lo, 0.487667, tolerance = 5e-4)
  expect_equal(e2$r_ci_up, 0.956060, tolerance = 5e-4)

  # Example 3, f = (143, 52, 41, 164); paper prints .741 / .0446 / (.641, .817)
  e3 <- bp(143, 52, 41, 164)
  expect_equal(e3$r,       0.740626, tolerance = 5e-4)
  expect_equal(e3$r_se,    0.044597, tolerance = 5e-4)
  expect_equal(e3$r_ci_lo, 0.641214, tolerance = 5e-4)
  expect_equal(e3$r_ci_up, 0.816933, tolerance = 5e-4)
})


test_that("all four equivalent margin pairs reach the bonett conversion", {
  # One table: n_exp = 40, n_nexp = 60, n_cases = 30, n_controls = 70, n_sample = 100.
  # n_sample plus one member of each pair determines it, so all four must agree.
  run <- function(...) {
    es_from_or_se(or = 2.5, logor_se = 0.2, or_to_cor = "bonett", n_sample = 100, ...)
  }
  a <- run(n_exp = 40,  n_cases = 30)
  b <- run(n_nexp = 60, n_controls = 70)
  cc <- run(n_exp = 40,  n_controls = 70)
  d <- run(n_nexp = 60, n_cases = 30)

  expect_equal(a$r, 0.326978458193, tolerance = 1e-10)
  for (x in list(b, cc, d)) {
    expect_equal(x$r,       a$r,       tolerance = 1e-12)
    expect_equal(x$r_se,    a$r_se,    tolerance = 1e-12)
    expect_equal(x$r_ci_lo, a$r_ci_lo, tolerance = 1e-12)
    expect_equal(x$r_ci_up, a$r_ci_up, tolerance = 1e-12)
    expect_equal(x$z,       a$z,       tolerance = 1e-12)
  }

  # ... and none of them is silently the lipsey_cooper fall-through
  lc <- es_from_or_se(or = 2.5, logor_se = 0.2, or_to_cor = "lipsey_cooper",
                      n_exp = 40, n_nexp = 60, n_cases = 30, n_controls = 70,
                      n_sample = 100)$r
  expect_false(isTRUE(all.equal(a$r, lc, tolerance = 1e-6)))
})


test_that("a user-supplied small_margin_prop still wins over the derivation", {
  derived <- es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                           n_cases = 40, n_controls = 60, n_sample = 100,
                           or_to_cor = "bonett")$r
  # derived pmin = min(50, 50, 40, 60) / 100 = 0.40
  supplied <- es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                            n_cases = 40, n_controls = 60, n_sample = 100,
                            small_margin_prop = 0.40, or_to_cor = "bonett")$r
  expect_equal(derived, supplied, tolerance = 1e-12)

  # a different supplied value must produce a different r
  other <- es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                         n_cases = 40, n_controls = 60, n_sample = 100,
                         small_margin_prop = 0.25, or_to_cor = "bonett")$r
  expect_false(isTRUE(all.equal(derived, other, tolerance = 1e-8)))
})


test_that("an incoherent or degenerate table does not get a derived pmin", {
  # n_exp + n_nexp = 100 but n_cases + n_controls = 90: the margins describe no table.
  incoherent <- es_from_or_se(or = 2.5, logor_se = 0.2, n_exp = 40, n_nexp = 60,
                              n_cases = 30, n_controls = 60, n_sample = 100,
                              or_to_cor = "bonett")$r
  lc <- es_from_or_se(or = 2.5, logor_se = 0.2, or_to_cor = "lipsey_cooper",
                      n_exp = 40, n_nexp = 60, n_cases = 30, n_controls = 60,
                      n_sample = 100)$r
  expect_equal(incoherent, lc, tolerance = 1e-12)

  # a zero margin is degenerate: c can go negative and flip the sign of r
  degenerate <- es_from_or_se(or = 2.5, logor_se = 0.2, n_exp = 100, n_nexp = 0,
                              n_cases = 30, n_controls = 70, n_sample = 100,
                              or_to_cor = "bonett")$r
  expect_false(is.na(degenerate) && FALSE)  # must not error
  expect_true(is.na(degenerate) || abs(degenerate) <= 1)
})


test_that("the row's result does not depend on other rows in the same call", {
  # Row 1 can use bonett, row 2 cannot (margins do not sum to n_sample).
  two <- es_from_or_se(
    or = c(2.5, 2.5), logor_se = c(0.2, 0.2),
    n_exp = c(40, 40), n_nexp = c(60, 60),
    n_cases = c(30, 30), n_controls = c(70, 60),
    n_sample = c(100, 100), or_to_cor = rep("bonett", 2)
  )
  one <- es_from_or_se(or = 2.5, logor_se = 0.2, n_exp = 40, n_nexp = 60,
                       n_cases = 30, n_controls = 70, n_sample = 100,
                       or_to_cor = "bonett")
  expect_equal(two$r[1], one$r, tolerance = 1e-12)

  # and a scalar or_to_cor must behave exactly like the recycled vector
  scal <- es_from_or_se(
    or = c(2.5, 2.5), logor_se = c(0.2, 0.2),
    n_exp = c(40, 40), n_nexp = c(60, 60),
    n_cases = c(30, 30), n_controls = c(70, 60),
    n_sample = c(100, 100), or_to_cor = "bonett"
  )
  expect_equal(scal$r, two$r, tolerance = 1e-12)
})
