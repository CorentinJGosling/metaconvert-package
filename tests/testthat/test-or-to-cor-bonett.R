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
  # The fallback target changed from lipsey_cooper to digby (roadmap 2.7); what this
  # test is about -- an incoherent table must not get a derived pmin, so bonett must
  # not run -- is unchanged.
  dig <- es_from_or_se(or = 2.5, logor_se = 0.2, or_to_cor = "digby",
                       n_exp = 40, n_nexp = 60, n_cases = 30, n_controls = 60,
                       n_sample = 100)$r
  expect_equal(incoherent, dig, tolerance = 1e-12)

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


## ---------------------------------------------------------------------------
## The fallback message used to offer 'small_margin_prop' as an alternative to the
## margins. It is not one: bonett's coefficient is
##   c = (1 - |n_exp/n_sample - n_cases/n_sample|/5 - (1/2 - small_margin_prop)^2) / 2
## which reads n_exp, n_cases AND n_sample as well, so the eligibility gate is a
## CONJUNCTION and the message was describing a disjunction. A user who followed it
## got no bonett estimate and no further explanation -- the row silently kept the
## stand-in computed earlier (lipsey_cooper before roadmap 2.7, digby since).
## ---------------------------------------------------------------------------

test_that("small_margin_prop alone does NOT make a row eligible for bonett", {
  # Everything bonett needs except n_sample and the case margin.
  got <- suppressMessages(es_from_or_se(
    or = 2.5, logor_se = 0.2, n_exp = 40, n_nexp = 60,
    small_margin_prop = 0.40, or_to_cor = "bonett")$r)
  dig <- es_from_or_se(or = 2.5, logor_se = 0.2, n_exp = 40, n_nexp = 60,
                       or_to_cor = "digby")$r
  # It fell back: the requested method did not run. (The fallback target became digby
  # in roadmap 2.7; the point here is that bonett did not run, not which stand-in ran.)
  expect_equal(got, dig, tolerance = 1e-12)

  # Completing the margin set is what actually enables it, and it moves the answer.
  ok <- suppressMessages(es_from_or_se(
    or = 2.5, logor_se = 0.2, n_exp = 40, n_nexp = 60,
    n_cases = 30, n_controls = 70, n_sample = 100,
    small_margin_prop = 0.40, or_to_cor = "bonett")$r)
  expect_false(isTRUE(all.equal(ok, dig, tolerance = 1e-8)))
})


test_that("the fallback message does not offer small_margin_prop as a substitute", {
  msg <- NULL
  withCallingHandlers(
    invisible(es_from_or_se(or = 2.5, logor_se = 0.2, n_exp = 40, n_nexp = 60,
                            small_margin_prop = 0.40, or_to_cor = "bonett")),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
  msg <- paste(msg, collapse = " ")
  expect_true(nzchar(msg))

  # The defect, stated exactly: the message used to read
  #   "For 'bonett', supply 'small_margin_prop', or 'n_sample' together with ..."
  # i.e. it named as sufficient the one input that cannot work on its own.
  expect_false(grepl("supply 'small_margin_prop'", msg, fixed = TRUE))

  # What it must name instead -- the conjunction the gate really tests.
  expect_true(grepl("n_sample", msg, fixed = TRUE))
  expect_true(grepl("n_exp", msg, fixed = TRUE))
  expect_true(grepl("n_cases", msg, fixed = TRUE))

  # And it must say plainly that small_margin_prop on its own is not enough,
  # otherwise a reader re-derives the same wrong conclusion from its presence.
  expect_true(grepl("small_margin_prop", msg, fixed = TRUE))
  expect_true(grepl("on its own", msg, fixed = TRUE))
})


## ---------------------------------------------------------------------------
## WHICH conversion an ineligible row falls back to (roadmap 2.7).
##
## It used to keep the "lipsey_cooper" R/Z that .es_from_d() had already computed --
## not a choice, just the value that happened to be there. Scored against the
## tetrachoric over 1176 coherent 2x2 tables that is the WORST of the four options
## (mean |error| 0.0950 vs digby 0.0169), and it needs MORE information than digby,
## which reads no margins at all.
## ---------------------------------------------------------------------------

test_that("an ineligible bonett row falls back to digby, not to lipsey_cooper", {
  got <- suppressMessages(es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                                        or_to_cor = "bonett"))
  dig <- es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                       or_to_cor = "digby")
  lc  <- es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                       or_to_cor = "lipsey_cooper")
  expect_equal(got$r, dig$r, tolerance = 1e-12)
  expect_equal(got$z, dig$z, tolerance = 1e-12)
  expect_equal(got$r_se, dig$r_se, tolerance = 1e-12)
  # and it is genuinely a different number from the old fallback
  expect_false(isTRUE(all.equal(got$r, lc$r, tolerance = 1e-8)))
})


test_that("the fallback now works on rows that used to come back empty", {
  # digby's coefficient c = 3/4 is a constant, so it reads no margins. lipsey_cooper
  # needs the arm sizes and returned NA without them, so a row carrying only an odds
  # ratio and its SE produced no correlation at all. This WIDENS coverage.
  got <- suppressMessages(es_from_or_se(or = 2, logor_se = 0.2, or_to_cor = "bonett"))
  expect_false(is.na(got$r))
  expect_equal(got$r, es_from_or_se(or = 2, logor_se = 0.2, or_to_cor = "digby")$r,
               tolerance = 1e-12)
})


test_that("a user who explicitly asks for lipsey_cooper still gets lipsey_cooper", {
  # "Do not remove what the user explicitly asked for" (roadmap 1.1). Rows whose
  # or_to_cor IS lipsey_cooper are excluded from the fallback set by construction.
  lc <- es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                      or_to_cor = "lipsey_cooper")
  expect_equal(round(lc$r, 8), 0.18768063)
  expect_silent(es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                              or_to_cor = "lipsey_cooper"))
})


test_that("the substitution message names the method actually used", {
  msg <- NULL
  withCallingHandlers(
    invisible(es_from_or_se(or = 2, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                            or_to_cor = "bonett")),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
  msg <- paste(msg, collapse = " ")
  expect_true(grepl("obtained with 'digby' instead", msg, fixed = TRUE))
  expect_false(grepl("obtained with 'lipsey_cooper' instead", msg, fixed = TRUE))
})
