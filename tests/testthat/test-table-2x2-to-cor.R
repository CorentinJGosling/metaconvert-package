# =============================================================================
# Roadmap item 1.3 -- table_2x2_to_cor, and the tetrachoric route around it.
#
# The argument was accepted and SILENTLY IGNORED by convert_df(): commented out of
# the per-row method-column loop and of all three es_from_2x2*() call sites, so
# convert_df(x, measure = "r", table_2x2_to_cor = "banana") ran without error and
# returned byte-identical output. The guard fired only on a direct es_from_2x2()
# call, and its message advertised three methods ('cooper_delta', 'cooper_std',
# 'lipsey') that no reachable code implemented.
#
# Only "tetrachoric" is offered, and that is BY DESIGN. phi -- the obvious
# alternative -- has an attainable range bounded by the margins, so studies of the
# same association with different event rates report different phi values and
# pooling them manufactures heterogeneity that is pure margin artefact. And the
# choice could not be offered responsibly anyway: a 2x2 with fixed n has three free
# parameters and the dichotomised-bivariate-normal family also has three, so the
# latent-normal model is saturated and admits no goodness-of-fit test.
#
# This file also pins the two changes that followed from that decision:
#   * the r-scale CI is now the tanh back-transform of the Fisher-z interval, so it
#     cannot leave the parameter space (measured escape 45.8% -> 0.0%);
#   * V35 discloses that the tetrachoric is weakly identified at extreme margins.
# =============================================================================

.tab <- function() data.frame(
  study_id = "s1",
  n_cases_exp = 30, n_controls_exp = 15,
  n_cases_nexp = 20, n_controls_nexp = 35
)

# --- validation ------------------------------------------------------------

test_that("convert_df() rejects an invalid table_2x2_to_cor instead of ignoring it", {
  expect_error(
    suppressMessages(convert_df(.tab(), measure = "r",
                                table_2x2_to_cor = "banana", verbose = FALSE)),
    "table_2x2_to_cor"
  )
  # The historical failure mode: no error at all.
  expect_error(
    suppressMessages(convert_df(.tab(), measure = "r",
                                table_2x2_to_cor = "lipsey", verbose = FALSE)),
    "tetrachoric"
  )
})

test_that("es_from_2x2() rejects invalid values with a well-formed message", {
  expect_error(es_from_2x2(30, 15, 20, 35, table_2x2_to_cor = "banana"), "banana")
  # The old message ended "Possible inputs are: 'tetrachoric', " -- a dangling comma
  # implying further options existed.
  msg <- tryCatch(es_from_2x2(30, 15, 20, 35, table_2x2_to_cor = "banana"),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("cooper_delta|cooper_std|lipsey", msg),
               info = "message must not advertise unimplemented methods")
  expect_false(grepl(",\\s*$", msg), info = "message must not end in a dangling comma")
})

test_that("the valid value still works on both paths", {
  expect_no_error(es_from_2x2(30, 15, 20, 35, table_2x2_to_cor = "tetrachoric"))
  expect_no_error(suppressMessages(
    convert_df(.tab(), measure = "r", table_2x2_to_cor = "tetrachoric", verbose = FALSE)))
})

# --- the r-scale CI cannot escape the parameter space ----------------------

test_that("the tetrachoric r CI stays inside [-1, 1]", {
  set.seed(4)
  n_checked <- 0; escaped <- 0
  for (i in 1:300) {
    N  <- sample(c(50, 100, 300), 1)
    pe <- runif(1, .3, .7); pc <- runif(1, .05, .5)
    n1 <- round(N * pe)
    a <- rbinom(1, n1, min(pc * 2.5, .95)); b <- n1 - a
    c <- rbinom(1, N - n1, pc); d <- (N - n1) - c
    if (min(a, b, c, d) < 1) next
    res <- es_from_2x2(a, b, c, d)
    if (is.na(res$r_ci_lo) || is.na(res$r_ci_up)) next
    n_checked <- n_checked + 1
    if (res$r_ci_lo < -1 || res$r_ci_up > 1) escaped <- escaped + 1
  }
  expect_gt(n_checked, 100)
  expect_equal(escaped, 0)
})

test_that("the r CI is the tanh back-transform of the z CI", {
  res <- es_from_2x2(30, 15, 20, 35)
  expect_equal(res$r_ci_lo, tanh(res$z_ci_lo), tolerance = 1e-10)
  expect_equal(res$r_ci_up, tanh(res$z_ci_up), tolerance = 1e-10)
  # ...and it still brackets the point estimate.
  expect_lt(res$r_ci_lo, res$r)
  expect_gt(res$r_ci_up, res$r)
})

test_that("reverse_2x2 reflects the CI rather than inverting it", {
  fwd <- es_from_2x2(30, 15, 20, 35)
  rev <- es_from_2x2(30, 15, 20, 35, reverse_2x2 = TRUE)
  expect_equal(rev$r, -fwd$r, tolerance = 1e-10)
  expect_equal(rev$r_ci_lo, -fwd$r_ci_up, tolerance = 1e-10)
  expect_equal(rev$r_ci_up, -fwd$r_ci_lo, tolerance = 1e-10)
  expect_lt(rev$r_ci_lo, rev$r_ci_up)   # never inverted
})

# --- V35: margin-heterogeneity disclosure ---------------------------------

test_that("V35 fires when the pool spans rare and common event rates", {
  d <- data.frame(
    study_id = c("rare", "common"),
    n_cases_exp = c(4, 60), n_controls_exp = c(96, 40),
    n_cases_nexp = c(2, 40), n_controls_nexp = c(98, 60)
  )
  s <- suppressMessages(summary(
    suppressMessages(convert_df(d, measure = "r", verbose = FALSE)), flags = TRUE))
  expect_true(all(grepl("Event rates across", s$flags_crude)),
              info = "every row in the pool should carry the disclosure")
  expect_true(all(grepl("^\\[INFO\\]|; \\[INFO\\]", s$flags_crude)))
})

test_that("V35 stays silent on a homogeneous mid-range pool", {
  d <- data.frame(
    study_id = c("a", "b"),
    n_cases_exp = c(50, 55), n_controls_exp = c(50, 45),
    n_cases_nexp = c(40, 42), n_controls_nexp = c(60, 58)
  )
  s <- suppressMessages(summary(
    suppressMessages(convert_df(d, measure = "r", verbose = FALSE)), flags = TRUE))
  expect_false(any(grepl("Event rates across", s$flags_crude)))
})

test_that("V35 is confined to correlation measures and honours its threshold", {
  d <- data.frame(
    study_id = c("rare", "common"),
    n_cases_exp = c(4, 60), n_controls_exp = c(96, 40),
    n_cases_nexp = c(2, 40), n_controls_nexp = c(98, 60)
  )
  # Not a correlation pool -> no disclosure (it is about the tetrachoric route).
  s_or <- suppressMessages(summary(
    suppressMessages(convert_df(d, measure = "logor", verbose = FALSE)), flags = TRUE))
  expect_false(any(grepl("Event rates across", s_or$flags_crude)))

  # Threshold configurable.
  s_off <- suppressMessages(summary(
    suppressMessages(convert_df(d, measure = "r", verbose = FALSE,
                                flag_options = list(margin_range_min = 0.99))),
    flags = TRUE))
  expect_false(any(grepl("Event rates across", s_off$flags_crude)))
})

# --- es_from_phi no longer swaps estimands silently ------------------------

test_that("es_from_phi() announces the estimand when the 2x2 cannot be rebuilt", {
  # .notify_once() is session-scoped, so clear the register first or an earlier test
  # in the same run will have consumed the one-time message.
  rm(list = ls(metaConvert:::.mcv_notices), envir = metaConvert:::.mcv_notices)
  expect_message(es_from_phi(phi = 0.2182, n_sample = 200),
                 "phi coefficient itself")
})

test_that("es_from_phi() warns loudly when ONE call mixes both estimands", {
  rm(list = ls(metaConvert:::.mcv_notices), envir = metaConvert:::.mcv_notices)
  expect_message(
    es_from_phi(phi = c(0.2182, 0.2182), n_sample = c(200, 200),
                n_cases = c(60, NA), n_exp = c(100, NA)),
    "TWO DIFFERENT correlation estimands"
  )
})

test_that("the two es_from_phi regimes really do differ, which is why it warns", {
  with_margins <- es_from_phi(phi = 0.2182, n_sample = 200, n_cases = 60, n_exp = 100)
  without      <- suppressMessages(es_from_phi(phi = 0.2182, n_sample = 200))
  expect_equal(without$r, 0.2182, tolerance = 1e-6)      # phi itself
  expect_gt(with_margins$r, 1.4 * without$r)             # tetrachoric, much larger
})
