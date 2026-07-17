# =============================================================================
# Cronbach's alpha (es_from_cronbach_alpha) — promoted from tests_save/checked/
# (P26, audit digest [16]).
#
# DE-CIRCULARIZED: the archived file re-typed the source formulas
# (log(1 - alpha), sqrt(2k/((k-1)(n-2)))) as its own expectations. Every
# numeric expectation below is anchored EXTERNALLY instead:
#   - metafor::escalc(measure = "ABT")  — Bonett (2002) transformed alpha.
#     NOTE the sign convention: metafor stores yi = -ln(1 - alpha) (higher =
#     more reliable) while metaConvert stores +ln(1 - alpha); the variances
#     are identical (a sign flip does not change the sampling variance).
#   - metafor::escalc(measure = "ARAW") — raw alpha with the Feldt et al.
#     (1987) large-sample variance.
#   - Feldt exact-F interval — an independent CI construction the Bonett
#     large-sample interval must approximate (coarse-tolerance sanity anchor).
# =============================================================================

test_that("es_from_cronbach_alpha (bonett) matches metafor escalc ABT (ES sign-flipped, SE exact)", {
  skip_if_not_installed("metafor")
  alphas <- c(0.85, 0.70, 0.95, 0.40)
  ns     <- c(200, 100, 50, 400)
  ks     <- c(10, 5, 20, 8)

  res <- es_from_cronbach_alpha(cronbach_alpha = alphas, n_sample = ns, n_items = ks)
  mfr <- metafor::escalc(measure = "ABT", ai = alphas, mi = ks, ni = ns)

  # metafor: yi = -ln(1 - alpha); metaConvert: +ln(1 - alpha)
  expect_equal(res$alpha, -as.numeric(mfr$yi), tolerance = 1e-10)
  expect_equal(res$alpha_se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
  expect_equal(res$info_used, rep("cronbach_alpha", 4))
})

test_that("es_from_cronbach_alpha (raw) matches metafor escalc ARAW (ES + SE exact)", {
  skip_if_not_installed("metafor")
  alphas <- c(0.85, 0.70, 0.95)
  ns     <- c(200, 100, 50)
  ks     <- c(10, 5, 20)

  res <- es_from_cronbach_alpha(cronbach_alpha = alphas, n_sample = ns, n_items = ks,
                                alpha_to_es = "raw")
  mfr <- metafor::escalc(measure = "ARAW", ai = alphas, mi = ks, ni = ns)

  expect_equal(res$alpha, as.numeric(mfr$yi), tolerance = 1e-10)
  expect_equal(res$alpha_se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-10)
})

test_that("back-transformed Bonett CI approximates the Feldt exact-F interval", {
  # Independent anchor: for a single-sample alpha, the exact interval is
  # 1 - (1 - alpha) * qf(p, n - 1, (n - 1)(k - 1))  (Feldt 1965; Feldt,
  # Woodruff, & Salih 1987). Bonett's large-sample interval is an
  # approximation to it, so the back-transformed bounds (with the swap:
  # upper alpha bound comes from the LOWER transformed bound) must sit close
  # at n = 200 — a coarse-tolerance sanity check, not a bit-exact identity.
  a <- 0.85; n <- 200; k <- 10
  res <- es_from_cronbach_alpha(cronbach_alpha = a, n_sample = n, n_items = k)

  # package CI is on the ln(1 - alpha) scale; back-transform with bound swap
  alpha_lo <- 1 - exp(res$alpha_ci_up)
  alpha_up <- 1 - exp(res$alpha_ci_lo)

  feldt_lo <- 1 - (1 - a) * qf(0.975, n - 1, (n - 1) * (k - 1))
  feldt_up <- 1 - (1 - a) * qf(0.025, n - 1, (n - 1) * (k - 1))
  # frozen values of the Feldt bounds (computed externally):
  expect_equal(feldt_lo, 0.8169948416075911, tolerance = 1e-12)
  expect_equal(feldt_up, 0.8791784502404276, tolerance = 1e-12)

  expect_lt(abs(alpha_lo - feldt_lo), 0.01)
  expect_lt(abs(alpha_up - feldt_up), 0.01)
  expect_true(alpha_lo < a && a < alpha_up)
})

test_that("es_from_cronbach_alpha validates alpha_to_es parameter", {
  expect_error(
    es_from_cronbach_alpha(cronbach_alpha = 0.85, n_sample = 200, n_items = 10,
                           alpha_to_es = "invalid"),
    "not in tolerated values"
  )
})

test_that("es_from_cronbach_alpha handles edge cases", {
  # Alpha = 0: ln(1 - 0) = 0 on the Bonett scale (hand identity, not a
  # re-typed source formula)
  res0 <- es_from_cronbach_alpha(cronbach_alpha = 0, n_sample = 100, n_items = 5)
  expect_equal(res0$alpha, 0, tolerance = 1e-10)

  # Alpha very close to 1 (should give large negative Bonett; ln(0.01) = -4.6)
  res_high <- es_from_cronbach_alpha(cronbach_alpha = 0.99, n_sample = 100, n_items = 10)
  expect_true(res_high$alpha < -4)

  # Alpha = 1 (boundary): under the Bonett transform ln(1 - 1) is undefined,
  # so the per-element guard (P16) now degrades the row to NA instead of
  # emitting -Inf (audit digest [17], fixed).
  res1 <- es_from_cronbach_alpha(cronbach_alpha = 1, n_sample = 100, n_items = 10)
  expect_true(is.na(res1$alpha))
  expect_true(is.na(res1$alpha_se))
})

test_that("es_from_cronbach_alpha handles vectorized input", {
  res <- es_from_cronbach_alpha(
    cronbach_alpha = c(0.80, 0.85, 0.90),
    n_sample = c(100, 200, 300),
    n_items = c(5, 10, 20)
  )
  expect_equal(nrow(res), 3)
  expect_equal(res$info_used, rep("cronbach_alpha", 3))
})

test_that("convert_df with measure='alpha' and default Bonett method matches metafor ABT", {
  skip_if_not_installed("metafor")
  dat <- data.frame(
    cronbach_alpha = c(0.85, 0.70),
    n_sample = c(200, 100),
    n_items = c(10, 5)
  )
  mc <- convert_df(dat, measure = "alpha", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  mfr <- metafor::escalc(measure = "ABT", ai = dat$cronbach_alpha,
                         mi = dat$n_items, ni = dat$n_sample)
  # summary() rounds to `digits` decimal places, so pipeline comparisons use
  # a slightly looser tolerance than the direct-call bit-match tests above
  expect_equal(s$es, -as.numeric(mfr$yi), tolerance = 1e-8)
  expect_equal(s$se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-8)
  expect_equal(s$info_used, rep("cronbach_alpha", 2))
})

test_that("convert_df with measure='alpha' and raw method matches metafor ARAW", {
  skip_if_not_installed("metafor")
  dat <- data.frame(
    cronbach_alpha = c(0.85, 0.70),
    n_sample = c(200, 100),
    n_items = c(10, 5)
  )
  mc <- convert_df(dat, measure = "alpha", alpha_to_es = "raw",
                   verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  mfr <- metafor::escalc(measure = "ARAW", ai = dat$cronbach_alpha,
                         mi = dat$n_items, ni = dat$n_sample)
  # summary() rounds to `digits` decimal places (see the bonett pipeline test)
  expect_equal(s$es, as.numeric(mfr$yi), tolerance = 1e-8)
  expect_equal(s$se, sqrt(as.numeric(mfr$vi)), tolerance = 1e-8)
})

test_that("alpha data_extraction_sheet works", {
  dat <- data_extraction_sheet(measure = "alpha", extension = "data.frame", verbose = FALSE)
  expect_true("cronbach_alpha" %in% colnames(dat))
  expect_true("n_items" %in% colnames(dat))
  expect_true("n_sample" %in% colnames(dat))
})
