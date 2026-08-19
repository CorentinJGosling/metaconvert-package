# Regression tests for the four blocking defects found while assessing the
# package for a Cronbach's alpha / McDonald's omega meta-analysis.
#
# Each block fails on the pre-fix tree. The failure modes share a signature --
# a plausible-looking number with no error, no warning and no NA -- so every
# assertion here checks the VALUE, not merely that the call succeeded.

# ---------------------------------------------------------------------------
# FIX 1 -- a scalar n_sample / n_items beside a vector of coefficients NA'd the
# standard error of every row after the first, while leaving the point
# estimates correct (so head() and any is.na(es) check both passed).
# ---------------------------------------------------------------------------
test_that("alpha: scalar n_sample/n_items recycle instead of NA-ing rows 2..k", {
  a <- c(0.72, 0.85, 0.90, 0.78)
  scal <- es_from_cronbach_alpha(cronbach_alpha = a, n_sample = 200, n_items = 10)
  vect <- es_from_cronbach_alpha(cronbach_alpha = a,
                                 n_sample = rep(200, 4), n_items = rep(10, 4))
  expect_equal(scal, vect)
  expect_false(anyNA(scal$alpha_se))
  # the pre-fix tree returned exactly one non-NA SE
  expect_equal(sum(!is.na(scal$alpha_se)), 4L)
})

test_that("icc: scalar n_sample/n_measurements recycle", {
  rho <- c(0.7, 0.8, 0.9)
  scal <- es_from_icc(icc = rho, n_sample = 50, n_measurements = 2)
  vect <- es_from_icc(icc = rho, n_sample = rep(50, 3), n_measurements = rep(2, 3))
  expect_equal(scal, vect)
  expect_false(anyNA(scal$icc_se))
})

test_that("reliability routes reject a genuinely incompatible length", {
  expect_error(es_from_cronbach_alpha(c(.7, .8, .9), c(100, 200), 10), "n_sample")
  expect_error(es_from_cronbach_alpha(c(.7, .8, .9), 100, c(10, 20)), "n_items")
  expect_error(es_from_icc(c(.7, .8, .9), c(50, 60), 2), "n_sample")
  expect_error(es_from_icc(c(.7, .8, .9), 50, c(2, 3)), "n_measurements")
})

test_that("omitted arguments still degrade to NA rather than erroring", {
  expect_equal(nrow(es_from_cronbach_alpha(c(0.7, 0.8))), 2L)
  expect_true(all(is.na(es_from_cronbach_alpha(c(0.7, 0.8))$alpha)))
})

# ---------------------------------------------------------------------------
# FIX 2 -- `digits` rounded es/se/CI IN the returned data.frame, not just the
# printout, so `se` (an inverse-variance weight downstream) lost precision.
# ---------------------------------------------------------------------------
test_that("summary() returns es/se/CI at full precision regardless of digits", {
  d <- data.frame(study_id = "A", cronbach_alpha = 0.995,
                  n_sample = 5000, n_items = 30)
  for (m in c("bonett", "raw")) {
    s3 <- summary(convert_df(d, measure = "alpha", alpha_to_es = m,
                             verbose = FALSE, split_adjusted = FALSE), digits = 3)
    s15 <- summary(convert_df(d, measure = "alpha", alpha_to_es = m,
                              verbose = FALSE, split_adjusted = FALSE), digits = 15)
    expect_identical(s3$es, s15$es, info = m)
    expect_identical(s3$se, s15$se, info = m)
    expect_identical(s3$es_ci_lo, s15$es_ci_lo, info = m)
    expect_identical(s3$es_ci_up, s15$es_ci_up, info = m)
  }
  # the raw-scale SE here is 1.017e-4: rounding to 3 dp made it exactly 0
  s_raw <- summary(convert_df(d, measure = "alpha", alpha_to_es = "raw",
                              verbose = FALSE, split_adjusted = FALSE))
  expect_gt(s_raw$se, 0)
})

test_that("digits still governs the human-readable summary string", {
  d <- data.frame(study_id = "A", cronbach_alpha = 0.85,
                  n_sample = 200, n_items = 10)
  s2 <- summary(convert_df(d, measure = "alpha", verbose = FALSE,
                           split_adjusted = FALSE), digits = 2)
  s6 <- summary(convert_df(d, measure = "alpha", verbose = FALSE,
                           split_adjusted = FALSE), digits = 6)
  expect_false(identical(s2$es_summary, s6$es_summary))
})

test_that("a near-degenerate SE pool now survives metafor::rma()", {
  skip_if_not_installed("metafor")
  d <- data.frame(study_id = paste0("S", 1:4),
                  cronbach_alpha = c(.99, .995, .98, .992),
                  n_sample = c(4000, 5000, 3000, 4500), n_items = 30)
  s <- summary(convert_df(d, measure = "alpha", alpha_to_es = "raw",
                          verbose = FALSE, split_adjusted = FALSE))
  expect_true(all(s$se > 0))
  m <- metafor::rma(yi = s$es, sei = s$se, method = "FE")
  expect_true(is.finite(as.numeric(m$beta[1])))
})

# ---------------------------------------------------------------------------
# FIX 3 -- the user_es_* passthrough wrote a RAW reliability into a
# Bonett-transformed column. An omega of 0.91 landed at +0.91 beside native
# alpha rows near -2.12, carrying ~19x the correct inverse-variance weight
# ((0.08635463 / 0.02)^2 = 18.64 on this fixture).
# ---------------------------------------------------------------------------
mk_user_dat <- function() {
  data.frame(
    study_id = c("native", "via_user"),
    cronbach_alpha = c(0.88, NA), n_sample = c(300, 300), n_items = c(10, 10),
    user_es_crude = c(NA, 0.91), user_se_crude = c(NA, 0.02),
    user_es_original_measure_crude = c(NA, "omega"),
    user_es_target_measure_crude = c(NA, "alpha")
  )
}

test_that("a raw coefficient is refused rather than written into a Bonett pool", {
  s <- suppressWarnings(summary(
    convert_df(mk_user_dat(), measure = "alpha", verbose = FALSE,
               split_adjusted = FALSE), digits = 15))
  expect_true(is.na(s$es[s$study_id == "via_user"]))
  expect_true(is.na(s$se[s$study_id == "via_user"]))
  # the legitimate native row is untouched
  expect_equal(s$es[s$study_id == "native"], log(1 - 0.88), tolerance = 1e-12)
  # and it must not be the old behaviour
  expect_false(isTRUE(all.equal(s$es[s$study_id == "via_user"], 0.91)))
})

test_that("the refusal explains the scale rather than blaming the extraction", {
  # two warnings fire: the pre-existing unknown-type notice, and the new scale
  # refusal. Capture both so neither escapes the test, and pin the wording of
  # the one that carries the diagnosis.
  w <- capture_warnings(
    convert_df(mk_user_dat(), measure = "alpha", verbose = FALSE,
               split_adjusted = FALSE))
  expect_true(any(grepl("analysed on the 'bonett' scale", w, fixed = TRUE)))
  expect_true(any(grepl("set to NA instead of being written onto the", w, fixed = TRUE)))
  # the escape hatch must be described as the dataset-wide switch it is
  expect_true(any(grepl("DATASET-WIDE change of analysis scale", w, fixed = TRUE)))
  expect_true(any(grepl("Unknown user_es_original_measure_crude", w, fixed = TRUE)))
})

test_that("under alpha_to_es = 'raw' the passthrough is legitimate and allowed", {
  s <- suppressWarnings(summary(
    convert_df(mk_user_dat(), measure = "alpha", alpha_to_es = "raw",
               verbose = FALSE, split_adjusted = FALSE), digits = 15))
  expect_equal(s$es[s$study_id == "via_user"], 0.91)
  expect_equal(s$es[s$study_id == "native"], 0.88)
})

test_that("measures whose scale is fixed by the measure name still pass through", {
  d <- data.frame(
    study_id = c("a", "b"), n_exp = 50, n_nexp = 50,
    user_es_crude = c(0.4, 0.5), user_se_crude = c(0.1, 0.1),
    user_es_original_measure_crude = c("unknown_thing", "unknown_thing"),
    user_es_target_measure_crude = c("g", "g")
  )
  s <- suppressWarnings(summary(
    convert_df(d, measure = "g", verbose = FALSE, split_adjusted = FALSE),
    digits = 15))
  expect_equal(s$es, c(0.4, 0.5))
})

test_that(".user_passthrough_blocked keys on the active transform", {
  expect_true(.user_passthrough_blocked("alpha", alpha_to_es = "bonett"))
  expect_false(.user_passthrough_blocked("alpha", alpha_to_es = "raw"))
  expect_true(.user_passthrough_blocked("icc", icc_to_es = "bonett"))
  expect_false(.user_passthrough_blocked("icc", icc_to_es = "raw"))
  expect_true(.user_passthrough_blocked("prop", prop_to_es = "logit"))
  expect_false(.user_passthrough_blocked("prop", prop_to_es = "raw"))
  expect_false(.user_passthrough_blocked("g"))
  expect_false(.user_passthrough_blocked(NA_character_))
})

# ---------------------------------------------------------------------------
# FIX 4 -- there was no exported back-transformation, and metafor's transf.iabt
# (the natural thing to reach for) silently maps the whole pool to 0, because
# metafor's ABT stores -log(1-a) where metaConvert stores +log(1-a).
# ---------------------------------------------------------------------------
test_that("reliability_backtransform round-trips the Bonett transform", {
  a <- c(0.78, 0.85, 0.91, 0.88)
  e <- es_from_cronbach_alpha(a, n_sample = 200, n_items = 10)
  expect_equal(reliability_backtransform(e$alpha)$reliability, a, tolerance = 1e-12)
})

test_that("reliability_backtransform swaps the CI bounds", {
  e <- es_from_cronbach_alpha(0.85, n_sample = 200, n_items = 10)
  bt <- reliability_backtransform(e$alpha, e$alpha_ci_lo, e$alpha_ci_up)
  # 1 - exp(x) is DECREASING: the transformed UPPER bound is the alpha LOWER bound
  expect_equal(bt$reliability_ci_lo, 1 - exp(e$alpha_ci_up), tolerance = 1e-12)
  expect_equal(bt$reliability_ci_up, 1 - exp(e$alpha_ci_lo), tolerance = 1e-12)
  expect_lt(bt$reliability_ci_lo, bt$reliability)
  expect_gt(bt$reliability_ci_up, bt$reliability)
})

test_that("reliability_backtransform consumes rma and predict.rma objects", {
  skip_if_not_installed("metafor")
  e <- es_from_cronbach_alpha(c(0.78, 0.85, 0.91, 0.88), n_sample = 200, n_items = 10)
  m <- metafor::rma(yi = e$alpha, sei = e$alpha_se, method = "REML")
  bt <- reliability_backtransform(m)
  expect_equal(bt$reliability, 1 - exp(as.numeric(m$beta[1])), tolerance = 1e-12)
  expect_equal(bt$reliability_ci_lo, 1 - exp(m$ci.ub), tolerance = 1e-12)
  expect_equal(bt$reliability_ci_up, 1 - exp(m$ci.lb), tolerance = 1e-12)
  # prediction interval present and wider than the CI on both sides
  expect_true(all(c("reliability_pi_lo", "reliability_pi_up") %in% names(bt)))
  expect_lt(bt$reliability_pi_lo, bt$reliability_ci_lo)
  expect_gt(bt$reliability_pi_up, bt$reliability_ci_up)
  expect_equal(reliability_backtransform(stats::predict(m))$reliability,
               bt$reliability, tolerance = 1e-12)
})

test_that("reliability_backtransform is NOT metafor::transf.iabt", {
  skip_if_not_installed("metafor")
  e <- es_from_cronbach_alpha(c(0.78, 0.85, 0.91, 0.88), n_sample = 200, n_items = 10)
  m <- metafor::rma(yi = e$alpha, sei = e$alpha_se, method = "REML")
  ours <- reliability_backtransform(m)$reliability
  theirs <- as.numeric(stats::predict(m, transf = metafor::transf.iabt)$pred)
  expect_gt(ours, 0.8)
  expect_equal(theirs, 0)          # documents the trap this helper exists to avoid
  expect_false(isTRUE(all.equal(ours, theirs)))
})

test_that("reliability_backtransform: raw scale passes through without swapping", {
  bt <- reliability_backtransform(0.85, 0.80, 0.90, method = "raw")
  expect_equal(bt$reliability, 0.85)
  expect_equal(bt$reliability_ci_lo, 0.80)
  expect_equal(bt$reliability_ci_up, 0.90)
})

test_that("reliability_backtransform does not clamp a negative reliability", {
  # transf.iabt floors at 0; a negative alpha is unusual but possible and the
  # user needs to see it
  expect_lt(reliability_backtransform(0.5)$reliability, 0)
  expect_equal(reliability_backtransform(0.5)$reliability, 1 - exp(0.5))
})

test_that("reliability_backtransform validates its arguments", {
  expect_error(reliability_backtransform(1, method = "hakstian_whalen"), "method")
  expect_error(reliability_backtransform(c(1, 2, 3), ci_lo = c(1, 2)), "ci_lo")
  expect_true(all(is.na(reliability_backtransform(c(-1, -2))$reliability_ci_lo)))
})

# --- object-shape hardening (found by adversarial review of the fix itself) ---

test_that("reliability_backtransform refuses a meta-regression instead of mangling it", {
  skip_if_not_installed("metafor")
  e <- es_from_cronbach_alpha(c(0.78, 0.85, 0.91, 0.88), n_sample = 200, n_items = 10)
  mods <- c(1, 1, 2, 2)
  # random-effects: used to abort inside `$<-.data.frame` (PI length k vs estimate length 1)
  mm <- metafor::rma(yi = e$alpha, sei = e$alpha_se, mods = ~ mods, method = "REML")
  expect_error(reliability_backtransform(mm), "meta-regression")
  # fixed-effect: no PI, so no length clash -- used to SILENTLY return the reference level
  mf <- metafor::rma(yi = e$alpha, sei = e$alpha_se, mods = ~ mods, method = "FE")
  expect_error(reliability_backtransform(mf), "meta-regression")
  # the route the error points at works and returns one row per fitted value
  expect_equal(nrow(reliability_backtransform(stats::predict(mm))), 4L)
})

test_that("reliability_backtransform reads a predict() table that went through as.data.frame()", {
  skip_if_not_installed("metafor")
  e <- es_from_cronbach_alpha(c(0.78, 0.85, 0.91, 0.88), n_sample = 200, n_items = 10)
  m <- metafor::rma(yi = e$alpha, sei = e$alpha_se, method = "REML")
  ref <- reliability_backtransform(stats::predict(m))
  # used to fall through to as.numeric(), FLATTENING the columns into six bogus
  # rows whose first value was coincidentally right
  got <- reliability_backtransform(as.data.frame(stats::predict(m)))
  expect_equal(nrow(got), 1L)
  expect_equal(got$reliability, ref$reliability, tolerance = 1e-12)
  expect_equal(got$reliability_ci_lo, ref$reliability_ci_lo, tolerance = 1e-12)
  expect_equal(got$reliability_pi_lo, ref$reliability_pi_lo, tolerance = 1e-12)
})

test_that("reliability_backtransform handles a list.rma carrying no CI (blup)", {
  skip_if_not_installed("metafor")
  e <- es_from_cronbach_alpha(c(0.78, 0.85, 0.91, 0.88), n_sample = 200, n_items = 10)
  m <- metafor::rma(yi = e$alpha, sei = e$alpha_se, method = "REML")
  # blup() has pi.lb/pi.ub but no ci.lb/ci.ub: the zero-length numeric used to
  # defeat the is.null() fallback and abort naming a ci_lo the caller never passed
  b <- reliability_backtransform(metafor::blup(m))
  expect_equal(nrow(b), 4L)
  expect_true(all(is.na(b$reliability_ci_lo)))
  expect_false(anyNA(b$reliability))
  expect_true(all(b$reliability_pi_lo < b$reliability_pi_up))
})
