# =============================================================================
# Audit fixes for the single-group proportion paths (R/es_from_PROP.R).
#
#   P9 — the two hard stop()s ("Proportions must be between 0 and 1",
#        "Number of cases cannot exceed sample size") aborted the ENTIRE
#        convert_df() pipeline under correct_inputs = FALSE, because
#        convert_df() calls every es_from_* function for every measure with
#        no tryCatch — one stray bad prop killed a whole measure = "d"
#        analysis too. The package convention is per-row NA + warning.
#
#   P5 — at p == 0/1 the continuity correction computed p_c = (x+0.5)/(n+1)
#        but the raw/logit SEs divided by the UNcorrected n, deviating from
#        metafor's PR/PLO convention (the roxygen's declared reference) by
#        exactly sqrt((n+1)/n). Boundary SEs must bit-match metafor::escalc.
#
# Oracle rule: expected values come from metafor::escalc computed at test
# time or from hand-derived binomial formulas written independently below —
# never from the package function under test.
# =============================================================================

# -----------------------------------------------------------------------------
# P9 — direct calls: per-row NA + warning instead of a hard stop
# -----------------------------------------------------------------------------

test_that("P9: out-of-range prop is set to NA per row with a warning, not an error", {
  expect_warning(
    res <- es_from_prop_single_group(prop = c(1.5, 0.3), n_sample = c(100, 100)),
    regexp = "outside \\[0, 1\\]"
  )
  expect_true(is.na(res$prop[1]))
  expect_true(is.na(res$prop_se[1]))
  # the valid row must be unaffected: hand-derived binomial SE sqrt(p(1-p)/n)
  expect_equal(res$prop[2], 0.3, tolerance = 1e-10)
  expect_equal(res$prop_se[2], sqrt(0.3 * 0.7 / 100), tolerance = 1e-10)
})

test_that("P9: negative prop is also NA'd per row, not an error", {
  expect_warning(
    res <- es_from_prop_single_group(prop = c(-0.2, 0.4), n_sample = c(50, 50)),
    regexp = "outside \\[0, 1\\]"
  )
  expect_true(is.na(res$prop[1]))
  expect_true(is.na(res$prop_se[1]))
  expect_equal(res$prop[2], 0.4, tolerance = 1e-10)
  expect_equal(res$prop_se[2], sqrt(0.4 * 0.6 / 50), tolerance = 1e-10)
})

test_that("P9: n_cases > n_sample is NA'd per pair with a warning, not an error", {
  expect_warning(
    res <- es_from_prop_single_group_counts(
      n_cases = c(110, 25), n_sample = c(100, 100)
    ),
    regexp = "n_cases > n_sample"
  )
  expect_true(is.na(res$prop[1]))
  expect_true(is.na(res$prop_se[1]))
  expect_equal(res$prop[2], 0.25, tolerance = 1e-10)
  expect_equal(res$prop_se[2], sqrt(0.25 * 0.75 / 100), tolerance = 1e-10)
})

test_that("P9: the unknown prop_to_es argument check remains a hard stop", {
  # a wrong prop_to_es string is a programming error, not bad data
  expect_error(
    es_from_prop_single_group(prop = 0.3, n_sample = 100, prop_to_es = "arcsine"),
    regexp = "not in tolerated values"
  )
})

# -----------------------------------------------------------------------------
# P9 — pipeline: convert_df(correct_inputs = FALSE) must not abort
# -----------------------------------------------------------------------------

test_that("P9: convert_df(measure='prop', correct_inputs=FALSE) survives an out-of-range prop", {
  dat <- data.frame(study_id = 1:2, prop = c(1.5, 0.3), n_sample = 100)
  expect_no_error(
    res <- suppressWarnings(
      convert_df(dat, measure = "prop", verbose = FALSE, correct_inputs = FALSE)
    )
  )
  smy <- suppressWarnings(summary(res))
  expect_true(is.na(smy$es_crude[smy$row_id == 1]))
  # summary() rounds for display, hence the loose tolerance on the valid row
  expect_equal(as.numeric(smy$es_crude[smy$row_id == 2]), 0.3, tolerance = 1e-3)
})

test_that("P9: the es_from_PROP.R stop no longer fires in the pipeline for n_cases > n_sample", {
  dat <- data.frame(study_id = 1:2, n_cases = c(110, 25), n_sample = 100)
  res <- tryCatch(
    suppressWarnings(
      convert_df(dat, measure = "prop", verbose = FALSE, correct_inputs = FALSE)
    ),
    error = function(e) conditionMessage(e)
  )
  # Pre-fix this call aborted with the es_from_PROP.R hard stop "Number of
  # cases cannot exceed sample size". That stop is gone. NOTE: a residual,
  # SEPARATE hard stop remains in R/es_from_PHI_CHISQ.R — metafor::conv.2x2
  # validates n_cases > n_sample even when phi/chisq is NA ("One or more
  # marginal counts are larger than the sample sizes.") — which is outside
  # the scope of this fix (file owned elsewhere). The conditional block below
  # verifies the full end-to-end behaviour automatically once that residual
  # stop is also converted to the per-row NA convention.
  expect_false(identical(res, "Number of cases cannot exceed sample size"))
  if (!is.character(res)) {
    smy <- suppressWarnings(summary(res))
    expect_true(is.na(smy$es_crude[smy$row_id == 1]))
    expect_equal(as.numeric(smy$es_crude[smy$row_id == 2]), 0.25, tolerance = 1e-3)
  }
})

test_that("P9: a stray invalid prop no longer kills a measure='d' analysis (correct_inputs=FALSE)", {
  dat <- data.frame(
    study_id = 1, mean_exp = 12, mean_nexp = 10,
    mean_sd_exp = 3, mean_sd_nexp = 3, n_exp = 40, n_nexp = 40,
    prop = 1.5
  )
  expect_no_error(
    res <- suppressWarnings(
      convert_df(dat, measure = "d", verbose = FALSE, correct_inputs = FALSE)
    )
  )
  smy <- suppressWarnings(summary(res))
  # hand-derived Cohen's d: pooled SD = sqrt((39*9 + 39*9)/78) = 3, d = (12-10)/3
  expect_equal(as.numeric(smy$es_crude[1]), 2 / 3, tolerance = 1e-3)
})

# -----------------------------------------------------------------------------
# P5 — boundary continuity correction must use the corrected denominator n + 1
# -----------------------------------------------------------------------------

test_that("P5: boundary (p = 0 or 1) raw/logit SEs bit-match metafor PR/PLO", {
  skip_if_not_installed("metafor")
  cases <- list(c(x = 0, n = 20), c(x = 20, n = 20), c(x = 0, n = 50))
  for (cs in cases) {
    x <- unname(cs["x"])
    n <- unname(cs["n"])
    mf_pr <- metafor::escalc(measure = "PR", xi = x, ni = n)
    mf_plo <- metafor::escalc(measure = "PLO", xi = x, ni = n)
    res_raw <- es_from_prop_single_group(
      prop = x / n, n_sample = n, prop_to_es = "raw"
    )
    res_logit <- es_from_prop_single_group(
      prop = x / n, n_sample = n, prop_to_es = "logit"
    )
    # point estimates were already bit-exact pre-fix: regression assertions
    expect_equal(res_raw$prop, as.numeric(mf_pr$yi), tolerance = 1e-10)
    expect_equal(res_logit$prop, as.numeric(mf_plo$yi), tolerance = 1e-10)
    # SEs: the fix — the corrected denominator n + 1 flows into the SE, so
    # sqrt(p_c(1-p_c)/(n+1)) [raw] and sqrt(1/((n+1) p_c (1-p_c))) [logit]
    expect_equal(res_raw$prop_se, sqrt(as.numeric(mf_pr$vi)), tolerance = 1e-10)
    expect_equal(res_logit$prop_se, sqrt(as.numeric(mf_plo$vi)), tolerance = 1e-10)
  }
})

test_that("P5: Freeman-Tukey boundary is untouched (already bit-exact vs metafor PFT)", {
  skip_if_not_installed("metafor")
  for (cs in list(c(0, 20), c(20, 20))) {
    x <- cs[1]
    n <- cs[2]
    mf_pft <- metafor::escalc(measure = "PFT", xi = x, ni = n)
    res_ft <- es_from_prop_single_group(
      prop = x / n, n_sample = n, prop_to_es = "freeman_tukey"
    )
    expect_equal(res_ft$prop, as.numeric(mf_pft$yi), tolerance = 1e-10)
    expect_equal(res_ft$prop_se, sqrt(as.numeric(mf_pft$vi)), tolerance = 1e-10)
  }
})

# -----------------------------------------------------------------------------
# Interior regression — nothing may move for non-boundary proportions
# -----------------------------------------------------------------------------

test_that("interior regression: p = 0.3, n = 100 bit-matches metafor PR/PLO/PFT (yi and SE)", {
  skip_if_not_installed("metafor")
  mf_pr <- metafor::escalc(measure = "PR", xi = 30, ni = 100)
  mf_plo <- metafor::escalc(measure = "PLO", xi = 30, ni = 100)
  mf_pft <- metafor::escalc(measure = "PFT", xi = 30, ni = 100)
  res_raw <- es_from_prop_single_group(prop = 0.3, n_sample = 100, prop_to_es = "raw")
  res_logit <- es_from_prop_single_group(prop = 0.3, n_sample = 100, prop_to_es = "logit")
  res_ft <- es_from_prop_single_group(prop = 0.3, n_sample = 100, prop_to_es = "freeman_tukey")
  expect_equal(res_raw$prop, as.numeric(mf_pr$yi), tolerance = 1e-10)
  expect_equal(res_raw$prop_se, sqrt(as.numeric(mf_pr$vi)), tolerance = 1e-10)
  expect_equal(res_logit$prop, as.numeric(mf_plo$yi), tolerance = 1e-10)
  expect_equal(res_logit$prop_se, sqrt(as.numeric(mf_plo$vi)), tolerance = 1e-10)
  expect_equal(res_ft$prop, as.numeric(mf_pft$yi), tolerance = 1e-10)
  expect_equal(res_ft$prop_se, sqrt(as.numeric(mf_pft$vi)), tolerance = 1e-10)
})

test_that("continuity correction still does not touch user-entered nonzero proportions", {
  # p = 0.05 with n = 20 (x = 1) is interior: no 0.5-correction may be applied
  res <- es_from_prop_single_group(prop = 0.05, n_sample = 20, prop_to_es = "raw")
  expect_identical(res$prop, 0.05)
  expect_equal(res$prop_se, sqrt(0.05 * 0.95 / 20), tolerance = 1e-14)
})
