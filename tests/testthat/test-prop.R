# =============================================================================
# Single-group proportions (es_from_PROP.R) — promoted from tests_save/checked/
# (P26, audit digest [25]).
#
# Adapted to the phase-1 fixes in R/es_from_PROP.R:
#   - P9: out-of-range proportions and n_cases > n_sample are now per-row
#     NA + warning instead of hard stop()s (the archived "ERROR PATHS"
#     section expected errors).
#   - P5: boundary (p = 0 or 1) raw/logit SEs now use the continuity-
#     corrected denominator n + 1 (metafor PR/PLO convention). The boundary
#     tests below only assert point estimates (unchanged); the boundary SE
#     oracle lives in test-EXTERNAL-PROP.R and test-prop-fixes.R.
# =============================================================================

# PROPORTION - SINGLE GROUP ----------------------------------------------------

test_that("es_from_prop_single_group computes raw proportion by default", {
  res <- es_from_prop_single_group(prop = 0.30, n_sample = 100)

  expect_equal(res$prop, 0.30, tolerance = 1e-10)
  expect_equal(res$prop_se, sqrt(0.30 * 0.70 / 100), tolerance = 1e-10)
  expect_equal(res$info_used, "prop_single_group")
})

test_that("es_from_prop_single_group with prop_to_es='logit'", {
  res <- es_from_prop_single_group(prop = 0.30, n_sample = 100, prop_to_es = "logit")

  expect_equal(res$prop, log(0.30 / 0.70), tolerance = 1e-10)
  expect_equal(res$prop_se, sqrt(1 / (100 * 0.30 * 0.70)), tolerance = 1e-10)
})

test_that("es_from_prop_single_group with prop_to_es='freeman_tukey'", {
  p <- 0.30; n <- 100; x <- p * n
  expected_ft <- 0.5 * (asin(sqrt(x / (n + 1))) + asin(sqrt((x + 1) / (n + 1))))
  expected_se <- 1 / sqrt(4 * n + 2)

  res <- es_from_prop_single_group(prop = p, n_sample = n, prop_to_es = "freeman_tukey")

  expect_equal(res$prop, expected_ft, tolerance = 1e-10)
  expect_equal(res$prop_se, expected_se, tolerance = 1e-10)
})

test_that("es_from_prop_single_group validates prop_to_es parameter", {
  expect_error(
    es_from_prop_single_group(prop = 0.30, n_sample = 100, prop_to_es = "invalid"),
    "not in tolerated values"
  )
})

test_that("es_from_prop_single_group handles edge case: prop = 0 with continuity correction", {
  res_raw <- es_from_prop_single_group(prop = 0, n_sample = 50, prop_to_es = "raw")
  expected_corrected <- 0.5 / 51
  expect_equal(res_raw$prop, expected_corrected, tolerance = 1e-10)

  res_logit <- es_from_prop_single_group(prop = 0, n_sample = 50, prop_to_es = "logit")
  expect_equal(res_logit$prop, log(expected_corrected / (1 - expected_corrected)), tolerance = 1e-10)

  # Freeman-Tukey handles zeros naturally (no continuity correction)
  res_ft <- es_from_prop_single_group(prop = 0, n_sample = 50, prop_to_es = "freeman_tukey")
  expect_equal(res_ft$prop, 0.5 * (asin(sqrt(0 / 51)) + asin(sqrt(1 / 51))), tolerance = 1e-10)
})

test_that("es_from_prop_single_group handles edge case: prop = 1 with continuity correction", {
  n <- 100
  expected_corrected <- (n + 0.5) / (n + 1)

  res_raw <- es_from_prop_single_group(prop = 1, n_sample = n, prop_to_es = "raw")
  expect_equal(res_raw$prop, expected_corrected, tolerance = 1e-10)

  res_logit <- es_from_prop_single_group(prop = 1, n_sample = n, prop_to_es = "logit")
  expect_equal(res_logit$prop, log(expected_corrected / (1 - expected_corrected)), tolerance = 1e-10)
})

test_that("es_from_prop_single_group handles vectorized input", {
  res <- es_from_prop_single_group(
    prop = c(0.30, 0.50, 0.75),
    n_sample = c(100, 200, 80)
  )
  expect_equal(nrow(res), 3)
  expect_equal(res$info_used, rep("prop_single_group", 3))
})

test_that("es_from_prop_single_group handles reversal", {
  res_fwd <- es_from_prop_single_group(prop = 0.30, n_sample = 100)
  res_rev <- es_from_prop_single_group(prop = 0.30, n_sample = 100, reverse_prop = TRUE)

  expect_equal(res_rev$prop, 1 - res_fwd$prop, tolerance = 1e-10)

  # Logit flips sign when reversed
  res_logit_fwd <- es_from_prop_single_group(prop = 0.30, n_sample = 100, prop_to_es = "logit")
  res_logit_rev <- es_from_prop_single_group(prop = 0.30, n_sample = 100, prop_to_es = "logit",
                                              reverse_prop = TRUE)
  expect_equal(res_logit_rev$prop, -res_logit_fwd$prop, tolerance = 1e-10)
})

test_that("es_from_prop_single_group raw CIs bounded [0, 1]", {
  res <- es_from_prop_single_group(prop = 0.05, n_sample = 20)
  expect_true(res$prop_ci_lo >= 0)
  expect_true(res$prop_ci_up <= 1)
})

# COUNT DATA ---------------------------------------------------------------

test_that("es_from_prop_single_group_counts delegates correctly", {
  res_counts <- es_from_prop_single_group_counts(n_cases = 30, n_sample = 100)
  res_prop <- es_from_prop_single_group(prop = 0.30, n_sample = 100)

  expect_equal(res_counts$prop, res_prop$prop, tolerance = 1e-10)
  expect_equal(res_counts$prop_se, res_prop$prop_se, tolerance = 1e-10)
  expect_equal(res_counts$info_used, "prop_single_group_counts")
})

test_that("es_from_prop_single_group_counts passes prop_to_es through", {
  res <- es_from_prop_single_group_counts(n_cases = 30, n_sample = 100, prop_to_es = "logit")
  expected <- es_from_prop_single_group(prop = 0.30, n_sample = 100, prop_to_es = "logit")

  expect_equal(res$prop, expected$prop, tolerance = 1e-10)
  expect_equal(res$prop_se, expected$prop_se, tolerance = 1e-10)
})

# CONVERT_DF INTEGRATION ---------------------------------------------------

test_that("convert_df with measure='prop' and default raw method", {
  dat <- data.frame(
    study_id = 1:3,
    prop = c(0.30, 0.50, 0.75),
    n_sample = c(100, 200, 80)
  )
  mc <- convert_df(dat, measure = "prop", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  expect_equal(s$es, dat$prop, tolerance = 1e-10)
  expect_equal(s$info_used, rep("prop_single_group", 3))
})

test_that("convert_df with measure='prop' and logit method", {
  dat <- data.frame(
    study_id = 1:3,
    prop = c(0.30, 0.50, 0.75),
    n_sample = c(100, 200, 80)
  )
  mc <- convert_df(dat, measure = "prop", prop_to_es = "logit",
                   verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  expected_logit <- log(dat$prop / (1 - dat$prop))
  expect_equal(s$es, expected_logit, tolerance = 1e-10)
})

test_that("convert_df with measure='prop' and freeman_tukey method", {
  dat <- data.frame(
    study_id = 1:3,
    prop = c(0.30, 0.50, 0.75),
    n_sample = c(100, 200, 80)
  )
  mc <- convert_df(dat, measure = "prop", prop_to_es = "freeman_tukey",
                   verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  x <- dat$prop * dat$n_sample
  expected_ft <- 0.5 * (asin(sqrt(x / (dat$n_sample + 1))) +
                         asin(sqrt((x + 1) / (dat$n_sample + 1))))
  expect_equal(s$es, expected_ft, tolerance = 1e-10)
})

test_that("convert_df with measure='prop' mixed input (prop + n_cases)", {
  dat <- data.frame(
    study_id = 1:4,
    prop = c(0.30, NA, 0.50, NA),
    n_cases = c(NA, 15, NA, 50),
    n_sample = c(100, 200, 150, 150)
  )
  mc <- convert_df(dat, measure = "prop", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  expect_equal(s$es, c(0.30, 15/200, 0.50, 50/150), tolerance = 1e-10)
  expect_equal(s$info_used,
               c("prop_single_group", "prop_single_group_counts",
                 "prop_single_group", "prop_single_group_counts"))
})

test_that("convert_df with measure='prop' hierarchy (prop prioritized over n_cases)", {
  dat <- data.frame(
    study_id = 1,
    prop = 0.35,
    n_cases = 30,
    n_sample = 100
  )
  mc <- convert_df(dat, measure = "prop", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  expect_equal(s$info_used, "prop_single_group")
  expect_equal(s$es, 0.35, tolerance = 1e-10)
})

test_that("prop data_extraction_sheet works", {
  dat <- data_extraction_sheet(measure = "prop", extension = "data.frame", verbose = FALSE)
  expect_true("prop" %in% colnames(dat))
  expect_true("n_cases" %in% colnames(dat))
  expect_true("n_sample" %in% colnames(dat))
})

# INVALID-INPUT PATHS (P9: per-row NA + warning, no longer hard stops) --------

test_that("es_from_prop_single_group warns and NAs rows with prop < 0 (was a stop)", {
  expect_warning(
    res <- es_from_prop_single_group(prop = c(-0.1, 0.30), n_sample = c(100, 100)),
    "outside \\[0, 1\\]"
  )
  expect_true(is.na(res$prop[1]))
  expect_true(is.na(res$prop_se[1]))
  # the valid row is unaffected
  expect_equal(res$prop[2], 0.30, tolerance = 1e-10)
  expect_equal(res$prop_se[2], sqrt(0.30 * 0.70 / 100), tolerance = 1e-10)
})

test_that("es_from_prop_single_group warns and NAs rows with prop > 1 (was a stop)", {
  expect_warning(
    res <- es_from_prop_single_group(prop = 1.5, n_sample = 100),
    "outside \\[0, 1\\]"
  )
  expect_true(is.na(res$prop))
  expect_true(is.na(res$prop_se))
})

test_that("es_from_prop_single_group_counts warns and NAs pairs with n_cases > n_sample (was a stop)", {
  expect_warning(
    res <- es_from_prop_single_group_counts(n_cases = c(110, 30), n_sample = c(100, 100)),
    "n_cases > n_sample"
  )
  expect_true(is.na(res$prop[1]))
  expect_true(is.na(res$prop_se[1]))
  expect_equal(res$prop[2], 0.30, tolerance = 1e-10)
})

# FREEMAN-TUKEY REVERSAL ----------------------------------------------------

test_that("Freeman-Tukey reversal is symmetric", {
  res_fwd <- es_from_prop_single_group(prop = 0.30, n_sample = 100,
                                        prop_to_es = "freeman_tukey")
  res_rev <- es_from_prop_single_group(prop = 0.30, n_sample = 100,
                                        prop_to_es = "freeman_tukey",
                                        reverse_prop = TRUE)
  # FT(0.30) + FT(0.70) should roughly sum to pi/2 (due to arcsin symmetry)
  # More precisely: FT(p) + FT(1-p) ~ pi/2 for large n
  res_70 <- es_from_prop_single_group(prop = 0.70, n_sample = 100,
                                       prop_to_es = "freeman_tukey")
  expect_equal(res_rev$prop, res_70$prop, tolerance = 1e-10)
})

# EDGE CASES ----------------------------------------------------------------

test_that("all-NA input returns all-NA output", {
  res <- es_from_prop_single_group(prop = c(NA, NA), n_sample = c(NA, NA))
  expect_true(all(is.na(res$prop)))
  expect_true(all(is.na(res$prop_se)))
})

test_that("n_sample = 1 produces valid output", {
  res <- es_from_prop_single_group(prop = 0.5, n_sample = 1)
  expect_true(is.finite(res$prop))
  expect_true(is.finite(res$prop_se))
})

test_that("extreme proportions with small n (Freeman-Tukey)", {
  res <- es_from_prop_single_group(prop = 0.001, n_sample = 5,
                                    prop_to_es = "freeman_tukey")
  expect_true(is.finite(res$prop))
  expect_true(is.finite(res$prop_se))
  expect_true(res$prop >= 0)
})
