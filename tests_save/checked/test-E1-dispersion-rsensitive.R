library(testthat)

# ==============================================================================
# Fix (b): E1 dispersion signal is now a k-robust max-abs-deviation-from-median
# (.dispersion_stat), replacing the raw sample SD.
# ==============================================================================

test_that(".dispersion_stat returns NA for fewer than two finite values", {
  expect_true(is.na(metaConvert:::.dispersion_stat(numeric(0))))
  expect_true(is.na(metaConvert:::.dispersion_stat(5)))
  expect_true(is.na(metaConvert:::.dispersion_stat(c(NA, 2))))
})

test_that(".dispersion_stat equals |a - b| / 2 for exactly two estimates", {
  expect_equal(metaConvert:::.dispersion_stat(c(0.5, 1.5)), 0.5)
  expect_equal(metaConvert:::.dispersion_stat(c(-1, 1)), 1)
})

test_that(".dispersion_stat is k-invariant to agreeing peers (unlike sd)", {
  # A lone outlier at distance 2 from a cluster: MaxAD stays 2 as agreeing
  # methods are added, whereas sd shrinks ~1/sqrt(k) and would fall under the
  # E1 threshold in exactly the data-rich rows where the outlier matters.
  expect_equal(metaConvert:::.dispersion_stat(c(1, 1, 3)), 2)
  expect_equal(metaConvert:::.dispersion_stat(c(1, 1, 1, 1, 3)), 2)
  expect_lt(sd(c(1, 1, 1, 1, 3)), 2)          # sd shrinks
  expect_lt(sd(c(1, 1, 1, 1, 3)), sd(c(1, 1, 3)))
})

test_that(".dispersion_stat ignores non-finite entries", {
  expect_equal(metaConvert:::.dispersion_stat(c(NA, 1, 3, Inf)), 1)  # median(1,3)=2
})

test_that("E1 message now reports max|ES - median|, not SD", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = 2, dispersion = 0.8, overlap = 1, diff_mm = 0.1,
    opts = opts, measure = "g"
  )
  expect_true(any(grepl("High dispersion across estimates: max\\|ES - median\\|", flags[[1]])))
})


# ==============================================================================
# Fix (a): SE/CI-driven DISCORDANT flags (A6, E2, E2b) are tagged r-sensitive
# when r_pre_post was defaulted. Point-estimate checks (E1, E3) are NOT tagged.
# ==============================================================================

test_that(".r_sensitive_note tags TRUE, is empty for FALSE, and has no ';'", {
  note <- metaConvert:::.r_sensitive_note(TRUE)
  expect_true(grepl("r-sensitive", note))
  expect_false(grepl(";", note))                 # must not create an untagged token
  expect_equal(metaConvert:::.r_sensitive_note(FALSE), "")
  expect_equal(metaConvert:::.r_sensitive_note(NA), "")
})

test_that("[E2b] low CI overlap is tagged r-sensitive only when r was defaulted", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = c(2, 2),
    dispersion    = c(0.1, 0.1),   # below E1 threshold -> isolate E2b
    overlap       = c(0.5, 0.5),   # 0.5 < overlap_min 0.85 -> E2b fires
    diff_mm       = c(0.1, 0.1),   # below E3 threshold
    opts = opts, measure = "g",
    r_defaulted = c(TRUE, FALSE)
  )
  expect_true(any(grepl("Low CI overlap", flags[[1]])))
  expect_true(any(grepl("r-sensitive", flags[[1]])))
  expect_true(any(grepl("Low CI overlap", flags[[2]])))
  expect_false(any(grepl("r-sensitive", flags[[2]])))
})

test_that("[E1] high dispersion is NOT tagged r-sensitive (point-estimate check)", {
  opts <- metaConvert:::.default_flag_options()
  flags <- metaConvert:::.flag_internal_consistency(
    n_estimations = 2, dispersion = 0.8, overlap = 1, diff_mm = 0.1,
    opts = opts, measure = "g", r_defaulted = TRUE
  )
  expect_true(any(grepl("High dispersion", flags[[1]])))
  expect_false(any(grepl("r-sensitive", flags[[1]])))
})

test_that("[A6] CI-width inconsistency is tagged r-sensitive only when r was defaulted", {
  flags <- metaConvert:::.flag_numeric_integrity(
    es       = c(0.5, 0.5),
    se       = c(0.1, 0.1),
    ci_lo    = c(0.4, 0.4),   # width 0.2 vs expected ~0.40 -> A6 fires
    ci_up    = c(0.6, 0.6),
    info_used = c("means_sd", "means_sd"),
    measure = "g",
    n_total = c(100, 100),
    r_defaulted = c(TRUE, FALSE)
  )
  expect_true(any(grepl("CI width inconsistent with SE", flags[[1]])))
  expect_true(any(grepl("r-sensitive", flags[[1]])))
  expect_true(any(grepl("CI width inconsistent with SE", flags[[2]])))
  expect_false(any(grepl("r-sensitive", flags[[2]])))
})

test_that("convert_df exposes a per-row r_defaulted attribute for pre/post data", {
  dat <- data.frame(
    n_exp = 40, n_nexp = 40,
    mean_pre_exp = 10, mean_pre_sd_exp = 3, mean_exp = 8, mean_sd_exp = 3,
    mean_pre_nexp = 10, mean_pre_sd_nexp = 3, mean_nexp = 9, mean_sd_nexp = 3
  )
  res <- suppressWarnings(suppressMessages(
    metaConvert::convert_df(dat, measure = "g", verbose = FALSE)
  ))
  rd <- attr(res, "r_defaulted")
  expect_type(rd, "logical")
  expect_true(any(rd))   # this row used the default r_pre_post
})
