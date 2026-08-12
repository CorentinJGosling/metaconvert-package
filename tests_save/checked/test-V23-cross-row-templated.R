library(testthat)

# ==============================================================================
# V23: Cross-row identical summary-statistic block (templated / non-independent
# data). Detects rows from DIFFERENT studies that share a byte-for-byte
# identical block of continuous summary statistics (means + SDs of both arms).
# This is the cross-study analogue of the duplicate-study check (Category H):
# H keys on identical study_id, V23 keys on identical *values* across different
# study_ids — the paper-mill / cloned-table / wrong-source signature.
#
# Motivating real case: Huang 2017 and Xu 2017 (different authors, sites, N,
# enrolment windows and journals) report byte-identical baseline CDR in both
# arms (2.20 +/- 0.83 / 2.19 +/- 0.81) — implausible for independent cohorts.
# ==============================================================================

V23 <- function(dat, ...) {
  metaConvert:::.validate_input_data(dat, verbose = FALSE, ...)$issues
}

# ------------------------------------------------------------------------------
# Positive detection
# ------------------------------------------------------------------------------

test_that("[V23] fires on identical baseline block across different study_ids (Huang/Xu)", {
  dat <- data.frame(
    study_id        = c("huang_2017", "xu_2017"),
    n_exp           = c(40, 36),
    n_nexp          = c(40, 36),
    mean_pre_exp    = c(2.20, 2.20),
    mean_pre_sd_exp = c(0.83, 0.83),
    mean_pre_nexp   = c(2.19, 2.19),
    mean_pre_sd_nexp = c(0.81, 0.81),
    stringsAsFactors = FALSE
  )
  iss <- V23(dat)
  expect_true(grepl("Identical baseline summary statistics", iss[1]))
  expect_true(grepl("Identical baseline summary statistics", iss[2]))
  # Each row names the OTHER row as "study_id (row N)"
  expect_true(grepl("xu_2017 (row 2)", iss[1], fixed = TRUE))
  expect_true(grepl("huang_2017 (row 1)", iss[2], fixed = TRUE))
  # [INFO] severity (preserved data, judgment-dependent)
  expect_true(grepl("^\\[INFO\\]", iss[1]))
})

test_that("[V23] fires on identical endpoint block", {
  dat <- data.frame(
    study_id     = c("p", "q"),
    n_exp        = c(50, 44),
    n_nexp       = c(50, 44),
    mean_exp     = c(12.34, 12.34),
    mean_sd_exp  = c(3.21, 3.21),
    mean_nexp    = c(9.87, 9.87),
    mean_sd_nexp = c(2.50, 2.50),
    stringsAsFactors = FALSE
  )
  iss <- V23(dat)
  expect_true(grepl("Identical endpoint summary statistics", iss[1]))
  expect_true(grepl("Identical endpoint summary statistics", iss[2]))
})

test_that("[V23] reports the shared values in the message", {
  dat <- data.frame(
    study_id        = c("a", "b"),
    mean_pre_exp    = c(2.20, 2.20),
    mean_pre_sd_exp = c(0.83, 0.83),
    mean_pre_nexp   = c(2.19, 2.19),
    mean_pre_sd_nexp = c(0.81, 0.81),
    stringsAsFactors = FALSE
  )
  iss <- V23(dat)
  expect_true(grepl("mean_pre_exp", iss[1]))
  expect_true(grepl("2.2", iss[1]))
})

test_that("[V23] fires across a duplicate group of 3 studies, naming all others", {
  dat <- data.frame(
    study_id        = c("a", "b", "c"),
    mean_pre_exp    = c(2.20, 2.20, 2.20),
    mean_pre_sd_exp = c(0.83, 0.83, 0.83),
    mean_pre_nexp   = c(2.19, 2.19, 2.19),
    mean_pre_sd_nexp = c(0.81, 0.81, 0.81),
    stringsAsFactors = FALSE
  )
  iss <- V23(dat)
  expect_equal(sum(grepl("Identical baseline", iss)), 3)
  # row 1 (study a) names studies b and c with their rows
  expect_true(grepl("b (row 2)", iss[1], fixed = TRUE))
  expect_true(grepl("c (row 3)", iss[1], fixed = TRUE))
})

# ------------------------------------------------------------------------------
# False-positive guards
# ------------------------------------------------------------------------------

test_that("[V23] does NOT fire when only one arm matches (shared-control multi-arm)", {
  # Drug-A-vs-placebo and Drug-B-vs-placebo: control arm repeats by design,
  # treatment arms differ. Requiring the FULL block to match protects this.
  dat <- data.frame(
    study_id     = c("trialA_arm1", "trialB_arm1"),
    mean_exp     = c(5.12, 8.44),   # treatment arms differ
    mean_sd_exp  = c(1.11, 1.22),
    mean_nexp    = c(2.19, 2.19),   # shared control arm
    mean_sd_nexp = c(0.81, 0.81),
    stringsAsFactors = FALSE
  )
  expect_false(any(grepl("Identical endpoint", V23(dat))))
})

test_that("[V23] does NOT fire when rows share the same study_id (Category H's job)", {
  # A 3-arm trial entered as two rows with identical pooled stats: same study_id
  # => skip (the duplicate-study check owns this; values legitimately repeat).
  dat <- data.frame(
    study_id     = c("trialX", "trialX"),
    mean_exp     = c(5.12, 5.12),
    mean_sd_exp  = c(1.11, 1.11),
    mean_nexp    = c(2.19, 2.19),
    mean_sd_nexp = c(0.81, 0.81),
    stringsAsFactors = FALSE
  )
  expect_false(any(grepl("Identical endpoint", V23(dat))))
})

test_that("[V23] does NOT fire on low-entropy all-integer blocks", {
  # Round integers collide too easily to be diagnostic; the non-integer gate
  # (templated_min_match) suppresses them.
  dat <- data.frame(
    study_id     = c("a", "b"),
    mean_exp     = c(10, 10),
    mean_sd_exp  = c(2, 2),
    mean_nexp    = c(8, 8),
    mean_sd_nexp = c(2, 2),
    stringsAsFactors = FALSE
  )
  expect_false(any(grepl("Identical endpoint", V23(dat))))
})

test_that("[V23] does NOT fire on distinct (non-identical) blocks", {
  dat <- data.frame(
    study_id     = c("a", "b"),
    mean_exp     = c(12.34, 12.35),  # differ in last decimal
    mean_sd_exp  = c(3.21, 3.21),
    mean_nexp    = c(9.87, 9.87),
    mean_sd_nexp = c(2.50, 2.50),
    stringsAsFactors = FALSE
  )
  expect_false(any(grepl("Identical endpoint", V23(dat))))
})

test_that("[V23] does NOT fire when the block is incomplete (any NA)", {
  dat <- data.frame(
    study_id     = c("a", "b"),
    mean_exp     = c(12.34, 12.34),
    mean_sd_exp  = c(3.21, 3.21),
    mean_nexp    = c(9.87, 9.87),
    mean_sd_nexp = c(NA, NA),  # block not fully present
    stringsAsFactors = FALSE
  )
  expect_false(any(grepl("Identical endpoint", V23(dat))))
})

# ------------------------------------------------------------------------------
# Configurability
# ------------------------------------------------------------------------------

test_that("[V23] respects the templated_min_match threshold", {
  # Only 2 of 4 values are non-integer and matching (mean_exp, mean_sd_exp).
  dat <- data.frame(
    study_id     = c("a", "b"),
    mean_exp     = c(12.34, 12.34),
    mean_sd_exp  = c(3.21, 3.21),
    mean_nexp    = c(9, 9),       # integer
    mean_sd_nexp = c(2, 2),       # integer
    stringsAsFactors = FALSE
  )
  expect_false(any(grepl("Identical endpoint", V23(dat, templated_min_match = 3))))
  expect_true(any(grepl("Identical endpoint", V23(dat, templated_min_match = 2))))
})

test_that("[V23] is disabled by enable_cross_row = FALSE", {
  dat <- data.frame(
    study_id        = c("a", "b"),
    mean_pre_exp    = c(2.20, 2.20),
    mean_pre_sd_exp = c(0.83, 0.83),
    mean_pre_nexp   = c(2.19, 2.19),
    mean_pre_sd_nexp = c(0.81, 0.81),
    stringsAsFactors = FALSE
  )
  expect_false(any(grepl("Identical baseline", V23(dat, enable_cross_row = FALSE))))
})

test_that("[V23] fires across different study_ids even when study_id column is absent", {
  # No study_id -> cannot rule out same-study, but identical multi-decimal block
  # still warrants a look. Message still emitted.
  dat <- data.frame(
    mean_pre_exp    = c(2.20, 2.20),
    mean_pre_sd_exp = c(0.83, 0.83),
    mean_pre_nexp   = c(2.19, 2.19),
    mean_pre_sd_nexp = c(0.81, 0.81),
    stringsAsFactors = FALSE
  )
  iss <- V23(dat)
  expect_true(any(grepl("Identical baseline", iss)))
  # No study_id -> reference falls back to bare "row N"
  expect_true(grepl("row 2", iss[1], fixed = TRUE))
  expect_false(grepl("(row 2)", iss[1], fixed = TRUE))
})

# ------------------------------------------------------------------------------
# Integration: V23 preserves data and surfaces through convert_df() -> summary()
# ------------------------------------------------------------------------------

test_that("[V23] preserves input data (warn only, no NA)", {
  dat <- data.frame(
    study_id        = c("a", "b"),
    mean_pre_exp    = c(2.20, 2.20),
    mean_pre_sd_exp = c(0.83, 0.83),
    mean_pre_nexp   = c(2.19, 2.19),
    mean_pre_sd_nexp = c(0.81, 0.81),
    stringsAsFactors = FALSE
  )
  result <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_equal(result$data$mean_pre_exp, c(2.20, 2.20))
  expect_equal(result$data$mean_pre_sd_nexp, c(0.81, 0.81))
})

test_that("[V23] surfaces in summary() flags_crude (crude-scope, not duplicated in adjusted)", {
  dat <- data.frame(
    study_id     = c("huang_2017", "xu_2017"),
    n_exp        = c(40, 36),
    n_nexp       = c(40, 36),
    mean_exp     = c(12.34, 12.34),
    mean_sd_exp  = c(3.21, 3.21),
    mean_nexp    = c(9.87, 9.87),
    mean_sd_nexp = c(2.50, 2.50),
    stringsAsFactors = FALSE
  )
  res <- convert_df(dat, measure = "g", verbose = FALSE)
  s <- summary(res, flags = TRUE)
  fc <- if ("flags_crude" %in% names(s)) s$flags_crude else s$flags
  expect_true(any(grepl("Identical endpoint summary statistics", fc)))
  if ("flags_adjusted" %in% names(s)) {
    expect_false(any(grepl("Identical endpoint summary statistics", s$flags_adjusted)))
  }
})
