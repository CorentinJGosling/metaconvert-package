# ==============================================================================
# Pin the tiered column layout introduced in metaConvert 1.1.0.
#
# `summary.metaConvert()` returns columns in named tiers so the eye walks
# the audit trail (provenance → quality → summary diagnostics) before
# landing on the chosen ES, then optionally drills into min/max details:
#
#   Identity   row_id, study_id, author, year, outcome (+ predictor,
#              info_expected, adjusted_input when present)
#
#   For each suffix (_crude, _adjusted, or "" when split_adjusted = FALSE):
#     Provenance  all_info, info_measure, es_guidance, es_selected
#     Quality     flags
#     Summary     es_summary, es_consistency, n_estimations, dispersion_es
#     Primary     es, se, es_ci_lo, es_ci_up, info_used, measure
#     Min/max     overlap_min_max, diff_min_max, min_info, max_info,
#                 min_es_value, min_es_se, min_es_ci_lo, min_es_ci_up,
#                 max_es_value, max_es_se, max_es_ci_lo, max_es_ci_up
#
#   Raw passthrough  appended at the end iff include_raw = TRUE.
#
# These tests fail loudly if anyone reorders or drops a tier column.
# ==============================================================================

library(testthat)
library(metaConvert)

# Helper: extract a vector of integer positions for a set of column names,
# returning -1 for any missing column so the test fails informatively.
.col_pos <- function(df, names) {
  out <- match(names, colnames(df))
  out[is.na(out)] <- -1L
  out
}

# Helper: assert that a sequence of named columns is contiguous and in order.
.expect_contiguous <- function(df, names, label) {
  pos <- .col_pos(df, names)
  expect_false(any(pos < 0),
               info = sprintf("%s — missing columns: %s", label,
                              paste(names[pos < 0], collapse = ", ")))
  expect_equal(diff(pos), rep(1L, length(pos) - 1L),
               info = sprintf("%s — columns not contiguous (positions: %s)",
                              label, paste(pos, collapse = ",")))
}

test_that("split_adjusted = TRUE produces tiered layout (default include_raw)", {
  obj <- convert_df(df.short, measure = "g", verbose = FALSE)
  res <- summary(obj, flags = TRUE, guidance = TRUE)

  # Identity tier appears first, in canonical order
  identity_present <- intersect(
    c("row_id", "study_id", "author", "year", "outcome",
      "predictor", "info_expected", "adjusted_input"),
    colnames(res)
  )
  .expect_contiguous(res, identity_present, "Identity tier")
  expect_equal(which(colnames(res) == identity_present[1]), 1L)

  # Per-suffix tier sequence: provenance → quality → summary → primary → minmax
  for (sfx in c("_crude", "_adjusted")) {
    .expect_contiguous(res, paste0(c("all_info", "info_measure",
                                     "es_guidance", "es_selected"), sfx),
                       paste0("Provenance", sfx))
    .expect_contiguous(res, paste0("flags", sfx), paste0("Quality", sfx))
    .expect_contiguous(res, paste0(c("es_summary", "es_consistency",
                                     "n_estimations", "dispersion_es"), sfx),
                       paste0("Summary", sfx))
    .expect_contiguous(res, paste0(c("es", "se", "es_ci_lo", "es_ci_up",
                                     "info_used", "measure"), sfx),
                       paste0("Primary", sfx))
    .expect_contiguous(res, paste0(c("overlap_min_max", "diff_min_max",
                                     "min_info", "max_info",
                                     "min_es_value", "min_es_se",
                                     "min_es_ci_lo", "min_es_ci_up",
                                     "max_es_value", "max_es_se",
                                     "max_es_ci_lo", "max_es_ci_up"), sfx),
                       paste0("Min/max", sfx))
  }

  # Provenance comes before Primary which comes before Min/max (per suffix)
  for (sfx in c("_crude", "_adjusted")) {
    expect_lt(which(colnames(res) == paste0("all_info", sfx)),
              which(colnames(res) == paste0("es", sfx)))
    expect_lt(which(colnames(res) == paste0("es", sfx)),
              which(colnames(res) == paste0("min_es_value", sfx)))
  }

  # _crude block comes before _adjusted block
  expect_lt(which(colnames(res) == "es_crude"),
            which(colnames(res) == "es_adjusted"))
})

test_that("include_raw = FALSE drops raw passthrough but keeps every tier column", {
  obj <- convert_df(df.short, measure = "g", verbose = FALSE)
  full <- summary(obj, flags = TRUE, guidance = TRUE, include_raw = TRUE)
  lean <- summary(obj, flags = TRUE, guidance = TRUE, include_raw = FALSE)

  # Lean is strictly smaller
  expect_lt(ncol(lean), ncol(full))
  # Lean is a strict prefix of full (tiered columns are at the front)
  expect_equal(colnames(lean), head(colnames(full), ncol(lean)))

  # Every tier column is present in the lean output
  tier_cols <- c(
    "row_id", "study_id", "author", "year", "outcome",
    paste0(c("all_info", "info_measure", "es_guidance", "es_selected",
             "flags",
             "es_summary", "es_consistency", "n_estimations", "dispersion_es",
             "es", "se", "es_ci_lo", "es_ci_up", "info_used", "measure",
             "overlap_min_max", "diff_min_max", "min_info", "max_info",
             "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
             "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up"),
           "_crude"),
    paste0(c("all_info", "info_measure", "es_guidance", "es_selected",
             "flags",
             "es_summary", "es_consistency", "n_estimations", "dispersion_es",
             "es", "se", "es_ci_lo", "es_ci_up", "info_used", "measure",
             "overlap_min_max", "diff_min_max", "min_info", "max_info",
             "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
             "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up"),
           "_adjusted")
  )
  missing_in_lean <- setdiff(tier_cols, colnames(lean))
  expect_equal(missing_in_lean, character(0),
               info = paste("Tier columns missing from lean output:",
                            paste(missing_in_lean, collapse = ", ")))

  # And NO raw input columns survive (mean_exp/n_exp are typical inputs)
  expect_false("mean_exp" %in% colnames(lean))
  expect_false("n_exp" %in% colnames(lean))
  expect_false("mean_sd_exp" %in% colnames(lean))
  # Pre-computed ES inputs are also dropped
  expect_false("cohen_d" %in% colnames(lean))
  expect_false("hedges_g" %in% colnames(lean))
})

test_that("split_adjusted = FALSE produces unsuffixed tiered layout", {
  obj <- convert_df(df.short, measure = "g", verbose = FALSE,
                    split_adjusted = FALSE)
  res <- summary(obj, flags = TRUE, guidance = TRUE, include_raw = FALSE)

  # All tier columns appear without _crude/_adjusted suffix
  expected <- c("row_id", "study_id", "author", "year", "outcome",
                "all_info", "info_measure", "es_guidance", "es_selected",
                "flags",
                "es_summary", "es_consistency", "n_estimations", "dispersion_es",
                "es", "se", "es_ci_lo", "es_ci_up", "info_used", "measure",
                "overlap_min_max", "diff_min_max", "min_info", "max_info",
                "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
                "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up")
  present <- intersect(expected, colnames(res))
  .expect_contiguous(res, present, "Unsplit tier sequence")
  # First column is row_id
  expect_equal(colnames(res)[1], "row_id")
})
