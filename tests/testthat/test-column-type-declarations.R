# Every column the extraction sheet DESCRIBES as character must actually be typed as
# character by .check_data().
#
# .check_data() types columns by a grepl() heuristic over their names, and anything it
# does not recognise falls through to numeric. A character column that is registered in
# expected_cols but missed by the heuristic therefore errors with "Non-numeric
# characters in column ..." the first time a user passes real values -- loudly, not
# silently (measured: a registered numeric column given text raises, it does not coerce
# to NA, contrary to what internal_check_data.R's own comment used to claim). Loud is
# good, but the error lands on the USER rather than on CI.
#
# icc_type, omega_type, omega_estimator, alpha_type and pool_side each needed an
# explicit override for exactly this reason, and each was discovered by someone running
# into it. This is a SCANNER, not a pin: it derives the expected type from the sheet's
# own "- character" description suffix, so a column added tomorrow is covered without
# anyone remembering to add it here.

test_that("every character-described column is typed as character by .check_data()", {
  sheet <- data_extraction_sheet(measure = "all")
  desc  <- trimws(as.character(unlist(sheet[1, ])))
  char_cols <- colnames(sheet)[grepl("- character$", desc)]
  expect_gt(length(char_cols), 5)   # the scan must not silently become vacuous

  # a plausible value per column: a bogus one could error for unrelated reasons, and
  # the point is the TYPE, not the semantics
  vals <- c(info_expected = "means_sd", unit_type = "raw_scale",
            alpha_type = "covariance", omega_type = "total",
            omega_estimator = "cfa_bifactor", icc_type = "agreement",
            user_es_original_measure_crude = "d", user_es_target_measure_crude = "d",
            user_es_original_measure_adj = "d", user_es_target_measure_adj = "d")

  base <- data.frame(study_id = paste0("S", 1:3), cronbach_alpha = c(.80, .85, .90),
                     n_sample = 200, n_items = 10, stringsAsFactors = FALSE)
  for (cc in char_cols) {
    d <- base
    d[[cc]] <- rep(if (cc %in% names(vals)) vals[[cc]] else "x", 3)
    err <- tryCatch({
      suppressMessages(suppressWarnings(convert_df(d, measure = "alpha", verbose = FALSE)))
      ""
    }, error = function(e) conditionMessage(e))
    # any unrelated error is fine; a TYPE error is the defect being scanned for
    expect_false(grepl("Non-numeric characters", err, fixed = TRUE),
                 info = paste0("column '", cc, "' is described as character in ",
                               "data_extraction_sheet() but .check_data() types it ",
                               "numeric -- add an override in internal_check_data.R"))
  }
})
