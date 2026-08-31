# =============================================================================
# Roadmap item 1.5 -- every column name data_extraction_sheet() hands to the user
# must be a name the pipeline actually reads.
#
# data_extraction_sheet() is documented as the definitive guide to the expected
# input format: the user fills it in and passes it to convert_df(). If the sheet
# emits a name no consumer reads, the user's value is silently discarded -- there
# is no error, because .check_data() ignores unrecognised columns.
#
# That happened: cols_req emitted "all_info_expected" (R/data_extraction.R:57),
# a string present at that one line and nowhere else in the package, while
# .check_data() (R/internal_check_data.R:47) and functions_summary.R:598 both read
# "info_expected". The documented info_expected output column therefore never
# populated for anyone who used the generated sheet.
#
# These tests are deliberately written against ALL measures rather than the
# default, because the sheet is measure-dependent and a phantom could hide in a
# single branch.
# =============================================================================

.sheet_measures <- function() {
  # Every value data_extraction_sheet() documents, read off the argument's own
  # choices vector so this helper cannot fall behind the function again: it used to
  # enumerate 16 of the 27, leaving the within-group (dw/gw/mdw), partial-correlation
  # (rp/zp), hazard-ratio, omega and log-scale branches unvisited by the two
  # assertions below.
  eval(formals(data_extraction_sheet)$measure)
}

test_that("every name the extraction sheet emits is recognised by .check_data()", {
  # .check_data()'s expected_cols is the authoritative set of names the pipeline
  # reads. A sheet name absent from it is a name whose value goes nowhere.
  recognised <- colnames(metaConvert:::.check_data(
    data.frame(study_id = "s1", n_exp = 10, n_nexp = 10,
               mean_exp = 1, mean_sd_exp = 1, mean_nexp = 0, mean_sd_nexp = 1),
    split_adjusted = TRUE, main_es = TRUE, format = "wide"
  ))

  for (m in .sheet_measures()) {
    nm <- unique(names(data_extraction_sheet(measure = m)))
    orphans <- setdiff(nm, recognised)
    expect_equal(
      orphans, character(0),
      info = paste0("measure = '", m, "': the sheet emits column(s) no consumer ",
                    "reads, so a user filling them in has the value silently ",
                    "discarded: ", paste(orphans, collapse = ", "))
    )
  }
})

test_that("the sheet emits info_expected, not the phantom all_info_expected", {
  for (m in .sheet_measures()) {
    nm <- names(data_extraction_sheet(measure = m))
    expect_true("info_expected" %in% nm,
                info = paste0("measure = '", m, "' is missing info_expected"))
    expect_false("all_info_expected" %in% nm,
                 info = paste0("measure = '", m, "' still emits the phantom ",
                               "all_info_expected"))
  }
})

test_that("a value filled into info_expected survives the round trip to summary()", {
  # The end-to-end behaviour the phantom broke: fill the sheet as generated, and
  # the value must reach the summary output rather than vanishing.
  d <- data.frame(
    study_id = "s1", info_expected = "means_sd",
    n_exp = 50, n_nexp = 50,
    mean_exp = 13, mean_sd_exp = 3, mean_nexp = 11, mean_sd_nexp = 3
  )
  res <- suppressMessages(convert_df(d, measure = "g", verbose = FALSE))
  out <- suppressMessages(summary(res))

  expect_true("info_expected" %in% names(out),
              info = "info_expected is absent from summary() output entirely.")
  expect_equal(as.character(out$info_expected), "means_sd",
               info = "info_expected reached summary() but lost its value.")
})

test_that("the identifier block round-trips as a whole", {
  # cols_req is spliced into every measure branch, so a defect there affects all of
  # them. Check the other five identifiers survive too, not just info_expected.
  d <- data.frame(
    study_id = "s1", author = "Smith", year = 2020,
    predictor = "drug", outcome = "depression", info_expected = "means_sd",
    n_exp = 50, n_nexp = 50,
    mean_exp = 13, mean_sd_exp = 3, mean_nexp = 11, mean_sd_nexp = 3
  )
  res <- suppressMessages(convert_df(d, measure = "g", verbose = FALSE))
  out <- suppressMessages(summary(res))

  for (col in c("study_id", "author", "year", "predictor", "info_expected")) {
    expect_true(col %in% names(out),
                info = paste0("identifier column '", col, "' lost by the pipeline"))
  }
  expect_equal(as.character(out$study_id), "s1")
  expect_equal(as.character(out$info_expected), "means_sd")
})

test_that("the sheet's header and its legend row stay the same length", {
  # cols_req and inf_req are parallel vectors (R/data_extraction.R:57-63). Editing
  # one without the other silently shifts every description by a column.
  for (m in .sheet_measures()) {
    s <- data_extraction_sheet(measure = m)
    expect_equal(ncol(s), length(names(s)),
                 info = paste0("measure = '", m, "'"))
    expect_gt(nrow(s), 0)
    # The legend row must describe the column it sits under: the info_expected
    # description is the one that moved when the name was wrong.
    idx <- which(names(s) == "info_expected")
    expect_length(idx, 1)
    expect_match(as.character(s[[idx]][1]), "expected input data",
                 info = paste0("measure = '", m, "': the legend under ",
                               "info_expected does not describe it"))
  }
})
