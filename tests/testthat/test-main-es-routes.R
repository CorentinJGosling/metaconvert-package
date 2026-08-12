## Regression tests for the route view (convert_df(main_es = FALSE)) and for the
## crashes / flag-scoping defects that went unnoticed because no test in the
## suite had ever executed that branch.

.routes_dat <- function(k = 6) {
  data.frame(
    study_id     = paste0("S", seq_len(k)),
    n_exp        = 50, n_nexp = 50,
    mean_exp     = 10 + seq_len(k) / 4,
    mean_nexp    = 8,
    mean_sd_exp  = 3, mean_sd_nexp = 3,
    mean_se_exp  = 3 / sqrt(50), mean_se_nexp = 3 / sqrt(50),
    student_t    = 3 + seq_len(k) / 10,
    stringsAsFactors = FALSE
  )
}

test_that("main_es = FALSE returns one row per estimation route", {
  d <- .routes_dat()
  s <- suppressMessages(summary(convert_df(d, measure = "g", main_es = FALSE,
                                           verbose = FALSE)))
  expect_gt(nrow(s), nrow(d))
  # every input comparison is represented
  expect_setequal(unique(s$row_id), seq_len(nrow(d)))
  # routes within a comparison are distinct methods
  per_row <- tapply(s$info_used, s$row_id, function(v) length(unique(v)))
  expect_true(all(per_row >= 2))
})

test_that("B1: es_selected = minimum/maximum is rejected with main_es = FALSE", {
  d <- .routes_dat()
  expect_error(
    convert_df(d, measure = "g", main_es = FALSE, es_selected = "minimum",
               verbose = FALSE),
    "incompatible with 'main_es = FALSE'"
  )
  expect_error(
    convert_df(d, measure = "g", main_es = FALSE, es_selected = "maximum",
               verbose = FALSE),
    "incompatible with 'main_es = FALSE'"
  )
  # still valid with main_es = TRUE
  expect_no_error(
    suppressMessages(summary(convert_df(d, measure = "g",
                                        es_selected = "minimum", verbose = FALSE)))
  )
})

test_that("B1: route rows keep distinct effect sizes (no min/max collapse)", {
  d <- .routes_dat()
  s <- suppressMessages(summary(convert_df(d, measure = "g", main_es = FALSE,
                                           verbose = FALSE)))
  n_distinct <- tapply(as.numeric(s$es), s$row_id, function(v) length(unique(v)))
  expect_true(any(n_distinct > 1))
})

test_that("B0: user columns named 'blank' or 'dat_long' do not break summary()", {
  d <- .routes_dat()
  d_blank <- d; d_blank$blank <- 1
  expect_no_error(suppressMessages(summary(convert_df(d_blank, measure = "g",
                                                      verbose = FALSE))))
  d_dl <- d; d_dl$dat_long <- 1
  expect_no_error(suppressMessages(summary(convert_df(d_dl, measure = "g",
                                                      main_es = FALSE, verbose = FALSE))))
  # and the numbers are unchanged by the presence of the extra column
  a <- suppressMessages(summary(convert_df(d, measure = "g", verbose = FALSE)))
  b <- suppressMessages(summary(convert_df(d_blank, measure = "g", verbose = FALSE)))
  expect_equal(as.numeric(a$es_crude), as.numeric(b$es_crude))
})

test_that("B2/B3: unestimable comparisons survive and receive guidance", {
  dz <- data.frame(study_id = paste0("Z", 1:4), n_exp = 30, n_nexp = 30)
  s <- suppressMessages(summary(convert_df(dz, measure = "g", main_es = FALSE,
                                           verbose = FALSE)))
  expect_equal(nrow(s), 4)
  expect_true(all(is.na(s$es)))
  expect_true(all(nchar(s$es_guidance) > 0))
})

test_that("B3: no comparison disappears from the route view", {
  s_wide <- suppressMessages(summary(convert_df(df.haza, measure = "g",
                                                verbose = FALSE)))
  s_long <- suppressMessages(summary(convert_df(df.haza, measure = "g",
                                                main_es = FALSE, verbose = FALSE)))
  expect_equal(length(unique(s_long$row_id)), length(unique(s_wide$row_id)))
  expect_equal(length(unique(s_long$row_id)), nrow(df.haza))
})

test_that("B4: distinct studies are not flagged as duplicates of themselves", {
  d <- .routes_dat(12)
  s <- suppressMessages(summary(convert_df(d, measure = "g", main_es = FALSE,
                                           verbose = FALSE)))
  expect_false(any(grepl("Duplicate study_id", s$flags, fixed = TRUE)))
})

test_that("B5: a cross-row ES outlier is still detected in the route view", {
  d <- data.frame(
    study_id     = c(paste0("H", 1:8), "OUT"),
    n_exp        = 60, n_nexp = 60,
    mean_exp     = c(10.0, 10.1, 9.9, 10.2, 9.8, 10.0, 10.1, 9.9, 25.0),
    mean_nexp    = 8,
    mean_sd_exp  = 3, mean_sd_nexp = 3,
    # only the outlier carries the extra routes -- the configuration that used to
    # let its own replicates inflate the IQR and mask it
    mean_se_exp  = c(rep(NA, 8), 3 / sqrt(60)),
    mean_se_nexp = c(rep(NA, 8), 3 / sqrt(60)),
    student_t    = c(rep(NA, 8), (25 - 8) / (3 * sqrt(2 / 60))),
    stringsAsFactors = FALSE
  )
  s_wide <- suppressMessages(summary(convert_df(d, measure = "g", verbose = FALSE)))
  s_long <- suppressMessages(summary(convert_df(d, measure = "g", main_es = FALSE,
                                                verbose = FALSE)))
  expect_true(any(grepl("outlier", s_wide$flags_crude, ignore.case = TRUE)))
  flagged <- unique(s_long$study_id[grepl("outlier", s_long$flags, ignore.case = TRUE)])
  expect_true("OUT" %in% flagged)
})

test_that("B8: split_adjusted = FALSE is honoured in the route view", {
  d <- data.frame(
    study_id = "A", n_exp = 50, n_nexp = 50,
    mean_exp = 10, mean_nexp = 8, mean_sd_exp = 3, mean_sd_nexp = 3,
    ancova_mean_exp = 10.1, ancova_mean_nexp = 8.1,
    ancova_mean_sd_exp = 2.8, ancova_mean_sd_nexp = 2.8,
    cov_outcome_r = 0.5, n_cov_ancova = 1
  )
  s_split <- suppressMessages(summary(convert_df(d, measure = "g", main_es = FALSE,
                                                 split_adjusted = TRUE, verbose = FALSE)))
  s_one   <- suppressMessages(summary(convert_df(d, measure = "g", main_es = FALSE,
                                                 split_adjusted = FALSE, verbose = FALSE)))
  expect_true("adjusted_input" %in% colnames(s_split))
  expect_false("adjusted_input" %in% colnames(s_one))
  # single pool: crude and adjusted routes are compared against one another
  expect_true(all(as.numeric(s_one$n_estimations) ==
                    length(unique(s_one$info_used))))
})

test_that("the default (main_es = TRUE) output is unaffected", {
  d <- .routes_dat()
  s <- suppressMessages(summary(convert_df(d, measure = "g", verbose = FALSE)))
  expect_equal(nrow(s), nrow(d))
  expect_true("es_crude" %in% colnames(s))
  expect_false(any(grepl("FORMULA", s$flags_crude, fixed = TRUE)))
})
