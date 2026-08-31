## V44 -- floor effect in the response distribution.
##
## The check exists because the OBVIOUS gate is the wrong one. Measured over 20
## response distributions (simulations/studies/11_reliability_se.R, block 11e),
## the number of response categories has Spearman rho = -0.16 with the variance
## error c, while the standardised position of the reported total-score mean has
## rho = -0.96. Three FIVE-POINT formats span c = 0.92 to 1.69, so a
## "dichotomous items are bad" rule both misses most of the problem and
## mis-attributes it.
##
## The single most important assertion in this file is the PHQ-9 one: a 0-based
## instrument read as 1-based shifts floor_position by 1/(C-1) ~ 0.25, which
## spans two bands, so defaulting `scale_min` would false-alarm on one of the
## most widely used instruments in the literature. V44 must stay SILENT when
## scale_min is absent rather than guess.

.v44_df <- function(scale_mean, scale_min = 1, n_items = 8, n_cat = 5, ...) {
  d <- data.frame(
    study_id = paste0("S", seq_along(scale_mean)),
    cronbach_alpha = rep(0.88, length(scale_mean)),
    n_sample = rep(300, length(scale_mean)),
    n_items = n_items,
    n_response_categories = n_cat,
    scale_mean = scale_mean,
    stringsAsFactors = FALSE
  )
  if (!is.null(scale_min)) d$scale_min <- scale_min
  extra <- list(...)
  for (nm in names(extra)) d[[nm]] <- extra[[nm]]
  d
}

.v44_flags <- function(d) {
  s <- suppressMessages(summary(
    convert_df(d, measure = "alpha", verbose = FALSE), flags = TRUE))
  fl <- if ("flags_crude" %in% names(s)) s$flags_crude else s$flags
  ifelse(is.na(fl), "", fl)
}

.v44_fired <- function(d) grepl("Floor effect", .v44_flags(d), fixed = TRUE)


test_that("V44 fires on a floored scale and stays silent on a symmetric one", {
  # item mean 1.62 on a 1-5 scale -> floor_position 0.16, the severe band
  expect_true(.v44_fired(.v44_df(13))[1])
  # item mean 3.0 -> floor_position 0.50, benign
  expect_false(.v44_fired(.v44_df(24))[1])
})


test_that("the severity bands follow the measured thresholds", {
  severe <- .v44_flags(.v44_df(13))[1]     # fp = 0.16
  mild   <- .v44_flags(.v44_df(16.7))[1]   # fp = 0.27
  expect_true(grepl("[UNUSUAL] Floor effect", severe, fixed = TRUE))
  expect_true(grepl("[INFO] Floor effect", mild, fixed = TRUE))

  # the band boundary is exclusive at the top: floor_position = 0.28 is NOT
  # flagged, so a row sitting exactly on the threshold stays quiet
  expect_false(.v44_fired(.v44_df(8 * (1 + 0.28 * 4)))[1])
})


test_that("a 0-based instrument is NOT read as a floored 1-based one", {
  # THE FALSE POSITIVE THIS CHECK EXISTS TO AVOID.
  # PHQ-9: 9 items scored 0-3, moderate sample, total mean 15 => item mean 1.67.
  # TRUE floor_position = 1.67/3 = 0.56 (benign).
  # Read as 1-based it would be (1.67-1)/3 = 0.22 -- the mild band, a false alarm.
  phq <- .v44_df(15, scale_min = 0, n_items = 9, n_cat = 4)
  expect_false(.v44_fired(phq)[1])

  # and the same instrument genuinely floored (mild sample) still fires
  phq_floor <- .v44_df(2.5, scale_min = 0, n_items = 9, n_cat = 4)  # item mean 0.28
  expect_true(.v44_fired(phq_floor)[1])
})


test_that("V44 stays silent rather than guessing when scale_min is absent", {
  # Same numbers that fire when scale_min = 1 is stated explicitly.
  d <- .v44_df(13, scale_min = NULL)
  expect_false("scale_min" %in% names(d))
  expect_false(.v44_fired(d)[1])

  # an NA scale_min is the same situation as an absent column
  d2 <- .v44_df(13, scale_min = NA_real_)
  expect_false(.v44_fired(d2)[1])
})


test_that("a study that reported its own SE is not flagged", {
  # V44 is about the CLOSED-FORM standard error. A study supplying its own SE or
  # CI is unaffected by any of this, so flagging it would be pure noise.
  expect_false(.v44_fired(.v44_df(13, cronbach_alpha_se = 0.02))[1])
  expect_false(.v44_fired(.v44_df(13, cronbach_alpha_ci_lo = 0.84,
                                  cronbach_alpha_ci_up = 0.92))[1])
  # but the same row without them does fire, so the guard is what silenced it
  expect_true(.v44_fired(.v44_df(13))[1])
})


test_that("a mean outside its own scale range is reported instead of banded", {
  # item mean 6.25 on a 1-5 scale: an unambiguous transcription error (an ITEM
  # mean entered where a TOTAL was asked for, or the wrong scale_min).
  fl <- .v44_flags(.v44_df(50))[1]
  expect_true(grepl("Scale mean outside its own range", fl, fixed = TRUE))
  # and it must NOT also emit a floor band for the same row
  expect_false(grepl("Floor effect", fl, fixed = TRUE))
})


test_that("V44 routes to the crude scope only", {
  # The leading quoted 'scale_mean' is what .v_flag_matches_scope() reads. A
  # message quoting no column matches EVERY scope and is duplicated into
  # flags_adjusted -- a scope these columns have nothing to do with.
  s <- suppressMessages(summary(
    convert_df(.v44_df(13), measure = "alpha", verbose = FALSE), flags = TRUE))
  skip_if_not("flags_adjusted" %in% names(s))
  expect_true(grepl("Floor effect", s$flags_crude[1], fixed = TRUE))
  expect_false(grepl("Floor effect", ifelse(is.na(s$flags_adjusted[1]), "",
                                            s$flags_adjusted[1]), fixed = TRUE))
})


test_that("the message contains no '; ' (the flag merge splits on it)", {
  fl <- .v44_flags(.v44_df(13))[1]
  parts <- unlist(strsplit(fl, "; ", fixed = TRUE))
  hit <- parts[grepl("Floor effect", parts, fixed = TRUE)]
  expect_length(hit, 1L)
  # the whole message survived as ONE token, so it still carries its severity tag
  expect_true(grepl("^\\[UNUSUAL\\]", hit))
  expect_true(grepl("I-squared", hit, fixed = TRUE))
})


test_that("the thresholds are configurable through flag_options", {
  # default: floor_position 0.50 is benign
  expect_false(.v44_fired(.v44_df(24))[1])
  # raise the mild band above it and the same row is flagged
  s <- suppressMessages(summary(
    convert_df(.v44_df(24), measure = "alpha", verbose = FALSE,
               flag_options = list(floor_position_mild = 0.60,
                                   floor_position_severe = 0.55)),
    flags = TRUE))
  fl <- if ("flags_crude" %in% names(s)) s$flags_crude else s$flags
  expect_true(grepl("Floor effect", ifelse(is.na(fl[1]), "", fl[1]), fixed = TRUE))
})


test_that("V44 is inert for reviews that do not carry the columns", {
  # An ordinary alpha sheet with none of the three columns must behave exactly as
  # before -- V44 is additive, not a new requirement.
  d <- data.frame(study_id = c("S1", "S2"), cronbach_alpha = c(0.88, 0.91),
                  n_sample = c(300, 250), n_items = 8)
  expect_false(any(grepl("Floor effect", .v44_flags(d), fixed = TRUE)))
  expect_false(any(grepl("Scale mean outside", .v44_flags(d), fixed = TRUE)))
})


test_that("the new columns are in the extraction sheet for alpha and omega only", {
  nu <- c("n_response_categories", "scale_mean", "scale_min")
  expect_true(all(nu %in% colnames(data_extraction_sheet(measure = "alpha"))))
  expect_true(all(nu %in% colnames(data_extraction_sheet(measure = "omega"))))
  # NOT in an SMD sheet: these are reliability-block columns, and mean_exp /
  # mean_nexp there are already spoken for by a known-groups validity contrast
  # that a measurement-properties review may carry on the very same row.
  expect_false(any(nu %in% colnames(data_extraction_sheet(measure = "d"))))
})
