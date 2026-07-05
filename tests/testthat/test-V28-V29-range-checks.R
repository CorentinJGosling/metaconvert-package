library(testthat)

# ==============================================================================
# V28 — Reported mean outside its own [min, max]  (.validate_input_data)
# min <= mean <= max is an assumption-free arithmetic identity; a violation is an
# unambiguous transcription/units error. [INVALID], warn-only (data preserved),
# rounding-aware so a mean that rounds to a bound is not falsely flagged.
# ==============================================================================

test_that("V28 stays silent when the mean lies inside its reported range", {
  dat <- data.frame(
    mean_exp = 10, min_exp = 2, max_exp = 20,
    mean_nexp = 9, min_nexp = 3, max_nexp = 18
  )
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("Mean outside its reported range", res$issues[1]))
})

test_that("V28 flags a mean below the reported min", {
  dat <- data.frame(mean_exp = 1, min_exp = 2, max_exp = 20)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("\\[INVALID\\] Mean outside its reported range for 'exp'", res$issues[1]))
  expect_true(grepl("mean = 1 not in \\[min = 2, max = 20\\]", res$issues[1]))
  # warn-only: value must NOT be set to NA
  expect_equal(res$data$mean_exp[1], 1)
})

test_that("V28 flags a mean above the reported max, on the nexp arm", {
  dat <- data.frame(mean_nexp = 25, min_nexp = 2, max_nexp = 20)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("\\[INVALID\\] Mean outside its reported range for 'nexp'", res$issues[1]))
})

test_that("V28 is rounding-aware: a mean that rounds to just below min is NOT flagged", {
  # min reported to 2 dp (2.34), mean reported to 1 dp (2.3). The true mean could
  # be 2.34, so 2.3 is within rounding of the min and must not be flagged.
  dat <- data.frame(mean_exp = 2.3, min_exp = 2.34, max_exp = 9.99)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("Mean outside its reported range", res$issues[1]))
})

test_that("V28 does not double-fire when min > max (that is V24's job)", {
  # Inverted range: V24 catches non-monotonicity; V28 requires lo <= max and stays quiet.
  dat <- data.frame(mean_exp = 5, min_exp = 20, max_exp = 2)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("Mean outside its reported range", res$issues[1]))
})


# ==============================================================================
# V29 — Range/SD (Wan 2014 / Hozo 2005) reverse consistency check
# (max - min)/SD should sit near xi(n) = 2*qnorm((n-0.375)/(n+0.25)); a ratio far
# outside 0.5*xi .. 2.5*xi signals SD-as-SE, variance-as-SD, or a range/IQR mix-up.
# [UNUSUAL], warn-only.
# ==============================================================================

test_that("V29 stays silent when range and SD are mutually consistent", {
  # n = 30 => xi ~ 4.08; range 40 / SD 10 = 4.0, ratio/xi ~ 0.98 (consistent)
  dat <- data.frame(min_exp = 5, max_exp = 45, mean_sd_exp = 10, n_exp = 30)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("Range/SD ratio implausible", res$issues[1]))
})

test_that("V29 flags an SD entered as an SE (range/SD explodes)", {
  # True SD 10 for n = 30, but SE = 10/sqrt(30) ~ 1.83 was entered instead.
  # range 40 / 1.83 ~ 21.9 >> 2.5*xi (~10.2)
  dat <- data.frame(min_exp = 5, max_exp = 45, mean_sd_exp = 1.83, n_exp = 30)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("\\[UNUSUAL\\] Range/SD ratio implausible for 'exp'", res$issues[1]))
  expect_true(grepl("SD-as-SE", res$issues[1]))
})

test_that("V29 flags a variance entered as an SD (range/SD collapses)", {
  # True SD 10, but the variance (100) was entered. range 40 / 100 = 0.4 < 0.5*xi (~2.04)
  dat <- data.frame(min_nexp = 5, max_nexp = 45, mean_sd_nexp = 100, n_nexp = 30)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("\\[UNUSUAL\\] Range/SD ratio implausible for 'nexp'", res$issues[1]))
  # warn-only: SD preserved
  expect_equal(res$data$mean_sd_nexp[1], 100)
})

test_that("V29 needs both a range and a directly reported SD (no false fire otherwise)", {
  # SD present but no range -> skipped
  dat <- data.frame(mean_sd_exp = 1.83, n_exp = 30)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_equal(res$issues[1], "")
})


# ==============================================================================
# V15 — n_exp + n_nexp vs n_sample. Retagged [INVALID] -> [UNUSUAL] because a
# mismatch is routinely legitimate (multi-arm trial, analysis subset).
# ==============================================================================

test_that("V15 flags a sample-size mismatch as [UNUSUAL], not [INVALID]", {
  dat <- data.frame(n_exp = 50, n_nexp = 50, n_sample = 120)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("\\[UNUSUAL\\] Sample size mismatch", res$issues[1]))
  expect_false(grepl("\\[INVALID\\] Sample size", res$issues[1]))
  expect_true(grepl("multi-arm trial", res$issues[1]))
})

test_that("V15 stays silent when the arm sizes add up", {
  dat <- data.frame(n_exp = 50, n_nexp = 50, n_sample = 100)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("Sample size mismatch", res$issues[1]))
})


# ==============================================================================
# V6 severity tag + global severity-prefix guard (via convert_df -> summary)
# Every flag emitted anywhere must begin with one of the four severity tokens.
# ==============================================================================

test_that("the V6 default-r note carries an [INFO] tag and every flag is tagged", {
  dat <- data.frame(
    n_exp = 40, n_nexp = 40,
    mean_pre_exp = 10, mean_pre_sd_exp = 3, mean_exp = 8, mean_sd_exp = 3,
    mean_pre_nexp = 10, mean_pre_sd_nexp = 3, mean_nexp = 9, mean_sd_nexp = 3
  )
  res <- suppressWarnings(suppressMessages(
    metaConvert::convert_df(dat, measure = "g", verbose = FALSE)
  ))
  s <- suppressWarnings(suppressMessages(summary(res, flags = TRUE)))

  flag_cols <- grep("flag", names(s), ignore.case = TRUE, value = TRUE)
  expect_true(length(flag_cols) > 0)

  tokens <- unlist(lapply(flag_cols, function(cc) {
    vals <- as.character(s[[cc]])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    unlist(strsplit(vals, "; "))
  }))
  tokens <- tokens[nzchar(tokens)]

  # The default-r note must now appear tagged [INFO]
  expect_true(any(grepl("^\\[INFO\\] Default r_pre_post", tokens)))

  # Global guard: every emitted flag begins with a known severity token
  bad <- tokens[!grepl("^\\[(INVALID|UNUSUAL|DISCORDANT|INFO)\\]", tokens)]
  expect_equal(bad, character(0))
})
