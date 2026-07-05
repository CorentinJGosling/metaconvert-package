# V24 — Median / quantile ordering check (.validate_input_data)
# A valid five-number summary requires min <= Q1 <= median <= Q3 <= max.
# Violations (e.g. Q3 < median, as in Rosen 2017's printed % change eGFR) are
# mathematically impossible and signal a primary-source error. The check is
# warn-only: values are preserved (faithful transcription), only flagged.

test_that("V24 stays silent on a valid five-number summary", {
  dat <- data.frame(
    min_exp = -68.6, q1_exp = -18.3, med_exp = -5.6, q3_exp = 8.3, max_exp = 56.7,
    min_nexp = -40.9, q1_nexp = -12.6, med_nexp = -5.1, q3_nexp = 7.9, max_nexp = 89.5
  )
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_false(grepl("Non-monotonic quantile summary", res$issues[1]))
})

test_that("V24 stays silent when only a valid median + IQR is given (no min/max)", {
  dat <- data.frame(
    q1_exp = -20.9, med_exp = -12.3, q3_exp = -0.3,
    q1_nexp = -19,  med_nexp = -9.8, q3_nexp = 0
  )
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_equal(res$issues[1], "")
})

test_that("V24 flags Q3 < median (the Rosen 2017 case) on both arms, data preserved", {
  dat <- data.frame(
    min_exp = -68.6, q1_exp = -18.3, med_exp = -5.6, q3_exp = -8.3, max_exp = -56.7,
    min_nexp = -40.9, q1_nexp = -12.6, med_nexp = -5.1, q3_nexp = -7.9, max_nexp = -89.5
  )
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("Non-monotonic quantile summary for 'exp'", res$issues[1]))
  expect_true(grepl("Q3 \\(-8.3\\) < median \\(-5.6\\)", res$issues[1]))
  expect_true(grepl("Non-monotonic quantile summary for 'nexp'", res$issues[1]))
  # warn-only: values must NOT be set to NA
  expect_equal(res$data$q3_exp[1], -8.3)
  expect_equal(res$data$max_exp[1], -56.7)
})

test_that("V24 flags Q1 > median and min > Q1", {
  dat <- data.frame(
    q1_exp = 5, med_exp = 2, q3_exp = 8        # Q1 (5) > median (2)
  )
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_true(grepl("median \\(2\\) < Q1 \\(5\\)", res$issues[1]))

  dat2 <- data.frame(
    min_nexp = 3, q1_nexp = 1, med_nexp = 2, q3_nexp = 4   # min (3) > Q1 (1)
  )
  res2 <- metaConvert:::.validate_input_data(dat2, verbose = FALSE)
  expect_true(grepl("Q1 \\(1\\) < min \\(3\\)", res2$issues[1]))
})

test_that("V24 handles NA gaps without false positives", {
  # median + max only, ascending: no flag
  dat <- data.frame(med_exp = -5, q3_exp = NA, max_exp = 10)
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  expect_equal(res$issues[1], "")
})
