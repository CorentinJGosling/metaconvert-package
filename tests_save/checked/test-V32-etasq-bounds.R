# V32: eta-squared (crude etasq and adjusted etasq_adj) must lie in [0, 1).
# A value >= 1 makes the Cohen's d conversion diverge (d = 2*sqrt(eta2/(1-eta2)),
# or the implied ANCOVA F = eta2*df/(1-eta2)) to Inf/NaN, so it is flagged INVALID
# and set to missing under correct_inputs = TRUE.

test_that("V32 flags etasq >= 1 and sets it to missing (correct_inputs = TRUE)", {
  df <- data.frame(etasq = c(0.28, 1, 1.5), n_exp = 20, n_nexp = 22)
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)

  expect_false(grepl("Out-of-range eta-squared", out$issues[1]),
               info = "A valid eta-squared (0.28) must not fire")
  expect_match(out$issues[2], "Out-of-range eta-squared: 'etasq' = 1",
               info = "etasq == 1 diverges and must be flagged")
  expect_match(out$issues[3], "Out-of-range eta-squared: 'etasq' = 1.5")
  expect_true(is.na(out$data$etasq[2]) && is.na(out$data$etasq[3]),
              info = "Out-of-range values are set to missing")
  expect_equal(out$data$etasq[1], 0.28)
})

test_that("V32 covers the adjusted etasq_adj column too", {
  df <- data.frame(etasq_adj = c(0.4, 1.2), n_exp = 20, n_nexp = 22)
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_false(grepl("Out-of-range eta-squared", out$issues[1]))
  expect_match(out$issues[2], "Out-of-range eta-squared: 'etasq_adj' = 1.2")
  expect_true(is.na(out$data$etasq_adj[2]))
})

test_that("V32 preserves the value under correct_inputs = FALSE but still flags", {
  df <- data.frame(etasq = 1.3, n_exp = 20, n_nexp = 22)
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE, correct_inputs = FALSE)
  expect_match(out$issues[1], "Out-of-range eta-squared: 'etasq' = 1.3")
  expect_equal(out$data$etasq[1], 1.3, info = "Value preserved when correct_inputs = FALSE")
})

test_that("V32 boundary: etasq just below 1 is valid, exactly 1 is invalid", {
  df <- data.frame(etasq = c(0.999, 1.0), n_exp = 20, n_nexp = 22)
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_false(grepl("Out-of-range eta-squared", out$issues[1]),
               info = "0.999 is a valid (if large) eta-squared")
  expect_match(out$issues[2], "Out-of-range eta-squared")
})
