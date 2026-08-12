# Regression test for the `discard` column feature (A4).
#
# .check_data() contained `situation <- situation[-removed_rows]` referencing an
# undefined variable, so ANY truthy value in the documented `discard` column aborted
# the whole convert_df() run with "object 'situation' not found". The row removal is
# done by the subset() call; the stray line was removed.

test_that("discard column drops flagged rows without crashing", {
  df <- data.frame(
    mean_exp = c(5, 6, 7), mean_nexp = c(4, 5, 6),
    mean_sd_exp = c(1, 1, 1), mean_sd_nexp = c(1, 1, 1),
    n_exp = c(20, 20, 20), n_nexp = c(20, 20, 20),
    discard = c("no", "yes", "no")
  )
  expect_no_error(
    s <- suppressWarnings(suppressMessages(summary(convert_df(df, measure = "d"))))
  )
  # rows 1 and 3 retained, row 2 discarded
  expect_equal(nrow(s), 2L)
})

test_that("discard accepts the documented spellings and the all-keep case is unchanged", {
  base <- data.frame(
    mean_exp = c(5, 6, 7), mean_nexp = c(4, 5, 6),
    mean_sd_exp = c(1, 1, 1), mean_sd_nexp = c(1, 1, 1),
    n_exp = c(20, 20, 20), n_nexp = c(20, 20, 20)
  )
  for (kw in c("yes", "Yes", "remove", "removed")) {
    df <- base; df$discard <- c("no", kw, "no")
    s <- suppressWarnings(suppressMessages(summary(convert_df(df, measure = "d"))))
    expect_equal(nrow(s), 2L, info = paste("discard keyword:", kw))
  }
  # no discards -> all three rows kept
  df0 <- base; df0$discard <- c("no", "no", "no")
  s0 <- suppressWarnings(suppressMessages(summary(convert_df(df0, measure = "d"))))
  expect_equal(nrow(s0), 3L)
})
