# V33: reported proportion vs its own counts (V33a) and GRIM-style
# achievable-fraction (V33b). Both are warn-only (data preserved) and
# rounding-aware via .count_decimals(); messages quote 'prop' so they route to
# the crude scope only.

v33 <- function(dat, ...) {
  metaConvert:::.validate_input_data(dat, verbose = FALSE, ...)
}

test_that("V33a fires when prop contradicts n_cases/n_sample beyond rounding", {
  dat <- data.frame(prop = 0.5, n_cases = 10, n_sample = 40)
  res <- v33(dat)
  expect_match(res$issues[1], "\\[INVALID\\] Proportion inconsistent with counts")
  # warn only: all three values preserved
  expect_equal(res$data$prop, 0.5)
  expect_equal(res$data$n_cases, 10)
  expect_equal(res$data$n_sample, 40)
})

test_that("V33a is rounding-aware and silent on consistent data", {
  # exact match
  expect_false(grepl("Proportion inconsistent",
    v33(data.frame(prop = 0.35, n_cases = 7, n_sample = 20))$issues[1]))
  # 0.3 reported at 1 dp is consistent with 7/20 = 0.35 (rounding error 0.05)
  expect_false(grepl("Proportion inconsistent",
    v33(data.frame(prop = 0.3, n_cases = 7, n_sample = 20))$issues[1]))
  # but 0.31 at 2 dp is NOT consistent with 7/20 = 0.35
  expect_match(
    v33(data.frame(prop = 0.31, n_cases = 7, n_sample = 20))$issues[1],
    "Proportion inconsistent")
})

test_that("V33b fires when prop is unachievable for any integer count", {
  # 0.33 on N = 8: nearest achievable 3/8 = 0.375 (or 2/8 = 0.25), gap >> 0.005
  res <- v33(data.frame(prop = 0.33, n_sample = 8))
  expect_match(res$issues[1], "\\[UNUSUAL\\] Proportion not achievable")
  expect_match(res$issues[1], "3/8", fixed = TRUE)
  expect_equal(res$data$prop, 0.33)  # preserved
})

test_that("V33b is silent when achievable, when N is large, or N non-integer", {
  expect_false(grepl("not achievable",
    v33(data.frame(prop = 0.375, n_sample = 8))$issues[1]))
  expect_false(grepl("not achievable",
    v33(data.frame(prop = 0.33, n_sample = 100))$issues[1]))
  # non-integer N: skipped
  expect_false(grepl("not achievable",
    v33(data.frame(prop = 0.33, n_sample = 8.5))$issues[1]))
  # exactly achievable: 3/10
  expect_false(grepl("not achievable",
    v33(data.frame(prop = 0.3, n_sample = 10))$issues[1]))
})

test_that("V33a takes precedence over V33b when counts are present", {
  # inconsistent AND unachievable: only the (stronger) V33a message is emitted
  res <- v33(data.frame(prop = 0.33, n_cases = 3, n_sample = 8))
  expect_match(res$issues[1], "Proportion inconsistent")
  expect_false(grepl("not achievable", res$issues[1]))
})

test_that("V33 defers to V5/V11 on their territory", {
  # n_cases > n_sample is V5's job (both set to NA under correct_inputs)
  res <- v33(data.frame(prop = 0.5, n_cases = 50, n_sample = 40))
  expect_match(res$issues[1], "Inconsistent input")
  expect_false(grepl("Proportion inconsistent", res$issues[1]))
  # out-of-range prop is V11's job
  res2 <- v33(data.frame(prop = 1.5, n_cases = 10, n_sample = 40))
  expect_match(res2$issues[1], "Out-of-range proportion")
  expect_false(grepl("Proportion inconsistent", res2$issues[1]))
})

test_that("V33 messages route to the crude scope and split cleanly", {
  res <- v33(data.frame(prop = 0.5, n_cases = 10, n_sample = 40))
  msg <- strsplit(res$issues[1], "; ", fixed = TRUE)[[1]]
  msg <- msg[grepl("Proportion inconsistent", msg)]
  expect_length(msg, 1)
  expect_true(metaConvert:::.v_flag_matches_scope(msg, "_crude"))
  expect_false(metaConvert:::.v_flag_matches_scope(msg, "_adjusted"))
  expect_true(grepl("^\\[INVALID\\]", msg))
})
