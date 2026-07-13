# V30: baseline statistic copied into the endpoint slot (Kanukula 2024 #9/#10).
# Fires on an exact baseline<->endpoint tie (a copy); silent on legitimately
# similar-but-not-identical timepoints (that is V18's territory).

test_that("V30 fires when both arms' endpoint SD equal their baseline SD (multi-decimal), and escalates when the means also match", {
  df <- data.frame(
    mean_pre_exp     = 24.5, mean_exp     = 24.5,
    mean_pre_nexp    = 25.3, mean_nexp    = 25.3,
    mean_pre_sd_exp  = 6.25, mean_sd_exp  = 6.25,
    mean_pre_sd_nexp = 5.80, mean_sd_nexp = 5.80
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_match(out$issues[1], "Baseline SD copied into endpoint slot",
               info = "V30 should fire when both arms' endpoint SD == baseline SD")
  expect_match(out$issues[1], "both group means are identical across timepoints",
               info = "V30 should escalate when the means are also copied")
})

test_that("V30 fires when both arms match with integer SDs, without the mean-copy escalation", {
  df <- data.frame(
    mean_pre_exp     = 24, mean_exp     = 21,
    mean_pre_nexp    = 25, mean_nexp    = 26,
    mean_pre_sd_exp  = 8,  mean_sd_exp  = 8,
    mean_pre_sd_nexp = 6,  mean_sd_nexp = 6
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_match(out$issues[1], "Baseline SD copied into endpoint slot",
               info = "V30 should fire on a both-arm integer SD tie")
  expect_false(grepl("both group means are identical", out$issues[1]),
               info = "No mean escalation when the endpoint means differ from baseline")
})

test_that("V30 fires on a lone arm only when the tied SD has >= 2 decimals", {
  # exp copied with 2 decimals (distinctive), nexp genuinely different
  df_fire <- data.frame(
    mean_pre_sd_exp  = 6.25, mean_sd_exp  = 6.25,
    mean_pre_sd_nexp = 6,    mean_sd_nexp = 9
  )
  out_fire <- metaConvert:::.validate_input_data(df_fire, verbose = FALSE)
  expect_match(out_fire$issues[1], "Baseline SD copied into endpoint slot",
               info = "A lone multi-decimal SD tie is distinctive enough to fire")

  # exp copied but only an integer (1 decimal) -> not distinctive enough
  df_silent <- data.frame(
    mean_pre_sd_exp  = 8, mean_sd_exp  = 8,
    mean_pre_sd_nexp = 6, mean_sd_nexp = 9
  )
  out_silent <- metaConvert:::.validate_input_data(df_silent, verbose = FALSE)
  expect_false(grepl("Baseline SD copied into endpoint slot", out_silent$issues[1]),
               info = "A lone integer SD tie collides by chance too often to flag")
})

test_that("V30 does NOT fire on legitimately similar-but-not-identical SDs", {
  df <- data.frame(
    mean_pre_sd_exp  = 8.00, mean_sd_exp  = 8.10,
    mean_pre_sd_nexp = 6.00, mean_sd_nexp = 6.20
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_false(grepl("Baseline SD copied into endpoint slot", out$issues[1]),
               info = "V30 must not fire when the SDs merely drift (that is V18's job)")
})

test_that("V30 is silent for an endpoint-only design (no baseline columns)", {
  df <- data.frame(
    mean_exp    = 24, mean_sd_exp  = 8,
    mean_nexp   = 26, mean_sd_nexp = 6
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_false(grepl("Baseline SD copied into endpoint slot", out$issues[1]),
               info = "Nothing to compare against when the baseline slots are empty")
})

test_that("V30 does not fire when a zero baseline SD coincides with a zero endpoint SD", {
  # zero SDs are V10's territory; V30's positive_only guard must exclude them
  df <- data.frame(
    mean_pre_sd_exp  = 0, mean_sd_exp  = 0,
    mean_pre_sd_nexp = 0, mean_sd_nexp = 0
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_false(grepl("Baseline SD copied into endpoint slot", out$issues[1]),
               info = "A 0 == 0 tie is a V10 zero-SD issue, not a V30 copy")
})
