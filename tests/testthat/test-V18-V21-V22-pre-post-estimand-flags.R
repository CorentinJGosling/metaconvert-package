test_that("V18 fires when SD_baseline / SD_endpoint < 0.7", {
  # PETRA-like: severity-screened baseline truncates SD_baseline
  df <- data.frame(
    n_exp            = 80,  n_nexp = 80,
    mean_pre_exp     = 30,  mean_pre_nexp     = 30,
    mean_pre_sd_exp  = 4,   mean_pre_sd_nexp  = 4,
    mean_exp         = 24,  mean_nexp         = 27,
    mean_sd_exp      = 8,   mean_sd_nexp      = 8
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_match(out$issues[1], "SD_baseline / SD_endpoint",
               info = "V18 should fire when sd_ratio = 0.5 < 0.7")
  expect_match(out$issues[1], "0.50",
               info = "V18 should report computed sd_ratio")
})

test_that("V18 does NOT fire when SDs are roughly equal", {
  df <- data.frame(
    n_exp            = 80,  n_nexp = 80,
    mean_pre_exp     = 30,  mean_pre_nexp     = 30,
    mean_pre_sd_exp  = 8,   mean_pre_sd_nexp  = 8,
    mean_exp         = 24,  mean_nexp         = 27,
    mean_sd_exp      = 8,   mean_sd_nexp      = 8
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_false(grepl("SD_baseline / SD_endpoint", out$issues[1]),
               info = "V18 should not fire at sd_ratio = 1.0")
})

test_that("V21 fires when |baseline imbalance| > 0.30", {
  df <- data.frame(
    n_exp            = 80,  n_nexp = 80,
    mean_pre_exp     = 32,  mean_pre_nexp     = 30,
    mean_pre_sd_exp  = 5,   mean_pre_sd_nexp  = 5,
    mean_exp         = 25,  mean_nexp         = 28,
    mean_sd_exp      = 5,   mean_sd_nexp      = 5
  )
  # imb = (32-30)/sqrt(((79*25)+(79*25))/158) = 2/5 = 0.40 > 0.30
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_match(out$issues[1], "standardised baseline imbalance",
               info = "V21 should fire at imb = 0.40 > 0.30")
  expect_match(out$issues[1], "0.40")
})

test_that("V21 does NOT fire when imbalance is small", {
  df <- data.frame(
    n_exp            = 80,  n_nexp = 80,
    mean_pre_exp     = 30,  mean_pre_nexp     = 30,
    mean_pre_sd_exp  = 5,   mean_pre_sd_nexp  = 5,
    mean_exp         = 25,  mean_nexp         = 28,
    mean_sd_exp      = 5,   mean_sd_nexp      = 5
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_false(grepl("standardised baseline imbalance", out$issues[1]),
               info = "V21 should not fire at imb = 0")
})

test_that("V22 fires when ANCOVA residual SD is provided without cov_outcome_r", {
  df <- data.frame(
    n_exp                = 80,  n_nexp = 80,
    ancova_mean_exp      = 24,  ancova_mean_nexp = 27,
    ancova_mean_sd_pooled = 6
    # cov_outcome_r missing
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_match(out$issues[1], "ANCOVA residual SD",
               info = "V22 should fire when ANCOVA residual SD lacks cov_outcome_r")
  expect_match(out$issues[1], "back-transformation")
})

test_that("V22 does NOT fire when cov_outcome_r is provided", {
  df <- data.frame(
    n_exp                = 80,  n_nexp = 80,
    ancova_mean_exp      = 24,  ancova_mean_nexp = 27,
    ancova_mean_sd_pooled = 6,
    cov_outcome_r        = 0.4
  )
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE)
  expect_false(grepl("ANCOVA residual SD", out$issues[1]),
               info = "V22 should not fire when cov_outcome_r is present and valid")
})

test_that("Configurable thresholds via flag_options work as expected", {
  df <- data.frame(
    n_exp            = 80,  n_nexp = 80,
    mean_pre_exp     = 30,  mean_pre_nexp     = 30,
    mean_pre_sd_exp  = 6,   mean_pre_sd_nexp  = 6,
    mean_exp         = 24,  mean_nexp         = 27,
    mean_sd_exp      = 8,   mean_sd_nexp      = 8
  )
  # sd_ratio = 6/8 = 0.75; below 0.9 threshold -> should fire
  out <- metaConvert:::.validate_input_data(df, verbose = FALSE,
                                            sd_ratio_bl_ep_min = 0.9)
  expect_match(out$issues[1], "SD_baseline / SD_endpoint",
               info = "Tightening V18 threshold to 0.9 should expose sd_ratio = 0.75")
  # sd_ratio = 0.75 above 0.6 threshold -> should NOT fire
  out2 <- metaConvert:::.validate_input_data(df, verbose = FALSE,
                                             sd_ratio_bl_ep_min = 0.6)
  expect_false(grepl("SD_baseline / SD_endpoint", out2$issues[1]),
               info = "Loosening V18 threshold to 0.6 should hide sd_ratio = 0.75")
})
