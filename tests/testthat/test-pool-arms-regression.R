# Regression tests for pool_arms() bugs fixed in main_pool_arms.R. Each targets a
# specific defect the original test-pool-arms.R could not catch (it always used a
# character study_id placed in the first column, with fully-tagged studies).

## ---- factor study_id: no phantom all-NA rows -------------------------------

test_that("pool: factor study_id with unused levels produces no phantom rows", {
  dat <- data.frame(
    study_id  = factor(c("S1", "S1", "S2"), levels = c("S1", "S2", "S3", "S4")),
    pool_side = c("exp", "exp", NA),
    n_exp = c(30, 25, 40), n_nexp = c(28, 28, 35),
    mean_exp = c(12, 14, 10), mean_sd_exp = c(3, 2.5, 4),
    mean_nexp = c(10, 10, 9.5), mean_sd_nexp = c(3, 3, 3.8)
  )
  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)
  expect_equal(nrow(res), 2)              # S1 pooled + S2 passthrough, no phantoms
  expect_false(any(is.na(res$study_id)))
  expect_equal(res$n_exp[res$study_id == "S1"], 55)
})

test_that("split: factor study_id with unused levels produces no phantom rows", {
  dat <- data.frame(
    study_id  = factor(c("S1", "S1", "S2"), levels = c("S1", "S2", "S3")),
    pool_side = c("exp", "exp", NA),
    n_exp = c(30, 25, 40), n_nexp = c(28, 28, 35)
  )
  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "split", verbose = FALSE)
  expect_equal(nrow(res), 3)
  expect_false(any(is.na(res$study_id)))
})

## ---- stale n_sample when the shared arm's n is unavailable ------------------

test_that("pool: n_sample is NA'd (not stale) when the shared-side n is absent", {
  dat <- data.frame(
    study_id = c("S1", "S1"), pool_side = c("exp", "exp"),
    n_exp = c(30, 25), n_sample = c(58, 53),   # no n_nexp column
    mean_exp = c(12.5, 14.2), mean_sd_exp = c(3.1, 2.8)
  )
  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)
  expect_equal(res$n_exp, 55)
  expect_true(is.na(res$n_sample))           # was 58 (stale first-row value)
})

## ---- split: warns on inconsistent shared-side values -----------------------

test_that("split: warns when a shared-side column differs across arms", {
  dat <- data.frame(
    study_id = c("S1", "S1"), pool_side = c("exp", "exp"),
    n_exp = c(30, 25), n_nexp = c(28, 40),          # shared n disagrees
    mean_nexp = c(10, 10), mean_sd_nexp = c(3, 3)   # (only n_nexp differs)
  )
  expect_warning(
    pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
              method = "split", verbose = TRUE),
    "differs across arms"
  )
})

## ---- mixed NA / non-NA pool_side within a study warns ----------------------

test_that("mixed tagged/untagged rows within one study emit a warning", {
  dat <- data.frame(
    study_id = c("S1", "S1", "S1"), pool_side = c("exp", "exp", NA),
    n_exp = c(30, 25, 20), n_nexp = c(28, 28, 28),
    mean_exp = c(12, 14, 9), mean_sd_exp = c(3, 2.5, 2),
    mean_nexp = c(10, 10, 10), mean_sd_nexp = c(3, 3, 3)
  )
  expect_warning(
    pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
              method = "pool", verbose = TRUE),
    "both tagged .* and untagged"
  )
  # a genuine passthrough (NA row in a DIFFERENT study) must NOT warn
  dat2 <- data.frame(
    study_id = c("S1", "S1", "S2"), pool_side = c("exp", "exp", NA),
    n_exp = c(30, 25, 40), n_nexp = c(28, 28, 35),
    mean_exp = c(12, 14, 10), mean_sd_exp = c(3, 2.5, 4),
    mean_nexp = c(10, 10, 9.5), mean_sd_nexp = c(3, 3, 3.8)
  )
  expect_no_warning(
    pool_arms(dat2, study_id = "study_id", pool_side = "pool_side",
              method = "pool", verbose = TRUE)
  )
})

## ---- warnings name the study_id, not the first column ----------------------

test_that("pool warnings name the study via study_id even when it is not column 1", {
  dat <- data.frame(
    note      = c("aaa", "bbb"),          # first column, NOT the id
    study_id  = c("RealStudy", "RealStudy"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25), n_nexp = c(28, 28),
    mean_exp = c(12, 14), mean_sd_exp = c(3, 2.5),
    mean_nexp = c(10, 11), mean_sd_nexp = c(3, 3)   # shared side differs -> warns
  )
  # capture warnings explicitly
  msgs <- character(0)
  withCallingHandlers(
    pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
              method = "pool", verbose = TRUE),
    warning = function(cnd) { msgs[[length(msgs) + 1]] <<- conditionMessage(cnd)
                              invokeRestart("muffleWarning") }
  )
  expect_true(any(grepl("RealStudy", msgs)))
  expect_false(any(grepl("aaa", msgs)))
})

## ---- user column named .orig_row is preserved ------------------------------

test_that("a user column named .orig_row is not clobbered", {
  dat <- data.frame(
    study_id = c("S1", "S1", "S2"), pool_side = c("exp", "exp", NA),
    n_exp = c(30, 25, 40), n_nexp = c(28, 28, 35),
    .orig_row = c(111, 222, 333),
    mean_exp = c(12, 14, 10), mean_sd_exp = c(3, 2.5, 4),
    mean_nexp = c(10, 10, 9.5), mean_sd_nexp = c(3, 3, 3.8)
  )
  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)
  expect_true(".orig_row" %in% colnames(res))
  expect_equal(res$.orig_row[res$study_id == "S2"], 333)
})

## ---- SE-error-bar-only studies keep a pooled interval ----------------------

test_that("pool: a study with only SE error bars is not dropped", {
  dat <- data.frame(
    study_id = c("S1", "S1"), pool_side = c("exp", "exp"),
    n_exp = c(30, 25), n_nexp = c(28, 28),
    plot_mean_exp = c(12, 14),
    plot_mean_se_lo_exp = c(11.5, 13.6),
    plot_mean_se_up_exp = c(12.5, 14.4),
    mean_nexp = c(10, 10), mean_sd_nexp = c(3, 3)
  )
  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  m_pooled <- (30 * 12 + 25 * 14) / 55
  expect_equal(res$plot_mean_exp, m_pooled)
  expect_false(is.na(res$plot_mean_se_lo_exp))
  expect_false(is.na(res$plot_mean_se_up_exp))
  # symmetric around the pooled mean
  expect_equal(res$plot_mean_se_up_exp - m_pooled,
               m_pooled - res$plot_mean_se_lo_exp, tolerance = 1e-10)

  # matches the intended derivation: SD_i = SE_bar_i * sqrt(n_i), Cochrane pool,
  # SE_pooled = SD_pooled / sqrt(n_pooled)
  sd1 <- 0.5 * sqrt(30); sd2 <- 0.4 * sqrt(25)
  sd_pool <- sqrt(((30 - 1) * sd1^2 + (25 - 1) * sd2^2 +
                     30 * 25 / 55 * (12^2 + 14^2 - 2 * 12 * 14)) / (55 - 1))
  expect_equal(res$plot_mean_se_up_exp - m_pooled, sd_pool / sqrt(55),
               tolerance = 1e-10)
})
