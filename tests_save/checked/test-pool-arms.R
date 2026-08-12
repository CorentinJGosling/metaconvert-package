test_that("pool method: two experimental arms pooled correctly (means + SDs)", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28),
    mean_exp = c(12.0, 14.0),
    mean_sd_exp = c(3.0, 2.5),
    mean_nexp = c(10.0, 10.0),
    mean_sd_nexp = c(3.0, 3.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(nrow(res), 1)
  # n_exp should be summed

  expect_equal(res$n_exp, 55)
  # n_nexp unchanged
  expect_equal(res$n_nexp, 28)
  # weighted mean: (30*12 + 25*14) / 55
  expect_equal(res$mean_exp, (30 * 12 + 25 * 14) / 55, tolerance = 1e-10)
  # nexp mean unchanged
  expect_equal(res$mean_nexp, 10.0)
  # nexp SD unchanged
  expect_equal(res$mean_sd_nexp, 3.0)
  # Cochrane pooled SD
  n1 <- 30; n2 <- 25; sd1 <- 3.0; sd2 <- 2.5; m1 <- 12.0; m2 <- 14.0
  expected_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2 +
                         n1 * n2 / (n1 + n2) * (m1^2 + m2^2 - 2 * m1 * m2)) /
                        (n1 + n2 - 1))
  expect_equal(res$mean_sd_exp, expected_sd, tolerance = 1e-10)
  # pool_side column should be removed
  expect_false("pool_side" %in% colnames(res))
})

test_that("pool method: nexp side pooled correctly", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("nexp", "nexp"),
    n_exp = c(40, 40),
    n_nexp = c(20, 30),
    mean_exp = c(5.0, 5.0),
    mean_sd_exp = c(2.0, 2.0),
    mean_nexp = c(7.0, 8.0),
    mean_sd_nexp = c(2.5, 3.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(nrow(res), 1)
  expect_equal(res$n_nexp, 50)
  expect_equal(res$n_exp, 40)
  expect_equal(res$mean_exp, 5.0)
  expect_equal(res$mean_sd_exp, 2.0)
  # weighted mean nexp: (20*7 + 30*8) / 50
  expect_equal(res$mean_nexp, (20 * 7 + 30 * 8) / 50, tolerance = 1e-10)
})

test_that("pool method: three arms pooled iteratively", {
  dat <- data.frame(
    study_id = c("S1", "S1", "S1"),
    pool_side = c("exp", "exp", "exp"),
    n_exp = c(20, 25, 30),
    n_nexp = c(28, 28, 28),
    mean_exp = c(10.0, 12.0, 14.0),
    mean_sd_exp = c(3.0, 2.5, 3.5),
    mean_nexp = c(9.0, 9.0, 9.0),
    mean_sd_nexp = c(2.8, 2.8, 2.8)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(nrow(res), 1)
  expect_equal(res$n_exp, 75)
  expect_equal(res$mean_exp, (20 * 10 + 25 * 12 + 30 * 14) / 75, tolerance = 1e-10)
  # SD should be a valid positive number
  expect_true(res$mean_sd_exp > 0)
})

test_that("pool method: 2x2 counts summed correctly", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(50, 60),
    n_nexp = c(55, 55),
    n_cases_exp = c(10, 15),
    n_controls_exp = c(40, 45),
    n_cases_nexp = c(8, 8),
    n_controls_nexp = c(47, 47)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(res$n_exp, 110)
  expect_equal(res$n_cases_exp, 25)
  expect_equal(res$n_controls_exp, 85)
  expect_equal(res$n_cases_nexp, 8)
  expect_equal(res$n_controls_nexp, 47)
})

test_that("pool method: proportions recomputed from pooled counts", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(50, 50),
    n_nexp = c(40, 40),
    n_cases_exp = c(10, 20),
    prop_cases_exp = c(0.2, 0.4),
    n_cases_nexp = c(5, 5),
    prop_cases_nexp = c(0.125, 0.125)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  # pooled prop = sum(cases) / sum(n) = 30 / 100
  expect_equal(res$prop_cases_exp, 0.3, tolerance = 1e-10)
})

test_that("pool method: person-time summed", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 40),
    n_nexp = c(35, 35),
    time_exp = c(100, 150),
    time_nexp = c(200, 200)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(res$time_exp, 250)
  expect_equal(res$time_nexp, 200)
})

test_that("pool method: CIs and medians set to NA", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28),
    mean_exp = c(12.0, 14.0),
    mean_sd_exp = c(3.0, 2.5),
    mean_ci_lo_exp = c(10.5, 12.8),
    mean_ci_up_exp = c(13.5, 15.2),
    med_exp = c(11.5, 13.8),
    mean_nexp = c(10.0, 10.0),
    mean_sd_nexp = c(3.0, 3.0)
  )

  res <- suppressWarnings(
    pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
              method = "pool", verbose = TRUE)
  )

  expect_true(is.na(res$mean_ci_lo_exp))
  expect_true(is.na(res$mean_ci_up_exp))
  expect_true(is.na(res$med_exp))
})

test_that("pool method: SE recomputed from pooled SD", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28),
    mean_exp = c(12.0, 14.0),
    mean_sd_exp = c(3.0, 2.5),
    mean_se_exp = c(0.55, 0.50),
    mean_nexp = c(10.0, 10.0),
    mean_sd_nexp = c(3.0, 3.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  # SE should be SD_pooled / sqrt(n_pooled)
  expected_se <- res$mean_sd_exp / sqrt(res$n_exp)
  expect_equal(res$mean_se_exp, expected_se, tolerance = 1e-10)
})

test_that("pool method: pass-through rows with NA pool_side", {
  dat <- data.frame(
    study_id = c("S1", "S1", "S2"),
    pool_side = c("exp", "exp", NA),
    n_exp = c(30, 25, 40),
    n_nexp = c(28, 28, 35),
    mean_exp = c(12.0, 14.0, 10.0),
    mean_sd_exp = c(3.0, 2.5, 4.0),
    mean_nexp = c(10.0, 10.0, 9.5),
    mean_sd_nexp = c(3.0, 3.0, 3.8)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(nrow(res), 2)
  # S2 row should be unchanged
  s2 <- res[res$study_id == "S2", ]
  expect_equal(s2$n_exp, 40)
  expect_equal(s2$mean_exp, 10.0)
})

test_that("pool method: single-row group warns and passes through", {
  dat <- data.frame(
    study_id = c("S1"),
    pool_side = c("exp"),
    n_exp = c(30),
    n_nexp = c(28),
    mean_exp = c(12.0),
    mean_sd_exp = c(3.0),
    mean_nexp = c(10.0),
    mean_sd_nexp = c(3.0)
  )

  expect_warning(
    res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                     method = "pool", verbose = TRUE),
    "only 1 row"
  )
  expect_equal(nrow(res), 1)
  expect_equal(res$n_exp, 30)
})

test_that("pool method: mixed pool_side within study errors", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "nexp"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28)
  )

  expect_error(
    pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
              method = "pool"),
    "mixed pool_side"
  )
})

test_that("pool method: invalid pool_side values error", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "both"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28)
  )

  expect_error(
    pool_arms(dat, study_id = "study_id", pool_side = "pool_side"),
    "invalid values"
  )
})

test_that("pool method: non-identical shared side warns", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28),
    mean_exp = c(12.0, 14.0),
    mean_sd_exp = c(3.0, 2.5),
    mean_nexp = c(10.0, 11.0),  # different!
    mean_sd_nexp = c(3.0, 3.0)
  )

  expect_warning(
    pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
              method = "pool", verbose = TRUE),
    "differs across arms"
  )
})

# ---- Split method tests ----

test_that("split method: shared arm n divided by k", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28),
    mean_exp = c(12.0, 14.0),
    mean_sd_exp = c(3.0, 2.5),
    mean_nexp = c(10.0, 10.0),
    mean_sd_nexp = c(3.0, 3.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "split", verbose = FALSE)

  expect_equal(nrow(res), 2)
  # Shared side (nexp) divided by 2
  expect_equal(res$n_nexp[1], 14)
  expect_equal(res$n_nexp[2], 14)
  # Pooled side (exp) unchanged
  expect_equal(res$n_exp[1], 30)
  expect_equal(res$n_exp[2], 25)
  # Means unchanged
  expect_equal(res$mean_exp[1], 12.0)
  expect_equal(res$mean_exp[2], 14.0)
})

test_that("split method: odd n handled with remainder distribution", {
  dat <- data.frame(
    study_id = c("S1", "S1", "S1"),
    pool_side = c("exp", "exp", "exp"),
    n_exp = c(20, 25, 30),
    n_nexp = c(29, 29, 29),
    mean_exp = c(10.0, 12.0, 14.0),
    mean_nexp = c(9.0, 9.0, 9.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "split", verbose = FALSE)

  expect_equal(nrow(res), 3)
  # 29 / 3 = 9 remainder 2 → first 2 get 10, last gets 9
  expect_equal(sum(res$n_nexp), 29)
  expect_true(all(res$n_nexp %in% c(9, 10)))
})

test_that("split method: 2x2 counts divided", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(50, 60),
    n_nexp = c(40, 40),
    n_cases_nexp = c(10, 10),
    n_controls_nexp = c(30, 30)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "split", verbose = FALSE)

  expect_equal(res$n_nexp[1], 20)
  expect_equal(res$n_nexp[2], 20)
  expect_equal(res$n_cases_nexp[1], 5)
  expect_equal(res$n_cases_nexp[2], 5)
  expect_equal(res$n_controls_nexp[1], 15)
  expect_equal(res$n_controls_nexp[2], 15)
})

test_that("split method: pass-through rows preserved", {
  dat <- data.frame(
    study_id = c("S1", "S1", "S2"),
    pool_side = c("exp", "exp", NA),
    n_exp = c(30, 25, 40),
    n_nexp = c(28, 28, 35)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "split", verbose = FALSE)

  expect_equal(nrow(res), 3)
  s2 <- res[res$study_id == "S2", ]
  expect_equal(s2$n_exp, 40)
  expect_equal(s2$n_nexp, 35)
})

test_that("split method: n_sample recomputed", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28),
    n_sample = c(58, 53)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "split", verbose = FALSE)

  expect_equal(res$n_sample[1], 30 + 14)
  expect_equal(res$n_sample[2], 25 + 14)
})

test_that("pool method: all NA pool_side returns data unchanged", {
  dat <- data.frame(
    study_id = c("S1", "S2"),
    pool_side = c(NA, NA),
    n_exp = c(30, 25),
    n_nexp = c(28, 28)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(nrow(res), 2)
  expect_equal(res$n_exp, c(30, 25))
  expect_false("pool_side" %in% colnames(res))
})

test_that("pool method: pre-post columns pooled correctly", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 20),
    n_nexp = c(25, 25),
    mean_pre_exp = c(8.0, 9.0),
    mean_pre_sd_exp = c(2.0, 2.5),
    mean_exp = c(12.0, 14.0),
    mean_sd_exp = c(3.0, 2.5),
    mean_pre_nexp = c(8.5, 8.5),
    mean_pre_sd_nexp = c(2.2, 2.2),
    mean_nexp = c(10.0, 10.0),
    mean_sd_nexp = c(3.0, 3.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(nrow(res), 1)
  # pre means: weighted average
  expect_equal(res$mean_pre_exp, (30 * 8 + 20 * 9) / 50, tolerance = 1e-10)
  # pre SDs: Cochrane formula
  expect_true(!is.na(res$mean_pre_sd_exp))
  expect_true(res$mean_pre_sd_exp > 0)
})

test_that("pool method: multiple studies processed independently", {
  dat <- data.frame(
    study_id = c("S1", "S1", "S2", "S2"),
    pool_side = c("exp", "exp", "nexp", "nexp"),
    n_exp = c(30, 25, 40, 40),
    n_nexp = c(28, 28, 20, 30),
    mean_exp = c(12.0, 14.0, 5.0, 5.0),
    mean_sd_exp = c(3.0, 2.5, 2.0, 2.0),
    mean_nexp = c(10.0, 10.0, 7.0, 8.0),
    mean_sd_nexp = c(3.0, 3.0, 2.5, 3.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(nrow(res), 2)
  # S1: exp pooled
  s1 <- res[res$study_id == "S1", ]
  expect_equal(s1$n_exp, 55)
  expect_equal(s1$n_nexp, 28)
  # S2: nexp pooled
  s2 <- res[res$study_id == "S2", ]
  expect_equal(s2$n_exp, 40)
  expect_equal(s2$n_nexp, 50)
})

test_that("pool method: row order preserved", {
  dat <- data.frame(
    study_id = c("S2", "S1", "S1"),
    pool_side = c(NA, "exp", "exp"),
    n_exp = c(40, 30, 25),
    n_nexp = c(35, 28, 28),
    mean_exp = c(10.0, 12.0, 14.0),
    mean_sd_exp = c(4.0, 3.0, 2.5),
    mean_nexp = c(9.5, 10.0, 10.0),
    mean_sd_nexp = c(3.8, 3.0, 3.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  # S2 was first in input, should still come first
  expect_equal(res$study_id[1], "S2")
  expect_equal(res$study_id[2], "S1")
})

test_that("missing column names error correctly", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25)
  )

  expect_error(
    pool_arms(dat, study_id = "bad_col", pool_side = "pool_side"),
    "not found"
  )
  expect_error(
    pool_arms(dat, study_id = "study_id", pool_side = "bad_col"),
    "not found"
  )
})

test_that("non-dataframe input errors", {
  expect_error(
    pool_arms(list(a = 1), study_id = "a", pool_side = "b"),
    "data.frame"
  )
})

test_that("pool method: n_sample recomputed after pooling", {
  dat <- data.frame(
    study_id = c("S1", "S1"),
    pool_side = c("exp", "exp"),
    n_exp = c(30, 25),
    n_nexp = c(28, 28),
    n_sample = c(58, 53),
    mean_exp = c(12.0, 14.0),
    mean_sd_exp = c(3.0, 2.5),
    mean_nexp = c(10.0, 10.0),
    mean_sd_nexp = c(3.0, 3.0)
  )

  res <- pool_arms(dat, study_id = "study_id", pool_side = "pool_side",
                   method = "pool", verbose = FALSE)

  expect_equal(res$n_sample, 55 + 28)
})
