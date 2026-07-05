test_that("es_from_spearman_rho converts correctly to Pearson r", {
  # Known conversion: r_p = 2 * sin(pi/6 * r_s)
  res <- es_from_spearman_rho(spearman_r = 0.55, n_sample = 100)
  expected_r <- 2 * sin(pi / 6 * 0.55)
  expect_equal(res$r, expected_r, tolerance = 1e-10)

  # r_s = 0 should give r_p = 0

  res0 <- es_from_spearman_rho(spearman_r = 0, n_sample = 50)
  expect_equal(res0$r, 0, tolerance = 1e-10)

  # r_s = 1 should give r_p = 2*sin(pi/6) = 1
  res1 <- es_from_spearman_rho(spearman_r = 1, n_sample = 50)
  expect_equal(res1$r, 1, tolerance = 1e-10)
})

test_that("es_from_spearman_rho SE uses delta method", {
  # SE via delta method: sqrt((pi/3 * cos(pi/6 * r_s))^2 * (1-r_s^2)^2 / (n-1))
  r_s <- 0.55
  n <- 100
  deriv <- (pi / 3) * cos(pi / 6 * r_s)
  var_spearman <- (1 - r_s^2)^2 / (n - 1)
  expected_se <- sqrt(deriv^2 * var_spearman)

  res <- es_from_spearman_rho(spearman_r = r_s, n_sample = n)
  expect_equal(res$r_se, expected_se, tolerance = 1e-10)
})

test_that("es_from_spearman_rho produces Fisher's z", {
  res <- es_from_spearman_rho(spearman_r = 0.40, n_sample = 80)
  expected_r <- 2 * sin(pi / 6 * 0.40)
  expected_z <- atanh(expected_r)
  expect_equal(res$z, expected_z, tolerance = 1e-10)
})

test_that("es_from_spearman_rho handles vectorized input", {
  res <- es_from_spearman_rho(
    spearman_r = c(0.3, 0.5, 0.7),
    n_sample = c(50, 100, 200)
  )
  expect_equal(nrow(res), 3)
  expect_equal(res$info_used, rep("spearman_r", 3))
})

test_that("es_from_spearman_rho handles reversal", {
  res_fwd <- es_from_spearman_rho(spearman_r = 0.50, n_sample = 100)
  res_rev <- es_from_spearman_rho(spearman_r = 0.50, n_sample = 100,
                                   reverse_spearman_r = TRUE)
  expect_equal(res_fwd$r, -res_rev$r, tolerance = 1e-10)
})

test_that("Spearman integrates into convert_df hierarchy for measure='r'", {
  dat <- data.frame(spearman_r = 0.55, n_sample = 100)
  mc <- convert_df(dat, measure = "r", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, flags = FALSE, guidance = FALSE)
  expect_equal(s$info_used, "spearman_r")
  expect_false(is.na(s$es))
})

test_that("Spearman feeds into d/g hierarchy", {
  dat <- data.frame(spearman_r = 0.55, n_sample = 100)
  mc <- convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, flags = FALSE, guidance = FALSE)
  expect_equal(s$info_used, "spearman_r")
  expect_false(is.na(s$es))
})
