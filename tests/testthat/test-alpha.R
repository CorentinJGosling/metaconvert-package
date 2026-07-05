test_that("es_from_cronbach_alpha computes Bonett transformation by default", {
  res <- es_from_cronbach_alpha(cronbach_alpha = 0.85, n_sample = 200, n_items = 10)

  # Default is Bonett: T = ln(1 - alpha)
  expect_equal(res$alpha, log(1 - 0.85), tolerance = 1e-10)
  # Bonett SE: sqrt(2k / ((k-1)(n-2))) per Bonett (2002, JEBS)
  expect_equal(res$alpha_se, sqrt(2 * 10 / (9 * 198)), tolerance = 1e-10)

  expect_equal(res$info_used, "cronbach_alpha")
})

test_that("es_from_cronbach_alpha with alpha_to_es='raw' returns raw alpha", {
  res <- es_from_cronbach_alpha(cronbach_alpha = 0.85, n_sample = 200, n_items = 10,
                                alpha_to_es = "raw")

  # Raw alpha value
  expect_equal(res$alpha, 0.85, tolerance = 1e-10)
  # Raw SE: (1 - alpha) * bonett_se
  bonett_se <- sqrt(2 * 10 / (9 * 198))
  expect_equal(res$alpha_se, (1 - 0.85) * bonett_se, tolerance = 1e-10)
})

test_that("es_from_cronbach_alpha validates alpha_to_es parameter", {
  expect_error(
    es_from_cronbach_alpha(cronbach_alpha = 0.85, n_sample = 200, n_items = 10,
                           alpha_to_es = "invalid"),
    "not in tolerated values"
  )
})

test_that("es_from_cronbach_alpha handles edge cases", {
  # Alpha = 0
  res0 <- es_from_cronbach_alpha(cronbach_alpha = 0, n_sample = 100, n_items = 5)
  expect_equal(res0$alpha, 0, tolerance = 1e-10)  # ln(1-0) = ln(1) = 0

  # Alpha very close to 1 (should give large negative Bonett)
  res_high <- es_from_cronbach_alpha(cronbach_alpha = 0.99, n_sample = 100, n_items = 10)
  expect_true(res_high$alpha < -4)  # ln(0.01) = -4.605

  # Alpha = 1 (boundary case - bonett = -Inf)
  res1 <- es_from_cronbach_alpha(cronbach_alpha = 1, n_sample = 100, n_items = 10)
  expect_equal(res1$alpha, -Inf)
})

test_that("es_from_cronbach_alpha handles vectorized input", {
  res <- es_from_cronbach_alpha(
    cronbach_alpha = c(0.80, 0.85, 0.90),
    n_sample = c(100, 200, 300),
    n_items = c(5, 10, 20)
  )
  expect_equal(nrow(res), 3)
  expect_equal(res$info_used, rep("cronbach_alpha", 3))
})

test_that("convert_df with measure='alpha' and default Bonett method", {
  dat <- data.frame(
    cronbach_alpha = c(0.85, 0.70),
    n_sample = c(200, 100),
    n_items = c(10, 5)
  )
  mc <- convert_df(dat, measure = "alpha", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, flags = FALSE, guidance = FALSE)

  # Default alpha_to_es = "bonett"
  expect_equal(s$es[1], log(1 - 0.85), tolerance = 1e-3)
  expect_equal(s$es[2], log(1 - 0.70), tolerance = 1e-3)
  expect_equal(s$info_used, rep("cronbach_alpha", 2))
})

test_that("convert_df with measure='alpha' and raw method", {
  dat <- data.frame(
    cronbach_alpha = c(0.85, 0.70),
    n_sample = c(200, 100),
    n_items = c(10, 5)
  )
  mc <- convert_df(dat, measure = "alpha", alpha_to_es = "raw",
                   verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, flags = FALSE, guidance = FALSE)

  # Raw alpha values
  expect_equal(s$es[1], 0.85, tolerance = 1e-3)
  expect_equal(s$es[2], 0.70, tolerance = 1e-3)
})

test_that("alpha data_extraction_sheet works", {
  dat <- data_extraction_sheet(measure = "alpha", extension = "data.frame", verbose = FALSE)
  expect_true("cronbach_alpha" %in% colnames(dat))
  expect_true("n_items" %in% colnames(dat))
  expect_true("n_sample" %in% colnames(dat))
})
