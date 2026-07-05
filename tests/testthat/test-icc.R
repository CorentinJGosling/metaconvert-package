test_that("es_from_icc computes Bonett transformation for ICC(2,1) by default", {
  res <- es_from_icc(icc = 0.80, n_sample = 50, n_measurements = 2, icc_type = "agreement")

  # Default is Bonett: T = ln(1 - ICC)
  expect_equal(res$icc, log(1 - 0.80), tolerance = 1e-10)

  # ICC(2,1) Bonett-transformed SE: sqrt(2*(1+(k-1)*rho)^2 / (k*(k-1)*(n-1)))
  # (1-rho)^2 cancels with delta method derivative 1/(1-rho)^2
  rho <- 0.80; k <- 2; n <- 50
  expected_se <- sqrt(2 * (1 + (k - 1) * rho)^2 / (k * (k - 1) * (n - 1)))
  expect_equal(res$icc_se, expected_se, tolerance = 1e-10)

  expect_equal(res$info_used, "icc")
})

test_that("es_from_icc computes Bonett transformation for ICC(3,1)", {
  res <- es_from_icc(icc = 0.75, n_sample = 100, n_measurements = 3, icc_type = "consistency")

  # ICC(3,1) Bonett-transformed SE: sqrt(2 / ((k-1)*(n-1)))
  # Fully variance-stabilized — no dependence on rho
  rho <- 0.75; k <- 3; n <- 100
  expected_se <- sqrt(2 / ((k - 1) * (n - 1)))
  expect_equal(res$icc_se, expected_se, tolerance = 1e-10)
})

test_that("es_from_icc with icc_to_es='raw' returns raw ICC", {
  res <- es_from_icc(icc = 0.80, n_sample = 50, n_measurements = 2, icc_to_es = "raw")

  # Raw ICC value
  expect_equal(res$icc, 0.80, tolerance = 1e-10)

  # Raw SE: (1 - rho) * bonett_transformed_se (inverse delta method)
  rho <- 0.80; k <- 2; n <- 50
  bonett_transformed_se <- sqrt(2 * (1 + (k - 1) * rho)^2 / (k * (k - 1) * (n - 1)))
  expect_equal(res$icc_se, (1 - rho) * bonett_transformed_se, tolerance = 1e-10)
})

test_that("es_from_icc validates icc_to_es parameter", {
  expect_error(
    es_from_icc(icc = 0.80, n_sample = 50, n_measurements = 2, icc_to_es = "invalid"),
    "not in tolerated values"
  )
})

test_that("es_from_icc handles vectorized input", {
  res <- es_from_icc(
    icc = c(0.70, 0.80, 0.90),
    n_sample = c(50, 100, 200),
    n_measurements = c(2, 3, 2),
    icc_type = c("agreement", "consistency", "agreement")
  )
  expect_equal(nrow(res), 3)
  expect_equal(res$info_used, rep("icc", 3))
})

test_that("es_from_icc defaults to agreement type", {
  res <- es_from_icc(icc = 0.80, n_sample = 50, n_measurements = 2)
  expect_equal(res$icc_type, "agreement")
})

test_that("convert_df with measure='icc' and default Bonett method", {
  dat <- data.frame(
    icc = c(0.80, 0.75),
    n_sample = c(50, 100),
    n_measurements = c(2, 3),
    icc_type = c("agreement", "consistency")
  )
  mc <- convert_df(dat, measure = "icc", verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, flags = FALSE, guidance = FALSE)

  # Default icc_to_es = "bonett"
  expect_equal(s$es[1], log(1 - 0.80), tolerance = 1e-3)
  expect_equal(s$es[2], log(1 - 0.75), tolerance = 1e-3)
  expect_equal(s$info_used, rep("icc", 2))
})

test_that("convert_df with measure='icc' and raw method", {
  dat <- data.frame(
    icc = c(0.80, 0.75),
    n_sample = c(50, 100),
    n_measurements = c(2, 3),
    icc_type = c("agreement", "consistency")
  )
  mc <- convert_df(dat, measure = "icc", icc_to_es = "raw",
                   verbose = FALSE, split_adjusted = FALSE)
  s <- summary(mc, flags = FALSE, guidance = FALSE)

  # Raw ICC values
  expect_equal(s$es[1], 0.80, tolerance = 1e-3)
  expect_equal(s$es[2], 0.75, tolerance = 1e-3)
})

test_that("icc data_extraction_sheet works", {
  dat <- data_extraction_sheet(measure = "icc", extension = "data.frame", verbose = FALSE)
  expect_true("icc" %in% colnames(dat))
  expect_true("n_measurements" %in% colnames(dat))
  expect_true("icc_type" %in% colnames(dat))
})
