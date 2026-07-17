# =============================================================================
# Intraclass correlation (es_from_icc) — promoted from tests_save/checked/
# (P26, audit digest [16]).
#
# DE-CIRCULARIZED: the archived file re-typed the source SE expression and
# asserted it was "confirmed by Monte Carlo" (no Monte Carlo exists in the
# repo). The SE expectations below are instead HAND-DERIVED through the
# independent F-route for the two-way consistency ICC(3,1):
#
#   ICC(3,1)-hat = (F0 - 1) / (F0 + k - 1),   F0 = MSR / MSE
#   =>  1 - ICC = k / (F0 + k - 1),  ln(1 - ICC) = ln k - ln(F0 + k - 1)
#   F0 ~ theta * F(n-1, (n-1)(k-1)),  theta = E[MSR]/E[MSE]
#                                           = (1 + (k-1) rho) / (1 - rho)
#   Var(ln F) ~ 2/(n-1) + 2/((n-1)(k-1))   (large-sample log-F variance)
#   delta method at F0 = theta, with theta + k - 1 = k / (1 - rho):
#   Var(ln(1 - ICC)) = theta^2 * Var(ln F) * (1 - rho)^2 / k^2
#
# The derivation never types the package's own closed form; the test computes
# the chain step by step.
#
# HONESTY NOTE (audit digest [14]): the package applies this SAME
# consistency-derived SE to icc_type = "agreement" (ICC(2,1)) rows. That is a
# leading-order proxy that assumes negligible between-rater variance; when
# rater variance > 0 it understates the ICC(2,1) sampling SD (the exact
# ICC(2,1) variance is not computable from the summary inputs icc/n/k). The
# tests below therefore anchor the SE through the ICC(3,1) derivation and
# only assert that the agreement path emits that same documented proxy.
# =============================================================================

# hand-derived F-route variance of ln(1 - ICC(3,1)) (see header)
.orc_icc31_var <- function(rho, k, n) {
  theta <- (1 + (k - 1) * rho) / (1 - rho)
  var_lnF <- 2 / (n - 1) + 2 / ((n - 1) * (k - 1))
  theta^2 * var_lnF * (1 - rho)^2 / k^2
}

test_that("es_from_icc (bonett, consistency) matches the independent F-route derivation", {
  rho <- 0.75; k <- 3; n <- 100
  res <- es_from_icc(icc = rho, n_sample = n, n_measurements = k,
                     icc_type = "consistency")

  # Bonett point estimate: T = ln(1 - ICC) (hand identity)
  expect_equal(res$icc, log(1 - rho), tolerance = 1e-10)
  expect_equal(res$icc_se, sqrt(.orc_icc31_var(rho, k, n)), tolerance = 1e-10)
  expect_equal(res$info_used, "icc")
})

test_that("es_from_icc (bonett, agreement default) emits the same consistency-derived proxy SE", {
  # The package uses one shared SE for both icc_type values (the ICC(3,1)
  # F-route value above). For agreement data with non-negligible rater
  # variance this is anti-conservative (digest [14]) — this test pins the
  # documented proxy behaviour, it does not certify the agreement variance.
  rho <- 0.80; k <- 2; n <- 50
  res_agr <- es_from_icc(icc = rho, n_sample = n, n_measurements = k,
                         icc_type = "agreement")
  res_con <- es_from_icc(icc = rho, n_sample = n, n_measurements = k,
                         icc_type = "consistency")

  expect_equal(res_agr$icc, log(1 - rho), tolerance = 1e-10)
  expect_equal(res_agr$icc_se, sqrt(.orc_icc31_var(rho, k, n)), tolerance = 1e-10)
  # icc_type has no computational effect on the SE (shared formula)
  expect_equal(res_agr$icc_se, res_con$icc_se, tolerance = 1e-15)
})

test_that("es_from_icc with icc_to_es='raw' returns raw ICC with delta back-transformed SE", {
  # If T = ln(1 - ICC) then ICC = 1 - e^T and |dICC/dT| = 1 - ICC, so the
  # raw-scale SE must be (1 - rho) times the transformed-scale SE (delta
  # method, derived here rather than re-typed from the source).
  rho <- 0.80; k <- 2; n <- 50
  res <- es_from_icc(icc = rho, n_sample = n, n_measurements = k, icc_to_es = "raw")

  expect_equal(res$icc, rho, tolerance = 1e-10)
  expect_equal(res$icc_se, (1 - rho) * sqrt(.orc_icc31_var(rho, k, n)),
               tolerance = 1e-10)
})

test_that("es_from_icc validates icc_to_es parameter", {
  expect_error(
    es_from_icc(icc = 0.80, n_sample = 50, n_measurements = 2, icc_to_es = "invalid"),
    "not in tolerated values"
  )
})

test_that("es_from_icc handles vectorized input (per-row icc_type)", {
  res <- es_from_icc(
    icc = c(0.70, 0.80, 0.90),
    n_sample = c(50, 100, 200),
    n_measurements = c(2, 3, 2),
    icc_type = c("agreement", "consistency", "agreement")
  )
  expect_equal(nrow(res), 3)
  expect_equal(res$info_used, rep("icc", 3))
  # each row's SE equals its own F-route value
  expect_equal(res$icc_se,
               sqrt(mapply(.orc_icc31_var, c(0.70, 0.80, 0.90), c(2, 3, 2),
                           c(50, 100, 200))),
               tolerance = 1e-10)
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
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  # Default icc_to_es = "bonett": T = ln(1 - ICC), SE from the F-route.
  # summary() rounds to `digits` decimal places, hence the 1e-8 tolerance.
  expect_equal(s$es, log(1 - dat$icc), tolerance = 1e-8)
  expect_equal(s$se,
               sqrt(mapply(.orc_icc31_var, dat$icc, dat$n_measurements, dat$n_sample)),
               tolerance = 1e-8)
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
  s <- summary(mc, digits = 11, flags = FALSE, guidance = FALSE)

  expect_equal(s$es, dat$icc, tolerance = 1e-8)
  expect_equal(s$se,
               (1 - dat$icc) *
                 sqrt(mapply(.orc_icc31_var, dat$icc, dat$n_measurements, dat$n_sample)),
               tolerance = 1e-8)
})

test_that("icc data_extraction_sheet works", {
  dat <- data_extraction_sheet(measure = "icc", extension = "data.frame", verbose = FALSE)
  expect_true("icc" %in% colnames(dat))
  expect_true("n_measurements" %in% colnames(dat))
  expect_true("icc_type" %in% colnames(dat))
})
