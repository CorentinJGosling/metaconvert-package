library(testthat)

# ==============================================================================
# F2 — measure-aware upper-SE ceiling (the former max(2, 5/sqrt(N)) collapsed to
# a flat 2 for every N >= 7, so its N-scaled branch was dead and it could not
# catch the SE ~ 1 that an SD-entered-as-SE produces for a standardised d/g).
# ==============================================================================

test_that("[F2] d/g SE ~1 at large N is now flagged (dead 5/sqrt(N) branch removed)", {
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(0.5, 0.5), se = c(1.5, 0.5), measure = "g", n_total = c(100, 100)
  )
  expect_true(any(grepl("Implausibly large SE", flags[[1]])))   # 1.5 > ceiling 1.0
  expect_false(any(grepl("Implausibly large SE", flags[[2]])))  # 0.5 < 1.0
})

test_that("[F2] ratio measures keep a permissive ceiling (rare-event SEs not flagged)", {
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(0.5, 0.5), se = c(1.5, 2.5), measure = "logor", n_total = c(100, 100)
  )
  expect_false(any(grepl("Implausibly large SE", flags[[1]])))  # 1.5 < ceiling 2.0
  expect_true(any(grepl("Implausibly large SE", flags[[2]])))   # 2.5 > 2.0
})

test_that("[F2] r/z use a tighter ceiling than the ratio measures", {
  flags <- metaConvert:::.flag_se_sample_size(
    es = c(0.3, 0.3), se = c(0.8, 0.4), measure = "r", n_total = c(100, 100)
  )
  expect_true(any(grepl("Implausibly large SE", flags[[1]])))   # 0.8 > 0.6
  expect_false(any(grepl("Implausibly large SE", flags[[2]])))  # 0.4 < 0.6
})


# ==============================================================================
# D2 — logOR/logRR sqrt(N) fallback is scaled by the balanced-event constant so
# fallback rows share the reconstruction pool's ~1 scale. Removes the old bimodal
# pool (recon ~1 vs raw se*sqrt(N) ~4-25) that flagged the minority mode wholesale,
# while KEEPING gross outliers detectable (covered by the existing SE=10 test).
# ==============================================================================

test_that("[D2] logOR sqrt(N)-fallback rows are not falsely flagged (scale-compatible)", {
  opts <- metaConvert:::.default_flag_options()
  # 8 reconstruction-normalised rows (n=100/100) + 2 forced-fallback rows
  # (n_exp > 5000 triggers the sqrt(N) fallback). All SEs consistent with N
  # (balanced ~50% events: SE ~ 4/sqrt(N)). Under the old code the 2 fallback
  # rows landed at ~4 vs a recon median ~1 and were flagged wholesale.
  n_exp  <- c(rep(100, 8), 6000, 6000)
  n_nexp <- n_exp
  N      <- n_exp + n_nexp
  se     <- c(0.27, 0.28, 0.29, 0.28, 0.27, 0.29, 0.28, 0.28,
              4 / sqrt(12000), 4 / sqrt(12000))
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = rep(0.5, 10), se = se, opts = opts, n_total = N,
    measure = "logor", n_exp = n_exp, n_nexp = n_nexp
  )
  for (i in seq_along(flags)) {
    expect_false(any(grepl("SE outlier", flags[[i]])),
      info = paste("row", i, "should not be flagged (fallback rows are scale-compatible)"))
  }
})

test_that("[D2] a gross SE outlier in the fallback pool is still flagged", {
  # Regression guard mirroring the existing SE=10 logOR case: the pathological
  # row's 2x2 reconstruction fails, so it goes to the (now /4-scaled) fallback,
  # but still lands far above the ~1 cluster and must fire.
  opts <- metaConvert:::.default_flag_options()
  es <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.7, 0.5)
  se <- c(0.28, 0.28, 0.29, 0.27, 0.28, 0.28, 0.28, 0.64, 10.0)
  flags <- metaConvert:::.flag_cross_row_outliers(
    es = es, se = se, opts = opts, n_total = rep(200, 9),
    measure = "logor", n_exp = rep(100, 9), n_nexp = rep(100, 9)
  )
  expect_true(any(grepl("SE outlier", flags[[9]])))
  expect_false(any(grepl("SE outlier", flags[[8]])))
})


# ==============================================================================
# Disattenuation — corrected-r CI is now the back-transformed Fisher-z interval,
# so it always respects (-1, 1); a symmetric Wald r-CI could exceed +/-1.
# ==============================================================================

test_that("es_disattenuate corrected-r CI respects (-1, 1) in the normal case", {
  res <- es_disattenuate(r = 0.6, r_se = 0.05,
                          reliability_x = 0.8, reliability_y = 0.8)
  expect_true(is.finite(res$r_corrected_ci_lo) && is.finite(res$r_corrected_ci_up))
  expect_true(res$r_corrected_ci_up <= 1 && res$r_corrected_ci_lo >= -1)
  expect_true(res$r_corrected_ci_lo < res$r_corrected && res$r_corrected < res$r_corrected_ci_up)
})

test_that("es_disattenuate leaves the r CI NA when the corrected r is clamped (>1)", {
  res <- suppressWarnings(es_disattenuate(
    r = 0.80, r_se = 0.03, reliability_x = 0.6, reliability_y = 0.6))
  expect_true(res$r_corrected > 1)                 # documented limitation, point unchanged
  expect_true(is.na(res$r_corrected_ci_lo))
  expect_true(is.na(res$r_corrected_ci_up))
  expect_true(is.finite(res$z_corrected))          # z stays finite (bounded)
})
