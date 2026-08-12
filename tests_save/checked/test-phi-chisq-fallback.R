# Regression test for the re-enabled phi/chisq -> R/Z/D/G fallback (A9).
#
# When the 2x2 table cannot be reconstructed (n_cases / n_exp missing), phi IS the Pearson
# correlation of the two binary variables, so R/Z/D/G are delivered via es_from_pearson_r()
# with STANDARD sampling variances (z_se = 1/sqrt(n-3), r_se = (1-r^2)/sqrt(n-1)) -- NOT
# the earlier ad hoc formula that inflated the SE by up to ~22% at large phi. OR/RR/NNT
# need the 2x2 margins and stay NA. The fallback must not override rows where the 2x2 was
# reconstructed.

test_that("phi-only row yields R/Z/D/G identical to es_from_pearson_r (standard SEs)", {
  e <- es_from_phi(phi = 0.3, n_sample = 120)
  pr <- es_from_pearson_r(pearson_r = 0.3, n_sample = 120)
  for (col in c("r", "r_se", "z", "z_se", "d", "d_se")) {
    expect_equal(e[[col]], pr[[col]], tolerance = 1e-9, info = col)
  }
  # standard closed-form SEs
  expect_equal(e$z_se, 1 / sqrt(120 - 3), tolerance = 1e-9)
  expect_equal(e$r_se, (1 - 0.3^2) / sqrt(120 - 1), tolerance = 1e-9)
  expect_equal(e$r, 0.3, tolerance = 1e-12)
  # OR/RR/NNT require the 2x2 margins -> NA in the phi-only case
  expect_true(is.na(e$logor))
})

test_that("full-info phi still uses the 2x2 route (fallback does not override)", {
  ef <- es_from_phi(phi = 0.3, n_sample = 120, n_cases = 20, n_exp = 40)
  expect_false(is.na(ef$logor))                       # OR available from the 2x2
  expect_false(isTRUE(all.equal(ef$r, 0.3)))          # r is the 2x2-derived value, not phi
})

test_that("chisq-only row uses r = sqrt(chisq/n) with standard SEs", {
  ec <- es_from_chisq(chisq = 4.21, n_sample = 78)
  expect_equal(ec$r, sqrt(4.21 / 78), tolerance = 1e-9)
  expect_equal(ec$z_se, 1 / sqrt(78 - 3), tolerance = 1e-9)
  expect_true(is.na(ec$logor))
})

test_that("reverse_phi flips the sign of the fallback estimates", {
  er <- es_from_phi(phi = 0.3, n_sample = 120, reverse_phi = TRUE)
  expect_equal(er$r, -0.3, tolerance = 1e-9)
  expect_equal(er$z, atanh(-0.3), tolerance = 1e-9)
})

test_that("convert_df delivers r for a phi-only study (pipeline)", {
  s <- suppressWarnings(suppressMessages(
    summary(convert_df(data.frame(phi = 0.3, n_sample = 120), measure = "r"))
  ))
  es_col <- grep("^es_crude$|^es$", colnames(s), value = TRUE)[1]
  expect_equal(s[[es_col]][1], 0.3, tolerance = 1e-3)
})
