# Regression test for the lipsey_cooper Fisher-z variance (A5).
#
# vz for z = atanh(r), r = d/sqrt(d^2 + a) with a = 1/(p*(1-p)) is vd/(d^2 + a) -- the
# exact Fisher transform of the (correct) r variance vr = a^2 vd/(d^2+a)^3. The code had
# vd/(vd + a) (the value esc::convert_d2r also returns), which is NOT Fisher-consistent.
# metaConvert intentionally uses the consistent delta value (it diverges from esc for
# large |d|; they agree as d -> 0). Only the z-scale SE/CI on the lipsey_cooper family
# are affected (r, r_se, r CI, z point estimate were already correct).

test_that("lipsey_cooper z_se is the Fisher transform of the r variance (vd/(d^2+a))", {
  d <- 1.5; n1 <- 100; n2 <- 100
  es <- es_from_cohen_d(cohen_d = d, n_exp = n1, n_nexp = n2, smd_to_cor = "lipsey_cooper")
  a  <- (n1 + n2)^2 / (n1 * n2)                 # == 1/(p*(1-p))
  vd <- 1 / n1 + 1 / n2 + d^2 / (2 * (n1 + n2)) # borenstein default variance
  r  <- d / sqrt(d^2 + a)
  vr <- a^2 * vd / (d^2 + a)^3

  expect_equal(es$z_se[1], sqrt(vd / (d^2 + a)), tolerance = 1e-6)
  # exact Fisher consistency with the r variance: vz = vr/(1-r^2)^2
  expect_equal(es$z_se[1], sqrt(vr) / (1 - r^2), tolerance = 1e-6)
  # r_se was always correct and is unchanged
  expect_equal(es$r_se[1], sqrt(vr), tolerance = 1e-6)
  # guard against a regression to the old esc-matching value vd/(vd+a)
  expect_lt(es$z_se[1], sqrt(vd / (vd + a)))
})

test_that("default viechtbauer z variance (1/(n-1)) is untouched", {
  es <- es_from_cohen_d(cohen_d = 1.5, n_exp = 100, n_nexp = 100)  # default smd_to_cor
  expect_equal(es$z_se[1], sqrt(1 / (100 + 100 - 1)), tolerance = 1e-6)
})

test_that("lipsey_cooper z fix propagates through the OR -> d -> z path", {
  es <- es_from_or_se(or = 15, logor_se = 0.3, n_exp = 100, n_nexp = 100,
                      or_to_cor = "lipsey_cooper")
  d <- es$d[1]; a <- (100 + 100)^2 / (100 * 100); vd <- (es$d_se[1])^2
  expect_equal(es$z_se[1], sqrt(vd / (d^2 + a)), tolerance = 1e-4)
})
