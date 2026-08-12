## `cor_to_smd` against its source paper.
##
## Mathur MB & VanderWeele TJ (2020) "A Simple, Interpretable Conversion from
## Pearson's Correlation to Cohen's d for Continuous Exposures", Epidemiology
## 31(2):e16-e18  (data-raw/R to SMD/R to SMD - Mathur 2020 Epidemio.pdf).
##
## The paper's two equations are two of the three shipped routes:
##   Eq (1.1)  d = 2r/sqrt(1-r^2)                          -> cor_to_smd = "cooper"
##   Eq (1.2)  d = r*Delta / (s_x*sqrt(1-r^2))              -> cor_to_smd = "mathur"
##             SE(d) = |d| * sqrt(1/(r^2*(N-3)) + 1/(2*(N-1)))
##
## Eq (1.1) is stated by the authors to be derived for a POINT-BISERIAL r (binary
## exposure); Eq (1.2) is their conversion for a continuous exposure, reported per
## Delta units of X. These are different estimands, not competing approximations.

test_that("mathur reproduces Mathur & VanderWeele Eq (1.2), estimate and SE", {
  for (r in c(-0.5, -0.15, 0.2, 0.35, 0.7)) {
    for (N in c(30, 120, 500)) {
      for (sx in c(0.8, 2.4)) {
        for (D in c(0.5, 1, 3)) {
          d_paper  <- r * D / (sx * sqrt(1 - r^2))
          se_paper <- abs(d_paper) * sqrt(1 / (r^2 * (N - 3)) + 1 / (2 * (N - 1)))
          got <- es_from_pearson_r(pearson_r = r, n_sample = N, sd_iv = sx,
                                   unit_increase_iv = D, unit_type = "value",
                                   cor_to_smd = "mathur")
          expect_equal(got$d, d_paper, tolerance = 1e-12,
                       info = sprintf("r=%s N=%s sx=%s D=%s", r, N, sx, D))
          expect_equal(got$d_se, se_paper, tolerance = 1e-12,
                       info = sprintf("r=%s N=%s sx=%s D=%s", r, N, sx, D))
        }
      }
    }
  }
})

test_that("cooper reproduces Eq (1.1)", {
  for (r in c(-0.5, 0.2, 0.35, 0.7)) {
    got <- es_from_pearson_r(pearson_r = r, n_sample = 120, cor_to_smd = "cooper")
    expect_equal(got$d, 2 * r / sqrt(1 - r^2), tolerance = 1e-12)
  }
})

test_that("Eq (1.1) equals Eq (1.2) at Delta = 2 s_x, but their SEs do not", {
  # The paper: applying Eq (1.1) to a continuous-X correlation "coincides with the
  # effect size associated with an increase in X of two standard deviations" --
  # and is equally explicit that even then "the standard error estimates in
  # Equations (1.1) and (1.2) will, in general, still not coincide".
  for (r in c(-0.4, 0.25, 0.6)) {
    for (N in c(40, 200)) {
      for (sx in c(0.8, 2.4, 7)) {
        cp <- es_from_pearson_r(pearson_r = r, n_sample = N, cor_to_smd = "cooper")
        mt <- es_from_pearson_r(pearson_r = r, n_sample = N, sd_iv = sx,
                                unit_increase_iv = 2, unit_type = "sd",
                                cor_to_smd = "mathur")
        expect_equal(mt$d, cp$d, tolerance = 1e-12,
                     info = sprintf("r=%s N=%s sx=%s", r, N, sx))
        # sanity: the equivalence must be free of s_x, since Delta = 2*s_x cancels
        expect_false(isTRUE(all.equal(mt$d_se, cp$d_se)),
                     info = sprintf("SEs should differ: r=%s N=%s sx=%s", r, N, sx))
      }
    }
  }
})

test_that("unit_type selects raw units vs SD units of the exposure", {
  # Documentation defect, not an arithmetic one: ?convert_df says unit_type must be
  # "sd" or "value", the shipped default is "raw_scale", and nothing validates the
  # argument -- any non-"sd" value means "Delta is in raw units of X", which IS the
  # paper's Delta. Pinned here so a future validation change is a deliberate one.
  r <- 0.4; N <- 100; sx <- 3; D <- 2
  raw <- es_from_pearson_r(pearson_r = r, n_sample = N, sd_iv = sx,
                           unit_increase_iv = D, unit_type = "value",
                           cor_to_smd = "mathur")
  dflt <- es_from_pearson_r(pearson_r = r, n_sample = N, sd_iv = sx,
                            unit_increase_iv = D, unit_type = "raw_scale",
                            cor_to_smd = "mathur")
  insd <- es_from_pearson_r(pearson_r = r, n_sample = N, sd_iv = sx,
                            unit_increase_iv = D, unit_type = "sd",
                            cor_to_smd = "mathur")
  expect_equal(dflt$d, raw$d, tolerance = 1e-12)          # default == raw units
  expect_equal(insd$d, raw$d * sx, tolerance = 1e-12)     # sd units scale by s_x
})
