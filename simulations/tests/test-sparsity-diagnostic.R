## =============================================================================
## Roadmap item 4.1 -- the se_ratio > 1 pattern in the binary studies is driven by
## SPARSITY, not by the +0.5 continuity correction the README suspected.
##
## The README carried this as "an unexplained pattern, not yet diagnosed", and
## named the continuity correction as the likely structural cause. The causal
## direction was backwards: the correction only ever applies in the sparse region,
## and removing it makes the ratio WORSE (1.20-1.71 against 1.08-1.32). What the
## split below shows is that where zero cells essentially never occur, every route
## is calibrated -- so the correction is suppressing part of the inflation, not
## producing it.
##
## These tests pin the diagnostic itself, so the claim in the README is
## reproducible rather than asserted.
## =============================================================================

test_that("zero_cell_prob() behaves like a probability of an empty cell", {
  # a big, balanced, common-outcome table: a zero cell is essentially impossible
  expect_lt(zero_cell_prob(n = 300, p_exp = 0.5, br = 0.50, rr = 1), 1e-10)
  # a small, rare-outcome table: near-certain
  expect_gt(zero_cell_prob(n = 50, p_exp = 0.5, br = 0.01, rr = 1), 0.5)
  # monotone in n at a fixed risk
  p <- vapply(c(20, 50, 100, 300, 1000),
              function(k) zero_cell_prob(k, 0.5, 0.05, 1), numeric(1))
  expect_true(all(diff(p) < 0))
  # rr * br >= 1 is not a probability; those conditions are dropped by the grid
  expect_true(is.na(zero_cell_prob(n = 100, p_exp = 0.5, br = 0.5, rr = 2)))
  # vectorised
  expect_length(zero_cell_prob(c(50, 300), c(0.5, 0.5), c(0.01, 0.5), c(1, 1)), 2)
  expect_true(all(zero_cell_prob(c(50, 300), 0.5, c(0.01, 0.5), 1) <= 1))
})

test_that("the dense region really is dense", {
  # If a condition with a material zero-cell probability leaked in, the whole
  # argument would be circular.
  agg <- utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                   "04_or_to_rr_nrep1000.csv"))
  agg <- agg[agg$target == "sample", ]
  pz <- zero_cell_prob(agg$n, agg$p_exp, agg$br, agg$rr)
  dense <- !is.na(pz) & pz < 1e-4
  expect_gt(sum(dense), 0)
  expect_lt(max(pz[dense]), 1e-4)
  expect_gt(min(pz[!dense & !is.na(pz)]), 1e-4)
})

test_that("study 04: every method is calibrated where zero cells do not occur", {
  agg <- utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                   "04_or_to_rr_nrep1000.csv"))
  r <- sparsity_split(agg)

  # This is the finding: dense ~ 1 for EVERY route, sparse inflated for every route.
  expect_true(all(abs(r$dense_mean - 1) < 0.02),
              info = paste(sprintf("%s=%.3f", r$method, r$dense_mean), collapse = " "))
  expect_true(all(r$sparse_mean > 1.2),
              info = paste(sprintf("%s=%.3f", r$method, r$sparse_mean), collapse = " "))
  # and the inflation is a real gap, not noise
  expect_true(all(r$sparse_mean - r$dense_mean > 0.2))
  expect_equal(unique(r$dense_n), 105L)
  expect_equal(unique(r$sparse_n), 255L)
})

test_that("transpose is included, which is what makes the correction the wrong suspect", {
  # `transpose` copies logor_se through untouched. If a per-method bug explained the
  # pattern, transpose would be exempt; it is not.
  agg <- utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                   "04_or_to_rr_nrep1000.csv"))
  r <- sparsity_split(agg)
  tr <- r[r$method == "transpose", ]
  expect_equal(nrow(tr), 1L)
  expect_lt(abs(tr$dense_mean - 1), 0.02)
  expect_gt(tr$sparse_mean, 1.2)
})

test_that("metafor_conv2x2 is a SELECTION artifact, not a miscalibration", {
  agg <- utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                   "04_or_to_rr_nrep1000.csv"))
  r <- sparsity_split(agg)
  mf <- r[r$method == "metafor_conv2x2", ]
  expect_equal(nrow(mf), 1L)
  # calibrated where it can reconstruct at all ...
  expect_lt(abs(mf$dense_mean - 1), 0.02)
  # ... but it fails on a large share of the sparse cells, so its full-grid figure
  # is computed on the survivors
  expect_gt(mf$nonest_max, 0.5)
  expect_gt(mf$nonest_mean, 0.05)
  # every other route reconstructs essentially always, which is why only this one
  # needs the caveat
  others <- r[r$method != "metafor_conv2x2", ]
  expect_true(all(others$nonest_mean < 0.01))
})

test_that("metaumbrella_exp is no longer an outlier (items 1.1 / 1.8 / 1.10)", {
  # The roadmap's original 4.1 text said "of the two methods that stay
  # miscalibrated, only metaumbrella_exp does". After the reconstruction fixes the
  # two variants are the same route.
  agg <- utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                   "04_or_to_rr_nrep1000.csv"))
  r <- sparsity_split(agg)
  a <- r[r$method == "metaumbrella_exp", ]
  b <- r[r$method == "metaumbrella_cases", ]
  expect_equal(a$dense_mean, b$dense_mean)
  expect_equal(a$sparse_mean, b$sparse_mean)
})

test_that("study 05 shows the same split, once the selection cases are set aside", {
  agg <- utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                   "05_rr_to_or_nrep1000.csv"))
  r <- sparsity_split(agg)

  # THE SELECTION CAVEAT APPLIES TO A SECOND METHOD, and this test is where that was
  # found. `grant` in study 05 has dense_mean 1.10, outside the calibrated band --
  # but its non-estimability reaches 1.000 (whole cells where rr * br_guess >= 1, and
  # item 1.4 now returns NA rather than an estimate with a NaN variance), so even its
  # dense-region figure is computed on survivors. It is the study-05 analogue of
  # metafor_conv2x2 in study 04, and must not be read as a calibration failure.
  clean <- r[r$nonest_mean < 0.01, ]
  selected <- r[r$nonest_mean >= 0.01, ]

  expect_gt(nrow(clean), 0)
  expect_true(all(abs(clean$dense_mean - 1) < 0.03),
              info = paste(sprintf("%s=%.3f", clean$method, clean$dense_mean), collapse = " "))
  expect_true(all(clean$sparse_mean > 1.1))

  # and the method failing the band is exactly the one carrying the selection
  expect_equal(selected$method, "grant")
  expect_gt(selected$nonest_max, 0.9)
})

test_that("exactly two routes across the two binary studies need the selection caveat", {
  # Stated as a fact about the shipped grids so a change is visible: a route whose
  # summary is conditional on having succeeded must be read with nonest_rate beside
  # it, never as a plain calibration number.
  s04 <- sparsity_split(utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                                  "04_or_to_rr_nrep1000.csv")))
  s05 <- sparsity_split(utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                                  "05_rr_to_or_nrep1000.csv")))
  flagged <- c(s04$method[s04$nonest_mean >= 0.05], s05$method[s05$nonest_mean >= 0.05])
  expect_setequal(flagged, c("metafor_conv2x2", "grant"))
})

test_that("sparsity_split() refuses an aggregate that is not a binary study", {
  agg <- utils::read.csv(file.path(.sim_root, "data", "aggregated",
                                   "01a_smd_to_cor_r_nrep1000.csv"))
  expect_error(sparsity_split(agg), "binary studies")
})
