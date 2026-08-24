## =============================================================================
## Roadmap item 3.1 -- the psychometric family had ZERO simulation coverage.
##
## Nine studies covered the conversion routes; grep -i "alpha|icc|omega" over
## studies/ and R/ returned nothing. That gap mattered more than a missing row in
## a coverage table, because the reliability audit produced several numbers that
## lived only in commit messages and code comments -- each a simulation run once
## and thrown away. A number recorded that way cannot be re-checked when the code
## moves underneath it, which is exactly what happened to the "coverage around
## 0.74-0.76" claim study 10a now replaces.
##
## These tests do NOT re-run the Monte Carlo (that is run_10's job and takes
## minutes). They assert the things that make the study trustworthy and that a
## later edit could quietly break:
##
##   * the study is registered and honours the runner contract;
##   * the DGP really produces the design it claims (a two-way model whose rater
##     variance is the requested share, and an ICC estimator that recovers the
##     requested value) -- otherwise the coverage table measures the wrong thing;
##   * rater_share = 0 is in the grid, because that control is what makes the
##     degradation attributable to rater variance rather than to a coding error;
##   * n spans the range over which coverage is claimed to WORSEN;
##   * the deterministic helpers return what the documentation cites them for.
## =============================================================================

.needs_pkg <- function() if (!METACONVERT_AVAILABLE) skip("metaConvert not loadable")


test_that("study 10 is registered and honours the runner contract", {
  expect_true(exists("run_10", envir = .GlobalEnv))
  fmls <- names(formals(get("run_10", envir = .GlobalEnv)))
  expect_true(all(c("nrep", "cores") %in% fmls))
  # run_everything() discovers by regex, so the NAME has to match too
  expect_true("run_10" %in% ls(envir = .GlobalEnv, pattern = "^run_[0-9]{2}$"))
})


test_that("the 10a grid keeps the control cell and the n range the claim needs", {
  g <- build_grid_10a()
  # rater_share = 0 is the one-way case the package formula is EXACT for. Without
  # it in the grid the study measures a degradation with nothing to compare it
  # against, and a reader cannot tell a real defect from a broken DGP.
  expect_true(0 %in% g$rater_share)
  expect_true(any(g$rater_share > 0))
  # the documented claim is that coverage WORSENS with n, so the grid has to span
  # enough n for that to be visible rather than asserted
  expect_true(all(c(20, 1000) %in% g$n))
  expect_gte(length(unique(g$n)), 4L)
  expect_true(all(g$icc > 0 & g$icc < 1))
})


test_that("gen_icc() produces the two-way design it claims", {
  .needs_pkg()
  set.seed(1)
  cond <- data.frame(n = 400, k = 4, icc = 0.6, rater_share = 0.5)
  d <- gen_icc(cond, nrep = 300)

  expect_equal(nrow(d), 300L)
  expect_true(all(c("icc", "n_sample", "n_measurements",
                    "theta_pop", "theta_sample") %in% names(d)))
  expect_true(all(d$n_sample == 400)); expect_true(all(d$n_measurements == 4))

  # the target is on the analysis scale the package reports, not the raw ICC
  expect_equal(unique(d$theta_pop), log(1 - 0.6), tolerance = 1e-12)

  # THE LOAD-BEARING CHECK: the ICC estimator must recover the requested value.
  # If it did not, the coverage table would be measuring estimator bias rather
  # than the standard error's behaviour, and the whole study would be wrong in a
  # way no amount of replication would reveal.
  expect_equal(mean(d$icc), 0.6, tolerance = 0.03)

  # and the estimates are not degenerate
  expect_true(all(is.finite(d$icc)))
  expect_true(all(abs(d$icc) < 1))
})


test_that("rater_share really controls the rater variance", {
  .needs_pkg()
  # At rater_share = 0 the design is one-way and the package SE is exact, so the
  # ICC estimator should be no more variable than at share > 0. If this ordering
  # broke, the control cell would stop being a control.
  set.seed(2)
  lo <- gen_icc(data.frame(n = 60, k = 2, icc = 0.8, rater_share = 0.0), 400)
  hi <- gen_icc(data.frame(n = 60, k = 2, icc = 0.8, rater_share = 0.8), 400)
  expect_lt(stats::sd(lo$icc), stats::sd(hi$icc))
  # both still target the same ICC
  expect_equal(mean(lo$icc), 0.8, tolerance = 0.05)
  expect_equal(mean(hi$icc), 0.8, tolerance = 0.08)
})


test_that("estimate_icc() returns the long frame the runner expects", {
  .needs_pkg()
  set.seed(3)
  cond <- data.frame(n = 50, k = 2, icc = 0.8, rater_share = 0.5)
  d <- gen_icc(cond, nrep = 40)
  e <- estimate_icc(d, cond)

  expect_true(all(c("method", "est", "se", "ci_lo", "ci_up") %in% names(e)))
  expect_equal(nrow(e), 80L)                       # 2 methods x 40 replications
  expect_setequal(unique(e$method), c("agreement", "consistency"))
  expect_true(all(is.finite(e$est)))
  expect_true(all(is.finite(e$se)))
})


test_that("icc_type is numerically inert -- the fact that killed remedy (a)", {
  .needs_pkg()
  # Flipping the default icc_type from "agreement" to "consistency" was one of the
  # two remedies roadmap 1.1 offered for the anti-conservative SE. It cannot be
  # one: both types share a single SE formula, so the change only relabels rows.
  # The difference is exactly 0, not merely small, and showing that beats
  # asserting it in a comment.
  z <- icc_type_is_inert()
  expect_gt(nrow(z), 10L)
  expect_equal(max(z$es_diff), 0)
  expect_equal(max(z$se_diff), 0)
})


test_that("icc_se_by_k() reproduces V42's premise", {
  .needs_pkg()
  # V42 exists because k moves the ICC SE far more than the 1.039x item-count
  # case V37 was written for -- but the figure is rho-dependent, which is why the
  # roadmap's single "1.435x" could not be reproduced as one number.
  z <- icc_se_by_k()
  expect_true(all(z$ratio_max_min > 1))
  # the SE falls as k rises, at every rho
  expect_true(all(z$se_at_k_max < z$se_at_k_min))
  # rho-dependence is the point: the ratio is not constant across rho
  expect_gt(diff(range(z$ratio_max_min)), 0.5)
  # and it comfortably exceeds the alpha item-count case V37 handles
  expect_true(all(z$ratio_max_min > 1.039))
})


test_that("v36_false_positive_rate() measures FALSE positives only", {
  .needs_pkg()
  # Every alpha is drawn independently, so there is no induction anywhere and
  # every flag is a false positive by construction. Small sweep: the published
  # table comes from run_10's own reps, not from here.
  z <- v36_false_positive_rate(k = c(10, 30), dp = c(2, 4),
                               n_round = TRUE, reps = 40L)
  expect_equal(nrow(z), 4L)
  expect_true(all(z$pool_fp >= 0 & z$pool_fp <= 1))
  expect_true(all(z$row_fp >= 0 & z$row_fp <= 1))
  # the pool rate rises with pool size at fixed precision -- the birthday
  # argument that refutes tightening the decimal gate
  for (d in unique(z$dp)) {
    s <- z[z$dp == d, ]
    s <- s[order(s$k), ]
    expect_gte(s$pool_fp[nrow(s)], s$pool_fp[1])
  }
})


test_that("the study documents why the Fisher-TF comparison is absent", {
  # Roadmap 3.1 also asks for the low-ICC Bonett-vs-Fisher-TF pooling comparison
  # (item 1.5). icc_to_es = "fisher_tf" is NOT IMPLEMENTED, so simulating it would
  # mean re-typing the candidate formula inline and scoring the package against a
  # private copy of a method it does not ship -- the practice this programme was
  # rebuilt to stop. The omission has to be explained in the file, or the next
  # reader records it as an oversight and adds exactly that.
  f <- file.path(.sim_root, "studies", "10_reliability.R")
  expect_true(file.exists(f))
  # collapse first: the explanation wraps across lines, so a per-line grep would
  # pass or fail on where the text happens to break
  src <- paste(readLines(f, warn = FALSE), collapse = " ")
  expect_true(grepl("fisher_tf", src, fixed = TRUE))
  expect_true(grepl("IS NOT IMPLEMENTED", src, fixed = TRUE))
  expect_true(grepl("1.5", src, fixed = TRUE))
})
