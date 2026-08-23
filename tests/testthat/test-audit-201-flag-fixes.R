# =============================================================================
# Audit 2.0.1 flag fixes: findings #9, #41, #42, #43.
#
# All four are diagnostic-only -- no effect size, standard error or confidence
# interval changes anywhere in this file. What changes is which rows get flagged
# and what the flag says, which is exactly what a reviewer acts on.
#
#  #9  A6 fired a false [DISCORDANT] CI-width flag on EVERY measure = "rp" row at
#      small n, because "rp" was missing from qt_measures and the check therefore
#      compared a t-based interval against the Wald-z width.
#  #41 A study_id containing "; " fragmented the V23 message and mis-routed the
#      leading fragment into both the crude and adjusted scopes.
#  #42 B1b and C4 were gated on measure == "r", so a partial correlation never
#      got a bounds or magnitude check at all.
#  #43 In a 1-vs-1 direction split, BOTH rows were told they opposed "the pool
#      majority" -- a self-contradictory pair of flags describing a pool that has
#      no majority.
# =============================================================================

fires <- function(f, pat) vapply(f, function(x) grepl(pat, x, fixed = TRUE), logical(1))

rp_run <- function(n_sample, n_covariates, linreg_t = 2.5) {
  summary(
    convert_df(data.frame(linreg_t = linreg_t, n_sample = n_sample,
                          n_covariates = n_covariates),
               measure = "rp", verbose = FALSE),
    flags = TRUE)
}


# --- #9 ----------------------------------------------------------------------

test_that("A6 does not fire on partial correlations built with the regression df", {
  # es_from_linreg_t() builds rp on qt(.975, n_sample - n_covariates - 2)
  # (R/es_from_REGRESSION.R:249, :268-269). Before the fix A6 expected the Wald-z
  # width and fired on every one of these rows.
  r <- rp_run(c(8, 10, 12, 15, 20, 23), n_covariates = 1)
  expect_false(any(fires(r$flags_crude, "CI width inconsistent")))
})

test_that("the A6 rp fix subtracts the covariates, not just the two", {
  # The crossover scales with n_covariates: measured last-firing N was 23 at k = 1,
  # 27 at k = 3 and 37 at k = 10. Adding "rp" to qt_measures WITHOUT subtracting
  # the covariate count would leave k >= 3 still false-firing.
  for (k in c(3, 10)) {
    r <- rp_run(c(15, 20, 27, 30, 37, 40), n_covariates = k)
    expect_false(any(fires(r$flags_crude, "CI width inconsistent")),
                 info = paste("n_covariates =", k))
  }
})

test_that("zp stays out of qt_measures: its interval really is built with qnorm", {
  # R/es_from_REGRESSION.R:276-277. Adding "zp" would introduce the mirror-image
  # false positive, so this pins the exclusion rather than assuming it.
  r <- summary(convert_df(data.frame(linreg_t = 2.5, n_sample = c(8, 10, 20),
                                     n_covariates = 1),
                          measure = "zp", verbose = FALSE), flags = TRUE)
  expect_false(any(fires(r$flags_crude, "CI width inconsistent")))
})


# --- #42 ---------------------------------------------------------------------

test_that("B1b flags a partial-correlation CI that escapes [-1, 1]", {
  r <- rp_run(c(8, 10, 12), n_covariates = 1)
  expect_true(all(r$es_ci_up_crude > 1))                       # the precondition
  expect_true(all(abs(r$es_crude) <= 1))                       # point estimate valid
  expect_true(all(fires(r$flags_crude, "escapes the parameter space")))
  # the message names the measure correctly
  expect_true(all(fires(r$flags_crude, "while partial r is valid")))
})

test_that("C4 flags a high partial correlation", {
  r <- rp_run(10, n_covariates = 1, linreg_t = 9)
  expect_gt(abs(r$es_crude[1]), 0.95)
  expect_true(fires(r$flags_crude, "High correlation: |partial r| ="))
})

test_that("B1b reports BOTH bounds when both escape, not only the upper one", {
  r <- summary(convert_df(data.frame(user_es_crude = 0, user_se_crude = 0.6,
                                     user_es_original_measure_crude = "r",
                                     n_sample = 50),
                          measure = "r", verbose = FALSE), flags = TRUE)
  expect_lt(r$es_ci_lo_crude[1], -1)
  expect_gt(r$es_ci_up_crude[1], 1)
  expect_true(fires(r$flags_crude, "lower bound is below -1"))
  expect_true(fires(r$flags_crude, "upper bound exceeds 1"))
})

test_that("the r-scale checks are unchanged for measure = 'r'", {
  # #42 widens a gate; it must not alter what the gate already did.
  r <- summary(convert_df(data.frame(user_es_crude = 0.98, user_se_crude = 0.05,
                                     user_es_original_measure_crude = "r",
                                     n_sample = 50),
                          measure = "r", verbose = FALSE), flags = TRUE)
  expect_true(fires(r$flags_crude, "escapes the parameter space"))
  expect_true(fires(r$flags_crude, "High correlation: |r| ="))
})


# --- #41 ---------------------------------------------------------------------

test_that("a study_id containing '; ' does not leak a V23 fragment into the adjusted scope", {
  # Two different studies sharing a byte-identical baseline block: V23's case.
  # The leading fragment carries no quoted column name, so .v_flag_matches_scope()
  # returned TRUE for every scope and the truncated text landed in both columns.
  d <- data.frame(
    study_id        = c("Huang; 2017", "Xu 2017"),
    mean_pre_exp    = c(2.20, 2.20), mean_pre_sd_exp  = c(0.83, 0.83),
    mean_pre_nexp   = c(2.19, 2.19), mean_pre_sd_nexp = c(0.81, 0.81),
    mean_exp        = c(1.5, 1.9),   mean_sd_exp      = c(0.70, 0.75),
    mean_nexp       = c(2.0, 2.1),   mean_sd_nexp     = c(0.80, 0.82),
    n_exp           = c(40, 50),     n_nexp           = c(40, 50))
  s <- summary(convert_df(d, measure = "g", verbose = FALSE), flags = TRUE)

  expect_true(fires(s$flags_crude[2], "Identical baseline summary statistics"))
  expect_false(fires(s$flags_adjusted[2], "Identical baseline"))
  expect_equal(s$flags_crude[2], gsub("Huang; 2017", "Huang, 2017", s$flags_crude[2],
                                      fixed = TRUE))
})


# --- #43 ---------------------------------------------------------------------

test_that("a 1-vs-1 direction split is not described as opposing a majority", {
  mk <- function(id, out, d) {
    data.frame(study_id = id, outcome = out, n_exp = 40, n_nexp = 40,
               mean_exp = d, mean_sd_exp = 1, mean_nexp = 0, mean_sd_nexp = 1)
  }
  rows <- do.call(rbind, c(
    lapply(1:9, function(k) rbind(mk(paste0("S", k, "a"), paste0("O", k),  0.6746),
                                  mk(paste0("S", k, "b"), paste0("O", k),  0.6746))),
    list(rbind(mk("S10a", "O10", 0.6746), mk("S10b", "O10", -0.6746)))))

  s <- summary(convert_df(rows, measure = "g", verbose = FALSE), flags = TRUE,
               flag_options = list(flag_group = "outcome"))
  hit <- grep("Direction", s$flags_crude)

  expect_length(hit, 2)                                     # both sides of the tie
  # both still warrant [UNUSUAL]: an evenly split pool is exactly as suspicious
  expect_true(all(fires(s$flags_crude[hit], "[UNUSUAL] Direction conflict")))
  expect_true(all(fires(s$flags_crude[hit], "the pool is evenly split")))
  expect_true(all(fires(s$flags_crude[hit], "with no majority direction")))
  # ...and neither is told it opposes a majority that does not exist
  expect_false(any(fires(s$flags_crude[hit], "opposite to the pool majority")))
})

test_that("a genuine majority still gets the minority/majority wording", {
  mk <- function(id, d) {
    data.frame(study_id = id, n_exp = 40, n_nexp = 40,
               mean_exp = d, mean_sd_exp = 1, mean_nexp = 0, mean_sd_nexp = 1)
  }
  rows <- rbind(mk("A", 0.6746), mk("B", 0.6746), mk("C", -0.6746), mk("D", -0.6746),
                mk("E", 0.6746), mk("F", 0.6746), mk("G", 0.6746))
  s <- summary(convert_df(rows, measure = "g", verbose = FALSE), flags = TRUE)
  hit <- grep("Direction", s$flags_crude)
  skip_if(length(hit) == 0, "direction-conflict thresholds not met in this pool")
  expect_true(any(fires(s$flags_crude[hit], "opposite to the pool majority")))
  expect_false(any(fires(s$flags_crude[hit], "the pool is evenly split")))
})


# --- rp follow-up to #42: the F family and the E1/E3 thresholds ----------------
#
# Same defect family as #42 (a gate that names 'r' and forgets 'rp'), found while
# fixing it and measured before being changed. Two distinct failures:
#
#  * .flag_se_sample_size() gated on a `standardized` vector that omitted rp/zp, so
#    it returned EARLY -- F1 and F2 never ran on a partial correlation at all.
#  * had it run, the switch fall-through would have given rp the ratio-measure
#    defaults: floor 0.1/sqrt(N) (flat, so it does not shrink with (1 - rp^2)) and
#    ceiling max(2, 8/sqrt(N)). Measured at N = 200, q = 2, rp = 0.95: the flat
#    floor 0.00707 sits ABOVE the true SE 0.00696, i.e. a false positive on a valid
#    row; and an SD-entered-as-SE error at SE 0.85-0.97 passes a ceiling of 2.0.

test_that('the F family runs on rp and zp at all', {
  f <- metaConvert:::.flag_se_sample_size
  # pre-fix this returned character(0) for every rp/zp row regardless of SE
  expect_length(f(es = 0.4, se = 0.970, measure = 'rp', n_total = 20)[[1]], 1)
  expect_length(f(es = 0.4, se = 0.970, measure = 'zp', n_total = 20)[[1]], 1)
})

test_that('rp and zp get the r/z SE ceiling, so an SD-as-SE error is caught', {
  f <- metaConvert:::.flag_se_sample_size
  for (N in c(20, 50, 200)) {
    bad <- (1 - 0.16) / sqrt(N - 5) * sqrt(N)          # SE inflated ~sqrt(N)
    ref <- f(es = 0.4, se = bad, measure = 'r',  n_total = N)[[1]]
    expect_true(any(grepl('Implausibly large SE', ref)))   # r catches it...
    for (m in c('rp', 'zp')) {
      got <- f(es = 0.4, se = bad, measure = m, n_total = N)[[1]]
      expect_true(any(grepl('Implausibly large SE', got)),
                  info = paste(m, 'N =', N, 'SE =', round(bad, 3)))
    }
  }
})

test_that('the rp SE floor scales with (1 - rp^2), so a valid high-rp row is silent', {
  f <- metaConvert:::.flag_se_sample_size
  for (N in c(50, 200, 500)) {
    se <- (1 - 0.95^2) / sqrt(N - 5)                   # the true SE(rp)
    got <- f(es = 0.95, se = se, measure = 'rp', n_total = N)[[1]]
    expect_length(got, 0)
    # the unscaled ratio-default floor would have flagged this at N = 200
    if (N == 200) expect_lt(se, 0.1 / sqrt(N))
  }
})

test_that('zp gets no (1 - es^2) scaling: it is on the unbounded Fisher scale', {
  f <- metaConvert:::.flag_se_sample_size
  # a zp of 2.0 is legitimate; scaling by (1 - 4) would make the floor negative
  expect_length(f(es = 2.0, se = 1 / sqrt(200 - 5), measure = 'zp', n_total = 200)[[1]], 0)
})

test_that('E1/E3 use the correlation thresholds for rp, not the SMD fallback', {
  o <- metaConvert:::.default_flag_options()
  expect_equal(o$dispersion_max_r, 0.15); expect_equal(o$diff_max_r, 0.30)
  expect_equal(o$dispersion_max_smd, 0.50); expect_equal(o$diff_max_smd, 1.00)

  # Two rp routes deliberately disagreeing: linreg_t = 2.5 against a
  # linreg_b/linreg_b_se implying t = 15. Dispersion 0.295 and min-max diff 0.59
  # sit BETWEEN the correlation thresholds (0.15 / 0.30) and the SMD fallback
  # (0.50 / 1.00), so this pool fires only if rp is on the r thresholds -- which
  # is precisely the discrimination this fix is about.
  s <- summary(convert_df(data.frame(linreg_t = 2.5, linreg_b = 3.0,
                                     linreg_b_se = 0.2, n_sample = 100,
                                     n_covariates = 1),
                          measure = 'rp', verbose = FALSE), flags = TRUE)
  expect_equal(s$n_estimations_crude[1], 2)
  expect_equal(round(s$diff_min_max_crude[1], 2), 0.59)
  expect_true(fires(s$flags_crude, 'High dispersion across estimates'))
  expect_true(fires(s$flags_crude, 'Large min-max difference'))
  # and the observed values really are below the fallback it used to get
  expect_lt(0.295, o$dispersion_max_smd)
  expect_lt(s$diff_min_max_crude[1], o$diff_max_smd)
})
