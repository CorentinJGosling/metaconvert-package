## Red-half tests for the 2026-08-28 audit, Tier-2 flag findings.
##
## Every test_that() below asserts the CORRECT behaviour and therefore FAILS at the
## time of writing. Do not "fix" these tests: fix R/internal_flags.R (and, for the
## last one, R/internal_generate_df.R).

# ---- helpers ---------------------------------------------------------------

# convert_df() + summary() both write progress to the console; keep the suite quiet.
audit_summary <- function(dat, measure, ...) {
  out <- NULL
  utils::capture.output(
    suppressMessages({
      res <- convert_df(dat, measure = measure, verbose = FALSE)
      out <- summary(res, flags = TRUE, guidance = FALSE, ...)
    })
  )
  as.data.frame(out)
}

# .flag_* helpers return a list (one character vector per row).
audit_fired <- function(flags, pattern) {
  vapply(flags, function(x) any(grepl(pattern, x, fixed = TRUE)), logical(1))
}

# ---------------------------------------------------------------------------
# AUDIT line 775 - B3 ("NNT magnitude < 1") is not rate-aware: it applies the
# risk-based bound |NNT| >= 1 to the person-time NNT that es_from_cases_time()
# itself computes. Per the verifier correction the same route-blindness hits
# B6/B6b, which judge the rate-scale rd (an incidence rate difference) against
# the risk-scale [-1, 1] bounds.
# ---------------------------------------------------------------------------
test_that("AUDIT-775: B3/B6 stand down on the rate-scale outputs of the cases_time route", {
  df <- data.frame(
    study_id = "B", n_cases_exp = 500, n_cases_nexp = 10,
    time_exp = 100, time_nexp = 100, baseline_rate = 0.1
  )

  nnt <- audit_summary(df, "nnt")
  # IRR = (500/100)/(10/100) = 50; IRD = baseline_rate * (1 - IRR) = 0.1 * (-49)
  # = -4.9 events per person-time; NNT = 1 / IRD = -0.2040816 person-time units.
  expect_identical(nnt$info_used_crude, "cases_time")
  expect_equal(nnt$es_crude, -1 / 4.9, tolerance = 1e-8)
  # A Mayne (2006) person-time NNT is a DURATION; |NNT| < 1 is entirely possible.
  expect_false(grepl("NNT magnitude < 1", nnt$flags_crude, fixed = TRUE))

  # Benefit side of the same route, baseline_rate auto-computed from the counts:
  # rate_nexp = 60/10 = 6/py, rate_exp = 30/10 = 3/py, IRD = 3/py, NNT = 1/3 py.
  df_benefit <- data.frame(
    study_id = "D", n_cases_exp = 30, n_cases_nexp = 60,
    time_exp = 10, time_nexp = 10
  )
  nnt_b <- audit_summary(df_benefit, "nnt")
  expect_equal(nnt_b$es_crude, 1 / 3, tolerance = 1e-8)
  expect_false(grepl("NNT magnitude < 1", nnt_b$flags_crude, fixed = TRUE))

  # Verifier correction: the same route writes the incidence rate difference into
  # the rd columns, so B6/B6b must not apply the risk-scale [-1, 1] bound either.
  rd <- audit_summary(df, "rd")
  expect_identical(rd$info_used_crude, "cases_time")
  expect_equal(rd$es_crude, -4.9, tolerance = 1e-8)
  expect_false(grepl("RD outside [-1, 1]", rd$flags_crude, fixed = TRUE))
  expect_false(grepl("RD CI bound outside [-1, 1]", rd$flags_crude, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# AUDIT line 445 - B5 ("minimum possible NNT") applies the benefit-side floor
# 1/baseline_risk to both signs. On the harm side the attainable floor is
# 1/(1 - baseline_risk), and for the person-time branch there is no floor at all
# (the incidence rate ratio is unbounded above).
# ---------------------------------------------------------------------------
test_that("AUDIT-445: B5 uses the harm-side bound 1/(1-BR), not 1/BR, for a negative NNT", {
  # Fully self-consistent 2x2: control 5/100, exposed 90/100, baseline_risk = 0.05.
  df <- data.frame(
    study_id = "A", n_cases_exp = 90, n_controls_exp = 10,
    n_cases_nexp = 5, n_controls_nexp = 95, baseline_risk = 0.05
  )
  nnt <- audit_summary(df, "nnt")
  # RD = 0.05 - 0.90 = -0.85  =>  NNH = 1 / -0.85 = -1.176471, computed by the
  # package from the counts. Harm-side floor = 1/(1 - 0.05) = 1.0526316, and
  # 1.176471 > 1.052632, so the value is attainable and must not be [INVALID].
  expect_equal(nnt$es_crude, -1 / 0.85, tolerance = 1e-8)
  expect_gt(abs(nnt$es_crude), 1 / (1 - 0.05))
  expect_false(grepl("Minimum possible NNT", nnt$flags_crude, fixed = TRUE))

  # Person-time harm side: IRR is unbounded above, so |IRD| is unbounded and no
  # positive lower bound on |NNT| exists.
  df2 <- data.frame(
    study_id = "B", n_cases_exp = 500, n_cases_nexp = 10,
    time_exp = 100, time_nexp = 100, baseline_rate = 0.1
  )
  nnt2 <- audit_summary(df2, "nnt")
  expect_lt(nnt2$es_crude, 0)
  expect_false(grepl("Minimum possible person-time NNT", nnt2$flags_crude, fixed = TRUE))

  # Regression guard: the benefit-side bound must survive the fix. Control 30/100
  # vs exposed 10/100 gives RD = 0.20, NNT = 5, below the declared floor 1/0.05 = 20.
  df3 <- data.frame(
    study_id = "C", n_cases_exp = 10, n_controls_exp = 90,
    n_cases_nexp = 30, n_controls_nexp = 70, baseline_risk = 0.05
  )
  nnt3 <- audit_summary(df3, "nnt")
  expect_equal(nnt3$es_crude, 5, tolerance = 1e-8)
  expect_true(grepl("Minimum possible NNT", nnt3$flags_crude, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# AUDIT line 555 - A6 treats "r" as a symmetric-t-interval measure, but the
# package builds the r CI by back-transforming a symmetric interval on the
# variance-stabilised z scale (smd_to_cor) or on the Fisher z scale (tetrachoric
# 2x2). The reported CI is correct; the flag is the error.
# ---------------------------------------------------------------------------
test_that("AUDIT-555: A6 does not flag a package-built r CI (back-transformed z interval)", {
  e <- es_from_means_sd(
    n_exp = 13, n_nexp = 13,
    mean_exp = 5, mean_sd_exp = 2, mean_nexp = 3, mean_sd_nexp = 2
  )

  # Reproduce the package's own construction: back-transform of the symmetric
  # interval on the variance-stabilised z scale (Viechtbauer, the default).
  p <- 0.5
  a <- sqrt(stats::dnorm(stats::qnorm(p))) / (p * (1 - p))^(1 / 4)
  bt <- function(zz) (1 / a) * ((exp(2 * zz / a) - 1) / (exp(2 * zz / a) + 1))
  z_lo <- e$z - stats::qnorm(.975) * e$z_se
  z_up <- e$z + stats::qnorm(.975) * e$z_se
  expect_equal(e$r_ci_lo, bt(z_lo), tolerance = 1e-10)
  expect_equal(e$r_ci_up, bt(z_up), tolerance = 1e-10)

  fl <- metaConvert:::.flag_numeric_integrity(
    e$r, e$r_se, e$r_ci_lo, e$r_ci_up,
    info_used = "means_sd", measure = "r", exp = FALSE,
    n_total = 26, n_exp = 13, n_nexp = 13
  )
  expect_false(any(grepl("CI width inconsistent with SE", fl[[1]], fixed = TRUE)))

  # Second back-transform in play: the tetrachoric route builds the r CI as
  # tanh() of the Fisher z interval, which is likewise asymmetric about r.
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    e2 <- es_from_2x2(
      n_cases_exp = 10, n_controls_exp = 10,
      n_cases_nexp = 5, n_controls_nexp = 15
    )
    expect_equal(e2$r_ci_lo, tanh(e2$z - stats::qnorm(.975) * e2$z_se), tolerance = 1e-10)
    expect_equal(e2$r_ci_up, tanh(e2$z + stats::qnorm(.975) * e2$z_se), tolerance = 1e-10)
    fl2 <- metaConvert:::.flag_numeric_integrity(
      e2$r, e2$r_se, e2$r_ci_lo, e2$r_ci_up,
      info_used = "2x2", measure = "r", exp = FALSE,
      n_total = 40, n_exp = 20, n_nexp = 20
    )
    expect_false(any(grepl("CI width inconsistent with SE", fl2[[1]], fixed = TRUE)))
  }

  # Regression guards: A6 must keep working where the r CI is genuinely wrong,
  # and must keep passing the route that really does build r +/- qt * se.
  fl_bad <- metaConvert:::.flag_numeric_integrity(
    e$r, e$r_se,
    e$r - (e$r - e$r_ci_lo) / 2, e$r + (e$r_ci_up - e$r) / 2,
    info_used = "means_sd", measure = "r", exp = FALSE,
    n_total = 26, n_exp = 13, n_nexp = 13
  )
  expect_true(any(grepl("CI width inconsistent with SE", fl_bad[[1]], fixed = TRUE)))

  pr <- es_from_pearson_r(pearson_r = 0.5, n_sample = 26)
  expect_equal(pr$r_ci_up - pr$r_ci_lo, 2 * stats::qt(.975, 24) * pr$r_se, tolerance = 1e-10)
  fl_pr <- metaConvert:::.flag_numeric_integrity(
    pr$r, pr$r_se, pr$r_ci_lo, pr$r_ci_up,
    info_used = "pearson_r", measure = "r", exp = FALSE,
    n_total = 26, n_exp = 13, n_nexp = 13
  )
  expect_false(any(grepl("CI width inconsistent with SE", fl_pr[[1]], fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# AUDIT line 566 - `skip_a6 <- identical(measure, "nnt") || exp` disables A6 for
# every exp-scale measure, so a user-entered SE given on the natural scale
# instead of the log scale has no Tier-2 detector. The check is well defined on
# the log scale: log(ci_up) - log(ci_lo) vs 2 * qnorm(.975) * se.
# ---------------------------------------------------------------------------
test_that("AUDIT-566: A6 runs on the log scale for exp-scale measures (or/rr/irr/hr)", {
  z <- stats::qnorm(.975)
  or <- 2.0
  se_log <- 0.20
  ci_lo <- or * exp(-z * se_log)
  ci_up <- or * exp(+z * se_log)
  # Observed log-scale width = 2 * z * 0.20 = 0.7839856.
  expect_equal(log(ci_up) - log(ci_lo), 2 * z * se_log, tolerance = 1e-10)

  # Wrong-scale SE: 0.4 is the NATURAL-scale SE (2.0 * 0.20). A6 should expect a
  # log width of 2 * z * 0.4 = 1.5679712 and find 0.7839856 - a 50% discrepancy,
  # far outside A6's own tolerance.
  fl_bad <- metaConvert:::.flag_numeric_integrity(
    or, 0.4, ci_lo, ci_up,
    info_used = "user_input_crude", measure = "logor", exp = TRUE,
    n_total = 100, n_exp = 50, n_nexp = 50
  )
  expect_true(any(grepl("CI width inconsistent with SE", fl_bad[[1]], fixed = TRUE)))

  # The correctly-scaled row must stay silent.
  fl_ok <- metaConvert:::.flag_numeric_integrity(
    or, se_log, ci_lo, ci_up,
    info_used = "user_input_crude", measure = "logor", exp = TRUE,
    n_total = 100, n_exp = 50, n_nexp = 50
  )
  expect_false(any(grepl("CI width inconsistent with SE", fl_ok[[1]], fixed = TRUE)))

  # End to end: six OR rows, row 1 declaring its SE on the natural scale.
  d <- data.frame(
    study_id = paste0("s", 1:6),
    user_es_crude = c(2.00, 1.90, 2.10, 2.00, 1.95, 2.05),
    user_es_original_measure_crude = "or",
    user_es_target_measure_crude = "or",
    stringsAsFactors = FALSE
  )
  d$user_ci_lo_crude <- d$user_es_crude * exp(-z * se_log)
  d$user_ci_up_crude <- d$user_es_crude * exp(+z * se_log)
  d$user_se_crude <- c(2.00 * se_log, rep(se_log, 5))
  res <- audit_summary(d, "or")
  expect_true(grepl("CI width inconsistent with SE", res$flags_crude[1], fixed = TRUE))
  expect_false(any(grepl("CI width inconsistent with SE", res$flags_crude[-1], fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# AUDIT line 577 - D2's SE normalisation for measure = "logrr" hands the RISK
# ratio (and Var(logRR)) to an OR-based 2x2 solver, so the normalised SEs no
# longer sit on a common ~1 scale and an ordinary small low-event-rate trial is
# flagged an SE outlier. The identical tables run as "or" raise nothing.
# ---------------------------------------------------------------------------
test_that("AUDIT-577: D2 does not raise an SE outlier on a correct 2x2 risk-ratio pool", {
  a  <- c(40, 10, 80, 20,  5, 50, 30, 12)
  b  <- c(20,  5, 60, 10,  2, 25, 15,  6)
  n1 <- c(100, 100, 100, 200, 50, 100, 100, 120)
  n2 <- n1
  d <- data.frame(
    study_id = paste0("s", seq_along(a)),
    n_cases_exp = a, n_cases_nexp = b,
    n_controls_exp = n1 - a, n_controls_nexp = n2 - b
  )

  rr <- audit_summary(d, "rr")
  # Study 5 is a=5/50 vs b=2/50: RR = 2.5 and
  # SE(logRR) = sqrt(1/5 - 1/50 + 1/2 - 1/50) = 0.8124038 - a correct table with
  # a correct SE, merely a small trial with a low event rate.
  expect_equal(rr$es_crude[5], 2.5, tolerance = 1e-8)
  expect_equal(rr$se_crude[5], sqrt(1 / 5 - 1 / 50 + 1 / 2 - 1 / 50), tolerance = 1e-8)
  expect_false(any(grepl("SE outlier", rr$flags_crude, fixed = TRUE)))

  # Control: the same tables as odds ratios, where the normalisation is exact.
  or <- audit_summary(d, "or")
  expect_false(any(grepl("SE outlier", or$flags_crude, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# AUDIT line 1009 - D1's `min_dev` floor has no entry for "rp"/"zp", so they fall
# to the 0 default reserved for md/nnt/prop and get outlier flags at deviations
# the same switch calls non-diagnostic for "r"/"z".
# ---------------------------------------------------------------------------
test_that("AUDIT-1009: D1's outlier_min_deviation floor covers rp and zp", {
  opts <- metaConvert:::.default_flag_options()
  es <- c(0.50, 0.50, 0.51, 0.51, 0.52, 0.52, 0.53, 0.53, 0.54, 0.70)
  # median = 0.52, IQR = 0.02, so row 10 clears the 3*IQR rule; its deviation
  # from the median is 0.18, well under the 0.3 floor "r" gets and the 0.5 "z" gets.
  run <- function(measure) {
    metaConvert:::.flag_cross_row_outliers(
      es, rep(0.2, 10), opts,
      info_used = rep("means_sd", 10), measure = measure,
      n_total = rep(100, 10), n_exp = rep(50, 10), n_nexp = rep(50, 10)
    )
  }
  # Reference behaviour of the families rp/zp belong to.
  expect_false(any(audit_fired(run("r"), "ES outlier")))
  expect_false(any(audit_fired(run("z"), "ES outlier")))
  # rp is bounded [-1, 1] exactly as r is; zp is on the same Fisher scale as z.
  expect_false(any(audit_fired(run("rp"), "ES outlier")))
  expect_false(any(audit_fired(run("zp"), "ES outlier")))
  # md / mdw are correctly floorless (arbitrary outcome units) and must stay so.
  expect_true(any(audit_fired(run("md"), "ES outlier")))
})

# ---------------------------------------------------------------------------
# AUDIT line 1035 - the same `min_dev` switch omits "dw"/"gw", which are
# standardised SMDs on exactly the same scale as d/g and should inherit the
# 1.0 SD floor.
# ---------------------------------------------------------------------------
test_that("AUDIT-1035: D1's outlier_min_deviation floor covers dw and gw", {
  opts <- metaConvert:::.default_flag_options()
  es <- c(0.50, 0.50, 0.51, 0.51, 0.52, 0.52, 0.53, 0.53, 0.54, 0.70)
  run <- function(measure) {
    metaConvert:::.flag_cross_row_outliers(
      es, rep(0.2, 10), opts,
      info_used = rep("means_sd", 10), measure = measure,
      n_total = rep(100, 10), n_exp = rep(50, 10), n_nexp = rep(50, 10)
    )
  }
  # Max deviation from the median is 0.18 SD, far under the 1.0 SD floor d/g get.
  expect_false(any(audit_fired(run("d"), "ES outlier")))
  expect_false(any(audit_fired(run("g"), "ES outlier")))
  expect_false(any(audit_fired(run("dw"), "ES outlier")))
  expect_false(any(audit_fired(run("gw"), "ES outlier")))
})

# ---------------------------------------------------------------------------
# AUDIT line 1048 - D3/V27 defines a row's spread as the pooled arm SD when both
# arm SDs are present and as se * sqrt(n_eff) otherwise. For a two-group md,
# SE_MD = SD * sqrt(1/n1 + 1/n2), so se * sqrt(N) = 2 * SD at balanced arms: the
# two definitions differ by a factor of 2 inside one pool. The implied SD must be
# se / sqrt(1/n_exp + 1/n_nexp).
# ---------------------------------------------------------------------------
test_that("AUDIT-1048: D3 puts arm-SD rows and SE-only rows on one spread scale", {
  opts <- metaConvert:::.default_flag_options()
  n1 <- 50; n2 <- 50
  run <- function(sd_true, has_arm_sd) {
    se <- sd_true * sqrt(1 / n1 + 1 / n2)
    metaConvert:::.flag_cross_row_outliers(
      rep(2, length(sd_true)), se, opts,
      info_used = rep("md_sd", length(sd_true)), measure = "md",
      n_total = rep(n1 + n2, length(sd_true)),
      n_exp = rep(n1, length(sd_true)), n_nexp = rep(n2, length(sd_true)),
      sd_exp = ifelse(has_arm_sd, sd_true, NA_real_),
      sd_nexp = ifelse(has_arm_sd, sd_true, NA_real_)
    )
  }

  # Ten studies; study 10 is genuinely lower-variance (3.5 vs a pool near 10),
  # i.e. only 2.9x below the pool - under the K = 5 gate on any single scale.
  sd_true <- c(10, 10, 11, 9, 10, 12, 9, 10, 10, 3.5)
  none <- rep(FALSE, 10)
  all_of <- rep(TRUE, 10)
  only10 <- c(rep(FALSE, 9), TRUE)

  expect_false(any(audit_fired(run(sd_true, none), "SD outlier")))
  expect_false(any(audit_fired(run(sd_true, all_of), "SD outlier")))
  # Mixed completeness - the normal case for an extraction sheet. The study is
  # identical in all three runs; only the presence of transcribed arm SDs changes.
  expect_false(any(audit_fired(run(sd_true, only10), "SD outlier")))

  # When D3 legitimately fires, it must report the pooled-SD scale. Here rows 1-9
  # have SD = 10 and row 10 has SD = 1 (a genuine SE-as-SD error, 10x too small);
  # no row reports arm SDs, so every spread comes from the SE. The implied SD is
  # se / sqrt(1/50 + 1/50) = 2 / 0.2 = 10 for rows 1-9 and 1 for row 10 - not the
  # doubled 20 / 2 the current se * sqrt(N) rule produces.
  msg <- run(c(rep(10, 9), 1), none)[[10]]
  hit <- msg[grepl("SD outlier", msg, fixed = TRUE)]
  expect_length(hit, 1L)
  spread_val <- suppressWarnings(as.numeric(sub(".*spread = ([0-9.]+) is.*", "\\1", hit)))
  median_val <- suppressWarnings(as.numeric(sub(".*pool median \\(([0-9.]+)\\).*", "\\1", hit)))
  expect_equal(spread_val, 1, tolerance = 1e-6)
  expect_equal(median_val, 10, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# AUDIT line 1061 - under exp = TRUE, E1 approximates the log-scale dispersion as
# dispersion / es[i], using the hierarchy-SELECTED estimate as the delta-method
# centre instead of the median the deviation was measured around. The printed
# number is not the quantity the message names, and whether E1 fires depends on
# which estimate was selected.
# ---------------------------------------------------------------------------
test_that("AUDIT-1061: E1's log-scale dispersion is exact and selection-invariant under exp=TRUE", {
  fe <- metaConvert:::.flag_internal_consistency
  ds <- metaConvert:::.dispersion_stat
  opts <- metaConvert:::.default_flag_options()

  run <- function(v, selected) {
    fe(2, ds(v), NA, max(v) - min(v), opts, "logor",
       min_info = "2x2", max_info = "means_sd",
       min_es = min(v), max_es = max(v), exp = TRUE, es = selected)[[1]]
  }

  # Case 1 - two estimates, OR 1.2 and OR 4.0. The true log-scale dispersion is
  # max|log(v) - median(log(v))| = (log(4) - log(1.2)) / 2 = 0.6019864, BELOW the
  # threshold of 1, so E1 must not fire whichever estimate the hierarchy picked.
  v1 <- c(1.2, 4.0)
  expect_equal(ds(log(v1)), 0.6019864, tolerance = 1e-6)
  expect_false(any(grepl("High dispersion", run(v1, 1.2), fixed = TRUE)))
  expect_false(any(grepl("High dispersion", run(v1, 4.0), fixed = TRUE)))

  # Case 2 - OR 1 and OR 10. True log-scale dispersion = log(10)/2 = 1.1512925,
  # ABOVE the threshold, so E1 must fire for either selection and must print the
  # true value, not dispersion/es (4.5 when 1 is selected, 0.45 when 10 is).
  v2 <- c(1, 10)
  expect_equal(ds(log(v2)), 1.1512925, tolerance = 1e-6)
  for (sel in v2) {
    msg <- run(v2, sel)
    expect_true(any(grepl("High dispersion", msg, fixed = TRUE)))
    expect_true(any(grepl("= 1.151 (log scale)", msg, fixed = TRUE)))
  }
})

# ---------------------------------------------------------------------------
# AUDIT line 1074 - .generate_df() treats col_min == col_max as a tie and
# re-derives both indices from the CI-lower matrix; when the tied routes have NA
# CI bounds (systematic for an NNT whose risk-difference CI crosses zero, per the
# Altman 1998 disjoint rule) which.min(all-NA)[1] is NA and every min_*/max_*
# diagnostic plus diff_min_max is blanked.
# ---------------------------------------------------------------------------
test_that("AUDIT-1074: min/max diagnostics survive tied routes with NA CI bounds", {
  mk <- function(a, b, c, d) {
    data.frame(
      study_id = "s", n_cases_exp = a, n_controls_exp = b,
      n_cases_nexp = c, n_controls_nexp = d, n_exp = a + b, n_nexp = c + d
    )
  }

  # exposed 30/100 vs control 20/100: RD = 0.20 - 0.30 = -0.10, NNT = -10, and the
  # RD CI crosses zero so both 2x2 routes return NA NNT bounds.
  crossing <- audit_summary(mk(30, 70, 20, 80), "nnt")
  expect_equal(crossing$n_estimations_crude, 2)
  expect_equal(crossing$es_crude, -10, tolerance = 1e-8)
  # Both tied routes return exactly -10, so the diagnostics are well defined.
  expect_equal(crossing$min_es_value_crude, -10, tolerance = 1e-8)
  expect_equal(crossing$max_es_value_crude, -10, tolerance = 1e-8)
  expect_equal(crossing$diff_min_max_crude, 0, tolerance = 1e-8)
  expect_false(is.na(crossing$min_info_crude))
  expect_false(is.na(crossing$max_info_crude))
  expect_true(crossing$min_info_crude %in% c("2x2", "2x2_sum"))
  expect_true(crossing$max_info_crude %in% c("2x2", "2x2_sum"))

  # Control: the same shape with a CI clear of zero already works today.
  clear <- audit_summary(mk(60, 40, 20, 80), "nnt")
  expect_equal(clear$es_crude, -2.5, tolerance = 1e-8)
  expect_equal(clear$diff_min_max_crude, 0, tolerance = 1e-8)
  expect_identical(clear$min_info_crude, "2x2")
})
