# Red-half regression tests for the OR/RR/RD/NNT findings of
# AUDIT-2026-08-28-findings.md. Every test in this file asserts the CORRECT
# behaviour, so each one FAILS until the corresponding defect is repaired.
#
# Where the audit's "Verifier correction" field disagrees with its "Fix" field,
# the correction governs (it is quoted in the comment above each block).

# --------------------------------------------------------------------------
# AUDIT line 365 (BLOCKER) --- .or_to_rr()'s metaumbrella_cases / _exp branches
# read logrr_se off the reconstructed 2x2 table and never consult the user's
# reported logor_se, so the returned RR precision is frozen: a 30x range in the
# reported OR standard error yields a byte-identical logrr_se.
#
# Correct behaviour: RR is a deterministic monotone function of (OR, margins),
# so the SE must be propagated by the delta method,
#     SE(logRR) = |d logRR / d logOR| * logor_se.
#
# The derivative is the ratio of the two MARGINAL standard errors at the
# reconstructed table -- Katz over Woolf -- because SE(logOR) and SE(logRR) describe
# the same product-binomial sampling model:
#     SE(logOR) = sqrt(1/a + 1/b + 1/c + 1/d)                     (Woolf)
#     SE(logRR) = sqrt(1/a - 1/n1 + 1/c - 1/n2)                   (Katz et al. 1978)
#     g = SE(logRR) / SE(logOR)
# The fixed-margins alternative g = (1/a+1/c)/(1/a+1/b+1/c+1/d) is the CONDITIONAL
# derivative: it holds the case margin fixed, but that margin is an observed
# statistic, not a design constant, and unlike for the odds ratio it is not ancillary
# for the risk ratio. A 4e5-replicate product-binomial Monte Carlo on 90/10 vs 60/40
# gives SD(logRR) = 0.088867, matching Katz (0.088192) and not the conditional form
# (0.070767).
#
# This block deliberately does NOT hard-code the derivative: an earlier version pinned
# 0.7630761098, the conditional value, and passed for two reasons that are worth
# recording. Its configuration is a 20%-event table, where the two forms agree to
# 0.26%; and the constant was transcribed from the same derivation the code used, so
# the test could only ever confirm that the code matched itself. The pin below is
# external instead: when the reported logor_se is exactly the one the reconstructed
# table implies, the propagated logrr_se must equal metafor's own RR standard error.
test_that("AUDIT-blocker: or_to_rr metaumbrella_* propagates the reported logor_se into logrr_se", {
  skip_if_not_installed("metafor")

  # A table whose event rates are high enough to separate the two candidate
  # derivatives (they differ by 24% here, against 0.3% at 20% events).
  a <- 90; b <- 10; cc <- 60; d <- 40
  n_exp <- a + b; n_nexp <- cc + d

  mf_or <- metafor::escalc(measure = "OR", ai = a, bi = b, ci = cc, di = d)
  mf_rr <- metafor::escalc(measure = "RR", ai = a, bi = b, ci = cc, di = d)
  or_hat     <- exp(as.numeric(mf_or$yi))
  logor_se   <- sqrt(as.numeric(mf_or$vi))
  target_se  <- sqrt(as.numeric(mf_rr$vi))

  for (method in c("metaumbrella_cases", "metaumbrella_exp")) {
    got <- suppressMessages(suppressWarnings(
      es_from_or_se(or = or_hat, logor_se = logor_se,
                    n_exp = n_exp, n_nexp = n_nexp,
                    n_cases = a + cc, n_controls = b + d, or_to_rr = method)
    ))
    # EXTERNAL pin: the reported SE is exactly the table's own, so the propagation
    # must reproduce metafor's RR standard error, not merely something proportional.
    expect_equal(got$logrr_se, target_se, tolerance = 1e-8,
                 info = paste("or_to_rr =", method))
    # ... and the same study entered as a 2x2 table must not disagree with itself.
    tab <- es_from_2x2(n_cases_exp = a, n_controls_exp = b,
                       n_cases_nexp = cc, n_controls_nexp = d)
    expect_equal(got$logrr_se, tab$logrr_se, tolerance = 1e-8,
                 info = paste("or_to_rr =", method))
  }

  # Same statement as a scale invariant, independent of any numeric constant:
  # logrr_se must be linear in logor_se with a zero intercept. This is the property
  # the propagation exists to provide -- before it, a 30x range in the reported
  # logor_se returned a byte-identical logrr_se.
  r1 <- suppressMessages(suppressWarnings(
    es_from_or_se(or = 2, logor_se = 0.05, n_exp = 300, n_nexp = 700,
                  n_cases = 200, n_controls = 800)))
  r2 <- suppressMessages(suppressWarnings(
    es_from_or_se(or = 2, logor_se = 1.5, n_exp = 300, n_nexp = 700,
                  n_cases = 200, n_controls = 800)))
  expect_equal(r2$logrr_se / r1$logrr_se, 1.5 / 0.05, tolerance = 1e-6)

  # The CI must be rebuilt from the propagated SE, not left at the crude-table one,
  # so a null OR cannot come back as a significant RR.
  expect_equal(r2$logrr_ci_up - r2$logrr_ci_lo, 2 * qnorm(.975) * r2$logrr_se,
               tolerance = 1e-8)
  expect_gt(r2$logrr_ci_up, 0)
  expect_lt(r2$logrr_ci_lo, 0)
})

# --------------------------------------------------------------------------
# AUDIT line 423 --- or_list_L2 (R/main_convert_df.R:1446) ranks es_odds_ratio
# (SE IMPUTED from the margins by .se_from_or) ABOVE es_odds_ratio_pval (SE
# recovered exactly as |logOR / qnorm(p/2)|), contradicting the @note added to
# es_from_or() in this release. On a row carrying or + or_pval + margins the
# pipeline selects the imputed SE, ~1.49x too wide, i.e. ~2.2x too little weight.
#
# Expected SE derived independently from the generating 2x2 (40,160,20,180):
#   SE(logOR) = sqrt(1/40 + 1/160 + 1/20 + 1/180) = 0.2946278255
# and the p-value was constructed from exactly that SE, so the p-value route
# recovers it to machine precision.
test_that("AUDIT-major: convert_df prefers the exact p-value OR SE over the imputed one", {
  a <- 40; b <- 160; cc <- 20; d <- 180
  or <- (a * d) / (b * cc)                       # 2.25
  se_true <- sqrt(1 / a + 1 / b + 1 / cc + 1 / d) # 0.2946278255
  pv <- 2 * pnorm(abs(log(or) / se_true), lower.tail = FALSE)

  x <- data.frame(study_id = "s1", or = or, or_pval = pv,
                  n_cases = 60, n_controls = 340)
  s <- suppressMessages(as.data.frame(
    summary(convert_df(x, measure = "logor", verbose = FALSE))))

  # Today: info_used_crude == "or", se_crude == 0.4390289 (ratio 1.4901 to truth).
  expect_identical(as.character(s$info_used_crude[1]), "or_pval")
  expect_equal(s$se_crude[1], se_true, tolerance = 1e-6)
  expect_equal(s$es_crude[1], log(or), tolerance = 1e-9)
})

# --------------------------------------------------------------------------
# AUDIT line 749 (+ Verifier correction at line 760) --- `or` receives no
# positivity guard on entry to the standalone ratio routes. A log OR typed into
# the adjacent `or` column (negative for any protective effect) self-blanks the
# log-scale outputs but leaves the Grant risk-difference block finite: it emits
# a fabricated rd, a NEGATIVE rd_se/nnt_se (the delta line omits abs() on the
# `or` factor) and a transposed CI. es_from_rr_se() has the identical defect.
# The verifier adds that or = 0 is the worse case (rd_se == 0, i.e. an infinite
# inverse-variance weight), which is why the guard must test `or <= 0` rather
# than the `< 0` used by .positive_columns().
test_that("AUDIT-minor: or <= 0 and rr <= 0 are guarded to NA before the Grant RD block", {
  neg <- suppressMessages(suppressWarnings(
    es_from_or_se(or = -0.7, logor_se = 0.2, baseline_risk = 0.3,
                  n_cases = 50, n_controls = 150)))
  # Today: rd = 0.7285714, rd_se = -0.122449, nnt = 1.372549, nnt_se = -0.2306805,
  #        rd_ci_lo = 0.968567 > rd_ci_up = 0.4885758 (transposed).
  expect_true(is.na(neg$rd))
  expect_true(is.na(neg$rd_se))
  expect_true(is.na(neg$nnt))
  expect_true(is.na(neg$nnt_se))

  zero <- suppressMessages(suppressWarnings(
    es_from_or_se(or = 0, logor_se = 0.2, baseline_risk = 0.3,
                  n_cases = 50, n_controls = 150)))
  # Today: rd = 0.3 with rd_se = 0 exactly (infinite meta-analytic weight).
  expect_true(is.na(zero$rd))
  expect_true(is.na(zero$rd_se))

  negrr <- suppressMessages(suppressWarnings(
    es_from_rr_se(rr = -0.5, logrr_se = 0.2, baseline_risk = 0.3,
                  n_exp = 50, n_nexp = 50)))
  # Today: rd = 0.45, rd_se = -0.03, nnt = 2.222222, nnt_se = -0.1481481.
  expect_true(is.na(negrr$rd))
  expect_true(is.na(negrr$rd_se))
  expect_true(is.na(negrr$nnt))
  expect_true(is.na(negrr$nnt_se))
})

# --------------------------------------------------------------------------
# AUDIT line 762 (+ Verifier correction at line 773) --- .se_from_or()
# (R/internal_multiple_formulas.R:696) returns NaN, not NA_real_, whenever its
# enumeration of candidate 2x2 tables comes back empty: the NaN is
# mean(numeric(0)) after the feasibility filter NA's every candidate. It
# propagates unguarded to logor_se / logor_ci_lo / logor_ci_up beside a finite
# logor, and .positive_or_na() cannot catch it (is.na(NaN) is TRUE, so the
# x <= 0 branch is skipped). The verifier widens the trigger beyond n_cases <= 1
# to any empty enumeration (n_controls <= 1; small margins with an extreme OR)
# and endorses `if (!is.finite(v_or_mean)) return(res)` as sufficient. The OR
# itself is the user's own input and must be kept (documented convention).
test_that("AUDIT-minor: .se_from_or returns NA rather than NaN on an empty enumeration", {
  # Internal helper: takes a length-3 vector c(or, n_cases, n_controls).
  for (v in list(c(2, 0, 50), c(2, 1, 50), c(2, 50, 1), c(100, 5, 5))) {
    got <- metaConvert:::.se_from_or(v)
    expect_false(is.nan(got$se),
                 info = paste("se_from_or(", paste(v, collapse = ", "), ")"))
    expect_false(is.nan(got$var),
                 info = paste("se_from_or(", paste(v, collapse = ", "), ")"))
    expect_true(is.na(got$se))
  }
  # A legitimate enumeration is untouched.
  ok <- metaConvert:::.se_from_or(c(2, 10, 50))
  expect_equal(ok$se, 0.8770824383, tolerance = 1e-8)

  # Exported route: NaN must not reach the log-OR columns, while the reported
  # OR is kept (es visible, se = NA -- the bare-omega convention).
  r <- suppressMessages(suppressWarnings(
    es_from_or(or = 2, n_cases = 1, n_controls = 50)))
  expect_equal(r$logor, log(2), tolerance = 1e-12)
  expect_false(is.nan(r$logor_se))
  expect_false(is.nan(r$logor_ci_lo))
  expect_false(is.nan(r$logor_ci_up))
  expect_true(is.na(r$logor_se))
})

# --------------------------------------------------------------------------
# AUDIT line 801 (+ Verifier correction at line 812) --- a derived rd_se of
# exactly 0 (a zero sampling variance, i.e. an infinite inverse-variance weight)
# reaches the output and is selected by summary(). At baseline_risk = 0 the
# Grant chain gives drd_dor = BR(1-BR)/denom^2 = 0 (R/es_from_stand_OR.R:559-561)
# and rd_se = BR*rr*logrr_se = 0 (R/es_from_stand_RR.R:206); es_from_2x2() has
# the same shape on a double-zero table, RD_SE = sqrt(pt(1-pt)/n1 + pc(1-pc)/n2).
#
# The verifier is explicit that ONLY the second fix option is admissible: pass
# the derived rd_se through .positive_or_na() at the point of emission. The
# first option (narrowing .baseline_risk_or_na() to x <= 0) must NOT be applied
# -- baseline_risk = 0 is the deliberately kept rare-disease limit, pinned by
# tests/testthat/test-baseline-risk-guard.R:90 -- so the log-scale outputs must
# survive untouched here.
test_that("AUDIT-minor: a derived rd_se of exactly zero is declined, not emitted", {
  a <- suppressMessages(suppressWarnings(
    es_from_or_se(or = 0.5, logor_se = 0.2, baseline_risk = 0)))
  expect_false(isTRUE(a$rd_se == 0))   # today: rd = 0, rd_se = 0
  expect_true(is.na(a$rd_se))
  # The rare-disease limit on the log scale must be preserved (BR = 0 stays legal).
  expect_equal(a$logor, log(0.5), tolerance = 1e-12)
  expect_equal(a$logor_se, 0.2, tolerance = 1e-12)

  b <- suppressMessages(suppressWarnings(
    es_from_rr_se(rr = 0.5, logrr_se = 0.2, baseline_risk = 0)))
  expect_false(isTRUE(b$rd_se == 0))
  expect_true(is.na(b$rd_se))
  expect_equal(b$logrr, log(0.5), tolerance = 1e-12)
  expect_equal(b$logrr_se, 0.2, tolerance = 1e-12)

  # Verifier-added scope: the 2x2 double-zero table has the identical failure.
  d <- suppressMessages(suppressWarnings(
    es_from_2x2(n_cases_exp = 0, n_controls_exp = 50,
                n_cases_nexp = 0, n_controls_nexp = 50)))
  expect_false(isTRUE(d$rd_se == 0))
  expect_true(is.na(d$rd_se))
})

# --------------------------------------------------------------------------
# AUDIT line 814 (+ Verifier correction at line 825) --- the Altman (1998)
# discontinuity test is written with strict inequalities
# (rd_ci_lo_raw < 0 & rd_ci_up_raw > 0), so a risk-difference bound landing
# EXACTLY on 0 skips the guard and the reciprocal yields a literal Inf, where
# the neighbouring input (bound epsilon below 0) correctly returns NA.
#
# The verifier records two corrections that this test encodes: (1) the boundary
# is the GENERIC outcome of the common reporting pattern "RD 0.20, 95% CI 0.00
# to 0.40" for es_from_rd_ci()/es_from_rd_pval(), not a measure-zero curiosity;
# (2) for a directly entered NEGATIVE rd the emitted bound is a wrong NUMBER,
# not just a wrong type -- nnt_ci_lo = 1/(+0) = +Inf where the confidence set is
# (-Inf, -2.5], giving an inverted interval (lo > up) with a sign-wrong bound.
test_that("AUDIT-minor: NNT CI bounds are NA, not Inf, when an RD CI limit is exactly zero", {
  se0 <- 0.2 / qnorm(.975)   # makes rd_ci_lo land exactly on 0 for rd = 0.2

  pos <- suppressMessages(suppressWarnings(es_from_rd_se(rd = 0.2, rd_se = se0)))
  expect_equal(pos$rd_ci_lo, 0, tolerance = 1e-12)  # confirms we are on the boundary
  expect_false(is.infinite(pos$nnt_ci_up))          # today: Inf
  expect_true(is.na(pos$nnt_ci_lo))
  expect_true(is.na(pos$nnt_ci_up))
  # Point estimate and SE are unaffected by the fix.
  expect_equal(pos$nnt, 5, tolerance = 1e-9)

  # Negative-RD sub-case: today nnt_ci_lo = Inf and nnt_ci_up = -2.5 (inverted,
  # sign-wrong). Both bounds must be NA.
  neg <- suppressMessages(suppressWarnings(es_from_rd_se(rd = -0.2, rd_se = se0)))
  expect_false(is.infinite(neg$nnt_ci_lo))
  expect_true(is.na(neg$nnt_ci_lo))
  expect_true(is.na(neg$nnt_ci_up))
  expect_equal(neg$nnt, -5, tolerance = 1e-9)

  # The reporting pattern that makes this generic rather than incidental.
  ci <- suppressMessages(suppressWarnings(
    es_from_rd_ci(rd = 0.2, rd_ci_lo = 0, rd_ci_up = 0.4)))
  expect_false(is.infinite(ci$nnt_ci_up))
  expect_true(is.na(ci$nnt_ci_lo))
  expect_true(is.na(ci$nnt_ci_up))

  # The neighbouring input (bound strictly below 0) already behaves correctly
  # and must keep doing so -- this half is a guard against over-correcting.
  nb <- suppressMessages(suppressWarnings(
    es_from_rd_se(rd = 0.2, rd_se = se0 * 1.0001)))
  expect_lt(nb$rd_ci_lo, 0)
  expect_true(is.na(nb$nnt_ci_lo))
  expect_true(is.na(nb$nnt_ci_up))
})
