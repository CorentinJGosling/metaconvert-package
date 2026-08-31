# Red-half regression tests for four defects of the user-input route recorded in
# AUDIT-2026-08-28-findings.md. Every assertion below states the CORRECT behaviour,
# so each test_that() fails on the code as it stands and passes once the defect is
# repaired. Nothing here is a snapshot of current output.

# Collect every warning an expression emits without letting it abort the call, so a
# test can assert on the ABSENCE of one particular warning rather than on silence.
.collect_warnings <- function(expr) {
  w <- character(0)
  suppressMessages(withCallingHandlers(
    force(expr),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  ))
  w
}

.quiet <- function(expr) suppressWarnings(suppressMessages(force(expr)))


# AUDIT-2026-08-28-findings.md:610-619
# Defect: for or/rr/irr/hr entered as a natural-scale CI *together with* a log-scale SE
# and no point estimate, the `ratio_from_ci` block is skipped (it is gated on
# is.na(user_se_crude)), so the ES is the ARITHMETIC midpoint of the natural CI
# ((1.5+4)/2 = 2.75) instead of the geometric one (sqrt(1.5*4) = 2.4494897).
test_that("AUDIT-ratio-ci-midpoint: a ratio ES recovered from a natural-scale CI is the geometric midpoint even when a log-scale SE is supplied", {
  ci_lo <- 1.5
  ci_up <- 4

  # Independently derived: on the log scale the estimate is the midpoint of the
  # reported bounds, i.e. the GEOMETRIC mean of the natural bounds.
  expected_log_es <- (log(ci_lo) + log(ci_up)) / 2      # 0.8958797346
  expected_se     <- (log(ci_up) - log(ci_lo)) / (2 * qnorm(.975))  # 0.2502161419

  expect_equal(exp(expected_log_es), sqrt(ci_lo * ci_up), tolerance = 1e-12)

  targets <- c(or = "logor", rr = "logrr", irr = "logirr", hr = "loghr")

  for (m in names(targets)) {
    col <- targets[[m]]

    res <- .quiet(es_from_user_crude(
      user_es_original_measure_crude = m,
      user_es_crude   = NA,
      user_se_crude   = expected_se,   # correct, on the log scale
      user_ci_lo_crude = ci_lo,
      user_ci_up_crude = ci_up,
      user_es_target_measure_crude = m,
      n_exp = 50, n_nexp = 50
    ))

    expect_equal(res[[col]], expected_log_es, tolerance = 1e-8,
                 info = paste0("crude route, measure = ", m))
    expect_equal(exp(res[[col]]), sqrt(ci_lo * ci_up), tolerance = 1e-8,
                 info = paste0("crude route, measure = ", m))
    # The reported CI must be preserved verbatim on the log scale, and the SE with it.
    expect_equal(res[[paste0(col, "_se")]], expected_se, tolerance = 1e-8)
    expect_equal(res[[paste0(col, "_ci_lo")]], log(ci_lo), tolerance = 1e-8)
    expect_equal(res[[paste0(col, "_ci_up")]], log(ci_up), tolerance = 1e-8)
  }

  # The adjusted route carries the identical block and the identical defect.
  res_adj <- .quiet(es_from_user_adj(
    user_es_original_measure_adj = "or",
    user_es_adj    = NA,
    user_se_adj    = expected_se,
    user_ci_lo_adj = ci_lo,
    user_ci_up_adj = ci_up,
    user_es_target_measure_adj = "or",
    n_exp = 50, n_nexp = 50
  ))
  expect_equal(res_adj$logor, expected_log_es, tolerance = 1e-8)

  # And the same row must not change value merely because the SE cell was filled in:
  # the CI-only entry (which already takes the geometric route) is the reference.
  res_ci_only <- .quiet(es_from_user_crude(
    user_es_original_measure_crude = "or",
    user_es_crude = NA, user_se_crude = NA,
    user_ci_lo_crude = ci_lo, user_ci_up_crude = ci_up,
    user_es_target_measure_crude = "or",
    n_exp = 50, n_nexp = 50
  ))
  res_with_se <- .quiet(es_from_user_crude(
    user_es_original_measure_crude = "or",
    user_es_crude = NA, user_se_crude = expected_se,
    user_ci_lo_crude = ci_lo, user_ci_up_crude = ci_up,
    user_es_target_measure_crude = "or",
    n_exp = 50, n_nexp = 50
  ))
  expect_equal(res_with_se$logor, res_ci_only$logor, tolerance = 1e-10)
})


# AUDIT-2026-08-28-findings.md:1126-1137 (see the Verifier correction, which OVERRIDES
# the filed Fix: the d<->g step for the within-group entered types must use the
# within-subject correction J(n - 1), NOT the pooled J(n_exp + n_nexp - 2) that
# routing through .es_from_d() would apply, and the r/z/OR family must stay refused).
# Defect: the "dw"/"gw" branches of .dispatch_user_conversion() never perform the
# Hedges step, so an entered dw yields gw = NA and an entered gw yields dw = NA -- the
# row is silently dropped from the pool, while the numerically identical d/g entry converts.
test_that("AUDIT-within-subject-user-es: an entered dw/gw converts to its g/d counterpart with the within-subject J(n-1) instead of being dropped", {
  n <- 50
  # Hedges' small-sample correction, closed form; within-subject df = n - 1.
  J <- function(df) exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
  j49 <- J(n - 1)                                  # 0.984602177236
  expect_equal(j49, 0.984602177236, tolerance = 1e-10)

  # -- entered dw, g-family target -------------------------------------------------
  res_dw <- .quiet(es_from_user_crude(
    user_es_original_measure_crude = "dw",
    user_es_crude = 0.5, user_se_crude = 0.2,
    user_es_target_measure_crude = "gw",
    n_exp = n, n_nexp = NA
  ))
  expect_false(is.na(res_dw$gw))
  expect_equal(res_dw$gw,    0.5 * j49, tolerance = 1e-8)   # 0.492301088618
  expect_equal(res_dw$gw_se, 0.2 * j49, tolerance = 1e-8)   # 0.196920435447
  # the entered value itself is untouched
  expect_equal(res_dw$dw, 0.5, tolerance = 1e-12)

  # -- entered gw, d-family target -------------------------------------------------
  res_gw <- .quiet(es_from_user_crude(
    user_es_original_measure_crude = "gw",
    user_es_crude = 0.5, user_se_crude = 0.2,
    user_es_target_measure_crude = "dw",
    n_exp = n, n_nexp = NA
  ))
  expect_false(is.na(res_gw$dw))
  expect_equal(res_gw$dw,    0.5 / j49, tolerance = 1e-8)   # 0.507819311759
  expect_equal(res_gw$dw_se, 0.2 / j49, tolerance = 1e-8)   # 0.203127724703
  expect_equal(res_gw$gw, 0.5, tolerance = 1e-12)

  # The pooled-df correction is the WRONG one here and must not be what appears.
  j98 <- J(2 * n - 2)                              # 0.992324087421
  expect_false(isTRUE(all.equal(res_dw$gw, 0.5 * j98, tolerance = 1e-8)))

  # Verifier correction: the r/z/OR family stays refused for a paired SMD (the
  # point-biserial and Cox maps assume two independent groups). This must remain NA.
  res_r <- .quiet(es_from_user_crude(
    user_es_original_measure_crude = "dw",
    user_es_crude = 0.5, user_se_crude = 0.2,
    user_es_target_measure_crude = "r",
    n_exp = n, n_nexp = NA
  ))
  expect_true(is.na(res_r$r))

  # End to end: neither within-subject row may drop out of its pool.
  dat <- data.frame(
    study_id = c("w_dw", "w_gw"),
    user_es_original_measure_crude = c("dw", "gw"),
    user_es_crude = 0.5, user_se_crude = 0.2,
    n_exp = n, n_nexp = NA
  )
  s_gw <- .quiet(summary(convert_df(dat, measure = "gw", verbose = FALSE)))
  s_gw <- s_gw[match(c("w_dw", "w_gw"), s_gw$study_id), ]
  expect_equal(s_gw$es_crude, c(0.5 * j49, 0.5), tolerance = 1e-8)

  s_dw <- .quiet(summary(convert_df(dat, measure = "dw", verbose = FALSE)))
  s_dw <- s_dw[match(c("w_dw", "w_gw"), s_dw$study_id), ]
  expect_equal(s_dw$es_crude, c(0.5, 0.5 / j49), tolerance = 1e-8)
})


# AUDIT-2026-08-28-findings.md:1100-1111 (Verifier correction CORRECTION 2 overrides the
# filed Fix: subsetting the merged column by the row indices the branch keeps -- not
# rep(..., length.out = length(nn)) -- is what is required; the filed recipe misaligns
# every row after a dropped one, which is strictly worse than the defect).
# Defect: convert_df() passes its SCALAR method arguments to the two user-input calls
# instead of the per-row columns it merged into x, so a user_input row silently gets the
# dataset default rather than the conversion its own cell requested.
test_that("AUDIT-per-row-methods-user-route: convert_df honours per-row smd_to_cor / or_to_cor / cor_to_smd columns on the user-input route", {
  # (a) per-row smd_to_cor, entered d -> target r ------------------------------------
  # Ground truth from the native route, which does honour the argument, plus the
  # closed-form point-biserial for the lipsey_cooper arm.
  r_viech  <- es_from_cohen_d(cohen_d = 0.5, n_exp = 50, n_nexp = 50,
                              smd_to_cor = "viechtbauer")$r          # 0.3068752868
  r_lipsey <- es_from_cohen_d(cohen_d = 0.5, n_exp = 50, n_nexp = 50,
                              smd_to_cor = "lipsey_cooper")$r        # 0.2425356250
  expect_equal(r_lipsey, 0.5 / sqrt(0.5^2 + 4), tolerance = 1e-10)   # d / sqrt(d^2 + 4)

  dat_smd <- data.frame(
    study_id = c("a", "b"),
    user_es_original_measure_crude = "d",
    user_es_crude = 0.5, user_se_crude = 0.2,
    n_exp = 50, n_nexp = 50,
    smd_to_cor = c("viechtbauer", "lipsey_cooper"),
    stringsAsFactors = FALSE
  )
  s_smd <- .quiet(summary(convert_df(dat_smd, measure = "r", verbose = FALSE)))
  s_smd <- s_smd[match(c("a", "b"), s_smd$study_id), ]
  expect_equal(s_smd$es_crude, c(r_viech, r_lipsey), tolerance = 1e-7)

  # (b) per-row or_to_cor, entered or -> target r ------------------------------------
  r_pearson <- .quiet(es_from_or_se(or = 2.5, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                                    or_to_cor = "pearson")$r)        # 0.3463354903
  r_digby   <- .quiet(es_from_or_se(or = 2.5, logor_se = 0.2, n_exp = 50, n_nexp = 50,
                                    or_to_cor = "digby")$r)          # 0.3306955648
  expect_false(isTRUE(all.equal(r_pearson, r_digby, tolerance = 1e-6)))

  dat_or <- data.frame(
    study_id = c("a", "b"),
    user_es_original_measure_crude = "or",
    user_es_crude = 2.5, user_se_crude = 0.2,
    n_exp = 50, n_nexp = 50,
    or_to_cor = c("pearson", "digby"),
    stringsAsFactors = FALSE
  )
  s_or <- .quiet(summary(convert_df(dat_or, measure = "r", verbose = FALSE)))
  s_or <- s_or[match(c("a", "b"), s_or$study_id), ]
  expect_equal(s_or$es_crude, c(r_pearson, r_digby), tolerance = 1e-7)

  # (c) per-row cor_to_smd, entered r -> target d, with a row the branch DROPS --------
  # Row "a" carries no SE and no CI, so the "r" branch of .dispatch_user_conversion()
  # excludes it. Rows "b" and "c" must still receive THEIR OWN cor_to_smd cells; a fix
  # that recycles the column onto the kept rows positionally would hand row "b" the
  # cell of row "a" and row "c" the cell of row "b".
  d_viech  <- es_from_pearson_r(pearson_r = 0.5, n_sample = 100, n_exp = 50, n_nexp = 50,
                                cor_to_smd = "viechtbauer")$d        # 0.8768565296
  d_cooper <- es_from_pearson_r(pearson_r = 0.5, n_sample = 100, n_exp = 50, n_nexp = 50,
                                cor_to_smd = "cooper")$d             # 1.1547005384
  expect_equal(d_cooper, 2 * 0.5 / sqrt(1 - 0.5^2), tolerance = 1e-10)  # 2r / sqrt(1 - r^2)

  dat_cor <- data.frame(
    study_id = c("a", "b", "c"),
    user_es_original_measure_crude = "r",
    user_es_crude = 0.5,
    user_se_crude = c(NA, 0.08, 0.08),
    n_exp = 50, n_nexp = 50, n_sample = 100,
    cor_to_smd = c("cooper", "viechtbauer", "cooper"),
    stringsAsFactors = FALSE
  )
  s_cor <- .quiet(summary(convert_df(dat_cor, measure = "d", verbose = FALSE)))
  s_cor <- s_cor[match(c("a", "b", "c"), s_cor$study_id), ]
  expect_true(is.na(s_cor$es_crude[1]))
  expect_equal(s_cor$es_crude[2:3], c(d_viech, d_cooper), tolerance = 1e-7)
})


# AUDIT-2026-08-28-findings.md:1372-1383 (Verifier correction: the fix is to add
# `user_es_crude > 1` to the predicate -- 0.00% false positives above 1, and the
# statistic is anti-informative below it).
# Defect: the scale-sanity check compares a LOG-scale SE against a NATURAL-scale ratio
# (`user_se_crude > user_es_crude`), so a perfectly correct protective OR/RR entry
# (OR = 0.20 with SE(logOR) = 0.30) is warned as probably mis-scaled.
test_that("AUDIT-protective-or-se-warning: the 'SE may be on the natural scale' warning does not fire for a valid protective OR/RR", {
  pattern <- "may be on the natural scale"

  w_or <- .collect_warnings(es_from_user_crude(
    user_es_original_measure_crude = "or",
    user_es_crude = 0.20, user_se_crude = 0.30,
    user_es_target_measure_crude = "or",
    n_exp = 50, n_nexp = 50
  ))
  expect_length(grep(pattern, w_or), 0)

  w_rr <- .collect_warnings(es_from_user_crude(
    user_es_original_measure_crude = "rr",
    user_es_crude = 0.20, user_se_crude = 0.30,
    user_es_target_measure_crude = "rr",
    n_exp = 50, n_nexp = 50
  ))
  expect_length(grep(pattern, w_rr), 0)

  # The adjusted route carries an identical copy of the predicate.
  w_adj <- .collect_warnings(es_from_user_adj(
    user_es_original_measure_adj = "or",
    user_es_adj = 0.20, user_se_adj = 0.30,
    user_es_target_measure_adj = "or",
    n_exp = 50, n_nexp = 50
  ))
  expect_length(grep(pattern, w_adj), 0)

  # The entry silenced above is correct, and must stay correct.
  res <- .quiet(es_from_user_crude(
    user_es_original_measure_crude = "or",
    user_es_crude = 0.20, user_se_crude = 0.30,
    user_es_target_measure_crude = "or",
    n_exp = 50, n_nexp = 50
  ))
  expect_equal(res$logor, log(0.20), tolerance = 1e-10)
  expect_equal(res$logor_se, 0.30, tolerance = 1e-10)

  # Anti-over-fix pin: above 1 the statistic does discriminate, so the warning must
  # survive there. (This half already holds today.)
  w_harmful <- .collect_warnings(es_from_user_crude(
    user_es_original_measure_crude = "or",
    user_es_crude = 1.5, user_se_crude = 2.0,
    user_es_target_measure_crude = "or",
    n_exp = 50, n_nexp = 50
  ))
  expect_gte(length(grep(pattern, w_harmful)), 1)
})
