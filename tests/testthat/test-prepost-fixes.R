# =============================================================================
# Coverage for the pre/post & change-score fixes and for the surface the audit
# found untested: pooled-vs-legacy standardizers, the corrected pooled morris_dz
# variance, degenerate inputs (SD = 0, |r| >= 1, tiny m), the per-row method
# vector, the convert_df pool_sd threading and coercion/r-imputation messages,
# and the E6/E7 estimand-mixing flags.
#
# Design rule for this file: assertions are anchored to an INDEPENDENT reference
# (metafor::escalc, or an algebraic identity), never to a re-typed copy of the
# package's own source expression.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Pooled morris_dz reproduces metafor::escalc(measure = "SMD") on change scores
# -----------------------------------------------------------------------------
test_that("pooled morris_dz matches escalc SMD on change scores (point + SE)", {
  skip_if_not_installed("metafor")
  set.seed(1)
  # deliberately UNEQUAL arm change-SDs: the regime where a per-arm standardizer
  # and a pooled standardizer diverge
  cases <- list(
    list(n1 = 29, n2 = 34, m1 = -12.4, s1 = 6.52, m2 = -3.5, s2 = 10.67),
    list(n1 = 39, n2 = 38, m1 = -14.9, s1 = 6.70, m2 = -7.4, s2 = 9.70),
    list(n1 = 21, n2 = 18, m1 = -5.8, s1 = 7.67, m2 = 1.9, s2 = 7.22),
    list(n1 = 12, n2 = 60, m1 = 3.0, s1 = 2.10, m2 = 0.4, s2 = 9.90)
  )
  for (cs in cases) {
    res <- es_from_mean_change_sd(
      n_exp = cs$n1, n_nexp = cs$n2,
      mean_change_exp = cs$m1, mean_change_sd_exp = cs$s1,
      mean_change_nexp = cs$m2, mean_change_sd_nexp = cs$s2,
      r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6,
      pre_post_to_smd = "morris_dz", pool_sd = TRUE
    )
    ref <- metafor::escalc(
      measure = "SMD",
      m1i = cs$m1, sd1i = cs$s1, n1i = cs$n1,
      m2i = cs$m2, sd2i = cs$s2, n2i = cs$n2
    )
    # point estimate: exact (same estimator)
    expect_equal(res$g, as.numeric(ref$yi), tolerance = 1e-6)
    # SE: escalc uses vi = 1/n1 + 1/n2 + g^2/(2N); the package uses the
    # J^2 (N/(n1 n2) + g^2/(2m)) convention. Equivalent to first order.
    # Under the OLD (defective) code this ratio was ~0.89 at r = 0.6, so 3%
    # is tight enough to catch a reintroduced 2(1-r) factor.
    expect_equal(res$g_se, sqrt(as.numeric(ref$vi)), tolerance = 3e-2)
  }
})

test_that("pooled morris_dz SE is INVARIANT to r_pre_post (no 2(1-r) factor)", {
  # The change-SD-standardized estimator conditions on the change SD, which is
  # supplied directly here; r plays no part in either the point estimate or the
  # variance. The pre-fix code made the SE swing with r -- this is the sharpest
  # regression guard against reintroducing that factor.
  mk <- function(r) es_from_mean_change_sd(
    n_exp = 29, n_nexp = 34,
    mean_change_exp = -12.4, mean_change_sd_exp = 6.52,
    mean_change_nexp = -3.5, mean_change_sd_nexp = 10.67,
    r_pre_post_exp = r, r_pre_post_nexp = r,
    pre_post_to_smd = "morris_dz", pool_sd = TRUE
  )
  a <- mk(0.2); b <- mk(0.6); cc <- mk(0.9)
  expect_equal(a$g, b$g, tolerance = 1e-12)
  expect_equal(a$g_se, b$g_se, tolerance = 1e-12)
  expect_equal(b$g_se, cc$g_se, tolerance = 1e-12)
})

test_that("pooled morris_drm = pooled morris_dz * sqrt(2(1-r)) in BOTH g and SE", {
  # Morris & DeShon (2002) metric transformation. Holds for the point estimate
  # by construction and for the SE because r is treated as a known constant.
  r <- 0.6
  args <- list(
    n_exp = 29, n_nexp = 34,
    mean_change_exp = -12.4, mean_change_sd_exp = 6.52,
    mean_change_nexp = -3.5, mean_change_sd_nexp = 10.67,
    r_pre_post_exp = r, r_pre_post_nexp = r, pool_sd = TRUE
  )
  dz <- do.call(es_from_mean_change_sd, c(args, pre_post_to_smd = "morris_dz"))
  drm <- do.call(es_from_mean_change_sd, c(args, pre_post_to_smd = "morris_drm"))
  k <- sqrt(2 * (1 - r))
  expect_equal(drm$g, dz$g * k, tolerance = 1e-10)
  expect_equal(drm$g_se, dz$g_se * k, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# 2. pool_sd default, and the legacy per-arm path
# -----------------------------------------------------------------------------
test_that("pool_sd defaults to TRUE in the two-group wrappers", {
  args <- list(
    n_exp = 29, n_nexp = 34,
    mean_change_exp = -12.4, mean_change_sd_exp = 6.52,
    mean_change_nexp = -3.5, mean_change_sd_nexp = 10.67,
    r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6
  )
  default <- do.call(es_from_mean_change_sd, args)
  pooled <- do.call(es_from_mean_change_sd, c(args, pool_sd = TRUE))
  legacy <- do.call(es_from_mean_change_sd, c(args, pool_sd = FALSE))

  expect_equal(default$g, pooled$g, tolerance = 1e-12)
  expect_equal(default$g_se, pooled$g_se, tolerance = 1e-12)
  # and the legacy path really is different on unequal arm SDs
  expect_false(isTRUE(all.equal(legacy$g, pooled$g, tolerance = 1e-3)))
})

test_that("legacy pool_sd = FALSE reproduces the per-arm difference of within-group g's", {
  # Regression test OF THE LEGACY PATH (retained for backward compatibility):
  # g = J(n1-1) * d_rm_exp - J(n2-1) * d_rm_nexp, each arm on its OWN change SD.
  n1 <- 29; n2 <- 34; r <- 0.6
  res <- es_from_mean_change_sd(
    n_exp = n1, n_nexp = n2,
    mean_change_exp = -12.4, mean_change_sd_exp = 6.52,
    mean_change_nexp = -3.5, mean_change_sd_nexp = 10.67,
    r_pre_post_exp = r, r_pre_post_nexp = r,
    pre_post_to_smd = "morris_drm", pool_sd = FALSE
  )
  k <- sqrt(2 * (1 - r))
  d_exp <- (-12.4 / 6.52) * k
  d_nexp <- (-3.5 / 10.67) * k
  Jf <- function(df) 1 - 3 / (4 * df - 1)
  g_legacy <- Jf(n1 - 1) * d_exp - Jf(n2 - 1) * d_nexp
  expect_equal(res$g, g_legacy, tolerance = 1e-3)
})

test_that("pooled and legacy coincide when the two arms' SDs and n are equal", {
  # The one regime in which the per-arm difference IS a valid between-group SMD.
  args <- list(
    n_exp = 40, n_nexp = 40,
    mean_change_exp = 6, mean_change_sd_exp = 8,
    mean_change_nexp = 2, mean_change_sd_nexp = 8,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "morris_drm"
  )
  pooled <- do.call(es_from_mean_change_sd, c(args, pool_sd = TRUE))
  legacy <- do.call(es_from_mean_change_sd, c(args, pool_sd = FALSE))
  # d is identical; g differs only through the J degrees of freedom (N-2 vs n-1)
  expect_equal(pooled$d, legacy$d, tolerance = 1e-10)
  expect_equal(pooled$g, legacy$g, tolerance = 1e-2)
})

# -----------------------------------------------------------------------------
# 3. Degenerate inputs return NA rather than Inf/NaN
# -----------------------------------------------------------------------------
test_that("a zero standardizing SD yields NA, not Inf", {
  z <- es_from_mean_change_sd(
    n_exp = 30, n_nexp = 30,
    mean_change_exp = 5, mean_change_sd_exp = 0,
    mean_change_nexp = 2, mean_change_sd_nexp = 4,
    r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6, pool_sd = FALSE
  )
  expect_true(is.na(z$d))
  expect_false(isTRUE(is.infinite(z$d)))

  zp <- es_from_mean_change_sd(
    n_exp = 30, n_nexp = 30,
    mean_change_exp = 5, mean_change_sd_exp = 0,
    mean_change_nexp = 2, mean_change_sd_nexp = 0,
    r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6, pool_sd = TRUE
  )
  expect_true(is.na(zp$d))
})

test_that("|r_pre_post| >= 1 yields NA, while r = 0.99 still computes", {
  base <- list(
    n_exp = 30, n_nexp = 30,
    mean_pre_exp = 10, mean_exp = 12, mean_pre_sd_exp = 3, mean_sd_exp = 3,
    mean_pre_nexp = 10, mean_nexp = 11, mean_pre_sd_nexp = 3, mean_sd_nexp = 3,
    pre_post_to_smd = "morris_dz"
  )
  # r = 1 with equal pre/post SDs collapses the change SD to 0
  bad <- do.call(es_from_means_sd_pre_post,
                 c(base, r_pre_post_exp = 1, r_pre_post_nexp = 1))
  expect_true(is.na(bad$d))
  expect_false(isTRUE(is.nan(bad$d)))

  ok <- do.call(es_from_means_sd_pre_post,
                c(base, r_pre_post_exp = 0.99, r_pre_post_nexp = 0.99))
  expect_true(is.finite(ok$d))
})

# -----------------------------------------------------------------------------
# 4. Per-row pre_post_to_smd vector (previously fell wholesale into drm)
# -----------------------------------------------------------------------------
test_that("a mixed pre_post_to_smd vector is applied per row, not wholesale", {
  mixed <- es_from_paired_t(
    paired_t_exp = c(3, 3), paired_t_nexp = c(1, 1),
    n_exp = c(50, 50), n_nexp = c(50, 50),
    r_pre_post_exp = c(0.6, 0.6), r_pre_post_nexp = c(0.6, 0.6),
    pre_post_to_smd = c("morris_dz", "morris_drm")
  )
  all_dz <- es_from_paired_t(
    paired_t_exp = 3, paired_t_nexp = 1, n_exp = 50, n_nexp = 50,
    r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6, pre_post_to_smd = "morris_dz"
  )
  all_drm <- es_from_paired_t(
    paired_t_exp = 3, paired_t_nexp = 1, n_exp = 50, n_nexp = 50,
    r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6, pre_post_to_smd = "morris_drm"
  )
  # row 1 must follow dz, row 2 must follow drm
  expect_equal(mixed$d[1], all_dz$d, tolerance = 1e-10)
  expect_equal(mixed$d[2], all_drm$d, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(mixed$d[1], mixed$d[2])))
})

test_that("a scalar pre_post_to_smd still vectorizes over rows", {
  res <- es_from_paired_t(
    paired_t_exp = c(3, 2), paired_t_nexp = c(1, 1),
    n_exp = c(50, 40), n_nexp = c(50, 40),
    r_pre_post_exp = c(0.6, 0.6), r_pre_post_nexp = c(0.6, 0.6)
  )
  expect_equal(nrow(res), 2)
  expect_true(all(is.finite(res$d)))
})

# -----------------------------------------------------------------------------
# 5. paired-t morris_dz variance follows the metafor SMCC convention
# -----------------------------------------------------------------------------
test_that("paired_t_single_group morris_dz matches escalc(SMCC, ti=)", {
  skip_if_not_installed("metafor")
  for (n in c(15, 50, 120)) {
    t_stat <- 3.2
    es <- es_from_paired_t_single_group(
      paired_t_exp = t_stat, n_exp = n, r_pre_post_exp = 0.6,
      pre_post_to_smd = "morris_dz"
    )
    ref <- metafor::escalc(measure = "SMCC", ti = t_stat, ni = n, ri = 0.6)
    expect_equal(es$g, as.numeric(ref$yi), tolerance = 1e-8)
    expect_equal(es$g_se, sqrt(as.numeric(ref$vi)), tolerance = 1e-8)
  }
})

test_that("the paired-t and mean-change routes agree on equivalent inputs (morris_dz)", {
  # Same study, two input formats. The SE convention must not depend on which
  # column the user happened to have.
  n <- 60; mean_change <- 4.5; sd_change <- 9.0
  t_stat <- mean_change / (sd_change / sqrt(n))

  via_t <- es_from_paired_t_single_group(
    paired_t_exp = t_stat, n_exp = n, r_pre_post_exp = 0.6,
    pre_post_to_smd = "morris_dz"
  )
  via_mc <- es_from_mean_change_sd_single_group(
    mean_change_exp = mean_change, mean_change_sd_exp = sd_change,
    n_exp = n, r_pre_post_exp = 0.6, pre_post_to_smd = "morris_dz"
  )
  expect_equal(via_t$g, via_mc$g, tolerance = 1e-8)
  expect_equal(via_t$g_se, via_mc$g_se, tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
# 6. convert_df: pool_sd threading, and the coercion / r-imputation messages
# -----------------------------------------------------------------------------
dat_change <- data.frame(
  study_id = c("s1", "s2"),
  n_exp = c(29, 39), n_nexp = c(34, 38),
  mean_change_exp = c(-12.4, -14.9), mean_change_sd_exp = c(6.52, 6.70),
  mean_change_nexp = c(-3.5, -7.4), mean_change_sd_nexp = c(10.67, 9.70),
  r_pre_post_exp = c(0.6, 0.6), r_pre_post_nexp = c(0.6, 0.6)
)

test_that("convert_df threads pool_sd to the mean-change route (not just accepted and dropped)", {
  pooled <- summary(convert_df(
    dat_change, measure = "g", main_es = TRUE, verbose = FALSE,
    hierarchy = "mean_change_sd", pool_sd = TRUE
  ), digits = 11)
  legacy <- summary(convert_df(
    dat_change, measure = "g", main_es = TRUE, verbose = FALSE,
    hierarchy = "mean_change_sd", pool_sd = FALSE
  ), digits = 11)

  expect_equal(unique(pooled$info_used_crude), "mean_change_sd")
  # arm SDs are unequal in both rows, so the two constructions MUST differ
  expect_false(isTRUE(all.equal(pooled$es_crude, legacy$es_crude, tolerance = 1e-3)))

  # and the pipeline result must equal the wrapper result (no silent drop)
  direct <- es_from_mean_change_sd(
    n_exp = dat_change$n_exp, n_nexp = dat_change$n_nexp,
    mean_change_exp = dat_change$mean_change_exp,
    mean_change_sd_exp = dat_change$mean_change_sd_exp,
    mean_change_nexp = dat_change$mean_change_nexp,
    mean_change_sd_nexp = dat_change$mean_change_sd_nexp,
    r_pre_post_exp = dat_change$r_pre_post_exp,
    r_pre_post_nexp = dat_change$r_pre_post_nexp,
    pre_post_to_smd = "cooper", pool_sd = TRUE
  )
  expect_equal(as.numeric(pooled$es_crude), direct$g, tolerance = 1e-8)
})

test_that("convert_df messages when bonett/morris_dav is coerced to cooper for change data", {
  expect_message(
    convert_df(dat_change, measure = "g", main_es = TRUE, verbose = TRUE,
               hierarchy = "mean_change_sd", pre_post_to_smd = "bonett"),
    "requires separate pre/post SDs"
  )
  # ... and stays silent when the requested method is actually usable
  msgs <- capture.output(
    invisible(convert_df(dat_change, measure = "g", main_es = TRUE, verbose = TRUE,
                         hierarchy = "mean_change_sd", pre_post_to_smd = "morris_dz")),
    type = "message"
  )
  expect_false(any(grepl("requires separate pre/post SDs", msgs)))
})

test_that("convert_df messages when r_pre_post is imputed for rows that use it", {
  no_r <- dat_change
  no_r$r_pre_post_exp <- NA
  no_r$r_pre_post_nexp <- NA
  expect_message(
    convert_df(no_r, measure = "g", main_es = TRUE, verbose = TRUE,
               hierarchy = "mean_change_sd"),
    "r_pre_post was not reported"
  )
})

test_that("the r-imputation note does NOT fire on data with no pre/post columns", {
  # Guard against the message becoming ambient noise: df.haza carries no
  # pre/post, mean-change or paired data, so no row can consume r_pre_post.
  msgs <- capture.output(
    invisible(convert_df(df.haza, measure = "g", main_es = TRUE, verbose = TRUE)),
    type = "message"
  )
  expect_false(any(grepl("r_pre_post was not reported", msgs)))
  expect_false(any(grepl("requires separate pre/post SDs", msgs)))
})

# -----------------------------------------------------------------------------
# 7. E6 / E7 estimand-mixing flags
# -----------------------------------------------------------------------------
test_that("E6 flags a pool mixing change-SD (morris_dz) rows with endpoint rows", {
  mixed <- data.frame(
    study_id = c("endpoint1", "endpoint2", "change1"),
    n_exp = c(30, 32, 29), n_nexp = c(30, 31, 34),
    # two endpoint (means_sd) rows -> raw-score-SD metric
    mean_exp = c(12, 13, NA), mean_sd_exp = c(4, 4.2, NA),
    mean_nexp = c(10, 10.5, NA), mean_sd_nexp = c(4.1, 4.0, NA),
    # one change-score row -> change-SD metric under morris_dz
    mean_change_exp = c(NA, NA, -12.4), mean_change_sd_exp = c(NA, NA, 6.52),
    mean_change_nexp = c(NA, NA, -3.5), mean_change_sd_nexp = c(NA, NA, 10.67),
    r_pre_post_exp = c(NA, NA, 0.6), r_pre_post_nexp = c(NA, NA, 0.6)
  )
  res <- summary(convert_df(
    mixed, measure = "g", main_es = TRUE, verbose = FALSE,
    hierarchy = "means_sd > mean_change_sd", pre_post_to_smd = "morris_dz"
  ), flags = TRUE)

  expect_true("flags_crude" %in% names(res))
  expect_true(any(grepl("Mixed SMD standardizers", res$flags_crude)))
  # the change-score row must be named
  chg <- which(res$info_used_crude == "mean_change_sd")
  expect_true(all(grepl("Mixed SMD standardizers", res$flags_crude[chg])))
})

test_that("E6 stays silent when every row is on the raw-score-SD metric", {
  # Same data, but morris_drm puts the change row back on the raw-score metric.
  mixed <- data.frame(
    study_id = c("endpoint1", "change1"),
    n_exp = c(30, 29), n_nexp = c(30, 34),
    mean_exp = c(12, NA), mean_sd_exp = c(4, NA),
    mean_nexp = c(10, NA), mean_sd_nexp = c(4.1, NA),
    mean_change_exp = c(NA, -12.4), mean_change_sd_exp = c(NA, 6.52),
    mean_change_nexp = c(NA, -3.5), mean_change_sd_nexp = c(NA, 10.67),
    r_pre_post_exp = c(NA, 0.6), r_pre_post_nexp = c(NA, 0.6)
  )
  res <- summary(convert_df(
    mixed, measure = "g", main_es = TRUE, verbose = FALSE,
    hierarchy = "means_sd > mean_change_sd", pre_post_to_smd = "morris_drm"
  ), flags = TRUE)
  expect_false(any(grepl("Mixed SMD standardizers", res$flags_crude)))
})

test_that("E7 flags paired-t rows sharing a pool with pooled-standardizer rows", {
  mixed <- data.frame(
    study_id = c("mc1", "pt1"),
    n_exp = c(29, 40), n_nexp = c(34, 40),
    mean_change_exp = c(-12.4, NA), mean_change_sd_exp = c(6.52, NA),
    mean_change_nexp = c(-3.5, NA), mean_change_sd_nexp = c(10.67, NA),
    paired_t_exp = c(NA, 3.1), paired_t_nexp = c(NA, 1.2),
    r_pre_post_exp = c(0.6, 0.6), r_pre_post_nexp = c(0.6, 0.6)
  )
  res <- summary(convert_df(
    mixed, measure = "g", main_es = TRUE, verbose = FALSE,
    hierarchy = "mean_change_sd > paired_t", pool_sd = TRUE,
    flag_options = list(enable_informational = TRUE)
  ), flags = TRUE)

  pt <- which(res$info_used_crude == "paired_t")
  expect_true(length(pt) > 0)
  expect_true(all(grepl("Per-arm standardizer", res$flags_crude[pt])))
})
