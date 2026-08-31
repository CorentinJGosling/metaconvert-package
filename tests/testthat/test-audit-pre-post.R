# Red tests for AUDIT-2026-08-28-findings.md -- the pre/post and paired family.
#
# Every test_that() below asserts the CORRECT behaviour, so it is expected to FAIL
# until the corresponding defect is repaired. Line numbers in the comments refer to
# AUDIT-2026-08-28-findings.md. Where a Verifier correction is present it overrides
# the filed Fix and the assertions follow the correction.

# --------------------------------------------------------------------------------
# AUDIT 456-465 -- es_from_means_sd_pre_post_single_group() computes the mdw variance
# inline from the RAW signed SDs (R/es_from_PAIRED_SINGLE_GROUP.R:146), bypassing the
# guard the SMD kernel applies. A negative SD sign-flips the -2*r*sd_pre*sd_post cross
# term, so dw/gw are correctly NA'd while mdw_se / mdw_ci are emitted finite and wrong
# (here 1.665833 instead of 2.876630, a ~3x inverse-variance weight error), and the CI
# flips the row from non-significant to significant.
# --------------------------------------------------------------------------------
test_that("AUDIT-456: the single-group mdw variance guards a negative SD instead of sign-flipping the cross term", {
  # Correct value, from the definition of Var(post - pre) for paired data:
  #   sqrt((sd_pre^2 + sd_post^2 - 2*r*sd_pre*sd_post) / n)
  clean_se <- sqrt((10^2 + 11^2 - 2 * (-0.5) * 10 * 11) / 40) # 2.8766300
  ok <- es_from_means_sd_pre_post_single_group(50, 55, 10, 11, 40, -0.5)
  expect_equal(ok$mdw, 5)
  expect_equal(ok$mdw_se, clean_se, tolerance = 1e-7)

  # A negative BASELINE SD must neutralise the mdw variance, exactly as
  # .single_group_pre_post_to_smd() already neutralises dw/gw.
  bad_pre <- es_from_means_sd_pre_post_single_group(50, 55, -10, 11, 40, -0.5)
  expect_true(is.na(bad_pre$d)) # already true today (kernel guard)
  expect_true(is.na(bad_pre$mdw_se))
  expect_true(is.na(bad_pre$mdw_ci_lo))
  expect_true(is.na(bad_pre$mdw_ci_up))

  # ... and so must a negative ENDPOINT SD.
  bad_post <- es_from_means_sd_pre_post_single_group(50, 55, 10, -11, 40, -0.5)
  expect_true(is.na(bad_post$mdw_se))
  expect_true(is.na(bad_post$mdw_ci_lo))
  expect_true(is.na(bad_post$mdw_ci_up))

  # Regression guard on the fix: the mean-change wrappers legitimately pass
  # mean_pre_sd_exp = 0, so the guard must be ">= 0", not .positive_or_na().
  mc <- es_from_mean_change_sd_single_group(
    mean_change_exp = 2, mean_change_sd_exp = 5, n_exp = 40, r_pre_post_exp = 0.5
  )
  expect_equal(mc$mdw_se, sqrt(25 / 40), tolerance = 1e-7) # 0.7905694
})

# --------------------------------------------------------------------------------
# AUDIT 866-877 (+ Verifier correction) -- four pre/post-means entry points bypass the
# R/internal_guards.R helpers: es_from_means_se_pre_post (sd <- se * sqrt(n), no
# .positive_or_na()), es_from_means_ci_pre_post (raw ci_up - ci_lo instead of
# .ci_width()), and the two single-group analogues. A negative SE or a transposed CI
# bound therefore inflates md_se / mdw_se by 2.24x / 3.00x with no NA and no warning.
# --------------------------------------------------------------------------------
test_that("AUDIT-866: pre/post SE and CI entry points neutralise a negative SE and a transposed CI", {
  se <- 2 / sqrt(30) # SD = 2 in every arm at every timepoint
  # md_se = sqrt(sum over arms of (sd_pre^2 + sd_post^2 - 2*r*sd_pre*sd_post)/n)
  clean_md_se <- sqrt(2 * (2^2 + 2^2 - 2 * 0.8 * 2 * 2) / 30) # 0.3265986

  ok <- es_from_means_se_pre_post(
    mean_pre_exp = 10, mean_exp = 13, mean_pre_se_exp = se, mean_se_exp = se,
    mean_pre_nexp = 10, mean_nexp = 11, mean_pre_se_nexp = se, mean_se_nexp = se,
    n_exp = 30, n_nexp = 30, r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8
  )
  expect_equal(ok$md_se, clean_md_se, tolerance = 1e-7)

  # A negative reported SE is arithmetically impossible: NA it, do not square it.
  bad_se <- es_from_means_se_pre_post(
    mean_pre_exp = 10, mean_exp = 13, mean_pre_se_exp = se, mean_se_exp = -se,
    mean_pre_nexp = 10, mean_nexp = 11, mean_pre_se_nexp = se, mean_se_nexp = se,
    n_exp = 30, n_nexp = 30, r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8
  )
  expect_true(is.na(bad_se$md_se))

  # A transposed CI is the SAME interval (.ci_width() = abs(up - lo)), so the route
  # must return the clean numbers, not a 2.24x-wide SE.
  tq <- qt(.975, 29)
  pre_lo <- 10 - tq * se; pre_up <- 10 + tq * se
  exp_lo <- 13 - tq * se; exp_up <- 13 + tq * se
  nex_lo <- 11 - tq * se; nex_up <- 11 + tq * se
  ci_args <- function(lo, up) {
    es_from_means_ci_pre_post(
      mean_pre_exp = 10, mean_pre_ci_lo_exp = pre_lo, mean_pre_ci_up_exp = pre_up,
      mean_exp = 13, mean_ci_lo_exp = lo, mean_ci_up_exp = up,
      mean_pre_nexp = 10, mean_pre_ci_lo_nexp = pre_lo, mean_pre_ci_up_nexp = pre_up,
      mean_nexp = 11, mean_ci_lo_nexp = nex_lo, mean_ci_up_nexp = nex_up,
      n_exp = 30, n_nexp = 30, r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8
    )
  }
  clean_ci <- ci_args(exp_lo, exp_up)
  transposed <- ci_args(exp_up, exp_lo)
  expect_equal(transposed$md_se, clean_md_se, tolerance = 1e-7)
  expect_equal(transposed$md_se, clean_ci$md_se, tolerance = 1e-10)
  expect_equal(transposed$md_ci_lo, clean_ci$md_ci_lo, tolerance = 1e-10)
  expect_equal(transposed$md_ci_up, clean_ci$md_ci_up, tolerance = 1e-10)

  # Single-group analogues (R/es_from_PAIRED_SINGLE_GROUP.R:249-250 and :354-355).
  clean_mdw_se <- sqrt((2^2 + 2^2 - 2 * 0.8 * 2 * 2) / 30) # 0.2309401
  sg_ok <- es_from_means_se_pre_post_single_group(10, 13, se, se, 30, 0.8)
  expect_equal(sg_ok$mdw_se, clean_mdw_se, tolerance = 1e-7)
  sg_bad <- es_from_means_se_pre_post_single_group(10, 13, se, -se, 30, 0.8)
  expect_true(is.na(sg_bad$mdw_se))

  sg_ci <- function(lo, up) {
    es_from_means_ci_pre_post_single_group(
      mean_pre_exp = 10, mean_pre_ci_lo_exp = pre_lo, mean_pre_ci_up_exp = pre_up,
      mean_exp = 13, mean_ci_lo_exp = lo, mean_ci_up_exp = up,
      n_exp = 30, r_pre_post_exp = 0.8
    )
  }
  expect_equal(sg_ci(exp_up, exp_lo)$mdw_se, clean_mdw_se, tolerance = 1e-7)

  # Regression guard on the fix (Verifier correction): es_from_mean_change_se() passes
  # mean_pre_se = 0 deliberately, so a blanket .positive_or_na() would wipe the whole
  # mean-change family. These must stay finite.
  mc <- es_from_mean_change_se(
    mean_change_exp = 3, mean_change_se_exp = se,
    mean_change_nexp = 1, mean_change_se_nexp = se,
    n_exp = 30, n_nexp = 30, r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8
  )
  expect_equal(mc$md_se, sqrt(2 * 2^2 / 30), tolerance = 1e-7) # 0.5163978
  expect_false(is.na(mc$d))
})

# --------------------------------------------------------------------------------
# AUDIT 892-903 -- the mdw block of es_from_means_sd_pre_post_single_group() indexes
# the RAW arguments with nn_miss, while the SMD block above it reads a data.frame()
# that already recycled length-1 columns. An explicitly scalar r_pre_post_exp (or
# n_exp / mean_pre_sd_exp / mean_sd_exp) therefore NAs mdw_se and the mdw CI on every
# row but the first, leaving a fully populated mdw column that rma() silently thins.
# --------------------------------------------------------------------------------
test_that("AUDIT-892: the single-group mdw block recycles length-1 arguments before indexing", {
  vec <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = c(10, 20, 30), mean_exp = c(12, 23, 34),
    mean_pre_sd_exp = c(5, 5, 5), mean_sd_exp = c(6, 6, 6),
    n_exp = c(40, 40, 40), r_pre_post_exp = rep(0.5, 3), pre_post_to_smd = "bonett"
  )
  # Closed form: sqrt((25 + 36 - 2*0.5*5*6) / 40) = sqrt(31/40) = 0.8803408 on all rows.
  expect_equal(vec$mdw_se, rep(sqrt((25 + 36 - 2 * 0.5 * 5 * 6) / 40), 3), tolerance = 1e-7)

  scalar_r <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = c(10, 20, 30), mean_exp = c(12, 23, 34),
    mean_pre_sd_exp = c(5, 5, 5), mean_sd_exp = c(6, 6, 6),
    n_exp = c(40, 40, 40), r_pre_post_exp = 0.5, pre_post_to_smd = "bonett"
  )
  expect_equal(scalar_r$mdw_se, vec$mdw_se, tolerance = 1e-10)
  expect_equal(scalar_r$mdw_ci_lo, vec$mdw_ci_lo, tolerance = 1e-10)
  expect_equal(scalar_r$mdw_ci_up, vec$mdw_ci_up, tolerance = 1e-10)

  scalar_n <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = c(10, 20, 30), mean_exp = c(12, 23, 34),
    mean_pre_sd_exp = c(5, 5, 5), mean_sd_exp = c(6, 6, 6),
    n_exp = 40, r_pre_post_exp = rep(0.5, 3), pre_post_to_smd = "bonett"
  )
  expect_equal(scalar_n$mdw_se, vec$mdw_se, tolerance = 1e-10)

  # It propagates through the mean-change wrapper: sqrt(25/40) = 0.7905694 on all rows.
  mc <- es_from_mean_change_sd_single_group(
    mean_change_exp = c(2, 3, 4), mean_change_sd_exp = c(5, 5, 5),
    n_exp = c(40, 40, 40), r_pre_post_exp = 0.5
  )
  expect_equal(mc$mdw_se, rep(sqrt(25 / 40), 3), tolerance = 1e-7)
})

# --------------------------------------------------------------------------------
# AUDIT 827-838 (regression from commit 1eae59c) -- the paired-t/F routes now reach
# .single_group_pre_post_to_smd() with mean_pre_sd = 0, whose morris_dz cross term
# 2*r*sd_pre*sd_post is identically zero but evaluates to NA once .guard_r_pre_post()
# has NA'd an |r| >= 1. d_z = t/sqrt(n) and its SMCC variance 1/n + g^2/(2n) are both
# r-free, so the whole row must be invariant to r; master returned the finite value.
# --------------------------------------------------------------------------------
test_that("AUDIT-827: the paired-t/F family under morris_dz is invariant to r_pre_post", {
  ref <- es_from_paired_t(5.1, 2.2, 30, 28, 0.6, 0.6, pre_post_to_smd = "morris_dz")
  # Independent closed form for d_z: t_exp/sqrt(n_exp) - t_nexp/sqrt(n_nexp).
  expect_equal(ref$d, 5.1 / sqrt(30) - 2.2 / sqrt(28), tolerance = 1e-7) # 0.5153673

  for (rv in c(1, -1, 1.5)) {
    got <- es_from_paired_t(5.1, 2.2, 30, 28, rv, rv, pre_post_to_smd = "morris_dz")
    expect_equal(got$d, 5.1 / sqrt(30) - 2.2 / sqrt(28), tolerance = 1e-7)
    expect_equal(got$d_se, ref$d_se, tolerance = 1e-10)
    expect_equal(got$g, ref$g, tolerance = 1e-10)
    expect_equal(got$g_se, ref$g_se, tolerance = 1e-10)
  }

  # Same delegation, same defect, via the F entry point.
  f_ref <- es_from_paired_f(5.1^2, 2.2^2, 30, 28, 0.6, 0.6, pre_post_to_smd = "morris_dz")
  f_bad <- es_from_paired_f(5.1^2, 2.2^2, 30, 28, 1, 1, pre_post_to_smd = "morris_dz")
  expect_equal(f_bad$d, f_ref$d, tolerance = 1e-10)
  expect_equal(f_bad$g_se, f_ref$g_se, tolerance = 1e-10)

  # ... and via the single-group sibling: d_z = t/sqrt(n) = 2.1/sqrt(20) = 0.4695743.
  sg_ref <- es_from_paired_t_single_group(2.1, 20, 0.8, pre_post_to_smd = "morris_dz")
  expect_equal(sg_ref$d, 2.1 / sqrt(20), tolerance = 1e-7)
  sg_bad <- es_from_paired_t_single_group(2.1, 20, 1, pre_post_to_smd = "morris_dz")
  expect_equal(sg_bad$d, 2.1 / sqrt(20), tolerance = 1e-7)
  expect_equal(sg_bad$d_se, sg_ref$d_se, tolerance = 1e-10)
})

# --------------------------------------------------------------------------------
# AUDIT 905-916 (Verifier correction: the fix belongs in the kernel, because the
# mean-change family reaches the same 0 * NA poisoning and does so on master too) --
# whenever a caller zeroes the pre slot, .single_group_pre_post_to_smd()'s morris_dz
# standardizer must reduce to sd_post for ANY r, including an NA'd out-of-range one.
# --------------------------------------------------------------------------------
test_that("AUDIT-905: the morris_dz kernel and the mean-change family are r-free when the pre SD is zero", {
  k_ref <- metaConvert:::.single_group_pre_post_to_smd(0, 2, 0, 1, 25, 0.5, "morris_dz")
  # Var(post - pre) = Var(post) when Var(pre) = 0, so d = (2 - 0)/1 = 2 for any r.
  expect_equal(unname(k_ref[1, "d"]), 2, tolerance = 1e-10)
  k_na <- metaConvert:::.single_group_pre_post_to_smd(0, 2, 0, 1, 25, NA_real_, "morris_dz")
  expect_equal(unname(k_na[1, ]), unname(k_ref[1, ]), tolerance = 1e-10)

  # Two-group mean-change route (pre-existing, not a 1eae59c regression).
  mc <- function(rv) {
    es_from_mean_change_sd(
      mean_change_exp = 5, mean_change_sd_exp = 8,
      mean_change_nexp = 2, mean_change_sd_nexp = 7,
      n_exp = 30, n_nexp = 32, r_pre_post_exp = rv, r_pre_post_nexp = rv,
      pre_post_to_smd = "morris_dz"
    )
  }
  mc_ref <- mc(0.6)
  expect_equal(mc_ref$d, 5 / 8 - 2 / 7, tolerance = 1e-7) # 0.3392857
  mc_bad <- mc(1)
  expect_equal(mc_bad$d, 5 / 8 - 2 / 7, tolerance = 1e-7)
  expect_equal(mc_bad$d_se, mc_ref$d_se, tolerance = 1e-10)
  expect_equal(mc_bad$g, mc_ref$g, tolerance = 1e-10)
  expect_equal(mc_bad$g_se, mc_ref$g_se, tolerance = 1e-10)

  # Single-group mean-change route: d_z = mean_change / sd_change = 5/8.
  sg_bad <- es_from_mean_change_sd_single_group(
    mean_change_exp = 5, mean_change_sd_exp = 8, n_exp = 30,
    r_pre_post_exp = 1, pre_post_to_smd = "morris_dz"
  )
  expect_equal(sg_bad$d, 5 / 8, tolerance = 1e-7)

  # Regression guard: morris_drm legitimately stays NA (its sqrt(2(1-r)) factor).
  drm_bad <- es_from_mean_change_sd(
    mean_change_exp = 5, mean_change_sd_exp = 8,
    mean_change_nexp = 2, mean_change_sd_nexp = 7,
    n_exp = 30, n_nexp = 32, r_pre_post_exp = 1, r_pre_post_nexp = 1,
    pre_post_to_smd = "morris_drm"
  )
  expect_true(is.na(drm_bad$d))
})

# --------------------------------------------------------------------------------
# AUDIT 853-864 (Verifier correction) -- an F statistic and a two-sided p-value carry
# no sign, so es_from_paired_f / _f_pval / _t_pval recover |t| for BOTH arms and then
# form d <- d_exp - d_nexp. When the control arm moved the other way the correct
# contrast is |d_exp| + |d_nexp|, and no per-arm sign column exists, so the direction
# cannot be expressed through these entry points at all. Assertions below pin the
# behaviourally testable half of the corrected fix (per-arm sign columns); the
# documentation half is not assertable from a portable test file.
# --------------------------------------------------------------------------------
test_that("AUDIT-853: the unsigned paired-F / paired-p routes can express opposite-direction arms", {
  n <- 30
  r <- 0.6
  sd_change <- 2
  sigma_raw <- sd_change / sqrt(2 * (1 - r)) # 2.236068
  t_exp <- 3 / (sd_change / sqrt(n)) # +8.215838 (improvement)
  t_nexp <- -1 / (sd_change / sqrt(n)) # -2.738613 (deterioration)
  true_drm <- (3 - (-1)) / sigma_raw # 1.788854

  # The signed entry point already gets it right -- regression guard.
  expect_equal(
    es_from_paired_t(t_exp, t_nexp, n, n, r, r, pre_post_to_smd = "cooper")$d,
    true_drm,
    tolerance = 1e-6
  )
  expect_equal(
    es_from_mean_change_pval(3, 2 * pt(-abs(t_exp), n - 1), -1,
                             2 * pt(-abs(t_nexp), n - 1), r, r, n, n,
                             pre_post_to_smd = "cooper")$d,
    true_drm,
    tolerance = 1e-6
  )

  p_exp <- 2 * pt(-abs(t_exp), n - 1)
  p_nexp <- 2 * pt(-abs(t_nexp), n - 1)

  expect_equal(
    es_from_paired_f(t_exp^2, t_nexp^2, n, n, r, r,
                     pre_post_to_smd = "cooper", reverse_paired_f_nexp = TRUE)$d,
    true_drm,
    tolerance = 1e-6
  )
  expect_equal(
    es_from_paired_f_pval(p_exp, p_nexp, n, n, r, r,
                          pre_post_to_smd = "cooper",
                          reverse_paired_f_pval_nexp = TRUE)$d,
    true_drm,
    tolerance = 1e-6
  )
  expect_equal(
    es_from_paired_t_pval(p_exp, p_nexp, n, n, r, r,
                          pre_post_to_smd = "cooper",
                          reverse_paired_t_pval_nexp = TRUE)$d,
    true_drm,
    tolerance = 1e-6
  )
})

# --------------------------------------------------------------------------------
# AUDIT 467-476 -- pre_post_to_smd_restricted (R/main_convert_df.R:857) is computed
# OUTSIDE with(), from the length-1 formal argument, so the per-row pre_post_to_smd
# COLUMN is discarded for all 8 mean-change routes and all 5 paired-t/F routes. Every
# row silently gets whichever method the scalar implies (by default bonett -> cooper).
# --------------------------------------------------------------------------------
test_that("AUDIT-467: convert_df honours a per-row pre_post_to_smd column on the paired-t and mean-change routes", {
  # Scalar reference runs (the two values the column asks for), taken from the
  # exported routes, which do honour a per-row method.
  ref_dz <- es_from_paired_t(rep(3.5, 2), rep(1.2, 2), rep(30, 2), rep(32, 2),
                             rep(0.6, 2), rep(0.6, 2),
                             pre_post_to_smd = c("morris_dz", "morris_dz"))
  ref_drm <- es_from_paired_t(rep(3.5, 2), rep(1.2, 2), rep(30, 2), rep(32, 2),
                              rep(0.6, 2), rep(0.6, 2),
                              pre_post_to_smd = c("morris_drm", "morris_drm"))
  expect_equal(ref_dz$g[1], 0.4153645, tolerance = 1e-6)
  expect_equal(ref_drm$g[1], 0.3715133, tolerance = 1e-6)

  dat <- data.frame(
    study_id = c("A", "B"), paired_t_exp = 3.5, paired_t_nexp = 1.2,
    n_exp = 30, n_nexp = 32, r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6,
    pre_post_to_smd = c("morris_dz", "morris_drm")
  )
  out <- NULL
  invisible(utils::capture.output(suppressMessages(
    out <- as.data.frame(summary(
      convert_df(dat, measure = "g", verbose = FALSE, split_adjusted = FALSE)
    ))
  )))
  a <- out[match("A", out$study_id), ]
  b <- out[match("B", out$study_id), ]
  expect_equal(a$es, ref_dz$g[1], tolerance = 1e-6)
  expect_equal(a$se, ref_dz$g_se[1], tolerance = 1e-6)
  expect_equal(b$es, ref_drm$g[1], tolerance = 1e-6)
  expect_equal(b$se, ref_drm$g_se[1], tolerance = 1e-6)

  # Same defect on the mean-change family.
  dat2 <- data.frame(
    study_id = c("A", "B"), mean_change_exp = 5, mean_change_sd_exp = 8,
    mean_change_nexp = 2, mean_change_sd_nexp = 7,
    n_exp = 30, n_nexp = 32, r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6,
    pre_post_to_smd = c("morris_dz", "morris_drm")
  )
  out2 <- NULL
  invisible(utils::capture.output(suppressMessages(
    out2 <- as.data.frame(summary(
      convert_df(dat2, measure = "g", verbose = FALSE, split_adjusted = FALSE)
    ))
  )))
  # Scalar reference values through the same pipeline.
  expect_equal(out2$es[match("A", out2$study_id)], 0.3299357, tolerance = 1e-6) # morris_dz
  expect_equal(out2$se[match("A", out2$study_id)], 0.2682760, tolerance = 1e-6)
  expect_equal(out2$es[match("B", out2$study_id)], 0.2951035, tolerance = 1e-6) # morris_drm
  expect_equal(out2$se[match("B", out2$study_id)], 0.2345127, tolerance = 1e-6)
})

# --------------------------------------------------------------------------------
# AUDIT 879-890 (Verifier correction: the same OR is duplicated at
# R/main_convert_df.R:409 for the verbose note, and the trigger is any row whose
# r-consuming data sits in one arm only) -- .r_defaulted ORs both arms, so a
# single-group row is always "defaulted" through its structurally absent nexp arm,
# even though r_pre_post_exp WAS supplied and IS what the arithmetic used.
# --------------------------------------------------------------------------------
test_that("AUDIT-879: V6 / r_defaulted does not fire on a single-group row that supplied r_pre_post_exp", {
  x <- data.frame(
    study_id = c("s1", "s2"),
    mean_pre_exp = c(50, 50), mean_exp = c(55, 56),
    mean_pre_sd_exp = c(10, 10), mean_sd_exp = c(11, 11),
    n_exp = c(40, 40), r_pre_post_exp = c(0.6, 0.6)
  )
  res <- convert_df(x, measure = "dw", verbose = FALSE)
  expect_equal(unname(attr(res, "r_defaulted")), c(FALSE, FALSE))

  s <- NULL
  invisible(utils::capture.output(suppressMessages(
    s <- as.data.frame(summary(res, flags = TRUE))
  )))
  expect_false(any(grepl("Default r_pre_post", s$flags_crude, fixed = TRUE)))

  # The verbose console note (R/main_convert_df.R:409) must agree with the flag.
  msgs <- utils::capture.output(
    invisible(convert_df(x, measure = "dw", verbose = TRUE)),
    type = "message"
  )
  expect_false(any(grepl("r_pre_post was not reported", msgs, fixed = TRUE)))

  # Positive control, unchanged: a row that genuinely omits r IS flagged.
  x_no_r <- x
  x_no_r$r_pre_post_exp <- NULL
  expect_true(all(attr(convert_df(x_no_r, measure = "dw", verbose = FALSE), "r_defaulted")))
})

# --------------------------------------------------------------------------------
# AUDIT 840-851 (Verifier correction: the `ok` widening is verified output-identical,
# and is the half that recovers the cost on datasets with no paired data) --
# .paired_t_to_smd()'s filter keys only on n >= 2, so every row enters the per-row
# mapply even when its paired statistic is NA. Asserted by counting kernel calls
# rather than by timing.
# --------------------------------------------------------------------------------
test_that("AUDIT-840: .paired_t_to_smd skips rows carrying no paired statistic", {
  real_kernel <- metaConvert:::.single_group_pre_post_to_smd
  calls <- 0L
  local_mocked_bindings(
    .single_group_pre_post_to_smd = function(...) {
      calls <<- calls + 1L
      real_kernel(...)
    },
    .package = "metaConvert"
  )

  k <- 40L
  na_stat <- rep(NA_real_, k)
  nvec <- rep(30, k)
  rvec <- rep(0.8, k)

  out <- es_from_paired_t(na_stat, na_stat, nvec, nvec, rvec, rvec)
  # No paired statistic anywhere: the kernel must never be entered.
  expect_identical(calls, 0L)
  # ... and the output must be unchanged by the skip (verified output-identical).
  expect_true(all(is.na(out$d)))
  expect_true(all(is.na(out$d_se)))
  expect_true(all(is.na(out$g)))

  # A pool that DOES carry paired data still reaches the kernel and is unaffected.
  calls <- 0L
  mixed <- es_from_paired_t(c(NA_real_, 2.5), c(NA_real_, 1.1), c(30, 30), c(30, 30),
                            c(0.8, 0.8), c(0.8, 0.8))
  expect_true(calls > 0L)
  expect_true(is.na(mixed$d[1]))
  expect_false(is.na(mixed$d[2]))

  # The .d_j() warning stream sits outside the filter and must stay unchanged.
  calls <- 0L
  seen <- character()
  withCallingHandlers(
    es_from_paired_t(na_stat, na_stat, c(1, rep(30, k - 1)), c(1, rep(30, k - 1)),
                     rvec, rvec),
    warning = function(cnd) {
      seen <<- c(seen, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("Hedges' J correction is undefined", seen, fixed = TRUE)))
})
