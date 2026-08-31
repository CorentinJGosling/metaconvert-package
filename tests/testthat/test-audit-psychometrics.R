# Red-half regression tests for the 2026-08-28 audit, psychometric-utils subsystem.
#
# Every test in this file is expected to FAIL until the corresponding defect in
# R/psychometric_utils.R is repaired. Each test asserts the CORRECT value,
# derived by hand from the formulas the function's own @details document, never
# merely "not the buggy value".

# Locate the package source root (needed by the documentation test only).
# Returns NA when the sources are not beside the tests (e.g. R CMD check, which
# copies tests/ but not R/ or man/), in which case that test skips.
.audit_pkg_root <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:6) {
    desc <- file.path(d, "DESCRIPTION")
    if (file.exists(desc) &&
        any(grepl("^Package:[[:space:]]*metaConvert", readLines(desc, warn = FALSE)))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  NA_character_
}


# ---------------------------------------------------------------------------
# AUDIT line 522-531 — compute_sem(): a length-1 `icc_se` (or any argument
# shorter than sd/icc) is not recycled, so the nested ifelse() collapses to
# length 1 and ROW 1's variance and CI are broadcast onto every row, producing
# 95% CIs that do not contain their own point estimate.
# ---------------------------------------------------------------------------
test_that("AUDIT-psychometrics: compute_sem recycles a scalar icc_se instead of broadcasting row 1", {

  # ---- regime A: external icc_se supplied (delta method) -------------------
  # Oracle, straight from the documented formula (?compute_sem @details case 1):
  #   Var(SEM) = SD^2 / (4 (1 - ICC)) * icc_se^2 + (1 - ICC) * SD^2 / (2 (n - 1))
  #   CI = SEM +/- qnorm(.975) * SE, truncated at 0
  # with sd = c(10, 20, 5), icc = c(.85, .60, .90), n = 100, icc_se = 0.03:
  #   row 1  sem  3.872983346  se 0.4751395329  ci [ 2.941726974,  4.804239718]
  #   row 2  sem 12.649110641  se 1.0164058284  ci [10.656991823, 14.641229458]
  #   row 3  sem  1.581138830  se 0.2624428750  ci [ 1.066760247,  2.095517413]
  res <- compute_sem(sd = c(10, 20, 5), icc = c(0.85, 0.60, 0.90),
                     n_sample = 100, icc_se = 0.03)

  expect_equal(res$sem, c(3.872983346, 12.649110641, 1.581138830), tolerance = 1e-7)
  expect_equal(res$sem_se, c(0.4751395329, 1.0164058284, 0.2624428750), tolerance = 1e-7)
  expect_equal(res$sem_ci_lo, c(2.941726974, 10.656991823, 1.066760247), tolerance = 1e-7)
  expect_equal(res$sem_ci_up, c(4.804239718, 14.641229458, 2.095517413), tolerance = 1e-7)

  # A 95% CI must contain its own point estimate. Today row 2's does not
  # (sem = 12.65 with a reported CI of [2.94, 4.80]).
  expect_true(all(res$sem >= res$sem_ci_lo & res$sem <= res$sem_ci_up))

  # Passing the SAME scalar as a full-length vector already gives the right
  # answer, so the defect is purely the recycling, not the formula: the two
  # calls must agree.
  res_vec <- compute_sem(sd = c(10, 20, 5), icc = c(0.85, 0.60, 0.90),
                         n_sample = 100, icc_se = rep(0.03, 3))
  expect_equal(res$sem_se, res_vec$sem_se, tolerance = 1e-10)
  expect_equal(res$sem_ci_lo, res_vec$sem_ci_lo, tolerance = 1e-10)
  expect_equal(res$sem_ci_up, res_vec$sem_ci_up, tolerance = 1e-10)

  # ---- regime B: same-sample chi-square, sd scalar shorter than icc --------
  # Oracle (?compute_sem @details case 2), df = (n - 1)(k - 1) = 99:
  #   SE = SEM / sqrt(2 df);  CI = SEM * sqrt(df / qchisq(c(.975, .025), df))
  # with sd = 10, icc = c(.85, .60, .90), n = 100, k = 2:
  #   row 1  sem 3.872983346  se 0.2752409413  ci [3.400505890, 4.499148917]
  #   row 2  sem 6.324555320  se 0.4494665750  ci [5.553002865, 7.347079416]
  #   row 3  sem 3.162277660  se 0.2247332875  ci [2.776501433, 3.673539708]
  res_ss <- compute_sem(sd = 10, icc = c(0.85, 0.60, 0.90), n_sample = 100)

  expect_equal(res_ss$sem, c(3.872983346, 6.324555320, 3.162277660), tolerance = 1e-7)
  expect_equal(res_ss$sem_se, c(0.2752409413, 0.4494665750, 0.2247332875), tolerance = 1e-7)
  expect_equal(res_ss$sem_ci_lo, c(3.400505890, 5.553002865, 2.776501433), tolerance = 1e-7)
  expect_equal(res_ss$sem_ci_up, c(4.499148917, 7.347079416, 3.673539708), tolerance = 1e-7)
  expect_true(all(res_ss$sem >= res_ss$sem_ci_lo & res_ss$sem <= res_ss$sem_ci_up))

  # ---- partial-length recycling: icc_se length 2 over 4 rows ---------------
  # This is a clean multiple, so nothing errors; today rows 3/4 silently
  # inherit rows 1/2's SE and CI. Oracle from the same case-1 formula with
  # icc_se recycled elementwise, i.e. c(0.03, 0.05, 0.03, 0.05).
  res_p <- compute_sem(sd = c(10, 20, 5, 8), icc = c(0.85, 0.60, 0.90, 0.70),
                       n_sample = 100, icc_se = c(0.03, 0.05))
  expect_equal(res_p$sem_se,
               c(0.4751395329, 1.1971135318, 0.2624428750, 0.4798989793),
               tolerance = 1e-7)
  expect_equal(res_p$sem_ci_lo,
               c(2.941726974, 10.302811233, 1.066760247, 3.441195744),
               tolerance = 1e-7)
  expect_equal(res_p$sem_ci_up,
               c(4.804239718, 14.995410048, 2.095517413, 5.322365176),
               tolerance = 1e-7)
  expect_true(all(res_p$sem >= res_p$sem_ci_lo & res_p$sem <= res_p$sem_ci_up))
})


# ---------------------------------------------------------------------------
# AUDIT line 533-542 — compute_sem(): the same-sample chi-square variance/CI is
# exact only for a consistency ICC(3,1); under the absolute-agreement ICC(2,1)
# that es_from_icc() returns BY DEFAULT the SE is up to 2.7x too small
# (measured coverage 0.80 at n = 30, 0.51 at n = 100) and neither ?compute_sem
# nor man/compute_sem.Rd mentions ICC type at all.
#
# AGREED REMEDY IS DOC-ONLY (audit conflict table, line 195: "Take the SEM
# finding's minimal option (@details + Rd caveat mirroring V31), keep the closed
# form computing by default"). This test therefore asserts (a) the caveat exists
# in both the roxygen block and the CRAN-visible Rd, and (b) the numbers are
# UNCHANGED — no NA on the default path, which would blank two vignette chunks.
# ---------------------------------------------------------------------------
test_that("AUDIT-psychometrics: compute_sem documents the ICC(2,1)-vs-ICC(3,1) caveat and keeps computing", {

  # (b) first: the closed form must keep returning numbers on the default path.
  # Oracle, df = (30 - 1)(2 - 1) = 29, sd = 10, icc = 0.85:
  #   sem = 10 * sqrt(0.15)            = 3.87298334621
  #   se  = sem / sqrt(2 * 29)         = 0.508547627716
  #   ci  = sem * sqrt(29 / qchisq(c(.975, .025), 29)) = [3.08447078739, 5.20651030125]
  res <- compute_sem(sd = 10, icc = 0.85, n_sample = 30, n_measurements = 2)
  expect_false(is.na(res$sem_se))
  expect_false(is.na(res$sem_ci_lo))
  expect_false(is.na(res$sem_ci_up))
  expect_equal(res$sem, 3.87298334621, tolerance = 1e-8)
  expect_equal(res$sem_se, 0.508547627716, tolerance = 1e-8)
  expect_equal(res$sem_ci_lo, 3.08447078739, tolerance = 1e-8)
  expect_equal(res$sem_ci_up, 5.20651030125, tolerance = 1e-8)

  # (a) the caveat. Sources are only present when the tests run beside them.
  root <- .audit_pkg_root()
  skip_if(is.na(root), "package sources not available (R CMD check copies tests/ only)")

  rd_path <- file.path(root, "man", "compute_sem.Rd")
  r_path  <- file.path(root, "R", "psychometric_utils.R")
  skip_if_not(file.exists(rd_path) && file.exists(r_path),
              "man/compute_sem.Rd or R/psychometric_utils.R not found")

  rd_txt <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")
  r_txt  <- paste(readLines(r_path, warn = FALSE), collapse = "\n")

  # The @details of compute_sem() only; do not let the ICC-type words of some
  # other function in the same file satisfy the check.
  sem_block <- sub(".*#' Compute Standard Error of Measurement", "", r_txt)
  sem_block <- sub("compute_sem <- function.*", "", sem_block)

  agree_re  <- "(?i)(absolute[- ]agreement|ICC\\(2, ?1\\))"
  consis_re <- "(?i)(consistency|ICC\\(3, ?1\\))"

  # CRAN-visible documentation
  expect_true(grepl(agree_re, rd_txt, perl = TRUE),
              info = "man/compute_sem.Rd must name the absolute-agreement ICC(2,1) case")
  expect_true(grepl(consis_re, rd_txt, perl = TRUE),
              info = "man/compute_sem.Rd must say the closed form is exact only for a consistency ICC(3,1)")

  # roxygen source the Rd is generated from
  expect_true(grepl(agree_re, sem_block, perl = TRUE),
              info = "compute_sem()'s @details must name the absolute-agreement ICC(2,1) case")
  expect_true(grepl(consis_re, sem_block, perl = TRUE),
              info = "compute_sem()'s @details must say the closed form is exact only for a consistency ICC(3,1)")
})


# ---------------------------------------------------------------------------
# AUDIT line 544-553 — compute_sdc(): SDC = c * SEM with c = 1.96 * sqrt(2) a
# known constant, so the SDC interval is c times the SEM interval. compute_sdc()
# instead discards the exact chi-square interval compute_sem() just returned and
# rebuilds sdc +/- 1.96 * sdc_se — reinstating exactly the Wald undercoverage
# (0.918 vs 0.95 at n = 10, k = 2) that ?compute_sem cites as its reason to
# abandon the Wald form one step upstream.
# ---------------------------------------------------------------------------
test_that("AUDIT-psychometrics: compute_sdc rescales the SEM interval instead of rebuilding a Wald CI", {

  s <- compute_sem(sd = 10, icc = 0.85, n_sample = 10, n_measurements = 2)

  # The upstream interval is the exact chi-square one, df = (10 - 1)(2 - 1) = 9.
  expect_equal(s$sem, 3.87298334621, tolerance = 1e-8)
  expect_equal(s$sem_ci_lo, 2.66397430067, tolerance = 1e-8)
  expect_equal(s$sem_ci_up, 7.07055783356, tolerance = 1e-8)

  multiplier <- stats::qnorm(0.975) * sqrt(2)   # 2.7718076487

  # Accept either shape of the agreed fix: explicit sem_ci_lo / sem_ci_up
  # arguments, or compute_sdc() taking the whole compute_sem() result.
  d <- tryCatch(
    compute_sdc(sem = s$sem, sem_se = s$sem_se,
                sem_ci_lo = s$sem_ci_lo, sem_ci_up = s$sem_ci_up),
    error = function(e) tryCatch(compute_sdc(s), error = function(e2) NULL))
  ok <- is.data.frame(d) && nrow(d) == length(s$sem) &&
    all(c("sdc", "sdc_se", "sdc_ci_lo", "sdc_ci_up") %in% names(d))
  if (!ok) {
    # Neither shape is available today: compute_sdc() has no way to receive the
    # SEM interval at all, so record the absence as NA and let the value
    # expectations below name the defect.
    d <- data.frame(sdc = NA_real_, sdc_se = NA_real_,
                    sdc_ci_lo = NA_real_, sdc_ci_up = NA_real_)
  }

  expect_equal(nrow(d), 1L)

  # Point estimate and SE are already right and must not move:
  #   sdc    = 2.7718076487 * 3.87298334621  = 10.7351648623
  #   sdc_se = 2.7718076487 * 0.912870929175 =  2.53030262376
  expect_equal(d$sdc, 10.7351648623, tolerance = 1e-8)
  expect_equal(d$sdc_se, 2.53030262376, tolerance = 1e-8)

  # The interval must be the rescaled SEM interval, not the Wald rebuild:
  #   lo = 2.7718076487 * 2.66397430067 =  7.38402434255  (package gives  5.775863)
  #   up = 2.7718076487 * 7.07055783356 = 19.5982262836   (package gives 15.694467)
  expect_equal(d$sdc_ci_lo, 7.38402434255, tolerance = 1e-8)
  expect_equal(d$sdc_ci_up, 19.5982262836, tolerance = 1e-8)

  # Stated as the scale-free identity, so it holds for any input, and it is
  # the property the Wald rebuild breaks: the SDC interval is a monotone
  # rescaling of the SEM interval.
  expect_equal(d$sdc_ci_lo, multiplier * s$sem_ci_lo, tolerance = 1e-10)
  expect_equal(d$sdc_ci_up, multiplier * s$sem_ci_up, tolerance = 1e-10)

  # The Wald fallback must survive for callers with no interval to rescale
  # (audit conflict table, line 195).
  d_wald <- compute_sdc(sem = s$sem, sem_se = s$sem_se)
  expect_equal(d_wald$sdc, 10.7351648623, tolerance = 1e-8)
  expect_equal(d_wald$sdc_se, 2.53030262376, tolerance = 1e-8)
})
