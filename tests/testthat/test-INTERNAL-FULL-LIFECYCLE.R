# Long-running file (> 25 s): executed locally and in CI (devtools::test and
# testthat set NOT_CRAN=true); skipped wholesale on CRAN to respect check-time
# limits. The fast pre/post files still run on CRAN.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
# ==============================================================================
# INTERNAL FULL LIFECYCLE TEST
# ==============================================================================
# Starting from a single real dataset (metaumbrella::df.SMC), derives ALL
# possible input formats (means+SD, means+SE, means+CI, mean change+SD/SE/CI/pval,
# paired t, paired t pval, paired F, paired F pval) and verifies that all
# convert_df() pathways produce identical effect sizes.
#
# This is the definitive internal consistency test for pre/post effect sizes.
#
# Dataset: metaumbrella::df.SMC (33 rows, 100% complete)
# ==============================================================================

library(testthat)
library(metaConvert)

has_metaumbrella <- requireNamespace("metaumbrella", quietly = TRUE)
if (!has_metaumbrella) {
  skip("metaumbrella package not available")
}

# ==============================================================================
# SECTION 0: DATA PREPARATION
# ==============================================================================
# Build a single wide dataframe with ALL derived columns from df.SMC

data("df.SMC", package = "metaumbrella")
dat <- df.SMC

# --- Map column names to metaConvert conventions ---
dat$n_exp <- dat$n_cases
dat$n_nexp <- dat$n_controls
dat$mean_pre_exp <- dat$mean_pre_cases
dat$mean_exp <- dat$mean_cases
dat$mean_pre_sd_exp <- dat$sd_pre_cases
dat$mean_sd_exp <- dat$sd_cases
dat$mean_pre_nexp <- dat$mean_pre_controls
dat$mean_nexp <- dat$mean_controls
dat$mean_pre_sd_nexp <- dat$sd_pre_controls
dat$mean_sd_nexp <- dat$sd_controls
dat$r_pre_post_exp <- dat$pre_post_cor
dat$r_pre_post_nexp <- dat$pre_post_cor

# --- Derive SE of pre/post means ---
dat$mean_pre_se_exp <- dat$mean_pre_sd_exp / sqrt(dat$n_exp)
dat$mean_se_exp <- dat$mean_sd_exp / sqrt(dat$n_exp)
dat$mean_pre_se_nexp <- dat$mean_pre_sd_nexp / sqrt(dat$n_nexp)
dat$mean_se_nexp <- dat$mean_sd_nexp / sqrt(dat$n_nexp)

# --- Derive CI of pre/post means ---
dat$mean_pre_ci_lo_exp <- dat$mean_pre_exp - qt(0.975, dat$n_exp - 1) * dat$mean_pre_se_exp
dat$mean_pre_ci_up_exp <- dat$mean_pre_exp + qt(0.975, dat$n_exp - 1) * dat$mean_pre_se_exp
dat$mean_ci_lo_exp <- dat$mean_exp - qt(0.975, dat$n_exp - 1) * dat$mean_se_exp
dat$mean_ci_up_exp <- dat$mean_exp + qt(0.975, dat$n_exp - 1) * dat$mean_se_exp
dat$mean_pre_ci_lo_nexp <- dat$mean_pre_nexp - qt(0.975, dat$n_nexp - 1) * dat$mean_pre_se_nexp
dat$mean_pre_ci_up_nexp <- dat$mean_pre_nexp + qt(0.975, dat$n_nexp - 1) * dat$mean_pre_se_nexp
dat$mean_ci_lo_nexp <- dat$mean_nexp - qt(0.975, dat$n_nexp - 1) * dat$mean_se_nexp
dat$mean_ci_up_nexp <- dat$mean_nexp + qt(0.975, dat$n_nexp - 1) * dat$mean_se_nexp

# --- Derive mean change statistics ---
dat$mean_change_exp <- dat$mean_exp - dat$mean_pre_exp
dat$mean_change_nexp <- dat$mean_nexp - dat$mean_pre_nexp

dat$mean_change_sd_exp <- sqrt(
  dat$mean_pre_sd_exp^2 + dat$mean_sd_exp^2 -
    2 * dat$r_pre_post_exp * dat$mean_pre_sd_exp * dat$mean_sd_exp
)
dat$mean_change_sd_nexp <- sqrt(
  dat$mean_pre_sd_nexp^2 + dat$mean_sd_nexp^2 -
    2 * dat$r_pre_post_nexp * dat$mean_pre_sd_nexp * dat$mean_sd_nexp
)

dat$mean_change_se_exp <- dat$mean_change_sd_exp / sqrt(dat$n_exp)
dat$mean_change_se_nexp <- dat$mean_change_sd_nexp / sqrt(dat$n_nexp)

dat$mean_change_ci_lo_exp <- dat$mean_change_exp - qt(0.975, dat$n_exp - 1) * dat$mean_change_se_exp
dat$mean_change_ci_up_exp <- dat$mean_change_exp + qt(0.975, dat$n_exp - 1) * dat$mean_change_se_exp
dat$mean_change_ci_lo_nexp <- dat$mean_change_nexp - qt(0.975, dat$n_nexp - 1) * dat$mean_change_se_nexp
dat$mean_change_ci_up_nexp <- dat$mean_change_nexp + qt(0.975, dat$n_nexp - 1) * dat$mean_change_se_nexp

# --- Derive paired t/F statistics ---
dat$paired_t_exp <- dat$mean_change_exp / dat$mean_change_se_exp
dat$paired_t_nexp <- dat$mean_change_nexp / dat$mean_change_se_nexp

dat$paired_f_exp <- dat$paired_t_exp^2
dat$paired_f_nexp <- dat$paired_t_nexp^2

dat$paired_t_pval_exp <- 2 * pt(-abs(dat$paired_t_exp), dat$n_exp - 1)
dat$paired_t_pval_nexp <- 2 * pt(-abs(dat$paired_t_nexp), dat$n_nexp - 1)

dat$paired_f_pval_exp <- dat$paired_t_pval_exp
dat$paired_f_pval_nexp <- dat$paired_t_pval_nexp

# --- Derive mean change p-values (same as paired t p-values) ---
dat$mean_change_pval_exp <- dat$paired_t_pval_exp
dat$mean_change_pval_nexp <- dat$paired_t_pval_nexp

# --- Remove original df.SMC columns to avoid interference ---
# Keep only metaConvert columns to ensure convert_df picks up our derived columns
cols_to_remove <- c("n_cases", "n_controls", "mean_pre_cases", "mean_cases",
                    "sd_pre_cases", "sd_cases", "mean_pre_controls", "mean_controls",
                    "sd_pre_controls", "sd_controls", "pre_post_cor",
                    "value", "se", "measure")
dat <- dat[, setdiff(names(dat), cols_to_remove)]

# --- Helper: run convert_df for a given hierarchy and extract results ---
#
# pool_sd defaults to TRUE, matching convert_df()'s own (current) default:
#   pool_sd = TRUE  -> between-group SMD = (difference in mean change) / (SD pooled
#                      across arms). Morris (2008)'s recommended common standardizer.
#   pool_sd = FALSE -> LEGACY per-arm construction: each arm is standardized by its
#                      OWN SD and the two within-group values are then subtracted.
#                      Valid only when the arms' SDs are equal; kept for backward
#                      compatibility and still exercised below.
run_pathway <- function(data, hierarchy, measure = "d",
                        pre_post_to_smd = "cooper", pool_sd = TRUE) {
  res <- convert_df(data,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = hierarchy,
    measure = measure,
    pre_post_to_smd = pre_post_to_smd,
    pool_sd = pool_sd
  )
  summ <- summary(res, digits = 11)
  list(
    es = summ$es_crude,
    se = summ$se_crude,
    info = summ$info_used_crude
  )
}

# IMPORTANT — why some comparisons below must pass pool_sd = FALSE.
#
# The paired-t / paired-F routes CANNOT pool. A paired t identifies each arm's
# mean_change and sd_change, but NOT the two arms' SD ratio, so the pooled
# standardizer is mathematically unrecoverable from it. Those routes therefore
# keep the per-arm construction and expose no pool_sd argument.
#
# Consequently a "paired_t == means/mean_change" equivalence assertion is only
# true when the means-side call is run on the same per-arm footing, i.e. with
# pool_sd = FALSE (or on a fixture whose arm SDs happen to be equal — df.SMC's
# are not). Every such comparison below passes pool_sd = FALSE on the means side
# and says so. All other comparisons run under the pooled default.


# ==============================================================================
# SECTION 1: TWO-GROUP — ALL PATHWAYS MUST MATCH (cooper/morris_drm)
# ==============================================================================
# These run under the CURRENT DEFAULT standardizer (pool_sd = TRUE): every
# means/mean-change route must agree once they share the pooled standardizer.

test_that("LIFECYCLE: means_sd vs means_se — cooper, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "cooper")
  alt <- run_pathway(dat, "means_se_pre_post", "d", "cooper")
  expect_equal(ref$info, rep("means_sd_pre_post", nrow(dat)))
  expect_equal(alt$info, rep("means_se_pre_post", nrow(dat)))
  expect_equal(ref$es, alt$es, tolerance = 1e-10)
  expect_equal(ref$se, alt$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: means_sd vs means_ci — cooper, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "cooper")
  alt <- run_pathway(dat, "means_ci_pre_post", "d", "cooper")
  expect_equal(alt$info, rep("means_ci_pre_post", nrow(dat)))
  expect_equal(ref$es, alt$es, tolerance = 1e-10)
  expect_equal(ref$se, alt$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: means_sd vs mean_change_sd — cooper, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "cooper")
  alt <- run_pathway(dat, "mean_change_sd", "d", "cooper")
  expect_equal(alt$info, rep("mean_change_sd", nrow(dat)))
  expect_equal(ref$es, alt$es, tolerance = 1e-10)
  expect_equal(ref$se, alt$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: means_sd vs mean_change_se — cooper, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "cooper")
  alt <- run_pathway(dat, "mean_change_se", "d", "cooper")
  expect_equal(alt$info, rep("mean_change_se", nrow(dat)))
  expect_equal(ref$es, alt$es, tolerance = 1e-10)
  expect_equal(ref$se, alt$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: means_sd vs mean_change_ci — cooper, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "cooper")
  alt <- run_pathway(dat, "mean_change_ci", "d", "cooper")
  expect_equal(alt$info, rep("mean_change_ci", nrow(dat)))
  expect_equal(ref$es, alt$es, tolerance = 1e-10)
  expect_equal(ref$se, alt$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: means_sd vs mean_change_pval — cooper, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "cooper")
  alt <- run_pathway(dat, "mean_change_pval", "d", "cooper")
  expect_equal(alt$info, rep("mean_change_pval", nrow(dat)))
  expect_equal(ref$es, alt$es, tolerance = 1e-6)
  expect_equal(ref$se, alt$se, tolerance = 1e-6)
})

test_that("LIFECYCLE: means_sd vs paired_t — cooper, d (legacy per-arm standardizer, pool_sd = FALSE)", {
  # paired_t cannot pool (it does not identify the arms' SD ratio), so it always
  # builds the per-arm standardizer. The means side is therefore run with
  # pool_sd = FALSE to put both routes on the same footing; this pins the LEGACY
  # per-arm path, which remains a supported (backward-compatible) construction.
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "cooper", pool_sd = FALSE)
  alt <- run_pathway(dat, "paired_t", "d", "cooper", pool_sd = FALSE)
  expect_equal(alt$info, rep("paired_t", nrow(dat)))
  expect_equal(ref$es, alt$es, tolerance = 1e-10)
  expect_equal(ref$se, alt$se, tolerance = 1e-10)
})

# NOTE: paired_t_pval, paired_f, paired_f_pval lose per-group sign information.
# For two-group designs, this means d = |d_exp| - |d_nexp| instead of d_exp - d_nexp,
# producing different magnitudes when individual group effects have different signs.
# So we compare these UNSIGNED pathways against EACH OTHER (not against means_sd).

test_that("LIFECYCLE: unsigned pathways match each other — cooper, d", {
  alt_tpval <- run_pathway(dat, "paired_t_pval", "d", "cooper")
  alt_f     <- run_pathway(dat, "paired_f", "d", "cooper")
  alt_fpval <- run_pathway(dat, "paired_f_pval", "d", "cooper")

  expect_equal(alt_tpval$info, rep("paired_t_pval", nrow(dat)))
  expect_equal(alt_f$info, rep("paired_f", nrow(dat)))
  expect_equal(alt_fpval$info, rep("paired_f_pval", nrow(dat)))

  # F = t^2, so paired_f and paired_t_pval should give identical abs(d)
  expect_equal(abs(alt_f$es), abs(alt_tpval$es), tolerance = 1e-10)
  expect_equal(alt_f$se, alt_tpval$se, tolerance = 1e-10)

  # F_pval and t_pval give same unsigned per-group d (pval roundtrip)
  expect_equal(abs(alt_fpval$es), abs(alt_tpval$es), tolerance = 1e-6)
  expect_equal(alt_fpval$se, alt_tpval$se, tolerance = 1e-6)

  # All produce valid (non-NA, finite) results
  expect_true(all(!is.na(alt_tpval$es) & is.finite(alt_tpval$es)))
  expect_true(all(!is.na(alt_f$es) & is.finite(alt_f$es)))
  expect_true(all(alt_tpval$se > 0))
})

# --- Same for g ---

test_that("LIFECYCLE: means_sd vs all pathways — cooper, g", {
  ref <- run_pathway(dat, "means_sd_pre_post", "g", "cooper")

  # SE-derived (exact)
  alt_se <- run_pathway(dat, "means_se_pre_post", "g", "cooper")
  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)

  # CI-derived (exact)
  alt_ci <- run_pathway(dat, "means_ci_pre_post", "g", "cooper")
  expect_equal(ref$es, alt_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_ci$se, tolerance = 1e-10)

  # Mean change SD (exact)
  alt_mc_sd <- run_pathway(dat, "mean_change_sd", "g", "cooper")
  expect_equal(ref$es, alt_mc_sd$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_sd$se, tolerance = 1e-10)

  # Mean change SE (exact)
  alt_mc_se <- run_pathway(dat, "mean_change_se", "g", "cooper")
  expect_equal(ref$es, alt_mc_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_se$se, tolerance = 1e-10)

  # Mean change CI (exact)
  alt_mc_ci <- run_pathway(dat, "mean_change_ci", "g", "cooper")
  expect_equal(ref$es, alt_mc_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_ci$se, tolerance = 1e-10)

  # Mean change pval (near-exact, sign from mean_change column)
  alt_mc_pval <- run_pathway(dat, "mean_change_pval", "g", "cooper")
  expect_equal(ref$es, alt_mc_pval$es, tolerance = 1e-6)
  expect_equal(ref$se, alt_mc_pval$se, tolerance = 1e-6)

  # Paired t (exact) — compared against the LEGACY per-arm means reference,
  # because paired_t cannot pool (see the note above run_pathway).
  ref_perarm <- run_pathway(dat, "means_sd_pre_post", "g", "cooper", pool_sd = FALSE)
  alt_t <- run_pathway(dat, "paired_t", "g", "cooper", pool_sd = FALSE)
  expect_equal(ref_perarm$es, alt_t$es, tolerance = 1e-10)
  expect_equal(ref_perarm$se, alt_t$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: unsigned pathways match each other — cooper, g", {
  # Unsigned pathways lose per-group sign; compare against each other only
  alt_tpval <- run_pathway(dat, "paired_t_pval", "g", "cooper")
  alt_f     <- run_pathway(dat, "paired_f", "g", "cooper")
  alt_fpval <- run_pathway(dat, "paired_f_pval", "g", "cooper")

  expect_equal(abs(alt_f$es), abs(alt_tpval$es), tolerance = 1e-10)
  expect_equal(alt_f$se, alt_tpval$se, tolerance = 1e-10)
  expect_equal(abs(alt_fpval$es), abs(alt_tpval$es), tolerance = 1e-6)
  expect_equal(alt_fpval$se, alt_tpval$se, tolerance = 1e-6)
})


# ==============================================================================
# SECTION 2: TWO-GROUP — ALL PATHWAYS WITH morris_dz
# ==============================================================================

test_that("LIFECYCLE: signed pathways match — morris_dz, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "morris_dz")

  # SE-derived (exact)
  alt_se <- run_pathway(dat, "means_se_pre_post", "d", "morris_dz")
  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)

  # CI-derived (exact)
  alt_ci <- run_pathway(dat, "means_ci_pre_post", "d", "morris_dz")
  expect_equal(ref$es, alt_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_ci$se, tolerance = 1e-10)

  # Mean change SD (exact)
  alt_mc_sd <- run_pathway(dat, "mean_change_sd", "d", "morris_dz")
  expect_equal(ref$es, alt_mc_sd$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_sd$se, tolerance = 1e-10)

  # Mean change SE (exact)
  alt_mc_se <- run_pathway(dat, "mean_change_se", "d", "morris_dz")
  expect_equal(ref$es, alt_mc_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_se$se, tolerance = 1e-10)

  # Mean change CI (exact)
  alt_mc_ci <- run_pathway(dat, "mean_change_ci", "d", "morris_dz")
  expect_equal(ref$es, alt_mc_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_ci$se, tolerance = 1e-10)

  # Mean change pval (near-exact, sign from mean_change column)
  alt_mc_pval <- run_pathway(dat, "mean_change_pval", "d", "morris_dz")
  expect_equal(ref$es, alt_mc_pval$es, tolerance = 1e-6)
  expect_equal(ref$se, alt_mc_pval$se, tolerance = 1e-6)

  # Paired t — compared against the LEGACY per-arm means reference, because
  # paired_t cannot pool (see the note above run_pathway).
  # BOTH d and SE now match EXACTLY. The SE tolerance used to be 0.02 to absorb a
  # variance-convention mismatch: the paired_t morris_dz variance was built as
  # J^2/n + g^2/(2n), while the means route used the metafor SMCC convention
  # (variance built from the CORRECTED g: 1/n + g^2/(2n)). The paired_t route now
  # uses the SMCC convention too, so the two agree to machine precision and the
  # tolerance is tightened from 0.02 to 1e-10.
  ref_perarm <- run_pathway(dat, "means_sd_pre_post", "d", "morris_dz", pool_sd = FALSE)
  alt_t <- run_pathway(dat, "paired_t", "d", "morris_dz", pool_sd = FALSE)
  expect_equal(ref_perarm$es, alt_t$es, tolerance = 1e-10)
  expect_equal(ref_perarm$se, alt_t$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: unsigned pathways match each other — morris_dz, d", {
  alt_tpval <- run_pathway(dat, "paired_t_pval", "d", "morris_dz")
  alt_f     <- run_pathway(dat, "paired_f", "d", "morris_dz")
  alt_fpval <- run_pathway(dat, "paired_f_pval", "d", "morris_dz")

  expect_equal(abs(alt_f$es), abs(alt_tpval$es), tolerance = 1e-10)
  expect_equal(alt_f$se, alt_tpval$se, tolerance = 1e-10)
  expect_equal(abs(alt_fpval$es), abs(alt_tpval$es), tolerance = 1e-6)
  expect_equal(alt_fpval$se, alt_tpval$se, tolerance = 1e-6)
})

test_that("LIFECYCLE: signed pathways match — morris_dz, g", {
  ref <- run_pathway(dat, "means_sd_pre_post", "g", "morris_dz")

  alt_se <- run_pathway(dat, "means_se_pre_post", "g", "morris_dz")
  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)

  alt_ci <- run_pathway(dat, "means_ci_pre_post", "g", "morris_dz")
  expect_equal(ref$es, alt_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_ci$se, tolerance = 1e-10)

  alt_mc_sd <- run_pathway(dat, "mean_change_sd", "g", "morris_dz")
  expect_equal(ref$es, alt_mc_sd$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_sd$se, tolerance = 1e-10)

  alt_mc_se <- run_pathway(dat, "mean_change_se", "g", "morris_dz")
  expect_equal(ref$es, alt_mc_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_se$se, tolerance = 1e-10)

  alt_mc_ci <- run_pathway(dat, "mean_change_ci", "g", "morris_dz")
  expect_equal(ref$es, alt_mc_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_ci$se, tolerance = 1e-10)

  alt_mc_pval <- run_pathway(dat, "mean_change_pval", "g", "morris_dz")
  expect_equal(ref$es, alt_mc_pval$es, tolerance = 1e-6)
  expect_equal(ref$se, alt_mc_pval$se, tolerance = 1e-6)

  # Paired t — LEGACY per-arm reference (paired_t cannot pool). Both g and SE now
  # match exactly; SE tolerance tightened 0.02 -> 1e-10 (see morris_dz d note).
  ref_perarm <- run_pathway(dat, "means_sd_pre_post", "g", "morris_dz", pool_sd = FALSE)
  alt_t <- run_pathway(dat, "paired_t", "g", "morris_dz", pool_sd = FALSE)
  expect_equal(ref_perarm$es, alt_t$es, tolerance = 1e-10)
  expect_equal(ref_perarm$se, alt_t$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: unsigned pathways match each other — morris_dz, g", {
  alt_tpval <- run_pathway(dat, "paired_t_pval", "g", "morris_dz")
  alt_f     <- run_pathway(dat, "paired_f", "g", "morris_dz")
  alt_fpval <- run_pathway(dat, "paired_f_pval", "g", "morris_dz")

  expect_equal(abs(alt_f$es), abs(alt_tpval$es), tolerance = 1e-10)
  expect_equal(alt_f$se, alt_tpval$se, tolerance = 1e-10)
  expect_equal(abs(alt_fpval$es), abs(alt_tpval$es), tolerance = 1e-6)
  expect_equal(alt_fpval$se, alt_tpval$se, tolerance = 1e-6)
})


# ==============================================================================
# SECTION 3: TWO-GROUP — BONETT AND MORRIS_DAV (pre-post means only)
# ==============================================================================
# These methods only work with pre-post means functions (not mean change or
# paired t/f, which require separate pre/post SDs unavailable from change scores)

test_that("LIFECYCLE: SD/SE/CI give identical results — bonett, g", {
  ref <- run_pathway(dat, "means_sd_pre_post", "g", "bonett")
  alt_se <- run_pathway(dat, "means_se_pre_post", "g", "bonett")
  alt_ci <- run_pathway(dat, "means_ci_pre_post", "g", "bonett")

  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)
  expect_equal(ref$es, alt_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_ci$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: SD/SE/CI give identical results — morris_dav, g", {
  ref <- run_pathway(dat, "means_sd_pre_post", "g", "morris_dav")
  alt_se <- run_pathway(dat, "means_se_pre_post", "g", "morris_dav")
  alt_ci <- run_pathway(dat, "means_ci_pre_post", "g", "morris_dav")

  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)
  expect_equal(ref$es, alt_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_ci$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: SD/SE/CI give identical results — bonett, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "bonett")
  alt_se <- run_pathway(dat, "means_se_pre_post", "d", "bonett")
  alt_ci <- run_pathway(dat, "means_ci_pre_post", "d", "bonett")

  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)
  expect_equal(ref$es, alt_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_ci$se, tolerance = 1e-10)
})

test_that("LIFECYCLE: SD/SE/CI give identical results — morris_dav, d", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "morris_dav")
  alt_se <- run_pathway(dat, "means_se_pre_post", "d", "morris_dav")
  alt_ci <- run_pathway(dat, "means_ci_pre_post", "d", "morris_dav")

  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)
  expect_equal(ref$es, alt_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_ci$se, tolerance = 1e-10)
})


# ==============================================================================
# SECTION 4: SINGLE-GROUP — ALL PATHWAYS MUST MATCH (cooper)
# ==============================================================================
# Single-group functions are called directly (not through convert_df)

test_that("LIFECYCLE: single-group — all pathways match, cooper", {
  for (i in 1:nrow(dat)) {
    row <- dat[i, ]

    # Reference: means + SD
    ref <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = row$mean_pre_exp, mean_exp = row$mean_exp,
      mean_pre_sd_exp = row$mean_pre_sd_exp, mean_sd_exp = row$mean_sd_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp,
      pre_post_to_smd = "cooper"
    )

    # Means + SE
    alt_se <- es_from_means_se_pre_post_single_group(
      mean_pre_exp = row$mean_pre_exp, mean_exp = row$mean_exp,
      mean_pre_se_exp = row$mean_pre_se_exp, mean_se_exp = row$mean_se_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp,
      pre_post_to_smd = "cooper"
    )
    expect_equal(ref$d, alt_se$d, tolerance = 1e-10,
                 info = paste("SE d mismatch, row", i))
    expect_equal(ref$d_se, alt_se$d_se, tolerance = 1e-10,
                 info = paste("SE d_se mismatch, row", i))

    # Means + CI
    alt_ci <- es_from_means_ci_pre_post_single_group(
      mean_pre_exp = row$mean_pre_exp, mean_exp = row$mean_exp,
      mean_pre_ci_lo_exp = row$mean_pre_ci_lo_exp,
      mean_pre_ci_up_exp = row$mean_pre_ci_up_exp,
      mean_ci_lo_exp = row$mean_ci_lo_exp,
      mean_ci_up_exp = row$mean_ci_up_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp,
      pre_post_to_smd = "cooper"
    )
    expect_equal(ref$d, alt_ci$d, tolerance = 1e-10,
                 info = paste("CI d mismatch, row", i))
    expect_equal(ref$d_se, alt_ci$d_se, tolerance = 1e-10,
                 info = paste("CI d_se mismatch, row", i))

    # Mean change + SD
    alt_mc_sd <- es_from_mean_change_sd_single_group(
      mean_change_exp = row$mean_change_exp,
      mean_change_sd_exp = row$mean_change_sd_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp
    )
    expect_equal(ref$d, alt_mc_sd$d, tolerance = 1e-10,
                 info = paste("MC-SD d mismatch, row", i))
    expect_equal(ref$d_se, alt_mc_sd$d_se, tolerance = 1e-10,
                 info = paste("MC-SD d_se mismatch, row", i))

    # Mean change + SE
    alt_mc_se <- es_from_mean_change_se_single_group(
      mean_change_exp = row$mean_change_exp,
      mean_change_se_exp = row$mean_change_se_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp
    )
    expect_equal(ref$d, alt_mc_se$d, tolerance = 1e-10,
                 info = paste("MC-SE d mismatch, row", i))
    expect_equal(ref$d_se, alt_mc_se$d_se, tolerance = 1e-10,
                 info = paste("MC-SE d_se mismatch, row", i))

    # Mean change + CI
    alt_mc_ci <- es_from_mean_change_ci_single_group(
      mean_change_exp = row$mean_change_exp,
      mean_change_ci_lo_exp = row$mean_change_ci_lo_exp,
      mean_change_ci_up_exp = row$mean_change_ci_up_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp
    )
    expect_equal(ref$d, alt_mc_ci$d, tolerance = 1e-10,
                 info = paste("MC-CI d mismatch, row", i))
    expect_equal(ref$d_se, alt_mc_ci$d_se, tolerance = 1e-10,
                 info = paste("MC-CI d_se mismatch, row", i))

    # Mean change + pval
    alt_mc_pval <- es_from_mean_change_pval_single_group(
      mean_change_exp = row$mean_change_exp,
      mean_change_pval_exp = row$paired_t_pval_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp
    )
    expect_equal(ref$d, alt_mc_pval$d, tolerance = 1e-6,
                 info = paste("MC-pval d mismatch, row", i))
    expect_equal(ref$d_se, alt_mc_pval$d_se, tolerance = 1e-6,
                 info = paste("MC-pval d_se mismatch, row", i))

    # Paired t
    alt_t <- es_from_paired_t_single_group(
      paired_t_exp = row$paired_t_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp
    )
    expect_equal(ref$d, alt_t$d, tolerance = 1e-10,
                 info = paste("paired_t d mismatch, row", i))
    expect_equal(ref$d_se, alt_t$d_se, tolerance = 1e-10,
                 info = paste("paired_t d_se mismatch, row", i))
  }
})


# ==============================================================================
# SECTION 5: CROSS-VALIDATION — TWO-GROUP d ~ SINGLE-GROUP exp - nexp
# ==============================================================================
# "two-group d == d_exp - d_nexp" IS the definition of the LEGACY per-arm
# standardizer, so this cross-validation is run with pool_sd = FALSE. It remains
# a valid regression test OF THAT PATH.
#
# Under the pooled default (pool_sd = TRUE) the identity does NOT hold, and must
# not: a common standardizer is a structurally different estimand. On this fixture
# the pooled d departs from (d_exp - d_nexp) by up to 0.46, which is exactly the
# per-arm bias the pooled default was introduced to remove. That contrast is
# pinned by the "pool_sd default" test in Section 6.

test_that("LIFECYCLE: two-group d == single_group_exp - single_group_nexp (bonett; legacy per-arm standardizer, pool_sd = FALSE)", {
  for (i in 1:nrow(dat)) {
    row <- dat[i, ]

    # Two-group, legacy per-arm standardizer
    two_grp <- run_pathway(dat[i, ], "means_sd_pre_post", "d", "bonett", pool_sd = FALSE)

    # Single-group exp
    sg_exp <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = row$mean_pre_exp, mean_exp = row$mean_exp,
      mean_pre_sd_exp = row$mean_pre_sd_exp, mean_sd_exp = row$mean_sd_exp,
      n_exp = row$n_exp, r_pre_post_exp = row$r_pre_post_exp,
      pre_post_to_smd = "bonett"
    )

    # Single-group nexp
    sg_nexp <- es_from_means_sd_pre_post_single_group(
      mean_pre_exp = row$mean_pre_nexp, mean_exp = row$mean_nexp,
      mean_pre_sd_exp = row$mean_pre_sd_nexp, mean_sd_exp = row$mean_sd_nexp,
      n_exp = row$n_nexp, r_pre_post_exp = row$r_pre_post_nexp,
      pre_post_to_smd = "bonett"
    )

    # Under the per-arm standardizer this is an exact identity, not an
    # approximation: max |diff| across the 33 rows is ~5e-12. Tolerance tightened
    # from the previous 0.1 (which was far looser than the identity warrants).
    expected_diff <- sg_exp$d - sg_nexp$d
    expect_equal(two_grp$es, expected_diff, tolerance = 1e-8,
                 info = paste("Two-group vs diff of single-group, row", i))
  }
})


# ==============================================================================
# SECTION 6: SANITY CHECKS ACROSS ALL 33 ROWS
# ==============================================================================

test_that("LIFECYCLE: all rows produce finite non-NA values", {
  for (method in c("cooper", "bonett", "morris_dav", "morris_dz")) {
    res <- run_pathway(dat, "means_sd_pre_post", "d", method)
    expect_true(all(!is.na(res$es)),
                info = paste(method, "has NA es values"))
    expect_true(all(is.finite(res$es)),
                info = paste(method, "has non-finite es values"))
    expect_true(all(!is.na(res$se)),
                info = paste(method, "has NA se values"))
    expect_true(all(is.finite(res$se)),
                info = paste(method, "has non-finite se values"))
  }
})

test_that("LIFECYCLE: g and d have the same sign for all methods", {
  for (method in c("cooper", "bonett", "morris_dav", "morris_dz")) {
    d_res <- run_pathway(dat, "means_sd_pre_post", "d", method)
    g_res <- run_pathway(dat, "means_sd_pre_post", "g", method)
    # Where d != 0, sign should match
    nonzero <- d_res$es != 0
    expect_true(all(sign(d_res$es[nonzero]) == sign(g_res$es[nonzero])),
                info = paste(method, "sign mismatch between d and g"))
  }
})

# NOTE: |g| <= |d| does NOT hold for two-group designs because
# g = J_exp * d_exp - J_nexp * d_nexp ≠ J * (d_exp - d_nexp).
# The bias correction is applied per-group BEFORE subtraction,
# so the two-group g can exceed d. This is tested in single-group
# Section 4 where the property does hold.

test_that("LIFECYCLE: SE is positive for all methods", {
  for (method in c("cooper", "bonett", "morris_dav", "morris_dz")) {
    d_res <- run_pathway(dat, "means_sd_pre_post", "d", method)
    expect_true(all(d_res$se > 0),
                info = paste(method, "has non-positive SE"))
  }
})

test_that("LIFECYCLE: different methods produce different values (not all identical)", {
  cooper_d <- run_pathway(dat, "means_sd_pre_post", "d", "cooper")
  bonett_d <- run_pathway(dat, "means_sd_pre_post", "d", "bonett")
  dz_d <- run_pathway(dat, "means_sd_pre_post", "d", "morris_dz")
  dav_d <- run_pathway(dat, "means_sd_pre_post", "d", "morris_dav")

  # Different standardization methods should NOT produce identical results
  # (unless there's a degenerate case, which doesn't apply to df.SMC)
  expect_false(all(abs(cooper_d$es - bonett_d$es) < 1e-10),
               info = "cooper and bonett should differ")
  expect_false(all(abs(cooper_d$es - dz_d$es) < 1e-10),
               info = "cooper and morris_dz should differ")
  expect_false(all(abs(bonett_d$es - dav_d$es) < 1e-10),
               info = "bonett and morris_dav should differ")
})

test_that("LIFECYCLE: convert_df defaults to the POOLED standardizer (pool_sd = TRUE)", {
  # Guards the default. If pool_sd ever silently reverts to the legacy per-arm
  # standardizer, the two-group SMDs change materially and this fails.
  for (method in c("cooper", "bonett", "morris_dav", "morris_dz")) {
    implicit <- run_pathway(dat, "means_sd_pre_post", "d", method)
    pooled   <- run_pathway(dat, "means_sd_pre_post", "d", method, pool_sd = TRUE)
    per_arm  <- run_pathway(dat, "means_sd_pre_post", "d", method, pool_sd = FALSE)

    # The default IS the pooled standardizer.
    expect_equal(implicit$es, pooled$es, tolerance = 1e-12,
                 info = paste(method, "default should equal pool_sd = TRUE"))
    expect_equal(implicit$se, pooled$se, tolerance = 1e-12,
                 info = paste(method, "default SE should equal pool_sd = TRUE"))

    # ...and it is genuinely a different estimand from the legacy per-arm path.
    # df.SMC's arms have unequal SDs, so the two must not coincide. (They would
    # coincide only on a fixture with equal arm SDs.)
    expect_false(all(abs(pooled$es - per_arm$es) < 1e-8),
                 info = paste(method, "pooled and per-arm standardizers should differ"))
  }
})


# ==============================================================================
# SECTION 7: RAW DATA → t.test → convert_df
# ==============================================================================
# Simulates the real-world workflow: researcher has raw paired data, runs
# t.test(paired=TRUE), extracts summary stats, and enters them into convert_df.
# Verifies all pathways produce identical results from a single ground truth.

# --- Generate raw paired data ---
set.seed(12345)
n_raw_exp <- 25
n_raw_nexp <- 30

# Experimental group: correlated pre-post with positive effect
pre_exp_raw <- rnorm(n_raw_exp, mean = 50, sd = 10)
post_exp_raw <- pre_exp_raw + rnorm(n_raw_exp, mean = 5, sd = 6)

# Control group: correlated pre-post with smaller effect
pre_nexp_raw <- rnorm(n_raw_nexp, mean = 50, sd = 10)
post_nexp_raw <- pre_nexp_raw + rnorm(n_raw_nexp, mean = 1, sd = 6)

# --- Run t.test for each group ---
tt_exp <- t.test(post_exp_raw, pre_exp_raw, paired = TRUE)
tt_nexp <- t.test(post_nexp_raw, pre_nexp_raw, paired = TRUE)

# --- Extract t.test outputs ---
# Experimental group
mc_exp <- as.numeric(tt_exp$estimate)
mc_se_exp <- tt_exp$stderr
mc_ci_lo_exp <- tt_exp$conf.int[1]
mc_ci_up_exp <- tt_exp$conf.int[2]
mc_pval_exp <- tt_exp$p.value
pt_exp <- as.numeric(tt_exp$statistic)

# Control group
mc_nexp <- as.numeric(tt_nexp$estimate)
mc_se_nexp <- tt_nexp$stderr
mc_ci_lo_nexp <- tt_nexp$conf.int[1]
mc_ci_up_nexp <- tt_nexp$conf.int[2]
mc_pval_nexp <- tt_nexp$p.value
pt_nexp <- as.numeric(tt_nexp$statistic)

# --- Compute additional stats from raw data ---
mc_sd_exp <- sd(post_exp_raw - pre_exp_raw)
mc_sd_nexp <- sd(post_nexp_raw - pre_nexp_raw)
r_exp <- cor(pre_exp_raw, post_exp_raw)
r_nexp <- cor(pre_nexp_raw, post_nexp_raw)

# Pre/post summary stats from raw data
mean_pre_exp_raw <- mean(pre_exp_raw)
mean_post_exp_raw <- mean(post_exp_raw)
sd_pre_exp_raw <- sd(pre_exp_raw)
sd_post_exp_raw <- sd(post_exp_raw)
se_pre_exp_raw <- sd_pre_exp_raw / sqrt(n_raw_exp)
se_post_exp_raw <- sd_post_exp_raw / sqrt(n_raw_exp)

mean_pre_nexp_raw <- mean(pre_nexp_raw)
mean_post_nexp_raw <- mean(post_nexp_raw)
sd_pre_nexp_raw <- sd(pre_nexp_raw)
sd_post_nexp_raw <- sd(post_nexp_raw)
se_pre_nexp_raw <- sd_pre_nexp_raw / sqrt(n_raw_nexp)
se_post_nexp_raw <- sd_post_nexp_raw / sqrt(n_raw_nexp)

# Pre/post CIs from raw data
ci_lo_pre_exp <- mean_pre_exp_raw - qt(0.975, n_raw_exp - 1) * se_pre_exp_raw
ci_up_pre_exp <- mean_pre_exp_raw + qt(0.975, n_raw_exp - 1) * se_pre_exp_raw
ci_lo_post_exp <- mean_post_exp_raw - qt(0.975, n_raw_exp - 1) * se_post_exp_raw
ci_up_post_exp <- mean_post_exp_raw + qt(0.975, n_raw_exp - 1) * se_post_exp_raw

ci_lo_pre_nexp <- mean_pre_nexp_raw - qt(0.975, n_raw_nexp - 1) * se_pre_nexp_raw
ci_up_pre_nexp <- mean_pre_nexp_raw + qt(0.975, n_raw_nexp - 1) * se_pre_nexp_raw
ci_lo_post_nexp <- mean_post_nexp_raw - qt(0.975, n_raw_nexp - 1) * se_post_nexp_raw
ci_up_post_nexp <- mean_post_nexp_raw + qt(0.975, n_raw_nexp - 1) * se_post_nexp_raw

# Derived: F and F pval
pf_exp <- pt_exp^2
pf_nexp <- pt_nexp^2

# --- Build 1-row dataframe for convert_df ---
dat_raw <- data.frame(
  # Sample sizes and correlations
  n_exp = n_raw_exp, n_nexp = n_raw_nexp,
  r_pre_post_exp = r_exp, r_pre_post_nexp = r_nexp,
  # Pre/post means and SDs
  mean_pre_exp = mean_pre_exp_raw, mean_exp = mean_post_exp_raw,
  mean_pre_sd_exp = sd_pre_exp_raw, mean_sd_exp = sd_post_exp_raw,
  mean_pre_nexp = mean_pre_nexp_raw, mean_nexp = mean_post_nexp_raw,
  mean_pre_sd_nexp = sd_pre_nexp_raw, mean_sd_nexp = sd_post_nexp_raw,
  # Pre/post SEs
  mean_pre_se_exp = se_pre_exp_raw, mean_se_exp = se_post_exp_raw,
  mean_pre_se_nexp = se_pre_nexp_raw, mean_se_nexp = se_post_nexp_raw,
  # Pre/post CIs
  mean_pre_ci_lo_exp = ci_lo_pre_exp, mean_pre_ci_up_exp = ci_up_pre_exp,
  mean_ci_lo_exp = ci_lo_post_exp, mean_ci_up_exp = ci_up_post_exp,
  mean_pre_ci_lo_nexp = ci_lo_pre_nexp, mean_pre_ci_up_nexp = ci_up_pre_nexp,
  mean_ci_lo_nexp = ci_lo_post_nexp, mean_ci_up_nexp = ci_up_post_nexp,
  # Mean change stats (from t.test)
  mean_change_exp = mc_exp, mean_change_nexp = mc_nexp,
  mean_change_sd_exp = mc_sd_exp, mean_change_sd_nexp = mc_sd_nexp,
  mean_change_se_exp = mc_se_exp, mean_change_se_nexp = mc_se_nexp,
  mean_change_ci_lo_exp = mc_ci_lo_exp, mean_change_ci_up_exp = mc_ci_up_exp,
  mean_change_ci_lo_nexp = mc_ci_lo_nexp, mean_change_ci_up_nexp = mc_ci_up_nexp,
  mean_change_pval_exp = mc_pval_exp, mean_change_pval_nexp = mc_pval_nexp,
  # Paired t/F (from t.test)
  paired_t_exp = pt_exp, paired_t_nexp = pt_nexp,
  paired_t_pval_exp = mc_pval_exp, paired_t_pval_nexp = mc_pval_nexp,
  paired_f_exp = pf_exp, paired_f_nexp = pf_nexp,
  paired_f_pval_exp = mc_pval_exp, paired_f_pval_nexp = mc_pval_nexp
)

# --- Test 7a: t.test outputs match manual calculations ---
test_that("LIFECYCLE-RAW: t.test outputs match manual calculations", {
  # estimate = mean(post) - mean(pre)
  expect_equal(mc_exp, mean(post_exp_raw) - mean(pre_exp_raw), tolerance = 1e-10)
  expect_equal(mc_nexp, mean(post_nexp_raw) - mean(pre_nexp_raw), tolerance = 1e-10)

  # stderr = sd(diff) / sqrt(n)
  expect_equal(mc_se_exp, sd(post_exp_raw - pre_exp_raw) / sqrt(n_raw_exp),
               tolerance = 1e-10)
  expect_equal(mc_se_nexp, sd(post_nexp_raw - pre_nexp_raw) / sqrt(n_raw_nexp),
               tolerance = 1e-10)

  # statistic = estimate / stderr
  expect_equal(pt_exp, mc_exp / mc_se_exp, tolerance = 1e-10)
  expect_equal(pt_nexp, mc_nexp / mc_se_nexp, tolerance = 1e-10)
})

# --- Test 7b: All signed pathways match (cooper, d) ---
test_that("LIFECYCLE-RAW: all signed pathways match — cooper, d", {
  ref <- run_pathway(dat_raw, "means_sd_pre_post", "d", "cooper")

  # means_se (exact)
  alt_se <- run_pathway(dat_raw, "means_se_pre_post", "d", "cooper")
  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)

  # mean_change_sd (exact)
  alt_mc_sd <- run_pathway(dat_raw, "mean_change_sd", "d", "cooper")
  expect_equal(ref$es, alt_mc_sd$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_sd$se, tolerance = 1e-10)

  # mean_change_se (exact)
  alt_mc_se <- run_pathway(dat_raw, "mean_change_se", "d", "cooper")
  expect_equal(ref$es, alt_mc_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_se$se, tolerance = 1e-10)

  # mean_change_ci (exact — t.test CI uses same qt)
  alt_mc_ci <- run_pathway(dat_raw, "mean_change_ci", "d", "cooper")
  expect_equal(ref$es, alt_mc_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_ci$se, tolerance = 1e-10)

  # mean_change_pval (near-exact — pval roundtrip)
  alt_mc_pval <- run_pathway(dat_raw, "mean_change_pval", "d", "cooper")
  expect_equal(ref$es, alt_mc_pval$es, tolerance = 1e-6)
  expect_equal(ref$se, alt_mc_pval$se, tolerance = 1e-6)

  # paired_t (exact) — against the LEGACY per-arm means reference, because
  # paired_t cannot pool (see the note above run_pathway). The simulated arms
  # here have unequal SDs, so the pooled and per-arm standardizers differ.
  ref_perarm <- run_pathway(dat_raw, "means_sd_pre_post", "d", "cooper", pool_sd = FALSE)
  alt_t <- run_pathway(dat_raw, "paired_t", "d", "cooper", pool_sd = FALSE)
  expect_equal(ref_perarm$es, alt_t$es, tolerance = 1e-10)
  expect_equal(ref_perarm$se, alt_t$se, tolerance = 1e-10)
})

# --- Test 7c: Unsigned pathways match each other ---
test_that("LIFECYCLE-RAW: unsigned pathways match each other — cooper, d", {
  alt_tpval <- run_pathway(dat_raw, "paired_t_pval", "d", "cooper")
  alt_f     <- run_pathway(dat_raw, "paired_f", "d", "cooper")
  alt_fpval <- run_pathway(dat_raw, "paired_f_pval", "d", "cooper")

  # F = t^2, so paired_f and paired_t_pval give identical abs(d)
  expect_equal(abs(alt_f$es), abs(alt_tpval$es), tolerance = 1e-10)
  expect_equal(alt_f$se, alt_tpval$se, tolerance = 1e-10)

  # F_pval and t_pval give same unsigned per-group d (pval roundtrip)
  expect_equal(abs(alt_fpval$es), abs(alt_tpval$es), tolerance = 1e-6)
  expect_equal(alt_fpval$se, alt_tpval$se, tolerance = 1e-6)
})

# --- Test 7d: All signed pathways match (cooper, g) ---
test_that("LIFECYCLE-RAW: all signed pathways match — cooper, g", {
  ref <- run_pathway(dat_raw, "means_sd_pre_post", "g", "cooper")

  alt_se <- run_pathway(dat_raw, "means_se_pre_post", "g", "cooper")
  expect_equal(ref$es, alt_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_se$se, tolerance = 1e-10)

  alt_mc_sd <- run_pathway(dat_raw, "mean_change_sd", "g", "cooper")
  expect_equal(ref$es, alt_mc_sd$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_sd$se, tolerance = 1e-10)

  alt_mc_se <- run_pathway(dat_raw, "mean_change_se", "g", "cooper")
  expect_equal(ref$es, alt_mc_se$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_se$se, tolerance = 1e-10)

  alt_mc_ci <- run_pathway(dat_raw, "mean_change_ci", "g", "cooper")
  expect_equal(ref$es, alt_mc_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_mc_ci$se, tolerance = 1e-10)

  alt_mc_pval <- run_pathway(dat_raw, "mean_change_pval", "g", "cooper")
  expect_equal(ref$es, alt_mc_pval$es, tolerance = 1e-6)
  expect_equal(ref$se, alt_mc_pval$se, tolerance = 1e-6)

  # paired_t — LEGACY per-arm reference (paired_t cannot pool).
  ref_perarm <- run_pathway(dat_raw, "means_sd_pre_post", "g", "cooper", pool_sd = FALSE)
  alt_t <- run_pathway(dat_raw, "paired_t", "g", "cooper", pool_sd = FALSE)
  expect_equal(ref_perarm$es, alt_t$es, tolerance = 1e-10)
  expect_equal(ref_perarm$se, alt_t$se, tolerance = 1e-10)
})

# --- Test 7e: means_ci pathway matches (from raw-derived pre/post CIs) ---
test_that("LIFECYCLE-RAW: means_ci_pre_post matches reference — cooper, d", {
  ref <- run_pathway(dat_raw, "means_sd_pre_post", "d", "cooper")
  alt_ci <- run_pathway(dat_raw, "means_ci_pre_post", "d", "cooper")
  expect_equal(ref$es, alt_ci$es, tolerance = 1e-10)
  expect_equal(ref$se, alt_ci$se, tolerance = 1e-10)
})

} # end NOT_CRAN gate
