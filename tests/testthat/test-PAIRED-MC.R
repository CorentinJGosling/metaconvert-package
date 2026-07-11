# Long-running file (> 25 s): executed locally and in CI (devtools::test and
# testthat set NOT_CRAN=true); skipped wholesale on CRAN to respect check-time
# limits. The fast pre/post files still run on CRAN.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
# ==============================================================================
# PAIRED MEAN CHANGE FORMAT EQUIVALENCE TESTS
# ==============================================================================
# Tests that mean_change_se, mean_change_ci, and mean_change_pval pathways
# produce the same results as paired_t (all derive from the same underlying
# statistics: mean change, SE of change, t = MC/SE, CI = MC ± t*SE, pval from t).
#
# SCOPE: these are equivalence tests of the LEGACY per-arm standardizer
# (pool_sd = FALSE), NOT of the current pooled default (pool_sd = TRUE).
#
# Why they must be: a paired t (or its CI/p-value re-expressions) identifies each
# arm's mean_change / sd_change, but it does NOT identify the two arms' SD ratio,
# so the pooled standardizer sqrt(((n1-1)*sd_c1^2 + (n2-1)*sd_c2^2)/(N-2)) is
# mathematically unrecoverable from a paired t. es_from_paired_t() therefore has
# no pool_sd argument and always uses the per-arm construction (standardize each
# arm by its OWN change SD, then subtract). The only pool_sd setting under which
# "mean_change_* == paired_t" can hold is FALSE, so the mean-change side is
# pinned to pool_sd = FALSE below. The pooled default (pool_sd = TRUE) is a
# different, non-equivalent estimand; it is covered by the two tests at the
# bottom of this file, which check the same format equivalence among the
# mean-change routes under the shipped default and assert that the two
# standardizers really do disagree on this data.
#
# Uses metaumbrella::df.SMC (33 rows of real data) instead of rnorm simulations.
# ==============================================================================

library(testthat)
library(metaConvert)

has_metaumbrella <- requireNamespace("metaumbrella", quietly = TRUE)
if (!has_metaumbrella) {
  skip("metaumbrella package not available")
}

# --- Data preparation (shared across all tests) ---
data("df.SMC", package = "metaumbrella")
dat <- df.SMC

# Map column names to metaConvert conventions
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

# Derive mean change statistics
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

dat$paired_t_exp <- dat$mean_change_exp / dat$mean_change_se_exp
dat$paired_t_nexp <- dat$mean_change_nexp / dat$mean_change_se_nexp

dat$mean_change_ci_lo_exp <- dat$mean_change_exp - qt(0.975, dat$n_exp - 1) * dat$mean_change_se_exp
dat$mean_change_ci_up_exp <- dat$mean_change_exp + qt(0.975, dat$n_exp - 1) * dat$mean_change_se_exp
dat$mean_change_ci_lo_nexp <- dat$mean_change_nexp - qt(0.975, dat$n_nexp - 1) * dat$mean_change_se_nexp
dat$mean_change_ci_up_nexp <- dat$mean_change_nexp + qt(0.975, dat$n_nexp - 1) * dat$mean_change_se_nexp

dat$mean_change_pval_exp <- 2 * pt(-abs(dat$paired_t_exp), dat$n_exp - 1)
dat$mean_change_pval_nexp <- 2 * pt(-abs(dat$paired_t_nexp), dat$n_nexp - 1)

# Remove original df.SMC columns to avoid interference
cols_to_remove <- c("n_cases", "n_controls", "mean_pre_cases", "mean_cases",
                    "sd_pre_cases", "sd_cases", "mean_pre_controls", "mean_controls",
                    "sd_pre_controls", "sd_controls", "pre_post_cor",
                    "value", "se", "measure")
dat <- dat[, setdiff(names(dat), cols_to_remove)]

# Helper function
#
# pool_sd = FALSE selects the legacy per-arm standardizer for the mean_change_*
# routes, which is the construction es_from_paired_t() is hard-wired to (see the
# header note). It is a no-op for the "paired_t" hierarchy itself, since that
# route takes no pool_sd argument -- passing it here simply keeps both sides of
# every comparison on one explicit, named construction.
#
# Results are memoized: every convert_df() call runs ALL ~76 methods over the 33
# rows, and the tests below request the same pathway several times (paired_t/d is
# wanted by three of them). The cache is keyed on the arguments that vary; every
# caller in this file passes the same `dat`, so `data` is deliberately not part
# of the key -- pass a different data frame and you must clear the cache.
.pathway_cache <- new.env(parent = emptyenv())

run_pathway <- function(data, hierarchy, measure = "d",
                        pre_post_to_smd = "cooper",
                        pool_sd = FALSE) {
  key <- paste(hierarchy, measure, pre_post_to_smd, pool_sd, sep = "|")
  hit <- .pathway_cache[[key]]
  if (!is.null(hit)) {
    return(hit)
  }

  res <- convert_df(data,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = hierarchy,
    measure = measure,
    pre_post_to_smd = pre_post_to_smd,
    pool_sd = pool_sd
  )
  out <- summary(res, digits = 11)
  .pathway_cache[[key]] <- out
  out
}


test_that("MC+SE vs paired_t — df.SMC (cooper, d and g; legacy per-arm standardizer, pool_sd = FALSE)", {
  ref_d <- run_pathway(dat, "paired_t", "d", "cooper")
  alt_d <- run_pathway(dat, "mean_change_se", "d", "cooper")
  expect_equal(unique(ref_d$info_used_crude), "paired_t")
  expect_equal(unique(alt_d$info_used_crude), "mean_change_se")
  expect_equal(ref_d$es_crude, alt_d$es_crude, tolerance = 1e-10)
  expect_equal(ref_d$se_crude, alt_d$se_crude, tolerance = 1e-10)

  ref_g <- run_pathway(dat, "paired_t", "g", "cooper")
  alt_g <- run_pathway(dat, "mean_change_se", "g", "cooper")
  expect_equal(ref_g$es_crude, alt_g$es_crude, tolerance = 1e-10)
  expect_equal(ref_g$se_crude, alt_g$se_crude, tolerance = 1e-10)
})

test_that("MC+CI vs paired_t — df.SMC (cooper, d and g; legacy per-arm standardizer, pool_sd = FALSE)", {
  ref_d <- run_pathway(dat, "paired_t", "d", "cooper")
  alt_d <- run_pathway(dat, "mean_change_ci", "d", "cooper")
  expect_equal(unique(alt_d$info_used_crude), "mean_change_ci")
  expect_equal(ref_d$es_crude, alt_d$es_crude, tolerance = 1e-10)
  expect_equal(ref_d$se_crude, alt_d$se_crude, tolerance = 1e-10)

  ref_g <- run_pathway(dat, "paired_t", "g", "cooper")
  alt_g <- run_pathway(dat, "mean_change_ci", "g", "cooper")
  expect_equal(ref_g$es_crude, alt_g$es_crude, tolerance = 1e-10)
  expect_equal(ref_g$se_crude, alt_g$se_crude, tolerance = 1e-10)
})

test_that("MC+pval vs paired_t — df.SMC (cooper, d and g; legacy per-arm standardizer, pool_sd = FALSE)", {
  ref_d <- run_pathway(dat, "paired_t", "d", "cooper")
  alt_d <- run_pathway(dat, "mean_change_pval", "d", "cooper")
  expect_equal(unique(alt_d$info_used_crude), "mean_change_pval")
  # pval roundtrip loses some precision
  expect_equal(ref_d$es_crude, alt_d$es_crude, tolerance = 1e-6)
  expect_equal(ref_d$se_crude, alt_d$se_crude, tolerance = 1e-6)

  ref_g <- run_pathway(dat, "paired_t", "g", "cooper")
  alt_g <- run_pathway(dat, "mean_change_pval", "g", "cooper")
  expect_equal(ref_g$es_crude, alt_g$es_crude, tolerance = 1e-6)
  expect_equal(ref_g$se_crude, alt_g$se_crude, tolerance = 1e-6)
})

# The three tests above pin the mean-change routes to pool_sd = FALSE, because
# that is the only construction a paired t can be compared against (see header).
# The format-equivalence property they check -- SD, SE, CI and p-value are four
# encodings of the same change statistics, so they must yield one ES -- is
# independent of the standardizer, and metaConvert now ships pool_sd = TRUE by
# default. Re-check the equivalence under the shipped default, using
# mean_change_sd (which the paired-t route cannot stand in for) as the reference.
test_that("MC format equivalence holds under the pooled default (pool_sd = TRUE)", {
  ref_d <- run_pathway(dat, "mean_change_sd", "d", "cooper", pool_sd = TRUE)
  expect_equal(unique(ref_d$info_used_crude), "mean_change_sd")

  for (m in c("d", "g")) {
    ref <- run_pathway(dat, "mean_change_sd", m, "cooper", pool_sd = TRUE)

    se_alt <- run_pathway(dat, "mean_change_se", m, "cooper", pool_sd = TRUE)
    expect_equal(unique(se_alt$info_used_crude), "mean_change_se")
    expect_equal(ref$es_crude, se_alt$es_crude, tolerance = 1e-10)
    expect_equal(ref$se_crude, se_alt$se_crude, tolerance = 1e-10)

    ci_alt <- run_pathway(dat, "mean_change_ci", m, "cooper", pool_sd = TRUE)
    expect_equal(unique(ci_alt$info_used_crude), "mean_change_ci")
    expect_equal(ref$es_crude, ci_alt$es_crude, tolerance = 1e-10)
    expect_equal(ref$se_crude, ci_alt$se_crude, tolerance = 1e-10)

    # pval roundtrip loses some precision
    pv_alt <- run_pathway(dat, "mean_change_pval", m, "cooper", pool_sd = TRUE)
    expect_equal(unique(pv_alt$info_used_crude), "mean_change_pval")
    expect_equal(ref$es_crude, pv_alt$es_crude, tolerance = 1e-6)
    expect_equal(ref$se_crude, pv_alt$se_crude, tolerance = 1e-6)
  }
})

# Guard the reason the tests above must pass pool_sd = FALSE: the pooled and
# per-arm standardizers are genuinely different estimands on df.SMC (whose arms
# have unequal change SDs), so the default must NOT silently reproduce the
# legacy per-arm numbers. If this ever starts failing, the pool_sd flip has been
# reverted or neutralised.
test_that("pooled and per-arm standardizers differ on df.SMC (unequal arm SDs)", {
  pooled  <- run_pathway(dat, "mean_change_sd", "g", "cooper", pool_sd = TRUE)
  perarm  <- run_pathway(dat, "mean_change_sd", "g", "cooper", pool_sd = FALSE)
  expect_false(isTRUE(all.equal(pooled$es_crude, perarm$es_crude, tolerance = 1e-6)))
})

} # end NOT_CRAN gate
