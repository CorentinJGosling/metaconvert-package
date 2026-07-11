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
run_pathway <- function(data, hierarchy, measure = "d",
                        pre_post_to_smd = "cooper") {
  res <- convert_df(data,
    verbose = FALSE,
    es_selected = "hierarchy",
    hierarchy = hierarchy,
    measure = measure,
    pre_post_to_smd = pre_post_to_smd
  )
  summary(res, digits = 11)
}


test_that("MC+SE vs paired_t — df.SMC (cooper, d and g)", {
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

test_that("MC+CI vs paired_t — df.SMC (cooper, d and g)", {
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

test_that("MC+pval vs paired_t — df.SMC (cooper, d and g)", {
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

} # end NOT_CRAN gate
