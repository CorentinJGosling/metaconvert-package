# ==============================================================================
# PAIRED STATISTICS FORMAT EQUIVALENCE AND REVERSE TESTS
# ==============================================================================
# Tests that paired_t, paired_t_pval, paired_f, paired_f_pval pathways
# produce consistent results and that reverse flags correctly flip signs.
#
# Uses metaumbrella::df.SMC (33 rows of real data).
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

# Map column names
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

# Derive mean change + paired statistics
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

dat$paired_f_exp <- dat$paired_t_exp^2
dat$paired_f_nexp <- dat$paired_t_nexp^2

dat$paired_t_pval_exp <- 2 * pt(-abs(dat$paired_t_exp), dat$n_exp - 1)
dat$paired_t_pval_nexp <- 2 * pt(-abs(dat$paired_t_nexp), dat$n_nexp - 1)

dat$paired_f_pval_exp <- dat$paired_t_pval_exp
dat$paired_f_pval_nexp <- dat$paired_t_pval_nexp

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


# ==============================================================================
# FORWARD TESTS: paired statistics vs means_sd (signed pathways)
# ==============================================================================

test_that("D/G — paired_t vs means_sd_pre_post — cooper", {
  ref <- run_pathway(dat, "means_sd_pre_post", "d", "cooper")
  alt <- run_pathway(dat, "paired_t", "d", "cooper")
  expect_equal(unique(ref$info_used_crude), "means_sd_pre_post")
  expect_equal(unique(alt$info_used_crude), "paired_t")
  expect_equal(ref$es_crude, alt$es_crude, tolerance = 1e-10)
  expect_equal(ref$se_crude, alt$se_crude, tolerance = 1e-10)

  ref_g <- run_pathway(dat, "means_sd_pre_post", "g", "cooper")
  alt_g <- run_pathway(dat, "paired_t", "g", "cooper")
  expect_equal(ref_g$es_crude, alt_g$es_crude, tolerance = 1e-10)
  expect_equal(ref_g$se_crude, alt_g$se_crude, tolerance = 1e-10)
})

# NOTE: paired_t_pval, paired_f, paired_f_pval lose per-group sign information.
# For two-group designs d = |d_exp| - |d_nexp| instead of d_exp - d_nexp,
# so we compare these unsigned pathways against EACH OTHER, not against means_sd.

test_that("D — unsigned pathways consistent — cooper", {
  alt_tpval <- run_pathway(dat, "paired_t_pval", "d", "cooper")
  alt_f     <- run_pathway(dat, "paired_f", "d", "cooper")
  alt_fpval <- run_pathway(dat, "paired_f_pval", "d", "cooper")
  expect_equal(unique(alt_tpval$info_used_crude), "paired_t_pval")
  expect_equal(unique(alt_f$info_used_crude), "paired_f")
  expect_equal(unique(alt_fpval$info_used_crude), "paired_f_pval")

  # F = t^2 → same |t| → same unsigned d
  expect_equal(abs(alt_f$es_crude), abs(alt_tpval$es_crude), tolerance = 1e-10)
  expect_equal(alt_f$se_crude, alt_tpval$se_crude, tolerance = 1e-10)

  # pval roundtrip
  expect_equal(abs(alt_fpval$es_crude), abs(alt_tpval$es_crude), tolerance = 1e-6)
  expect_equal(alt_fpval$se_crude, alt_tpval$se_crude, tolerance = 1e-6)
})


# ==============================================================================
# REVERSE TESTS: setting reverse flag flips d/g/r/z/logor, preserves SE
# ==============================================================================

# Helper: run reverse test for a given pathway and reverse column name
run_reverse_test <- function(data, hierarchy, reverse_col) {
  measures <- c("d", "g", "r", "z", "logor")

  # Forward (no reverse)
  fwd <- lapply(measures, function(m) run_pathway(data, hierarchy, m, "cooper"))
  names(fwd) <- measures

  # Add reverse flag
  data_rv <- data
  data_rv[[reverse_col]] <- TRUE

  # Reversed
  rev <- lapply(measures, function(m) run_pathway(data_rv, hierarchy, m, "cooper"))
  names(rev) <- measures

  list(fwd = fwd, rev = rev)
}

test_that("REVERSE — paired_t: es flips, SE preserved", {
  rv <- run_reverse_test(dat, "paired_t", "reverse_paired_t")
  expect_equal(rv$fwd$d$info_used_crude, rv$rev$d$info_used_crude)
  for (m in c("d", "g", "r", "z", "logor")) {
    expect_equal(rv$fwd[[m]]$es_crude, -rv$rev[[m]]$es_crude, tolerance = 1e-10,
                 info = paste("paired_t reverse", m, "es"))
    expect_equal(rv$fwd[[m]]$se_crude, rv$rev[[m]]$se_crude, tolerance = 1e-10,
                 info = paste("paired_t reverse", m, "se"))
  }
})

test_that("REVERSE — paired_t_pval: es flips, SE preserved", {
  rv <- run_reverse_test(dat, "paired_t_pval", "reverse_paired_t_pval")
  expect_equal(unique(rv$fwd$d$info_used_crude), "paired_t_pval")
  for (m in c("d", "g", "r", "z", "logor")) {
    expect_equal(rv$fwd[[m]]$es_crude, -rv$rev[[m]]$es_crude, tolerance = 1e-10,
                 info = paste("paired_t_pval reverse", m, "es"))
    expect_equal(rv$fwd[[m]]$se_crude, rv$rev[[m]]$se_crude, tolerance = 1e-10,
                 info = paste("paired_t_pval reverse", m, "se"))
  }
})

test_that("REVERSE — paired_f: es flips, SE preserved", {
  rv <- run_reverse_test(dat, "paired_f", "reverse_paired_f")
  expect_equal(unique(rv$fwd$d$info_used_crude), "paired_f")
  for (m in c("d", "g", "r", "z", "logor")) {
    expect_equal(rv$fwd[[m]]$es_crude, -rv$rev[[m]]$es_crude, tolerance = 1e-10,
                 info = paste("paired_f reverse", m, "es"))
    expect_equal(rv$fwd[[m]]$se_crude, rv$rev[[m]]$se_crude, tolerance = 1e-10,
                 info = paste("paired_f reverse", m, "se"))
  }
})

test_that("REVERSE — paired_f_pval: es flips, SE preserved", {
  rv <- run_reverse_test(dat, "paired_f_pval", "reverse_paired_f_pval")
  expect_equal(unique(rv$fwd$d$info_used_crude), "paired_f_pval")
  for (m in c("d", "g", "r", "z", "logor")) {
    expect_equal(rv$fwd[[m]]$es_crude, -rv$rev[[m]]$es_crude, tolerance = 1e-10,
                 info = paste("paired_f_pval reverse", m, "es"))
    expect_equal(rv$fwd[[m]]$se_crude, rv$rev[[m]]$se_crude, tolerance = 1e-10,
                 info = paste("paired_f_pval reverse", m, "se"))
  }
})
