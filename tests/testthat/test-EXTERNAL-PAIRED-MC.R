# ==============================================================================
# EXTERNAL VALIDATION: es_from_PAIRED_MC.R
# ==============================================================================
# Validates all functions in es_from_PAIRED_MC.R against external packages:
#   - metafor (MD/MC for raw mean difference)
#   - Manual formulas (Morris & DeShon 2002 for cooper/morris_drm and morris_dz)
#
# Tests cover:
#   1. es_from_mean_change_sd() - cooper and morris_dz methods
#   2. es_from_mean_change_se() - cooper and morris_dz methods
#   3. es_from_mean_change_ci() - cooper and morris_dz methods
#   4. es_from_mean_change_pval() - cooper and morris_dz methods
#
# Dataset: metaumbrella::df.SMC (compute mean change from pre/post means)
# ==============================================================================

library(testthat)
library(metaConvert)

# Check for required packages
has_metafor <- requireNamespace("metafor", quietly = TRUE)
has_metaumbrella <- requireNamespace("metaumbrella", quietly = TRUE)

# Skip all tests if metaumbrella not available
if (!has_metaumbrella) {
  skip("metaumbrella package not available - skipping all EXTERNAL-PAIRED-MC tests")
}

# Load dataset and prepare mean change data
data("df.SMC", package = "metaumbrella")

# ==============================================================================
# SECTION 1: es_from_mean_change_sd() - TWO-GROUP
# ==============================================================================

## 1.1 COOPER vs Morris & DeShon (2002) Manual Formula ====
test_that("MC-SD: cooper two-group matches Morris & DeShon (2002) (d + d_se)", {
  dat <- df.SMC[1:5, ]
  dat$n_exp <- dat$n_cases
  dat$n_nexp <- dat$n_controls

  # Calculate mean change and SD of change
  dat$mean_change_exp <- dat$mean_cases - dat$mean_pre_cases
  dat$mean_change_nexp <- dat$mean_controls - dat$mean_pre_controls

  dat$r_pre_post_exp <- 0.6
  dat$r_pre_post_nexp <- 0.6

  dat$mean_change_sd_exp <- sqrt(dat$sd_pre_cases^2 + dat$sd_cases^2 -
                                  2 * dat$r_pre_post_exp * dat$sd_pre_cases * dat$sd_cases)
  dat$mean_change_sd_nexp <- sqrt(dat$sd_pre_controls^2 + dat$sd_controls^2 -
                                   2 * dat$r_pre_post_nexp * dat$sd_pre_controls * dat$sd_controls)

  for (i in 1:nrow(dat)) {
    # metaConvert
    mc <- es_from_mean_change_sd(
      mean_change_exp = dat$mean_change_exp[i],
      mean_change_sd_exp = dat$mean_change_sd_exp[i],
      n_exp = dat$n_exp[i],
      mean_change_nexp = dat$mean_change_nexp[i],
      mean_change_sd_nexp = dat$mean_change_sd_nexp[i],
      n_nexp = dat$n_nexp[i],
      r_pre_post_exp = dat$r_pre_post_exp[i],
      r_pre_post_nexp = dat$r_pre_post_nexp[i],
      pre_post_to_smd = "cooper"
    )

    # Manual Morris & DeShon (2002) calculation
    # Experimental group
    d_rm_exp <- dat$mean_change_exp[i] / dat$mean_change_sd_exp[i] *
                sqrt(2 * (1 - dat$r_pre_post_exp[i]))
    var_rm_exp <- 2 * (1 - dat$r_pre_post_exp[i]) / dat$n_exp[i] +
                  d_rm_exp^2 / (2 * dat$n_exp[i])

    # Control group
    d_rm_nexp <- dat$mean_change_nexp[i] / dat$mean_change_sd_nexp[i] *
                 sqrt(2 * (1 - dat$r_pre_post_nexp[i]))
    var_rm_nexp <- 2 * (1 - dat$r_pre_post_nexp[i]) / dat$n_nexp[i] +
                   d_rm_nexp^2 / (2 * dat$n_nexp[i])

    # Difference
    d_expected <- d_rm_exp - d_rm_nexp
    se_expected <- sqrt(var_rm_exp + var_rm_nexp)

    expect_equal(mc$d, d_expected, tolerance = 1e-10,
                 label = paste0("Study ", i, ": cooper d matches Morris & DeShon"))
    expect_equal(mc$d_se, se_expected, tolerance = 1e-10,
                 label = paste0("Study ", i, ": cooper d_se matches Morris & DeShon"))
  }
})

## 1.2 MORRIS_DZ vs metafor SMCC Manual Formula ====
test_that("MC-SD: morris_dz two-group matches metafor SMCC formula (d + d_se)", {
  dat <- df.SMC[1:5, ]
  dat$n_exp <- dat$n_cases
  dat$n_nexp <- dat$n_controls

  # Calculate mean change and SD of change
  dat$mean_change_exp <- dat$mean_cases - dat$mean_pre_cases
  dat$mean_change_nexp <- dat$mean_controls - dat$mean_pre_controls

  dat$r_pre_post_exp <- 0.55
  dat$r_pre_post_nexp <- 0.55

  dat$mean_change_sd_exp <- sqrt(dat$sd_pre_cases^2 + dat$sd_cases^2 -
                                  2 * dat$r_pre_post_exp * dat$sd_pre_cases * dat$sd_cases)
  dat$mean_change_sd_nexp <- sqrt(dat$sd_pre_controls^2 + dat$sd_controls^2 -
                                   2 * dat$r_pre_post_nexp * dat$sd_pre_controls * dat$sd_controls)

  .d_j <- function(df) exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2))

  for (i in 1:nrow(dat)) {
    # metaConvert
    mc <- es_from_mean_change_sd(
      mean_change_exp = dat$mean_change_exp[i],
      mean_change_sd_exp = dat$mean_change_sd_exp[i],
      n_exp = dat$n_exp[i],
      mean_change_nexp = dat$mean_change_nexp[i],
      mean_change_sd_nexp = dat$mean_change_sd_nexp[i],
      n_nexp = dat$n_nexp[i],
      r_pre_post_exp = dat$r_pre_post_exp[i],
      r_pre_post_nexp = dat$r_pre_post_nexp[i],
      pre_post_to_smd = "morris_dz"
    )

    # Manual metafor SMCC formula
    # Experimental group
    J_exp <- .d_j(dat$n_exp[i] - 1)
    d_z_exp <- dat$mean_change_exp[i] / dat$mean_change_sd_exp[i]
    g_z_exp <- d_z_exp * J_exp
    var_g_z_exp <- 1 / dat$n_exp[i] + g_z_exp^2 / (2 * dat$n_exp[i])
    var_d_z_exp <- var_g_z_exp / (J_exp^2)

    # Control group
    J_nexp <- .d_j(dat$n_nexp[i] - 1)
    d_z_nexp <- dat$mean_change_nexp[i] / dat$mean_change_sd_nexp[i]
    g_z_nexp <- d_z_nexp * J_nexp
    var_g_z_nexp <- 1 / dat$n_nexp[i] + g_z_nexp^2 / (2 * dat$n_nexp[i])
    var_d_z_nexp <- var_g_z_nexp / (J_nexp^2)

    # Difference
    d_expected <- d_z_exp - d_z_nexp
    se_expected <- sqrt(var_d_z_exp + var_d_z_nexp)

    expect_equal(mc$d, d_expected, tolerance = 1e-10,
                 label = paste0("Study ", i, ": morris_dz d matches SMCC formula"))
    expect_equal(mc$d_se, se_expected, tolerance = 1e-10,
                 label = paste0("Study ", i, ": morris_dz d_se matches SMCC formula"))
  }
})

## 1.3 RAW MD vs metafor MD ====
test_that("MC-SD: raw MD two-group matches metafor MD (md + md_se)", {
  skip_if_not(has_metafor, "metafor not available")

  dat <- df.SMC[1:5, ]
  dat$n_exp <- dat$n_cases
  dat$n_nexp <- dat$n_controls

  # Calculate mean change and SD of change
  dat$mean_change_exp <- dat$mean_cases - dat$mean_pre_cases
  dat$mean_change_nexp <- dat$mean_controls - dat$mean_pre_controls

  dat$r_pre_post_exp <- 0.7
  dat$r_pre_post_nexp <- 0.7

  dat$mean_change_sd_exp <- sqrt(dat$sd_pre_cases^2 + dat$sd_cases^2 -
                                  2 * dat$r_pre_post_exp * dat$sd_pre_cases * dat$sd_cases)
  dat$mean_change_sd_nexp <- sqrt(dat$sd_pre_controls^2 + dat$sd_controls^2 -
                                   2 * dat$r_pre_post_nexp * dat$sd_pre_controls * dat$sd_controls)

  for (i in 1:nrow(dat)) {
    # metaConvert
    mc <- es_from_mean_change_sd(
      mean_change_exp = dat$mean_change_exp[i],
      mean_change_sd_exp = dat$mean_change_sd_exp[i],
      n_exp = dat$n_exp[i],
      mean_change_nexp = dat$mean_change_nexp[i],
      mean_change_sd_nexp = dat$mean_change_sd_nexp[i],
      n_nexp = dat$n_nexp[i],
      r_pre_post_exp = dat$r_pre_post_exp[i],
      r_pre_post_nexp = dat$r_pre_post_nexp[i],
      pre_post_to_smd = "cooper" # Method doesn't affect MD
    )

    # metafor MD (two independent groups with mean difference as outcome)
    mf_data <- data.frame(
      m1i = dat$mean_change_exp[i],
      m2i = dat$mean_change_nexp[i],
      sd1i = dat$mean_change_sd_exp[i],
      sd2i = dat$mean_change_sd_nexp[i],
      n1i = dat$n_exp[i],
      n2i = dat$n_nexp[i]
    )
    mf <- metafor::escalc(measure = "MD", m1i = m1i, m2i = m2i,
                          sd1i = sd1i, sd2i = sd2i,
                          n1i = n1i, n2i = n2i, data = mf_data)

    expect_equal(mc$md, as.numeric(mf$yi), tolerance = 1e-6,
                 label = paste0("Study ", i, ": raw MD matches metafor MD"))
    expect_equal(mc$md_se, sqrt(as.numeric(mf$vi)), tolerance = 1e-6,
                 label = paste0("Study ", i, ": raw MD SE matches metafor MD"))
  }
})

## 1.4 EQUIVALENCE: Mean change ≡ Pre/Post with mean_pre=0 ====
test_that("MC-SD: mean change EXACTLY equals means_sd_pre_post(mean_pre=0)", {
  dat <- df.SMC[1:5, ]
  dat$n_exp <- dat$n_cases
  dat$n_nexp <- dat$n_controls

  # Calculate mean change and SD of change
  dat$mean_change_exp <- dat$mean_cases - dat$mean_pre_cases
  dat$mean_change_nexp <- dat$mean_controls - dat$mean_pre_controls

  dat$r_pre_post_exp <- 0.65
  dat$r_pre_post_nexp <- 0.65

  dat$mean_change_sd_exp <- sqrt(dat$sd_pre_cases^2 + dat$sd_cases^2 -
                                  2 * dat$r_pre_post_exp * dat$sd_pre_cases * dat$sd_cases)
  dat$mean_change_sd_nexp <- sqrt(dat$sd_pre_controls^2 + dat$sd_controls^2 -
                                   2 * dat$r_pre_post_nexp * dat$sd_pre_controls * dat$sd_controls)

  methods <- c("cooper", "morris_dz")

  for (method in methods) {
    for (i in 1:nrow(dat)) {
      # From mean change
      mc_change <- es_from_mean_change_sd(
        mean_change_exp = dat$mean_change_exp[i],
        mean_change_sd_exp = dat$mean_change_sd_exp[i],
        n_exp = dat$n_exp[i],
        mean_change_nexp = dat$mean_change_nexp[i],
        mean_change_sd_nexp = dat$mean_change_sd_nexp[i],
        n_nexp = dat$n_nexp[i],
        r_pre_post_exp = dat$r_pre_post_exp[i],
        r_pre_post_nexp = dat$r_pre_post_nexp[i],
        pre_post_to_smd = method
      )

      # From pre/post with mean_pre=0
      mc_prepost <- es_from_means_sd_pre_post(
        mean_pre_exp = 0,
        mean_exp = dat$mean_change_exp[i],
        mean_pre_sd_exp = 0,
        mean_sd_exp = dat$mean_change_sd_exp[i],
        mean_pre_nexp = 0,
        mean_nexp = dat$mean_change_nexp[i],
        mean_pre_sd_nexp = 0,
        mean_sd_nexp = dat$mean_change_sd_nexp[i],
        n_exp = dat$n_exp[i],
        n_nexp = dat$n_nexp[i],
        r_pre_post_exp = dat$r_pre_post_exp[i],
        r_pre_post_nexp = dat$r_pre_post_nexp[i],
        pre_post_to_smd = method
      )

      expect_equal(mc_change$d, mc_prepost$d, tolerance = 1e-15,
                   label = paste0(method, " study ", i, ": mean_change ≡ pre/post (d)"))
      expect_equal(mc_change$d_se, mc_prepost$d_se, tolerance = 1e-15,
                   label = paste0(method, " study ", i, ": mean_change ≡ pre/post (d_se)"))
      expect_equal(mc_change$g, mc_prepost$g, tolerance = 1e-15,
                   label = paste0(method, " study ", i, ": mean_change ≡ pre/post (g)"))
      expect_equal(mc_change$g_se, mc_prepost$g_se, tolerance = 1e-15,
                   label = paste0(method, " study ", i, ": mean_change ≡ pre/post (g_se)"))
    }
  }
})

# ==============================================================================
# SECTION 2: es_from_mean_change_se() - TWO-GROUP
# ==============================================================================

test_that("MC-SE: all methods match SD version after SE→SD conversion", {
  dat <- df.SMC[1:3, ]
  dat$n_exp <- dat$n_cases
  dat$n_nexp <- dat$n_controls

  # Calculate mean change and SD of change
  dat$mean_change_exp <- dat$mean_cases - dat$mean_pre_cases
  dat$mean_change_nexp <- dat$mean_controls - dat$mean_pre_controls

  dat$r_pre_post_exp <- 0.6
  dat$r_pre_post_nexp <- 0.6

  dat$mean_change_sd_exp <- sqrt(dat$sd_pre_cases^2 + dat$sd_cases^2 -
                                  2 * dat$r_pre_post_exp * dat$sd_pre_cases * dat$sd_cases)
  dat$mean_change_sd_nexp <- sqrt(dat$sd_pre_controls^2 + dat$sd_controls^2 -
                                   2 * dat$r_pre_post_nexp * dat$sd_pre_controls * dat$sd_controls)

  # Convert SD to SE
  dat$mean_change_se_exp <- dat$mean_change_sd_exp / sqrt(dat$n_exp)
  dat$mean_change_se_nexp <- dat$mean_change_sd_nexp / sqrt(dat$n_nexp)

  methods <- c("cooper", "morris_dz")

  for (method in methods) {
    for (i in 1:nrow(dat)) {
      # From SD
      mc_sd <- es_from_mean_change_sd(
        mean_change_exp = dat$mean_change_exp[i],
        mean_change_sd_exp = dat$mean_change_sd_exp[i],
        n_exp = dat$n_exp[i],
        mean_change_nexp = dat$mean_change_nexp[i],
        mean_change_sd_nexp = dat$mean_change_sd_nexp[i],
        n_nexp = dat$n_nexp[i],
        r_pre_post_exp = dat$r_pre_post_exp[i],
        r_pre_post_nexp = dat$r_pre_post_nexp[i],
        pre_post_to_smd = method
      )

      # From SE
      mc_se <- es_from_mean_change_se(
        mean_change_exp = dat$mean_change_exp[i],
        mean_change_se_exp = dat$mean_change_se_exp[i],
        n_exp = dat$n_exp[i],
        mean_change_nexp = dat$mean_change_nexp[i],
        mean_change_se_nexp = dat$mean_change_se_nexp[i],
        n_nexp = dat$n_nexp[i],
        r_pre_post_exp = dat$r_pre_post_exp[i],
        r_pre_post_nexp = dat$r_pre_post_nexp[i],
        pre_post_to_smd = method
      )

      expect_equal(mc_se$d, mc_sd$d, tolerance = 1e-10,
                   label = paste0(method, " study ", i, ": SE→SD gives same d"))
      expect_equal(mc_se$d_se, mc_sd$d_se, tolerance = 1e-10,
                   label = paste0(method, " study ", i, ": SE→SD gives same d_se"))
    }
  }
})

# ==============================================================================
# SECTION 3: es_from_mean_change_ci() - TWO-GROUP
# ==============================================================================

test_that("MC-CI: all methods match SD version after CI→SD conversion", {
  dat <- df.SMC[1:3, ]
  dat$n_exp <- dat$n_cases
  dat$n_nexp <- dat$n_controls

  # Calculate mean change and SD of change
  dat$mean_change_exp <- dat$mean_cases - dat$mean_pre_cases
  dat$mean_change_nexp <- dat$mean_controls - dat$mean_pre_controls

  dat$r_pre_post_exp <- 0.6
  dat$r_pre_post_nexp <- 0.6

  dat$mean_change_sd_exp <- sqrt(dat$sd_pre_cases^2 + dat$sd_cases^2 -
                                  2 * dat$r_pre_post_exp * dat$sd_pre_cases * dat$sd_cases)
  dat$mean_change_sd_nexp <- sqrt(dat$sd_pre_controls^2 + dat$sd_controls^2 -
                                   2 * dat$r_pre_post_nexp * dat$sd_pre_controls * dat$sd_controls)

  # Convert SD to CI
  dat$mean_change_se_exp <- dat$mean_change_sd_exp / sqrt(dat$n_exp)
  dat$mean_change_se_nexp <- dat$mean_change_sd_nexp / sqrt(dat$n_nexp)

  dat$mean_change_ci_lo_exp <- dat$mean_change_exp - qt(0.975, dat$n_exp - 1) * dat$mean_change_se_exp
  dat$mean_change_ci_up_exp <- dat$mean_change_exp + qt(0.975, dat$n_exp - 1) * dat$mean_change_se_exp
  dat$mean_change_ci_lo_nexp <- dat$mean_change_nexp - qt(0.975, dat$n_nexp - 1) * dat$mean_change_se_nexp
  dat$mean_change_ci_up_nexp <- dat$mean_change_nexp + qt(0.975, dat$n_nexp - 1) * dat$mean_change_se_nexp

  methods <- c("cooper", "morris_dz")

  for (method in methods) {
    for (i in 1:nrow(dat)) {
      # From SD
      mc_sd <- es_from_mean_change_sd(
        mean_change_exp = dat$mean_change_exp[i],
        mean_change_sd_exp = dat$mean_change_sd_exp[i],
        n_exp = dat$n_exp[i],
        mean_change_nexp = dat$mean_change_nexp[i],
        mean_change_sd_nexp = dat$mean_change_sd_nexp[i],
        n_nexp = dat$n_nexp[i],
        r_pre_post_exp = dat$r_pre_post_exp[i],
        r_pre_post_nexp = dat$r_pre_post_nexp[i],
        pre_post_to_smd = method
      )

      # From CI
      mc_ci <- es_from_mean_change_ci(
        mean_change_exp = dat$mean_change_exp[i],
        mean_change_ci_lo_exp = dat$mean_change_ci_lo_exp[i],
        mean_change_ci_up_exp = dat$mean_change_ci_up_exp[i],
        mean_change_nexp = dat$mean_change_nexp[i],
        mean_change_ci_lo_nexp = dat$mean_change_ci_lo_nexp[i],
        mean_change_ci_up_nexp = dat$mean_change_ci_up_nexp[i],
        n_exp = dat$n_exp[i],
        n_nexp = dat$n_nexp[i],
        r_pre_post_exp = dat$r_pre_post_exp[i],
        r_pre_post_nexp = dat$r_pre_post_nexp[i],
        pre_post_to_smd = method
      )

      expect_equal(mc_ci$d, mc_sd$d, tolerance = 1e-9,
                   label = paste0(method, " study ", i, ": CI→SD gives same d"))
      expect_equal(mc_ci$d_se, mc_sd$d_se, tolerance = 1e-9,
                   label = paste0(method, " study ", i, ": CI→SD gives same d_se"))
    }
  }
})

# ==============================================================================
# SECTION 4: es_from_mean_change_pval() - TWO-GROUP
# ==============================================================================

test_that("MC-PVAL: all methods match SD version after pval→t→SE→SD conversion", {
  dat <- df.SMC[1:3, ]
  dat$n_exp <- dat$n_cases
  dat$n_nexp <- dat$n_controls

  # Calculate mean change and SD of change
  dat$mean_change_exp <- dat$mean_cases - dat$mean_pre_cases
  dat$mean_change_nexp <- dat$mean_controls - dat$mean_pre_controls

  dat$r_pre_post_exp <- 0.6
  dat$r_pre_post_nexp <- 0.6

  dat$mean_change_sd_exp <- sqrt(dat$sd_pre_cases^2 + dat$sd_cases^2 -
                                  2 * dat$r_pre_post_exp * dat$sd_pre_cases * dat$sd_cases)
  dat$mean_change_sd_nexp <- sqrt(dat$sd_pre_controls^2 + dat$sd_controls^2 -
                                   2 * dat$r_pre_post_nexp * dat$sd_pre_controls * dat$sd_controls)

  # Calculate t-statistics and p-values
  dat$t_exp <- dat$mean_change_exp / (dat$mean_change_sd_exp / sqrt(dat$n_exp))
  dat$t_nexp <- dat$mean_change_nexp / (dat$mean_change_sd_nexp / sqrt(dat$n_nexp))
  dat$pval_exp <- 2 * pt(abs(dat$t_exp), dat$n_exp - 1, lower.tail = FALSE)
  dat$pval_nexp <- 2 * pt(abs(dat$t_nexp), dat$n_nexp - 1, lower.tail = FALSE)

  methods <- c("cooper", "morris_dz")

  for (method in methods) {
    for (i in 1:nrow(dat)) {
      # From SD
      mc_sd <- es_from_mean_change_sd(
        mean_change_exp = dat$mean_change_exp[i],
        mean_change_sd_exp = dat$mean_change_sd_exp[i],
        n_exp = dat$n_exp[i],
        mean_change_nexp = dat$mean_change_nexp[i],
        mean_change_sd_nexp = dat$mean_change_sd_nexp[i],
        n_nexp = dat$n_nexp[i],
        r_pre_post_exp = dat$r_pre_post_exp[i],
        r_pre_post_nexp = dat$r_pre_post_nexp[i],
        pre_post_to_smd = method
      )

      # From p-value
      mc_pval <- es_from_mean_change_pval(
        mean_change_exp = dat$mean_change_exp[i],
        mean_change_pval_exp = dat$pval_exp[i],
        mean_change_nexp = dat$mean_change_nexp[i],
        mean_change_pval_nexp = dat$pval_nexp[i],
        n_exp = dat$n_exp[i],
        n_nexp = dat$n_nexp[i],
        r_pre_post_exp = dat$r_pre_post_exp[i],
        r_pre_post_nexp = dat$r_pre_post_nexp[i],
        pre_post_to_smd = method
      )

      expect_equal(mc_pval$d, mc_sd$d, tolerance = 1e-9,
                   label = paste0(method, " study ", i, ": pval→SD gives same d"))
      expect_equal(mc_pval$d_se, mc_sd$d_se, tolerance = 1e-9,
                   label = paste0(method, " study ", i, ": pval→SD gives same d_se"))
    }
  }
})
