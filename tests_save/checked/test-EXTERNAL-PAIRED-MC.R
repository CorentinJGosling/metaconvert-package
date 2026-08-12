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

## 1.1 COOPER vs Morris & DeShon (2002) Manual Formula (DEFAULT per-arm path) ====
# NOTE: this test pins the per-arm standardizer (pool_sd = FALSE), which is the package
# DEFAULT: each arm is standardized by its OWN change SD, the two within-group d_rm values
# are subtracted, and their sampling variances are ADDED (the arms are independent). This is
# Morris (2008) d_ppc1 / Becker (1988) — the same construction metafor users implement as
# escalc(measure = "SMCR") per arm followed by `yi = yT - yC; vi = vT + vC`. The Morris &
# DeShon (2002) comparator below is built exactly that way (d_rm_exp - d_rm_nexp,
# var = var_rm_exp + var_rm_nexp), so it encodes the per-arm construction. pool_sd = FALSE is
# passed explicitly here for self-documentation (it is also the default; test 1.2 asserts that).
# The OPT-IN pool_sd = TRUE path (a single SD pooled across arms, Morris 2008 d_ppc2) is a
# different estimand and is covered separately in test 1.2b.
test_that("MC-SD: cooper two-group matches Morris & DeShon (2002) (d + d_se) (per-arm standardizer, pool_sd = FALSE)", {
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
      pre_post_to_smd = "cooper",
      # per-arm standardizer (the default): required for the comparator below, which
      # subtracts two separately-standardized within-group d_rm values
      pool_sd = FALSE
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

## 1.2 MORRIS_DZ vs metafor SMCC Manual Formula (DEFAULT per-arm path) ====
# NOTE: this test pins the per-arm standardizer (pool_sd = FALSE), which is the package
# DEFAULT. The comparator is metafor's SMCC (change-score) formula applied SEPARATELY TO EACH
# ARM, with the two within-group d_z values then SUBTRACTED (var = var_exp + var_nexp) — i.e.
# it encodes the per-arm construction. The per-arm SMCC variance (built from the corrected g:
# var_g = 1/n + g^2/(2n), var_d = var_g / J^2) is unchanged.
#
# This test ALSO regression-guards the DEFAULT itself: the wrapper is called a second time
# with pool_sd left unspecified, and the result must be identical to the explicit
# pool_sd = FALSE call. (pool_sd = TRUE, the arm-pooled Morris 2008 d_ppc2 estimand, is
# opt-in and is covered separately in test 1.2b.)
test_that("MC-SD: morris_dz two-group matches per-arm metafor SMCC formula (d + d_se); pool_sd defaults to FALSE (per-arm d_ppc1, Becker 1988 / Morris d_ppc1)", {
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
      pre_post_to_smd = "morris_dz",
      # per-arm standardizer (the default): required for the comparator below, which applies
      # SMCC to each arm separately and subtracts the two within-group d_z values
      pool_sd = FALSE
    )

    # DEFAULT GUARD: calling the wrapper WITHOUT pool_sd must reproduce the per-arm
    # (pool_sd = FALSE) result exactly — i.e. the default is the per-arm d_ppc1 path.
    mc_default <- es_from_mean_change_sd(
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
    expect_equal(mc_default$d, mc$d, tolerance = 1e-15,
                 label = paste0("Study ", i, ": default pool_sd is the per-arm path (d)"))
    expect_equal(mc_default$d_se, mc$d_se, tolerance = 1e-15,
                 label = paste0("Study ", i, ": default pool_sd is the per-arm path (d_se)"))

    # Manual per-arm metafor SMCC formula
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

## 1.2b MORRIS_DZ pooled (OPT-IN pool_sd = TRUE) vs live metafor::escalc(SMD) ====
# Covers the OPT-IN pooled path (pool_sd = TRUE, passed EXPLICITLY — it is NOT the default):
# a single standardizing SD pooled across the two arms, so the between-group estimate is
# (mean change difference) / (pooled change SD). This is Morris (2008) d_ppc2, and it is
# exactly a two-sample SMD computed on change scores — so it can be anchored against a LIVE
# metafor::escalc(measure = "SMD") call rather than a (circular) transcription of the source.
#
# The pooled variance formulas are correct and validated; they are simply no longer applied
# silently, because d_ppc2 assumes the two arms' true standardizing SDs are equal, whereas the
# default per-arm d_ppc1 (tests 1.1 / 1.2) does not. Choosing between them is a genuine
# analytic decision, so it is now the caller's.
#
#   - POINT ESTIMATE: exact. metaConvert's g equals metafor's yi to machine precision.
#   - SE: exact against metafor's default vtype = "LS" variance.
test_that("MC-SD: morris_dz pooled (opt-in pool_sd = TRUE) matches metafor::escalc(SMD) on change scores", {
  skip_if_not_installed("metafor")

  dat <- df.SMC[1:5, ]
  dat$n_exp <- dat$n_cases
  dat$n_nexp <- dat$n_controls

  dat$mean_change_exp <- dat$mean_cases - dat$mean_pre_cases
  dat$mean_change_nexp <- dat$mean_controls - dat$mean_pre_controls

  dat$r_pre_post_exp <- 0.55
  dat$r_pre_post_nexp <- 0.55

  dat$mean_change_sd_exp <- sqrt(dat$sd_pre_cases^2 + dat$sd_cases^2 -
                                  2 * dat$r_pre_post_exp * dat$sd_pre_cases * dat$sd_cases)
  dat$mean_change_sd_nexp <- sqrt(dat$sd_pre_controls^2 + dat$sd_controls^2 -
                                   2 * dat$r_pre_post_nexp * dat$sd_pre_controls * dat$sd_controls)

  .d_j <- function(df) exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))

  for (i in 1:nrow(dat)) {
    # metaConvert, OPT-IN pooled standardizer -- pool_sd = TRUE passed EXPLICITLY
    # (the default is pool_sd = FALSE, the per-arm path guarded in test 1.2).
    mc <- es_from_mean_change_sd(
      mean_change_exp = dat$mean_change_exp[i],
      mean_change_sd_exp = dat$mean_change_sd_exp[i],
      n_exp = dat$n_exp[i],
      mean_change_nexp = dat$mean_change_nexp[i],
      mean_change_sd_nexp = dat$mean_change_sd_nexp[i],
      n_nexp = dat$n_nexp[i],
      r_pre_post_exp = dat$r_pre_post_exp[i],
      r_pre_post_nexp = dat$r_pre_post_nexp[i],
      pre_post_to_smd = "morris_dz",
      pool_sd = TRUE
    )

    # Live external reference: two-sample SMD on the change scores, metafor's
    # DEFAULT vtype = "LS". Once the change SD is pooled across arms the pooled
    # morris_dz estimator IS an independent-groups Hedges g on the change scores,
    # so agreement is now EXACT in both the point estimate and the variance --
    # not merely close. (It was previously pinned to vtype = "LS2", the convention
    # that carried a spurious leading J^2 and understated vi by up to 9%.)
    mf <- metafor::escalc(
      measure = "SMD", vtype = "LS",
      m1i = dat$mean_change_exp[i], m2i = dat$mean_change_nexp[i],
      sd1i = dat$mean_change_sd_exp[i], sd2i = dat$mean_change_sd_nexp[i],
      n1i = dat$n_exp[i], n2i = dat$n_nexp[i]
    )
    g_ref <- as.numeric(mf$yi)
    se_ref <- sqrt(as.numeric(mf$vi))

    N <- dat$n_exp[i] + dat$n_nexp[i]
    J <- .d_j(N - 2)

    # Point estimates: EXACT agreement with metafor.
    expect_equal(mc$g, g_ref, tolerance = 1e-10,
                 label = paste0("Study ", i, ": pooled morris_dz g matches metafor SMD yi"))
    expect_equal(mc$d, g_ref / J, tolerance = 1e-10,
                 label = paste0("Study ", i, ": pooled morris_dz d matches metafor SMD yi / J"))

    # SE: BIT-EXACT agreement with metafor's default LS variance.
    expect_equal(mc$g_se, se_ref, tolerance = 1e-12,
                 label = paste0("Study ", i, ": pooled morris_dz g_se matches metafor SMD (LS)"))

    # Anti-regression: the OLD variance carried a spurious 2*(1-r_avg) factor on the
    # N/(n_exp*n_nexp) term. Assert we are NOT back on that formula (it is several percent off).
    m <- N - 2
    r_avg <- (dat$n_exp[i] * dat$r_pre_post_exp[i] +
                dat$n_nexp[i] * dat$r_pre_post_nexp[i]) / N
    old_g_se <- sqrt(J^2 * (2 * (1 - r_avg) * N / (dat$n_exp[i] * dat$n_nexp[i]) +
                              mc$g^2 / (2 * m)))
    expect_gt(abs(mc$g_se - old_g_se) / se_ref, 1e-2)
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
