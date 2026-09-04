# smd_var = "hedges_olkin" on every pre/post call below (2.1.0). These tests pin
# agreement with metafor's SMCC/SMCR/SMCRH/SMCRPH variances and the Bonett (2008) /
# Morris & DeShon (2002) formulas, which are the Hedges-Olkin ("LS") form. Since 2.1.0
# the pre/post routes honour smd_var like the two-group routes, and the package
# default "borenstein" applies J^2 to the d-scale variance instead; the metafor form is
# reached with smd_var = "hedges_olkin", which is what these pins now request. The
# exception is the cooper / morris_drm blocks, whose Morris & DeShon (2002) oracle is
# the d-scale variance 2(1-r)/n + d^2/(2n): that is the J^2-outside ("borenstein")
# convention, so those calls request the default explicitly.

# Formula Verification Tests Against Primary Literature
# Validates that metaConvert implementations match published formulas
# References: Bonett (2008), Cooper (2019), Morris & DeShon (2002)

library(testthat)
library(metaConvert)

# Helper function to calculate Hedges' correction factor J
.calculate_J <- function(df) {
  exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
}

# Test 1: Bonett (2008) Formula Verification - Single Group ====
test_that("bonett formula matches Bonett (2008) for single-group", {
  # Variance formula (SMCRH): Bonett, D. G. (2008). Confidence intervals for
  # standardized linear contrasts of means. Psychological Methods, 13(2),
  # 99-109. https://doi.org/10.1037/1082-989X.13.2.99
  # Pretest-posttest-control design context: Morris, S. B. (2008). Estimating
  # effect sizes from pretest-posttest-control group designs. Organizational
  # Research Methods, 11(2), 364-386. https://doi.org/10.1177/1094428106291059

  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 12
  n <- 30
  r <- 0.6

  # Manual calculation per Bonett (2008)
  # Point estimate: d = (mean_post - mean_pre) / sd_pre
  d_bonett_manual <- (mean_post - mean_pre) / sd_pre

  # Calculate sd_change (used in variance formula)
  sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)

  # Hedges' correction factor
  J <- .calculate_J(n - 1)
  g_bonett_manual <- d_bonett_manual * J

  # Variance formula from Bonett (2008):
  # Implementation calculates var_g first, then var_d = var_g / J²
  # var(g) = sd_change² / (sd_pre² × (n-1)) + g² / (2 × (n-1))
  var_g_manual <- sd_change^2 / (sd_pre^2 * (n - 1)) +
                  g_bonett_manual^2 / (2 * (n - 1))
  var_bonett_manual <- var_g_manual / (J^2)
  se_bonett_manual <- sqrt(var_bonett_manual)

  # metaConvert calculation
  result <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre_exp = mean_pre,
    mean_exp = mean_post,
    mean_pre_sd_exp = sd_pre,
    mean_sd_exp = sd_post,
    n_exp = n,
    r_pre_post_exp = r,
    pre_post_to_smd = "bonett"
  )

  # Verify against manual calculation
  expect_equal(result$d, d_bonett_manual, tolerance = 1e-10,
               label = "Cohen's d should match Bonett (2008) formula")
  expect_equal(result$d_se, se_bonett_manual, tolerance = 1e-10,
               label = "SE should match Bonett (2008) variance formula")

  # Verify Hedges' g
  expect_equal(result$g, g_bonett_manual, tolerance = 1e-10,
               label = "Hedges' g should equal d × J")
})

test_that("bonett formula with different r values", {
  # Test that formula works correctly across range of r values
  mean_pre <- 45
  mean_post <- 52
  sd_pre <- 8
  sd_post <- 9
  n <- 25

  r_values <- c(0, 0.3, 0.5, 0.7, 0.9)

  for (r in r_values) {
    # Manual calculation
    d_manual <- (mean_post - mean_pre) / sd_pre
    sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
    J <- .calculate_J(n - 1)
    g_manual <- d_manual * J
    var_g_manual <- sd_change^2 / (sd_pre^2 * (n - 1)) + g_manual^2 / (2 * (n - 1))
    var_manual <- var_g_manual / (J^2)

    # metaConvert
    result <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
      mean_pre, mean_post, sd_pre, sd_post, n, r,
      pre_post_to_smd = "bonett"
    )

    expect_equal(result$d, d_manual, tolerance = 1e-10,
                 label = paste("d correct for r =", r))
    expect_equal(result$d_se^2, var_manual, tolerance = 1e-10,
                 label = paste("Variance correct for r =", r))
  }
})

test_that("bonett formula with variance heterogeneity", {
  # Test that formula handles different sd_pre vs sd_post
  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 8
  sd_post <- 16  # 2x larger
  n <- 35
  r <- 0.65

  # Manual calculation
  d_manual <- (mean_post - mean_pre) / sd_pre  # Uses sd_pre only
  sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
  J <- .calculate_J(n - 1)
  g_manual <- d_manual * J
  var_g_manual <- sd_change^2 / (sd_pre^2 * (n - 1)) + g_manual^2 / (2 * (n - 1))
  var_manual <- var_g_manual / (J^2)

  result <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  expect_equal(result$d, d_manual, tolerance = 1e-10)
  expect_equal(result$d_se^2, var_manual, tolerance = 1e-10)
})

# Test 2: Bonett (2008) Formula Verification - Two Groups ====
test_that("bonett two-group matches Bonett (2008) (per-arm standardizer, the pool_sd = FALSE default)", {
  # PER-ARM PATH (pool_sd = FALSE) -- this is the DEFAULT.
  # The manual calculation below standardizes EACH ARM by ITS OWN baseline SD, then
  # subtracts the two within-group values and ADDS their sampling variances (the arms
  # are independent). That is Morris (2008) d_ppc1 / Becker (1988), and it is exactly
  # what pool_sd = FALSE does. pool_sd is passed explicitly below for clarity, but it
  # is the default. The opt-in alternative (pool_sd = TRUE, a single standardizing SD
  # pooled across arms = Morris 2008 d_ppc2) is covered by the "pooled standardizer"
  # tests further down; it is more efficient but assumes the two arms' true SDs are
  # equal, so the package does not make that assumption silently.
  #
  # Two-group pre-post design
  mean_pre_exp <- 50
  mean_exp <- 60
  sd_pre_exp <- 10
  sd_post_exp <- 12
  n_exp <- 30
  r_exp <- 0.6

  mean_pre_nexp <- 51
  mean_post_nexp <- 53
  sd_pre_nexp <- 9
  sd_post_nexp <- 11
  n_nexp <- 28
  r_nexp <- 0.65

  # Manual calculation for experimental group
  d_exp_manual <- (mean_exp - mean_pre_exp) / sd_pre_exp
  sd_change_exp <- sqrt(sd_pre_exp^2 + sd_post_exp^2 - 2 * r_exp * sd_pre_exp * sd_post_exp)
  J_exp <- .calculate_J(n_exp - 1)
  g_exp_manual <- d_exp_manual * J_exp
  var_g_exp <- sd_change_exp^2 / (sd_pre_exp^2 * (n_exp - 1)) +
               g_exp_manual^2 / (2 * (n_exp - 1))
  var_exp_manual <- var_g_exp / (J_exp^2)

  # Manual calculation for control group
  d_nexp_manual <- (mean_post_nexp - mean_pre_nexp) / sd_pre_nexp
  sd_change_nexp <- sqrt(sd_pre_nexp^2 + sd_post_nexp^2 - 2 * r_nexp * sd_pre_nexp * sd_post_nexp)
  J_nexp <- .calculate_J(n_nexp - 1)
  g_nexp_manual <- d_nexp_manual * J_nexp
  var_g_nexp <- sd_change_nexp^2 / (sd_pre_nexp^2 * (n_nexp - 1)) +
                g_nexp_manual^2 / (2 * (n_nexp - 1))
  var_nexp_manual <- var_g_nexp / (J_nexp^2)

  # Combined effect (difference of the two SEPARATELY standardized within-group d's)
  d_combined_manual <- d_exp_manual - d_nexp_manual
  se_combined_manual <- sqrt(var_exp_manual + var_nexp_manual)

  # metaConvert calculation -- pool_sd = FALSE (the default) selects the per-arm rule
  # that the manual calculation above replicates.
  result <- es_from_means_sd_pre_post(smd_var = "hedges_olkin", 
    mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r_exp, r_nexp,
    pre_post_to_smd = "bonett",
    pool_sd = FALSE
  )

  expect_equal(result$d, d_combined_manual, tolerance = 1e-10,
               label = "Two-group d should equal difference of within-group d's")
  expect_equal(result$d_se, se_combined_manual, tolerance = 1e-10,
               label = "Two-group SE should be sqrt of sum of variances")
})

# Test 3: Cooper (2019) / Morris & DeShon (2002) Formula - Single Group ====
test_that("cooper formula matches Morris & DeShon (2002) for single-group", {
  # Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates
  # in meta-analysis with repeated measures and independent-groups designs.
  # Psychological Methods, 7(1), 105-125.
  # https://doi.org/10.1037/1082-989X.7.1.105
  #
  # Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019).
  # The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
  # (Table 12.2)

  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 12
  n <- 30
  r <- 0.6

  # Manual calculation per Morris & DeShon (2002) Equation 8
  # d_rm = (mean_diff / sd_change) × sqrt(2(1-r))
  sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
  d_rm_manual <- (mean_post - mean_pre) / sd_change * sqrt(2 * (1 - r))

  # Variance formula (Morris & DeShon 2002)
  # var(d_rm) = 2(1-r)/n + d_rm²/(2n)
  var_rm_manual <- 2 * (1 - r) / n + d_rm_manual^2 / (2 * n)
  se_rm_manual <- sqrt(var_rm_manual)

  # metaConvert calculation
  result <- es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "cooper"
  )

  expect_equal(result$d, d_rm_manual, tolerance = 1e-10,
               label = "Cooper d_rm should match Morris & DeShon (2002) Eq 8")
  expect_equal(result$d_se, se_rm_manual, tolerance = 1e-10,
               label = "Cooper SE should match Morris & DeShon variance formula")
})

test_that("cooper formula verified across r values", {
  # Verify that formula is r-dependent as specified
  mean_pre <- 45
  mean_post <- 52
  sd_pre <- 8
  sd_post <- 9
  n <- 25

  r_values <- c(0, 0.3, 0.5, 0.7, 0.9)

  for (r in r_values) {
    # Manual calculation
    sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
    d_rm_manual <- (mean_post - mean_pre) / sd_change * sqrt(2 * (1 - r))
    var_rm_manual <- 2 * (1 - r) / n + d_rm_manual^2 / (2 * n)

    # metaConvert
    result <- es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
      mean_pre, mean_post, sd_pre, sd_post, n, r,
      pre_post_to_smd = "cooper"
    )

    expect_equal(result$d, d_rm_manual, tolerance = 1e-10,
                 label = paste("Cooper d_rm correct for r =", r))
    expect_equal(result$d_se^2, var_rm_manual, tolerance = 1e-10,
                 label = paste("Cooper variance correct for r =", r))
  }
})

test_that("cooper r-dependency verification", {
  # Verify variance properties for cooper/morris_drm standardizer
  # Formula: d_rm = (mean_diff / sd_change) × sqrt(2(1-r))
  # where sd_change = sqrt(sd_pre² + sd_post² - 2r×sd_pre×sd_post)
  # The relationship between r and d_rm is complex due to both terms

  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 10
  sd_post <- 10
  n <- 30

  r_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)

  results <- lapply(r_values, function(r) {
    es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
      mean_pre, mean_post, sd_pre, sd_post, n, r,
      pre_post_to_smd = "cooper"
    )
  })

  # SE should decrease monotonically with r (from variance formula)
  se_values <- sapply(results, function(x) x$d_se)
  expect_true(all(diff(se_values) < 0),
              label = "Cooper SE should decrease as r increases")

  # All results should be finite and positive
  d_values <- sapply(results, function(x) x$d)
  expect_true(all(is.finite(d_values)),
              label = "All d values should be finite")
  expect_true(all(d_values > 0),
              label = "All d values should be positive for positive mean change")
})

# Test 4: Cooper Formula - Two Groups ====
test_that("cooper two-group matches Morris & DeShon (2002) (per-arm standardizer, the pool_sd = FALSE default)", {
  # PER-ARM PATH (pool_sd = FALSE) -- the DEFAULT. As in Test 2, the manual
  # calculation computes a d_rm PER ARM (each divided by its OWN change SD),
  # subtracts them and adds their variances. That is the pool_sd = FALSE
  # construction. The opt-in pool_sd = TRUE alternative divides the between-arm
  # change difference by a SINGLE change SD pooled across arms (covered by the
  # "pooled standardizer" tests below).
  mean_pre_exp <- 30
  mean_exp <- 45
  sd_pre_exp <- 10
  sd_post_exp <- 11
  n_exp <- 25
  r_exp <- 0.7

  mean_pre_nexp <- 32
  mean_post_nexp <- 34
  sd_pre_nexp <- 9
  sd_post_nexp <- 10
  n_nexp <- 25
  r_nexp <- 0.7

  # Manual calculation for each group
  sd_change_exp <- sqrt(sd_pre_exp^2 + sd_post_exp^2 - 2 * r_exp * sd_pre_exp * sd_post_exp)
  d_rm_exp <- (mean_exp - mean_pre_exp) / sd_change_exp * sqrt(2 * (1 - r_exp))
  var_rm_exp <- 2 * (1 - r_exp) / n_exp + d_rm_exp^2 / (2 * n_exp)

  sd_change_nexp <- sqrt(sd_pre_nexp^2 + sd_post_nexp^2 - 2 * r_nexp * sd_pre_nexp * sd_post_nexp)
  d_rm_nexp <- (mean_post_nexp - mean_pre_nexp) / sd_change_nexp * sqrt(2 * (1 - r_nexp))
  var_rm_nexp <- 2 * (1 - r_nexp) / n_nexp + d_rm_nexp^2 / (2 * n_nexp)

  # Combined (difference of the two SEPARATELY standardized within-group d_rm's)
  d_combined <- d_rm_exp - d_rm_nexp
  se_combined <- sqrt(var_rm_exp + var_rm_nexp)

  # metaConvert -- pool_sd = FALSE (the default) selects the per-arm rule replicated above.
  result <- es_from_means_sd_pre_post(smd_var = "borenstein", 
    mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r_exp, r_nexp,
    pre_post_to_smd = "cooper",
    pool_sd = FALSE
  )

  expect_equal(result$d, d_combined, tolerance = 1e-10)
  expect_equal(result$d_se, se_combined, tolerance = 1e-10)
})

# Test 4b: POOLED standardizer (OPT-IN, pool_sd = TRUE) ====
# Morris (2008) d_ppc2: divide the between-arm difference in mean change by a SINGLE
# standardizing SD pooled across arms, rather than standardizing each arm by its own SD
# and subtracting (d_ppc1, the pool_sd = FALSE DEFAULT). d_ppc2 is more efficient but
# assumes the two arms' true standardizing SDs are equal, so it must be requested
# explicitly with pool_sd = TRUE.

test_that("pool_sd defaults to FALSE (per-arm d_ppc1, Becker 1988 / Morris d_ppc1)", {
  mean_pre_exp <- 30; mean_exp <- 45; sd_pre_exp <- 10; sd_post_exp <- 11
  mean_pre_nexp <- 32; mean_post_nexp <- 34; sd_pre_nexp <- 9; sd_post_nexp <- 10
  n_exp <- 25; n_nexp <- 25; r_exp <- 0.7; r_nexp <- 0.7

  call_pp <- function(...) {
    es_from_means_sd_pre_post(smd_var = "hedges_olkin", 
      mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
      mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
      n_exp, n_nexp, r_exp, r_nexp, ...
    )
  }

  for (method in c("bonett", "cooper", "morris_drm", "morris_dz", "morris_dav")) {
    res_default <- call_pp(pre_post_to_smd = method)
    res_pooled <- call_pp(pre_post_to_smd = method, pool_sd = TRUE)
    res_perarm <- call_pp(pre_post_to_smd = method, pool_sd = FALSE)

    # The default must BE the per-arm (d_ppc1) path...
    expect_equal(res_default$d, res_perarm$d, tolerance = 1e-12,
                 label = paste("default == pool_sd=FALSE for", method))
    expect_equal(res_default$d_se, res_perarm$d_se, tolerance = 1e-12,
                 label = paste("default SE == pool_sd=FALSE SE for", method))

    # ...and must NOT silently be the pooled (d_ppc2) path. The arms' SDs differ here
    # (10/11 vs 9/10), so the two constructions are genuinely distinct.
    expect_false(isTRUE(all.equal(res_default$d, res_pooled$d, tolerance = 1e-6)),
                 label = paste("default differs from pool_sd=TRUE for", method))
  }
})

test_that("opt-in pooled standardizer (pool_sd = TRUE) divides the between-arm change difference by ONE pooled SD", {
  # Definitional check of the OPT-IN pooled estimator (plain arithmetic, not a source
  # transcription): d = (change_exp - change_nexp) / sd_pooled, where sd_pooled is
  # the df-weighted pool of the per-arm standardizing SDs. pool_sd = TRUE is passed
  # explicitly: it is NOT the default (the default is the per-arm d_ppc1 rule).
  mean_pre_exp <- 30; mean_exp <- 45; sd_pre_exp <- 10; sd_post_exp <- 11
  mean_pre_nexp <- 32; mean_post_nexp <- 34; sd_pre_nexp <- 9; sd_post_nexp <- 10
  n_exp <- 25; n_nexp <- 25; r_exp <- 0.7; r_nexp <- 0.7

  N <- n_exp + n_nexp
  m <- N - 2                     # pooled degrees of freedom
  r_avg <- (n_exp * r_exp + n_nexp * r_nexp) / N
  mean_diff <- (mean_exp - mean_pre_exp) - (mean_post_nexp - mean_pre_nexp)

  pool2 <- function(sd_e, sd_n) {
    sqrt(((n_exp - 1) * sd_e^2 + (n_nexp - 1) * sd_n^2) / m)
  }
  sd_change_exp <- sqrt(sd_pre_exp^2 + sd_post_exp^2 - 2 * r_exp * sd_pre_exp * sd_post_exp)
  sd_change_nexp <- sqrt(sd_pre_nexp^2 + sd_post_nexp^2 - 2 * r_nexp * sd_pre_nexp * sd_post_nexp)

  call_pp <- function(method) {
    es_from_means_sd_pre_post(smd_var = "hedges_olkin", 
      mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
      mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
      n_exp, n_nexp, r_exp, r_nexp, pre_post_to_smd = method,
      pool_sd = TRUE
    )
  }

  # bonett: standardizer = pooled BASELINE SD
  expect_equal(call_pp("bonett")$d,
               mean_diff / pool2(sd_pre_exp, sd_pre_nexp),
               tolerance = 1e-10,
               label = "pooled bonett divides by the pooled baseline SD")

  # morris_dz: standardizer = pooled CHANGE-SCORE SD
  expect_equal(call_pp("morris_dz")$d,
               mean_diff / pool2(sd_change_exp, sd_change_nexp),
               tolerance = 1e-10,
               label = "pooled morris_dz divides by the pooled change SD")

  # morris_drm / cooper: pooled change SD, rescaled to the raw-score metric
  expect_equal(call_pp("cooper")$d,
               mean_diff / pool2(sd_change_exp, sd_change_nexp) * sqrt(2 * (1 - r_avg)),
               tolerance = 1e-10,
               label = "pooled cooper rescales by sqrt(2(1-r))")

  # morris_dav: standardizer = pooled AVERAGE of pre/post SDs
  sd_av_exp <- sqrt((sd_pre_exp^2 + sd_post_exp^2) / 2)
  sd_av_nexp <- sqrt((sd_pre_nexp^2 + sd_post_nexp^2) / 2)
  expect_equal(call_pp("morris_dav")$d,
               mean_diff / pool2(sd_av_exp, sd_av_nexp),
               tolerance = 1e-10,
               label = "pooled morris_dav divides by the pooled average SD")
})

test_that("opt-in pooled morris_dz (pool_sd = TRUE) reproduces metafor::escalc(measure='SMD') on change scores", {
  # EXTERNAL reference (not a source transcription): with the change scores treated
  # as the outcome, the pooled morris_dz estimator IS the ordinary two-sample
  # Hedges' g, so metafor must reproduce its point estimate exactly.
  # This covers the OPT-IN pooled path -- pool_sd = TRUE is passed explicitly below,
  # because the default is the per-arm d_ppc1 rule (whose estimand is different).
  skip_if_not_installed("metafor")

  mean_pre_exp <- 30; mean_exp <- 45; sd_pre_exp <- 10; sd_post_exp <- 11
  mean_pre_nexp <- 32; mean_post_nexp <- 34; sd_pre_nexp <- 9; sd_post_nexp <- 10
  n_exp <- 25; n_nexp <- 25; r_exp <- 0.7; r_nexp <- 0.7

  sd_change_exp <- sqrt(sd_pre_exp^2 + sd_post_exp^2 - 2 * r_exp * sd_pre_exp * sd_post_exp)
  sd_change_nexp <- sqrt(sd_pre_nexp^2 + sd_post_nexp^2 - 2 * r_nexp * sd_pre_nexp * sd_post_nexp)

  ref <- metafor::escalc(
    measure = "SMD",
    m1i = mean_exp - mean_pre_exp, m2i = mean_post_nexp - mean_pre_nexp,
    sd1i = sd_change_exp, sd2i = sd_change_nexp,
    n1i = n_exp, n2i = n_nexp
  )

  result <- es_from_means_sd_pre_post(smd_var = "hedges_olkin", 
    mean_pre_exp, mean_exp, sd_pre_exp, sd_post_exp,
    mean_pre_nexp, mean_post_nexp, sd_pre_nexp, sd_post_nexp,
    n_exp, n_nexp, r_exp, r_nexp,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE                  # opt-in: NOT the default
  )

  expect_equal(result$g, as.numeric(ref$yi), tolerance = 1e-8,
               label = "pooled morris_dz g == metafor Hedges' g on change scores")

  # The dz variance must NOT carry a 2*(1-r) factor. Here 2*(1-0.7) = 0.6, so the
  # old (defective) variance was ~40% too small; requiring agreement with metafor's
  # variance to within 5% pins the fix without re-transcribing the source formula.
  # (The residual few-% gap is the known J^2 / (2m vs 2N) convention difference.)
  expect_equal(result$g_se^2 / as.numeric(ref$vi), 1, tolerance = 0.05,
               label = "pooled morris_dz variance agrees with metafor (no spurious 2(1-r))")
})

# Test 5: Morris_drm alias verification ====
test_that("morris_drm is exact alias of cooper", {
  # Verify that morris_drm and cooper produce identical results
  # They should be the same formula, just different names

  mean_pre <- 55
  mean_post <- 65
  sd_pre <- 12
  sd_post <- 14
  n <- 40
  r <- 0.75

  result_cooper <- es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "cooper"
  )

  result_morris_drm <- es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "morris_drm"
  )

  # Should be identical to machine precision
  expect_equal(result_cooper$d, result_morris_drm$d, tolerance = 1e-15,
               label = "morris_drm should be exact alias of cooper")
  expect_equal(result_cooper$d_se, result_morris_drm$d_se, tolerance = 1e-15)
  expect_equal(result_cooper$g, result_morris_drm$g, tolerance = 1e-15)
})

# Test 6: Mean Difference (MDw) Formula Verification ====
test_that("within-group mean difference formula is straightforward difference", {
  # MDw (within-group mean difference) should simply be mean_post - mean_pre
  # var(MDw) = sd_change² / n

  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 12
  n <- 30
  r <- 0.6

  # Manual calculation
  mdw_manual <- mean_post - mean_pre
  sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
  var_mdw_manual <- sd_change^2 / n
  se_mdw_manual <- sqrt(var_mdw_manual)

  # metaConvert (works for any standardizer)
  result <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  expect_equal(result$mdw, mdw_manual, tolerance = 1e-10,
               label = "MDw should equal mean difference")
  expect_equal(result$mdw_se, se_mdw_manual, tolerance = 1e-10,
               label = "MDw SE should follow standard formula")
})

# Test 7: Special Cases and Boundary Conditions ====
test_that("bonett equals cooper when sd_pre = sd_post and specific r", {
  # When sd_pre = sd_post (homogeneous variance), different formulas converge
  mean_pre <- 40
  mean_post <- 50
  sd <- 10  # Equal SD
  n <- 30
  r <- 0.5

  result_bonett <- es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
    mean_pre, mean_post, sd, sd, n, r,
    pre_post_to_smd = "bonett"
  )

  result_cooper <- es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
    mean_pre, mean_post, sd, sd, n, r,
    pre_post_to_smd = "cooper"
  )

  # When sd_pre = sd_post:
  # bonett: d = mean_diff / sd
  # cooper: d = (mean_diff / sd_change) * sqrt(2(1-r))
  # sd_change = sqrt(2*sd² - 2*r*sd²) = sd*sqrt(2(1-r))
  # So cooper = mean_diff / (sd*sqrt(2(1-r))) * sqrt(2(1-r)) = mean_diff / sd = bonett

  expect_equal(result_bonett$d, result_cooper$d, tolerance = 1e-10,
               label = "Bonett = Cooper when sd_pre = sd_post")
})

test_that("formulas handle r=0 correctly", {
  # When r=0 (no correlation), formulas should still work
  mean_pre <- 60
  mean_post <- 70
  sd_pre <- 12
  sd_post <- 14
  n <- 50
  r <- 0

  # Bonett formula
  result_bonett <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  d_bonett_manual <- (mean_post - mean_pre) / sd_pre
  expect_equal(result_bonett$d, d_bonett_manual, tolerance = 1e-10)

  # Cooper formula
  result_cooper <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "cooper"
  )

  sd_change <- sqrt(sd_pre^2 + sd_post^2)  # r=0, so no covariance term
  d_cooper_manual <- (mean_post - mean_pre) / sd_change * sqrt(2)
  expect_equal(result_cooper$d, d_cooper_manual, tolerance = 1e-10)
})

test_that("formulas handle r approaching 1", {
  # When r approaches 1, variance should approach 0
  mean_pre <- 55
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 10
  n <- 40
  r <- 0.99

  result_cooper <- es_from_means_sd_pre_post_single_group(smd_var = "borenstein", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "cooper"
  )

  # Manual calculation
  sd_change <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
  d_manual <- (mean_post - mean_pre) / sd_change * sqrt(2 * (1 - r))
  var_manual <- 2 * (1 - r) / n + d_manual^2 / (2 * n)

  expect_equal(result_cooper$d, d_manual, tolerance = 1e-10)
  expect_equal(result_cooper$d_se^2, var_manual, tolerance = 1e-10)

  # Variance should be very small
  expect_true(result_cooper$d_se < 0.1,
              label = "SE should be very small when r approaches 1")
})

# Test 8: Hedges' correction factor verification ====
test_that("Hedges g correction factor matches standard formula", {
  # J(df) = Γ(df/2) / (sqrt(df/2) × Γ((df-1)/2))
  # where Γ is the gamma function

  n <- 25
  df <- n - 1

  # Manual calculation
  J_manual <- exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))

  # Get result from metaConvert
  mean_pre <- 40
  mean_post <- 50
  sd_pre <- 10
  sd_post <- 10
  r <- 0.6

  result <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "bonett"
  )

  # g should equal d × J
  g_expected <- result$d * J_manual
  expect_equal(result$g, g_expected, tolerance = 1e-10,
               label = "Hedges g should equal d × J")

  # Variance of g should equal J² × variance of d
  var_g_expected <- result$d_se^2 * J_manual^2
  expect_equal(result$g_se^2, var_g_expected, tolerance = 1e-10,
               label = "var(g) should equal J² × var(d)")
})

test_that("Hedges correction more important for small samples", {
  # J approaches 1 as n increases, so correction matters more for small n

  mean_pre <- 50
  mean_post <- 60
  sd_pre <- 10
  sd_post <- 10
  r <- 0.6

  # Small sample
  n_small <- 10
  result_small <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n_small, r,
    pre_post_to_smd = "bonett"
  )

  J_small <- .calculate_J(n_small - 1)
  correction_small <- abs(1 - J_small)

  # Large sample
  n_large <- 200
  result_large <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n_large, r,
    pre_post_to_smd = "bonett"
  )

  J_large <- .calculate_J(n_large - 1)
  correction_large <- abs(1 - J_large)

  # Correction should be larger for small sample
  expect_true(correction_small > correction_large,
              label = "Hedges correction more important for small n")

  # For large n, J should be very close to 1
  expect_true(J_large > 0.99,
              label = "J should approach 1 for large n")
})

# Test 9: Confidence Interval Formula Verification ====
test_that("confidence intervals use correct t-distribution", {
  mean_pre <- 45
  mean_post <- 52
  sd_pre <- 8
  sd_post <- 9
  n <- 25
  r <- 0.7

  result <- es_from_means_sd_pre_post_single_group(smd_var = "hedges_olkin", 
    mean_pre, mean_post, sd_pre, sd_post, n, r,
    pre_post_to_smd = "cooper"
  )

  # Manual CI calculation
  t_crit <- qt(0.975, df = n - 1)
  ci_lo_manual <- result$d - t_crit * result$d_se
  ci_up_manual <- result$d + t_crit * result$d_se

  expect_equal(result$d_ci_lo, ci_lo_manual, tolerance = 1e-10,
               label = "Lower CI should use t-distribution")
  expect_equal(result$d_ci_up, ci_up_manual, tolerance = 1e-10,
               label = "Upper CI should use t-distribution")

  # CI should be symmetric around point estimate
  expect_equal(result$d - result$d_ci_lo,
               result$d_ci_up - result$d,
               tolerance = 1e-10,
               label = "CI should be symmetric")
})
