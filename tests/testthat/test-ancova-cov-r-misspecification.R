# =============================================================================
# Roadmap item 2.1 -- pin HOW cov_outcome_r enters each ANCOVA-family route.
#
# The @param documentation makes a factual claim about what happens when a user
# mis-specifies cov_outcome_r (supplying the total-sample correlation, or the
# square root of the whole model R-squared, instead of the pooled within-group
# correlation). That claim splits the family in two, and the split was documented
# wrongly:
#
#   RESIDUAL-SD routes (14): cov_outcome_r rescales the residual SD to the
#     marginal scale AND enters the sampling variance. Estimate and SE move by the
#     same factor, so the t statistic -- and the p-value -- are unchanged. A wrong
#     value is genuinely invisible. This is what the docs say.
#
#   MARGINAL-SD routes (2): es_from_ancova_means_sd_pooled_crude and
#     es_from_cohen_d_adj receive an ALREADY-marginal SD, so cov_outcome_r enters
#     the sampling variance ONLY, via the (1 - R^2) factor of Cooper eq. 12.26.
#     Nothing cancels: the estimate is untouched, the SE shrinks, the p-value
#     MOVES, and the error is visible. es_from_ancova_means_sd_pooled_crude
#     nevertheless carried the residual-SD wording verbatim -- the exact inverse of
#     its own behaviour -- telling users a maximally-detectable error was
#     undetectable.
#
# These tests assert the behaviour, so the documentation cannot silently drift
# from it again. They are behavioural, not textual: they would fail if a future
# change altered which routes rescale the standardizer.
# =============================================================================

# A single ANCOVA scenario expressed in every input form the family accepts, so
# each route sees the same underlying study.
.anc <- local({
  n1 <- 60; n2 <- 60; N <- n1 + n2
  s_res <- 3; md <- 2
  se_md <- s_res * sqrt(1 / n1 + 1 / n2)
  tstat <- md / se_md
  list(
    n1 = n1, n2 = n2, s_res = s_res, md = md, se_md = se_md, tstat = tstat,
    fstat = tstat^2,
    pval = 2 * stats::pt(-abs(tstat), N - 3),
    eta = tstat^2 / (tstat^2 + (N - 3)),
    ci_lo = md - stats::qt(.975, N - 3) * se_md,
    ci_up = md + stats::qt(.975, N - 3) * se_md,
    se_arm = s_res / sqrt(n1)
  )
})

# Every route, evaluated at one value of cov_outcome_r.
.anc_routes <- function(r) {
  a <- .anc
  list(
    ancova_means_sd = es_from_ancova_means_sd(
      ancova_mean_exp = 12, ancova_mean_nexp = 10,
      ancova_mean_sd_exp = a$s_res, ancova_mean_sd_nexp = a$s_res,
      cov_outcome_r = r, n_cov_ancova = 1, n_exp = a$n1, n_nexp = a$n2),
    ancova_means_se = es_from_ancova_means_se(
      ancova_mean_exp = 12, ancova_mean_nexp = 10,
      ancova_mean_se_exp = a$se_arm, ancova_mean_se_nexp = a$se_arm,
      cov_outcome_r = r, n_cov_ancova = 1, n_exp = a$n1, n_nexp = a$n2),
    ancova_md_sd = es_from_ancova_md_sd(
      ancova_md = a$md, ancova_md_sd = a$s_res,
      cov_outcome_r = r, n_cov_ancova = 1, n_exp = a$n1, n_nexp = a$n2),
    ancova_md_se = es_from_ancova_md_se(
      ancova_md = a$md, ancova_md_se = a$se_md,
      cov_outcome_r = r, n_cov_ancova = 1, n_exp = a$n1, n_nexp = a$n2),
    ancova_md_ci = es_from_ancova_md_ci(
      ancova_md = a$md, ancova_md_ci_lo = a$ci_lo, ancova_md_ci_up = a$ci_up,
      cov_outcome_r = r, n_cov_ancova = 1, n_exp = a$n1, n_nexp = a$n2),
    ancova_md_pval = es_from_ancova_md_pval(
      ancova_md = a$md, ancova_md_pval = a$pval,
      cov_outcome_r = r, n_cov_ancova = 1, n_exp = a$n1, n_nexp = a$n2),
    ancova_t = es_from_ancova_t(
      ancova_t = a$tstat, cov_outcome_r = r, n_cov_ancova = 1,
      n_exp = a$n1, n_nexp = a$n2),
    ancova_t_pval = es_from_ancova_t_pval(
      ancova_t_pval = a$pval, cov_outcome_r = r, n_cov_ancova = 1,
      n_exp = a$n1, n_nexp = a$n2),
    ancova_f = es_from_ancova_f(
      ancova_f = a$fstat, cov_outcome_r = r, n_cov_ancova = 1,
      n_exp = a$n1, n_nexp = a$n2),
    ancova_f_pval = es_from_ancova_f_pval(
      ancova_f_pval = a$pval, cov_outcome_r = r, n_cov_ancova = 1,
      n_exp = a$n1, n_nexp = a$n2),
    etasq_adj = es_from_etasq_adj(
      etasq_adj = a$eta, cov_outcome_r = r, n_cov_ancova = 1,
      n_exp = a$n1, n_nexp = a$n2),
    pooled_adj = es_from_ancova_means_sd_pooled_adj(
      ancova_mean_exp = 12, ancova_mean_nexp = 10, ancova_mean_sd_pooled = a$s_res,
      cov_outcome_r = r, n_cov_ancova = 1, n_exp = a$n1, n_nexp = a$n2),
    # --- the two marginal-SD routes ---
    pooled_crude = es_from_ancova_means_sd_pooled_crude(
      ancova_mean_exp = 12, ancova_mean_nexp = 10, mean_sd_pooled = a$s_res,
      cov_outcome_r = r, n_cov_ancova = 1, n_exp = a$n1, n_nexp = a$n2),
    cohen_d_adj = es_from_cohen_d_adj(
      cohen_d_adj = a$md / a$s_res, cov_outcome_r = r, n_cov_ancova = 1,
      n_exp = a$n1, n_nexp = a$n2)
  )
}

.MARGINAL_SD_ROUTES <- c("pooled_crude", "cohen_d_adj")

test_that("on residual-SD routes a wrong cov_outcome_r cancels: the t statistic is unchanged", {
  lo <- .anc_routes(0.0)
  hi <- .anc_routes(0.9)
  residual <- setdiff(names(lo), .MARGINAL_SD_ROUTES)

  for (nm in residual) {
    d_ratio  <- hi[[nm]]$d / lo[[nm]]$d
    se_ratio <- hi[[nm]]$d_se / lo[[nm]]$d_se
    # Both move...
    expect_false(isTRUE(all.equal(d_ratio, 1)),
                 info = paste0(nm, ": the effect size should depend on cov_outcome_r"))
    # ...by the same factor, so the t ratio is 1 and the p-value does not move.
    expect_equal(d_ratio / se_ratio, 1, tolerance = 1e-8,
                 info = paste0(nm, ": estimate and SE must move by the SAME factor, ",
                               "which is what makes the error undetectable here."))
  }
})

test_that("on the two marginal-SD routes the estimate is INVARIANT and only the SE moves", {
  # This is the fact the pooled_crude documentation used to deny.
  lo <- .anc_routes(0.0)
  hi <- .anc_routes(0.9)

  for (nm in .MARGINAL_SD_ROUTES) {
    expect_equal(hi[[nm]]$d, lo[[nm]]$d, tolerance = 1e-12,
                 info = paste0(nm, ": the effect size must NOT depend on cov_outcome_r ",
                               "-- the supplied SD is already marginal."))
    expect_lt(hi[[nm]]$d_se, lo[[nm]]$d_se)
  }
})

test_that("on the two marginal-SD routes the p-value DOES move (the error is visible)", {
  lo <- .anc_routes(0.0)
  hi <- .anc_routes(0.9)

  for (nm in .MARGINAL_SD_ROUTES) {
    t_ratio <- (hi[[nm]]$d / hi[[nm]]$d_se) / (lo[[nm]]$d / lo[[nm]]$d_se)
    expect_gt(t_ratio, 2)   # measured 2.07 at this scenario
    # Documented consequence: a ~4x inflation of the inverse-variance weight.
    w_ratio <- (lo[[nm]]$d_se / hi[[nm]]$d_se)^2
    expect_gt(w_ratio, 4)
  }
})

test_that("the SE deflation on marginal-SD routes tracks Cooper eq. 12.26's (1 - R^2)", {
  # The leading variance term carries (1 - R^2); the d^2/(2N) term does not, so the
  # realised ratio sits slightly above sqrt(1 - R^2) rather than equalling it.
  for (r in c(0.3, 0.5, 0.9)) {
    lo <- .anc_routes(0)[["pooled_crude"]]
    hi <- .anc_routes(r)[["pooled_crude"]]
    ratio <- hi$d_se / lo$d_se
    expect_gt(ratio, sqrt(1 - r^2))
    expect_lt(ratio, 1)
  }
})

test_that("exactly two exported routes have the marginal-SD structure", {
  # Guards the documentation's "the only other route with this structure" claim.
  lo <- .anc_routes(0.0)
  hi <- .anc_routes(0.9)
  invariant <- names(lo)[vapply(names(lo),
    function(nm) isTRUE(all.equal(hi[[nm]]$d, lo[[nm]]$d, tolerance = 1e-12)),
    logical(1))]
  expect_setequal(invariant, .MARGINAL_SD_ROUTES)
})

test_that("V22 does not fire on a WRONG cov_outcome_r (documented gap)", {
  # The @param text tells users V22 covers 'missing', not 'wrong'. Pin that, so the
  # claim stays true (or the test fails when someone extends V22).
  d <- data.frame(
    study_id = "s1", n_exp = 60, n_nexp = 60,
    ancova_mean_exp = 12, ancova_mean_nexp = 10, mean_sd_pooled = 3,
    cov_outcome_r = 0.9, n_cov_ancova = 1
  )
  res <- suppressMessages(convert_df(d, measure = "g", verbose = FALSE))
  s <- suppressMessages(summary(res, flags = TRUE))
  fl <- unlist(s[, grep("^flags", names(s)), drop = FALSE])
  expect_false(any(grepl("residual SD without", fl, fixed = TRUE)),
               info = "V22 fired on a present-but-wrong cov_outcome_r; update the docs.")
})
