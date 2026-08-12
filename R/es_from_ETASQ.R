#' Convert an eta-squared value to various effect size measures
#'
#' @param etasq an eta-squared value (binary predictor, ANOVA model), defined as
#'   \eqn{SS_{effect} / SS_{total}}. For a one-way, two-group ANOVA this equals the
#'   partial eta-squared \eqn{SS_{effect} / (SS_{effect} + SS_{error})}, so the two
#'   coincide and either label is correct here (unlike the ANCOVA case handled by
#'   \code{\link{es_from_etasq_adj}()}, which requires the partial form specifically).
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation.
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_etasq a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first computes a Cohen's d (D) and Hedges' g (G)
#' from the eta squared of a binary predictor (ANOVA model).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a Cohen's d** the following formula is used (Cohen, 1988):
#' \deqn{d = 2 * \sqrt{\frac{etasq}{1 - etasq}}}
#'
#' Note that this closed form is the large-sample, **equal-groups** (\eqn{n\_exp = n\_nexp})
#' limit of the exact F-based conversion. For markedly unequal group sizes it is
#' approximate; supplying the ANOVA F (\code{\link{es_from_anova_f}()}), which
#' carries the per-group sample sizes, is more accurate.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 11. ANOVA statistics, Student's t-test, or point-bis correlation'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences. Routledge.
#'
#' @export es_from_etasq
#'
#' @md
#'
#' @examples
#' es_from_etasq(etasq = 0.28, n_exp = 20, n_nexp = 22)
es_from_etasq <- function(etasq, n_exp, n_nexp, smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_etasq) {
  if (missing(reverse_etasq)) reverse_etasq <- rep(FALSE, length(etasq))
  reverse_etasq[is.na(reverse_etasq)] <- FALSE


  d <- 2 * (sqrt(etasq / (1 - etasq)))

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse = reverse_etasq
  )

  es$info_used <- "etasq"

  return(es)
}

#' Convert an adjusted eta-squared value (i.e., from an ANCOVA) to various effect size measures
#'
#' @param etasq_adj an adjusted eta-squared value obtained from an ANCOVA model. This
#'   must be the **partial** eta-squared of the group effect,
#'   \eqn{\eta_p^2 = SS_{group} / (SS_{group} + SS_{error})} — the quantity that
#'   most software (SPSS, \code{car}, \code{effectsize}) prints for an ANCOVA term.
#'   Do **not** supply the classical eta-squared \eqn{SS_{group} / SS_{total}}: in an
#'   ANCOVA the total sum of squares also contains the covariate variance, so the
#'   classical value is smaller than the partial one and would be inverted into too
#'   small an F, biasing the effect size toward zero. (In a covariate-free one-way
#'   ANOVA the two definitions coincide; see \code{\link{es_from_etasq}()}.)
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param cov_outcome_r correlation between the outcome and covariate (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_etasq a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the adjusted (partial) eta-squared of a binary predictor
#' (ANCOVA model) into the ANCOVA F-statistic it implies, and then relies on the
#' calculations of the \code{\link{es_from_ancova_f}()} function. The returned
#' Cohen's d (D) and Hedges' g (G) are therefore expressed on the
#' **marginal (unadjusted) SD scale** — consistent with the other
#' \code{es_from_ancova_*} functions — not on the residual (covariate-adjusted)
#' SD scale on which a partial eta-squared is natively defined.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To convert the adjusted eta-squared into an ANCOVA F-statistic**,
#' the partial-eta-squared identity \eqn{\eta_p^2 = F / (F + df)} is inverted
#' (this is why a *partial* eta-squared, not a classical one, is required):
#' \deqn{df = n\_exp + n\_nexp - 2 - n\_cov\_ancova}
#' \deqn{ancova\_f = \frac{etasq\_adj * df}{1 - etasq\_adj}}
#'
#' **To estimate a Cohen's d** the formula used is (table 12.3 in Cooper):
#' \deqn{cohen\_d = \sqrt{ancova\_f * \frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp}} * \sqrt{1 - cov\_out\_cor^2}}
#'
#' Note that the back-transformation to the marginal SD scale requires
#' \code{cov_outcome_r}; when it is missing, the effect size estimates are
#' returned as NA (the marginal scale is not identified from the adjusted
#' eta-squared alone).
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d_adj}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 18. Adjusted: ANCOVA statistics, eta-squared'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @note
#' Two cautions apply. (1) \code{etasq_adj} must be the *partial* eta-squared (see
#' the parameter description); a classical eta-squared biases the effect size
#' toward zero. (2) The sampling variance is rebuilt from Cooper's eq. 12.26, which
#' sets the covariate-imbalance ("leverage") term
#' \eqn{D = (\bar{x}_{exp} - \bar{x}_{nexp})^2 / SS_x} of the exact ANCOVA variance
#' to zero (Lai & Kelley, 2012, eq. 5); \eqn{D} is not recoverable from summary
#' statistics. Because partial eta-squared is an exact function of the t statistic
#' (\eqn{\eta^2_p = t^2/(t^2 + df)}), this route inherits the behaviour of
#' \code{\link{es_from_ancova_t}()}: \strong{the point estimate is NOT unbiased under
#' covariate imbalance}. Both \code{d} and its standard error are attenuated by
#' \eqn{1/\sqrt{1 + \delta_x^2/4}} for equal arms (about 3% at
#' \eqn{\delta_x = 0.5}, 11% at \eqn{\delta_x = 1}, negligible when randomised).
#' The p-value is unaffected because the two shrink together, but the magnitude is
#' understated and the understatement does not average out on pooling. See
#' \code{\link{es_from_ancova_t}()} for the full statement, and prefer
#' \code{\link{es_from_ancova_means_sd}()} when adjusted means are available.
#'
#' @references
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence intervals. British Journal of Mathematical and Statistical Psychology, 65(2), 350-370.
#'
#' @export es_from_etasq_adj
#'
#' @md
#'
#' @examples
#' es_from_etasq_adj(etasq_adj = 0.28, n_cov_ancova = 3, cov_outcome_r = 0.2, n_exp = 20, n_nexp = 22)
es_from_etasq_adj <- function(etasq_adj, n_exp, n_nexp, n_cov_ancova, cov_outcome_r,
                              smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_etasq) {
  if (missing(reverse_etasq)) reverse_etasq <- rep(FALSE, length(etasq_adj))
  reverse_etasq[is.na(reverse_etasq)] <- FALSE

  # A partial eta-squared from an ANCOVA is natively defined on the RESIDUAL
  # (covariate-adjusted) SD scale. Convert it back to the ANCOVA F-statistic it
  # implies (eta_p^2 = F / (F + df_err), so F = eta_p^2 * df_err / (1 - eta_p^2))
  # and delegate to es_from_ancova_f(), which back-transforms the point estimate
  # to the MARGINAL (unadjusted-SD) scale via sqrt(1 - cov_outcome_r^2)
  # (table 12.3 in Cooper). This keeps a single source of truth for the ANCOVA
  # family: the point estimate, the (1 - cov_outcome_r^2)-shrunk sampling
  # variance and the df = N - 2 - n_cov_ancova CIs are mutually consistent, and
  # etasq_adj agrees exactly with es_from_ancova_f fed the algebraically
  # equivalent statistic. When cov_outcome_r is NA the marginal scale is not
  # identified and the output degrades to NA (rather than silently emitting a
  # residual-scale value).
  df_err <- n_exp + n_nexp - 2 - n_cov_ancova
  f_implied <- etasq_adj * df_err / (1 - etasq_adj)

  es <- es_from_ancova_f(
    ancova_f = f_implied, cov_outcome_r = cov_outcome_r,
    n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_ancova_f = reverse_etasq
  )

  es$info_used <- "etasq_adj"

  return(es)
}
