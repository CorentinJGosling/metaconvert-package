#' Convert an eta-squared value to various effect size measures
#'
#' @param etasq an eta-squared value (binary predictor, ANOVA model))
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation.
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
es_from_etasq <- function(etasq, n_exp, n_nexp, smd_to_cor = "viechtbauer", reverse_etasq) {
  if (missing(reverse_etasq)) reverse_etasq <- rep(FALSE, length(etasq))
  reverse_etasq[is.na(reverse_etasq)] <- FALSE


  d <- 2 * (sqrt(etasq / (1 - etasq)))

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_etasq
  )

  es$info_used <- "etasq"

  return(es)
}

#' Convert an adjusted eta-squared value (i.e., from an ANCOVA) to various effect size measures
#'
#' @param etasq_adj an adjusted eta-squared value (i.e., obtained from an ANCOVA model)
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param cov_outcome_r correlation between the outcome and covariate (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_etasq a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first computes an adjusted Cohen's d (D) and Hedges' g (G)
#' from the adjusted eta squared of a binary predictor (ANCOVA model).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a Cohen's d** the following formula is used (Cohen, 1988):
#' \deqn{d\_adj = 2 * \sqrt{\frac{etasq\_adj}{1 - etasq\_adj}}}
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
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences. Routledge.
#'
#' @export es_from_etasq_adj
#'
#' @md
#'
#' @examples
#' es_from_etasq_adj(etasq = 0.28, n_cov_ancova = 3, cov_outcome_r = 0.2, n_exp = 20, n_nexp = 22)
es_from_etasq_adj <- function(etasq_adj, n_exp, n_nexp, n_cov_ancova, cov_outcome_r,
                              smd_to_cor = "viechtbauer", reverse_etasq) {
  if (missing(reverse_etasq)) reverse_etasq <- rep(FALSE, length(etasq_adj))
  reverse_etasq[is.na(reverse_etasq)] <- FALSE

  d <- 2 * (sqrt(etasq_adj / (1 - etasq_adj)))

  es <- .es_from_d(
    d = d, n_cov_ancova = n_cov_ancova, cov_outcome_r = cov_outcome_r,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_etasq
  )

  es$info_used <- "etasq_adj"

  return(es)
}
