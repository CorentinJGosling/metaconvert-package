#' Convert a t-statistic obtained from an ANCOVA model into several effect size measures.
#'
#' @param ancova_t a t-statistic from an ANCOVA (binary predictor)
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_ancova_t a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first computes an "adjusted" Cohen's d (D), and
#' Hedges' g (G) from the t-value of an ANCOVA (binary predictor).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a Cohen's d** the formula used is (table 12.3 in Cooper):
#' \deqn{cohen\_d =  ancova\_t* \sqrt{\frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp}}\sqrt{1 - cov\_out\_cor^2}}
#'
#' **To estimate other effect size measures**,
#' Calculations of the \code{\link{es_from_cohen_d_adj}()} are applied.
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
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_ancova_t
#'
#' @md
#'
#' @examples
#' es_from_ancova_t(ancova_t = 2, cov_outcome_r = 0.2, n_cov_ancova = 3, n_exp = 20, n_nexp = 20)
es_from_ancova_t <- function(ancova_t, cov_outcome_r, n_cov_ancova, n_exp, n_nexp,
                             smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_ancova_t) {
  if (missing(reverse_ancova_t)) reverse_ancova_t <- rep(FALSE, length(ancova_t))
  reverse_ancova_t[is.na(reverse_ancova_t)] <- FALSE

  d <- ancova_t * sqrt(1/n_exp + 1/n_nexp) * sqrt(1 - cov_outcome_r^2)

  es <- .es_from_d(
    d = d, adjusted = TRUE, cov_outcome_r = cov_outcome_r,
    n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse = reverse_ancova_t
  )

  es$info_used <- "ancova_t"

  return(es)
}

#' Convert a F-statistic obtained from an ANCOVA model into several effect size measures.
#'
#' @param ancova_f a F-statistic from an ANCOVA (binary predictor)
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_ancova_f a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first computes an "adjusted" Cohen's d (D), and
#' Hedges' g (G) from the F-value of an ANCOVA (binary predictor).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a Cohen's d** the formula used is (table 12.3 in Cooper):
#' \deqn{cohen\_d =  \sqrt{ancova\_f * \frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp}} * \sqrt{1 - cov\_out\_cor^2}}
#'
#' **To estimate other effect size measures**,
#' Calculations of the \code{\link{es_from_cohen_d_adj}()} are applied.
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
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_ancova_f
#'
#' @md
#'
#' @examples
#' es_from_ancova_f(ancova_f = 4, cov_outcome_r = 0.2, n_cov_ancova = 3, n_exp = 20, n_nexp = 20)
es_from_ancova_f <- function(ancova_f, cov_outcome_r, n_cov_ancova, n_exp, n_nexp,
                             smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_ancova_f) {
  if (missing(reverse_ancova_f)) reverse_ancova_f <- rep(FALSE, length(ancova_f))
  reverse_ancova_f[is.na(reverse_ancova_f)] <- FALSE

  t <- sqrt(ancova_f)

  es <- es_from_ancova_t(
    ancova_t = t, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_ancova_t = reverse_ancova_f
  )

  es$info_used <- "ancova_f"

  return(es)
}

#' Convert a two-tailed p-value of an ANCOVA t-test into several effect size measures.
#'
#' @param ancova_t_pval a two-tailed p-value of a t-test in an ANCOVA (binary predictor)
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_ancova_t_pval a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the p-value of an ANCOVA (binary predictor) into a t value,
#' and then relies on the calculations of the \code{\link{es_from_ancova_t}()} function.
#'
#' **To convert the p-value into a t-value,** the following formula is used (table 12.3 in Cooper):
#' \deqn{df = n\_exp + n\_nexp - 2 - n\_cov\_ancova}
#' \deqn{t = | qt(ancova\_t\_pval/2, df = df) |}
#' Then, calculations of the \code{\link{es_from_ancova_t}()} are applied.
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
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
#' @export es_from_ancova_t_pval
#'
#' @md
#'
#' @examples
#' es_from_ancova_t_pval(
#'   ancova_t_pval = 0.05, cov_outcome_r = 0.2,
#'   n_cov_ancova = 3, n_exp = 20, n_nexp = 20
#' )
es_from_ancova_t_pval <- function(ancova_t_pval, cov_outcome_r, n_cov_ancova, n_exp, n_nexp,
                                  smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_ancova_t_pval) {
  if (missing(reverse_ancova_t_pval)) reverse_ancova_t_pval <- rep(FALSE, length(ancova_t_pval))
  reverse_ancova_t_pval[is.na(reverse_ancova_t_pval)] <- FALSE

  t_inv <- abs(qt(
    p = ancova_t_pval / 2,
    df = n_exp + n_nexp - 2 - n_cov_ancova,
    lower.tail = FALSE
  ))

  es <- es_from_ancova_t(
    ancova_t = t_inv, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_ancova_t = reverse_ancova_t_pval
  )

  es$info_used <- "ancova_t_pval"

  return(es)
}

#' Convert a two-tailed p-value of an ANCOVA t-test into several effect size measures.
#'
#' @param ancova_f_pval a two-tailed p-value of an F-test in an ANCOVA (binary predictor)
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_ancova_f_pval a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the p-value of an ANCOVA (binary predictor) into a t value,
#' and then relies on the calculations of the \code{\link{es_from_ancova_t}()} function.
#'
#' **To convert the p-value into a t-value,** the following formula is used (table 12.3 in Cooper):
#' \deqn{df = n\_exp + n\_nexp - 2 - n\_cov\_ancova}
#' \deqn{t = | qt(ancova\_f\_pval/2, df = df) |}
#' Then, calculations of the \code{\link{es_from_ancova_t}()} are applied.
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
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
#' @export es_from_ancova_f_pval
#'
#' @md
#'
#' @examples
#' es_from_ancova_f_pval(
#'   ancova_f_pval = 0.05, cov_outcome_r = 0.2,
#'   n_cov_ancova = 3, n_exp = 20, n_nexp = 20
#' )
es_from_ancova_f_pval <- function(ancova_f_pval, cov_outcome_r, n_cov_ancova, n_exp, n_nexp,
                                  smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_ancova_f_pval) {
  if (missing(reverse_ancova_f_pval)) reverse_ancova_f_pval <- rep(FALSE, length(ancova_f_pval))
  reverse_ancova_f_pval[is.na(reverse_ancova_f_pval)] <- FALSE

  t_inv <- abs(qt(
    p = ancova_f_pval / 2,
    df = n_exp + n_nexp - 2 - n_cov_ancova,
    lower.tail = FALSE
  ))

  es <- es_from_ancova_t(
    ancova_t = t_inv, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_ancova_t = reverse_ancova_f_pval
  )

  es$info_used <- "ancova_f_pval"

  return(es)
}
