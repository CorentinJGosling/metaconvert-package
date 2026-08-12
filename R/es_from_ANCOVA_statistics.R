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
#' @note
#' The reported ANCOVA t already embeds the exact standard error
#' \eqn{s_{res}\sqrt{1/n_{exp} + 1/n_{nexp} + D}} with the covariate-imbalance
#' ("leverage") term \eqn{D = (\bar{x}_{exp} - \bar{x}_{nexp})^2 / SS_x}
#' (Lai & Kelley, 2012, eq. 5). The back-transformation used here rebuilds the
#' variance from Cooper's eq. 12.26, which sets \eqn{D = 0}: it assumes a covariate
#' balanced across groups and treats \code{cov_outcome_r} as known. \eqn{D} cannot
#' be recovered from summary statistics, since it needs the covariate group means
#' and the within-group sum of squares.
#'
#' \strong{Unlike the adjusted-means routes, the point estimate is therefore NOT
#' unbiased under covariate imbalance.} Because \code{d} is recovered by inverting
#' that same \eqn{D = 0} expression, both \code{d} and its standard error are
#' attenuated by the identical factor
#' \eqn{1/\sqrt{1 + D/(1/n_{exp} + 1/n_{nexp})}}, which for equal arms is
#' \eqn{1/\sqrt{1 + \delta_x^2/4}} where \eqn{\delta_x} is the standardised
#' covariate imbalance: about 3% at \eqn{\delta_x = 0.5} and 11% at
#' \eqn{\delta_x = 1}, and negligible in randomised designs.
#'
#' Because the two shrink together, the implied test statistic, the p-value and
#' whether the confidence interval excludes zero are all unaffected. What is
#' affected is the \emph{magnitude}, which is understated and does not average out
#' on pooling (fixed-effect pooled bias about -11% at \eqn{\delta_x = 1}). Prefer
#' \code{\link{es_from_ancova_means_sd}()} or \code{\link{es_from_ancova_md_sd}()}
#' when the study reports adjusted means, or an adjusted mean difference together
#' with the residual SD.
#'
#' The standard error is a separate matter, and switching route does \emph{not}
#' fix it: because Cooper's eq. 12.26 omits \eqn{D} on \emph{every} ANCOVA route,
#' all of them understate the exact ANCOVA sampling variance by the same factor
#' (about 26% too much inverse-variance weight at \eqn{\delta_x = 1}, on the
#' adjusted-means routes as much as here). Relative to one another the routes
#' differ by under 1.5% in weight, so the point-estimate attenuation above passes
#' through to the pooled estimate essentially undamped.
#'
#' The size of the attenuation depends on how much covariate leverage the supplied
#' input carries, which splits the affected routes into two tiers. Routes derived
#' from the \emph{combined} mean-difference standard error --- \code{ancova_t},
#' \code{ancova_f}, \code{ancova_t_pval}, \code{etasq_adj}, \code{ancova_md_se},
#' \code{ancova_md_ci}, \code{ancova_md_pval} --- carry the full between-arm
#' \eqn{D} and take the factor above. \code{\link{es_from_ancova_means_se}()} sees
#' only each arm's own leverage \eqn{(\bar{x}_j - \bar{x})^2 / SS_x} and is
#' therefore attenuated \emph{less}: at \eqn{\delta_x = 1} with n = 60/60 the
#' measured ratios are 0.942 and 0.893 respectively. Both invert with
#' \eqn{SD = SE\sqrt{n}} (i.e. leverage = 0), so each is deflated by exactly the
#' leverage its input contained; for equal arms that gives
#' \eqn{1/\sqrt{1 + \delta_x^2/8}}, exactly half the combined-SE routes'
#' \eqn{D} contribution. At exact covariate balance every ANCOVA route agrees to
#' ~4e-16.
#'
#' Note this is the standard aggregate-data conversion, shared with other
#' implementations of Cooper's table 12.3; only the unbiasedness of the point
#' estimate is at issue.
#'
#' @references
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence intervals. British Journal of Mathematical and Statistical Psychology, 65(2), 350-370.
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
#' @note
#' Like \code{\link{es_from_ancova_t}()}, this route recovers \code{d} by inverting an
#' expression that sets the covariate-imbalance ("leverage") term \eqn{D} of the exact
#' ANCOVA variance to zero, so \strong{the point estimate is not unbiased under
#' covariate imbalance}: \code{d} and its standard error are attenuated by the same
#' factor, \eqn{1/\sqrt{1 + \delta_x^2/4}} for equal arms (about 3% at
#' \eqn{\delta_x = 0.5}, 11% at \eqn{\delta_x = 1}, negligible when randomised). The
#' p-value is unaffected. See \code{\link{es_from_ancova_t}()} for the full statement.
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
#' @note
#' Like \code{\link{es_from_ancova_t}()}, this route recovers \code{d} by inverting an
#' expression that sets the covariate-imbalance ("leverage") term \eqn{D} of the exact
#' ANCOVA variance to zero, so \strong{the point estimate is not unbiased under
#' covariate imbalance}: \code{d} and its standard error are attenuated by the same
#' factor, \eqn{1/\sqrt{1 + \delta_x^2/4}} for equal arms (about 3% at
#' \eqn{\delta_x = 0.5}, 11% at \eqn{\delta_x = 1}, negligible when randomised). The
#' p-value is unaffected. See \code{\link{es_from_ancova_t}()} for the full statement.
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
#' @note
#' Like \code{\link{es_from_ancova_t}()}, this route recovers \code{d} by inverting an
#' expression that sets the covariate-imbalance ("leverage") term \eqn{D} of the exact
#' ANCOVA variance to zero, so \strong{the point estimate is not unbiased under
#' covariate imbalance}: \code{d} and its standard error are attenuated by the same
#' factor, \eqn{1/\sqrt{1 + \delta_x^2/4}} for equal arms (about 3% at
#' \eqn{\delta_x = 0.5}, 11% at \eqn{\delta_x = 1}, negligible when randomised). The
#' p-value is unaffected. See \code{\link{es_from_ancova_t}()} for the full statement.
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
