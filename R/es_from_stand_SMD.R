#' Convert a Cohen's d value to several effect size measures
#'
#' @param cohen_d Cohen's d (i.e., standardized mean difference) value.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_d a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function estimates the standard error of a Cohen's d value and computes a Hedges' g (G).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the standard error of Cohen's d**, the following formula is used (formula 12.13 in Cooper):
#' \deqn{cohen\_d\_se =  \sqrt{\frac{n\_exp+n\_nexp}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2*(n\_exp+n\_nexp)}}}
#' \deqn{cohen\_d\_ci\_lo = cohen\_d - cohen\_d\_se * qt(.975, df = n\_exp+n\_nexp-2)}
#' \deqn{cohen\_d\_ci\_up = cohen\_d + cohen\_d\_se * qt(.975, df = n\_exp+n\_nexp-2)}
#'
#' **To estimate the Hedges' g and its standard error**, the following formulas are used (Hedges, 1981):
#' \deqn{df = n\_exp + n\_nexp - 2}
#' \deqn{J = exp(\log_{gamma}(\frac{df}{2}) - 0.5 * \log(\frac{df}{2}) - \log_{gamma}(\frac{df - 1}{2}))}
#' \deqn{hedges\_g = cohen\_d * J}
#' \deqn{hedges\_g\_se = \sqrt{cohen\_d\_se^2 * J^2}}
#' \deqn{hedges\_g\_ci\_lo = hedges\_g - hedges\_g\_se * qt(.975, df = n\_exp+n\_nexp-2)}
#' \deqn{hedges\_g\_ci\_up = hedges\_g + hedges\_g\_se * qt(.975, df = n\_exp+n\_nexp-2)}
#'
#' **To estimate the log odds ratio and its standard error**, the following formulas are used (formulas 12.34-12.35 in Cooper):
#' \deqn{logor = \frac{cohen\_d * \pi}{\sqrt{3}}}
#' \deqn{logor\_se = \sqrt{\frac{cohen\_d\_se^2 * \pi^2}{3}}}
#' \deqn{logor\_lo = logor - logor\_se * qnorm(.975)}
#' \deqn{logor\_up = logor + logor\_se * qnorm(.975)}
#' Note that this conversion assumes that responses within the two groups follow logistic distributions.
#'
#' **To estimate the correlation coefficient and its standard error**, various formulas can be used.
#'
#' **A.** To estimate the 'biserial' correlation (\code{smd_to_cor="viechtbauer"}), the following formulas are used (formulas 5, 8, 13, 17, 18, 19 in Viechtbauer):
#' \deqn{h = \frac{n\_exp + n\_nexp}{n\_exp} + \frac{n\_exp + n\_nexp}{n\_nexp}}
#' \deqn{r.pb = \frac{cohen\_d}{\sqrt{cohen\_d^2 + h}}}
#' \deqn{p = \frac{n\_exp}{n\_exp + n\_nexp}}
#' \deqn{q = 1 - p}
#' \deqn{R = \frac{\sqrt{p*q}}{dnorm(qnorm(1-p)) * r.pb}}
#' \deqn{R\_var = \frac{1}{n\_exp + n\_nexp - 1} * (\frac{\sqrt{p*q}}{dnorm(qnorm(1-p))} - R^2)^2}
#' \deqn{R\_se = \sqrt{R\_var}}
#' \deqn{a = \frac{\sqrt{dnorm(qnorm(1-p))}}{(p*q)^\frac{1}{4}}}
#' \deqn{Z = \frac{a}{2} * \log(\frac{1+a*R}{1-a*R})}
#' \deqn{Z\_var = \frac{1}{n - 1}}
#' \deqn{Z\_se = \sqrt{Z\_var}}
#' \deqn{Z\_ci\_lo = Z - qnorm(.975) * Z\_se}
#' \deqn{Z\_ci\_up = Z + qnorm(.975) * Z\_se}
#' \deqn{R\_ci\_lo = tanh(Z\_lo)}
#' \deqn{R\_ci\_up = tanh(Z\_up)}
#'
#' **B.** To estimate the correlation coefficient according to Cooper et al. (2019) (formulas 12.40-42)
#' and Borenstein et al. (2009) (formulas 54-56),
#' the following formulas are used (\code{smd_to_cor="lipsey_cooper"}):
#' \deqn{p = \frac{n\_exp}{n\_exp + n\_nexp}}
#' \deqn{R = \frac{cohen\_d}{\sqrt{cohen\_d^2 + 1 / (p * (1 - p))}}}
#' \deqn{a = \frac{(n\_exp + n\_nexp)^2}{(n\_exp*n\_nexp)}}
#' \deqn{var\_R = \frac{a^2 * cohen\_d\_se^2}{(cohen\_d^2 + a)^3}}
#' \deqn{R\_se = \sqrt{R\_var}}
#' \deqn{R\_ci\_lo = R - qt(.975, n\_exp+n\_nexp- 2) * R\_se}
#' \deqn{R\_ci\_up = R + qt(.975, n\_exp+n\_nexp- 2) * R\_se}
#' \deqn{Z = atanh(R)}
#' \deqn{Z\_var = \frac{cohen\_d\_se^2}{cohen\_d\_se^2 + (1 / p*(1-p))}}
#' \deqn{Z\_se = \sqrt{Z\_var}}
#' \deqn{Z\_ci\_lo = Z - qnorm(.975) * Z\_se}
#' \deqn{Z\_ci\_up = Z + qnorm(.975) * Z\_se}
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2021). Introduction to meta-analysis. John Wiley & Sons.
#'
#' Hedges LV (1981): Distribution theory for Glass's estimator of effect size and related estimators. Journal of Educational and Behavioral Statistics, 6, 107-28
#'
#' Jacobs, P., & Viechtbauer, W. (2017). Estimation of the biserial correlation and its sampling variance for use in meta-analysis. Research synthesis methods, 8(2), 161-180.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 1. Cohen's d or Hedges' g'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_cohen_d
#'
#' @md
#'
#' @examples
#' es_from_cohen_d(cohen_d = 1, n_exp = 20, n_nexp = 20)
es_from_cohen_d <- function(cohen_d, n_exp, n_nexp, smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_d) {
  if (missing(reverse_d)) reverse_d <- rep(FALSE, length(n_exp))
  reverse_d[is.na(reverse_d)] <- FALSE


  es <- .es_from_d(
    d = cohen_d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse = reverse_d
  )

  es$info_used <- "cohen_d"
  return(es)
}

#' Convert an adjusted Cohen's d value to several effect size measures
#'
#' @param cohen_d_adj Adjusted Cohen's d (i.e., standardized mean difference) value.
#' @param n_cov_ancova number of covariates
#' @param cov_outcome_r covariate-outcome correlation (in case of multiple covariates, the multiple correlation)
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_d a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function estimates the standard error of an adjusted Cohen's d value and Hedges' g (G), and
#' converts an odds ratio (OR) and correlation coefficients (R/Z).
#'
#' **To estimate the standard error of Cohen's d**, the following formula is used (table 12.3 in Cooper):
#' \deqn{d\_se = \sqrt{\frac{n\_exp+n\_nexp}{n\_exp*n\_nexp} * (1 - cov\_outcome\_r^2) + \frac{cohen\_d\_adj^2}{2*(n\_exp+n\_nexp)}}}
#'
#' **To estimate other effect size measures**, calculations of the
#'  \code{\link{es_from_cohen_d}()} function are used (with the exception of the degree of freedom
#'  that is estimated as \code{df = n_exp + n_nexp - 2 - n_cov_ancova}).
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
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 1. Cohen's d or Hedges' g'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_cohen_d_adj
#'
#' @md
#'
#' @examples
#' es_from_cohen_d_adj(cohen_d_adj = 1, n_cov_ancova = 4, cov_outcome_r = .30, n_exp = 20, n_nexp = 20)
es_from_cohen_d_adj <- function(cohen_d_adj, n_cov_ancova, cov_outcome_r, n_exp, n_nexp,
                                smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_d) {
  if (missing(reverse_d)) reverse_d <- rep(FALSE, length(n_exp))
  reverse_d[is.na(reverse_d)] <- FALSE

  es <- .es_from_d(
    d = cohen_d_adj, adjusted = TRUE, n_cov_ancova = n_cov_ancova,
    cov_outcome_r = cov_outcome_r,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse = reverse_d
  )

  es$info_used <- "cohen_d_adj"
  return(es)
}

#' Convert a Hedges' g value to other effect size measures (G, OR, COR)
#'
#' @param hedges_g Hedges' g value
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{hedges_g} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_g a logical value indicating whether the direction of the \code{hedges_g} value should be flipped.
#'
#' @details
#' This function estimates the standard error of the Hedges' g and the Cohen's d (D).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate standard error of Hedges'g**, the following formula is used (Hedges, 1981):
#' \deqn{df = n\_exp + n\_nexp - 2}
#' \deqn{hedges\_g\_se = \sqrt{cohen\_d\_se^2 * J^2}}
#' \deqn{hedges\_g\_ci\_lo = hedges\_g - hedges\_g\_se * qt(.975, df = n\_exp+n\_nexp-2)}
#' \deqn{hedges\_g\_ci\_up = hedges\_g + hedges\_g\_se * qt(.975, df = n\_exp+n\_nexp-2)}
#'
#' **To estimate the Cohen's d value**, the following formula is used (Hedges, 1981):
#' \deqn{J = exp(\log_{gamma}(\frac{df}{2}) - 0.5 * \log(\frac{df}{2}) - \log_{gamma}(\frac{df - 1}{2}))}
#' \deqn{cohen\_d = \frac{hedges\_g}{J}}
#' \deqn{cohen\_d\_se = \sqrt{(\frac{n\_exp+n\_nexp}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2*(n\_exp+n\_nexp)})}}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @references
#' Hedges LV (1981): Distribution theory for Glass's estimator of effect size and related estimators. Journal of Educational and Behavioral Statistics, 6, 107-28
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 1. Cohen's d or Hedges' g'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_hedges_g
#'
#' @md
#'
#' @examples
#' es_from_hedges_g(hedges_g = 0.243, n_exp = 20, n_nexp = 20)
es_from_hedges_g <- function(hedges_g, n_exp, n_nexp, smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_g) {
  if (missing(reverse_g)) reverse_g <- rep(FALSE, length(n_exp))
  reverse_g[is.na(reverse_g)] <- FALSE
  if (length(reverse_g) == 1) reverse_g = c(rep(reverse_g, length(hedges_g)))
  if (length(reverse_g) != length(hedges_g)) stop("The length of the 'reverse_g' argument is incorrectly specified.")

  df <- n_exp + n_nexp - 2

  hedges_g <- ifelse(reverse_g, -hedges_g, hedges_g)

  J <- .d_j(df)
  d <- hedges_g / J

  es <- .es_from_d(d = d, n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor, smd_var = smd_var)

  es$info_used <- "hedges_g"
  return(es)
}
