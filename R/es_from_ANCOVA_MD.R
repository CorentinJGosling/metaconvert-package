#' Convert an adjusted mean difference and adjusted standard deviation between two independent groups obtained from an ANCOVA model into several effect size measures
#'
#' @param ancova_md adjusted mean difference between two independent groups
#' @param ancova_md_sd covariate-adjusted pooled within-group standard deviation
#' @param cov_outcome_r correlation between the outcome and covariate (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_ancova_md a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first computes an "adjusted" Cohen's d (D), Hedges' g (G)
#' from the adjusted mean difference (MD).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the SE of the adjusted MD**:
#' \deqn{md\_se = ancova\_md\_sd * \sqrt{\frac{1}{n\_exp} + \frac{1}{n\_nexp}}}
#' \deqn{md\_lo = md - md\_se * qt(.975, n\_exp + n\_nexp-2-n\_cov\_ancova)}
#' \deqn{md\_up = md + md\_se * qt(.975, n\_exp + n\_nexp-2-n\_cov\_ancova)}
#'
#' **To recover the unadjusted pooled SD** (for Cohen's d, table 12.3 in Cooper):
#' \deqn{md\_sd = \frac{ancova\_md\_sd}{\sqrt{1 - cor\_outcome\_r^2}}}
#'
#' **To estimate the Cohen's d** (table 12.3 in Cooper):
#' \deqn{d = \frac{ancova\_md}{md\_sd}}
#'
#' **To estimate other effect size measures**,
#' Calculations of the \code{\link{es_from_cohen_d_adj}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 20. Adjusted: Mean difference and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_ancova_md_sd
#'
#' @md
#'
#' @examples
#' es_from_ancova_md_sd(
#'   ancova_md = 4, ancova_md_sd = 2,
#'   cov_outcome_r = 0.5, n_cov_ancova = 5,
#'   n_exp = 20, n_nexp = 22
#' )
es_from_ancova_md_sd <- function(ancova_md, ancova_md_sd,
                                 cov_outcome_r, n_cov_ancova,
                                 n_exp, n_nexp,
                                 smd_to_cor = "viechtbauer", smd_var = "borenstein",
                                 reverse_ancova_md) {
  if (missing(reverse_ancova_md)) reverse_ancova_md <- rep(FALSE, length(ancova_md))
  reverse_ancova_md[is.na(reverse_ancova_md)] <- FALSE
  if (length(reverse_ancova_md) == 1) reverse_ancova_md = c(rep(reverse_ancova_md, length(ancova_md)))
  if (length(reverse_ancova_md) != length(ancova_md)) stop("The length of the 'reverse_ancova_md' argument is incorrectly specified.")

  md_sd <- ancova_md_sd / sqrt(1 - cov_outcome_r^2)

  d <- ancova_md / md_sd

  es <- .es_from_d(
    d = d, cov_outcome_r = cov_outcome_r,
    adjusted = TRUE,
    n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse = reverse_ancova_md
  )

  es$info_used <- "ancova_md_sd"

  es$md <- ifelse(reverse_ancova_md, -ancova_md, ancova_md)
  es$md_se <- ancova_md_sd * sqrt(1 / n_exp + 1 / n_nexp)
  es$md_ci_lo <- es$md - qt(.975, n_exp + n_nexp - 2 - n_cov_ancova) * es$md_se
  es$md_ci_up <- es$md + qt(.975, n_exp + n_nexp - 2 - n_cov_ancova) * es$md_se

  return(es)
}


#' Convert an adjusted mean difference and standard error between two independent groups obtained from an ANCOVA model into several effect size measures
#'
#' @param ancova_md adjusted mean difference between two independent groups
#' @param ancova_md_se covariate-adjusted standard error of the mean difference
#' @param cov_outcome_r correlation between the outcome and covariate (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_ancova_md a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean difference (MD) standard error into a standard deviation,
#' and then relies on the calculations of the \code{\link{es_from_ancova_md_sd}} function.
#'
#' **To convert the standard error into a standard deviation**, the following formula is used.
#' \deqn{ancova\_md\_sd = \frac{ancova\_md\_se}{\sqrt{1 / n_exp + 1 / n_nexp}}}
#' Calculations of the \code{\link{es_from_ancova_md_sd}()} are then applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 20. Adjusted: Mean difference and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_ancova_md_se
#'
#' @md
#'
#' @examples
#' es_from_ancova_md_se(
#'   ancova_md = 4, ancova_md_se = 2,
#'   cov_outcome_r = 0.5, n_cov_ancova = 5,
#'   n_exp = 20, n_nexp = 22
#' )
es_from_ancova_md_se <- function(ancova_md, ancova_md_se,
                                 cov_outcome_r, n_cov_ancova,
                                 n_exp, n_nexp,
                                 smd_to_cor = "viechtbauer", smd_var = "borenstein",
                                 reverse_ancova_md) {
  if (missing(reverse_ancova_md)) reverse_ancova_md <- rep(FALSE, length(ancova_md))
  reverse_ancova_md[is.na(reverse_ancova_md)] <- FALSE
  if (length(reverse_ancova_md) == 1) reverse_ancova_md = c(rep(reverse_ancova_md, length(ancova_md)))
  if (length(reverse_ancova_md) != length(ancova_md)) stop("The length of the 'reverse_ancova_md' argument is incorrectly specified.")

  ancova_md_sd <- ancova_md_se / sqrt(1/n_exp + 1/n_nexp)

  es <- es_from_ancova_md_sd(ancova_md = ancova_md,
                             ancova_md_sd = ancova_md_sd,
                             cov_outcome_r = cov_outcome_r,
                             n_cov_ancova = n_cov_ancova,
                             n_exp = n_exp, n_nexp = n_nexp,
                             smd_to_cor = smd_to_cor, smd_var = smd_var,
                             reverse_ancova_md = reverse_ancova_md)

  es$info_used <- "ancova_md_se"

  return(es)
}


#' Convert an adjusted mean difference and adjusted standard deviation between two independent groups obtained from an ANCOVA model into several effect size measures
#'
#' @param ancova_md adjusted mean difference between two independent groups
#' @param ancova_md_ci_lo lower bound of the covariate-adjusted 95% CI of the mean difference
#' @param ancova_md_ci_up upper bound of the covariate-adjusted 95% CI of the mean difference
#' @param cov_outcome_r correlation between the outcome and covariate (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_ancova_md a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean difference (MD) 95% CI into a standard error,
#' and then relies on the calculations of the \code{\link{es_from_ancova_md_se}} function.
#'
#' **To convert the 95% CI into a standard error,** the following formula is used (table 12.3 in Cooper):
#' \deqn{md\_se = \frac{ancova\_md\_ci\_up - ancova\_md\_ci\_lo}{(2 * qt(0.975, n\_exp + n\_nexp - 2 - n\_cov\_ancova))}}
#' Calculations of the \code{\link{es_from_ancova_md_se}()} are then applied.
#'
#' @export es_from_ancova_md_ci
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 20. Adjusted: Mean difference and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_ancova_md_ci(
#'   ancova_md = 4, ancova_md_ci_lo = 2,
#'   ancova_md_ci_up = 6,
#'   cov_outcome_r = 0.5, n_cov_ancova = 5,
#'   n_exp = 20, n_nexp = 22
#' )
es_from_ancova_md_ci <- function(ancova_md, ancova_md_ci_lo, ancova_md_ci_up,
                                 cov_outcome_r, n_cov_ancova,
                                 n_exp, n_nexp,
                                 smd_to_cor = "viechtbauer", smd_var = "borenstein",
                                 reverse_ancova_md) {
  if (missing(reverse_ancova_md)) reverse_ancova_md <- rep(FALSE, length(ancova_md))
  reverse_ancova_md[is.na(reverse_ancova_md)] <- FALSE
  if (length(reverse_ancova_md) == 1) reverse_ancova_md = c(rep(reverse_ancova_md, length(ancova_md)))
  if (length(reverse_ancova_md) != length(ancova_md)) stop("The length of the 'reverse_ancova_md' argument is incorrectly specified.")

  ancova_md_se <- (ancova_md_ci_up - ancova_md_ci_lo) / (2 * qt(0.975, n_exp + n_nexp - 2 - n_cov_ancova))

  es <- es_from_ancova_md_se(ancova_md = ancova_md,
                             ancova_md_se = ancova_md_se,
                             cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
                             n_exp = n_exp, n_nexp = n_nexp,
                             smd_to_cor = smd_to_cor, smd_var = smd_var,
                             reverse_ancova_md = reverse_ancova_md)

  es$info_used <- "ancova_md_ci"

  return(es)
}

#' Convert an adjusted mean difference and adjusted standard deviation between two independent groups obtained from an ANCOVA model into several effect size measures
#'
#' @param ancova_md adjusted mean difference between two independent groups
#' @param ancova_md_pval p-value (two-tailed) of the adjusted mean difference
#' @param cov_outcome_r correlation between the outcome and covariate (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples.
#' @param reverse_ancova_md a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean difference (MD) p-value into a standard error,
#' and then relies on the calculations of the \code{\link{es_from_ancova_md_se}()} function.
#'
#' **To convert the p-value into a standard error,** the following formula is used (table 12.3 in Cooper):
#' \deqn{t = qt(p = \frac{ancova\_md\_pval}{2}, df = n\_exp + n\_nexp - 2 - n\_cov\_ancova)}
#' \deqn{ancova\_md\_se = | \frac{ancova\_md}{t} |}
#' Calculations of the \code{\link{es_from_ancova_md_se}()} are then applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 20. Adjusted: Mean difference and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_ancova_md_pval
#'
#' @md
#'
#' @examples
#' es_from_ancova_md_pval(
#'   ancova_md = 4, ancova_md_pval = 0.05,
#'   cov_outcome_r = 0.5, n_cov_ancova = 5,
#'   n_exp = 20, n_nexp = 22
#' )
es_from_ancova_md_pval <- function(ancova_md, ancova_md_pval,
                                   cov_outcome_r, n_cov_ancova,
                                   n_exp, n_nexp,
                                   smd_to_cor = "viechtbauer", smd_var = "borenstein",
                                   reverse_ancova_md) {
  if (missing(reverse_ancova_md)) reverse_ancova_md <- rep(FALSE, length(ancova_md))
  reverse_ancova_md[is.na(reverse_ancova_md)] <- FALSE
  if (length(reverse_ancova_md) == 1) reverse_ancova_md = c(rep(reverse_ancova_md, length(ancova_md)))
  if (length(reverse_ancova_md) != length(ancova_md)) stop("The length of the 'reverse_ancova_md' argument is incorrectly specified.")

  t <- qt(p = ancova_md_pval / 2,
          df = n_exp + n_nexp - 2 - n_cov_ancova,
          lower.tail = FALSE)

  ancova_md_se <- abs(ancova_md / t)

  es <- es_from_ancova_md_se(ancova_md = ancova_md,
                             ancova_md_se = ancova_md_se,
                             cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
                             n_exp = n_exp, n_nexp = n_nexp,
                             smd_to_cor = smd_to_cor, smd_var = smd_var,
                             reverse_ancova_md = reverse_ancova_md)

  es$info_used <- "ancova_md_pval"

  return(es)
}
