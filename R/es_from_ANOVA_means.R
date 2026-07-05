#' Convert means and standard deviations of two independent groups into several effect size measures
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_sd_exp standard deviation of participants in the experimental/exposed group.
#' @param mean_sd_nexp standard deviation of participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first computes a Cohen's d (D), Hedges' g (G) and mean difference (MD)
#' from the means and standard deviations of two independent groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a mean difference**  (formulas 12.1-12.6 in Cooper):
#' \deqn{md = mean\_exp - mean\_nexp}
#' \deqn{md\_se = \sqrt{\frac{mean\_sd\_exp^2}{n\_exp} + \frac{mean\_sd\_nexp^2}{n\_nexp}}}
#' \deqn{md\_ci\_lo = md - md\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#' \deqn{md\_ci\_up = md + md\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#'
#' **To estimate a Cohen's d** the following formulas are used (formulas 12.10-12.18 in Cooper):
#' \deqn{mean\_sd\_pooled = \sqrt{\frac{(n\_exp - 1) * sd\_exp^2 + (n\_nexp - 1) * sd\_nexp^2}{n\_exp+n\_nexp-2}}}
#' \deqn{cohen\_d =  \frac{mean\_exp - mean\_nexp}{mean\_sd\_pooled}}
#' \deqn{cohen\_d\_se = \frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2(n\_exp+n\_nexp)}}
#' \deqn{cohen\_d\_ci\_lo = cohen\_d - cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#' \deqn{cohen\_d\_ci\_up = cohen\_d + cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_means_sd
#'
#' @md
#'
#' @examples
#' es_from_means_sd(
#'   n_exp = 55, n_nexp = 55,
#'   mean_exp = 2.3, mean_sd_exp = 1.2,
#'   mean_nexp = 1.9, mean_sd_nexp = 0.9
#' )
es_from_means_sd <- function(mean_exp, mean_sd_exp, mean_nexp, mean_sd_nexp, n_exp, n_nexp,
                             smd_to_cor = "viechtbauer", reverse_means) {
  if (missing(reverse_means)) reverse_means <- rep(FALSE, length(mean_exp))
  reverse_means[is.na(reverse_means)] <- FALSE
  if (length(reverse_means) == 1) reverse_means = c(rep(reverse_means, length(mean_exp)))
  if (length(reverse_means) != length(mean_exp)) stop("The length of the 'reverse_means' argument is incorrectly specified.")

  pooled_sd <- sqrt(((n_exp - 1) * mean_sd_exp^2 + (n_nexp - 1) * mean_sd_nexp^2) / (n_exp + n_nexp - 2))

  d <- (mean_exp - mean_nexp) / pooled_sd

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_means
  )

  es$info_used <- "means_sd"

  es$md <- ifelse(reverse_means, mean_nexp - mean_exp, mean_exp - mean_nexp)
  es$md <- ifelse(!is.na(es$md) & !is.na(mean_sd_exp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp),
  es$md, NA
  )
  es$md_se <- ifelse(!is.na(es$md), sqrt(mean_sd_exp^2 / n_exp + mean_sd_nexp^2 / n_nexp), NA)
  es$md_ci_lo <- ifelse(!is.na(es$md), es$md - qt(.975, n_exp + n_nexp - 2) * es$md_se, NA)
  es$md_ci_up <- ifelse(!is.na(es$md), es$md + qt(.975, n_exp + n_nexp - 2) * es$md_se, NA)

  return(es)
}

#' Convert means and standard errors of two independent groups several effect size measures
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_se_exp standard error of participants in the experimental/exposed group.
#' @param mean_se_nexp standard error of participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the standard errors of two independent groups into standard deviations,
#' and then relies on the calculations of the \code{\link{es_from_means_sd}()} function.
#'
#' **To convert the standard errors into standard deviations**, the following formula is used.
#' \deqn{mean\_sd\_exp = mean\_se\_exp * \sqrt{n\_exp}}
#' \deqn{mean\_sd\_nexp = mean\_se\_nexp * \sqrt{n\_nexp}}
#' Then, calculations of the \code{\link{es_from_means_sd}()} are applied.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @export es_from_means_se
#'
#' @md
#'
#' @examples
#' es_from_means_se(
#'   mean_exp = 42, mean_se_exp = 11,
#'   mean_nexp = 42, mean_se_nexp = 15,
#'   n_exp = 43, n_nexp = 34
#' )
es_from_means_se <- function(mean_exp, mean_se_exp, mean_nexp, mean_se_nexp, n_exp, n_nexp,
                             smd_to_cor = "viechtbauer", reverse_means) {
  if (missing(reverse_means)) reverse_means <- rep(FALSE, length(mean_exp))
  reverse_means[is.na(reverse_means)] <- FALSE

  sd_exp <- mean_se_exp * sqrt(n_exp)
  sd_nexp <- mean_se_nexp * sqrt(n_nexp)

  es <- es_from_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = mean_exp, mean_sd_exp = sd_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = sd_nexp,
    smd_to_cor = smd_to_cor, reverse_means = reverse_means
  )

  es$info_used <- "means_se"

  return(es)
}

#' Convert means of two groups and the pooled standard deviation into several effect size measures
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_sd_pooled pooled standard deviation across both groups.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first computes a Cohen's d (D), Hedges' g (G) and mean difference (MD)
#' from the means of two independent groups and the pooled standard deviation across the groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a mean difference**  (formulas 12.1-12.6 in Cooper):
#' \deqn{md = mean\_exp - mean\_nexp}
#' \deqn{md\_se = \sqrt{\frac{n_exp+n_nexp}{n_exp*n_nexp} * mean_sd_pooled^2}}
#' \deqn{md\_ci\_lo = md - md\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#' \deqn{md\_ci\_up = md + md\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#'
#' **To estimate a Cohen's d** the following formulas are used (formulas 12.10-12.18 in Cooper):
#' \deqn{cohen\_d =  \frac{mean\_exp - mean\_nexp}{means\_sd\_pooled}}
#' \deqn{cohen\_d\_se = \frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2(n\_exp+n\_nexp)}}
#' \deqn{cohen\_d\_ci\_lo = cohen\_d - cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#' \deqn{cohen\_d\_ci\_up = cohen\_d + cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_means_sd_pooled
#'
#' @md
#'
#' @examples
#' es_from_means_sd_pooled(
#'   n_exp = 55, n_nexp = 55,
#'   mean_exp = 2.3, mean_nexp = 1.9,
#'   mean_sd_pooled = 0.9
#' )
es_from_means_sd_pooled <- function(mean_exp, mean_nexp, mean_sd_pooled, n_exp, n_nexp,
                                    smd_to_cor = "viechtbauer", reverse_means) {
  if (missing(reverse_means)) reverse_means <- rep(FALSE, length(mean_exp))
  reverse_means[is.na(reverse_means)] <- FALSE
  if (length(reverse_means) == 1) reverse_means = c(rep(reverse_means, length(mean_exp)))
  if (length(reverse_means) != length(mean_exp)) stop("The length of the 'reverse_means' argument is incorrectly specified.")

  d <- (mean_exp - mean_nexp) / mean_sd_pooled

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_means
  )

  es$info_used <- "means_sd_pooled"

  es$md <- ifelse(reverse_means, mean_nexp - mean_exp, mean_exp - mean_nexp)
  es$md <- ifelse(!is.na(es$md) & !is.na(mean_sd_pooled) &
    !is.na(n_exp) & !is.na(n_nexp),
  es$md, NA
  )
  es$md_se <- ifelse(!is.na(es$md), sqrt((n_exp + n_nexp) / (n_exp * n_nexp) * mean_sd_pooled^2), NA)
  es$md_ci_lo <- ifelse(!is.na(es$md), es$md - qt(.975, n_exp + n_nexp - 2) * es$md_se, NA)
  es$md_ci_up <- ifelse(!is.na(es$md), es$md + qt(.975, n_exp + n_nexp - 2) * es$md_se, NA)

  return(es)
}

#' Convert means and 95% CI of two independent groups several effect size measures
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_ci_lo_exp lower bound of the 95% CI of the mean of the experimental/exposed group
#' @param mean_ci_up_exp upper bound of the 95% CI of the mean of the experimental/exposed group
#' @param mean_ci_lo_nexp lower bound of the 95% CI of the mean of the non-experimental/non-exposed group.
#' @param mean_ci_up_nexp upper bound of the 95% CI of the mean of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param reverse_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the 95% CI of two independent groups into a standard error,
#' and then relies on the calculations of the \code{\link{es_from_means_se}()} function.
#'
#' **To convert the 95% CIs into standard errors,** the following formula is used (table 12.3 in Cooper):
#' \deqn{mean\_se\_exp = \frac{mean\_ci\_up\_exp - mean\_ci\_lo\_exp}{2 * qt{(0.975, df = n\_exp - 1)}}}
#' \deqn{mean\_se\_nexp = \frac{mean\_ci\_up\_nexp - mean\_ci\_lo\_nexp}{2 * qt{(0.975, df = n\_nexp - 1)}}}
#' Calculations of the \code{\link{es_from_means_se}()} are then applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @export es_from_means_ci
#'
#' @md
#'
#' @examples
#' es_from_means_ci(
#'   n_exp = 55, n_nexp = 55,
#'   mean_exp = 25, mean_ci_lo_exp = 15, mean_ci_up_exp = 35,
#'   mean_nexp = 18, mean_ci_lo_nexp = 12, mean_ci_up_nexp = 24
#' )
es_from_means_ci <- function(mean_exp, mean_ci_lo_exp, mean_ci_up_exp,
                             mean_nexp, mean_ci_lo_nexp, mean_ci_up_nexp,
                             n_exp, n_nexp, smd_to_cor = "viechtbauer",
                             max_asymmetry = 10,
                             reverse_means) {
  if (missing(reverse_means)) reverse_means <- rep(FALSE, length(mean_exp))
  reverse_means[is.na(reverse_means)] <- FALSE

  df_exp <- n_exp - 1
  df_nexp <- n_nexp - 1

  se_exp <- (mean_ci_up_exp - mean_ci_lo_exp) / (2 * qt(0.975, df_exp))
  se_nexp <- (mean_ci_up_nexp - mean_ci_lo_nexp) / (2 * qt(0.975, df_nexp))

  es <- es_from_means_se(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = mean_exp,
    mean_se_exp = se_exp,
    mean_nexp = mean_nexp, mean_se_nexp = se_nexp,
    smd_to_cor = smd_to_cor, reverse_means = reverse_means
  )

  es$info_used <- "means_ci"

  return(es)
}
