#' Convert a mean difference between two independent groups and standard deviation into several effect size measures
#'
#' @param md mean difference between two independent groups
#' @param md_sd standard deviation of the mean difference
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_md a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean difference and 95% CI into a Cohen's d (D) and Hedges' g (G).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **The formula used to obtain the Cohen's d** is:
#' \deqn{d = \frac{md}{md\_sd}}
#' Note that this formula is perfectly accurate only if the \code{md_sd} has been estimated by assuming that the variance of the two groups is equal.
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
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 10. Mean difference and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_md_sd
#'
#' @md
#'
#' @examples
#' es_from_md_sd(md = 4, md_sd = 2, n_exp = 20, n_nexp = 22)
es_from_md_sd <- function(md, md_sd, n_exp, n_nexp, smd_to_cor = "viechtbauer", reverse_md) {
  if (missing(reverse_md)) reverse_md <- rep(FALSE, length(md))
  reverse_md[is.na(reverse_md)] <- FALSE
  if (length(reverse_md) == 1) reverse_md = c(rep(reverse_md, length(md)))
  if (length(reverse_md) != length(md)) stop("The length of the 'reverse_md' argument is incorrectly specified.")

  d <- md / md_sd

  d_se <- sqrt((n_exp + n_nexp) / (n_exp * n_nexp) + (d^2) / (2 * (n_exp + n_nexp)))

  es <- .es_from_d(
    d = d, d_se = d_se, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_md
  )

  es$info_used <- "md_sd"

  es$md <- ifelse(reverse_md, -md, md)
  es$md_se <- ifelse(!is.na(es$md), md_sd * sqrt(1 / n_exp + 1 / n_nexp), NA)
  es$md_ci_lo <- ifelse(!is.na(es$md), es$md - qt(.975, n_exp + n_nexp - 2) * es$md_se, NA)
  es$md_ci_up <- ifelse(!is.na(es$md), es$md + qt(.975, n_exp + n_nexp - 2) * es$md_se, NA)

  return(es)
}

#' Convert a mean difference between two independent groups and its standard error into several effect size measures
#'
#' @param md mean difference between two independent groups
#' @param md_se standard error of the mean difference
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_md a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function the standard error of a mean difference into a
#' standard deviation:
#' \deqn{inv\_n = \frac{1}{n\_exp} + \frac{1}{n\_nexp}}
#' \deqn{md\_sd = \frac{md\_se}{\sqrt{inv\_n}}}
#'
#' Calculations of the \code{\link{es_from_md_sd}} function are then used to estimate
#' the Cohen's d and other effect size measures.
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 10. Mean difference and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_md_se
#'
#' @md
#'
#' @examples
#' es_from_md_se(md = 4, md_se = 2, n_exp = 20, n_nexp = 22)
es_from_md_se <- function(md, md_se, n_exp, n_nexp, smd_to_cor = "viechtbauer", reverse_md) {
  if (missing(reverse_md)) reverse_md <- rep(FALSE, length(md))
  reverse_md[is.na(reverse_md)] <- FALSE

  md_sd <- md_se / sqrt(1 / n_exp + 1 / n_nexp)

  es <- es_from_md_sd(
    md = md, md_sd = md_sd,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_md = reverse_md
  )

  es$info_used <- "md_se"

  return(es)
}

#' Convert a mean difference between two independent groups and 95% CI into several effect size measures
#'
#' @param md mean difference between two independent groups
#' @param md_ci_lo lower bound of the 95% CI of the mean difference
#' @param md_ci_up upper bound of the 95% CI of the mean difference
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param reverse_md a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts 95% CI of a mean difference into a standard error
#' (Cochrane Handbook section 6.5.2.3):
#' \deqn{md\_se = \frac{md\_ci\_up - md\_ci\_lo}{2 * qt(0.975, df = n\_exp + n\_nexp - 2)}}
#'
#' Calculations of the \code{\link{es_from_md_se}()} function are
#' then used to estimate the Cohen's d and other effect size measures.
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 10. Mean difference and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_md_ci
#'
#' @md
#'
#' @examples
#' es_from_md_ci(md = 4, md_ci_lo = 2, md_ci_up = 6, n_exp = 20, n_nexp = 22)
es_from_md_ci <- function(md, md_ci_lo, md_ci_up, n_exp, n_nexp,
                          smd_to_cor = "viechtbauer", max_asymmetry = 10, reverse_md) {
  if (missing(reverse_md)) reverse_md <- rep(FALSE, length(md))
  reverse_md[is.na(reverse_md)] <- FALSE



  md_se <- (md_ci_up - md_ci_lo) / (2 * qt(0.975, n_exp + n_nexp - 2))

  es <- es_from_md_se(
    md = md, md_se = md_se,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_md = reverse_md
  )

  es$info_used <- "md_ci"

  return(es)
}

#' Convert a mean difference between two independent groups and its p-value into several effect size measures
#'
#' @param md mean difference between two independent groups
#' @param md_pval p-value of the mean difference
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_md a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the p-value of a mean difference
#' into a standard error (Cochrane Handbook section 6.5.2.3):
#' \deqn{t = qt(\frac{md\_pval}{2}, df = n\_exp + n\_nexp - 2)}
#' \deqn{md\_se = |\frac{md}{t}|}
#' Calculations of the \code{\link{es_from_md_se}} function are then used to estimate the Cohen's d and other effect size measures.
#'
#' @references
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 10. Mean difference and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_md_pval
#'
#' @md
#'
#' @examples
#' es_from_md_pval(md = 4, md_pval = 0.024, n_exp = 20, n_nexp = 22)
es_from_md_pval <- function(md, md_pval, n_exp, n_nexp, smd_to_cor = "viechtbauer", reverse_md) {
  if (missing(reverse_md)) reverse_md <- rep(FALSE, length(md))
  reverse_md[is.na(reverse_md)] <- FALSE

  t <- qt(p = md_pval / 2, df = n_exp + n_nexp - 2, lower.tail = FALSE)

  md_se <- abs(md / t)

  es <- es_from_md_se(
    md = md, md_se = md_se,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_md = reverse_md
  )

  es$info_used <- "md_pval"

  return(es)
}

