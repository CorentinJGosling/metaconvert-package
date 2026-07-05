#' Convert mean changes and standard deviations of two independent groups into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the experimental/exposed group.
#' @param mean_change_sd_exp standard deviation of the mean change (i.e., SD of the difference scores)
#'   for participants in the experimental/exposed group.
#' @param mean_change_nexp mean change of participants in the non-experimental/non-exposed group.
#' @param mean_change_sd_nexp standard deviation of the mean change for participants in the
#'   non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group (only used with \code{pre_post_to_smd = "morris_drm"}, see details).
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group (only used with \code{pre_post_to_smd = "morris_drm"}, see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the mean change into a SMD ("morris_drm" or "morris_dz", see details).
#' @param pool_sd a logical value indicating whether the SD used to standardize the effect size should be pooled across the two groups.
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function computes a Cohen's d (D) and Hedges' g (G) from the mean change and standard deviation
#' of change scores of two independent groups. Odds ratio (OR) and correlation coefficients (R/Z) are
#' then converted from Cohen's d.
#'
#' Two formulas can be used to obtain the SMD (Morris & DeShon, 2002):
#' \deqn{d_{z} = \frac{mean\_change}{sd\_change}}
#' \deqn{d_{rm} = d_{z} \times \sqrt{2 \times (1 - r\_pre\_post)}}
#' The 'morris_drm' formula (default, alias 'cooper') requires the pre-post correlation while
#' 'morris_dz' does not. Note that d_rm and d_z are not expressed on the same scale and should not be
#' combined in a same meta-analysis (Morris & DeShon, 2002). Under 'morris_drm', the resulting effect
#' size directly depends on the r_pre_post value entered; if r is unknown, the 'morris_dz' formula
#' can be used instead. The 'bonett' and 'morris_dav' formulas require separate pre/post SDs and are
#' thus not available for mean change data.
#'
#' This function simply internally calls the \code{\link{es_from_means_sd_pre_post}} function but setting:
#' \deqn{mean\_pre\_exp = 0}
#' \deqn{mean\_pre\_sd\_exp = 0}
#' \deqn{mean\_exp = mean\_change\_exp}
#' \deqn{mean\_sd\_exp = mean\_change\_sd\_exp}
#' \deqn{mean\_pre\_nexp = 0}
#' \deqn{mean\_pre\_sd\_nexp = 0}
#' \deqn{mean\_nexp = mean\_change\_nexp}
#' \deqn{mean\_sd\_nexp = mean\_change\_sd\_nexp}
#'
#' To know more about the calculations, see \code{\link{es_from_means_sd_pre_post}} function.
#'
#' @references
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 14. Paired: mean change, and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_mean_change_sd
#'
#' @md
#'
#' @examples
#' es_from_mean_change_sd(
#'   n_exp = 36, n_nexp = 35,
#'   mean_change_exp = 8.4, mean_change_sd_exp = 9.13,
#'   mean_change_nexp = 2.43, mean_change_sd_nexp = 6.61,
#'   r_pre_post_exp = 0.2, r_pre_post_nexp = 0.2
#' )
es_from_mean_change_sd <- function(mean_change_exp, mean_change_sd_exp,
                                   mean_change_nexp, mean_change_sd_nexp,
                                   r_pre_post_exp, r_pre_post_nexp,
                                   n_exp, n_nexp,
                                   smd_to_cor = "viechtbauer",
                                   pre_post_to_smd = "cooper",
                                   pool_sd = FALSE,
                                   reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(mean_change_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "mean_change",
    func_name = "es_from_mean_change_sd"
  )

  es <- es_from_means_sd_pre_post(
    mean_pre_exp = rep(0, length(mean_change_exp)),
    mean_exp = mean_change_exp,
    mean_pre_sd_exp = rep(0, length(mean_change_exp)),
    mean_sd_exp = mean_change_sd_exp,
    mean_pre_nexp = rep(0, length(mean_change_exp)),
    mean_nexp = mean_change_nexp,
    mean_pre_sd_nexp = rep(0, length(mean_change_exp)),
    mean_sd_nexp = mean_change_sd_nexp,
    n_exp = n_exp,
    n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    pool_sd = pool_sd,
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_sd"

  return(es)
}
#' Convert mean changes and standard errors of two independent groups into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the experimental/exposed group.
#' @param mean_change_se_exp standard error of the mean change for participants in the experimental/exposed group.
#' @param mean_change_nexp mean change of participants in the non-experimental/non-exposed group.
#' @param mean_change_se_nexp standard error of the mean change for participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group (only used with \code{pre_post_to_smd = "morris_drm"}, see details).
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group (only used with \code{pre_post_to_smd = "morris_drm"}, see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the mean change into a SMD ("morris_drm" or "morris_dz", see details).
#' @param pool_sd a logical value indicating whether the SD used to standardize the effect size should be pooled across the two groups.
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean change and standard errors of two independent groups into a Cohen's d.
#' The Cohen's d is then converted to other effect size measures.
#'
#' This function simply internally calls the \code{\link{es_from_means_se_pre_post}} function but setting:
#' \deqn{mean\_pre\_exp = 0}
#' \deqn{mean\_pre\_se\_exp = 0}
#' \deqn{mean\_exp = mean\_change\_exp}
#' \deqn{mean\_se\_exp = mean\_change\_se\_exp}
#' \deqn{mean\_pre\_nexp = 0}
#' \deqn{mean\_pre\_se\_nexp = 0}
#' \deqn{mean\_nexp = mean\_change\_nexp}
#' \deqn{mean\_se\_nexp = mean\_change\_se\_nexp}
#'
#' To know more about the calculations, see \code{\link{es_from_means_se_pre_post}} function.
#'
#' @references
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 14. Paired: mean change, and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_mean_change_se
#'
#' @md
#'
#' @examples
#' es_from_mean_change_se(
#'   n_exp = 36, n_nexp = 35,
#'   mean_change_exp = 8.4, mean_change_se_exp = 9.13,
#'   mean_change_nexp = 2.43, mean_change_se_nexp = 6.61,
#'   r_pre_post_exp = 0.2, r_pre_post_nexp = 0.2
#' )
es_from_mean_change_se <- function(mean_change_exp, mean_change_se_exp,
                                   mean_change_nexp, mean_change_se_nexp,
                                   r_pre_post_exp, r_pre_post_nexp,
                                   n_exp, n_nexp,
                                   smd_to_cor = "viechtbauer",
                                   pre_post_to_smd = "cooper",
                                   pool_sd = FALSE,
                                   reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(mean_change_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "mean_change",
    func_name = "es_from_mean_change_se"
  )

  es <- es_from_means_se_pre_post(
    mean_pre_exp = rep(0, length(mean_change_exp)),
    mean_exp = mean_change_exp,
    mean_pre_se_exp = rep(0, length(mean_change_exp)),
    mean_se_exp = mean_change_se_exp,
    mean_pre_nexp = rep(0, length(mean_change_exp)),
    mean_nexp = mean_change_nexp,
    mean_pre_se_nexp = rep(0, length(mean_change_exp)),
    mean_se_nexp = mean_change_se_nexp,
    n_exp = n_exp,
    n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    pool_sd = pool_sd,
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_se"

  return(es)
}
#' Convert mean changes and 95% CI of two independent groups into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the experimental/exposed group.
#' @param mean_change_ci_lo_exp lower bound of the 95% CI around the mean change of the experimental/exposed group.
#' @param mean_change_ci_up_exp upper bound of the 95% CI around the mean change of the experimental/exposed group.
#' @param mean_change_nexp mean change of participants in the non-experimental/non-exposed group.
#' @param mean_change_ci_lo_nexp lower bound of the 95% CI around the mean change of the non-experimental/non-exposed group.
#' @param mean_change_ci_up_nexp upper bound of the 95% CI around the mean change of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group (only used with \code{pre_post_to_smd = "morris_drm"}, see details).
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group (only used with \code{pre_post_to_smd = "morris_drm"}, see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param pool_sd a logical value indicating whether the SD used to standardize the effect size should be pooled across the two groups.
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#' @param pre_post_to_smd formula used to convert the mean change into a SMD ("morris_drm" or "morris_dz", see details).
#'
#' @details
#' This function converts the mean change and 95% CI of two independent groups into a Cohen's d.
#' The Cohen's d is then converted to other effect size measures.
#'
#' This function simply internally calls the \code{\link{es_from_means_ci_pre_post}} function but setting:
#' \deqn{mean\_pre\_exp = 0}
#' \deqn{mean\_pre\_ci\_lo\_exp = 0}
#' \deqn{mean\_pre\_ci\_up\_exp = 0}
#'
#' \deqn{mean\_exp = mean\_change\_exp}
#' \deqn{mean\_ci\_lo\_exp = mean\_change\_ci\_lo\_exp}
#' \deqn{mean\_ci\_up\_exp = mean\_change\_ci\_up\_exp}
#'
#' \deqn{mean\_pre\_nexp = 0}
#' \deqn{mean\_pre\_ci\_lo\_nexp = 0}
#' \deqn{mean\_pre\_ci\_up\_nexp = 0}
#'
#' \deqn{mean\_nexp = mean\_change\_nexp}
#' \deqn{mean\_ci\_lo\_nexp = mean\_change\_ci\_lo\_nexp}
#' \deqn{mean\_ci\_up\_nexp = mean\_change\_ci\_up\_nexp}
#'
#' To know more about the calculations, see \code{\link{es_from_means_sd_pre_post}} function.
#'
#' @references
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @export es_from_mean_change_ci
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 14. Paired: mean change, and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_mean_change_ci(
#'   n_exp = 36, n_nexp = 35,
#'   mean_change_exp = 8.4,
#'   mean_change_ci_lo_exp = 6.4, mean_change_ci_up_exp = 10.4,
#'   mean_change_nexp = 2.43,
#'   mean_change_ci_lo_nexp = 1.43, mean_change_ci_up_nexp = 3.43,
#'   r_pre_post_exp = 0.2, r_pre_post_nexp = 0.2
#' )
es_from_mean_change_ci <- function(mean_change_exp,
                                   mean_change_ci_lo_exp, mean_change_ci_up_exp,
                                   mean_change_nexp,
                                   mean_change_ci_lo_nexp, mean_change_ci_up_nexp,
                                   r_pre_post_exp, r_pre_post_nexp,
                                   n_exp, n_nexp, max_asymmetry = 10,
                                   smd_to_cor = "viechtbauer", pre_post_to_smd = "cooper",
                                   pool_sd = FALSE, reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(mean_change_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "mean_change",
    func_name = "es_from_mean_change_ci"
  )

  es <- es_from_means_ci_pre_post(
    mean_pre_exp = rep(0, length(mean_change_exp)),
    mean_pre_ci_lo_exp = rep(0, length(mean_change_exp)),
    mean_pre_ci_up_exp = rep(0, length(mean_change_exp)),

    mean_pre_nexp = rep(0, length(mean_change_exp)),
    mean_pre_ci_lo_nexp = rep(0, length(mean_change_exp)),
    mean_pre_ci_up_nexp = rep(0, length(mean_change_exp)),

    mean_exp = mean_change_exp, mean_ci_lo_exp = mean_change_ci_lo_exp, mean_ci_up_exp = mean_change_ci_up_exp,
    mean_nexp = mean_change_nexp, mean_ci_lo_nexp = mean_change_ci_lo_nexp, mean_ci_up_nexp = mean_change_ci_up_nexp,
    n_exp = n_exp,
    n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    pool_sd = pool_sd,
    max_asymmetry = max_asymmetry,
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_ci"

  return(es)
}


#' Convert mean changes and p-values of two independent groups into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the experimental/exposed group.
#' @param mean_change_pval_exp p-value of the mean change for participants in the experimental/exposed group.
#' @param mean_change_nexp mean change of participants in the non-experimental/non-exposed group.
#' @param mean_change_pval_nexp p-value of the mean change for participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group (only used with \code{pre_post_to_smd = "morris_drm"}, see details).
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group (only used with \code{pre_post_to_smd = "morris_drm"}, see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pool_sd a logical value indicating whether the SD used to standardize the effect size should be pooled across the two groups.
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#' @param pre_post_to_smd formula used to convert the mean change into a SMD ("morris_drm" or "morris_dz", see details).
#'
#' @details
#' This function converts the mean change and associated p-values of two independent groups into a Cohen's d.
#' The Cohen's d is then converted to other effect size measures.
#'
#' To start, this function estimates the mean change standard errors from the p-values:
#' \deqn{t\_exp <- qt(p = mean\_change\_pval\_exp / 2, df = n\_exp - 1, lower.tail = FALSE)}
#' \deqn{t\_nexp <- qt(p = mean\_change\_pval\_nexp / 2, df = n\_nexp - 1, lower.tail = FALSE)}
#' \deqn{mean\_change\_se\_exp <- |\frac{mean\_change\_exp}{t\_exp}|}
#' \deqn{mean\_change\_se\_nexp <- |\frac{mean\_change\_nexp}{t\_nexp}|}
#'
#' Then, this function simply internally calls the \code{\link{es_from_means_se_pre_post}} function but setting:
#' \deqn{mean\_pre\_exp = 0}
#' \deqn{mean\_pre\_se\_exp = 0}
#' \deqn{mean\_exp = mean\_change\_exp}
#' \deqn{mean\_se\_exp = mean\_change\_se\_exp}
#' \deqn{mean\_pre\_nexp = 0}
#' \deqn{mean\_pre\_se\_nexp = 0}
#' \deqn{mean\_nexp = mean\_change\_nexp}
#' \deqn{mean\_se\_nexp = mean\_change\_se\_nexp}
#'
#' To know more about other calculations, see \code{\link{es_from_means_sd_pre_post}} function.
#'
#' @references
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @export es_from_mean_change_pval
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 14. Paired: mean change, and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_mean_change_pval(
#'   n_exp = 36, n_nexp = 35,
#'   mean_change_exp = 8.4, mean_change_pval_exp = 0.13,
#'   mean_change_nexp = 2.43, mean_change_pval_nexp = 0.61,
#'   r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8
#' )
es_from_mean_change_pval <- function(mean_change_exp, mean_change_pval_exp,
                                   mean_change_nexp, mean_change_pval_nexp,
                                   r_pre_post_exp, r_pre_post_nexp,
                                   n_exp, n_nexp,
                                   smd_to_cor = "viechtbauer", pre_post_to_smd = "cooper",
                                   pool_sd = FALSE, reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(mean_change_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "mean_change",
    func_name = "es_from_mean_change_pval"
  )

  t_exp <- qt(p = mean_change_pval_exp / 2, df = n_exp - 1, lower.tail = FALSE)
  t_nexp <- qt(p = mean_change_pval_nexp / 2, df = n_nexp - 1, lower.tail = FALSE)

  mean_change_se_exp <- abs(mean_change_exp / t_exp)
  mean_change_se_nexp <- abs(mean_change_nexp / t_nexp)

  es <- es_from_means_se_pre_post(
    mean_pre_exp = rep(0, length(mean_change_exp)),
    mean_exp = mean_change_exp,
    mean_pre_se_exp = rep(0, length(mean_change_exp)),
    mean_se_exp = mean_change_se_exp,
    mean_pre_nexp = rep(0, length(mean_change_exp)),
    mean_nexp = mean_change_nexp,
    mean_pre_se_nexp = rep(0, length(mean_change_exp)),
    mean_se_nexp = mean_change_se_nexp,
    n_exp = n_exp,
    n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    pool_sd = pool_sd,
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_pval"

  return(es)
}
