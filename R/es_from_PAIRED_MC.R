#' Convert mean changes and standard deviations of two independent groups into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the experimental/exposed group.
#' @param mean_change_sd_exp standard deviation of the mean change for participants in the experimental/exposed group.
#' @param mean_change_nexp mean change of participants in the non-experimental/non-exposed group.
#' @param mean_change_sd_nexp standard deviation of the mean change for participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first computes a Cohen's d (D), Hedges' g (G)
#' from the mean change (MC) and standard deviations of two independent groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' This function simply internally calls the \code{\link{es_from_means_sd_pre_post}} function but setting:
#' \deqn{mean\_pre\_exp = mean\_change\_exp}
#' \deqn{mean\_pre\_sd\_exp = mean\_change\_sd\_exp}
#' \deqn{mean\_exp = 0}
#' \deqn{mean\_sd\_exp = 0}
#' \deqn{mean\_pre\_nexp = mean\_change\_nexp}
#' \deqn{mean\_pre\_sd\_nexp = mean\_change\_sd\_nexp}
#' \deqn{mean\_nexp = 0}
#' \deqn{mean\_sd\_nexp = 0}
#'
#' To know more about the calculations, see \code{\link{es_from_means_sd_pre_post}} function.
#'
#' @references
#' Bonett, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364–386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
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
                                   smd_to_cor = "viechtbauer", reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.5, length(mean_change_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.5
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.5, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.5

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       mean_change_sd_exp, mean_change_sd_nexp,
                       error_message = paste0("The number of people exposed/non-exposed and standard deviations ",
                                              "should be >0."),
                       func = "es_from_mean_change_sd")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  es <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_change_exp,
    mean_exp = 0,
    mean_pre_sd_exp = mean_change_sd_exp,
    mean_sd_exp = 0,
    mean_pre_nexp = mean_change_nexp,
    mean_nexp = 0,
    mean_pre_sd_nexp = mean_change_sd_nexp,
    mean_sd_nexp = 0,
    n_exp = n_exp,
    n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = "cooper",
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
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean change and standard errors of two independent groups
#' into a Cohen's d. The Cohen's d is then converted to other effect size measures.
#'
#' This function simply internally calls the \code{\link{es_from_means_se_pre_post}} function but setting:
#' \deqn{mean\_pre\_exp = mean\_change\_exp}
#' \deqn{mean\_pre\_se\_exp = mean\_change\_se\_exp}
#' \deqn{mean\_exp = 0}
#' \deqn{mean\_se\_exp = 0}
#' \deqn{mean\_pre\_nexp = mean\_change\_nexp}
#' \deqn{mean\_pre\_se\_nexp = mean\_change\_se\_nexp}
#' \deqn{mean\_nexp = 0}
#' \deqn{mean\_se\_nexp = 0}
#'
#' To know more about the calculations, see \code{\link{es_from_means_se_pre_post}} function.
#'
#' @references
#' Bonett, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364–386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
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
                                   smd_to_cor = "viechtbauer", reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.5, length(mean_change_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.5
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.5, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.5

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       mean_change_se_exp, mean_change_se_nexp,
                       error_message = paste0("The number of people exposed/non-exposed and standard errors ",
                                              "should be >0."),
                       func = "es_from_mean_change_se")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  es <- es_from_means_se_pre_post(
    mean_pre_exp = mean_change_exp,
    mean_exp = 0,
    mean_pre_se_exp = mean_change_se_exp,
    mean_se_exp = 0,
    mean_pre_nexp = mean_change_nexp,
    mean_nexp = 0,
    mean_pre_se_nexp = mean_change_se_nexp,
    mean_se_nexp = 0,
    n_exp = n_exp,
    n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = "cooper",
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_se"

  return(es)
}
#' Convert mean changes and standard deviations of two independent groups into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the experimental/exposed group.
#' @param mean_change_ci_lo_exp lower bound of the 95% CI around the mean change of the experimental/exposed group.
#' @param mean_change_ci_up_exp upper bound of the 95% CI around the mean change of the experimental/exposed group.
#' @param mean_change_nexp mean change of participants in the non-experimental/non-exposed group.
#' @param mean_change_ci_lo_nexp lower bound of the 95% CI around the mean change of the non-experimental/non-exposed group.
#' @param mean_change_ci_up_nexp upper bound of the 95% CI around the mean change of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean change and 95% CI of two independent groups
#' into a Cohen's d. The Cohen's d is then converted to other effect size measures.
#'
#' This function simply internally calls the \code{\link{es_from_means_ci_pre_post}} function but setting:
#' \deqn{mean\_pre\_exp = mean\_change\_exp}
#' \deqn{mean\_pre\_ci\_lo\_exp = mean\_change\_ci\_lo\_exp}
#' \deqn{mean\_pre\_ci\_up\_exp = mean\_change\_ci\_up\_exp}
#'
#' \deqn{mean\_exp = 0}
#' \deqn{mean\_ci\_lo\_exp = 0}
#' \deqn{mean\_ci\_up\_exp = 0}
#'
#' \deqn{mean\_pre\_nexp = mean\_change\_nexp}
#' \deqn{mean\_pre\_ci\_lo\_nexp = mean\_change\_ci\_lo\_nexp}
#' \deqn{mean\_pre\_ci\_up\_nexp = mean\_change\_ci\_up\_nexp}
#'
#' \deqn{mean\_nexp = 0}
#' \deqn{mean\_ci\_lo\_nexp = 0}
#' \deqn{mean\_ci\_up\_nexp = 0}
#'
#' To know more about the calculations, see \code{\link{es_from_means_sd_pre_post}} function.
#'
#' @references
#' Bonett, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364–386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
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
                                   smd_to_cor = "viechtbauer", reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.5, length(mean_change_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.5
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.5, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.5

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       error_message = paste0("The number of people exposed/non-exposed ",
                                              "should be >0."),
                       func = "es_from_mean_change_ci")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })
  tryCatch({
    .validate_ci_symmetry(mean_change_exp, mean_change_ci_lo_exp, mean_change_ci_up_exp,
                          func = "es_from_mean_change_ci",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })
  tryCatch({
    .validate_ci_symmetry(mean_change_nexp, mean_change_ci_lo_nexp, mean_change_ci_up_nexp,
                          func = "es_from_mean_change_ci",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })

  es <- es_from_means_ci_pre_post(
    mean_pre_exp = mean_change_exp,
    mean_pre_ci_lo_exp = mean_change_ci_lo_exp,
    mean_pre_ci_up_exp = mean_change_ci_up_exp,

    mean_pre_nexp = mean_change_nexp,
    mean_pre_ci_lo_nexp = mean_change_ci_lo_nexp,
    mean_pre_ci_up_nexp = mean_change_ci_up_nexp,

    mean_exp = 0, mean_ci_lo_exp = 0, mean_ci_up_exp = 0,
    mean_nexp = 0, mean_ci_lo_nexp = 0, mean_ci_up_nexp = 0,
    n_exp = n_exp,
    n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = "cooper",
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_ci"

  return(es)
}


#' Convert mean changes and standard deviations of two independent groups into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the experimental/exposed group.
#' @param mean_change_pval_exp p-value of the mean change for participants in the experimental/exposed group.
#' @param mean_change_nexp mean change of participants in the non-experimental/non-exposed group.
#' @param mean_change_pval_nexp p-value of the mean change for participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean change and associated p-values of two independent groups
#' into a Cohen's d. The Cohen's d is then converted to other effect size measures.
#'
#' To start, this function estimates the mean change standard errors from the p-values:
#' \deqn{t\_exp <- qt(p = mean\_change\_pval\_exp / 2, df = n\_exp - 1, lower.tail = FALSE)}
#' \deqn{t\_nexp <- qt(p = mean\_change\_pval\_nexp / 2, df = n\_nexp - 1, lower.tail = FALSE)}
#' \deqn{mean\_change\_se\_exp <- |\frac{mean\_change\_exp}{t\_exp}|}
#' \deqn{mean\_change\_se\_nexp <- |\frac{mean\_change\_nexp}{t\_nexp}|}

#' Then, this function simply internally calls the \code{\link{es_from_means_se_pre_post}} function but setting:
#' \deqn{mean\_pre\_exp = mean\_change\_exp}
#' \deqn{mean\_pre\_se\_exp = mean\_change\_se\_exp}
#' \deqn{mean\_exp = 0}
#' \deqn{mean\_se\_exp = 0}
#' \deqn{mean\_pre\_nexp = mean\_change\_nexp}
#' \deqn{mean\_pre\_se\_nexp = mean\_change\_se\_nexp}
#' \deqn{mean\_nexp = 0}
#' \deqn{mean\_se\_nexp = 0}
#'
#' To know more about other calculations, see \code{\link{es_from_means_sd_pre_post}} function.
#'
#' @references
#' Bonett, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364–386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
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
                                   smd_to_cor = "viechtbauer", reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.5, length(mean_change_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.5
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.5, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.5

  tryCatch({
    .validate_positive(n_exp, n_nexp, mean_change_pval_exp, mean_change_pval_nexp,
                       error_message = paste0("The number of people exposed/non-exposed and p-values ",
                                              "should be >0."),
                       func = "es_from_mean_change_pval")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  t_exp <- qt(p = mean_change_pval_exp / 2, df = n_exp - 1, lower.tail = FALSE)
  t_nexp <- qt(p = mean_change_pval_nexp / 2, df = n_nexp - 1, lower.tail = FALSE)

  mean_change_se_exp <- abs(mean_change_exp / t_exp)
  mean_change_se_nexp <- abs(mean_change_nexp / t_nexp)

  es <- es_from_means_se_pre_post(
    mean_pre_exp = mean_change_exp,
    mean_exp = 0,
    mean_pre_se_exp = mean_change_se_exp,
    mean_se_exp = 0,
    mean_pre_nexp = mean_change_nexp,
    mean_nexp = 0,
    mean_pre_se_nexp = mean_change_se_nexp,
    mean_se_nexp = 0,
    n_exp = n_exp,
    n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = "cooper",
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_pval"

  return(es)
}
