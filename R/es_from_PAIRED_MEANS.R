#' Convert pre-post means of two independent groups into various effect size measures
#'
#' @param mean_pre_exp mean of the experimental/exposed group at baseline
#' @param mean_exp mean of the experimental/exposed group at follow up
#' @param mean_pre_sd_exp standard deviation of the experimental/exposed group at baseline
#' @param mean_sd_exp standard deviation of the experimental/exposed group at follow up
#' @param mean_pre_nexp mean of the non-experimental/non-exposed group at baseline
#' @param mean_nexp mean of the non-experimental/non-exposed group at follow up
#' @param mean_pre_sd_nexp standard deviation of the non-experimental/non-exposed group at baseline
#' @param mean_sd_nexp standard deviation of the non-experimental/non-exposed group at follow up
#' @param n_exp number of the experimental/exposed group
#' @param n_nexp number of the non-experimental/non-exposed group
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param pre_post_to_smd formula used to convert the pre and post means/SD into a SMD (see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_means_pre_post a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts pre-post means of two independent groups into a Cohen's d (D) and Hedges' g (G).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **Two approaches can be used to compute the Cohen's d.**
#'
#' In these two approaches, the standard deviation of the difference within each group first needs to be obtained:
#' \deqn{adj\_exp = 2*r\_pre\_post\_exp*mean\_pre\_sd\_exp*mean\_sd\_exp}
#' \deqn{sd\_change\_exp = \sqrt{mean\_pre\_sd\_exp^2 + mean\_sd\_exp^2 - adj\_exp}}
#' \deqn{adj\_nexp = 2*r\_pre\_post\_nexp*mean\_pre\_sd\_nexp*mean\_sd\_nexp}
#' \deqn{sd\_change\_nexp = \sqrt{mean\_pre\_sd\_nexp^2 + mean\_sd\_nexp^2 - adj\_nexp}}
#'
#' 1. In the approach described by Bonett (\code{pre_post_to_smd = "bonett"}), one Cohen's d per group is obtained by
#' standardizing the pre-post mean difference by the standard deviation at baseline (Bonett, 2008):
#' \deqn{cohen\_d\_exp = \frac{mean\_pre\_exp - mean\_exp}{mean\_pre\_sd\_exp}}
#' \deqn{cohen\_d\_nexp = \frac{mean\_pre\_nexp - mean\_nexp}{mean\_pre\_sd\_nexp}}
#' \deqn{cohen\_d\_se\_exp = \sqrt{\frac{sd\_change\_exp^2}{mean\_pre\_sd\_exp^2 * (n\_exp - 1) + g\_exp^2 / (2 * (n\_exp - 1))}}}
#' \deqn{cohen\_d\_se\_nexp = \sqrt{\frac{sd\_change\_nexp^2}{mean\_pre\_sd\_nexp^2 * (n\_nexp - 1) + g\_nexp^2 / (2 * (n\_nexp - 1))}}}
#'
#' 2. In the approach described by Cooper (\code{pre_post_to_smd = "cooper"}), the following formulas are used:
#' \deqn{cohen\_d\_exp = \frac{mean\_pre\_exp - mean\_exp}{sd\_change\_exp} * \sqrt{2 * (1 - r\_pre\_post\_exp)}}
#' \deqn{cohen\_d\_nexp = \frac{mean\_pre\_nexp - mean\_nexp}{sd\_change\_nexp} * \sqrt{2 * (1 - r\_pre\_post\_nexp)}}
#' \deqn{cohen\_d\_se\_exp = \frac{2 * (1 - r\_pre\_post\_exp)}{n\_exp} + \frac{cohen\_d\_exp^2}{2 * n\_exp}}
#' \deqn{cohen\_d\_se\_nexp = \frac{2 * (1 - r\_pre\_post\_nexp)}{n\_nexp} + \frac{cohen\_d\_nexp^2}{2 * n\_nexp}}
#'
#' Last, the Cohen's d reflecting the within-group change from baseline to follow-up are combined into one Cohen's d:
#' \deqn{cohen\_d = d\_exp - d\_nexp}
#' \deqn{cohen\_d\_se = \sqrt{cohen\_d\_se\_exp^2 + cohen\_d\_se\_nexp^2}}
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
#'  \code{required input data} \tab See 'Section 15. Paired: pre-post means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Bonett, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364–386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_means_sd_pre_post
#'
#' @md
#'
#' @examples
#' es_from_means_sd_pre_post(
#'   n_exp = 36, n_nexp = 35,
#'   mean_pre_exp = 98, mean_exp = 102,
#'   mean_pre_sd_exp = 16, mean_sd_exp = 17,
#'   mean_pre_nexp = 96, mean_nexp = 102,
#'   mean_pre_sd_nexp = 14, mean_sd_nexp = 15,
#'   r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8
#' )
es_from_means_sd_pre_post <- function(mean_pre_exp, mean_exp, mean_pre_sd_exp, mean_sd_exp,
                                      mean_pre_nexp, mean_nexp, mean_pre_sd_nexp, mean_sd_nexp,
                                      n_exp, n_nexp, r_pre_post_exp, r_pre_post_nexp,
                                      smd_to_cor = "viechtbauer",
                                      pre_post_to_smd = "bonett",
                                      reverse_means_pre_post) {
  if (!all(pre_post_to_smd %in% c("bonett", "cooper"))) {
    stop(paste0(
      "'", paste(unique(pre_post_to_smd), sep=" / "), "' not in tolerated values for the 'pre_post_to_smd' argument.",
      "Possible inputs are: 'bonett' or 'cooper'"
    ))
  }

  tryCatch({
    .validate_positive(n_exp, n_nexp, mean_pre_sd_exp, mean_pre_sd_nexp, mean_sd_exp, mean_sd_nexp,
                       error_message = paste0("The number of people exposed/non-exposed and standard deviations ",
                                              "should be >0."),
                       func = "es_from_means_sd_pre_post")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  if (missing(reverse_means_pre_post)) reverse_means_pre_post <- rep(FALSE, length(mean_pre_exp))
  reverse_means_pre_post[is.na(reverse_means_pre_post)] <- FALSE
  if (length(reverse_means_pre_post) == 1) reverse_means_pre_post = c(rep(reverse_means_pre_post, length(mean_pre_exp)))
  if (length(reverse_means_pre_post) != length(mean_pre_exp)) stop("The length of the 'reverse_means_pre_post' argument of incorrectly specified.")

  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.5, length(mean_pre_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.5
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.5, length(mean_pre_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.5

  g <- g_se <- g_ci_lo <- g_ci_up <-
    d <- d_se <- d_ci_lo <- d_ci_up <- rep(NA, length(mean_pre_exp))

  nn_miss <- which(!is.na(mean_pre_exp) & !is.na(mean_pre_sd_exp) &
    !is.na(mean_pre_nexp) & !is.na(mean_pre_sd_nexp) &
    !is.na(mean_exp) & !is.na(mean_sd_exp) &
    !is.na(mean_nexp) & !is.na(mean_sd_nexp) &
    !is.na(r_pre_post_exp) & !is.na(r_pre_post_nexp) &
    !is.na(n_exp) & !is.na(n_nexp))

  dat_smd_pre_post <- data.frame(
    mean_pre_exp = mean_pre_exp, mean_pre_sd_exp = mean_pre_sd_exp,
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp,
    mean_pre_nexp = mean_pre_nexp, mean_pre_sd_nexp = mean_pre_sd_nexp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,
    pre_post_to_smd = pre_post_to_smd
  )
  if (length(nn_miss) != 0) {
    smd_pp <- t(mapply(.pre_post_to_smd,
      mean_pre_exp = dat_smd_pre_post$mean_pre_exp[nn_miss],
      mean_pre_sd_exp = dat_smd_pre_post$mean_pre_sd_exp[nn_miss],
      mean_exp = dat_smd_pre_post$mean_exp[nn_miss],
      mean_sd_exp = dat_smd_pre_post$mean_sd_exp[nn_miss],
      mean_pre_nexp = dat_smd_pre_post$mean_pre_nexp[nn_miss],
      mean_pre_sd_nexp = dat_smd_pre_post$mean_pre_sd_nexp[nn_miss],
      mean_nexp = dat_smd_pre_post$mean_nexp[nn_miss],
      mean_sd_nexp = dat_smd_pre_post$mean_sd_nexp[nn_miss],
      n_exp = dat_smd_pre_post$n_exp[nn_miss],
      n_nexp = dat_smd_pre_post$n_nexp[nn_miss],
      r_pre_post_exp = dat_smd_pre_post$r_pre_post_exp[nn_miss],
      r_pre_post_nexp = dat_smd_pre_post$r_pre_post_nexp[nn_miss],
      pre_post_to_smd = dat_smd_pre_post$pre_post_to_smd[nn_miss]
    ))

    d[nn_miss] <- smd_pp[, 1]
    d_se[nn_miss] <- sqrt(smd_pp[, 2])
    d_ci_lo[nn_miss] <- smd_pp[, 3]
    d_ci_up[nn_miss] <- smd_pp[, 4]
    g[nn_miss] <- smd_pp[, 5]
    g_se[nn_miss] <- sqrt(smd_pp[, 6])
    g_ci_lo[nn_miss] <- smd_pp[, 7]
    g_ci_up[nn_miss] <- smd_pp[, 8]
  }

  es <- .es_from_d(
    d = d, d_se = d_se, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_means_pre_post
  )

  es$g <- ifelse(reverse_means_pre_post, -g, g)
  es$g_se <- g_se
  es$g_ci_lo <- g_ci_lo
  es$g_ci_up <- g_ci_up

  md_exp <- mean_pre_exp - mean_exp
  md_var_exp <- (mean_pre_sd_exp^2 + mean_sd_exp^2 -
    2 * r_pre_post_exp * mean_pre_sd_exp * mean_sd_exp) / n_exp

  md_nexp <- mean_pre_nexp - mean_nexp
  md_var_nexp <- (mean_pre_sd_nexp^2 + mean_sd_nexp^2 -
    2 * r_pre_post_nexp * mean_pre_sd_nexp * mean_sd_nexp) / n_nexp

  es$md <- ifelse(reverse_means_pre_post, md_nexp - md_exp, md_exp - md_nexp)
  es$md_se <- sqrt(md_var_exp + md_var_nexp)
  es$md_ci_lo <- es$md - es$md_se * qt(.975, n_exp + n_nexp - 2)
  es$md_ci_up <- es$md + es$md_se * qt(.975, n_exp + n_nexp - 2)

  es$info_used <- "means_sd_pre_post"

  return(es)
}

#' Convert pre-post means of two independent groups into various effect size measures
#'
#' @param mean_pre_exp mean of the experimental/exposed group at baseline
#' @param mean_pre_se_exp standard error of the experimental/exposed group at baseline
#' @param mean_exp mean of the experimental/exposed group at follow up
#' @param mean_se_exp standard error of the experimental/exposed group at follow up
#' @param mean_pre_nexp mean of the non-experimental/non-exposed group at baseline
#' @param mean_pre_se_nexp standard error of the non-experimental/non-exposed group at baseline
#' @param mean_nexp mean of the non-experimental/non-exposed group at follow up
#' @param mean_se_nexp standard error of the non-experimental/non-exposed group at follow up
#' @param n_exp number of the experimental/exposed group
#' @param n_nexp number of the non-experimental/non-exposed group
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param pre_post_to_smd formula used to convert the pre and post means/SD into a SMD (see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_means_pre_post a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the pre/post standard errors of two independent groups into standard deviations (Section 6.5.2.2 in the Cochrane Handbook).
#' \deqn{mean\_pre\_sd\_exp = mean\_pre\_se\_exp * \sqrt{n\_exp}}
#' \deqn{mean\_pre\_sd\_nexp = mean\_pre\_se\_nexp * \sqrt{n\_nexp}}
#' \deqn{mean\_sd\_exp = mean\_se\_exp * \sqrt{n\_exp}}
#' \deqn{mean\_sd\_nexp = mean\_se\_nexp * \sqrt{n\_nexp}}
#'
#' Then, calculations of the \code{\link{es_from_means_sd_pre_post}()} are applied.
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
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 15. Paired: pre-post means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_means_se_pre_post
#'
#' @md
#'
#' @examples
#' es_from_means_sd_pre_post(
#'   n_exp = 36, n_nexp = 35,
#'   mean_pre_exp = 98, mean_exp = 102,
#'   mean_pre_sd_exp = 16, mean_sd_exp = 17,
#'   mean_pre_nexp = 96, mean_nexp = 102,
#'   mean_pre_sd_nexp = 14, mean_sd_nexp = 15,
#'   r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8
#' )
es_from_means_se_pre_post <- function(mean_pre_exp, mean_exp, mean_pre_se_exp, mean_se_exp,
                                      mean_pre_nexp, mean_nexp, mean_pre_se_nexp, mean_se_nexp,
                                      n_exp, n_nexp, r_pre_post_exp, r_pre_post_nexp,
                                      smd_to_cor = "viechtbauer",
                                      pre_post_to_smd = "bonett",
                                      reverse_means_pre_post) {
  if (!all(pre_post_to_smd %in% c("bonett", "cooper"))) {
    stop(paste0(
      "'",
      unique(pre_post_to_smd[!pre_post_to_smd %in% c("bonett", "cooper")]),
      "' not in tolerated values for the 'pre_post_to_smd' argument.",
      "Possible inputs are: 'bonett' or 'cooper'"
    ))
  }

  tryCatch({
    .validate_positive(n_exp, n_nexp, mean_pre_se_exp, mean_pre_se_nexp, mean_se_exp, mean_se_nexp,
                       error_message = paste0("The number of people exposed/non-exposed and standard errors ",
                                              "should be >0."),
                       func = "es_from_means_se_pre_post")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })


  if (missing(reverse_means_pre_post)) reverse_means_pre_post <- rep(FALSE, length(mean_pre_exp))
  reverse_means_pre_post[is.na(reverse_means_pre_post)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.5, length(mean_pre_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.5
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.5, length(mean_pre_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.5

  sd_pre_exp <- mean_pre_se_exp * sqrt(n_exp)
  sd_exp <- mean_se_exp * sqrt(n_exp)
  sd_pre_nexp <- mean_pre_se_nexp * sqrt(n_nexp)
  sd_nexp <- mean_se_nexp * sqrt(n_nexp)

  es <- es_from_means_sd_pre_post(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_sd_exp = sd_pre_exp, mean_sd_exp = sd_exp,
    mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_nexp,
    mean_pre_sd_nexp = sd_pre_nexp, mean_sd_nexp = sd_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_means_pre_post = reverse_means_pre_post
  )

  es$info_used <- "means_se_pre_post"

  return(es)
}

#' Convert pre-post means of two independent groups into various effect size measures
#'
#' @param mean_pre_exp mean of the experimental/exposed group at baseline
#' @param mean_exp mean of the experimental/exposed group at follow up
#' @param mean_pre_nexp mean of the non-experimental/non-exposed group at baseline
#' @param mean_nexp mean of the non-experimental/non-exposed group at follow up
#' @param mean_pre_ci_lo_exp lower bound of the 95% CI of the mean of the experimental/exposed group at baseline
#' @param mean_pre_ci_up_exp upper bound of the 95% CI of the mean of the experimental/exposed group at baseline
#' @param mean_ci_lo_exp lower bound of the 95% CI of the mean of the experimental/exposed group at follow up
#' @param mean_ci_up_exp upper bound of the 95% CI of the mean of the experimental/exposed group at follow up
#' @param mean_pre_ci_lo_nexp lower bound of the 95% CI of the mean of the non-experimental/non-exposed group at baseline
#' @param mean_pre_ci_up_nexp upper bound of the 95% CI of the mean of the non-experimental/non-exposed group at baseline
#' @param mean_ci_lo_nexp lower bound of the 95% CI of the mean of the non-experimental/non-exposed group at follow up
#' @param mean_ci_up_nexp upper bound of the 95% CI of the mean of the non-experimental/non-exposed group at follow up
#' @param n_exp number of the experimental/exposed group
#' @param n_nexp number of the non-experimental/non-exposed group
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param pre_post_to_smd formula used to convert the pre and post means/SD into a SMD (see details).
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_means_pre_post a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the bounds of the 95% CI of the pre/post means of two independent groups into standard errors (Section 6.3.1 in the Cochrane Handbook).
#' \deqn{mean\_pre\_se\_exp = \frac{mean\_pre\_ci\_up\_exp - mean\_pre\_ci\_lo\_exp}{2 * qt{(0.975, df = n\_exp - 1)}}}
#' \deqn{mean\_pre\_se\_nexp = \frac{mean\_pre\_ci\_up\_nexp - mean\_pre\_ci\_lo\_nexp}{2 * qt{(0.975, df = n\_nexp - 1)}}}
#' \deqn{mean\_se\_exp = \frac{mean\_ci\_up\_exp - mean\_ci\_lo\_exp}{2 * qt{(0.975, df = n\_exp - 1)}}}
#' \deqn{mean\_se\_nexp = \frac{mean\_ci\_up\_nexp - mean\_ci\_lo\_nexp}{2 * qt{(0.975, df = n\_nexp - 1)}}}
#'
#' Then, calculations of the \code{\link{es_from_means_se_pre_post}} are applied.
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @export es_from_means_ci_pre_post
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 15. Paired: pre-post means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_means_ci_pre_post(
#'   n_exp = 36, n_nexp = 35,
#'   mean_pre_exp = 98,
#'   mean_pre_ci_lo_exp = 88,
#'   mean_pre_ci_up_exp = 108,
#'   mean_exp = 102,
#'   mean_ci_lo_exp = 92,
#'   mean_ci_up_exp = 112,
#'   mean_pre_nexp = 96,
#'   mean_pre_ci_lo_nexp = 86,
#'   mean_pre_ci_up_nexp = 106,
#'   mean_nexp = 102,
#'   mean_ci_lo_nexp = 92,
#'   mean_ci_up_nexp = 112,
#'   r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8
#' )
es_from_means_ci_pre_post <- function(mean_pre_exp, mean_exp,
                                      mean_pre_ci_lo_exp, mean_pre_ci_up_exp,
                                      mean_ci_lo_exp, mean_ci_up_exp,
                                      mean_pre_nexp, mean_nexp,
                                      mean_pre_ci_lo_nexp, mean_pre_ci_up_nexp,
                                      mean_ci_lo_nexp, mean_ci_up_nexp,
                                      n_exp, n_nexp, r_pre_post_exp, r_pre_post_nexp,
                                      smd_to_cor = "viechtbauer",
                                      pre_post_to_smd = "bonett", max_asymmetry = 10,
                                      reverse_means_pre_post) {

  if (missing(reverse_means_pre_post)) reverse_means_pre_post <- rep(FALSE, length(mean_pre_exp))
  reverse_means_pre_post[is.na(reverse_means_pre_post)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.5, length(mean_pre_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.5
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.5, length(mean_pre_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.5

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       error_message = paste0("The number of people exposed/non-exposed ",
                                              "should be >0."),
                       func = "es_from_means_ci_pre_post")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  tryCatch({
    .validate_ci_symmetry(mean_exp, mean_ci_lo_exp, mean_ci_up_exp,
                          func = "es_from_means_ci_pre_post",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })
  tryCatch({
    .validate_ci_symmetry(mean_nexp, mean_ci_lo_nexp, mean_ci_up_nexp,
                          func = "es_from_means_ci_pre_post",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })

  tryCatch({
    .validate_ci_symmetry(mean_pre_exp, mean_pre_ci_lo_exp, mean_pre_ci_up_exp,
                          func = "es_from_means_ci_pre_post",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })
  tryCatch({
    .validate_ci_symmetry(mean_pre_nexp, mean_pre_ci_lo_nexp, mean_pre_ci_up_nexp,
                          func = "es_from_means_ci_pre_post",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })

  df_exp <- n_exp - 1
  df_nexp <- n_nexp - 1

  se_pre_exp <- (mean_pre_ci_up_exp - mean_pre_ci_lo_exp) / (2 * qt(0.975, df_exp))
  se_exp <- (mean_ci_up_exp - mean_ci_lo_exp) / (2 * qt(0.975, df_exp))
  se_pre_nexp <- (mean_pre_ci_up_nexp - mean_pre_ci_lo_nexp) / (2 * qt(0.975, df_nexp))
  se_nexp <- (mean_ci_up_nexp - mean_ci_lo_nexp) / (2 * qt(0.975, df_nexp))

  es <- es_from_means_se_pre_post(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_se_exp = se_pre_exp, mean_se_exp = se_exp,
    mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_nexp,
    mean_pre_se_nexp = se_pre_nexp, mean_se_nexp = se_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_means_pre_post = reverse_means_pre_post
  )

  es$info_used <- "means_ci_pre_post"

  return(es)
}
