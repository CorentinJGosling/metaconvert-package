#' Convert means and standard deviations of two independent groups obtained from an ANCOVA model into several effect size measures
#'
#' @param ancova_mean_exp adjusted mean of participants in the experimental/exposed group.
#' @param ancova_mean_nexp adjusted mean of participants in the non-experimental/non-exposed group.
#' @param ancova_mean_sd_exp adjusted standard deviation of participants in the experimental/exposed group.
#' @param ancova_mean_sd_nexp adjusted standard deviation of participants in the non-experimental/non-exposed group.
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_ancova_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first computes an "adjusted" mean difference (MD),
#' Cohen's d (D) and Hedges' g (G) from the adjusted means and standard deviations.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' This function start by estimating the non-adjusted standard deviation of the two groups (formula 12.24 in Cooper);
#' \deqn{mean\_sd\_exp = \frac{ancova\_mean\_sd\_exp}{\sqrt{1 - cov\_outcome\_r^2}}}
#' \deqn{mean\_sd\_nexp = \frac{ancova\_mean\_sd\_nexp}{\sqrt{1 - cov\_outcome\_r^2}}}
#'
#' **To obtain the mean difference**, the following formulas are used (authors calculations):
#' \deqn{md = ancova\_mean\_exp - ancova\_mean\_nexp}
#' \deqn{md\_se = \sqrt{\frac{mean\_sd\_exp^2}{n\_exp} + \frac{mean\_sd\_nexp^2}{n\_nexp}}}
#' \deqn{md\_ci\_lo = md - md\_se * qt(.975, n\_exp+n\_nexp-2-n\_cov\_ancova)}
#' \deqn{md\_ci\_up = md + md\_se * qt(.975, n\_exp+n\_nexp-2-n\_cov\_ancova)}
#'
#' **To obtain the Cohen's d**, the following formulas are used (table 12.3 in Cooper):
#' \deqn{mean\_sd\_pooled = \sqrt{\frac{(n\_exp - 1) * ancova\_mean\_exp^2 + (n\_nexp - 1) * ancova\_mean\_nexp^2}{n\_exp+n\_nexp-2}}}
#' \deqn{cohen\_d =  \frac{ancova\_mean\_exp - ancova\_mean\_nexp}{mean\_sd\_pooled}}
#' \deqn{cohen\_d\_se = \frac{(n\_exp+n\_nexp)*(1-cov\_outcome\_r^2)}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2(n\_exp+n\_nexp)}}
#' \deqn{cohen\_d\_ci\_lo = cohen\_d - cohen\_d\_se * qt(.975, n\_exp + n\_nexp - 2 - n\_cov\_ancova)}
#' \deqn{cohen\_d\_ci\_up = cohen\_d + cohen\_d\_se * qt(.975, n\_exp + n\_nexp - 2 - n\_cov\_ancova)}
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
#'  \code{required input data} \tab See 'Section 19. Adjusted: Means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_ancova_means_sd
#'
#' @md
#'
#' @examples
#' es_from_ancova_means_sd(
#'   n_exp = 55, n_nexp = 55,
#'   ancova_mean_exp = 2.3, ancova_mean_sd_exp = 1.2,
#'   ancova_mean_nexp = 1.9, ancova_mean_sd_nexp = 0.9,
#'   cov_outcome_r = 0.2, n_cov_ancova = 3
#' )
#'
es_from_ancova_means_sd <- function(n_exp, n_nexp,
                                    ancova_mean_exp, ancova_mean_nexp,
                                    ancova_mean_sd_exp, ancova_mean_sd_nexp,
                                    cov_outcome_r, n_cov_ancova,
                                    smd_to_cor = "viechtbauer", reverse_ancova_means) {
  if (missing(reverse_ancova_means)) reverse_ancova_means <- rep(FALSE, length(ancova_mean_exp))
  reverse_ancova_means[is.na(reverse_ancova_means)] <- FALSE
  if (length(reverse_ancova_means) == 1) reverse_ancova_means = c(rep(reverse_ancova_means, length(ancova_mean_exp)))
  if (length(reverse_ancova_means) != length(ancova_mean_exp)) stop("The length of the 'reverse_ancova_means' argument of incorrectly specified.")

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       ancova_mean_sd_exp, ancova_mean_sd_nexp,
                       cov_outcome_r, n_cov_ancova,
                       error_message = paste0("The number of people exposed/non-exposed, ANCOVA SDs ",
                                              "as well as the correlation and number of covariates in ANCOVA ",
                                              "should be >0."),
                       func = "es_from_ancova_means_sd")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  df <- n_exp + n_nexp - 2 - n_cov_ancova
  ancova_mean_sd_pooled <- sqrt(((n_exp - 1) * ancova_mean_sd_exp^2 +
    (n_nexp - 1) * ancova_mean_sd_nexp^2) / (n_exp + n_nexp - 2))

  es <- es_from_ancova_means_sd_pooled_adj(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_sd_pooled = ancova_mean_sd_pooled, cov_outcome_r = cov_outcome_r,
    n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  )


  es$md <- ifelse(reverse_ancova_means,
                  ancova_mean_nexp - ancova_mean_exp,
                  ancova_mean_exp - ancova_mean_nexp)
  mean_sd_exp <- ancova_mean_sd_exp / sqrt(1 - cov_outcome_r^2)
  mean_sd_nexp <- ancova_mean_sd_nexp / sqrt(1 - cov_outcome_r^2)
  es$md_se <- sqrt(mean_sd_exp^2/n_exp + mean_sd_nexp^2/n_nexp)
  es$md_ci_lo <- es$md - qt(.975, df) * es$md_se
  es$md_ci_up <- es$md + qt(.975, df) * es$md_se

  es$info_used <- "ancova_means_sd"

  return(es)
}
#' Convert means and standard errors of two independent groups obtained from an ANCOVA model into several effect size measures
#'
#' @param ancova_mean_exp adjusted mean of participants in the experimental/exposed group.
#' @param ancova_mean_nexp adjusted mean of participants in the non-experimental/non-exposed group.
#' @param ancova_mean_se_exp adjusted standard error of participants in the experimental/exposed group.
#' @param ancova_mean_se_nexp adjusted standard error of participants in the non-experimental/non-exposed group.
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_ancova_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the adjusted means standard errors of two independent groups
#' into standard deviations,
#' and then relies on the calculations of the \code{\link{es_from_ancova_means_sd}} function.
#'
#' **To convert the standard errors into standard deviations**, the following formula is used.
#' \deqn{ancova\_mean\_sd\_exp = ancova\_mean\_se\_exp * \sqrt{n\_exp}}
#' \deqn{ancova\_mean\_sd\_nexp = ancova\_mean\_se\_nexp * \sqrt{n\_nexp}}
#' Calculations of the \code{\link{es_from_ancova_means_sd}()} are then applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 19. Adjusted: Means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @export es_from_ancova_means_se
#'
#' @md
#'
#' @examples
#' es_from_ancova_means_se(
#'   n_exp = 55, n_nexp = 55,
#'   ancova_mean_exp = 2.3, ancova_mean_se_exp = 1.2,
#'   ancova_mean_nexp = 1.9, ancova_mean_se_nexp = 0.9,
#'   cov_outcome_r = 0.2, n_cov_ancova = 3
#' )
es_from_ancova_means_se <- function(n_exp, n_nexp,
                                    ancova_mean_exp, ancova_mean_nexp,
                                    ancova_mean_se_exp, ancova_mean_se_nexp,
                                    cov_outcome_r, n_cov_ancova,
                                    smd_to_cor = "viechtbauer", reverse_ancova_means) {
  if (missing(reverse_ancova_means)) reverse_ancova_means <- rep(FALSE, length(ancova_mean_exp))
  reverse_ancova_means[is.na(reverse_ancova_means)] <- FALSE


  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       ancova_mean_se_exp, ancova_mean_se_nexp,
                       cov_outcome_r, n_cov_ancova,
                       error_message = paste0("The number of people exposed/non-exposed, ANCOVA SEs ",
                                              "as well as the correlation and number of covariates in ANCOVA ",
                                              "should be >0."),
                       func = "es_from_ancova_means_se")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })


  ancova_mean_sd_exp <- ancova_mean_se_exp * sqrt(n_exp)
  ancova_mean_sd_nexp <- ancova_mean_se_nexp * sqrt(n_nexp)


  es <- es_from_ancova_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_sd_exp = ancova_mean_sd_exp, ancova_mean_sd_nexp = ancova_mean_sd_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  )

  es$info_used <- "ancova_means_se"

  return(es)
}
#' Convert means and 95% CIs of two independent groups obtained from an ANCOVA model into several effect size measures
#'
#' @param ancova_mean_exp adjusted mean of participants in the experimental/exposed group.
#' @param ancova_mean_nexp adjusted mean of participants in the non-experimental/non-exposed group.
#' @param ancova_mean_ci_lo_exp lower bound of the adjusted 95% CI of the mean of the experimental/exposed group
#' @param ancova_mean_ci_up_exp upper bound of the adjusted 95% CI of the mean of the experimental/exposed group
#' @param ancova_mean_ci_lo_nexp lower bound of the adjusted 95% CI of the mean of the non-experimental/non-exposed group.
#' @param ancova_mean_ci_up_nexp upper bound of the adjusted 95% CI of the mean of the non-experimental/non-exposed group.
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param reverse_ancova_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the adjusted means 95% CI of two independent groups
#' into a standard error,
#' and then relies on the calculations of the \code{\link{es_from_ancova_means_se}()} function.
#'
#' **To convert the 95% CIs into standard errors,** the following formula is used (table 12.3 in Cooper):
#' \deqn{ancova\_mean\_se\_exp = \frac{ancova\_mean\_ci\_up\_exp - ancova\_mean\_ci\_lo\_exp}{2 * qt(0.975, df = n\_exp - 1)}}
#' \deqn{ancova\_mean\_se\_nexp = \frac{ancova\_mean\_ci\_up\_nexp - ancova\_mean\_ci\_lo\_nexp}{2 * qt(0.975, df = n\_nexp - 1)}}
#' Calculations of the \code{\link{es_from_ancova_means_se}()} are then applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 19. Adjusted: Means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_ancova_means_ci
#'
#' @md
#'
#' @examples
#' es_from_ancova_means_ci(
#'   n_exp = 55, n_nexp = 55, cov_outcome_r = 0.5, n_cov_ancova = 4,
#'   ancova_mean_exp = 25, ancova_mean_ci_lo_exp = 15, ancova_mean_ci_up_exp = 35,
#'   ancova_mean_nexp = 18, ancova_mean_ci_lo_nexp = 12, ancova_mean_ci_up_nexp = 24
#' )
es_from_ancova_means_ci <- function(n_exp, n_nexp,
                                    ancova_mean_exp, ancova_mean_ci_lo_exp, ancova_mean_ci_up_exp,
                                    ancova_mean_nexp, ancova_mean_ci_lo_nexp, ancova_mean_ci_up_nexp,
                                    cov_outcome_r, n_cov_ancova,
                                    max_asymmetry = 10,
                                    smd_to_cor = "viechtbauer", reverse_ancova_means) {
  if (missing(reverse_ancova_means)) reverse_ancova_means <- rep(FALSE, length(ancova_mean_exp))
  reverse_ancova_means[is.na(reverse_ancova_means)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       cov_outcome_r, n_cov_ancova,
                       error_message = paste0("The number of people exposed/non-exposed ",
                                              "as well as the correlation and number of covariates in ANCOVA ",
                                              "should be >0."),
                       func = "es_from_ancova_means_ci")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })
  tryCatch({
    .validate_ci_symmetry(ancova_mean_exp, ancova_mean_ci_lo_exp, ancova_mean_ci_up_exp,
                          func = "es_from_ancova_means_ci",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })
  tryCatch({
    .validate_ci_symmetry(ancova_mean_nexp, ancova_mean_ci_lo_nexp, ancova_mean_ci_up_nexp,
                          func = "es_from_ancova_means_ci",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })

  df_exp <- n_exp - 1
  df_nexp <- n_nexp - 1

  ancova_se_exp <- (ancova_mean_ci_up_exp - ancova_mean_ci_lo_exp) / (2 * qt(0.975, df_exp))
  ancova_se_nexp <- (ancova_mean_ci_up_nexp - ancova_mean_ci_lo_nexp) / (2 * qt(0.975, df_nexp))

  es <- es_from_ancova_means_se(
    n_exp = n_exp, n_nexp = n_nexp,
    ancova_mean_exp = ancova_mean_exp,
    ancova_mean_se_exp = ancova_se_exp,
    ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_se_nexp = ancova_se_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  )

  es$info_used <- "ancova_means_ci"

  return(es)
}
#' Convert means and adjusted pooled standard deviation of two independent groups obtained from an ANCOVA model into several effect size measures
#'
#' @param ancova_mean_exp adjusted mean of participants in the experimental/exposed group.
#' @param ancova_mean_nexp adjusted mean of participants in the non-experimental/non-exposed group.
#' @param ancova_mean_sd_pooled adjusted pooled standard deviation.
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_ancova_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the adjusted pooled standard deviations of two independent groups
#' into a crude pooled standard deviation.
#' and then relies on the calculations of the \code{\link{es_from_ancova_means_sd_pooled_crude}()} function.
#'
#' **To convert the adjusted pooled SD into a crude pooled SD** (table 12.3 in Cooper):
#' \deqn{mean\_sd\_pooled =  \frac{ancova\_mean\_sd\_pooled}{\sqrt{1 - cov\_outcome\_r^2}}}
#'
#' Calculations of the \code{\link{es_from_ancova_means_sd_pooled_crude}()} are then applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 19. Adjusted: Means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_ancova_means_sd_pooled_adj
#'
#' @md
#'
#' @examples
#' es_from_ancova_means_sd_pooled_adj(
#'   ancova_mean_exp = 98, ancova_mean_nexp = 87,
#'   ancova_mean_sd_pooled = 17, cov_outcome_r = 0.2,
#'   n_cov_ancova = 3, n_exp = 20, n_nexp = 20
#' )
es_from_ancova_means_sd_pooled_adj <- function(ancova_mean_exp, ancova_mean_nexp,
                                               ancova_mean_sd_pooled, cov_outcome_r,
                                               n_cov_ancova, n_exp, n_nexp,
                                               smd_to_cor = "viechtbauer", reverse_ancova_means) {
  if (missing(reverse_ancova_means)) reverse_ancova_means <- rep(FALSE, length(ancova_mean_exp))
  reverse_ancova_means[is.na(reverse_ancova_means)] <- FALSE
  if (length(reverse_ancova_means) == 1) reverse_ancova_means = c(rep(reverse_ancova_means, length(ancova_mean_exp)))
  if (length(reverse_ancova_means) != length(ancova_mean_exp)) stop("The length of the 'reverse_ancova_means' argument of incorrectly specified.")

  tryCatch({
    .validate_positive(n_exp, n_nexp, ancova_mean_sd_pooled,
                       cov_outcome_r, n_cov_ancova,
                       error_message = paste0("The number of people exposed/non-exposed, ANCOVA pooled SD ",
                                              "as well as the correlation and number of covariates in ANCOVA ",
                                              "should be >0."),
                       func = "es_from_ancova_means_sd_pooled_adj")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  sd_within <- ancova_mean_sd_pooled / sqrt(1 - cov_outcome_r^2)

  d <- (ancova_mean_exp - ancova_mean_nexp) / sd_within

  es <- .es_from_d(
    d = d, cov_outcome_r = cov_outcome_r,
    adjusted = TRUE,
    n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_ancova_means
  )

  es$info_used <- "ancova_means_sd_pooled_adj"

  es$md <- ifelse(reverse_ancova_means,
    ancova_mean_nexp - ancova_mean_exp,
    ancova_mean_exp - ancova_mean_nexp
  )
  es$md_se <- sqrt((n_exp + n_nexp) / (n_exp * n_nexp) * (1 - cov_outcome_r^2) * sd_within^2)
  es$md_ci_lo <- es$md - qt(.975, n_exp + n_nexp - 2 - n_cov_ancova) * es$md_se
  es$md_ci_up <- es$md + qt(.975, n_exp + n_nexp - 2 - n_cov_ancova) * es$md_se

  return(es)
}

#' Convert adjusted means obtained from an ANCOVA model and crude pooled standard deviation of two independent groups into several effect size measures
#'
#' @param ancova_mean_exp adjusted mean of participants in the experimental/exposed group.
#' @param ancova_mean_nexp adjusted mean of participants in the non-experimental/non-exposed group.
#' @param mean_sd_pooled crude pooled standard deviation.
#' @param cov_outcome_r correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the adjusted \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_ancova_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first computes an "adjusted" mean difference (MD) and Cohen's d (D)
#' from the adjusted means and crude pooled standard deviation of two independent groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the Cohen's d**:
#' \deqn{d =  \frac{ancova\_mean\_exp - ancova\_mean\_nexp\_adj}{mean\_sd\_pooled}}
#'
#' **To estimate the mean difference**:
#' \deqn{md =  ancova\_mean\_exp - ancova\_mean\_nexp\_adj}
#' \deqn{md\_se =  \sqrt{\frac{n\_exp + n\_nexp}{n\_exp * n\_nexp} * (1 - cov\_outcome\_r^2) * mean\_sd\_pooled^2}}
#'
#' Then, calculations of the \code{\link{es_from_ancova_means_sd}()} and \code{\link{es_from_cohen_d_adj}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 19. Adjusted: Means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_ancova_means_sd_pooled_crude
#'
#' @md
#'
#' @examples
#' es_from_ancova_means_sd_pooled_crude(
#'   ancova_mean_exp = 29, ancova_mean_nexp = 34,
#'   mean_sd_pooled = 7, cov_outcome_r = 0.2,
#'   n_cov_ancova = 3, n_exp = 20, n_nexp = 20
#' )
es_from_ancova_means_sd_pooled_crude <- function(ancova_mean_exp, ancova_mean_nexp,
                                                 mean_sd_pooled, cov_outcome_r, n_cov_ancova,
                                                 n_exp, n_nexp, smd_to_cor = "viechtbauer", reverse_ancova_means) {
  if (missing(reverse_ancova_means)) reverse_ancova_means <- rep(FALSE, length(ancova_mean_exp))
  reverse_ancova_means[is.na(reverse_ancova_means)] <- FALSE
  if (length(reverse_ancova_means) == 1) reverse_ancova_means = c(rep(reverse_ancova_means, length(ancova_mean_exp)))
  if (length(reverse_ancova_means) != length(ancova_mean_exp)) stop("The length of the 'reverse_ancova_means' argument of incorrectly specified.")

  tryCatch({
    .validate_positive(n_exp, n_nexp, mean_sd_pooled,
                       cov_outcome_r, n_cov_ancova,
                       error_message = paste0("The number of people exposed/non-exposed, pooled SD ",
                                              "as well as the correlation and number of covariates in ANCOVA ",
                                              "should be >0."),
                       func = "es_from_ancova_means_sd_pooled_crude")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  d <- (ancova_mean_exp - ancova_mean_nexp) / mean_sd_pooled

  es <- .es_from_d(
    d = d, cov_outcome_r = cov_outcome_r,
    adjusted = TRUE,
    n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_ancova_means
  )

  es$info_used <- "ancova_means_sd_pooled"

  es$md <- ifelse(reverse_ancova_means,
    ancova_mean_nexp - ancova_mean_exp,
    ancova_mean_exp - ancova_mean_nexp
  )

  es$md_se <- sqrt((n_exp + n_nexp) / (n_exp * n_nexp) * (1 - cov_outcome_r^2) * mean_sd_pooled^2)
  es$md_ci_lo <- es$md - qt(.975, n_exp + n_nexp - 2 - n_cov_ancova) * es$md_se
  es$md_ci_up <- es$md + qt(.975, n_exp + n_nexp - 2 - n_cov_ancova) * es$md_se

  return(es)
}
