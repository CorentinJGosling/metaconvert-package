#' Convert a standardized regression coefficient and the standard deviation of the dependent variable
#' into several effect size measures
#'
#' @param beta_std a standardized regression coefficient value (binary predictor, no other covariables in the model)
#' @param sd_dv standard deviation of the dependent variable
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_beta_std a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts a standardized linear regression coefficient
#' (coming from a model with only one binary predictor), into an
#' unstandardized linear regression coefficient.
#'
#' \deqn{sd\_dummy = \sqrt{\frac{n_exp - (n_exp^2 / (n_exp + n_nexp))}{(n_exp + n_nexp - 1)}}}
#' \deqn{unstd\_beta = beta\_std * \frac{sd\_dv}{sd\_dummy}}
#'
#' Calculations of the \code{\link{es_from_beta_unstd}} functions are then used.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 13. (Un-)Standardized regression coefficient'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' @export es_from_beta_std
#'
#' @md
#'
#' @examples
#' es_from_beta_std(beta_std = 2.1, sd_dv = 0.98, n_exp = 20, n_nexp = 22)
es_from_beta_std <- function(beta_std, sd_dv, n_exp, n_nexp,
                             smd_to_cor = "viechtbauer", reverse_beta_std) {
  if (missing(reverse_beta_std)) reverse_beta_std <- rep(FALSE, length(beta_std))
  reverse_beta_std[is.na(reverse_beta_std)] <- FALSE
  if (length(reverse_beta_std) == 1) reverse_beta_std = c(rep(reverse_beta_std, length(beta_std)))
  if (length(reverse_beta_std) != length(beta_std)) stop("The length of the 'reverse_beta_std' argument is incorrectly specified.")

  beta_std <- ifelse(reverse_beta_std, -beta_std, beta_std)

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       error_message = paste0("The number of people exposed/non-exposed, and standard deviation of the outcome ",
                                              "should be >0."),
                       func = "es_from_beta_std")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })


  sd_dummy <- sqrt((n_exp - (n_exp^2 / (n_exp + n_nexp))) / (n_exp + n_nexp - 1))

  unstd_beta <- beta_std * (sd_dv / sd_dummy)

  es <- es_from_beta_unstd(
    beta_unstd = unstd_beta, sd_dv = sd_dv,
    n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor
  )

  es$info_used <- "beta_std"

  return(es)
}

#' Convert an unstandardized regression coefficient and the standard deviation of the dependent variable
#' into several effect size measures
#'
#' @param beta_unstd an unstandardized regression coefficient value (binary predictor, no other covariables in the model)
#' @param sd_dv standard deviation of the dependent variable
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_beta_unstd a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function estimates a Cohen's d (D) and Hedges' g (G) from an unstandardized linear regression coefficient (coming from a model with only one binary predictor),
#' and the standard deviation of the dependent variable.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **The formula used to obtain the Cohen's d is**:
#' \deqn{N = n\_exp + n\_nexp}
#' \deqn{sd\_pooled = \sqrt{\frac{sd\_dv^2 * (N - 1) - unstd\_beta^2 * \frac{n\_exp * n\_nexp}{N}}{N - 2}}}
#' \deqn{cohen\_d = \frac{unstd\_beta}{sd\_pooled}}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 13. (Un-)Standardized regression coefficient'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' @export es_from_beta_unstd
#'
#' @md
#'
#' @examples
#' es_from_beta_unstd(beta_unstd = 2.1, sd_dv = 0.98, n_exp = 20, n_nexp = 22)
es_from_beta_unstd <- function(beta_unstd, sd_dv, n_exp, n_nexp,
                               smd_to_cor = "viechtbauer", reverse_beta_unstd) {
  if (missing(reverse_beta_unstd)) reverse_beta_unstd <- rep(FALSE, length(beta_unstd))
  reverse_beta_unstd[is.na(reverse_beta_unstd)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       error_message = paste0("The number of people exposed/non-exposed, and standard deviation of the outcome ",
                                              "should be >0."),
                       func = "es_from_beta_unstd")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })


  sd_pooled <- suppressWarnings(
    sqrt(abs(((sd_dv^2 * (n_exp + n_nexp - 1)) - (beta_unstd^2 * ((n_exp * n_nexp) / (n_exp + n_nexp)))) /
      (n_exp + n_nexp - 2)))
  )

  d <- beta_unstd / sd_pooled

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_beta_unstd
  )

  es$info_used <- "beta_unstd"

  return(es)
}
