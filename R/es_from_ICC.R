#' Convert an intraclass correlation coefficient (ICC) into an effect size measure
#'
#' @param icc intraclass correlation coefficient value
#' @param n_sample total sample size (number of subjects)
#' @param n_measurements number of measurements or raters
#' @param icc_type ICC type: \code{"agreement"} for ICC(2,1) (two-way random, absolute agreement)
#'   or \code{"consistency"} for ICC(3,1) (two-way mixed, consistency). Default is \code{"agreement"}.
#' @param icc_to_es method used to compute the effect size from ICC.
#'   Must be either \code{"bonett"} or \code{"raw"}.
#'
#' @details
#' This function computes an effect size from an ICC.
#'
#' 1. When \code{icc_to_es = "bonett"} (default), the Bonett (2002) transformation
#' is applied:
#' \deqn{T(ICC) = \ln(1 - ICC)}
#' For ICC(2,1) (\code{icc_type = "agreement"}):
#' \deqn{T\_se = \sqrt{\frac{2 (1 + (k-1) ICC)^2}{k (k - 1)(n - 1)}}}
#' For ICC(3,1) (\code{icc_type = "consistency"}):
#' \deqn{T\_se = \sqrt{\frac{2}{(k - 1)(n - 1)}}}
#'
#' 2. When \code{icc_to_es = "raw"}, the raw ICC is used and its standard error
#' is obtained by the delta method.
#'
#' The ICC(2,1) SE is approximated using the one-way random model formula
#' (Bonett, 2002; Donner & Eliasziw, 1987).
#'
#' @export es_from_icc
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for estimating intraclass correlations
#' with desired precision. Statistics in Medicine, 21(9), 1331-1335.
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: uses in assessing rater
#' reliability. Psychological Bulletin, 86(2), 420-428.
#'
#' @md
#'
#' @return
#' This function estimates the standard error of the ICC.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab icc\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab N/A\cr
#'  \tab \cr
#'  \code{required input data} \tab icc + n_sample + n_measurements\cr
#'  \tab \cr
#' }
#'
#' @examples
#' es_from_icc(
#'   icc = 0.80, n_sample = 50, n_measurements = 2, icc_type = "agreement"
#' )
es_from_icc <- function(icc, n_sample, n_measurements, icc_type = "agreement",
                        icc_to_es = "bonett") {

  if (missing(n_sample)) n_sample <- rep(NA, length(icc))
  if (missing(n_measurements)) n_measurements <- rep(NA, length(icc))
  if (missing(icc_type)) icc_type <- rep("agreement", length(icc))
  icc_type[is.na(icc_type)] <- "agreement"
  if (length(icc_type) == 1) icc_type <- rep(icc_type, length(icc))

  if (!all(icc_type %in% c("agreement", "consistency"))) {
    stop(paste0("'", unique(icc_type[!icc_type %in% c("agreement", "consistency")]),
                "' not in tolerated values for the 'icc_type' argument. ",
                "Possible inputs are: 'agreement', 'consistency'"))
  }

  if (!icc_to_es %in% c("bonett", "raw")) {
    stop(paste0("'", icc_to_es, "' not in tolerated values for the 'icc_to_es' argument. ",
                "Possible inputs are: 'bonett', 'raw'"))
  }

  nn_miss <- which(!is.na(icc) & !is.na(n_sample) & !is.na(n_measurements))

  n <- length(icc)
  icc_es <- rep(NA_real_, n)
  icc_es_se <- rep(NA_real_, n)

  if (length(nn_miss) != 0) {
    rho <- icc[nn_miss]
    ns <- n_sample[nn_miss]
    k <- n_measurements[nn_miss]
    type <- icc_type[nn_miss]

    bonett_transformed_se <- ifelse(
      type == "agreement",
      sqrt(2 * (1 + (k - 1) * rho)^2 / (k * (k - 1) * (ns - 1))),
      sqrt(2 / ((k - 1) * (ns - 1)))
    )

    if (icc_to_es == "bonett") {
      icc_es[nn_miss] <- log(1 - rho)
      icc_es_se[nn_miss] <- bonett_transformed_se
    } else {
      icc_es[nn_miss] <- rho
      icc_es_se[nn_miss] <- (1 - rho) * bonett_transformed_se
    }
  }

  icc_ci_lo <- icc_es - qnorm(0.975) * icc_es_se
  icc_ci_up <- icc_es + qnorm(0.975) * icc_es_se

  result <- data.frame(
    icc = icc_es,
    icc_se = icc_es_se,
    icc_ci_lo = icc_ci_lo,
    icc_ci_up = icc_ci_up,
    n_sample = n_sample,
    n_measurements = n_measurements,
    icc_type = icc_type,
    info_used = "icc"
  )

  return(result)
}
