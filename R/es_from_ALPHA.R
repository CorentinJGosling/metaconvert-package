#' Convert Cronbach's alpha into an effect size measure
#'
#' @param cronbach_alpha Cronbach's alpha reliability coefficient
#' @param n_sample the total number of participants that completed the scale
#' @param n_items number of items in the scale
#' @param alpha_to_es method used to compute the effect size from Cronbach's alpha.
#'   Must be either \code{"bonett"} or \code{"raw"}.
#'
#' @details
#' This function computes an effect size from a Cronbach's alpha.
#'
#' 1. When \code{alpha_to_es = "bonett"} (default), the Bonett (2002) transformation
#' is applied:
#' \deqn{T(\alpha) = \ln(1 - \alpha)}
#' \deqn{T\_se = \sqrt{\frac{2k}{(k - 1)(n - 2)}}}
#'
#' 2. When \code{alpha_to_es = "raw"}, the raw alpha is used. Its standard error is
#' the delta-method back-transform of the Bonett (2002) transformed variance
#' (equivalently the van Zyl, Neudecker & Nel, 2000, asymptotic variance; this is
#' the form implemented by \code{metafor}'s \code{measure = "ARAW"}):
#' \deqn{\alpha\_se = (1 - \alpha) \sqrt{\frac{2k}{(k - 1)(n - 2)}}}
#' Note that Feldt et al.'s (1987) classical asymptotic variance instead uses an
#' \eqn{(n - 1)} denominator; the \eqn{(n - 2)} form above follows Bonett (2002)
#' for consistency with the Bonett transformation used in method 1.
#'
#' The Bonett transformation stabilizes the variance and is recommended for meta-analysis.
#'
#' @export es_from_cronbach_alpha
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient alpha.
#' Journal of Educational and Behavioral Statistics, 27(4), 335-340.
#'
#' Feldt, L. S., Woodruff, D. J., & Salih, F. A. (1987). Statistical inference for coefficient alpha.
#' Applied Psychological Measurement, 11(1), 93-103.
#'
#' van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution of the maximum likelihood
#' estimator of Cronbach's alpha. Psychometrika, 65(3), 271-280.
#'
#' @md
#'
#' @return
#' This function estimates the standard error of the Cronbach's alpha.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab alpha\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab N/A\cr
#'  \tab \cr
#'  \code{required input data} \tab cronbach_alpha + n_sample + n_items\cr
#'  \tab \cr
#' }
#'
#' @examples
#' es_from_cronbach_alpha(
#'   cronbach_alpha = 0.85, n_sample = 200, n_items = 10
#' )
es_from_cronbach_alpha <- function(cronbach_alpha, n_sample, n_items,
                                   alpha_to_es = "bonett") {

  if (missing(n_sample)) n_sample <- rep(NA, length(cronbach_alpha))
  if (missing(n_items)) n_items <- rep(NA, length(cronbach_alpha))

  if (!alpha_to_es %in% c("bonett", "raw")) {
    stop(paste0("'", alpha_to_es, "' not in tolerated values for the 'alpha_to_es' argument. ",
                "Possible inputs are: 'bonett', 'raw'"))
  }

  nn_miss <- which(!is.na(cronbach_alpha) & !is.na(n_sample) & !is.na(n_items))

  n <- length(cronbach_alpha)
  alpha_es <- rep(NA_real_, n)
  alpha_es_se <- rep(NA_real_, n)

  if (length(nn_miss) != 0) {
    a <- cronbach_alpha[nn_miss]
    ns <- n_sample[nn_miss]
    k <- n_items[nn_miss]

    if (alpha_to_es == "bonett") {
      alpha_es[nn_miss] <- log(1 - a)
      alpha_es_se[nn_miss] <- sqrt(2 * k / ((k - 1) * (ns - 2)))
    } else {
      alpha_es[nn_miss] <- a
      alpha_es_se[nn_miss] <- (1 - a) * sqrt(2 * k / ((k - 1) * (ns - 2)))
    }
  }

  alpha_ci_lo <- alpha_es - qnorm(0.975) * alpha_es_se
  alpha_ci_up <- alpha_es + qnorm(0.975) * alpha_es_se

  result <- data.frame(
    alpha = alpha_es,
    alpha_se = alpha_es_se,
    alpha_ci_lo = alpha_ci_lo,
    alpha_ci_up = alpha_ci_up,
    n_sample = n_sample,
    n_items = n_items,
    info_used = "cronbach_alpha"
  )

  return(result)
}
