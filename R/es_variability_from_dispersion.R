#' Convert means and/or standard deviations of two independent groups into two effect measures (VR/CVR)
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_sd_exp standard deviation of participants in the experimental/exposed group.
#' @param mean_sd_nexp standard deviation of participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param reverse_means_variability a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the means and standard deviations of two independent groups into a log variability ratio (VR) and a log coefficient of variation ratio (CVR).
#'
#' The formulas used to obtain the log VR are (formulas 5 and 15, Senior et al. 2020):
#' \deqn{logvr = log(\frac{mean\_sd\_exp}{mean\_sd\_nexp}) + \frac{1}{2 * (n\_exp - 1)} - \frac{1}{2 * (n\_nexp - 1)}}
#' \deqn{logvr\_se = \sqrt{\frac{1}{2} * (\frac{n\_nexp}{(n\_nexp - 1)^2} + \frac{n\_exp}{(n\_exp - 1)^2})}}
#' \deqn{logvr\_ci\_lo = logvr - qnorm(.975) * logvr\_se}
#' \deqn{logvr\_ci\_up = logvr + qnorm(.975) * logvr\_se}
#'
#' The formulas used to obtain the log CVR are (formulas 6 and 16, Senior et al. 2020):
#' \deqn{cvt = mean\_sd\_exp / mean\_exp}
#' \deqn{cvc = mean\_sd\_nexp / mean\_nexp}
#' \deqn{logcvr = log(\frac{cvt}{cvc}) + \frac{1}{2} * (\frac{1}{n\_exp - 1} - \frac{1}{n\_nexp - 1}) + \frac{1}{2} * (\frac{mean\_sd\_nexp^2}{n\_nexp * mean\_nexp^2} - \frac{mean\_sd\_exp^2}{n\_exp * mean\_exp^2})}
#' \deqn{vt\_exp = \frac{mean\_sd\_exp^2}{n\_exp * mean\_exp^2} + \frac{mean\_sd\_exp^4}{2 * n\_exp^2 * mean\_exp^4} + \frac{n\_exp}{2 * (n\_exp - 1)^2}}
#' \deqn{vt\_nexp = \frac{mean\_sd\_nexp^2}{n\_nexp * mean\_nexp^2} + \frac{mean\_sd\_nexp^4}{2 * n\_nexp^2 * mean\_nexp^4} + \frac{n\_nexp}{2 * (n\_nexp - 1)^2}}
#' \deqn{logcvr\_se = \sqrt{vt\_exp + vt\_nexp}}
#' \deqn{logcvr\_ci\_lo = logcvr - qnorm(.975) * logcvr\_se}
#' \deqn{logcvr\_ci\_up = logcvr + qnorm(.975) * logcvr\_se}
#'
#' @references
#' Senior, A. M., Viechtbauer, W., & Nakagawa, S. (2020). Revisiting and expanding the meta-analysis of variation: The log coefficient of variation ratio. Research Synthesis Methods, 11(4), 553-567. https://doi.org/10.1002/jrsm.1423
#'
#' @return
#' This function estimates VR and CVR
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab VR + CVR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab No conversion performed\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 6. Variability indices'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_variab_from_means_sd
#'
#' @md
#'
#' @examples
#' es_variab_from_means_sd(
#'   n_exp = 55, n_nexp = 55,
#'   mean_exp = 2.3, mean_sd_exp = 1.2,
#'   mean_nexp = 1.9, mean_sd_nexp = 0.9
#' )
es_variab_from_means_sd <- function(mean_exp, mean_nexp, mean_sd_exp, mean_sd_nexp,
                                    n_exp, n_nexp, reverse_means_variability) {
  if (missing(reverse_means_variability)) reverse_means_variability <- rep(FALSE, length(mean_sd_exp))
  reverse_means_variability[is.na(reverse_means_variability)] <- FALSE
  if (length(reverse_means_variability) == 1) reverse_means_variability = c(rep(reverse_means_variability, length(mean_sd_exp)))
  if (length(reverse_means_variability) != length(mean_sd_exp)) stop("The length of the 'reverse_means_variability' argument is incorrectly specified.")

  es <- data.frame(logcvr = rep(NA, length(mean_sd_exp)))
  # es$logcvr = log(mean_sd_exp/mean_exp) + 1/(2 * (n_exp - 1)) - log(mean_sd_nexp/mean_nexp) - 1/(2 * (n_nexp - 1))
  # yi <- log(sd1i/m1i) + 1/(2 * (n1i - 1)) - log(sd2i/m2i) -
  #   1/(2 * (n2i - 1))

  # lnCVRIND2 (equation 6) Nakagawa 202
  CVT = mean_sd_exp / mean_exp
  CVC = mean_sd_nexp / mean_nexp
  logcvr <- suppressWarnings(
    log(CVT / CVC) +
      1/2 * (1/(n_exp-1) - 1/(n_nexp-1)) +
      1/2 * (mean_sd_nexp^2 / (n_nexp * mean_nexp^2) - mean_sd_exp^2 / (n_exp * mean_exp^2))
  )
  es$logcvr <- ifelse(reverse_means_variability, -logcvr, logcvr)

  # s2IND2 (equation 16) Senior et al. 2020.
  # Note the 1/2 factor on the Var(ln s) term: for a normal sample
  # Var(ln s) is about trigamma((n-1)/2)/4, which is about 0.5 * n/(n-1)^2, and the
  # same quantity in logvr_se below is written 1/2 * n/(n-1)^2. Without the 1/2 the
  # lnCVR variance is about twice too large, making the SE about 40% too wide and the
  # inverse-variance weights about 2x off.
  es$logcvr_se <- sqrt(
    mean_sd_exp^2 / (n_exp * mean_exp^2) +
    mean_sd_exp^4 / (2 * n_exp^2 * mean_exp^4) +
    0.5 * n_exp / ((n_exp - 1)^2) +
    mean_sd_nexp^2 / (n_nexp * mean_nexp^2) +
    mean_sd_nexp^4 / (2 * n_nexp^2 * mean_nexp^4) +
    0.5 * n_nexp / ((n_nexp - 1)^2))

  es$logcvr_ci_lo <- es$logcvr - qnorm(.975) * es$logcvr_se
  es$logcvr_ci_up <- es$logcvr + qnorm(.975) * es$logcvr_se

  logvr <- suppressWarnings(
    log(mean_sd_exp / mean_sd_nexp) + 1 / (2 * (n_exp - 1)) - 1 / (2 * (n_nexp - 1))
  )
  es$logvr <- ifelse(reverse_means_variability, -logvr, logvr)
  es$logvr_se <- sqrt(1/2 * (n_nexp/((n_nexp - 1)^2) + n_exp/((n_exp - 1)^2)))
  # es$logvr_se <- sqrt(1 / (2 * (n_exp - 1)) + 1 / (2 * (n_nexp - 1)))
  es$logvr_ci_lo <- es$logvr - qnorm(.975) * es$logvr_se
  es$logvr_ci_up <- es$logvr + qnorm(.975) * es$logvr_se

  es$info_used <- "variability_means_sd"

  return(es)
}

#' Convert means and/or standard errors of two independent groups into two effect measures (VR/CVR)
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_se_exp standard error of participants in the experimental/exposed group.
#' @param mean_se_nexp standard error of participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param reverse_means_variability a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the standard errors into standard deviations
#' (\eqn{SE = SD / \sqrt{n}}, so \eqn{SD = SE \times \sqrt{n}}).
#' \deqn{mean\_sd\_exp = mean\_se\_exp * \sqrt{n\_exp}}
#' \deqn{mean\_sd\_nexp = mean\_se\_nexp * \sqrt{n\_nexp}}
#'
#' Then, calculations of the \code{\link{es_variab_from_means_sd}} are applied.
#'
#' @return
#' This function estimates VR and CVR
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab VR + CVR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab No conversion performed\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 6. Variability indices'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Senior, A. M., Viechtbauer, W., & Nakagawa, S. (2020). Revisiting and expanding the meta-analysis of variation: The log coefficient of variation ratio. Research Synthesis Methods, 11(4), 553-567. https://doi.org/10.1002/jrsm.1423
#'
#' @export es_variab_from_means_se
#'
#' @md
#'
#' @examples
#' es_variab_from_means_se(
#'   mean_exp = 42, mean_se_exp = 11,
#'   mean_nexp = 42, mean_se_nexp = 15,
#'   n_exp = 43, n_nexp = 34
#' )
es_variab_from_means_se <- function(mean_exp, mean_nexp, mean_se_exp, mean_se_nexp,
                                    n_exp, n_nexp, reverse_means_variability) {
  if (missing(reverse_means_variability)) reverse_means_variability <- rep(FALSE, length(mean_se_exp))
  reverse_means_variability[is.na(reverse_means_variability)] <- FALSE

  sd_exp <- mean_se_exp * sqrt(n_exp)
  sd_nexp <- mean_se_nexp * sqrt(n_nexp)

  es <- es_variab_from_means_sd(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_sd_exp = sd_exp, mean_sd_nexp = sd_nexp,
    n_exp = n_exp, n_nexp = n_nexp, reverse_means_variability = reverse_means_variability
  )
  es$info_used <- "variability_means_se"

  return(es)
}

#' Convert means and/or 95% CIs of two independent groups into two effect measures (VR/CVR)
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_ci_lo_exp lower bound of the 95% CI of the mean of the experimental/exposed group
#' @param mean_ci_up_exp upper bound of the 95% CI of the mean of the experimental/exposed group
#' @param mean_ci_lo_nexp lower bound of the 95% CI of the mean of the non-experimental/non-exposed group.
#' @param mean_ci_up_nexp upper bound of the 95% CI of the mean of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param reverse_means_variability a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the bounds of the 95% CI of the means of two independent groups into standard errors.
#' \deqn{mean\_se\_exp = \frac{mean\_ci\_up\_exp - mean\_ci\_lo\_exp}{2 * qt{(0.975, df = n\_exp - 1)}}}
#' \deqn{mean\_se\_nexp = \frac{mean\_ci\_up\_nexp - mean\_ci\_lo\_nexp}{2 * qt{(0.975, df = n\_nexp - 1)}}}
#'
#' Then, calculations of the \code{\link{es_variab_from_means_se}} are applied.
#'
#' @references
#' Senior, A. M., Viechtbauer, W., & Nakagawa, S. (2020). Revisiting and expanding the meta-analysis of variation: The log coefficient of variation ratio. Research Synthesis Methods, 11(4), 553-567. https://doi.org/10.1002/jrsm.1423
#'
#' @export es_variab_from_means_ci
#'
#' @return
#' This function estimates VR and CVR
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab VR + CVR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab No conversion performed\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 6. Variability indices'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_variab_from_means_ci(
#'   mean_exp = 42, mean_ci_lo_exp = 32, mean_ci_up_exp = 52,
#'   mean_nexp = 42, mean_ci_lo_nexp = 37, mean_ci_up_nexp = 47,
#'   n_exp = 43, n_nexp = 34
#' )
es_variab_from_means_ci <- function(mean_exp, mean_nexp,
                                    mean_ci_lo_exp, mean_ci_up_exp,
                                    mean_ci_lo_nexp, mean_ci_up_nexp,
                                    n_exp, n_nexp, reverse_means_variability) {
  if (missing(reverse_means_variability)) reverse_means_variability <- rep(FALSE, length(mean_ci_up_exp))
  reverse_means_variability[is.na(reverse_means_variability)] <- FALSE

  se_exp <- (mean_ci_up_exp - mean_ci_lo_exp) / (2 * qt(0.975, n_exp - 1))
  se_nexp <- (mean_ci_up_nexp - mean_ci_lo_nexp) / (2 * qt(0.975, n_nexp - 1))

  es <- es_variab_from_means_se(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_se_exp = se_exp, mean_se_nexp = se_nexp,
    n_exp = n_exp, n_nexp = n_nexp, reverse_means_variability = reverse_means_variability
  )
  es$info_used <- "variability_means_ci"

  return(es)
}
