#' Convert the number of cases and the person-time of disease-free observation in two independent groups into an incidence rate ratio (IRR)
#'
#' @param n_cases_exp number of cases in the exposed group
#' @param n_cases_nexp number of cases in the non-exposed group
#' @param time_exp person-time of disease-free observation in the exposed group
#' @param time_nexp person-time of disease-free observation in the non-exposed group
#' @param baseline_rate incidence rate of events (per person-time) in the non-exposed group (n_cases_nexp / time_nexp is used when missing)
#' @param reverse_irr a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function estimates the incidence rate ratio from the number of cases and
#' the person-time of disease-free observation in two independent groups.
#'
#' **The formula used to obtain the IRR and its standard error** are (Cochrane Handbook (section 6.7.1):
#' \deqn{logirr = log(\frac{n\_cases\_exp / time\_exp}{n\_cases\_nexp / time\_nexp)}}
#' \deqn{logirr\_se = \sqrt{\frac{1}{n\_cases\_exp} + \frac{1}{n\_cases\_nexp}}}
#'
#' **To estimate a person-time NNT** (Mayne et al., 2006), the following formulas are used:
#' \deqn{ird = baseline\_rate \times (1 - irr)}
#' \deqn{nnt = \frac{1}{ird}}
#' where \code{ird} is the incidence rate difference and \code{baseline_rate} is the
#' incidence rate in the control group.
#'
#' **To estimate the standard error of the IRD**, two formulas are used.
#' When \code{baseline_rate} is missing:
#' \deqn{ird\_se = \sqrt{\frac{n\_cases\_exp}{time\_exp^2} + \frac{n\_cases\_nexp}{time\_nexp^2}}}
#' When \code{baseline_rate} is entered by users, the delta method is used:
#' \deqn{ird\_se = baseline\_rate \times IRR \times logirr\_se}
#'
#' @return
#' This function estimates IRR and, when baseline rate information is available, NNT.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab IRR + NNT\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab N/A\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 5. Incidence Ratio Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_cases_time
#'
#' @references
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' Mayne, T. J., Whalen, E., & Rost, K. (2006). Annualized was found better than absolute risk reduction in the calculation of number needed to treat in chronic conditions. Journal of clinical epidemiology, 59(3), 217-223.
#'
#' @md
#'
#' @examples
#' es_from_cases_time(
#'   n_cases_exp = 241, n_cases_nexp = 554,
#'   time_exp = 12.764, time_nexp = 19.743
#' )
es_from_cases_time <- function(n_cases_exp, n_cases_nexp, time_exp, time_nexp,
                                baseline_rate, reverse_irr) {
  if (missing(reverse_irr)) reverse_irr <- rep(FALSE, length(n_cases_exp))
  reverse_irr[is.na(reverse_irr)] <- FALSE
  if (length(reverse_irr) == 1) reverse_irr = c(rep(reverse_irr, length(n_cases_exp)))
  if (length(reverse_irr) != length(n_cases_exp)) stop("The length of the 'reverse_irr' argument is incorrectly specified.")
  if (missing(baseline_rate)) baseline_rate <- rep(NA_real_, length(n_cases_exp))

  # IRR (log scale)
  logirr_raw <- log((n_cases_exp / time_exp) / (n_cases_nexp / time_nexp))
  logirr_se <- sqrt(1 / n_cases_exp + 1 / n_cases_nexp)

  es <- data.frame(
    logirr = ifelse(reverse_irr, -logirr_raw, logirr_raw),
    logirr_se = logirr_se
  )

  es$logirr_ci_lo <- es$logirr - es$logirr_se * qnorm(.975)
  es$logirr_ci_up <- es$logirr + es$logirr_se * qnorm(.975)

  # Person-time NNT (Mayne et al. 2006)
  provided_baseline <- !is.na(baseline_rate)
  baseline_rate <- ifelse(is.na(baseline_rate), n_cases_nexp / time_nexp, baseline_rate)

  irr <- exp(logirr_raw)

  ird <- baseline_rate * (1 - irr)

  es$rd <- ifelse(reverse_irr, -ird, ird)
  # delta method if baseline_rate given, poisson variance otherwise
  ird_se <- ifelse(
    provided_baseline,
    baseline_rate * irr * logirr_se,
    sqrt(n_cases_exp / time_exp^2 + n_cases_nexp / time_nexp^2)
  )
  es$rd_se <- ird_se
  es$rd_ci_lo <- es$rd - qnorm(.975) * ird_se
  es$rd_ci_up <- es$rd + qnorm(.975) * ird_se

  es$nnt <- ifelse(ird == 0, NA, 1 / ird)
  es$nnt <- ifelse(reverse_irr, -es$nnt, es$nnt)
  es$nnt_se <- ifelse(ird == 0, NA, ird_se / ird^2)
  ird_ci_lo <- ird - qnorm(.975) * ird_se
  ird_ci_up <- ird + qnorm(.975) * ird_se
  crosses_zero <- (ird_ci_lo < 0 & ird_ci_up > 0) | ird == 0
  es$nnt_ci_lo <- ifelse(crosses_zero, NA,
                          ifelse(reverse_irr, -1 / ird_ci_lo, 1 / ird_ci_up))
  es$nnt_ci_up <- ifelse(crosses_zero, NA,
                          ifelse(reverse_irr, -1 / ird_ci_up, 1 / ird_ci_lo))

  es$info_used <- "cases_time"
  return(es)
}
