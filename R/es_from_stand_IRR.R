#' Convert the number of cases and the person-time of disease-free observation in two independent groups into an incidence rate ratio (IRR)
#'
#' @param n_cases_exp number of cases in the exposed group
#' @param n_cases_nexp number of cases in the non-exposed group
#' @param time_exp person-time of disease-free observation in the exposed group
#' @param time_nexp person-time of disease-free observation in the non-exposed group
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
#' @return
#' This function estimates IRR.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab IRR\cr
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
#' @md
#'
#' @examples
#' es_from_cases_time(
#'   n_cases_exp = 241, n_cases_nexp = 554,
#'   time_exp = 12.764, time_nexp = 19.743
#' )
es_from_cases_time <- function(n_cases_exp, n_cases_nexp, time_exp, time_nexp, reverse_irr) {
  if (missing(reverse_irr)) reverse_irr <- rep(FALSE, length(n_cases_exp))
  reverse_irr[is.na(reverse_irr)] <- FALSE
  if (length(reverse_irr) == 1) reverse_irr = c(rep(reverse_irr, length(n_cases_exp)))
  if (length(reverse_irr) != length(n_cases_exp)) stop("The length of the 'reverse_irr' argument is incorrectly specified.")

  tryCatch({
    .validate_positive(n_cases_exp, n_cases_nexp, time_exp, time_nexp,
                       error_message = paste0("The number of people exposed/non-exposed and time of disease-free observation time ",
                                              "should be >0."),
                       func = "es_from_cases_time")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  es <- data.frame(
    logirr = log((n_cases_exp / time_exp) / (n_cases_nexp / time_nexp)),
    logirr_se = sqrt(1 / n_cases_exp + 1 / n_cases_nexp)
  )

  es$logirr <- ifelse(reverse_irr, -es$logirr, es$logirr)

  es$logirr_ci_lo <- es$logirr - es$logirr_se * qnorm(.975)
  es$logirr_ci_up <- es$logirr + es$logirr_se * qnorm(.975)

  es$info_used <- "cases_time"
  return(es)
}
