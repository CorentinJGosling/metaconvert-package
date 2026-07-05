#' Disattenuate (correct for unreliability) a correlation coefficient
#'
#' @param r observed correlation coefficient
#' @param r_se standard error of the observed correlation
#' @param reliability_x reliability of measure X (e.g., target PROM)
#' @param reliability_y reliability of measure Y (e.g., comparator instrument)
#' @param n_sample sample size (not used in the calculations)
#'
#' @details
#' Corrects an observed correlation for attenuation due to measurement error
#' in both measures, using the classical disattenuation formula
#' (Hunter & Schmidt, 2004; Spearman, 1904):
#'
#' \deqn{r_c = \frac{r_{obs}}{\sqrt{rel_x \times rel_y}}}
#'
#' where \eqn{rel_x} and \eqn{rel_y} are reliability coefficients for the
#' two measures (e.g., Cronbach's alpha, test-retest ICC).
#'
#' The standard error of the corrected correlation is approximated as:
#' \deqn{r_c\_se = \frac{r\_se}{\sqrt{rel_x \times rel_y}}}
#'
#' This approximation treats the reliabilities as known constants, so the SE
#' and CI are a lower bound when the reliabilities are themselves estimated
#' (Hunter & Schmidt, 2004, Ch. 3).
#'
#' The function also provides Fisher's z transformation of the corrected
#' correlation for use in meta-analysis, with
#' \eqn{SE(z_c) = SE(r_c) / (1 - r_c^2)} (delta method). The corrected-r
#' confidence interval is obtained by back-transforming the Fisher-z interval
#' (\eqn{\tanh}), so it always lies within (-1, 1).
#'
#' This function is typically applied to the results of
#' \code{summary(convert_df(..., measure = "r"))}:
#'
#' \preformatted{
#' res <- summary(convert_df(my_data, measure = "r"))
#' corrected <- es_disattenuate(
#'   r = res$es, r_se = res$se,
#'   reliability_x = my_data$rel_target,
#'   reliability_y = my_data$rel_comparator,
#'   n_sample = my_data$n_sample
#' )
#' }
#'
#' @export es_disattenuate
#'
#' @references
#' Hunter, J. E., & Schmidt, F. L. (2004). Methods of Meta-Analysis: Correcting
#' Error and Bias in Research Findings (2nd ed.). Sage Publications.
#'
#' Spearman, C. (1904). The proof and measurement of association between two things.
#' The American Journal of Psychology, 15(1), 72-101.
#'
#' @md
#'
#' @return
#' A data.frame with the corrected correlation and its Fisher's z
#' transformation (with their standard errors and 95% CIs), and the
#' attenuation factor.
#'
#' @examples
#' es_disattenuate(r = 0.50, r_se = 0.05,
#'                 reliability_x = 0.85, reliability_y = 0.80,
#'                 n_sample = 100)
#'
#' # Only correct for one measure's unreliability (set other to 1)
#' es_disattenuate(r = 0.50, r_se = 0.05,
#'                 reliability_x = 0.85, reliability_y = 1.0,
#'                 n_sample = 100)
es_disattenuate <- function(r, r_se, reliability_x, reliability_y, n_sample) {

  if (missing(n_sample)) n_sample <- rep(NA, length(r))

  if (any(!is.na(reliability_x) & (reliability_x <= 0 | reliability_x > 1), na.rm = TRUE)) {
    warning("reliability_x outside (0, 1] detected; corrected r may be invalid or Inf")
  }
  if (any(!is.na(reliability_y) & (reliability_y <= 0 | reliability_y > 1), na.rm = TRUE)) {
    warning("reliability_y outside (0, 1] detected; corrected r may be invalid or Inf")
  }

  A <- sqrt(reliability_x * reliability_y)

  r_corrected <- r / A

  # first-order approx
  r_corrected_se <- r_se / A

  extreme <- which(!is.na(r_corrected) & abs(r_corrected) > 0.999)
  if (length(extreme) > 0) {
    warning(sprintf(
      "Corrected r exceeds 0.999 in %d row(s); values clamped to +/-0.9999 before Fisher's z. Check reliability inputs - this often reflects a reliability value that is too low for the observed correlation.",
      length(extreme)
    ), call. = FALSE)
  }
  # bound corrected r to avoid Inf in atanh
  r_bounded <- pmin(pmax(r_corrected, -0.9999), 0.9999)
  z_corrected <- atanh(r_bounded)
  # dz/dr = 1/(1-r^2)
  z_corrected_se <- r_corrected_se / (1 - r_bounded^2)
  z_corrected_ci_lo <- z_corrected - qnorm(0.975) * z_corrected_se
  z_corrected_ci_up <- z_corrected + qnorm(0.975) * z_corrected_se

  # CI from back-transformed z interval
  r_corrected_ci_lo <- tanh(z_corrected_ci_lo)
  r_corrected_ci_up <- tanh(z_corrected_ci_up)
  # NA when r was clamped
  if (length(extreme) > 0) {
    r_corrected_ci_lo[extreme] <- NA_real_
    r_corrected_ci_up[extreme] <- NA_real_
  }

  result <- data.frame(
    r_corrected = r_corrected,
    r_corrected_se = r_corrected_se,
    r_corrected_ci_lo = r_corrected_ci_lo,
    r_corrected_ci_up = r_corrected_ci_up,
    z_corrected = z_corrected,
    z_corrected_se = z_corrected_se,
    z_corrected_ci_lo = z_corrected_ci_lo,
    z_corrected_ci_up = z_corrected_ci_up,
    attenuation_factor = A
  )

  return(result)
}
