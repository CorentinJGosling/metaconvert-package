#' Disattenuate (correct for unreliability) a correlation coefficient
#'
#' @param r observed correlation coefficient (values outside \eqn{[-1, 1]}
#'   trigger a warning; the arithmetic is still applied so the inputs can be
#'   inspected)
#' @param r_se standard error of the observed correlation. Optional: when omitted,
#'   it is derived from \code{n_sample} (see Details).
#' @param reliability_x reliability of measure X (e.g., target PROM)
#' @param reliability_y reliability of measure Y (e.g., comparator instrument)
#' @param n_sample sample size. Used only to derive \code{r_se} when \code{r_se}
#'   is not supplied; otherwise unused.
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
#' When \code{r_se} is not supplied, it is first obtained from the sample size as
#' the large-sample Pearson correlation SE \eqn{\sqrt{(1 - r^2)^2 / (n - 1)}}
#' (Cooper et al., 2019).
#'
#' The function also provides Fisher's z transformation of the corrected
#' correlation for use in meta-analysis, with
#' \eqn{SE(z_c) = SE(r_c) / (1 - r_c^2)} (delta method). The corrected-r
#' confidence interval is obtained by back-transforming the Fisher-z interval
#' (\eqn{\tanh}), so it always lies within (-1, 1).
#'
#' When the corrected correlation is extreme (\eqn{|r_c| > 0.999}, including
#' the mathematically impossible \eqn{|r_c| \ge 1} that inconsistent inputs
#' produce), no meaningful Fisher's z or CI exists: the corrected-r CI and all
#' Fisher's z outputs (\code{z_corrected}, \code{z_corrected_se} and its CI
#' bounds) are set to \code{NA} and a warning is emitted. The corrected point
#' estimate and its first-order standard error are always returned as
#' computed, so the offending inputs can be inspected.
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
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of
#' research synthesis and meta-analysis. Russell Sage Foundation.
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

  if (missing(n_sample)) n_sample <- rep(NA_real_, length(r))
  if (missing(r_se)) r_se <- rep(NA_real_, length(r))
  # Recycle scalar inputs to the length of r: the masked indexing below
  # (n_sample[need_se]) NA-pads when the vector is shorter than the mask,
  # which used to silently yield NA SEs for every row after the first.
  if (length(r_se) == 1) r_se <- rep(r_se, length(r))
  if (length(n_sample) == 1) n_sample <- rep(n_sample, length(r))
  if (length(reliability_x) == 1) reliability_x <- rep(reliability_x, length(r))
  if (length(reliability_y) == 1) reliability_y <- rep(reliability_y, length(r))

  # When r_se is not supplied, derive it from the sample size using the
  # large-sample Pearson correlation variance (Cooper et al., 2019):
  # SE(r) = sqrt((1 - r^2)^2 / (n - 1)).
  need_se <- is.na(r_se) & !is.na(n_sample) & !is.na(r)
  r_se[need_se] <- sqrt((1 - r[need_se]^2)^2 / (n_sample[need_se] - 1))

  if (any(!is.na(r) & abs(r) > 1, na.rm = TRUE)) {
    warning("Observed r outside [-1, 1] detected; the disattenuated r is invalid")
  }
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

  # Single constant for the extreme trigger AND the atanh bound: a corrected r
  # beyond it has no meaningful Fisher's z, so every derived-scale output is
  # suppressed below (the bounded value never surfaces in the output).
  extreme_threshold <- 0.999
  extreme <- which(!is.na(r_corrected) & abs(r_corrected) > extreme_threshold)
  if (length(extreme) > 0) {
    warning(sprintf(
      "Corrected r exceeds %s in %d row(s); the CI and Fisher's z outputs are set to NA as unreliable. Check the reliability inputs (and the observed r) - this often reflects a reliability value that is too low for the observed correlation.",
      extreme_threshold, length(extreme)
    ), call. = FALSE)
  }
  # bound corrected r to avoid Inf in atanh (extreme rows are NA'd below)
  r_bounded <- pmin(pmax(r_corrected, -extreme_threshold), extreme_threshold)
  z_corrected <- atanh(r_bounded)
  # dz/dr = 1/(1-r^2)
  z_corrected_se <- r_corrected_se / (1 - r_bounded^2)
  z_corrected_ci_lo <- z_corrected - qnorm(0.975) * z_corrected_se
  z_corrected_ci_up <- z_corrected + qnorm(0.975) * z_corrected_se

  # CI from back-transformed z interval
  r_corrected_ci_lo <- tanh(z_corrected_ci_lo)
  r_corrected_ci_up <- tanh(z_corrected_ci_up)
  # For extreme rows, NA every derived-scale output (r CI AND the z columns):
  # a bounded atanh would be a pure artifact (identical whatever the input),
  # and r_c >= 1 has no valid Fisher's z or SE. Only the point estimate and
  # its first-order SE are kept.
  if (length(extreme) > 0) {
    r_corrected_ci_lo[extreme] <- NA_real_
    r_corrected_ci_up[extreme] <- NA_real_
    z_corrected[extreme] <- NA_real_
    z_corrected_se[extreme] <- NA_real_
    z_corrected_ci_lo[extreme] <- NA_real_
    z_corrected_ci_up[extreme] <- NA_real_
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
