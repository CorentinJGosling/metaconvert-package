#' Convert a Spearman's rank correlation coefficient to several effect size measures
#'
#' @param spearman_r a Spearman's rank correlation coefficient value
#' @param n_sample the total number of participants
#' @param sd_iv the standard deviation of the independent variable
#' @param unit_increase_iv a value of the independent variable that will be used to estimate the Cohen's d (see details).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Must be either "sd" or "value"
#' @param n_exp number of the experimental/exposed group
#' @param n_nexp number of the non-experimental/non-exposed group
#' @param cor_to_smd formula used to convert the derived Pearson's r value into a SMD.
#' @param reverse_spearman_r a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first converts a Spearman's rank correlation coefficient to an approximate Pearson's
#' correlation coefficient using the formula proposed by Rupinski & Dunlap (1996):
#' \deqn{r_p = 2 \sin(\pi / 6 \times r_s)}
#'
#' The standard error of the resulting Pearson's r is derived via the delta method:
#' \deqn{r_p\_se = \sqrt{\left(\frac{\pi}{3} \cos\left(\frac{\pi}{6} r_s\right)\right)^2 \times \frac{(1 - r_s^2)^2}{n - 1}}}
#'
#' The converted Pearson's r is then further converted to a Fisher's z, Cohen's d, Hedges' g, and
#' odds ratio using the same formulas as \code{\link{es_from_pearson_r}()}.
#'
#' @export es_from_spearman_rho
#'
#' @references
#' Rupinski, M. T., & Dunlap, W. P. (1996). Approximating Pearson product-moment correlations from
#' Kendall's tau and Spearman's rho. Educational and Psychological Measurement, 56(3), 419-429.
#'
#' @md
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab R + Z\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + OR\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 4. Pearson's r or Fisher's z'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @examples
#' es_from_spearman_rho(
#'   spearman_r = .55, n_sample = 100
#' )
es_from_spearman_rho <- function(spearman_r, n_sample,
                                 n_exp, n_nexp, cor_to_smd = "viechtbauer",
                                 sd_iv, unit_increase_iv, unit_type = "raw_scale",
                                 reverse_spearman_r) {
  if (missing(reverse_spearman_r)) {
    reverse_spearman_r <- rep(FALSE, length(spearman_r))
  }
  reverse_spearman_r[is.na(reverse_spearman_r)] <- FALSE
  if (length(reverse_spearman_r) == 1) reverse_spearman_r = c(rep(reverse_spearman_r, length(spearman_r)))
  if (length(reverse_spearman_r) != length(spearman_r)) stop("The length of the 'reverse_spearman_r' argument is incorrectly specified.")

  if (missing(n_exp)) {
    n_exp <- rep(NA, length(spearman_r))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA, length(spearman_r))
  }
  if (missing(n_sample)) {
    n_sample <- rep(NA, length(spearman_r))
  }
  if (missing(sd_iv)) {
    sd_iv <- rep(NA, length(spearman_r))
  }
  if (missing(unit_increase_iv)) {
    unit_increase_iv <- rep(NA, length(spearman_r))
  }
  if (missing(unit_type)) {
    unit_type <- rep(NA, length(spearman_r))
  }

  # Rupinski & Dunlap (1996)
  r <- 2 * sin(pi / 6 * spearman_r)

  es <- es_from_pearson_r(
    pearson_r = r, sd_iv = sd_iv, n_sample = n_sample,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    reverse_pearson_r = reverse_spearman_r
  )

  # delta method: d/dr_s 2*sin(pi/6*r_s) = (pi/3)*cos(pi/6*r_s); var(r_s) = (1-r_s^2)^2/(n-1)
  r_applied <- ifelse(reverse_spearman_r, -spearman_r, spearman_r)
  deriv <- (pi / 3) * cos(pi / 6 * r_applied)
  var_spearman <- (1 - r_applied^2)^2 / (n_sample - 1)
  r_se_delta <- sqrt(deriv^2 * var_spearman)

  es$r_se <- r_se_delta
  es$r_ci_lo <- es$r - qt(.975, n_sample - 2) * r_se_delta
  es$r_ci_up <- es$r + qt(.975, n_sample - 2) * r_se_delta

  # dz/dr = 1/(1-r^2)
  r_bounded <- pmin(pmax(es$r, -0.9999), 0.9999)
  es$z_se <- r_se_delta / (1 - r_bounded^2)
  es$z_ci_lo <- es$z - qnorm(.975) * es$z_se
  es$z_ci_up <- es$z + qnorm(.975) * es$z_se

  es$info_used <- "spearman_r"
  return(es)
}
