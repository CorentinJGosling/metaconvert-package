#' Convert a Spearman's rank correlation coefficient to several effect size measures
#'
#' @param spearman_r a Spearman's rank correlation coefficient value
#' @param n_sample the total number of participants
#' @param sd_iv the standard deviation of the independent variable
#' @param unit_increase_iv a value of the independent variable that will be used to estimate the Cohen's d (see details).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Use '"sd"' when the increase is expressed in standard deviations of the independent variable, and '"raw_scale"' when it is in the raw units of that variable. '"value"' and '"raw_data"' are accepted synonyms of '"raw_scale"'. Defaults to '"raw_scale"'. Read only by \code{cor_to_smd = "mathur"}; any other value is rejected rather than silently treated as raw units.
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
#' The standard error of the resulting Pearson's r is derived via the delta method.
#' The sampling variance of the Spearman rho under bivariate normality carries an
#' extra \eqn{(1 + r_s^2 / 2)} factor relative to the Pearson variance
#' (Bonett & Wright, 2000); omitting it would understate the SE, increasingly so
#' for strong correlations:
#' \deqn{r_p\_se = \sqrt{\left(\frac{\pi}{3} \cos\left(\frac{\pi}{6} r_s\right)\right)^2 \times \left(1 + \frac{r_s^2}{2}\right)\frac{(1 - r_s^2)^2}{n - 1}}}
#'
#' The converted Pearson's r is then further converted to a Fisher's z, Cohen's d, Hedges' g, and
#' odds ratio using the same formulas as \code{\link{es_from_pearson_r}()}. This
#' delta-method SE is propagated to the Fisher's z and (for the
#' \code{cor_to_smd = "viechtbauer"} and \code{cor_to_smd = "cooper"} paths,
#' whose SMD standard errors are proportional to the input r SE) to the Cohen's d,
#' Hedges' g and odds ratio, so all measures share the corrected r variance.
#' The \code{cor_to_smd = "mathur"} SMD variance contains no r SE term, so its
#' d/g/OR standard errors are left as computed by \code{\link{es_from_pearson_r}()}.
#'
#' @export es_from_spearman_rho
#'
#' @references
#' Rupinski, M. T., & Dunlap, W. P. (1996). Approximating Pearson product-moment correlations from
#' Kendall's tau and Spearman's rho. Educational and Psychological Measurement, 56(3), 419-429.
#'
#' Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for estimating Pearson, Kendall
#' and Spearman correlations. Psychometrika, 65(1), 23-28.
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
  .validate_unit_type(unit_type)

  # Mirror the sample-size fallback es_from_pearson_r applies internally, so a
  # direct call supplying only n_exp/n_nexp gets the same delta-method r/z SEs
  # (and corrected d/g/OR SEs) as the equivalent n_sample call, instead of NA
  # r_se/z_se next to uncorrected SMD SEs.
  n_sample <- ifelse(is.na(n_sample), n_exp + n_nexp, n_sample)

  # Rupinski & Dunlap (1996)
  r <- 2 * sin(pi / 6 * spearman_r)

  es <- es_from_pearson_r(
    pearson_r = r, sd_iv = sd_iv, n_sample = n_sample,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    reverse_pearson_r = reverse_spearman_r
  )

  # SE of the converted Pearson r_p via the delta method on the Spearman -> Pearson
  # transform r_p = 2*sin(pi/6*r_s):  d r_p / d r_s = (pi/3)*cos(pi/6*r_s).
  # The sampling variance of the Spearman rho under bivariate normality carries an
  # extra (1 + r_s^2/2) factor relative to the Pearson variance (Bonett & Wright,
  # 2000); using the plain Pearson form would understate the SE (increasingly so
  # for strong correlations).
  r_applied <- ifelse(reverse_spearman_r, -spearman_r, spearman_r)
  deriv <- (pi / 3) * cos(pi / 6 * r_applied)
  var_spearman <- (1 + r_applied^2 / 2) * (1 - r_applied^2)^2 / (n_sample - 1)
  r_se_delta <- sqrt(deriv^2 * var_spearman)

  # Reuse the (correct) converted Pearson r and z that es_from_pearson_r computed.
  r_p_applied <- es$r
  z_applied <- es$z

  # Propagate the corrected r_p SE into the converted SMD/OR measures. For the
  # (default) viechtbauer path, d_se is proportional to the input r_se
  # (d_se = |d(rtod)/dr_p| * r_se), and the cooper path is likewise exactly
  # linear in r_se (d_se = sqrt(4 * r_se^2 / (1 - r^2)^3)), so scaling d_se by
  # r_se_delta / r_se_pearson yields delta-consistent d/g/OR SEs on both paths;
  # g, OR and all CIs are then rebuilt by .es_from_d exactly as elsewhere in the
  # package. Only the mathur SMD variance genuinely does not route through r_se
  # (d_se = |d| * sqrt(1/(r^2 (n-3)) + 1/(2 (n-1))) has no r_se term), so those
  # rows are left unchanged (scale = 1).
  r_se_pearson <- es$r_se
  scale <- ifelse(cor_to_smd %in% c("viechtbauer", "cooper") &
                    !is.na(r_se_pearson) &
                    r_se_pearson > 0 & is.finite(r_se_delta),
                  r_se_delta / r_se_pearson, 1)
  es <- .es_from_d(d = es$d, d_se = es$d_se * scale,
                   n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample)

  # Overwrite the round-trip r/z from .es_from_d with the direct Spearman->Pearson
  # values and their delta-method SEs.
  es$r <- r_p_applied
  es$r_se <- r_se_delta
  es$r_ci_lo <- r_p_applied - qt(.975, n_sample - 2) * r_se_delta
  es$r_ci_up <- r_p_applied + qt(.975, n_sample - 2) * r_se_delta

  # dz/dr = 1/(1-r^2)
  es$z <- z_applied
  r_bounded <- pmin(pmax(r_p_applied, -0.9999), 0.9999)
  es$z_se <- r_se_delta / (1 - r_bounded^2)
  es$z_ci_lo <- es$z - qnorm(.975) * es$z_se
  es$z_ci_up <- es$z + qnorm(.975) * es$z_se

  es$info_used <- "spearman_r"
  return(es)
}
