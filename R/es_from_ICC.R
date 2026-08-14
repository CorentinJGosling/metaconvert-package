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
#' with the one-way random-model / two-way-consistency leading-order sampling
#' variance (Bonett, 2002; Donner & Eliasziw, 1987):
#' \deqn{T\_se = \sqrt{\frac{2 (1 + (k-1) ICC)^2}{k (k - 1)(n - 1)}}}
#'
#' 2. When \code{icc_to_es = "raw"}, the raw ICC is used and its standard error
#' is obtained by the delta method (\eqn{(1 - ICC)} times the transformed-scale SE).
#'
#' **Scope of the SE formula.** For the two-way consistency ICC(3,1)
#' (\code{icc_type = "consistency"}) this SE is exact at leading order: deriving
#' it from \eqn{F_0 = MSR/MSE} (with degrees of freedom \eqn{n-1} and
#' \eqn{(n-1)(k-1)}) reduces to the same expression, and it still depends on the
#' ICC value (it is \emph{not} \eqn{\rho}-free). For the two-way
#' absolute-agreement ICC(2,1) (\code{icc_type = "agreement"}, the default) the
#' same formula is only a **one-way approximation that assumes negligible
#' between-rater variance**: when raters differ systematically
#' (\eqn{\sigma^2_{rater} > 0}), the ICC(2,1) estimator depends on the
#' between-rater mean square, which has only \eqn{k - 1} degrees of freedom, so
#' its true sampling variance does not shrink at the \eqn{1/n} rate this formula
#' assumes and the reported SE/CI can be markedly anti-conservative (simulation:
#' 95% CI coverage around 0.74-0.76 with moderate rater variance, degrading as
#' \eqn{n} grows). The exact ICC(2,1) variance requires the rater-variance
#' component, which summary data do not report; a per-row informational flag
#' (V31) marks agreement-type rows for this reason. If the raters are known to
#' be exchangeable (negligible rater variance), the approximation is accurate.
#'
#' **Scale note.** Under the default \code{icc_to_es = "bonett"} the returned
#' \code{icc} column and its \code{icc_se} are both on the \eqn{\ln(1 - ICC)} scale:
#' the example below returns \code{icc = -1.609} for an input ICC of 0.80. To recover
#' the raw-scale SE - which is what the \code{icc_se} argument of
#' \code{\link{compute_sem}} expects - multiply \code{icc_se} by \eqn{1 - ICC} using
#' the ICC value you supplied, not the returned \code{icc} column. Calling
#' \code{es_from_icc} with \code{icc_to_es = "raw"} returns both quantities on the
#' raw scale directly.
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
#'  \code{natural effect size measure} \tab icc, returned as ln(1 - icc) under the default \code{icc_to_es = "bonett"}\cr
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

  # P16: per-element guards -- |icc| > 1 is impossible (V11 bounds), icc = 1
  # has no Bonett transform (log(0)) and a degenerate raw SE, n <= 1 makes the
  # SE denominator (n - 1) non-positive, and k < 2 leaves the ICC undefined.
  # Direct calls degrade to NA instead of emitting Inf/NaN.
  invalid <- (!is.na(icc) & (icc > 1 | icc < -1)) |
    (!is.na(icc) & icc == 1) |
    (!is.na(n_sample) & n_sample <= 1) |
    (!is.na(n_measurements) & n_measurements < 2)

  nn_miss <- which(!is.na(icc) & !is.na(n_sample) & !is.na(n_measurements) &
                     !invalid)

  n <- length(icc)
  icc_es <- rep(NA_real_, n)
  icc_es_se <- rep(NA_real_, n)

  if (length(nn_miss) != 0) {
    rho <- icc[nn_miss]
    ns <- n_sample[nn_miss]
    k <- n_measurements[nn_miss]

    # Bonett (2002) variance-stabilised SE of ln(1 - ICC) for a single-measure ICC.
    # For the two-way consistency ICC(3,1), deriving Var(ln(1 - ICC)) directly from
    # F0 = MSR / MSE (df = n - 1 and (n - 1)(k - 1)) gives
    # [(1 + (k - 1) * rho) / k]^2 * 2k / ((k - 1)(n - 1)), which reduces to the
    # expression below (confirmed by simulation for the consistency case) -- it is
    # NOT independent of rho. For the two-way absolute-agreement ICC(2,1) the same
    # expression is only a one-way approximation that assumes NEGLIGIBLE
    # between-rater variance; with sigma^2_rater > 0 the ICC(2,1) estimator depends
    # on MSC (k - 1 df only) and this SE is anti-conservative (see the roxygen
    # details and the V31 informational flag). The exact ICC(2,1) variance needs
    # the rater-variance component, which summary data do not carry.
    bonett_transformed_se <- sqrt(
      2 * (1 + (k - 1) * rho)^2 / (k * (k - 1) * (ns - 1))
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
