#' Convert McDonald's omega into an effect size measure
#'
#' @param omega McDonald's omega reliability coefficient (natural scale, usually between 0 and 1)
#' @param omega_se standard error of omega, \strong{on the natural (coefficient) scale}
#' @param omega_ci_lo lower bound of the 95% confidence interval of omega (natural scale)
#' @param omega_ci_up upper bound of the 95% confidence interval of omega (natural scale)
#' @param n_sample the total number of participants that completed the scale. Not used to
#'   compute the standard error (see \strong{Details}); carried through for reporting and
#'   moderator analyses.
#' @param n_items number of items in the scale. Carried through for reporting; not used in
#'   the standard error.
#' @param omega_type which omega is reported: \code{"total"} (default), \code{"hierarchical"},
#'   \code{"asymptotic"} or \code{"subscale"}. These are \strong{different estimands} and must
#'   not be pooled together (see \strong{Details}).
#' @param omega_to_es method used to compute the effect size from omega. One of
#'   \code{"bonett"} (default), \code{"raw"} or \code{"hakstian_whalen"}, matching
#'   \code{\link{es_from_cronbach_alpha}}.
#'
#' @details
#' \strong{Why this function needs a standard error or a confidence interval, and
#' \code{es_from_cronbach_alpha()} does not.} Coefficient alpha has a closed-form sampling
#' variance in \eqn{(n, k)} alone -- Bonett's \eqn{2k/((k-1)(n-2))} and the
#' Hakstian-Whalen variance both descend from the Feldt (1965) / Kristof (1963) result
#' that \eqn{(1 - r)/(1 - \alpha)} is distributed as \eqn{F}, which holds under
#' \strong{essential tau-equivalence} (equal true-score loadings), compound symmetry and
#' multivariate normality. Omega exists precisely to \emph{drop} tau-equivalence: it is
#' defined for a congeneric model with freely estimated loadings. Borrowing either alpha
#' variance for omega would therefore be a category error, and there is no (n, k)-only
#' replacement: the asymptotic variance of omega depends on the full covariance matrix of
#' the estimated loadings and error variances (Raykov, 2002), which primary studies do not
#' print, and the recommended interval is a bootstrap (BCa) one computed from raw data
#' (Kelley & Pornprasertmanit, 2016).
#'
#' The standard error is therefore taken from what the primary study reported, in this
#' order:
#'
#' \enumerate{
#'   \item \code{omega_se}, interpreted on the \strong{natural} scale and mapped to the
#'     analysis scale by the delta method:
#'     \deqn{SE_T = SE_\omega / (1 - \omega)}  (bonett)
#'     \deqn{SE_T = SE_\omega / (3 (1 - \omega)^{2/3})}  (hakstian_whalen)
#'   \item otherwise the reported confidence interval, whose \strong{bounds} are transformed
#'     and differenced. Transforming the bounds rather than the point estimate is what makes
#'     this correct for the asymmetric bootstrap intervals omega usually carries, since a
#'     BCa interval is transformation-respecting.
#'   \item otherwise the effect size is returned with \code{omega_se = NA}. The row is then
#'     visible and countable in the pipeline (it appears in the extraction sheet and in the
#'     reporting-rate accounting) but is dropped by \code{summary()} and cannot be pooled --
#'     which is the honest outcome for a bare omega with no uncertainty attached.
#' }
#'
#' \strong{Lossiness.} Step 2 collapses a possibly asymmetric interval into a single
#' standard error, from which a symmetric Wald interval is rebuilt (the convention used by
#' every other \code{*_ci} route in this package). The point estimate and the interval
#' width are preserved; the asymmetry is not.
#'
#' \strong{omega_type is an estimand, not a label.} \code{"total"} is the proportion of
#' total score variance attributable to all common factors; \code{"hierarchical"} is the
#' proportion attributable to the \emph{general} factor only, and is systematically smaller.
#' They answer different questions and must not share a pool. \code{convert_df()} flags a
#' pool that mixes them.
#'
#' \strong{Do not pool omega with alpha either.} Under a congeneric model
#' \eqn{\omega_{total} \ge \alpha}, with equality only under exact tau-equivalence, so the
#' two are different quantities and their average is neither. There is a second, subtler
#' reason specific to this package: alpha's standard errors here are analytic functions of
#' \eqn{(n, k)} while omega's come from reported (often bootstrap) intervals, so the
#' \emph{weighting scheme} would differ systematically by metric, making metric perfectly
#' confounded with precision. Run two pools and report them side by side; if a study
#' reports both, estimate the gap as a within-study paired contrast rather than by
#' differencing two separately pooled means (the two coefficients are computed on the same
#' data and are near-perfectly dependent).
#'
#' \strong{Heywood cases.} An estimated \eqn{\omega \ge 1} arises from an improper solution
#' (a negative estimated error variance) and is not simply a transcription error. It has no
#' Bonett or Hakstian-Whalen transform, so those scales return \code{NA}; \code{"raw"}
#' returns the value so it stays inspectable.
#'
#' @export es_from_omega
#'
#' @references
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for population
#' reliability coefficients: Evaluation of methods, recommendations, and software for
#' composite measures. Psychological Methods, 21(1), 69-92.
#'
#' McDonald, R. P. (1999). Test theory: A unified treatment. Lawrence Erlbaum.
#'
#' Raykov, T. (2002). Analytic estimation of standard error and confidence interval for
#' scale reliability. Multivariate Behavioral Research, 37(1), 89-103.
#'
#' @md
#'
#' @return
#' This function estimates the standard error of McDonald's omega.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab omega\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab N/A\cr
#'  \tab \cr
#'  \code{required input data} \tab omega + (omega_se OR omega_ci_lo & omega_ci_up)\cr
#'  \tab \cr
#' }
#'
#' @examples
#' # from a reported standard error (natural scale)
#' es_from_omega(omega = 0.86, omega_se = 0.021, n_sample = 300, n_items = 12)
#'
#' # from a reported (possibly asymmetric) confidence interval
#' es_from_omega(omega = 0.86, omega_ci_lo = 0.81, omega_ci_up = 0.90,
#'               n_sample = 300, n_items = 12)
#'
#' # a bare omega: kept and visible, but with no SE it cannot be pooled
#' es_from_omega(omega = 0.86)
es_from_omega <- function(omega, omega_se, omega_ci_lo, omega_ci_up,
                          n_sample, n_items, omega_type = "total",
                          omega_to_es = "bonett") {

  len <- length(omega)
  if (missing(omega_se)) omega_se <- rep(NA_real_, len)
  if (missing(omega_ci_lo)) omega_ci_lo <- rep(NA_real_, len)
  if (missing(omega_ci_up)) omega_ci_up <- rep(NA_real_, len)
  if (missing(n_sample)) n_sample <- rep(NA_real_, len)
  if (missing(n_items)) n_items <- rep(NA_real_, len)
  if (missing(omega_type)) omega_type <- rep("total", len)

  # Same recycling contract as es_from_cronbach_alpha(): a length-1 argument beside a
  # vector of coefficients is the natural way to write an extraction sheet, and without
  # this the row-index subsetting below would silently NA every row after the first.
  rec <- function(v, nm) {
    if (length(v) == 1) v <- rep(v, len)
    if (length(v) != len)
      stop(paste0("The length of the '", nm, "' argument is incorrectly specified."))
    v
  }
  omega_se <- rec(omega_se, "omega_se")
  omega_ci_lo <- rec(omega_ci_lo, "omega_ci_lo")
  omega_ci_up <- rec(omega_ci_up, "omega_ci_up")
  n_sample <- rec(n_sample, "n_sample")
  n_items <- rec(n_items, "n_items")
  omega_type <- rec(omega_type, "omega_type")
  omega_type[is.na(omega_type)] <- "total"

  if (!all(omega_type %in% c("total", "hierarchical", "asymptotic", "subscale"))) {
    stop(paste0("'", paste(unique(omega_type[!omega_type %in%
                c("total", "hierarchical", "asymptotic", "subscale")]), collapse = "', '"),
                "' not in tolerated values for the 'omega_type' argument. ",
                "Possible inputs are: 'total', 'hierarchical', 'asymptotic', 'subscale'"))
  }
  if (!omega_to_es %in% c("bonett", "raw", "hakstian_whalen")) {
    stop(paste0("'", omega_to_es, "' not in tolerated values for the 'omega_to_es' argument. ",
                "Possible inputs are: 'bonett', 'raw', 'hakstian_whalen'"))
  }

  # A non-positive dispersion cannot be a standard error; a transposed interval is the
  # same interval (R/internal_guards.R).
  omega_se <- .positive_or_na(omega_se)
  ci_lo <- .ci_lower(omega_ci_lo, omega_ci_up)
  ci_up <- .ci_upper(omega_ci_lo, omega_ci_up)

  # forward transform and its derivative-based SE map
  fwd <- function(w) {
    if (identical(omega_to_es, "bonett")) return(log(1 - w))
    if (identical(omega_to_es, "hakstian_whalen")) return(1 - (1 - w)^(1 / 3))
    w
  }
  se_map <- function(w, se) {
    if (identical(omega_to_es, "bonett")) return(se / (1 - w))
    if (identical(omega_to_es, "hakstian_whalen")) return(se / (3 * (1 - w)^(2 / 3)))
    se
  }

  # omega >= 1 (a Heywood case) has no bonett/hakstian_whalen transform; keep it on the
  # raw scale so it stays inspectable rather than vanishing.
  transformable <- if (identical(omega_to_es, "raw")) {
    !is.na(omega)
  } else {
    !is.na(omega) & omega < 1
  }

  es <- rep(NA_real_, len)
  se <- rep(NA_real_, len)
  es[transformable] <- fwd(omega[transformable])

  # 1. reported SE (natural scale) -> delta method
  from_se <- which(transformable & !is.na(omega_se))
  if (length(from_se)) se[from_se] <- se_map(omega[from_se], omega_se[from_se])

  # 2. otherwise the reported CI, transformed at the BOUNDS so an asymmetric bootstrap
  #    interval maps correctly (BCa is transformation-respecting)
  from_ci <- which(transformable & is.na(se) &
                     !is.na(ci_lo) & !is.na(ci_up) & ci_lo < 1 & ci_up < 1)
  if (length(from_ci)) {
    t_lo <- fwd(ci_lo[from_ci])
    t_up <- fwd(ci_up[from_ci])
    se[from_ci] <- abs(t_up - t_lo) / (2 * qnorm(0.975))
  }
  se <- .positive_or_na(se)

  omega_ci_lo_out <- es - qnorm(0.975) * se
  omega_ci_up_out <- es + qnorm(0.975) * se

  data.frame(
    omega = es,
    omega_se = se,
    omega_ci_lo = omega_ci_lo_out,
    omega_ci_up = omega_ci_up_out,
    n_sample = n_sample,
    n_items = n_items,
    omega_type = omega_type,
    info_used = "omega"
  )
}
