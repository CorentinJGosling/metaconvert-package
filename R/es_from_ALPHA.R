#' Convert Cronbach's alpha into an effect size measure
#'
#' @param cronbach_alpha Cronbach's alpha reliability coefficient
#' @param n_sample the total number of participants that completed the scale
#' @param n_items number of items in the scale
#' @param cronbach_alpha_se standard error of alpha, \strong{on the natural (coefficient)
#'   scale}. Optional. When supplied it takes precedence over the closed-form
#'   \eqn{(n, k)} standard error.
#' @param cronbach_alpha_ci_lo lower bound of the 95% confidence interval of alpha (natural scale)
#' @param cronbach_alpha_ci_up upper bound of the 95% confidence interval of alpha (natural scale)
#' @param alpha_to_es method used to compute the effect size from Cronbach's alpha.
#'   Must be one of \code{"bonett"} (default), \code{"raw"} or \code{"hakstian_whalen"}.
#'
#' @details
#'
#' \strong{Where the standard error comes from.} Three sources, in this order (the same
#' precedence \code{\link{es_from_omega}} and \code{\link{es_from_icc}} use):
#' \enumerate{
#'   \item \code{cronbach_alpha_se}, read on the natural scale and delta-mapped onto the
#'     analysis scale;
#'   \item otherwise \code{cronbach_alpha_ci_lo} / \code{cronbach_alpha_ci_up}, transformed at
#'     the \strong{bounds} so an asymmetric interval maps correctly instead of being
#'     symmetrised first;
#'   \item otherwise the closed-form \eqn{(n, k)} standard error below.
#' }
#' \code{n_sample} and \code{n_items} gate only the third route, not the point estimate, so a
#' study reporting \emph{"alpha = .88, 95% CI .85 to .91"} without an item count is usable
#' rather than being dropped whole. A row with none of the three keeps its effect size and
#' gets \code{alpha_se = NA}, so it stays visible and countable but is left out of the pool
#' by \code{summary()}.
#'
#' Unlike omega, alpha does have a closed-form variance in \eqn{(n, k)}, so the third route
#' is the usual one. The first two matter because a reported interval is sometimes all a
#' study gives -- and because they are the only way to enter an alpha whose item count was
#' never reported.
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
#' 3. When \code{alpha_to_es = "hakstian_whalen"}, the Hakstian & Whalen (1976)
#' cube-root normalising transformation is applied. This is the transformation
#' Rodriguez & Maeda (2006) recommend for the meta-analysis of coefficient alpha,
#' and the one used by most published reliability-generalization syntheses:
#' \deqn{T(\alpha) = 1 - (1 - \alpha)^{1/3}}
#' \deqn{T\_var = \frac{18 k (n - 1) (1 - \alpha)^{2/3}}{(k - 1)(9n - 11)^2}}
#' It derives from Paulson's (1942) normalising transformation of the \eqn{F}
#' distribution applied to the Feldt (1965) / Kristof (1963) result that
#' \eqn{(1 - r)/(1 - \alpha)} is distributed as \eqn{F}.
#'
#' \strong{Orientation.} Rodriguez & Maeda write \eqn{T = (1 - \alpha)^{1/3}},
#' which \emph{decreases} in alpha; this function instead stores
#' \eqn{1 - (1 - \alpha)^{1/3}}, which \emph{increases} in alpha, matching
#' \code{metafor}'s \code{measure = "AHW"} bit for bit (both \code{yi} and
#' \code{vi}). The two differ by a constant, so the sampling variance is
#' identical, but the increasing orientation means \code{metafor::transf.iahw()}
#' back-transforms this output correctly and moderator coefficients read in the
#' same direction as alpha. Contrast the Bonett scale, where metaConvert and
#' metafor deliberately differ in sign (see \code{\link{reliability_backtransform}}).
#'
#' \strong{Choosing between them.} Bonett's standard error does not involve alpha
#' at all, so two studies with the same \eqn{n} and \eqn{k} receive identical
#' weight whether their alpha is .85 or .98; the Hakstian-Whalen variance does
#' involve alpha and separates them (a 3.8:1 weight ratio at \eqn{n = 200},
#' \eqn{k = 10}). Use \code{"hakstian_whalen"} to reproduce or audit the
#' reliability-generalization method literature, and \code{"bonett"} (the
#' default) otherwise.
#'
#' Both the Bonett and Hakstian-Whalen transformations stabilise the variance and
#' are suitable for meta-analysis; the raw scale is not recommended.
#'
#' @export es_from_cronbach_alpha
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient alpha.
#' Journal of Educational and Behavioral Statistics, 27(4), 335-340.
#'
#' Hakstian, A. R., & Whalen, T. E. (1976). A k-sample significance test for independent
#' alpha coefficients. Psychometrika, 41(2), 219-231.
#'
#' Rodriguez, M. C., & Maeda, Y. (2006). Meta-analysis of coefficient alpha.
#' Psychological Methods, 11(3), 306-322.
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
                                   cronbach_alpha_se, cronbach_alpha_ci_lo,
                                   cronbach_alpha_ci_up,
                                   alpha_to_es = "bonett") {

  len <- length(cronbach_alpha)
  if (missing(n_sample)) n_sample <- rep(NA, len)
  if (missing(n_items)) n_items <- rep(NA, len)
  if (missing(cronbach_alpha_se)) cronbach_alpha_se <- rep(NA_real_, len)
  if (missing(cronbach_alpha_ci_lo)) cronbach_alpha_ci_lo <- rep(NA_real_, len)
  if (missing(cronbach_alpha_ci_up)) cronbach_alpha_ci_up <- rep(NA_real_, len)

  # A reliability-generalization sheet fixes the instrument, so `n_items` (and
  # often `n_sample`) is naturally passed as ONE number beside a vector of
  # coefficients. Without this recycling the subsetting below
  # (`n_sample[nn_miss]`) indexes a length-1 vector at positions 2..k and
  # returns NA, so every row after the first got an NA standard error while its
  # point estimate stayed correct -- a fully populated `alpha` column hiding a
  # pool that rma() then silently drops. Same idiom as es_from_pearson_r().
  if (length(n_sample) == 1) n_sample <- rep(n_sample, length(cronbach_alpha))
  if (length(n_items) == 1) n_items <- rep(n_items, length(cronbach_alpha))
  if (length(n_sample) != length(cronbach_alpha)) stop("The length of the 'n_sample' argument is incorrectly specified.")
  if (length(n_items) != length(cronbach_alpha)) stop("The length of the 'n_items' argument is incorrectly specified.")
  rec <- function(v, nm) {
    if (length(v) == 1) return(rep(v, len))
    if (length(v) != len) {
      stop(paste0("The length of the '", nm, "' argument is incorrectly specified."))
    }
    v
  }
  cronbach_alpha_se <- rec(cronbach_alpha_se, "cronbach_alpha_se")
  cronbach_alpha_ci_lo <- rec(cronbach_alpha_ci_lo, "cronbach_alpha_ci_lo")
  cronbach_alpha_ci_up <- rec(cronbach_alpha_ci_up, "cronbach_alpha_ci_up")

  # a non-positive dispersion cannot be a standard error; a transposed interval is the
  # same interval (R/internal_guards.R)
  cronbach_alpha_se <- .positive_or_na(cronbach_alpha_se)
  ci_lo <- .ci_lower(cronbach_alpha_ci_lo, cronbach_alpha_ci_up)
  ci_up <- .ci_upper(cronbach_alpha_ci_lo, cronbach_alpha_ci_up)

  if (!alpha_to_es %in% c("bonett", "raw", "hakstian_whalen")) {
    stop(paste0("'", alpha_to_es, "' not in tolerated values for the 'alpha_to_es' argument. ",
                "Possible inputs are: 'bonett', 'raw', 'hakstian_whalen'"))
  }

  # P16: per-element guards -- alpha > 1 is impossible (V11 upper bound),
  # n <= 2 makes the SE denominator (n - 2) non-positive, and a single item
  # (k < 2) leaves alpha undefined. alpha = 1 exactly is a valid boundary on
  # the raw scale (es = 1, se = 0) but has no Bonett transform (log(0)).
  # Direct calls degrade to NA instead of emitting Inf/NaN (the pipeline
  # already NAs the impossible values via Tier-1 validation).
  # Forward transform and the delta map that carries a NATURAL-scale reported SE onto
  # it -- the same pair es_from_omega() uses, since alpha and omega share these scales.
  fwd <- function(a) {
    if (identical(alpha_to_es, "bonett")) return(log(1 - a))
    if (identical(alpha_to_es, "hakstian_whalen")) return(1 - (1 - a)^(1 / 3))
    a
  }
  se_map <- function(a, se) {
    if (identical(alpha_to_es, "bonett")) return(se / (1 - a))
    if (identical(alpha_to_es, "hakstian_whalen")) return(se / (3 * (1 - a)^(2 / 3)))
    se
  }

  # alpha > 1 is impossible (V11 upper bound); alpha = 1 exactly is a valid boundary on
  # the raw scale (es = 1, se = 0) but has no Bonett transform (log(0)). n_sample and
  # n_items are NOT part of this test any more: they gate only the COMPUTED standard
  # error (route 3 below), not the point estimate, so a study reporting
  # "alpha = .88, 95% CI [.85, .91]" without an item count is now usable instead of
  # being dropped whole -- the gap roadmap item 2.2 exists for.
  valid_es <- !is.na(cronbach_alpha) & cronbach_alpha <= 1 &
    !(cronbach_alpha == 1 & alpha_to_es == "bonett")

  n <- length(cronbach_alpha)
  alpha_es <- rep(NA_real_, n)
  alpha_es_se <- rep(NA_real_, n)
  alpha_es[valid_es] <- fwd(cronbach_alpha[valid_es])

  # --- standard error, in order of preference -------------------------------
  # 1. reported SE, read on the NATURAL (coefficient) scale and delta-mapped
  from_se <- which(valid_es & !is.na(cronbach_alpha_se))
  if (length(from_se)) {
    alpha_es_se[from_se] <- se_map(cronbach_alpha[from_se], cronbach_alpha_se[from_se])
  }

  # 2. otherwise the reported CI, transformed at the BOUNDS so an asymmetric interval
  #    maps correctly rather than being symmetrised first. A bound of exactly 1 has no
  #    log(1 - .) or (1 - .)^(1/3), so it is excluded on those scales only.
  ci_usable <- if (identical(alpha_to_es, "raw")) rep(TRUE, n) else (ci_lo < 1 & ci_up < 1)
  ci_usable[is.na(ci_usable)] <- FALSE
  from_ci <- which(valid_es & is.na(alpha_es_se) &
                     !is.na(ci_lo) & !is.na(ci_up) & ci_usable)
  if (length(from_ci)) {
    alpha_es_se[from_ci] <-
      abs(fwd(ci_up[from_ci]) - fwd(ci_lo[from_ci])) / (2 * qnorm(0.975))
  }

  # 3. otherwise the closed-form (n, k) SE. n <= 2 makes the denominator (n - 2)
  #    non-positive and a single item leaves alpha undefined, so both gate this route.
  nn_miss <- which(valid_es & is.na(alpha_es_se) &
                     !is.na(n_sample) & n_sample > 2 &
                     !is.na(n_items) & n_items >= 2)

  if (length(nn_miss) != 0) {
    a <- cronbach_alpha[nn_miss]
    ns <- n_sample[nn_miss]
    k <- n_items[nn_miss]

    if (alpha_to_es == "bonett") {
      alpha_es_se[nn_miss] <- sqrt(2 * k / ((k - 1) * (ns - 2)))
    } else if (alpha_to_es == "hakstian_whalen") {
      # Hakstian & Whalen (1976), the transform Rodriguez & Maeda (2006)
      # recommend for meta-analysis of alpha. Their T is (1 - alpha)^(1/3),
      # which DECREASES in alpha; we store metafor's orientation
      # 1 - (1 - alpha)^(1/3), which increases, so that metafor::transf.iahw()
      # inverts our output correctly. Shifting by a constant leaves the
      # variance untouched, so the Rodriguez & Maeda variance below applies to
      # either orientation unchanged (verified bit-exact against
      # metafor::escalc(measure = "AHW")).
      alpha_es_se[nn_miss] <- sqrt(
        18 * k * (ns - 1) * (1 - a)^(2 / 3) / ((k - 1) * (9 * ns - 11)^2)
      )
    } else {
      alpha_es_se[nn_miss] <- (1 - a) * sqrt(2 * k / ((k - 1) * (ns - 2)))
    }
  }
  # NO .positive_or_na() on the result. The REPORTED se is already guarded at the top,
  # and on the raw scale alpha = 1 legitimately yields se = 0 -- a documented boundary
  # (es = 1, se = 0) that V40 deliberately stays silent about. Guarding here NA'd it,
  # which drops the row from summary() entirely.

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
