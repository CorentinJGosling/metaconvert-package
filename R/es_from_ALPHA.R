#' Convert Cronbach's alpha into an effect size measure
#'
#' @param cronbach_alpha Cronbach's alpha reliability coefficient
#' @param n_sample the total number of participants that completed the scale
#' @param n_items number of items in the scale
#' @param alpha_to_es method used to compute the effect size from Cronbach's alpha.
#'   Must be one of \code{"bonett"} (default), \code{"raw"} or \code{"hakstian_whalen"}.
#'
#' @details
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
#' \deqn{T(lpha) = 1 - (1 - lpha)^{1/3}}
#' \deqn{T\_var = rac{18 k (n - 1) (1 - lpha)^{2/3}}{(k - 1)(9n - 11)^2}}
#' It derives from Paulson's (1942) normalising transformation of the \eqn{F}
#' distribution applied to the Feldt (1965) / Kristof (1963) result that
#' \eqn{(1 - r)/(1 - lpha)} is distributed as \eqn{F}.
#'
#' \strong{Orientation.} Rodriguez & Maeda write \eqn{T = (1 - lpha)^{1/3}},
#' which \emph{decreases} in alpha; this function instead stores
#' \eqn{1 - (1 - lpha)^{1/3}}, which \emph{increases} in alpha, matching
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
                                   alpha_to_es = "bonett") {

  if (missing(n_sample)) n_sample <- rep(NA, length(cronbach_alpha))
  if (missing(n_items)) n_items <- rep(NA, length(cronbach_alpha))

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
  invalid <- (!is.na(cronbach_alpha) & cronbach_alpha > 1) |
    (!is.na(cronbach_alpha) & cronbach_alpha == 1 & alpha_to_es == "bonett") |
    (!is.na(n_sample) & n_sample <= 2) |
    (!is.na(n_items) & n_items < 2)

  nn_miss <- which(!is.na(cronbach_alpha) & !is.na(n_sample) & !is.na(n_items) &
                     !invalid)

  n <- length(cronbach_alpha)
  alpha_es <- rep(NA_real_, n)
  alpha_es_se <- rep(NA_real_, n)

  if (length(nn_miss) != 0) {
    a <- cronbach_alpha[nn_miss]
    ns <- n_sample[nn_miss]
    k <- n_items[nn_miss]

    if (alpha_to_es == "bonett") {
      alpha_es[nn_miss] <- log(1 - a)
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
      alpha_es[nn_miss] <- 1 - (1 - a)^(1 / 3)
      alpha_es_se[nn_miss] <- sqrt(
        18 * k * (ns - 1) * (1 - a)^(2 / 3) / ((k - 1) * (9 * ns - 11)^2)
      )
    } else {
      alpha_es[nn_miss] <- a
      alpha_es_se[nn_miss] <- (1 - a) * sqrt(2 * k / ((k - 1) * (ns - 2)))
    }
  }

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
