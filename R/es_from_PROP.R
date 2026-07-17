#' Convert single-group proportion into an effect size measure
#'
#' @param prop proportion of cases in the single group (0-1)
#' @param n_sample total sample size
#' @param prop_to_es method used to compute the effect size from the proportion.
#'   Must be \code{"raw"} (default), \code{"logit"}, or \code{"freeman_tukey"}.
#' @param reverse_prop a logical value indicating whether the direction of the proportion should be flipped.
#'
#' @details
#' This function computes an effect size from a single-group proportion.
#'
#' 1. **When \code{prop_to_es = "raw"}** (default), the proportion is used as the effect size
#' and its standard error is:
#' \deqn{prop\_se = \sqrt{\frac{prop \times (1 - prop)}{n}}}
#'
#' 2. **When \code{prop_to_es = "logit"}**, the proportion is converted to a log odds:
#' \deqn{logit = \log\left(\frac{prop}{1 - prop}\right)}
#' \deqn{logit\_se = \sqrt{\frac{1}{n \times prop \times (1 - prop)}}}
#'
#' 3. **When \code{prop_to_es = "freeman_tukey"}**, the Freeman-Tukey double arcsine
#' transformation is applied (Barendregt et al., 2013):
#' \deqn{FT = \frac{1}{2}\left[\arcsin\left(\sqrt{\frac{x}{n+1}}\right) + \arcsin\left(\sqrt{\frac{x+1}{n+1}}\right)\right]}
#' with \eqn{x = n \times prop}, and:
#' \deqn{FT\_se = \frac{1}{\sqrt{4n + 2}}}
#'
#' When prop is equal to 0 or 1, a 0.5 correction is applied for the raw and logit methods:
#' \eqn{prop\_corrected = \frac{x + 0.5}{n + 1}}. For the raw method this shifts the
#' reported point estimate itself (a boundary proportion of exactly 0 or 1 is
#' returned as \eqn{(x + 0.5)/(n + 1)}, i.e. nudged toward the interior), matching
#' the convention of \code{metafor}'s \code{measure = "PR"}. The corrected
#' denominator \eqn{n + 1} also replaces \eqn{n} in the raw and logit standard
#' errors for these boundary rows (i.e.,
#' \eqn{\sqrt{p_c (1 - p_c) / (n + 1)}} and \eqn{\sqrt{1 / ((n + 1) \times p_c \times (1 - p_c))}}),
#' matching \code{metafor}'s \code{"PR"}/\code{"PLO"} convention throughout.
#'
#' Proportions outside \eqn{[0, 1]} are set to NA with a warning (the affected
#' rows return NA, but do not abort the run).
#'
#' @references
#' Barendregt, J. J., Doi, S. A., Lee, Y. Y., Norman, R. E., & Vos, T. (2013).
#' Meta-analysis of prevalence. Journal of Epidemiology and Community Health, 67(11), 974-978.
#'
#' Miller, J. J. (1978). The inverse of the Freeman-Tukey double arcsine transformation.
#' The American Statistician, 32(4), 138-138. (Back-transformation left to the
#' user; not applied here: Freeman-Tukey results stay on the transformed scale.)
#'
#' @return
#' This function estimates a single-group proportion.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab PROP\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab N/A\cr
#'  \tab \cr
#'  \code{required input data} \tab prop + n_sample\cr
#'  \tab \cr
#' }
#'
#' Note: when \code{prop_to_es = "logit"} or \code{"freeman_tukey"}, the returned
#' \code{prop}, \code{prop_se}, \code{prop_ci_lo} and \code{prop_ci_up} columns are
#' on the transformed (log-odds / Freeman-Tukey double-arcsine) scale, \strong{not}
#' the `[0, 1]` proportion scale, and are not back-transformed. Only
#' \code{prop_to_es = "raw"} returns values on the proportion scale.
#'
#' @export es_from_prop_single_group
#'
#' @md
#'
#' @examples
#' es_from_prop_single_group(prop = 0.30, n_sample = 100)
es_from_prop_single_group <- function(prop, n_sample, prop_to_es = "raw", reverse_prop) {
  if (missing(reverse_prop)) reverse_prop <- rep(FALSE, length(prop))
  reverse_prop[is.na(reverse_prop)] <- FALSE
  if (length(reverse_prop) == 1) reverse_prop <- rep(reverse_prop, length(prop))
  if (length(reverse_prop) != length(prop)) {
    stop("The length of the 'reverse_prop' argument is incorrectly specified.")
  }

  if (!prop_to_es %in% c("raw", "logit", "freeman_tukey")) {
    stop(paste0("'", prop_to_es, "' not in tolerated values for the 'prop_to_es' argument. ",
                "Possible inputs are: 'raw', 'logit', 'freeman_tukey'"))
  }

  # P9: per-row NA + warning instead of a hard stop. convert_df() calls this
  # function for every measure with no tryCatch, so a stop() here would abort
  # an entire pipeline run (including unrelated measures) because of one bad
  # cell. Tier-1 validation reports the issue in both correct_inputs modes.
  out_of_range <- !is.na(prop) & (prop < 0 | prop > 1)
  if (any(out_of_range)) {
    warning(sprintf(
      "es_from_prop_single_group: %d row(s) had a proportion outside [0, 1]; these proportions were set to NA.",
      sum(out_of_range)
    ), call. = FALSE)
    prop[out_of_range] <- NA_real_
  }

  nn_miss <- which(!is.na(prop) & !is.na(n_sample))

  n <- length(prop)
  prop_es <- rep(NA_real_, n)
  prop_es_se <- rep(NA_real_, n)

  if (length(nn_miss) != 0) {
    p <- prop[nn_miss]
    ns <- n_sample[nn_miss]

    p <- ifelse(reverse_prop[nn_miss], 1 - p, p)

    # correcting by 0.5 when prop = 0 or 1
    n_cases <- p * ns
    needs_correction <- (p == 0 | p == 1)
    p_corrected <- p
    p_corrected[needs_correction] <- (n_cases[needs_correction] + 0.5) / (ns[needs_correction] + 1)
    # P5: the continuity-corrected denominator n + 1 must also flow into the
    # SE for boundary rows (metafor PR/PLO convention: vi = 1/x_c + 1/(n_c - x_c)
    # with x_c = x + 0.5, n_c = n + 1); interior rows keep the raw n.
    ns_corrected <- ns
    ns_corrected[needs_correction] <- ns[needs_correction] + 1

    if (prop_to_es == "raw") {
      prop_es[nn_miss] <- p_corrected
      prop_es_se[nn_miss] <- sqrt(p_corrected * (1 - p_corrected) / ns_corrected)

    } else if (prop_to_es == "logit") {
      prop_es[nn_miss] <- log(p_corrected / (1 - p_corrected))
      prop_es_se[nn_miss] <- sqrt(1 / (ns_corrected * p_corrected * (1 - p_corrected)))

    } else if (prop_to_es == "freeman_tukey") {
      # no correction for FT
      x_ft <- ns * p
      prop_es[nn_miss] <- 0.5 * (asin(sqrt(x_ft / (ns + 1))) +
                                  asin(sqrt((x_ft + 1) / (ns + 1))))
      prop_es_se[nn_miss] <- 1 / sqrt(4 * ns + 2)
    }
  }

  # CIs
  prop_ci_lo <- prop_es - qnorm(0.975) * prop_es_se
  prop_ci_up <- prop_es + qnorm(0.975) * prop_es_se

  # bound to [0,1]
  if (prop_to_es == "raw") {
    prop_ci_lo <- pmax(0, prop_ci_lo)
    prop_ci_up <- pmin(1, prop_ci_up)
  }

  result <- data.frame(
    prop = prop_es,
    prop_se = prop_es_se,
    prop_ci_lo = prop_ci_lo,
    prop_ci_up = prop_ci_up,
    n_sample = n_sample,
    info_used = "prop_single_group"
  )

  return(result)
}


#' Convert single-group case counts into a proportion effect size measure
#'
#' @param n_cases number of cases in the single group
#' @param n_sample total sample size
#' @param prop_to_es method used to compute the effect size from the proportion.
#'   Must be \code{"raw"} (default), \code{"logit"}, or \code{"freeman_tukey"}.
#' @param reverse_prop a logical value indicating whether the direction of the proportion should be flipped.
#'
#' @details
#' This is a convenience function that converts case counts to proportions
#' and then calls \code{\link{es_from_prop_single_group}()}.
#'
#' The proportion is calculated as:
#' \deqn{prop = \frac{n\_cases}{n\_sample}}
#'
#' Then, calculations of the \code{\link{es_from_prop_single_group}()} are applied.
#'
#' Rows where \code{n_cases > n_sample} are set to NA (both values) with a
#' warning; the affected rows return NA, but do not abort the run.
#'
#' @references
#' Barendregt, J. J., Doi, S. A., Lee, Y. Y., Norman, R. E., & Vos, T. (2013).
#' Meta-analysis of prevalence. Journal of Epidemiology and Community Health, 67(11), 974-978.
#'
#' @return
#' This function returns the same output as \code{\link{es_from_prop_single_group}()}.
#' See that function's documentation for details.
#'
#' @export es_from_prop_single_group_counts
#'
#' @md
#'
#' @examples
#' es_from_prop_single_group_counts(n_cases = 25, n_sample = 100)
es_from_prop_single_group_counts <- function(n_cases, n_sample, prop_to_es = "raw", reverse_prop) {
  if (missing(reverse_prop)) reverse_prop <- rep(FALSE, length(n_cases))
  reverse_prop[is.na(reverse_prop)] <- FALSE

  # P9: per-row NA + warning instead of a hard stop (see
  # es_from_prop_single_group). Both members of an inconsistent pair are set
  # to NA, mirroring the Tier-1 part-total convention.
  inconsistent <- !is.na(n_cases) & !is.na(n_sample) & n_cases > n_sample
  if (any(inconsistent)) {
    warning(sprintf(
      "es_from_prop_single_group_counts: %d row(s) had n_cases > n_sample; both values were set to NA for these rows.",
      sum(inconsistent)
    ), call. = FALSE)
    n_cases[inconsistent] <- NA_real_
    n_sample[inconsistent] <- NA_real_
  }

  prop <- n_cases / n_sample

  result <- es_from_prop_single_group(prop = prop, n_sample = n_sample,
                                      prop_to_es = prop_to_es, reverse_prop = reverse_prop)

  result$info_used <- "prop_single_group_counts"

  return(result)
}
