#' Convert means and standard deviations of two independent groups into several effect size measures
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_sd_exp standard deviation of participants in the experimental/exposed group.
#' @param mean_sd_nexp standard deviation of participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
#' @param smd_denom standardizer for the standardized mean difference. "pooled" (default) uses the pooled endpoint SD (Cohen's d / Hedges' g); "glass" (alias "control") uses the control (non-experimental) endpoint SD (Glass's delta); "glass_robust" (alias "control_robust") is Glass's delta with a heteroscedasticity-consistent sampling variance.
#' @param reverse_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first computes a Cohen's d (D), Hedges' g (G) and mean difference (MD)
#' from the means and standard deviations of two independent groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a mean difference**  (formulas 12.1-12.6 in Cooper):
#' \deqn{md = mean\_exp - mean\_nexp}
#' \deqn{md\_se = \sqrt{\frac{mean\_sd\_exp^2}{n\_exp} + \frac{mean\_sd\_nexp^2}{n\_nexp}}}
#' The confidence interval uses the unequal-variance (Welch-Satterthwaite) degrees of freedom
#' that match this standard error (as in \code{stats::t.test(var.equal = FALSE)}):
#' \deqn{df = \frac{(s_e^2/n\_exp + s_n^2/n\_nexp)^2}{(s_e^2/n\_exp)^2/(n\_exp - 1) + (s_n^2/n\_nexp)^2/(n\_nexp - 1)}}
#' \deqn{md\_ci\_lo = md - md\_se * qt(.975, df)}
#' \deqn{md\_ci\_up = md + md\_se * qt(.975, df)}
#'
#' **To estimate a Cohen's d** the following formulas are used (formulas 12.10-12.18 in Cooper):
#' \deqn{mean\_sd\_pooled = \sqrt{\frac{(n\_exp - 1) * sd\_exp^2 + (n\_nexp - 1) * sd\_nexp^2}{n\_exp+n\_nexp-2}}}
#' \deqn{cohen\_d =  \frac{mean\_exp - mean\_nexp}{mean\_sd\_pooled}}
#' \deqn{cohen\_d\_se = \sqrt{\frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2(n\_exp+n\_nexp)}}}
#' \deqn{cohen\_d\_ci\_lo = cohen\_d - cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#' \deqn{cohen\_d\_ci\_up = cohen\_d + cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_means_sd
#'
#' @md
#'
#' @examples
#' es_from_means_sd(
#'   n_exp = 55, n_nexp = 55,
#'   mean_exp = 2.3, mean_sd_exp = 1.2,
#'   mean_nexp = 1.9, mean_sd_nexp = 0.9
#' )
es_from_means_sd <- function(mean_exp, mean_sd_exp, mean_nexp, mean_sd_nexp, n_exp, n_nexp,
                             smd_to_cor = "viechtbauer", smd_var = "borenstein",
                             smd_denom = "pooled", reverse_means) {
  n_row <- length(mean_exp)
  if (missing(reverse_means)) reverse_means <- rep(FALSE, n_row)
  reverse_means[is.na(reverse_means)] <- FALSE
  if (length(reverse_means) == 1) reverse_means = c(rep(reverse_means, n_row))
  if (length(reverse_means) != n_row) stop("The length of the 'reverse_means' argument is incorrectly specified.")

  # smd_denom may be a single value or specified per row (like smd_var / smd_to_cor).
  # A blank/NA cell defaults to "pooled". Each row is standardised by its OWN choice, so a
  # mixed column can no longer silently collapse to the first row's value (previously the
  # match.arg(...[1]) took the first element and applied it to every row).
  smd_denom_in <- as.character(smd_denom)
  if (length(smd_denom_in) == 1) smd_denom_in <- rep(smd_denom_in, n_row)
  if (length(smd_denom_in) != n_row) stop("The length of the 'smd_denom' argument is incorrectly specified.")
  smd_denom_in[is.na(smd_denom_in)] <- "pooled"
  smd_denom <- .normalize_smd_denom(smd_denom_in) # "glass" -> "control", etc.
  if (any(is.na(smd_denom))) {
    stop(paste0(
      "'", paste(unique(smd_denom_in[is.na(smd_denom)]), collapse = "', '"),
      "' not in tolerated values for the 'smd_denom' argument. Possible inputs are: ",
      "'pooled' (the default), 'glass' (alias 'control') or 'glass_robust' (alias 'control_robust')."
    ), call. = FALSE)
  }

  # recycle the other per-row method arguments so they can be subset alongside the rows
  if (length(smd_var) == 1) smd_var <- rep(smd_var, n_row)
  if (length(smd_to_cor) == 1) smd_to_cor <- rep(smd_to_cor, n_row)

  es <- NULL
  for (dn in c("pooled", "control", "control_robust")) {
    sel <- which(smd_denom == dn)
    if (length(sel) == 0) next
    if (dn == "pooled") {
      pooled_sd <- sqrt(((n_exp[sel] - 1) * mean_sd_exp[sel]^2 +
        (n_nexp[sel] - 1) * mean_sd_nexp[sel]^2) / (n_exp[sel] + n_nexp[sel] - 2))
      d <- (mean_exp[sel] - mean_nexp[sel]) / pooled_sd
      es_sel <- .es_from_d(
        d = d, n_exp = n_exp[sel], n_nexp = n_nexp[sel],
        smd_to_cor = smd_to_cor[sel], smd_var = smd_var[sel], reverse = reverse_means[sel]
      )
    } else {
      # Glass's delta: standardise by the control (non-experimental) endpoint SD.
      es_sel <- .glass_from_means(
        mean_exp = mean_exp[sel], mean_sd_exp = mean_sd_exp[sel],
        mean_nexp = mean_nexp[sel], mean_sd_nexp = mean_sd_nexp[sel],
        n_exp = n_exp[sel], n_nexp = n_nexp[sel],
        smd_to_cor = smd_to_cor[sel], smd_var = smd_var[sel],
        reverse = reverse_means[sel], robust = (dn == "control_robust")
      )
    }
    es_sel$.row_order <- sel
    es <- if (is.null(es)) es_sel else rbind(es, es_sel)
  }
  es <- es[order(es$.row_order), , drop = FALSE]
  es$.row_order <- NULL
  rownames(es) <- NULL

  es$info_used <- "means_sd"

  es$md <- ifelse(reverse_means, mean_nexp - mean_exp, mean_exp - mean_nexp)
  es$md <- ifelse(!is.na(es$md) & !is.na(mean_sd_exp) & !is.na(mean_sd_nexp) &
    !is.na(n_exp) & !is.na(n_nexp),
  es$md, NA
  )
  # Unpooled (Welch) SE for the mean difference; pair it with the Welch-Satterthwaite df
  # (not the pooled n1 + n2 - 2) so the CI is coherent with the unequal-variance SE. This
  # matches stats::t.test(var.equal = FALSE). Falls back to the pooled df when the
  # Welch df is undefined (e.g. both arm SDs zero, or an arm n < 2).
  se_e2 <- mean_sd_exp^2 / n_exp
  se_n2 <- mean_sd_nexp^2 / n_nexp
  md_df <- (se_e2 + se_n2)^2 / (se_e2^2 / (n_exp - 1) + se_n2^2 / (n_nexp - 1))
  md_df <- ifelse(is.finite(md_df) & md_df > 0, md_df, n_exp + n_nexp - 2)
  es$md_se <- ifelse(!is.na(es$md), sqrt(se_e2 + se_n2), NA)
  es$md_ci_lo <- ifelse(!is.na(es$md), es$md - qt(.975, md_df) * es$md_se, NA)
  es$md_ci_up <- ifelse(!is.na(es$md), es$md + qt(.975, md_df) * es$md_se, NA)

  return(es)
}

#' Convert means and standard errors of two independent groups into several effect size measures
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_se_exp standard error of participants in the experimental/exposed group.
#' @param mean_se_nexp standard error of participants in the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
#' @param smd_denom standardizer for the standardized mean difference. "pooled" (default) uses the pooled endpoint SD (Cohen's d / Hedges' g); "glass" (alias "control") uses the control (non-experimental) endpoint SD (Glass's delta); "glass_robust" (alias "control_robust") is Glass's delta with a heteroscedasticity-consistent sampling variance.
#' @param reverse_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the standard errors of two independent groups into standard deviations,
#' and then relies on the calculations of the \code{\link{es_from_means_sd}()} function.
#'
#' **To convert the standard errors into standard deviations**, the following formula is used.
#' \deqn{mean\_sd\_exp = mean\_se\_exp * \sqrt{n\_exp}}
#' \deqn{mean\_sd\_nexp = mean\_se\_nexp * \sqrt{n\_nexp}}
#' Then, calculations of the \code{\link{es_from_means_sd}()} are applied.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @export es_from_means_se
#'
#' @md
#'
#' @examples
#' es_from_means_se(
#'   mean_exp = 42, mean_se_exp = 11,
#'   mean_nexp = 42, mean_se_nexp = 15,
#'   n_exp = 43, n_nexp = 34
#' )
es_from_means_se <- function(mean_exp, mean_se_exp, mean_nexp, mean_se_nexp, n_exp, n_nexp,
                             smd_to_cor = "viechtbauer", smd_var = "borenstein",
                             smd_denom = "pooled", reverse_means) {
  if (missing(reverse_means)) reverse_means <- rep(FALSE, length(mean_exp))
  reverse_means[is.na(reverse_means)] <- FALSE

  sd_exp <- mean_se_exp * sqrt(n_exp)
  sd_nexp <- mean_se_nexp * sqrt(n_nexp)

  es <- es_from_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = mean_exp, mean_sd_exp = sd_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = sd_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, smd_denom = smd_denom,
    reverse_means = reverse_means
  )

  es$info_used <- "means_se"

  return(es)
}

#' Convert means of two groups and the pooled standard deviation into several effect size measures
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_sd_pooled pooled standard deviation across both groups.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
#' @param reverse_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function first computes a Cohen's d (D), Hedges' g (G) and mean difference (MD)
#' from the means of two independent groups and the pooled standard deviation across the groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a mean difference**  (formulas 12.1-12.6 in Cooper):
#' \deqn{md = mean\_exp - mean\_nexp}
#' \deqn{md\_se = \sqrt{\frac{n\_exp+n\_nexp}{n\_exp*n\_nexp} * mean\_sd\_pooled^2}}
#' \deqn{md\_ci\_lo = md - md\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#' \deqn{md\_ci\_up = md + md\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#'
#' **To estimate a Cohen's d** the following formulas are used (formulas 12.10-12.18 in Cooper):
#' \deqn{cohen\_d =  \frac{mean\_exp - mean\_nexp}{mean\_sd\_pooled}}
#' \deqn{cohen\_d\_se = \sqrt{\frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2(n\_exp+n\_nexp)}}}
#' \deqn{cohen\_d\_ci\_lo = cohen\_d - cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#' \deqn{cohen\_d\_ci\_up = cohen\_d + cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_means_sd_pooled
#'
#' @md
#'
#' @examples
#' es_from_means_sd_pooled(
#'   n_exp = 55, n_nexp = 55,
#'   mean_exp = 2.3, mean_nexp = 1.9,
#'   mean_sd_pooled = 0.9
#' )
es_from_means_sd_pooled <- function(mean_exp, mean_nexp, mean_sd_pooled, n_exp, n_nexp,
                                    smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_means) {
  if (missing(reverse_means)) reverse_means <- rep(FALSE, length(mean_exp))
  reverse_means[is.na(reverse_means)] <- FALSE
  if (length(reverse_means) == 1) reverse_means = c(rep(reverse_means, length(mean_exp)))
  if (length(reverse_means) != length(mean_exp)) stop("The length of the 'reverse_means' argument is incorrectly specified.")

  # Unlike es_from_means_sd(), which pools the two arm SDs and so squares them on
  # the way, this route divides by the supplied SD directly: a negative one
  # sign-flips d while md keeps its own sign. See R/internal_guards.R.
  mean_sd_pooled <- .positive_or_na(mean_sd_pooled)

  d <- (mean_exp - mean_nexp) / mean_sd_pooled

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse = reverse_means
  )

  es$info_used <- "means_sd_pooled"

  es$md <- ifelse(reverse_means, mean_nexp - mean_exp, mean_exp - mean_nexp)
  es$md <- ifelse(!is.na(es$md) & !is.na(mean_sd_pooled) &
    !is.na(n_exp) & !is.na(n_nexp),
  es$md, NA
  )
  es$md_se <- ifelse(!is.na(es$md), sqrt((n_exp + n_nexp) / (n_exp * n_nexp) * mean_sd_pooled^2), NA)
  es$md_ci_lo <- ifelse(!is.na(es$md), es$md - qt(.975, n_exp + n_nexp - 2) * es$md_se, NA)
  es$md_ci_up <- ifelse(!is.na(es$md), es$md + qt(.975, n_exp + n_nexp - 2) * es$md_se, NA)

  return(es)
}

#' Convert means and 95% CI of two independent groups into several effect size measures
#'
#' @param mean_exp mean of participants in the experimental/exposed group.
#' @param mean_nexp mean of participants in the non-experimental/non-exposed group.
#' @param mean_ci_lo_exp lower bound of the 95% CI of the mean of the experimental/exposed group
#' @param mean_ci_up_exp upper bound of the 95% CI of the mean of the experimental/exposed group
#' @param mean_ci_lo_nexp lower bound of the 95% CI of the mean of the non-experimental/non-exposed group.
#' @param mean_ci_up_nexp upper bound of the 95% CI of the mean of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
#' @param smd_denom standardizer for the standardized mean difference. "pooled" (default) uses the pooled endpoint SD (Cohen's d / Hedges' g); "glass" (alias "control") uses the control (non-experimental) endpoint SD (Glass's delta); "glass_robust" (alias "control_robust") is Glass's delta with a heteroscedasticity-consistent sampling variance.
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param reverse_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the 95% CI of two independent groups into a standard error,
#' and then relies on the calculations of the \code{\link{es_from_means_se}()} function.
#'
#' **To convert the 95% CIs into standard errors,** the following formula is used (table 12.3 in Cooper):
#' \deqn{mean\_se\_exp = \frac{mean\_ci\_up\_exp - mean\_ci\_lo\_exp}{2 * qt{(0.975, df = n\_exp - 1)}}}
#' \deqn{mean\_se\_nexp = \frac{mean\_ci\_up\_nexp - mean\_ci\_lo\_nexp}{2 * qt{(0.975, df = n\_nexp - 1)}}}
#' Calculations of the \code{\link{es_from_means_se}()} are then applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @export es_from_means_ci
#'
#' @md
#'
#' @examples
#' es_from_means_ci(
#'   n_exp = 55, n_nexp = 55,
#'   mean_exp = 25, mean_ci_lo_exp = 15, mean_ci_up_exp = 35,
#'   mean_nexp = 18, mean_ci_lo_nexp = 12, mean_ci_up_nexp = 24
#' )
es_from_means_ci <- function(mean_exp, mean_ci_lo_exp, mean_ci_up_exp,
                             mean_nexp, mean_ci_lo_nexp, mean_ci_up_nexp,
                             n_exp, n_nexp, smd_to_cor = "viechtbauer",
                             smd_var = "borenstein", smd_denom = "pooled",
                             max_asymmetry = 10,
                             reverse_means) {
  if (missing(reverse_means)) reverse_means <- rep(FALSE, length(mean_exp))
  reverse_means[is.na(reverse_means)] <- FALSE

  df_exp <- n_exp - 1
  df_nexp <- n_nexp - 1

  se_exp <- (mean_ci_up_exp - mean_ci_lo_exp) / (2 * qt(0.975, df_exp))
  se_nexp <- (mean_ci_up_nexp - mean_ci_lo_nexp) / (2 * qt(0.975, df_nexp))

  es <- es_from_means_se(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = mean_exp,
    mean_se_exp = se_exp,
    mean_nexp = mean_nexp, mean_se_nexp = se_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, smd_denom = smd_denom,
    reverse_means = reverse_means
  )

  es$info_used <- "means_ci"

  return(es)
}
