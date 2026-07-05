#' Convert a risk difference value and its standard error into several effect size measures
#'
#' @param rd risk difference value (control risk minus treatment risk)
#' @param rd_se standard error of the risk difference
#' @param n_exp number of participants in the exposed/treatment group
#' @param n_nexp number of participants in the non-exposed/control group
#' @param n_cases number of cases/events across both groups
#' @param n_controls number of controls/no-event across both groups
#' @param n_sample total number of participants in the sample
#' @param baseline_risk proportion of cases in the non-exposed/control group.
#'   Required for converting RD to OR, RR, and SMD measures.
#' @param reverse_rd a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts a risk difference (RD) and its standard error into an odds ratio (OR),
#' risk ratio (RR), number needed to treat (NNT), Cohen's d (D), Hedges' g (G),
#' and correlation coefficients (R/Z).
#'
#' **NNT is always computed from RD:**
#' \deqn{nnt = \frac{1}{rd}}
#' \deqn{nnt\_se = \frac{rd\_se}{rd^2}}
#'
#' **When \code{baseline_risk} is available, the following conversions are performed.**
#'
#' **To estimate the odds ratio:**
#' Let \eqn{pt = baseline\_risk - rd} (treatment group risk). Then:
#' \deqn{or = \frac{pt \times (1 - baseline\_risk)}{baseline\_risk \times (1 - pt)}}
#' \deqn{logor\_se = \frac{rd\_se}{pt \times (1 - pt)}}
#' where the SE is derived via the delta method from \eqn{\frac{d(\log OR)}{d(RD)} = \frac{-1}{pt \times (1 - pt)}}.
#'
#' **To estimate the risk ratio:**
#' \deqn{rr = 1 - \frac{rd}{baseline\_risk}}
#' \deqn{logrr\_se = \frac{rd\_se}{|baseline\_risk - rd|}}
#' where the SE is derived via the delta method from \eqn{\frac{d(\log RR)}{d(RD)} = \frac{-1}{baseline\_risk - rd}}.
#'
#' **To estimate the Cohen's d and Hedges' g:**
#' The OR is first converted to Cohen's d using the formulas of Cooper et al. (2019):
#' \deqn{d = \log(or) \times \frac{\sqrt{3}}{\pi}}
#' \deqn{d\_se = \sqrt{\frac{logor\_se^2 \times 3}{\pi^2}}}
#' Then, the standard conversion to Hedges' g and correlation coefficients is applied.
#'
#' Note that the conversions to OR and RR assume the baseline risk is a fixed constant.
#'
#' @references
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Deeks, J.J. (2002). Issues in the selection of a summary statistic for meta-analysis of clinical trials with binary outcomes. Statistics in Medicine, 21(11), 1575-1600.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab RD\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + RR + NNT\cr
#'  \code{} \tab D + G + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab rd + rd_se\cr
#'  \tab \cr
#' }
#'
#' @export es_from_rd_se
#'
#' @md
#'
#' @examples
#' es_from_rd_se(rd = 0.15, rd_se = 0.05,
#'               baseline_risk = 0.30, n_exp = 100, n_nexp = 100)
es_from_rd_se <- function(rd, rd_se,
                          baseline_risk,
                          n_exp, n_nexp, n_cases, n_controls, n_sample,
                          reverse_rd) {

  if (missing(rd)) rd <- rep(NA_real_, length(rd_se))
  if (missing(rd_se)) rd_se <- rep(NA_real_, length(rd))
  if (missing(baseline_risk)) baseline_risk <- rep(NA_real_, length(rd))
  if (missing(n_exp)) n_exp <- rep(NA_real_, length(rd))
  if (missing(n_nexp)) n_nexp <- rep(NA_real_, length(rd))
  if (missing(n_cases)) n_cases <- rep(NA_real_, length(rd))
  if (missing(n_controls)) n_controls <- rep(NA_real_, length(rd))
  if (missing(n_sample)) n_sample <- rep(NA_real_, length(rd))
  if (missing(reverse_rd)) reverse_rd <- rep(FALSE, length(rd))
  reverse_rd[is.na(reverse_rd)] <- FALSE
  if (length(reverse_rd) == 1) reverse_rd <- rep(reverse_rd, length(rd))
  if (length(reverse_rd) != length(rd)) {
    stop("The length of the 'reverse_rd' argument is incorrectly specified.")
  }

  n_sample <- ifelse(is.na(n_sample),
    ifelse(!is.na(n_cases) & !is.na(n_controls),
           n_cases + n_controls,
           n_exp + n_nexp),
    n_sample
  )

  # OR -------
  treatment_risk <- baseline_risk - rd
  valid_pt <- !is.na(treatment_risk) & treatment_risk > 0 & treatment_risk < 1 &
              !is.na(baseline_risk) & baseline_risk > 0 & baseline_risk < 1

  or <- ifelse(valid_pt,
    (treatment_risk * (1 - baseline_risk)) / (baseline_risk * (1 - treatment_risk)),
    NA_real_
  )
  logOR <- suppressWarnings(log(or))
  # delta method
  logor_se <- ifelse(valid_pt, rd_se / abs(treatment_risk * (1 - treatment_risk)), NA_real_)

  d <- logOR * sqrt(3) / pi
  d_se <- sqrt(logor_se^2 * 3 / (pi^2))

  es <- .es_from_d(
    d = d, d_se = d_se, n_exp = n_exp, n_nexp = n_nexp,
    n_sample = n_sample, reverse = reverse_rd,
    smd_to_cor = rep("lipsey_cooper", length(d))
  )

  row_miss <- which(is.na(d_se))
  es$d[row_miss] <- es$d_se[row_miss] <-
    es$d_ci_lo[row_miss] <- es$d_ci_up[row_miss] <-
    es$g[row_miss] <- es$g_se[row_miss] <-
    es$g_ci_lo[row_miss] <- es$g_ci_up[row_miss] <- NA
  es$r[row_miss] <- es$r_se[row_miss] <-
    es$r_ci_lo[row_miss] <- es$r_ci_up[row_miss] <-
    es$z[row_miss] <- es$z_se[row_miss] <-
    es$z_ci_lo[row_miss] <- es$z_ci_up[row_miss] <- NA

  es$logor <- ifelse(reverse_rd, -logOR, logOR)
  es$logor_se <- logor_se
  es$logor_ci_lo <- es$logor - qnorm(.975) * logor_se
  es$logor_ci_up <- es$logor + qnorm(.975) * logor_se

  # RR -------
  valid_pt_rr <- !is.na(treatment_risk) & treatment_risk >= 0 & treatment_risk <= 1 &
                 !is.na(baseline_risk) & baseline_risk > 0 & baseline_risk <= 1
  rr <- ifelse(valid_pt_rr,
    treatment_risk / baseline_risk,
    NA_real_
  )
  rr <- ifelse(!is.na(rr) & rr > 0, rr, NA_real_)
  logRR <- suppressWarnings(log(rr))
  logrr_se <- ifelse(!is.na(rr) & !is.na(treatment_risk) & treatment_risk > 0,
    rd_se / abs(treatment_risk),
    NA_real_
  )

  es$logrr <- ifelse(reverse_rd, -logRR, logRR)
  es$logrr_se <- logrr_se
  es$logrr_ci_lo <- es$logrr - qnorm(.975) * logrr_se
  es$logrr_ci_up <- es$logrr + qnorm(.975) * logrr_se

  # RD -------
  es$rd <- ifelse(reverse_rd, -rd, rd)
  es$rd_se <- rd_se
  es$rd_ci_lo <- es$rd - qnorm(.975) * rd_se
  es$rd_ci_up <- es$rd + qnorm(.975) * rd_se

  # NNT -------
  es$nnt <- ifelse(rd == 0, NA, 1 / rd)
  es$nnt <- ifelse(reverse_rd, -es$nnt, es$nnt)
  es$nnt_se <- ifelse(rd == 0, NA, rd_se / rd^2)
  rd_ci_lo_raw <- rd - qnorm(.975) * rd_se
  rd_ci_up_raw <- rd + qnorm(.975) * rd_se
  crosses_zero <- (rd_ci_lo_raw < 0 & rd_ci_up_raw > 0) | rd == 0
  es$nnt_ci_lo <- ifelse(crosses_zero, NA,
                          ifelse(reverse_rd, -1 / rd_ci_lo_raw, 1 / rd_ci_up_raw))
  es$nnt_ci_up <- ifelse(crosses_zero, NA,
                          ifelse(reverse_rd, -1 / rd_ci_up_raw, 1 / rd_ci_lo_raw))

  es$info_used <- "rd_se"
  return(es)
}


#' Convert a risk difference value and its 95% confidence interval to several effect size measures
#'
#' @param rd risk difference value
#' @param rd_ci_lo lower bound of the 95% CI around the risk difference
#' @param rd_ci_up upper bound of the 95% CI around the risk difference
#' @param n_exp number of participants in the exposed/treatment group
#' @param n_nexp number of participants in the non-exposed/control group
#' @param n_cases number of cases/events across both groups
#' @param n_controls number of controls/no-event across both groups
#' @param n_sample total number of participants in the sample
#' @param baseline_risk proportion of cases in the non-exposed/control group
#' @param reverse_rd a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function computes the standard error of the risk difference from its 95% CI
#' (Section 6.5.2.2 in the Cochrane Handbook):
#' \deqn{rd\_se = \frac{rd\_ci\_up - rd\_ci\_lo}{2 \times qnorm(.975)}}
#'
#' Then, calculations of \code{\link{es_from_rd_se}()} are applied.
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab RD\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + RR + NNT\cr
#'  \code{} \tab D + G + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab rd + rd_ci_lo + rd_ci_up\cr
#'  \tab \cr
#' }
#'
#' @export es_from_rd_ci
#'
#' @md
#'
#' @examples
#' es_from_rd_ci(
#'   rd = 0.15, rd_ci_lo = 0.05, rd_ci_up = 0.25,
#'   baseline_risk = 0.30, n_exp = 100, n_nexp = 100
#' )
es_from_rd_ci <- function(rd, rd_ci_lo, rd_ci_up,
                          baseline_risk,
                          n_exp, n_nexp, n_cases, n_controls, n_sample,
                          reverse_rd) {

  if (missing(rd)) rd <- rep(NA_real_, length(rd_ci_lo))
  if (missing(rd_ci_lo)) rd_ci_lo <- rep(NA_real_, length(rd))
  if (missing(rd_ci_up)) rd_ci_up <- rep(NA_real_, length(rd))
  if (missing(baseline_risk)) baseline_risk <- rep(NA_real_, length(rd))
  if (missing(n_exp)) n_exp <- rep(NA_real_, length(rd))
  if (missing(n_nexp)) n_nexp <- rep(NA_real_, length(rd))
  if (missing(n_cases)) n_cases <- rep(NA_real_, length(rd))
  if (missing(n_controls)) n_controls <- rep(NA_real_, length(rd))
  if (missing(n_sample)) n_sample <- rep(NA_real_, length(rd))
  if (missing(reverse_rd)) reverse_rd <- rep(FALSE, length(rd))
  reverse_rd[is.na(reverse_rd)] <- FALSE

  rd_se <- (rd_ci_up - rd_ci_lo) / (2 * qnorm(.975))

  es <- es_from_rd_se(
    rd = rd, rd_se = rd_se,
    baseline_risk = baseline_risk,
    n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls,
    n_sample = n_sample,
    reverse_rd = reverse_rd
  )

  es$info_used <- "rd_ci"
  return(es)
}


#' Convert a risk difference value and its p-value to several effect size measures
#'
#' @param rd risk difference value
#' @param rd_pval p-value of the risk difference
#' @param n_exp number of participants in the exposed/treatment group
#' @param n_nexp number of participants in the non-exposed/control group
#' @param n_cases number of cases/events across both groups
#' @param n_controls number of controls/no-event across both groups
#' @param n_sample total number of participants in the sample
#' @param baseline_risk proportion of cases in the non-exposed/control group
#' @param reverse_rd_pval a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function uses the p-value of the risk difference to obtain the standard error
#' (Section 6.3.2 in the Cochrane Handbook):
#' \deqn{z = qnorm(rd\_pval / 2, lower.tail = FALSE)}
#' \deqn{rd\_se = |\frac{rd}{z}|}
#'
#' Then, calculations of \code{\link{es_from_rd_se}()} are applied.
#'
#' @references
#' Higgins, J. P., Thomas, J., Chandler, J., Cumpston, M., Li, T., Page, M. J., & Welch, V. A. (Eds.). (2019). Cochrane handbook for systematic reviews of interventions. John Wiley & Sons.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab RD\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + RR + NNT\cr
#'  \code{} \tab D + G + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab rd + rd_pval\cr
#'  \tab \cr
#' }
#'
#' @export es_from_rd_pval
#'
#' @md
#'
#' @examples
#' es_from_rd_pval(
#'   rd = 0.15, rd_pval = 0.01,
#'   baseline_risk = 0.30, n_exp = 100, n_nexp = 100
#' )
es_from_rd_pval <- function(rd, rd_pval,
                            baseline_risk,
                            n_exp, n_nexp, n_cases, n_controls, n_sample,
                            reverse_rd_pval) {

  if (missing(rd)) rd <- rep(NA_real_, length(rd_pval))
  if (missing(rd_pval)) rd_pval <- rep(NA_real_, length(rd))
  if (missing(baseline_risk)) baseline_risk <- rep(NA_real_, length(rd))
  if (missing(n_exp)) n_exp <- rep(NA_real_, length(rd))
  if (missing(n_nexp)) n_nexp <- rep(NA_real_, length(rd))
  if (missing(n_cases)) n_cases <- rep(NA_real_, length(rd))
  if (missing(n_controls)) n_controls <- rep(NA_real_, length(rd))
  if (missing(n_sample)) n_sample <- rep(NA_real_, length(rd))
  if (missing(reverse_rd_pval)) reverse_rd_pval <- rep(FALSE, length(rd))
  reverse_rd_pval[is.na(reverse_rd_pval)] <- FALSE

  z_rd <- qnorm(rd_pval / 2, lower.tail = FALSE)
  rd_se <- abs(rd / z_rd)

  es <- es_from_rd_se(
    rd = rd, rd_se = rd_se,
    baseline_risk = baseline_risk,
    n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls,
    n_sample = n_sample,
    reverse_rd = reverse_rd_pval
  )

  es$info_used <- "rd_pval"
  return(es)
}
