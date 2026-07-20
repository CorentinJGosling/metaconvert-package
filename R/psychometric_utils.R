#' Compute change-score reliability
#'
#' @param reliability single-occasion reliability coefficient (e.g., Cronbach's alpha or ICC)
#' @param r_pre_post pre-post correlation of observed scores
#'
#' @details
#' Computes the reliability of a change (difference) score from a single-occasion
#' reliability coefficient and the pre-post correlation of observed scores,
#' assuming equal variances AND equal reliabilities at pre and post
#' (Lord, 1963; Cronbach & Furby, 1970):
#'
#' \deqn{rel_{change} = \frac{rel - r_{12}}{1 - r_{12}}}
#'
#' This is the equal-variance/equal-reliability special case of the general
#' Lord (1963) formula; when the pre and post variances (or reliabilities)
#' differ materially, the general formula should be used instead.
#'
#' This is useful for disattenuating change-score correlations, where the
#' reliability of the change score is needed rather than the single-occasion reliability.
#'
#' A negative result (which occurs when \code{r_pre_post} exceeds
#' \code{reliability}) is population-impossible under classical test theory
#' and signals inconsistent inputs; it is returned as computed but triggers a
#' warning, since a negative value fed to \code{\link{es_disattenuate}} as a
#' reliability would produce \code{NaN} (square root of a negative number).
#'
#' @export reliability_change_score
#'
#' @references
#' Lord, F. M. (1963). Elementary models for measuring change.
#' In C. W. Harris (Ed.), Problems in measuring change. University of Wisconsin Press.
#'
#' Cronbach, L. J., & Furby, L. (1970). How we should measure "change" - or should we?
#' Psychological Bulletin, 74(1), 68-80.
#'
#' @md
#'
#' @return
#' A data.frame containing the change-score reliability (\code{rel_change}).
#'
#' @examples
#' reliability_change_score(reliability = 0.85, r_pre_post = 0.60)
reliability_change_score <- function(reliability, r_pre_post) {

  if (any(!is.na(r_pre_post) & r_pre_post >= 1, na.rm = TRUE)) {
    warning("r_pre_post >= 1 will produce Inf/NaN (division by zero in Lord formula)")
  }
  if (any(!is.na(reliability) & (reliability < 0 | reliability > 1), na.rm = TRUE)) {
    warning("reliability outside [0, 1] detected; results may be invalid")
  }
  if (any(!is.na(r_pre_post) & (r_pre_post < -1 | r_pre_post > 1), na.rm = TRUE)) {
    warning("r_pre_post outside [-1, 1] detected; results may be invalid")
  }

  rel_change <- (reliability - r_pre_post) / (1 - r_pre_post)

  # r_pre_post > reliability is population-impossible under classical test
  # theory (it would require a true-score correlation > 1); the negative
  # output is preserved but flagged, since it will produce NaN downstream if
  # used as a reliability in es_disattenuate().
  if (any(!is.na(rel_change) & rel_change < 0, na.rm = TRUE)) {
    warning("Negative change-score reliability produced (r_pre_post > reliability): inputs are inconsistent under classical test theory, and a negative value used as a reliability in es_disattenuate() will produce NaN")
  }

  result <- data.frame(
    rel_change = rel_change
  )

  return(result)
}


#' Compute Standard Error of Measurement (SEM) from reliability and SD
#'
#' @param sd standard deviation of the PROM/test scores
#' @param icc intraclass correlation coefficient (test-retest reliability)
#' @param n_sample sample size (used for the SEM sampling-variance estimation).
#'   Required in BOTH variance regimes: it feeds \eqn{Var(SD)} in the
#'   external-\code{icc_se} delta method and the degrees of freedom of the
#'   same-sample chi-square variance, so when it is omitted the SEM standard
#'   error and CI are \code{NA}.
#' @param icc_se standard error of the ICC (optional), on the RAW ICC scale.
#'   Note that \code{\link{es_from_icc}} under its default
#'   \code{icc_to_es = "bonett"} returns a column named \code{icc_se} on the
#'   \eqn{\ln(1 - ICC)} scale: that value must be converted first
#'   (\code{raw_se = bonett_se * (1 - icc)}) before being passed here (or use
#'   \code{icc_to_es = "raw"}, which returns the raw-scale SE directly).
#'   When supplied, the ICC is treated as an independent external estimate;
#'   when omitted, the ICC and SD are assumed to come from the same sample
#'   (see Details).
#' @param n_measurements number of measurement occasions or raters (k) used to
#'   estimate the ICC; default 2 (test-retest); must be >= 2. Only used when
#'   \code{icc_se} is not supplied.
#'
#' @details
#' Computes the standard error of measurement (SEM) from a reliability
#' coefficient and the standard deviation of scores:
#'
#' \deqn{SEM = SD \times \sqrt{1 - ICC}}
#'
#' The sampling variance of SEM is computed in one of two ways:
#'
#' 1. **When \code{icc_se} is supplied**, the ICC is treated as an external
#' estimate, independent of the within-sample SD, and the bivariate delta method
#' (with zero covariance) is used:
#' \deqn{Var(SEM) = \frac{SD^2}{4(1 - ICC)} \times Var(ICC) + (1 - ICC) \times Var(SD)}
#' with \eqn{Var(ICC) = icc\_se^2} and \eqn{Var(SD) \approx SD^2 / (2(n-1))}.
#'
#' 2. **When \code{icc_se} is not supplied**, the ICC and SD are assumed to come
#' from the same reliability sample, where they are strongly positively
#' correlated. Because \eqn{SEM = SD\sqrt{1 - ICC}} equals the within-subject
#' residual root-mean-square \eqn{\sqrt{MSE}}, and \eqn{MSE / \sigma_e^2} follows
#' a scaled chi-square with \eqn{(n-1)(k-1)} degrees of freedom, the exact
#' sampling variance is
#' \deqn{Var(SEM) = \frac{SEM^2}{2(n-1)(k-1)}}
#' where \eqn{k} is the number of measurement occasions/raters
#' (\code{n_measurements}, default 2 for test-retest). Treating the same-sample SD
#' and ICC as independent (the delta method of case 1) would overestimate this
#' variance by roughly 2-6x, so the exact form is used instead.
#'
#' The 95% CI is likewise regime-specific. In the same-sample regime (case 2)
#' the degrees of freedom \eqn{df = (n - 1)(k - 1)} are known, so the exact
#' chi-square interval
#' \deqn{[SEM \sqrt{df / \chi^2_{0.975, df}},\; SEM \sqrt{df / \chi^2_{0.025, df}}]}
#' is returned; the symmetric Wald interval undercovers at small df (89.7%
#' instead of 95% at n = 10, k = 2). When \code{icc_se} is supplied (case 1)
#' the df behind the external estimate is unknown, so the Wald interval
#' \eqn{SEM \pm 1.96 \times SE}, truncated at 0, is kept.
#'
#' Degenerate inputs: \code{n_measurements} = 1 makes the same-sample df zero,
#' so the SE and CI are returned as \code{NA} with a warning; an ICC > 1 is
#' impossible and yields \code{NA} for the SEM, its SE and CI (with a warning).
#'
#' @export compute_sem
#'
#' @references
#' Weir, J. P. (2005). Quantifying test-retest reliability using the intraclass correlation
#' coefficient and the SEM. Journal of Strength and Conditioning Research, 19(1), 231-240.
#'
#' @md
#'
#' @return
#' A data.frame with the SEM, its standard error and 95% CI.
#'
#' @examples
#' compute_sem(sd = 10, icc = 0.85, n_sample = 100)
compute_sem <- function(sd, icc, n_sample, icc_se, n_measurements = 2) {

  if (missing(icc_se)) icc_se <- rep(NA_real_, length(sd))
  if (missing(n_sample)) n_sample <- rep(NA_real_, length(sd))
  if (missing(n_measurements)) n_measurements <- rep(2, length(sd))
  n_measurements[is.na(n_measurements)] <- 2
  if (length(n_measurements) == 1) n_measurements <- rep(n_measurements, length(sd))

  if (any(!is.na(sd) & sd < 0, na.rm = TRUE)) {
    warning("Negative SD detected; SEM will be invalid")
  }
  if (any(!is.na(icc) & (icc < 0 | icc > 1), na.rm = TRUE)) {
    warning("ICC outside [0, 1] detected; SEM is invalid (set to NA when ICC > 1)")
  }
  # icc_se must be the RAW-scale SE of the ICC. A raw-scale SE larger than
  # (1 - icc) is implausible and is the signature of es_from_icc()'s default
  # Bonett output (the SE of ln(1 - ICC)) being passed through unchanged.
  # Warn only - the value is used exactly as supplied.
  scale_suspect <- !is.na(icc_se) & !is.na(icc) & icc < 1 & icc_se > (1 - icc)
  if (any(scale_suspect, na.rm = TRUE)) {
    warning("icc_se looks like a transformed-scale (Bonett, ln(1 - ICC)) standard error (icc_se > 1 - icc); compute_sem() needs the RAW-scale ICC standard error. If it comes from es_from_icc() with the default icc_to_es = 'bonett', convert it first: raw_se = bonett_se * (1 - icc). The value is used as supplied.")
  }

  # ICC > 1 is impossible: sqrt(1 - icc) would be NaN, and no variance is
  # defensible - return NA (not NaN paired with a zero SE).
  sem <- sd * sqrt(pmax(1 - icc, 0))
  sem[!is.na(icc) & icc > 1] <- NA_real_
  k <- n_measurements

  # Same-sample df = (n - 1)(k - 1); k = 1 (or n = 1) makes it zero, so the
  # chi-square variance/CI below would divide by zero - invalidate the df and
  # warn where the same-sample regime actually needs it.
  df_sem <- (n_sample - 1) * (k - 1)
  df_sem[!is.na(df_sem) & df_sem < 1] <- NA_real_
  if (any(is.na(icc_se) & !is.na(k) & k < 2, na.rm = TRUE)) {
    warning("n_measurements must be >= 2 to estimate the same-sample SEM sampling variance; SE and CI set to NA for the affected row(s)")
  }

  # Sampling variance of SEM, in two regimes:
  #
  # (A) icc_se SUPPLIED -- the ICC is an external estimate, independent of the
  #     within-sample SD, so the bivariate delta method with Cov(SD, ICC) = 0
  #     applies:
  #       Var(SEM) = (SD^2 / (4(1 - ICC))) * Var(ICC) + (1 - ICC) * Var(SD)
  #     with Var(ICC) = icc_se^2 and Var(SD) = SD^2 / (2(n - 1)).
  #
  # (B) icc_se NOT supplied -- ICC and SD come from the SAME reliability sample and
  #     are strongly positively correlated; treating them as independent (case A)
  #     overestimates Var(SEM) by ~2-6x. Since SEM = SD * sqrt(1 - ICC) = sqrt(MSE)
  #     (the within-subject residual RMS), and MSE / sigma_e^2 ~ chi-square / df
  #     with df = (n - 1)(k - 1), the exact sampling variance is
  #       Var(SEM) = SEM^2 / (2 (n - 1)(k - 1)),   k = n_measurements (default 2).
  #     Verified by simulation to within ~1%.
  var_sd <- sd^2 / (2 * (n_sample - 1))
  var_sem_delta <- (sd^2 / (4 * (1 - icc))) * icc_se^2 + (1 - icc) * var_sd
  var_sem_exact <- sem^2 / (2 * df_sem)

  # icc > 1: no defensible variance (NA, not the old maximally-confident 0);
  # icc == 1 exactly: the coherent degenerate case sem = 0, se = 0.
  var_sem <- ifelse(!is.na(icc) & icc > 1, NA_real_,
                    ifelse(!is.na(icc) & icc == 1, 0,
                           ifelse(!is.na(icc_se), var_sem_delta, var_sem_exact)))
  sem_se <- sqrt(var_sem)

  # CI, regime-specific:
  # - same-sample (icc_se missing): df is known, use the exact chi-square
  #   interval for sigma_e (the Wald interval undercovers at small df,
  #   89.7% at n = 10, k = 2);
  # - external icc_se: df unknown, keep the truncated Wald interval.
  same_sample <- is.na(icc_se)
  sem_ci_lo <- ifelse(same_sample,
                      sem * sqrt(df_sem / stats::qchisq(0.975, df_sem)),
                      pmax(0, sem - qnorm(0.975) * sem_se))
  sem_ci_up <- ifelse(same_sample,
                      sem * sqrt(df_sem / stats::qchisq(0.025, df_sem)),
                      sem + qnorm(0.975) * sem_se)

  result <- data.frame(
    sem = sem,
    sem_se = sem_se,
    sem_ci_lo = sem_ci_lo,
    sem_ci_up = sem_ci_up
  )

  return(result)
}


#' Compute Smallest Detectable Change (SDC) from SEM
#'
#' @param sem standard error of measurement
#' @param sem_se standard error of the SEM (optional; for propagating uncertainty)
#'
#' @details
#' Computes the smallest detectable change (SDC) from the standard error
#' of measurement:
#'
#' \deqn{SDC = 1.96 \times \sqrt{2} \times SEM}
#'
#' The SDC represents the smallest change in scores that can be detected
#' beyond measurement error with 95% confidence. If \code{sem_se} is provided,
#' the uncertainty in SDC is propagated:
#'
#' \deqn{SDC\_se = 1.96 \times \sqrt{2} \times SEM\_se}
#'
#' @export compute_sdc
#'
#' @references
#' de Vet, H. C. W., Terwee, C. B., Mokkink, L. B., & Knol, D. L. (2011).
#' Measurement in Medicine: A Practical Guide. Cambridge University Press.
#'
#' @md
#'
#' @return
#' A data.frame with the SDC, its standard error and 95% CI.
#'
#' @examples
#' sem_res <- compute_sem(sd = 10, icc = 0.85, n_sample = 100)
#' compute_sdc(sem = sem_res$sem, sem_se = sem_res$sem_se)
compute_sdc <- function(sem, sem_se) {

  if (missing(sem_se)) sem_se <- rep(NA, length(sem))

  multiplier <- qnorm(0.975) * sqrt(2)

  sdc <- multiplier * sem
  sdc_se <- multiplier * sem_se

  sdc_ci_lo <- pmax(0, sdc - qnorm(0.975) * sdc_se)
  sdc_ci_up <- sdc + qnorm(0.975) * sdc_se

  result <- data.frame(
    sdc = sdc,
    sdc_se = sdc_se,
    sdc_ci_lo = sdc_ci_lo,
    sdc_ci_up = sdc_ci_up
  )

  return(result)
}
