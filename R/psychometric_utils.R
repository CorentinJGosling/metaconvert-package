#' Compute change-score reliability
#'
#' @param reliability single-occasion reliability coefficient (e.g., Cronbach's alpha or ICC)
#' @param r_pre_post pre-post correlation of observed scores
#'
#' @details
#' Computes the reliability of a change (difference) score from a single-occasion
#' reliability coefficient and the pre-post correlation of observed scores,
#' assuming equal variances at pre and post (Lord, 1963; Cronbach & Furby, 1970):
#'
#' \deqn{rel_{change} = \frac{rel - r_{12}}{1 - r_{12}}}
#'
#' This is useful for disattenuating change-score correlations, where the
#' reliability of the change score is needed rather than the single-occasion reliability.
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

  result <- data.frame(
    rel_change = rel_change
  )

  return(result)
}


#' Compute Standard Error of Measurement (SEM) from reliability and SD
#'
#' @param sd standard deviation of the PROM/test scores
#' @param icc intraclass correlation coefficient (test-retest reliability)
#' @param n_sample sample size (used for delta-method variance estimation)
#' @param icc_se standard error of the ICC (optional; if not provided, a rough approximation is used)
#'
#' @details
#' Computes the standard error of measurement (SEM) from a reliability
#' coefficient and the standard deviation of scores:
#'
#' \deqn{SEM = SD \times \sqrt{1 - ICC}}
#'
#' The sampling variance of SEM is derived via the delta method:
#' \deqn{Var(SEM) = \frac{SD^2}{4(1 - ICC)} \times Var(ICC) + (1 - ICC) \times Var(SD)}
#'
#' where \eqn{Var(SD) \approx SD^2 / (2(n-1))}.
#'
#' If \code{icc_se} is not provided, the variance of ICC is approximated as
#' \eqn{Var(ICC) \approx (1 - ICC^2)^2 / (n - 1)}, which is a rough approximation.
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
compute_sem <- function(sd, icc, n_sample, icc_se) {

  if (missing(icc_se)) icc_se <- rep(NA, length(sd))
  if (missing(n_sample)) n_sample <- rep(NA_real_, length(sd))

  if (any(!is.na(sd) & sd < 0, na.rm = TRUE)) {
    warning("Negative SD detected; SEM will be invalid")
  }
  if (any(!is.na(icc) & (icc < 0 | icc > 1), na.rm = TRUE)) {
    warning("ICC outside [0, 1] detected; SEM will be invalid (NaN or negative)")
  }

  sem <- sd * sqrt(1 - icc)

  # d(SEM)/d(ICC) = -SD / (2*sqrt(1-ICC))
  # d(SEM)/d(SD) = sqrt(1-ICC) = SEM/SD
  var_sd <- sd^2 / (2 * (n_sample - 1))
  var_icc <- ifelse(!is.na(icc_se), icc_se^2,
                    (1 - icc^2)^2 / (n_sample - 1))  # rough approximation

  var_sem <- ifelse(icc == 1, 0,
                    (sd^2 / (4 * (1 - icc))) * var_icc + (1 - icc) * var_sd)
  sem_se <- sqrt(var_sem)

  sem_ci_lo <- pmax(0, sem - qnorm(0.975) * sem_se)
  sem_ci_up <- sem + qnorm(0.975) * sem_se

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
