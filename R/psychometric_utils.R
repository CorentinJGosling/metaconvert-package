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

#' Back-transform a pooled reliability estimate to the coefficient scale
#'
#' @param x a numeric vector of estimates on the analysis scale, OR an
#'   \code{rma}/\code{rma.mv} object returned by \code{metafor}, OR the output of
#'   \code{metafor::predict.rma()}.
#' @param ci_lo lower confidence bound(s) on the analysis scale. Ignored when
#'   \code{x} is an \code{rma} or \code{predict.rma} object.
#' @param ci_up upper confidence bound(s) on the analysis scale. Ignored when
#'   \code{x} is an \code{rma} or \code{predict.rma} object.
#' @param method the transformation that produced \code{x}, i.e. the value passed
#'   as \code{alpha_to_es} / \code{icc_to_es} in \code{\link{convert_df}}. Must be
#'   either \code{"bonett"} (default) or \code{"raw"}.
#'
#' @details
#' Under the default \code{"bonett"} scale, \code{\link{convert_df}} returns
#' \eqn{T = \ln(1 - \rho)}, so the coefficient is recovered with
#'
#' \deqn{\rho = 1 - \exp(T)}
#'
#' \strong{This map is monotone decreasing, so the confidence bounds swap}: the
#' upper bound on the \eqn{T} scale is the \emph{lower} bound on the reliability
#' scale. Getting that swap wrong is the single most-botched step in published
#' reliability-generalization work, and it is silent -- the interval still looks
#' like an interval. This helper performs the swap for you.
#'
#' It also exists to keep users away from \code{metafor::transf.iabt}, which is
#' \strong{not} a valid back-transformation for \code{metaConvert} output.
#' \code{transf.iabt} implements \eqn{1 - \exp(-x)} and then \emph{clamps
#' negative inputs to zero}, because \code{metafor}'s \code{measure = "ABT"}
#' stores \eqn{-\ln(1 - \rho)} where \code{metaConvert} follows Bonett (2002) and
#' stores \eqn{+\ln(1 - \rho)}. The standard errors are identical, but every sign
#' is opposite, so \code{metaConvert}'s values are always negative and
#' \code{transf.iabt} maps the whole pool to exactly 0 with no error, no warning
#' and no \code{NA}.
#'
#' A negative back-transformed reliability is returned as computed rather than
#' clamped: a negative Cronbach's alpha is unusual but mathematically possible
#' (negative average inter-item covariance), and silently flooring it at zero
#' would hide the very extraction problem worth seeing.
#'
#' Under \code{method = "raw"} the analysis scale already is the coefficient
#' scale, so values pass through unchanged and the bounds are not swapped. Both
#' methods are accepted so that the same reporting code works whichever
#' \code{alpha_to_es} / \code{icc_to_es} was used.
#'
#' @export reliability_backtransform
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for testing and estimating
#' coefficient alpha. Journal of Educational and Behavioral Statistics, 27(4),
#' 335-340.
#'
#' @md
#'
#' @return
#' A data.frame with the back-transformed reliability and its interval
#' (\code{reliability}, \code{reliability_ci_lo}, \code{reliability_ci_up}).
#' When \code{x} is an \code{rma} object or \code{predict.rma} output that
#' carries a prediction interval, \code{reliability_pi_lo} and
#' \code{reliability_pi_up} are added.
#'
#' @examples
#' # A per-study summary() column on the Bonett scale
#' reliability_backtransform(c(-2.12, -1.90, -1.61))
#'
#' # A pooled metafor model: CI and prediction interval, bounds swapped for you
#' \donttest{
#' es <- es_from_cronbach_alpha(
#'   cronbach_alpha = c(0.78, 0.85, 0.91, 0.88),
#'   n_sample = 200, n_items = 10
#' )
#' m <- metafor::rma(yi = es$alpha, sei = es$alpha_se, method = "REML")
#' reliability_backtransform(m)
#' }
reliability_backtransform <- function(x, ci_lo, ci_up, method = "bonett") {

  if (!method %in% c("bonett", "raw")) {
    stop(paste0("'", method, "' not in tolerated values for the 'method' argument. ",
                "Possible inputs are: 'bonett', 'raw'"))
  }

  pi_lo <- NULL
  pi_up <- NULL

  # metafor objects carry their own bounds; taking them from the object rather
  # than from the user is the point of accepting them here, since hand-copying
  # ci.lb/ci.ub (or pi.lb/pi.ub) into the wrong slot is exactly the error this
  # function exists to prevent.
  if (inherits(x, "rma")) {
    fit <- x
    # A moderator model has no single pooled reliability: its coefficients are
    # contrasts on the transformed scale, and 1 - exp() of a contrast is not a
    # reliability of anything. Refuse explicitly and name the route that works,
    # rather than returning the reference level (fixed-effect fits, which carry
    # no prediction interval, would otherwise do exactly that, silently).
    n_coef <- if (!is.null(fit$X)) ncol(fit$X) else length(fit$beta)
    if (isTRUE(n_coef > 1L)) {
      stop(paste0(
        "'x' is a meta-regression with ", n_coef, " coefficients, which has no ",
        "single pooled reliability to back-transform. Pass the fitted values ",
        "instead -- reliability_backtransform(predict(x)) -- or, for named ",
        "moderator levels, predict(x, newmods = ...)."
      ))
    }
    x <- as.numeric(fit$beta[1])
    ci_lo <- as.numeric(fit$ci.lb[1])
    ci_up <- as.numeric(fit$ci.ub[1])
    pred <- try(stats::predict(fit), silent = TRUE)
    if (!inherits(pred, "try-error") && !is.null(pred$pi.lb)) {
      pi_lo <- as.numeric(pred$pi.lb)
      pi_up <- as.numeric(pred$pi.ub)
    }
  } else if (inherits(x, "list.rma") ||
             (is.list(x) && !is.null(x[["pred"]]))) {
    # A prediction table is a prediction table whether or not it has been through
    # as.data.frame() / tibble(); excluding data.frames here sent that (routine)
    # shape to as.numeric(), which FLATTENED the columns and returned one bogus
    # "reliability" per column with row 1 coincidentally correct. Read members by
    # [[ ]] so partial matching cannot pick up an unrelated `pred_type` column.
    pred <- x
    x <- as.numeric(pred[["pred"]])
    # blup() and some predict() shapes carry pi.* but no ci.*; a zero-length
    # numeric here would defeat the is.null() fallback below and abort with a
    # complaint about a `ci_lo` the caller never passed.
    ci_lo <- if (is.null(pred[["ci.lb"]])) rep(NA_real_, length(x)) else as.numeric(pred[["ci.lb"]])
    ci_up <- if (is.null(pred[["ci.ub"]])) rep(NA_real_, length(x)) else as.numeric(pred[["ci.ub"]])
    if (!is.null(pred[["pi.lb"]])) {
      pi_lo <- as.numeric(pred[["pi.lb"]])
      pi_up <- as.numeric(pred[["pi.ub"]])
    }
  }

  x <- as.numeric(x)
  if (missing(ci_lo) || is.null(ci_lo)) ci_lo <- rep(NA_real_, length(x))
  if (missing(ci_up) || is.null(ci_up)) ci_up <- rep(NA_real_, length(x))
  ci_lo <- as.numeric(ci_lo)
  ci_up <- as.numeric(ci_up)

  if (length(ci_lo) == 1) ci_lo <- rep(ci_lo, length(x))
  if (length(ci_up) == 1) ci_up <- rep(ci_up, length(x))
  if (length(ci_lo) != length(x)) stop("The length of the 'ci_lo' argument is incorrectly specified.")
  if (length(ci_up) != length(x)) stop("The length of the 'ci_up' argument is incorrectly specified.")

  # bonett: 1 - exp(t) is DECREASING, so the transformed-scale upper bound is the
  # reliability-scale LOWER bound. raw: identity, bounds keep their roles.
  if (method == "bonett") {
    bt <- function(t) 1 - exp(t)
    out_lo <- bt(ci_up)
    out_up <- bt(ci_lo)
  } else {
    bt <- function(t) t
    out_lo <- ci_lo
    out_up <- ci_up
  }

  result <- data.frame(
    reliability = bt(x),
    reliability_ci_lo = out_lo,
    reliability_ci_up = out_up
  )

  # defensive: never recycle a prediction interval of a different length onto the
  # estimate (that mismatch is what made the moderator case abort inside
  # `$<-.data.frame` with a message naming no argument of this function).
  if (!is.null(pi_lo) && length(pi_lo) != length(x)) { pi_lo <- NULL; pi_up <- NULL }

  if (!is.null(pi_lo)) {
    if (method == "bonett") {
      result$reliability_pi_lo <- bt(pi_up)
      result$reliability_pi_up <- bt(pi_lo)
    } else {
      result$reliability_pi_lo <- pi_lo
      result$reliability_pi_up <- pi_up
    }
  }

  return(result)
}
