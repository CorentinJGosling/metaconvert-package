#' Convert an intraclass correlation coefficient (ICC) into an effect size measure
#'
#' @param icc intraclass correlation coefficient value
#' @param n_sample total sample size (number of subjects)
#' @param n_measurements number of measurements or raters
#' @param icc_type ICC type: \code{"agreement"} for ICC(2,1) (two-way random, absolute agreement)
#'   or \code{"consistency"} for ICC(3,1) (two-way mixed, consistency). Default is \code{"agreement"}.
#'   Matching is case- and punctuation-insensitive, so \code{"ICC(2,1)"}, \code{"absolute agreement"}
#'   and \code{"two-way random"} are all accepted for the first, and \code{"ICC(3,1)"} /
#'   \code{"two-way mixed"} for the second.
#'
#'   \strong{Average-measures ICCs} - \code{"average"}, \code{"ICC(2,k)"}, \code{"ICC(3,k)"} - are
#'   also recognised, and are \emph{not} the same estimand: they describe the mean of
#'   \code{n_measurements} measurements rather than one. Such a value is stepped down to the
#'   single-measurement ICC with the inverse Spearman-Brown formula
#'   \deqn{ICC_1 = \frac{ICC_k}{k - (k - 1) ICC_k}}
#'   so that it is on the same scale as the rest of the pool, and the row is flagged
#'   \code{[INFO]} recording both values. If \code{n_measurements} is missing the step-down
#'   cannot be performed, and the row is set to \code{NA} with an \code{[INVALID]} flag rather
#'   than pooled as if it were single-measures - entering an ICC(2,k) as if it were ICC(2,1)
#'   overstates the reliability by up to 1.9 log units (k = 10, ICC = 0.95) while leaving the
#'   standard error almost unchanged, so nothing downstream reveals the error.
#' @param icc_se standard error of the ICC, \strong{on the natural (raw ICC) scale}.
#'   Optional. When supplied it takes precedence over the closed-form \code{(n, k)}
#'   standard error.
#' @param icc_ci_lo lower bound of the 95% confidence interval of the ICC (natural scale)
#' @param icc_ci_up upper bound of the 95% confidence interval of the ICC (natural scale)
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
#' \strong{Where the standard error comes from.} Three sources, in this order
#' (the same precedence \code{\link{es_from_omega}} uses):
#' \enumerate{
#'   \item \code{icc_se}, read on the natural scale and delta-mapped onto the
#'     analysis scale;
#'   \item otherwise \code{icc_ci_lo} / \code{icc_ci_up}, transformed at the
#'     \strong{bounds} so an asymmetric (e.g. bootstrap) interval maps correctly
#'     instead of being symmetrised first;
#'   \item otherwise the closed-form \eqn{(n, k)} standard error above.
#' }
#' The ordering matters because of the coverage problem documented below: for an
#' agreement-type ICC the computed standard error is a one-way approximation, so
#' wherever the study reported its own uncertainty that is the better source. The
#' ICC literature reports intervals routinely (\emph{"ICC 0.94, 95\% CI 0.86 to
#' 0.98"}), and such a row is now usable even when no sample size is given -
#' \code{n_sample} and \code{n_measurements} gate only source 3, not the point
#' estimate. A row with none of the three keeps its effect size and gets
#' \code{icc_se = NA}, so it stays visible and countable but is dropped from the
#' pool by \code{summary()}.
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
                        icc_se, icc_ci_lo, icc_ci_up,
                        icc_to_es = "bonett") {

  len <- length(icc)
  if (missing(n_sample)) n_sample <- rep(NA_real_, len)
  if (missing(n_measurements)) n_measurements <- rep(NA_real_, len)
  if (missing(icc_type)) icc_type <- rep("agreement", len)
  if (missing(icc_se)) icc_se <- rep(NA_real_, len)
  if (missing(icc_ci_lo)) icc_ci_lo <- rep(NA_real_, len)
  if (missing(icc_ci_up)) icc_ci_up <- rep(NA_real_, len)
  icc_type[is.na(icc_type)] <- "agreement"

  # See the note in es_from_cronbach_alpha(): `icc_type` was already recycled
  # here, but the numeric arguments were not, so a scalar `n_sample` or
  # `n_measurements` beside a vector of ICCs returned an NA standard error for
  # every row after the first while leaving the point estimates intact.
  rec <- function(v, nm) {
    if (length(v) == 1) return(rep(v, len))
    if (length(v) != len) {
      stop(paste0("The length of the '", nm, "' argument is incorrectly specified."))
    }
    v
  }
  n_sample       <- rec(n_sample, "n_sample")
  n_measurements <- rec(n_measurements, "n_measurements")
  icc_type       <- rec(icc_type, "icc_type")
  icc_se         <- rec(icc_se, "icc_se")
  icc_ci_lo      <- rec(icc_ci_lo, "icc_ci_lo")
  icc_ci_up      <- rec(icc_ci_up, "icc_ci_up")

  # A non-positive dispersion cannot be a standard error; a transposed interval is
  # the same interval (R/internal_guards.R), as for omega.
  icc_se <- .positive_or_na(icc_se)
  ci_lo  <- .ci_lower(icc_ci_lo, icc_ci_up)
  ci_up  <- .ci_upper(icc_ci_lo, icc_ci_up)

  icc_type <- .normalise_icc_type(icc_type)

  # Roadmap 1.2. An average-measures ICC describes the mean of k measurements, so
  # it is a DIFFERENT estimand from the single-measurement ICC every other row in
  # the pool carries. Step it down with Spearman-Brown where k is known; where it
  # is not, leave the value alone and let the n_measurements guard below set the
  # row to NA (a step-down cannot be guessed, and pooling the un-stepped value is
  # the silent error this exists to stop). .validate_input_data() raises [INFO] on
  # the stepped rows and [INVALID] on the ones dropped for want of k.
  avg <- .icc_is_average(icc_type)
  can_step <- avg & !is.na(icc) & !is.na(n_measurements) &
    is.finite(n_measurements) & n_measurements >= 2 & !is.na(icc) & abs(icc) <= 1
  if (any(can_step)) {
    k_s   <- n_measurements[can_step]
    rho_k <- icc[can_step]
    # The REPORTED uncertainty is on the average-measures scale too, so it has to
    # travel with the point estimate or the row ends up with an ICC(1) estimate
    # carrying an ICC(k) interval. Spearman-Brown is monotone increasing
    # (f'(x) = k / (k - (k-1)x)^2 > 0), so the CI bounds map straight through and
    # keep their order; a reported SE maps by that same derivative.
    icc_se[can_step] <- icc_se[can_step] * k_s / (k_s - (k_s - 1) * rho_k)^2
    ci_lo[can_step]  <- .icc_step_down(ci_lo[can_step], k_s)
    ci_up[can_step]  <- .icc_step_down(ci_up[can_step], k_s)
    icc[can_step]    <- .icc_step_down(rho_k, k_s)
  }
  # Report the estimand actually computed. A row that could not be stepped down is
  # NA'd below, so no un-stepped average value ever reaches the output.
  icc_type[avg] <- sub("_average$", "", icc_type[avg])
  if (any(avg & !can_step)) icc[avg & !can_step] <- NA_real_

  if (!icc_to_es %in% c("bonett", "raw")) {
    stop(paste0("'", icc_to_es, "' not in tolerated values for the 'icc_to_es' argument. ",
                "Possible inputs are: 'bonett', 'raw'"))
  }

  # Forward transform and the delta map that carries a RAW-scale reported SE onto it.
  fwd    <- function(r) if (identical(icc_to_es, "bonett")) log(1 - r) else r
  se_map <- function(r, s) if (identical(icc_to_es, "bonett")) s / (1 - r) else s

  # P16: per-element guards -- |icc| > 1 is impossible (V11 bounds) and icc = 1
  # has no Bonett transform (log(0)) and a degenerate raw SE. Direct calls degrade
  # to NA instead of emitting Inf/NaN. n_sample/n_measurements are NOT part of this
  # test any more: they gate only the COMPUTED standard error (route 3 below), not
  # the point estimate, so a study reporting "ICC 0.94 (95% CI 0.86-0.98)" with no
  # sample size is now usable instead of being dropped whole.
  valid_es <- !is.na(icc) & icc >= -1 & icc < 1

  n <- length(icc)
  icc_es    <- rep(NA_real_, n)
  icc_es_se <- rep(NA_real_, n)
  icc_es[valid_es] <- fwd(icc[valid_es])

  # --- standard error, in order of preference -------------------------------
  # Mirrors es_from_omega(). The ordering is deliberate: for an agreement-type
  # ICC the computed (n, k) SE is a one-way approximation whose coverage runs
  # 0.82 at n = 20 down to 0.14 at n = 1000 (see the @details), so wherever the
  # study reported its own uncertainty that is the better source. For a
  # consistency ICC the computed SE is exact at leading order and route 3 is
  # reached whenever nothing was reported.

  # 1. reported SE, interpreted on the RAW ICC scale and delta-mapped
  from_se <- which(valid_es & !is.na(icc_se))
  if (length(from_se)) icc_es_se[from_se] <- se_map(icc[from_se], icc_se[from_se])

  # 2. otherwise the reported CI, transformed at the BOUNDS so an asymmetric
  #    (e.g. bootstrap) interval maps correctly rather than being symmetrised
  #    first. A bound of exactly 1 has no log(1 - .) and is excluded on the
  #    Bonett scale only.
  ci_usable <- if (identical(icc_to_es, "raw")) {
    rep(TRUE, n)
  } else {
    ci_lo < 1 & ci_up < 1
  }
  ci_usable[is.na(ci_usable)] <- FALSE
  from_ci <- which(valid_es & is.na(icc_es_se) &
                     !is.na(ci_lo) & !is.na(ci_up) & ci_usable)
  if (length(from_ci)) {
    icc_es_se[from_ci] <-
      abs(fwd(ci_up[from_ci]) - fwd(ci_lo[from_ci])) / (2 * qnorm(0.975))
  }

  # 3. otherwise the closed-form (n, k) SE.
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
  from_nk <- which(valid_es & is.na(icc_es_se) &
                     !is.na(n_sample) & n_sample > 1 &
                     !is.na(n_measurements) & n_measurements >= 2)
  if (length(from_nk)) {
    rho <- icc[from_nk]
    ns  <- n_sample[from_nk]
    k   <- n_measurements[from_nk]
    bonett_transformed_se <- sqrt(
      2 * (1 + (k - 1) * rho)^2 / (k * (k - 1) * (ns - 1))
    )
    icc_es_se[from_nk] <- if (identical(icc_to_es, "bonett")) {
      bonett_transformed_se
    } else {
      (1 - rho) * bonett_transformed_se
    }
  }
  icc_es_se <- .positive_or_na(icc_es_se)

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


# Normalise icc_type to one of 'agreement' / 'consistency'.
#
# ONE definition, called from three places -- es_from_icc() (the route),
# convert_df() (which normalises the stored column before validation) and the V31
# check in internal_flags.R. Before this, es_from_icc() normalised while V31 used
# bare string equality on the raw cell, so the rows that GET the anti-conservative
# agreement SE were largely the rows that did not get warned about it: of ten
# spellings that resolve to 'agreement', only the two already spelled exactly
# 'agreement' fired V31. Two places independently deciding what 'agreement' means
# is the same rot recorded for the hardcoded route list in roadmap 1.2.
#
# The key keeps DIGITS ([^a-z0-9], not [^a-z]). The old inline version stripped
# them, which silently made the 'icc21' and 'icc31' map entries unreachable: every
# ICC(2,1) / ICC(3,1) / icc21 spelling collapsed to the key 'icc', missed the map,
# and fell through to the 'agreement' default -- so a CONSISTENCY ICC entered as
# 'ICC(3,1)' was relabelled as agreement, with a warning that named it unrecognised.
# The shared SE is unaffected (it is one formula for both types), so the cost was a
# wrong estimand label and a spurious V31 on consistency rows, not a wrong number.
#
# Idempotent, like .normalise_omega_type(): convert_df() normalises the column and
# es_from_icc() normalises again, so every output must also be a valid input.
#
# AVERAGE-MEASURES levels are a separate axis, not aliases. icc_type carries two
# independent facts -- the model (absolute agreement vs consistency) and the UNIT
# (a single measurement vs the mean of k). Before roadmap 1.2 the second was not
# representable at all: 'average', 'ICC(2,k)' and 'icc2k' were unrecognised, warned,
# and fell back to 'agreement', so an average-measures value was computed as if it
# were single-measures. The SE ratio stays near 1 (1.05-1.39x), so nothing
# downstream looks wrong, while the POINT ESTIMATE is off by up to 1.9 log units:
#
#   k    ICC_avg   true ICC_1   es (wrong)   es (right)   error
#   2    0.90      0.8182       -2.3026      -1.7047      -0.598
#   5    0.90      0.6429       -2.3026      -1.0296      -1.273
#   10   0.95      0.6552       -2.9957      -1.0647      -1.931
#
# The four *_average levels are resolved by es_from_icc(), which steps them down
# with Spearman-Brown and reports the single-measures level it actually computed.
.normalise_icc_type <- function(x, warn = TRUE) {
  raw <- as.character(x)
  key <- gsub('[^a-z0-9]', '', tolower(trimws(raw)))
  map <- c(
    agreement = 'agreement', absoluteagreement = 'agreement',
    icc21 = 'agreement', twowayrandom = 'agreement', absolute = 'agreement',
    agree = 'agreement', icc2 = 'agreement',
    consistency = 'consistency', icc31 = 'consistency',
    twowaymixed = 'consistency', consistent = 'consistency', icc3 = 'consistency',
    # average-measures. A bare 'average' / 'ICC(k)' does not say which model, so it
    # takes the package default (agreement), exactly as a bare NA does.
    average = 'agreement_average', averagemeasures = 'agreement_average',
    averagemeasure = 'agreement_average', avg = 'agreement_average',
    mean = 'agreement_average', icc2k = 'agreement_average',
    agreementaverage = 'agreement_average',
    absoluteagreementaverage = 'agreement_average',
    twowayrandomaverage = 'agreement_average',
    icc3k = 'consistency_average', consistencyaverage = 'consistency_average',
    twowaymixedaverage = 'consistency_average'
  )
  out <- unname(map[key])
  bad <- which(is.na(out) & !is.na(raw) & nzchar(key))
  if (length(bad) && isTRUE(warn)) {
    warning(paste0(
      "Unrecognised icc_type value(s): '", paste(unique(raw[bad]), collapse = "', '"),
      "'. Treated as 'agreement'. Recognised values are 'agreement' and ",
      "'consistency' (case-insensitive; 'ICC(2,1)' / 'ICC(3,1)' / 'absolute ",
      "agreement' / 'two-way random' / 'two-way mixed' are also accepted), plus ",
      "the average-measures forms 'average' / 'ICC(2,k)' / 'ICC(3,k)'."))
  }
  out[is.na(out)] <- 'agreement'
  out
}


# TRUE for the absolute-agreement family, single- or average-measures. V31 keys on
# this: after a Spearman-Brown step-down an ICC(2,k) row is computed with the SAME
# one-way agreement SE approximation, so it needs the same note.
.icc_is_agreement <- function(type) {
  type %in% c('agreement', 'agreement_average')
}


# TRUE for the average-measures forms, i.e. the value describes the mean of k
# measurements rather than one.
.icc_is_average <- function(type) {
  type %in% c('agreement_average', 'consistency_average')
}


# Spearman-Brown, inverted: recover the single-measurement ICC from an
# average-of-k ICC.  ICC_k = k*ICC_1 / (1 + (k-1)*ICC_1)  =>
.icc_step_down <- function(rho_k, k) {
  # Denominator k - (k-1)*rho_k is >= 1 for any rho_k <= 1 and k >= 2, so it cannot
  # vanish for an in-range ICC (V11 already bounds icc to [-1, 1]).
  rho_k / (k - (k - 1) * rho_k)
}
