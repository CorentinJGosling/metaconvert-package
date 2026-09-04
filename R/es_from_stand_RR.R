#' Convert a risk ratio value and standard error to various effect size measures
#'
#' @param rr risk ratio value
#' @param logrr log risk ratio value
#' @param logrr_se standard error of the log risk ratio
#' @param n_cases number of cases/events
#' @param n_controls number of controls/no-event
#' @param n_exp number of participants in the exposed group
#' @param n_nexp number of participants in the non-exposed group
#' @param baseline_risk proportion of cases in the non-exposed group (required for the \code{rr_to_or = "grant"} argument, and for the risk difference and the NNT). It is never derived from other columns.
#' @param reverse_rr a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param rr_to_or formula used to convert the \code{rr} value into an odds ratio (see details).
#'
#' @details
#' This function converts the (log) risk ratio (RR) value and its standard error
#' to odds ratio (OR) and number needed to treat.
#'
#' **To estimate the odds ratio and its standard error**, various formulas can be used.
#'
#' **A.** First, the approach identified by the \code{rr_to_or = "grant"} argument value
#' can be used. The estimator is the inverse of the OR-to-RR conversion due to Zhang and
#' Yu (1998); Grant (2014) restates the same expression and is cited here for the
#' communication framing, not as the source of the formula. The argument value
#' \code{"grant"} is retained for backward compatibility.
#' Neither paper gives a variance, standard error or confidence interval for the converted
#' value. To derive the variance, we used this formula to convert the bounds of the 95% CI,
#' which were then used to obtain the variance.
#'
#' This value requires (rr + logrr_se + baseline_risk) to generate an OR.
#' The 95% CI of the RR is reconstructed internally from \code{rr} and \code{logrr_se}.
#' The following formulas are used (br = baseline_risk):
#' \deqn{or = \frac{rr * (1 - br)}{1 - rr * br}}
#' \deqn{or\_ci\_lo = \frac{rr\_ci\_lo * (1 - br)}{1 - rr\_ci\_lo * br}}
#' \deqn{or\_ci\_up = \frac{rr\_ci\_up * (1 - br)}{1 - rr\_ci\_up * br}}
#' \deqn{logor\_se = \frac{log(or\_ci\_up) - log(or\_ci\_lo)}{2 * qnorm(.975)}}
#'
#' **B.** Second, the formulas implemented in the metaumbrella package can be used (\code{rr_to_or = "metaumbrella"}, the default).
#' This value requires (rr + logrr_se + n_cases + n_controls) to generate an OR.
#' More precisely, we previously developed functions that simulate all combinations of the possible number of cases and controls
#' in the exposed and non-exposed groups compatible with the actual value of the RR.
#' Then, the functions select the contingency table whose standard error coincides best with the standard error reported.
#' The RR value and its standard are obtained from this estimated contingency table.
#'
#' **C.** Third, it is possible to transpose the RR to a OR (\code{rr_to_or = "transpose"}).
#' This value requires (rr + logrr_se) to generate an OR.
#' It is known that OR and RR are similar when the baseline risk is small.
#' Therefore, users can request to simply transpose the RR value & standard error into a OR value & standard error.
#' \deqn{or = rr}
#' \deqn{logor\_se = logrr\_se}
#'
#' **D.** Fourth, it is possible to recreate the 2x2 table using the dipietrantonj's formulas (\code{rr_to_or = "dipietrantonj"}).
#' This value requires (rr + logrr_se + n_exp + n_nexp) to generate an OR. Information on this approach can be retrieved in
#' Di Pietrantonj (2006).
#'
#' **To estimate the NNT**, the formulas used are :
#' \deqn{nnt = \frac{1}{br * (1 - rr)}}
#'
#' @references
#' Di Pietrantonj C. (2006). Four-fold table cell frequencies imputation in meta analysis. Statistics in medicine, 25(13), 2299-2322. https://doi.org/10.1002/sim.2287
#'
#' Gosling, C. J., Solanes, A., Fusar-Poli, P., & Radua, J. (2023). metaumbrella: the first comprehensive suite to perform data analysis in umbrella reviews with stratification of the evidence. BMJ mental health, 26(1), e300534. https://doi.org/10.1136/bmjment-2022-300534
#'
#' Grant R. L. (2014). Converting an odds ratio to a range of plausible relative risks for better communication of research findings. BMJ (Clinical research ed.), 348, f7450. https://doi.org/10.1136/bmj.f7450
#'
#' Veroniki, A. A., Pavlides, M., Patsopoulos, N. A., & Salanti, G. (2013). Reconstructing 2x2 contingency tables from odds ratios using the Di Pietrantonj method: difficulties, constraints and impact in meta-analysis results. Research synthesis methods, 4(1), 78-94. https://doi.org/10.1002/jrsm.1061
#'
#' Zhang, J., & Yu, K. F. (1998). What's the relative risk? A method of correcting the odds ratio in cohort studies of common outcomes. JAMA, 280(19), 1690-1691. https://doi.org/10.1001/jama.280.19.1690
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab RR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + NNT + RD\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 3. Risk Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_rr_se
#'
#' @md
#'
#' @examples
#' es_from_rr_se(rr = 2.12, logrr_se = 0.242, n_cases = 120, n_controls = 44,
#'               baseline_risk = 0.20)
es_from_rr_se <- function(rr, logrr, logrr_se, baseline_risk,
                          n_exp, n_nexp, n_cases, n_controls,
                          rr_to_or = "metaumbrella",
                          reverse_rr) {

  if (missing(rr)) rr <- rep(NA_real_, length(logrr))
  if (missing(logrr)) logrr <- rep(NA_real_, length(rr))
  if (missing(logrr_se)) logrr_se <- rep(NA_real_, length(rr))
  # Propagated unchanged into logrr_se and the OR/NNT conversions; a non-positive
  # value becomes a negative sampling variance. See R/internal_guards.R.
  logrr_se <- .positive_or_na(logrr_se)
  if (missing(baseline_risk)) {
    baseline_risk <- rep(NA_real_, length(rr))
  }
  # Divided by, not merely reported: outside [0, 1) both (1 - BR) and (1 - RR*BR)
  # change sign together, so the grant conversion stays finite and positive while
  # returning a negative standard error and a transposed interval. See
  # R/internal_guards.R.
  baseline_risk <- .baseline_risk_or_na(baseline_risk)
  # A ratio must be strictly positive. The log-scale outputs self-blank on rr <= 0,
  # but the risk-difference block does not: it returns a fabricated rd with a negative
  # rd_se and a transposed CI (and rd_se exactly 0 at rr = 0).
  rr <- .ratio_or_na(rr)
  if (missing(n_exp)) {
    n_exp <- rep(NA_real_, length(rr))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA_real_, length(rr))
  }
  if (missing(n_cases)) {
    n_cases <- rep(NA_real_, length(rr))
  }
  if (missing(n_controls)) {
    n_controls <- rep(NA_real_, length(rr))
  }
  if (missing(reverse_rr)) reverse_rr <- rep(FALSE, length(rr))
  reverse_rr[is.na(reverse_rr)] <- FALSE
  if (length(reverse_rr) == 1) reverse_rr = c(rep(reverse_rr, length(rr)))
  if (length(reverse_rr) != length(rr)) stop("The length of the 'reverse_rr' argument is incorrectly specified.")

  if (!all(rr_to_or %in% c("metaumbrella", "transpose",
                           "grant", "dipietrantonj"))) {
    stop(paste0("'",
                unique(rr_to_or[!rr_to_or %in% c("metaumbrella", "transpose", "grant", "dipietrantonj")]),
                "' not in tolerated values for the 'rr_to_or' argument.
                Possible inputs are: 'metaumbrella', 'transpose', 'grant', 'dipietrantonj'"))
  }


  rr <- ifelse(is.na(rr) & !is.na(logrr), exp(logrr), rr)
  # rr <- ifelse(reverse_rr, 1 / rr, rr)

  # RR -------
  es <- data.frame(
    logrr = ifelse(reverse_rr, -log(rr), log(rr)),
    logrr_se = logrr_se,
    logrr_ci_lo = ifelse(reverse_rr, -log(rr) - qnorm(.975) * logrr_se,
                                      log(rr) - qnorm(.975) * logrr_se),
    logrr_ci_up = ifelse(reverse_rr, -log(rr) + qnorm(.975) * logrr_se,
                                      log(rr) + qnorm(.975) * logrr_se)
  )
  rr_ci_lo <- exp(log(rr) - qnorm(.975) * logrr_se)
  rr_ci_up <- exp(log(rr) + qnorm(.975) * logrr_se)

  # OR -------
  es$logor <- es$logor_se <- es$logor_ci_lo <- es$logor_ci_up <- NA

  dat_or <- data.frame(
    rr = rr, logrr_se = logrr_se, rr_ci_lo = rr_ci_lo, rr_ci_up = rr_ci_up,
    n_cases = n_cases, n_controls = n_controls, n_exp = n_exp, n_nexp = n_nexp,
    baseline_risk = baseline_risk, rr_to_or = rr_to_or
  )

  nn_miss <- with(dat_or, which(
      (rr_to_or == "grant" & !is.na(rr) & !is.na(baseline_risk) & !is.na(rr_ci_lo) & !is.na(rr_ci_up)) |
      (rr_to_or == "metaumbrella" & !is.na(rr) & !is.na(logrr_se) & !is.na(n_cases) & !is.na(n_controls)) |
      (rr_to_or == "transpose" & !is.na(rr) & !is.na(logrr_se)) |
      (rr_to_or == "dipietrantonj" & !is.na(rr) & !is.na(rr_ci_lo) & !is.na(rr_ci_up) & !is.na(n_exp) & !is.na(n_nexp))
  ))


  if (length(nn_miss) != 0) {
    res_or <- t(mapply(.rr_to_or,
      rr = dat_or$rr[nn_miss],
      logrr_se = dat_or$logrr_se[nn_miss],
      rr_ci_lo = dat_or$rr_ci_lo[nn_miss],
      rr_ci_up = dat_or$rr_ci_up[nn_miss],
      n_cases = dat_or$n_cases[nn_miss],
      n_controls = dat_or$n_controls[nn_miss],
      n_exp = dat_or$n_exp[nn_miss],
      n_nexp = dat_or$n_nexp[nn_miss],
      baseline_risk = dat_or$baseline_risk[nn_miss],
      rr_to_or = dat_or$rr_to_or[nn_miss]
    ))

    es$logor[nn_miss] <- ifelse(reverse_rr[nn_miss], -res_or[, 1], res_or[, 1])
    es$logor_se[nn_miss] <- res_or[, 2]
    # On reverse the point estimate is negated, so the symmetric, log-scale CI must be
    # negated and swapped: new_lo = -old_up, new_up = -old_lo. Swapping alone leaves a
    # wrong-signed, inverted interval (lo > up) that does not bracket the reversed
    # estimate.
    es$logor_ci_lo[nn_miss] <- ifelse(reverse_rr[nn_miss], -res_or[, 4], res_or[, 3])
    es$logor_ci_up[nn_miss] <- ifelse(reverse_rr[nn_miss], -res_or[, 3], res_or[, 4])
  }

  # Risk difference from RR + baseline_risk
  # The risk difference is baseline_risk - (rr * baseline_risk), so it presumes the
  # implied exposed risk rr * baseline_risk is itself a probability. When it exceeds 1
  # the pair of risks does not exist and the RD and NNT describe nothing: at rr = 2 with
  # baseline_risk = 0.9 the implied exposed risk is 1.8, and the route returns
  # rd = -0.900 and nnt = -1.111 while correctly returning the odds ratio as NA. Neither
  # value trips a flag, since B6 tests |rd| > 1 and |-0.9| is inside the bound, so the
  # impossible pair would reach the output unremarked.
  #
  # NB the OR path already declines here (Grant's transform needs rr * baseline_risk < 1
  # and returns NA outside it), so this closes the one route that did not.
  .risk_pair_ok <- is.na(rr) | is.na(baseline_risk) | (rr * baseline_risk <= 1)
  rd <- ifelse(.risk_pair_ok, baseline_risk * (1 - rr), NA_real_)

  # delta method
  rd_se <- ifelse(.risk_pair_ok, baseline_risk * rr * logrr_se, NA_real_)

  # A derived risk-difference SE of exactly 0 is a zero sampling variance -- an infinite
  # inverse-variance weight, or an rma() abort. It arises at baseline_risk = 0, where
  # rd_se = BR * rr * logrr_se collapses. Given the delta method own assumption that
  # baseline_risk is a known constant, rd = 0 IS the correct conditional point estimate
  # there; what is wrong is shipping it as a poolable row, so the row declines both.
  # The guard is keyed on the SE and never on rd itself: a genuinely null risk
  # difference has rd = 0 with a perfectly good standard error and must survive.
  # baseline_risk = 0 stays legal on the log scale (the rare-disease limit, where the
  # grant conversions become the identity). See R/internal_guards.R.
  degenerate_rd <- !is.na(rd_se) & rd_se <= 0
  rd_se <- .positive_or_na(rd_se)
  rd <- ifelse(degenerate_rd, NA_real_, rd)

  es$rd <- ifelse(reverse_rr, -rd, rd)
  es$rd_se <- rd_se
  es$rd_ci_lo <- es$rd - qnorm(.975) * rd_se
  es$rd_ci_up <- es$rd + qnorm(.975) * rd_se

  es$nnt <- ifelse(rd == 0, NA, 1 / rd)
  es$nnt <- ifelse(reverse_rr, -es$nnt, es$nnt)
  es$nnt_se <- ifelse(rd == 0, NA, rd_se / rd^2)
  rd_ci_lo_raw <- rd - qnorm(.975) * rd_se
  rd_ci_up_raw <- rd + qnorm(.975) * rd_se
  # Non-strict: a bound landing exactly on 0 is still the Altman discontinuity, and
  # the reciprocal of +0 is a literal Inf where the neighbouring input returns NA.
  crosses_zero <- (rd_ci_lo_raw <= 0 & rd_ci_up_raw >= 0) | rd == 0
  es$nnt_ci_lo <- ifelse(crosses_zero, NA,
                          ifelse(reverse_rr, -1 / rd_ci_lo_raw, 1 / rd_ci_up_raw))
  es$nnt_ci_up <- ifelse(crosses_zero, NA,
                          ifelse(reverse_rr, -1 / rd_ci_up_raw, 1 / rd_ci_lo_raw))

  es$info_used <- "rr_se"
  return(es)
}


#' Convert a risk ratio value and 95% confidence interval to various effect size measures
#'
#' @param rr risk ratio value
#' @param logrr log risk ratio value
#' @param rr_ci_lo lower bound of the 95% CI around the risk ratio value
#' @param rr_ci_up upper bound of the 95% CI around the risk ratio value
#' @param logrr_ci_lo lower bound of the 95% CI around the log risk ratio value
#' @param logrr_ci_up upper bound of the 95% CI around the log risk ratio value
#' @param n_cases number of cases/events (required for \code{rr_to_or = "metaumbrella"}, the default)
#' @param n_controls number of controls/no-event (required for \code{rr_to_or = "metaumbrella"}, the default)
#' @param n_exp number of participants in the exposed group (required for \code{rr_to_or = "dipietrantonj"})
#' @param n_nexp number of participants in the non-exposed group (required for \code{rr_to_or = "dipietrantonj"})
#' @param baseline_risk proportion of cases in the non-exposed group (required for \code{rr_to_or = "grant"}, and for the risk difference and the NNT)
#' @param reverse_rr a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param rr_to_or formula used to convert the \code{rr} value into an odds ratio (see details).
#'
#' @details
#' This function uses the 95% CI of the (log) risk ratio to obtain the standard error (Section 6.5.2.2 in the Cochrane Handbook).
#' \deqn{logrr\_se = \frac{\log{rr\_ci\_up} - \log{rr\_ci\_lo}}{2 * qnorm(.975)}}
#'
#' Then, calculations of the \code{\link{es_from_rr_se}()} are applied.
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab RR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + NNT + RD\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 3. Risk Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_rr_ci
#'
#' @md
#'
#' @examples
#' es_from_rr_ci(
#'   rr = 1, rr_ci_lo = 0.5, rr_ci_up = 2,
#'   n_cases = 42, n_controls = 38, baseline_risk = 0.08
#' )
es_from_rr_ci <- function(rr, rr_ci_lo, rr_ci_up, logrr, logrr_ci_lo, logrr_ci_up, baseline_risk,
                          n_exp, n_nexp, n_cases, n_controls, rr_to_or = "metaumbrella",
                          max_asymmetry = 10, reverse_rr) {
  if (missing(rr)) {
    rr <- rep(NA_real_, length(logrr))
  }
  if (missing(logrr)) {
    logrr <- rep(NA_real_, length(rr))
  }
  if (missing(baseline_risk)) {
    baseline_risk <- rep(NA_real_, length(rr))
  }
  # Divided by, not merely reported: outside [0, 1) both (1 - BR) and (1 - RR*BR)
  # change sign together, so the grant conversion stays finite and positive while
  # returning a negative standard error and a transposed interval. See
  # R/internal_guards.R.
  baseline_risk <- .baseline_risk_or_na(baseline_risk)
  # A ratio must be strictly positive. The log-scale outputs self-blank on rr <= 0,
  # but the risk-difference block does not: it returns a fabricated rd with a negative
  # rd_se and a transposed CI (and rd_se exactly 0 at rr = 0).
  rr <- .ratio_or_na(rr)
  if (missing(n_exp)) {
    n_exp <- rep(NA_real_, length(rr))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA_real_, length(rr))
  }
  if (missing(n_cases)) {
    n_cases <- rep(NA_real_, length(rr))
  }
  if (missing(n_controls)) {
    n_controls <- rep(NA_real_, length(rr))
  }
  if (missing(rr_ci_lo)) {
    rr_ci_lo <- rep(NA_real_, length(rr))
  }
  if (missing(rr_ci_up)) {
    rr_ci_up <- rep(NA_real_, length(rr))
  }
  if (missing(logrr_ci_lo)) {
    logrr_ci_lo <- rep(NA_real_, length(rr))
  }
  if (missing(logrr_ci_up)) {
    logrr_ci_up <- rep(NA_real_, length(rr))
  }
  if (missing(reverse_rr)) {
    reverse_rr <- rep(FALSE, length(rr))
  }
  reverse_rr[is.na(reverse_rr)] <- FALSE

  rr <- ifelse(is.na(rr) & !is.na(logrr), exp(logrr), rr)
  logrr_ci_lo <- ifelse(is.na(logrr_ci_lo) & !is.na(rr_ci_lo), log(rr_ci_lo), logrr_ci_lo)
  logrr_ci_up <- ifelse(is.na(logrr_ci_up) & !is.na(rr_ci_up), log(rr_ci_up), logrr_ci_up)

  logrr_se <- .ci_width(logrr_ci_lo, logrr_ci_up) / (2 * qnorm(.975))
  es <- es_from_rr_se(
    rr = rr, logrr_se = logrr_se,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls,
    rr_to_or = rr_to_or, reverse_rr = reverse_rr
  )

  es$info_used <- "rr_ci"

  return(es)
}

#' Convert a risk ratio value and its p-value to various effect size measures
#'
#' @param rr risk ratio value
#' @param logrr log risk ratio value
#' @param rr_pval p-value of the risk ratio
#' @param n_cases number of cases/events (required for \code{rr_to_or = "metaumbrella"}, the default)
#' @param n_controls number of controls/no-event (required for \code{rr_to_or = "metaumbrella"}, the default)
#' @param n_exp number of participants in the exposed group (required for \code{rr_to_or = "dipietrantonj"})
#' @param n_nexp number of participants in the non-exposed group (required for \code{rr_to_or = "dipietrantonj"})
#' @param baseline_risk proportion of cases in the non-exposed group (required for \code{rr_to_or = "grant"}, and for the risk difference and the NNT)
#' @param reverse_rr_pval a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param rr_to_or formula used to convert the \code{rr} value into an odds ratio (see details).
#'
#' @details
#' This function uses the p-value of the (log) risk ratio to obtain the standard error (Section 6.3.2 in the Cochrane Handbook).
#' \deqn{logrr\_z = qnorm(rr_pval/2, lower.tail=FALSE)}
#' \deqn{logrr\_se = |\frac{\log(rr)}{logrr\_z}|}
#'
#' Then, calculations of \code{\link{es_from_rr_se}} are applied.
#'
#' @references
#' Higgins, J. P., Thomas, J., Chandler, J., Cumpston, M., Li, T., Page, M. J., & Welch, V. A. (Eds.). (2019). Cochrane handbook for systematic reviews of interventions. John Wiley & Sons.
#'
#' @export es_from_rr_pval
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab RR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + NNT + RD\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 3. Risk Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_rr <- es_from_rr_pval(
#'   rr = 3.51, rr_pval = 0.001,
#'   n_cases = 12, n_controls = 68
#' )
es_from_rr_pval <- function(rr, logrr, rr_pval, baseline_risk,
                            n_exp, n_nexp, n_cases, n_controls,
                            rr_to_or = "metaumbrella",
                            reverse_rr_pval) {
  if (missing(rr)) rr <- rep(NA_real_, length(logrr))
  if (missing(logrr)) logrr <- rep(NA_real_, length(rr))
  if (missing(baseline_risk)) {
    baseline_risk <- rep(NA_real_, length(rr))
  }
  # Divided by, not merely reported: outside [0, 1) both (1 - BR) and (1 - RR*BR)
  # change sign together, so the grant conversion stays finite and positive while
  # returning a negative standard error and a transposed interval. See
  # R/internal_guards.R.
  baseline_risk <- .baseline_risk_or_na(baseline_risk)
  # A ratio must be strictly positive. The log-scale outputs self-blank on rr <= 0,
  # but the risk-difference block does not: it returns a fabricated rd with a negative
  # rd_se and a transposed CI (and rd_se exactly 0 at rr = 0).
  rr <- .ratio_or_na(rr)
  if (missing(n_exp)) {
    n_exp <- rep(NA_real_, length(rr))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA_real_, length(rr))
  }
  if (missing(n_cases)) {
    n_cases <- rep(NA_real_, length(rr))
  }
  if (missing(n_controls)) {
    n_controls <- rep(NA_real_, length(rr))
  }
  if (missing(reverse_rr_pval)) reverse_rr_pval <- rep(FALSE, length(rr))
  reverse_rr_pval[is.na(reverse_rr_pval)] <- FALSE

  rr <- ifelse(is.na(rr) & !is.na(logrr), exp(logrr), rr)
  # Mirror es_from_or_pval: take the positive critical value and the absolute
  # magnitude, so rr == 1 (log(rr) = 0) evaluates to 0 here rather than to the 0/0 =
  # NaN a sign(log(rr)) denominator would produce. The direction of the effect is
  # carried by rr itself, so the SE only needs its magnitude.
  #
  # What this function returns at rr == 1 is NA rather than 0; the 0 above describes
  # the intermediate value. es_from_rr_se() passes it through .positive_or_na(), and a
  # zero SE is non-positive: it would be an infinite inverse-variance weight. NA is the
  # intended outcome, since rr = 1 beside a large p-value carries no information about
  # the standard error at all.
  # Measured and pinned in tests/testthat/test-unverified-routes-external.R, which
  # also checks that es_from_or_pval behaves identically, as "mirror" claims.
  rr_pval <- .pval_or_na(rr_pval)  # p <= 0 or p >= 1: no finite se, see .pval_or_na()
  z_rr <- qnorm(rr_pval / 2, lower.tail = FALSE)
  logrr_se <- abs(log(rr) / z_rr)

  es <- es_from_rr_se(
    rr = rr, logrr_se = logrr_se, n_cases = n_cases,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_controls = n_controls,
    rr_to_or = rr_to_or, reverse_rr = reverse_rr_pval
  )

  es$info_used <- "rr_pval"

  return(es)
}
