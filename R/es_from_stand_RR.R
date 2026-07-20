#' Convert a risk ratio value and standard error to various effect size measures
#'
#' @param rr risk ratio value
#' @param logrr log risk ratio value
#' @param logrr_se standard error of the log risk ratio
#' @param n_cases number of cases/events
#' @param n_controls number of controls/no-event
#' @param n_exp number of participants in the exposed group
#' @param n_nexp number of participants in the non-exposed group
#' @param baseline_risk proportion of cases in the non-exposed group (n_cases_nexp / n_nexp is used when missing)
#' @param reverse_rr a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param rr_to_or formula used to convert the \code{rr} value into an odds ratio (see details).
#'
#' @details
#' This function converts the (log) risk ratio (RR) value and its standard error
#' to odds ratio (OR) and number needed to treat.
#'
#' **To estimate the odds ratio and its standard error**, various formulas can be used.
#'
#' **A.** First, the approach described in Grant (2014) can be used.
#' However, in the paper, only the formula to convert an RR value to a OR value
#' is described. To derive the variance, we used this formula to convert the bounds of the 95% CI, which
#' were then used to obtain the variance.
#'
#' This argument requires (rr + baseline_risk + rr_ci_lo + rr_ci_up) to generate a RR.
#' The following formulas are used (br = baseline_risk):
#' \deqn{or = \frac{rr * (1 - br)}{1 - rr * br}}
#' \deqn{or\_ci\_lo = \frac{rr\_ci\_lo * (1 - br)}{1 - rr\_ci\_lo * br}}
#' \deqn{or\_ci\_up = \frac{rr\_ci\_up * (1 - br)}{1 - rr\_ci\_up * br}}
#' \deqn{logor\_se = \frac{log(or\_ci\_up) - log(or\_ci\_lo)}{2 * qnorm(.975)}}
#'
#' **B.** Second, the formulas implemented in the metaumbrella package can be used (\code{or_to_rr = "metaumbrella_exp"}).
#' This argument requires (rr + logrr_se + n_exp + n_nexp) to generate a OR.
#' More precisely, we previously developed functions that simulate all combinations of the possible number of cases and controls
#' in the exposed and non-exposed groups compatible with the actual value of the RR.
#' Then, the functions select the contingency table whose standard error coincides best with the standard error reported.
#' The RR value and its standard are obtained from this estimated contingency table.
#'
#' **C.** Third, it is possible to transpose the RR to a OR (\code{rr_to_or = "transpose"}).
#' This argument requires (rr + logrr_se) to generate a OR.
#' It is known that OR and RR are similar when the baseline risk is small.
#' Therefore, users can request to simply transpose the RR value & standard error into a OR value & standard error.
#' \deqn{or = rr}
#' \deqn{logor\_se = logrr\_se}
#'
#' **D.** Fourth, it is possible to recreate the 2x2 table using the dipietrantonj's formulas (\code{rr_to_or = "dipietrantonj"}).
#' This argument requires (rr + logrr_ci_lo + logrr_ci_lo) to generate a OR. Information on this approach can be retrieved in
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
#' es_from_rr_se(rr = 2.12, logrr_se = 0.242, n_exp = 120, n_nexp = 44)
es_from_rr_se <- function(rr, logrr, logrr_se, baseline_risk,
                          n_exp, n_nexp, n_cases, n_controls,
                          rr_to_or = "metaumbrella",
                          reverse_rr) {

  if (missing(rr)) rr <- rep(NA_real_, length(logrr))
  if (missing(logrr)) logrr <- rep(NA_real_, length(rr))
  if (missing(logrr_se)) logrr_se <- rep(NA_real_, length(rr))
  if (missing(baseline_risk)) {
    baseline_risk <- rep(NA_real_, length(rr))
  }
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
    # On reverse the point estimate is negated, so the (symmetric, log-scale) CI must be
    # NEGATED AND SWAPPED: new_lo = -old_up, new_up = -old_lo. Swapping alone left a
    # wrong-signed, inverted interval (lo > up) that did not bracket the reversed estimate.
    es$logor_ci_lo[nn_miss] <- ifelse(reverse_rr[nn_miss], -res_or[, 4], res_or[, 3])
    es$logor_ci_up[nn_miss] <- ifelse(reverse_rr[nn_miss], -res_or[, 3], res_or[, 4])
  }

  # Risk difference from RR + baseline_risk
  rd <- baseline_risk * (1 - rr)

  es$rd <- ifelse(reverse_rr, -rd, rd)
  # delta method
  rd_se <- baseline_risk * rr * logrr_se
  es$rd_se <- rd_se
  es$rd_ci_lo <- es$rd - qnorm(.975) * rd_se
  es$rd_ci_up <- es$rd + qnorm(.975) * rd_se

  es$nnt <- ifelse(rd == 0, NA, 1 / rd)
  es$nnt <- ifelse(reverse_rr, -es$nnt, es$nnt)
  es$nnt_se <- ifelse(rd == 0, NA, rd_se / rd^2)
  rd_ci_lo_raw <- rd - qnorm(.975) * rd_se
  rd_ci_up_raw <- rd + qnorm(.975) * rd_se
  crosses_zero <- (rd_ci_lo_raw < 0 & rd_ci_up_raw > 0) | rd == 0
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
#' @param n_cases number of cases/events
#' @param n_controls number of controls/no-event
#' @param n_exp number of participants in the exposed group (only required for the \code{rr_to_or = "grant_CI"}, \code{rr_to_or = "grant_2x2"} arguments).
#' @param n_nexp number of participants in the non-exposed group (only required for the \code{rr_to_or = "grant_CI"}, \code{rr_to_or = "grant_2x2"} arguments).
#' @param baseline_risk proportion of cases in the non-exposed group (only required for the \code{rr_to_or = "grant_CI"} and \code{rr_to_or = "grant_2x2"} arguments).
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
#'  \code{converted effect size measure} \tab OR + NNT\cr
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

  logrr_se <- (logrr_ci_up - logrr_ci_lo) / (2 * qnorm(.975))
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
#' @param n_cases number of cases/events
#' @param n_controls number of controls/no-event
#' @param n_exp number of participants in the exposed group (only required for the \code{rr_to_or = "grant_CI"}, \code{rr_to_or = "grant_2x2"} arguments).
#' @param n_nexp number of participants in the non-exposed group (only required for the \code{rr_to_or = "grant_CI"}, \code{rr_to_or = "grant_2x2"} arguments).
#' @param baseline_risk proportion of cases in the non-exposed group (only required for the \code{rr_to_or = "grant_CI"} and \code{rr_to_or = "grant_2x2"} arguments).
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
#'  \code{converted effect size measure} \tab OR + NNT\cr
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
  # magnitude, so rr == 1 (log(rr) = 0) gives logrr_se = 0 (finite) rather than the
  # 0/0 = NaN produced by the earlier sign(log(rr)) denominator. The direction of the
  # effect is carried by rr itself, so the SE only needs its magnitude.
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
