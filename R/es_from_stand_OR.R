#' Convert an odds ratio value and its standard error into several effect size measures
#'
#' @param or odds ratio value
#' @param logor log odds ratio value
#' @param logor_se the standard error of the log odds ratio
#' @param n_cases number of cases/events across exposed/non-exposed groups
#' @param n_controls number of controls/no-event across exposed/non-exposed groups
#' @param n_exp number of participants in the exposed group
#' @param n_nexp number of participants in the non-exposed group
#' @param n_sample total number of participants in the sample
#' @param baseline_risk proportion of cases in the non-exposed group
#' @param small_margin_prop smallest margin proportion of cases/events in the underlying 2x2 table
#' @param reverse_or a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param or_to_rr formula used to convert the \code{or} value into a risk ratio (see details).
#' @param or_to_cor formula used to convert the \code{or} value into a correlation coefficient (see details).
#'
#' @details
#' This function converts the log odds ratio into a Risk ratio (RR), Cohen's d (D), Hedges' g (G)
#' and correlation coefficients (R/Z).
#'
#' **To estimate the Cohen's d value and its standard error**
#' The following formulas are used (Cooper et al., 2019):
#' \deqn{d = \log(or) * \frac{\sqrt{3}}{\pi}}
#' \deqn{d\_se = \sqrt{\frac{logor\_se^2 * 3}{\pi^2}}}
#'
#' **To estimate the risk ratio and its standard error, various formulas can be used.**
#'
#' **A.** First, the approach described in Grant (2014) can be used.
#' However, in the paper, only the formula to convert an OR value to a RR value
#' is described.
#' To derive the variance, we used this formula to convert the bounds of the 95% CI, which
#' were then used to obtain the variance.
#'
#' This argument requires (or + baseline_risk + or_ci_lo + or_ci_up) to generate a RR.
#' The following formulas are used (br = baseline_risk):
#' \deqn{rr = \frac{or}{1 - br + br*or}}
#' \deqn{rr\_ci\_lo = \frac{or\_ci\_lo}{1 - br + br*or\_ci\_lo}}
#' \deqn{rr\_ci\_up = \frac{or\_ci\_up}{1 - br + br*or\_ci\_up}}
#' \deqn{logrr\_se = \frac{log(rr\_ci\_up) - log(rr\_ci\_lo)}{2 * qnorm(.975)}}
#'
#'
#' **B.** Second, the formulas implemented in the metaumbrella package can be used
#' (\code{or_to_rr = "metaumbrella_cases"} or \code{or_to_rr = "metaumbrella_exp"}).
#' This argument requires (or + logor_se + n_cases + n_controls) or (or + logor_se + n_exp + n_nexp)
#' to generate a RR.
#' More precisely, when the OR value and its standard error, plus either
#' (i) the number of cases and controls or
#' (ii) the number of participants in the exposed and non-exposed groups,
#' are available, we previously developed functions that simulate all combinations of the possible
#' number of cases and controls
#' in the exposed and non-exposed groups compatible with the actual value of the OR.
#' Then, the functions select the contingency table whose standard error coincides best with
#' the standard error reported.
#' The RR value and its standard are obtained from this estimated contingency table.
#'
#' **C.** Third, it is possible to transpose the RR to a OR (\code{or_to_rr = "transpose"}).
#' This argument requires (or + logor_se) to generate a OR.
#' It is known that OR and RR are similar when the baseline risk is small.
#' Therefore, users can request to simply transpose the OR value & standard error into a RR value & standard error.
#' \deqn{rr = or}
#' \deqn{logrr\_se = logor\_se}
#'
#' **D.** Fourth, it is possible to recreate the 2x2 table using the dipietrantonj's formulas (\code{or_to_rr = "dipietrantonj"}).
#' This argument requires (or + logor_ci_lo + logor_ci_lo) to generate a RR. Information on this approach can be retrieved in
#' Di Pietrantonj (2006).
#'
#' **To estimate the NNT, the formulas used are :**
#' \deqn{\frac{(1 - br * (1 - or))}{(1 - br) * (br * (1 - or))}}
#'
#' **To estimate a correlation coefficient, various formulas can be used.**
#'
#' **A.** First, the approach described in Pearson (1900) can be used (\code{or_to_cor = "pearson"}).
#' This argument requires (or + logor_se) to generate a R/Z.
#' It converts the OR value and its standard error to a tetrachoric correlation.
#' Note that the formula assumes that each cell of the 2x2 used to estimate the OR has been added 1/2 before estimating the OR value and its standard error.
#' If it is not the case, formulas can produce slightly less accurate results.
#'
#' \deqn{c = \frac{1}{2}}
#' \deqn{r = \cos{\frac{\pi}{1+or^c}}}
#' \deqn{r\_se = logor\_se * ((\pi * c * or^c) * \frac{\sin(\pi / (1+or^c))}{1+or^c})^2}
#' \deqn{or\_ci\_lo = exp(log(or) - qnorm(.975)*logor\_se)}
#' \deqn{or\_ci\_up = exp(log(or) + qnorm(.975)*logor\_se)}
#' \deqn{r\_ci\_lo = cos(\frac{\pi}{1 + or\_ci\_lo^c})}
#' \deqn{r\_ci\_up = cos(\frac{\pi}{1 + or\_ci\_up^c})}
#' \deqn{z = atanh(r)}
#' \deqn{z\_se = \sqrt{\frac{r\_se^2}{(1 - r^2)^2}}}
#' \deqn{z\_ci\_lo = atanh(r\_lo)}
#' \deqn{z\_ci\_up = atanh(r\_up)}
#'
#' **B.** Second, the approach described in Digby (1983) can be used (\code{or_to_cor = "digby"}).
#' This argument requires (or + logor_se) to generate a R/Z.
#' It converts the OR value and its standard error to a tetrachoric correlation.
#' Note that the formula assumes that each cell of the 2x2 used to estimate the OR has been added 1/2 before estimating the OR value and its standard error.
#' If it is not the case, formulas can produce slightly less accurate results.
#'
#' \deqn{c = \frac{3}{4}}
#' \deqn{r = \frac{or^c - 1}{or^c + 1}}
#' \deqn{r\_se = \sqrt{\frac{c^2}{4} * (1 - r^2)^2 * logor\_se}}
#' \deqn{z = atanh(r)}
#' \deqn{z\_se = \sqrt{\frac{r\_se^2}{(1 - r^2)^2}}}
#' \deqn{z\_ci\_lo = z - qnorm(.975)*\sqrt{\frac{c^2}{4} * logor\_se}}
#' \deqn{z\_ci\_up = z + qnorm(.975)*\sqrt{\frac{c^2}{4} * logor\_se}}
#' \deqn{r\_ci\_lo = tanh(z\_lo)}
#' \deqn{r\_ci\_up = tanh(z\_up)}
#'
#' **C.** Third, the approach described in Bonett (2005) can be used (\code{or_to_cor = "bonett"}).
#' This argument requires (or + logor_se + n_cases + n_exp + small_margin_prop) to generate a R/Z.
#' Note that the formula assumes that each cell of the 2x2 used to estimate the OR has been added 1/2 before estimating the OR value and its standard error.
#' If it is not the case, formulas can produce slightly less accurate results.
#'
#' \deqn{c = \frac{\frac{1 - |n\_exp - n\_cases|}{5} - (0.5 - small\_margin\_prop)^2}{2}}
#' \deqn{r = \cos{\frac{\pi}{1+or^c}}}
#' \deqn{r\_se = logor\_se * ((\pi * c * or^c) * \frac{\sin(\frac{\pi}{1+or^c})}{1+or^c})^2}
#' \deqn{or\_ci\_lo = exp(log(or) - qnorm(.975)*logor\_se)}
#' \deqn{or\_ci\_up = exp(log(or) + qnorm(.975)*logor\_se)}
#' \deqn{r\_ci\_lo = cos(\frac{\pi}{1 + or\_ci\_lo^c})}
#' \deqn{r\_ci\_up = cos(\frac{\pi}{1 + or\_ci\_up^c})}
#' \deqn{z = atanh(r)}
#' \deqn{z\_se = \sqrt{\frac{r\_se^2}{(1 - r^2)^2}} }
#' \deqn{z\_ci\_lo = atanh(r\_lo)}
#' \deqn{z\_ci\_up = atanh(r\_up)}
#'
#' **D.** Last, the approach described in Cooper et al. (2019) can be used (\code{or_to_cor = "lipsey_cooper"}).
#' This argument requires (or + logor_se + n_exp + n_nexp) to generate a R/Z.
#' As shown above, the function starts to estimate a SMD from the OR.
#' Then, as described in \code{\link{es_from_cohen_d}}, it converts this Cohen's d value into a correlation
#' coefficient using the \code{"lipsey_cooper"} formulas.
#'
#'
#' @export es_from_or_se
#'
#' @md
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab RR + NNT\cr
#'  \code{} \tab D + G + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 2. Odds Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Bonett, Douglas G. and Robert M. Price. (2005). Inferential Methods for the Tetrachoric Correlation Coefficient. Journal of Educational and Behavioral Statistics 30:213-25.
#'
#' Bonett, D. G., & Price, R. M. (2007). Statistical inference for generalized Yule coefficients in 2× 2 contingency tables. Sociological methods & research, 35(3), 429-446.
#'
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Di Pietrantonj C. (2006). Four-fold table cell frequencies imputation in meta analysis. Statistics in medicine, 25(13), 2299–2322. https://doi.org/10.1002/sim.2287
#'
#' Digby, Peter G. N. (1983). Approximating the Tetrachoric Correlation Coefficient. Biometrics 39:753-7.
#'
#' Gosling, C. J., Solanes, A., Fusar-Poli, P., & Radua, J. (2023). metaumbrella: the first comprehensive suite to perform data analysis in umbrella reviews with stratification of the evidence. BMJ mental health, 26(1), e300534. https://doi.org/10.1136/bmjment-2022-300534
#'
#' Grant R. L. (2014). Converting an odds ratio to a range of plausible relative risks for better communication of research findings. BMJ (Clinical research ed.), 348, f7450. https://doi.org/10.1136/bmj.f7450
#'
#' Pearson, K. (1900). Mathematical Contributions to the Theory of Evolution. VII: On the Correlation of Characters Not Quantitatively Measurable. Philosophical Transactions of the Royal Statistical Society of London, Series A 19:1-47
#'
#' Veroniki, A. A., Pavlides, M., Patsopoulos, N. A., & Salanti, G. (2013). Reconstructing 2x2 contingency tables from odds ratios using the Di Pietrantonj method: difficulties, constraints and impact in meta-analysis results. Research synthesis methods, 4(1), 78–94. https://doi.org/10.1002/jrsm.1061
#'
#' @examples
#' es_from_or_se(or = 2.12, logor_se = 0.242, n_exp = 120, n_nexp = 44)
es_from_or_se <- function(or, logor, logor_se, baseline_risk,
                          small_margin_prop,
                          n_exp, n_nexp, n_cases, n_controls, n_sample,
                          or_to_rr = "metaumbrella_cases",
                          or_to_cor = "pearson", reverse_or) {

  if (!all(or_to_cor %in% c("pearson", "digby", "bonett", "lipsey_cooper"))) {
    stop(paste0("'",
                unique(or_to_cor[!or_to_cor %in% c("pearson", "digby", "bonett", "lipsey_cooper")]),
                "' not in tolerated values for the 'or_to_cor' argument.
                Possible inputs are: 'pearson', 'digby', 'bonett', 'lipsey_cooper'"))

  } else if (!all(or_to_rr %in% c("metaumbrella_cases", "metaumbrella_exp", "transpose", "grant", "dipietrantonj"))) {
    stop(paste0("'",
                unique(or_to_rr[!or_to_rr %in% c("metaumbrella_cases", "metaumbrella_exp", "transpose", "grant", "dipietrantonj")]),
                "' not in tolerated values for the 'or_to_rr' argument.
                Possible inputs are: 'metaumbrella_cases', 'metaumbrella_exp', 'transpose', 'grant', 'dipietrantonj'"))
  }


  # or = dat$or; logor = dat$logor; logor_se = dat$logor_se;
  # baseline_risk = dat$baseline_risk;
  # small_margin_prop = dat$small_margin_prop;
  # n_exp = dat$n_exp; n_nexp = dat$n_nexp; n_cases = dat$n_cases; n_controls = dat$n_controls;
  # or_to_rr = "metaumbrella_cases";
  # or_to_cor = "cooper_delta";
  # reverse_or <- rep(FALSE, length(or))
  # or_to_rr = "dipietrantonj"
  if (missing(or)) {
    or <- rep(NA_real_, length(logor))
  }
  if (missing(logor)) {
    logor <- rep(NA_real_, length(or))
  }
  if (missing(baseline_risk)) {
    baseline_risk <- rep(NA_real_, length(or))
  }
  if (missing(small_margin_prop)) {
    small_margin_prop <- rep(NA_real_, length(or))
  }
  if (missing(n_exp)) {
    n_exp <- rep(NA_real_, length(or))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA_real_, length(or))
  }
  if (missing(n_cases)) {
    n_cases <- rep(NA_real_, length(or))
  }
  if (missing(n_controls)) {
    n_controls <- rep(NA_real_, length(or))
  }
  if (missing(reverse_or)) {
    reverse_or <- rep(FALSE, length(or))
  }
  if (missing(n_sample)) {
    n_sample <- rep(NA_real_, length(or))
  }
  reverse_or[is.na(reverse_or)] <- FALSE
  if (length(reverse_or) == 1) reverse_or = c(rep(reverse_or, length(or)))
  if (length(reverse_or) != length(or)) stop("The length of the 'reverse_or' argument is incorrectly specified.")

  tryCatch({
    .validate_positive(logor_se, baseline_risk, small_margin_prop, n_exp, n_nexp,
                       n_cases, n_controls, n_sample,
                       error_message = paste0("The number of people exposed/non-exposed, cases/controls, total sample size, ",
                                              "baseline risk, standard error of the logOR",
                                              "should be >0."),
                       func = "es_from_or_se")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  or <- ifelse(is.na(or) & !is.na(logor), exp(logor), or)

  logOR <- suppressWarnings(log(or))

  d <- logOR * sqrt(3) / pi
  d_se <- sqrt(logor_se^2 * 3 / (pi^2))

  n_sample <- ifelse(is.na(n_sample),
    ifelse(!is.na(n_cases) & !is.na(n_controls),
           n_cases + n_controls,
           n_exp + n_nexp),
    n_sample
  )

  es <- .es_from_d(
    d = d, d_se = d_se, n_exp = n_exp, n_nexp = n_nexp,
    n_sample = n_sample, reverse = reverse_or,
    smd_to_cor = rep("lipsey_cooper", length(d))
  )

  row_miss = which(is.na(d_se))
  es$d[row_miss] <- es$d_se[row_miss] <-
    es$d_ci_lo[row_miss] <- es$d_ci_up[row_miss] <-
    es$g[row_miss] <- es$g_se[row_miss] <-
    es$g_ci_lo[row_miss] <- es$g_ci_up[row_miss] <- NA
  es$r[row_miss] <- es$r_se[row_miss] <-
    es$r_ci_lo[row_miss] <- es$r_ci_up[row_miss] <-
    es$z[row_miss] <- es$z_se[row_miss] <-
    es$z_ci_lo[row_miss] <- es$z_ci_up[row_miss] <- NA


  # OR -------
  es$logor <- ifelse(reverse_or, -logOR, logOR)
  es$logor_se <- logor_se
  es$logor_ci_lo <- es$logor - es$logor_se * qnorm(.975)
  es$logor_ci_up <- es$logor + es$logor_se * qnorm(.975)
  or_ci_lo <- exp(log(or) - logor_se * qnorm(.975))
  or_ci_up <- exp(log(or) + logor_se * qnorm(.975))

  # RR -------
  es$logrr <- es$logrr_se <- es$logrr_ci_lo <- es$logrr_ci_up <- NA

  dat_rr <- data.frame(
    or = or, logor_se = logor_se, or_ci_lo = or_ci_lo, or_ci_up = or_ci_up,
    n_cases = n_cases, n_controls = n_controls, n_exp = n_exp, n_nexp = n_nexp,
    baseline_risk = baseline_risk, or_to_rr = or_to_rr
  )

  # print(paste0(or, ", ", or_ci_lo, ", ", or_ci_up, ", ",n_exp, ", ",n_nexp))
  nn_miss <- with(dat_rr, which(
      (or_to_rr == "grant" & !is.na(or) & !is.na(baseline_risk)) |
      (or_to_rr == "metaumbrella_cases" & !is.na(or) & !is.na(logor_se) &
         !is.na(n_cases) & !is.na(n_controls)) |
      (or_to_rr == "metaumbrella_exp" & !is.na(or) & !is.na(logor_se) &
         !is.na(n_exp) & !is.na(n_nexp)) |
      (or_to_rr == "transpose" & !is.na(or) & !is.na(logor_se)) |
      (or_to_rr == "dipietrantonj" & !is.na(or) & !is.na(or_ci_lo) & !is.na(or_ci_up) &
         !is.na(n_exp) & !is.na(n_nexp))
  ))


  if (length(nn_miss) != 0) {
    res_rr <- data.frame(t(mapply(.or_to_rr,
      or = dat_rr$or[nn_miss],
      logor_se = dat_rr$logor_se[nn_miss],
      or_ci_lo = dat_rr$or_ci_lo[nn_miss],
      or_ci_up = dat_rr$or_ci_up[nn_miss],
      n_cases = dat_rr$n_cases[nn_miss],
      n_controls = dat_rr$n_controls[nn_miss],
      n_exp = dat_rr$n_exp[nn_miss],
      n_nexp = dat_rr$n_nexp[nn_miss],
      baseline_risk = dat_rr$baseline_risk[nn_miss],
      or_to_rr = dat_rr$or_to_rr[nn_miss]
    )))
    res_rr$logrr = unlist(res_rr$logrr)
    res_rr$logrr_se = unlist(res_rr$logrr_se)
    res_rr$logrr_ci_lo = unlist(res_rr$logrr_ci_lo)
    res_rr$logrr_ci_up = unlist(res_rr$logrr_ci_up)
    es$logrr[nn_miss] <- ifelse(reverse_or[nn_miss], -res_rr[, 1], res_rr[, 1])
    es$logrr_se[nn_miss] <- res_rr[, 2]
    es$logrr_ci_lo[nn_miss] <- ifelse(reverse_or[nn_miss], res_rr[, 4], res_rr[, 3])
    es$logrr_ci_up[nn_miss] <- ifelse(reverse_or[nn_miss], res_rr[, 3], res_rr[, 4])
  }

  # COR -------
  dat_cor <- data.frame(
    or = or, logor_se = logor_se,
    n_cases = n_cases, n_exp = n_exp,
    n_sample = n_sample,
    small_margin_prop = small_margin_prop,
    or_to_cor = or_to_cor
  )

  nn_miss <- with(dat_cor, which(
    (or_to_cor == "bonett" & !is.na(or) & !is.na(logor_se) &
       !is.na(small_margin_prop) & !is.na(n_sample) &
      !is.na(n_exp) & !is.na(n_cases)) |
      (or_to_cor == "pearson" & !is.na(or) & !is.na(logor_se)) |
      (or_to_cor == "digby" & !is.na(or) & !is.na(logor_se))
  ))

  if (length(nn_miss) != 0) {
    no_cooper = which(or_to_cor != "lipsey_cooper")
    es$r[no_cooper] <- es$r_se[no_cooper] <- es$r_ci_lo[no_cooper] <- es$r_ci_up[no_cooper] <-
      es$z[no_cooper] <- es$z_se[no_cooper] <- es$z_ci_lo[no_cooper] <- es$z_ci_up[no_cooper] <- NA

    res_cor <- t(mapply(.or_to_cor,
      or = dat_cor$or[nn_miss],
      logor_se = dat_cor$logor_se[nn_miss],
      n_cases = dat_cor$n_cases[nn_miss],
      n_exp = dat_cor$n_exp[nn_miss],
      n_sample = dat_cor$n_sample[nn_miss],
      small_margin_prop = dat_cor$small_margin_prop[nn_miss],
      or_to_cor = dat_cor$or_to_cor[nn_miss]
    ))

    es$r[nn_miss] <- ifelse(reverse_or[nn_miss], -res_cor[, 1], res_cor[, 1]) # res_cor[, 1]
    es$r_se[nn_miss] <- res_cor[, 2]
    es$r_ci_lo[nn_miss] <- ifelse(reverse_or[nn_miss], res_cor[, 4], res_cor[, 3]) # res_cor[, 3]
    es$r_ci_up[nn_miss] <- ifelse(reverse_or[nn_miss], res_cor[, 3], res_cor[, 4])
    es$z[nn_miss] <- ifelse(reverse_or[nn_miss], -res_cor[, 5], res_cor[, 5]) # res_cor[, 5]
    es$z_se[nn_miss] <- res_cor[, 6]
    es$z_ci_lo[nn_miss] <- ifelse(reverse_or[nn_miss], res_cor[, 8], res_cor[, 7]) # res_cor[, 7]
    es$z_ci_up[nn_miss] <- ifelse(reverse_or[nn_miss], res_cor[, 7], res_cor[, 8]) # res_cor[, 8]
  }

  es$nnt <- (1 - baseline_risk * (1 - or)) /
    ((1 - baseline_risk) * (baseline_risk * (1 - or)))

  es$nnt <- ifelse(reverse_or, -es$nnt, es$nnt)

  es$info_used <- "or_se"
  return(es)
}

#' Convert an odds ratio value to several effect size measures
#'
#' @param or odds ratio value
#' @param logor log odds ratio value
#' @param n_cases number of cases/events
#' @param n_controls number of controls/no-event
#' @param n_exp number of participants in the exposed group
#' @param n_nexp number of participants in the non-exposed group
#' @param n_sample total number of participants in the sample
#' @param baseline_risk proportion of cases in the non-exposed group
#' @param small_margin_prop smallest margin proportion of the underlying 2x2 table
#' @param reverse_or a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param or_to_rr formula used to convert the \code{or} value into a risk ratio (see details).
#' @param or_to_cor formula used to convert the \code{or} value into a correlation coefficient (see details).
#'
#' @details
#' This function computes the standard error of the log odds ratio.
#' Risk ratio (RR), Cohen's d (D), Hedges' g (G) and correlation coefficients (R/Z),
#' are converted from the odds ratio value.
#'
#' **Estimation of the standard error of the log OR.**
#' This function generates the standard error of an odds ratio (OR) based on the OR value and the number of cases and controls.
#' More precisely, this function simulates all combinations of the possible number of cases and controls in the exposed and non-exposed groups
#' compatible with the reported OR value and with the overall number of cases and controls. Then, our function assumes that the variance of the OR
#' is equal to the mean of the standard error of all possible situations. This estimation thus necessarily comes with some imprecision and should not
#' be used before having requested the value (or raw data) to authors of the original report.
#'
#' **Conversion of other effect size measures.**
#' Calculations of \code{\link{es_from_or_se}()} are then applied to
#' estimate the other effect size measures
#'
#' @references
#' Gosling, C. J., Solanes, A., Fusar-Poli, P., & Radua, J. (2023). metaumbrella: the first comprehensive suite to perform data analysis in umbrella reviews with stratification of the evidence. BMJ mental health, 26(1), e300534. https://doi.org/10.1136/bmjment-2022-300534
#'
#' @export es_from_or
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab N/A\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + RR + NNT\cr
#'  \code{} \tab D + G + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 2. Odds Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_or_guess <- es_from_or(or = 0.5, n_cases = 210, n_controls = 220)
#' es_or <- es_from_or_se(or = 0.5, logor_se = 0.4, n_cases = 210, n_controls = 220)
#' round(es_or_guess$logor_se, 0.10) == round(es_or$logor_se, 0.10)
es_from_or <- function(or, logor, n_cases, n_controls, n_sample,
                       small_margin_prop, baseline_risk,
                       n_exp, n_nexp,
                       or_to_cor = "bonett", or_to_rr = "metaumbrella_cases",
                       reverse_or) {
  if (missing(or)) {
    or <- rep(NA_real_, length(logor))
  }
  if (missing(logor)) {
    logor <- rep(NA_real_, length(or))
  }
  if (missing(baseline_risk)) {
    baseline_risk <- rep(NA_real_, length(or))
  }
  if (missing(small_margin_prop)) {
    small_margin_prop <- rep(NA_real_, length(or))
  }
  if (missing(n_exp)) {
    n_exp <- rep(NA_real_, length(or))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA_real_, length(or))
  }
  if (missing(n_cases)) {
    n_cases <- rep(NA_real_, length(or))
  }
  if (missing(n_controls)) {
    n_controls <- rep(NA_real_, length(or))
  }
  if (missing(reverse_or)) {
    reverse_or <- rep(FALSE, length(or))
  }
  if (missing(n_sample)) {
    n_sample <- rep(NA_real_, length(or))
  }
  reverse_or[is.na(reverse_or)] <- FALSE

  tryCatch({
    .validate_positive(baseline_risk, small_margin_prop, n_exp, n_nexp,
                       n_cases, n_controls, n_sample,
                       error_message = paste0("The number of people exposed/non-exposed, cases/controls, total sample size, ",
                                              "baseline risk,",
                                              "should be >0."),
                       func = "es_from_or")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  or <- ifelse(is.na(or) & !is.na(logor), exp(logor), or)

  dat_or_se <- data.frame(or, n_cases, n_controls)
  log_or_se <- do.call(rbind, apply(dat_or_se, 1, .se_from_or))

  es <- es_from_or_se(
    or = or, logor_se = log_or_se$se,
    n_cases = n_cases, n_controls = n_controls,
    small_margin_prop = small_margin_prop,
    baseline_risk = baseline_risk, n_sample=n_sample,
    n_exp = n_exp, n_nexp = n_nexp,
    or_to_cor = or_to_cor, or_to_rr = or_to_rr,
    reverse_or = reverse_or
  )

  es$info_used <- "or"
  return(es)
}

#' Convert an odds ratio value and its 95% confidence interval to several effect size measures
#'
#' @param or odds ratio value
#' @param logor log odds ratio value
#' @param or_ci_lo lower bound of the 95% CI around the odds ratio value
#' @param or_ci_up upper bound of the 95% CI around the odds ratio value
#' @param logor_ci_lo lower bound of the 95% CI around the log odds ratio value
#' @param logor_ci_up upper bound of the 95% CI around the log odds ratio value
#' @param n_cases number of cases/events
#' @param n_controls number of controls/no-event
#' @param n_exp number of participants in the exposed group (only required for the \code{or_to_rr = "grant"}, and \code{or_to_rr = "metaumbrella_exp"} arguments)
#' @param n_nexp number of participants in the non-exposed group (only required for the \code{or_to_rr = "grant"}, and \code{or_to_rr = "metaumbrella_exp"} arguments)
#' @param n_sample total number of participants in the sample
#' @param baseline_risk proportion of cases in the non-exposed group (only required for the \code{or_to_rr = "grant"} argument).
#' @param small_margin_prop smallest margin proportion of the underlying 2x2 table
#' @param reverse_or a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param or_to_cor formula used to convert the \code{or} value into a correlation coefficient (see details).
#' @param or_to_rr formula used to convert the \code{or} value into a risk ratio (see details).
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#'
#' @details
#' This function computes the standard error of the (log) odds ratio
#' into a standard error (Section 6.5.2.2 in the Cochrane Handbook).
#' \deqn{logor\_se = \frac{\log{or\_ci\_up} - \log{or\_ci\_lo}}{2 * qnorm(.975)}}
#'
#' Then, calculations of \code{\link{es_from_or_se}} are applied.
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @export es_from_or_ci
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab RR + NNT\cr
#'  \code{} \tab D + G + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 2. Odds Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_or <- es_from_or_ci(
#'   or = 1, or_ci_lo = 0.5, or_ci_up = 2,
#'   n_cases = 42, n_controls = 38, baseline_risk = 0.08,
#'   or_to_rr = "grant"
#' )
es_from_or_ci <- function(or, or_ci_lo, or_ci_up, logor, logor_ci_lo, logor_ci_up,
                          baseline_risk, small_margin_prop, n_exp, n_nexp,
                          n_cases, n_controls, n_sample, max_asymmetry = 10,
                          or_to_cor = "bonett", or_to_rr = "metaumbrella_cases",
                          reverse_or) {
  if (missing(or)) {
    or <- rep(NA_real_, length(logor))
  }
  if (missing(logor)) {
    logor <- rep(NA_real_, length(or))
  }
  if (missing(baseline_risk)) {
    baseline_risk <- rep(NA_real_, length(or))
  }
  if (missing(small_margin_prop)) {
    small_margin_prop <- rep(NA_real_, length(or))
  }
  if (missing(n_exp)) {
    n_exp <- rep(NA_real_, length(or))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA_real_, length(or))
  }
  if (missing(n_cases)) {
    n_cases <- rep(NA_real_, length(or))
  }
  if (missing(n_controls)) {
    n_controls <- rep(NA_real_, length(or))
  }
  if (missing(or_ci_lo)) {
    or_ci_lo <- rep(NA_real_, length(or))
  }
  if (missing(logor_ci_lo)) {
    logor_ci_lo <- rep(NA_real_, length(or))
  }
  if (missing(or_ci_up)) {
    or_ci_up <- rep(NA_real_, length(or))
  }
  if (missing(logor_ci_up)) {
    logor_ci_up <- rep(NA_real_, length(or))
  }
  if (missing(reverse_or)) {
    reverse_or <- rep(FALSE, length(or))
  }
  if (missing(n_sample)) {
    n_sample <- rep(FALSE, length(or))
  }
  reverse_or[is.na(reverse_or)] <- FALSE

  or <- ifelse(is.na(or) & !is.na(logor), exp(logor), or)

  tryCatch({
    .validate_positive(baseline_risk, small_margin_prop, n_exp, n_nexp,
                       n_cases, n_controls, n_sample,
                       error_message = paste0("The number of people exposed/non-exposed, cases/controls, total sample size, ",
                                              "baseline risk, standard error of the logOR",
                                              "should be >0."),
                       func = "es_from_or_se")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  tryCatch({
    .validate_ci_symmetry(logor, logor_ci_lo, logor_ci_up,
                          func = "es_from_or_ci",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })
  tryCatch({
    .validate_ci_symmetry(log(or), log(or_ci_lo), log(or_ci_up),
                          func = "es_from_or_ci",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })


  logor_ci_lo <- ifelse(is.na(logor_ci_lo) & !is.na(or_ci_lo), log(or_ci_lo), logor_ci_lo)
  logor_ci_up <- ifelse(is.na(logor_ci_up) & !is.na(or_ci_up), log(or_ci_up), logor_ci_up)
  logor_se <- (logor_ci_up - logor_ci_lo) / (2 * qnorm(.975))

  es <- es_from_or_se(
    or = or, logor_se = logor_se,
    n_cases = n_cases, n_controls = n_controls,
    baseline_risk = baseline_risk, small_margin_prop = small_margin_prop,
    n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample,
    or_to_cor = or_to_cor, or_to_rr = or_to_rr,
    reverse_or = reverse_or
  )

  es$info_used <- "or_ci"

  return(es)
}

#' Convert an odds ratio value and its standard error to several effect size measures
#'
#' @param or odds ratio value
#' @param logor log odds ratio value
#' @param or_pval p-value of the (log) odds ratio
#' @param n_cases number of cases/events
#' @param n_controls number of controls/no-event
#' @param n_exp number of participants in the exposed group (only required for the \code{or_to_rr = "grant"}, and \code{or_to_rr = "metaumbrella_exp"} arguments)
#' @param n_nexp number of participants in the non-exposed group (only required for the \code{or_to_rr = "grant"}, and \code{or_to_rr = "metaumbrella_exp"} arguments)
#' @param n_sample total number of participants in the sample
#' @param baseline_risk proportion of cases in the non-exposed group (only required for the \code{or_to_rr = "grant"} argument).
#' @param small_margin_prop smallest margin proportion of the underlying 2x2 table
#' @param reverse_or_pval a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param or_to_cor formula used to convert the \code{or} value into a correlation coefficient (see details).
#' @param or_to_rr formula used to convert the \code{or} value into a risk ratio (see details).
#'
#' @details
#' This function computes the standard error of the (log) odds ratio into
#' from a p-value (Section 6.3.2 in the Cochrane Handbook).
#' \deqn{logor\_z = qnorm(or_pval/2, lower.tail=FALSE)}
#' \deqn{logor\_se = |\frac{\log(or)}{logor\_z}|}
#'
#' Then, calculations of \code{\link{es_from_or_se}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab RR + NNT\cr
#'  \code{} \tab D + G + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 2. Odds Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_or_pval
#'
#' @references
#' Higgins, J. P., Thomas, J., Chandler, J., Cumpston, M., Li, T., Page, M. J., & Welch, V. A. (Eds.). (2019). Cochrane handbook for systematic reviews of interventions. John Wiley & Sons.
#'
#' @md
#'
#' @examples
#' es_or <- es_from_or_pval(
#'   or = 3.51, or_pval = 0.001,
#'   n_cases = 12, n_controls = 68
#' )
es_from_or_pval <- function(or, logor, or_pval, baseline_risk, small_margin_prop,
                            n_exp, n_nexp, n_cases, n_controls, n_sample,
                            or_to_rr = "metaumbrella_cases",
                            or_to_cor = "bonett", reverse_or_pval) {
  if (missing(or)) {
    or <- rep(NA_real_, length(logor))
  }
  if (missing(logor)) {
    logor <- rep(NA_real_, length(or))
  }
  if (missing(baseline_risk)) {
    baseline_risk <- rep(NA_real_, length(or))
  }
  if (missing(small_margin_prop)) {
    small_margin_prop <- rep(NA_real_, length(or))
  }
  if (missing(n_exp)) {
    n_exp <- rep(NA_real_, length(or))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA_real_, length(or))
  }
  if (missing(n_cases)) {
    n_cases <- rep(NA_real_, length(or))
  }
  if (missing(n_controls)) {
    n_controls <- rep(NA_real_, length(or))
  }
  if (missing(reverse_or_pval)) {
    reverse_or_pval <- rep(FALSE, length(or))
  }
  if (missing(n_sample)) {
    n_sample <- rep(FALSE, length(or))
  }
  reverse_or_pval[is.na(reverse_or_pval)] <- FALSE

  tryCatch({
    .validate_positive(es_from_or_pval, baseline_risk,
                       small_margin_prop, n_exp, n_nexp, n_cases,
                       n_controls, n_sample,
                       error_message = paste0("The number of people exposed/non-exposed, cases/controls, total sample size, ",
                                              "baseline risk, p-value of the OR",
                                              "should be >0."),
                       func = "es_from_or_pval")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  or <- ifelse(is.na(or) & !is.na(logor), exp(logor), or)
  logOR <- suppressWarnings(log(or))

  z_or <- qnorm(or_pval / 2, lower.tail = FALSE)
  logor_se <- abs(logOR / z_or)

  es <- es_from_or_se(
    or = or, logor_se = logor_se,
    baseline_risk = baseline_risk, small_margin_prop = small_margin_prop,
    n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample,
    n_cases = n_cases, n_controls = n_controls,
    or_to_cor = or_to_cor, or_to_rr = or_to_rr,
    reverse_or = reverse_or_pval
  )

  es$info_used <- "or_pval"

  return(es)
}
