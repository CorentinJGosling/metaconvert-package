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
#' @param small_margin_prop smallest margin proportion of cases/events in the underlying 2x2 table (a proportion in (0, 0.5])
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
#' **A.** First, the approach identified by the \code{or_to_rr = "grant"} argument value
#' can be used. The estimator itself is due to Zhang and Yu (1998); Grant (2014) restates
#' the same expression and is cited here for the communication framing (translating one OR
#' into a range of plausible RRs), not as the source of the formula. The argument value
#' \code{"grant"} is retained for backward compatibility.
#' Neither paper gives a variance, standard error or confidence interval for the converted
#' value. To derive the variance, we used this formula to convert the bounds of the 95% CI,
#' which were then used to obtain the variance.
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
#' **C.** Third, it is possible to transpose the OR to a RR (\code{or_to_rr = "transpose"}).
#' This argument requires (or + logor_se) to generate a RR.
#' It is known that OR and RR are similar when the baseline risk is small.
#' Therefore, users can request to simply transpose the OR value & standard error into a RR value & standard error.
#' \deqn{rr = or}
#' \deqn{logrr\_se = logor\_se}
#'
#' **D.** Fourth, it is possible to recreate the 2x2 table using the dipietrantonj's formulas (\code{or_to_rr = "dipietrantonj"}).
#' This argument requires (or + logor_se + n_exp + n_nexp) to generate a RR. The 95% CI of
#' the OR, from which the 2x2 table is reconstructed, is obtained internally from \code{or}
#' and \code{logor_se}. Information on this approach can be retrieved in
#' Di Pietrantonj (2006).
#'
#' **To estimate the NNT, the formulas used are :**
#' \deqn{treatment\_risk = \frac{or \times br}{1 - br + or \times br}}
#' \deqn{rd = br - treatment\_risk}
#' \deqn{nnt = \frac{1}{rd} = \frac{1 - br \times (1 - or)}{br \times (1 - br) \times (1 - or)}}
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
#' \deqn{r\_se = logor\_se * (\pi * c * or^c) * \frac{\sin(\pi / (1+or^c))}{(1+or^c)^2}}
#' \deqn{or\_ci\_lo = exp(log(or) - qnorm(.975)*logor\_se)}
#' \deqn{or\_ci\_up = exp(log(or) + qnorm(.975)*logor\_se)}
#' \deqn{r\_ci\_lo = cos(\frac{\pi}{1 + or\_ci\_lo^c})}
#' \deqn{r\_ci\_up = cos(\frac{\pi}{1 + or\_ci\_up^c})}
#' \deqn{z = atanh(r)}
#' \deqn{z\_se = \sqrt{\frac{r\_se^2}{(1 - r^2)^2}}}
#' \deqn{z\_ci\_lo = atanh(r\_ci\_lo)}
#' \deqn{z\_ci\_up = atanh(r\_ci\_up)}
#'
#' **B.** Second, the approach described in Digby (1983) can be used (\code{or_to_cor = "digby"}).
#' This argument requires (or + logor_se) to generate a R/Z.
#' It converts the OR value and its standard error to a tetrachoric correlation.
#' Note that the formula assumes that each cell of the 2x2 used to estimate the OR has been added 1/2 before estimating the OR value and its standard error.
#' If it is not the case, formulas can produce slightly less accurate results.
#'
#' \deqn{c = \frac{3}{4}}
#' \deqn{r = \frac{or^c - 1}{or^c + 1}}
#' \deqn{r\_se = \sqrt{\frac{c^2}{4} * (1 - r^2)^2 * logor\_se^2}}
#' \deqn{z = atanh(r)}
#' \deqn{z\_se = \sqrt{\frac{r\_se^2}{(1 - r^2)^2}}}
#' \deqn{z\_ci\_lo = z - qnorm(.975)*\sqrt{\frac{c^2}{4} * logor\_se^2}}
#' \deqn{z\_ci\_up = z + qnorm(.975)*\sqrt{\frac{c^2}{4} * logor\_se^2}}
#' \deqn{r\_ci\_lo = tanh(z\_ci\_lo)}
#' \deqn{r\_ci\_up = tanh(z\_ci\_up)}
#'
#' **C.** Third, the approach described in Bonett (2005) can be used (\code{or_to_cor = "bonett"}).
#' This argument requires (or + logor_se + n_sample + the 2x2 margins) to generate a R/Z.
#' Note that the formula assumes that each cell of the 2x2 used to estimate the OR has been added 1/2 before estimating the OR value and its standard error.
#' If it is not the case, formulas can produce slightly less accurate results.
#'
#' \code{small_margin_prop} is Bonett and Price's \eqn{p_{min}}, the smallest of the four
#' marginal proportions of the underlying 2x2 table. When it is left blank it is derived
#' as \eqn{min(n\_exp, n\_nexp, n\_cases, n\_controls) / n\_sample}, so the method runs on
#' the margins alone. The two margin pairs (\code{n_exp}/\code{n_nexp} and
#' \code{n_cases}/\code{n_controls}) each sum to \code{n_sample}, so it is enough to supply
#' \code{n_sample} plus \emph{one member of each pair}: any of the four combinations
#' (\code{n_exp} or \code{n_nexp}) x (\code{n_cases} or \code{n_controls}) works, and the
#' missing margins are back-filled. A value supplied by the user always takes precedence;
#' it must be a proportion in (0, 0.5], and a value outside that range is not validated
#' and can drive \eqn{c} negative, which flips the sign of \code{r}.
#' The derivation is skipped when the implied margins are not all strictly positive or do
#' not sum to \code{n_sample}. A row that cannot use the Bonett conversion keeps the
#' \code{"lipsey_cooper"} R/Z computed from the Cohen's d, and a message reports the
#' substitution. This decision is taken row by row and does not depend on the other rows
#' of the call.
#'
#' \deqn{c = \frac{1 - \frac{|\frac{n\_exp}{n\_sample} - \frac{n\_cases}{n\_sample}|}{5} - (0.5 - small\_margin\_prop)^2}{2}}
#' \deqn{r = \cos{\frac{\pi}{1+or^c}}}
#' \deqn{r\_se = logor\_se * (\pi * c * or^c) * \frac{\sin(\frac{\pi}{1+or^c})}{(1+or^c)^2}}
#' \deqn{or\_ci\_lo = exp(log(or) - qnorm(.975)*logor\_se)}
#' \deqn{or\_ci\_up = exp(log(or) + qnorm(.975)*logor\_se)}
#' \deqn{r\_ci\_lo = cos(\frac{\pi}{1 + or\_ci\_lo^c})}
#' \deqn{r\_ci\_up = cos(\frac{\pi}{1 + or\_ci\_up^c})}
#' \deqn{z = atanh(r)}
#' \deqn{z\_se = \sqrt{\frac{r\_se^2}{(1 - r^2)^2}} }
#' \deqn{z\_ci\_lo = atanh(r\_ci\_lo)}
#' \deqn{z\_ci\_up = atanh(r\_ci\_up)}
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
#'  \code{converted effect size measure} \tab RR + NNT + RD\cr
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
#' Bonett, D. G., & Price, R. M. (2007). Statistical inference for generalized Yule coefficients in 2* 2 contingency tables. Sociological methods & research, 35(3), 429-446.
#'
#' Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Di Pietrantonj C. (2006). Four-fold table cell frequencies imputation in meta analysis. Statistics in medicine, 25(13), 2299-2322. https://doi.org/10.1002/sim.2287
#'
#' Digby, Peter G. N. (1983). Approximating the Tetrachoric Correlation Coefficient. Biometrics 39:753-7.
#'
#' Gosling, C. J., Solanes, A., Fusar-Poli, P., & Radua, J. (2023). metaumbrella: the first comprehensive suite to perform data analysis in umbrella reviews with stratification of the evidence. BMJ mental health, 26(1), e300534. https://doi.org/10.1136/bmjment-2022-300534
#'
#' Grant R. L. (2014). Converting an odds ratio to a range of plausible relative risks for better communication of research findings. BMJ (Clinical research ed.), 348, f7450. https://doi.org/10.1136/bmj.f7450
#'
#' Pearson, K. (1900). Mathematical Contributions to the Theory of Evolution. VII: On the Correlation of Characters Not Quantitatively Measurable. Philosophical Transactions of the Royal Statistical Society of London, Series A 19:1-47
#'
#' Veroniki, A. A., Pavlides, M., Patsopoulos, N. A., & Salanti, G. (2013). Reconstructing 2x2 contingency tables from odds ratios using the Di Pietrantonj method: difficulties, constraints and impact in meta-analysis results. Research synthesis methods, 4(1), 78-94. https://doi.org/10.1002/jrsm.1061
#'
#' Zhang, J., & Yu, K. F. (1998). What's the relative risk? A method of correcting the odds ratio in cohort studies of common outcomes. JAMA, 280(19), 1690-1691. https://doi.org/10.1001/jama.280.19.1690
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

  or <- ifelse(is.na(or) & !is.na(logor), exp(logor), or)

  # A non-positive standard error would be carried unchanged into logor_se and
  # from there into a negative sampling variance. See R/internal_guards.R.
  logor_se <- .positive_or_na(logor_se)

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
    # Columns come out of t(mapply(...)) positionally as logrr / logrr_se /
    # logrr_ci_lo / logrr_ci_up. .mapply_col() forces each to a plain numeric
    # vector of the right length; the previous unlist()-based extraction was a
    # no-op on the numeric-matrix path (the data.frame columns are named X1..X4,
    # not logrr*) and silently recycled values on the list-matrix path.
    rr_es    <- .mapply_col(res_rr, 1)
    rr_se    <- .mapply_col(res_rr, 2)
    rr_ci_lo <- .mapply_col(res_rr, 3)
    rr_ci_up <- .mapply_col(res_rr, 4)
    es$logrr[nn_miss] <- ifelse(reverse_or[nn_miss], -rr_es, rr_es)
    es$logrr_se[nn_miss] <- rr_se
    # On reverse: negate AND swap the CI bounds (new_lo = -old_up, new_up = -old_lo);
    # swapping alone left a wrong-signed, inverted interval that did not bracket -logrr.
    es$logrr_ci_lo[nn_miss] <- ifelse(reverse_or[nn_miss], -rr_ci_up, rr_ci_lo)
    es$logrr_ci_up[nn_miss] <- ifelse(reverse_or[nn_miss], -rr_ci_lo, rr_ci_up)
  }

  # COR -------
  # Derive small_margin_prop when the user left it blank.
  #
  # Bonett & Price (2005, p. 216) define pmin as "the smallest marginal proportion" of
  # the 2x2 table, i.e. min(p1+, p2+, p+1, p+2). Those four margins are exactly
  # n_exp, n_nexp, n_cases and n_controls over n_sample. The two margin PAIRS
  # (n_exp/n_nexp and n_cases/n_controls) each sum to n_sample, so n_sample plus one
  # member of each pair determines the whole table -- there are four equivalent ways to
  # describe it, and all four must work.
  #
  # The back-fill below is therefore SYMMETRIC. An earlier version derived only
  # n_nexp from n_exp and n_controls from n_cases, so a user who entered the other
  # member of either pair fell out of the gate and silently kept the lipsey_cooper
  # value: for one table (n_exp 40, n_nexp 60, n_cases 30, n_controls 70, n_sample 100,
  # or 2.5, logor_se 0.2) only n_exp + n_cases gave the bonett r = 0.326978458193; the
  # other three gave 0.243943602135 / 0.243470956833 / 0.243943602135, bit-identical to
  # or_to_cor = "lipsey_cooper" on the same inputs.
  #
  # n_exp and n_cases are then the ones .or_to_cor() consumes (its c formula names them
  # directly), so the back-filled values, not the raw arguments, go into dat_cor.
  #
  # Without this, small_margin_prop was a user-entered column with no derivation
  # anywhere, so a blank cell -- the normal case -- dropped the row out of the gate and
  # silently left the lipsey_cooper r computed by .es_from_d() above in place. Since
  # convert_df() defaults to or_to_cor = "bonett", the default method usually did not
  # run: or = 2, logor_se = 0.2, n_exp = n_nexp = 50, n_cases = 40, n_controls = 60
  # returned r = 0.1876806337 (lipsey_cooper) where bonett gives 0.2586008221.
  # Worse, the fall-through was not even stable: because the block below NA's every
  # non-lipsey_cooper row before refilling only the gated ones, a blank-pmin row
  # returned NA rather than 0.1876806337 whenever some OTHER row in the same data
  # frame did carry small_margin_prop -- so a row's own ES depended on another row.
  #
  # Margins are taken UNCORRECTED (no +0.5 per cell). That reproduces all three of
  # Bonett & Price's worked examples as printed: Ex1 f=(203,186,167,374) -> .333/.048/
  # (.237,.424); Ex2 f=(4,6,1,89) -> .831/.108/(.488,.956); Ex3 f=(143,52,41,164) ->
  # .741/.0446/(.641,.817). The paper's own prose ("computed ... after 0.5 has been
  # added to each cell frequency") would give .835 on Ex2, so its printed numbers use
  # the uncorrected margins and this convention is the faithful one.
  #
  # A user-supplied small_margin_prop always wins. The derivation is skipped unless the
  # four implied margins are strictly positive and both margin pairs sum to n_sample:
  # a degenerate or self-contradictory table has no coherent pmin (and can send c
  # negative, flipping the sign of r), so such a row keeps the old fall-through.
  n_exp_marg <- ifelse(is.na(n_exp), n_sample - n_nexp, n_exp)
  n_cases_marg <- ifelse(is.na(n_cases), n_sample - n_controls, n_cases)
  n_nexp_marg <- ifelse(is.na(n_nexp), n_sample - n_exp_marg, n_nexp)
  n_controls_marg <- ifelse(is.na(n_controls), n_sample - n_cases_marg, n_controls)
  margins_ok <- !is.na(n_sample) &
    !is.na(n_exp_marg) & !is.na(n_cases_marg) &
    !is.na(n_nexp_marg) & !is.na(n_controls_marg) &
    n_exp_marg > 0 & n_nexp_marg > 0 & n_cases_marg > 0 & n_controls_marg > 0 &
    abs(n_exp_marg + n_nexp_marg - n_sample) <= 1e-8 * pmax(1, abs(n_sample)) &
    abs(n_cases_marg + n_controls_marg - n_sample) <= 1e-8 * pmax(1, abs(n_sample))
  small_margin_prop <- ifelse(
    is.na(small_margin_prop) & margins_ok,
    pmin(n_exp_marg, n_nexp_marg, n_cases_marg, n_controls_marg) / n_sample,
    small_margin_prop
  )

  # Only hand .or_to_cor() a back-filled margin when the implied table is coherent;
  # otherwise keep the raw value so an incoherent row still fails the gate below.
  dat_cor <- data.frame(
    or = or, logor_se = logor_se,
    n_cases = ifelse(margins_ok, n_cases_marg, n_cases),
    n_exp = ifelse(margins_ok, n_exp_marg, n_exp),
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

  # A row that requested a conversion other than "lipsey_cooper" but does not carry that
  # method's inputs keeps the "lipsey_cooper" R/Z computed earlier, and the substitution
  # is reported. Such a row used to be blanked to NA, which had two consequences. First,
  # activating the small_margin_prop derivation switched the blanking on for the whole
  # call, so margin-poor rows that the previous release estimated came back as NA and the
  # study was lost. Second, the decision was taken INSIDE the `if (length(nn_miss) != 0)`
  # block, i.e. only when some OTHER row qualified, so one and the same row returned a
  # number when called alone and NA when called beside an eligible peer.
  # Note this is read from the RECYCLED column, not the raw argument: with a scalar
  # or_to_cor, which(or_to_cor != "lipsey_cooper") has length 1 and names only row 1.
  # Restricted to rows that actually have an odds ratio to convert: a row with no OR has
  # no correlation by any method, so it has not "fallen back" to anything.
  fallback <- setdiff(
    which(dat_cor$or_to_cor != "lipsey_cooper" &
            !is.na(dat_cor$or) & !is.na(dat_cor$logor_se)),
    nn_miss
  )
  if (length(fallback) != 0) {
    asked <- unique(as.character(dat_cor$or_to_cor[fallback]))
    message("or_to_cor = '", paste(asked, collapse = "' / '"), "': row(s) ",
            paste(fallback, collapse = ", "), " do not carry the inputs this conversion ",
            "requires, so their R and Z were obtained with 'lipsey_cooper' instead. ",
            "For 'bonett', supply 'small_margin_prop', or 'n_sample' together with one of ",
            "('n_exp', 'n_nexp') and one of ('n_cases', 'n_controls').")
  }

  if (length(nn_miss) != 0) {
    res_cor <- .mapply_memo(.or_to_cor,
      or = dat_cor$or[nn_miss],
      logor_se = dat_cor$logor_se[nn_miss],
      n_cases = dat_cor$n_cases[nn_miss],
      n_exp = dat_cor$n_exp[nn_miss],
      n_sample = dat_cor$n_sample[nn_miss],
      small_margin_prop = dat_cor$small_margin_prop[nn_miss],
      or_to_cor = dat_cor$or_to_cor[nn_miss]
    )

    # On reverse: negate AND swap each CI (new_lo = -old_up, new_up = -old_lo). For r,
    # tanh is odd so negate-and-swap of the r bounds is correct even for the asymmetric,
    # z-back-transformed r interval. Swapping alone left inverted, wrong-signed bounds.
    cor_v <- lapply(1:8, function(j) .mapply_col(res_cor, j))
    es$r[nn_miss] <- ifelse(reverse_or[nn_miss], -cor_v[[1]], cor_v[[1]])
    es$r_se[nn_miss] <- cor_v[[2]]
    es$r_ci_lo[nn_miss] <- ifelse(reverse_or[nn_miss], -cor_v[[4]], cor_v[[3]])
    es$r_ci_up[nn_miss] <- ifelse(reverse_or[nn_miss], -cor_v[[3]], cor_v[[4]])
    es$z[nn_miss] <- ifelse(reverse_or[nn_miss], -cor_v[[5]], cor_v[[5]])
    es$z_se[nn_miss] <- cor_v[[6]]
    es$z_ci_lo[nn_miss] <- ifelse(reverse_or[nn_miss], -cor_v[[8]], cor_v[[7]])
    es$z_ci_up[nn_miss] <- ifelse(reverse_or[nn_miss], -cor_v[[7]], cor_v[[8]])
  }

  # Risk difference from OR + baseline_risk (Grant 2014)
  treatment_risk <- (or * baseline_risk) / (1 - baseline_risk + or * baseline_risk)
  rd <- baseline_risk - treatment_risk

  es$rd <- ifelse(reverse_or, -rd, rd)
  # delta method
  drd_dor <- baseline_risk * (1 - baseline_risk) / (1 - baseline_risk + or * baseline_risk)^2
  rd_se <- abs(drd_dor) * or * logor_se
  es$rd_se <- rd_se
  es$rd_ci_lo <- es$rd - qnorm(.975) * rd_se
  es$rd_ci_up <- es$rd + qnorm(.975) * rd_se

  es$nnt <- ifelse(rd == 0, NA, 1 / rd)
  es$nnt <- ifelse(reverse_or, -es$nnt, es$nnt)
  es$nnt_se <- ifelse(rd == 0, NA, rd_se / rd^2)
  rd_ci_lo_raw <- rd - qnorm(.975) * rd_se
  rd_ci_up_raw <- rd + qnorm(.975) * rd_se
  crosses_zero <- (rd_ci_lo_raw < 0 & rd_ci_up_raw > 0) | rd == 0
  es$nnt_ci_lo <- ifelse(crosses_zero, NA,
                          ifelse(reverse_or, -1 / rd_ci_lo_raw, 1 / rd_ci_up_raw))
  es$nnt_ci_up <- ifelse(crosses_zero, NA,
                          ifelse(reverse_or, -1 / rd_ci_up_raw, 1 / rd_ci_lo_raw))

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
#' @param baseline_risk proportion of cases in the non-exposed group (n_cases_nexp / n_nexp is used when missing)
#' @param small_margin_prop smallest margin proportion of the underlying 2x2 table (a proportion in (0, 0.5])
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
#' compatible with the reported OR value and with the overall number of cases and controls. Then, our function assumes that the variance of the log OR
#' is equal to the mean of the variance of all possible situations (the reported standard error is the square root of that mean variance).
#' This estimation thus necessarily comes with some imprecision and should not
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
#'  \code{converted effect size measure} \tab OR + RR + NNT + RD\cr
#'  \code{} \tab D + G + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 2. Odds Ratio'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @note
#' **The imputed standard error is systematically too wide, and increasingly so in
#' larger studies.** With no reported uncertainty to work from, this function
#' enumerates every 2x2 table compatible with the odds ratio and the case/control
#' margins and averages the log-OR *variance*. Variance is convex in the cell counts
#' and diverges as any cell approaches zero, so the average is pulled upward by
#' near-degenerate tables that carry the same weight as plausible ones.
#'
#' Measured against the standard error the realised table would have given, over
#' simulated studies crossing exposure prevalence, baseline risk, odds ratio and
#' sample size, the median ratio of imputed to true SE is about **1.4**, and it is not
#' a constant offset: because the flat enumeration weight puts \eqn{1/(n\_cases - 1)}
#' on the table with one case in the exposed arm, whose \eqn{1/a} term is 1, the
#' inflation **grows with study size** -- roughly 1.2x in the smallest studies to 2.0x
#' in the largest, at fixed exposure prevalence.
#'
#' Two consequences worth knowing. An SE 1.4x too wide gives that study about half its
#' correct inverse-variance weight; and because the inflation is size-dependent rather
#' than uniform, it also distorts weights *between* imputed studies, which shifts the
#' pooled point estimate rather than merely widening its interval.
#'
#' The remedy is to supply the uncertainty rather than have it imputed: prefer
#' \code{\link{es_from_or_se}()}, \code{\link{es_from_or_ci}()} or
#' \code{\link{es_from_or_pval}()} whenever the source reports a standard error, a
#' confidence interval or a p-value. \code{\link{convert_df}()} does this
#' automatically -- when a row reports an SE or a CI, this route is skipped entirely
#' rather than contributing a dominated second estimate.
#'
#' @examples
#' es_or_guess <- es_from_or(or = 0.5, n_cases = 210, n_controls = 220)
#' es_or <- es_from_or_se(or = 0.5, logor_se = 0.4, n_cases = 210, n_controls = 220)
#' es_or_guess$logor_se # standard error simulated from the margins (~1.4x too wide)
#' es_or$logor_se # standard error as reported by the source study
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
#' @param small_margin_prop smallest margin proportion of the underlying 2x2 table (a proportion in (0, 0.5])
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
#'  \code{converted effect size measure} \tab RR + NNT + RD\cr
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
    n_sample <- rep(NA_real_, length(or))
  }
  reverse_or[is.na(reverse_or)] <- FALSE

  or <- ifelse(is.na(or) & !is.na(logor), exp(logor), or)


  logor_ci_lo <- ifelse(is.na(logor_ci_lo) & !is.na(or_ci_lo), log(or_ci_lo), logor_ci_lo)
  logor_ci_up <- ifelse(is.na(logor_ci_up) & !is.na(or_ci_up), log(or_ci_up), logor_ci_up)
  logor_se <- .ci_width(logor_ci_lo, logor_ci_up) / (2 * qnorm(.975))

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

#' Convert an odds ratio value and its p-value to several effect size measures
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
#' @param small_margin_prop smallest margin proportion of the underlying 2x2 table (a proportion in (0, 0.5])
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
#'  \code{converted effect size measure} \tab RR + NNT + RD\cr
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
    n_sample <- rep(NA_real_, length(or))
  }
  reverse_or_pval[is.na(reverse_or_pval)] <- FALSE

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

#' Convert an odds ratio or risk ratio and a Wald t-statistic from a regression model into several effect size measures
#'
#' @param or odds ratio value (from logistic regression)
#' @param logor log odds ratio value
#' @param rr risk ratio value (from log-binomial or modified Poisson regression)
#' @param logrr log risk ratio value
#' @param logreg_t a t-statistic (Wald statistic) from a logistic, log-binomial, or modified Poisson regression model
#' @param n_cases number of cases/events across exposed/non-exposed groups
#' @param n_controls number of controls/no-event across exposed/non-exposed groups
#' @param n_exp number of participants in the exposed group
#' @param n_nexp number of participants in the non-exposed group
#' @param n_sample total number of participants in the sample
#' @param baseline_risk proportion of cases in the non-exposed group
#' @param small_margin_prop smallest margin proportion of cases/events in the underlying 2x2 table (a proportion in (0, 0.5])
#' @param reverse_logreg_t a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param or_to_rr formula used to convert the \code{or} value into a risk ratio (see details).
#' @param or_to_cor formula used to convert the \code{or} value into a correlation coefficient (see details).
#' @param rr_to_or formula used to convert the \code{rr} value into an odds ratio (see details).
#'
#' @details
#' This function derives the standard error of the log odds ratio (or log risk ratio) from a
#' Wald t-statistic reported in a regression model.
#'
#' **To estimate the standard error of the log OR (or log RR)**, the formulas used are:
#' \deqn{t = \frac{\beta}{SE(\beta)}}
#' \deqn{SE(\beta) = \frac{|\beta|}{|t|}}
#' where \eqn{\beta} is \eqn{\log(OR)} or \eqn{\log(RR)} depending on the model.
#'
#' Then, if an OR (or logOR) is entered, calculations of \code{\link{es_from_or_se}()} are applied.
#' If a RR (or logRR) is entered, calculations of \code{\link{es_from_rr_se}()} are applied.
#'
#' Note that the standardized-mean-difference and correlation conversions (D, G, R, Z) are
#' produced for **OR inputs only**. RR inputs are treated as a ratio measure and yield
#' RR + OR + NNT + RD: RR is not converted to a standardized mean difference or correlation,
#' because that would require going through the OR and the baseline risk (see
#' \code{\link{es_from_rr_se}}). To obtain a SMD or correlation from an RR, convert it to an
#' OR first (supplying the baseline risk) and then use the OR path.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR \cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR inputs: D + G + R + Z (+ RR + NNT + RD); RR inputs: OR + NNT + RD \cr
#' }
#'
#' @references
#' Sanchez-Meca, J., Marin-Martinez, F., & Chacon-Moscoso, S. (2003). Effect-size indices for
#' dichotomized outcomes in meta-analysis. \emph{Psychological Methods}, 8(4), 448--467.
#'
#' @md
#'
#' @export es_from_logreg_t
#'
#' @examples
#' es_or <- es_from_logreg_t(
#'   or = 2.12, logreg_t = 3.21,
#'   n_cases = 50, n_controls = 150
#' )
#'
#' es_rr <- es_from_logreg_t(
#'   rr = 1.5, logreg_t = 2.8,
#'   n_exp = 100, n_nexp = 100
#' )
es_from_logreg_t <- function(or, logor, rr, logrr, logreg_t,
                         baseline_risk, small_margin_prop,
                         n_exp, n_nexp, n_cases, n_controls, n_sample,
                         or_to_rr = "metaumbrella_cases",
                         or_to_cor = "bonett",
                         rr_to_or = "metaumbrella",
                         reverse_logreg_t) {

  len <- if (!missing(or)) length(or) else if (!missing(logor)) length(logor) else if (!missing(rr)) length(rr) else if (!missing(logrr)) length(logrr) else length(logreg_t)

  if (missing(or)) or <- rep(NA_real_, len)
  if (missing(logor)) logor <- rep(NA_real_, len)
  if (missing(rr)) rr <- rep(NA_real_, len)
  if (missing(logrr)) logrr <- rep(NA_real_, len)
  if (missing(baseline_risk)) baseline_risk <- rep(NA_real_, len)
  if (missing(small_margin_prop)) small_margin_prop <- rep(NA_real_, len)
  if (missing(n_exp)) n_exp <- rep(NA_real_, len)
  if (missing(n_nexp)) n_nexp <- rep(NA_real_, len)
  if (missing(n_cases)) n_cases <- rep(NA_real_, len)
  if (missing(n_controls)) n_controls <- rep(NA_real_, len)
  if (missing(n_sample)) n_sample <- rep(NA_real_, len)
  if (missing(logreg_t)) logreg_t <- rep(NA_real_, len)
  if (missing(reverse_logreg_t)) reverse_logreg_t <- rep(FALSE, len)
  reverse_logreg_t[is.na(reverse_logreg_t)] <- FALSE

  has_or <- !is.na(or) | !is.na(logor)
  has_rr <- !is.na(rr) | !is.na(logrr)

  # OR -------
  or <- ifelse(is.na(or) & !is.na(logor), exp(logor), or)
  logOR <- suppressWarnings(log(or))
  logor_se <- abs(logOR / logreg_t)

  es_or <- es_from_or_se(
    or = or, logor_se = logor_se,
    baseline_risk = baseline_risk, small_margin_prop = small_margin_prop,
    n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample,
    n_cases = n_cases, n_controls = n_controls,
    or_to_cor = or_to_cor, or_to_rr = or_to_rr,
    reverse_or = reverse_logreg_t
  )

  # RR -------
  rr <- ifelse(is.na(rr) & !is.na(logrr), exp(logrr), rr)
  logRR <- suppressWarnings(log(rr))
  logrr_se <- abs(logRR / logreg_t)

  es_rr <- es_from_rr_se(
    rr = rr, logrr_se = logrr_se,
    baseline_risk = baseline_risk,
    n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls,
    rr_to_or = rr_to_or,
    reverse_rr = reverse_logreg_t
  )

  use_rr <- has_rr & !has_or

  if (!any(use_rr)) {
    es_or$info_used <- "logreg_t"
    return(es_or)
  } else if (all(use_rr)) {
    es_rr$info_used <- "logreg_t"
    return(es_rr)
  }

  all_cols <- union(names(es_or), names(es_rr))
  es <- data.frame(matrix(NA_real_, nrow = len, ncol = 0))
  for (col in all_cols) {
    if (col == "info_used") next
    or_val <- if (col %in% names(es_or)) es_or[[col]] else rep(NA_real_, len)
    rr_val <- if (col %in% names(es_rr)) es_rr[[col]] else rep(NA_real_, len)
    es[[col]] <- ifelse(use_rr, rr_val, or_val)
  }
  es$info_used <- "logreg_t"

  return(es)
}
