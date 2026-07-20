#' Convert a 2x2 table into several effect size measures
#'
#' @param n_cases_exp number of cases/events in the exposed group
#' @param n_cases_nexp number of cases/events in the non exposed group
#' @param n_controls_exp number of controls/no-event in the exposed group
#' @param n_controls_nexp number of controls/no-event in the non exposed group
#' @param table_2x2_to_cor formula used to obtain a correlation coefficient from the contingency table (see details).
#' @param reverse_2x2 a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first computes (log) odds ratio (OR), (log) risk ratio (RR) and number needed to treat (NNT)
#' from the 2x2 table. Note that if a cell is equal to 0, we applied the typical adjustment (add 0.5) to all cells.
#' This adjustment is used for the OR/RR only; the RD and NNT are obtained from the raw cell counts.
#' Cohen's d (D), Hedges' g (G) and correlation coefficients (R/Z) are then estimated from the OR.
#'
#' **To estimate an OR**, the formulas used (Box 6.4.a in the Cochrane Handbook) are:
#' \deqn{logor = log(\frac{n\_cases\_exp / n\_cases\_nexp}{n\_controls\_exp / n\_controls\_nexp})}
#' \deqn{logor\_se = \sqrt{\frac{1}{n\_cases\_exp} + \frac{1}{n\_cases\_nexp} + \frac{1}{n\_controls\_exp} + \frac{1}{n\_controls\_nexp}}}
#'
#' **To estimate an RR**, the formulas used (Box 6.4.a in the Cochrane Handbook) are:
#' \deqn{logrr = log(\frac{n\_cases\_exp / n\_exp}{n\_cases\_nexp / n\_nexp})}
#' \deqn{logrr\_se = \sqrt{\frac{1}{n\_cases\_exp} - \frac{1}{n\_exp} + \frac{1}{n\_cases\_nexp} - \frac{1}{n\_nexp}}}
#'
#' **To estimate a risk difference (RD) and NNT**, the formulas used are (Wen et al., 2005; Altman, 1998):
#' \deqn{pt = \frac{n\_cases\_exp}{n\_cases\_exp + n\_controls\_exp}}
#' \deqn{pc = \frac{n\_cases\_nexp}{n\_cases\_nexp + n\_controls\_nexp}}
#' \deqn{rd = pc - pt}
#' \deqn{rd\_se = \sqrt{\frac{pt(1-pt)}{n\_exp} + \frac{pc(1-pc)}{n\_nexp}}}
#' \deqn{nnt = \frac{1}{rd}}
#' \deqn{nnt\_se = \frac{rd\_se}{rd^2}}
#' Note that NNT confidence intervals are set to NA when the RD confidence interval crosses zero
#' (discontinuous CI; Altman, 1998).
#'
#' **Direction convention.** The risk difference is defined as \eqn{rd = pc - pt}
#' (control risk minus exposed risk), so a POSITIVE RD means the control group has the
#' higher risk. This is the OPPOSITE direction to the OR and RR produced from the same
#' 2x2 table, which are exposed-over-non-exposed (an OR/RR \eqn{> 1} means the exposed
#' group has the higher risk). Consequently, for the same table, a protective exposure
#' yields \eqn{OR < 1}, \eqn{RR < 1} but \eqn{RD > 0}; keep this in mind when pooling RD
#' alongside OR/RR, and use \code{reverse_2x2} if you need to align the directions.
#'
#' **To convert the 2x2 table into a SMD**,
#' the function estimates an OR value from the 2x2 table (formula above)
#' that is then converted to a SMD
#' (see formula in \code{\link{es_from_or_se}()}).
#'
#' **To convert the 2x2 table into a correlation coefficient**,
#' For now, only the tetrachoric correlation is currently proposed
#' - \code{table_2x2_to_cor = "tetrachoric"}.
#' Given the heavy calculations required for this effect size measure,
#' we relied on the implementation of the formulas of the 'metafor' package. More
#' information can be retrieved here
#' (https://wviechtb.github.io/metafor/reference/escalc.html#-b-measures-for-two-dichotomous-variables).
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT + RD\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 7. Contingency (2x2) table or proportions'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @export es_from_2x2
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Available from www.training.cochrane.org/handbook.
#'
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' Sedgwick, P. (2013). What is number needed to treat (NNT)? Bmj, 347.
#'
#' Altman, D. G. (1998). Confidence intervals for the number needed to treat. BMJ, 317(7168), 1309-1312.
#'
#' Wen, S., Zhang, L., & Yang, B. (2005). Two approaches to incorporate clinical data uncertainty into number needed to treat. Journal of Clinical Pharmacy and Therapeutics, 30(2), 105-109.
#'
#' @examples
#' es_from_2x2(n_cases_exp = 467, n_cases_nexp = 22087, n_controls_exp = 261, n_controls_nexp = 8761)
es_from_2x2 <- function(n_cases_exp, n_cases_nexp,
                        n_controls_exp, n_controls_nexp,
                        table_2x2_to_cor = "tetrachoric", reverse_2x2) {
  if (missing(reverse_2x2)) reverse_2x2 <- rep(FALSE, length(n_cases_exp))
  reverse_2x2[is.na(reverse_2x2)] <- FALSE
  if (length(reverse_2x2) == 1) reverse_2x2 = c(rep(reverse_2x2, length(n_cases_exp)))
  if (length(reverse_2x2) != length(n_cases_exp)) stop("The length of the 'reverse_2x2' argument is incorrectly specified.")

  if (!all(table_2x2_to_cor %in% c("tetrachoric"))) { # , "cooper", "lipsey"
    stop(paste0(
      "'",
      unique(table_2x2_to_cor[!table_2x2_to_cor %in% c("tetrachoric")]), #, "cooper", "lipsey"
      "' not in tolerated values for the 'table_2x2_to_cor' argument.",
      " Possible inputs are: 'tetrachoric', "#'cooper', 'lipsey'
    ))
  }

  # rd from raw counts
  n_exp_raw <- n_cases_exp + n_controls_exp
  n_nexp_raw <- n_cases_nexp + n_controls_nexp
  pc_raw <- n_cases_nexp / n_nexp_raw
  pt_raw <- n_cases_exp / n_exp_raw
  rd <- pc_raw - pt_raw
  rd_se <- sqrt(pt_raw * (1 - pt_raw) / n_exp_raw + pc_raw * (1 - pc_raw) / n_nexp_raw)

  # 0.5 correction for or/rr
  zero <- which(n_cases_exp == 0 | n_cases_nexp == 0 | n_controls_exp == 0 | n_controls_nexp == 0)
  n_cases_exp[zero] <- n_cases_exp[zero] + 0.5
  n_cases_nexp[zero] <- n_cases_nexp[zero] + 0.5
  n_controls_exp[zero] <- n_controls_exp[zero] + 0.5
  n_controls_nexp[zero] <- n_controls_nexp[zero] + 0.5
  n_exp <- n_cases_exp + n_controls_exp
  n_nexp <- n_cases_nexp + n_controls_nexp

  or_raw <- suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp))
  or <- ifelse(reverse_2x2, 1 / or_raw, or_raw)
  se_or <- suppressWarnings(sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp))

  rr_raw <- suppressWarnings((n_cases_exp / n_exp) / (n_cases_nexp / n_nexp))
  rr <- ifelse(reverse_2x2, 1 / rr_raw, rr_raw)
  se_rr <- suppressWarnings(sqrt(1 / n_cases_exp - 1 / n_exp + 1 / n_cases_nexp - 1 / n_nexp))

  logOR <- suppressWarnings(log(or))
  d <- logOR * sqrt(3) / pi
  d_se <- sqrt(se_or^2 * 3 / (pi^2))

  es <- .es_from_d(
    d = d, d_se = d_se,
    n_exp = n_exp,
    n_nexp = n_nexp,
    smd_to_cor = rep("lipsey_cooper", length(d))
  )

  es$logor <- logOR
  es$logor_se <- se_or
  es$logor_ci_lo <- es$logor - qnorm(.975) * es$logor_se
  es$logor_ci_up <- es$logor + qnorm(.975) * es$logor_se

  es$logrr <- log(rr)
  es$logrr_se <- se_rr
  es$logrr_ci_lo <- es$logrr - qnorm(.975) * es$logrr_se
  es$logrr_ci_up <- es$logrr + qnorm(.975) * es$logrr_se

  dat2x2 <- data.frame(
    n_cases_exp = n_cases_exp, n_controls_exp = n_controls_exp,
    n_cases_nexp = n_cases_nexp, n_controls_nexp = n_controls_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    table_2x2_to_cor = table_2x2_to_cor, reverse_2x2 = reverse_2x2
  )

  nn_miss <- which(
    (dat2x2$table_2x2_to_cor == "lipsey" |
       dat2x2$table_2x2_to_cor == "tetrachoric") &
      !is.na(dat2x2$n_cases_exp) & !is.na(dat2x2$n_controls_exp) &
      !is.na(dat2x2$n_cases_nexp) & !is.na(dat2x2$n_controls_nexp)
  )

  if (length(nn_miss) != 0) {
    # no_cooper = which(dat2x2$table_2x2_to_cor != "cooper")
    #
    # es$r[no_cooper] <- es$r_se[no_cooper] <-
    #   es$r_ci_lo[no_cooper] <- es$r_ci_up[no_cooper] <-
    #   es$z[no_cooper] <- es$z_se[no_cooper] <-
    #   es$z_ci_lo[no_cooper] <- es$z_ci_up[no_cooper] <- NA

    res_tet <- suppressWarnings(t(mapply(.contingency_to_cor,
      n_cases_exp = dat2x2$n_cases_exp[nn_miss],
      n_controls_exp = dat2x2$n_controls_exp[nn_miss],
      n_cases_nexp = dat2x2$n_cases_nexp[nn_miss],
      n_controls_nexp = dat2x2$n_controls_nexp[nn_miss],
      table_2x2_to_cor = dat2x2$table_2x2_to_cor[nn_miss],
      reverse_2x2 = dat2x2$reverse_2x2[nn_miss]
    )))

    es$r[nn_miss] <- res_tet[, 1]
    es$r_se[nn_miss] <- suppressWarnings(sqrt(res_tet[, 2]))
    es$r_ci_lo[nn_miss] <- res_tet[, 3]
    es$r_ci_up[nn_miss] <- res_tet[, 4]

    es$z[nn_miss] <- res_tet[, 5]
    es$z_se[nn_miss] <- suppressWarnings(sqrt(res_tet[, 6]))
    es$z_ci_lo[nn_miss] <- res_tet[, 7]
    es$z_ci_up[nn_miss] <- res_tet[, 8]
  }

  es$rd <- ifelse(reverse_2x2, -rd, rd)
  es$rd_se <- rd_se
  es$rd_ci_lo <- es$rd - qnorm(.975) * rd_se
  es$rd_ci_up <- es$rd + qnorm(.975) * rd_se

  es$nnt <- ifelse(rd == 0, NA, ifelse(reverse_2x2, -1 / rd, 1 / rd))
  es$nnt_se <- ifelse(rd == 0, NA, rd_se / rd^2)
  rd_ci_lo_raw <- rd - qnorm(.975) * rd_se
  rd_ci_up_raw <- rd + qnorm(.975) * rd_se
  crosses_zero <- (rd_ci_lo_raw < 0 & rd_ci_up_raw > 0)
  es$nnt_ci_lo <- ifelse(crosses_zero | rd == 0, NA,
                          ifelse(reverse_2x2, -1 / rd_ci_lo_raw, 1 / rd_ci_up_raw))
  es$nnt_ci_up <- ifelse(crosses_zero | rd == 0, NA,
                          ifelse(reverse_2x2, -1 / rd_ci_up_raw, 1 / rd_ci_lo_raw))

  es$info_used <- "2x2"
  return(es)
}

#' Convert a table with the number of cases and row marginal sums into several effect size measures
#'
#' @param n_cases_exp number of cases/events in the exposed group
#' @param n_cases_nexp number of cases/events in the non exposed group
#' @param n_exp total number of participants in the exposed group
#' @param n_nexp total number of participants in the non exposed group
#' @param table_2x2_to_cor formula used to obtain a correlation coefficient from the contingency table (see details).
#' @param reverse_2x2 a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function uses the number of cases in both the exposed
#' and non-exposed groups and the total number of participants exposed and non-exposed
#' to recreate a 2x2 table.
#' Then relies on the calculations of the \code{\link{es_from_2x2}} function.
#' \deqn{n\_controls\_exp = n\_exp - n\_cases\_exp}
#' \deqn{n\_controls\_nexp = n\_nexp - n\_cases\_nexp}
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT + RD\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 7. Contingency (2x2) table or proportions'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#'
#' @md
#'
#' @export es_from_2x2_sum
#'
#' @examples
#' es_from_2x2_sum(n_cases_exp = 10, n_exp = 40, n_cases_nexp = 25, n_nexp = 47)
es_from_2x2_sum <- function(n_cases_exp, n_exp, n_cases_nexp, n_nexp,
                            table_2x2_to_cor = "tetrachoric", reverse_2x2) {
  if (missing(reverse_2x2)) reverse_2x2 <- rep(FALSE, length(n_cases_exp))
  reverse_2x2[is.na(reverse_2x2)] <- FALSE

  es <- es_from_2x2(
    n_cases_exp = n_cases_exp,
    n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_exp - n_cases_exp,
    n_controls_nexp = n_nexp - n_cases_nexp,
    table_2x2_to_cor = table_2x2_to_cor,
    reverse_2x2 = reverse_2x2
  )

  es$info_used <- "2x2_sum"

  return(es)
}

#' Convert the proportion of occurrence of a binary event in two independent groups into several effect size measures
#'
#' @param prop_cases_exp proportion of cases/events in the exposed group (ranging from 0 to 1)
#' @param prop_cases_nexp proportion of cases/events in the non-exposed group (ranging from 0 to 1)
#' @param n_exp total number of participants in the exposed group
#' @param n_nexp total number of participants in the non exposed group
#' @param table_2x2_to_cor formula used to obtain a correlation coefficient from the contingency table (see details).
#' @param reverse_prop a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#'
#' @details
#' This function uses the proportions and sample size to
#' recreate the 2x2 table, and
#' then relies on the calculations of the \code{\link{es_from_2x2_sum}()} function.
#'
#' The formulas used is to obtain the 2x2 table are
#' \deqn{n\_cases\_exp = prop\_cases\_exp * n\_exp}
#' \deqn{n\_cases\_nexp = prop\_cases\_nexp * n\_nexp}
#' \deqn{n\_controls\_exp = (1 - prop\_cases\_exp) * n\_exp}
#' \deqn{n\_controls\_nexp = (1 - prop\_cases\_nexp) * n\_nexp}
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT + RD\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 7. Contingency (2x2) table or proportions'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_2x2_prop
#'
#'
#' @md
#'
#' @examples
#' es_from_2x2_prop(prop_cases_exp = 0.80, prop_cases_nexp = 0.60, n_exp = 10, n_nexp = 20)
es_from_2x2_prop <- function(prop_cases_exp, prop_cases_nexp, n_exp, n_nexp,
                             table_2x2_to_cor = "tetrachoric", reverse_prop) {
  if (missing(reverse_prop)) reverse_prop <- rep(FALSE, length(prop_cases_exp))
  reverse_prop[is.na(reverse_prop)] <- FALSE


  n_cases_exp <- round(prop_cases_exp * n_exp)
  n_cases_nexp <- round(prop_cases_nexp * n_nexp)
  n_controls_exp <- n_exp - n_cases_exp
  n_controls_nexp <- n_nexp - n_cases_nexp

  es <- es_from_2x2(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp,
    table_2x2_to_cor = table_2x2_to_cor, reverse_2x2 = reverse_prop
  )

  es$info_used <- "2x2_prop"

  return(es)
}


