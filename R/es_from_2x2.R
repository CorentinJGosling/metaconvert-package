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
#' This adjustment feeds the OR, the RR, the SMD (D/G) converted from the OR, **and the correlation
#' coefficients (R/Z)**, which are solved from the corrected table; only the RD and the NNT are
#' obtained from the raw cell counts.
#' Cohen's d (D) and Hedges' g (G) are then estimated from the OR. The correlation coefficients
#' (R/Z) are **not** obtained from the OR: they come from the tetrachoric solve on the 2x2 table
#' itself (see below), which is a different route and returns a different value.
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
#' (control risk minus exposed risk), so a positive RD means the control group has the
#' higher risk. This is the opposite direction to the OR and RR produced from the same
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
#' **To convert the 2x2 table into a correlation coefficient**, only the tetrachoric
#' correlation is currently available:
#' - \code{table_2x2_to_cor = "tetrachoric"}.
#' Given the heavy calculations required for this effect size measure,
#' we relied on the implementation of the formulas of the 'metafor' package. More
#' information can be retrieved here
#' (https://wviechtb.github.io/metafor/reference/escalc.html#-b-measures-for-two-dichotomous-variables).
#'
#' That parity holds **on a table with no zero cell**. On a table with a zero cell it does not, and
#' the difference is large. 'metafor' deliberately does not apply its \code{add}/\code{to} continuity
#' correction to \code{"RTET"}, so it returns the boundary maximum-likelihood estimate
#' (\eqn{r = \pm 1}) with an enormous variance -- its way of reporting that the correlation is not
#' identified. \emph{metaConvert} solves the +0.5-corrected table instead, so it returns an interior
#' estimate with an ordinary-looking, finite standard error. That standard error is the sampling SE
#' \emph{of the shrunk table} and does not account for the shrinkage, so it is anti-conservative for
#' the value actually reported: on \eqn{a/b/c/d = 0/20/10/10} \emph{metaConvert} gives
#' \eqn{r = -0.856} with \eqn{SE = 0.113}, against \eqn{r = -1} with \eqn{SE = 149} from 'metafor' --
#' an inverse-variance weight about 1300 times larger. This is a deliberate choice, not an oversight:
#' the boundary estimate cannot be pooled. But treat the correlation from any zero-cell table as a
#' shrunk value whose precision is overstated, and consider excluding such rows or handling them in a
#' sensitivity analysis.
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

  if (!all(table_2x2_to_cor %in% c("tetrachoric"))) {
    stop(paste0(
      "'",
      paste(unique(table_2x2_to_cor[!table_2x2_to_cor %in% c("tetrachoric")]),
            collapse = "', '"),
      "' not in tolerated values for the 'table_2x2_to_cor' argument. ",
      "The only possible input is 'tetrachoric'.\n",
      "  A 2x2 table is converted to a correlation by assuming both binary variables ",
      "are dichotomised continua, and estimating the correlation of those latent ",
      "variables. This is a deliberate restriction, not a temporary one -- see ",
      "?convert_df. If your variables are genuinely dichotomous rather than ",
      "dichotomised, prefer a binary effect size (measure = 'logor', 'rr' or 'rd')."
    ), call. = FALSE)
  }

  # A cell count cannot be negative. The containment lives here rather than in the two
  # wrappers because es_from_2x2() is itself exported and equally exposed, and because
  # es_from_2x2_sum() (n_controls_exp = n_exp - n_cases_exp) and es_from_2x2_prop()
  # (the same after round(prop * n)) manufacture the impossible cell inside the route,
  # where no column-keyed Tier-1 check in convert_df() can see it. Guarding the parent
  # covers all three from the single audited place R/internal_guards.R calls for.
  #
  # Nothing else declines the row: the OR arm divides by the negative cell and goes
  # NaN, but the RR arm reads only the row margins -- rr = (a/n_exp)/(c/n_nexp) -- so it
  # emits a finite, correctly-signed, plausible-looking log risk ratio with a finite SE
  # and CI from a table that cannot exist, and the RD/NNT arm survives every mild
  # overshoot (41 cases of 40 -> rd = -0.525 with se 0.0749, inside B6 bound). The
  # whole table is neutralised, not just the offending cell: at least one of the four
  # is wrong and which one cannot be inferred, exactly as V7 treats an inconsistent
  # additive triple.
  bad_cell <- (!is.na(n_cases_exp) & n_cases_exp < 0) |
    (!is.na(n_cases_nexp) & n_cases_nexp < 0) |
    (!is.na(n_controls_exp) & n_controls_exp < 0) |
    (!is.na(n_controls_nexp) & n_controls_nexp < 0)
  n_cases_exp <- ifelse(bad_cell, NA_real_, n_cases_exp)
  n_cases_nexp <- ifelse(bad_cell, NA_real_, n_cases_nexp)
  n_controls_exp <- ifelse(bad_cell, NA_real_, n_controls_exp)
  n_controls_nexp <- ifelse(bad_cell, NA_real_, n_controls_nexp)

  # rd from raw counts
  n_exp_raw <- n_cases_exp + n_controls_exp
  n_nexp_raw <- n_cases_nexp + n_controls_nexp
  pc_raw <- n_cases_nexp / n_nexp_raw
  pt_raw <- n_cases_exp / n_exp_raw
  rd <- pc_raw - pt_raw
  rd_se <- suppressWarnings(sqrt(pt_raw * (1 - pt_raw) / n_exp_raw + pc_raw * (1 - pc_raw) / n_nexp_raw))

  # A risk-difference SE of exactly 0 is a zero sampling variance -- an infinite
  # inverse-variance weight, or an rma() abort. On a 2x2 it means a double-zero (or
  # double-full) table, where both arm variances vanish. metafor never emits such a row
  # either: escalc(measure = "RD") continuity-corrects it, or returns NA under
  # drop00 = TRUE. Keyed on the SE and never on rd itself -- a balanced 10/50 vs 10/50
  # table has rd = 0 with a perfectly good SE and must survive. See R/internal_guards.R.
  degenerate_rd <- !is.na(rd_se) & rd_se <= 0
  rd_se <- .positive_or_na(rd_se)
  rd <- ifelse(degenerate_rd, NA_real_, rd)

  # The tetrachoric solve below must see the RAW table, so keep a copy before the cells
  # are corrected in place. The +0.5 is a ratio-measure device: it exists because a zero
  # cell makes the odds ratio and the risk ratio undefined, and @details documents it as
  # applying to the OR/RR only. metafor deliberately does NOT apply it to measure =
  # "RTET" (not even under to = "all"), because the tetrachoric is not undefined on such
  # a table -- its ML lands on the boundary r = +-1 with an enormous variance, which is
  # how the estimator says the correlation is not identified by these counts. Shrinking
  # the table first replaces that honest refusal with an interior estimate carrying an
  # ordinary-looking SE that does not account for the shrinkage: on 0/20/10/10 it
  # reported r = -0.856 with se = 0.113 against metafor's -1 with se = 149.0, i.e. a
  # 1313x inverse-variance weight on a study whose correlation the data cannot pin down.
  raw_cases_exp     <- n_cases_exp
  raw_controls_exp  <- n_controls_exp
  raw_cases_nexp    <- n_cases_nexp
  raw_controls_nexp <- n_controls_nexp

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

  # Raw cells, not the +0.5-corrected ones -- see the note beside the correction above.
  dat2x2 <- data.frame(
    n_cases_exp = raw_cases_exp, n_controls_exp = raw_controls_exp,
    n_cases_nexp = raw_cases_nexp, n_controls_nexp = raw_controls_nexp,
    n_exp = raw_cases_exp + raw_controls_exp,
    n_nexp = raw_cases_nexp + raw_controls_nexp,
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

    res_tet <- suppressWarnings(.mapply_memo(.contingency_to_cor,
      n_cases_exp = dat2x2$n_cases_exp[nn_miss],
      n_controls_exp = dat2x2$n_controls_exp[nn_miss],
      n_cases_nexp = dat2x2$n_cases_nexp[nn_miss],
      n_controls_nexp = dat2x2$n_controls_nexp[nn_miss],
      table_2x2_to_cor = dat2x2$table_2x2_to_cor[nn_miss],
      reverse_2x2 = dat2x2$reverse_2x2[nn_miss]
    ))

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
  # Non-strict, with the rd == 0 clause folded in as on the four sibling routes: a
  # bound landing exactly on 0 is still the Altman discontinuity, and the reciprocal of
  # +0 is a literal Inf where the neighbouring input returns NA.
  crosses_zero <- (rd_ci_lo_raw <= 0 & rd_ci_up_raw >= 0) | rd == 0
  es$nnt_ci_lo <- ifelse(crosses_zero, NA,
                          ifelse(reverse_2x2, -1 / rd_ci_lo_raw, 1 / rd_ci_up_raw))
  es$nnt_ci_up <- ifelse(crosses_zero, NA,
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
#' then relies on the calculations of the \code{\link{es_from_2x2}()} function.
#'
#' The formulas used to obtain the 2x2 table are
#' \deqn{n\_cases\_exp = round(prop\_cases\_exp * n\_exp)}
#' \deqn{n\_cases\_nexp = round(prop\_cases\_nexp * n\_nexp)}
#' \deqn{n\_controls\_exp = n\_exp - n\_cases\_exp}
#' \deqn{n\_controls\_nexp = n\_nexp - n\_cases\_nexp}
#' The numbers of cases are rounded to integers so that the reconstructed table is a
#' genuine contingency table. The resulting effect sizes therefore differ slightly from
#' those obtained from the unrounded cell counts.
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


  # A proportion outside [0, 1] is neutralised per row, following the convention of
  # es_from_prop_single_group(): a warning rather than a silent NA or a hard stop, since
  # one bad cell must never abort a whole convert_df() run. The @param already says
  # "ranging from 0 to 1"; until now that was documented and unenforced, and the trigger
  # is the same percentage mix-up baseline_risk is already guarded for. Left unchecked,
  # round(1.2 * 100) = 120 cases out of 100 manufactures a negative control cell that
  # only es_from_2x2() own guard would see, and the row would reach the pool as
  # es_crude = 24 with no flag.
  bad_prop <- (!is.na(prop_cases_exp) & (prop_cases_exp < 0 | prop_cases_exp > 1)) |
    (!is.na(prop_cases_nexp) & (prop_cases_nexp < 0 | prop_cases_nexp > 1))
  if (any(bad_prop)) {
    warning(
      "Proportions outside [0, 1] were passed to es_from_2x2_prop() (row(s) ",
      paste(which(bad_prop), collapse = ", "),
      "). A proportion is not a percentage: these rows are set to NA.",
      call. = FALSE
    )
    prop_cases_exp <- ifelse(bad_prop, NA_real_, prop_cases_exp)
    prop_cases_nexp <- ifelse(bad_prop, NA_real_, prop_cases_nexp)
  }

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


