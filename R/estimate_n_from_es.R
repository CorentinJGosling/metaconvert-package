#' Estimate sample size and time for the number of cases exposed and non-exposed from IRR, variance and overall number of cases
#'
#' @param irr IRR
#' @param var variance
#' @param n_cases number of cases
#' @param time time
#' @param time_exp time of exposed group
#' @param time_nexp time of non exposed group
#'
#' @noRd
.estimate_n_from_irr <- function(irr, var, n_cases, time = NA, time_exp = NA, time_nexp = NA) {
  # (re)-estimate n and time from irr and var
  n_cases_exp <- optimize(function(n_cases_exp, var, n_cases) {
    (1 / n_cases_exp + 1 / (n_cases - n_cases_exp) - var)^2
  }, c(0, n_cases), var, n_cases)$minimum
  n_cases_nexp <- n_cases - n_cases_exp
  if (is.na(time) && is.na(time_exp) && is.na(time_nexp)) {
    time <- 1
  }
  if (!is.na(time)) {
    time_exp <- n_cases_exp / (n_cases_exp + irr * n_cases_nexp) * time
    time_nexp <- time - time_exp
  } else if (!is.na(time_exp)) {
    time <- time_exp * (n_cases_exp + irr * n_cases_nexp) / n_cases_exp
    time_nexp <- time - time_exp
  } else {
    time_exp <- (n_cases_exp * time_nexp) / (irr * n_cases_nexp)
    time <- time_exp + time_nexp
  }
  return(data.frame(n_cases_exp, time_exp, n_cases_nexp, time_nexp))
}

# NB two further helpers -- .estimate_n_from_or_and_n_cases() and
# .estimate_n_from_or_and_n_exp() -- were defined here as well until they were
# removed. They were SHADOWED DEAD CODE: identically-named definitions live in
# R/internal_multiple_formulas.R, and because DESCRIPTION has no Collate field the
# files are sourced alphabetically, so internal_multiple_formulas.R ("i") overwrote
# this file ("e"). Patching the copies here had no effect, and by the end they had
# drifted: the live versions solve the 2x2 exactly from a second margin, while these
# still only searched by matching the reported variance. Edit the live copies.
#
# The two helpers below are NOT duplicated and ARE live -- do not remove this file.

#' Estimate the number of cases and controls exposed and non-exposed given the risk ratio
#'
#' @param rr RR
#' @param var variance
#' @param n_cases number of cases
#' @param n_controls number of controls
#'
#' @noRd
.estimate_n_from_rr <- function(rr, var, n_cases, n_controls) {
  # uncorrected
  n_cases_nexp_sim1 <- 0:n_cases
  n_controls_nexp_sim1 <- round(n_cases_nexp_sim1 * ((rr * (n_cases + n_controls)) / (n_cases + (rr - 1) * n_cases_nexp_sim1) - 1))
  n_cases_exp_sim1 <- n_cases - n_cases_nexp_sim1
  n_controls_exp_sim1 <- n_controls - n_controls_nexp_sim1

  # we take only positives (no-zero)
  idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
  n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
  n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
  n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
  n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

  # corregint 0.5
  n_cases_nexp_sim2 <- 0:n_cases
  n_controls_nexp_sim2 <- ((n_cases + n_controls - n_cases_nexp_sim2 + 1) - (n_cases + n_controls + 2) * (n_cases - n_cases_nexp_sim2 + 0.5) / ((n_cases_nexp_sim2 + 0.5) * rr + n_cases - n_cases_nexp_sim2 + 0.5))
  n_cases_exp_sim2 <- n_cases - n_cases_nexp_sim2
  n_controls_exp_sim2 <- n_controls - n_controls_nexp_sim2

  # we take the ones with some 0 but non negative
  idx_some_zero <- which(
    (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
      (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
  )
  n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
  n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
  n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
  n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

  #########
  n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)
  n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
  n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)
  n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)

  some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
  # var(log RR) = 1/a - 1/n1 + 1/c - 1/n2 (arm totals SUBTRACTED); mirrors the fixed
  # live copy .metaumbrella_rr_se_to_or in R/internal_multiple_formulas.R. This helper
  # is currently unused, but is kept in sync so a future load-order change cannot
  # silently reintroduce the added-sign bug. (Also drops the old unconditional
  # var_sim overwrite that had discarded the continuity-corrected branch.)
  var_sim <- ifelse(some_zero,
    1 / ((n_cases + 1) - (n_cases_nexp_sim + 0.5)) - 1 / ((n_cases + 1) + (n_controls + 1) - ((n_cases_nexp_sim + 0.5) + (n_controls_nexp_sim + 0.5))) +
      1 / (n_cases_nexp_sim + 0.5) - 1 / ((n_cases_nexp_sim + 0.5) + (n_controls_nexp_sim + 0.5)),
    1 / (n_cases - n_cases_nexp_sim) - 1 / (n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim)) +
      1 / n_cases_nexp_sim - 1 / (n_cases_nexp_sim + n_controls_nexp_sim)
  )

  best <- order((var_sim - var)^2)[1]
  n_cases_nexp <- n_cases_nexp_sim[best]
  n_controls_nexp <- n_controls_nexp_sim[best]
  n_cases_exp <- n_cases - n_cases_nexp
  n_controls_exp <- n_controls - n_controls_nexp
  n_exp <- n_cases_exp + n_controls_exp
  n_nexp <- n_cases_nexp + n_controls_nexp

  return(data.frame(n_cases_exp, n_exp, n_cases_nexp, n_nexp))
}
