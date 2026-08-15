library(MonteCarlo); library(tidyverse); library(ggplot2)
.estimate_n_from_or_and_n_cases = function (or, var, n_cases, n_controls) {

  res = data.frame(n_cases_exp = NA, n_cases_nexp = NA, n_controls_exp = NA, n_controls_nexp = NA)

  if (!is.na(or) & !is.na(var) & !is.na(n_cases) & !is.na(n_controls)) {

    # Create all possibilites of n
    n_cases_nexp_sim1 = 0:n_cases
    n_controls_nexp_sim1 = round(n_controls * (1 - (n_cases - n_cases_nexp_sim1) / (n_cases + (or - 1) * n_cases_nexp_sim1)))
    n_cases_exp_sim1 = n_cases - n_cases_nexp_sim1
    n_controls_exp_sim1 = n_controls - n_controls_nexp_sim1
    # sim1: possiblities with strictly positive n
    idx_non_zero <- which(
      n_cases_nexp_sim1 > 0 &
        n_controls_nexp_sim1 > 0 &
        n_cases_exp_sim1 > 0 &
        n_controls_exp_sim1 > 0
    )
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # sim2: possiblities with positive n and at least one zero (Add 0.5 to the possiblities with any 0)
    n_cases_nexp_sim2 = 0:n_cases
    n_controls_nexp_sim2 = round((n_controls + 0.5) - (n_controls + 1) * (n_cases - n_cases_nexp_sim2 + 0.5) / ((n_cases_nexp_sim2 + 0.5) * or + n_cases - n_cases_nexp_sim2 + 0.5))
    n_cases_exp_sim2 = n_cases - n_cases_nexp_sim2
    n_controls_exp_sim2 = n_controls - n_controls_nexp_sim2
    # (n_cases_exp_sim2 + 0.5) / (n_cases_nexp_sim2 + 0.5) / (n_controls_exp_sim2 + 0.5) * (n_controls_nexp_sim2 + 0.5)
    # select the ones with some 0 but non-negative
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    # join both previous vectors
    n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)
    n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)
    n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)

    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0

    var_sim <- ifelse(some_zero,
                      1 / ((n_cases+ 1 ) - (n_cases_nexp_sim + 0.5)) + 1 / ((n_controls + 1) - (n_controls_nexp_sim + 0.5)) + 1 / (n_cases_nexp_sim + 0.5) + 1 / (n_controls_nexp_sim + 0.5),
                      1 / (n_cases - n_cases_nexp_sim) + 1 / (n_controls - n_controls_nexp_sim) + 1 / n_cases_nexp_sim + 1 / n_controls_nexp_sim)

    #var_sim2 = 1 / ((n_cases+1) - (n_cases_nexp_sim+0.5)) + 1 / ((n_controls+1) - (n_controls_nexp_sim+0.5)) + 1 / (n_cases_nexp_sim+0.5) + 1 / (n_controls_nexp_sim+0.5)

    best = order((var_sim - var)^2)[1]

    res$n_cases_nexp = n_cases_nexp_sim[best]
    res$n_controls_nexp = n_controls_nexp_sim[best]
    res$n_cases_exp = n_cases - res$n_cases_nexp
    res$n_controls_exp = n_controls - res$n_controls_nexp

  }

  return(res)
}

#' Estimate the n, using the variance, the number of exposed and non-exposed subjects
#'
#' @param or OR
#' @param var variance
#' @param n_exp number of exposed participants
#' @param n_nexp number of non exposed participants
#'
#' @noRd
.estimate_n_from_or_and_n_exp = function (or, var, n_exp, n_nexp) {

  res = data.frame(n_cases_exp = NA, n_cases_nexp = NA, n_controls_exp = NA, n_controls_nexp = NA)

  if (!is.na(or) & !is.na(var) & !is.na(n_exp) & !is.na(n_nexp)) {
    # first: uncorrected values with 0
    n_controls_exp_sim1 = 0:n_exp
    n_controls_nexp_sim1 = round(n_nexp / (1 + (n_exp - n_controls_exp_sim1) / (or * n_controls_exp_sim1)))
    n_cases_exp_sim1 = n_exp - n_controls_exp_sim1
    n_cases_nexp_sim1 = n_nexp - n_controls_nexp_sim1
    # we take the ones without 0 and non-negative
    idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # correcting by 0.5
    n_controls_exp_sim2 = 0:n_exp
    n_controls_nexp_sim2 = round((n_nexp + 0.5) - ((n_nexp + 1)*(n_exp - n_controls_exp_sim2 + 0.5)) / ((n_controls_exp_sim2 + 0.5) * or + n_exp - n_controls_exp_sim2 + 0.5 ))
    n_cases_exp_sim2 = n_exp - n_controls_exp_sim2
    n_cases_nexp_sim2 = n_nexp - n_controls_nexp_sim2

    #SELECT THE ONES THAT HAS SOME 0 BUT NO NEGATIVE ONES
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    n_controls_exp_sim = append(n_controls_exp_sim1, n_controls_exp_sim2)
    n_controls_nexp_sim = append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_exp_sim = append(n_cases_exp_sim1, n_cases_exp_sim2)
    n_cases_nexp_sim = append(n_cases_nexp_sim1, n_cases_nexp_sim2)


    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
    var_sim <- ifelse(some_zero,
                      1 / ((n_exp+1) - (n_controls_exp_sim+0.5) + 1/(n_controls_exp_sim+0.5) + 1/((n_nexp+1) - (n_controls_nexp_sim+0.5)) + 1/(n_controls_nexp_sim+0.5)),
                      1 / (n_exp - n_controls_exp_sim) + 1 / n_controls_exp_sim + 1 / (n_nexp - n_controls_nexp_sim) + 1 / n_controls_nexp_sim
    )

    #var_sim = 1 / (n_exp - n_controls_exp_sim) + 1 / n_controls_exp_sim + 1 / (n_nexp - n_controls_nexp_sim) + 1 / n_controls_nexp_sim
    best = order((var_sim - var)^2)[1]
    res$n_controls_exp = n_controls_exp_sim[best]
    res$n_controls_nexp = n_controls_nexp_sim[best]
    res$n_cases_exp = n_exp - res$n_controls_exp
    res$n_cases_nexp = n_nexp - res$n_controls_nexp
  }
  return(res)
}

########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
OR_to_RR <- function(rr, br, n, p, br_guess) {

  n_exp = sum(rbinom(n, 1, p))
  n_exp = ifelse(n_exp == 0, sum(rbinom(n, 1, p)), n_exp)

  # n_exp = p * n        # a + b
  n_cases_exp = (rr * br) * n_exp
  n_controls_exp = n_exp - n_cases_exp
  n_nexp = (1 - p) * n  # c + d
  n_cases_nexp = br * n_nexp
  n_controls_nexp = n_nexp - n_cases_nexp
  n_cases = n_cases_exp + n_cases_nexp
  n_controls = n_controls_exp + n_controls_nexp

  # ==============================================================
  # Estimate OR/RR from raw data =================================
  # ==============================================================

  zero = which(n_cases_exp == 0 | n_cases_nexp == 0 | n_controls_exp == 0 | n_controls_nexp == 0)

  n_cases_exp[zero] = n_cases_exp[zero] + 0.5
  n_cases_nexp[zero] = n_cases_nexp[zero] + 0.5
  n_controls_exp[zero] = n_controls_exp[zero] + 0.5
  n_controls_nexp[zero] = n_controls_nexp[zero] + 0.5

  or_raw_data = suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp))
  or_se_raw_data = suppressWarnings(sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp))
  or_pval_raw_data = 1 - 2 * abs(pnorm(log(or_raw_data)/or_se_raw_data) - 0.5)

  rr_raw_data = suppressWarnings((n_cases_exp / n_exp) / (n_cases_nexp / n_nexp))
  rr_se_raw_data = suppressWarnings(sqrt(1 / n_cases_exp - 1 / n_exp + 1 / n_cases_nexp - 1 / n_nexp))
  rr_pval_raw_data = 1 - 2 * abs(pnorm(log(rr_raw_data)/rr_se_raw_data) - 0.5)

  br_raw_data = n_cases_nexp / n_nexp

  # =========================================================================
  # 1. Convert OR to RR according to Grant_2x2 ==============================
  # =========================================================================

  logrr_grant_2x2 = suppressWarnings(log(or_raw_data / (1 - br_guess + (br_guess * or_raw_data))))

  a_est = exp(logrr_grant_2x2) * br_guess * n_exp
  b_est = (1 - exp(logrr_grant_2x2) * br_guess) * n_exp
  c_est = br_guess * n_nexp
  d_est = (1 - br_guess) * n_nexp

  logrr_se_grant_2x2 = suppressWarnings(sqrt(1/a_est - 1/n_exp + 1/c_est - 1/n_nexp))
  logrr_pval_grant_2x2 = 1 - 2 * abs(pnorm(log(rr_raw_data)/rr_se_raw_data) - 0.5)
  logrr_ci_lo_grant_2x2 = logrr_grant_2x2 - qnorm(.975) * logrr_se_grant_2x2
  logrr_ci_up_grant_2x2 = logrr_grant_2x2 + qnorm(.975) * logrr_se_grant_2x2

  # =========================================================================
  # 2. Convert OR to RR according to Grant_delta ============================
  # =========================================================================
  logrr_grant_delta = logrr_grant_2x2

  logrr_se_grant_delta = sqrt(do.call(numDeriv::grad,
      list(x = log(or_raw_data),
           br_guess=br_guess,
           func = function(or_raw_data, br_guess) {
             exp(or_raw_data) / (1 - br_guess + (br_guess * exp(or_raw_data)))
             }))^2 * or_se_raw_data^2)
  logrr_pval_grant_delta = 1 - 2 * abs(pnorm(logrr_grant_delta/logrr_se_grant_delta) - 0.5)
  logrr_ci_lo_grant_delta = logrr_grant_delta - qnorm(.975) * logrr_se_grant_delta
  logrr_ci_up_grant_delta = logrr_grant_delta + qnorm(.975) * logrr_se_grant_delta

  # =========================================================================
  # 3. Convert OR to RR according to Grant_CI ===============================
  # =========================================================================
  logrr_grant_CI = logrr_grant_2x2
  or_ci_lo = exp(log(or_raw_data) - qnorm(.975) * or_se_raw_data)
  or_ci_up = exp(log(or_raw_data) + qnorm(.975) * or_se_raw_data)

  logrr_ci_lo_grant_CI = suppressWarnings(log(or_ci_lo / (1 - br_guess + (br_guess * or_ci_lo))))
  logrr_ci_up_grant_CI = suppressWarnings(log(or_ci_up / (1 - br_guess + (br_guess * or_ci_up))))
  logrr_se_grant_CI = (logrr_ci_up_grant_CI - logrr_ci_lo_grant_CI)/(2 * qnorm(.975))
  logrr_pval_grant_CI = 1 - 2 * abs(pnorm(logrr_grant_CI/logrr_se_grant_CI) - 0.5)
  logrr_ci_lo_grant_CI = logrr_grant_CI - qnorm(.975) * logrr_se_grant_CI
  logrr_ci_up_grant_CI = logrr_grant_CI + qnorm(.975) * logrr_se_grant_CI

  # =========================================================================
  # 4. Convert OR to RR according to us (cases) =============================
  # =========================================================================
  contingency_meta_cases = .estimate_n_from_or_and_n_cases(
    or = or_raw_data, var = or_se_raw_data^2,
    n_cases = n_cases, n_controls = n_controls)

    logrr_meta_cases = log(with(contingency_meta_cases, suppressWarnings((n_cases_exp / (n_cases_exp + n_controls_exp)) / (n_cases_nexp / (n_cases_nexp + n_controls_nexp)))))
    logrr_se_meta_cases = with(contingency_meta_cases, suppressWarnings(sqrt(1 / n_cases_exp - 1 / (n_cases_exp + n_controls_exp) + 1 / n_cases_nexp - 1 / (n_cases_nexp + n_controls_nexp))))
    logrr_pval_meta_cases = 1 - 2 * abs(pnorm(logrr_meta_cases/logrr_se_meta_cases) - 0.5)
    logrr_ci_lo_meta_cases = logrr_meta_cases - qnorm(.975) * logrr_se_meta_cases
    logrr_ci_up_meta_cases = logrr_meta_cases + qnorm(.975) * logrr_se_meta_cases


    # =========================================================================
    # 5. Convert OR to RR according to us (exp) ===============================
    # =========================================================================
    contingency_meta_exp = .estimate_n_from_or_and_n_exp(
      or = or_raw_data, var = or_se_raw_data^2, n_exp = n_exp, n_nexp = n_nexp)

    logrr_meta_exp = log(with(contingency_meta_exp, suppressWarnings((n_cases_exp / (n_cases_exp + n_controls_exp)) / (n_cases_nexp / (n_cases_nexp + n_controls_nexp)))))
    logrr_se_meta_exp = with(contingency_meta_exp, suppressWarnings(sqrt(1 / n_cases_exp - 1 / (n_cases_exp + n_controls_exp) + 1 / n_cases_nexp - 1 / (n_cases_nexp + n_controls_nexp))))
    logrr_pval_meta_exp = 1 - 2 * abs(pnorm(logrr_meta_exp/logrr_se_meta_exp) - 0.5)
    logrr_ci_lo_meta_exp = logrr_meta_exp - qnorm(.975) * logrr_se_meta_exp
    logrr_ci_up_meta_exp = logrr_meta_exp + qnorm(.975) * logrr_se_meta_exp

    # =========================================================================
    # 6. Convert OR to RR - transpose =========================================
    # =========================================================================
    logrr_transpose = log(or_raw_data)
    logrr_se_transpose = or_se_raw_data
    logrr_pval_transpose = 1 - 2 * abs(pnorm(logrr_transpose/logrr_se_transpose) - 0.5)
    logrr_ci_lo_transpose = logrr_transpose - qnorm(.975) * logrr_se_transpose
    logrr_ci_up_transpose = logrr_transpose + qnorm(.975) * logrr_se_transpose

    # =========================================================================
    # 7. Convert OR to RR - dipietrantonj =========================================
    # =========================================================================

    n_dec = max(
      nchar(gsub("^.+[.]", "", or_raw_data)),
      nchar(gsub("^.+[.]", "", or_ci_lo)),
      nchar(gsub("^.+[.]", "", or_ci_up))
    )

    estim = estimraw::estim_raw(
      es=or_raw_data, lb=or_ci_lo, ub=or_ci_up,
      m1=n_exp, m2=n_nexp, dec = n_dec, measure = "or")

    if (length(estim) != 4) {
      logrr_dipietrantonj = suppressWarnings(log((estim[[1]]$a[1] / (estim[[1]]$a[1] + estim[[1]]$b[1])) / (estim[[1]]$c[1] / (estim[[1]]$c[1] + estim[[1]]$d[1]))))
      logrr_se_dipietrantonj = suppressWarnings(sqrt(1 / estim[[1]]$a[1] - 1 / (estim[[1]]$a[1] + estim[[1]]$b[1]) + 1 / estim[[1]]$c[1] - 1 / (estim[[1]]$c[1] + estim[[1]]$d[1])))

    } else {
      logrr_dipietrantonj = log(suppressWarnings((estim$a[1] / (estim$a[1] + estim$b[1])) / (estim$c[1] / (estim$c[1] + estim$d[1]))))
      logrr_se_dipietrantonj = suppressWarnings(sqrt(1 / estim$a[1] - 1 / (estim$a[1] + estim$b[1]) + 1 / estim$c[1] - 1 / (estim$c[1] + estim$d[1])))
    }
    logrr_pval_dipietrantonj = 1 - 2 * abs(pnorm(logrr_dipietrantonj/logrr_se_dipietrantonj) - 0.5)
    logrr_ci_lo_dipietrantonj = logrr_dipietrantonj - qnorm(.975) * logrr_se_dipietrantonj
    logrr_ci_up_dipietrantonj = logrr_dipietrantonj + qnorm(.975) * logrr_se_dipietrantonj

  return(list(
    "or_raw_data" = log(or_raw_data),
    "br_raw_data" = br_raw_data,
    "rr_raw_data" = log(rr_raw_data),
    "rr_se_raw" = rr_se_raw_data,
    "rr_pval_raw" = rr_pval_raw_data,
    # metaumbrella cases
    "rr_metaumbrella_cases" = logrr_meta_cases,
    "rr_se_metaumbrella_cases" = logrr_se_meta_cases,
    "rr_pval_metaumbrella_cases" = logrr_pval_meta_cases,
    "rr_ci_lo_metaumbrella_cases" = logrr_ci_lo_meta_cases,
    "rr_ci_up_metaumbrella_cases" = logrr_ci_up_meta_cases,

    # metaumbrella exp
    "rr_metaumbrella_exp" = logrr_meta_exp,
    "rr_se_metaumbrella_exp" = logrr_se_meta_exp,
    "rr_pval_metaumbrella_exp" = logrr_pval_meta_exp,
    "rr_ci_lo_metaumbrella_exp" = logrr_ci_lo_meta_exp,
    "rr_ci_up_metaumbrella_exp" = logrr_ci_up_meta_exp,

    # grant delta
    "rr_grant_delta" = logrr_grant_delta,
    "rr_se_grant_delta" = logrr_se_grant_delta,
    "rr_pval_grant_delta" = logrr_pval_grant_delta,
    "rr_ci_lo_grant_delta" = logrr_ci_lo_grant_delta,
    "rr_ci_up_grant_delta" = logrr_ci_up_grant_delta,

    # grant 2x2
    "rr_grant_2x2" = logrr_grant_2x2,
    "rr_se_grant_2x2" = logrr_se_grant_2x2,
    "rr_pval_grant_2x2" = logrr_pval_grant_2x2,
    "rr_ci_lo_grant_2x2" = logrr_ci_lo_grant_2x2,
    "rr_ci_up_grant_2x2" = logrr_ci_up_grant_2x2,

    # grant CI
    "rr_grant_CI" = logrr_grant_CI,
    "rr_se_grant_CI" = logrr_se_grant_CI,
    "rr_pval_grant_CI" = logrr_pval_grant_CI,
    "rr_ci_lo_grant_CI" = logrr_ci_lo_grant_CI,
    "rr_ci_up_grant_CI" = logrr_ci_up_grant_CI,

    # transpose
    "rr_transpose" = logrr_transpose,
    "rr_se_transpose" = logrr_se_transpose,
    "rr_pval_transpose" = logrr_pval_transpose,
    "rr_ci_lo_transpose" = logrr_ci_lo_transpose,
    "rr_ci_up_transpose" = logrr_ci_up_transpose,

    # dipietrantonj
    "rr_dipietrantonj" = logrr_dipietrantonj,
    "rr_se_dipietrantonj" = logrr_se_dipietrantonj,
    "rr_pval_dipietrantonj" = logrr_pval_dipietrantonj,
    "rr_ci_lo_dipietrantonj" = logrr_ci_lo_dipietrantonj,
    "rr_ci_up_dipietrantonj" = logrr_ci_up_dipietrantonj
  ))
}

rr_grid <- c(0.25, 0.5, 0.75, 1)
n_grid <- c(25, 50, 75, 100, 300)
br_grid <- c(0.1, 0.25, 0.5)
br_guess_grid <- c(0.1, 0.25, 0.5)
p_grid <- c(0.3, 0.5, 0.7)

param_list = list("rr" = rr_grid, "n" = n_grid, "br" = br_grid, "p" = p_grid, "br_guess" = br_guess_grid)
MC_result <- MonteCarlo(func = OR_to_RR, nrep = 50,
                        max_grid = 10000, param_list = param_list)
data_sim <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data")
rio::export(data_sim, "OR_to_RR_sim.txt")



res_or = read.delim("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/OR_to_RR_sim.txt") %>%
  mutate_if(is.character, as.numeric)

res_or[res_or==Inf|res_or==-Inf] <- NA
mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}
var_na = function(x) {
  res = var(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}

checks = res_or %>%
  group_by(rr, br) %>%
  summarise(rr_raw = exp(mean_na(rr_raw_data)),
            br_raw = mean_na(br_raw_data))

method = c("metaumbrella_cases", "metaumbrella_exp", "grant_CI", "grant_2x2", "grant_delta", "dipietrantonj", "transpose")

for (met in method) {
  print(met)
  res_or[, paste0("rr_var_", met)] <- res_or[, paste0("rr_se_", met)]^2
  res_or[, paste0("ci_cov_rr_", met)] <- res_or[, paste0("rr_ci_lo_", met)] <= res_or$rr_raw_data & res_or[, paste0("rr_ci_up_", met)] >= res_or$rr_raw_data
}
res = res_or %>%
  group_by(rr, n, br, br_guess) %>%
  summarise(
    n_sim = n(),

    var_rr_metaumbrella_cases = var_na(rr_metaumbrella_cases),
    rr_metaumbrella_cases = mean_na(rr_metaumbrella_cases),
    mean_rr_var_metaumbrella_cases = mean_na(rr_var_metaumbrella_cases),
    bias_ci_rr_metaumbrella_cases = mean_na(ci_cov_rr_metaumbrella_cases),

    var_rr_metaumbrella_exp = var_na(rr_metaumbrella_exp),
    rr_metaumbrella_exp = mean_na(rr_metaumbrella_exp),
    mean_rr_var_metaumbrella_exp = mean_na(rr_var_metaumbrella_exp),
    bias_ci_rr_metaumbrella_exp = mean_na(ci_cov_rr_metaumbrella_exp),

    var_rr_dipietrantonj = var_na(rr_dipietrantonj),
    rr_dipietrantonj = mean_na(rr_dipietrantonj),
    mean_rr_var_dipietrantonj = mean_na(rr_var_dipietrantonj),
    bias_ci_rr_dipietrantonj = mean_na(ci_cov_rr_dipietrantonj),

    var_rr_transpose = var_na(rr_transpose),
    rr_transpose = mean_na(rr_transpose),
    mean_rr_var_transpose = mean_na(rr_var_transpose),
    bias_ci_rr_transpose = mean_na(ci_cov_rr_transpose),

    var_rr_grant_CI = var_na(rr_grant_CI),
    rr_grant_CI = mean_na(rr_grant_CI),
    mean_rr_var_grant_CI = mean_na(rr_var_grant_CI),
    bias_ci_rr_grant_CI = mean_na(ci_cov_rr_grant_CI),

    var_rr_grant_delta = var_na(rr_grant_delta),
    rr_grant_delta = mean_na(rr_grant_delta),
    mean_rr_var_grant_delta = mean_na(rr_var_grant_delta),
    bias_ci_rr_grant_delta = mean_na(ci_cov_rr_grant_delta),

    var_rr_grant_2x2 = var_na(rr_grant_2x2),
    rr_grant_2x2 = mean_na(rr_grant_2x2),
    mean_rr_var_grant_2x2 = mean_na(rr_var_grant_2x2),
    bias_ci_rr_grant_2x2 = mean_na(ci_cov_rr_grant_2x2)
  )

res = res_or %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("acc_es") |
                 starts_with("acc_se") |
                 starts_with("acc_pval"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(rr, n, br, br_guess, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("acc_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("acc_se", method_long, fixed = TRUE) ~ "_se",
                           grepl("acc_pval", method_long, fixed = TRUE) ~ "_pval"),
         method = case_when(grepl("metaumbrella_cases", method_long, fixed = TRUE) ~ "metaumbrella_cases",
                            grepl("metaumbrella_exp", method_long, fixed = TRUE) ~ "metaumbrella_exp",
                            grepl("grant_CI", method_long, fixed = TRUE) ~ "grant_CI",
                            grepl("grant_2x2", method_long, fixed = TRUE) ~ "grant_2x2",
                            grepl("grant_delta", method_long, fixed = TRUE) ~ "grant_delta",
                            grepl("dipietrantonj", method_long, fixed = TRUE) ~ "dipietrantonj",
                            grepl("transpose", method_long, fixed = TRUE) ~ "transpose")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")
err_plot = res
err_plot$rr <- paste0("RR=", err_plot$rr)
err_plot$br <- paste0("br=", err_plot$br)
err_plot$br_guess <- paste0("guessed br=", err_plot$br_guess)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "200"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=200"))

err_plot$method[err_plot$method == "dipietrantonj"] <- "dipietrantonj"
err_plot$method[err_plot$method == "transpose"] <- "RR = OR"
err_plot$method[err_plot$method == "metaumbrella_cases"] <- "metaumbrella (cases)"
err_plot$method[err_plot$method == "metaumbrella_exp"] <- "metaumbrella (exp)"
err_plot$method[err_plot$method == "grant_CI"] <- "Grant (CI)"
err_plot$method[err_plot$method == "grant_2x2"] <- "Grant (2x2)"
err_plot$method[err_plot$method == "grant_delta"] <- "Grant (delta)"

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/OR_to_RR_AGG.txt")
