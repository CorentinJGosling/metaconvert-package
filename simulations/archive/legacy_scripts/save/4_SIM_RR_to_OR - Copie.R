library(MonteCarlo); library(tidyverse); library(ggplot2)
.metaumbrella_rr_se_to_or <- function(rr, logrr, logrr_se, n_cases, n_controls) {

  es = data.frame(value = NA, se = NA)

  if (!is.na(rr) & !is.na(logrr_se) & !is.na(n_cases) & !is.na(n_controls)) {

    # uncorrected
    n_cases_nexp_sim1 = 0:n_cases
    n_controls_nexp_sim1 = round(n_cases_nexp_sim1 * ((rr * (n_cases + n_controls)) / (n_cases + (rr - 1) * n_cases_nexp_sim1) - 1))
    n_cases_exp_sim1 = n_cases - n_cases_nexp_sim1
    n_controls_exp_sim1 = n_controls - n_controls_nexp_sim1

    # we take only positives (no-zero)
    idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # corregint 0.5
    n_cases_nexp_sim2 = 0:n_cases
    n_controls_nexp_sim2 = ((n_cases + n_controls - n_cases_nexp_sim2 + 1) - (n_cases + n_controls + 2) * (n_cases - n_cases_nexp_sim2 + 0.5) / ((n_cases_nexp_sim2 + 0.5) * rr + n_cases - n_cases_nexp_sim2 + 0.5))
    n_cases_exp_sim2 = n_cases - n_cases_nexp_sim2
    n_controls_exp_sim2 = n_controls - n_controls_nexp_sim2

    # we take the ones with some 0 but non negative
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    #
    n_controls_exp_sim = append(n_controls_exp_sim1, n_controls_exp_sim2)
    n_controls_nexp_sim = append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_nexp_sim = append(n_cases_nexp_sim1, n_cases_nexp_sim2)
    n_cases_exp_sim = append(n_cases_exp_sim1, n_cases_exp_sim2)

    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
    var_sim = ifelse(some_zero,
                     1 / ((n_cases+1) - (n_cases_nexp_sim+0.5)) + 1 / ((n_cases+1) + (n_controls+1) - ((n_cases_nexp_sim+0.5) + (n_controls_nexp_sim+0.5))) +
                       1 / (n_cases_nexp_sim+0.5) + 1 / ((n_cases_nexp_sim+0.5) + (n_controls_nexp_sim+0.5)),
                     1 / (n_cases - n_cases_nexp_sim) + 1 / (n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim)) +
                       1 / n_cases_nexp_sim + 1 / (n_cases_nexp_sim + n_controls_nexp_sim)
    )

    var_sim = 1 / (n_cases - n_cases_nexp_sim) + 1 / (n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim)) +
      1 / n_cases_nexp_sim + 1 / (n_cases_nexp_sim + n_controls_nexp_sim)

    best = order((var_sim - logrr_se^2)^2)[1]
    n_cases_nexp = n_cases_nexp_sim[best]
    n_controls_nexp = n_controls_nexp_sim[best]
    n_cases_exp = n_cases - n_cases_nexp
    n_controls_exp = n_controls - n_controls_nexp


    # es$n_cases_nexp = n_cases_nexp
    # es$n_controls_nexp = n_controls_nexp
    # es$n_cases_exp = n_cases_exp
    # es$n_controls_exp = n_controls_exp

    es$value = (n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp)
    es$se = sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp)
  }
  return(es)
}


# or = 0.25; br = 0.1; br_guess=0.25; n=25; p=0.5
########################################################################################
RR_to_OR <- function(or, br, br_guess, n, p) {
  # n1 = round(p * n)
  n1 = sum(rbinom(n, 1, p))
  n1 = ifelse(n1 == 0, sum(rbinom(n, 1, p)), n1)

  n2 = n - n1
  n_cases_nexp = round(br * n2)
  n_cases_nexp = ifelse(n_cases_nexp < 0, 0, n_cases_nexp)
  n_controls_nexp = round(n2 - n_cases_nexp)
  n_controls_nexp = ifelse(n_controls_nexp < 0, 0, n_controls_nexp)
  rr = n2 / (n_cases_nexp + n_controls_nexp / or)
  n_cases_exp = round((rr * br) * n1)
  n_cases_exp = ifelse(n_cases_exp < 0, 0, n_cases_exp)
  n_controls_exp = round(n1 - n_cases_exp)
  n_controls_exp = ifelse(n_controls_exp < 0, 0, n_controls_exp)
  n_cases = n_cases_exp + n_cases_nexp
  n_controls = n_controls_exp + n_controls_nexp
  n_exp = n_cases_exp + n_controls_exp
  n_nexp = n_cases_nexp + n_controls_nexp

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
  rr_ci_lo_raw_data = exp(log(rr_raw_data) - qnorm(.975) * rr_se_raw_data)
  rr_ci_up_raw_data = exp(log(rr_raw_data) + qnorm(.975) * rr_se_raw_data)

  br_raw_data = n_cases_nexp / n_nexp

  # =========================================================================
  # 1. Convert RR to OR according to Grant_2x2 ==============================
  # =========================================================================
  logor_grant_2x2 = suppressWarnings(log(rr_raw_data * (1 - br_guess) / (1 - rr_raw_data * br_guess)))
  a_est = rr_raw_data * br_guess * n_exp
  b_est = (1 - rr_raw_data * br_guess) * n_exp
  c_est = br_guess * n_nexp
  d_est = (1 - br_guess) * n_nexp
  logor_se_grant_2x2 = suppressWarnings(sqrt(1 / a_est + 1 / b_est + 1 / c_est + 1 / d_est))
  logor_pval_grant_2x2 = 1 - 2 * abs(pnorm(logor_grant_2x2/logor_se_grant_2x2) - 0.5)

  # =========================================================================
  # 2. Convert RR to OR according to Grant_delta ============================
  # =========================================================================
  logor_grant_delta = logor_grant_2x2
  logor_se_grant_delta = sqrt(do.call(numDeriv::grad, list(
    x = log(rr_raw_data),
    baseline_risk=br_guess,
    func = function(rr, baseline_risk) {
      exp(rr) * (1 - baseline_risk) / (1 - rr * baseline_risk)
      }))^2 * rr_se_raw_data^2)
  logor_pval_grant_delta = 1 - 2 * abs(pnorm(logor_grant_delta/logor_se_grant_delta) - 0.5)

  # =========================================================================
  # 3. Convert RR to OR according to Grant_CI ===============================
  # =========================================================================
  logor_grant_CI = logor_grant_2x2
  logor_ci_lo_grant = suppressWarnings(log(rr_ci_lo_raw_data * (1 - br_guess) / (1 - rr_ci_lo_raw_data * br_guess)))
  logor_ci_up_grant = suppressWarnings(log(rr_ci_up_raw_data * (1 - br_guess) / (1 - rr_ci_up_raw_data * br_guess)))
  logor_se_grant_CI = (logor_ci_up_grant - logor_ci_lo_grant)/(2 * qnorm(.975))
  logor_pval_grant_CI = 1 - 2 * abs(pnorm(logor_grant_CI/logor_se_grant_CI) - 0.5)

  # =========================================================================
  # 4. Convert RR to OR - transpose =========================================
  # =========================================================================
  logor_transpose = log(rr_raw_data)
  logor_se_transpose = rr_se_raw_data
  logor_pval_transpose = 1 - 2 * abs(pnorm(logor_transpose / logor_se_transpose) - 0.5)

  # =========================================================================
  # 5. Convert RR to OR according to us (cases) =============================
  # =========================================================================
  raw_res = .metaumbrella_rr_se_to_or(rr = rr_raw_data, logrr_se = rr_se_raw_data,
                                      n_cases = n_cases, n_controls = n_controls)
  logor_metaumbrella = log(raw_res$value)
  logor_se_metaumbrella = raw_res$se
  logor_pval_metaumbrella = 1 - 2 * abs(pnorm(logor_metaumbrella / logor_se_metaumbrella) - 0.5)

  # =========================================================================
  # 6. Convert RR to OR according to us (cases) =============================
  # =========================================================================
  n_dec = max(
    nchar(gsub("^.+[.]", "", rr_raw_data)),
    nchar(gsub("^.+[.]", "", rr_ci_lo_raw_data)),
    nchar(gsub("^.+[.]", "", rr_ci_up_raw_data))
  )

  estim = estimraw::estim_raw(
    es=rr_raw_data, lb=rr_ci_lo_raw_data, ub=rr_ci_up_raw_data,
    m1=n_exp, m2=n_nexp, dec = n_dec, measure = "rr")

  if (length(estim) != 4) {
    logor_dipietrantonj = suppressWarnings(
      suppressWarnings(log((estim[[1]]$a[1] * estim[[1]]$d[1]) / (estim[[1]]$b[1] * estim[[1]]$c[1])))
      )
    logor_se_dipietrantonj = suppressWarnings(sqrt(1/estim[[1]]$a[1] + 1/estim[[1]]$b[1] + 1/estim[[1]]$c[1] + 1/estim[[1]]$d[1]))

  } else {
    logor_dipietrantonj = suppressWarnings(log((estim$a[1] * estim$d[1]) / (estim$b[1] * estim$c[1])))
    logor_se_dipietrantonj = suppressWarnings(sqrt(1/estim$a[1] + 1/estim$b[1] + 1/estim$c[1] + 1/estim$d[1]))
  }
  logor_pval_dipietrantonj = 1 - 2 * abs(pnorm(logor_dipietrantonj/logor_se_dipietrantonj) - 0.5)

  # return result:
  return(list(
    "rr_raw_data" = log(rr_raw_data),
    "br_raw_data" = br_raw_data,
    "or_raw_data" = log(or_raw_data),
    "or_se_raw_data" = or_se_raw_data,
    "or_pval_raw_data" = or_pval_raw_data,
    # metaumbella
    "or_metaumbrella" = logor_metaumbrella,
    "or_se_metaumbrella" = logor_se_metaumbrella,
    "or_pval_metaumbrella" = logor_pval_metaumbrella,
    # grant delta
    "or_grant_delta" = logor_grant_delta,
    "or_se_grant_delta" = logor_se_grant_delta,
    "or_pval_grant_delta" = logor_pval_grant_delta,
    # grant 2x2
    "or_grant_2x2" = logor_grant_2x2,
    "or_se_grant_2x2" = logor_se_grant_2x2,
    "or_pval_grant_2x2" = logor_pval_grant_2x2,
    # grant CI
    "or_grant_CI" = logor_grant_CI,
    "or_se_grant_CI" = logor_se_grant_CI,
    "or_pval_grant_CI" = logor_pval_grant_CI,
    # transpose
    "or_transpose" = logor_transpose,
    "or_se_transpose" = logor_se_transpose,
    "or_pval_transpose" = logor_pval_transpose,
    # dipietrantonj
    "or_dipietrantonj" = logor_dipietrantonj,
    "or_se_dipietrantonj" = logor_se_dipietrantonj,
    "or_pval_dipietrantonj" = logor_pval_dipietrantonj

  ))
}
or_grid <- c(0.25, 0.5, 0.75, 1)
n_grid <- c(25, 50, 75, 100, 200)
br_grid <- c(0.1, 0.25, 0.5)
br_guess_grid <- c(0.1, 0.25, 0.5)
p_grid <- c(0.3, 0.5, 0.7)
# collect parameter grids in list:
param_list = list("or" = or_grid, "n" = n_grid, "br" = br_grid,
                  "br_guess" = br_guess_grid, "p" = p_grid)
MC_result <- MonteCarlo(func = RR_to_OR, nrep = 500, param_list = param_list, max_grid = 10000)
df <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data")
rio::export(df, "RR_to_OR_sim.txt")

res_rr = read.delim("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/RR_to_OR_sim.txt") %>%
  mutate_if(is.character, as.numeric)

res_rr[res_rr==Inf|res_rr==-Inf] <- NA
mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}

checks = res_rr %>%
  group_by(or, br) %>%
  summarise(or_raw = exp(mean_na(or_raw_data)),
            br_raw = mean_na(br_raw_data))

method = c("metaumbrella", "grant_CI", "grant_2x2", "grant_delta", "dipietrantonj", "transpose")
for (met in method) {
  print(met)
  res_rr[, paste0("acc_es_", met)] <- NA
  res_rr[, paste0("acc_se_", met)] <- NA
  res_rr[, paste0("acc_pval_", met)] <- NA
  res_rr[, paste0("acc_es_", met)] <- abs(res_rr[, paste0("or_", met)] - res_rr$or_raw_data)
  res_rr[, paste0("acc_se_", met)] <- abs(res_rr[, paste0("or_se_", met)] - res_rr$or_se_raw_data)
  res_rr[, paste0("acc_pval_", met)] <- abs(res_rr[, paste0("or_pval_", met)] - res_rr$or_pval_raw_data)
}
res = res_rr %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("acc_es") |
                 starts_with("acc_se") |
                 starts_with("acc_pval"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(or, n, br, br_guess, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("acc_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("acc_se", method_long, fixed = TRUE) ~ "_se",
                           grepl("acc_pval", method_long, fixed = TRUE) ~ "_pval"),
         method = case_when(grepl("metaumbrella", method_long, fixed = TRUE) ~ "metaumbrella",
                            grepl("grant_CI", method_long, fixed = TRUE) ~ "grant_CI",
                            grepl("grant_2x2", method_long, fixed = TRUE) ~ "grant_2x2",
                            grepl("grant_delta", method_long, fixed = TRUE) ~ "grant_delta",
                            grepl("dipietrantonj", method_long, fixed = TRUE) ~ "dipietrantonj",
                            grepl("transpose", method_long, fixed = TRUE) ~ "transpose")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")
err_plot = res
err_plot$rr <- paste0("OR=", err_plot$or)
err_plot$br <- paste0("br=", err_plot$br)
err_plot$br_guess <- paste0("guessed br=", err_plot$br_guess)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "200"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=200"))

err_plot$method[err_plot$method == "dipietrantonj"] <- "dipietrantonj"
err_plot$method[err_plot$method == "transpose"] <- "RR = OR"
err_plot$method[err_plot$method == "metaumbrella"] <- "metaumbrella"
err_plot$method[err_plot$method == "grant_CI"] <- "Grant (CI)"
err_plot$method[err_plot$method == "grant_2x2"] <- "Grant (2x2)"
err_plot$method[err_plot$method == "grant_delta"] <- "Grant (delta)"

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/RR_to_OR_AGG.txt")

