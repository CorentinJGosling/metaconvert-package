library(MonteCarlo); library(tidyverse); library(ggplot2)
.se_from_or_sim = function (or, n_cases, n_controls) {

    ca_ex = 1:(n_cases - 1)
    ca_ne = n_cases - ca_ex
    co_ex = round(n_controls / (1 + ca_ne * or / ca_ex))
    co_ex[which(co_ex < 1 | co_ex > n_controls - 1)] = NA
    co_ne = n_controls - co_ex
    v_or_mean = mean(1/ca_ex + 1/ca_ne + 1/co_ex + 1/co_ne, na.rm = TRUE)

    res = data.frame(
      value = rep(NA, length(or)),
      var = rep(NA, length(or)),
      se = rep(NA, length(or))
    )
    res$value = or
    res$var = v_or_mean
    res$se = sqrt(v_or_mean)


  return(res)
}
# or=1.2;
# br=0.3;
# n=50;
# p=0.4;
########################################################################################
SE_OR <- function(or, br, n, p) {
  # n1 = round(p * n)
  # n2 = n - n1
  # c = br * n2
  # d = n2 - c
  # rr = n2 / (c + d / or)
  # a = (rr * br) * n1
  # b = n1 - a
  n1 = sum(rbinom(n, 1, p))
  n1 = ifelse(n1 == 0, sum(rbinom(n, 1, p)), n1)

  # n1 = round(p * n)
  n2 = n - n1
  n_cases_nexp = br * n2
  n_controls_nexp = n2 - n_cases_nexp
  rr = n2 / (n_cases_nexp + n_controls_nexp / or)
  n_cases_exp = (rr * br) * n1
  n_controls_exp = n1 - n_cases_exp
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
  pval_raw_data = 1 - 2 * abs(pnorm(log(or_raw_data)/or_se_raw_data) - 0.5)

  br_raw_data = n_cases_nexp / n_nexp


  # =========================================================================
  # 1. Obtain SE from OR + n_cases/exp ======================================
  # =========================================================================

  or_se_metaumbrella = .se_from_or_sim(or = or_raw_data, n_cases = n_cases,
                           n_controls = n_controls)$se
  pval_metaumbrella = 1 - 2 * abs(pnorm(log(or_raw_data) / or_se_metaumbrella) - 0.5)

  # return result:
  return(list(
    "or_raw_data" = or_raw_data,
    "br_raw_data" = br_raw_data,
    "or_se_raw_data" = or_se_raw_data,
    "pval_raw_data" = pval_raw_data,
    "or_metaumbrella" = or_raw_data,
    "or_se_metaumbrella" = or_se_metaumbrella,
    "pval_metaumbrella" = pval_metaumbrella
  ))
}

or_grid <- c(0.25, 0.5, 0.75, 1)
n_grid <- c(25, 50, 75, 100, 200)
br_grid <- c(0.1, 0.25, 0.5)
p_grid <- c(0.3, 0.5, 0.7)

param_list = list("or" = or_grid, "n" = n_grid, "br" = br_grid, "p" = p_grid)
MC_result <- MonteCarlo(func = SE_OR, nrep = 5000, param_list = param_list)
df <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/")
rio::export(df, "OR_SE_sim.txt")

res_seor = read.delim("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/OR_SE_sim.txt")
mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}

checks = res_seor %>%
  group_by(or, br) %>%
  summarise(or_raw_data = exp(mean_na(log(or_raw_data))),
            br_raw_data = mean_na(br_raw_data))

method = c("metaumbrella")
for (met in method) {
  res_seor[, paste0("acc_es", met)] <- abs(log(res_seor[, paste0("or_", met)]) -
                                               log(res_seor$or_raw_data))
  res_seor[, paste0("acc_se", met)] <- abs(log(res_seor[, paste0("or_", met)]) -
                                                     log(res_seor$or_se_raw_data))
  res_seor[, paste0("acc_pval", met)] <- abs(res_seor[, paste0("pval_", met)] -
                                             res_seor$pval_raw_data)
}
res = res_seor %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols =  starts_with("acc_"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(or, br, n, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("acc_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("acc_se", method_long, fixed = TRUE) ~ "_se",
                           grepl("acc_pval", method_long, fixed = TRUE) ~ "_pval"),
         method = "metaumbrella") %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot = res

err_plot$or <- paste0("OR=", err_plot$or)
err_plot$br <- paste0("br=", err_plot$br)
err_plot$n <- paste0("n=", err_plot$n)

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/OR_SE_AGG.txt")

