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

  n_exp = ifelse(runif(1) > p * n - floor(p * n), floor(p * n), ceiling(p * n))
  n_nexp = n - n_exp

  rr = 1 / ((1 - br) / or + br)

  n_cases_exp = sum(rbinom(1, n_exp, rr * br))
  n_controls_exp = n_exp - n_cases_exp

  n_cases_nexp = sum(rbinom(n_nexp, 1, br))
  n_controls_nexp = n_nexp - n_cases_nexp

  n_cases = n_cases_exp + n_cases_nexp
  n_controls = n_controls_exp + n_controls_nexp
  n_exp = n_cases_exp + n_controls_exp
  n_nexp = n_cases_nexp + n_controls_nexp
  # Estimate OR/RR from raw data =================================
  zero = which(n_cases_exp == 0 | n_cases_nexp == 0 | n_controls_exp == 0 | n_controls_nexp == 0)

  n_cases_exp[zero] = n_cases_exp[zero] + 0.5
  n_cases_nexp[zero] = n_cases_nexp[zero] + 0.5
  n_controls_exp[zero] = n_controls_exp[zero] + 0.5
  n_controls_nexp[zero] = n_controls_nexp[zero] + 0.5

  or_raw_data = suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp))
  or_se_raw_data = suppressWarnings(sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp))
  br_raw_data = n_cases_nexp / n_nexp


  # 1. Obtain SE from OR + n_cases/exp ======================================

  res = .se_from_or_sim(or = or_raw_data, n_cases = n_cases,
                        n_controls = n_controls)
  or_metaumbrella = log(res$value)
  or_se_metaumbrella = res$se
  or_ci_lo_metaumbrella = or_metaumbrella - qnorm(.975) * or_se_metaumbrella
  or_ci_up_metaumbrella = or_metaumbrella + qnorm(.975) * or_se_metaumbrella

  return(list(
    "or_raw_data" = log(or_raw_data),
    "br_raw_data" = br_raw_data,
    "or_se_raw_data" = or_se_raw_data,
    "or_metaumbrella" = or_metaumbrella,
    "or_se_metaumbrella" = or_se_metaumbrella,
    "or_ci_lo_metaumbrella" = or_ci_lo_metaumbrella,
    "or_ci_up_metaumbrella" = or_ci_up_metaumbrella
  ))
}

or_grid <- c(2, 1, 0.75, 0.5, 0.25)
n_grid <- c(25, 50, 75, 100, 300)
br_grid <- c(0.5, 0.15, 0.05, 0.01)
br_guess_grid <- c(0.5, 0.15, 0.05, 0.01)
p_grid <- c(0.3, 0.5, 0.7)

param_list = list("or" = or_grid, "n" = n_grid, "br" = br_grid, "p" = p_grid)
MC_result <- MonteCarlo(func = SE_OR, nrep = 10000,
                        time_n_test = TRUE,
                        param_list = param_list)
data_sim <- MakeFrame(MC_result)

rio::export(data_sim, "D:/simulations/data/OR_SE_sim10000.txt")

res_seor = read.delim("D:/simulations/data/OR_SE_sim10000.txt")
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

checks = res_seor %>%
  group_by(or) %>%
  summarise(or_raw_data = exp(mean_na(log(or_raw_data))),
            br_raw_data = mean_na(br_raw_data),
            n=n())

method = c("metaumbrella")
for (met in method) {

  res_seor[, paste0("bias_es_or_", met)] <- res_seor[, paste0("or_", met)] - res_seor$or_raw_data
  res_seor[, paste0("or_var_", met)] = res_seor[, paste0("or_se_", met)]^2

  res_seor[, paste0("ci_cov_or_", met)] = res_seor[, paste0("or_ci_lo_", met)] <= res_seor$or_raw_data &
    res_seor[, paste0("or_ci_up_", met)] >= res_seor$or_raw_data
}

res = res_seor %>%
  group_by(or, n, br) %>%
  summarise(
    n_sim = n(),
    bias_es_or_metaumbrella = mean_na(bias_es_or_metaumbrella),
    var_or_metaumbrella = var_na(or_metaumbrella),
    mean_or_var_metaumbrella = mean_na(or_var_metaumbrella),
    bias_ci_or_metaumbrella = mean_na(ci_cov_or_metaumbrella)
    )

for (met in method) {
  print(met)
  res[, paste0("bias_var_or_", met)] <- res[, paste0("mean_or_var_", met)] / res[, paste0("var_or_", met)]
}

err_plot = res %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("bias_es") |
                 starts_with("bias_var") |
                 starts_with("bias_ci"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(or, n, br, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("bias_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("bias_var", method_long, fixed = TRUE) ~ "_var",
                           grepl("bias_ci", method_long, fixed = TRUE) ~ "_ci"),
         method = case_when(grepl("metaumbrella", method_long, fixed = TRUE) ~ "metaumbrella")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot$or <- paste0("OR=", err_plot$or)
err_plot$br <- paste0("br=", err_plot$br)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "300"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=300"))

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/metaConvert/simulations/data_agg/OR_SE_AGG10000.txt")

