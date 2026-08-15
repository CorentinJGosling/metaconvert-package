library(MonteCarlo); library(tidyverse); library(ggplot2)

md_to_d <- function(d, n, p) {

  # generate sample scores
  n_exp = round(p * n)
  n_nexp = n - n_exp
  scores_grp1 <- rnorm(n_exp, d, 1)
  scores_grp2 <- rnorm(n_nexp, 0, 1)
  mean_exp <- mean(scores_grp1)
  mean_nexp <- mean(scores_grp2)
  mean_sd_exp <- sd(scores_grp1)
  mean_sd_nexp <- sd(scores_grp2)

  # SMD
  pooled_sd <- sqrt(((n_exp - 1)*mean_sd_exp^2 + (n_nexp - 1)*mean_sd_nexp^2) / (n_exp + n_nexp - 2))
  d_raw_data <- (mean_exp - mean_nexp) / pooled_sd
  d_se_raw_data <- sqrt(1/n_exp  + 1/n_nexp + (d_raw_data^2)/(2*(n_exp+n_nexp)))
  pval_raw_data = 1 - 2 * abs(pt(d_raw_data/d_se_raw_data, n_exp+n_nexp-2) - 0.5)

  # MD - assumption equal variance
  md_equal = mean_exp - mean_nexp
  md_se_equal = sqrt((1/n_exp + 1/n_nexp) * pooled_sd^2)
  md_sd_equal = md_se_equal / sqrt(1/n_exp + 1/n_nexp)

  d_equal = md_equal / md_sd_equal
  d_se_equal <- sqrt(1/n_exp  + 1/n_nexp + (d_equal^2)/(2*(n_exp+n_nexp)))
  pval_equal = 1 - 2 * abs(pt(d_equal/d_se_equal, n_exp+n_nexp-2) - 0.5)

  # MD - assumption non-equal variance
  md_non_equal = mean_exp - mean_nexp
  md_se_non_equal = sqrt(mean_sd_exp^2/n_exp + mean_sd_nexp^2/n_nexp)
  md_sd_non_equal = md_se_non_equal / sqrt(1/n_exp + 1/n_nexp)

  d_non_equal = md_non_equal / md_sd_non_equal
  d_se_non_equal <- sqrt(1/n_exp  + 1/n_nexp + (d_non_equal^2)/(2*(n_exp+n_nexp)))
  pval_non_equal = 1 - 2 * abs(pt(d_non_equal/d_se_non_equal, n_exp+n_nexp-2) - 0.5)

  # SMD convert

  return(list(
    "d_raw_data" = d_raw_data,
    "d_se_raw_data" = d_se_raw_data,
    "pval_raw_data" = pval_raw_data,
    "d_equal" = d_equal,
    "d_se_equal" = d_se_equal,
    "pval_equal" = pval_equal,
    "d_non_equal" = d_non_equal,
    "d_se_non_equal" = d_se_non_equal,
    "pval_non_equal" = pval_non_equal))
}

d_grid <- c(0, 0.2, 0.5, 0.8)
n_grid <- c(25, 50, 75, 100, 200)
p_grid <- c(0.3, 0.5, 0.7)

# Collect MCMC results ---------------
param_list = list("d" = d_grid, "n" = n_grid, "p" = p_grid)
MC_result <- MonteCarlo(func = md_to_d, nrep = 5000, param_list = param_list)
res <- MakeFrame(MC_result)

# Export the raw dataset ---------------
setwd("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data")
rio::export(res, "MD_to_SMD_sim.txt")

res = read.delim("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/MD_to_SMD_sim.txt")
# Check accuracy dataset creation ---------------
checks = res %>%
  group_by(d) %>%
  summarise(d_raw_data = mean(d_raw_data))

res_long = res %>%
  mutate(acc_es_d_equal = abs(d_equal - d_raw_data),
         acc_se_d_equal = abs(d_se_equal - d_se_raw_data),
         acc_pval_equal = abs(pval_equal - pval_raw_data),
         acc_es_d_non_equal = abs(d_non_equal - d_raw_data),
         acc_se_d_non_equal = abs(d_se_non_equal - d_se_raw_data),
         acc_pval_d_non_equal = abs(pval_non_equal - pval_raw_data)) %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("acc_es") |
                 starts_with("acc_se") |
                 starts_with("acc_pval"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(d, n, p, method_long) %>%
  summarise(bias = mean(bias),
            n_sim = n()) %>%
  mutate(pivot = case_when(grepl("acc_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("acc_se", method_long, fixed = TRUE) ~ "_se",
                           grepl("acc_pval", method_long, fixed = TRUE) ~ "_pval"),
        method = case_when(grepl("non_equal", method_long, fixed = TRUE) ~ "metaumbrella_non_equal_var",
                           !grepl("non_equal", method_long, fixed = TRUE) ~ "metaumbrella_equal_var")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")
res_long$n = paste0("n=", res_long$n)
res_long$d = paste0("d=", res_long$d)

rio::export(res_long, "MD_to_SMD_AGG.txt")
