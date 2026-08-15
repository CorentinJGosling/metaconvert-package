library(MonteCarlo); library(tidyverse)
cont_to_z <- function(r, n, p1, p2) {

  # generate sample scores
  sig <- rbind(c(1, r), c(r, 1))

  # create the mean vector
  mu <- c(0, 0)
  # generate the multivariate normal distribution
  res <- as.data.frame(MASS::mvrnorm(n=n, mu=mu, Sigma=sig))

  # ==============================================================
  # Estimate COR from raw data =================================
  # ==============================================================
  r_raw_data = as.numeric(cor.test(~res$V1 + res$V2)$estimate)
  r_se_raw_data =  sqrt((1-r_raw_data^2)^2 / (n-1))
  r_pval_raw_data = 1 - 2 * abs(pt(r_raw_data / r_se_raw_data, n - 2) - 0.5)

  z_raw_data = atanh(r_raw_data)
  z_se_raw_data = sqrt(1/(n - 3))
  z_pval_raw_data = 1 - 2 * abs(pnorm(z_raw_data / z_se_raw_data) - 0.5)

  # ==============================================================
  # Obtain 2x2 table  ============================================
  # ==============================================================
  x1 = 1-p1
  x2 = 1-p2
  cut1 = quantile(res$V1, x1)
  cut2 = quantile(res$V2, x2)
  n_cases_exp = nrow(subset(res, V1 >= cut1 & V2 >= cut2))
  n_controls_exp = nrow(subset(res, V1 < cut1 & V2 >= cut2))
  n_cases_nexp = nrow(subset(res, V1 >= cut1 & V2 < cut2))
  n_controls_nexp = nrow(subset(res, V1 < cut1 & V2 < cut2))
  n_exp = n_cases_exp + n_controls_exp
  n_nexp =  n_cases_nexp + n_controls_nexp
  n_cases = n_cases_exp + n_cases_nexp
  n_controls = n_controls_exp + n_controls_nexp
  n_sample = n_exp + n_nexp

  # ==============================================================
  # Estimate OR from raw data =================================
  # ==============================================================
  zero = which(n_cases_exp == 0 | n_cases_nexp == 0 | n_controls_exp == 0 | n_controls_nexp == 0)

  n_cases_exp[zero] = n_cases_exp[zero] + 0.5
  n_cases_nexp[zero] = n_cases_nexp[zero] + 0.5
  n_controls_exp[zero] = n_controls_exp[zero] + 0.5
  n_controls_nexp[zero] = n_controls_nexp[zero] + 0.5

  or_raw_data = suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp))
  or_se_raw_data = suppressWarnings(sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp))

  # ==============================================================
  # Conversion 1. Cooper 2x2 => OR => SMD => R ===================
  # ==============================================================
  d = log(or_raw_data) * sqrt(3) / pi
  vd =  sqrt(or_se_raw_data^2 * 3 / (pi^2))

  a = ((n_exp + n_nexp)^2) / (n_exp*n_nexp)
  p = n_exp / (n_exp + n_nexp)

  r_cooper = d/sqrt(d^2 + 1/(p * (1 - p)))
  r_se_cooper = sqrt(a^2 * vd / ((d^2 + a)^3))

  z_cooper = atanh(r_cooper)
  z_se_cooper = sqrt(vd / (vd + 1/(p*(1-p))))

  r_pval_cooper = 1 - 2 * abs(pt(r_cooper / r_se_cooper, n_sample - 2) - 0.5)
  z_pval_cooper = 1 - 2 * abs(pnorm(z_cooper / z_se_cooper) - 0.5)

  # ==============================================================
  # Conversion 2. bonett =========================================
  # ==============================================================
  small_margin_prop = min(
    n_cases/n_sample,
    n_controls/n_sample,
    n_exp/n_sample,
    n_nexp/n_sample
  )

  c = (1 - abs(n_exp - n_cases)/5 - (1/2 - small_margin_prop)^2)/2

  r_bonett = cos(pi/(1+or_se_raw_data^c))
  r_se_bonett = or_se_raw_data * (pi*c*or_se_raw_data^c) * sin(pi/(1+or_se_raw_data^c)) / (1+or_se_raw_data^c)^2

  # or_ci_lo = exp(log(or_se_raw_data) - qnorm(.975)*or_se_raw_data)
  # or_ci_up = exp(log(or_se_raw_data) + qnorm(.975)*or_se_raw_data)
  # r_ci_lo_bonett = cos(pi/(1 + or_ci_lo^c))
  # r_ci_up_bonett = cos(pi/(1 + or_ci_up^c))

  z_bonett = atanh(r_bonett)
  z_se_bonett = sqrt(r_se_bonett^2 / ((1 - r_bonett^2)^2)) # delta method
  # z_ci_lo_bonett = atanh(r_ci_lo_bonett)
  # z_ci_lo_bonett = atanh(r_ci_up_bonett)

  r_pval_bonett = 1 - 2 * abs(pt(r_bonett / r_se_bonett, n_sample - 2) - 0.5)
  z_pval_bonett = 1 - 2 * abs(pnorm(z_bonett / z_se_bonett) - 0.5)

  # ==============================================================
  # Conversion 3. digby =========================================
  # ==============================================================
  c = 1/2

  r_pearson = cos(pi/(1+or_se_raw_data^c))
  r_se_pearson = or_se_raw_data * (pi*c*or_se_raw_data^c) * sin(pi/(1+or_se_raw_data^c)) / (1+or_se_raw_data^c)^2

  # r_ci_lo_pearson = cos(pi/(1 + or_ci_lo^c))
  # r_ci_up_pearson = cos(pi/(1 + or_ci_up^c))

  z_pearson = atanh(r_pearson)
  z_se_pearson = sqrt(r_se_pearson^2 / ((1 - r_pearson^2)^2)) # delta method
  # z_ci_lo_pearson = atanh(r_ci_lo_pearson)
  # z_ci_lo_pearson = atanh(r_ci_up_pearson)
  r_pval_pearson = 1 - 2 * abs(pt(r_pearson / r_se_pearson, n_sample - 2) - 0.5)
  z_pval_pearson = 1 - 2 * abs(pnorm(z_pearson / z_se_pearson) - 0.5)


  # ==============================================================
  # Conversion 3. digby  =======================================
  # ==============================================================
  c = 3/4

  r_digby = (or_se_raw_data^c - 1)/(or_se_raw_data^c + 1)
  r_se_digby = sqrt((c^2 / 4) * (1 - r_digby^2)^2 * or_se_raw_data^2)

  z_digby = atanh(r_digby)
  z_se_digby = sqrt(r_se_digby^2 / ((1 - r_digby^2)^2)) # delta method
  # z_ci_lo_digby = z_digby - qnorm(.975) * sqrt(c^2/4 * or_se_raw_data^2)
  # z_ci_up_digby = z_digby + qnorm(.975) * sqrt(c^2/4 * or_se_raw_data^2)

  # r_ci_lo_digby = tanh(z_ci_lo_digby)
  # r_ci_up_digby = tanh(z_ci_up_digby)
  r_pval_digby = 1 - 2 * abs(pt(r_digby / r_se_digby, n_sample - 2) - 0.5)
  z_pval_digby = 1 - 2 * abs(pnorm(z_digby / z_se_digby) - 0.5)

  # return result:
  return(list(
    "ratio_cases" = n_cases/n_sample,
    "ratio_exp" = n_exp/n_sample,
    "or_raw_data" = or_raw_data,
    "or_se_raw_data" = or_se_raw_data,
    "br_raw_data" = n_cases_nexp / n_nexp,
    "r_raw_data" = r_raw_data,
    "r_se_raw_data" = r_se_raw_data,
    "r_pval_raw_data" = r_pval_raw_data,
    "z_raw_data" = z_raw_data,
    "z_se_raw_data" = z_se_raw_data,
    "z_pval_raw_data" = z_pval_raw_data,
    # es
    "r_cooper"           = r_cooper,
    "z_cooper"           = z_cooper,
    "r_bonett"           = r_bonett,
    "z_bonett"           = z_bonett,
    "r_pearson"          = r_pearson,
    "z_pearson"          = z_pearson,
    "r_digby"            = r_digby,
    "z_digby"            = z_digby,
    # se
    "r_se_cooper"         = r_se_cooper,
    "z_se_cooper"         = z_se_cooper,
    "r_se_bonett"         = r_se_bonett,
    "z_se_bonett"         = z_se_bonett,
    "r_se_pearson"        = r_se_pearson,
    "z_se_pearson"        = z_se_pearson,
    "r_se_digby"          = r_se_digby,
    "z_se_digby"          = z_se_digby,
    # pval
    "r_pval_cooper"       = r_pval_cooper,
    "z_pval_cooper"       = z_pval_cooper,
    "r_pval_bonett"       = r_pval_bonett,
    "z_pval_bonett"       = z_pval_bonett,
    "r_pval_pearson"      = r_pval_pearson,
    "z_pval_pearson"      = z_pval_pearson,
    "r_pval_digby"        = r_pval_digby,
    "z_pval_digby"        = z_pval_digby
  ))
}
r_grid <- c(0, 0.1, 0.2, 0.3, 0.5)
n_grid <- c(25, 50, 75, 100, 200)
p1_grid = c(0.1, 0.3, 0.5)
p2_grid = c(0.1, 0.3, 0.5)
# p = expand.grid(p_gen, p_gen)
# p1_grid <- p[,1]
# p2_grid <- p[,2]

# collect parameter grids in list:
param_list = list("r" = r_grid, "n" = n_grid, "p1" = p1_grid, "p2" = p1_grid)
MC_result <- MonteCarlo(func = cont_to_z, nrep = 5000, param_list = param_list, max_grid = 10000)
df <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data")
rio::export(df, "OR_to_COR_sim.txt")


res_cont = read.delim("OR_to_COR_sim.txt")
res_cont[res_cont==Inf|res_cont==-Inf] <- NA
mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}

checks = res_cont %>%
  group_by(r) %>%
  summarise(r_raw_data = mean_na(r_raw_data))

method = c("pearson", "digby", "bonett", "cooper")
for (met in method) {
  res_cont[, paste0("acc_es_r_", met)] <- abs(res_cont[, paste0("r_", met)] - res_cont$r_raw_data)
  res_cont[, paste0("acc_se_r_", met)] <- abs(res_cont[, paste0("r_se_", met)] - res_cont$r_se_raw_data)
  res_cont[, paste0("acc_pval_r_", met)] <- abs(res_cont[, paste0("r_pval_", met)] - res_cont$r_pval_raw)
  res_cont[, paste0("acc_es_z_", met)] <- abs(res_cont[, paste0("z_", met)] - res_cont$z_raw_data)
  res_cont[, paste0("acc_se_z_", met)] <- abs(res_cont[, paste0("z_se_", met)] - res_cont$z_se_raw_data)
  res_cont[, paste0("acc_pval_z_", met)] <- abs(res_cont[, paste0("z_pval_", met)] - res_cont$z_pval_raw)
}
res = res_cont %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("acc_es") |
                 starts_with("acc_se") |
                 starts_with("acc_pval"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(r, n, method_long) %>% #p1, p2,
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("acc_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("acc_se", method_long, fixed = TRUE) ~ "_se",
                           grepl("acc_pval", method_long, fixed = TRUE) ~ "_pval"),
          method = case_when(grepl("r_digby", method_long, fixed = TRUE) ~ "r_digby",
                             grepl("z_digby", method_long, fixed = TRUE) ~ "z_digby",
                             grepl("r_bonett", method_long, fixed = TRUE) ~ "r_bonett",
                             grepl("z_bonett", method_long, fixed = TRUE) ~ "z_bonett",
                             grepl("r_pearson", method_long, fixed = TRUE) ~ "r_pearson",
                             grepl("z_pearson", method_long, fixed = TRUE) ~ "z_pearson",
                             grepl("r_cooper", method_long, fixed = TRUE) ~ "r_cooper",
                             grepl("z_cooper", method_long, fixed = TRUE) ~ "z_cooper")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot = res

err_plot$method[err_plot$method == "r_cooper"] <- "Cooper (r)"
err_plot$method[err_plot$method == "z_cooper"] <- "Cooper (z)"
err_plot$method[err_plot$method == "r_bonett"] <- "Bonett (r)"
err_plot$method[err_plot$method == "z_bonett"] <- "Bonett (z)"
err_plot$method[err_plot$method == "r_pearson"] <- "Pearson (r)"
err_plot$method[err_plot$method == "z_pearson"] <- "Pearson (z)"
err_plot$method[err_plot$method == "r_digby"] <- "Digby (r)"
err_plot$method[err_plot$method == "z_digby"] <- "Digby (z)"

err_plot$r <- paste0("r=", err_plot$r)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "200"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/OR_to_COR_AGG.txt")
