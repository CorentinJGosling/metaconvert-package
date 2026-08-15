library(MonteCarlo); library(tidyverse)

r=0.2; n=200; p1=0.3; p2=0.5
or_to_z <- function(r, n, p1, p2) {

  sig <- rbind(c(1, r), c(r, 1))

  mu <- c(0, 0)
  res <- as.data.frame(MASS::mvrnorm(n=n, mu=mu, Sigma=sig))

  # Estimate COR from raw data =================================
  r_raw_data = as.numeric(cor.test(~res$V1 + res$V2)$estimate)
  z_raw_data = atanh(r_raw_data)

  x1 = 1-p1
  x2 = 1-p2
  cut1 = quantile(res$V1, x1)
  cut2 = quantile(res$V2, x2)

  n_cases_exp = nrow(subset(res, V1 >= cut1 & V2 >= cut2))
  n_controls_exp = nrow(subset(res, V1 >= cut1 & V2 < cut2))
  n_cases_nexp = nrow(subset(res, V1 < cut1 & V2 >= cut2))
  n_controls_nexp = nrow(subset(res, V1 < cut1 & V2 < cut2))
  cbind(r_raw_data, psych::phi(matrix(c(n_cases_exp, n_cases_nexp,
                                        n_controls_exp, n_controls_nexp), nrow=2)))

  n_exp = n_cases_exp + n_controls_exp
  n_nexp =  n_cases_nexp + n_controls_nexp
  n_cases = n_cases_exp + n_cases_nexp
  n_controls = n_controls_exp + n_controls_nexp
  n_sample = n_exp + n_nexp

  zero = which(n_cases_exp == 0 | n_cases_nexp == 0 | n_controls_exp == 0 | n_controls_nexp == 0)
  n_cases_exp[zero] = n_cases_exp[zero] + 0.5
  n_cases_nexp[zero] = n_cases_nexp[zero] + 0.5
  n_controls_exp[zero] = n_controls_exp[zero] + 0.5
  n_controls_nexp[zero] = n_controls_nexp[zero] + 0.5

  or_raw_data = suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp))
  or_se_raw_data = suppressWarnings(sqrt(1/n_cases_exp + 1/n_controls_exp +
                                           1/n_cases_nexp + 1/n_controls_nexp))
  or_ci_lo = exp(log(or_raw_data) - qnorm(.975) * or_se_raw_data)
  or_ci_up = exp(log(or_raw_data) + qnorm(.975) * or_se_raw_data)

  # Conversion 1. Cooper 2x2 => OR => SMD => R ===================
  d = log(or_raw_data) * sqrt(3) / pi
  vd =  sqrt(or_se_raw_data^2 * 3 / (pi^2))

  a = ((n_exp + n_nexp)^2) / (n_exp*n_nexp)
  p = n_exp / (n_exp + n_nexp)

  r_cooper = d/sqrt(d^2 + 1/(p * (1 - p)))
  r_se_cooper = sqrt(a^2 * vd / ((d^2 + a)^3))

  z_cooper = atanh(r_cooper)
  z_se_cooper = sqrt(vd / (vd + 1/(p*(1-p))))

  r_ci_lo_cooper <- r_cooper - qt(.975, df = n_exp + n_nexp - 2) * r_se_cooper
  r_ci_up_cooper <- r_cooper + qt(.975, df = n_exp + n_nexp - 2) * r_se_cooper
  z_ci_lo_cooper <- z_cooper - qnorm(.975) * z_se_cooper
  z_ci_up_cooper <- z_cooper + qnorm(.975) * z_se_cooper

  # Conversion 2. bonett =========================================
  small_margin_prop = min(
    n_cases/n_sample,
    n_controls/n_sample,
    n_exp/n_sample,
    n_nexp/n_sample
  )

  c = (1 - abs(n_exp/n_sample - n_cases/n_sample)/5 -
         (1/2 - small_margin_prop)^2)/2

  r_bonett = cos(pi/(1+or_raw_data^c))
  r_se_bonett = or_se_raw_data * (pi*c*or_raw_data^c) * sin(pi/(1+or_raw_data^c)) / (1+or_raw_data^c)^2

  r_ci_lo_bonett <- cos(pi / (1 + or_ci_lo^c))
  r_ci_up_bonett <- cos(pi / (1 + or_ci_up^c))

  z_bonett = atanh(r_bonett)
  z_se_bonett <- sqrt(r_se_bonett^2 / ((1 - r_bonett^2)^2))
  z_ci_lo_bonett <- atanh(r_ci_lo_bonett)
  z_ci_up_bonett <- atanh(r_ci_up_bonett)

  # Conversion 3. pearson =========================================
  c = 1/2

  r_pearson = cos(pi/(1+or_raw_data^c))
  r_se_pearson = or_se_raw_data * (pi*c*or_raw_data^c) * sin(pi/(1+or_raw_data^c)) / (1+or_raw_data^c)^2

  z_pearson = atanh(r_pearson)
  z_se_pearson = sqrt(r_se_pearson^2 / ((1 - r_pearson^2)^2))
  r_ci_lo_pearson <- cos(pi / (1 + or_ci_lo^c))
  r_ci_up_pearson <- cos(pi / (1 + or_ci_up^c))
  z_ci_lo_pearson = atanh(r_ci_lo_pearson)
  z_ci_up_pearson = atanh(r_ci_up_pearson)

  # Conversion 3. digby  =======================================
  c = 3/4

  r_digby = (or_raw_data^c - 1)/(or_raw_data^c + 1)
  r_se_digby = sqrt((c^2 / 4) * (1 - r_digby^2)^2 * or_se_raw_data^2)

  z_digby = atanh(r_digby)
  z_se_digby = sqrt(r_se_digby^2 / ((1 - r_digby^2)^2))
  z_ci_lo_digby = z_digby - qnorm(.975) * sqrt(c^2/4 * or_se_raw_data^2)
  z_ci_up_digby = z_digby + qnorm(.975) * sqrt(c^2/4 * or_se_raw_data^2)

  r_ci_lo_digby = tanh(z_ci_lo_digby)
  r_ci_up_digby = tanh(z_ci_up_digby)

  # return result:
  return(list(
    "or_raw_data" = or_raw_data,
    "br_raw_data" = n_cases_nexp / n_nexp,
    "r_raw_data" = r_raw_data,
    "z_raw_data" = z_raw_data,
    # es
    "r_pearson"           = r_pearson,
    "z_pearson"           = z_pearson,
    "r_cooper"           = r_cooper,
    "z_cooper"           = z_cooper,
    "r_digby"      = r_digby,
    "z_digby"      = z_digby,
    "r_bonett"           = r_bonett,
    "z_bonett"           = z_bonett,
    # se
    "r_se_pearson"         = r_se_pearson,
    "z_se_pearson"         = z_se_pearson,
    "r_se_cooper"         = r_se_cooper,
    "z_se_cooper"         = z_se_cooper,
    "r_se_digby"    = r_se_digby,
    "z_se_digby"    = z_se_digby,
    "r_se_bonett"         = r_se_bonett,
    "z_se_bonett"         = z_se_bonett,
    # se
    "r_ci_lo_cooper"      = r_ci_lo_cooper,
    "z_ci_lo_cooper"      = z_ci_lo_cooper,
    "r_ci_lo_pearson"      = r_ci_lo_pearson,
    "z_ci_lo_pearson"      = z_ci_lo_pearson,
    "r_ci_lo_digby"    = r_ci_lo_digby,
    "z_ci_lo_digby"    = z_ci_lo_digby,
    "r_ci_lo_bonett"         = r_ci_lo_bonett,
    "z_ci_lo_bonett"         = z_ci_lo_bonett,
    # se
    "r_ci_up_pearson"      = r_ci_up_pearson,
    "z_ci_up_pearson"      = z_ci_up_pearson,
    "r_ci_up_cooper"      = r_ci_up_cooper,
    "z_ci_up_cooper"      = z_ci_up_cooper,
    "r_ci_up_digby"    = r_ci_up_digby,
    "z_ci_up_digby"    = z_ci_up_digby,
    "r_ci_up_bonett"         = r_ci_up_bonett,
    "z_ci_up_bonett"         = z_ci_up_bonett
  ))
}
r_grid <- c(0, 0.1, 0.2, 0.3, 0.5)
n_grid <- c(25, 50, 75, 100, 300)
p1_grid = c(0.3, 0.5, 0.7)
p2_grid = c(0.3, 0.5, 0.7)

# collect parameter grids in list:
param_list = list("r" = r_grid, "n" = n_grid, "p1" = p1_grid, "p2" = p1_grid)

start_time <- Sys.time()
MC_result <- MonteCarlo(func = or_to_z, nrep = 10000, param_list = param_list, max_grid = 10000)
end_time <- Sys.time()
end_time - start_time

data_sim <- MakeFrame(MC_result)
rio::export(data_sim, "D:/simulations/data/OR_to_COR_sim.txt")

res_cont = read.delim("D:/simulations/data/OR_to_COR_sim.txt")
res_cont[res_cont==Inf|res_cont==-Inf] <- NA
checks = res_cont %>%
  group_by(r) %>%
  summarise(r_raw_data = mean(r_raw_data))
dat = res_cont
mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}

method = c("pearson", "digby", "bonett", "cooper")
for (met in method) {
  dat[, paste0("bias_es_r_", met)] <- dat[, paste0("r_", met)] - dat$r_raw_data
  dat[, paste0("r_var_", met)] = dat[, paste0("r_se_", met)]^2

  dat[, paste0("ci_cov_r_", met)] = dat[, paste0("r_ci_lo_", met)] <= dat$r_raw_data &
    dat[, paste0("r_ci_up_", met)] >= dat$r_raw_data

  dat[, paste0("bias_es_z_", met)] <- dat[, paste0("z_", met)] - dat$z_raw_data
  dat[, paste0("z_var_", met)] = dat[, paste0("z_se_", met)]^2

  dat[, paste0("ci_cov_z_", met)] = dat[, paste0("z_ci_lo_", met)] <= dat$z_raw_data &
    dat[, paste0("z_ci_up_", met)] >= dat$r_raw_data
}

res = dat %>%
  group_by(r, n) %>%
  summarise(
    n_sim = n(),
    bias_es_r_pearson = mean_na(bias_es_r_pearson),
    var_r_pearson = var(r_pearson),
    mean_r_var_pearson = mean_na(r_var_pearson),
    bias_ci_r_pearson = mean_na(ci_cov_r_pearson),

    bias_es_r_digby = mean_na(bias_es_r_digby),
    var_r_digby = var(r_digby),
    mean_r_var_digby = mean_na(r_var_digby),
    bias_ci_r_digby = mean_na(ci_cov_r_digby),

    bias_es_r_bonett = mean_na(bias_es_r_bonett),
    var_r_bonett = var(r_bonett),
    mean_r_var_bonett = mean_na(r_var_bonett),
    bias_ci_r_bonett = mean_na(ci_cov_r_bonett),

    bias_es_r_cooper = mean_na(bias_es_r_cooper),
    var_r_cooper = var(r_cooper),
    mean_r_var_cooper = mean_na(r_var_cooper),
    bias_ci_r_cooper = mean_na(ci_cov_r_cooper),

    bias_es_z_pearson = mean_na(bias_es_z_pearson),
    var_z_pearson = var(r_pearson),
    mean_z_var_pearson = mean_na(r_var_pearson),
    bias_ci_z_pearson = mean_na(ci_cov_z_pearson),

    bias_es_z_digby = mean_na(bias_es_z_digby),
    var_z_digby = var(r_digby),
    mean_z_var_digby = mean_na(r_var_digby),
    bias_ci_z_digby = mean_na(ci_cov_z_digby),

    bias_es_z_bonett = mean_na(bias_es_z_bonett),
    var_z_bonett = var(r_bonett),
    mean_z_var_bonett = mean_na(r_var_bonett),
    bias_ci_z_bonett = mean_na(ci_cov_z_bonett),

    bias_es_z_cooper = mean_na(bias_es_z_cooper),
    var_z_cooper = var(r_cooper),
    mean_z_var_cooper = mean_na(r_var_cooper),
    bias_ci_z_cooper = mean_na(ci_cov_z_cooper)

  )

for (met in method) {
  print(met)
  res[, paste0("bias_var_r_", met)] <- res[, paste0("mean_r_var_", met)] / res[, paste0("var_r_", met)]
  res[, paste0("bias_var_z_", met)] <- res[, paste0("mean_z_var_", met)] / res[, paste0("var_z_", met)]
}

err_plot = res %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("bias_es") |
                 starts_with("bias_var") |
                 starts_with("bias_ci"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(r, n, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("bias_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("bias_var", method_long, fixed = TRUE) ~ "_var",
                           grepl("bias_ci", method_long, fixed = TRUE) ~ "_ci"),
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

err_plot$method[err_plot$method == "r_cooper"] <- "Cooper (r)"
err_plot$method[err_plot$method == "z_cooper"] <- "Cooper (z)"
err_plot$method[err_plot$method == "r_bonett"] <- "Bonett (r)"
err_plot$method[err_plot$method == "z_bonett"] <- "Bonett (z)"
err_plot$method[err_plot$method == "r_pearson"] <- "Pearson (r)"
err_plot$method[err_plot$method == "z_pearson"] <- "Pearson (z)"
err_plot$method[err_plot$method == "r_digby"] <- "Digby (r)"
err_plot$method[err_plot$method == "z_digby"] <- "Digby (z)"

err_plot$r <- paste0("r=", err_plot$r)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "300"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=300"))
rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/metaConvert/simulations/data_agg/OR_to_COR_AGG.txt")
nrow(read.delim("D:/simulations/data/SMD_to_COR_sim.txt"))
