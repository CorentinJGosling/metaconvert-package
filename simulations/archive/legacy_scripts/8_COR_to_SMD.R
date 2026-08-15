library(MonteCarlo); library(tidyverse); library(ggplot2)
sd_pooled = function(sd1, sd2, n1, n2) {
  sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))
}
.d_j <- function(x) {
  j <- ifelse(x <= 1, NA, 1) * exp(lgamma(x / 2) - 0.5 * log(x / 2) - lgamma((x - 1) / 2))
  return(j)
}

means_adj_to_smd <- function(d, n, p) {
  n_exp = round(n*p)
  n_nexp = n - n_exp

  grp_exp = rnorm(n_exp, d, 1)
  grp_nexp = rnorm(n_nexp, 0, 1)

  mean_exp <- mean(grp_exp)
  mean_sd_exp <- sd(grp_exp)
  mean_nexp <- mean(grp_nexp)
  mean_sd_nexp <- sd(grp_nexp)
  pooled_sd_crude = sd_pooled(sd1 = mean_sd_exp,  n1 = n_exp,
                              sd2 = mean_sd_nexp, n2 = n_nexp)
  d_raw_data = (mean_exp - mean_nexp) / pooled_sd_crude

  cor_grp = cor.test(~grp_exp + grp_nexp)
  r_raw_data = as.numeric(cor_grp$estimate)

  # viechtbauer ===========================
    res_g <- metafor::conv.delta(
      yi = r, vi = r_se^2, transf = metafor::transf.rtod, var.names = c("g", "g_var")
    )
    J <- .d_j(n_exp+n_nexp - 2)
    res_g$d <- res_g$g / J
    res_g$d_se <- sqrt(res_g$g_var / J^2)

  # Cooper ===========================
  d <- 2 * r / sqrt(1 - r^2)
  d_se <- sqrt(4 * sqrt(r_se) / ((1 - r^2)^3))

  # Mathur  =======================

  increase <- ifelse(unit_type == "sd",
                     unit_increase_iv * sd_iv,
                     unit_increase_iv)

  d <- r * increase / (sd_iv * sqrt((1 - r^2)))
  d_se <- abs(d) * sqrt(1 / (r^2 * (n_sample - 3)) + 1 / (2 * (n_sample - 1)))


  sd_sigma = mod_wo$sigma / sqrt(1 - r_cov_guess^2)
  sd_pooled_est = sd_pooled(sd1 = ancova_mean_sd_exp_model, n1 = n_exp,
                            sd2 = ancova_mean_sd_nexp_model, n2 = n_nexp) / sqrt(1 - r_cov_guess^2)
  sd_min = min(sds) / sqrt(1 - r_cov_guess^2)
  sd_max = max(sds) / sqrt(1 - r_cov_guess^2)

  d_crude_est = (ancova_mean_exp - ancova_mean_nexp) / sd_pooled_est
  d_se_crude_est = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_crude_est^2/(2*(n_sample)))
  d_ci_lo_crude_est <- d_crude_est - qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_crude_est
  d_ci_up_crude_est <- d_crude_est + qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_crude_est

  d_sigma = (ancova_mean_exp - ancova_mean_nexp) / sd_sigma
  d_se_sigma = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_sigma^2/(2*(n_sample)))
  d_ci_lo_sigma <- d_sigma - qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_sigma
  d_ci_up_sigma <- d_sigma + qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_sigma

  d_min = (ancova_mean_exp - ancova_mean_nexp) / sd_min
  d_se_min = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_min^2/(2*(n_sample)))
  d_ci_lo_min <- d_min - qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_min
  d_ci_up_min <- d_min + qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_min

  d_max = (ancova_mean_exp - ancova_mean_nexp) / sd_max
  d_se_max = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_max^2/(2*(n_sample)))
  d_ci_lo_max <- d_max - qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_max
  d_ci_up_max <- d_max + qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_max

  d_ttest <- -coef(summary(mod_int))[2,3] * sqrt(1/n_exp + 1/n_nexp) * sqrt(1 - r_cov_guess^2)
  d_se_ttest = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_ttest^2/(2*(n_sample)))
  d_ci_lo_ttest <- d_ttest - qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_ttest
  d_ci_up_ttest <- d_ttest + qt(.975, df = n_exp + n_nexp - 2 - 1) * d_se_ttest

  return(list(
    "r_raw_data" = r_raw_data,

    "d_raw_data" = d_raw_data,
    "d_se_raw_data" = d_se_raw_data,
    "d_ci_lo_raw_data" = d_ci_lo_raw_data,
    "d_ci_up_raw_data" = d_ci_up_raw_data,

    "d_ttest" = d_ttest,
    "d_se_ttest" = d_se_ttest,
    "d_ci_lo_ttest" = d_ci_lo_ttest,
    "d_ci_up_ttest" = d_ci_up_ttest,

    "d_crude_est" = d_crude_est,
    "d_se_crude_est" = d_se_crude_est,
    "d_ci_lo_crude_est" = d_ci_lo_crude_est,
    "d_ci_up_crude_est" = d_ci_up_crude_est,

    "d_sigma" = d_sigma,
    "d_se_sigma" = d_se_sigma,
    "d_ci_lo_sigma" = d_ci_lo_sigma,
    "d_ci_up_sigma" = d_ci_up_sigma,

    "d_min" = d_min,
    "d_se_min" = d_se_min,
    "d_ci_lo_min" = d_ci_lo_min,
    "d_ci_up_min" = d_ci_up_min,

    "d_max" = d_max,
    "d_se_max" = d_se_max,
    "d_ci_lo_max" = d_ci_lo_max,
    "d_ci_up_max" = d_ci_up_max

    ))
}

r_cov_guess_grid <- c(0.1, 0.3, 0.5)
d_grid <- c(0, 0.2, 0.5, 0.8)
r_grid <- c(0.1, 0.3, 0.5)
n_grid <- c(25, 50, 75, 100, 300)
p_grid <- c(0.3, 0.5, 0.7)
# collect parameter grids in list:
param_list = list("d" = d_grid, "r_cov" = r_grid, "n" = n_grid, "p" = p_grid, "r_cov_guess" = r_cov_guess_grid)
MC_result <- MonteCarlo(func = means_adj_to_smd, nrep = 2500,
                        time_n_test = TRUE,
                        param_list = param_list, max_grid = 10000)
datsim <- MakeFrame(MC_result)

rio::export(datsim, "D:/simulations/data/MEANS_ADJ_to_SMD_sim2500d.txt")

dat = rbind(
  read.delim("D:/simulations/data/MEANS_ADJ_to_SMD_sim2500a.txt"),
  read.delim("D:/simulations/data/MEANS_ADJ_to_SMD_sim2500b.txt"),
  read.delim("D:/simulations/data/MEANS_ADJ_to_SMD_sim2500c.txt"),
  read.delim("D:/simulations/data/MEANS_ADJ_to_SMD_sim2500d.txt"))

mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}

checks = dat %>%
  group_by(d, r_cov) %>%
  summarise(d_raw_data = mean(d_raw_data),
            r_raw_data = mean(r_raw_data))

method = c("crude_est", "ttest", "sigma", "min", "max")
for (met in method) {
  dat[, paste0("bias_es_d_", met)] <- dat[, paste0("d_", met)] - dat$d_raw_data
  dat[, paste0("d_var_", met)] = dat[, paste0("d_se_", met)]^2

  dat[, paste0("ci_cov_d_", met)] = dat[, paste0("d_ci_lo_", met)] <= dat$d_raw_data &
                                    dat[, paste0("d_ci_up_", met)] >= dat$d_raw_data
}

res = dat %>%
  group_by(d, n, r_cov, p, r_cov_guess) %>%
  summarise(
    n_sim = n(),
    bias_es_d_crude_est = mean(bias_es_d_crude_est),
    var_d_crude_est = var(d_crude_est),
    mean_d_var_crude_est = mean(d_var_crude_est),
    bias_ci_d_crude_est = mean(ci_cov_d_crude_est),

    bias_es_d_sigma = mean(bias_es_d_sigma),
    var_d_sigma = var(d_sigma),
    mean_d_var_sigma = mean(d_var_sigma),
    bias_ci_d_sigma = mean(ci_cov_d_sigma),

    bias_es_d_ttest = mean(bias_es_d_ttest),
    var_d_ttest = var(d_ttest),
    mean_d_var_ttest = mean(d_var_ttest),
    bias_ci_d_ttest = mean(ci_cov_d_ttest),

    bias_es_d_min = mean(bias_es_d_min),
    var_d_min = var(d_min),
    mean_d_var_min = mean(d_var_min),
    bias_ci_d_min = mean(ci_cov_d_min),

    bias_es_d_max = mean(bias_es_d_max),
    var_d_max = var(d_max),
    mean_d_var_max = mean(d_var_max),
    bias_ci_d_max = mean(ci_cov_d_max)
  )

for (met in method) {
  print(met)
  res[, paste0("bias_var_d_", met)] <- res[, paste0("mean_d_var_", met)] / res[, paste0("var_d_", met)]
}

err_plot = res %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("bias_es") |
                 starts_with("bias_var") |
                 starts_with("bias_ci"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(d, n, r_cov, p, r_cov_guess, method_long) %>%
  summarise(bias = mean(bias)) %>%
  mutate(pivot = case_when(grepl("bias_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("bias_var", method_long, fixed = TRUE) ~ "_var",
                           grepl("bias_ci", method_long, fixed = TRUE) ~ "_ci"),
         method = case_when(grepl("crude_est", method_long, fixed = TRUE) ~ "crude_est",
                            grepl("sigma", method_long, fixed = TRUE) ~ "sigma",
                            grepl("ttest", method_long, fixed = TRUE) ~ "ttest",
                            grepl("min", method_long, fixed = TRUE) ~ "min",
                            grepl("max", method_long, fixed = TRUE) ~ "max")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot$d <- paste0("d=", err_plot$d)
err_plot$r_cov <- paste0("r_cov=", err_plot$r_cov)
err_plot$r_cov_guess <- paste0("guessed r_cov=", err_plot$r_cov_guess)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "300"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=300"))
err_plot$method[err_plot$method == "crude_est"] <- "SD = pooled SDs"
err_plot$method[err_plot$method == "sigma"] <- "SD = sigma(model)"
err_plot$method[err_plot$method == "min"] <- "SD = min(SDs)"
err_plot$method[err_plot$method == "max"] <- "SD = max(SDs)"
err_plot$method[err_plot$method == "ttest"] <- "Comparison (from t-test)"

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/metaConvert/simulations/data_agg/MEANS_ADJ_to_SMD_AGG10000.txt")
