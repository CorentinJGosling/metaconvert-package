library(MonteCarlo); library(tidyverse); library(ggplot2)
library(emmeans)
sd_pooled = function(sd1, sd2, n1, n2) {
  sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))
}



# time4 <- Sys.time()
# time5 <- Sys.time()
# time6 <- Sys.time()
# end_time - start_time

#
d = 0.5
r_cov = 0.1
n=1000000
p=0.3
r_cov_guess=0.1

means_adj_to_smd <- function(d, r_cov, r_cov_guess, n, p) {
  # start_time <- Sys.time()
  n_exp = round(n*p)
  n_nexp = n - n_exp

  # generate sample scores
  sig <- rbind(c(1, r_cov), c(r_cov, 1))
  mu_exp <- c(d, -0.1)
  mu_nexp <- c(0, 0.2)

  # generate the VD + set up correlation between cov & VD
  # (this information is needed in the formulas)
  res_exp <- as.data.frame(MASS::mvrnorm(n=n_exp, mu=mu_exp, Sigma=sig))
  res_nexp <- as.data.frame(MASS::mvrnorm(n=n_nexp, mu=mu_nexp, Sigma=sig))

  # dataset ready
  # the dataset is composed of 3 variables : V1=DV, V2=covariate, group
  res = rbind(cbind(res_exp, group="g1"),
              cbind(res_nexp, group="g2"))

  # checking the r is correct
  mod_cov = summary(lm(res$V1 ~ res$V2))
  r_raw_data = sqrt(mod_cov$r.squared)

  # variables for calculations
  mean_exp <- mean(res_exp$V1)
  mean_sd_exp <- sd(res_exp$V1)
  mean_nexp <- mean(res_nexp$V1)
  mean_sd_nexp <- sd(res_nexp$V1)
  mod_int = lm(V1 ~ group + V2, res)
  mod_wo = summary(lm(V1 ~ 0 + group + V2, res))

  # this is taking 2/3 time of the simulations.
  # Do you know another way to estimate adjusted means?
  res_ancova = data.frame(emmeans(mod_int, ~group))

  # variables
  n_sample = n_exp+n_nexp
  ancova_mean_exp = res_ancova$emmean[1]
  ancova_mean_nexp = res_ancova$emmean[2]
  ancova_mean_sd_exp_model = coef(mod_wo)[1,2] * sqrt(n_exp)
  ancova_mean_sd_nexp_model = coef(mod_wo)[2,2] * sqrt(n_nexp)
  sds = c(ancova_mean_sd_exp_model, ancova_mean_sd_nexp_model)

  # ============================================================
  # Determine true SMD from raw data ===========================
  # ============================================================

  # true raw SMD
  pooled_sd_crude = sd_pooled(sd1 = mean_sd_exp,  n1 = n_exp,
                              sd2 = mean_sd_nexp, n2 = n_nexp)
  d_crude = (mean_exp - mean_nexp) / pooled_sd_crude
  d_se_crude = sqrt( (1/n_exp + 1/n_nexp) + d^2/(2*(n_sample)) )
  p_value_crude = 1 - 2 * abs(pt(d_crude / d_se_crude, n_sample-2) - 0.5)

  # true adjusted SMD
  d_ancova = (ancova_mean_exp - ancova_mean_nexp) / pooled_sd_crude
  d_se_ancova = sqrt( (1/n_exp + 1/n_nexp) * (1-r_raw_data^2) + d_ancova^2/(2*(n_sample)) )
  p_value_ancova = 1 - 2 * abs(pt(d_ancova / d_se_ancova,
                                  n_exp + n_nexp - 2 - 1) - 0.5)

  # time3 <- Sys.time()
  # ============================================================
  # Estimate SMD from various input data =======================
  # ============================================================

  # 0. estimated SMD via standard pooled
  sd_std_pooled = sd_pooled(sd1 = ancova_mean_sd_exp_model, n1 = n_exp,
                            sd2 = ancova_mean_sd_nexp_model, n2 = n_nexp) / sqrt(1 - r_cov_guess^2)

  d_pooled = (ancova_mean_exp - ancova_mean_nexp) / sd_std_pooled
  d_se_pooled = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_pooled^2/(2*(n_sample)))
  p_value_pooled = 1 - 2 * abs(pt(d_pooled / d_se_pooled,
                                 n_exp + n_nexp - 2 - 1) - 0.5)

  # 1. estimated SMD via SIGMA
  sd_sigma = mod_wo$sigma / sqrt(1 - r_cov_guess^2)
  d_sigma = (ancova_mean_exp - ancova_mean_nexp) / sd_sigma
  d_se_sigma = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_sigma^2/(2*(n_sample)))
  p_value_sigma = 1 - 2 * abs(pt(d_sigma / d_se_sigma,
                                 n_exp + n_nexp - 2 - 1) - 0.5)

  # 2. estimated SMD via standard mean
  sd_mean = mean(sds) / sqrt(1 - r_cov_guess^2)

  d_mean = (ancova_mean_exp - ancova_mean_nexp) / sd_mean
  d_se_mean = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_mean^2/(2*(n_sample)))
  p_value_mean = 1 - 2 * abs(pt(d_mean / d_se_mean,
                                 n_exp + n_nexp - 2 - 1) - 0.5)

  # 3. estimated SMD via weighted mean
  # not sure why there is a -1 here, but I strictly followed your formula
  sd_w_mean = sqrt(weighted.mean(sds^2,
                                table(res$group) - 1)) / sqrt(1 - r_cov_guess^2)

  d_w_mean = (ancova_mean_exp - ancova_mean_nexp) / sd_w_mean
  d_se_w_mean = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_w_mean^2/(2*(n_sample)))
  p_value_w_mean = 1 - 2 * abs(pt(d_w_mean / d_se_w_mean,
                                 n_exp + n_nexp - 2 - 1) - 0.5)
  # 4. estimated SMD via SIGMA
  sd_min = min(sds) / sqrt(1 - r_cov_guess^2)

  d_min = (ancova_mean_exp - ancova_mean_nexp) / sd_min
  d_se_min = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_min^2/(2*(n_sample)))
  p_value_min = 1 - 2 * abs(pt(d_min / d_se_min,
                                 n_exp + n_nexp - 2 - 1) - 0.5)
  # 5. estimated SMD via max
  sd_max = max(sds) / sqrt(1 - r_cov_guess^2)

  d_max = (ancova_mean_exp - ancova_mean_nexp) / sd_max
  d_se_max = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_max^2/(2*(n_sample)))
  p_value_max = 1 - 2 * abs(pt(d_max / d_se_max,
                                 n_exp + n_nexp - 2 - 1) - 0.5)
  # 6. estimated SMD via EMMEANS
  ancova_mean_sd_exp_emm = res_ancova$SE[1] * sqrt(n_exp)
  ancova_mean_sd_nexp_emm = res_ancova$SE[2] * sqrt(n_nexp)
  sd_emmeans = sd_pooled(sd1 = ancova_mean_sd_exp_emm,  n1 = n_exp,
                         sd2 = ancova_mean_sd_nexp_emm, n2 = n_nexp) / sqrt(1 - r_cov_guess^2)

  d_emmeans = (ancova_mean_exp - ancova_mean_nexp) / sd_emmeans
  d_se_emmeans = sqrt(((1/n_exp + 1/n_nexp) * (1-r_cov_guess^2)) + d_emmeans^2/(2*(n_sample)))
  p_value_emmeans = 1 - 2 * abs(pt(d_emmeans / d_se_emmeans,
                                 n_exp + n_nexp - 2 - 1) - 0.5)
  # ============================================================
  # output =====================================================
  # ============================================================

  return(list(
    "d_raw_data" = d_crude,
    "r_raw_data" = r_raw_data,

    "d_ancova" = d_ancova,
    "d_se_ancova" = d_se_ancova,
    "p_value_ancova" = p_value_ancova,

    "d_pooled" = d_pooled,
    "d_se_pooled" = d_se_pooled,
    "p_value_pooled" = p_value_pooled,

    "d_sigma" = d_sigma,
    "d_se_sigma" = d_se_sigma,
    "p_value_sigma" = p_value_sigma,

    "d_mean" = d_mean,
    "d_se_mean" = d_se_mean,
    "p_value_mean" = p_value_mean,

    "d_w_mean" = d_w_mean,
    "d_se_w_mean" = d_se_w_mean,
    "p_value_w_mean" = p_value_w_mean,

    "d_min" = d_min,
    "d_se_min" = d_se_min,
    "p_value_min" = p_value_min,

    "d_emmeans" = d_emmeans,
    "d_se_emmeans" = d_se_emmeans,
    "p_value_emmeans" = p_value_emmeans,

    "d_max" = d_max,
    "d_se_max" = d_se_max,
    "p_value_max" = p_value_max
  ))
}

r_cov_guess_grid <- c(0.1, 0.3, 0.5)
d_grid <- c(0, 0.2, 0.5, 0.8)
r_grid <- c(0.1, 0.3, 0.5)
n_grid <- c(25, 50, 75, 100, 200)
p_grid <- c(0.3, 0.5, 0.7)
# collect parameter grids in list:
param_list = list("d" = d_grid, "r_cov" = r_grid, "n" = n_grid, "p" = p_grid, "r_cov_guess" = r_cov_guess_grid)
MC_result <- MonteCarlo(func = means_adj_to_smd, nrep = 10,
                        ncpus = 1,
                        param_list = param_list, max_grid = 10000)
res <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/")
rio::export(res, "MEANS_ADJ_to_SMD_sim.txt")



dat = read.delim("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/MEANS_ADJ_to_SMD_sim.txt")
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
  summarise(d_raw_data = mean_na(d_raw_data),
            d_ancova = mean_na(d_ancova),
            r_raw_data = mean_na(r_raw_data))

method = c("pooled", "sigma", "mean", "w_mean", "min", "max", "emmeans")
for (met in method) {
  dat[, paste0("acc_es_", met)] <- abs(dat[, paste0("d_", met)] -
                                        dat$d_ancova)
  dat[, paste0("acc_se_", met)] <- abs(dat[, paste0("d_se_", met)] -
                                        dat$d_se_ancova)
  dat[, paste0("acc_pval_", met)] <- abs(dat[, paste0("p_value_", met)] -
                                          dat$p_value_ancova)
}
res = dat %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("acc_es") |
                 starts_with("acc_se") |
                 starts_with("acc_pval"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(d, n, r_cov, r_cov_guess, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("acc_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("acc_se", method_long, fixed = TRUE) ~ "_se",
                           grepl("acc_pval", method_long, fixed = TRUE) ~ "_pval"),
         method = case_when(grepl("sigma", method_long, fixed = TRUE) ~ "sigma",
                            grepl("mean", method_long) &
                              !grepl("w_mean", method_long, fixed = TRUE) &
                              !grepl("emmeans", method_long, fixed = TRUE) ~ "mean",
                            grepl("pooled", method_long, fixed = TRUE) ~ "pooled",
                            grepl("w_mean", method_long, fixed = TRUE) ~ "w_mean",
                            grepl("min", method_long, fixed = TRUE) ~ "min",
                            grepl("max", method_long, fixed = TRUE) ~ "max",
                            grepl("emmeans", method_long, fixed = TRUE) ~ "emmeans")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot = res

err_plot$d <- paste0("SMD=", err_plot$d)
err_plot$r_cov <- paste0("r_cov=", err_plot$r_cov)
err_plot$r_cov_guess <- paste0("guessed r_cov=", err_plot$r_cov_guess)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "200"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
err_plot$method[err_plot$method == "pooled"] <- "SD = pooled SDs"
err_plot$method[err_plot$method == "sigma"] <- "SD = sigma(model)"
err_plot$method[err_plot$method == "mean"] <- "SD = mean SDs"
err_plot$method[err_plot$method == "w_mean"] <- "SD = weighted mean SDs"
err_plot$method[err_plot$method == "min"] <- "SD = min(SDs)"
err_plot$method[err_plot$method == "max"] <- "SD = max(SDs)"
err_plot$method[err_plot$method == "emmeans"] <- "SD = emmean(model)"

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/MEANS_ADJ_to_SMD_AGG.txt")
