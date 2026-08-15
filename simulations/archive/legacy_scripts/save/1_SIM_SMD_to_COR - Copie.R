library(MonteCarlo); library(tidyverse); library(ggplot2)

smd_to_cor <- function(r, n, p) {

  # generate sample scores
  sig <- rbind(c(1, r), c(r, 1))

  # create the mean vector
  mu <- c(0, 0)

  # generate the multivariate normal distribution
  res <- as.data.frame(MASS::mvrnorm(n=n, mu=mu, Sigma=sig))
  r_raw_data = as.numeric(cor.test(~res$V1 + res$V2)$estimate)
  v.r_raw_data =  (1-r_raw_data^2)^2 / (n-1)
  r_raw_data.pval = 1 - 2 * abs(pt(r_raw_data / sqrt(v.r_raw_data), n - 2) - 0.5)
  z_raw_data = atanh(r_raw_data)
  v.z_raw_data = 1/(n - 3)
  z_raw_data.pval = 1 - 2 * abs(pnorm(z_raw_data / sqrt(v.z_raw_data)) - 0.5)

  # obtain SMD from categorisation of the V1
  x = 1-p
  cut = quantile(res$V1, x)
  n_exp = length(res$V1[res$V1 >= cut])
  n_nexp = length(res$V1[res$V1 < cut])
  scores_grp1 <- res$V2[res$V1 >= cut]
  scores_grp2 <- res$V2[res$V1 < cut]
  mean_exp <- mean(scores_grp1)
  mean_nexp <- mean(scores_grp2)
  mean_sd_exp <- sd(scores_grp1)
  mean_sd_nexp <- sd(scores_grp2)
  pooled_sd <- sqrt(((n_exp - 1) * mean_sd_exp^2 + (n_nexp - 1) * mean_sd_nexp^2) / (n_exp + n_nexp - 2))

  # ==========================================================
  # Estimate d from raw data =================================
  # ==========================================================
  d <- (mean_exp - mean_nexp) / pooled_sd
  vd = sqrt((n_exp+n_nexp)/(n_exp*n_nexp) + d^2/(2*(n_exp+n_nexp)))
  pval.d = 1 - 2 * abs(pt(d / sqrt(vd), n_exp + n_nexp - 2) - 0.5)

  # =========================================================================
  # 1. Convert d to r according to VIECHTBAUER ==============================
  # =========================================================================
  df =  n_exp + n_nexp - 2
  h <- df/n_exp + df/n_nexp
  p <- n_exp / (n_exp + n_nexp)
  q <- n_nexp / (n_exp + n_nexp)
  r_pb <- d / sqrt(d^2 + h)
  f  <- dnorm(qnorm(p, lower.tail=FALSE))
  rb  <- sqrt(p*q) / f * r_pb
  r_viechtbauer <- ifelse(rb > 1, 1, ifelse(rb < -1, -1, rb))
  vr_viechtbauer = 1/(n_exp + n_nexp - 1) * (p*q/f^2 - (3/2 + (1 - p*qnorm(p, lower.tail=FALSE)/f)*(1 + q*qnorm(p, lower.tail=FALSE)/f)) * r_viechtbauer^2 + r_viechtbauer^4)
  fzp <- dnorm(qnorm(p))
  a_viechtbauer <- sqrt(fzp) / (p*(1-p))^(1/4)
  zb = (a_viechtbauer/2) * log((1+a_viechtbauer*r_viechtbauer)/(1-a_viechtbauer*r_viechtbauer))

  z_viechtbauer <- ifelse(zb > 1, 1, ifelse(zb < -1, -1, zb))
  vz_viechtbauer = 1/(n_exp + n_nexp - 1)

  # ===========================================================================
  # 2. Convert d to r according to COOPER/LIPSEY ==============================
  # ===========================================================================
  a = ((n_exp + n_nexp)^2) / (n_exp*n_nexp)
  p = n_exp / (n_exp + n_nexp)

  r_lipsey = d/sqrt(d^2 + 1/(p * (1 - p)))
  vr_lipsey = a^2 * vd / ((d^2 + a)^3)

  z_lipsey = atanh(r_lipsey)
  vz_lipsey = vd / (vd + 1/(p*(1-p)))

# return result:
return(list(
  "d_raw_data" = d,
  "r_raw_data" = r_raw_data,
  "z_raw_data" = z_raw_data,
  "r_se_raw_data" = sqrt(v.r_raw_data),
  "z_se_raw_data" = sqrt(v.z_raw_data),
  "r_pval_raw_data" = r_raw_data.pval,
  "z_pval_raw_data" = z_raw_data.pval,
  # es
  "r_lipsey_cooper"     = r_lipsey,
  "z_lipsey_cooper"     = z_lipsey,
  "r_viecht"            = r_viechtbauer,
  "z_viecht"            = z_viechtbauer,
  # se
  "r_se_lipsey_cooper"     = sqrt(vr_lipsey),
  "z_se_lipsey_cooper"     = sqrt(vz_lipsey),
  "r_se_viecht"            = sqrt(vr_viechtbauer),
  "z_se_viecht"            = sqrt(vz_viechtbauer),
  # pval
  "r_pval_lipsey_cooper"     = 1 - 2 * abs(pt(r_lipsey / sqrt(vr_lipsey), n_exp + n_nexp - 2) - 0.5),
  "z_pval_lipsey_cooper"     = 1 - 2 * abs(pnorm(z_lipsey / sqrt(vz_lipsey)) - 0.5),
  "r_pval_viecht"            = 1 - 2 * abs(pt(r_viechtbauer / sqrt(vr_viechtbauer), n_exp + n_nexp - 2) - 0.5),
  "z_pval_viecht"            = 1 - 2 * abs(pnorm(z_viechtbauer / sqrt(vz_viechtbauer)) - 0.5)))
}

r_grid <- c(0, 0.1, 0.2, 0.3, 0.5)
n_grid <- c(25, 50, 75, 100, 200)
p_grid <- c(0.3, 0.5, 0.7)

# collect parameter grids in list:
param_list = list("r" = r_grid, "n" = n_grid, "p" = p_grid)
MC_result <- MonteCarlo(func = smd_to_cor, nrep = 500, param_list = param_list)
df <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data")
rio::export(df, "SMD_to_COR_sim.txt")

#=====================================================================#
#==========================AGG========================================#
#========================DATASET======================================#
#=====================================================================#
res_smd = read.delim("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/SMD_to_COR_sim.txt") %>%
  mutate_if(is.character, as.numeric)

checks = res_smd %>%
  group_by(r) %>%
  summarise(sd_r_raw_data = sd(r_raw_data),
            r_raw_data = mean(r_raw_data),
            sd_r_viecht = sd(r_viecht),
            r_viecht = mean(r_viecht))

method = c("viecht", "lipsey_cooper")
for (met in method) {
  res_smd[, paste0("acc_es_r_", met)] <- abs(res_smd[, paste0("r_", met)] -
                                               res_smd$r_raw_data)
  res_smd[, paste0("acc_se_r_", met)] <- abs(res_smd[, paste0("r_se_", met)] -
                                               res_smd$r_se_raw_data)
  res_smd[, paste0("acc_pval_r_", met)] <- abs(res_smd[, paste0("r_pval_", met)] -
                                                 res_smd$r_pval_raw)
  res_smd[, paste0("acc_es_z_", met)] <- abs(res_smd[, paste0("z_", met)] -
                                               res_smd$z_raw_data)
  res_smd[, paste0("acc_se_z_", met)] <- abs(res_smd[, paste0("z_se_", met)] -
                                               res_smd$z_se_raw_data)
  res_smd[, paste0("acc_pval_z_", met)] <- abs(res_smd[, paste0("z_pval_", met)] -
                                                 res_smd$z_pval_raw)
}
res = res_smd %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("acc_es") |
                 starts_with("acc_se") |
                 starts_with("acc_pval"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(r, n, p, method_long) %>%
  summarise(bias = mean_na(bias),
            n_sim = n()) %>%
  mutate(pivot = case_when(grepl("acc_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("acc_se", method_long, fixed = TRUE) ~ "_se",
                           grepl("acc_pval", method_long, fixed = TRUE) ~ "_pval"),
  method = case_when(grepl("r_lipsey", method_long, fixed = TRUE) ~ "r_lipsey",
                     grepl("z_lipsey", method_long, fixed = TRUE) ~ "z_lipsey",
                     grepl("r_viecht", method_long, fixed = TRUE) ~ "r_viecht",
                     grepl("z_viecht", method_long, fixed = TRUE) ~ "z_viecht")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot = res

err_plot$method[err_plot$method == "r_lipsey"] <- "Lipsey-Cooper (r)"
err_plot$method[err_plot$method == "z_lipsey"] <- "Lipsey-Cooper (z)"
err_plot$method[err_plot$method == "r_viecht"] <- "Viechtbauer (r)"
err_plot$method[err_plot$method == "z_viecht"] <- "Viechtbauer (z)"

err_plot$r <- paste0("r=", err_plot$r)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "200"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=200"))

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/SMD_to_COR_AGG.txt")
