library(MonteCarlo); library(tidyverse); library(ggplot2)

smd_to_cor <- function(r, n, p) {

  sig <- rbind(c(1, r), c(r, 1))

  mu <- c(0, 0)

  res <- as.data.frame(MASS::mvrnorm(n=n, mu=mu, Sigma=sig))
  r_raw_data = as.numeric(cor.test(~res$V1 + res$V2)$estimate)
  z_raw_data = atanh(r_raw_data)

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
  df =  n_exp + n_nexp - 2
  pooled_sd <- sqrt(((n_exp - 1) * mean_sd_exp^2 + (n_nexp - 1) * mean_sd_nexp^2) / df)

  d_raw_data = d = (mean_exp - mean_nexp) / pooled_sd
  d_se_raw_data = sqrt((n_exp+n_nexp)/(n_exp*n_nexp) + d^2/(2*(n_exp+n_nexp)))
  vd = d_se_raw_data^2
  h <- df/n_exp + df/n_nexp
  p <- n_exp / (n_exp + n_nexp)
  q <- n_nexp / (n_exp + n_nexp)
  r_pb <- d / sqrt(d^2 + h)
  f  <- dnorm(qnorm(p, lower.tail=FALSE))
  r_viechtbauer <- sqrt(p * q) / f * r_pb
  r_trunc = ifelse(r_viechtbauer > 1, 1, ifelse(r_viechtbauer < -1, -1, r_viechtbauer))
  r_se_viechtbauer <- sqrt(1 / (n_exp + n_nexp - 1) *
    (p * q / f^2 - (3 / 2 + (1 - p * qnorm(p, lower.tail = FALSE) / f) *
                      (1 + q * qnorm(p, lower.tail = FALSE) / f)) *
       r_trunc^2 + r_trunc^4))

  fzp <- dnorm(qnorm(p))
  a_viechtbauer <- sqrt(fzp) / (p*(1-p))^(1/4)
  z_viechtbauer = (a_viechtbauer/2) * log((1+a_viechtbauer*r_trunc)/
                                          (1-a_viechtbauer*r_trunc))
  z_se_viechtbauer = sqrt(1/(n_exp + n_nexp - 1))
  z_ci_lo_viechtbauer <- z_viechtbauer - qnorm(.975) * z_se_viechtbauer
  z_ci_up_viechtbauer <- z_viechtbauer + qnorm(.975) * z_se_viechtbauer
  r_ci_lo_viechtbauer <- (1/a_viechtbauer) * ((exp(2*z_ci_lo_viechtbauer / a_viechtbauer) - 1) / (exp(2*z_ci_lo_viechtbauer / a_viechtbauer) + 1))
  r_ci_up_viechtbauer <- (1/a_viechtbauer) * ((exp(2*z_ci_up_viechtbauer / a_viechtbauer) - 1) / (exp(2*z_ci_up_viechtbauer / a_viechtbauer) + 1))

  a <- ((n_exp + n_nexp)^2) / (n_exp * n_nexp)
  p <- n_exp / (n_exp + n_nexp)
  r_lipsey <- d / sqrt(d^2 + 1 / (p * (1 - p)))
  r_se_lipsey <- sqrt(a^2 * vd / ((d^2 + a)^3))
  z_lipsey <- atanh(r_lipsey)
  z_se_lipsey <- sqrt(vd / (vd + 1 / (p * (1 - p))))
  r_ci_lo_lipsey <- r_lipsey - qt(.975, df = n_exp + n_nexp - 2) * r_se_lipsey
  r_ci_up_lipsey <- r_lipsey + qt(.975, df = n_exp + n_nexp - 2) * r_se_lipsey
  z_ci_lo_lipsey <- z_lipsey - qnorm(.975) * z_se_lipsey
  z_ci_up_lipsey <- z_lipsey + qnorm(.975) * z_se_lipsey

return(list(
  "d_raw_data"    = d,
  "r_raw_data"    = r_raw_data,
  "z_raw_data"    = z_raw_data,

  "r_lipsey_cooper"           = r_lipsey,
  "r_se_lipsey_cooper"        = r_se_lipsey,
  "r_ci_lo_lipsey_cooper"     = r_ci_lo_lipsey,
  "r_ci_up_lipsey_cooper"     = r_ci_up_lipsey,

  "z_lipsey_cooper"           = z_lipsey,
  "z_se_lipsey_cooper"        = z_se_lipsey,
  "z_ci_lo_lipsey_cooper"     = z_ci_lo_lipsey,
  "z_ci_up_lipsey_cooper"     = z_ci_up_lipsey,

  "r_viechtbauer"             = r_viechtbauer,
  "r_se_viechtbauer"          = r_se_viechtbauer,
  "r_ci_lo_viechtbauer"       = r_ci_lo_viechtbauer,
  "r_ci_up_viechtbauer"       = r_ci_up_viechtbauer,

  "z_viechtbauer"             = z_viechtbauer,
  "z_se_viechtbauer"          = z_se_viechtbauer,
  "z_ci_lo_viechtbauer"       = z_ci_lo_viechtbauer,
  "z_ci_up_viechtbauer"       = z_ci_up_viechtbauer))
}

r_grid <- c(0, 0.1, 0.25, 0.5, 0.75)
n_grid <- c(25, 50, 75, 100, 300)
p_grid <- c(0.3, 0.5, 0.7)

# collect parameter grids in list:
param_list = list("r" = r_grid, "n" = n_grid, "p" = p_grid)
MC_result <- MonteCarlo(func = smd_to_cor, nrep = 10000, param_list = param_list)
datasim <- MakeFrame(MC_result)

rio::export(datasim, "D:/simulations/data/SMD_to_COR_sim.txt")

#=====================================================================#
#==========================AGG========================================#
#========================DATASET======================================#
#=====================================================================#
res_smd = read.delim("D:/simulations/data/SMD_to_COR_sim.txt") %>%
  mutate_if(is.character, as.numeric)
mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}

breks = res_smd %>%
  group_by(r) %>%
  summarise(r_raw_data = mean_na(r_raw_data))
ggplot(res_smd, aes(x = r_raw_data,
                   group=as.character(r),
                   fill=as.character(r))) +
  geom_density(alpha=0.5) +
  geom_vline(data=breks, aes(xintercept = r_raw_data),
             color="red") +
  scale_x_continuous(name="r",
                     breaks=breks$r_raw_data,labels=paste0(round(breks$r_raw_data,2)))

checks = res_smd %>%
  group_by(r) %>%
  summarise(sd_r_raw_data = sd(r_raw_data),
            r_raw_data = mean(r_raw_data))

for (met in c("lipsey_cooper", "viechtbauer")) {

  res_smd[, paste0("bias_es_r_", met)] <- res_smd[, paste0("r_", met)] - res_smd$r_raw_data
  res_smd[, paste0("bias_es_z_", met)] <- res_smd[, paste0("z_", met)] - res_smd$z_raw_data

  res_smd[, paste0("r_var_", met)] = res_smd[, paste0("r_se_", met)]^2
  res_smd[, paste0("z_var_", met)] = res_smd[, paste0("z_se_", met)]^2

  res_smd[, paste0("ci_cov_r_", met)] = res_smd[, paste0("r_ci_lo_", met)] <= res_smd$r_raw_data &
                                        res_smd[, paste0("r_ci_up_", met)] >= res_smd$r_raw_data
  res_smd[, paste0("ci_cov_z_", met)] = res_smd[, paste0("z_ci_lo_", met)] <= res_smd$z_raw_data &
                                        res_smd[, paste0("z_ci_up_", met)] >= res_smd$z_raw_data
}

res = res_smd %>%
  group_by(r, p, n) %>%
  summarise(
    n_sim = n(),
    bias_es_r_viechtbauer = mean(bias_es_r_viechtbauer),
    var_r_viechtbauer = var(r_viechtbauer),
    mean_r_var_viechtbauer = mean(r_var_viechtbauer),
    bias_ci_r_viechtbauer = mean(ci_cov_r_viechtbauer),

    bias_es_z_viechtbauer = mean(bias_es_z_viechtbauer),
    var_z_viechtbauer = var(z_viechtbauer),
    mean_z_var_viechtbauer = mean(z_var_viechtbauer),
    bias_ci_z_viechtbauer = mean(ci_cov_z_viechtbauer),

    bias_es_r_lipsey_cooper = mean(bias_es_r_lipsey_cooper),
    bias_es_z_lipsey_cooper = mean(bias_es_z_lipsey_cooper),

    var_r_lipsey_cooper = var(r_lipsey_cooper),
    mean_r_var_lipsey_cooper = mean(r_var_lipsey_cooper),
    var_z_lipsey_cooper = var(z_lipsey_cooper),
    mean_z_var_lipsey_cooper = mean(z_var_lipsey_cooper),

    bias_ci_r_lipsey_cooper = mean(ci_cov_r_lipsey_cooper),
    bias_ci_z_lipsey_cooper = mean(ci_cov_z_lipsey_cooper)
  )

for (met in c("viechtbauer", "lipsey_cooper")) {
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
  group_by(r, n, p, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("bias_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("bias_var", method_long, fixed = TRUE) ~ "_var",
                           grepl("bias_ci", method_long, fixed = TRUE) ~ "_ci"),
  method = case_when(grepl("r_lipsey", method_long, fixed = TRUE) ~ "r_lipsey",
                     grepl("z_lipsey", method_long, fixed = TRUE) ~ "z_lipsey",
                     grepl("r_viecht", method_long, fixed = TRUE) ~ "r_viechtbauer",
                     grepl("z_viecht", method_long, fixed = TRUE) ~ "z_viechtbauer")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot$method[err_plot$method == "r_lipsey"] <- "Lipsey-Cooper (r)"
err_plot$method[err_plot$method == "z_lipsey"] <- "Lipsey-Cooper (z)"
err_plot$method[err_plot$method == "r_viechtbauer"] <- "Viechtbauer (r)"
err_plot$method[err_plot$method == "z_viechtbauer"] <- "Viechtbauer (z)"

err_plot$r <- paste0("r=", err_plot$r)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "300"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=300"))

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/metaConvert/simulations/data_agg/SMD_to_COR_AGG.txt")


View(read.delim("C:/Users/Corentin Gosling/drive_gmail/Recherche/metaConvert/simulations/data_agg/SMD_to_COR_AGG.txt"))









