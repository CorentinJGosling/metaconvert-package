library(MonteCarlo); library(tidyverse); library(ggplot2)
d=0
r=0
r_guess=0
n=25
p=0.3
means_adj_to_smd <- function(d, r, r_guess, n, p) {
  # print(paste0(d, "_", r, "_", r_guess, "_", n, "_",p))
  n1 = round(n*p)
  n2 = n - n1
  vd1=rnorm(n1, d, 1)
  vd2=rnorm(n2, 0, 1)
  grp = factor(rep(c(0,1), c(n1,n2)))
  vd = c(vd1, vd2)
  # length(vd)
  cor = 1/r
  cov1 = vd + rnorm(n, 0, cor)
  # cov2 = vd + rnorm(n, 0, 1/r)
  # cov3 = vd + rnorm(n, 0, 1/r)
  # cov4 = vd + rnorm(n, 0, 1/r)
  #
  # mod = summary(lm(vd~cov1+cov2+cov3+cov4))
  #
  m1 <- mean(vd1)
  m2 <- mean(vd2)
  sd1 <- sd(vd1)
  sd2 <- sd(vd2)

  mod = lm(vd~cov1)
  mod_adj = lm(vd~grp+cov1)
  r_raw_data = sqrt(summary(mod)$r.squared)
  p_value_raw = summary(mod_adj)$coefficients[2,4]

  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  d_crude <- (m1 - m2) / pooled_sd

  res = data.frame(emmeans::emmeans(mod_adj, ~grp))
  em_mean1 = res$emmean[1]
  em_mean2 = res$emmean[2]
  sd_adj1 = res$SE[1] * sqrt(n1 - 1) / sqrt(1 - r_guess^2)
  sd_adj2 = res$SE[2] * sqrt(n2 - 1) / sqrt(1 - r_guess^2)
  sd_adj3 = res$SE[1] * sqrt(n1 - 1)
  sd_adj4 = res$SE[2] * sqrt(n2 - 1)

  #  raw MD
  md = em_mean1 - em_mean2
  md_se = sqrt(sd1/n1 + sd2/n2)
  p_value = 1 - 2 * abs(pt(md / md_se, n1+n2-2-1) - 0.5)

  #  estimated SD
  md_est1 = em_mean1 - em_mean2
  md_se_est1 = sqrt(sd_adj1/n1 * (1 - r_guess^2) + sd_adj2/n2 * (1 - r_guess^2))
  p_value_est1 = 1 - 2 * abs(pt(md_est1 / md_se_est1, n1+n2-2-1) - 0.5)

  #  estimated SD
  md_est2 = em_mean1 - em_mean2
  md_se_est2 = sqrt(sd_adj3/n1 * (1 - r_guess^2) + sd_adj4/n2 * (1 - r_guess^2))
  p_value_est2 = 1 - 2 * abs(pt(md_est2 / md_se_est2, n1+n2-2-1) - 0.5)

  # return result:
  return(list(
    "p_value_ancova" = p_value_raw,
    "d_raw_data" = d_crude,
    "md_raw_data" = md,
    "md_se_raw_data" = md_se,
    "p_value_raw_data" = p_value,
    "r_raw_data" = r_raw_data,
    "md_est1" = md_est1,
    "md_se_est1" = md_se_est1,
    "p_value_est1" = p_value_est1,
    "md_est2" = md_est2,
    "md_se_est2" = md_se_est2,
    "p_value_est2" = p_value_est2

  ))
}

r_guess_grid <- c(0.1, 0.3, 0.5)
d_grid <- c(0, 0.2, 0.5, 0.8)
r_grid <- c(0.1, 0.3, 0.5)
n_grid <- c(25, 50, 100, 200)
p_grid <- c(0.3, 0.5, 0.7)
# collect parameter grids in list:
param_list = list("d" = d_grid,"r" = r_grid, "n" = n_grid, "p" = p_grid, "r_guess" = r_guess_grid)
MC_result <- MonteCarlo(func = means_adj_to_smd, nrep = 10, param_list = param_list, max_grid = 10000)
df <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/Documents/es.utils/simulations/data")
rio::export(df, "MEANS_ADJ_to_MD_sim.txt")
