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
  # Conversion 2. Lipsey =========================================
  # ==============================================================
  r_lipsey = (n_cases_exp*n_controls_nexp - n_controls_exp*n_cases_nexp) /
    sqrt((n_cases_exp+ n_controls_exp)*(n_cases_nexp+n_controls_nexp)*
           (n_cases_exp+n_cases_nexp)*(n_controls_exp+n_cases_nexp))
  z_lipsey = atanh(r_lipsey)
  z_se_lipsey = sqrt(or_se_raw_data^2 * (z_lipsey^2) / (log(or_raw_data)^2))

  z_lo_lipsey = z_lipsey - qnorm(.975)*z_se_lipsey
  z_up_lipsey = z_lipsey + qnorm(.975)*z_se_lipsey
  r_lo_lipsey = tanh(z_lo_lipsey)
  r_up_lipsey = tanh(z_up_lipsey)

  r_se_lipsey = (r_up_lipsey - r_lo_lipsey)/(2*qt(.975, n_sample - 2))
  r_pval_lipsey = 1 - 2 * abs(pt(r_lipsey / r_se_lipsey, n_sample - 2) - 0.5)
  z_pval_lipsey = 1 - 2 * abs(pnorm(z_lipsey / z_se_lipsey) - 0.5)


  # ==============================================================
  # Conversion 3. Tetrachoric correlation  =======================
  # ==============================================================
  r_tetrachoric = v.r_tetrachoric = NA
  tet = metafor:::.rtet(ai = n_cases_exp, bi= n_controls_exp,
                        ci = n_cases_nexp, di= n_controls_nexp,
                        maxcor = 0.9999)
  r_tetrachoric = as.numeric(tet$yi)
  r_se_tetrachoric = sqrt(as.numeric(tet$vi))
  z_tetrachoric = atanh(r_tetrachoric)
  z_se_tetrachoric = sqrt(r_se_tetrachoric^2 / ((1 - r_tetrachoric^2)^2))
  r_pval_tetrachoric = 1 - 2 * abs(pt(r_tetrachoric / r_se_tetrachoric, n_sample - 2) - 0.5)
  z_pval_tetrachoric = 1 - 2 * abs(pnorm(-z_tetrachoric / z_se_tetrachoric) - 0.5)

  # return result:
  return(list(
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
    "r_lipsey"           = r_lipsey,
    "z_lipsey"           = z_lipsey,
    "r_tetrachoric"      = r_tetrachoric,
    "z_tetrachoric"      = z_tetrachoric,
    "r_cooper"           = r_cooper,
    "z_cooper"           = z_cooper,
    # se
    "r_se_lipsey"         = r_se_lipsey,
    "z_se_lipsey"         = z_se_lipsey,
    "r_se_tetrachoric"    = r_se_tetrachoric,
    "z_se_tetrachoric"    = z_se_tetrachoric,
    "r_se_cooper"         = r_se_cooper,
    "z_se_cooper"         = z_se_cooper,
    # pval
    "r_pval_cooper"       = r_pval_cooper,
    "z_pval_cooper"       = z_pval_cooper,
    "r_pval_lipsey"       = r_pval_lipsey,
    "z_pval_lipsey"       = z_pval_lipsey,
    "r_pval_tetrachoric"  = r_pval_tetrachoric,
    "z_pval_tetrachoric"  = z_pval_tetrachoric
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
MC_result <- MonteCarlo(func = cont_to_z, nrep = 50, param_list = param_list, max_grid = 10000)
df <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data")
rio::export(df, "2x2_to_COR_sim.txt")


res_cont = read.delim("2x2_to_COR_sim.txt")
res_cont[res_cont==Inf|res_cont==-Inf] <- NA

checks = res_cont %>%
  group_by(r) %>%
  summarise(r_raw_data = mean_na(r_raw_data))

method = c("lipsey", "tetrachoric", "cooper")
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
          method = case_when(grepl("r_lipsey", method_long, fixed = TRUE) ~ "r_lipsey",
                             grepl("z_lipsey", method_long, fixed = TRUE) ~ "z_lipsey",
                             grepl("r_tetrachoric", method_long, fixed = TRUE) ~ "r_tetrachoric",
                             grepl("z_tetrachoric", method_long, fixed = TRUE) ~ "z_tetrachoric",
                             grepl("r_cooper", method_long, fixed = TRUE) ~ "r_cooper",
                             grepl("z_cooper", method_long, fixed = TRUE) ~ "z_cooper")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot = res

err_plot$method[err_plot$method == "r_cooper"] <- "Cooper (r)"
err_plot$method[err_plot$method == "z_cooper"] <- "Cooper (z)"
err_plot$method[err_plot$method == "r_lipsey"] <- "Lipsey (r)"
err_plot$method[err_plot$method == "z_lipsey"] <- "Lipsey (z)"
err_plot$method[err_plot$method == "r_tetrachoric"] <- "Tetrachoric (r)"
err_plot$method[err_plot$method == "z_tetrachoric"] <- "Tetrachoric (z)"

err_plot$r <- paste0("r=", err_plot$r)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "200"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=200"))
rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/simulations/data/2x2_to_COR_AGG.txt")
