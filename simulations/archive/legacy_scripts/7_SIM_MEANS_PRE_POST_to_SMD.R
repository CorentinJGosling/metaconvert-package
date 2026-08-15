library(MonteCarlo); library(tidyverse); library(ggplot2)
sd_pooled = function(sd1, sd2, n1, n2) {
  sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))
}
.d_j <- function(x) {
  j <- ifelse(x <= 1, NA, 1) * exp(lgamma(x / 2) - 0.5 * log(x / 2) - lgamma((x - 1) / 2))
  return(j)
}
p=0.5;n=50;r_pre_post=0.5;d_posttest=1;d_pre=0.5;d_time_common=0.5; r_pre_post_guess=0.5
# https://cran.r-project.org/web/packages/TOSTER/vignettes/SMD_calcs.html#cohens-dz-change-scores

means_pre_post_to_smd <- function(d_posttest, d_pre,
                                  d_time_common,
                                  r_pre_post, r_pre_post_guess, n, p) {
  n_exp = round(n*p)
  n_nexp = n - n_exp
  sig <- rbind(c(1, r_pre_post),
               c(r_pre_post, 1))

  mu_exp <- c(d_pre, d_pre + d_time_common + d_posttest)
  mu_nexp <- c(0, d_time_common)
  # mu_exp <- c(d_pre, d_pre+d_posttest)
  # mu_nexp <- c(0, d_pre)

  res_exp <- as.data.frame(MASS::mvrnorm(n=n_exp, mu=mu_exp, Sig=sig))
  res_nexp <- as.data.frame(MASS::mvrnorm(n=n_nexp, mu=mu_nexp, Sig=sig))

  res = rbind(cbind(res_exp, group="g1"),
              cbind(res_nexp, group="g2"))

  mean_pre_exp <- mean(res_exp$V1)
  mean_pre_sd_exp <- sd(res_exp$V1)
  mean_exp <- mean(res_exp$V2)
  mean_sd_exp <- sd(res_exp$V2)

  mean_pre_nexp <- mean(res_nexp$V1)
  mean_pre_sd_nexp <- sd(res_nexp$V1)
  mean_nexp <- mean(res_nexp$V2)
  mean_sd_nexp <- sd(res_nexp$V2)

  r_raw_data = as.numeric(cor.test(~res$V1 + res$V2)$estimate)

  # Determine true SMD from raw data ===========================
  mod_wo = summary(lm(V2 ~ 0 + group + V1, res))
  sd_post = mod_wo$sigma #/ sqrt(1 - r_pre_post^2)
  coeff = coef(mod_wo)
  ancova_mean_exp = coeff[1,1]
  ancova_mean_nexp = coeff[2,1]
  d_raw_data = (ancova_mean_exp - ancova_mean_nexp) / sd_post
  J_exp <- .d_j(n_exp - 1)
  J_nexp <- .d_j(n_nexp - 1)

  # bonett
  var_exp_1 <- mean_pre_sd_exp^2 + mean_sd_exp^2 -
    2 * r_pre_post_guess * mean_pre_sd_exp * mean_sd_exp
  d_exp_1 <- (mean_pre_exp - mean_exp) / mean_pre_sd_exp
  g_exp_1 <- d_exp_1 * J_exp

  var_d_exp_1 <- var_exp_1 / (mean_pre_sd_exp^2 * (n_exp - 1)) + g_exp_1^2 / (2 * (n_exp - 1))

  d_nexp_1 <- (mean_pre_nexp - mean_nexp) / mean_pre_sd_nexp
  g_nexp_1 <- d_nexp_1 * J_nexp

  var_nexp_1 <- mean_pre_sd_nexp^2 + mean_sd_nexp^2 - 2 * r_pre_post_guess * mean_pre_sd_nexp * mean_sd_nexp
  var_d_nexp_1 <- var_nexp_1 / (mean_pre_sd_nexp^2 * (n_nexp - 1)) + g_nexp_1^2 / (2 * (n_nexp - 1))

  d_bonett <- d_exp_1 - d_nexp_1
  d_se_bonett <- sqrt(var_d_exp_1 + var_d_nexp_1)
  d_ci_lo_bonett <- d_bonett - d_se_bonett * qt(.975, n_exp + n_nexp - 2)
  d_ci_up_bonett <- d_bonett + d_se_bonett * qt(.975, n_exp + n_nexp - 2)

  # cooper
  sd_diff_exp <- sqrt(mean_pre_sd_exp^2 + mean_sd_exp^2 -
                        (2 * r_pre_post_guess * mean_pre_sd_exp * mean_sd_exp))

  sd_diff_nexp <- sqrt(mean_pre_sd_nexp^2 + mean_sd_nexp^2 -
                         (2 * r_pre_post_guess * mean_pre_sd_nexp * mean_sd_nexp))
  d_exp_2 <- (mean_pre_exp - mean_exp) / sd_diff_exp * sqrt(2 * (1 - r_pre_post_guess))

  d_nexp_2 <- (mean_pre_nexp - mean_nexp) / sd_diff_nexp * sqrt(2 * (1 - r_pre_post_guess))

  d_var_exp_2 <- 2 * (1 - r_pre_post_guess)/n_exp + d_exp_2^2 / (2*n_exp)
  d_var_nexp_2 <- 2 * (1 - r_pre_post_guess)/n_nexp + d_nexp_2^2 / (2*n_nexp)

  d_cooper <- d_exp_2 - d_nexp_2
  d_se_cooper <- sqrt(d_var_exp_2 + d_var_nexp_2)
  d_ci_lo_cooper <- d_cooper - d_se_cooper * qt(.975, n_exp + n_nexp - 2)
  d_ci_up_cooper <- d_cooper + d_se_cooper * qt(.975, n_exp + n_nexp - 2)

  # dz
  d_z_exp = (mean_pre_exp - mean_exp) / sd_diff_exp
  d_var_z_exp = 1/n_exp + (d_z_exp^2)/(2*n_exp)

  d_z_nexp = (mean_pre_nexp - mean_nexp) / sd_diff_nexp
  d_var_z_nexp = 1/n_nexp + (d_z_nexp^2)/(2*n_nexp)

  d_z = d_z_exp - d_z_nexp
  d_se_z = sqrt(d_var_z_exp + d_var_z_nexp)
  d_ci_lo_z <- d_z - qt(.975, df = n_exp + n_nexp - 2) * d_se_z
  d_ci_up_z <- d_z + qt(.975, df = n_exp + n_nexp - 2) * d_se_z

  # dav
  d_av_exp = (mean_pre_exp - mean_exp) / mean(c(mean_sd_exp, mean_pre_sd_exp))
  d_var_av_exp = 1/n_exp + (d_av_exp^2)/(2*n_exp)

  d_av_nexp = (mean_pre_nexp - mean_nexp) / mean(c(mean_sd_nexp, mean_pre_sd_nexp))
  d_var_av_nexp = 1/n_nexp + (d_av_nexp^2)/(2*n_nexp)

  d_av = d_av_exp - d_av_nexp
  d_se_av = sqrt(d_var_av_exp + d_var_av_nexp)
  d_ci_lo_av <- d_av - qt(.975, df = n_exp + n_nexp - 2) * d_se_av
  d_ci_up_av <- d_av + qt(.975, df = n_exp + n_nexp - 2) * d_se_av

  # dpost
  d_post <- (mean_exp - mean_nexp) / sqrt(((n_exp-1)*mean_sd_exp^2 + (n_nexp-1)*mean_sd_nexp^2) / (n_exp + n_nexp - 2))
  d_se_post = sqrt( (1/n_exp + 1/n_nexp)  + d_post^2/(2*(n_exp + n_nexp)) )
  d_ci_lo_post <- d_post - qt(.975, df = n_exp + n_nexp - 2) * d_se_post
  d_ci_up_post <- d_post + qt(.975, df = n_exp + n_nexp - 2) * d_se_post


  return(list(
    "r_raw_data" = r_raw_data,

    "d_raw_data" = d_raw_data,

    "d_bonett" = d_bonett,
    "d_se_bonett" = d_se_bonett,
    "d_ci_lo_bonett" = d_ci_lo_bonett,
    "d_ci_up_bonett" = d_ci_up_bonett,

    "d_z" = d_z,
    "d_se_z" = d_se_z,
    "d_ci_lo_z" = d_ci_lo_z,
    "d_ci_up_z" = d_ci_up_z,

    "d_av" = d_av,
    "d_se_av" = d_se_av,
    "d_ci_lo_av" = d_ci_lo_av,
    "d_ci_up_av" = d_ci_up_av,

    "d_cooper" = d_cooper,
    "d_se_cooper" = d_se_cooper,
    "d_ci_lo_cooper" = d_ci_lo_cooper,
    "d_ci_up_cooper" = d_ci_up_cooper,

    "d_post" = d_post,
    "d_se_post" = d_se_post,
    "d_ci_lo_post" = d_ci_lo_post,
    "d_ci_up_post" = d_ci_up_post

    ))
}
r_pre_post_grid <- c(0.3, 0.5, 0.7)
r_pre_post_guess_grid <- c(0.3, 0.5, 0.7)
d_pre_grid <- c(0, 0.4, 0.8)
d_time_common_grid <- c(0, 0.4, 0.8)
d_post_grid <- c(0, 0.2, 0.5, 0.8)
n_grid <- c(25, 50, 75, 100, 300)
p_grid <- c(0.3, 0.5, 0.7)
# collect parameter grids in list:
param_list = list("r_pre_post" = r_pre_post_grid,
                  "r_pre_post_guess" = r_pre_post_guess_grid,
                  "d_pre" = d_pre_grid,
                  "d_posttest" = d_post_grid,
                  "d_time_common" = d_time_common_grid,
                  "n" = n_grid,
                  "p" = p_grid)

MC_result <- MonteCarlo(func = means_pre_post_to_smd,
                        nrep = 100,
                        time_n_test = TRUE,
                        param_list = param_list,
                        max_grid = 10000)
datsim <- MakeFrame(MC_result)

rio::export(datsim, "D:/simulations/data/MEANS_PRE_POST_to_SMD_sim1000a.txt")

dat = rbind(
  read.delim("D:/simulations/data/MEANS_PRE_POST_to_SMD_sim100a.txt"))
  # read.delim("D:/simulations/data/MEANS_ADJ_to_SMD_sim2500b.txt"),
  # read.delim("D:/simulations/data/MEANS_ADJ_to_SMD_sim2500c.txt"),
  # read.delim("D:/simulations/data/MEANS_ADJ_to_SMD_sim2500d.txt"))

mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}

checks = dat %>%
  group_by(d_poscooper, r_pre_post) %>%
  summarise(d_raw_data = mean(d_raw_data),
            r_raw_data = mean(r_raw_data))

method = c("bonett", "cooper", "post", "z", "av")
for (met in method) {
  dat[, paste0("bias_es_d_", met)] <- dat[, paste0("d_", met)] - dat$d_raw_data
  dat[, paste0("d_var_", met)] = dat[, paste0("d_se_", met)]^2

  dat[, paste0("ci_cov_d_", met)] = dat[, paste0("d_ci_lo_", met)] <= dat$d_raw_data &
                                    dat[, paste0("d_ci_up_", met)] >= dat$d_raw_data
}

res = dat %>%
  group_by(d_posttest, d_pre, n, r_pre_post, p, r_pre_post_guess) %>%
  summarise(
    n_sim = n(),
    bias_es_d_bonett = mean(bias_es_d_bonett),
    var_d_bonett = var(d_bonett),
    mean_d_var_bonett = mean(d_var_bonett),
    bias_ci_d_bonett = mean(ci_cov_d_bonett),

    bias_es_d_post = mean(bias_es_d_post),
    var_d_post = var(d_post),
    mean_d_var_post = mean(d_var_post),
    bias_ci_d_post = mean(ci_cov_d_post),

    bias_es_d_z = mean(bias_es_d_z),
    var_d_z = var(d_z),
    mean_d_var_z = mean(d_var_z),
    bias_ci_d_z = mean(ci_cov_d_z),

    bias_es_d_cooper = mean(bias_es_d_cooper),
    var_d_cooper = var(d_cooper),
    mean_d_var_cooper = mean(d_var_cooper),
    bias_ci_d_cooper = mean(ci_cov_d_cooper)
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
  group_by(d_posttest, d_pre, n, r_pre_post, p, r_pre_post_guess, method_long) %>%
  summarise(bias = mean(bias)) %>%
  mutate(pivot = case_when(grepl("bias_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("bias_var", method_long, fixed = TRUE) ~ "_var",
                           grepl("bias_ci", method_long, fixed = TRUE) ~ "_ci"),
         method = case_when(grepl("bonett", method_long, fixed = TRUE) ~ "bonett",
                            grepl("post", method_long, fixed = TRUE) ~ "post",
                            grepl("z", method_long, fixed = TRUE) ~ "z",
                            grepl("cooper", method_long, fixed = TRUE) ~ "cooper")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot$d <- paste0("d=", err_plot$d_posttest)
err_plot$p <- paste0("% of exposed=", err_plot$p)
err_plot$r_pre_post <- paste0("r_pre_post=", err_plot$r_pre_post)
err_plot$r_pre_post_guess <- paste0("guessed r_pre_post=", err_plot$r_pre_post_guess)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "300"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=300"))
err_plot$method[err_plot$method == "bonett"] <- "Bonett"
err_plot$method[err_plot$method == "post"] <- "Post-test"
err_plot$method[err_plot$method == "cooper"] <- "Cooper (Cohen's drm)"
err_plot$method[err_plot$method == "z"] <- "Cohen's dz"
err_plot$method[err_plot$method == "z"] <- "Cumming (Cohen's av)"
err_plot$bi <- paste0("bi (d=", err_plot$d_pre, ")")

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/metaConvert/simulations/data_agg/MEANS_PRE_POST_to_SMD_AGG2500a.txt")
