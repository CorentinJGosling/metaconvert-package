library(MonteCarlo); library(tidyverse); library(ggplot2)

# estim2 = estimraw::estim_raw(
#   es=rr_raw_data, lb=rr_ci_lo_raw_data, ub=rr_ci_up_raw_data,
#   m1=n_exp, m2=n_nexp, dec = 3, e1 = n_cases, measure = "rr")
#
# logor_dipietrantonj_wo_cases =  suppressWarnings(log((estim2$a[1] * estim2$d[1]) / (estim2$b[1] * estim2$c[1])))
# logor_se_dipietrantonj_wo_cases = suppressWarnings(sqrt(1/estim2$a[1] + 1/estim2$b[1] + 1/estim2$c[1] + 1/estim2$d[1]))
# logor_dipietrantonj_wo_cases, logor_se_dipietrantonj_wo_cases

.estim_metaC <- function(rr_raw_data, rr_ci_lo_raw_data, rr_ci_up_raw_data,
                         n_exp, n_nexp) {
  tryCatch(
    expr = {
      estim = estimraw::estim_raw(
        es=rr_raw_data, lb=rr_ci_lo_raw_data, ub=rr_ci_up_raw_data,
        m1=n_exp, m2=n_nexp, dec = 3, measure = "rr")

      logor =  suppressWarnings(log((estim$a[1] * estim$d[1]) / (estim$b[1] * estim$c[1])))
      logor_se = suppressWarnings(sqrt(1/estim$a[1] + 1/estim$b[1] + 1/estim$c[1] + 1/estim$d[1]))

      dat <- cbind(logor, logor_se)
      return(dat)
    },
    error = function(e) {
      dat <- cbind(NA, NA)
      return(dat)
    },
    warning = function(w) {
      dat <- cbind(NA, NA)
      return(dat)
    },
    finally = {
    }
  )
}
.metaumbrella_rr_se_to_or <- function(rr, logrr, logrr_se, n_cases, n_controls) {

  es = data.frame(value = NA, se = NA)

  if (!is.na(rr) & !is.na(logrr_se) & !is.na(n_cases) & !is.na(n_controls)) {

    n_cases_nexp_sim1 = 0:n_cases
    n_controls_nexp_sim1 = round(n_cases_nexp_sim1 * ((rr * (n_cases + n_controls)) / (n_cases + (rr - 1) * n_cases_nexp_sim1) - 1))
    n_cases_exp_sim1 = n_cases - n_cases_nexp_sim1
    n_controls_exp_sim1 = n_controls - n_controls_nexp_sim1

    idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    n_cases_nexp_sim2 = 0:n_cases
    n_controls_nexp_sim2 = ((n_cases + n_controls - n_cases_nexp_sim2 + 1) - (n_cases + n_controls + 2) * (n_cases - n_cases_nexp_sim2 + 0.5) / ((n_cases_nexp_sim2 + 0.5) * rr + n_cases - n_cases_nexp_sim2 + 0.5))
    n_cases_exp_sim2 = n_cases - n_cases_nexp_sim2
    n_controls_exp_sim2 = n_controls - n_controls_nexp_sim2

    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    n_controls_exp_sim = append(n_controls_exp_sim1, n_controls_exp_sim2)
    n_controls_nexp_sim = append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_nexp_sim = append(n_cases_nexp_sim1, n_cases_nexp_sim2)
    n_cases_exp_sim = append(n_cases_exp_sim1, n_cases_exp_sim2)

    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
    var_sim = ifelse(some_zero,
                     1 / ((n_cases+1) - (n_cases_nexp_sim+0.5)) + 1 / ((n_cases+1) + (n_controls+1) - ((n_cases_nexp_sim+0.5) + (n_controls_nexp_sim+0.5))) +
                       1 / (n_cases_nexp_sim+0.5) + 1 / ((n_cases_nexp_sim+0.5) + (n_controls_nexp_sim+0.5)),
                     1 / (n_cases - n_cases_nexp_sim) + 1 / (n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim)) +
                       1 / n_cases_nexp_sim + 1 / (n_cases_nexp_sim + n_controls_nexp_sim)
    )

    var_sim = 1 / (n_cases - n_cases_nexp_sim) + 1 / (n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim)) +
      1 / n_cases_nexp_sim + 1 / (n_cases_nexp_sim + n_controls_nexp_sim)

    best = order((var_sim - logrr_se^2)^2)[1]
    n_cases_nexp = n_cases_nexp_sim[best]
    n_controls_nexp = n_controls_nexp_sim[best]
    n_cases_exp = n_cases - n_cases_nexp
    n_controls_exp = n_controls - n_controls_nexp

    es$value = (n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp)
    es$se = sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp)
  }
  return(es)
}


# salanti
# n_cases_nexp = sum(rbinom(n_nexp, 1, br))
# n_controls_nexp = n_nexp - n_cases_nexp
#
# t = or * n_cases_nexp / (or * n_cases_nexp + n_controls_nexp)
# n_cases_exp = sum(rbinom(n_exp, 1, t))
# n_controls_exp = n_exp - n_cases_exp
#
# n_cases = n_cases_exp + n_cases_nexp
# n_controls = n_controls_exp + n_controls_nexp
# n_exp = n_cases_exp + n_controls_exp
# n_nexp = n_cases_nexp + n_controls_nexp

# or = 0.25; br = 0.1; br_guess=0.25; n=25; p=0.5
or=2; br=0.5; br_guess=0.5; n=25; p=0.7
#######################################################################################
RR_to_OR <- function(or, br, br_guess, n, p) {

  n_exp = ifelse(runif(1) > p * n - floor(p * n), floor(p * n), ceiling(p * n))
  n_nexp = n - n_exp

  rr = 1 / ((1 - br) / or + br)

  n_cases_exp = sum(rbinom(1, n_exp, rr * br))
  n_controls_exp = n_exp - n_cases_exp

  n_cases_nexp = sum(rbinom(n_nexp, 1, br))
  n_controls_nexp = n_nexp - n_cases_nexp

  n_cases = n_cases_exp + n_cases_nexp
  n_controls = n_controls_exp + n_controls_nexp
  n_exp = n_cases_exp + n_controls_exp
  n_nexp = n_cases_nexp + n_controls_nexp
  br_raw_data = n_cases_nexp / n_nexp

  # Estimate OR/RR from raw data =================================
  zero = which(n_cases_exp == 0 | n_cases_nexp == 0 | n_controls_exp == 0 | n_controls_nexp == 0)

  n_cases_exp[zero] = n_cases_exp[zero] + 0.5
  n_cases_nexp[zero] = n_cases_nexp[zero] + 0.5
  n_controls_exp[zero] = n_controls_exp[zero] + 0.5
  n_controls_nexp[zero] = n_controls_nexp[zero] + 0.5

  br_raw_data = n_cases_nexp / n_nexp
  or_raw_data = suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp))
  or_se_raw_data = suppressWarnings(sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp))

  rr_raw_data = suppressWarnings((n_cases_exp / n_exp) / (n_cases_nexp / n_nexp))
  rr_se_raw_data = suppressWarnings(sqrt(1 / n_cases_exp - 1 / n_exp + 1 / n_cases_nexp - 1 / n_nexp))

  n_cases_exp = n_cases_nexp = n_controls_exp = n_controls_nexp = NULL

  # 3. Convert RR to OR according to grant ===============================
  logor_grant = suppressWarnings(log(rr_raw_data * (1 - br_guess) / (1 - rr_raw_data * br_guess)))
  rr_ci_lo_raw_data = exp(log(rr_raw_data) - qnorm(.975) * rr_se_raw_data)
  rr_ci_up_raw_data = exp(log(rr_raw_data) + qnorm(.975) * rr_se_raw_data)

  logor_ci_lo_grant = suppressWarnings(
    log(rr_ci_lo_raw_data * (1 - br_guess) / (1 - rr_ci_lo_raw_data * br_guess)))

  logor_ci_up_grant = suppressWarnings(
    log(rr_ci_up_raw_data * (1 - br_guess) / (1 - rr_ci_up_raw_data * br_guess)))

  logor_se_grant = (logor_ci_up_grant - logor_ci_lo_grant)/(2 * qnorm(.975))

  # 4. Convert RR to OR - transpose =========================================
  logor_transpose = log(rr_raw_data)
  logor_se_transpose = rr_se_raw_data
  logor_ci_lo_transpose = logor_transpose - qnorm(.975) * logor_se_transpose
  logor_ci_up_transpose = logor_transpose + qnorm(.975) * logor_se_transpose

  # 5. Convert RR to OR according to us (cases) =============================
  raw_res = .metaumbrella_rr_se_to_or(rr = rr_raw_data, logrr_se = rr_se_raw_data,
                                      n_cases = n_cases, n_controls = n_controls)
  logor_metaumbrella = log(raw_res$value)
  logor_se_metaumbrella = raw_res$se
  logor_ci_lo_metaumbrella = logor_metaumbrella - qnorm(.975) * logor_se_metaumbrella
  logor_ci_up_metaumbrella = logor_metaumbrella + qnorm(.975) * logor_se_metaumbrella

  # 6. Convert RR to OR according to (dipie) =============================

  estim = .estim_metaC(rr_raw_data, rr_ci_lo_raw_data, rr_ci_up_raw_data,
                       n_exp, n_nexp)

  logor_dipietrantonj = as.numeric(estim[,1])
  logor_se_dipietrantonj = as.numeric(estim[,2])

  logor_ci_lo_dipietrantonj = logor_dipietrantonj - qnorm(.975) * logor_se_dipietrantonj
  logor_ci_up_dipietrantonj = logor_dipietrantonj + qnorm(.975) * logor_se_dipietrantonj

  return(list(
    "rr_raw_data" = log(rr_raw_data),
    "br_raw_data" = br_raw_data,
    "or_raw_data" = log(or_raw_data),
    "or_se_raw_data" = or_se_raw_data,

    "or_metaumbrella" = logor_metaumbrella,
    "or_se_metaumbrella" = logor_se_metaumbrella,
    "or_ci_lo_metaumbrella" = logor_ci_lo_metaumbrella,
    "or_ci_up_metaumbrella" = logor_ci_up_metaumbrella,

    "or_grant" = logor_grant,
    "or_se_grant" = logor_se_grant,
    "or_ci_lo_grant" = logor_ci_lo_grant,
    "or_ci_up_grant" = logor_ci_up_grant,

    "or_transpose" = logor_transpose,
    "or_se_transpose" = logor_se_transpose,
    "or_ci_lo_transpose" = logor_ci_lo_transpose,
    "or_ci_up_transpose" = logor_ci_up_transpose,

    "or_dipietrantonj" = logor_dipietrantonj,
    "or_se_dipietrantonj" = logor_se_dipietrantonj,
    "or_ci_lo_dipietrantonj" = logor_ci_lo_dipietrantonj,
    "or_ci_up_dipietrantonj" = logor_ci_up_dipietrantonj

  ))
}
or_grid <- c(2, 1, 0.75, 0.5, 0.25)
n_grid <- c(25, 50, 75, 100, 300)
br_grid <- c(0.7, 0.4, 0.15, 0.05, 0.01)
br_guess_grid <- c(0.7, 0.4, 0.15, 0.05, 0.01)
p_grid <- c(0.3, 0.5, 0.7)

param_list = list("or" = or_grid, "n" = n_grid, "br" = br_grid,
                  "br_guess" = br_guess_grid, "p" = p_grid)
MC_result <- MonteCarlo(func = RR_to_OR, nrep = 2500,
                        time_n_test = TRUE,
                        param_list = param_list, max_grid = 10000)
dat_rr <- MakeFrame(MC_result)

rio::export(dat_rr, "D:/simulations/data/RR_to_OR_sim2500d.txt")

res_rr = bind_rows(
  read.delim("D:/simulations/data/RR_to_OR_sim2500a.txt"),
  read.delim("D:/simulations/data/RR_to_OR_sim2500b.txt"),
  read.delim("D:/simulations/data/RR_to_OR_sim2500c.txt"),
  read.delim("D:/simulations/data/RR_to_OR_sim2500d.txt")) %>%
  mutate_if(is.character, as.numeric)


res_rr[res_rr==Inf|res_rr==-Inf] <- NA

method=c("metaumbrella", "grant",
         "dipietrantonj", "transpose")
for (met in method) {
  row = NA
  row = which(is.na(res_rr[, paste0("or_", met)]))
  res_rr[row, paste0("or_se_", met)] <- NA
  res_rr[row, paste0("or_ci_lo_", met)] <- NA
  res_rr[row, paste0("or_ci_up_", met)] <- NA
}


mean_na = function(x) {
  res = mean(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}
var_na = function(x) {
  res = var(x, na.rm = TRUE)
  if (sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    return(res)
  }
}
breks = res_rr %>%
  group_by(or) %>%
  summarise(or_raw = mean_na(or_raw_data),
            br_raw_data = mean_na(br_raw_data),
            n=n())
res_rr$logor = round(log(res_rr$or),2)
ggplot(res_rr, aes(x = or_raw_data,
                   group=as.character(logor),
                   fill=as.character(logor))) +
  geom_density(alpha=0.5) +
  geom_vline(data=breks, aes(xintercept = or_raw)) +
  scale_x_continuous(name="or",
                     breaks=breks$or_raw,labels=paste0(round(breks$or_raw,2)))

checks = res_rr %>%
  group_by(or, br) %>%
  summarise(n =n(),
            or_raw = exp(mean_na(or_raw_data)),
            br_raw = mean_na(br_raw_data))

for (met in method) {
  res_rr[, paste0("bias_es_or_", met)] <- res_rr[, paste0("or_", met)] - res_rr$or_raw_data
  res_rr[, paste0("or_var_", met)] = res_rr[, paste0("or_se_", met)]^2

  res_rr[, paste0("ci_cov_or_", met)] = res_rr[, paste0("or_ci_lo_", met)] <= res_rr$or_raw_data &
    res_rr[, paste0("or_ci_up_", met)] >= res_rr$or_raw_data
}


res = res_rr %>%
  group_by(or, n, br, br_guess) %>%
  summarise(
    n_sim = n(),
    bias_es_or_metaumbrella = mean_na(bias_es_or_metaumbrella),
    var_or_metaumbrella = var_na(or_metaumbrella),
    mean_or_var_metaumbrella = mean_na(or_var_metaumbrella),
    bias_ci_or_metaumbrella = mean_na(ci_cov_or_metaumbrella),

    bias_es_or_dipietrantonj = mean_na(bias_es_or_dipietrantonj),
    var_or_dipietrantonj = var_na(or_dipietrantonj),
    mean_or_var_dipietrantonj = mean_na(or_var_dipietrantonj),
    bias_ci_or_dipietrantonj = mean_na(ci_cov_or_dipietrantonj),

    bias_es_or_transpose = mean_na(bias_es_or_transpose),
    var_or_transpose = var_na(or_transpose),
    mean_or_var_transpose = mean_na(or_var_transpose),
    bias_ci_or_transpose = mean_na(ci_cov_or_transpose),

    bias_es_or_grant = mean_na(bias_es_or_grant),
    var_or_grant = var_na(or_grant),
    mean_or_var_grant = mean_na(or_var_grant),
    bias_ci_or_grant = mean_na(ci_cov_or_grant)
  )

for (met in method) {
  print(met)
  res[, paste0("bias_var_or_", met)] <- res[, paste0("mean_or_var_", met)] / res[, paste0("var_or_", met)]
}

err_plot = res %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("bias_es") |
                 starts_with("bias_var") |
                 starts_with("bias_ci"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(or, n, br, br_guess, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("bias_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("bias_var", method_long, fixed = TRUE) ~ "_var",
                           grepl("bias_ci", method_long, fixed = TRUE) ~ "_ci"),
         method = case_when(grepl("metaumbrella", method_long, fixed = TRUE) ~ "metaumbrella",
                            grepl("grant", method_long, fixed = TRUE) ~ "grant",
                            grepl("dipietrantonj", method_long, fixed = TRUE) ~ "dipietrantonj",
                            grepl("transpose", method_long, fixed = TRUE) ~ "transpose")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot$or <- paste0("or=", err_plot$or)
err_plot$br <- paste0("br=", err_plot$br)
err_plot$br_guess <- paste0("guessed br=", err_plot$br_guess)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "300"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=300"))

err_plot$method[err_plot$method == "dipietrantonj"] <- "Di Pietrantonj"
err_plot$method[err_plot$method == "transpose"] <- "or = rr"
err_plot$method[err_plot$method == "metaumbrella"] <- "metaumbrella"
err_plot$method[err_plot$method == "grant"] <- "Grant"

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/metaConvert/simulations/data_agg/RR_to_OR_AGG10000.txt")

# 1. Convert RR to OR according to Grant_2x2 ==============================
# logor_grant_2x2 = suppressWarnings(log(rr_raw_data * (1 - br_guess) / (1 - rr_raw_data * br_guess)))
# a_est = rr_raw_data * br_guess * n_exp
# b_est = (1 - rr_raw_data * br_guess) * n_exp
# c_est = br_guess * n_nexp
# d_est = (1 - br_guess) * n_nexp
# if (a_est==0|b_est==0|c_est==0|d_est==0) {
#   logor_se_grant_2x2 = logor_ci_lo_grant_2x2 = logor_ci_up_grant_2x2 = NA
# } else {
#   logor_se_grant_2x2 = sqrt(1/a_est + 1/b_est + 1/c_est + 1/d_est)
#   logor_ci_lo_grant_2x2 = logor_grant_2x2 - qnorm(.975) * logor_se_grant_2x2
#   logor_ci_up_grant_2x2 = logor_grant_2x2 + qnorm(.975) * logor_se_grant_2x2
# }
#
# # 2. Convert RR to OR according to Grant_delta ============================
# logor_grant_delta = logor_grant_2x2
# logor_se_grant_delta = sqrt(do.call(numDeriv::grad, list(
#   x = log(rr_raw_data),
#   baseline_risk=br_guess,
#   func = function(rr, baseline_risk) {
#     exp(rr) * (1 - baseline_risk) / (1 -  exp(rr) * baseline_risk)
#   }))^2 * rr_se_raw_data^2)
# logor_ci_lo_grant_delta = logor_grant_delta - qnorm(.975) * logor_se_grant_delta
# logor_ci_up_grant_delta = logor_grant_delta + qnorm(.975) * logor_se_grant_delta
