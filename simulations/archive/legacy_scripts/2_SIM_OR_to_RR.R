library(MonteCarlo); library(tidyverse); library(ggplot2)
.estim_metaRR <- function(or_raw_data, or_ci_lo_raw_data, or_ci_up_raw_data,
                         n_exp, n_nexp, br_raw_data) {
  tryCatch(
    expr = {
      estim = estimraw::estim_raw(
        es=or_raw_data, lb=or_ci_lo_raw_data, ub=or_ci_up_raw_data,
        m1=n_exp, m2=n_nexp, dec = 3, measure = "or")

      if (length(estim) != 4) {
        numb = which.min(
          c(
            abs(estim[[1]]$c[1] / (estim[[1]]$c[1] + estim[[1]]$d[1]) - br_raw_data),
            abs(estim[[2]]$c[1] / (estim[[2]]$c[1] + estim[[2]]$d[1]) - br_raw_data)
          )
        )

        rr_dipietrantonj_rand = suppressWarnings((estim[[1]]$a[1] / (estim[[1]]$a[1]+estim[[1]]$b[1])) /
                                        (estim[[1]]$c[1] / (estim[[1]]$c[1]+estim[[1]]$d[1])))
        logrr_se_dipietrantonj_rand = suppressWarnings(sqrt(1/estim[[1]]$a[1] - 1 / (estim[[1]]$a[1]+estim[[1]]$b[1]) +
                                              1/estim[[1]]$c[1] - 1 / (estim[[1]]$c[1]+estim[[1]]$d[1])))

        rr_dipietrantonj_br = suppressWarnings((estim[[numb]]$a[1] / (estim[[numb]]$a[1]+estim[[numb]]$b[1])) /
                                        (estim[[numb]]$c[1] / (estim[[1]]$c[1]+estim[[numb]]$d[1])))
        logrr_se_dipietrantonj_br = suppressWarnings(sqrt(1/estim[[numb]]$a[1] - 1 / (estim[[numb]]$a[numb]+estim[[numb]]$b[1]) +
                                            1/estim[[numb]]$c[1] - 1 / (estim[[numb]]$c[numb]+estim[[numb]]$d[1])))

      } else {
        rr_dipietrantonj_rand = rr_dipietrantonj_br = suppressWarnings((estim[[1]]$a[1] / (estim[[1]]$a[1]+estim[[1]]$b[1])) /
                                                        (estim[[1]]$c[1] / (estim[[1]]$c[1]+estim[[1]]$d[1])))
        logrr_se_dipietrantonj_rand = logrr_se_dipietrantonj_br = suppressWarnings(sqrt(1/estim[[1]]$a[1] - 1 / (estim[[1]]$a[1]+estim[[1]]$b[1]) +
                                                1/estim[[1]]$c[1] - 1 / (estim[[1]]$c[1]+estim[[1]]$d[1])))

      }
      dat <- cbind(logrr_br = log(rr_dipietrantonj_br),
                   logrr_se_br = logrr_se_dipietrantonj_br,
                   logrr_rand = log(rr_dipietrantonj_rand),
                   logrr_se_rand = logrr_se_dipietrantonj_rand)
      return(dat)
    },
    error = function(e) {
      dat <- cbind(NA, NA, NA, NA)
      return(dat)
    },
    warning = function(w) {
      dat <- cbind(NA, NA, NA, NA)
      return(dat)
    },
    finally = {
    }
  )
}

.estimate_n_from_or_and_n_cases = function (or, var, n_cases, n_controls) {

  res = data.frame(n_cases_exp = NA, n_cases_nexp = NA, n_controls_exp = NA, n_controls_nexp = NA)

  if (!is.na(or) & !is.na(var) & !is.na(n_cases) & !is.na(n_controls)) {

    n_cases_nexp_sim1 = 0:n_cases
    n_controls_nexp_sim1 = round(n_controls * (1 - (n_cases - n_cases_nexp_sim1) /
                                                 (n_cases + (or - 1) * n_cases_nexp_sim1)))
    n_cases_exp_sim1 = n_cases - n_cases_nexp_sim1
    n_controls_exp_sim1 = n_controls - n_controls_nexp_sim1
    idx_non_zero <- which(
      n_cases_nexp_sim1 > 0 &
        n_controls_nexp_sim1 > 0 &
        n_cases_exp_sim1 > 0 &
        n_controls_exp_sim1 > 0
    )
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    n_cases_nexp_sim2 = 0:n_cases
    n_controls_nexp_sim2 = round((n_controls + 0.5) - (n_controls + 1) * (n_cases - n_cases_nexp_sim2 + 0.5) /
                                   ((n_cases_nexp_sim2 + 0.5) * or + n_cases - n_cases_nexp_sim2 + 0.5))
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

    n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)
    n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)
    n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)

    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0

    var_sim <- ifelse(some_zero,
                      1 / ((n_cases+ 1 ) - (n_cases_nexp_sim + 0.5)) + 1 / ((n_controls + 1) - (n_controls_nexp_sim + 0.5)) + 1 / (n_cases_nexp_sim + 0.5) + 1 / (n_controls_nexp_sim + 0.5),
                      1 / (n_cases - n_cases_nexp_sim) + 1 / (n_controls - n_controls_nexp_sim) + 1 / n_cases_nexp_sim + 1 / n_controls_nexp_sim)

    best = order((var_sim - var)^2)[1]

    res$n_cases_nexp = n_cases_nexp_sim[best]
    res$n_controls_nexp = n_controls_nexp_sim[best]
    res$n_cases_exp = n_cases - res$n_cases_nexp
    res$n_controls_exp = n_controls - res$n_controls_nexp

  }

  return(res)
}

#' Estimate the n, using the variance, the number of exposed and non-exposed subjects
#'
#' @param or OR
#' @param var variance
#' @param n_exp number of exposed participants
#' @param n_nexp number of non exposed participants
#'
#' @noRd
.estimate_n_from_or_and_n_exp = function (or, var, n_exp, n_nexp) {

  res = data.frame(n_cases_exp = NA, n_cases_nexp = NA, n_controls_exp = NA, n_controls_nexp = NA)

  if (!is.na(or) & !is.na(var) & !is.na(n_exp) & !is.na(n_nexp)) {
    n_controls_exp_sim1 = 0:n_exp
    n_controls_nexp_sim1 = round(n_nexp / (1 + (n_exp - n_controls_exp_sim1) / (or * n_controls_exp_sim1)))
    n_cases_exp_sim1 = n_exp - n_controls_exp_sim1
    n_cases_nexp_sim1 = n_nexp - n_controls_nexp_sim1
    idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    n_controls_exp_sim2 = 0:n_exp
    n_controls_nexp_sim2 = round((n_nexp + 0.5) - ((n_nexp + 1)*(n_exp - n_controls_exp_sim2 + 0.5)) / ((n_controls_exp_sim2 + 0.5) * or + n_exp - n_controls_exp_sim2 + 0.5 ))
    n_cases_exp_sim2 = n_exp - n_controls_exp_sim2
    n_cases_nexp_sim2 = n_nexp - n_controls_nexp_sim2

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
    n_cases_exp_sim = append(n_cases_exp_sim1, n_cases_exp_sim2)
    n_cases_nexp_sim = append(n_cases_nexp_sim1, n_cases_nexp_sim2)


    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
    var_sim <- ifelse(some_zero,
                      1 / ((n_exp+1) - (n_controls_exp_sim+0.5) + 1/(n_controls_exp_sim+0.5) + 1/((n_nexp+1) - (n_controls_nexp_sim+0.5)) + 1/(n_controls_nexp_sim+0.5)),
                      1 / (n_exp - n_controls_exp_sim) + 1 / n_controls_exp_sim + 1 / (n_nexp - n_controls_nexp_sim) + 1 / n_controls_nexp_sim
    )

    best = order((var_sim - var)^2)[1]
    res$n_controls_exp = n_controls_exp_sim[best]
    res$n_controls_nexp = n_controls_nexp_sim[best]
    res$n_cases_exp = n_exp - res$n_controls_exp
    res$n_cases_nexp = n_nexp - res$n_controls_nexp
  }
  return(res)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
rr=2; br=0.5; n=25; p=0.3; br_guess=0.5

OR_to_RR <- function(rr, br, n, p, br_guess) {

  n_exp = ifelse(runif(1) > p * n - floor(p * n), floor(p * n), ceiling(p * n))
  n_nexp = n - n_exp

  n_cases_exp = sum(rbinom(1, n_exp, rr * br))
  n_controls_exp = n_exp - n_cases_exp

  n_cases_nexp = sum(rbinom(n_nexp, 1, br))
  n_controls_nexp = n_nexp - n_cases_nexp

  n_cases = n_cases_exp + n_cases_nexp
  n_controls = n_controls_exp + n_controls_nexp
  n_exp = n_cases_exp + n_controls_exp
  n_nexp = n_cases_nexp + n_controls_nexp

  zero = which(n_cases_exp == 0 | n_cases_nexp == 0 | n_controls_exp == 0 | n_controls_nexp == 0)

  n_cases_exp[zero] = n_cases_exp[zero] + 0.5
  n_cases_nexp[zero] = n_cases_nexp[zero] + 0.5
  n_controls_exp[zero] = n_controls_exp[zero] + 0.5
  n_controls_nexp[zero] = n_controls_nexp[zero] + 0.5

  br_raw_data = n_cases_nexp / n_nexp
  or_raw_data = suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp))
  or_se_raw_data = suppressWarnings(sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp))
  or_ci_lo_raw_data = exp(log(or_raw_data) - qnorm(.975) * or_se_raw_data)
  or_ci_up_raw_data = exp(log(or_raw_data) + qnorm(.975) * or_se_raw_data)

  rr_raw_data = suppressWarnings((n_cases_exp / n_exp) / (n_cases_nexp / n_nexp))

  logrr_grant = suppressWarnings(log(or_raw_data / (1 - br_guess + (br_guess * or_raw_data))))

  logrr_ci_lo_grant = suppressWarnings(log(or_ci_lo_raw_data / (1 - br_guess + (br_guess * or_ci_lo_raw_data))))
  logrr_ci_up_grant = suppressWarnings(log(or_ci_up_raw_data / (1 - br_guess + (br_guess * or_ci_up_raw_data))))
  logrr_se_grant = (logrr_ci_up_grant - logrr_ci_lo_grant)/(2 * qnorm(.975))
  logrr_ci_lo_grant = logrr_grant - qnorm(.975) * logrr_se_grant
  logrr_ci_up_grant = logrr_grant + qnorm(.975) * logrr_se_grant

  contingency_meta_cases = .estimate_n_from_or_and_n_cases(
    or = or_raw_data, var = or_se_raw_data^2,
    n_cases = n_cases, n_controls = n_controls)

    logrr_meta_cases = log(with(contingency_meta_cases, suppressWarnings(
      (n_cases_exp / (n_cases_exp + n_controls_exp)) /
      (n_cases_nexp / (n_cases_nexp + n_controls_nexp)))))
    logrr_se_meta_cases = with(contingency_meta_cases, suppressWarnings(
      sqrt(1 / n_cases_exp - 1 / (n_cases_exp + n_controls_exp) +
           1 / n_cases_nexp - 1 / (n_cases_nexp + n_controls_nexp))))
    logrr_ci_lo_meta_cases = logrr_meta_cases - qnorm(.975) * logrr_se_meta_cases
    logrr_ci_up_meta_cases = logrr_meta_cases + qnorm(.975) * logrr_se_meta_cases


    contingency_meta_exp = .estimate_n_from_or_and_n_exp(
      or = or_raw_data, var = or_se_raw_data^2, n_exp = n_exp, n_nexp = n_nexp)

    logrr_meta_exp = log(with(contingency_meta_exp, suppressWarnings(
      (n_cases_exp / (n_cases_exp + n_controls_exp)) /
      (n_cases_nexp / (n_cases_nexp + n_controls_nexp)))))
    logrr_se_meta_exp = with(contingency_meta_exp, suppressWarnings(
      sqrt(1 / n_cases_exp - 1 / (n_cases_exp + n_controls_exp) +
           1 / n_cases_nexp - 1 / (n_cases_nexp + n_controls_nexp))))
    logrr_ci_lo_meta_exp = logrr_meta_exp - qnorm(.975) * logrr_se_meta_exp
    logrr_ci_up_meta_exp = logrr_meta_exp + qnorm(.975) * logrr_se_meta_exp

    logrr_transpose = log(or_raw_data)
    logrr_se_transpose = or_se_raw_data
    logrr_ci_lo_transpose = logrr_transpose - qnorm(.975) * logrr_se_transpose
    logrr_ci_up_transpose = logrr_transpose + qnorm(.975) * logrr_se_transpose

    estim = .estim_metaRR(or_raw_data, or_ci_lo_raw_data, or_ci_up_raw_data,
                         n_exp, n_nexp, br_raw_data)

    logrr_dipietrantonj_br = as.numeric(estim[,1])
    logrr_se_dipietrantonj_br = as.numeric(estim[,2])
    logrr_dipietrantonj_rand = as.numeric(estim[,3])
    logrr_se_dipietrantonj_rand = as.numeric(estim[,4])

    logrr_ci_lo_dipietrantonj_rand = logrr_dipietrantonj_rand - qnorm(.975) * logrr_se_dipietrantonj_rand
    logrr_ci_up_dipietrantonj_rand = logrr_dipietrantonj_rand + qnorm(.975) * logrr_se_dipietrantonj_rand

    logrr_ci_lo_dipietrantonj_br = logrr_dipietrantonj_br - qnorm(.975) * logrr_se_dipietrantonj_br
    logrr_ci_up_dipietrantonj_br = logrr_dipietrantonj_br + qnorm(.975) * logrr_se_dipietrantonj_br

  return(list(
    "or_raw_data" = log(or_raw_data),
    "br_raw_data" = br_raw_data,
    "rr_raw_data" = log(rr_raw_data),

    "rr_metaumbrella_cases" = logrr_meta_cases,
    "rr_se_metaumbrella_cases" = logrr_se_meta_cases,
    "rr_ci_lo_metaumbrella_cases" = logrr_ci_lo_meta_cases,
    "rr_ci_up_metaumbrella_cases" = logrr_ci_up_meta_cases,

    "rr_metaumbrella_exp" = logrr_meta_exp,
    "rr_se_metaumbrella_exp" = logrr_se_meta_exp,
    "rr_ci_lo_metaumbrella_exp" = logrr_ci_lo_meta_exp,
    "rr_ci_up_metaumbrella_exp" = logrr_ci_up_meta_exp,

    "rr_grant" = logrr_grant,
    "rr_se_grant" = logrr_se_grant,
    "rr_ci_lo_grant" = logrr_ci_lo_grant,
    "rr_ci_up_grant" = logrr_ci_up_grant,

    "rr_transpose" = logrr_transpose,
    "rr_se_transpose" = logrr_se_transpose,
    "rr_ci_lo_transpose" = logrr_ci_lo_transpose,
    "rr_ci_up_transpose" = logrr_ci_up_transpose,

    "rr_dipietrantonj_br" = logrr_dipietrantonj_br,
    "rr_se_dipietrantonj_br" = logrr_se_dipietrantonj_br,
    "rr_ci_lo_dipietrantonj_br" = logrr_ci_lo_dipietrantonj_br,
    "rr_ci_up_dipietrantonj_br" = logrr_ci_up_dipietrantonj_br,

    "rr_dipietrantonj_rand" = logrr_dipietrantonj_rand,
    "rr_se_dipietrantonj_rand" = logrr_se_dipietrantonj_rand,
    "rr_ci_lo_dipietrantonj_rand" = logrr_ci_lo_dipietrantonj_rand,
    "rr_ci_up_dipietrantonj_rand" = logrr_ci_up_dipietrantonj_rand

  ))
}

rr_grid <- c(2, 1, 0.75, 0.5, 0.25)
n_grid <- c(25, 50, 75, 100, 300)
br_grid <- c(0.5, 0.15, 0.05, 0.01)
br_guess_grid <- c(0.5, 0.15, 0.05, 0.01)
p_grid <- c(0.3, 0.5, 0.7)

param_list = list("rr" = rr_grid, "n" = n_grid, "br" = br_grid,
                  "p" = p_grid, "br_guess" = br_guess_grid)
start_time <- Sys.time()
MC_result <- MonteCarlo(func = OR_to_RR, nrep = 7500,
                        time_n_test = TRUE,
                        max_grid = 10000, param_list = param_list)
end_time <- Sys.time()
end_time - start_time
data_sim <- MakeFrame(MC_result)

rio::export(data_sim, "D:/simulations/data/OR_to_RR_sim7500a.txt")

res_or = bind_rows(
  read.delim("D:/simulations/data/OR_to_RR_sim2500a.txt"),
  read.delim("D:/simulations/data/OR_to_RR_sim7500a.txt")) %>%
  mutate_if(is.character, as.numeric)

res_or[res_or==Inf|res_or==-Inf] <- NA

method=c("metaumbrella_cases", "metaumbrella_exp", "grant",
         "dipietrantonj_rand",
         "dipietrantonj_br", "transpose")

for (met in method) {
  row = NA
  row = which(is.na(res_or[, paste0("rr_", met)]))
  res_or[row, paste0("rr_se_", met)] <- NA
  res_or[row, paste0("rr_ci_lo_", met)] <- NA
  res_or[row, paste0("rr_ci_up_", met)] <- NA
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
breks = res_or %>%
  group_by(rr) %>%
  summarise(rr_raw = mean_na(rr_raw_data),
            br_raw_data = mean_na(br_raw_data),
            n=n())
res_or$logrr = round(log(res_or$rr),2)
ggplot(res_or, aes(x = rr_raw_data,
                   group=as.character(logrr),
                   fill=as.character(logrr))) +
  geom_density(alpha=0.5) +
  geom_vline(data=breks, aes(xintercept = rr_raw),
             color="red") +
  scale_x_continuous(name="rr",
                   breaks=breks$rr_raw,labels=paste0(round(breks$rr_raw,2)))

ggplot(res_or, aes(x = br_raw_data,
                   group=as.character(br),
                   fill=as.character(br))) +
  geom_density(alpha=0.5) +
  geom_vline(data=breks, aes(xintercept = br_raw_data),
             color="red") +
  scale_x_continuous(name="br")

checks = res_or %>%
  group_by(rr, br) %>%
  summarise(n =n(),
            rr_raw = exp(mean_na(rr_raw_data)),
            br_raw = mean_na(br_raw_data))
View(checks)
for (met in method) {

  res_or[, paste0("bias_es_rr_", met)] <- res_or[, paste0("rr_", met)] - res_or$rr_raw_data
  res_or[, paste0("rr_var_", met)] = res_or[, paste0("rr_se_", met)]^2

  res_or[, paste0("ci_cov_rr_", met)] = res_or[, paste0("rr_ci_lo_", met)] <= res_or$rr_raw_data &
                                        res_or[, paste0("rr_ci_up_", met)] >= res_or$rr_raw_data
}


res = res_or %>%
  group_by(rr, n, br, br_guess) %>%
  summarise(
    n_sim = n(),
    bias_es_rr_metaumbrella_cases = mean_na(bias_es_rr_metaumbrella_cases),
    var_rr_metaumbrella_cases = var_na(rr_metaumbrella_cases),
    mean_rr_var_metaumbrella_cases = mean_na(rr_var_metaumbrella_cases),
    bias_ci_rr_metaumbrella_cases = mean_na(ci_cov_rr_metaumbrella_cases),

    bias_es_rr_metaumbrella_exp = mean_na(bias_es_rr_metaumbrella_exp),
    var_rr_metaumbrella_exp = var_na(rr_metaumbrella_exp),
    mean_rr_var_metaumbrella_exp = mean_na(rr_var_metaumbrella_exp),
    bias_ci_rr_metaumbrella_exp = mean_na(ci_cov_rr_metaumbrella_exp),

    bias_es_rr_dipietrantonj_rand = mean_na(bias_es_rr_dipietrantonj_rand),
    var_rr_dipietrantonj_rand = var_na(rr_dipietrantonj_rand),
    mean_rr_var_dipietrantonj_rand = mean_na(rr_var_dipietrantonj_rand),
    bias_ci_rr_dipietrantonj_rand = mean_na(ci_cov_rr_dipietrantonj_rand),

    bias_es_rr_dipietrantonj_br = mean_na(bias_es_rr_dipietrantonj_br),
    var_rr_dipietrantonj_br = var_na(rr_dipietrantonj_br),
    mean_rr_var_dipietrantonj_br = mean_na(rr_var_dipietrantonj_br),
    bias_ci_rr_dipietrantonj_br = mean_na(ci_cov_rr_dipietrantonj_br),

    bias_es_rr_transpose = mean_na(bias_es_rr_transpose),
    var_rr_transpose = var_na(rr_transpose),
    mean_rr_var_transpose = mean_na(rr_var_transpose),
    bias_ci_rr_transpose = mean_na(ci_cov_rr_transpose),

    bias_es_rr_grant = mean_na(bias_es_rr_grant),
    var_rr_grant = var_na(rr_grant),
    mean_rr_var_grant = mean_na(rr_var_grant),
    bias_ci_rr_grant = mean_na(ci_cov_rr_grant)
  )

for (met in method) {
  print(met)
  res[, paste0("bias_var_rr_", met)] <- res[, paste0("mean_rr_var_", met)] / res[, paste0("var_rr_", met)]
}

err_plot = res %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(cols = starts_with("bias_es") |
                 starts_with("bias_var") |
                 starts_with("bias_ci"),
               names_to="method_long",
               values_to="bias") %>%
  group_by(rr, n, br, br_guess, method_long) %>%
  summarise(bias = mean_na(bias)) %>%
  mutate(pivot = case_when(grepl("bias_es", method_long, fixed = TRUE) ~ "_es",
                           grepl("bias_var", method_long, fixed = TRUE) ~ "_var",
                           grepl("bias_ci", method_long, fixed = TRUE) ~ "_ci"),
         method = case_when(grepl("metaumbrella_cases", method_long, fixed = TRUE) ~ "metaumbrella_cases",
                            grepl("metaumbrella_exp", method_long, fixed = TRUE) ~ "metaumbrella_exp",
                            grepl("grant", method_long, fixed = TRUE) ~ "grant",
                            grepl("grant_2x2", method_long, fixed = TRUE) ~ "grant_2x2",
                            grepl("grant_delta", method_long, fixed = TRUE) ~ "grant_delta",
                            grepl("dipietrantonj_br", method_long, fixed = TRUE) ~ "dipietrantonj_br",
                            grepl("dipietrantonj_rand", method_long, fixed = TRUE) ~ "dipietrantonj_rand",
                            grepl("transpose", method_long, fixed = TRUE) ~ "transpose")) %>%
  select(-c(method_long)) %>%
  pivot_wider(names_from="pivot", values_from="bias", names_prefix="bias")

err_plot$rr <- paste0("RR=", err_plot$rr)
err_plot$br <- paste0("br=", err_plot$br)
err_plot$br_guess <- paste0("guessed br=", err_plot$br_guess)
err_plot$n <- factor(err_plot$n, levels = c("25", "50", "75", "100", "300"),
                     labels = c("n=25", "n=50", "n=75", "n=100", "n=300"))

err_plot$method[err_plot$method == "dipietrantonj_br"] <- "dipietrantonj (br)"
err_plot$method[err_plot$method == "dipietrantonj_rand"] <- "dipietrantonj (random)"
err_plot$method[err_plot$method == "transpose"] <- "OR = RR"
err_plot$method[err_plot$method == "metaumbrella_cases"] <- "metaumbrella (cases)"
err_plot$method[err_plot$method == "metaumbrella_exp"] <- "metaumbrella (exp)"
err_plot$method[err_plot$method == "grant"] <- "Grant"

rio::export(err_plot, "C:/Users/Corentin Gosling/drive_gmail/Recherche/metaConvert/simulations/data_agg/OR_to_RR_AGG10000.txt")
