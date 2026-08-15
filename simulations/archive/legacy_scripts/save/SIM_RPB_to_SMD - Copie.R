library(MonteCarlo); library(tidyverse); library(ggplot2)

rpb_to_d <- function(d, n, p) {
  # d = 0.8
  # n = 40
  # p = 0.5
  # generate sample scores
  n1 = round(p * n)
  n2 = n - n1
  scores_grp1 <- rnorm(n1, d, 1)
  scores_grp2 <- rnorm(n2, 0, 1)
  m1 <- mean(scores_grp1)
  m2 <- mean(scores_grp2)
  sd1 <- sd(scores_grp1)
  sd2 <- sd(scores_grp2)
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

  d_raw_data <- (m1 - m2) / pooled_sd
  d_se_raw_data <- sqrt(1/n1 + 1/n2 + d^2/(2*(n1+n2)))
  pval_raw_data = 1 - 2 * abs(pt(d_raw_data / d_se_raw_data, n1+n2-2) - 0.5)

  cor = cor.test(~c(scores_grp1, scores_grp2) + rep(c(0, 1), c(n1, n2)))

  r_pb_raw_data_escalc <- metafor::escalc(vtype="ST", measure = "RPB", m1i = m1, m2i = m2, sd1i = sd1, sd2i = sd2, n1i = n1, n2i = n2)
  r_pb_raw_data <- as.numeric(as.character(r_pb_raw_data_escalc$yi))
  se.r_pb_raw_data <- sqrt(as.numeric(as.character(r_pb_raw_data_escalc$vi)))
  # pval.r_pb_raw = 1 - 2 * abs(pt(r_pb_raw_data * sqrt((n1+n2-2)/(1-r_pb_raw_data^2)), n1+n2-2) - 0.5)
  # pval.r_pb_raw = 1 - 2 * abs(pt(r_pb_raw_data/sqrt((n1+n2-2)/(1-r_pb_raw_data^2)), n1+n2-2) - 0.5)
  pval.r_pb_raw = 2 * pnorm(-r_pb_raw_data / se.r_pb_raw_data)

  # Convert rpb to d according to Lispey
  d.lipsey = r_pb_raw_data / (sqrt((1-r_pb_raw_data^2) * (n1/(n1+n2)) * (1-(n1/(n1+n2)))))
  v.d = 1/n1 + 1/n2 + d^2/(2*(n1+n2))
  # pval.d_lipsey = 2 * pnorm(-d.lipsey / sqrt(v.d))
  pval.d_lipsey = 1 - 2 * abs(pt(d.lipsey / sqrt(v.d), n1+n2-2) - 0.5)

  # Convert rpb to d according to Wolfgang
  df = n1 + n2 - 2
  h = df/n1 + df/n2
  d.viecht = r_pb_raw_data * sqrt(h) / sqrt(1 - r_pb_raw_data^2)
  pval.d_viecht = 1 - 2 * abs(pt(d.viecht / sqrt(v.d), n1+n2-2) - 0.5)

  # return result:
  return(list(
    "d_raw_data" = d_raw_data,
    "d_se_raw_data" = d_se_raw_data,
    "pval_raw_data" = pval_raw_data,
    "pval_lipsey" = pval.d_lipsey,
    "pval_viecht" = pval.d_viecht,
    "d_viecht" = d.viecht,
    "d_lipsey" = d.lipsey,
    "d_se_viecht" = sqrt(v.d),
    "d_se_lipsey" = sqrt(v.d)))
}
d_grid <- c(0, 0.2, 0.5, 0.8)
n_grid <- c(25, 50, 75, 100, 200)
p_grid <- c(0.3, 0.5, 0.7)

# collect parameter grids in list:
param_list = list("d" = d_grid, "n" = n_grid, "p" = p_grid)
MC_result <- MonteCarlo(func = rpb_to_d, nrep = 500, param_list = param_list)
res <- MakeFrame(MC_result)

setwd("C:/Users/Corentin Gosling/Documents/metaConvert/simulations")
rio::export(res, "RPB_to_SMD.txt")
