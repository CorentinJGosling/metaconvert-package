#' Convert a phi value to several effect size measures
#'
#' @param phi phi value
#' @param n_sample total number of participants in the sample
#' @param n_cases total number of cases/events
#' @param n_exp total number of participants in the exposed group
#' @param reverse_phi a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' The functions computes an odds ratio (OR), risk ratio (RR), and number needed to treat (NNT)
#' from the the phi coefficient, the total number of participants,
#' the total number of cases and the total number of people exposed.
#' Cohen's d (D) and Hedges' g (G) are tried to be obtained from the OR, or are converted
#' using the approach by Lipsey et al. (2001).
#' The correlation coefficients (R/Z) are converted by assuming that the phi coefficient
#' is equal to a R, and the variances of R and Z are obtained using the approach proposed
#' by Lipsey et al. (2001) as well as by our own calculations.
#'
#' **To estimate the OR, RR, NNT,**,
#' this function reconstructs a 2x2 table (using the approach proposed by Viechtbauer, 2023).
#'
#' Then, the calculations of the \code{\link{es_from_2x2}()} function are applied.
#'
#' **To estimate the Cohen's d (D) and Hedges' g (G)**, the function first tries to convert it from
#' the OR obtained using the approach described above. If not possible (e.g., the number of cases and exposed are missing)
#' the function converts the Cohen's d from the Phi coefficient using the approach proposed by Lipsey et al. (2001):
#' \deqn{d = \frac{2 * phi}{\sqrt{1 - phi^2}}}
#' \deqn{d\_se = \sqrt{\frac{d}{phi^2 * n\_sample}}}
#'
#' **To estimate the correlation coefficients (R/Z)**, this function assumes that the
#' phi coefficient is equal to a correlation coefficient, and then obtains the variance using the
#' formula proposed by Lipsey et al. (2001):
#' \deqn{r = phi}
#' \deqn{z = atanh(r)}
#' \deqn{z\_se = \frac{z^2}{phi^2 * n\_sample}}
#' \deqn{effective\_n = \frac{1}{z\_se + 3}}
#' \deqn{r\_se = \sqrt{\frac{(1 - r^2)^2}{effective\_n - 1}}}
#'
#' Note that the approach to determine the standard error of R was developed by our team.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 8. Phi or chi-square'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Viechtbauer (2023). Accessed at https://wviechtb.github.io/metafor/reference/conv.2x2.html.
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' @export es_from_phi
#'
#' @md
#'
#' @examples
#' es_from_phi(phi = 0.3, n_sample = 120, n_cases = 20, n_exp = 40)
es_from_phi <- function(phi, n_cases, n_exp,
                        n_sample,
                        reverse_phi) {
  if (missing(reverse_phi)) reverse_phi <- rep(FALSE, length(phi))
  reverse_phi[is.na(reverse_phi)] <- FALSE

  if (missing(n_sample)) n_sample <- rep(NA, length(phi))
  if (missing(n_exp)) n_exp <- rep(NA, length(phi))
  if (missing(n_cases)) n_cases <- rep(NA, length(phi))

  if (length(reverse_phi) == 1) reverse_phi = c(rep(reverse_phi, length(phi)))
  if (length(reverse_phi) != length(phi)) stop("The length of the 'reverse_phi' argument of incorrectly specified.")

  tryCatch({
    .validate_positive(n_cases, n_exp, n_sample,
                       error_message = paste0("The number of cases, people exposed, total sample ",
                                              "should be >0."),
                       func = "es_from_phi")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  cont_table <- suppressWarnings(
    metafor::conv.2x2(
      ri = phi, ni = n_sample,
      n1i = n_exp, n2i = n_cases
    )
  )
  es <- es_from_2x2(
    n_cases_exp = cont_table$ai,
    n_controls_exp = cont_table$bi,
    n_cases_nexp = cont_table$ci,
    n_controls_nexp = cont_table$di,
    reverse_2x2 = reverse_phi
  )

  # miss_d = which(is.na(es$d))
  #
  # es$d[miss_d] <- (2 * phi[miss_d]) / sqrt(1 - phi[miss_d]^2)
  # X2 <- phi^2 * n_sample
  # es$d_se[miss_d] = sqrt(es$d[miss_d]^2 / X2[miss_d])
  # res_d <- .es_from_d(d = es$d, d_se = es$d_se, n_sample = n_sample, reverse = reverse_phi)
  # cols_d <- c("d", "d_se", "d_ci_lo", "d_ci_up",
  #            "g", "g_se", "g_ci_lo", "g_ci_up")
  # es[miss_d, cols_d] <- res_d[miss_d, cols_d]
  # es$r <- ifelse(reverse_phi, -phi, phi)
  # es$z <- atanh(es$r)
  # es$z_se <- sqrt(es$z^2 / X2)
  # es$z_ci_lo <- es$z - qnorm(.975) * es$z_se
  # es$z_ci_up <- es$z + qnorm(.975) * es$z_se
  # es$r_ci_lo <- tanh(es$z_ci_lo)
  # es$r_ci_up <- tanh(es$z_ci_up)
  # effective_n = 1/(es$z_se^2) + 3
  # es$r_se = sqrt((1 - es$r^2)^2 / (effective_n - 1))

  es$info_used <- "phi"

  return(es)
}

#' Convert a chi-square value to several effect size measures
#'
#' @param chisq value of the chi-squared
#' @param n_sample total number of participants in the sample
#' @param n_cases total number of cases/events
#' @param n_exp total number of participants in the exposed group
#' @param yates_chisq a logical value indicating whether the Chi square has been performed using Yate's correction for continuity.
#' @param reverse_chisq a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts a chi-square value (with one degree of freedom)
#' into a phi coefficient (Lipsey et al. 2001):
#' \deqn{phi = \sqrt{\frac{chisq^2}{n\_sample}}}.
#'
#' Note that if \code{yates_chisq = "TRUE"}, a small correction is added.
#'
#' Then, the phi coefficient is converted to other effect size measures (see \code{\link{es_from_phi}}).
#'
#' @references
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 8. Phi or chi-square'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_chisq
#'
#' @md
#'
#' @examples
#' es_from_chisq(chisq = 4.21, n_sample = 78, n_cases = 51, n_exp = 50)
es_from_chisq <- function(chisq, n_sample, n_cases, n_exp,
                          yates_chisq = FALSE,
                          reverse_chisq) {
  if (missing(reverse_chisq)) reverse_chisq <- rep(FALSE, length(chisq))
  if (missing(n_sample)) n_sample <- rep(NA, length(chisq))
  if (missing(n_exp)) n_exp <- rep(NA, length(chisq))
  if (missing(n_cases)) n_cases <- rep(NA, length(chisq))
  reverse_chisq[is.na(reverse_chisq)] <- FALSE

  if (length(reverse_chisq) == 1) reverse_chisq = c(rep(reverse_chisq, length(chisq)))
  if (length(reverse_chisq) != length(chisq)) stop("The length of the 'reverse_chisq' argument of incorrectly specified.")

  tryCatch({
    .validate_positive(n_cases, n_exp, n_sample, chisq,
                       error_message = paste0("The chi-square, number of cases, people exposed, total sample ",
                                              "should be >0."),
                       func = "es_from_chisq")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  cont_table <- suppressWarnings(
    metafor::conv.2x2(
      x2i = chisq, ni = n_sample,
      n1i = n_exp, n2i = n_cases,
      correct=yates_chisq
    )
  )

  es <- es_from_2x2(
    n_cases_exp = cont_table$ai,
    n_controls_exp = cont_table$bi,
    n_cases_nexp = cont_table$ci,
    n_controls_nexp = cont_table$di,
    reverse_2x2 = reverse_chisq
  )

  # miss_d = which(is.na(es$d))
  #
  # es$d[miss_d] <- 2 * sqrt(chisq[miss_d] / (n_sample[miss_d] - chisq[miss_d]))
  # es$d_se[miss_d] = sqrt(es$d[miss_d]^2 / chisq[miss_d])
  # res_d <- .es_from_d(d = es$d, d_se = es$d_se,
  #                     n_sample = n_sample, reverse = reverse_chisq)
  # cols_d <- c("d", "d_se", "d_ci_lo", "d_ci_up",
  #             "g", "g_se", "g_ci_lo", "g_ci_up")
  # es[miss_d, cols_d] <- res_d[miss_d, cols_d]
  # es$r <- sqrt(chisq/n_sample)
  # es$r <- ifelse(reverse_chisq, -es$r, es$r)
  # es$z <- atanh(es$r)
  # es$z_se <- sqrt(es$z^2 / chisq)
  # es$z_ci_lo <- es$z - qnorm(.975) * es$z_se
  # es$z_ci_up <- es$z + qnorm(.975) * es$z_se
  # es$r_ci_lo <- tanh(es$z_ci_lo)
  # es$r_ci_lo <- tanh(es$z_ci_up)
  # effective_n = 1/(es$z_se^2) + 3
  # es$r_se = sqrt((1 - es$r^2)^2 / (effective_n - 1))

  es$info_used <- "chisq"

  return(es)
}

#' Convert a p-value of a chi-square to several effect size measures
#'
#' @param chisq_pval p-value of a chi-square coefficient
#' @param n_sample total number of participants in the sample
#' @param n_cases total number of cases/events
#' @param n_exp total number of participants in the exposed group
#' @param yates_chisq a logical value indicating whether the Chi square has been performed using Yate's correction for continuity.
#' @param reverse_chisq_pval a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts a chi-square value (with one degree of freedom)
#' into a chi-square coefficient (Section 3.12 in Lipsey et al., 2001):
#' \deqn{chisq = qchisq(chisq\_pval, df = 1, lower.tail = FALSE)}
#'
#' Note that if \code{yates_chisq = "TRUE"}, a small correction is added.
#'
#' Then, the chisq coefficient is converted to other effect size measures (see \code{\link{es_from_chisq}}).
#'
#' @references
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 8. Phi or chi-square'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_chisq_pval
#'
#' @md
#'
#' @examples
#' es_from_chisq_pval(chisq_pval = 0.2, n_sample = 42, n_exp = 25, n_cases = 13)
es_from_chisq_pval <- function(chisq_pval, n_sample, n_cases, n_exp,
                               yates_chisq = FALSE,
                               reverse_chisq_pval) {
  if (missing(reverse_chisq_pval)) reverse_chisq_pval <- rep(FALSE, length(chisq_pval))
  reverse_chisq_pval[is.na(reverse_chisq_pval)] <- FALSE

  tryCatch({
    .validate_positive(n_cases, n_exp, n_sample, chisq_pval,
                       error_message = paste0("The chi-square p-value, number of cases, people exposed, total sample ",
                                              "should be >0."),
                       func = "es_from_chisq_pval")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  chisq <- stats::qchisq(p = chisq_pval, df = 1, lower.tail = FALSE)

  es <- es_from_chisq(
    chisq = chisq, n_sample = n_sample,
    yates_chisq = yates_chisq,
    n_cases = n_cases,
    n_exp = n_exp,
    reverse_chisq = reverse_chisq_pval
  )

  es$info_used <- "chisq_pval"
  return(es)
}

# es_d <- .es_from_d(d = d, d_se = d_se, n_sample = n_sample,
#                    reverse = reverse_chisq)
# OR
# px = n_exp/n_sample
# py = n_cases/n_sample
# n_cases_exp = suppressWarnings(
#   n_sample * px * py *
#     sqrt((chisq * px * py * (1-px) * (1-py))/n_sample))
#
# n_controls_exp = n_exp - n_cases_exp
# n_cases_nexp = n_cases - n_cases_exp
# n_controls_nexp = n_sample - (n_cases_exp + n_controls_exp + n_cases_nexp)
#
# es = es_from_2x2(n_cases_exp = n_cases_exp, n_controls_exp = n_controls_exp,
#                  n_cases_nexp = n_cases_nexp, n_controls_nexp = n_controls_nexp,
#                  reverse_2x2 = reverse_chisq)
#
# col.es = which(colnames(es) %in% c("d", "d_se", "d_ci_lo", "d_ci_up",
#                                    "g", "g_se", "g_ci_lo", "g_ci_up"))
# col.esd = which(colnames(es_d) %in% c("d", "d_se", "d_ci_lo", "d_ci_up",
#                                    "g", "g_se", "g_ci_lo", "g_ci_up"))
#
# es[, col.es] <- es_d[, col.esd]
#
# #COR
# r = r_se = r_ci_lo = r_ci_up =
#   z = z_se = z_ci_lo = z_ci_up = rep(NA, length(chisq))
#
# dat_cor = data.frame(n_cases_exp = n_cases_exp, n_controls_exp = n_controls_exp,
#                      n_cases_nexp = n_cases_nexp, n_controls_nexp = n_controls_nexp,
#                      n_sample = n_sample, reverse_chisq = reverse_chisq,
#                      chisq = chisq, chisq_to_cor = chisq_to_cor)
#
# nn_miss = which(
#   (!is.na(n_cases_exp) & !is.na(n_controls_exp) & !is.na(n_cases_nexp) & !is.na(n_controls_nexp)) |
#     (!is.na(chisq) & !is.na(n_sample))
# )
#
# if (length(nn_miss) != 0) {
#   cor = t(mapply(.chi_to_cor,
#                  chisq = dat_cor$chisq[nn_miss],
#                  n_sample = dat_cor$n_sample[nn_miss],
#                  n_cases_exp = dat_cor$n_cases_exp[nn_miss],
#                  n_controls_exp = dat_cor$n_controls_exp[nn_miss],
#                  n_cases_nexp = dat_cor$n_cases_nexp[nn_miss],
#                  n_controls_nexp = dat_cor$n_controls_nexp[nn_miss],
#                  reverse_chisq = dat_cor$reverse_chisq[nn_miss],
#                  chisq_to_cor = dat_cor$chisq_to_cor[nn_miss]))
#
#   r[nn_miss] = cor[, 1]
#   r_se[nn_miss] = sqrt(cor[, 2])
#   r_ci_lo[nn_miss] = cor[, 3]
#   r_ci_up[nn_miss] = cor[, 4]
#   z[nn_miss] = cor[, 5]
#   z_se[nn_miss] = sqrt(cor[, 6])
#   z_ci_lo[nn_miss] = cor[, 7]
#   z_ci_up[nn_miss] = cor[, 8]
# }
#
# col = which(colnames(es) %in% c("r", "r_se", "r_ci_lo", "r_ci_up",
#                                 "z", "z_se", "z_ci_lo", "z_ci_up"))
# es[, col] <- cbind(r, r_se, r_ci_lo, r_ci_up,
#                   z, z_se, z_ci_lo, z_ci_up)







# es_d <- .es_from_d(d = d, d_se = d_se, n_sample = n_sample, reverse = reverse_phi)

# px = n_exp/n_sample
# py = n_cases/n_sample
# n_cases_exp = n_sample * suppressWarnings(px*py + phi*sqrt(px*py*(1-px)*(1-py)))
# n_controls_exp = n_exp - n_cases_exp
# n_cases_nexp = n_cases - n_cases_exp
# n_controls_nexp = n_sample - (n_cases_exp + n_controls_exp + n_cases_nexp)
#
# es = es_from_2x2(n_cases_exp = n_cases_exp, n_controls_exp = n_controls_exp,
#                  n_cases_nexp = n_cases_nexp, n_controls_nexp = n_controls_nexp,
#                  reverse_2x2 = reverse_phi)
# col.es = which(colnames(es) %in% c("d", "d_se", "d_ci_lo", "d_ci_up",
#                                    "g", "g_se", "g_ci_lo", "g_ci_up"))
# col.esd = which(colnames(es_d) %in% c("d", "d_se", "d_ci_lo", "d_ci_up",
#                                       "g", "g_se", "g_ci_lo", "g_ci_up"))
#
# es[, col.es] <- es_d[, col.esd]
# #COR
# r = r_se = r_ci_lo = r_ci_up =
# z = z_se = z_ci_lo = z_ci_up = rep(NA, length(phi))
#
# dat_cor = data.frame(n_cases_exp = n_cases_exp, n_controls_exp = n_controls_exp,
#                      n_cases_nexp = n_cases_nexp, n_controls_nexp = n_controls_nexp,
#                      n_sample = n_sample, reverse_phi = reverse_phi,
#                      phi = phi, phi_to_cor = phi_to_cor)
#
# nn_miss = which(
#   (!is.na(n_cases_exp) & !is.na(n_controls_exp) & !is.na(n_cases_nexp) & !is.na(n_controls_nexp)) |
#     (!is.na(phi) & !is.na(n_sample))
# )
#
# if (length(nn_miss) != 0) {
#   cor = t(mapply(.phi_to_cor,
#                  phi = dat_cor$phi[nn_miss],
#                  n_sample = dat_cor$n_sample[nn_miss],
#                  n_cases_exp = dat_cor$n_cases_exp[nn_miss],
#                  n_controls_exp = dat_cor$n_controls_exp[nn_miss],
#                  n_cases_nexp = dat_cor$n_cases_nexp[nn_miss],
#                  n_controls_nexp = dat_cor$n_controls_nexp[nn_miss],
#                  reverse_phi = dat_cor$reverse_phi[nn_miss],
#                  phi_to_cor = dat_cor$phi_to_cor[nn_miss]))
#
#   r[nn_miss] = cor[, 1]
#   r_se[nn_miss] = sqrt(cor[, 2])
#   r_ci_lo[nn_miss] = cor[, 3]
#   r_ci_up[nn_miss] = cor[, 4]
#   z[nn_miss] = cor[, 5]
#   z_se[nn_miss] = sqrt(cor[, 6])
#   z_ci_lo[nn_miss] = cor[, 7]
#   z_ci_up[nn_miss] = cor[, 8]
# }
#
# col = which(colnames(es) %in% c("r", "r_se", "r_ci_lo", "r_ci_up",
#                                 "z", "z_se", "z_ci_lo", "z_ci_up"))
# es[, col] <- cbind(r, r_se, r_ci_lo, r_ci_up,
#                    z, z_se, z_ci_lo, z_ci_up)

