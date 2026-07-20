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
#' **To estimate D, G and the correlation coefficients (R/Z) when the 2x2 table cannot be
#' reconstructed** (e.g., the number of cases or exposed participants is missing), the phi
#' coefficient -- which is the Pearson correlation of the two binary variables -- is
#' treated directly as a correlation coefficient and passed to
#' \code{\link{es_from_pearson_r}()}. This yields R and Z with their standard large-sample
#' sampling variances,
#' \deqn{r = phi, \quad r\_se = \frac{1 - r^2}{\sqrt{n\_sample - 1}}}
#' \deqn{z = atanh(r), \quad z\_se = \frac{1}{\sqrt{n\_sample - 3}}}
#' and D and G are then obtained from R exactly as in \code{\link{es_from_pearson_r}()}
#' (the R to D step assumes balanced groups). In this situation the OR, RR and NNT require
#' the 2x2 margins and are returned as NA. When the 2x2 table can be reconstructed, D, G,
#' R and Z are instead derived from it (see \code{\link{es_from_2x2}()}).
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
  if (length(reverse_phi) != length(phi)) stop("The length of the 'reverse_phi' argument is incorrectly specified.")

  # A phi coefficient is a correlation and must lie in [-1, 1]; metafor::conv.2x2()
  # raises a HARD ERROR (not a warning) on |phi| > 1, which would abort an entire
  # convert_df() run over the other (valid) rows. Guard here so a single out-of-range
  # value degrades to an all-NA row (the package's one-bad-cell-cannot-abort contract).
  # convert_df()'s Tier-1 validation already NA's such values under correct_inputs =
  # TRUE; this is the backstop for direct es_from_phi() calls.
  phi_valid <- phi
  phi_valid[!is.na(phi_valid) & abs(phi_valid) > 1] <- NA_real_

  cont_table <- suppressWarnings(
    metafor::conv.2x2(
      ri = phi_valid, ni = n_sample,
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

  # Fallback for rows where the 2x2 table could NOT be reconstructed (n_cases / n_exp
  # missing, so conv.2x2 returns NA cells and OR/RR/NNT are NA): a phi coefficient IS the
  # Pearson correlation of two binary variables, so deliver R/Z/D/G by treating phi as a
  # Pearson r and delegating to the vetted es_from_pearson_r() engine. This gives the
  # STANDARD sampling variances (z_se = 1/sqrt(n-3), r_se = (1-r^2)/sqrt(n-1), matching
  # metafor) rather than the earlier ad hoc formula that inflated the SE by up to ~22% at
  # large phi. A phi-only row therefore yields the SAME R/Z/D/G as an equivalent
  # pearson_r row (internal consistency). NB: the R -> D step assumes balanced groups;
  # OR/RR/NNT require the 2x2 margins and remain NA when they are missing.
  miss <- which(is.na(es$r) & !is.na(phi) & !is.na(n_sample) & abs(phi) < 1)
  if (length(miss) > 0) {
    pr <- es_from_pearson_r(
      pearson_r = phi[miss], n_sample = n_sample[miss],
      reverse_pearson_r = reverse_phi[miss]
    )
    fb_cols <- intersect(
      c("d", "d_se", "d_ci_lo", "d_ci_up", "g", "g_se", "g_ci_lo", "g_ci_up",
        "r", "r_se", "r_ci_lo", "r_ci_up", "z", "z_se", "z_ci_lo", "z_ci_up"),
      intersect(colnames(es), colnames(pr))
    )
    es[miss, fb_cols] <- pr[, fb_cols]
  }

  es$info_used <- "phi"

  return(es)
}

#' Convert a chi-square value to several effect size measures
#'
#' @param chisq value of the chi-squared
#' @param n_sample total number of participants in the sample
#' @param n_cases total number of cases/events
#' @param n_exp total number of participants in the exposed group
#' @param yates_chisq logical value (or vector of length \code{length(chisq)})
#'   indicating whether the chi-square has been performed using Yates' correction
#'   for continuity. When a vector is supplied, each row is back-transformed using
#'   its own setting. Missing values are treated as \code{FALSE}.
#' @param reverse_chisq a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts a chi-square value (with one degree of freedom)
#' into a phi coefficient (Lipsey et al. 2001):
#' \deqn{phi = \sqrt{\frac{chisq}{n\_sample}}}
#' and then converts it to other effect size measures exactly as in
#' \code{\link{es_from_phi}()} (including the correlation-based R/Z/D/G fallback with
#' standard sampling variances when the 2x2 table cannot be reconstructed).
#'
#' Note that if \code{yates_chisq = TRUE}, the chi-square value is interpreted
#' as Yates-corrected when back-transforming to a 2x2 contingency table; this is
#' propagated row by row when a vector is supplied.
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
  if (length(reverse_chisq) != length(chisq)) stop("The length of the 'reverse_chisq' argument is incorrectly specified.")

  if (length(yates_chisq) == 1) yates_chisq <- rep(yates_chisq, length(chisq))
  if (length(yates_chisq) != length(chisq)) stop("The length of the 'yates_chisq' argument is incorrectly specified.")
  yates_chisq[is.na(yates_chisq)] <- FALSE

  cont_table <- suppressWarnings(
    metafor::conv.2x2(
      x2i = chisq, ni = n_sample,
      n1i = n_exp, n2i = n_cases,
      correct = yates_chisq
    )
  )

  es <- es_from_2x2(
    n_cases_exp = cont_table$ai,
    n_controls_exp = cont_table$bi,
    n_cases_nexp = cont_table$ci,
    n_controls_nexp = cont_table$di,
    reverse_2x2 = reverse_chisq
  )

  # Fallback for rows where the 2x2 table could NOT be reconstructed (n_cases / n_exp
  # missing). For a 1-df chi-square, phi = sqrt(chisq / n) is the Pearson correlation of
  # the two binary variables (its sign is not identified by chisq alone, so it is taken
  # positive and flipped by reverse_chisq). Deliver R/Z/D/G by treating phi as a Pearson
  # r via the vetted es_from_pearson_r() engine, giving the STANDARD SEs (matching
  # metafor) rather than the earlier ad hoc, SE-inflating formula. OR/RR/NNT require the
  # 2x2 margins and remain NA; the R -> D step assumes balanced groups.
  r_chi <- suppressWarnings(sqrt(chisq / n_sample))
  miss <- which(is.na(es$r) & !is.na(r_chi) & is.finite(r_chi) & r_chi < 1 & !is.na(n_sample))
  if (length(miss) > 0) {
    pr <- es_from_pearson_r(
      pearson_r = r_chi[miss], n_sample = n_sample[miss],
      reverse_pearson_r = reverse_chisq[miss]
    )
    fb_cols <- intersect(
      c("d", "d_se", "d_ci_lo", "d_ci_up", "g", "g_se", "g_ci_lo", "g_ci_up",
        "r", "r_se", "r_ci_lo", "r_ci_up", "z", "z_se", "z_ci_lo", "z_ci_up"),
      intersect(colnames(es), colnames(pr))
    )
    es[miss, fb_cols] <- pr[, fb_cols]
  }

  es$info_used <- "chisq"

  return(es)
}

#' Convert a p-value of a chi-square to several effect size measures
#'
#' @param chisq_pval p-value of a chi-square coefficient
#' @param n_sample total number of participants in the sample
#' @param n_cases total number of cases/events
#' @param n_exp total number of participants in the exposed group
#' @param yates_chisq logical value (or vector of length \code{length(chisq_pval)})
#'   indicating whether the chi-square has been performed using Yates' correction
#'   for continuity. When a vector is supplied, each row is back-transformed using
#'   its own setting. Missing values are treated as \code{FALSE}.
#' @param reverse_chisq_pval a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts a chi-square value (with one degree of freedom)
#' into a chi-square coefficient (Section 3.12 in Lipsey et al., 2001):
#' \deqn{chisq = qchisq(chisq\_pval, df = 1, lower.tail = FALSE)}
#'
#' Note that if \code{yates_chisq = TRUE}, the chi-square value is interpreted
#' as Yates-corrected when back-transforming to a 2x2 contingency table; this is
#' propagated row by row when a vector is supplied.
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

