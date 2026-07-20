#' Convert two paired t-test value of two independent groups into several effect size measures
#'
#' @param paired_t_exp Paired t-test value of the experimental/exposed group.
#' @param paired_t_nexp Paired t-test value of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the paired t-test value into a SMD (see details).
#' @param reverse_paired_t a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts paired t-tests of two independent groups value into a Cohen's d (D) and Hedges' g (G) (table 12.2 in Cooper).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the Cohen's d,** the following formulas are used (Cooper et al., 2019):
#' \deqn{cohen\_d\_exp = paired\_t\_exp * \sqrt{\frac{2 * (1 - r\_pre\_post\_exp)}{n\_exp}}}
#' \deqn{cohen\_d\_nexp = paired\_t\_nexp * \sqrt{\frac{2 * (1 - r\_pre\_post\_nexp)}{n\_nexp}}}
#' \deqn{cohen\_d\_se\_exp = \sqrt{\frac{2 * (1 - r\_pre\_post\_exp)}{n\_exp} + \frac{d\_exp^2}{2 * n\_exp}}}
#' \deqn{cohen\_d\_se\_nexp = \sqrt{\frac{2 * (1 - r\_pre\_post\_nexp)}{n\_nexp} + \frac{d\_nexp^2}{2 * n\_nexp}}}
#' \deqn{cohen\_d = d\_exp - d\_nexp}
#' \deqn{d\_se = \sqrt{cohen\_d\_se\_exp^2 + cohen\_d\_se\_nexp^2}}
#'
#' When \code{pre_post_to_smd = "morris_dz"}, the mean difference is standardized by the standard deviation
#' of the change score and the pre-post correlation is no longer involved (Morris & DeShon, 2002):
#' \deqn{cohen\_d\_exp = \frac{paired\_t\_exp}{\sqrt{n\_exp}}}
#' \deqn{cohen\_d\_nexp = \frac{paired\_t\_nexp}{\sqrt{n\_nexp}}}
#'
#' Note that the Cohen's d obtained from a paired t-test strongly depends on the pre-post correlation.
#' When \code{r_pre_post_exp} / \code{r_pre_post_nexp} are not indicated, a value of 0.8 is assumed and
#' users should conduct sensitivity analyses with other plausible values.
#'
#' **No \code{pool_sd} argument.** Unlike the means and mean-change converters, this function cannot
#' pool the standardizing SD across arms: a paired t-statistic identifies each arm's
#' \eqn{mean\_change / sd\_change} ratio but NOT the ratio of the two arms' SDs, so the pooled
#' standardizer is not recoverable from the reported statistic. Each arm is therefore standardized by
#' its own SD and the two within-group values are subtracted -- a construction that coincides with the
#' pooled one only when the arms' SDs are equal. When rows from this route are combined in one pool
#' with rows that DO use a pooled standardizer, \code{summary(..., flags = TRUE)} raises an
#' informational flag. If the arm SDs are reported, prefer \code{\link{es_from_means_sd_pre_post}} or
#' \code{\link{es_from_mean_change_sd}}.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 16. Paired: Paired F- or t-test'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @export es_from_paired_t
#'
#' @md
#'
#' @examples
#' es_from_paired_t(paired_t_exp = 2.1, paired_t_nexp = 4.2, n_exp = 20, n_nexp = 22)
es_from_paired_t <- function(paired_t_exp, paired_t_nexp, n_exp, n_nexp,
                             r_pre_post_exp, r_pre_post_nexp,
                             smd_to_cor = "viechtbauer",
                             pre_post_to_smd = "cooper",
                             reverse_paired_t) {
  if (missing(reverse_paired_t)) reverse_paired_t <- rep(FALSE, length(paired_t_exp))
  reverse_paired_t[is.na(reverse_paired_t)] <- FALSE
  if (length(reverse_paired_t) == 1) reverse_paired_t = c(rep(reverse_paired_t, length(paired_t_exp)))
  if (length(reverse_paired_t) != length(paired_t_exp)) stop("The length of the 'reverse_paired_t' argument is incorrectly specified.")

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "paired_t",
    func_name = "es_from_paired_t"
  )

  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(paired_t_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(paired_t_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  r_pre_post_exp <- .guard_r_pre_post(r_pre_post_exp)
  r_pre_post_nexp <- .guard_r_pre_post(r_pre_post_nexp)

  J_exp <- .d_j(n_exp - 1)
  J_nexp <- .d_j(n_nexp - 1)

  # per-row method selection: a vector mixing "morris_dz" and "morris_drm" rows
  # previously fell wholesale into the drm branch (all(...) is a single logical).
  # Recycled to the input length so a scalar method still yields one value per row.
  dz <- rep(pre_post_to_smd == "morris_dz", length.out = length(paired_t_exp))

  # morris dz. metafor SMCC convention (variance built from the corrected g),
  # matching .single_group_pre_post_to_smd so the paired-t and mean-change
  # routes return identical SEs on equivalent inputs
  d_exp_dz <- paired_t_exp / sqrt(n_exp)
  d_nexp_dz <- paired_t_nexp / sqrt(n_nexp)
  d_var_exp_dz <- (1 / n_exp + (J_exp * d_exp_dz)^2 / (2 * n_exp)) / J_exp^2
  d_var_nexp_dz <- (1 / n_nexp + (J_nexp * d_nexp_dz)^2 / (2 * n_nexp)) / J_nexp^2

  # morris drm
  d_exp_drm <- paired_t_exp * sqrt((2 * (1 - r_pre_post_exp)) / n_exp)
  d_nexp_drm <- paired_t_nexp * sqrt((2 * (1 - r_pre_post_nexp)) / n_nexp)
  d_var_exp_drm <- 2 * (1 - r_pre_post_exp) / n_exp + d_exp_drm^2 / (2 * n_exp)
  d_var_nexp_drm <- 2 * (1 - r_pre_post_nexp) / n_nexp + d_nexp_drm^2 / (2 * n_nexp)

  d_exp <- ifelse(dz, d_exp_dz, d_exp_drm)
  d_nexp <- ifelse(dz, d_nexp_dz, d_nexp_drm)
  d_var_exp <- ifelse(dz, d_var_exp_dz, d_var_exp_drm)
  d_var_nexp <- ifelse(dz, d_var_nexp_dz, d_var_nexp_drm)

  # A within-subject arm needs n >= 2 to have any estimable variance. The morris_drm
  # branch carries no J(n-1) factor, so unlike morris_dz it stays finite at n = 1 (and
  # blows up to Inf at n = 0), leaking a spurious d/g/logOR/r downstream. NA such arms
  # so this route matches the mean-change route (which NAs via the n >= 2 kernel guard).
  bad_exp <- !(is.finite(n_exp) & n_exp >= 2)
  bad_nexp <- !(is.finite(n_nexp) & n_nexp >= 2)
  d_exp[bad_exp] <- NA_real_; d_var_exp[bad_exp] <- NA_real_
  d_nexp[bad_nexp] <- NA_real_; d_var_nexp[bad_nexp] <- NA_real_

  g_exp <- J_exp * d_exp
  g_nexp <- J_nexp * d_nexp

  g_var_exp <- J_exp^2 * d_var_exp
  g_var_nexp <- J_nexp^2 * d_var_nexp

  d <- d_exp - d_nexp
  d_se <- sqrt(d_var_exp + d_var_nexp)

  es <- .es_from_d(
    d = d, d_se = d_se, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_paired_t
  )

  es$g <- ifelse(reverse_paired_t, g_nexp - g_exp, g_exp - g_nexp)
  es$g_se <- sqrt(g_var_exp + g_var_nexp)
  es$g_ci_lo <- es$g - qt(.975, n_exp + n_nexp - 2) * es$g_se
  es$g_ci_up <- es$g + qt(.975, n_exp + n_nexp - 2) * es$g_se

  es$info_used <- "paired_t"
  return(es)
}

#' Convert two paired t-test p-value obtained from two independent groups into several effect size measures
#'
#' @param paired_t_pval_exp P-value of the paired t-test value of the experimental/exposed group.
#' @param paired_t_pval_nexp P-value of the paired t-test value of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the paired t-test value into a SMD (see details).
#' @param reverse_paired_t_pval a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the p-values of two paired t-test obtained from two independent groups value into a Cohen's d (D) and Hedges' g (G) (table 12.2 in Cooper).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the Cohen's d,** the following formulas are used (Cooper et al., 2019):
#' This function converts a Student's t-test value into a Cohen's d (table 12.2 in Cooper).
#' \deqn{paired\_t\_exp = qt(\frac{paired\_t\_pval\_exp}{2}, df = n\_exp - 1) * \sqrt{\frac{2 * (1 - r\_pre\_post\_exp)}{n\_exp}}}
#' \deqn{paired\_t\_nexp = qt(\frac{paired\_t\_pval\_nexp}{2}, df = n\_nexp - 1) * \sqrt{\frac{2 * (1 - r\_pre\_post\_nexp)}{n\_nexp}}}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 16. Paired: Paired F- or t-test'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_paired_t_pval
#'
#' @md
#'
#' @examples
#' es_from_paired_t_pval(paired_t_pval_exp = 0.4, paired_t_pval_nexp = 0.01, n_exp = 19, n_nexp = 22)
es_from_paired_t_pval <- function(paired_t_pval_exp, paired_t_pval_nexp, n_exp, n_nexp,
                                  r_pre_post_exp, r_pre_post_nexp,
                                  smd_to_cor = "viechtbauer",
                                  pre_post_to_smd = "cooper",
                                  reverse_paired_t_pval) {
  if (missing(reverse_paired_t_pval)) reverse_paired_t_pval <- rep(FALSE, length(paired_t_pval_exp))
  reverse_paired_t_pval[is.na(reverse_paired_t_pval)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(paired_t_pval_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(paired_t_pval_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  paired_t_exp <- qt(p = paired_t_pval_exp / 2, df = n_exp - 1, lower.tail = FALSE)

  paired_t_nexp <- qt(p = paired_t_pval_nexp / 2, df = n_nexp - 1, lower.tail = FALSE)

  es <- es_from_paired_t(
    paired_t_exp = paired_t_exp, paired_t_nexp = paired_t_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_paired_t = reverse_paired_t_pval
  )

  es$info_used <- "paired_t_pval"
  return(es)
}

#' Convert two paired ANOVA f value of two independent groups into several effect size measures
#'
#' @param paired_f_exp Paired ANOVA F value of the experimental/exposed group.
#' @param paired_f_nexp Paired ANOVA F value of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the paired t-test value into a SMD (see details).
#' @param reverse_paired_f a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the paired F-test obtained from two independent groups value into a Cohen's d (D) and Hedges' g (G) (table 12.2 in Cooper).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the Cohen's d,** the following formulas are used (Cooper et al., 2019):
#' This function converts a Student's t-test value into a Cohen's d (table 12.2 in Cooper).
#' \deqn{paired\_t\_exp = \sqrt{paired\_f\_exp}}
#' \deqn{paired\_t\_nexp = \sqrt{paired\_f\_nexp}}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_paired_t}()} are applied.
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 16. Paired: Paired F- or t-test'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_paired_f
#'
#' @md
#'
#' @examples
#' es_from_paired_f(paired_f_exp = 2.1, paired_f_nexp = 4.2, n_exp = 20, n_nexp = 22)
es_from_paired_f <- function(paired_f_exp, paired_f_nexp, n_exp, n_nexp,
                             r_pre_post_exp, r_pre_post_nexp,
                             smd_to_cor = "viechtbauer",
                             pre_post_to_smd = "cooper",
                             reverse_paired_f) {
  if (missing(reverse_paired_f)) reverse_paired_f <- rep(FALSE, length(paired_f_exp))
  reverse_paired_f[is.na(reverse_paired_f)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(paired_f_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(paired_f_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8


  paired_t_exp <- sqrt(paired_f_exp)

  paired_t_nexp <- sqrt(paired_f_nexp)

  es <- es_from_paired_t(
    paired_t_exp = paired_t_exp, paired_t_nexp = paired_t_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_paired_t = reverse_paired_f
  )

  es$info_used <- "paired_f"
  return(es)
}

#' Convert two paired ANOVA f p-value of two independent groups into several effect size measures
#'
#' @param paired_f_pval_exp P-value of the paired ANOVA F of the experimental/exposed group.
#' @param paired_f_pval_nexp P-value of the paired ANOVA F of the non-experimental/non-exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param r_pre_post_exp pre-post correlation in the experimental/exposed group
#' @param r_pre_post_nexp pre-post correlation in the non-experimental/non-exposed group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the paired t-test value into a SMD (see details).
#' @param reverse_paired_f_pval a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the p-values of two paired F-test obtained from two independent groups value into a Cohen's d (D) and Hedges' g (G) (table 12.2 in Cooper).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the Cohen's d,** the following formulas are used (Cooper et al., 2019):
#' This function converts a Student's t-test value into a Cohen's d (table 12.2 in Cooper).
#' \deqn{paired\_t\_exp = qt(\frac{paired\_f\_pval\_exp}{2}, df = n\_exp - 1) * \sqrt{\frac{2 * (1 - r_pre_post_exp)}{n_exp}}}
#' \deqn{paired\_t\_nexp = qt(\frac{paired\_f\_pval\_nexp}{2}, df = n\_nexp - 1) * \sqrt{\frac{2 * (1 - r_pre_post_nexp)}{n_nexp}}}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_paired_t}()} are applied.
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 16. Paired: Paired F- or t-test'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_paired_f_pval
#'
#' @md
#'
#' @examples
#' es_from_paired_f_pval(paired_f_pval_exp = 0.4, paired_f_pval_nexp = 0.01, n_exp = 19, n_nexp = 22)
es_from_paired_f_pval <- function(paired_f_pval_exp, paired_f_pval_nexp, n_exp, n_nexp,
                                  r_pre_post_exp, r_pre_post_nexp,
                                  smd_to_cor = "viechtbauer",
                                  pre_post_to_smd = "cooper",
                                  reverse_paired_f_pval) {
  if (missing(reverse_paired_f_pval)) reverse_paired_f_pval <- rep(FALSE, length(paired_f_pval_exp))
  reverse_paired_f_pval[is.na(reverse_paired_f_pval)] <- FALSE
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(paired_f_pval_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(paired_f_pval_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  es <- es_from_paired_t_pval(
    paired_t_pval_exp = paired_f_pval_exp,
    paired_t_pval_nexp = paired_f_pval_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp,
    r_pre_post_nexp = r_pre_post_nexp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_paired_t_pval = reverse_paired_f_pval
  )

  es$info_used <- "paired_f_pval"
  return(es)
}
