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
#' \eqn{mean\_change / sd\_change} ratio but not the ratio of the two arms' SDs, so the pooled
#' standardizer is not recoverable from the reported statistic. Each arm is therefore standardized by
#' its own SD and the two within-group values are subtracted, a construction that coincides with the
#' pooled one only when the arms' SDs are equal. When rows from this route are combined in one pool
#' with rows that do use a pooled standardizer, \code{summary(..., flags = TRUE)} raises an
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

  # The morris_dz and morris_drm arithmetic is delegated to
  # .single_group_pre_post_to_smd() rather than written out again here, so that the
  # paired-t and mean-change routes return identical SEs by construction rather than by
  # both copies being kept in step.
  res_exp  <- .paired_t_to_smd(paired_t_exp,  n_exp,  r_pre_post_exp,  pre_post_to_smd)
  res_nexp <- .paired_t_to_smd(paired_t_nexp, n_nexp, r_pre_post_nexp, pre_post_to_smd)

  d_exp <- res_exp[, "d"];   d_var_exp <- res_exp[, "var_d"]
  d_nexp <- res_nexp[, "d"]; d_var_nexp <- res_nexp[, "var_d"]
  g_exp <- res_exp[, "g"];   g_var_exp <- res_exp[, "var_g"]
  g_nexp <- res_nexp[, "g"]; g_var_nexp <- res_nexp[, "var_g"]

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
#' @param reverse_paired_t_pval_exp a logical value indicating that the experimental/exposed group changed in the
#'   NEGATIVE direction. A two-sided p-value is unsigned, so this is the only way to state that arm's own
#'   direction (see details).
#' @param reverse_paired_t_pval_nexp a logical value indicating that the non-experimental/non-exposed group
#'   changed in the NEGATIVE direction (see details).
#'
#' @details
#' This function converts the p-values of two paired t-test obtained from two independent groups value into a Cohen's d (D) and Hedges' g (G) (table 12.2 in Cooper).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the Cohen's d,** the p-values are first converted into paired t-test values:
#' \deqn{paired\_t\_exp = qt(1 - \frac{paired\_t\_pval\_exp}{2}, df = n\_exp - 1)}
#' \deqn{paired\_t\_nexp = qt(1 - \frac{paired\_t\_pval\_nexp}{2}, df = n\_nexp - 1)}
#'
#' which are then converted into a Cohen's d (Cooper et al., 2019):
#' \deqn{cohen\_d\_exp = paired\_t\_exp * \sqrt{\frac{2 * (1 - r\_pre\_post\_exp)}{n\_exp}}}
#' \deqn{cohen\_d\_nexp = paired\_t\_nexp * \sqrt{\frac{2 * (1 - r\_pre\_post\_nexp)}{n\_nexp}}}
#'
#' **A two-sided p-value carries no sign.** The quantile above returns \eqn{|t|}, so only the
#' MAGNITUDE of each arm's change is recoverable and, unless told otherwise, this function assumes
#' that both arms moved in the SAME direction. The between-group estimate is formed as
#' \eqn{cohen\_d = d\_exp - d\_nexp}; when the control arm actually moved the other way -- it
#' deteriorated while the treatment arm improved, an ordinary trial result -- the correct contrast
#' is \eqn{|d\_exp| + |d\_nexp|}, so the assumption biases the estimate by \eqn{2|d\_nexp|}.
#' "Biased by 2|d_nexp|" is not the same as "halved": what the function returns is
#' \eqn{\frac{|d\_exp| - |d\_nexp|}{|d\_exp| + |d\_nexp|}}{(|d_exp| - |d_nexp|)/(|d_exp| + |d_nexp|)}
#' of the true contrast. It therefore collapses to exactly ZERO when the two arms' magnitudes are
#' equal, is halved only in the particular case \eqn{|d\_exp| = 3|d\_nexp|}, and reverses sign when
#' the control arm moved more than the treatment arm. Measured on es_from_paired_f() at
#' n_exp = n_nexp = 50 and r_pre_post = 0.8: two equal F values of 9 return d = 0 where the true
#' contrast is 0.5367; F_exp = 81 against F_nexp = 9 returns 0.5367 against a true 1.0733 (the exact
#' halving); and F_exp = 9 against F_nexp = 81 returns -0.5367 against a true +1.0733. The
#' assumption cannot be verified from the reported statistic, so no quality flag fires on such a row.
#'
#' Mark the arm that moved down with \code{reverse_paired_t_pval_exp} /
#' \code{reverse_paired_t_pval_nexp}, which flip the sign of that arm's recovered t before the
#' contrast is formed. \code{reverse_paired_t_pval} flips the whole contrast instead, so it cannot
#' express two arms moving apart. Both per-arm flags are registered input columns, so
#' they are usable from \code{\link{convert_df}()} as well as on a direct call. When the signed t values are reported,
#' \code{\link{es_from_paired_t}()} needs no such assumption, and the
#' \code{\link{es_from_mean_change_pval}()} family carries the direction in the sign of the mean
#' change.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @note
#' Both arms are assumed to have changed in the same direction, and a two-sided p-value can neither
#' confirm nor refute that. Use \code{\link{es_from_paired_t}()} with signed t values, the
#' \code{es_from_mean_change_*} family, or the per-arm \code{reverse_paired_t_pval_exp} /
#' \code{reverse_paired_t_pval_nexp} flags when the arms diverge.
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
                                  reverse_paired_t_pval,
                                  reverse_paired_t_pval_exp,
                                  reverse_paired_t_pval_nexp) {
  if (missing(reverse_paired_t_pval)) reverse_paired_t_pval <- rep(FALSE, length(paired_t_pval_exp))
  reverse_paired_t_pval[is.na(reverse_paired_t_pval)] <- FALSE
  reverse_paired_t_pval_exp <- .arm_direction_flag(
    if (missing(reverse_paired_t_pval_exp)) NULL else reverse_paired_t_pval_exp,
    length(paired_t_pval_exp), "reverse_paired_t_pval_exp"
  )
  reverse_paired_t_pval_nexp <- .arm_direction_flag(
    if (missing(reverse_paired_t_pval_nexp)) NULL else reverse_paired_t_pval_nexp,
    length(paired_t_pval_exp), "reverse_paired_t_pval_nexp"
  )
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(paired_t_pval_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(paired_t_pval_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  paired_t_exp <- qt(p = paired_t_pval_exp / 2, df = n_exp - 1, lower.tail = FALSE)

  paired_t_nexp <- qt(p = paired_t_pval_nexp / 2, df = n_nexp - 1, lower.tail = FALSE)

  # A two-sided p-value is unsigned, so both quantiles above are |t|. The per-arm flags
  # are the only way to say that an arm moved down; without them the contrast
  # d_exp - d_nexp assumes the two arms moved the same way (see the details section).
  paired_t_exp <- ifelse(reverse_paired_t_pval_exp, -paired_t_exp, paired_t_exp)
  paired_t_nexp <- ifelse(reverse_paired_t_pval_nexp, -paired_t_nexp, paired_t_nexp)

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
#' @param reverse_paired_f_exp a logical value indicating that the experimental/exposed group changed in the
#'   NEGATIVE direction. A paired F value is unsigned, so this is the only way to state that arm's own
#'   direction (see details).
#' @param reverse_paired_f_nexp a logical value indicating that the non-experimental/non-exposed group
#'   changed in the NEGATIVE direction (see details).
#'
#' @details
#' This function converts the paired F-test obtained from two independent groups value into a Cohen's d (D) and Hedges' g (G) (table 12.2 in Cooper).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the Cohen's d,** the paired F values are first converted into paired t-test values:
#' \deqn{paired\_t\_exp = \sqrt{paired\_f\_exp}}
#' \deqn{paired\_t\_nexp = \sqrt{paired\_f\_nexp}}
#'
#' **A paired F value carries no sign.** It is the square of the paired t, so the square roots above
#' return \eqn{|t|} and only the MAGNITUDE of each arm's change is recoverable. Unless told otherwise,
#' this function therefore assumes that both arms moved in the SAME direction. The between-group
#' estimate is formed as \eqn{cohen\_d = d\_exp - d\_nexp}; when the control arm actually moved the
#' other way -- it deteriorated while the treatment arm improved, an ordinary trial result -- the
#' correct contrast is \eqn{|d\_exp| + |d\_nexp|}, so the assumption biases the estimate by
#' \eqn{2|d\_nexp|}.
#' "Biased by 2|d_nexp|" is not the same as "halved": what the function returns is
#' \eqn{\frac{|d\_exp| - |d\_nexp|}{|d\_exp| + |d\_nexp|}}{(|d_exp| - |d_nexp|)/(|d_exp| + |d_nexp|)}
#' of the true contrast. It therefore collapses to exactly ZERO when the two arms' magnitudes are
#' equal, is halved only in the particular case \eqn{|d\_exp| = 3|d\_nexp|}, and reverses sign when
#' the control arm moved more than the treatment arm. Measured on es_from_paired_f() at
#' n_exp = n_nexp = 50 and r_pre_post = 0.8: two equal F values of 9 return d = 0 where the true
#' contrast is 0.5367; F_exp = 81 against F_nexp = 9 returns 0.5367 against a true 1.0733 (the exact
#' halving); and F_exp = 9 against F_nexp = 81 returns -0.5367 against a true +1.0733. The
#' assumption cannot be verified from the reported statistic, so no quality flag fires on such a row.
#'
#' Mark the arm that moved down with \code{reverse_paired_f_exp} / \code{reverse_paired_f_nexp},
#' which flip the sign of that arm's recovered t before the contrast is formed.
#' \code{reverse_paired_f} flips the whole contrast instead, so it cannot express two arms moving
#' apart. Both per-arm flags are registered input columns, so
#' they are usable from \code{\link{convert_df}()} as well as on a direct call. When the signed t values are reported, \code{\link{es_from_paired_t}()} needs no such
#' assumption, and the \code{es_from_mean_change_*} family carries the direction in the sign of the
#' mean change.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_paired_t}()} are applied.
#'
#' @note
#' Both arms are assumed to have changed in the same direction, and a paired F value can neither
#' confirm nor refute that. Use \code{\link{es_from_paired_t}()} with signed t values, the
#' \code{es_from_mean_change_*} family, or the per-arm \code{reverse_paired_f_exp} /
#' \code{reverse_paired_f_nexp} flags when the arms diverge.
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
                             reverse_paired_f,
                             reverse_paired_f_exp,
                             reverse_paired_f_nexp) {
  if (missing(reverse_paired_f)) reverse_paired_f <- rep(FALSE, length(paired_f_exp))
  reverse_paired_f[is.na(reverse_paired_f)] <- FALSE
  reverse_paired_f_exp <- .arm_direction_flag(
    if (missing(reverse_paired_f_exp)) NULL else reverse_paired_f_exp,
    length(paired_f_exp), "reverse_paired_f_exp"
  )
  reverse_paired_f_nexp <- .arm_direction_flag(
    if (missing(reverse_paired_f_nexp)) NULL else reverse_paired_f_nexp,
    length(paired_f_exp), "reverse_paired_f_nexp"
  )
  if (missing(r_pre_post_nexp)) r_pre_post_nexp <- rep(0.8, length(paired_f_exp))
  r_pre_post_nexp[is.na(r_pre_post_nexp)] <- 0.8
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(paired_f_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8


  paired_t_exp <- sqrt(paired_f_exp)

  paired_t_nexp <- sqrt(paired_f_nexp)

  # An F is the square of the t, so both square roots above are |t|. The per-arm flags
  # are the only way to say that an arm moved down; without them the contrast
  # d_exp - d_nexp assumes the two arms moved the same way (see the details section).
  paired_t_exp <- ifelse(reverse_paired_f_exp, -paired_t_exp, paired_t_exp)
  paired_t_nexp <- ifelse(reverse_paired_f_nexp, -paired_t_nexp, paired_t_nexp)

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
#' @param reverse_paired_f_pval_exp a logical value indicating that the experimental/exposed group changed in the
#'   NEGATIVE direction. A two-sided p-value is unsigned, so this is the only way to state that arm's own
#'   direction (see details).
#' @param reverse_paired_f_pval_nexp a logical value indicating that the non-experimental/non-exposed group
#'   changed in the NEGATIVE direction (see details).
#'
#' @details
#' This function converts the p-values of two paired F-test obtained from two independent groups value into a Cohen's d (D) and Hedges' g (G) (table 12.2 in Cooper).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate the Cohen's d,** the p-values are first converted into paired t-test values:
#' \deqn{paired\_t\_exp = qt(1 - \frac{paired\_f\_pval\_exp}{2}, df = n\_exp - 1)}
#' \deqn{paired\_t\_nexp = qt(1 - \frac{paired\_f\_pval\_nexp}{2}, df = n\_nexp - 1)}
#'
#' which are then converted into a Cohen's d (Cooper et al., 2019):
#' \deqn{cohen\_d\_exp = paired\_t\_exp * \sqrt{\frac{2 * (1 - r\_pre\_post\_exp)}{n\_exp}}}
#' \deqn{cohen\_d\_nexp = paired\_t\_nexp * \sqrt{\frac{2 * (1 - r\_pre\_post\_nexp)}{n\_nexp}}}
#'
#' **A two-sided p-value carries no sign.** The quantile above returns \eqn{|t|}, so only the
#' MAGNITUDE of each arm's change is recoverable and, unless told otherwise, this function assumes
#' that both arms moved in the SAME direction. The between-group estimate is formed as
#' \eqn{cohen\_d = d\_exp - d\_nexp}; when the control arm actually moved the other way -- it
#' deteriorated while the treatment arm improved, an ordinary trial result -- the correct contrast
#' is \eqn{|d\_exp| + |d\_nexp|}, so the assumption biases the estimate by \eqn{2|d\_nexp|}.
#' "Biased by 2|d_nexp|" is not the same as "halved": what the function returns is
#' \eqn{\frac{|d\_exp| - |d\_nexp|}{|d\_exp| + |d\_nexp|}}{(|d_exp| - |d_nexp|)/(|d_exp| + |d_nexp|)}
#' of the true contrast. It therefore collapses to exactly ZERO when the two arms' magnitudes are
#' equal, is halved only in the particular case \eqn{|d\_exp| = 3|d\_nexp|}, and reverses sign when
#' the control arm moved more than the treatment arm. Measured on es_from_paired_f() at
#' n_exp = n_nexp = 50 and r_pre_post = 0.8: two equal F values of 9 return d = 0 where the true
#' contrast is 0.5367; F_exp = 81 against F_nexp = 9 returns 0.5367 against a true 1.0733 (the exact
#' halving); and F_exp = 9 against F_nexp = 81 returns -0.5367 against a true +1.0733. The
#' assumption cannot be verified from the reported statistic, so no quality flag fires on such a row.
#'
#' Mark the arm that moved down with \code{reverse_paired_f_pval_exp} /
#' \code{reverse_paired_f_pval_nexp}, which flip the sign of that arm's recovered t before the
#' contrast is formed. \code{reverse_paired_f_pval} flips the whole contrast instead, so it cannot
#' express two arms moving apart. Both per-arm flags are registered input columns, so
#' they are usable from \code{\link{convert_df}()} as well as on a direct call. When the signed t values are reported,
#' \code{\link{es_from_paired_t}()} needs no such assumption, and the
#' \code{\link{es_from_mean_change_pval}()} family carries the direction in the sign of the mean
#' change.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_paired_t}()} are applied.
#'
#' @note
#' Both arms are assumed to have changed in the same direction, and a two-sided p-value can neither
#' confirm nor refute that. Use \code{\link{es_from_paired_t}()} with signed t values, the
#' \code{es_from_mean_change_*} family, or the per-arm \code{reverse_paired_f_pval_exp} /
#' \code{reverse_paired_f_pval_nexp} flags when the arms diverge.
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
                                  reverse_paired_f_pval,
                                  reverse_paired_f_pval_exp,
                                  reverse_paired_f_pval_nexp) {
  if (missing(reverse_paired_f_pval)) reverse_paired_f_pval <- rep(FALSE, length(paired_f_pval_exp))
  reverse_paired_f_pval[is.na(reverse_paired_f_pval)] <- FALSE
  # A two-sided p-value is unsigned, so the per-arm flags are the only way to say that an
  # arm moved down. They are handed straight to es_from_paired_t_pval(), which flips that
  # arm's recovered t before the contrast is formed (see the details section).
  reverse_paired_f_pval_exp <- .arm_direction_flag(
    if (missing(reverse_paired_f_pval_exp)) NULL else reverse_paired_f_pval_exp,
    length(paired_f_pval_exp), "reverse_paired_f_pval_exp"
  )
  reverse_paired_f_pval_nexp <- .arm_direction_flag(
    if (missing(reverse_paired_f_pval_nexp)) NULL else reverse_paired_f_pval_nexp,
    length(paired_f_pval_exp), "reverse_paired_f_pval_nexp"
  )
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
    reverse_paired_t_pval = reverse_paired_f_pval,
    reverse_paired_t_pval_exp = reverse_paired_f_pval_exp,
    reverse_paired_t_pval_nexp = reverse_paired_f_pval_nexp
  )

  es$info_used <- "paired_f_pval"
  return(es)
}
