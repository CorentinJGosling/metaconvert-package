#' Convert pre-post means of a single group into standard effect size measures
#'
#' @param mean_pre_exp mean of the group at baseline
#' @param mean_exp mean of the group at follow up
#' @param mean_pre_sd_exp standard deviation of the group at baseline
#' @param mean_sd_exp standard deviation of the group at follow up
#' @param n_exp number of participants in the group
#' @param r_pre_post_exp pre-post correlation within the group
#' @param pre_post_to_smd formula used to convert the pre and post means/SD into a SMD (see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_means_pre_post a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts pre-post means of a single group into within-group Cohen's d (dw), Hedges' g (gw), and mean difference (mdw).
#' These within-group effect sizes quantify the standardized change from baseline to follow-up within one group.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from Cohen's d.
#'
#' Four formulas can be used to convert the pre/post means into a SMD (\code{pre_post_to_smd} argument):
#'
#' 1. \code{"bonett"}: baseline SD standardizer (Bonett, 2008; equivalent to the SMCRH measure in metafor)
#' \deqn{d\_w = \frac{mean\_post\_exp - mean\_pre\_exp}{mean\_pre\_sd\_exp}}
#' \deqn{var(g\_w) = \frac{sd\_change^2}{sd\_pre^2(n-1)} + \frac{g\_w^2}{2(n-1)}}
#' where \deqn{sd\_change = \sqrt{sd\_pre^2 + sd\_post^2 - 2r\,sd\_pre\,sd\_post}}
#'
#' 2. \code{"cooper"} (alias: "morris_drm"): raw score standardizer (Cooper 2019, Morris & DeShon 2002)
#' \deqn{d\_rm = \frac{mean\_post\_exp - mean\_pre\_exp}{sd\_change} * \sqrt{2 * (1 - r\_pre\_post\_exp)}}
#' \deqn{var(d\_rm) = \frac{2 * (1 - r\_pre\_post\_exp)}{n\_exp} + \frac{d\_rm^2}{2 * n\_exp}}
#'
#' 3. \code{"morris_dz"}: change score standardizer (Morris & DeShon, 2002; equivalent to the SMCC measure in metafor)
#' \deqn{d\_z = \frac{mean\_post\_exp - mean\_pre\_exp}{sd\_change}}
#' \deqn{var(g\_z) = \frac{1}{n} + \frac{g\_z^2}{2n}}
#'
#' 4. \code{"morris_dav"}: average SD standardizer (Bonett, 2008, eq. 10; equivalent to the SMCRPH measure in metafor)
#' \deqn{sd\_av = \sqrt{(sd\_pre^2 + sd\_post^2)/2}}
#' \deqn{d\_av = \frac{mean\_post\_exp - mean\_pre\_exp}{sd\_av}}
#' \deqn{g\_av = d\_av \times J(mi), \quad mi = \frac{2(n\_exp - 1)}{1 + r\_pre\_post\_exp^2}}
#' \deqn{var(g\_av) = \frac{sd\_change^2}{sd\_av^2 (n - 1)} + \frac{g\_av^2 (sd\_pre^4 + sd\_post^4 + 2 r^2 sd\_pre^2 sd\_post^2)}{8\, sd\_av^4 (n - 1)}}
#'
#' This heteroscedasticity-robust variance is used in preference to the homoscedastic
#' form \eqn{2(1 - r)/n + g\_av^2(1 + r^2)/(4n)} (metafor's SMCRP), which assumes
#' \eqn{sd\_pre = sd\_post} and understates \eqn{var(g\_av)} when the two differ.
#'
#' The within-group Hedges' g is obtained by applying a bias correction factor J to d:
#' \deqn{g\_w = d\_w * J(n\_exp-1)}
#'
#' The within-group mean difference is simply:
#' \deqn{md\_w = mean\_post\_exp - mean\_pre\_exp}
#' \deqn{var(md\_w) = \frac{sd\_change^2}{n\_exp}}
#'
#' To estimate other effect size measures,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied, treating the single group as a matched design.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MDw + Dw + Gw\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 15. Paired: pre-post means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. Psychological Methods, 13(2), 99-109.
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386. https://doi.org/10.1177/1094428106291059
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125. https://doi.org/10.1037/1082-989X.7.1.105
#'
#' @export es_from_means_sd_pre_post_single_group
#'
#' @md
#'
#' @examples
#' es_from_means_sd_pre_post_single_group(
#'   n_exp = 36,
#'   mean_pre_exp = 98, mean_exp = 102,
#'   mean_pre_sd_exp = 16, mean_sd_exp = 17,
#'   r_pre_post_exp = 0.8
#' )
es_from_means_sd_pre_post_single_group <- function(mean_pre_exp, mean_exp,
                                                    mean_pre_sd_exp, mean_sd_exp,
                                                    n_exp, r_pre_post_exp = 0.8,
                                                    pre_post_to_smd = "bonett",
                                                    smd_to_cor = "viechtbauer",
                                                    reverse_means_pre_post) {
  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("bonett", "morris_drm", "morris_dz", "morris_dav"),
    context = "pre_post_means",
    func_name = "es_from_means_sd_pre_post_single_group"
  )

  if (missing(reverse_means_pre_post)) reverse_means_pre_post <- rep(FALSE, length(mean_pre_exp))
  reverse_means_pre_post[is.na(reverse_means_pre_post)] <- FALSE
  if (length(reverse_means_pre_post) == 1) reverse_means_pre_post <- c(rep(reverse_means_pre_post, length(mean_pre_exp)))
  if (length(reverse_means_pre_post) != length(mean_pre_exp)) stop("The length of the 'reverse_means_pre_post' argument is incorrectly specified.")

  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_pre_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  g <- g_se <- g_var <- g_ci_lo <- g_ci_up <-
    d <- d_se <- d_var <- d_ci_lo <- d_ci_up <-
    mdw <- mdw_se <- mdw_ci_lo <- mdw_ci_up <- rep(NA, length(mean_pre_exp))

  nn_miss <- which(!is.na(mean_pre_exp) & !is.na(mean_pre_sd_exp) &
    !is.na(mean_exp) & !is.na(mean_sd_exp) &
    !is.na(r_pre_post_exp) & !is.na(n_exp))

  dat_smd_pre_post <- data.frame(
    mean_pre = mean_pre_exp, mean_pre_sd = mean_pre_sd_exp,
    mean_post = mean_exp, mean_post_sd = mean_sd_exp,
    n = n_exp,
    r_pre_post = r_pre_post_exp,
    pre_post_to_smd = pre_post_to_smd
  )

  if (length(nn_miss) != 0) {
    smd_pp <- t(mapply(.single_group_pre_post_to_smd,
      mean_pre = dat_smd_pre_post$mean_pre[nn_miss],
      mean_post = dat_smd_pre_post$mean_post[nn_miss],
      mean_pre_sd = dat_smd_pre_post$mean_pre_sd[nn_miss],
      mean_post_sd = dat_smd_pre_post$mean_post_sd[nn_miss],
      n = dat_smd_pre_post$n[nn_miss],
      r_pre_post = dat_smd_pre_post$r_pre_post[nn_miss],
      pre_post_to_smd = dat_smd_pre_post$pre_post_to_smd[nn_miss]
    ))

    d[nn_miss] <- smd_pp[, 1]
    d_var[nn_miss] <- smd_pp[, 2]
    d_se[nn_miss] <- sqrt(smd_pp[, 2])
    d_ci_lo[nn_miss] <- smd_pp[, 3]
    d_ci_up[nn_miss] <- smd_pp[, 4]
    g[nn_miss] <- smd_pp[, 5]
    g_var[nn_miss] <- smd_pp[, 6]
    g_se[nn_miss] <- sqrt(smd_pp[, 6])
    g_ci_lo[nn_miss] <- smd_pp[, 7]
    g_ci_up[nn_miss] <- smd_pp[, 8]

    mdw[nn_miss] <- mean_exp[nn_miss] - mean_pre_exp[nn_miss]
    mdw_var <- (mean_pre_sd_exp[nn_miss]^2 + mean_sd_exp[nn_miss]^2 -
                2 * r_pre_post_exp[nn_miss] * mean_pre_sd_exp[nn_miss] * mean_sd_exp[nn_miss]) / n_exp[nn_miss]
    mdw_se[nn_miss] <- sqrt(mdw_var)
    mdw_ci_lo[nn_miss] <- mdw[nn_miss] - qt(0.975, n_exp[nn_miss] - 1) * mdw_se[nn_miss]
    mdw_ci_up[nn_miss] <- mdw[nn_miss] + qt(0.975, n_exp[nn_miss] - 1) * mdw_se[nn_miss]
  }

  es <- .es_from_d(
    d = d, d_se = d_se, n_exp = n_exp, n_nexp = n_exp,
    smd_to_cor = smd_to_cor, reverse = reverse_means_pre_post
  )

  # replace d/g from .es_from_d by the pre/post values
  es$d <- ifelse(reverse_means_pre_post, -d, d)
  es$d_se <- d_se
  es$d_ci_lo <- ifelse(reverse_means_pre_post, -d_ci_up, d_ci_lo)
  es$d_ci_up <- ifelse(reverse_means_pre_post, -d_ci_lo, d_ci_up)

  es$g <- ifelse(reverse_means_pre_post, -g, g)
  es$g_se <- g_se
  es$g_ci_lo <- ifelse(reverse_means_pre_post, -g_ci_up, g_ci_lo)
  es$g_ci_up <- ifelse(reverse_means_pre_post, -g_ci_lo, g_ci_up)

  es$mdw <- ifelse(reverse_means_pre_post, -mdw, mdw)
  es$mdw_se <- mdw_se
  es$mdw_ci_lo <- ifelse(reverse_means_pre_post, -mdw_ci_up, mdw_ci_lo)
  es$mdw_ci_up <- ifelse(reverse_means_pre_post, -mdw_ci_lo, mdw_ci_up)

  es$info_used <- "means_sd_pre_post_single_group"

  return(es)
}

#' Convert pre-post means and standard errors of a single group into standard effect size measures
#'
#' @param mean_pre_exp mean of the group at baseline
#' @param mean_exp mean of the group at follow up
#' @param mean_pre_se_exp standard error of the group at baseline
#' @param mean_se_exp standard error of the group at follow up
#' @param n_exp number of participants in the group
#' @param r_pre_post_exp pre-post correlation within the group
#' @param pre_post_to_smd formula used to convert the pre and post means/SD into a SMD (see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_means_pre_post a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the pre/post standard errors of a single group into standard deviations (Section 6.5.2.2 in the Cochrane Handbook).
#' \deqn{mean\_pre\_sd\_exp = mean\_pre\_se\_exp * \sqrt{n\_exp}}
#' \deqn{mean\_post\_sd\_exp = mean\_post\_se\_exp * \sqrt{n\_exp}}
#'
#' Then, calculations of the \code{\link{es_from_means_sd_pre_post_single_group}()} are applied.
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. Psychological Methods, 13(2), 99-109.
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386.
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MDw + Dw + Gw\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 15. Paired: pre-post means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_means_se_pre_post_single_group
#'
#' @md
#'
#' @examples
#' es_from_means_se_pre_post_single_group(
#'   n_exp = 36,
#'   mean_pre_exp = 98, mean_exp = 102,
#'   mean_pre_se_exp = 2.67, mean_se_exp = 2.83,
#'   r_pre_post_exp = 0.8
#' )
es_from_means_se_pre_post_single_group <- function(mean_pre_exp, mean_exp,
                                                    mean_pre_se_exp, mean_se_exp,
                                                    n_exp, r_pre_post_exp = 0.8,
                                                    pre_post_to_smd = "bonett",
                                                    smd_to_cor = "viechtbauer",
                                                    reverse_means_pre_post) {
  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("bonett", "morris_drm", "morris_dz", "morris_dav"),
    context = "pre_post_means",
    func_name = "es_from_means_se_pre_post_single_group"
  )

  if (missing(reverse_means_pre_post)) reverse_means_pre_post <- rep(FALSE, length(mean_pre_exp))
  reverse_means_pre_post[is.na(reverse_means_pre_post)] <- FALSE
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_pre_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  sd_pre <- mean_pre_se_exp * sqrt(n_exp)
  sd_post <- mean_se_exp * sqrt(n_exp)

  es <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_sd_exp = sd_pre, mean_sd_exp = sd_post,
    n_exp = n_exp,
    r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_means_pre_post = reverse_means_pre_post
  )

  es$info_used <- "means_se_pre_post_single_group"

  return(es)
}

#' Convert pre-post means and 95% CI of a single group into standard effect size measures
#'
#' @param mean_pre_exp mean of the group at baseline
#' @param mean_exp mean of the group at follow up
#' @param mean_pre_ci_lo_exp lower bound of the 95% CI of the mean at baseline
#' @param mean_pre_ci_up_exp upper bound of the 95% CI of the mean at baseline
#' @param mean_ci_lo_exp lower bound of the 95% CI of the mean at follow up
#' @param mean_ci_up_exp upper bound of the 95% CI of the mean at follow up
#' @param n_exp number of participants in the group
#' @param r_pre_post_exp pre-post correlation within the group
#' @param pre_post_to_smd formula used to convert the pre and post means/SD into a SMD (see details).
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param max_asymmetry percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param reverse_means_pre_post a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the bounds of the 95% CI of the pre/post means of a single group into standard errors (Section 6.3.1 in the Cochrane Handbook).
#' \deqn{mean\_pre\_se\_exp = \frac{mean\_pre\_ci\_up\_exp - mean\_pre\_ci\_lo\_exp}{2 * qt{(0.975, df = n\_exp - 1)}}}
#' \deqn{mean\_post\_se\_exp = \frac{mean\_post\_ci\_up\_exp - mean\_post\_ci\_lo\_exp}{2 * qt{(0.975, df = n\_exp - 1)}}}
#'
#' Then, calculations of the \code{\link{es_from_means_se_pre_post_single_group}()} are applied.
#'
#' @references
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022.
#'
#' Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. Psychological Methods, 13(2), 99-109.
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386.
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_means_ci_pre_post_single_group
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MDw + Dw + Gw\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 15. Paired: pre-post means and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_means_ci_pre_post_single_group(
#'   n_exp = 36,
#'   mean_pre_exp = 98,
#'   mean_pre_ci_lo_exp = 88, mean_pre_ci_up_exp = 108,
#'   mean_exp = 102,
#'   mean_ci_lo_exp = 92, mean_ci_up_exp = 112,
#'   r_pre_post_exp = 0.8
#' )
es_from_means_ci_pre_post_single_group <- function(mean_pre_exp, mean_exp,
                                                    mean_pre_ci_lo_exp, mean_pre_ci_up_exp,
                                                    mean_ci_lo_exp, mean_ci_up_exp,
                                                    n_exp, r_pre_post_exp = 0.8,
                                                    pre_post_to_smd = "bonett",
                                                    smd_to_cor = "viechtbauer",
                                                    max_asymmetry = 10,
                                                    reverse_means_pre_post) {

  if (missing(reverse_means_pre_post)) reverse_means_pre_post <- rep(FALSE, length(mean_pre_exp))
  reverse_means_pre_post[is.na(reverse_means_pre_post)] <- FALSE
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_pre_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  .check_ci_asymmetry <- function(val, lo, up, label) {
    half_up <- abs(up - val)
    half_lo <- abs(val - lo)
    denom <- pmax(half_up, half_lo)
    asym_pct <- ifelse(denom == 0, 0, abs(half_up - half_lo) / denom * 100)
    bad <- which(!is.na(asym_pct) & asym_pct > max_asymmetry)
    if (length(bad) > 0) {
      warning(sprintf("Asymmetric CI detected for '%s' in row(s) %s (asymmetry > %s%%)",
                      label, paste(bad, collapse = ", "), max_asymmetry), call. = FALSE)
    }
  }
  .check_ci_asymmetry(mean_pre_exp, mean_pre_ci_lo_exp, mean_pre_ci_up_exp, "mean_pre_exp")
  .check_ci_asymmetry(mean_exp, mean_ci_lo_exp, mean_ci_up_exp, "mean_exp")

  df <- n_exp - 1

  se_pre <- (mean_pre_ci_up_exp - mean_pre_ci_lo_exp) / (2 * qt(0.975, df))
  se_post <- (mean_ci_up_exp - mean_ci_lo_exp) / (2 * qt(0.975, df))

  es <- es_from_means_se_pre_post_single_group(
    mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
    mean_pre_se_exp = se_pre, mean_se_exp = se_post,
    n_exp = n_exp,
    r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_means_pre_post = reverse_means_pre_post
  )

  es$info_used <- "means_ci_pre_post_single_group"

  return(es)
}

#' Convert mean change and standard deviation of a single group into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the group.
#' @param mean_change_sd_exp standard deviation of the mean change for participants in the group.
#' @param n_exp number of participants in the group.
#' @param r_pre_post_exp pre-post correlation within the group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the mean change into a SMD (\code{"cooper"} by default; \code{"morris_dz"} also accepted, see details).
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function computes within-group Cohen's d (dw), Hedges' g (gw), and mean difference (mdw)
#' from the mean change (MC) and standard deviation of a single group.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from Cohen's d.
#'
#' The effect size (d_rm) depends linearly on the pre-post correlation:
#' \deqn{d\_rm = (mean\_change\_exp / sd\_change) \times \sqrt{2 \times (1 - r\_pre\_post\_exp)}}
#' If \code{r_pre_post_exp} is not provided, a value of 0.8 is assumed. A misspecified correlation
#' biases the point estimate, not only its precision; when the correlation is not reported,
#' a sensitivity analysis using several plausible values (e.g., 0.3, 0.5, 0.7) is advisable.
#'
#' This function simply internally calls the \code{\link{es_from_means_sd_pre_post_single_group}} function but setting:
#' \deqn{mean\_pre\_exp = 0}
#' \deqn{mean\_pre\_sd\_exp = 0}
#' \deqn{mean\_exp = mean\_change\_exp}
#' \deqn{mean\_sd\_exp = mean\_change\_sd\_exp}
#'
#' To know more about the calculations, see \code{\link{es_from_means_sd_pre_post_single_group}} function.
#'
#' @references
#' Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. Psychological Methods, 13(2), 99-109.
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386.
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MDw + Dw + Gw\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 14. Paired: mean change, and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_mean_change_sd_single_group
#'
#' @md
#'
#' @examples
#' es_from_mean_change_sd_single_group(
#'   n_exp = 36,
#'   mean_change_exp = 8.4, mean_change_sd_exp = 9.13,
#'   r_pre_post_exp = 0.2
#' )
es_from_mean_change_sd_single_group <- function(mean_change_exp, mean_change_sd_exp,
                                                 n_exp, r_pre_post_exp = 0.8,
                                                 smd_to_cor = "viechtbauer",
                                                 pre_post_to_smd = "cooper",
                                                 reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "mean_change",
    func_name = "es_from_mean_change_sd_single_group"
  )

  es <- es_from_means_sd_pre_post_single_group(
    mean_pre_exp = rep(0, length(mean_change_exp)),
    mean_exp = mean_change_exp,
    mean_pre_sd_exp = rep(0, length(mean_change_exp)),
    mean_sd_exp = mean_change_sd_exp,
    n_exp = n_exp,
    r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_sd_single_group"

  return(es)
}

#' Convert mean change and standard error of a single group into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the group.
#' @param mean_change_se_exp standard error of the mean change for participants in the group.
#' @param n_exp number of participants in the group.
#' @param r_pre_post_exp pre-post correlation within the group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the mean change into a SMD (\code{"cooper"} by default; \code{"morris_dz"} also accepted, see details).
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean change and standard error of a single group
#' into within-group Cohen's d, Hedges' g, and mean difference. The effect sizes are then converted to other measures.
#'
#' The effect size (d_rm) depends linearly on the pre-post correlation (0.8 assumed when
#' \code{r_pre_post_exp} is missing); see \code{\link{es_from_mean_change_sd_single_group}}.
#'
#' This function simply internally calls the \code{\link{es_from_means_se_pre_post_single_group}} function but setting:
#' \deqn{mean\_pre\_exp = 0}
#' \deqn{mean\_pre\_se\_exp = 0}
#' \deqn{mean\_exp = mean\_change\_exp}
#' \deqn{mean\_se\_exp = mean\_change\_se\_exp}
#'
#' To know more about the calculations, see \code{\link{es_from_means_se_pre_post_single_group}} function.
#'
#' @references
#' Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. Psychological Methods, 13(2), 99-109.
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386.
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MDw + Dw + Gw\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 14. Paired: mean change, and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_mean_change_se_single_group
#'
#' @md
#'
#' @examples
#' es_from_mean_change_se_single_group(
#'   n_exp = 36,
#'   mean_change_exp = 8.4, mean_change_se_exp = 1.52,
#'   r_pre_post_exp = 0.2
#' )
es_from_mean_change_se_single_group <- function(mean_change_exp, mean_change_se_exp,
                                                 n_exp, r_pre_post_exp = 0.8,
                                                 smd_to_cor = "viechtbauer",
                                                 pre_post_to_smd = "cooper",
                                                 reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "mean_change",
    func_name = "es_from_mean_change_se_single_group"
  )

  es <- es_from_means_se_pre_post_single_group(
    mean_pre_exp = rep(0, length(mean_change_exp)),
    mean_exp = mean_change_exp,
    mean_pre_se_exp = rep(0, length(mean_change_exp)),
    mean_se_exp = mean_change_se_exp,
    n_exp = n_exp,
    r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_se_single_group"

  return(es)
}

#' Convert mean change and 95% CI of a single group into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the group.
#' @param mean_change_ci_lo_exp lower bound of the 95% CI around the mean change.
#' @param mean_change_ci_up_exp upper bound of the 95% CI around the mean change.
#' @param n_exp number of participants in the group.
#' @param r_pre_post_exp pre-post correlation within the group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the mean change into a SMD (\code{"cooper"} by default; \code{"morris_dz"} also accepted, see details).
#' @param max_asymmetry percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean change and 95% CI of a single group
#' into within-group Cohen's d, Hedges' g, and mean difference. The effect sizes are then converted to other measures.
#'
#' The effect size (d_rm) depends linearly on the pre-post correlation (0.8 assumed when
#' \code{r_pre_post_exp} is missing); see \code{\link{es_from_mean_change_sd_single_group}}.
#'
#' This function simply internally calls the \code{\link{es_from_means_ci_pre_post_single_group}} function but setting:
#' \deqn{mean\_pre\_exp = 0}
#' \deqn{mean\_pre\_ci\_lo\_exp = 0}
#' \deqn{mean\_pre\_ci\_up\_exp = 0}
#' \deqn{mean\_exp = mean\_change\_exp}
#' \deqn{mean\_ci\_lo\_exp = mean\_change\_ci\_lo\_exp}
#' \deqn{mean\_ci\_up\_exp = mean\_change\_ci\_up\_exp}
#'
#' To know more about the calculations, see \code{\link{es_from_means_sd_pre_post_single_group}} function.
#'
#' @references
#' Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. Psychological Methods, 13(2), 99-109.
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386.
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_mean_change_ci_single_group
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MDw + Dw + Gw\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 14. Paired: mean change, and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_mean_change_ci_single_group(
#'   n_exp = 36,
#'   mean_change_exp = 8.4,
#'   mean_change_ci_lo_exp = 6.4, mean_change_ci_up_exp = 10.4,
#'   r_pre_post_exp = 0.2
#' )
es_from_mean_change_ci_single_group <- function(mean_change_exp,
                                                 mean_change_ci_lo_exp, mean_change_ci_up_exp,
                                                 n_exp, r_pre_post_exp = 0.8,
                                                 max_asymmetry = 10,
                                                 smd_to_cor = "viechtbauer",
                                                 pre_post_to_smd = "cooper",
                                                 reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "mean_change",
    func_name = "es_from_mean_change_ci_single_group"
  )

  es <- es_from_means_ci_pre_post_single_group(
    mean_pre_exp = rep(0, length(mean_change_exp)),
    mean_pre_ci_lo_exp = rep(0, length(mean_change_exp)),
    mean_pre_ci_up_exp = rep(0, length(mean_change_exp)),
    mean_exp = mean_change_exp,
    mean_ci_lo_exp = mean_change_ci_lo_exp,
    mean_ci_up_exp = mean_change_ci_up_exp,
    n_exp = n_exp,
    r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    max_asymmetry = max_asymmetry,
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_ci_single_group"

  return(es)
}

#' Convert mean change and p-value of a single group into standard effect size measures
#'
#' @param mean_change_exp mean change of participants in the group.
#' @param mean_change_pval_exp p-value of the mean change for participants in the group.
#' @param n_exp number of participants in the group.
#' @param r_pre_post_exp pre-post correlation within the group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the mean change into a SMD (\code{"cooper"} by default; \code{"morris_dz"} also accepted, see details).
#' @param reverse_mean_change a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the mean change and associated p-value of a single group
#' into within-group Cohen's d, Hedges' g, and mean difference. The effect sizes are then converted to other measures.
#'
#' The effect size (d_rm) depends linearly on the pre-post correlation (0.8 assumed when
#' \code{r_pre_post_exp} is missing); see \code{\link{es_from_mean_change_sd_single_group}}.
#'
#' To start, this function estimates the mean change standard error from the p-value:
#' \deqn{t <- qt(p = mean\_change\_pval\_exp / 2, df = n\_exp - 1, lower.tail = FALSE)}
#' \deqn{mean\_change\_se\_exp <- |\frac{mean\_change\_exp}{t}|}
#'
#' Then, this function simply internally calls the \code{\link{es_from_means_se_pre_post_single_group}} function but setting:
#' \deqn{mean\_pre\_exp = 0}
#' \deqn{mean\_pre\_se\_exp = 0}
#' \deqn{mean\_exp = mean\_change\_exp}
#' \deqn{mean\_se\_exp = mean\_change\_se\_exp}
#'
#' To know more about other calculations, see \code{\link{es_from_means_sd_pre_post_single_group}} function.
#'
#' @references
#' Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. Psychological Methods, 13(2), 99-109.
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control group designs. Organizational Research Methods, 11(2), 364-386.
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_mean_change_pval_single_group
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MDw + Dw + Gw\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 14. Paired: mean change, and dispersion'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_mean_change_pval_single_group(
#'   n_exp = 36,
#'   mean_change_exp = 8.4, mean_change_pval_exp = 0.13,
#'   r_pre_post_exp = 0.8
#' )
es_from_mean_change_pval_single_group <- function(mean_change_exp, mean_change_pval_exp,
                                                   n_exp, r_pre_post_exp = 0.8,
                                                   smd_to_cor = "viechtbauer",
                                                   pre_post_to_smd = "cooper",
                                                   reverse_mean_change) {
  if (missing(reverse_mean_change)) reverse_mean_change <- rep(FALSE, length(mean_change_exp))
  reverse_mean_change[is.na(reverse_mean_change)] <- FALSE
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(mean_change_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "mean_change",
    func_name = "es_from_mean_change_pval_single_group"
  )

  t_stat <- qt(p = mean_change_pval_exp / 2, df = n_exp - 1, lower.tail = FALSE)
  mean_change_se <- abs(mean_change_exp / t_stat)

  es <- es_from_means_se_pre_post_single_group(
    mean_pre_exp = rep(0, length(mean_change_exp)),
    mean_exp = mean_change_exp,
    mean_pre_se_exp = rep(0, length(mean_change_exp)),
    mean_se_exp = mean_change_se,
    n_exp = n_exp,
    r_pre_post_exp = r_pre_post_exp,
    smd_to_cor = smd_to_cor,
    pre_post_to_smd = pre_post_to_smd,
    reverse_means_pre_post = reverse_mean_change
  )

  es$info_used <- "mean_change_pval_single_group"

  return(es)
}

#' Convert paired t-test statistic from a single group into standard effect size measures
#'
#' @param paired_t_exp paired t-test statistic (pre-post comparison) for the group.
#' @param n_exp number of participants in the group.
#' @param r_pre_post_exp pre-post correlation within the group
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param pre_post_to_smd formula used to convert the paired t statistic into a SMD (\code{"cooper"} by default; \code{"morris_dz"} also accepted, see details).
#' @param reverse_paired_t a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts a paired t-test statistic from a single group
#' into within-group Cohen's d, Hedges' g, and mean difference using the raw score standardizer (morris_drm/cooper).
#'
#' The paired t-test statistic is related to the mean difference and standard error of difference:
#' \deqn{t\_paired = \frac{mean\_diff}{SE\_diff}}
#'
#' where:
#' \deqn{SE\_diff = \frac{SD\_diff}{\sqrt{n}}}
#' \deqn{SD\_diff = \sqrt{SD\_pre^2 + SD\_post^2 - 2 \times r \times SD\_pre \times SD\_post}}
#'
#' This function calculates the raw score standardized mean difference (d_rm, Morris & DeShon 2002):
#' \deqn{d\_rm = t\_paired \times \sqrt{\frac{2 \times (1 - r)}{n}}}
#'
#' The variance is calculated as:
#' \deqn{var(d\_rm) = \frac{2 \times (1 - r)}{n} + \frac{d\_rm^2}{2 \times n}}
#'
#' The result depends on the pre-post correlation (0.8 assumed when \code{r_pre_post_exp}
#' is missing); see \code{\link{es_from_mean_change_sd_single_group}}.
#'
#' The within-group Hedges' g is obtained by applying a bias correction factor J to d.
#'
#' To estimate other effect size measures, calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @references
#' Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. Psychological Methods, 7(1), 105-125.
#'
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MDw + Dw + Gw\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 15. Paired: pre-post paired t-test'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_paired_t_single_group
#'
#' @md
#'
#' @examples
#' es_from_paired_t_single_group(
#'   paired_t_exp = 2.5,
#'   n_exp = 20,
#'   r_pre_post_exp = 0.7
#' )
es_from_paired_t_single_group <- function(paired_t_exp, n_exp, r_pre_post_exp = 0.8,
                                          smd_to_cor = "viechtbauer",
                                          pre_post_to_smd = "cooper",
                                          reverse_paired_t) {
  if (missing(reverse_paired_t)) reverse_paired_t <- rep(FALSE, length(paired_t_exp))
  reverse_paired_t[is.na(reverse_paired_t)] <- FALSE
  if (missing(r_pre_post_exp)) r_pre_post_exp <- rep(0.8, length(paired_t_exp))
  r_pre_post_exp[is.na(r_pre_post_exp)] <- 0.8

  pre_post_to_smd <- .validate_pre_post_to_smd(
    pre_post_to_smd,
    allowed_methods = c("morris_drm", "morris_dz"),
    context = "paired_t",
    func_name = "es_from_paired_t_single_group"
  )

  r_pre_post_exp <- .guard_r_pre_post(r_pre_post_exp)

  J <- .d_j(n_exp - 1)

  # per-row method selection (a mixed vector previously fell wholesale into drm).
  # Recycled to the input length so a scalar method still yields one value per row.
  dz <- rep(pre_post_to_smd == "morris_dz", length.out = length(paired_t_exp))

  # d_z. metafor SMCC convention (variance built from the corrected g), matching
  # .single_group_pre_post_to_smd so the paired-t and mean-change routes return
  # identical SEs on equivalent inputs
  d_dz <- paired_t_exp / sqrt(n_exp)
  d_var_dz <- (1 / n_exp + (J * d_dz)^2 / (2 * n_exp)) / J^2

  # d_rm
  d_drm <- paired_t_exp * sqrt(2 * (1 - r_pre_post_exp) / n_exp)
  d_var_drm <- 2 * (1 - r_pre_post_exp) / n_exp + d_drm^2 / (2 * n_exp)

  d <- ifelse(dz, d_dz, d_drm)
  d_var <- ifelse(dz, d_var_dz, d_var_drm)

  # A single paired observation (n < 2) has no estimable within-subject variance. The
  # morris_drm branch carries no J(n-1) factor, so it stays finite at n = 1 (and Inf at
  # n = 0), leaking a spurious d/g/logOR/r (and an impossible |r| > 1) downstream. NA
  # such rows so this route matches the mean-change route (which NAs via its n >= 2
  # kernel guard); the qt(.975, n-1) CI is already NaN there.
  bad_n <- !(is.finite(n_exp) & n_exp >= 2)
  d[bad_n] <- NA_real_; d_var[bad_n] <- NA_real_

  d_se <- sqrt(d_var)

  d_ci_lo <- d - qt(0.975, n_exp - 1) * d_se
  d_ci_up <- d + qt(0.975, n_exp - 1) * d_se

  g <- d * J
  g_se <- d_se * J
  g_ci_lo <- d_ci_lo * J
  g_ci_up <- d_ci_up * J

  # md not recoverable from the t statistic alone
  md <- rep(NA_real_, length(paired_t_exp))
  md_se <- rep(NA_real_, length(paired_t_exp))
  md_ci_lo <- rep(NA_real_, length(paired_t_exp))
  md_ci_up <- rep(NA_real_, length(paired_t_exp))

  es <- .es_from_d(
    d = d, d_se = d_se, n_exp = n_exp, n_nexp = n_exp,
    smd_to_cor = smd_to_cor, reverse = reverse_paired_t
  )

  es$d <- ifelse(reverse_paired_t, -d, d)
  es$d_se <- d_se
  es$d_ci_lo <- ifelse(reverse_paired_t, -d_ci_up, d_ci_lo)
  es$d_ci_up <- ifelse(reverse_paired_t, -d_ci_lo, d_ci_up)

  es$g <- ifelse(reverse_paired_t, -g, g)
  es$g_se <- g_se
  es$g_ci_lo <- ifelse(reverse_paired_t, -g_ci_up, g_ci_lo)
  es$g_ci_up <- ifelse(reverse_paired_t, -g_ci_lo, g_ci_up)

  es$md <- ifelse(reverse_paired_t, -md, md)
  es$md_se <- md_se
  es$md_ci_lo <- ifelse(reverse_paired_t, -md_ci_up, md_ci_lo)
  es$md_ci_up <- ifelse(reverse_paired_t, -md_ci_lo, md_ci_up)

  es$info_used <- "paired_t_single_group"

  return(es)
}
