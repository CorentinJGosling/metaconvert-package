#' Convert a standardized regression coefficient and the standard deviation of the dependent variable
#' into several effect size measures
#'
#' @param beta_std a standardized regression coefficient value (binary predictor, no other covariables in the model)
#' @param sd_dv standard deviation of the dependent variable
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_beta_std a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts a standardized linear regression coefficient
#' (coming from a model with only one binary predictor), into an
#' unstandardized linear regression coefficient.
#'
#' \deqn{sd\_dummy = \sqrt{\frac{n_exp - (n_exp^2 / (n_exp + n_nexp))}{(n_exp + n_nexp - 1)}}}
#' \deqn{unstd\_beta = beta\_std * \frac{sd\_dv}{sd\_dummy}}
#'
#' Calculations of the \code{\link{es_from_beta_unstd}} functions are then used.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 13. (Un-)Standardized regression coefficient'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' @export es_from_beta_std
#'
#' @md
#'
#' @examples
#' es_from_beta_std(beta_std = 0.35, sd_dv = 0.98, n_exp = 20, n_nexp = 22)
es_from_beta_std <- function(beta_std, sd_dv, n_exp, n_nexp,
                             smd_to_cor = "viechtbauer", reverse_beta_std) {
  if (missing(reverse_beta_std)) reverse_beta_std <- rep(FALSE, length(beta_std))
  reverse_beta_std[is.na(reverse_beta_std)] <- FALSE
  if (length(reverse_beta_std) == 1) reverse_beta_std = c(rep(reverse_beta_std, length(beta_std)))
  if (length(reverse_beta_std) != length(beta_std)) stop("The length of the 'reverse_beta_std' argument is incorrectly specified.")

  beta_std <- ifelse(reverse_beta_std, -beta_std, beta_std)


  sd_dummy <- sqrt((n_exp - (n_exp^2 / (n_exp + n_nexp))) / (n_exp + n_nexp - 1))

  unstd_beta <- beta_std * (sd_dv / sd_dummy)

  es <- es_from_beta_unstd(
    beta_unstd = unstd_beta, sd_dv = sd_dv,
    n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor
  )

  es$info_used <- "beta_std"

  return(es)
}

#' Convert an unstandardized regression coefficient and the standard deviation of the dependent variable
#' into several effect size measures
#'
#' @param beta_unstd an unstandardized regression coefficient value (binary predictor, no other covariables in the model)
#' @param sd_dv standard deviation of the dependent variable
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_beta_unstd a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function estimates a Cohen's d (D) and Hedges' g (G) from an unstandardized linear regression coefficient (coming from a model with only one binary predictor),
#' and the standard deviation of the dependent variable.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **The formula used to obtain the Cohen's d is**:
#' \deqn{N = n\_exp + n\_nexp}
#' \deqn{sd\_pooled = \sqrt{\frac{sd\_dv^2 * (N - 1) - unstd\_beta^2 * \frac{n\_exp * n\_nexp}{N}}{N - 2}}}
#' \deqn{cohen\_d = \frac{unstd\_beta}{sd\_pooled}}
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
#'  \code{required input data} \tab See 'Section 13. (Un-)Standardized regression coefficient'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' @export es_from_beta_unstd
#'
#' @md
#'
#' @examples
#' es_from_beta_unstd(beta_unstd = 0.7, sd_dv = 0.98, n_exp = 20, n_nexp = 22)
es_from_beta_unstd <- function(beta_unstd, sd_dv, n_exp, n_nexp,
                               smd_to_cor = "viechtbauer", reverse_beta_unstd) {
  if (missing(reverse_beta_unstd)) reverse_beta_unstd <- rep(FALSE, length(beta_unstd))
  reverse_beta_unstd[is.na(reverse_beta_unstd)] <- FALSE


  # Within-group variance from the ANOVA SS decomposition
  # SS_total = SS_between + SS_within, with SS_between = beta_unstd^2 * n_exp*n_nexp/N.
  # This is >= 0 for any genuine binary-predictor OLS; it goes negative only when the
  # reported sd_dv is too small to be consistent with beta_unstd (a mathematically
  # impossible input, e.g. |standardized beta| > 1, or the residual SD entered as sd_dv).
  # Such rows are set to NA rather than masked with abs(), which would fabricate a
  # plausible-looking effect size.
  within_var <- ((sd_dv^2 * (n_exp + n_nexp - 1)) -
                 (beta_unstd^2 * ((n_exp * n_nexp) / (n_exp + n_nexp)))) /
                (n_exp + n_nexp - 2)

  n_impossible <- sum(within_var < 0, na.rm = TRUE)
  if (n_impossible > 0) {
    warning(sprintf(
      "es_from_beta_unstd: %d row(s) had a reported sd_dv too small to be consistent with beta_unstd (implied within-group variance < 0). These inputs are mathematically impossible (e.g. |standardized beta| > 1); the effect sizes are set to NA for those rows.",
      n_impossible
    ), call. = FALSE)
  }

  within_var <- ifelse(within_var < 0, NA_real_, within_var)
  sd_pooled <- sqrt(within_var) # NA where the inputs were impossible

  d <- beta_unstd / sd_pooled

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_beta_unstd
  )

  es$info_used <- "beta_unstd"

  return(es)
}


#' Convert a t-statistic from a linear regression model to several effect size measures
#'
#' @param linreg_t a t-statistic from a linear regression model
#' @param n_sample the total number of participants
#' @param n_covariates the number of covariates in the model (excluding the predictor of interest).
#' @param sd_iv the standard deviation of the independent variable (optional, see details)
#' @param unit_increase_iv a value of the independent variable that will be used to estimate the Cohen's d (optional, see details).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Must be either "sd" or "raw_scale"
#' @param n_exp number of the experimental/exposed group (optional)
#' @param n_nexp number of the non-experimental/non-exposed group (optional)
#' @param cor_to_smd formula used to convert a \code{pearson_r} or \code{fisher_z} value into a SMD.
#' @param reverse_linreg_t a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts a t-statistic from a linear regression model into a
#' partial correlation coefficient (Rp) and its Fisher's z transformation (Zp).
#'
#' The partial correlation is obtained as (Aloe & Thompson, 2013, Eq. 2; Gustafson, 1961):
#' \deqn{r_p = \frac{t}{\sqrt{t^2 + df}}}
#' where \eqn{df = n\_sample - n\_covariates - 2} is the residual degrees of freedom
#' of the regression model (i.e. \eqn{n} minus the intercept, the focal predictor,
#' and the \code{n_covariates} covariates).
#'
#' Its sampling variance is estimated as recommended by van Aert & Goos (2023, Eq. 5):
#' \deqn{var(r_p) = \frac{(1 - r_p^2)^2}{df}}
#'
#' The Fisher's z transformation and its variance are:
#' \deqn{z_p = atanh(r_p)}
#' \deqn{var(z_p) = \frac{1}{n - n\_covariates - 3}}
#'
#' Cohen's d, Hedges' g and odds ratio are then converted from the partial correlation
#' using the calculations of the \code{\link{es_from_pearson_r}} function.
#' A partial correlation controls for covariates and thus targets a different estimand
#' than a two-group comparison or a bivariate correlation (Aloe & Thompson, 2013);
#' the converted d/g/OR values should not be pooled with such effect sizes.
#' For meta-analyses of regression results, use \code{measure = "rp"} or
#' \code{measure = "zp"} in \code{\link{convert_df}}.
#'
#' @export es_from_linreg_t
#'
#' @references
#' Aloe, A. M., & Thompson, C. G. (2013). The synthesis of partial effect sizes.
#' \emph{Journal of the Society for Social Work and Research}, 4(4), 390--405.
#'
#' Gustafson, R. L. (1961). Partial correlations in regression computations.
#' \emph{Journal of the American Statistical Association}, 56, 363--367.
#'
#' Mathur, M. B., & VanderWeele, T. J. (2020). A simple, interpretable conversion from
#' Pearson's correlation to Cohen's d for continuous exposures.
#' \emph{Epidemiology}, 31(2), e16--e18.
#'
#' van Aert, R. C. M., & Goos, C. (2023). A critical reflection on computing the sampling
#' variance of the partial correlation coefficient. \emph{Research Synthesis Methods},
#' 14(3), 520--525.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab Rp (partial correlation) + Zp (Fisher's z of partial r)\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + OR (approximate, see details)\cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_linreg_t(linreg_t = 2.5, n_sample = 100, n_covariates = 2)
#'
#' # Reproduce Aloe & Thompson (2013) Table 1, Cole et al. (2004)
#' es_from_linreg_t(linreg_t = 6.19, n_sample = 232, n_covariates = 6)
es_from_linreg_t <- function(linreg_t, n_sample, n_covariates,
                             sd_iv, unit_increase_iv, unit_type = "raw_scale",
                             n_exp, n_nexp, cor_to_smd = "viechtbauer",
                             reverse_linreg_t) {

  if (missing(reverse_linreg_t)) reverse_linreg_t <- rep(FALSE, length(linreg_t))
  reverse_linreg_t[is.na(reverse_linreg_t)] <- FALSE
  if (length(reverse_linreg_t) == 1) reverse_linreg_t <- rep(reverse_linreg_t, length(linreg_t))
  if (length(reverse_linreg_t) != length(linreg_t)) stop("The length of the 'reverse_linreg_t' argument is incorrectly specified.")

  if (missing(n_exp)) n_exp <- rep(NA, length(linreg_t))
  if (missing(n_nexp)) n_nexp <- rep(NA, length(linreg_t))
  if (missing(sd_iv)) sd_iv <- rep(NA, length(linreg_t))
  if (missing(unit_increase_iv)) unit_increase_iv <- rep(NA, length(linreg_t))
  if (missing(unit_type)) unit_type <- rep(NA, length(linreg_t))

  if (!all(cor_to_smd %in% c("cooper", "mathur", "viechtbauer"))) {
    stop(paste0("'",
                unique(cor_to_smd[!cor_to_smd %in% c("cooper", "mathur", "viechtbauer")]),
                "' not in tolerated values for the 'cor_to_smd' argument.
                Possible inputs are: 'cooper', 'mathur', 'viechtbauer'"))
  }

  df <- n_sample - n_covariates - 2

  if (any(df <= 0, na.rm = TRUE)) {
    stop("Degrees of freedom must be positive. Check n_sample and n_covariates.")
  }

  rp <- linreg_t / sqrt(linreg_t^2 + df)
  rp <- ifelse(reverse_linreg_t, -rp, rp)

  rp_se <- sqrt((1 - rp^2)^2 / df)
  rp_ci_lo <- rp - qt(.975, df) * rp_se
  rp_ci_up <- rp + qt(.975, df) * rp_se

  zp <- atanh(rp)
  # target Fisher-z variance of a partial correlation controlling for
  # n_covariates variables is 1 / (n - n_covariates - 3); with the residual
  # df = n - n_covariates - 2 this is 1 / (df - 1).
  zp_se <- sqrt(1 / (df - 1))
  zp_ci_lo <- zp - qnorm(.975) * zp_se
  zp_ci_up <- zp + qnorm(.975) * zp_se

  n_sample_calc <- ifelse(is.na(n_sample), n_exp + n_nexp, n_sample)
  n_exp <- ifelse(is.na(n_exp), n_sample_calc / 2, n_exp)
  n_nexp <- ifelse(is.na(n_nexp), n_sample_calc / 2, n_nexp)

  dat_cor <- data.frame(
    r = rp, r_se = rp_se,
    sd_iv = sd_iv, n_sample = n_sample_calc,
    unit_increase_iv = unit_increase_iv,
    unit_type = unit_type,
    cor_to_smd = cor_to_smd
  )

  nn_miss <- with(dat_cor, which(
    (cor_to_smd == "mathur" & !is.na(r) & !is.na(sd_iv) &
      !is.na(n_sample) & !is.na(unit_increase_iv) & !is.na(unit_type)) |
      (cor_to_smd == "cooper" & !is.na(r) & !is.na(r_se)) |
      (cor_to_smd == "viechtbauer" & !is.na(r) & !is.na(r_se) & !is.na(n_sample))
  ))

  es <- data.frame(
    d = rep(NA, nrow(dat_cor)),
    d_se = rep(NA, nrow(dat_cor))
  )

  if (length(nn_miss) != 0) {
    res_d <- t(mapply(.cor_to_smd,
      r = dat_cor$r[nn_miss],
      r_se = dat_cor$r_se[nn_miss],
      n_sample = dat_cor$n_sample[nn_miss],
      sd_iv = dat_cor$sd_iv[nn_miss],
      unit_increase_iv = dat_cor$unit_increase_iv[nn_miss],
      unit_type = dat_cor$unit_type[nn_miss],
      cor_to_smd = dat_cor$cor_to_smd[nn_miss]
    ))

    es$d[nn_miss] <- unlist(res_d[, 1])
    es$d_se[nn_miss] <- unlist(res_d[, 2])
  }

  es <- .es_from_d(d = es$d, d_se = es$d_se, n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample_calc)

  es$rp <- rp
  es$rp_se <- rp_se
  es$rp_ci_lo <- rp_ci_lo
  es$rp_ci_up <- rp_ci_up

  es$zp <- zp
  es$zp_se <- zp_se
  es$zp_ci_lo <- zp_ci_lo
  es$zp_ci_up <- zp_ci_up

  es$info_used <- "linreg_t"

  return(es)
}


#' Convert an unstandardized regression coefficient and its standard error
#' into several effect size measures
#'
#' @param linreg_b unstandardized regression coefficient from a linear regression model
#' @param linreg_b_se standard error of the regression coefficient
#' @param n_sample the total number of participants
#' @param n_covariates the number of covariates in the model (excluding the predictor of interest).
#' @param sd_iv the standard deviation of the independent variable (optional, see details)
#' @param unit_increase_iv a value of the independent variable that will be used to estimate the Cohen's d (optional, see details).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Must be either "sd" or "raw_scale"
#' @param n_exp number of the experimental/exposed group (optional)
#' @param n_nexp number of the non-experimental/non-exposed group (optional)
#' @param cor_to_smd formula used to convert a \code{pearson_r} or \code{fisher_z} value into a SMD.
#' @param reverse_linreg_b a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function derives a t-statistic from a regression coefficient and its standard error,
#' then converts it into a partial correlation and other effect size measures.
#'
#' The Wald t-statistic is computed as:
#' \deqn{t = \frac{b}{SE(b)}}
#'
#' Once the t-statistic is obtained, all subsequent conversions follow the same formulas as
#' in \code{\link{es_from_linreg_t}}.
#'
#' For binary predictors without covariates, use \code{\link{es_from_beta_unstd}} instead,
#' which converts directly to Cohen's d.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab Rp (partial correlation) + Zp (Fisher's z of partial r)\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + OR (approximate, see \code{\link{es_from_linreg_t}})\cr
#' }
#'
#' @references
#' Aloe, A. M., & Thompson, C. G. (2013). The synthesis of partial effect sizes.
#' \emph{Journal of the Society for Social Work and Research}, 4(4), 390--405.
#'
#' van Aert, R. C. M., & Goos, C. (2023). A critical reflection on computing the sampling
#' variance of the partial correlation coefficient. \emph{Research Synthesis Methods},
#' 14(3), 520--525.
#'
#' @export es_from_linreg_b_se
#'
#' @md
#'
#' @examples
#' es_from_linreg_b_se(linreg_b = 1.5, linreg_b_se = 0.6, n_sample = 100, n_covariates = 2)
es_from_linreg_b_se <- function(linreg_b, linreg_b_se, n_sample, n_covariates,
                                sd_iv, unit_increase_iv, unit_type = "raw_scale",
                                n_exp, n_nexp, cor_to_smd = "viechtbauer",
                                reverse_linreg_b) {

  if (missing(reverse_linreg_b)) reverse_linreg_b <- rep(FALSE, length(linreg_b))
  reverse_linreg_b[is.na(reverse_linreg_b)] <- FALSE
  if (length(reverse_linreg_b) == 1) reverse_linreg_b <- rep(reverse_linreg_b, length(linreg_b))
  if (length(reverse_linreg_b) != length(linreg_b)) stop("The length of the 'reverse_linreg_b' argument is incorrectly specified.")

  if (missing(n_exp)) n_exp <- rep(NA, length(linreg_b))
  if (missing(n_nexp)) n_nexp <- rep(NA, length(linreg_b))
  if (missing(sd_iv)) sd_iv <- rep(NA, length(linreg_b))
  if (missing(unit_increase_iv)) unit_increase_iv <- rep(NA, length(linreg_b))
  if (missing(unit_type)) unit_type <- rep(NA, length(linreg_b))

  linreg_t <- linreg_b / linreg_b_se

  es <- es_from_linreg_t(
    linreg_t = linreg_t, n_sample = n_sample, n_covariates = n_covariates,
    sd_iv = sd_iv, unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_linreg_t = reverse_linreg_b
  )

  es$info_used <- "linreg_b_se"

  return(es)
}


#' Convert an unstandardized regression coefficient and its confidence interval
#' into several effect size measures
#'
#' @param linreg_b unstandardized regression coefficient from a linear regression model
#' @param linreg_b_ci_lo lower bound of the 95% confidence interval of the regression coefficient
#' @param linreg_b_ci_up upper bound of the 95% confidence interval of the regression coefficient
#' @param n_sample the total number of participants
#' @param n_covariates the number of covariates in the model (excluding the predictor of interest).
#' @param sd_iv the standard deviation of the independent variable (optional, see details)
#' @param unit_increase_iv a value of the independent variable that will be used to estimate the Cohen's d (optional, see details).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Must be either "sd" or "raw_scale"
#' @param n_exp number of the experimental/exposed group (optional)
#' @param n_nexp number of the non-experimental/non-exposed group (optional)
#' @param cor_to_smd formula used to convert a \code{pearson_r} or \code{fisher_z} value into a SMD.
#' @param reverse_linreg_b a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function derives the standard error from the 95% confidence interval using a
#' t-distribution with \eqn{df = n\_sample - n\_covariates - 1} degrees of freedom:
#' \deqn{SE(b) = \frac{ci\_up - ci\_lo}{2 \times qt(.975, df)}}
#'
#' Then, calculations of the \code{\link{es_from_linreg_b_se}} function are applied.
#' For binary predictors without covariates, use \code{\link{es_from_beta_unstd}} instead.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab Rp (partial correlation) + Zp (Fisher's z of partial r)\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + OR (approximate, see \code{\link{es_from_linreg_t}})\cr
#' }
#'
#' @references
#' Aloe, A. M., & Thompson, C. G. (2013). The synthesis of partial effect sizes.
#' \emph{Journal of the Society for Social Work and Research}, 4(4), 390--405.
#'
#' @export es_from_linreg_b_ci
#'
#' @md
#'
#' @examples
#' es_from_linreg_b_ci(
#'   linreg_b = 1.5, linreg_b_ci_lo = 0.3, linreg_b_ci_up = 2.7,
#'   n_sample = 100, n_covariates = 2
#' )
es_from_linreg_b_ci <- function(linreg_b, linreg_b_ci_lo, linreg_b_ci_up,
                                n_sample, n_covariates,
                                sd_iv, unit_increase_iv, unit_type = "raw_scale",
                                n_exp, n_nexp, cor_to_smd = "viechtbauer",
                                reverse_linreg_b) {

  if (missing(reverse_linreg_b)) reverse_linreg_b <- rep(FALSE, length(linreg_b))
  reverse_linreg_b[is.na(reverse_linreg_b)] <- FALSE
  if (length(reverse_linreg_b) == 1) reverse_linreg_b <- rep(reverse_linreg_b, length(linreg_b))
  if (length(reverse_linreg_b) != length(linreg_b)) stop("The length of the 'reverse_linreg_b' argument is incorrectly specified.")

  if (missing(n_exp)) n_exp <- rep(NA, length(linreg_b))
  if (missing(n_nexp)) n_nexp <- rep(NA, length(linreg_b))
  if (missing(sd_iv)) sd_iv <- rep(NA, length(linreg_b))
  if (missing(unit_increase_iv)) unit_increase_iv <- rep(NA, length(linreg_b))
  if (missing(unit_type)) unit_type <- rep(NA, length(linreg_b))

  df <- n_sample - n_covariates - 2
  linreg_b_se <- (linreg_b_ci_up - linreg_b_ci_lo) / (2 * qt(.975, df))

  es <- es_from_linreg_b_se(
    linreg_b = linreg_b, linreg_b_se = linreg_b_se,
    n_sample = n_sample, n_covariates = n_covariates,
    sd_iv = sd_iv, unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_linreg_b = reverse_linreg_b
  )

  es$info_used <- "linreg_b_ci"

  return(es)
}


#' Convert an unstandardized regression coefficient and its p-value
#' into several effect size measures
#'
#' @param linreg_b unstandardized regression coefficient from a linear regression model
#' @param linreg_b_pval two-sided p-value of the regression coefficient
#' @param n_sample the total number of participants
#' @param n_covariates the number of covariates in the model (excluding the predictor of interest).
#' @param sd_iv the standard deviation of the independent variable (optional, see details)
#' @param unit_increase_iv a value of the independent variable that will be used to estimate the Cohen's d (optional, see details).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Must be either "sd" or "raw_scale"
#' @param n_exp number of the experimental/exposed group (optional)
#' @param n_nexp number of the non-experimental/non-exposed group (optional)
#' @param cor_to_smd formula used to convert a \code{pearson_r} or \code{fisher_z} value into a SMD.
#' @param reverse_linreg_b_pval a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function recovers the t-statistic from the two-sided p-value:
#' \deqn{t = qt(1 - pval/2, df) \times sign(b)}
#' where \eqn{df = n\_sample - n\_covariates - 1}.
#'
#' The sign of the regression coefficient gives the direction, which the two-tailed
#' p-value does not carry.
#' Then, calculations of the \code{\link{es_from_linreg_t}} function are applied.
#' For binary predictors without covariates, use \code{\link{es_from_beta_unstd}} instead.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab Rp (partial correlation) + Zp (Fisher's z of partial r)\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + OR (approximate, see \code{\link{es_from_linreg_t}})\cr
#' }
#'
#' @references
#' Aloe, A. M., & Thompson, C. G. (2013). The synthesis of partial effect sizes.
#' \emph{Journal of the Society for Social Work and Research}, 4(4), 390--405.
#'
#' @export es_from_linreg_b_pval
#'
#' @md
#'
#' @examples
#' es_from_linreg_b_pval(
#'   linreg_b = 1.5, linreg_b_pval = 0.01,
#'   n_sample = 100, n_covariates = 2
#' )
es_from_linreg_b_pval <- function(linreg_b, linreg_b_pval,
                                  n_sample, n_covariates,
                                  sd_iv, unit_increase_iv, unit_type = "raw_scale",
                                  n_exp, n_nexp, cor_to_smd = "viechtbauer",
                                  reverse_linreg_b_pval) {

  if (missing(reverse_linreg_b_pval)) reverse_linreg_b_pval <- rep(FALSE, length(linreg_b))
  reverse_linreg_b_pval[is.na(reverse_linreg_b_pval)] <- FALSE
  if (length(reverse_linreg_b_pval) == 1) reverse_linreg_b_pval <- rep(reverse_linreg_b_pval, length(linreg_b))
  if (length(reverse_linreg_b_pval) != length(linreg_b)) stop("The length of the 'reverse_linreg_b_pval' argument is incorrectly specified.")

  if (missing(n_exp)) n_exp <- rep(NA, length(linreg_b))
  if (missing(n_nexp)) n_nexp <- rep(NA, length(linreg_b))
  if (missing(sd_iv)) sd_iv <- rep(NA, length(linreg_b))
  if (missing(unit_increase_iv)) unit_increase_iv <- rep(NA, length(linreg_b))
  if (missing(unit_type)) unit_type <- rep(NA, length(linreg_b))

  df <- n_sample - n_covariates - 2
  linreg_t <- qt(1 - linreg_b_pval / 2, df) * sign(linreg_b)

  es <- es_from_linreg_t(
    linreg_t = linreg_t, n_sample = n_sample, n_covariates = n_covariates,
    sd_iv = sd_iv, unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_linreg_t = reverse_linreg_b_pval
  )

  es$info_used <- "linreg_b_pval"

  return(es)
}
