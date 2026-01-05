#' Converts the means and bounds of an error bar (generally extracted from a plot) into four effect measures (SMD, MD, OR, COR)
#'
#' @param plot_mean_exp mean of participants in the experimental/exposed group (extracted from a plot).
#' @param plot_mean_nexp mean of participants in the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_mean_sd_lo_exp lower bound of an error bar depicting -1 SD from the mean of the experimental/exposed group (extracted from a plot).
#' @param plot_mean_sd_lo_nexp lower bound of an error bar depicting -1 SD from the mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_mean_sd_up_exp upper bound of an error bar depicting +1 SD from the mean of the experimental/exposed group (extracted from a plot).
#' @param plot_mean_sd_up_nexp upper bound of an error bar depicting +1 SD from the mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_mean_se_lo_exp lower bound of an error bar depicting -1 SE from the mean of the experimental/exposed group (extracted from a plot).
#' @param plot_mean_se_lo_nexp lower bound of an error bar depicting -1 SE from the mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_mean_se_up_exp upper bound of an error bar depicting +1 SE from the mean of the experimental/exposed group (extracted from a plot).
#' @param plot_mean_se_up_nexp upper bound of an error bar depicting +1 SE from the mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_mean_ci_lo_exp lower bound of an error bar depicting the 95% CI of the mean of the experimental/exposed group (extracted from a plot).
#' @param plot_mean_ci_lo_nexp lower bound of an error bar depicting the 95% CI of the mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_mean_ci_up_exp upper bound of an error bar depicting the 95% CI of the mean of the experimental/exposed group (extracted from a plot).
#' @param plot_mean_ci_up_nexp upper bound of an error bar depicting the 95% CI of the mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_plot_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function uses the bounds of an error bar of a mean obtained from a plot into a standard deviation.
#' Then, a mean difference (MD), Cohen's d (D), and Hedges' g (G) are estimated.
#' Odds ratio (OR), risk ratio (RR) and correlation coefficients (R/Z) are converted from the Cohen's d value.
#'
#' **To convert the bound of an error bar into a standard deviation,**
#' this function always prioritizes information from the \code{plot_mean_sd_*} arguments,
#' then those from the \code{plot_mean_se_*} arguments,
#' then those from the \code{plot_mean_ci_*} arguments.
#'
#' 1. If the bounds of the standard deviations are provided, the following formulas are used:
#' \deqn{mean\_sd\_lo\_exp = plot\_mean\_exp - plot\_mean\_sd\_lo\_exp}
#' \deqn{mean\_sd\_up\_exp = plot\_mean\_sd\_up\_exp - plot\_mean\_exp}
#' \deqn{mean\_sd\_exp = \frac{mean\_sd\_lo\_exp + mean\_sd\_up\_exp}{2}}
#'
#' \deqn{mean\_sd\_lo\_nexp = plot\_mean\_nexp - plot\_mean\_sd\_lo\_nexp}
#' \deqn{mean\_sd\_up\_nexp = plot\_mean\_sd\_up\_nexp - plot\_mean\_nexp}
#' \deqn{mean\_sd\_nexp = \frac{mean\_sd\_lo\_nexp + mean\_sd\_up\_nexp}{2}}
#'
#' Note that if only one bound (e.g., the upper bound) is provided, it will be the
#' only information used to estimate the standard deviation value.
#'
#' Then, calculations of the \code{\link{es_from_means_sd}} are used.
#'
#' 2. If the bounds of the standard errors are provided, the following formulas are used:
#' \deqn{mean\_se\_lo\_exp = plot\_mean\_exp - plot\_mean\_se\_lo\_exp}
#' \deqn{mean\_se\_up\_exp = plot\_mean\_se\_up\_exp - plot\_mean\_exp}
#' \deqn{mean\_se\_exp = \frac{mean\_se\_lo\_exp + mean\_se\_up\_exp}{2}}
#'
#' \deqn{mean\_se\_lo\_nexp = plot\_mean\_nexp - plot\_mean\_se\_lo\_nexp}
#' \deqn{mean\_se\_up\_nexp = plot\_mean\_se\_up\_nexp - plot\_mean\_nexp}
#' \deqn{mean\_se\_nexp = \frac{mean\_se\_lo\_nexp + mean\_se\_up\_nexp}{2}}
#'
#' Note that if only one bound (e.g., the upper bound) is provided, it will be the
#' only information used to estimate the standard error value.
#'
#' Then, calculations of the \code{\link{es_from_means_se}()} are used.
#'
#' 3. If the bounds of the 95% confidence intervals are provided, the calculations
#' of the \code{\link{es_from_means_ci}} are used.
#'
#' @export es_from_plot_means
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 21. From plot: means and dispersion (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_plot_means(
#'   n_exp = 35, n_nexp = 35,
#'   plot_mean_exp = 89, plot_mean_nexp = 104,
#'   plot_mean_sd_lo_exp = 69, plot_mean_sd_lo_nexp = 83,
#'   plot_mean_sd_up_exp = 109, plot_mean_sd_up_nexp = 125
#' )
es_from_plot_means <- function(n_exp, n_nexp,
                               plot_mean_exp, plot_mean_nexp,
                               plot_mean_sd_lo_exp, plot_mean_sd_lo_nexp,
                               plot_mean_sd_up_exp, plot_mean_sd_up_nexp,
                               plot_mean_se_lo_exp, plot_mean_se_lo_nexp,
                               plot_mean_se_up_exp, plot_mean_se_up_nexp,
                               plot_mean_ci_lo_exp, plot_mean_ci_lo_nexp,
                               plot_mean_ci_up_exp, plot_mean_ci_up_nexp,
                               smd_to_cor = "viechtbauer", reverse_plot_means) {
  if (missing(plot_mean_sd_lo_exp)) plot_mean_sd_lo_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_sd_lo_nexp)) plot_mean_sd_lo_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_sd_up_exp)) plot_mean_sd_up_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_sd_up_nexp)) plot_mean_sd_up_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_se_lo_exp)) plot_mean_se_lo_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_se_lo_nexp)) plot_mean_se_lo_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_se_up_exp)) plot_mean_se_up_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_se_up_nexp)) plot_mean_se_up_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_ci_lo_exp)) plot_mean_ci_lo_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_ci_lo_nexp)) plot_mean_ci_lo_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_ci_up_exp)) plot_mean_ci_up_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_mean_ci_up_nexp)) plot_mean_ci_up_nexp <- rep(NA_real_, length(n_exp))

  if (missing(reverse_plot_means)) reverse_plot_means <- rep(FALSE, length(n_exp))
  reverse_plot_means[is.na(reverse_plot_means)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       error_message = paste0("The number of people exposed/non-exposed ",
                                              "should be >0."),
                       func = "es_from_plot_means")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  ## SD ------
  #### exp
  mean_sd_exp_lo <- plot_mean_exp - plot_mean_sd_lo_exp
  mean_sd_exp_up <- plot_mean_sd_up_exp - plot_mean_exp
  mean_sd_exp_transit <- apply(cbind(mean_sd_exp_lo, mean_sd_exp_up), 1, mean, na.rm = TRUE)
  mean_sd_exp <- ifelse(is.na(mean_sd_exp_lo) &
                          is.na(mean_sd_exp_up), NA_real_,
                        mean_sd_exp_transit)
  #### nexp
  mean_sd_nexp_lo <- plot_mean_nexp - plot_mean_sd_lo_nexp
  mean_sd_nexp_up <- plot_mean_sd_up_nexp - plot_mean_nexp
  mean_sd_nexp_transit <- apply(cbind(mean_sd_nexp_lo, mean_sd_nexp_up), 1, mean, na.rm = TRUE)
  mean_sd_nexp <- ifelse(is.na(mean_sd_nexp_lo) &
                           is.na(mean_sd_nexp_up), NA_real_,
                         mean_sd_nexp_transit)

  ## SE ------
  #### exp

  mean_se_exp_lo <- plot_mean_exp - plot_mean_se_lo_exp
  mean_se_exp_up <- plot_mean_se_up_exp - plot_mean_exp
  mean_se_exp_transit <- apply(cbind(mean_se_exp_lo, mean_se_exp_up), 1, mean, na.rm = TRUE)
  mean_se_exp <- ifelse(is.na(mean_se_exp_lo) & is.na(mean_se_exp_up),
                        NA_real_, mean_se_exp_transit)
  #### nexp
  mean_se_nexp_lo <- plot_mean_nexp - plot_mean_se_lo_nexp
  mean_se_nexp_up <- plot_mean_se_up_nexp - plot_mean_nexp
  mean_se_nexp_transit <- apply(cbind(mean_se_nexp_lo, mean_se_nexp_up), 1, mean, na.rm = TRUE)
  mean_se_nexp <- ifelse(is.na(mean_se_nexp_lo) & is.na(mean_se_nexp_up), NA_real_,
                         mean_se_nexp_transit)

  es <- es_from_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = plot_mean_exp, mean_sd_exp = mean_sd_exp,
    mean_nexp = plot_mean_nexp, mean_sd_nexp = mean_sd_nexp,
    smd_to_cor = smd_to_cor, reverse_means = reverse_plot_means
  )

  es_1 <- es_from_means_se(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = plot_mean_exp, mean_se_exp = mean_se_exp,
    mean_nexp = plot_mean_nexp, mean_se_nexp = mean_se_nexp,
    smd_to_cor = smd_to_cor, reverse_means = reverse_plot_means
  )

  es_2 <- es_from_means_ci(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = plot_mean_exp, mean_nexp = plot_mean_nexp,
    mean_ci_lo_exp = plot_mean_ci_lo_exp, mean_ci_lo_nexp = plot_mean_ci_lo_nexp,
    mean_ci_up_exp = plot_mean_ci_up_exp, mean_ci_up_nexp = plot_mean_ci_up_nexp,
    smd_to_cor = smd_to_cor, reverse_means = reverse_plot_means
  )

  row_miss <- which(is.na(es$d) & is.na(es$d_se))
  if (length(row_miss) > 0) {
    es[row_miss, ] <- es_1[row_miss, ]
    row_miss_2 <- which(is.na(es$d) & is.na(es$d_se))
    if (length(row_miss_2) > 0) {
      es[row_miss_2, ] <- es_2[row_miss_2, ]
    }
  }

  es$info_used <- "means_plot"

  return(es)
}

#' Converts the means and bounds of an error bar (generally extracted from a plot) into four effect measures (SMD, MD, OR, COR)
#'
#' @param plot_ancova_mean_exp ancova_mean of participants in the experimental/exposed group (extracted from a plot).
#' @param plot_ancova_mean_nexp ancova_mean of participants in the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_ancova_mean_sd_lo_exp lower bound of an error bar depicting -1 SD from the ancova_mean of the experimental/exposed group (extracted from a plot).
#' @param plot_ancova_mean_sd_lo_nexp lower bound of an error bar depicting -1 SD from the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_ancova_mean_sd_up_exp upper bound of an error bar depicting +1 SD from the ancova_mean of the experimental/exposed group (extracted from a plot).
#' @param plot_ancova_mean_sd_up_nexp upper bound of an error bar depicting +1 SD from the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_ancova_mean_se_lo_exp lower bound of an error bar depicting -1 SE from the ancova_mean of the experimental/exposed group (extracted from a plot).
#' @param plot_ancova_mean_se_lo_nexp lower bound of an error bar depicting -1 SE from the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_ancova_mean_se_up_exp upper bound of an error bar depicting +1 SE from the ancova_mean of the experimental/exposed group (extracted from a plot).
#' @param plot_ancova_mean_se_up_nexp upper bound of an error bar depicting +1 SE from the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_ancova_mean_ci_lo_exp lower bound of an error bar depicting the 95% CI of the ancova_mean of the experimental/exposed group (extracted from a plot).
#' @param plot_ancova_mean_ci_lo_nexp lower bound of an error bar depicting the 95% CI of the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param plot_ancova_mean_ci_up_exp upper bound of an error bar depicting the 95% CI of the ancova_mean of the experimental/exposed group (extracted from a plot).
#' @param plot_ancova_mean_ci_up_nexp upper bound of an error bar depicting the 95% CI of the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).
#' @param cov_outcome_r correlation between the outcome and covariate (multiple correlation when multiple covariates are included in the ANCOVA model).
#' @param n_cov_ancova number of covariates in the ANCOVA model.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).
#' @param reverse_plot_ancova_means a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function uses the bounds of an error bar of a mean obtained from a plot into a standard deviation.
#' Then, a mean difference (MD), Cohen's d (D), and Hedges' g (G) are estimated.
#' Odds ratio (OR), risk ratio (RR) and correlation coefficients (R/Z) are converted from the Cohen's d value.
#'
#' **To convert the bound of an error bar into a standard deviation,**
#' this function always prioritizes information from the \code{plot_ancova_mean_sd_*} arguments,
#' then those from the \code{plot_ancova_mean_se_*} arguments,
#' then those from the \code{plot_ancova_mean_ci_*} arguments.
#'
#' 1. If the bounds of the standard deviations are provided, the following formulas are used:
#' \deqn{ancova\_mean\_sd\_lo\_exp = plot\_ancova\_mean\_exp - plot\_ancova\_mean\_sd\_lo\_exp}
#' \deqn{ancova\_mean\_sd\_up\_exp = plot\_ancova\_mean\_sd\_up\_exp - plot\_ancova\_mean\_exp}
#' \deqn{ancova\_mean\_sd\_exp = \frac{ancova\_mean\_sd\_lo\_exp + ancova\_mean\_sd\_up\_exp}{2}}
#'
#' \deqn{mean\_sd\_lo\_nexp = plot\_ancova\_mean\_nexp - plot\_ancova\_mean\_sd\_lo\_nexp}
#' \deqn{mean\_sd\_up\_nexp = plot\_ancova\_mean\_sd\_up\_nexp - plot\_ancova\_mean\_nexp}
#' \deqn{mean\_sd\_nexp = \frac{mean\_sd\_lo\_nexp + mean\_sd\_up\_nexp}{2}}
#'
#' Then, calculations of the \code{\link{es_from_ancova_means_sd}} are used.
#'
#' 2. If the bounds of the standard errors are provided, the following formulas are used:
#' \deqn{ancova\_mean\_se\_lo\_exp = plot\_ancova\_mean\_exp - plot\_ancova\_mean\_se\_lo\_exp}
#' \deqn{ancova\_mean\_se\_up\_exp = plot\_ancova\_mean\_se\_up\_exp - plot\_ancova\_mean\_exp}
#' \deqn{ancova\_mean\_se\_exp = \frac{ancova\_mean\_se\_lo\_exp + ancova\_mean\_se\_up\_exp}{2}}
#'
#' \deqn{mean\_se\_lo\_nexp = plot\_ancova\_mean\_nexp - plot\_ancova\_mean\_se\_lo\_nexp}
#' \deqn{mean\_se\_up\_nexp = plot\_ancova\_mean\_se\_up\_nexp - plot\_ancova\_mean\_nexp}
#' \deqn{mean\_se\_nexp = \frac{mean\_se\_lo\_nexp + mean\_se\_up\_nexp}{2}}
#'
#' Then, calculations of the \code{\link{es_from_ancova_means_se}} are used.
#'
#' 3. If the bounds of the 95% confidence intervals are provided, the calculations
#' of the \code{\link{es_from_ancova_means_ci}()} are used.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab MD + D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 22. From plot: adjusted means and dispersion (adjusted)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_plot_ancova_means
#'
#' @md
#'
#' @examples
#' es_from_plot_ancova_means(
#'   n_exp = 35, n_nexp = 35,
#'   cov_outcome_r = 0.2, n_cov_ancova = 4,
#'   plot_ancova_mean_exp = 89, plot_ancova_mean_nexp = 104,
#'   plot_ancova_mean_sd_lo_exp = 69, plot_ancova_mean_sd_lo_nexp = 83,
#'   plot_ancova_mean_sd_up_exp = 109, plot_ancova_mean_sd_up_nexp = 125
#' )
es_from_plot_ancova_means <- function(n_exp, n_nexp,
                                      plot_ancova_mean_exp, plot_ancova_mean_nexp,
                                      plot_ancova_mean_sd_lo_exp, plot_ancova_mean_sd_lo_nexp,
                                      plot_ancova_mean_sd_up_exp, plot_ancova_mean_sd_up_nexp,
                                      plot_ancova_mean_se_lo_exp, plot_ancova_mean_se_lo_nexp,
                                      plot_ancova_mean_se_up_exp, plot_ancova_mean_se_up_nexp,
                                      plot_ancova_mean_ci_lo_exp, plot_ancova_mean_ci_lo_nexp,
                                      plot_ancova_mean_ci_up_exp, plot_ancova_mean_ci_up_nexp,
                                      cov_outcome_r, n_cov_ancova,
                                      smd_to_cor = "viechtbauer", reverse_plot_ancova_means) {
  if (missing(plot_ancova_mean_sd_lo_exp)) plot_ancova_mean_sd_lo_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_sd_lo_nexp)) plot_ancova_mean_sd_lo_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_sd_up_exp)) plot_ancova_mean_sd_up_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_sd_up_nexp)) plot_ancova_mean_sd_up_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_se_lo_exp)) plot_ancova_mean_se_lo_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_se_lo_nexp)) plot_ancova_mean_se_lo_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_se_up_exp)) plot_ancova_mean_se_up_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_se_up_nexp)) plot_ancova_mean_se_up_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_ci_lo_exp)) plot_ancova_mean_ci_lo_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_ci_lo_nexp)) plot_ancova_mean_ci_lo_nexp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_ci_up_exp)) plot_ancova_mean_ci_up_exp <- rep(NA_real_, length(n_exp))
  if (missing(plot_ancova_mean_ci_up_nexp)) plot_ancova_mean_ci_up_nexp <- rep(NA_real_, length(n_exp))

  if (missing(reverse_plot_ancova_means)) reverse_plot_ancova_means <- rep(FALSE, length(n_exp))
  reverse_plot_ancova_means[is.na(reverse_plot_ancova_means)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       cov_outcome_r, n_cov_ancova,
                       error_message = paste0("The number of people exposed/non-exposed ",
                                              "as well as the correlation and number of covariates in ANCOVA ",
                                              "should be >0."),
                       func = "es_from_plot_ancova_means")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  ## SD ------
  #### exp
  ancova_mean_sd_exp_lo <- plot_ancova_mean_exp - plot_ancova_mean_sd_lo_exp
  ancova_mean_sd_exp_up <- plot_ancova_mean_sd_up_exp - plot_ancova_mean_exp
  ancova_mean_sd_exp_transit <- apply(cbind(ancova_mean_sd_exp_lo, ancova_mean_sd_exp_up), 1, mean, na.rm = TRUE)
  ancova_mean_sd_exp <- ifelse(is.na(ancova_mean_sd_exp_lo) & is.na(ancova_mean_sd_exp_up), NA_real_, ancova_mean_sd_exp_transit)
  #### nexp
  ancova_mean_sd_nexp_lo <- plot_ancova_mean_nexp - plot_ancova_mean_sd_lo_nexp
  ancova_mean_sd_nexp_up <- plot_ancova_mean_sd_up_nexp - plot_ancova_mean_nexp
  ancova_mean_sd_nexp_transit <- apply(cbind(ancova_mean_sd_nexp_lo, ancova_mean_sd_nexp_up), 1, mean, na.rm = TRUE)
  ancova_mean_sd_nexp <- ifelse(is.na(ancova_mean_sd_nexp_lo) & is.na(ancova_mean_sd_nexp_up), NA_real_, ancova_mean_sd_nexp_transit)

  ## SE ------
  #### exp
  ancova_mean_se_exp_lo <- plot_ancova_mean_exp - plot_ancova_mean_se_lo_exp
  ancova_mean_se_exp_up <- plot_ancova_mean_se_up_exp - plot_ancova_mean_exp
  ancova_mean_se_exp_transit <- apply(cbind(ancova_mean_se_exp_lo, ancova_mean_se_exp_up), 1, mean, na.rm = TRUE)
  ancova_mean_se_exp <- ifelse(is.na(ancova_mean_se_exp_lo) & is.na(ancova_mean_se_exp_up), NA_real_, ancova_mean_se_exp_transit)
  # ancova_mean_sd_conv_exp = ancova_mean_se_exp * sqrt(n_exp)
  # ancova_mean_sd_conv_exp = ancova_mean_se_exp * sqrt(n_exp - 1)
  #### nexp
  ancova_mean_se_nexp_lo <- plot_ancova_mean_nexp - plot_ancova_mean_se_lo_nexp
  ancova_mean_se_nexp_up <- plot_ancova_mean_se_up_nexp - plot_ancova_mean_nexp
  ancova_mean_se_nexp_transit <- apply(cbind(ancova_mean_se_nexp_lo, ancova_mean_se_nexp_up), 1, mean, na.rm = TRUE)
  ancova_mean_se_nexp <- ifelse(is.na(ancova_mean_se_nexp_lo) & is.na(ancova_mean_se_nexp_up), NA_real_, ancova_mean_se_nexp_transit)
  # ancova_mean_sd_conv_nexp = ancova_mean_se_nexp * sqrt(n_nexp - 1)
  # ancova_mean_sd_conv_nexp = ancova_mean_se_nexp * sqrt(n_nexp)

  es <- es_from_ancova_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    ancova_mean_exp = plot_ancova_mean_exp,
    ancova_mean_sd_exp = ancova_mean_sd_exp,
    ancova_mean_nexp = plot_ancova_mean_nexp,
    ancova_mean_sd_nexp = ancova_mean_sd_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_plot_ancova_means
  )
  es_1 <- es_from_ancova_means_se(
    n_exp = n_exp, n_nexp = n_nexp,
    ancova_mean_exp = plot_ancova_mean_exp,
    ancova_mean_se_exp = ancova_mean_se_exp,
    ancova_mean_nexp = plot_ancova_mean_nexp,
    ancova_mean_se_nexp = ancova_mean_se_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_plot_ancova_means
  )
  es_2 <- es_from_ancova_means_ci(
    n_exp = n_exp, n_nexp = n_nexp,
    ancova_mean_exp = plot_ancova_mean_exp,
    ancova_mean_ci_lo_exp = plot_ancova_mean_ci_lo_exp,
    ancova_mean_ci_up_exp = plot_ancova_mean_ci_up_exp,
    ancova_mean_nexp = plot_ancova_mean_nexp,
    ancova_mean_ci_lo_nexp = plot_ancova_mean_ci_lo_nexp,
    ancova_mean_ci_up_nexp = plot_ancova_mean_ci_up_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_plot_ancova_means
  )
  # ## CI ------
  # #### exp
  # ancova_mean_se_exp_conv = (plot_ancova_mean_ci_up_exp - plot_ancova_mean_ci_lo_exp) / (2 * qt(0.975, n_exp - 1))
  # ancova_mean_sd_conv_exp_ci = ancova_mean_se_exp_conv * sqrt(n_exp)
  # # ancova_mean_sd_conv_exp_ci = ancova_mean_se_exp_conv * sqrt(n_exp - 1)
  # #### nexp
  # ancova_mean_se_nexp_conv = (plot_ancova_mean_ci_up_nexp - plot_ancova_mean_ci_lo_nexp) / (2 * qt(0.975, n_nexp - 1))
  # # ancova_mean_sd_conv_nexp_ci = ancova_mean_se_nexp_conv * sqrt(n_nexp - 1)
  # ancova_mean_sd_conv_nexp_ci = ancova_mean_se_nexp_conv * sqrt(n_nexp)
  #
  #
  # es_1 = es_from_ancova_means_sd(n_exp = n_exp, n_nexp = n_nexp,
  #                             ancova_mean_exp = plot_ancova_mean_exp, ancova_mean_sd_exp = ancova_mean_sd_conv_exp,
  #                             ancova_mean_nexp = plot_ancova_mean_nexp, ancova_mean_sd_nexp = ancova_mean_sd_conv_nexp,
  #                             cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
  #                             smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_plot_ancova_means)
  #
  # es_2 = es_from_ancova_means_sd(n_exp = n_exp, n_nexp = n_nexp,
  #                             ancova_mean_exp = plot_ancova_mean_exp, ancova_mean_sd_exp = ancova_mean_sd_conv_exp_ci,
  #                             ancova_mean_nexp = plot_ancova_mean_nexp, ancova_mean_sd_nexp = ancova_mean_sd_conv_nexp_ci,
  #                             cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
  #                             smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_plot_ancova_means)

  row_miss <- which(is.na(es$d) & is.na(es$d_se))
  if (length(row_miss) > 0) {
    es[row_miss, ] <- es_1[row_miss, ]
    row_miss_2 <- which(is.na(es$d) & is.na(es$d_se))
    if (length(row_miss_2) > 0) {
      es[(row_miss_2), ] <- es_2[(row_miss_2), ]
    }
  }
  es$info_used <- "ancova_means_plot"

  return(es)
}

# m1 = rnorm(40, 50, 10)
# m2 = rnorm(40, 40, 10)
# sd = c(rnorm(15, 10, 2), rep(NA, 25))
# se = c(rep(NA, 15), rnorm(15, 8, 2), rep(NA, 10))
# ci = c(rep(NA, 30), abs(rnorm(10, 4, 1)))
# n_exp = runif(40, 10, 50)
# n_nexp = runif(40, 10, 50)
# reverse_means = rep(NA, 40)
# dat = data.frame(
#   n_exp = runif(40, 10, 50),
#   n_nexp = runif(40, 10, 50),
#   plot_mean_exp = m1,
#   plot_mean_nexp = m2,
#   plot_mean_sd_lo_exp = c(m1 - sd),
#   plot_mean_sd_lo_nexp = c(m2 - sd),
#   plot_mean_sd_up_exp = c(m1 + sd),
#   plot_mean_sd_up_nexp = c(m2 + sd),
#   plot_mean_se_lo_exp = c( m1 - se),
#   plot_mean_se_lo_nexp = c(m2 - se),
#   plot_mean_se_up_exp = c(m1 + se),
#   plot_mean_se_up_nexp = c(m2 + se),
#   plot_mean_ci_lo_exp = c(m1 - ci),
#   plot_mean_ci_lo_nexp = m2 - ci,
#   plot_mean_ci_up_exp = m1 + ci,
#   plot_mean_ci_up_nexp = m2 + ci)
#
# m1 = rnorm(40, 50, 10)
# m2 = rnorm(40, 40, 10)
# sd = c(rnorm(15, 10, 2), rep(NA, 25))
# se = c(rep(NA, 15), rnorm(15, 8, 2), rep(NA, 10))
# ci = c(rep(NA, 30), abs(rnorm(10, 4, 1)))
#
#   n_exp = runif(40, 10, 50)
#   n_nexp = runif(40, 10, 50)
#   plot_mean_exp = m1
#   plot_mean_nexp = m2
#   plot_mean_sd_lo_exp = c(m1 - sd)
#   plot_mean_sd_lo_nexp = c(m2 - sd)
#   plot_mean_sd_up_exp = c(m1 + sd)
#   plot_mean_sd_up_nexp = c(m2 + sd)
#   plot_mean_se_lo_exp = c( m1 - se)
#   plot_mean_se_lo_nexp = c(m2 - se)
#   plot_mean_se_up_exp = c(m1 + se)
#   plot_mean_se_up_nexp = c(m2 + se)
#   plot_mean_ci_lo_exp = c(m1 - ci)
#   plot_mean_ci_lo_nexp = m2 - ci
#   plot_mean_ci_up_exp = m1 + ci
#   plot_mean_ci_up_nexp = m2 + ci
#
#   res = es_from_plot_means(n_exp = dat$n_exp, n_nexp=dat$n_nexp,
#                      plot_mean_exp=dat$plot_mean_exp,
#                      plot_mean_nexp = dat$plot_mean_nexp,
#                      plot_mean_sd_lo_exp=dat$plot_mean_sd_lo_exp, plot_mean_sd_lo_nexp=dat$plot_mean_sd_lo_nexp,
#                      plot_mean_sd_up_exp=dat$plot_mean_sd_up_exp, plot_mean_sd_up_nexp=dat$plot_mean_sd_up_nexp,
#                      plot_mean_se_lo_exp=dat$plot_mean_se_lo_exp, plot_mean_se_lo_nexp=dat$plot_mean_se_lo_nexp,
#                      plot_mean_se_up_exp=dat$plot_mean_se_up_exp, plot_mean_se_up_nexp=dat$plot_mean_se_up_nexp,
#                      plot_mean_ci_lo_exp=dat$plot_mean_ci_lo_exp, plot_mean_ci_lo_nexp=dat$plot_mean_ci_lo_nexp,
#                      plot_mean_ci_up_exp=dat$plot_mean_ci_up_exp, plot_mean_ci_up_nexp=dat$plot_mean_ci_up_nexp,
#                      reverse_plot_means)
