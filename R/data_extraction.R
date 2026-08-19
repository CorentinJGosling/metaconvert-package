#' Data extraction sheet generator
#'
#' @param name Name of the file created
#' @param measure Target effect size measure (one of the 14 available in metaConvert). Default is "all".
#' @param extension Extension of the file created. Most common are ".xlsx", ".csv" or ".txt". It is also possible to generate an R dataframe object by using the "data.frame" extension.
#' @param type_of_measure One of "natural+converted" or "natural" (see details).
#' @param verbose logical variable indicating whether some information should be printed (e.g., the location where the sheet is created when using ".xlsx", ".csv" or ".txt" extensions)
#'
#' @details
#' This function generates, on your computer, a data extraction sheet that contains the name of columns
#' that can be used by our tools to estimate various effect size measures.
#'
#' If you select a specific measure (e.g., \code{measure = "g"}), you will be presented only with most common
#' information allowing to estimate this measure (e.g., you will not be provided with columns for contingency
#' tables if you request a data extraction sheet for \code{measure = "g"}).
#'
#'
#' ## Measure
#' You can specify a specific effect size measures (among those available in the \code{\link{metaConvert-package}}).
#' Doing this, the data extraction sheet will contain only the columns of the input data allowing **a natural estimation
#' of the effect size measure**. For example, if you request \code{measure="d"} the data extraction
#' sheet will not contain the columns for the contingency table since, although the \code{\link{convert_df}}
#' function allows you to **convert** a contingency table into a "d", this requires to convert the "OR"
#' that is naturally estimated from the contingency table into a "d".
#'
#'
#' This table is designed to be used in combination with tables showing the combination of input data
#' leading to estimate each of the effect size measures (https://metaconvert.org/input.html)
#'
#' ## Extension
#' You can export a file in various formats outside R
#' (by indicating, for example, ".txt", ".xlsx", or ".csv") in the \code{extension} argument.
#' You can also visualise this dataset directly in R by setting \code{extension = "data.frame"}.
#'
#' @return
#' This function returns a data extraction sheet that contains all the information
#' necessary to estimate any effect size using the metaConvert tools.
#'
#' @md
#'
#' @export data_extraction_sheet
#'
#' @examples
#' data_extraction_sheet(measure = "md", extension = "data.frame")
data_extraction_sheet <- function(measure = c("d", "g", "md", "dw", "gw", "mdw",
                                              "or", "rr", "nnt", "rd",
                                              "r", "z", "logvr", "logcvr", "irr",
                                              "prop"),
                                  type_of_measure = c("natural", "natural+converted"),
                                  name = "mcv_data_extraction",
                                  extension = c("data.frame", ".txt", ".csv", ".xlsx"),
                                  verbose = TRUE
                                  ) {
  type_of_measure = type_of_measure[1]
  extension = extension[1]
  measure = measure[1]
  # NB the sixth name is "info_expected", NOT "all_info_expected". This line used to
  # emit "all_info_expected", a string that appeared nowhere else in the
  # package: .check_data() (R/internal_check_data.R:47) and summary() (
  # R/functions_summary.R:598) both read "info_expected". A user who filled the sheet
  # exactly as generated therefore had the value silently discarded -- .check_data()
  # manufactured info_expected all-NA and functions_summary.R dropped it, so the
  # documented column (man/summary.metaConvert.Rd) never populated for anyone using
  # the generated sheet. Pinned by tests/testthat/test-extraction-sheet-column-names.R.
  cols_req = c("study_id", "author", "year", "predictor", "outcome", "info_expected")
  inf_req = c("unique ID for each study - numeric/character",
              "first author name - numeric/character",
              "year of publication - numeric/character",
              "intervention/exposure - numeric/character",
              "outcome - numeric/character",
              "expected input data leading to an effect size estimate - character")

  cols_sample = "n_sample"
  inf_sample = c("total number of participants - numeric")

  cols_exposed = c("n_exp", "n_nexp")
  inf_exposed = c("number of participants in the exposed/experimental group - numeric",
                  "number of participants in the non-exposed/non-experimental group - numeric")

  cols_cases = c("n_cases", "n_controls")
  inf_cases = c("number of cases/events across exposed/non-exposed groups - numeric",
                "number of controls/no-event across exposed/non-exposed groups - numeric")

  cols_input_crude = c("user_es_original_measure_crude",
    "user_es_target_measure_crude",
    "user_es_crude",
    "user_se_crude",
    "user_ci_lo_crude",
    "user_ci_up_crude")
  inf_input_crude = c("type of ES provided: d, g, or, logor, rr, logrr, r, z, rd, etc. - character",
                      "target ES for conversion (defaults to measure param) - character",
                      "value of an effect size - numeric",
                      "standard error of the effect size - numeric",
                      "lower bound of the 95% CI of the effect size measure - numeric",
                      "upper bound of the 95% CI of the effect size measure - numeric")

  cols_input_adjusted = c(
    "user_es_original_measure_adj",
    "user_es_target_measure_adj",
    "user_es_adj",
    "user_se_adj",
    "user_ci_lo_adj",
    "user_ci_up_adj")
  inf_input_adjusted = c("type of adjusted ES provided: d, g, or, logor, rr, logrr, r, z, rd, etc. - character",
                      "target ES for conversion (defaults to measure param) - character",
                      "value of the adjusted effect size - numeric",
                      "standard error of the effect size - numeric",
                      "adjusted lower bound of the 95% CI of the effect size measure - numeric",
                      "adjusted upper bound of the 95% CI of the effect size measure - numeric")

  cols_means_post = c(
    # means post
    "reverse_means",
    "mean_exp", "mean_sd_exp", "mean_se_exp", "mean_ci_lo_exp", "mean_ci_up_exp",
    "mean_nexp", "mean_sd_nexp", "mean_se_nexp", "mean_ci_lo_nexp", "mean_ci_up_nexp",
    "mean_sd_pooled",

    # plot
    "reverse_plot_means",
    "plot_mean_exp", "plot_mean_nexp",
    "plot_mean_sd_lo_exp", "plot_mean_sd_lo_nexp",
    "plot_mean_sd_up_exp", "plot_mean_sd_up_nexp",
    "plot_mean_se_lo_exp", "plot_mean_se_lo_nexp",
    "plot_mean_se_up_exp", "plot_mean_se_up_nexp",
    "plot_mean_ci_lo_exp", "plot_mean_ci_lo_nexp",
    "plot_mean_ci_up_exp", "plot_mean_ci_up_nexp",
    "reverse_means_variability")

  # md
  cols_med = c(    # quarts
    "reverse_med", "min_exp", "q1_exp", "med_exp", "q3_exp", "max_exp",
    "min_nexp", "q1_nexp", "med_nexp", "q3_nexp", "max_nexp"
  )
  inf_med = c(
    "whether the direction of the effect size generated from the medians/quartiles/ranges should be flipped - logical",
    "minimum value of the experimental/exposed group - numeric",
    "first quartile of the experimental/exposed group - numeric",
    "median value of the experimental/exposed group - numeric",
    "third quartile of the experimental/exposed group - numeric",
    "maximum value of the experimental/exposed group - numeric",
    "minimum value of the non-experimental/non-exposed group - numeric",
    "first quartile of the non-experimental/non-exposed group - numeric",
    "median value of the non-experimental/non-exposed group - numeric",
    "third quartile of the non-experimental/non-exposed group - numeric",
    "maximum value of the non-experimental/non-exposed  group - numeric"
  )
  # means post
  inf_means_post = c("whether the direction of the effect size generated from the means should be flipped - logical",
    "mean of participants in the experimental/exposed group - numeric",
    "standard deviation of participants in the experimental/exposed group - numeric",
    "standard error of participants in the experimental/exposed group - numeric",
    "lower bound of the 95% CI of the mean of the experimental/exposed group - numeric",
    "upper bound of the 95% CI of the mean of the experimental/exposed group - numeric",
    "mean of participants in the non-experimental/non-exposed group - numeric",
    "standard deviation of participants in the non-experimental/non-exposed group - numeric",
    "standard error of participants in the non-experimental/non-exposed group - numeric",
    "lower bound of the 95% CI of the mean of the non-experimental/non-exposed group - numeric",
    "upper bound of the 95% CI of the mean of the non-experimental/non-exposed group - numeric",
    "pooled standard deviation across both groups - numeric",

    # plot
    "whether the direction of the effect size generated from the means extracted from a plot should be flipped - logical",

    "mean of participants in the experimental/exposed group (extracted from a plot) - numeric",
    "mean of participants in the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting -1 SD from the mean of the experimental/exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting -1 SD from the mean of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting +1 SD from the mean of the experimental/exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting +1 SD from the mean of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting -1 SE from the mean of the experimental/exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting -1 SE from the mean of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting +1 SE from the mean of the experimental/exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting +1 SE from the mean of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting the 95% CI of the mean of the experimental/exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting the 95% CI of the mean of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting the 95% CI of the mean of the experimental/exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting the 95% CI of the mean of the non-experimental/non-exposed group (extracted from a plot) - numeric",

    "whether the direction of the effect size generated from the means variability should be flipped - logical")

  cols_md =  c(
    "reverse_md", "md", "md_sd", "md_se", "md_ci_lo", "md_ci_up",
   "md_pval")

  inf_md = c("whether the direction of the effect size generated from the mean difference should be flipped - logical",
  "mean difference - numeric",
  "standard deviation of the mean difference - numeric",
  "standard error of the mean difference - numeric",
  "lower bound of the 95% CI of the mean difference - numeric",
  "upper bound of the 95% CI of the mean difference - numeric",
  "p-value of the mean difference - numeric")

  cols_anova = c(
    # smd/g
    "reverse_d", "cohen_d", "reverse_g", "hedges_g",
    # t f eta
    "reverse_student_t", "student_t",
    "reverse_student_t_pval", "student_t_pval",
    "reverse_anova_f", "anova_f",
    "reverse_anova_f_pval", "anova_f_pval",
    "reverse_etasq", "etasq",
    "reverse_pt_bis_r", "pt_bis_r", "reverse_pt_bis_r_pval", "pt_bis_r_pval")

  inf_anova = c(
    # smd/g
    "whether the direction of the cohen_d value should be flipped - logical",
    "value of the Cohen's d - numeric",
    "whether the direction of the hedges_g value should be flipped - logical",
    "value of the Hedges' g - numeric",

    # t f
    "whether the direction of the effect size generated from the Student's t-test should be flipped - logical",
    "Student's t-test value - numeric",
    "whether the direction of the effect size generated from the Student's t-test p-value should be flipped - logical",
    "two-tailed p-value of a Student's t-test - numeric",
    "whether the direction of the effect size generated from the ANOVA F should be flipped - logical",
    "ANOVA F-test value - numeric",
    "whether the direction of the effect size generated from the ANOVA p-value should be flipped - logical",
    "two-tailed p-value of an ANOVA - numeric",
    "whether the direction of the effect size generated from the eta-squared should be flipped - logical",
    "eta-squared value - numeric",
    "whether the direction of the effect size generated from the point-biserial correlation should be flipped - logical",
    "point-biserial correlation coefficient value - numeric",
    "whether the direction of the effect size generated from the point-biserial correlation p-value should be flipped - logical",
    "p-value of a point-biserial correlation - numeric")

  cols_regression = c("reverse_beta_std", "beta_std", "reverse_beta_unstd", "beta_unstd", "sd_dv")
  inf_regression = c("whether the direction of the effect size generated from the standardized beta should be flipped - logical",
    "a standardized regression coefficient value (the predictor of interest must be binary) - numeric",
    "whether the direction of the effect size generated from the unstandardized beta should be flipped - logical",
    "a non-standardized regression coefficient value (the predictor of interest must be binary) - numeric",
    "standard deviation of the dependent variable (not the standard deviation of the regression coefficient) - numeric")

  cols_linreg = c("reverse_linreg_t", "linreg_t", "n_covariates")
  inf_linreg = c(
    "whether the direction of the effect size generated from the regression t should be flipped - logical",
    "t-statistic from a linear regression model - numeric",
    "number of covariates in the regression model (excluding predictor of interest) - numeric")

  cols_linreg_b = c("reverse_linreg_b", "linreg_b", "linreg_b_se",
                    "linreg_b_ci_lo", "linreg_b_ci_up",
                    "reverse_linreg_b_pval", "linreg_b_pval")
  inf_linreg_b = c(
    "whether the direction of the effect size generated from the regression coefficient should be flipped - logical",
    "unstandardized regression coefficient from a linear regression model - numeric",
    "standard error of the regression coefficient - numeric",
    "lower bound of the 95% CI of the regression coefficient - numeric",
    "upper bound of the 95% CI of the regression coefficient - numeric",
    "whether the direction of the effect size generated from the regression coefficient p-value should be flipped - logical",
    "p-value of the regression coefficient - numeric")

  cols_pre_means = c("reverse_means_pre_post", "mean_pre_exp", "mean_pre_sd_exp", "mean_pre_se_exp", "mean_pre_se_nexp", "mean_pre_ci_lo_exp", "mean_pre_ci_up_exp",
                     "mean_pre_nexp", "mean_pre_sd_nexp", "mean_pre_ci_lo_nexp", "mean_pre_ci_up_nexp",
                     "reverse_mean_change", "mean_change_exp", "mean_change_nexp",
                     "mean_change_sd_exp", "mean_change_sd_nexp",
                     "mean_change_se_exp", "mean_change_se_nexp",
                     "mean_change_ci_lo_exp", "mean_change_ci_up_exp",
                     "mean_change_ci_lo_nexp", "mean_change_ci_up_nexp",
                     "mean_change_pval_exp", "mean_change_pval_nexp",
                     "r_pre_post_exp", "r_pre_post_nexp")

  inf_pre_means = c("whether the direction of the effect size generated from the pre/post means should be flipped - logical",
  "mean of the experimental/exposed group at baseline - numeric",
  "standard deviation of the experimental/exposed group at baseline - numeric",
  "standard error of the experimental/exposed group at baseline - numeric",
  "lower bound of the 95% CI of the mean of the experimental/exposed group at baseline - numeric",
  "upper bound of the 95% CI of the mean of the experimental/exposed group at baseline - numeric",
  "mean of the non-experimental/non-exposed group at baseline - numeric",
  "standard deviation of the non-experimental/non-exposed group at baseline - numeric",
  "standard error of the non-experimental/non-exposed group at baseline - numeric",
  "lower bound of the 95% CI of the mean of the non-experimental/non-exposed group at baseline - numeric",
  "upper bound of the 95% CI of the mean of the non-experimental/non-exposed group at baseline - numeric",

  "whether the direction of the effect size generated from the mean change should be flipped - logical",
  "mean change of participants in the experimental/exposed group - numeric",
  "mean change of participants in the non-experimental/non-exposed group - numeric",
  "standard deviation of the mean change for participants in the experimental/exposed group - numeric",
  "standard deviation of the mean change for participants in the non-experimental/non-exposed group - numeric",
  "standard error of the mean change for participants in the experimental/exposed group - numeric",
  "standard error of the mean change for participants in the non-experimental/non-exposed group - numeric",
  "lower bound of the 95% CI of the mean change for the experimental/exposed group - numeric",
  "upper bound of the 95% CI of the mean change for the experimental/exposed group - numeric",
  "lower bound of the 95% CI of the mean change for the non-experimental/non-exposed group - numeric",
  "upper bound of the 95% CI of the mean change for the non-experimental/non-exposed group - numeric",

  "p-value of the mean change for the experimental/exposed group - numeric",
  "p-value of the mean change for the non-experimental/non-exposed group - numeric",

  "pre-post correlation in the experimental/exposed group - numeric",
  "pre-post correlation in the non-experimental/non-exposed group - numeric")

  cols_paired_statistics = c("reverse_paired_t", "paired_t_exp", "paired_t_nexp",
                             "reverse_paired_t_pval", "paired_t_pval_exp", "paired_t_pval_nexp",
                             "reverse_paired_f", "paired_f_exp", "paired_f_nexp",
                             "reverse_paired_f_pval",  "paired_f_pval_exp", "paired_f_pval_nexp")
  inf_paired_statistics = c(
  "whether the direction of the effect size generated from the paired t-tests should be flipped - logical",
  "paired t-test value of the experimental/exposed group - numeric",
  "paired t-test value of the non-experimental/non-exposed group - numeric",

  "whether the direction of the effect size generated from the paired t-test p-values should be flipped - logical",
  "p-value of the paired t-test value of the experimental/exposed group - numeric",
  "p-value of the paired t-test value of the non-experimental/non-exposed group - numeric",

  "whether the direction of the effect size generated from the paired F-tests should be flipped - logical",
  "paired ANOVA F value of the experimental/exposed group - numeric",
  "paired ANOVA F value of the non-experimental/non-exposed group - numeric",

  "whether the direction of the effect size generated from the paired F-tests p-values should be flipped - logical",
  "p-value of the paired ANOVA-F of the experimental/exposed group - numeric",
  "p-value of the paired ANOVA-F of the non-experimental/non-exposed group - numeric")

  cols_prop_single = c("reverse_prop", "prop", "n_cases")
  inf_prop_single = c(
    "whether the direction of the proportion should be flipped - logical",
    "proportion in single group (0-1) - numeric",
    "number of cases (alternative to proportion) - numeric")

  cols_ancova_mean = c(
    "reverse_ancova_means", "ancova_mean_sd_pooled", "cov_outcome_r", "n_cov_ancova",
     "ancova_mean_exp", "ancova_mean_nexp",
     "ancova_mean_sd_exp", "ancova_mean_sd_nexp",
     "ancova_mean_se_exp", "ancova_mean_se_nexp",
     "ancova_mean_ci_lo_exp", "ancova_mean_ci_up_exp",
     "ancova_mean_ci_lo_nexp", "ancova_mean_ci_up_nexp",

     "reverse_plot_ancova_means",
     "plot_ancova_mean_exp",
     "plot_ancova_mean_sd_lo_exp", "plot_ancova_mean_sd_up_exp",
     "plot_ancova_mean_se_lo_exp", "plot_ancova_mean_se_up_exp",
     "plot_ancova_mean_ci_lo_exp",  "plot_ancova_mean_ci_up_exp",

     "plot_ancova_mean_nexp",
     "plot_ancova_mean_sd_lo_nexp", "plot_ancova_mean_sd_up_nexp",
     "plot_ancova_mean_se_lo_nexp", "plot_ancova_mean_se_up_nexp",
     "plot_ancova_mean_ci_lo_nexp", "plot_ancova_mean_ci_up_nexp")

  inf_ancova_mean = c(
    "whether the direction of the effect size generated from the ANCOVA means should be flipped - logical",
    "adjusted pooled standard deviation - numeric",
    "correlation between the outcome and covariate(s) (multiple correlation when multiple covariates are included in the ANCOVA model) - numeric",
    "number of covariates in the ANCOVA model - numeric",

    "adjusted mean of participants in the experimental/exposed group - numeric",
    "adjusted standard deviation of participants in the experimental/exposed group - numeric",
    "adjusted standard error of participants in the experimental/exposed group - numeric",
    "lower bound of the adjusted 95% CI of the mean of the experimental/exposed group - numeric",
    "upper bound of the adjusted 95% CI of the mean of the experimental/exposed group - numeric",

    "adjusted mean of participants in the non-experimental/non-exposed group - numeric",
    "adjusted standard deviation of participants in the non-experimental/non-exposed group - numeric",
    "adjusted standard error of participants in the non-experimental/non-exposed group - numeric",
    "lower bound of the adjusted 95% CI of the mean of the non-experimental/non-exposed group - numeric",
    "upper bound of the adjusted 95% CI of the mean of the non-experimental/non-exposed group - numeric",

    "whether the direction of the effect size generated from the ANCOVA means extracted from a plot should be flipped - logical",
    "adjusted mean (from ANCOVA) of participants in the experimental/exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting -1 SD from the adjusted mean (from ANCOVA) of the experimental/exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting +1 SD from the adjusted mean (from ANCOVA) of the experimental/exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting -1 SE from the adjusted mean (from ANCOVA) of the experimental/exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting +1 SE from the adjusted mean (from ANCOVA) of the experimental/exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting the 95% CI of the adjusted mean (from ANCOVA) of the experimental/exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting the 95% CI of the adjusted mean (from ANCOVA) of the experimental/exposed group (extracted from a plot) - numeric",

    "adjusted mean (from ANCOVA) of participants in the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting -1 SD from the adjusted mean (from ANCOVA) of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting +1 SD from the adjusted mean (from ANCOVA) of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting -1 SE from the adjusted mean (from ANCOVA) of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting +1 SE from the adjusted mean (from ANCOVA) of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "lower bound of an error bar depicting the 95% CI of the adjusted mean (from ANCOVA) of the non-experimental/non-exposed group (extracted from a plot) - numeric",
    "upper bound of an error bar depicting the 95% CI of the adjusted mean (from ANCOVA) of the non-experimental/non-exposed group (extracted from a plot) - numeric"
  )
    # ANCOVA
  cols_ancova_stat = c("reverse_ancova_t", "ancova_t",
    "reverse_ancova_f", "ancova_f",
    "reverse_ancova_t_pval", "ancova_t_pval",
    "reverse_ancova_f_pval", "ancova_f_pval",
    "reverse_ancova_md", "ancova_md", "ancova_md_sd",
    "ancova_md_se", "ancova_md_ci_lo", "ancova_md_ci_up",
    "ancova_md_pval",
    "cohen_d_adj",
    "etasq_adj", "reverse_etasq_adj")

  inf_ancova_stat =  c("whether the direction of the effect size generated from the ANCOVA t-test should be flipped - logical",
    "a t-statistic from an ANCOVA - numeric",
    "whether the direction of the effect size generated from the ANCOVA F-test should be flipped - logical",
    "a F-statistic from an ANCOVA - numeric",
    "whether the direction of the effect size generated from the ANCOVA t-test p-value should be flipped - logical",
    "a two-tailed p-value from an ANCOVA - numeric",
    "whether the direction of the effect size generated from the ANCOVA F-test p-value should be flipped - logical",
    "a two-tailed p-value from an ANCOVA - numeric",

    "whether the direction of the effect size generated from the ANCOVA mean difference should be flipped - logical",
    "adjusted mean difference between two independent groups - numeric",
    "standard deviation of the adjusted mean difference - numeric",
    "standard error of the adjusted mean difference - numeric",
    "lower bound of the 95% CI around the adjusted mean difference - numeric",
    "upper bound of the 95% CI around the adjusted mean difference - numeric",

    "p-value of the adjusted mean difference - numeric",

    "value of an adjusted Cohen's d - numeric",
    "value of an adjusted eta-square - numeric",
    "whether the direction of the effect size generated from the adjusted eta-square should be flipped - logical")


    # 2x2 table
    cols_2x2 = c(
    "reverse_2x2", "baseline_risk", "small_margin_prop",
    "n_cases_exp", "n_cases_nexp", "n_controls_exp", "n_controls_nexp",
    "reverse_prop", "prop_cases_exp", "prop_cases_nexp",
    "reverse_chisq", "chisq", "reverse_chisq_pval", "chisq_pval",
    "yates_chisq",
    "reverse_phi", "phi")

    inf_2x2 = c(
      "whether the direction of the effect size generated from the contingency table should be flipped - logical",
      "proportion of cases/events in the non-exposed group - numeric",
      "smallest margin proportion of cases/events in the underlying 2x2 table - numeric",
      "number of cases/events in the exposed group - numeric",
      "number of cases/events in the non exposed group - numeric",
      "number of controls/no-event in the exposed group - numeric",
      "number of controls/no-event in the non exposed group - numeric",
      "whether the direction of the cohen_d value should be flipped - logical",
      "proportion of cases/events in the exposed group (ranging from 0 to 1) - numeric",
      "proportion of cases/events in the non-exposed group (ranging from 0 to 1) - numeric",
      "whether the direction of the effect size generated from the chi-square value should be flipped - logical",
      "chi-square value - numeric",
      "whether the direction of the effect size generated from the chi-square p-value should be flipped - logical",
      "chi-square p-value - numeric",
      "whether the chi-square statistic/p-value reported in the primary was computed with Yates' continuity correction - logical (default FALSE)",
      "whether the direction of the effect size generated from the phi value should be flipped - logical",
      "phi value - numeric"
    )

    cols_or = c(
    # or
    "reverse_or", "or", "logor", "logor_se", "or_ci_lo", "or_ci_up", "logor_ci_lo", "logor_ci_up",
    "reverse_or_pval", "or_pval",
    "reverse_logreg_t", "logreg_t")

    inf_or = c(
      "whether the direction of the effect size generated from the odds ratio should be flipped - logical",
      "odds ratio value - numeric",
      "log odds ratio value - numeric",
      "standard error of the log odds ratio - numeric",
      "lower bound of the 95% CI of the risk ratio - numeric",
      "upper bound of the 95% CI of the log risk ratio - numeric",
      "lower bound of the 95% CI of the risk ratio - numeric",
      "upper bound of the 95% CI of the log risk ratio - numeric",
      "whether the direction of the cohen_d value should be flipped - logical",
      "p-value of an odds ratio - numeric",
      "whether the direction of the effect size generated from the logistic regression t should be flipped - logical",
      "Wald t-statistic from a logistic regression model - numeric"
    )
    cols_rr = c(
    # rr
    "reverse_rr", "rr", "logrr", "logrr_se",
    "rr_ci_lo", "rr_ci_up", "logrr_ci_lo", "logrr_ci_up", "reverse_rr_pval", "rr_pval")
    inf_rr = c(
      "whether the direction of the effect size generated from the risk ratio should be flipped - logical",
      "risk ratio value - numeric",
      "log risk ratio value - numeric",
      "log risk ratio standard error - numeric",
      "lower bound of the 95% CI of the risk ratio - numeric",
      "upper bound of the 95% CI of the risk ratio - numeric",
      "lower bound of the 95% CI of the log risk ratio - numeric",
      "upper bound of the 95% CI of the log risk ratio - numeric",
      "whether the direction of the risk ratio generated from a p-value should be flipped - logical",
      "p-value of a risk ratio -  numeric"
    )

    cols_rd = c(
    # rd (risk difference)
    "reverse_rd", "rd", "rd_se",
    "rd_ci_lo", "rd_ci_up",
    "reverse_rd_pval", "rd_pval")

    inf_rd = c(
      "whether the direction of the effect size generated from the risk difference should be flipped - logical",
      "risk difference value (control risk minus treatment risk) - numeric",
      "standard error of the risk difference - numeric",
      "lower bound of the 95% CI of the risk difference - numeric",
      "upper bound of the 95% CI of the risk difference - numeric",
      "whether the direction of the effect size generated from the risk difference p-value should be flipped - logical",
      "p-value of the risk difference - numeric"
    )

    cols_r = c(
    "reverse_pearson_r", "pearson_r", "reverse_fisher_z",
    "fisher_z", "reverse_spearman_r", "spearman_r",
    "unit_increase_iv", "unit_type", "sd_iv")

    inf_r = c(
      "whether the direction of the effect size generated from the Pearson's correlation should be flipped - logical",
      "Pearson's correlation coefficient value  - numeric",
      "whether the direction of the effect size generated from the Fisher's z should be flipped - logical",
      "Fisher's r-to-z transformed correlation coefficient - numeric",
      "whether the direction of the effect size generated from the Spearman's correlation should be flipped - logical",
      "Spearman's rank correlation coefficient - numeric",
      "a value of the independent variable that will be used to estimate the Cohen's d - numeric",
      "type of unit for the unit_increase_iv variable. Must be 'sd' or 'raw_scale' ('value' and 'raw_data' are synonyms of 'raw_scale') - character",
      "standard deviation of the independent variable - numeric"
      )

    cols_alpha = c("cronbach_alpha", "n_items")
    inf_alpha = c(
      "Cronbach's alpha reliability coefficient - numeric",
      "number of items in the scale - numeric")

    cols_omega = c("omega", "omega_se", "omega_ci_lo", "omega_ci_up", "omega_type",
                   "omega_estimator")
    inf_omega = c(
      "McDonald's omega reliability coefficient (natural scale) - numeric",
      "standard error of omega, ON THE NATURAL SCALE - numeric",
      "lower bound of the 95% CI of omega (natural scale) - numeric",
      "upper bound of the 95% CI of omega (natural scale) - numeric",
      "which omega: 'total' (default), 'hierarchical', 'asymptotic' or 'subscale' - character",
      "how omega was estimated: 'cfa_bifactor', 'cfa_1factor', 'efa_schmid_leiman', 'first_pc', 'first_pf' or 'unspecified' - character")

    cols_icc = c("icc", "n_measurements", "icc_type")
    inf_icc = c(
      "intraclass correlation coefficient - numeric",
      "number of measurements or raters - numeric",
      "ICC type: 'agreement' (ICC(2,1)) or 'consistency' (ICC(3,1)) - character")

    cols_irr = c(
      # survival
    "reverse_irr", "time_exp", "time_nexp", "baseline_rate")

    inf_irr = c(
      "whether the direction of the effect size generated from the number of cases & times should be flipped - logical",
      "person-time of disease-free observation in the exposed group - numeric",
      "person-time of disease-free observation in the non-exposed group - numeric",
      "incidence rate of events in the non-exposed/control group (events per person-time) - numeric"
    )


    if (measure %in% c("d", "g")) {
      if (type_of_measure == "natural") {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_means_post,
                             inf_md, inf_anova, inf_regression,
                             inf_pre_means, inf_paired_statistics,
                             inf_ancova_mean, inf_ancova_stat, inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_means_post,
                           cols_md, cols_anova, cols_regression,
                           cols_pre_means, cols_paired_statistics,
                           cols_ancova_mean, cols_ancova_stat, cols_input_crude, cols_input_adjusted)

      } else {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_means_post,
                             inf_md, inf_anova, inf_regression,
                             inf_pre_means, inf_med,
                             inf_paired_statistics,
                             inf_ancova_mean, inf_ancova_stat,
                             inf_or, inf_2x2, inf_cases,
                             inf_sample, inf_r,
                             inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_means_post,
                           cols_md, cols_anova, cols_regression,
                           cols_pre_means, cols_med,
                           cols_paired_statistics,
                           cols_ancova_mean, cols_ancova_stat,
                           cols_or ,cols_2x2, cols_cases,
                           cols_sample, cols_r,
                           cols_input_crude, cols_input_adjusted)

      }
    } else if (measure == "md") {
      if (type_of_measure == "natural") {

      dat = data.frame(t(c(inf_req, inf_exposed, inf_md, inf_means_post,
                           inf_pre_means, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_exposed, cols_md, cols_means_post,
                         cols_pre_means, cols_input_crude, cols_input_adjusted)
      } else {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_md, inf_means_post,
                             inf_pre_means, inf_med, inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_md, cols_means_post,
                           cols_pre_means, cols_med, cols_input_crude, cols_input_adjusted)

      }
    } else if (measure %in% c("dw", "gw", "mdw")) {
      if (type_of_measure == "natural") {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_pre_means, inf_paired_statistics,
                             inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_pre_means, cols_paired_statistics,
                           cols_input_crude, cols_input_adjusted)
      } else {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_pre_means, inf_paired_statistics,
                             inf_means_post, inf_md, inf_anova, inf_regression,
                             inf_ancova_mean, inf_ancova_stat, inf_med,
                             inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_pre_means, cols_paired_statistics,
                           cols_means_post, cols_md, cols_anova, cols_regression,
                           cols_ancova_mean, cols_ancova_stat, cols_med,
                           cols_input_crude, cols_input_adjusted)
      }
    } else if (measure %in% c("or")) {
      if (type_of_measure == "natural") {
      dat = data.frame(t(c(inf_req, inf_exposed, inf_cases, inf_2x2,
                           inf_or, inf_rr, inf_rd,
                           inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_exposed, cols_cases, cols_2x2,
                         cols_or, cols_rr, cols_rd,
                         cols_input_crude, cols_input_adjusted)
      } else {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_cases, inf_2x2,
                             inf_or, inf_rr, inf_rd,
                             inf_means_post,
                             inf_md, inf_anova, inf_regression,
                             inf_pre_means, inf_paired_statistics,
                             inf_ancova_mean, inf_ancova_stat, inf_med,
                             inf_sample, inf_r,
                             inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_cases, cols_2x2,
                           cols_or, cols_rr, cols_rd,
                           cols_means_post,
                           cols_md, cols_anova, cols_regression,
                           cols_pre_means, cols_paired_statistics,
                           cols_ancova_mean, cols_ancova_stat, cols_med,
                           cols_sample, cols_r,
                           cols_input_crude, cols_input_adjusted)

      }
    } else if (measure %in% c("rr", "nnt", "rd")) {
      if (type_of_measure == "natural") {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_cases, inf_2x2,
                             inf_rr, inf_rd, inf_irr,
                             inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_cases, cols_2x2,
                           cols_rr, cols_rd, cols_irr,
                           cols_input_crude, cols_input_adjusted)
      } else {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_cases, inf_2x2,
                             inf_or, inf_rr, inf_rd, inf_irr,
                             inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_cases, cols_2x2,
                           cols_or, cols_rr, cols_rd, cols_irr,
                           cols_input_crude, cols_input_adjusted)

      }
    } else if (measure %in% c("r", "z")) {
      if (type_of_measure == "natural") {

      dat = data.frame(t(c(inf_req, inf_sample, inf_r, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_sample, cols_r, cols_input_crude, cols_input_adjusted)
      } else {
        dat = data.frame(t(c(inf_req, inf_sample, inf_r,
                             inf_cases, inf_2x2,
                             inf_or,
                             inf_means_post,
                             inf_md, inf_anova, inf_regression,
                             inf_pre_means, inf_paired_statistics,
                             inf_ancova_mean, inf_ancova_stat, inf_med,
                             inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_sample, cols_r,
                           cols_cases, cols_2x2,
                           cols_or,
                           cols_means_post,
                           cols_md, cols_anova, cols_regression,
                           cols_pre_means, cols_paired_statistics,
                           cols_ancova_mean, cols_ancova_stat, cols_med,
                           cols_input_crude, cols_input_adjusted)

      }
    } else if (measure %in% c("rp", "zp")) {

      dat = data.frame(t(c(inf_req, inf_sample, inf_linreg, inf_linreg_b, inf_r, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_sample, cols_linreg, cols_linreg_b, cols_r, cols_input_crude, cols_input_adjusted)

    } else if (measure %in% c("logvr", "logcvr")) {
      if (type_of_measure == "natural") {

      dat = data.frame(t(c(inf_req, inf_exposed, inf_means_post, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_exposed, cols_means_post, cols_input_crude, cols_input_adjusted)
      } else {
        dat = data.frame(t(c(inf_req, inf_exposed, inf_means_post,
                             inf_med,
                             inf_input_crude, inf_input_adjusted)))
        colnames(dat) <- c(cols_req, cols_exposed, cols_means_post,
                           cols_med,
                           cols_input_crude, cols_input_adjusted)

      }
    } else if (measure %in% c("irr")) {

      dat = data.frame(t(c(inf_req, inf_exposed, inf_irr, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_exposed, cols_irr, cols_input_crude, cols_input_adjusted)

    } else if (measure %in% c("hr")) {

      # HR can only enter via user_input_crude / user_input_adj (no raw-data
      # methods, because HR requires per-subject time-to-event data).
      dat = data.frame(t(c(inf_req, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_input_crude, cols_input_adjusted)

    } else if (measure %in% c("prop")) {

      dat = data.frame(t(c(inf_req, inf_sample, inf_prop_single, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_sample, cols_prop_single, cols_input_crude, cols_input_adjusted)

    } else if (measure == "alpha") {

      dat = data.frame(t(c(inf_req, inf_sample, inf_alpha, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_sample, cols_alpha, cols_input_crude, cols_input_adjusted)

    } else if (measure == "omega") {

      dat = data.frame(t(c(inf_req, inf_sample, inf_omega, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_sample, cols_omega, cols_input_crude, cols_input_adjusted)

    } else if (measure == "icc") {

      dat = data.frame(t(c(inf_req, inf_sample, inf_icc, inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_sample, cols_icc, cols_input_crude, cols_input_adjusted)

    } else if (measure %in% c("all")) {
      dat = data.frame(t(c(inf_req, inf_exposed, inf_means_post,
                           inf_md, inf_anova, inf_regression,
                           inf_linreg, inf_linreg_b,
                           inf_med,
                           inf_pre_means, inf_paired_statistics,
                           inf_ancova_mean, inf_ancova_stat,
                           inf_cases, inf_2x2,
                           inf_or, inf_rr, inf_rd,
                           inf_sample, inf_r,
                           inf_irr,
                           inf_prop_single,
                           inf_alpha, inf_omega, inf_icc,
                           inf_input_crude, inf_input_adjusted)))
      colnames(dat) <- c(cols_req, cols_exposed, cols_means_post,
                         cols_md, cols_anova, cols_regression,
                         cols_linreg, cols_linreg_b,
                         cols_med,
                         cols_pre_means, cols_paired_statistics,
                         cols_ancova_mean, cols_ancova_stat,
                         cols_cases, cols_2x2,
                         cols_or, cols_rr, cols_rd,
                         cols_sample, cols_r,
                         cols_irr,
                         cols_prop_single,
                         cols_alpha, cols_omega, cols_icc,
                         cols_input_crude, cols_input_adjusted)

    }

  if (missing(name)) name = "data_extract_metaConvert"
  if (extension == "data.frame") {
    dat
  } else if (extension != ".xlsx") {
    if (verbose) {
      message(paste0("Sheet created at location:\n\n", getwd(), "/",
                       name, extension,
                       "\n\n Good luck with your research!"))
    }
    rio::export(dat, paste0(name, extension))
  } else {
    if (verbose) {
      message(paste0("Sheet created at location:\n\n", getwd(), "/",
                       name, extension,
                       "\n\n Good luck with your research!"))
    }
    rio::export(dat, paste0(name, extension), overwrite = TRUE)
  }
}
#' Overview of effect size measures generated from each type of input data
#'
#' @param name Name of the file created
#' @param measure Target effect size measure. One of "all", "d", "g", "md", "or", "rr",
#'   "nnt", "r", "z", "irr", "logvr" or "logcvr". Default is "all". The risk difference and
#'   the psychometric measures (rd, alpha, icc, prop, hr) are not covered by this table;
#'   use \code{\link{data_extraction_sheet}()} for those.
#' @param extension Extension of the file created. Most common are ".xlsx", ".csv" or ".txt". It is also possible to generate an R dataframe object by using the "data.frame" extension.
#' @param type_of_measure One of "natural+converted" or "natural" (see details).
#' @param verbose logical variable indicating whether some information should be printed (e.g., the location where the sheet is created when using ".xlsx", ".csv" or ".txt" extensions)
#'
#' @details
#' This function generates, on your computer on in the console,
#' a dataset showing each effect size measure computed from each type of input data.
#' The exact combination and names of input data required are available in the links.
#'
#' The \code{measure} argument allows to filter the dataset created.
#' Only the input data allowing to estimate the selected effect size measure will be shown. Default is "all".

#' The \code{type_of_measure} argument allows to filter the dataset created.
#' - If "natural+converted" is selected, the dataset will contain all input data allowing to naturally estimate
#' and to convert the selected effect size measure
#' - If "natural" is selected, the dataset will contain all input data allowing to naturally estimate
#' the selected effect size measure
#'
#' ## Extension
#' You can export a file in various formats outside R
#' (by indicating, for example, ".txt", ".xlsx", or ".csv") in the \code{extension} argument.
#' You can also visualise this dataset directly in R by setting
#' \code{extension = "data.frame"}.
#'
#' This table is designed to be used in combination with tables showing the combination of input data
#' leading to estimate each of the effect size measures (https://metaconvert.org/input.html)
#'
#' @return
#' This function returns a table dataset presenting the input data enabling to
#' compute each effect size measure.
#'
#' @md
#'
#' @export see_input_data
#'
#' @examples
#' see_input_data(measure = "md", extension = "data.frame")

see_input_data <- function(measure = c("all", "d", "g", "md", "or", "rr", "nnt",
                                            "r", "z", "logvr", "logcvr", "irr"),
                          type_of_measure = c("natural", "natural+converted"),
                          name = "mcv_input_data",
                          extension = c("data.frame", ".txt", ".csv", ".xlsx"),
                          verbose = TRUE) {

  type_of_measure = type_of_measure[1]
  extension = extension[1]
  measure = measure[1]

  # The switch() blocks below have no default, so an unmatched value silently produced a
  # NULL 'dat' and the subsequent order() failed with "argument 1 is not a vector". The
  # measures this table does not cover are rejected here, by name, instead.
  supported_measures = c("all", "d", "g", "md", "or", "rr", "nnt",
                         "r", "z", "Z", "logvr", "logcvr", "irr")
  if (!measure %in% supported_measures) {
    stop("see_input_data(): unsupported 'measure' value \"", measure, "\". Supported ",
         "values are ", paste0("\"", setdiff(supported_measures, "Z"), "\"", collapse = ", "),
         ". The risk difference and the psychometric measures (\"rd\", \"alpha\", \"icc\", ",
         "\"prop\", \"hr\") are not covered by this table; use data_extraction_sheet() ",
         "for those.", call. = FALSE)
  }

  hier_list = convert_df(df.haza[1, ], measure = "d", verbose=FALSE)
  hier = as.character(sapply(hier_list, function(x) unique(x[, "info_used"])))

  dat = data.frame(
    hierarch_name = hier,
    name_input_data = rep(NA, length(hier)),
    list_input_data = rep(NA, length(hier)),
    corresponding_R_function = rep(NA, length(hier)),
    natural_effect_size_measure = rep(NA, length(hier)),
    converted_effect_size_measure = rep(NA, length(hier)),
    adjusted_input_data = rep(NA, length(hier))
  )

  dat[dat$hierarch_name == "user_input_crude",
      2:7] <- c("Crude effect size measure indicated by the user",
                "Section 23. https://metaconvert.org/input.html",
                "es_from_user_crude()",
                "Any effect size measure",
                "N/A",
                "Non-adjusted")
  dat[dat$hierarch_name == "user_input_adj",
      2:7] <- c("Adjusted effect size measure indicated by the user",
                "Section 24. https://metaconvert.org/input.html",
                "es_from_user_adj()",
                "Any adjusted effect size measure", "N/A",
                "Adjusted")
  dat[dat$hierarch_name == "cohen_d",
      2:7] <- c("Cohen's d value",
                "Section 1. https://metaconvert.org/input.html",
                "es_from_cohen_d()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "cohen_d_adj",
      2:7] <- c("Section 17. Adjusted Cohen's d value",
                "https://metaconvert.org/input.html",
                "es_from_cohen_d_adj()",
                "D+G", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "hedges_g",
      2:7] <- c("Hedges'g value",
                "Section 1. https://metaconvert.org/input.html",
                "es_from_hedges_g()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "or",
      2:7] <- c("Odds ratio value",
                "Section 2. https://metaconvert.org/input.html",
                "es_from_or()",
                "N/A", "OR+RR+NNT+D+G+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "or_se",
      2:7] <- c("Odds ratio value + standard error",
                "Section 2. https://metaconvert.org/input.html",
                "es_from_or_se()",
                "OR", "RR+NNT+D+G+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "or_ci",
      2:7] <- c("Odds ratio value + 95% confidence interval",
                "Section 2. https://metaconvert.org/input.html",
                "es_from_or_ci()",
                "OR", "RR+NNT+D+G+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "or_pval",
      2:7] <- c("Odds ratio value + p-value",
                "Section 2. https://metaconvert.org/input.html",
                "es_from_or_pval()",
                "OR", "RR+NNT+D+G+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "logreg_t",
      2:7] <- c("OR or RR + Wald t-statistic from logistic/log-binomial regression",
                "Section 2. https://metaconvert.org/input.html",
                "es_from_logreg_t()",
                "OR+RR", "NNT+D+G+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "pearson_r",
      2:7] <- c("Pearson's correlation coefficient",
                "Section 4. https://metaconvert.org/input.html",
                "es_from_pearson_r()",
                "R+Z", "D+G+OR",
                "Non-adjusted")
  dat[dat$hierarch_name == "fisher_z",
      2:7] <- c("Fisher's r-to-z correlation coefficient",
                "Section 4. https://metaconvert.org/input.html",
                "es_from_fisher_z()",
                "R+Z", "D+G+OR",
                "Non-adjusted")
  dat[dat$hierarch_name == "spearman_r",
      2:7] <- c("Spearman's rank correlation coefficient",
                "Section 4. https://metaconvert.org/input.html",
                "es_from_spearman_rho()",
                "R+Z", "D+G+OR",
                "Non-adjusted")

  dat[dat$hierarch_name == "means_sd",
      2:7] <- c("Means and standard devations of two independent groups",
                "Section 9. https://metaconvert.org/input.html",
                "es_from_means_sd()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "means_sd_pooled",
      2:7] <- c("Means and pooled standard devation of two independent groups",
                "Section 9. https://metaconvert.org/input.html",
                "es_from_means_sd_pooled()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "means_se",
      2:7] <- c("Means and standard errors of two independent groups",
                "Section 9. https://metaconvert.org/input.html",
                "es_from_means_se()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "means_ci",
      2:7] <- c("Means and 95% confidence intervals of two independent groups",
                "Section 9. https://metaconvert.org/input.html",
                "es_from_means_ci()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")


  dat[dat$hierarch_name == "mean_change_sd",
      2:7] <- c("Mean changes from pre- to post-test and standard devations of two independent groups",
                "Section 15. https://metaconvert.org/input.html",
                "es_from_mean_change_sd()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "mean_change_se",
      2:7] <- c("Mean changes from pre- to post-test and standard errors of two independent groups",
                "Section 15. https://metaconvert.org/input.html",
                "es_from_mean_change_se()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "mean_change_ci",
      2:7] <- c("Mean changes from pre- to post-test and 95% confidence intervals of two independent groups",
                "Section 15. https://metaconvert.org/input.html",
                "es_from_mean_change_ci()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "mean_change_pval",
      2:7] <- c("Mean changes from pre- to post-test and p-values of two independent groups",
                "Section 15. https://metaconvert.org/input.html",
                "es_from_mean_change_pval()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "means_sd_pre_post",
      2:7] <- c("Means and standard devations at pre- and post-test of two independent groups",
                "Section 14. https://metaconvert.org/input.html",
                "es_from_means_sd_pre_post()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "means_se_pre_post",
      2:7] <- c("Means and standard errors at pre- and post-test of two independent groups",
                "Section 14. https://metaconvert.org/input.html",
                "es_from_means_se_pre_post()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "means_ci_pre_post",
      2:7] <- c("Means and 95% confidence intervals at pre- and post-test of two independent groups",
                "Section 14. https://metaconvert.org/input.html",
                "es_from_means_ci_pre_post()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "paired_t",
      2:7] <- c("Paired t-tests of two independent groups",
                "Section 16. https://metaconvert.org/input.html",
                "es_from_paired_t()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "paired_t_pval",
      2:7] <- c("Paired t-test p-values of two independent groups",
                "Section 16. https://metaconvert.org/input.html",
                "es_from_paired_t_pval()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "paired_f",
      2:7] <- c("Paired F-tests of two independent groups",
                "Section 16. https://metaconvert.org/input.html",
                "es_from_paired_f()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "paired_f_pval",
      2:7] <- c("Paired F-test p-values of two independent groups",
                "Section 16. https://metaconvert.org/input.html",
                "es_from_paired_f_pval()",
                "D+G", "OR+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "student_t",
      2:7] <- c("Student t-test value",
                "Section 11. https://metaconvert.org/input.html",
                "es_from_student_t()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "student_t_pval",
      2:7] <- c("Student t-test p-values of two independent groups",
                "Section 11. https://metaconvert.org/input.html",
                "es_from_student_t_pval()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "anova_f",
      2:7] <- c("F-value from a one-way ANOVA",
                "Section 11. https://metaconvert.org/input.html",
                "es_from_anova_f()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "anova_f_pval",
      2:7] <- c("p-value from a one-way ANOVA",
                "Section 11. https://metaconvert.org/input.html",
                "es_from_anova_pval()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "etasq",
      2:7] <- c("Eta-squared value from an ANOVA model",
                "Section 11. https://metaconvert.org/input.html",
                "es_from_etasq",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "pt_bis_r",
      2:7] <- c("Point-biserial correlation coefficient value",
                "Section 11. https://metaconvert.org/input.html",
                "es_from_pt_bis_r",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "pt_bis_r_pval",
      2:7] <- c("P-value of a point-biserial correlation",
                "Section 11. https://metaconvert.org/input.html",
                "es_from_pt_bis_r_pval",
                "D+G", "OR+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "md_sd",
      2:7] <- c("Mean difference value between two independent groups and standard deviation",
                "Section 10. https://metaconvert.org/input.html",
                "es_from_md_sd()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "md_se",
      2:7] <- c("Mean difference value between two independent groups and standard error",
                "Section 10. https://metaconvert.org/input.html",
                "es_from_md_se()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "md_ci",
      2:7] <- c("Mean difference value between two independent groups and 95% confidence interval",
                "Section 10. https://metaconvert.org/input.html",
                "es_from_md_ci()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "md_pval",
      2:7] <- c("Mean difference value between two independent groups and p-value",
                "Section 10. https://metaconvert.org/input.html",
                "es_from_md_pval()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "med_quarts",
      2:7] <- c("Medians and quartiles of two independent groups",
                "Section 12. https://metaconvert.org/input.html",
                "es_from_med_quarts()",
                "N/A", "D+G+MD+OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "med_min_max",
      2:7] <- c("Medians and minimum+maximum values of two independent groups",
                "Section 12. https://metaconvert.org/input.html",
                "es_from_med_min_max()",
                "N/A", "D+G+MD+OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "med_min_max_quarts",
      2:7] <- c("Medians, quartiles and minimum+maximum values of two independent groups",
                "Section 12. https://metaconvert.org/input.html",
                "es_from_med_min_max_quarts()",
                "N/A", "D+G+MD+OR+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "means_plot",
      2:7] <- c("Means and bounds of standard deviations, standard errors, and/or 95% CIs of two independent groups extracted from a plot",
                "Section 21. https://metaconvert.org/input.html",
                "es_from_plot_means()",
                "D+G+MD", "OR+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "ancova_means_plot",
      2:7] <- c("Adjusted means and bounds of standard deviations, standard errors, and/or 95% CIs of two independent groups extracted from a plot",
                "Section 22. https://metaconvert.org/input.html",
                "es_from_plot_ancova_means()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")


  dat[dat$hierarch_name == "ancova_means_sd",
      2:7] <- c("Adjusted means and standard devations (from an ANCOVA) of two independent groups",
                "Section 19. https://metaconvert.org/input.html",
                "es_from_ancova_means_sd()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_means_se",
      2:7] <- c("Adjusted means and standard errors (from an ANCOVA) of two independent groups",
                "Section 19. https://metaconvert.org/input.html",
                "es_from_ancova_means_se()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_means_ci",
      2:7] <- c("Adjusted means and 95% confidence intervals (from an ANCOVA) of two independent groups",
                "Section 19. https://metaconvert.org/input.html",
                "es_from_ancova_means_ci()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_means_sd_pooled",
      2:7] <- c("Adjusted means (from an ANCOVA) and crude pooled standard deviation of two independent groups",
                "Section 19. https://metaconvert.org/input.html",
                "es_from_ancova_means_sd_pooled_crude()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_means_sd_pooled_adj",
      2:7] <- c("Adjusted means and pooled standard deviation (from an ANCOVA) of two independent groups",
                "Section 19. https://metaconvert.org/input.html",
                "es_from_ancova_means_sd_pooled_adj()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")

  dat[dat$hierarch_name == "ancova_t",
      2:7] <- c("T-test value from an ANCOVA model",
                "Section 18. https://metaconvert.org/input.html",
                "es_from_ancova_t()",
                "D+G", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_f",
      2:7] <- c("F-value from an ANCOVA model",
                "Section 18. https://metaconvert.org/input.html",
                "es_from_ancova_f()",
                "D+G", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_t_pval",
      2:7] <- c("p-value of a t-test from an ANCOVA model",
                "Section 18. https://metaconvert.org/input.html",
                "es_from_ancova_t_pval()",
                "D+G", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_f_pval",
      2:7] <- c("p-value of an F-test from an ANCOVA model",
                "Section 18. https://metaconvert.org/input.html",
                "es_from_ancova_f_pval()",
                "D+G", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "etasq_adj",
      2:7] <- c("Eta-squared value from an ANOVA model",
                "https://metaconvert.org/input.html#ancova_statistics",
                "es_from_etasq_adj()",
                "D+G", "OR+R+Z",
                "Adjusted")

  dat[dat$hierarch_name == "ancova_md_sd",
      2:7] <- c("Adjusted mean difference value between two independent groups and standard deviation (from an ANCOVA model)",
                "Section 20. https://metaconvert.org/input.html",
                "es_from_ancova_md_sd()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_md_se",
      2:7] <- c("Adjusted mean difference value between two independent groups and standard error (from an ANCOVA model)",
                "Section 20. https://metaconvert.org/input.html",
                "es_from_ancova_md_se()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_md_ci",
      2:7] <- c("Adjusted mean difference value between two independent groups and 95% confidence interval (from an ANCOVA model)",
                "Section 20. https://metaconvert.org/input.html",
                "es_from_ancova_md_ci()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")
  dat[dat$hierarch_name == "ancova_md_pval",
      2:7] <- c("Adjusted mean difference value between two independent groups and p-value (from an ANCOVA model)",
                "Section 20. https://metaconvert.org/input.html",
                "es_from_ancova_md_pval()",
                "D+G+MD", "OR+R+Z",
                "Adjusted")

  dat[dat$hierarch_name == "chisq",
      2:7] <- c("Chi-square value (obtained from a 2x2 table)",
                "Section 8. https://metaconvert.org/input.html",
                "es_from_chisq()",
                "OR+RR+R+Z", "D+G",
                "Non-adjusted")
  dat[dat$hierarch_name == "chisq_pval",
      2:7] <- c("P-value of a Chi-square (obtained from a 2x2 table)",
                "Section 8. https://metaconvert.org/input.html",
                "es_from_chisq_pval()",
                "OR+RR+R+Z", "D+G",
                "Non-adjusted")
  dat[dat$hierarch_name == "phi",
      2:7] <- c("Phi coefficient",
                "Section 8. https://metaconvert.org/input.html",
                "es_from_phi()",
                "OR+RR+R+Z", "D+G",
                "Non-adjusted")


  dat[dat$hierarch_name == "2x2",
      2:7] <- c("Number of participants in each cell of a 2x2 table",
                "Section 7. https://metaconvert.org/input.html",
                "es_from_2x2()",
                "OR+RR+NNT", "D+G+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "2x2_sum",
      2:7] <- c("Number of participants 'developping the disease' and row marginal sums",
                "Section 7. https://metaconvert.org/input.html",
                "es_from_2x2_sum()",
                "OR+RR+NNT", "D+G+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "2x2_prop",
      2:7] <- c("Proportion of participants 'developping the disease' and row marginal proportions",
                "Section 7. https://metaconvert.org/input.html",
                "es_from_2x2_prop()",
                "OR+RR+NNT", "D+G+R+Z",
                "Non-adjusted")


  dat[dat$hierarch_name == "rr_se",
      2:7] <- c("Risk ratio value + standard error",
                "Section 3. https://metaconvert.org/input.html",
                "es_from_rr_se()",
                "RR", "OR+NNT",
                "Non-adjusted")
  dat[dat$hierarch_name == "rr_ci",
      2:7] <- c("Risk ratio value + 95% confidence interval",
                "Section 3. https://metaconvert.org/input.html",
                "es_from_rr_ci()",
                "RR", "OR+NNT",
                "Non-adjusted")
  dat[dat$hierarch_name == "rr_pval",
      2:7] <- c("Risk ratio value + p-value",
                "Section 3. https://metaconvert.org/input.html",
                "es_from_rr_pval",
                "RR", "OR+NNT",
                "Non-adjusted")

  dat[dat$hierarch_name == "beta_std",
      2:7] <- c("Standardized regression estimate (with a dichotomous predictor)",
                "Section 13. https://metaconvert.org/input.html",
                "es_from_beta_std()",
                "D+G", "OR+R+Z",
                "Non-adjusted")
  dat[dat$hierarch_name == "beta_unstd",
      2:7] <- c("Unstandardized regression estimate (with a dichotomous predictor)",
                "Section 13. https://metaconvert.org/input.html",
                "es_from_beta_unstd()",
                "D+G", "OR+R+Z",
                "Non-adjusted")

  dat[dat$hierarch_name == "cases_time",
      2:7] <- c("Number of cases and time of disee-free observation time in two independent groups",
                "Section 5. https://metaconvert.org/input.html",
                "es_from_cases_time()",
                "IRR", "N/A",
                "Non-adjusted")

  dat[dat$hierarch_name == "variability_means_sd",
      2:7] <- c("Standard devations (and means) of two independent groups",
                "Section 6. https://metaconvert.org/input.html",
                "es_variab_from_means_sd()",
                "VR+CVR", "N/A",
                "Non-adjusted")
  dat[dat$hierarch_name == "variability_means_ci",
      2:7] <- c("95% CIs around means (and means) of two independent groups",
                "Section 6. https://metaconvert.org/input.html",
                "es_variab_from_means_ci()",
                "VR+CVR", "N/A",
                "Non-adjusted")
  dat[dat$hierarch_name == "variability_means_se",
      2:7] <- c("Standard errors (and means) of two independent groups",
                "Section 6. https://metaconvert.org/input.html",
                "es_variab_from_means_se()",
                "VR+CVR", "N/A",
                "Non-adjusted")

  if (type_of_measure == "natural+converted") {
    dat = switch(measure,
                 "d" = subset(dat, grepl("D+G", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "g" = subset(dat, grepl("G", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "md" = subset(dat, grepl("MD", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "or" = subset(dat, grepl("OR", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "rr" = subset(dat, grepl("RR", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "nnt" = subset(dat, grepl("NNT", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "r" = subset(dat, grepl("R+Z", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "z" = ,
                 "Z" = subset(dat, grepl("R+Z", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "irr" = subset(dat, grepl("IRR", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "logvr" = subset(dat, grepl("VR", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "logcvr" = subset(dat, grepl("CVR", paste0(dat[, 6], dat[, 5]), fixed = TRUE)),
                 "all" = dat
    )
  } else if (type_of_measure == "natural") {
    dat = switch(measure,
                 "d" = subset(dat, grepl("D+G", paste0(dat[, 5]), fixed = TRUE)),
                 "g" = subset(dat, grepl("G", paste0(dat[, 5]), fixed = TRUE)),
                 "md" = subset(dat, grepl("MD", paste0(dat[, 5]), fixed = TRUE)),
                 "or" = subset(dat, grepl("OR", paste0(dat[, 5]), fixed = TRUE)),
                 "rr" = subset(dat, grepl("RR", paste0(dat[, 5]), fixed = TRUE)),
                 "nnt" = subset(dat, grepl("NNT", paste0(dat[, 5]), fixed = TRUE)),
                 "r" = subset(dat, grepl("R+Z", paste0(dat[, 5]), fixed = TRUE)),
                 "z" = ,
                 "Z" = subset(dat, grepl("R+Z", paste0(dat[, 5]), fixed = TRUE)),
                 "irr" = subset(dat, grepl("IRR", paste0(dat[, 5]), fixed = TRUE)),
                 "logvr" = subset(dat, grepl("VR", paste0(dat[, 5]), fixed = TRUE)),
                 "logcvr" = subset(dat, grepl("CVR", paste0(dat[, 5]), fixed = TRUE)),
                 "all" = subset(dat, dat[, 5] != "N/A")
    )
  } else {
    stop(paste0("Incorrect type_of_measure: ", type_of_measure))
  }

  dat = dat[order(dat$list_input_data), ]
  dat = dat[order(dat$adjusted_input_data, decreasing=TRUE), ]

  if (missing(name)) name = "data_extract_metaConvert"
  if (extension == "data.frame") {
    dat
  } else if (extension != ".xlsx") {
    rio::export(dat, paste0(name, extension))
    if (verbose) {
      message(paste0("Sheet created at location:\n\n", getwd(), "/",
                         name, extension,
                         "\n\n Good luck with your research!"))
    }
  } else {
    rio::export(dat, paste0(name, extension), overwrite = TRUE)
    if (verbose) {
      message(paste0("Sheet created at location:\n\n", getwd(), "/",
                         name, extension,
                         "\n\n Good luck with your research!"))
    }
  }
}


