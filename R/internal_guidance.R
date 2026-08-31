################### es guidance ###################

#' mapping of methods to their required input columns
#' @param omega_se_source the active `omega_se_source`. Only the omega entry depends on
#'   it: the closed form needs `n_sample` and `n_items`, the reported route cannot use them.
#' @return named list of lists with $specific and $shared
#' @noRd
.method_required_columns <- function(omega_se_source = "reported") {
  list(
    # smd
    cohen_d = list(
      specific = c("cohen_d"),
      shared   = c("n_exp", "n_nexp")
    ),
    hedges_g = list(
      specific = c("hedges_g"),
      shared   = c("n_exp", "n_nexp")
    ),
    cohen_d_adj = list(
      specific = c("cohen_d_adj"),
      shared   = c("n_exp", "n_nexp", "n_cov_ancova", "cov_outcome_r")
    ),

    # or
    or = list(
      specific = list(c("or", "logor")),
      shared   = c("n_exp", "n_nexp")
    ),
    or_se = list(
      specific = list(c("or", "logor"), "logor_se"),
      shared   = character(0)
    ),
    or_ci = list(
      specific = list(c("or", "logor"), c("or_ci_lo", "logor_ci_lo"), c("or_ci_up", "logor_ci_up")),
      shared   = character(0)
    ),
    or_pval = list(
      specific = list(c("or", "logor"), "or_pval"),
      shared   = c("n_exp", "n_nexp")
    ),
    logreg_t = list(
      specific = list(c("or", "logor", "rr", "logrr"), "logreg_t"),
      shared   = c("n_exp", "n_nexp")
    ),

    # correlation
    pearson_r = list(
      specific = c("pearson_r"),
      shared   = c("n_sample")
    ),
    fisher_z = list(
      specific = c("fisher_z"),
      shared   = c("n_sample")
    ),

    # means
    means_sd = list(
      specific = c("mean_exp", "mean_sd_exp", "mean_nexp", "mean_sd_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    means_se = list(
      specific = c("mean_exp", "mean_se_exp", "mean_nexp", "mean_se_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    means_ci = list(
      specific = c("mean_exp", "mean_ci_lo_exp", "mean_ci_up_exp",
                    "mean_nexp", "mean_ci_lo_nexp", "mean_ci_up_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    means_sd_pooled = list(
      specific = c("mean_exp", "mean_nexp", "mean_sd_pooled"),
      shared   = c("n_exp", "n_nexp")
    ),

    # means from plot
    means_plot = list(
      specific = c("plot_mean_exp", "plot_mean_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),

    # pre-post means
    means_sd_pre_post = list(
      specific = c("mean_pre_exp", "mean_exp", "mean_pre_sd_exp", "mean_sd_exp",
                    "mean_pre_nexp", "mean_nexp", "mean_pre_sd_nexp", "mean_sd_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    means_se_pre_post = list(
      specific = c("mean_pre_exp", "mean_exp", "mean_pre_se_exp", "mean_se_exp",
                    "mean_pre_nexp", "mean_nexp", "mean_pre_se_nexp", "mean_se_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    means_ci_pre_post = list(
      specific = c("mean_pre_exp", "mean_exp",
                    "mean_pre_ci_lo_exp", "mean_pre_ci_up_exp",
                    "mean_ci_lo_exp", "mean_ci_up_exp",
                    "mean_pre_nexp", "mean_nexp",
                    "mean_pre_ci_lo_nexp", "mean_pre_ci_up_nexp",
                    "mean_ci_lo_nexp", "mean_ci_up_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),

    # mean change
    mean_change_sd = list(
      specific = c("mean_change_exp", "mean_change_sd_exp",
                    "mean_change_nexp", "mean_change_sd_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    mean_change_se = list(
      specific = c("mean_change_exp", "mean_change_se_exp",
                    "mean_change_nexp", "mean_change_se_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    mean_change_ci = list(
      specific = c("mean_change_exp",
                    "mean_change_ci_lo_exp", "mean_change_ci_up_exp",
                    "mean_change_nexp",
                    "mean_change_ci_lo_nexp", "mean_change_ci_up_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    mean_change_pval = list(
      specific = c("mean_change_exp", "mean_change_pval_exp",
                    "mean_change_nexp", "mean_change_pval_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),

    # paired statistics
    paired_t = list(
      specific = c("paired_t_exp", "paired_t_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    paired_t_pval = list(
      specific = c("paired_t_pval_exp", "paired_t_pval_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    paired_f = list(
      specific = c("paired_f_exp", "paired_f_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    paired_f_pval = list(
      specific = c("paired_f_pval_exp", "paired_f_pval_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),

    # anova / student t
    student_t = list(
      specific = c("student_t"),
      shared   = c("n_exp", "n_nexp")
    ),
    student_t_pval = list(
      specific = c("student_t_pval"),
      shared   = c("n_exp", "n_nexp")
    ),
    anova_f = list(
      specific = c("anova_f"),
      shared   = c("n_exp", "n_nexp")
    ),
    anova_f_pval = list(
      specific = c("anova_f_pval"),
      shared   = c("n_exp", "n_nexp")
    ),
    etasq = list(
      specific = c("etasq"),
      shared   = c("n_exp", "n_nexp")
    ),
    etasq_adj = list(
      specific = c("etasq_adj"),
      shared   = c("n_exp", "n_nexp", "n_cov_ancova", "cov_outcome_r")
    ),
    pt_bis_r = list(
      specific = c("pt_bis_r"),
      shared   = c("n_exp", "n_nexp")
    ),
    pt_bis_r_pval = list(
      specific = c("pt_bis_r_pval"),
      shared   = c("n_exp", "n_nexp")
    ),

    # mean difference
    md_sd = list(
      specific = c("md", "md_sd"),
      shared   = c("n_exp", "n_nexp")
    ),
    md_se = list(
      specific = c("md", "md_se"),
      shared   = c("n_exp", "n_nexp")
    ),
    md_ci = list(
      specific = c("md", "md_ci_lo", "md_ci_up"),
      shared   = c("n_exp", "n_nexp")
    ),
    md_pval = list(
      specific = c("md", "md_pval"),
      shared   = c("n_exp", "n_nexp")
    ),

    # medians
    med_quarts = list(
      specific = c("q1_exp", "med_exp", "q3_exp",
                    "q1_nexp", "med_nexp", "q3_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    med_min_max = list(
      specific = c("min_exp", "med_exp", "max_exp",
                    "min_nexp", "med_nexp", "max_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    med_min_max_quarts = list(
      specific = c("min_exp", "q1_exp", "med_exp", "q3_exp", "max_exp",
                    "min_nexp", "q1_nexp", "med_nexp", "q3_nexp", "max_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),

    # regression
    beta_std = list(
      specific = c("beta_std", "sd_dv"),
      shared   = c("n_exp", "n_nexp")
    ),
    beta_unstd = list(
      specific = c("beta_unstd", "sd_dv"),
      shared   = c("n_exp", "n_nexp")
    ),
    linreg_t = list(
      specific = c("linreg_t"),
      shared   = c("n_sample", "n_covariates")
    ),
    linreg_b_se = list(
      specific = c("linreg_b", "linreg_b_se"),
      shared   = c("n_sample", "n_covariates")
    ),
    linreg_b_ci = list(
      specific = c("linreg_b", "linreg_b_ci_lo", "linreg_b_ci_up"),
      shared   = c("n_sample", "n_covariates")
    ),
    linreg_b_pval = list(
      specific = c("linreg_b", "linreg_b_pval"),
      shared   = c("n_sample", "n_covariates")
    ),

    # 2x2 table
    "2x2" = list(
      specific = c("n_cases_exp", "n_cases_nexp",
                    "n_controls_exp", "n_controls_nexp"),
      shared   = character(0)
    ),
    "2x2_sum" = list(
      specific = c("n_cases_exp", "n_cases_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    "2x2_prop" = list(
      specific = c("prop_cases_exp", "prop_cases_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),

    # chi-square / phi
    chisq = list(
      specific = c("chisq"),
      shared   = c("n_sample", "n_cases", "n_exp")
    ),
    chisq_pval = list(
      specific = c("chisq_pval"),
      shared   = c("n_sample", "n_cases", "n_exp")
    ),
    phi = list(
      specific = c("phi"),
      shared   = c("n_sample", "n_cases", "n_exp")
    ),

    # risk ratio
    rr_se = list(
      specific = list(c("rr", "logrr"), "logrr_se"),
      shared   = character(0)
    ),
    rr_ci = list(
      specific = list(c("rr", "logrr"), c("rr_ci_lo", "logrr_ci_lo"), c("rr_ci_up", "logrr_ci_up")),
      shared   = character(0)
    ),
    rr_pval = list(
      specific = list(c("rr", "logrr"), "rr_pval"),
      shared   = c("n_exp", "n_nexp")
    ),

    # risk difference
    rd_se = list(
      specific = c("rd", "rd_se"),
      shared   = character(0)
    ),
    rd_ci = list(
      specific = c("rd", "rd_ci_lo", "rd_ci_up"),
      shared   = character(0)
    ),
    rd_pval = list(
      specific = c("rd", "rd_pval"),
      shared   = character(0)
    ),

    # ancova means
    ancova_means_sd = list(
      specific = c("ancova_mean_exp", "ancova_mean_nexp",
                    "ancova_mean_sd_exp", "ancova_mean_sd_nexp"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_means_se = list(
      specific = c("ancova_mean_exp", "ancova_mean_nexp",
                    "ancova_mean_se_exp", "ancova_mean_se_nexp"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_means_ci = list(
      specific = c("ancova_mean_exp", "ancova_mean_nexp",
                    "ancova_mean_ci_lo_exp", "ancova_mean_ci_up_exp",
                    "ancova_mean_ci_lo_nexp", "ancova_mean_ci_up_nexp"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_means_sd_pooled = list(
      specific = c("ancova_mean_exp", "ancova_mean_nexp", "mean_sd_pooled"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_means_sd_pooled_adj = list(
      specific = c("ancova_mean_exp", "ancova_mean_nexp", "ancova_mean_sd_pooled"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_means_plot = list(
      specific = c("plot_ancova_mean_exp", "plot_ancova_mean_nexp"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),

    # ancova statistics
    ancova_t = list(
      specific = c("ancova_t"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_f = list(
      specific = c("ancova_f"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_t_pval = list(
      specific = c("ancova_t_pval"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_f_pval = list(
      specific = c("ancova_f_pval"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),

    # ancova mean difference
    ancova_md_sd = list(
      specific = c("ancova_md", "ancova_md_sd"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_md_se = list(
      specific = c("ancova_md", "ancova_md_se"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_md_ci = list(
      specific = c("ancova_md", "ancova_md_ci_lo", "ancova_md_ci_up"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),
    ancova_md_pval = list(
      specific = c("ancova_md", "ancova_md_pval"),
      shared   = c("n_exp", "n_nexp", "cov_outcome_r", "n_cov_ancova")
    ),

    # user input
    user_input_crude = list(
      specific = c("user_es_crude", "user_se_crude"),
      shared   = c("user_es_original_measure_crude")
    ),
    user_input_adj = list(
      specific = c("user_es_adj", "user_se_adj"),
      shared   = c("user_es_original_measure_adj")
    ),

    # irr
    cases_time = list(
      specific = c("n_cases_exp", "n_cases_nexp", "time_exp", "time_nexp"),
      shared   = character(0)
    ),

    # variability ratios
    variability_means_sd = list(
      specific = c("mean_exp", "mean_nexp", "mean_sd_exp", "mean_sd_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    variability_means_se = list(
      specific = c("mean_exp", "mean_nexp", "mean_se_exp", "mean_se_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),
    variability_means_ci = list(
      specific = c("mean_exp", "mean_nexp",
                    "mean_ci_lo_exp", "mean_ci_up_exp",
                    "mean_ci_lo_nexp", "mean_ci_up_nexp"),
      shared   = c("n_exp", "n_nexp")
    ),

    # single group pre-post means
    means_sd_pre_post_single_group = list(
      specific = c("mean_pre_exp", "mean_exp",
                    "mean_pre_sd_exp", "mean_sd_exp"),
      shared   = c("n_exp")
    ),
    means_se_pre_post_single_group = list(
      specific = c("mean_pre_exp", "mean_exp",
                    "mean_pre_se_exp", "mean_se_exp"),
      shared   = c("n_exp")
    ),
    means_ci_pre_post_single_group = list(
      specific = c("mean_pre_exp", "mean_exp",
                    "mean_pre_ci_lo_exp", "mean_pre_ci_up_exp",
                    "mean_ci_lo_exp", "mean_ci_up_exp"),
      shared   = c("n_exp")
    ),

    # single group mean change
    mean_change_sd_single_group = list(
      specific = c("mean_change_exp", "mean_change_sd_exp"),
      shared   = c("n_exp")
    ),
    mean_change_se_single_group = list(
      specific = c("mean_change_exp", "mean_change_se_exp"),
      shared   = c("n_exp")
    ),
    mean_change_ci_single_group = list(
      specific = c("mean_change_exp",
                    "mean_change_ci_lo_exp", "mean_change_ci_up_exp"),
      shared   = c("n_exp")
    ),
    mean_change_pval_single_group = list(
      specific = c("mean_change_exp", "mean_change_pval_exp"),
      shared   = c("n_exp")
    ),

    # single group paired t
    paired_t_single_group = list(
      specific = c("paired_t_exp"),
      shared   = c("n_exp")
    ),

    # single group proportions
    prop_single_group = list(
      specific = c("prop"),
      shared   = c("n_sample")
    ),
    prop_single_group_counts = list(
      specific = c("n_cases"),
      shared   = c("n_sample")
    ),

    # psychometric
    # omega alone is not enough: without a reported SE (or CI) there is no
    # aggregate-data sampling variance for omega, so the row cannot be pooled.
    # Listing omega_se as "shared" makes a bare-omega row a near miss, so the
    # user is told what to add instead of getting the generic
    # "No partial input data found" fallback.
    # Under omega_se_source = "closed_form" the SE is computed from n_sample and
    # n_items instead, so a row lacking them is dropped for a reason the "reported"
    # wording cannot state. Naming them is conditional rather than unconditional
    # because on the default "reported" route they cannot help, and advising a user
    # to hunt for a sample size that the active setting will ignore is worse than
    # saying nothing.
    omega = list(
      specific = c("omega"),
      shared = if (identical(omega_se_source, "closed_form")) {
        c("omega_se", "n_sample", "n_items")
      } else {
        c("omega_se")
      }
    ),
    # An alpha row needs a variance from somewhere: the study's own SE or CI, or
    # n_sample + n_items for the closed-form (n, k) SE. Listing all four as "shared"
    # makes a bare-alpha row a near miss that names every route out of it.
    cronbach_alpha = list(
      specific = c("cronbach_alpha"),
      shared   = c("n_sample", "n_items", "cronbach_alpha_se",
                   "cronbach_alpha_ci_lo", "cronbach_alpha_ci_up")
    ),
    # An ICC row needs a variance from somewhere: the study's own icc_se or CI,
    # or n_sample + n_measurements for the closed-form (n, k) SE. Listing all four
    # as "shared" makes a bare-ICC row a near miss naming every route out of it.
    icc = list(
      specific = c("icc"),
      shared   = c("n_sample", "n_measurements", "icc_se", "icc_ci_lo", "icc_ci_up")
    ),
    spearman_r = list(
      specific = c("spearman_r"),
      shared   = c("n_sample")
    )
  )
}


#' human-readable method descriptions
#' @return named character vector
#' @noRd
.method_descriptions <- function() {
  c(
    cohen_d                = "Cohen's d (direct input)",
    hedges_g               = "Hedges' g (direct input)",
    cohen_d_adj            = "Adjusted Cohen's d (direct input)",

    or                     = "Odds ratio (direct input)",
    or_se                  = "Odds ratio + standard error",
    or_ci                  = "Odds ratio + confidence interval",
    or_pval                = "Odds ratio + p-value",
    logreg_t               = "Logistic regression t-value",

    pearson_r              = "Pearson correlation (direct input)",
    fisher_z               = "Fisher's z (direct input)",

    means_sd               = "Means + SDs of two groups",
    means_se               = "Means + SEs of two groups",
    means_ci               = "Means + CIs of two groups",
    means_sd_pooled        = "Means + pooled SD",
    means_plot             = "Means from plot",

    means_sd_pre_post      = "Pre/post means + SDs of two groups",
    means_se_pre_post      = "Pre/post means + SEs of two groups",
    means_ci_pre_post      = "Pre/post means + CIs of two groups",

    mean_change_sd         = "Mean change + SDs of two groups",
    mean_change_se         = "Mean change + SEs of two groups",
    mean_change_ci         = "Mean change + CIs of two groups",
    mean_change_pval       = "Mean change + p-values of two groups",

    paired_t               = "Paired t-values of two groups",
    paired_t_pval          = "Paired t p-values of two groups",
    paired_f               = "Paired F-values of two groups",
    paired_f_pval          = "Paired F p-values of two groups",

    student_t              = "Student's t-value",
    student_t_pval         = "Student's t p-value",
    anova_f                = "ANOVA F-value",
    anova_f_pval           = "ANOVA F p-value",
    etasq                  = "Eta-squared",
    etasq_adj              = "Adjusted eta-squared",
    pt_bis_r               = "Point-biserial correlation",
    pt_bis_r_pval          = "Point-biserial correlation p-value",

    md_sd                  = "Mean difference + SD",
    md_se                  = "Mean difference + SE",
    md_ci                  = "Mean difference + CI",
    md_pval                = "Mean difference + p-value",

    med_quarts             = "Medians + quartiles",
    med_min_max            = "Medians + min/max",
    med_min_max_quarts     = "Medians + quartiles + min/max",

    beta_std               = "Standardized regression coefficient",
    beta_unstd             = "Unstandardized regression coefficient",
    linreg_t               = "Linear regression t-value",
    linreg_b_se            = "Regression coefficient + standard error",
    linreg_b_ci            = "Regression coefficient + confidence interval",
    linreg_b_pval          = "Regression coefficient + p-value",

    "2x2"                  = "2x2 contingency table (full)",
    "2x2_sum"              = "2x2 table (cases + group totals)",
    "2x2_prop"             = "2x2 table (proportions)",

    chisq                  = "Chi-squared value",
    chisq_pval             = "Chi-squared p-value",
    phi                    = "Phi coefficient",

    rr_se                  = "Risk ratio + standard error",
    rr_ci                  = "Risk ratio + confidence interval",
    rr_pval                = "Risk ratio + p-value",

    rd_se                  = "Risk difference + standard error",
    rd_ci                  = "Risk difference + confidence interval",
    rd_pval                = "Risk difference + p-value",

    ancova_means_sd        = "ANCOVA adjusted means + SDs",
    ancova_means_se        = "ANCOVA adjusted means + SEs",
    ancova_means_ci        = "ANCOVA adjusted means + CIs",
    ancova_means_sd_pooled = "ANCOVA adjusted means + pooled SD",
    ancova_means_sd_pooled_adj = "ANCOVA adjusted means + adjusted pooled SD",
    ancova_means_plot      = "ANCOVA adjusted means from plot",

    ancova_t               = "ANCOVA t-value",
    ancova_f               = "ANCOVA F-value",
    ancova_t_pval          = "ANCOVA t p-value",
    ancova_f_pval          = "ANCOVA F p-value",

    ancova_md_sd           = "ANCOVA mean difference + SD",
    ancova_md_se           = "ANCOVA mean difference + SE",
    ancova_md_ci           = "ANCOVA mean difference + CI",
    ancova_md_pval         = "ANCOVA mean difference + p-value",

    user_input_crude       = "User-provided ES (crude) - set user_es_original_measure for conversion",
    user_input_adj         = "User-provided ES (adjusted) - set user_es_original_measure for conversion",

    cases_time             = "Cases + person-time (IRR)",

    variability_means_sd   = "Variability ratio from means + SDs",
    variability_means_se   = "Variability ratio from means + SEs",
    variability_means_ci   = "Variability ratio from means + CIs",

    means_sd_pre_post_single_group  = "Pre/post means + SDs (single group)",
    means_se_pre_post_single_group  = "Pre/post means + SEs (single group)",
    means_ci_pre_post_single_group  = "Pre/post means + CIs (single group)",

    mean_change_sd_single_group     = "Mean change + SD (single group)",
    mean_change_se_single_group     = "Mean change + SE (single group)",
    mean_change_ci_single_group     = "Mean change + CI (single group)",
    mean_change_pval_single_group   = "Mean change + p-value (single group)",

    paired_t_single_group           = "Paired t-value (single group)",

    prop_single_group               = "Proportion (direct input)",
    prop_single_group_counts        = "Cases + total (counts)",

    cronbach_alpha                  = "Cronbach's alpha (needs n_sample + n_items, or a reported cronbach_alpha_se / cronbach_alpha_ci_lo + cronbach_alpha_ci_up)",
    omega                           = "McDonald's omega (needs omega_se on the natural scale, or omega_ci_lo + omega_ci_up)",
    icc                             = "ICC (needs n_sample + n_measurements, or a reported icc_se / icc_ci_lo + icc_ci_up)",
    spearman_r                      = "Spearman correlation (converted to Pearson)"
  )
}


#' methods applicable to a given measure (crude/adjusted)
#' @param measure the effect size measure
#' @param suffix "" or "_crude" or "_adjusted"
#' @return character vector
#' @noRd
.get_applicable_methods <- function(measure, suffix) {

  smd_post <- c(
    "cohen_d", "hedges_g",
    "means_sd", "means_se", "means_ci", "means_sd_pooled", "means_plot",
    "student_t", "student_t_pval", "anova_f", "anova_f_pval",
    "etasq", "pt_bis_r", "pt_bis_r_pval",
    "md_sd", "md_se", "md_ci", "md_pval",
    "med_quarts", "med_min_max", "med_min_max_quarts",
    "beta_std", "beta_unstd"
  )

  smd_paired <- c(
    "means_sd_pre_post", "means_se_pre_post", "means_ci_pre_post",
    "mean_change_sd", "mean_change_se", "mean_change_ci", "mean_change_pval",
    "paired_t", "paired_t_pval", "paired_f", "paired_f_pval"
  )

  smd_adj <- c(
    "cohen_d_adj", "etasq_adj",
    "ancova_means_sd", "ancova_means_se", "ancova_means_ci",
    "ancova_means_sd_pooled", "ancova_means_sd_pooled_adj", "ancova_means_plot",
    "ancova_t", "ancova_f", "ancova_t_pval", "ancova_f_pval",
    "ancova_md_sd", "ancova_md_se", "ancova_md_ci", "ancova_md_pval"
  )

  or_cat  <- c("or", "or_se", "or_ci", "or_pval", "logreg_t")
  cor_cat <- c("pearson_r", "fisher_z", "spearman_r")
  cont_cat <- c("2x2", "2x2_sum", "2x2_prop")
  phi_cat <- c("chisq", "chisq_pval", "phi")
  rr_cat  <- c("rr_se", "rr_ci", "rr_pval")
  rd_cat  <- c("rd_se", "rd_ci", "rd_pval")
  irr_cat <- c("cases_time")
  var_cat <- c("variability_means_sd", "variability_means_se", "variability_means_ci")

  within_group <- c(
    "means_sd_pre_post_single_group", "means_se_pre_post_single_group",
    "means_ci_pre_post_single_group",
    "mean_change_sd_single_group", "mean_change_se_single_group",
    "mean_change_ci_single_group", "mean_change_pval_single_group",
    "paired_t_single_group"
  )

  prop_sg     <- c("prop_single_group", "prop_single_group_counts")
  partial_cor <- c("linreg_t", "linreg_b_se", "linreg_b_ci", "linreg_b_pval")

  if (measure %in% c("dw", "gw", "mdw")) {
    return(c(within_group, "user_input_crude"))
  } else if (measure == "prop") {
    return(c(prop_sg, "user_input_crude"))
  } else if (measure == "alpha") {
    return(c("cronbach_alpha", "user_input_crude"))
  } else if (measure == "omega") {
    return(c("omega", "user_input_crude"))
  } else if (measure == "icc") {
    return(c("icc", "user_input_crude"))
  } else if (measure %in% c("hr", "loghr")) {
    # HR needs per-subject time-to-event data the wide extraction sheet cannot hold,
    # so convert_df() runs the user-input methods alone (hierarchy c(USER_crude,
    # USER_adjusted)). Without this case HR falls through to the switch default and
    # guidance offers about 55 means/SD/variability methods that can never yield a hazard
    # ratio. Respect the crude/adjusted split so es_guidance_crude/_adjusted are right.
    if (suffix == "_crude") return("user_input_crude")
    if (suffix == "_adjusted") return("user_input_adj")
    return(c("user_input_crude", "user_input_adj"))
  } else if (measure %in% c("rp", "zp")) {
    return(c(partial_cor, "user_input_crude"))
  }

  crude_cats <- switch(
    measure,
    "d" =, "g" =, "md" =
      c(smd_post, smd_paired, or_cat, cont_cat, cor_cat, phi_cat),
    "logor" =, "or" =
      c(or_cat, cont_cat, rr_cat, rd_cat, phi_cat, cor_cat, smd_post),
    "logrr" =, "rr" =
      c(rr_cat, cont_cat, or_cat, rd_cat, phi_cat),
    "logirr" =, "irr" =
      irr_cat,
    "nnt" =
      c(rd_cat, cont_cat, or_cat, rr_cat, irr_cat, phi_cat),
    "rd" =
      c(rd_cat, cont_cat, or_cat, rr_cat, irr_cat, phi_cat),
    "r" =, "z" =
      c(cor_cat, cont_cat, or_cat, phi_cat, smd_post, smd_paired),
    "logvr" =, "logcvr" =
      c(var_cat, smd_post, smd_paired),
    c(smd_post, smd_paired, or_cat, cont_cat, cor_cat,
      phi_cat, rr_cat, irr_cat, var_cat)
  )

  adj_cats <- switch(
    measure,
    "d" =, "g" =, "md" = smd_adj,
    # ancova methods not used for other measures
    character(0)
  )

  crude_all <- unique(c(crude_cats, "user_input_crude"))
  adj_all   <- unique(c(adj_cats, "user_input_adj"))

  if (suffix == "_crude") {
    return(crude_all)
  } else if (suffix == "_adjusted") {
    return(adj_all)
  } else {
    return(unique(c(crude_all, adj_all)))
  }
}


#' true if at least one column of the group is non-NA
#' @param spec_entry column name(s)
#' @param row_data one row of data
#' @return logical
#' @noRd
.has_any_value <- function(spec_entry, row_data) {
  for (col in spec_entry) {
    if (col %in% names(row_data) && !is.na(row_data[[col]])) return(TRUE)
  }
  return(FALSE)
}


#' near-miss methods for one row
#' @param row_data one row of raw_data
#' @param methods_to_check method names
#' @param method_cols output of .method_required_columns()
#' @param descriptions output of .method_descriptions()
#' @param max_suggestions max near-misses shown
#' @return character
#' @noRd
.diagnose_missing_data <- function(row_data, methods_to_check, method_cols,
                                    descriptions = .method_descriptions(),
                                    max_suggestions = 3) {
  near_misses <- list()

  for (method in methods_to_check) {
    if (!method %in% names(method_cols)) next

    mc <- method_cols[[method]]
    spec <- mc$specific
    shared <- mc$shared

    all_specific <- character(0)
    has_any_specific <- FALSE

    if (is.list(spec)) {
      for (entry in spec) {
        if (.has_any_value(entry, row_data)) {
          has_any_specific <- TRUE
        }
        all_specific <- c(all_specific, entry)
      }
    } else {
      all_specific <- spec
      for (col in spec) {
        if (col %in% names(row_data) && !is.na(row_data[[col]])) {
          has_any_specific <- TRUE
          break
        }
      }
    }

    if (!has_any_specific) next

    all_required <- unique(c(all_specific, shared))
    missing_cols <- character(0)

    if (is.list(spec)) {
      for (entry in spec) {
        if (length(entry) > 1) {
          # missing only if none of the alternatives present
          if (!.has_any_value(entry, row_data)) {
            missing_cols <- c(missing_cols, entry[1])
          }
        } else {
          if (!.has_any_value(entry, row_data)) {
            missing_cols <- c(missing_cols, entry)
          }
        }
      }
    } else {
      for (col in spec) {
        if (!(col %in% names(row_data)) || is.na(row_data[[col]])) {
          missing_cols <- c(missing_cols, col)
        }
      }
    }

    for (col in shared) {
      if (!(col %in% names(row_data)) || is.na(row_data[[col]])) {
        missing_cols <- c(missing_cols, col)
      }
    }

    missing_specific <- character(0)
    missing_shared   <- character(0)
    for (col in missing_cols) {
      if (col %in% shared) {
        missing_shared <- c(missing_shared, col)
      } else {
        missing_specific <- c(missing_specific, col)
      }
    }

    if (length(missing_cols) > 0) {
      near_misses[[length(near_misses) + 1]] <- list(
        method  = method,
        missing = missing_cols,
        missing_specific = missing_specific,
        missing_shared   = missing_shared,
        only_shared = length(missing_specific) == 0,
        n_missing = length(missing_cols)
      )
    }
  }

  if (length(near_misses) == 0) {
    return("No partial input data found. See data_extraction_sheet() for required columns.")
  }

  # methods missing only shared cols first
  only_shared <- vapply(near_misses, function(x) x$only_shared, logical(1))
  n_miss <- vapply(near_misses, function(x) x$n_missing, integer(1))
  near_misses <- near_misses[order(!only_shared, n_miss)]

  near_misses <- near_misses[seq_len(min(length(near_misses), max_suggestions))]

  msgs <- vapply(near_misses, function(nm) {
    label <- if (nm$method %in% names(descriptions)) {
      descriptions[[nm$method]]
    } else {
      nm$method
    }
    if (nm$only_shared) {
      paste0(label, " [complete input]: only add '",
             paste(nm$missing_shared, collapse = "' + '"), "'")
    } else {
      paste0(label, ": add '", paste(nm$missing, collapse = "' + '"), "'")
    }
  }, character(1))

  paste(msgs, collapse = "; ")
}


#' add the es_guidance column to the summary result
#' @param res the assembled summary result
#' @param suffix "" or "_crude" or "_adjusted"
#' @param raw_data original input data
#' @param measure the effect size measure
#' @param object the metaConvert object
#' @return data.frame
#' @noRd
.add_es_guidance <- function(res, suffix, raw_data, measure, object) {
  n <- nrow(res)
  guidance_col <- paste0("es_guidance", suffix)
  es_col <- paste0("es", suffix)

  method_cols  <- .method_required_columns(
    omega_se_source = if (is.null(attr(object, "omega_se_source"))) "reported"
                      else attr(object, "omega_se_source")
  )
  descriptions <- .method_descriptions()

  methods_to_check <- .get_applicable_methods(measure, suffix)

  # Do not advertise a route that the active analysis scale forbids. For
  # alpha/icc/prop the user-input route is refused whenever *_to_es is not the
  # identity (see .user_passthrough_blocked() in R/es_from_USER.R), so telling a
  # user to complete their user_es_* columns would send them to a dead end --
  # the row would come back NA with a warning either way.
  if (!is.null(object) && measure %in% c("alpha", "omega", "icc", "prop")) {
    blocked <- .user_passthrough_blocked(
      measure,
      alpha_to_es = if (is.null(attr(object, "alpha_to_es"))) "bonett" else attr(object, "alpha_to_es"),
      icc_to_es   = if (is.null(attr(object, "icc_to_es")))   "bonett" else attr(object, "icc_to_es"),
      prop_to_es  = if (is.null(attr(object, "prop_to_es")))  "raw"    else attr(object, "prop_to_es"))
    if (isTRUE(blocked))
      methods_to_check <- setdiff(methods_to_check, c("user_input_crude", "user_input_adj"))
  }

  res[[guidance_col]] <- rep("", n)

  # only rows with NA es
  for (i in seq_len(n)) {
    es_val <- suppressWarnings(as.numeric(as.character(res[[es_col]][i])))
    if (is.na(es_val)) {
      row_id <- res$row_id[i]
      raw_row_idx <- which(raw_data$row_id == row_id)

      if (length(raw_row_idx) == 1) {
        row_data <- as.list(raw_data[raw_row_idx, ])
        res[[guidance_col]][i] <- .diagnose_missing_data(
          row_data, methods_to_check, method_cols, descriptions
        )
      }
    }
  }

  return(res)
}
