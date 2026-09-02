#' This function checks the input data to have the correct format and that it has no incongruencies
#'
#' @param x a data.frame passed to convert_df() to be checked for correct format and consistency
#'
#' @noRd
.check_data <- function(x, split_adjusted = TRUE, format = "wide", main_es = TRUE) {
  #### Check input type ----------------

  # we check if the input object is a dataframe. if it is not, we display an error message
  if (!("data.frame" %in% class(x))) {
    stop("The object passed to the convert_df() function is not a data.frame object.")
  } else if (length(class(x)) > 1 & "data.frame" %in% class(x)) { # if it is not only a dataframe, we convert it
    x <- data.frame(x)
  }

  if (length(colnames(x)) == 0) {
    stop("Dataframe passed to the convert_df() function has no column. Check format of the dataset.")
  } else if (nrow(x) == 0) {
    stop("Dataframe passed to the convert_df() function has no row. Check format of the dataset.")
  }


  #### Initialize some settings ------
  status <- "OK" # this is the general status of the whole analysis
  error_msgs <- c() # this will contain all errors and warnings for each study
  row_wrong <- c() # row of error
  col_wrong <- c() # col of error
  column_errors <- rep("", nrow(x))
  column_type_errors <- rep("", nrow(x))
  user_input_es <- user_input_se <-
    user_input_es_crd <- user_input_se_crd <-
    user_input_es_adj <- user_input_se_adj <- rep(NA_real_, nrow(x))
  x$row_id <- 1:nrow(x)

  # Backwards compatibility: map old user_es_measure column names to new names
  if ("user_es_measure_crude" %in% colnames(x) && !"user_es_original_measure_crude" %in% colnames(x)) {
    x$user_es_original_measure_crude <- x$user_es_measure_crude
  }
  if ("user_es_measure_adj" %in% colnames(x) && !"user_es_original_measure_adj" %in% colnames(x)) {
    x$user_es_original_measure_adj <- x$user_es_measure_adj
  }

  returned <- paste0("Checkings finished.") # object returned


  expected_cols <- c(
    "study_id", "author", "year", "predictor", "outcome", "info_expected",
    # sample size
    "n_sample", "n_exp", "n_nexp",
    # CRUDE input
    "user_es_original_measure_crude",
    "user_es_target_measure_crude",
    "user_es_crude",
    "user_se_crude",
    "user_ci_lo_crude",
    "user_ci_up_crude",
    # smd/g
    "reverse_d", "cohen_d", "reverse_g", "hedges_g",
    # regression
    "reverse_beta_std", "beta_std", "reverse_beta_unstd", "beta_unstd", "sd_dv",
    # regression t-statistic (partial correlation)
    "reverse_linreg_t", "linreg_t", "n_covariates",
    # regression coefficient (partial correlation)
    "reverse_linreg_b", "linreg_b", "linreg_b_se",
    "linreg_b_ci_lo", "linreg_b_ci_up",
    "reverse_linreg_b_pval", "linreg_b_pval",
    # md
    "reverse_md", "md", "md_sd", "md_se", "md_ci_lo", "md_ci_up",
    "md_pval",
    # quarts
    "reverse_med", "min_exp", "q1_exp", "med_exp", "q3_exp", "max_exp",
    "min_nexp", "q1_nexp", "med_nexp", "q3_nexp", "max_nexp",
    # means post
    "reverse_means", "mean_sd_pooled",
    "mean_exp", "mean_sd_exp", "mean_se_exp", "mean_ci_lo_exp", "mean_ci_up_exp",
    "mean_nexp", "mean_sd_nexp", "mean_se_nexp", "mean_ci_lo_nexp", "mean_ci_up_nexp",
    # t f eta
    "reverse_student_t", "student_t",
    "reverse_student_t_pval", "student_t_pval",
    "reverse_anova_f", "anova_f",
    "reverse_anova_f_pval", "anova_f_pval",
    "reverse_etasq", "etasq",
    "reverse_pt_bis_r", "pt_bis_r", "reverse_pt_bis_r_pval", "pt_bis_r_pval",

    # pre/post
    "reverse_means_pre_post", "mean_pre_exp", "mean_pre_nexp", "mean_pre_sd_exp", "mean_pre_sd_nexp",
    "mean_pre_se_exp", "mean_pre_se_nexp", "mean_pre_ci_lo_exp", "mean_pre_ci_up_exp", "mean_pre_ci_lo_nexp", "mean_pre_ci_up_nexp",
    "reverse_mean_change", "mean_change_exp", "mean_change_nexp",
    "mean_change_sd_exp", "mean_change_sd_nexp",
    "mean_change_se_exp", "mean_change_se_nexp",
    "mean_change_ci_lo_exp", "mean_change_ci_up_exp",
    "mean_change_ci_lo_nexp", "mean_change_ci_up_nexp",
    "mean_change_pval_exp", "mean_change_pval_exp",
    "mean_change_pval_nexp", "mean_change_pval_nexp",

    "r_pre_post_exp", "r_pre_post_nexp",
    "paired_t_exp", "paired_t_nexp", "reverse_paired_t",
    "paired_t_pval_exp", "paired_t_pval_nexp", "reverse_paired_t_pval",
    "reverse_paired_t_pval_exp", "reverse_paired_t_pval_nexp",
    "paired_f_exp", "paired_f_nexp", "reverse_paired_f",
    "reverse_paired_f_exp", "reverse_paired_f_nexp",
    "paired_f_pval_exp", "paired_f_pval_nexp", "reverse_paired_f_pval",
    "reverse_paired_f_pval_exp", "reverse_paired_f_pval_nexp",
    # plot
    "reverse_plot_means",
    "plot_mean_exp", "plot_mean_nexp",
    "plot_mean_sd_lo_exp", "plot_mean_sd_lo_nexp",
    "plot_mean_sd_up_exp", "plot_mean_sd_up_nexp",
    "plot_mean_se_lo_exp", "plot_mean_se_lo_nexp",
    "plot_mean_se_up_exp", "plot_mean_se_up_nexp",
    "plot_mean_ci_lo_exp", "plot_mean_ci_lo_nexp",
    "plot_mean_ci_up_exp", "plot_mean_ci_up_nexp",
    # ADJUSTED input
    "user_es_original_measure_adj",
    "user_es_target_measure_adj",
    "user_es_adj",
    "user_se_adj",
    "user_ci_lo_adj",
    "user_ci_up_adj",
    "cohen_d_adj",
    "etasq_adj", "reverse_etasq_adj",
    # ANCOVA
    "reverse_ancova_means", "ancova_mean_sd_pooled", "cov_outcome_r", "n_cov_ancova",
    "ancova_mean_exp", "ancova_mean_nexp",
    "ancova_mean_sd_exp", "ancova_mean_sd_nexp",
    "ancova_mean_se_exp", "ancova_mean_se_nexp",
    "ancova_mean_ci_lo_exp", "ancova_mean_ci_up_exp",
    "ancova_mean_ci_lo_nexp", "ancova_mean_ci_up_nexp",
    "reverse_ancova_t", "ancova_t",
    "reverse_ancova_f", "ancova_f",
    "reverse_ancova_t_pval", "ancova_t_pval",
    "reverse_ancova_f_pval", "ancova_f_pval",
    "reverse_ancova_md", "ancova_md", "ancova_md_sd", "ancova_md_se",
    "ancova_md_ci_lo", "ancova_md_ci_up",
    "ancova_md_pval",
    # plot
    "reverse_plot_ancova_means",
    "plot_ancova_mean_exp", "plot_ancova_mean_nexp",
    "plot_ancova_mean_sd_lo_exp", "plot_ancova_mean_sd_lo_nexp",
    "plot_ancova_mean_sd_up_exp", "plot_ancova_mean_sd_up_nexp",
    "plot_ancova_mean_se_lo_exp", "plot_ancova_mean_se_lo_nexp",
    "plot_ancova_mean_se_up_exp", "plot_ancova_mean_se_up_nexp",
    "plot_ancova_mean_ci_lo_exp", "plot_ancova_mean_ci_lo_nexp",
    "plot_ancova_mean_ci_up_exp", "plot_ancova_mean_ci_up_nexp",

    # variability
    "reverse_means_variability", "reverse_means_change_variability",

    # 2x2 table
    "reverse_2x2", "baseline_risk", "baseline_rate", "small_margin_prop",
    "n_cases_exp", "n_cases_nexp", "n_controls_exp", "n_controls_nexp",
    "reverse_prop", "prop_cases_exp", "prop_cases_nexp",
    "n_cases", "n_controls",
    # single-group proportion
    "prop",
    # or
    "reverse_or", "or", "logor", "logor_se", "or_ci_lo", "or_ci_up", "logor_ci_lo", "logor_ci_up",
    "reverse_or_pval", "or_pval",
    "reverse_logreg_t", "logreg_t",
    # rr
    "reverse_rr", "rr", "logrr", "logrr_se",
    "rr_ci_lo", "rr_ci_up", "logrr_ci_lo", "logrr_ci_up",
    "reverse_rr_pval", "rr_pval",
    # rd (risk difference)
    "reverse_rd", "rd", "rd_se",
    "rd_ci_lo", "rd_ci_up",
    "reverse_rd_pval", "rd_pval",
    # X2, PHI, COR PB
    "reverse_chisq", "chisq", "reverse_chisq_pval", "chisq_pval",
    "yates_chisq",
    "reverse_phi", "phi", #"reverse_phi_pval", "phi_pval",
    # r
    "reverse_pearson_r", "pearson_r", "reverse_fisher_z", "fisher_z",
    "unit_increase_iv",
    "unit_type",
    "sd_iv",
    # survival
    "time_exp", "time_nexp", "reverse_irr",
    # correlation (Spearman)
    "reverse_spearman_r", "spearman_r",
    # psychometrics
    "cronbach_alpha", "cronbach_alpha_se", "cronbach_alpha_ci_lo",
    "cronbach_alpha_ci_up", "n_items",
    "n_response_categories", "scale_mean", "scale_min",
    "omega", "omega_se", "omega_ci_lo", "omega_ci_up", "omega_type", "omega_estimator",
    "icc", "icc_se", "icc_ci_lo", "icc_ci_up", "n_measurements", "icc_type",
    "discard",
    # multi-arm trial pooling
    "pool_side"
  )
  # res = data.frame(Comment = rep("N/A", 5))
  # res[, c(expected_cols)] <- NA
  # res$study_id = paste("study ", 1:5)
  # res$author = c("Bern", "Smith", "Jonh", "Doe", "Bakari")
  # res$year = c(2018, 2011, 2023, 2001, 2015)
  # res$factor = "Early Intensive Behavioral Intervention"
  # res$outcome = "Vineland-III"
  # rio::export(res, "C:/Users/Corentin Gosling/drive_gmail/Recherche/es.utils/web/background/extraction_sheet.xlsx", overwrite=TRUE)

  expected_cols_type <- ifelse(grepl("study_id", expected_cols, fixed = TRUE) |
    grepl("info_expected", expected_cols, fixed = TRUE) |
    grepl("author", expected_cols, fixed = TRUE) |
    grepl("year", expected_cols, fixed = TRUE) |
    grepl("predictor", expected_cols, fixed = TRUE) |
    grepl("outcome", expected_cols, fixed = TRUE) |
    grepl("measure", expected_cols, fixed = TRUE) |
    grepl("unit_type", expected_cols, fixed = TRUE) |
    grepl("icc_type", expected_cols, fixed = TRUE) |
    grepl("discard", expected_cols, fixed = TRUE),
  "char", "numeric"
  )

  # Set reverse_* and discard columns to logical type
  expected_cols_type[grepl("reverse", expected_cols, fixed = TRUE)] <- "logical"
  expected_cols_type[grepl("discard", expected_cols, fixed = TRUE)] <- "logical"

  expected_cols_type[which(expected_cols == "cov_outcome_r")] <- "numeric"
  # n_measurements contains "measure", so the grepl above would otherwise type it
  # "char"; force it numeric. n_items matches no char-triggering grepl, so its override
  # is a harmless guard kept for symmetry.
  expected_cols_type[which(expected_cols == "n_measurements")] <- "numeric"
  expected_cols_type[which(expected_cols == "n_items")] <- "numeric"
  # V44's two inputs. Neither matches a char-triggering grepl, so these are guards like
  # n_items, but they are the columns the floor-effect flag would go dark on, without
  # saying so, if a future pattern ever started matching them.
  expected_cols_type[which(expected_cols == "n_response_categories")] <- "numeric"
  expected_cols_type[which(expected_cols == "scale_mean")] <- "numeric"
  expected_cols_type[which(expected_cols == "scale_min")] <- "numeric"
  # omega_type is a character estimand label; without this override the automatic
  # grepl() detection types it numeric and every value becomes NA (same reason
  # n_items / n_measurements need explicit overrides).
  expected_cols_type[which(expected_cols == "omega_type")] <- "char"
  expected_cols_type[which(expected_cols == "omega_estimator")] <- "char"
  expected_cols_type[which(expected_cols == "pool_side")] <- "char"
  # per-row Yates flag for chi-square back-derivation
  expected_cols_type[which(expected_cols == "yates_chisq")] <- "logical"
  # View(cbind(expected_cols, expected_cols_type))
  # colnames(x) <- tolower(colnames(x))

  # remove all rows that should be discarded from analyses
  if (any(x$discard %in% c("yes", "Yes", "remove", "removed", "TRUE", TRUE))) {
    removed_rows <- which(x$discard %in% c("yes", "Yes", "remove", "removed", "TRUE"))
    status <- ifelse(status == "ERROR", "ERROR", "WARNING")
    error_msgs <- append(error_msgs, paste0("Some rows of the original dataset have been removed based on the 'discard' column inputs: rows = ", paste(removed_rows, collapse = ", "), " (only a warning, not an error)."))
    column_errors <- column_errors[-removed_rows]
    column_type_errors <- column_type_errors[-removed_rows]
    x <- subset(x, !discard %in% c("yes", "Yes", "remove", "removed", "TRUE", TRUE))
    # NB: a stray `situation <- situation[-removed_rows]` was removed here. `situation`
    # was never defined anywhere in this function, so any truthy value in the documented
    # `discard` column aborted the whole convert_df() run with "object 'situation' not
    # found". The actual row removal is performed by the subset() call above.
    # row_id was assigned as 1:nrow() before the discard drop, so it is now
    # non-contiguous (e.g. c(1, 3)). Downstream code (the dispersion merge in
    # .generate_df) requires row_id == seq_len(nrow(x)); re-sequence it here so the
    # retained rows are renumbered and the invariant holds.
    x$row_id <- seq_len(nrow(x))
  }

  # set as NA columns not included in the dataset
  for (col in expected_cols) {
    if (!(col %in% colnames(x))) {
      x[, col] <- NA
    }
  }

  #### set "na", "inf" or blank as NA
  x[x == "" | x == " " | x == "na" | x == "NA" | x == "n/a" | x == "N/A" |
      x == "inf" | x == "infinity" | x == "Infinity" | x == "INFINITY" | x == "INF" |
      x == "Inf"] <- NA

  #### Convert numerical and logical columns to correct type
  for (j in 1:length(expected_cols)) {
    idx <- which(colnames(x) == expected_cols[j])
    x[, idx] <- as.character(x[, idx])
    if (expected_cols_type[j] == "numeric") {
      # check the presence of . and , in numeric columns
      if (any(grepl("\\.", x[, idx]) & grepl(",", x[, idx]))) {
        commas_and_points <- which(grepl("\\.", x[, idx]) & grepl(",", x[, idx]))
        status <- "ERROR"
        # column_errors[commas_and_points] = paste(column_errors[commas_and_points], "Column '", paste(colnames(x)[idx]), "' contains both '.' and ','. //")
        stop(paste("Column '", paste(colnames(x)[idx]), "' contains both '.' and ','. //"))
        # column_type_errors[commas_and_points] = "ERROR"
      } else if (any(grepl(",", x[, idx]))) {
        commas <- which(grepl(",", x[, idx]))
        status <- ifelse(status == "ERROR", "ERROR", "WARNING")
        column_errors[commas] <- paste(column_errors[commas], "A ',' has been converted to a '.' in column '", paste(colnames(x)[idx]), "' (this is only a warning, not an error). //")
        warning("A ',' has been converted to a '.' in column '", paste(colnames(x)[idx]), "' (this is only a warning, not an error).")
        x[, idx] <- gsub(",", ".", x[, idx])
      }
      # check the presence of non-numeric characters in numeric columns
      if (any(suppressWarnings(is.na(as.numeric(as.character(na.omit(x[, idx]))))))) {
        not_num <- x[which(!is.na(x[, idx])), ][suppressWarnings(is.na(as.numeric(as.character(na.omit(x[, idx]))))), ]$row_id
        stop("Non-numeric characters ('", paste(unique(x[not_num, idx]), collapse = " ' / ' "), "') in column '", paste(colnames(x)[idx], "'. //"))
      }
      # convert numeric columns to numeric format
      x[, idx] <- as.numeric(as.character(x[, idx]))
    } else if (expected_cols_type[j] == "logical") {
      # convert logical columns to logical format
      x[, idx] <- as.logical(x[, idx])
      # Set NA values to FALSE for logical columns - except yates_chisq,
      # where NA means "use the scalar yates_chisq default supplied to convert_df()"
      # (handled downstream in main_convert_df.R).
      if (expected_cols[j] != "yates_chisq") {
        x[is.na(x[, idx]), idx] <- FALSE
      }
    }
  }



  x$n_cases <- with(x, ifelse(is.na(n_cases) & !is.na(n_cases_exp) & !is.na(n_cases_nexp),
    n_cases_exp + n_cases_nexp,
    n_cases
  ))
  x$n_controls <- with(x, ifelse(is.na(n_controls) & !is.na(n_controls_exp) & !is.na(n_controls_nexp),
    n_controls_exp + n_controls_nexp,
    n_controls
  ))
  x$n_exp <- with(x, ifelse(is.na(n_exp) & !is.na(n_cases_exp) & !is.na(n_controls_exp),
    n_cases_exp + n_controls_exp,
    n_exp
  ))
  x$n_nexp <- with(x, ifelse(is.na(n_nexp) & !is.na(n_cases_nexp) & !is.na(n_controls_nexp),
    n_cases_nexp + n_controls_nexp,
    n_nexp
  ))
  x$n_sample <- with(x, ifelse(is.na(n_sample) & !is.na(n_exp) & !is.na(n_nexp),
    n_exp + n_nexp,
    ifelse(is.na(n_sample) & !is.na(n_cases) & !is.na(n_controls),
      n_cases + n_controls, n_sample
    )
  ))

  col_names <- c(
    "all_info", "measure", "info_measure", "n_estimations", "es_selected", "info_used",
    "es", "se", "es_ci_lo", "es_ci_up",
    "min_info", "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
    "max_info", "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up",
    "diff_min_max", "overlap_min_max", "dispersion_es"
  )

  if (split_adjusted == TRUE & format == "wide" & main_es == TRUE) {
    res_crude <- data.frame(matrix(NA, nrow = nrow(x), ncol = length(col_names)))
    colnames(res_crude) <- paste0(col_names, "_crude")

    res_adjusted <- data.frame(matrix(NA, nrow = nrow(x), ncol = length(col_names)))
    colnames(res_adjusted) <- paste0(col_names, "_adjusted")

    res <- cbind(
      x[, c("row_id", "study_id", "author", "year", "predictor", "outcome", "info_expected")],
      res_crude, res_adjusted,
      subset(x, select = -c(row_id, study_id, author, year, predictor, outcome, info_expected))
    )
  } else {#if ((split_adjusted == TRUE & format == "long") | main_es == FALSE) {
    res_crude <- data.frame(matrix(NA, nrow = nrow(x), ncol = length(col_names)))
    colnames(res_crude) <- col_names
    res <- cbind(
      x[, c("row_id", "study_id", "author", "year", "predictor", "outcome", "info_expected")],
      adjusted_input = NA,
      res_crude,
      subset(x, select = -c(row_id, study_id, author, year, predictor, outcome, info_expected))
    )
  }

  return(res)
}
# rio::export(data.frame(do.call(rbind, lapply(object, function(x) x[1, "info_used"]))),
#             "possible_inputs.xlsx", overwrite = TRUE)

# if (es_adj == "separate" & !is.na(es_adj)) {
#
#   # we extract information on ES already given by the users
#   if ("es_crude" %in% colnames(x) & "se_crude" %in% colnames(x)) {
#     user_input_es_crd = x$es_crude
#     user_input_se_crd = x$se_crude
#   } else {
#     es_crude = se_crude = rep(NA_real_, nrow(x))
#   }
#
#   if ("es_adjusted" %in% colnames(x) & "se_adjusted" %in% colnames(x)) {
#     user_input_es_adj = x$es_adjusted
#     user_input_se_adj = x$se_adjusted
#   } else {
#     es_adjusted = se_adjusted = rep(NA_real_, nrow(x))
#   }
#
#   res_return_crude <- data.frame(matrix(NA, nrow = nrow(x), ncol = length(col_names)))
#   colnames(res_return_crude) <- paste0(col_names, "_crude")
#
#   res_return_adjusted <- data.frame(matrix(NA, nrow = nrow(x), ncol = length(col_names)))
#   colnames(res_return_adjusted) <- paste0(col_names, "_adjusted")
#
#   res_return = cbind(res_return_crude, res_return_adjusted)
#
# } else {
#
#   # we extract information on ES already given by the users
#   if ("es" %in% colnames(x) & "se" %in% colnames(x)) {
#     user_input_es = x$es
#     user_input_se = x$se
#   } else {
#     es_crude = se_crude = rep(NA_real_, nrow(x))
#   }
#   res_return <- data.frame(matrix(NA, nrow = nrow(x), ncol = length(col_names)))
#   colnames(res_return) <- col_names
# }
#
# res = cbind(res_return, subset(x, select = -c(row_id)))
# attr(res, "user_input_es") <- user_input_es
# attr(res, "user_input_se") <- user_input_se
# attr(res, "user_input_es_crd") <- user_input_es_crd
# attr(res, "user_input_se_crd") <- user_input_se_crd
# attr(res, "user_input_es_adj") <- user_input_es_adj
# attr(res, "user_input_se_adj") <- user_input_se_adj
# attr(res, "user_input_se_adj") <- user_input_se_adj
