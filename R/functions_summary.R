#' Synthesize information of an object of class \dQuote{metaConvert} into a dataframe
#'
#' @param object an object of class \dQuote{metaConvert}
#' @param digits an integer value specifying the number of decimal places for the rounding of numeric values. Default is 3.
#' @param ... other arguments that can be passed to the function
#'
#' @details
#' Summary method for objects of class \dQuote{metaConvert} produced by the \code{\link{convert_df}}
#' function. This function automatically:
#' 1. computes all effect sizes from all available input data
#' 1. selects, if requested, a  main effect size for each association/comparison using the information passed by
#' the user in the \code{es_selected} argument of the \code{convert_df} function
#' 1. identifies the smallest and largest effect size for each association/comparison
#' 1. estimates the absolute difference between the smallest and largest effect size for each
#' association/comparison
#' 1. estimates the percentage of overlap between the 95% confidence intervals of the smallest and
#' largest effect size for each association/comparison
#'
#' @return
#' This function returns a dataframe with many columns. We present below the information stored in each column of the returned dataframe
#'
#' **1. Raw user information.**
#' The first columns placed at the left of the returned dataset are simply information provided
#' by the users to facilitate the identification of each row.
#' If the following columns are missing in the original dataset, these columns will not appear in
#' the returned dataset.
#'
#' \tabular{ll}{
#'  \code{row_id} \tab Row number in the original dataset.\cr
#'  \tab \cr
#'  \code{study_id} \tab Identifier of the study.\cr
#'  \tab \cr
#'  \code{author} \tab Name of the author of the study.\cr
#'  \tab \cr
#'  \code{year} \tab Year of publication of the study.\cr
#'  \tab \cr
#'  \code{predictor} \tab Name of the predictor (intervention, risk factor, etc.).\cr
#'  \tab \cr
#'  \code{outcome} \tab Name of the outcome.\cr
#'  \tab \cr
#'  \code{info_expected} \tab Types of input data users expect to be used to estimate their effect size measure.\cr
#'  \tab \cr
#' }
#'
#' **2. Information on generated effect sizes.**
#' Then, the function returns information on calculations. For example, users can retrieve
#' the effect size measure estimated, the number and type(s) of input data allowing to estimate the
#' chosen effect size measure, and the method used to obtain a unique effect size if overlapping
#' input data were available.
#' These columns could have several suffix.
#' * If users requested to separate crude and adjusted estimates,
#' then the following columns will be presented with both a "_crude" suffix and a "_adjusted" suffix.
#' * If users did not request to separate the presentation of crude and adjusted estimates, the following columns
#' will have no suffix.
#'
#' For example, let's take column "all_info". It can be "all_info_crude" (all input data used to estimate any crude effect size),
#' "all_info_adjusted" (all input data leading to estimate any adjusted effect size),
#' or "all_info" (all input data leading to estimate any crude or adjusted effect sizes).
#'
#' To facilitate the presentation, we thus refer to these columns as \code{name_of_the_column*},
#' the \code{*} meaning that it could end by _crude, _adjusted or "".
#'
#' \tabular{ll}{
#'  \code{all_info*} \tab list of input data available in the dataset that was used to estimate any effect size measure.\cr
#'  \tab \cr
#'  \code{measure*} \tab effect size measure requested by the user.\cr
#'  \tab \cr
#'  \code{info_measure*} \tab input data available to estimate the requested effect size measure.\cr
#'  \tab \cr
#'  \code{n_estimations*} \tab number of input data available to estimate the requested effect size measure.\cr
#'  \tab \cr
#'  \code{es_selected*} \tab method chosen by users to estimate the main effect size when overlapping data are present.\cr
#'  \tab \cr
#'  \code{info_used*} \tab type of input data used to estimate the main effect size.\cr
#'  \tab \cr
#' }
#'
#' **3. Main effect size.**
#' The following columns contain the key information, namely, the main effect size + standard error + 95% CI.
#'
#' Again, the suffix of these columns can vary depending on the separation of effect sizes
#' estimated from crude and adjusted input data.
#'
#' \tabular{ll}{
#'  \code{es*} \tab main effect size value.\cr
#'  \tab \cr
#'  \code{se*} \tab standard error of the effect size.\cr
#'  \tab \cr
#'  \code{es_ci_lo*} \tab lower bound of the 95% CI around the effect size.\cr
#'  \tab \cr
#'  \code{es_ci_up*} \tab upper bound of the 95% CI around the effect size.\cr
#'  \tab \cr
#' }
#'
#' **4. Overlapping effect sizes**
#' These columns are useful ONLY if a given comparison (i.e., row) has multiple input data
#' enabling to compute the requested effect size measure.
#'
#' These columns identify the smallest/largest effect size per comparison,
#' and some indicators of consistency.
#'
#' Again, the suffix of these columns can vary depending on the separation of effect sizes
#' estimated from crude and adjusted input data.
#'
#' \tabular{ll}{
#'  \code{min_info*} \tab type of input data leading to the smallest effect size for the comparison.\cr
#'  \tab \cr
#'  \code{min_es_value*} \tab smallest effect size value for the comparison.\cr
#'  \tab \cr
#'  \code{min_es_se*} \tab standard error of the smallest effect size for the comparison.\cr
#'  \tab \cr
#'  \code{min_es_ci_lo*} \tab lower bound of the 95% CI of the smallest effect size for the comparison.\cr
#'  \tab \cr
#'  \code{min_es_ci_up*} \tab upper bound of the 95% CI of the smallest effect size for the comparison.\cr
#'  \tab \cr
#' }
#'
#' \tabular{ll}{
#'  \code{max_info*} \tab type of input data leading to the largest effect size for the comparison.\cr
#'  \tab \cr
#'  \code{max_es_value*} \tab largest effect size value for the comparison.\cr
#'  \tab \cr
#'  \code{max_es_se*} \tab standard error of the largest effect size for the comparison.\cr
#'  \tab \cr
#'  \code{max_es_ci_lo*} \tab lower bound of the 95% CI of the largest effect size for the comparison.\cr
#'  \tab \cr
#'  \code{max_es_ci_up*} \tab upper bound of the 95% CI of the largest effect size for the comparison.\cr
#'  \tab \cr
#' }
#'
#' \tabular{ll}{
#'  \code{diff_min_max*} \tab difference between the smallest and largest effect size for the comparison.\cr
#'  \tab \cr
#'  \code{overlap_min_max*} \tab % of overlap between the 95% CIs of the largest/smallest effect sizes for the comparison.\cr
#'  \tab \cr
#'  \code{dispersion_es*} \tab standard deviation of all effect sizes for the comparison.\cr
#'  \tab \cr
#' }
#'
#' @seealso
#' \code{\link{metaConvert-package}} for the formatting of well-formatted datasets\cr
#' \code{\link{convert_df}} for estimating effect sizes from a dataset\cr
#'
#' @exportS3Method
#' @export summary.metaConvert
#'
#' @md
#' @examples
#' ### generate a summary of the results of an umbrella object
#' summary(
#'   convert_df(df.haza, measure = "g"),
#'   digits = 5)
summary.metaConvert <- function(object, digits = 3, ...) {
  # object = convert_df(dat, verbose = FALSE,
  #                     or_to_rr = "dipietrantonj", measure="nnt")
  # digits = 3
  # object = convert_df(subset(df.haza, !is.na(n_cases_exp) & !is.na(n_cases_nexp) &
  #                              !is.na(n_controls_exp) & !is.na(n_controls_nexp)),
  #                     verbose = FALSE, hierarchy = "2x2", measure = "nnt")
  # rio::export(data.frame(do.call(rbind, lapply(object, function(x) x[1, "info_used"]))),
  #             "possible_inputs.xlsx", overwrite = TRUE)
  raw_data <- attr(object, "raw_data")
  hierarchy <- attr(object, "hierarchy")
  exp <- attr(object, "exp")
  measure <- attr(object, "measure")
  split_adjusted <- attr(object, "split_adjusted")
  es_selected <- attr(object, "es_selected")
  format <- attr(object, "format_adjusted")
  main_es <- attr(object, "main_es")

  # add missing columns to each dataset (used to easily detect information leading to an ES)
  list_df_es_enh <- lapply(object, .add_columns, y = c(
    "d", "d_se",
    "g", "g_se",
    "md", "md_se",
    "r", "r_se",
    "z", "z_se",
    "logor", "logor_se",
    "logrr", "logrr_se",
    "logirr", "logirr_se",
    "logcvr", "logcvr_se",
    "logvr", "logvr_se",
    "nnt"
  ))

  # extract the values for the correct effect measure
  df_es <- lapply(object, .extract_es, measure = measure, exp = exp)

  # name of all information that could lead to an ES
  list_order <- as.character(sapply(df_es, function(x) x$info_used[1]))
  # rio::export(list_order, "list_possible_inputs.xlsx")

  # hierarchy indicated by user
  if (es_selected == "auto") {
    ordering_full <- list_order
  } else {
    ordering_raw <- gsub(" ", "", unlist(strsplit(hierarchy, split = ">", fixed = TRUE)))
    ordering_full <- append(ordering_raw, list_order[!list_order %in% ordering_raw])
  }

  # append the user hierarchy with all information

  # warn users if wrong inputs have been indicated in the hierarchy
  if (hierarchy == "hierarchy" & length(ordering_full[!ordering_full %in% list_order] > 0)) {
    stop(paste0(
      "Watch out! You have indicated elements that do not exist in the hierarchy. Please discard and replace the following elements: ",
      paste(ordering_full[!ordering_full %in% list_order], collapse = " + ")
    ))
  }

  # final list of hierarchy properly organized
  ordering_tot <- ordering_full[ordering_full %in% list_order]

  adj_list <- ordering_tot[which(grepl("adj", ordering_tot, fixed = TRUE) |
    grepl("ancova", ordering_tot, fixed = TRUE))]

  ordering_crude <- as.character(ordering_tot[-c(which(ordering_tot %in% adj_list))])
  ordering_adj <- as.character(ordering_tot[which(ordering_tot %in% adj_list)])

  # -----------------------------------------------------------------------

  if (split_adjusted == TRUE & format == "wide" & main_es == TRUE) {
    res1 <- .generate_df(
      x = raw_data, list_df = df_es, ordering = ordering_crude, exp = exp,
      digits = digits, suffix = "_crude", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh
    )
    # x = raw_data; list_df = df_es; ordering = ordering_crude; exp = exp;
    # digits = 3; suffix = "_crude"; measure = measure; main_es=TRUE
    res <- .generate_df(
      x = res1, list_df = df_es, ordering = ordering_adj, exp = exp,
      digits = digits, suffix = "_adjusted", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh
    )
    # x = res1; list_df = df_es; ordering = ordering_adj; exp = exp;
    # digits = 3; suffix = "_adjusted"; measure = measure
  } else if (split_adjusted == TRUE & format == "long" & main_es == TRUE) {
    res1 <- .generate_df(
      x = raw_data, list_df = df_es, ordering = ordering_crude, exp = exp,
      digits = digits, suffix = "", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh,
      main_es = TRUE
    )
    res1$adjusted_input = FALSE
    # x = raw_data; list_df = df_es; ordering = ordering_crude; exp = exp;
    # digits = 3; suffix = ""; measure = measure; es_selected = "hierarchy"; list_df_es_enh = list_df_es_enh

    res2 <- .generate_df(
      x = raw_data, list_df = df_es, ordering = ordering_adj, exp = exp,
      digits = digits, suffix = "", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh,
      main_es = TRUE
    )
    res2$adjusted_input = TRUE

    res1_sub = subset(res1, !is.na(res1$measure))
    res2_sub = subset(res2, !is.na(res2$measure))
    res_transit <- rbind(res1_sub, res2_sub)
    res <- res_transit[order(res_transit$row_id), ]

  } else if (split_adjusted == FALSE & main_es == TRUE) {
    res <- .generate_df(
      x = raw_data, list_df = df_es, ordering = ordering_tot, exp = exp,
      digits = digits, suffix = "", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh
    )

    res = subset(res, select = -c(adjusted_input))
  } else if (main_es == FALSE) {
    res1 <- .generate_df(
      x = raw_data, list_df = df_es, ordering = ordering_crude, exp = exp,
      digits = digits, suffix = "", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh,
      main_es = FALSE
    )
    res1$adjusted_input = FALSE
    # x = raw_data; list_df = df_es; ordering = ordering_crude; exp = exp;
    # digits = 3; suffix = ""; measure = measure; es_selected = "hierarchy"; list_df_es_enh = list_df_es_enh

    res2 <- .generate_df(
      x = raw_data, list_df = df_es, ordering = ordering_adj, exp = exp,
      digits = digits, suffix = "", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh,
      main_es = FALSE
    )
    res2$adjusted_input = TRUE

    res1_sub = subset(res1, !is.na(res1$measure))
    res2_sub = subset(res2, !is.na(res2$measure))
    res_transit <- rbind(res1_sub, res2_sub)
    res <- res_transit[order(res_transit$row_id), ]
    res$es_selected = "no selection"
  } else {
    stop("The combination of 'main_es', 'format_adjusted' and 'split_adjusted' is incorrect. Check documentation for more info.")
  }

  if (all(is.na(res$author))) res = subset(res, select = -c(author))
  if (all(is.na(res$year))) res = subset(res, select = -c(year))
  if (all(is.na(res$predictor))) res = subset(res, select = -c(predictor))
  if (all(is.na(res$outcome))) res = subset(res, select = -c(outcome))
  if (all(is.na(res$info_expected))) res = subset(res, select = -c(info_expected))

  main_cols <- which(colnames(res) %in% c(
    "es", "se", "es_ci_lo", "es_ci_up",
    "es_crude", "se_crude", "es_ci_lo_crude", "es_ci_up_crude",
    "es_adjusted", "se_adjusted", "es_ci_lo_adjusted", "es_ci_up_adjusted"
  ))

  res[, main_cols] <- lapply(res[, main_cols], function(x) as.numeric(as.character(x)))
  return(res)
}


#' Print a summary of an object of class \dQuote{metaConvert}
#'
#' @param x an object of class \dQuote{metaConvert}
#' @param ... other arguments that can be passed to the function
#'
#' @details
#' Summary method for objects of class \dQuote{metaConvert}.
#'
#' @return
#' Implicitly calls the \code{\link{summary.metaConvert}} function.
#'
#' @export
#'
#' @md
#'
#' @seealso
#' \code{\link{summary.metaConvert}}
#'
#' @examples
#' ### print the results of an object of class metaConvert
#' convert_df(df.haza, measure = "g")
print.metaConvert <- function(x, ...) {
  y <- summary.metaConvert(x, digits = 3, ...)
  print(y)
}
