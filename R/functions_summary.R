#' Keep one placeholder row for comparisons no estimation route could reach.
#'
#' The long-format branches build the crude and adjusted pools with two separate
#' \code{.generate_df()} calls and then drop every row whose \code{measure} is NA.
#' That silently deleted whole comparisons from the output: the row count no
#' longer matched the input, the console banner reported a vacuous
#' "ES estimated: n/n (100%)", and -- because \code{.add_es_guidance()} only
#' writes where \code{es} is NA -- the entire \code{es_guidance} feature became
#' unreachable in these modes. Re-attach the untouched crude row for any
#' comparison that survived in neither pool.
#'
#' @noRd
.retain_unestimable <- function(res_transit, res_all) {
  missing_ids <- setdiff(res_all$row_id, res_transit$row_id)
  if (length(missing_ids) == 0) return(res_transit)
  filler <- res_all[res_all$row_id %in% missing_ids, , drop = FALSE]
  rbind(res_transit, filler)
}

#' @noRd
.compact_measure_label <- function(measure, exp) {
  if (exp && measure == "logor") return("OR")
  if (exp && measure == "logrr") return("RR")
  if (exp && measure == "logirr") return("IRR")
  if (exp && measure == "loghr") return("HR")
  labels <- c(d = "d", g = "g", md = "MD", r = "r", z = "z",
              logor = "logOR", logrr = "logRR", logirr = "logIRR",
              loghr = "logHR",
              logvr = "logVR", logcvr = "logCVR",
              nnt = "NNT", rd = "RD", dw = "d_w", gw = "g_w", mdw = "MD_w",
              rp = "r_p", zp = "z_p",
              prop = "prop", alpha = "alpha", icc = "ICC")
  if (measure %in% names(labels)) labels[measure] else measure
}

.make_es_summary <- function(res, measure, exp, digits, suffix) {
  es_col <- paste0("es", suffix)
  ci_lo_col <- paste0("es_ci_lo", suffix)
  ci_up_col <- paste0("es_ci_up", suffix)
  out_col <- paste0("es_summary", suffix)
  label <- .compact_measure_label(measure, exp)

  res[[out_col]] <- ifelse(
    is.na(res[[es_col]]), NA_character_,
    paste0(label, " = ",
           round(as.numeric(as.character(res[[es_col]])), digits),
           " [",
           round(as.numeric(as.character(res[[ci_lo_col]])), digits),
           ", ",
           round(as.numeric(as.character(res[[ci_up_col]])), digits),
           "]")
  )
  res
}

.make_es_consistency <- function(res, digits, suffix) {
  es_col <- paste0("es", suffix)
  n_est_col <- paste0("n_estimations", suffix)
  min_info_col <- paste0("min_info", suffix)
  min_val_col <- paste0("min_es_value", suffix)
  max_info_col <- paste0("max_info", suffix)
  max_val_col <- paste0("max_es_value", suffix)
  overlap_col <- paste0("overlap_min_max", suffix)
  disp_col <- paste0("dispersion_es", suffix)
  out_col <- paste0("es_consistency", suffix)

  n_est <- suppressWarnings(as.numeric(as.character(res[[n_est_col]])))
  overlap_raw <- suppressWarnings(as.numeric(as.character(res[[overlap_col]])))
  disp_raw <- suppressWarnings(as.numeric(as.character(res[[disp_col]])))
  min_val <- suppressWarnings(as.numeric(as.character(res[[min_val_col]])))
  max_val <- suppressWarnings(as.numeric(as.character(res[[max_val_col]])))

  res[[out_col]] <- ifelse(
    is.na(res[[es_col]]), NA_character_,
    ifelse(!is.na(n_est) & n_est > 1,
      paste0("min: ", res[[min_info_col]], " = ", round(min_val, digits),
             ", max: ", res[[max_info_col]], " = ", round(max_val, digits),
             ", overlap: ", round(overlap_raw * 100), "%",
             # 'dispersion_es' is .dispersion_stat() = max|v - median(v)|, NOT a
             # sample SD (see the rationale at internal_generate_df.R). Label it
             # for what it is -- this string is the one place a user reads it.
             ", max dev: ", round(disp_raw, digits)),
      "")
  )
  res
}

# round ES/SE/CI columns (flags need full precision)
.round_numeric_cols <- function(res, digits, suffix) {
  cols_to_round <- paste0(
    c("es", "se", "es_ci_lo", "es_ci_up",
      "overlap_min_max", "diff_min_max", "dispersion_es",
      "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
      "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up"),
    suffix
  )
  for (col in cols_to_round) {
    if (!col %in% colnames(res)) next
    vals <- res[[col]]
    numeric_rows <- which(vals != "< 2 types of input data available")
    if (length(numeric_rows) > 0) {
      res[numeric_rows, col] <- round(
        as.numeric(as.character(vals[numeric_rows])), digits)
    }
  }
  res
}

#' Synthesize information of an object of class \dQuote{metaConvert} into a dataframe
#'
#' @param object an object of class \dQuote{metaConvert}
#' @param digits an integer value specifying the number of decimal places for the rounding of numeric values. Default is 3.
#' @param flags a logical value indicating whether quality/plausibility flags should be generated. Default is TRUE. The \code{flags} column (\code{flags_crude} / \code{flags_adjusted} when \code{split_adjusted = TRUE} and \code{format = "wide"}) includes both input validation flags (generated during \code{\link{convert_df}}) and post-computation quality flags.
#' @param flag_options a named list of thresholds overriding the defaults (and any options set in \code{\link{convert_df}}).
#' Available options:
#' \itemize{
#'   \item \code{smd_max} (default 3): |SMD| above this value is flagged
#'   \item \code{r_max} (default 0.95): same, for |r|
#'   \item \code{log_or_max} (default 5): same, for |logOR| and |logRR|
#'   \item \code{n_min} (default 10): flag sample sizes below this value
#'   \item \code{iqr_mult} (default 3): IQR multiplier for the cross-row outlier detection
#'   \item \code{dispersion_max_smd} / \code{dispersion_max_r} / \code{dispersion_max_logor} (defaults 0.5 / 0.15 / 1.0): maximum absolute deviation of the per-method estimates from their median, for the SMD family (d, g, dw, gw, md, mdw), for r, and for the log-ratio family (logOR, logRR, logIRR, logHR) respectively. Other measures use the SMD value
#'   \item \code{diff_max_smd} / \code{diff_max_r} / \code{diff_max_logor} (defaults 1.0 / 0.3 / 2.0): maximum min-max ES difference across estimation methods
#'   \item \code{overlap_min} (default 0.80): minimum CI overlap between the min/max estimates (0-1 scale)
#'   \item \code{enable_cross_row} (default TRUE): enable/disable the cross-row checks
#'   \item \code{direction_conflict_min} (default 2): number of significantly-positive and significantly-negative studies needed to raise the direction conflict flag (G1)
#'   \item \code{flag_group} (default NULL): one or more input-column names (e.g. \code{"outcome"} or \code{c("outcome", "subgroup")}) used to scope the CROSS-ROW checks (ES/SE/SD outliers, direction conflict, study-duplication, NNT-type and standardizer mixing). With the default \code{NULL} the whole dataset is one pool. Set it for multivariate / multi-outcome data, where a deviation or a repeated \code{study_id} is only meaningful WITHIN a group of comparable rows; a study contributing several outcomes is then no longer flagged as a duplicate, and outliers are judged against same-group peers. Per-row and cross-method checks, and the byte-identical templated-data check (V23), are unaffected.
#' }
#' @param guidance a logical value indicating whether missing data guidance should be generated for rows where the effect size is NA. Default is TRUE. When enabled, a column \code{es_guidance} (or \code{es_guidance_crude}/\code{es_guidance_adjusted}) is appended, listing the closest estimation methods and which specific columns are missing.
#' @param include_raw a logical value indicating whether the raw input columns should be appended after the effect size columns in the returned dataframe. Default is TRUE.
#' @param formulas a logical value indicating whether the sensitivity of each effect size to the choice of conversion formula (\code{or_to_rr}, \code{cor_to_smd}, \code{pre_post_to_smd}, ...) should be reported. Default is FALSE. When TRUE, a concise \code{[FORMULA]} message is appended to the flags column for every comparison whose estimate depends on such a choice, and the complete table of alternative estimates (as returned by \code{\link{es_formulas}}) is attached as the \code{"formulas"} attribute. \code{\link{es_formulas}} re-estimates one effect size per comparison, so the message is appended only where the reported \code{info_used} is the type of input data that estimate was derived from: an estimate obtained from other input data is left unannotated, even when it belongs to the same comparison. This reports methodological uncertainty rather than a data quality problem: unlike the \code{[DISCORDANT]} flags, a discrepancy between conversion formulae reflects differing statistical assumptions rather than an extraction error. Each alternative formula requires one additional evaluation of \code{\link{convert_df}}.
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
#'  \code{dispersion_es*} \tab maximum absolute deviation of the effect sizes from their median for the comparison (used by the cross-method discordance check).\cr
#'  \tab \cr
#' }
#'
#' **5. Quality/plausibility flags**
#' When \code{flags = TRUE} (the default), a flags* column is added, containing
#' the semicolon-separated list of issues detected for the row (empty when none).
#'
#' \tabular{ll}{
#'  \code{flags*} \tab quality/plausibility flags for the effect size in this row.\cr
#'  \tab \cr
#' }
#'
#' Flag categories:
#' \itemize{
#'   \item \strong{A (Numeric integrity)}: Inf/NaN values, negative SE, inverted CI, ES outside CI
#'   \item \strong{B (Bounds violations)}: correlation outside \eqn{[-1,1]}, non-positive OR/RR, NNT in (-1,1), proportion outside \eqn{[0,1]}, RD outside \eqn{[-1,1]}
#'   \item \strong{C (Plausibility)}: unusually large SMD, zero SE, high correlation, large logOR/logRR, small sample size
#'   \item \strong{D (Cross-row outliers)}: ES or SE is an outlier relative to other rows (IQR method)
#'   \item \strong{E (Internal consistency)}: high dispersion across estimation methods, low/zero CI overlap, large min-max difference
#' }
#'
#' **6. Missing data guidance**
#' When \code{guidance = TRUE} (the default), an es_guidance* column indicates,
#' for rows without an effect size, the closest estimation methods and which
#' input columns are missing (empty otherwise).
#'
#' \tabular{ll}{
#'  \code{es_guidance*} \tab guidance on which columns to add to obtain an effect size.\cr
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
#' ### generate a summary of the results of a metaConvert object
#' summary(
#'   convert_df(df.haza, measure = "g"),
#'   digits = 5)
summary.metaConvert <- function(object, digits = 3, flags = TRUE, flag_options = list(), guidance = TRUE, include_raw = TRUE, formulas = FALSE, ...) {
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
    "rp", "rp_se",
    "zp", "zp_se",
    "nnt", "nnt_se", "nnt_ci_lo", "nnt_ci_up",
    "rd", "rd_se", "rd_ci_lo", "rd_ci_up",
    "prop", "prop_se", "prop_ci_lo", "prop_ci_up",
    "alpha", "alpha_se", "alpha_ci_lo", "alpha_ci_up",
    "icc", "icc_se", "icc_ci_lo", "icc_ci_up"
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

  idx_adj <- which(ordering_tot %in% adj_list)
  ordering_crude <- if (length(idx_adj) > 0) as.character(ordering_tot[-idx_adj]) else as.character(ordering_tot)
  ordering_adj <- as.character(ordering_tot[idx_adj])

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
    res_transit <- .retain_unestimable(rbind(res1_sub, res2_sub), res1)
    res <- res_transit[order(res_transit$row_id), ]

  } else if (split_adjusted == FALSE & main_es == TRUE) {
    res <- .generate_df(
      x = raw_data, list_df = df_es, ordering = ordering_tot, exp = exp,
      digits = digits, suffix = "", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh
    )

    res = subset(res, select = -c(adjusted_input))
  } else if (main_es == FALSE & split_adjusted == FALSE) {
    # split_adjusted = FALSE was previously ignored in the route view: crude and
    # adjusted routes were always built as two separate pools, so a comparison
    # offering both never had them compared against each other. Honour it here by
    # running the routes through a single ordering, as main_es = TRUE does.
    res <- .generate_df(
      x = raw_data, list_df = df_es, ordering = ordering_tot, exp = exp,
      digits = digits, suffix = "", measure = measure,
      es_selected = es_selected, list_df_es_enh = list_df_es_enh,
      main_es = FALSE
    )
    res <- .retain_unestimable(subset(res, !is.na(res$measure)), res)
    res <- res[order(res$row_id), ]
    res = subset(res, select = -c(adjusted_input))
    res$es_selected = rep("no selection", nrow(res))
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
    res_transit <- .retain_unestimable(rbind(res1_sub, res2_sub), res1)
    res <- res_transit[order(res_transit$row_id), ]
    # rep() rather than a scalar: assigning a length-1 value to a 0-row
    # data.frame is an error, and res can legitimately be empty when the input
    # contains no estimable comparison at all.
    res$es_selected = rep("no selection", nrow(res))
  } else {
    stop("The combination of 'main_es', 'format_adjusted' and 'split_adjusted' is incorrect. Check documentation for more info.")
  }

  # compact summary columns
  if (split_adjusted == TRUE & format == "wide" & main_es == TRUE) {
    res <- .make_es_summary(res, measure, exp, digits, "_crude")
    res <- .make_es_consistency(res, digits, "_crude")
    res <- .make_es_summary(res, measure, exp, digits, "_adjusted")
    res <- .make_es_consistency(res, digits, "_adjusted")
  } else {
    res <- .make_es_summary(res, measure, exp, digits, "")
    res <- .make_es_consistency(res, digits, "")
  }

  # quality flags
  if (flags) {
    stored_opts <- attr(object, "flag_options")
    if (is.null(stored_opts)) stored_opts <- .default_flag_options()
    opts <- stored_opts
    opts[names(flag_options)] <- flag_options
    input_val <- attr(object, "input_validation")
    r_def <- attr(object, "r_defaulted")
    alpha_method <- attr(object, "alpha_to_es")
    if (is.null(alpha_method)) alpha_method <- "bonett"
    icc_method <- attr(object, "icc_to_es")
    if (is.null(icc_method)) icc_method <- "bonett"
    prop_method <- attr(object, "prop_to_es")
    if (is.null(prop_method)) prop_method <- "raw"
    pp_method <- attr(object, "pre_post_to_smd")
    if (is.null(pp_method)) pp_method <- "bonett"
    pool_sd_used <- attr(object, "pool_sd")
    if (is.null(pool_sd_used)) pool_sd_used <- FALSE
    smd_denom_used <- attr(object, "smd_denom_used")
    # Which transform each route put in the z column (E8). Derived in convert_df()
    # from the frames themselves; see .z_transform_by_route().
    z_transform_used <- attr(object, "z_transform")

    # info_used values ranked by the hierarchy actually in force -- the same
    # ordering_crude / ordering_adj this function already uses to select the reported
    # effect size. The cross-row checks need it to pick, as a comparison's
    # representative, the route that would actually be selected (see .flag_es_quality).
    # NB the element order of the metaConvert list is fixed and does NOT track the
    # hierarchy, so it cannot be used for this.

    if (split_adjusted == TRUE & format == "wide" & main_es == TRUE) {
      res <- .flag_es_quality(res, measure, exp, "_crude", raw_data, opts, input_val,
                              alpha_to_es = alpha_method, icc_to_es = icc_method,
                              prop_to_es = prop_method, pre_post_to_smd = pp_method,
                              pool_sd = pool_sd_used, r_defaulted = r_def,
                              smd_denom = smd_denom_used, es_order = ordering_crude,
                              z_transform = z_transform_used)
      res <- .flag_es_quality(res, measure, exp, "_adjusted", raw_data, opts, input_val,
                              alpha_to_es = alpha_method, icc_to_es = icc_method,
                              prop_to_es = prop_method, pre_post_to_smd = pp_method,
                              pool_sd = pool_sd_used, r_defaulted = r_def,
                              smd_denom = smd_denom_used, es_order = ordering_adj,
                              z_transform = z_transform_used)
    } else {
      res <- .flag_es_quality(res, measure, exp, "", raw_data, opts, input_val,
                              alpha_to_es = alpha_method, icc_to_es = icc_method,
                              prop_to_es = prop_method, pre_post_to_smd = pp_method,
                              pool_sd = pool_sd_used, r_defaulted = r_def,
                              smd_denom = smd_denom_used, es_order = ordering_tot,
                              z_transform = z_transform_used)
    }
  }

  # Sensitivity to the conversion formula: append a concise [FORMULA] message to
  # the flags column and attach the complete table of alternative estimates.
  # Disabled by default, since each alternative formula requires one additional
  # evaluation of convert_df(), and since this reports methodological uncertainty
  # rather than a data quality problem: it must not appear unrequested alongside
  # the metaDETECT flags.
  fx_store <- NULL
  if (isTRUE(formulas)) {
    fx <- try(es_formulas(object, digits = digits, verbose = FALSE), silent = TRUE)
    if (!inherits(fx, "try-error") && nrow(fx) > 0) {
      tokens <- .formula_tokens(fx, digits)
      if (!is.null(tokens)) {
        # es_formulas() re-runs convert_df() with main_es = TRUE and
        # split_adjusted = FALSE, so a token describes one estimate only: the
        # one derived from the type of input data named in fx$info_used. Attach
        # it where that input data is the one reported, or the disclosure lands
        # on estimates it does not describe -- the adjusted column of a split
        # run (whose estimate can be invariant to the parameter quoted, or
        # missing altogether), or a route the parameter cannot reach.
        described <- unique(paste(fx$row_id, as.character(fx$info_used), sep = "\r"))
        for (sfx in if (split_adjusted == TRUE & format == "wide" & main_es == TRUE)
                      c("_crude", "_adjusted") else "") {
          fcol <- paste0("flags", sfx)
          icol <- paste0("info_used", sfx)
          if (!fcol %in% colnames(res) || !icol %in% colnames(res)) next
          info_row <- as.character(res[[icol]])
          tok <- tokens[as.character(res$row_id)]
          add <- !is.na(tok) & !is.na(info_row) &
            paste(res$row_id, info_row, sep = "\r") %in% described
          cur <- res[[fcol]]
          cur[is.na(cur)] <- ""
          res[[fcol]][add] <- ifelse(nchar(cur[add]) > 0,
                                     paste0(cur[add], "; ", tok[add]), tok[add])
        }
      }
      # attached at the very end: `[.data.frame` drops non-standard attributes,
      # and res is still reordered/subset below.
      fx_store <- fx
    }
  }

  # missing data guidance
  if (guidance) {
    if (split_adjusted == TRUE & format == "wide" & main_es == TRUE) {
      res <- .add_es_guidance(res, "_crude", raw_data, measure, object)
      res <- .add_es_guidance(res, "_adjusted", raw_data, measure, object)
    } else {
      res <- .add_es_guidance(res, "", raw_data, measure, object)
    }
  }

  # rounding (flags above use full precision)
  if (split_adjusted == TRUE & format == "wide" & main_es == TRUE) {
    res <- .round_numeric_cols(res, digits, "_crude")
    res <- .round_numeric_cols(res, digits, "_adjusted")
  } else {
    res <- .round_numeric_cols(res, digits, "")
  }

  # reorder columns
  all_cols <- colnames(res)

  identity_names   <- c("row_id", "study_id", "author", "year", "outcome",
                        "predictor", "info_expected", "adjusted_input")
  provenance_names <- c("all_info", "info_measure", "es_guidance", "es_selected")
  quality_names    <- c("flags")
  summary_names    <- c("es_summary", "es_consistency", "n_estimations",
                        "dispersion_es")
  primary_names    <- c("es", "se", "es_ci_lo", "es_ci_up",
                        "info_used", "measure")
  minmax_names     <- c("overlap_min_max", "diff_min_max",
                        "min_info", "max_info",
                        "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
                        "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up")

  if (split_adjusted == TRUE & format == "wide" & main_es == TRUE) {
    suffixes <- c("_crude", "_adjusted")
  } else {
    suffixes <- ""
  }

  ordered <- character(0)
  ordered <- c(ordered, identity_names[identity_names %in% all_cols])
  for (sfx in suffixes) {
    for (block in list(provenance_names, quality_names, summary_names,
                       primary_names, minmax_names)) {
      cols <- paste0(block, sfx)
      ordered <- c(ordered, cols[cols %in% all_cols])
    }
  }
  remaining <- all_cols[!all_cols %in% ordered]
  if (isTRUE(include_raw)) {
    ordered <- c(ordered, remaining)
  }
  res <- res[, ordered, drop = FALSE]

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

  # diagnostic overview
  .measure_label <- function(m) {
    labels <- c(d = "Cohen's d", g = "Hedges' g", md = "Mean difference",
                r = "Pearson r", z = "Fisher's z", or = "Odds ratio",
                rr = "Risk ratio", irr = "Incidence rate ratio",
                logor = "log(OR)", logrr = "log(RR)", logirr = "log(IRR)",
                logvr = "log(VR)", logcvr = "log(CVR)",
                nnt = "NNT", rd = "Risk difference",
                dw = "Within-group d", gw = "Within-group g",
                mdw = "Within-group MD",
                rp = "Partial r", zp = "Fisher's z of partial r",
                prop = "Proportion", alpha = "Cronbach's alpha",
                icc = "ICC")
    if (m %in% names(labels)) labels[m] else m
  }

  .summary_for_suffix <- function(res, suffix) {
    es_col <- paste0("es", suffix)
    info_col <- paste0("info_used", suffix)
    flags_col <- paste0("flags", suffix)
    if (!es_col %in% colnames(res)) return(NULL)

    n_total <- nrow(res)
    es_vals <- res[[es_col]]
    n_estimated <- sum(!is.na(es_vals))
    # floor, not round: 288/289 rounds up to a self-contradicting "100%"
    pct <- floor(100 * n_estimated / n_total)

    info_vals <- res[[info_col]]
    info_vals <- info_vals[!is.na(info_vals) & info_vals != ""]
    method_freq <- sort(table(info_vals), decreasing = TRUE)

    n_flags <- 0
    if (flags_col %in% colnames(res)) {
      flag_vals <- res[[flags_col]]
      n_flags <- sum(!is.na(flag_vals) & nchar(flag_vals) > 0)
    }

    # the input row (row_id), not the position in res: the route view and the
    # long format hold several rows per comparison, so a position identifies no
    # study the reader can look up
    missing_rows <- unique(res$row_id[is.na(es_vals)])

    list(n_total = n_total, n_estimated = n_estimated, pct = pct,
         method_freq = method_freq, n_flags = n_flags,
         missing_rows = missing_rows)
  }

  is_wide_split <- split_adjusted == TRUE & format == "wide" & main_es == TRUE

  if (is_wide_split) {
    info_crude <- .summary_for_suffix(res, "_crude")
    info_adj   <- .summary_for_suffix(res, "_adjusted")

    message("\n-- metaConvert summary --")
    message("Measure: ", .measure_label(measure), "  |  ", info_crude$n_total, " studies")

    # Crude
    message("\nCrude estimates: ", info_crude$n_estimated, "/", info_crude$n_total,
            " (", info_crude$pct, "%)")
    if (length(info_crude$method_freq) > 0) {
      top <- utils::head(info_crude$method_freq, 5)
      for (i in seq_along(top)) {
        message("  ", format(names(top)[i], width = 25), " ", top[i])
      }
    }
    if (info_crude$n_flags > 0) message(info_crude$n_flags, " quality flag(s) raised")
    if (length(info_crude$missing_rows) > 0 && length(info_crude$missing_rows) <= 10) {
      message("Missing ES in input rows: ", paste(info_crude$missing_rows, collapse = ", "))
    } else if (length(info_crude$missing_rows) > 10) {
      message(length(info_crude$missing_rows), " input rows with missing ES")
    }

    # Adjusted
    if (!is.null(info_adj)) {
      message("\nAdjusted estimates: ", info_adj$n_estimated, "/", info_adj$n_total,
              " (", info_adj$pct, "%)")
      if (length(info_adj$method_freq) > 0) {
        top <- utils::head(info_adj$method_freq, 5)
        for (i in seq_along(top)) {
          message("  ", format(names(top)[i], width = 25), " ", top[i])
        }
      }
      if (info_adj$n_flags > 0) message(info_adj$n_flags, " quality flag(s) raised")
    }
  } else {
    info <- .summary_for_suffix(res, "")
    if (!is.null(info)) {
      # In the route view a comparison occupies one row per estimation route,
      # and in the long format one row per populated scope, so nrow(res) counts
      # neither in studies. Report both, and do not claim a method was
      # "selected" when main_es = FALSE selects nothing.
      n_cmp <- length(unique(res$row_id))
      route_view <- isFALSE(main_es)
      message("\n-- metaConvert summary --")
      if (route_view || n_cmp != info$n_total) {
        unit <- if (route_view) "estimation routes" else "crude/adjusted rows"
        message("Measure: ", .measure_label(measure), "  |  ", n_cmp,
                " comparisons  |  ", info$n_total, " ", unit)
      } else {
        message("Measure: ", .measure_label(measure), "  |  ", info$n_total, " studies")
      }
      message("ES estimated: ", info$n_estimated, "/", info$n_total,
              " (", info$pct, "%)")
      if (length(info$method_freq) > 0) {
        message(if (route_view) "\nRoutes available:" else "\nMethods selected:")
        top <- utils::head(info$method_freq, 5)
        for (i in seq_along(top)) {
          message("  ", format(names(top)[i], width = 25), " ", top[i])
        }
      }
      if (info$n_flags > 0) message(info$n_flags, " quality flag(s) raised")
      if (length(info$missing_rows) > 0 && length(info$missing_rows) <= 10) {
        message("Missing ES in input rows: ", paste(info$missing_rows, collapse = ", "))
      } else if (length(info$missing_rows) > 10) {
        message(length(info$missing_rows), " input rows with missing ES")
      }
    }
  }

  if (!is.null(fx_store)) attr(res, "formulas") <- fx_store

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
#' \code{x}, invisibly. Called for its side effect: the summary produced by
#' \code{\link{summary.metaConvert}} is printed.
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
  invisible(x)
}

#' Convert a \dQuote{metaConvert} object to a dataframe
#'
#' @param x an object of class \dQuote{metaConvert}
#' @param ... other arguments passed to \code{\link{summary.metaConvert}}
#'
#' @details
#' Convenience wrapper around \code{summary.metaConvert(x, ...)}.
#'
#' @return
#' A dataframe containing the effect size results, diagnostics, and all input data columns.
#'
#' @seealso
#' \code{\link{summary.metaConvert}}
#'
#' @exportS3Method
#'
#' @md
#'
#' @examples
#' ### get the full output as a dataframe
#' as.data.frame(convert_df(df.haza, measure = "g"))
as.data.frame.metaConvert <- function(x, ...) {
  summary.metaConvert(x, ...)
}
