#' default flag options
#' @noRd
.default_flag_options <- function() {
  list(
    smd_max = 3, # C1
    r_max = 0.95, # C4
    log_or_max = 5, # C5
    n_min = 10, # C6
    iqr_mult = 3, # D1/D2
    dispersion_max_smd = 0.5, # E1
    dispersion_max_r = 0.15, # E1
    dispersion_max_logor = 1.0, # E1
    diff_max_smd = 1.0, # E3
    diff_max_r = 0.3, # E3
    diff_max_logor = 2.0, # E3
    overlap_min = 0.85, # E2b
    enable_cross_row = TRUE, # D1/D2
    # One or more raw-data column names to scope the CROSS-ROW checks by (D1/D2/D3,
    # G, H, E4/E6/E7). NOT V23, which is a Tier-1 input-data property computed in
    # .validate_input_data() and is deliberately left dataset-wide.
    # NULL (default) = whole dataset is one pool, as before.
    # Set it (e.g. "outcome", or c("outcome","subgroup")) for multivariate / multi-
    # outcome extraction sheets, where deviations are only meaningful WITHIN a group
    # of comparable rows. Per-row and cross-method checks are unaffected.
    flag_group = NULL,
    alpha_max = 0.99, # C7
    icc_max = 0.99, # C8
    group_ratio_max = 10, # C10
    direction_conflict_min = 2, # G1
    direction_conflict_pct = 10, # G1
    direction_conflict_lone_pct = 20, # G1, Inf to disable the lone tier
    direction_conflict_info = TRUE, # G1
    direction_disagreement_min = 0.1, # E5, unit-dependent for md/mdw
    enable_informational = FALSE, # V8
    ignore_sign_reversal = FALSE, # E1-E3
    outlier_min_deviation = NULL, # D1
    se_outlier_min_ratio = 3, # D2
    se_ratio_extreme = 10, # D2b
    sd_ratio_max = 10, # V12
    sd_ratio_bl_ep_min = 0.7, # V18
    baseline_imbalance_max = 0.30, # V21
    templated_min_match = 3, # V23
    paired_as_indep_tol = 0.07, # V26
    sd_outlier_ratio = 5, # V27/D3
    sd_outlier_smd = FALSE # V27/D3, off: smd pools mix instrument scales
  )
}


#' Build a per-row grouping key for the CROSS-ROW checks from opts$flag_group.
#'
#' Returns a length-n character vector used to scope D1/D2/D3, G, H and E4/E6/E7.
#' V23 is NOT scoped: it is a Tier-1 input-data property evaluated in
#' .validate_input_data(), which never sees this key.
#' When flag_group is NULL / absent / names no present column, every row gets
#' the same key ("__all__") -- i.e. one pool, the pre-existing whole-dataset
#' behaviour. Rows with an NA in a grouping column fall into a shared "<NA>" level.
#' The key doubles as the human-readable group label shown in flag messages.
#'
#' @param data data.frame the group column(s) are read from (the raw input data)
#' @param flag_group NULL, or one or more column names in `data`
#' @param idx integer row map (length n) from result rows to rows of `data`; NULL
#'   means `data` is already row-aligned
#' @param n number of rows the key must cover (defaults to length(idx) / nrow(data))
#' @noRd
.build_group_key <- function(data, flag_group, idx = NULL, n = NULL) {
  if (is.null(n)) n <- if (!is.null(idx)) length(idx) else nrow(data)
  if (is.null(flag_group) || length(flag_group) == 0) return(rep("__all__", n))
  cols <- intersect(as.character(flag_group), colnames(data))
  if (length(cols) == 0) return(rep("__all__", n))
  parts <- lapply(cols, function(cn) {
    v <- as.character(data[[cn]])
    if (!is.null(idx)) v <- v[idx]
    v[is.na(v)] <- "<NA>"
    v
  })
  do.call(paste, c(parts, sep = " / "))
}


#' Split rows by group_key, run a cross-row flag function within each group, and
#' reassemble a length-n list of flag vectors (global row order preserved).
#'
#' `fn` receives a vector of GLOBAL row indices for one group and must return a
#' list of character vectors the same length, in that index order. When more than
#' one group is present, each emitted flag gets a " (within group '<g>')" suffix
#' (no "; ", so the downstream flag merge is unaffected) so the reviewer knows the
#' comparison set. Requires self-contained messages (no absolute row-number
#' references) -- true for the D-family and G checks.
#'
#' @param group_key length-n character vector from .build_group_key()
#' @param n total number of rows
#' @param fn function(idx) -> list of length(idx) character vectors
#' @noRd
.by_group <- function(group_key, n, fn) {
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)
  groups <- split(seq_len(n), group_key)
  multi <- length(groups) > 1L
  # Iterate positionally, NOT over names(groups): a group whose key is the empty
  # string cannot be reached with groups[[""]] (it returns NULL) and its rows would
  # silently lose every cross-row flag.
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    sub <- fn(idx)
    for (j in seq_along(idx)) {
      msgs <- sub[[j]]
      if (length(msgs) && multi) {
        msgs <- paste0(msgs, " (within group '", .group_label(names(groups)[k]), "')")
      }
      flags[[idx[j]]] <- msgs
    }
  }
  flags
}


#' Make a grouping-column value safe to interpolate into a flag message.
#'
#' The merged flag string is split on "; ", so a group label containing that
#' sequence (e.g. an outcome named "HAM-D; total score") would create an untagged
#' token downstream. Replace it with a comma.
#'
#' @param g character scalar, the raw group key
#' @noRd
.group_label <- function(g) {
  gsub("; ", ", ", as.character(g), fixed = TRUE)
}


#' max absolute deviation from the median (used by E1)
#' for two values equals |a - b| / 2, so E1 == E3 at k = 2
#' @param v numeric vector of per-method effect-size estimates for one row
#' @return numeric scalar max(|v - median(v)|), or NA if < 2 finite values
#' @noRd
.dispersion_stat <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2) return(NA_real_)
  max(abs(v - stats::median(v)))
}

#' note for r-sensitive flags (A6, E2, E2b), which depend on a defaulted r_pre_post
#' @param is_def logical scalar (TRUE when this row's r_pre_post was defaulted)
#' @return character scalar (empty string when not r-sensitive)
#' @noRd
.r_sensitive_note <- function(is_def) {
  if (isTRUE(is_def)) {
    # no "; " here, the flag merge splits on it
    " (r-sensitive: depends on the default r_pre_post - re-check with the true pre-post correlation)"
  } else {
    ""
  }
}


#########

#' Count the number of decimal places in a numeric value
#'
#' Used by the CI asymmetry check to estimate rounding precision.
#' @param x numeric scalar
#' @return integer, number of decimal places (0 for integers, max 10)
#' @noRd
.count_decimals <- function(x) {
  if (is.na(x) || !is.finite(x)) return(2L)
  s <- format(x, scientific = FALSE, drop0trailing = TRUE)
  if (!grepl("\\.", s)) return(1L)  # integers assumed at least 1dp precision
  dp <- nchar(sub(".*\\.", "", s))
  min(dp, 10L)
}


#' Columns that must be non-negative in the input data
#' @return character vector of column names
#' @noRd
.positive_columns <- function() {
  c(
    # Sample sizes
    "n_exp", "n_nexp", "n_sample", "n_cases", "n_controls",
    "n_cases_exp", "n_cases_nexp", "n_controls_exp", "n_controls_nexp",
    # SDs
    "mean_sd_exp", "mean_sd_nexp", "mean_sd_pooled",
    "mean_pre_sd_exp", "mean_pre_sd_nexp",
    "mean_change_sd_exp", "mean_change_sd_nexp",
    "ancova_mean_sd_exp", "ancova_mean_sd_nexp", "ancova_mean_sd_pooled",
    "md_sd", "ancova_md_sd",
    # SEs
    "mean_se_exp", "mean_se_nexp", "mean_pre_se_exp", "mean_pre_se_nexp",
    "mean_change_se_exp", "mean_change_se_nexp",
    "ancova_mean_se_exp", "ancova_mean_se_nexp",
    "md_se", "ancova_md_se", "logor_se", "logrr_se", "logirr_se",
    "rd_se", "linreg_b_se",
    "user_se_crude", "user_se_adj",
    # Psychometric counts
    "n_items", "n_measurements",
    # Statistics & other
    "anova_f", "chisq", "etasq", "etasq_adj",
    "time_exp", "time_nexp", "baseline_rate"
    # p-values checked in .bounded_columns
  )
}

#' CI triplets to check for consistency in the input data
#'
#' Each triplet has val/lo/up column names and a scale indicator: "additive"
#' (symmetric CIs), "exp" (natural-scale OR/RR, checked on the log scale), or
#' "user" (detected at runtime from the original measure).
#'
#' @return list of lists, each with val/lo/up/scale
#' @noRd
.ci_triplets <- function() {
  list(
    # Means (additive scale)
    list(val = "mean_exp",         lo = "mean_ci_lo_exp",         up = "mean_ci_up_exp",         scale = "additive"),
    list(val = "mean_nexp",        lo = "mean_ci_lo_nexp",        up = "mean_ci_up_nexp",        scale = "additive"),
    list(val = "mean_pre_exp",     lo = "mean_pre_ci_lo_exp",     up = "mean_pre_ci_up_exp",     scale = "additive"),
    list(val = "mean_pre_nexp",    lo = "mean_pre_ci_lo_nexp",    up = "mean_pre_ci_up_nexp",    scale = "additive"),
    list(val = "mean_change_exp",  lo = "mean_change_ci_lo_exp",  up = "mean_change_ci_up_exp",  scale = "additive"),
    list(val = "mean_change_nexp", lo = "mean_change_ci_lo_nexp", up = "mean_change_ci_up_nexp", scale = "additive"),
    list(val = "ancova_mean_exp",  lo = "ancova_mean_ci_lo_exp",  up = "ancova_mean_ci_up_exp",  scale = "additive"),
    list(val = "ancova_mean_nexp", lo = "ancova_mean_ci_lo_nexp", up = "ancova_mean_ci_up_nexp", scale = "additive"),
    # Mean differences (additive scale)
    list(val = "md",               lo = "md_ci_lo",               up = "md_ci_up",               scale = "additive"),
    list(val = "ancova_md",        lo = "ancova_md_ci_lo",        up = "ancova_md_ci_up",        scale = "additive"),
    # User input (scale depends on user_es_original_measure - detected at runtime)
    list(val = "user_es_crude",    lo = "user_ci_lo_crude",       up = "user_ci_up_crude",       scale = "user"),
    list(val = "user_es_adj",      lo = "user_ci_lo_adj",         up = "user_ci_up_adj",         scale = "user"),
    # Log-scale ratios (additive on log scale)
    list(val = "logor",            lo = "logor_ci_lo",            up = "logor_ci_up",            scale = "additive"),
    list(val = "logrr",            lo = "logrr_ci_lo",            up = "logrr_ci_up",            scale = "additive"),
    list(val = "logirr",           lo = "logirr_ci_lo",           up = "logirr_ci_up",           scale = "additive"),
    list(val = "loghr",            lo = "loghr_ci_lo",            up = "loghr_ci_up",            scale = "additive"),
    # Risk difference (additive scale)
    list(val = "rd",        lo = "rd_ci_lo",        up = "rd_ci_up",        scale = "additive"),
    # Regression coefficient (additive scale)
    list(val = "linreg_b",  lo = "linreg_b_ci_lo",  up = "linreg_b_ci_up",  scale = "additive"),
    # Natural-scale ratios (exp scale - inherently asymmetric CIs)
    list(val = "or",               lo = "or_ci_lo",               up = "or_ci_up",               scale = "exp"),
    list(val = "rr",               lo = "rr_ci_lo",               up = "rr_ci_up",               scale = "exp"),
    # Plot-extracted means (additive scale)
    list(val = "plot_mean_exp",         lo = "plot_mean_ci_lo_exp",         up = "plot_mean_ci_up_exp",         scale = "additive"),
    list(val = "plot_mean_nexp",        lo = "plot_mean_ci_lo_nexp",        up = "plot_mean_ci_up_nexp",        scale = "additive"),
    list(val = "plot_ancova_mean_exp",  lo = "plot_ancova_mean_ci_lo_exp",  up = "plot_ancova_mean_ci_up_exp",  scale = "additive"),
    list(val = "plot_ancova_mean_nexp", lo = "plot_ancova_mean_ci_lo_nexp", up = "plot_ancova_mean_ci_up_nexp", scale = "additive")
  )
}

#' Validate input data before effect size computation
#'
#' Sets invalid values to NA and records per-row issues. When
#' \code{correct_inputs = FALSE}, invalid values are flagged but preserved.
#'
#' @param x data.frame, the input data after .check_data()
#' @param max_asymmetry numeric, percentage threshold for CI asymmetry detection (default 10)
#' @param verbose logical, if TRUE emit console messages for detected issues (default TRUE)
#' @param correct_inputs logical, if TRUE (default) invalid values are set to missing;
#'   if FALSE, values are flagged but preserved in the returned data.
#' @return list with $data (cleaned data.frame) and $issues (character vector, one per row)
#' @noRd
.validate_input_data <- function(x, max_asymmetry = 10, verbose = TRUE,
                                  enable_informational = FALSE,
                                  sd_ratio_max = 10,
                                  sd_ratio_bl_ep_min = 0.7,
                                  baseline_imbalance_max = 0.30,
                                  enable_cross_row = TRUE,
                                  templated_min_match = 3,
                                  measure = NULL,
                                  paired_as_indep_tol = 0.07,
                                  correct_inputs = TRUE) {
  n <- nrow(x)
  row_issues <- vector("list", n)
  for (i in seq_len(n)) row_issues[[i]] <- character(0)

  # Check non-negative columns
  for (col in .positive_columns()) {
    if (!col %in% colnames(x)) next
    bad <- which(x[[col]] < 0 & !is.na(x[[col]]))
    if (length(bad) > 0) {
      if (correct_inputs) {
        x[[col]][bad] <- NA_real_
        for (i in bad) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Negative input: '%s' set to missing", col))
        }
      } else {
        for (i in bad) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Negative input: '%s' = %g", col, x[[col]][i]))
        }
      }
    }
  }

  # V10: Zero SD/SE (impossible for real data, indicates data entry error)
  .sd_columns <- c("mean_sd_exp", "mean_sd_nexp", "mean_sd_pooled",
                    "mean_pre_sd_exp", "mean_pre_sd_nexp",
                    "mean_change_sd_exp", "mean_change_sd_nexp",
                    "ancova_mean_sd_exp", "ancova_mean_sd_nexp", "ancova_mean_sd_pooled",
                    "md_sd", "ancova_md_sd")
  .se_columns <- c("mean_se_exp", "mean_se_nexp",
                    "mean_change_se_exp", "mean_change_se_nexp",
                    "ancova_mean_se_exp", "ancova_mean_se_nexp",
                    "ancova_md_se")
  for (col in c(.sd_columns, .se_columns)) {
    if (!col %in% colnames(x)) next
    zeros <- which(x[[col]] == 0 & !is.na(x[[col]]))
    if (length(zeros) > 0) {
      type_label <- if (col %in% .sd_columns) "SD" else "SE"
      if (correct_inputs) {
        x[[col]][zeros] <- NA_real_
        for (i in zeros) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Zero %s: '%s' = 0 (the effect size or its precision cannot be validly computed)",
                    type_label, col))
        }
      } else {
        for (i in zeros) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Zero %s: '%s' = 0 (the effect size or its precision cannot be validly computed)",
                    type_label, col))
        }
      }
    }
  }

  # V11: Bounded columns - values outside theoretical range set to missing
  # Correlations must be in [-1, 1]; reliability/proportion must be in [0, 1]
  .bounded_columns <- list(
    # alpha can be negative, only > 1 impossible (V25 warns on negatives)
    list(col = "cronbach_alpha",  lo = -Inf, up = 1, label = "alpha"),
    # icc bounded below by -1/(k-1) >= -1
    list(col = "icc",             lo = -1, up = 1, label = "ICC"),
    list(col = "prop",            lo = 0, up = 1, label = "proportion"),
    list(col = "pearson_r",       lo = -1, up = 1, label = "correlation"),
    list(col = "spearman_r",      lo = -1, up = 1, label = "correlation"),
    # phi is a 2x2 correlation coefficient in [-1, 1]; an out-of-range value
    # otherwise reaches metafor::conv.2x2(), which raises a HARD ERROR that would
    # abort the whole convert_df() run (violating the one-bad-cell-cannot-abort rule)
    list(col = "phi",             lo = -1, up = 1, label = "phi (correlation)"),
    # p-values must be in [0, 1]
    list(col = "student_t_pval",        lo = 0, up = 1, label = "p-value"),
    list(col = "anova_f_pval",          lo = 0, up = 1, label = "p-value"),
    list(col = "ancova_t_pval",         lo = 0, up = 1, label = "p-value"),
    list(col = "ancova_f_pval",         lo = 0, up = 1, label = "p-value"),
    list(col = "chisq_pval",            lo = 0, up = 1, label = "p-value"),
    list(col = "paired_t_pval_exp",     lo = 0, up = 1, label = "p-value"),
    list(col = "paired_t_pval_nexp",    lo = 0, up = 1, label = "p-value"),
    list(col = "paired_f_pval_exp",     lo = 0, up = 1, label = "p-value"),
    list(col = "paired_f_pval_nexp",    lo = 0, up = 1, label = "p-value"),
    list(col = "mean_change_pval_exp",  lo = 0, up = 1, label = "p-value"),
    list(col = "mean_change_pval_nexp", lo = 0, up = 1, label = "p-value"),
    list(col = "rr_pval",               lo = 0, up = 1, label = "p-value"),
    list(col = "rd_pval",               lo = 0, up = 1, label = "p-value"),
    list(col = "linreg_b_pval",         lo = 0, up = 1, label = "p-value"),
    list(col = "or_pval",               lo = 0, up = 1, label = "p-value"),
    list(col = "md_pval",               lo = 0, up = 1, label = "p-value"),
    list(col = "ancova_md_pval",        lo = 0, up = 1, label = "p-value"),
    list(col = "pt_bis_r_pval",         lo = 0, up = 1, label = "p-value"),
    # Baseline risk must be a proportion (0-1), not a percentage
    list(col = "baseline_risk",          lo = 0, up = 1, label = "baseline risk (proportion)"),
    # Covariate-outcome (multiple) correlation of an ANCOVA. |R| = 1 is handled
    # separately just below, because it is in range but still degenerate.
    list(col = "cov_outcome_r",   lo = -1, up = 1, label = "covariate-outcome correlation")
  )
  for (bc in .bounded_columns) {
    if (!bc$col %in% colnames(x)) next
    vals <- x[[bc$col]]
    bad <- which(!is.na(vals) & (vals < bc$lo | vals > bc$up))
    if (length(bad) > 0) {
      # Describe the valid range; one-sided when a bound is infinite
      range_txt <- if (is.infinite(bc$lo)) {
        sprintf("<= %g", bc$up)
      } else if (is.infinite(bc$up)) {
        sprintf(">= %g", bc$lo)
      } else {
        sprintf("[%g, %g]", bc$lo, bc$up)
      }
      if (correct_inputs) {
        x[[bc$col]][bad] <- NA_real_
        for (i in bad) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Out-of-range %s: '%s' = %g (valid: %s), set to missing",
                    bc$label, bc$col, vals[i], range_txt))
        }
      } else {
        for (i in bad) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Out-of-range %s: '%s' = %g (valid: %s)",
                    bc$label, bc$col, vals[i], range_txt))
        }
      }
    }
  }

  # V34: |cov_outcome_r| == 1 is inside the [-1, 1] bound checked above but is still
  # a non-identified input for every ANCOVA route, exactly as eta-squared = 1 is for
  # V32. Cooper eq. 12.24 divides the residual SD by sqrt(1 - R^2) and eq. 12.26
  # multiplies the leading variance term by (1 - R^2), so at |R| = 1 the marginal SD
  # is infinite and the sampling variance collapses: 13 of the 14 ANCOVA routes return
  # ES = 0 with SE = 0 (infinite meta-analytic weight), and the two that build the
  # SMD from an already-marginal SD (ancova_means_sd_pooled_crude, cohen_d_adj) return
  # a perfectly plausible ES with a silently deflated SE and no other flag. Handled
  # here rather than by widening the bound so the message names the real problem.
  if ("cov_outcome_r" %in% colnames(x)) {
    vals <- x$cov_outcome_r
    bad <- which(!is.na(vals) & abs(vals) == 1)
    if (length(bad) > 0) {
      if (correct_inputs) x$cov_outcome_r[bad] <- NA_real_
      for (i in bad) {
        row_issues[[i]] <- c(row_issues[[i]], sprintf(
          "[INVALID] Degenerate covariate-outcome correlation: 'cov_outcome_r' = %g implies a residual SD of zero, so the marginal-scale back-transformation 1/sqrt(1 - R^2) diverges and the sampling variance collapses to zero%s",
          vals[i], if (correct_inputs) ", set to missing" else ""))
      }
    }
  }

  # V32: eta-squared must lie in [0, 1). A value >= 1 makes the Cohen's d
  # conversion diverge -- d = 2*sqrt(eta2 / (1 - eta2)) for a raw eta-squared,
  # and the implied ANCOVA F = eta2 * df / (1 - eta2) for an adjusted one -- so
  # such a row is set to missing (correct_inputs = TRUE) or flagged (FALSE).
  # Negatives are already caught by the non-negative-column check above. A
  # large-but-valid eta-squared that implies a huge SMD is surfaced downstream
  # by the post-computation Large-SMD check, so no separate plausibility tier is
  # added here.
  for (col in c("etasq", "etasq_adj")) {
    if (!col %in% colnames(x)) next
    vals <- x[[col]]
    bad <- which(!is.na(vals) & vals >= 1)
    if (length(bad) == 0) next
    if (correct_inputs) {
      x[[col]][bad] <- NA_real_
      for (i in bad) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INVALID] Out-of-range eta-squared: '%s' = %g (valid: [0, 1)), set to missing",
                  col, vals[i]))
      }
    } else {
      for (i in bad) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INVALID] Out-of-range eta-squared: '%s' = %g (valid: [0, 1))",
                  col, vals[i]))
      }
    }
  }

  # V25: negative alpha/icc kept, warn only
  for (rel in list(c("cronbach_alpha", "alpha"), c("icc", "ICC"))) {
    if (!rel[1] %in% colnames(x)) next
    neg <- which(!is.na(x[[rel[1]]]) & x[[rel[1]]] < 0)
    for (i in neg) {
      row_issues[[i]] <- c(row_issues[[i]],
        sprintf("[UNUSUAL] Negative %s: '%s' = %g. Mathematically possible but indicates serious measurement problems or an extraction error - verify against the primary report",
                rel[2], rel[1], x[[rel[1]]][i]))
    }
  }

  # V31: agreement-type ICC rows use a one-way variance approximation that
  # assumes negligible between-rater variance -- the SE is anti-conservative
  # when raters differ systematically (informational, always active, like
  # V18/V21). icc_type = NA or an absent column resolves to the package
  # default "agreement".
  if ("icc" %in% colnames(x)) {
    icc_type_eff <- if ("icc_type" %in% colnames(x)) {
      tt <- as.character(x[["icc_type"]])
      tt[is.na(tt)] <- "agreement"
      tt
    } else {
      rep("agreement", nrow(x))
    }
    agr <- which(!is.na(x[["icc"]]) & icc_type_eff == "agreement")
    for (i in agr) {
      row_issues[[i]] <- c(row_issues[[i]],
        "[INFO] ICC agreement-type SE assumes negligible between-rater variance and can be anti-conservative when raters differ systematically (see the es_from_icc documentation)")
    }
  }

  # Check CI consistency
  for (tri in .ci_triplets()) {
    if (!all(c(tri$val, tri$lo, tri$up) %in% colnames(x))) next
    val <- x[[tri$val]]; lo <- x[[tri$lo]]; up <- x[[tri$up]]

    # V2: Inverted CI (lo > up) -- rounding-aware: only flag when the inversion
    # exceeds the combined rounding error of the two bounds (mirrors V4/V28), so a
    # borderline pair like lo = 0.15, up = 0.14 reported at different precisions is
    # preserved rather than destroyed.
    both_lu <- which(!is.na(lo) & !is.na(up) & lo > up)
    inverted <- both_lu[vapply(both_lu, function(i) {
      pad <- 0.5 * 10^(-.count_decimals(lo[i])) + 0.5 * 10^(-.count_decimals(up[i]))
      lo[i] > up[i] + pad
    }, logical(1))]
    if (length(inverted) > 0) {
      if (correct_inputs) {
        x[[tri$val]][inverted] <- NA_real_
        x[[tri$lo]][inverted] <- NA_real_
        x[[tri$up]][inverted] <- NA_real_
        for (i in inverted) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Inverted CI for '%s' (lo > up) set to missing", tri$val))
        }
      } else {
        for (i in inverted) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Inverted CI for '%s': lo = %g > up = %g",
                    tri$val, lo[i], up[i]))
        }
      }
    }

    # V3: Value outside CI (re-read after V2 cleanup) -- rounding-aware. The point
    # value and the CI bounds may be reported to different decimal precisions, so a
    # value that merely rounds to a bound (e.g. md = 0.1 [1 dp] vs lo = 0.15 [2 dp],
    # both consistent with a true md in [0.145, 0.155)) must NOT be destroyed. Flag
    # (and, under correct_inputs, NA) only when the value lies outside by more than
    # the combined rounding error of the value and the offending bound.
    val <- x[[tri$val]]; lo <- x[[tri$lo]]; up <- x[[tri$up]]
    cand <- which(!is.na(val) & !is.na(lo) & !is.na(up) & lo <= up &
                    (val < lo | val > up))
    outside <- cand[vapply(cand, function(i) {
      re_val <- 0.5 * 10^(-.count_decimals(val[i]))
      pad_lo <- re_val + 0.5 * 10^(-.count_decimals(lo[i]))
      pad_up <- re_val + 0.5 * 10^(-.count_decimals(up[i]))
      (val[i] < lo[i] - pad_lo) || (val[i] > up[i] + pad_up)
    }, logical(1))]
    if (length(outside) > 0) {
      if (correct_inputs) {
        x[[tri$val]][outside] <- NA_real_
        x[[tri$lo]][outside] <- NA_real_
        x[[tri$up]][outside] <- NA_real_
        for (i in outside) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Value outside CI for '%s' set to missing", tri$val))
        }
      } else {
        for (i in outside) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Value outside CI for '%s': %g not in [%g, %g]",
                    tri$val, val[i], lo[i], up[i]))
        }
      }
    }
  }

  # V4: ci asymmetry, warn only (rounding); or/rr checked on log scale
  for (tri in .ci_triplets()) {
    if (!all(c(tri$val, tri$lo, tri$up) %in% colnames(x))) next
    val <- x[[tri$val]]; lo <- x[[tri$lo]]; up <- x[[tri$up]]

    has_all <- which(!is.na(val) & !is.na(lo) & !is.na(up) & lo <= up)
    if (length(has_all) == 0) next

    # For exp-scale CIs, convert to log scale before checking asymmetry
    if (identical(tri$scale, "exp")) {
      # Need all values > 0 for log transform
      has_all <- has_all[val[has_all] > 0 & lo[has_all] > 0 & up[has_all] > 0]
      if (length(has_all) == 0) next
      check_val <- log(val[has_all])
      check_lo <- log(lo[has_all])
      check_up <- log(up[has_all])
      log_transformed <- rep(TRUE, length(has_all))
    } else if (identical(tri$scale, "user")) {
      # User-entered data: detect scale from user_es_original_measure column.
      # Ratio types (OR, RR, HR, IRR) need log-transform; others are additive.
      measure_col <- if (grepl("crude$", tri$val)) {
        "user_es_original_measure_crude"
      } else {
        "user_es_original_measure_adj"
      }
      if (!measure_col %in% colnames(x)) next
      ratio_types <- c("or", "rr", "hr", "irr")
      is_ratio <- tolower(as.character(x[[measure_col]])) %in% ratio_types

      ratio_idx <- has_all[is_ratio[has_all]]
      additive_idx <- has_all[!is_ratio[has_all]]

      if (length(ratio_idx) > 0) {
        ratio_idx <- ratio_idx[val[ratio_idx] > 0 & lo[ratio_idx] > 0 & up[ratio_idx] > 0]
      }
      # Combine: ratio on log scale, additive on original scale
      check_val <- numeric(length(has_all))
      check_lo <- numeric(length(has_all))
      check_up <- numeric(length(has_all))
      # Map back to positions in has_all
      for (j in seq_along(has_all)) {
        idx <- has_all[j]
        if (idx %in% ratio_idx) {
          check_val[j] <- log(val[idx])
          check_lo[j] <- log(lo[idx])
          check_up[j] <- log(up[idx])
        } else {
          check_val[j] <- val[idx]
          check_lo[j] <- lo[idx]
          check_up[j] <- up[idx]
        }
      }
      # Remove rows that were ratio but had non-positive values (couldn't log-transform)
      keep <- rep(TRUE, length(has_all))
      for (j in seq_along(has_all)) {
        idx <- has_all[j]
        if (is_ratio[idx] && !(idx %in% ratio_idx)) keep[j] <- FALSE
      }
      has_all <- has_all[keep]
      check_val <- check_val[keep]
      check_lo <- check_lo[keep]
      check_up <- check_up[keep]
      if (length(has_all) == 0) next
      # Track which rows were log-transformed (for rounding computation)
      log_transformed <- has_all %in% ratio_idx
    } else {
      check_val <- val[has_all]
      check_lo <- lo[has_all]
      check_up <- up[has_all]
      log_transformed <- rep(FALSE, length(has_all))
    }

    lower_dist <- check_val - check_lo
    upper_dist <- check_up - check_val
    avg_dist <- (lower_dist + upper_dist) / 2

    # Only check where avg_dist > 0 (non-zero CI width)
    pos <- avg_dist > 0
    if (!any(pos)) next

    checkable <- has_all[pos]
    is_log <- log_transformed[pos]
    asym_pct <- (abs(upper_dist[pos] - lower_dist[pos]) / avg_dist[pos]) * 100

    # per-row threshold: max asymmetry explainable by rounding at the reported precision
    # max rounding asymmetry = sum of the three rounding errors / avg half-width
    rounding_threshold <- rep(max_asymmetry, sum(pos))
    orig_val <- val[checkable]
    orig_lo  <- lo[checkable]
    orig_up  <- up[checkable]
    for (j in seq_along(checkable)) {
      dp_val <- .count_decimals(orig_val[j])
      dp_lo  <- .count_decimals(orig_lo[j])
      dp_up  <- .count_decimals(orig_up[j])
      re_val <- 0.5 * 10^(-dp_val)
      re_lo  <- 0.5 * 10^(-dp_lo)
      re_up  <- 0.5 * 10^(-dp_up)
      if (identical(tri$scale, "exp") || is_log[j]) {
        # On log scale, rounding error ~ re / value (first-order)
        if (orig_val[j] > 0 && orig_lo[j] > 0 && orig_up[j] > 0) {
          re_val <- re_val / orig_val[j]
          re_lo  <- re_lo  / orig_lo[j]
          re_up  <- re_up  / orig_up[j]
        }
      }
      max_rounding_asym <- (re_val + re_lo + re_up) / avg_dist[pos][j] * 100
      rounding_threshold[j] <- max(max_asymmetry, max_rounding_asym * 2)
    }

    asymmetric <- checkable[asym_pct > rounding_threshold]
    if (length(asymmetric) > 0) {
      for (i in asymmetric) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INFO] Asymmetric CI for '%s' (>%g%% asymmetry)",
                  tri$val, max_asymmetry))
      }
    }
  }

  # V5: Logical consistency - cases cannot exceed total sample size
  # Each pair is (part, total): n_cases_exp <= n_exp, etc.
  part_total_pairs <- list(
    c("n_cases_exp",  "n_exp"),
    c("n_cases_nexp", "n_nexp"),
    c("n_controls_exp",  "n_exp"),
    c("n_controls_nexp", "n_nexp"),
    c("n_cases",   "n_sample"),
    c("n_controls", "n_sample")
  )
  for (pair in part_total_pairs) {
    part_col <- pair[1]; total_col <- pair[2]
    if (!part_col %in% colnames(x) || !total_col %in% colnames(x)) next
    bad <- which(!is.na(x[[part_col]]) & !is.na(x[[total_col]]) &
                 x[[part_col]] > x[[total_col]])
    if (length(bad) > 0) {
      if (correct_inputs) {
        x[[part_col]][bad] <- NA_real_
        x[[total_col]][bad] <- NA_real_
        for (i in bad) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Inconsistent input: '%s' > '%s', both set to missing",
                    part_col, total_col))
        }
      } else {
        for (i in bad) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Inconsistent input: '%s' = %g > '%s' = %g",
                    part_col, x[[part_col]][i], total_col, x[[total_col]][i]))
        }
      }
    }
  }

  # V15: warn only, mismatch is often legitimate (multi-arm trial, analysis subset)
  if (all(c("n_exp", "n_nexp", "n_sample") %in% colnames(x))) {
    ne <- x[["n_exp"]]; nn <- x[["n_nexp"]]; ns <- x[["n_sample"]]
    mismatch <- which(!is.na(ne) & !is.na(nn) & !is.na(ns) &
                      (ne + nn) != ns)
    if (length(mismatch) > 0) {
      for (i in mismatch) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[UNUSUAL] Sample size mismatch: n_exp + n_nexp = %g != n_sample = %g. Often legitimate (multi-arm trial, or n_sample is an analysis subset) - verify these refer to the same pool",
                  ne[i] + nn[i], ns[i]))
      }
    }
  }

  # V7: a + b must equal the margin; all three set to NA (cannot tell which is wrong)
  additive_triples <- list(
    list(a = "n_cases_exp",  b = "n_controls_exp",  total = "n_exp"),
    list(a = "n_cases_nexp", b = "n_controls_nexp", total = "n_nexp"),
    list(a = "n_cases",      b = "n_controls",      total = "n_sample"),
    # row margins
    list(a = "n_cases_exp",  b = "n_cases_nexp",    total = "n_cases"),
    list(a = "n_controls_exp", b = "n_controls_nexp", total = "n_controls")
  )
  for (triple in additive_triples) {
    if (!all(c(triple$a, triple$b, triple$total) %in% colnames(x))) next
    a <- x[[triple$a]]; b <- x[[triple$b]]; total <- x[[triple$total]]
    mismatch <- which(!is.na(a) & !is.na(b) & !is.na(total) &
                      (a + b) != total)
    if (length(mismatch) > 0) {
      if (correct_inputs) {
        x[[triple$a]][mismatch] <- NA_real_
        x[[triple$b]][mismatch] <- NA_real_
        x[[triple$total]][mismatch] <- NA_real_
        for (i in mismatch) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] 2x2 inconsistency: %s + %s != %s",
                    triple$a, triple$b, triple$total))
        }
      } else {
        for (i in mismatch) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] 2x2 inconsistency: %s + %s != %s (sum=%g, total=%g)",
                    triple$a, triple$b, triple$total,
                    a[i] + b[i], total[i]))
        }
      }
    }
  }

  # V33: reported proportion vs its own counts (rounding-aware, warn only, data
  # preserved -- two independently transcribed values disagree and neither can be
  # inferred to be the wrong one; both prop routes still compute and the generic
  # cross-method discordance thresholds are far too coarse on the [0,1] scale to
  # catch this). V33a: prop != n_cases/n_sample beyond the rounding error of the
  # reported prop. V33b (GRIM-style, only when n_cases is absent): prop is not
  # achievable as k/n_sample for ANY integer k, beyond the reported rounding --
  # the count analogue of the mean-granularity (GRIM) test. Messages carry the
  # quoted 'prop' so .v_flag_matches_scope routes them to the crude scope only.
  if (all(c("prop", "n_sample") %in% colnames(x))) {
    p_v  <- suppressWarnings(as.numeric(x[["prop"]]))
    ns_v <- suppressWarnings(as.numeric(x[["n_sample"]]))
    nc_v <- if ("n_cases" %in% colnames(x)) {
      suppressWarnings(as.numeric(x[["n_cases"]]))
    } else {
      rep(NA_real_, n)
    }
    for (i in seq_len(n)) {
      if (is.na(p_v[i]) || is.na(ns_v[i]) || ns_v[i] <= 0) next
      if (p_v[i] < 0 || p_v[i] > 1) next  # V11's job
      tol <- 0.5 * 10^(-.count_decimals(p_v[i])) + 1e-12
      if (!is.na(nc_v[i])) {
        # V33a: direct consistency against the reported counts
        if (nc_v[i] > ns_v[i]) next  # V5's job
        computed <- nc_v[i] / ns_v[i]
        if (abs(p_v[i] - computed) > tol) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INVALID] Proportion inconsistent with counts: 'prop' = %g but n_cases/n_sample = %g/%g = %.4g. One of the three values is wrong - verify against the primary report",
                    p_v[i], nc_v[i], ns_v[i], computed))
        }
      } else {
        # V33b: achievable-fraction (GRIM-style) check, integer N only
        if (ns_v[i] < 2 || ns_v[i] != round(ns_v[i])) next
        k <- max(0, min(ns_v[i], round(p_v[i] * ns_v[i])))
        nearest <- k / ns_v[i]
        if (abs(p_v[i] - nearest) > tol) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[UNUSUAL] Proportion not achievable for any integer count: 'prop' = %g with n_sample = %d (nearest achievable: %d/%d = %.4g). The proportion, the sample size, or both may be wrong - verify against the primary report",
                    p_v[i], as.integer(ns_v[i]), as.integer(k),
                    as.integer(ns_v[i]), nearest))
        }
      }
    }
  }

  # V8: Zero cell in 2x2 table (warn only - continuity correction will be applied)
  # Controlled by enable_informational (default FALSE: suppressed)
  cell_cols <- c("n_cases_exp", "n_cases_nexp", "n_controls_exp", "n_controls_nexp")
  present_cells <- cell_cols[cell_cols %in% colnames(x)]
  if (length(present_cells) == 4) {
    # V8b: Double-zero studies - both arms have zero events (degenerate OR/RR)
    cases_exp <- x[["n_cases_exp"]]; cases_nexp <- x[["n_cases_nexp"]]
    double_zero <- which(!is.na(cases_exp) & !is.na(cases_nexp) &
                         cases_exp == 0 & cases_nexp == 0)
    if (length(double_zero) > 0) {
      for (i in double_zero) {
        row_issues[[i]] <- c(row_issues[[i]],
          "[INFO] Double-zero study: both arms have zero events (n_cases_exp = 0, n_cases_nexp = 0). OR/RR are degenerate even with continuity correction")
      }
    }
    # Also check double-zero on the controls side
    if (all(c("n_controls_exp", "n_controls_nexp") %in% colnames(x))) {
      ctrl_exp <- x[["n_controls_exp"]]; ctrl_nexp <- x[["n_controls_nexp"]]
      double_zero_ctrl <- which(!is.na(ctrl_exp) & !is.na(ctrl_nexp) &
                                ctrl_exp == 0 & ctrl_nexp == 0)
      if (length(double_zero_ctrl) > 0) {
        for (i in double_zero_ctrl) {
          row_issues[[i]] <- c(row_issues[[i]],
            "[INFO] Double-zero study: both arms have zero non-events (n_controls_exp = 0, n_controls_nexp = 0). OR/RR are degenerate even with continuity correction")
        }
      }
    }

    # V8: Individual zero cells (informational)
    if (isTRUE(enable_informational)) {
      for (col in present_cells) {
        zero_rows <- which(!is.na(x[[col]]) & x[[col]] == 0)
        if (length(zero_rows) > 0) {
          for (i in zero_rows) {
            row_issues[[i]] <- c(row_issues[[i]],
              sprintf("[INFO] Zero cell in 2x2 table ('%s' = 0): OR/RR use +0.5 continuity correction",
                      col))
          }
        }
      }
    }
  }

  # V16: n = 1 per arm (SMD undefined - pooled SD has df = 0)
  for (ncol in c("n_exp", "n_nexp")) {
    if (!ncol %in% colnames(x)) next
    one_rows <- which(!is.na(x[[ncol]]) & x[[ncol]] == 1)
    if (length(one_rows) > 0) {
      for (i in one_rows) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INVALID] %s = 1: within-group variance is undefined (df = 0), SMD cannot be estimated",
                  ncol))
      }
    }
  }

  # V24: min <= q1 <= median <= q3 <= max, warn only
  # reverse_med flips sign at computation, raw cells stay ascending
  .quantile_groups <- list(
    list(suffix = "exp",
         cols  = c("min_exp",  "q1_exp",  "med_exp",  "q3_exp",  "max_exp")),
    list(suffix = "nexp",
         cols  = c("min_nexp", "q1_nexp", "med_nexp", "q3_nexp", "max_nexp"))
  )
  .quantile_labels <- c(min_exp = "min", q1_exp = "Q1", med_exp = "median",
                        q3_exp = "Q3", max_exp = "max",
                        min_nexp = "min", q1_nexp = "Q1", med_nexp = "median",
                        q3_nexp = "Q3", max_nexp = "max")
  for (qg in .quantile_groups) {
    present <- qg$cols[qg$cols %in% colnames(x)]
    if (length(present) < 2) next
    for (i in seq_len(n)) {
      vals <- vapply(present, function(cc) {
        v <- x[[cc]][i]; if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v))
      }, numeric(1))
      ok <- !is.na(vals)
      if (sum(ok) < 2) next
      pv <- vals[ok]; pn <- present[ok]   # canonical min->max order preserved
      bad <- which(diff(pv) < 0)
      if (length(bad) > 0) {
        bp <- bad[1]
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INVALID] Non-monotonic quantile summary for '%s' (need min<=Q1<=median<=Q3<=max): %s (%g) < %s (%g)",
                  qg$suffix,
                  .quantile_labels[[pn[bp + 1]]], pv[bp + 1],
                  .quantile_labels[[pn[bp]]],     pv[bp]))
      }
    }
  }

  # arms used by V28 and V29 (endpoint only, no min_pre/max_pre columns)
  .range_arms <- list(
    list(sfx = "exp",  mean = "mean_exp",  min = "min_exp",  max = "max_exp",
         sd = "mean_sd_exp",  n = "n_exp"),
    list(sfx = "nexp", mean = "mean_nexp", min = "min_nexp", max = "max_nexp",
         sd = "mean_sd_nexp", n = "n_nexp")
  )

  # V28: mean must lie in [min, max], warn only; tolerance for rounding of both values
  for (arm in .range_arms) {
    if (!all(c(arm$mean, arm$min, arm$max) %in% colnames(x))) next
    mv <- suppressWarnings(as.numeric(x[[arm$mean]]))
    lo <- suppressWarnings(as.numeric(x[[arm$min]]))
    hi <- suppressWarnings(as.numeric(x[[arm$max]]))
    chk <- which(!is.na(mv) & !is.na(lo) & !is.na(hi) & lo <= hi)
    for (i in chk) {
      re_mean <- 0.5 * 10^(-.count_decimals(mv[i]))
      pad_lo <- re_mean + 0.5 * 10^(-.count_decimals(lo[i]))
      pad_hi <- re_mean + 0.5 * 10^(-.count_decimals(hi[i]))
      if (mv[i] < lo[i] - pad_lo || mv[i] > hi[i] + pad_hi) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INVALID] Mean outside its reported range for '%s': mean = %g not in [min = %g, max = %g]",
                  arm$sfx, mv[i], lo[i], hi[i]))
      }
    }
  }

  # V29: (max - min)/SD should be near xi(n) = 2*qnorm((n - .375)/(n + .25))
  # (Wan 2014), warn only
  range_sd_lo_mult <- 0.5
  range_sd_hi_mult <- 2.5
  for (arm in .range_arms) {
    if (!all(c(arm$min, arm$max, arm$sd, arm$n) %in% colnames(x))) next
    lo  <- suppressWarnings(as.numeric(x[[arm$min]]))
    hi  <- suppressWarnings(as.numeric(x[[arm$max]]))
    sdv <- suppressWarnings(as.numeric(x[[arm$sd]]))
    nv  <- suppressWarnings(as.numeric(x[[arm$n]]))
    chk <- which(!is.na(lo) & !is.na(hi) & !is.na(sdv) & !is.na(nv) &
                 hi > lo & sdv > 0 & nv >= 2)
    for (i in chk) {
      xi <- 2 * stats::qnorm((nv[i] - 0.375) / (nv[i] + 0.25))
      if (!is.finite(xi) || xi <= 0) next
      r_obs <- (hi[i] - lo[i]) / sdv[i]
      if (r_obs < range_sd_lo_mult * xi || r_obs > range_sd_hi_mult * xi) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[UNUSUAL] Range/SD ratio implausible for '%s': (max - min)/SD = %.2f (expected ~%.1f for n = %d). Possible SD-as-SE, variance-as-SD, or range/IQR mix-up",
                  arm$sfx, r_obs, xi, as.integer(nv[i])))
      }
    }
  }

  # V26: within-subject measure but reported CI matches the independent-groups SE
  # implied critical value compared to z and t; no assumed r needed. warn only
  within_meas <- !is.null(measure) && tolower(measure) %in% c("mdw", "dw", "gw")
  if (within_meas) {
    z_crit <- stats::qnorm(0.975)
    for (sfx in c("crude", "adj")) {
      lo_col <- paste0("user_ci_lo_", sfx)
      up_col <- paste0("user_ci_up_", sfx)
      need26 <- c(lo_col, up_col, "mean_sd_exp", "mean_sd_nexp", "n_exp", "n_nexp")
      if (!all(need26 %in% colnames(x))) next
      for (i in seq_len(n)) {
        lo  <- suppressWarnings(as.numeric(x[[lo_col]][i]))
        up  <- suppressWarnings(as.numeric(x[[up_col]][i]))
        sde <- suppressWarnings(as.numeric(x[["mean_sd_exp"]][i]))
        sdn <- suppressWarnings(as.numeric(x[["mean_sd_nexp"]][i]))
        ne  <- suppressWarnings(as.numeric(x[["n_exp"]][i]))
        nn  <- suppressWarnings(as.numeric(x[["n_nexp"]][i]))
        if (any(is.na(c(lo, up, sde, sdn, ne, nn)))) next
        # within-subject signature: equal arm n; need usable SDs and a real CI
        if (ne != nn || ne < 2 || sde <= 0 || sdn <= 0 || up <= lo) next
        indep_se <- sqrt(sde^2 / ne + sdn^2 / nn)
        if (!is.finite(indep_se) || indep_se <= 0) next
        implied_crit <- ((up - lo) / 2) / indep_se
        t_crit <- stats::qt(0.975, ne + nn - 2)
        matches <- (abs(implied_crit - z_crit) / z_crit <= paired_as_indep_tol) ||
                   (abs(implied_crit - t_crit) / t_crit <= paired_as_indep_tol)
        if (isTRUE(matches)) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[UNUSUAL] Within-subject design (%s) but the reported CI matches an independent-groups SE (n_exp = n_nexp = %g, implied critical value %.2f). Paired pre/post data analysed as two independent groups overestimates the variance - verify the SE formula (a paired/within-subject SE should be used)",
                    tolower(measure), ne, implied_crit))
        }
      }
    }
  }

  # V9: possible change score entered as endpoint (small SD/|mean|)
  if (isTRUE(enable_informational)) {
    m_exp_col <- "mean_exp"; sd_exp_col <- "mean_sd_exp"
    m_nexp_col <- "mean_nexp"; sd_nexp_col <- "mean_sd_nexp"
    has_exp <- m_exp_col %in% colnames(x) && sd_exp_col %in% colnames(x)
    has_nexp <- m_nexp_col %in% colnames(x) && sd_nexp_col %in% colnames(x)

    for (i in seq_len(n)) {
      m_e <- if (has_exp) x[[m_exp_col]][i] else NA_real_
      s_e <- if (has_exp) x[[sd_exp_col]][i] else NA_real_
      m_n <- if (has_nexp) x[[m_nexp_col]][i] else NA_real_
      s_n <- if (has_nexp) x[[sd_nexp_col]][i] else NA_real_

      cv_e <- if (!is.na(m_e) && !is.na(s_e) && s_e > 0 && abs(m_e) > 0) s_e / abs(m_e) else NA_real_
      cv_n <- if (!is.na(m_n) && !is.na(s_n) && s_n > 0 && abs(m_n) > 0) s_n / abs(m_n) else NA_real_

      # negative mean, small cv
      if (!is.na(cv_e) && !is.na(m_e) && m_e < 0 && cv_e < 0.3) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INFO] Possible change score as endpoint: '%s' = %g with small SD = %g (SD/|mean| = %s)",
                  m_exp_col, m_e, s_e, round(cv_e, 2)))
      }
      if (!is.na(cv_n) && !is.na(m_n) && m_n < 0 && cv_n < 0.3) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INFO] Possible change score as endpoint: '%s' = %g with small SD = %g (SD/|mean| = %s)",
                  m_nexp_col, m_n, s_n, round(cv_n, 2)))
      }

      # positive means: both cvs must be very small
      if (!is.na(cv_e) && !is.na(cv_n) &&
          !is.na(m_e) && m_e > 0 && !is.na(m_n) && m_n > 0 &&
          cv_e < 0.1 && cv_n < 0.1) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INFO] Possible change scores as endpoint: both group means positive with very small SD/|mean| (exp: %s, nexp: %s)",
                  round(cv_e, 2), round(cv_n, 2)))
      }
    }
  }

  # V12: extreme sd ratio between arms (Kanukula 2024 #5)
  sd_pairs <- list(
    c("mean_sd_exp", "mean_sd_nexp"),
    c("mean_pre_sd_exp", "mean_pre_sd_nexp"),
    c("mean_change_sd_exp", "mean_change_sd_nexp"),
    c("ancova_mean_sd_exp", "ancova_mean_sd_nexp")
  )
  for (pair in sd_pairs) {
    if (!all(pair %in% colnames(x))) next
    sd1 <- x[[pair[1]]]; sd2 <- x[[pair[2]]]
    checkable <- which(!is.na(sd1) & !is.na(sd2) & sd1 > 0 & sd2 > 0)
    if (length(checkable) == 0) next
    ratio <- pmax(sd1[checkable], sd2[checkable]) / pmin(sd1[checkable], sd2[checkable])
    flagged <- checkable[ratio > sd_ratio_max]
    if (length(flagged) > 0) {
      for (i in flagged) {
        r <- max(sd1[i], sd2[i]) / min(sd1[i], sd2[i])
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[UNUSUAL] Extreme SD ratio between arms: %s = %g, %s = %g (ratio: %.1f, threshold: %g). Possible variance/SD confusion or extraction error",
                  pair[1], sd1[i], pair[2], sd2[i], r, sd_ratio_max))
      }
    }
  }

  # V14: identical sds across arms (Kanukula 2024 #22)
  if (isTRUE(enable_informational)) {
    for (pair in sd_pairs) {
      if (!all(pair %in% colnames(x))) next
      sd1 <- x[[pair[1]]]; sd2 <- x[[pair[2]]]
      identical_rows <- which(!is.na(sd1) & !is.na(sd2) & sd1 > 0 & sd2 > 0 &
                              sd1 == sd2)
      if (length(identical_rows) > 0) {
        for (i in identical_rows) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INFO] Identical SDs across arms: %s = %s = %g. Verify this is not the pooled or between-group SD",
                    pair[1], pair[2], sd1[i]))
        }
      }
    }
  }

  # V17: identical means across arms
  if (isTRUE(enable_informational)) {
    mean_pairs <- list(
      c("mean_exp", "mean_nexp"),
      c("mean_pre_exp", "mean_pre_nexp"),
      c("mean_change_exp", "mean_change_nexp"),
      c("ancova_mean_exp", "ancova_mean_nexp")
    )
    for (pair in mean_pairs) {
      if (!all(pair %in% colnames(x))) next
      m1 <- x[[pair[1]]]; m2 <- x[[pair[2]]]
      identical_rows <- which(!is.na(m1) & !is.na(m2) & m1 == m2)
      if (length(identical_rows) > 0) {
        for (i in identical_rows) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INFO] Identical means across arms: %s = %s = %g. Verify this is not a copy-paste error",
                    pair[1], pair[2], m1[i]))
        }
      }
    }
  }

  # V18: sd_baseline / sd_endpoint below threshold (df-weighted pooled sds)
  if (all(c("mean_pre_sd_exp", "mean_pre_sd_nexp",
            "mean_sd_exp",     "mean_sd_nexp",
            "n_exp", "n_nexp") %in% colnames(x))) {
    n_e <- x$n_exp; n_n <- x$n_nexp
    spre_e <- x$mean_pre_sd_exp; spre_n <- x$mean_pre_sd_nexp
    spost_e <- x$mean_sd_exp;    spost_n <- x$mean_sd_nexp
    checkable <- which(!is.na(spre_e) & !is.na(spre_n) &
                        !is.na(spost_e) & !is.na(spost_n) &
                        !is.na(n_e) & !is.na(n_n) &
                        spre_e > 0 & spre_n > 0 &
                        spost_e > 0 & spost_n > 0 &
                        n_e + n_n > 2)
    for (i in checkable) {
      sp_pre  <- sqrt(((n_e[i] - 1) * spre_e[i]^2 +
                       (n_n[i] - 1) * spre_n[i]^2) / (n_e[i] + n_n[i] - 2))
      sp_post <- sqrt(((n_e[i] - 1) * spost_e[i]^2 +
                       (n_n[i] - 1) * spost_n[i]^2) / (n_e[i] + n_n[i] - 2))
      q <- sp_pre / sp_post
      if (is.finite(q) && q < sd_ratio_bl_ep_min) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INFO] SD_baseline / SD_endpoint = %.2f (below %.2f). Pre-post SMD formulas (Bonett, d_rm, d_av) target a baseline-SD-standardised estimand structurally different from the ANCOVA-on-endpoint-SD target at this variance regime - choice of estimand non-trivial",
                  q, sd_ratio_bl_ep_min))
      }
    }
  }

  # V30: baseline statistic copied into the endpoint slot (Kanukula 2024 #9/#10).
  # Two timepoints from the same patients essentially never reproduce an
  # identical SD (or mean); an exact baseline<->endpoint tie signals a copy.
  # Complementary to V18, which tests the sd_baseline/sd_endpoint ratio and by
  # design ignores the ratio = 1.00 a copy produces. Warn only, data preserved,
  # always active (a high-specificity exact-equality check, like V24/V28).
  .v30_equal <- function(pre_col, post_col, positive_only) {
    if (!all(c(pre_col, post_col) %in% colnames(x))) return(rep(NA_real_, n))
    pre  <- suppressWarnings(as.numeric(x[[pre_col]]))
    post <- suppressWarnings(as.numeric(x[[post_col]]))
    ok <- !is.na(pre) & !is.na(post) & pre == post
    if (positive_only) ok <- ok & pre > 0
    ifelse(ok, post, NA_real_)
  }
  v30_sd_e <- .v30_equal("mean_pre_sd_exp",  "mean_sd_exp",  TRUE)
  v30_sd_n <- .v30_equal("mean_pre_sd_nexp", "mean_sd_nexp", TRUE)
  v30_mn_e <- .v30_equal("mean_pre_exp",     "mean_exp",     FALSE)
  v30_mn_n <- .v30_equal("mean_pre_nexp",    "mean_nexp",    FALSE)
  for (i in seq_len(n)) {
    sd_e <- v30_sd_e[i]; sd_n <- v30_sd_n[i]
    has_e <- !is.na(sd_e); has_n <- !is.na(sd_n)
    both_sd <- has_e && has_n
    # a lone-arm SD tie needs >= 2 decimals to be distinctive (integer SDs
    # collide by chance too often to flag on a single arm)
    lone_val <- if (has_e && !has_n) sd_e else if (has_n && !has_e) sd_n else NA_real_
    fire <- both_sd || (!is.na(lone_val) && .count_decimals(lone_val) >= 2)
    if (!isTRUE(fire)) next
    # NB: no "; " anywhere in this message - the flag merge in .flag_es_quality
    # splits on "; ", so an internal one would shatter and mis-route the flag
    parts <- character(0)
    if (has_e) parts <- c(parts, sprintf("'mean_sd_exp' = 'mean_pre_sd_exp' = %s",
                                          format(sd_e, scientific = FALSE)))
    if (has_n) parts <- c(parts, sprintf("'mean_sd_nexp' = 'mean_pre_sd_nexp' = %s",
                                          format(sd_n, scientific = FALSE)))
    mean_note <- if (!is.na(v30_mn_e[i]) && !is.na(v30_mn_n[i])) {
      " (both group means are identical across timepoints too - almost certainly a copied baseline block)"
    } else {
      ""
    }
    row_issues[[i]] <- c(row_issues[[i]], sprintf(
      "[UNUSUAL] Baseline SD copied into endpoint slot (%s). Two timepoints from the same patients virtually never reproduce an identical SD - verify the endpoint (post-treatment) SD was not overwritten with the baseline SD%s",
      paste(parts, collapse = " and "), mean_note))
  }

  # V21: standardised baseline imbalance above threshold
  if (all(c("mean_pre_exp", "mean_pre_nexp",
            "mean_pre_sd_exp", "mean_pre_sd_nexp",
            "n_exp", "n_nexp") %in% colnames(x))) {
    checkable <- which(!is.na(x$mean_pre_exp) & !is.na(x$mean_pre_nexp) &
                        !is.na(x$mean_pre_sd_exp) & !is.na(x$mean_pre_sd_nexp) &
                        !is.na(x$n_exp) & !is.na(x$n_nexp) &
                        x$mean_pre_sd_exp > 0 & x$mean_pre_sd_nexp > 0 &
                        x$n_exp + x$n_nexp > 2)
    for (i in checkable) {
      sp_pre <- sqrt(((x$n_exp[i] - 1) * x$mean_pre_sd_exp[i]^2 +
                      (x$n_nexp[i] - 1) * x$mean_pre_sd_nexp[i]^2) /
                     (x$n_exp[i] + x$n_nexp[i] - 2))
      if (is.finite(sp_pre) && sp_pre > 0) {
        imb <- abs(x$mean_pre_exp[i] - x$mean_pre_nexp[i]) / sp_pre
        if (is.finite(imb) && imb > baseline_imbalance_max) {
          row_issues[[i]] <- c(row_issues[[i]],
            sprintf("[INFO] |standardised baseline imbalance| = %.2f (above %.2f). Endpoint-only SMD formulas inherit ANCOVA omitted-variable bias of magnitude r * imb_std - ANCOVA-adjusted statistics recommended at this imbalance level",
                    imb, baseline_imbalance_max))
        }
      }
    }
  }

  # V22: ancova residual sd without cov_outcome_r; marginal back-transformation
  # 1/sqrt(1 - R^2) impossible
  ancova_resid_cols <- c("ancova_mean_sd_pooled",
                          "ancova_mean_sd_exp", "ancova_mean_sd_nexp",
                          "ancova_md_sd")
  any_ancova_sd <- rep(FALSE, n)
  for (col_name in ancova_resid_cols) {
    if (col_name %in% colnames(x)) {
      vals <- x[[col_name]]
      any_ancova_sd <- any_ancova_sd | (!is.na(vals) & vals > 0)
    }
  }
  cov_r_present <- if ("cov_outcome_r" %in% colnames(x)) {
    # back-transformation uses R^2, so negative r is fine; only |r| >= 1 invalid
    !is.na(x$cov_outcome_r) & abs(x$cov_outcome_r) < 1
  } else {
    rep(FALSE, n)
  }
  # present-but-unusable is a different user error from never-supplied, and the two
  # need different advice: one is "add the column", the other is "the value you gave
  # cannot be a correlation"
  cov_r_supplied <- if ("cov_outcome_r" %in% colnames(x)) !is.na(x$cov_outcome_r) else rep(FALSE, n)
  flagged_v22 <- which(any_ancova_sd & !cov_r_present)
  for (i in flagged_v22) {
    why <- if (isTRUE(cov_r_supplied[i])) {
      "cov_outcome_r is outside (-1, 1) and cannot be used"
    } else {
      "cov_outcome_r is missing"
    }
    row_issues[[i]] <- c(row_issues[[i]], paste0(
      "[UNUSUAL] ANCOVA residual SD provided but ", why,
      ". Residual SD is conditional on the covariate(s) - without a usable cov_outcome_r the",
      " marginal-scale back-transformation 1/sqrt(1 - R^2) cannot be applied and the",
      " standardized effect size is returned as NA. Supply the pooled within-group",
      " correlation between the outcome and the covariate(s)",
      " (Cochrane Handbook section 6.5.2.5 does not flag this step)"))
  }

  # V23: identical stat block across different study_id (templated data);
  # same study_id skipped, full block required
  if (isTRUE(enable_cross_row) && n >= 2) {
    stat_blocks <- list(
      baseline = c("mean_pre_exp", "mean_pre_sd_exp", "mean_pre_nexp", "mean_pre_sd_nexp"),
      endpoint = c("mean_exp", "mean_sd_exp", "mean_nexp", "mean_sd_nexp"),
      change   = c("mean_change_exp", "mean_change_sd_exp", "mean_change_nexp", "mean_change_sd_nexp")
    )
    sid <- if ("study_id" %in% colnames(x)) as.character(x[["study_id"]]) else rep(NA_character_, n)
    for (bl_name in names(stat_blocks)) {
      cols <- stat_blocks[[bl_name]]
      cols <- cols[cols %in% colnames(x)]
      if (length(cols) < 3) next
      mat <- suppressWarnings(matrix(as.numeric(as.matrix(x[, cols, drop = FALSE])),
                                     nrow = n, dimnames = list(NULL, cols)))
      complete <- which(rowSums(!is.finite(mat)) == 0)
      if (length(complete) < 2) next
      fp <- vapply(complete, function(i)
        paste(format(mat[i, ], digits = 15, scientific = FALSE), collapse = "|"),
        character(1))
      for (key in unique(fp[duplicated(fp)])) {
        grp  <- complete[fp == key]
        vals <- mat[grp[1], ]
        # require >= templated_min_match non-integer values
        if (sum(vals != round(vals)) < templated_min_match) next
        for (i in grp) {
          others <- setdiff(grp, i)
          others <- others[is.na(sid[i]) | is.na(sid[others]) | sid[others] != sid[i]]
          if (length(others) == 0) next
          row_issues[[i]] <- c(row_issues[[i]], sprintf(
            paste0("[INFO] Identical %s summary statistics shared with %s (%s). ",
                   "Independent studies rarely reproduce multi-decimal means/SDs exactly - ",
                   "verify these rows are not templated, cloned, or extracted from the same source"),
            bl_name, paste(.row_ref(others, sid), collapse = ", "),
            paste(sprintf("'%s'=%s", cols, format(vals, scientific = FALSE)), collapse = ", ")))
        }
      }
    }
  }

  issues <- vapply(row_issues, function(v) {
    if (length(v) == 0) return("")
    paste(v, collapse = "; ")
  }, character(1))

  if (verbose) {
    any_bad <- which(nchar(issues) > 0)
    if (length(any_bad) > 0) {
      msgs <- paste0("  Row ", any_bad, ": ", issues[any_bad])
      message("Input data validation:\n", paste(msgs, collapse = "\n"))
    }
  }

  list(data = x, issues = issues)
}


#########

#' Format a method suffix like " (from means_sd)"
#' @param method character scalar or NA/NULL
#' @return character scalar (empty string if method is unavailable)
#' @noRd
.method_suffix <- function(method) {
  if (is.null(method) || is.na(method) || !nzchar(method)) return("")
  paste0(" (from ", method, ")")
}

#' Format a min/max method suffix for cross-method flags
#' @param min_info character scalar or NA/NULL
#' @param max_info character scalar or NA/NULL
#' @return character scalar like " - min: means_sd, max: student_t"
#' @noRd
.minmax_suffix <- function(min_info, max_info) {
  parts <- character(0)
  if (!is.null(min_info) && !is.na(min_info) && nzchar(min_info))
    parts <- c(parts, paste0("min: ", min_info))
  if (!is.null(max_info) && !is.na(max_info) && nzchar(max_info))
    parts <- c(parts, paste0("max: ", max_info))
  if (length(parts) == 0) return("")
  paste0(" - ", paste(parts, collapse = ", "))
}

#' Format a cross-row reference as "study_id (row N)" or "row N"
#'
#' Used by cross-row flags (V23, Category H) so a referenced study can be
#' identified by its study_id without counting rows. Falls back to the bare
#' row number when no study_id is available for that row.
#'
#' @param idx integer vector of row indices to reference
#' @param study_id character vector of per-row study IDs (or NULL), indexed by idx
#' @return character vector of formatted references, same length as idx
#' @noRd
.row_ref <- function(idx, study_id = NULL) {
  vapply(idx, function(i) {
    sid_i <- if (!is.null(study_id) && !is.na(i) &&
                 i >= 1 && i <= length(study_id)) {
      as.character(study_id[i])
    } else {
      NA_character_
    }
    if (!is.na(sid_i) && nzchar(sid_i)) {
      sprintf("%s (row %d)", sid_i, as.integer(i))
    } else {
      sprintf("row %d", as.integer(i))
    }
  }, character(1))
}

#' TRUE when every multi-method row has min ~ -max (sign-reversal convention,
#' min ~ 1/max on the natural ratio scale)
#'
#' @param n_estimations numeric vector, number of estimation methods per row
#' @param min_es,max_es numeric vectors of min/max cross-method ES
#' @param on_exp_scale logical, TRUE when min/max are natural-scale ratios
#' @return logical scalar: TRUE if all multi-method rows are reversals
#' @noRd
.rows_all_reversed <- function(n_estimations, min_es, max_es, on_exp_scale = FALSE) {
  if (is.null(min_es) || is.null(max_es)) return(FALSE)
  multi <- which(!is.na(n_estimations) & n_estimations >= 2 &
                 !is.na(min_es) & !is.na(max_es) &
                 is.finite(min_es) & is.finite(max_es))
  if (length(multi) == 0) return(FALSE)
  rev <- vapply(multi, function(i) {
    lo <- min_es[i]; hi <- max_es[i]
    if (on_exp_scale) {
      if (lo <= 0 || hi <= 0) return(FALSE)
      lo <- log(lo); hi <- log(hi)
    }
    mag <- max(abs(lo), abs(hi))
    mag > 0 && abs(lo + hi) < 0.01 * mag
  }, logical(1))
  all(rev)
}


#########

#' Check for numeric integrity issues (Inf, NaN, negative SE, CI problems)
#' @param es numeric vector of effect sizes
#' @param se numeric vector of standard errors
#' @param ci_lo numeric vector of CI lower bounds
#' @param ci_up numeric vector of CI upper bounds
#' @param info_used character vector of method names (optional, for traceability)
#' @return list of character vectors (one per row)
#' @noRd
.flag_numeric_integrity <- function(es, se, ci_lo, ci_up, info_used = NULL,
                                     measure = NULL, exp = FALSE,
                                     n_total = NULL, n_exp = NULL, n_nexp = NULL,
                                     enable_informational = FALSE,
                                     r_defaulted = NULL,
                                     prop_to_es = "raw",
                                     smd_denom = NULL,
                                     n_cov_ancova = NULL) {
  n <- length(es)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  for (i in seq_len(n)) {
    msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""
    # r-sensitivity note for the SE-driven A6 CI-width check (see .r_sensitive_note)
    rsuf <- if (!is.null(r_defaulted) && length(r_defaulted) >= i) {
      .r_sensitive_note(r_defaulted[i])
    } else {
      ""
    }
    # A1: es non-finite
    if (is.nan(es[i]) || (!is.na(es[i]) && is.infinite(es[i]))) {
      flags[[i]] <- c(flags[[i]], paste0("[INVALID] ES is non-finite or undefined", msuf))
    }
    # A2: SE is Inf or NaN
    if (is.nan(se[i]) || (!is.na(se[i]) && is.infinite(se[i]))) {
      flags[[i]] <- c(flags[[i]], paste0("[INVALID] SE is non-finite or undefined", msuf))
    }
    # A3: SE is negative
    if (!is.na(se[i]) && is.finite(se[i]) && se[i] < 0) {
      flags[[i]] <- c(flags[[i]], paste0("[INVALID] Negative SE: ", round(se[i], 3), msuf))
    }
    # A4: CI lower > CI upper
    if (!is.na(ci_lo[i]) && !is.na(ci_up[i]) &&
        is.finite(ci_lo[i]) && is.finite(ci_up[i]) &&
        ci_lo[i] > ci_up[i]) {
      flags[[i]] <- c(flags[[i]], paste0("[INVALID] Inverted CI: lower > upper", msuf))
    }
    # A5: ES not within CI bounds
    if (!is.na(es[i]) && !is.na(ci_lo[i]) && !is.na(ci_up[i]) &&
        is.finite(es[i]) && is.finite(ci_lo[i]) && is.finite(ci_up[i]) &&
        ci_lo[i] <= ci_up[i]) {
      if (es[i] < ci_lo[i] || es[i] > ci_up[i]) {
        flags[[i]] <- c(flags[[i]], paste0("[INVALID] ES outside its CI bounds", msuf))
      }
    }

    # A6: ci width vs se
    # skipped for nnt (discontinuous CI) and exp=TRUE (se stays on the log scale)
    skip_a6 <- identical(measure, "nnt") || exp
    if (!skip_a6 &&
        !is.na(es[i]) && !is.na(se[i]) && !is.na(ci_lo[i]) && !is.na(ci_up[i]) &&
        is.finite(es[i]) && is.finite(se[i]) && is.finite(ci_lo[i]) && is.finite(ci_up[i]) &&
        se[i] > 0 && ci_lo[i] <= ci_up[i]) {
      ci_width <- ci_up[i] - ci_lo[i]
      # package CIs: qt for d/g/md/r and dw/gw/mdw, qnorm otherwise
      # user CIs: checked against z width, t width also accepted
      qt_measures <- c("d", "g", "dw", "gw", "md", "mdw", "r")
      n_i <- if (!is.null(n_total)) n_total[i] else NA_real_
      # paired measures: df = n_pairs - 1 (n_total double counts the subjects)
      within_measures <- c("dw", "gw", "mdw")
      n_exp_i <- if (!is.null(n_exp)) n_exp[i] else NA_real_
      n_nexp_i <- if (!is.null(n_nexp)) n_nexp[i] else NA_real_
      if (!is.null(measure) && measure %in% within_measures &&
          !is.na(n_exp_i) && is.finite(n_exp_i) && n_exp_i > 2) {
        t_df <- n_exp_i - 1
        n_eff <- n_exp_i
      } else if (!is.na(n_i) && n_i > 4) {
        t_df <- n_i - 2
        n_eff <- n_i
      } else {
        t_df <- NA_real_
        n_eff <- n_i
      }
      is_user_es <- !is.null(info_used) && !is.na(info_used[i]) &&
                    grepl("^user_(es|input)", info_used[i])
      # P7: Glass rows (smd_denom = control/control_robust) build their d/g CIs
      # on the control-arm df (SMD1 convention), not the pooled N - 2. The Glass df
      # applies ONLY when the SELECTED estimate actually came from the endpoint-means
      # family -- the only methods that honour smd_denom. When the hierarchy picked a
      # non-means method (cohen_d, etasq, t/F, ANCOVA, medians...), the CI is built on
      # the pooled N - 2 regardless of the smd_denom input column, so keying the Glass
      # df off the raw input alone would false-fire the A6 CI-width flag on valid rows.
      glass_row <- !is_user_es && !is.null(smd_denom) &&
        length(smd_denom) >= i && !is.na(smd_denom[i]) &&
        smd_denom[i] %in% c("control", "control_robust") &&
        !is.null(measure) && measure %in% c("d", "g") &&
        !is.null(info_used) && !is.na(info_used[i]) &&
        info_used[i] %in% c("means_sd", "means_se", "means_ci") &&
        !is.na(n_nexp_i) && is.finite(n_nexp_i) && n_nexp_i > 2
      if (glass_row) t_df <- n_nexp_i - 1
      # P9: an ADJUSTED row's d/g interval is built by .es_from_d() on
      # qt(.975, n_exp + n_nexp - 2 - n_cov_ancova), so A6 must expect that df too --
      # otherwise it compares the row's own CI against a wider-df expectation and
      # fires [DISCORDANT] on a perfectly consistent interval. The gap only clears the
      # tolerance at small n (17% against a 10% tolerance at N = 9, q = 3), but it is a
      # pure false positive when it does. The adjusted routes are exactly those whose
      # info_used starts with "ancova" or ends in "_adj" (see the es_from_ancova_*,
      # es_from_cohen_d_adj and es_from_etasq_adj entry points); user-entered rows are
      # excluded because the source's own CI construction is unknown and is checked
      # against the Wald-z width instead.
      if (!is_user_es && !is.null(info_used) && !is.na(info_used[i]) &&
          grepl("^ancova|_adj$", info_used[i]) &&
          !is.null(n_cov_ancova) && length(n_cov_ancova) >= i &&
          !is.na(n_cov_ancova[i]) && is.finite(n_cov_ancova[i]) &&
          n_cov_ancova[i] > 0 && !is.na(t_df) && t_df - n_cov_ancova[i] > 0) {
        t_df <- t_df - n_cov_ancova[i]
      }
      if (is_user_es) {
        z_crit <- stats::qnorm(0.975)
      } else if (!is.null(measure) && measure %in% qt_measures &&
                 !is.na(t_df)) {
        z_crit <- stats::qt(0.975, t_df)
      } else {
        z_crit <- stats::qnorm(0.975)
      }
      expected_width <- 2 * z_crit * se[i]
      # ~10% at N=9, ~5% at N=36
      tol <- if (!is.na(n_eff) && n_eff > 0) max(0.05, 0.30 / sqrt(n_eff)) else 0.10
      # P8: raw-scale proportion CIs are deliberately clamped to [0, 1] by
      # es_from_prop_single_group(); when a bound sits at the clamp, compare
      # against the [0, 1]-clipped expected width instead of 2 * z * se.
      prop_clamped <- !is_user_es && identical(measure, "prop") &&
        identical(prop_to_es, "raw") &&
        (isTRUE(all.equal(ci_lo[i], 0)) || isTRUE(all.equal(ci_up[i], 1)))
      if (prop_clamped) {
        expected_width <- min(es[i] + z_crit * se[i], 1) -
          max(es[i] - z_crit * se[i], 0)
      }
      fires_a6 <- expected_width > 0 &&
        abs(ci_width - expected_width) / expected_width > tol
      # P6: package-computed md CIs legitimately use any df between the Welch
      # minimum (min(n1, n2) - 1; es_from_means_sd) and the pooled N - 2
      # (es_from_md_*, ANCOVA md) -- accept the whole band before flagging.
      if (fires_a6 && !is_user_es && identical(measure, "md") &&
          !is.na(n_exp_i) && !is.na(n_nexp_i) &&
          is.finite(n_exp_i) && is.finite(n_nexp_i) &&
          min(n_exp_i, n_nexp_i) > 2 && !is.na(t_df)) {
        w_pooled <- 2 * stats::qt(0.975, n_exp_i + n_nexp_i - 2) * se[i]
        w_welch_min <- 2 * stats::qt(0.975, min(n_exp_i, n_nexp_i) - 1) * se[i]
        fires_a6 <- ci_width < (1 - tol) * w_pooled ||
          ci_width > (1 + tol) * w_welch_min
        # report the closest band edge in the message
        expected_width <- if (ci_width < w_pooled) w_pooled else w_welch_min
      }
      if (fires_a6) {
        matches_t <- FALSE
        if (is_user_es && !is.null(measure) && measure %in% qt_measures &&
            !is.na(t_df)) {
          t_width <- 2 * stats::qt(0.975, t_df) * se[i]
          matches_t <- t_width > 0 &&
            abs(ci_width - t_width) / t_width <= tol
        }
        if (matches_t) {
          if (enable_informational) {
            flags[[i]] <- c(flags[[i]],
              sprintf("[INFO] CI width matches a t-based CI rather than Wald-z: z-expected %.3f, actual %.3f. Benign if the source built t CIs. If the source used Wald CIs, check the SE formula (e.g. paired vs independent)%s",
                      expected_width, ci_width, msuf))
          }
        } else {
          flags[[i]] <- c(flags[[i]],
            sprintf("[DISCORDANT] CI width inconsistent with SE: expected %.3f, actual %.3f%s%s",
                    expected_width, ci_width, msuf, rsuf))
        }
      }
    }
  }
  return(flags)
}


#########

#' Check for theoretically impossible values based on the effect size measure
#' @param es numeric vector of effect sizes
#' @param se numeric vector of standard errors
#' @param ci_lo numeric vector of CI lower bounds
#' @param ci_up numeric vector of CI upper bounds
#' @param measure character, the effect size measure
#' @param exp logical, whether OR/RR/IRR are on the exp (natural) scale
#' @param info_used character vector of method names (optional, for traceability)
#' @return list of character vectors (one per row)
#' @noRd
.flag_bounds_violations <- function(es, se, ci_lo, ci_up, measure, exp,
                                     info_used = NULL,
                                     baseline_risk = NULL,
                                     baseline_rate = NULL,
                                     alpha_to_es = "bonett",
                                     icc_to_es = "bonett",
                                     prop_to_es = "raw") {
  n <- length(es)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  for (i in seq_len(n)) {
    if (is.na(es[i])) next
    msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""

    # B1: Correlation outside [-1, 1]
    if (measure == "r" && is.finite(es[i]) && abs(es[i]) > 1) {
      flags[[i]] <- c(flags[[i]],
        paste0("[INVALID] r outside [-1, 1]: r = ", round(es[i], 3), msuf))
    }

    # B1b: raw-scale Wald CI escaping [-1, 1] while the point estimate is valid.
    # Same geometry as B7b/B8b for alpha and ICC: r is bounded identically and its
    # interval is built identically (es +/- z*se), so a symmetric Wald interval can
    # leave the parameter space near |r| = 1 even when r itself is fine. Common on
    # the tetrachoric route with small off-diagonal cells; B1 above tests only the
    # point estimate, so without this the bound passes through silently.
    #
    # [INFO], not [UNUSUAL], deliberately. The data, the table and the point
    # estimate are all sound -- there is nothing to verify and nothing extracted
    # wrongly. This is a property of the interval metaConvert itself constructs, and
    # the user has no corrective action available (table_2x2_to_cor accepts only
    # "tetrachoric"). That is V31's situation, not B7b/B8b's: B7b/B8b are [UNUSUAL]
    # because they can point the user at alpha_to_es/icc_to_es = "bonett", which
    # cannot overshoot. Here there is no such escape, so the flag is a disclosure.
    if (measure == "r" && is.finite(es[i]) && abs(es[i]) <= 1) {
      lo_out <- !is.na(ci_lo[i]) && is.finite(ci_lo[i]) && ci_lo[i] < -1
      up_out <- !is.na(ci_up[i]) && is.finite(ci_up[i]) && ci_up[i] > 1
      if (lo_out || up_out) {
        bound <- if (up_out) ci_up[i] else ci_lo[i]
        side <- if (up_out) "upper bound exceeds 1" else "lower bound is below -1"
        flags[[i]] <- c(flags[[i]],
          paste0("[INFO] Correlation CI ", side, " (", round(bound, 3),
                 ") while r is valid - the symmetric Wald interval escapes the parameter space, interpret the interval with caution", msuf))
      }
    }

    # B2: OR/RR/IRR non-positive on natural scale
    if (exp && measure %in% c("logor", "logrr", "logirr", "loghr") &&
        is.finite(es[i]) && es[i] <= 0) {
      flags[[i]] <- c(flags[[i]],
        paste0("[INVALID] ", toupper(sub("log", "", measure)),
               " is non-positive: ", round(es[i], 3), msuf))
    }

    # B2b: CI bound non-positive on natural scale (Riley 2023)
    if (exp && measure %in% c("logor", "logrr", "logirr", "loghr")) {
      bound_problems <- character(0)
      if (!is.na(ci_lo[i]) && is.finite(ci_lo[i]) && ci_lo[i] <= 0) {
        bound_problems <- c(bound_problems,
                            paste0("CI lower = ", round(ci_lo[i], 3)))
      }
      if (!is.na(ci_up[i]) && is.finite(ci_up[i]) && ci_up[i] <= 0) {
        bound_problems <- c(bound_problems,
                            paste0("CI upper = ", round(ci_up[i], 3)))
      }
      if (length(bound_problems) > 0) {
        flags[[i]] <- c(flags[[i]],
          paste0("[INVALID] ", toupper(sub("log", "", measure)),
                 " CI bound non-positive: ",
                 paste(bound_problems, collapse = ", "), msuf))
      }
    }

    # B3: NNT magnitude < 1 (mathematically impossible: cannot treat fewer
    # than one person to prevent one event)
    if (measure == "nnt" && is.finite(es[i]) &&
        abs(es[i]) < 1 && es[i] != 0) {
      flags[[i]] <- c(flags[[i]],
        paste0("[INVALID] NNT magnitude < 1: NNT = ", round(es[i], 3), msuf))
    }

    # B4: Proportion outside [0, 1] (raw scale only - transformed scales are unbounded)
    if (measure == "prop" && prop_to_es == "raw" && is.finite(es[i]) &&
        (es[i] < 0 || es[i] > 1)) {
      flags[[i]] <- c(flags[[i]],
        paste0("[INVALID] Proportion outside [0, 1]: prop = ", round(es[i], 3), msuf))
    }

    # B5: NNT below minimum possible given baseline risk
    if (measure == "nnt" && is.finite(es[i]) && abs(es[i]) > 0) {
      br <- if (!is.null(baseline_risk) && length(baseline_risk) >= i) baseline_risk[i] else NA
      brate <- if (!is.null(baseline_rate) && length(baseline_rate) >= i) baseline_rate[i] else NA

      if (!is.na(br) && br > 0 && br < 1) {
        min_nnt <- 1 / br
        if (abs(es[i]) < min_nnt) {
          flags[[i]] <- c(flags[[i]],
            paste0("[INVALID] Minimum possible NNT = ", round(min_nnt, 1),
                   " (baseline_risk = ", round(br * 100, 1), "%)",
                   " - reported |NNT| = ", round(abs(es[i]), 3), msuf))
        }
      }
      if (!is.na(brate) && brate > 0) {
        min_nnt_pt <- 1 / brate
        if (abs(es[i]) < min_nnt_pt) {
          flags[[i]] <- c(flags[[i]],
            paste0("[INVALID] Minimum possible person-time NNT = ", round(min_nnt_pt, 1),
                   " (baseline_rate = ", round(brate, 4), ")",
                   " - reported |NNT| = ", round(abs(es[i]), 3), msuf))
        }
      }
    }

    # B6: RD outside [-1, 1]
    if (measure == "rd" && is.finite(es[i]) && abs(es[i]) > 1) {
      flags[[i]] <- c(flags[[i]],
        paste0("[INVALID] RD outside [-1, 1]: rd = ", round(es[i], 3), msuf))
    }

    # B6b: RD CI bound outside [-1, 1]
    if (measure == "rd") {
      bound_problems <- character(0)
      if (!is.na(ci_lo[i]) && is.finite(ci_lo[i]) && abs(ci_lo[i]) > 1) {
        bound_problems <- c(bound_problems,
                            paste0("CI lower = ", round(ci_lo[i], 3)))
      }
      if (!is.na(ci_up[i]) && is.finite(ci_up[i]) && abs(ci_up[i]) > 1) {
        bound_problems <- c(bound_problems,
                            paste0("CI upper = ", round(ci_up[i], 3)))
      }
      if (length(bound_problems) > 0) {
        flags[[i]] <- c(flags[[i]],
          paste0("[INVALID] RD CI bound outside [-1, 1]: ",
                 paste(bound_problems, collapse = ", "), msuf))
      }
    }

    # B7: alpha bounds
    # bonett scale: es > 0 means negative alpha, possible but suspect
    if (measure == "alpha" && is.finite(es[i])) {
      if (identical(alpha_to_es, "bonett")) {
        if (es[i] > 0) {
          flags[[i]] <- c(flags[[i]],
            paste0("[UNUSUAL] Bonett alpha ln(1-alpha) = ", round(es[i], 3),
                   " > 0 implies a negative Cronbach's alpha - possible but indicates serious measurement problems, verify extraction", msuf))
        }
      } else if (es[i] > 1) {
        flags[[i]] <- c(flags[[i]],
          paste0("[INVALID] Alpha exceeds 1: alpha = ", round(es[i], 3), msuf))
      }
      # B7b: raw-scale Wald CI escaping the upper bound (the point estimate is
      # valid, the symmetric interval is not; the bonett scale cannot overshoot)
      if (!identical(alpha_to_es, "bonett") &&
          !is.na(ci_up[i]) && is.finite(ci_up[i]) && ci_up[i] > 1 &&
          es[i] <= 1) {
        flags[[i]] <- c(flags[[i]],
          paste0("[UNUSUAL] Raw-scale alpha CI upper bound exceeds 1 (",
                 round(ci_up[i], 3),
                 ") - the symmetric Wald interval escapes the parameter space, consider alpha_to_es = 'bonett'", msuf))
      }
    }

    # B8: ICC bounds
    if (measure == "icc" && is.finite(es[i])) {
      if (identical(icc_to_es, "bonett")) {
        if (es[i] > 0) {
          flags[[i]] <- c(flags[[i]],
            paste0("[UNUSUAL] Bonett ICC ln(1-ICC) = ", round(es[i], 3),
                   " > 0 implies a negative ICC - possible but indicates serious measurement problems, verify extraction", msuf))
        }
      } else if (es[i] > 1) {
        flags[[i]] <- c(flags[[i]],
          paste0("[INVALID] ICC exceeds 1: icc = ", round(es[i], 3), msuf))
      }
      # B8b: raw-scale Wald CI escaping the parameter space [-1, 1]
      if (!identical(icc_to_es, "bonett") &&
          !is.na(es[i]) && es[i] <= 1 && es[i] >= -1) {
        icc_bound_problems <- character(0)
        if (!is.na(ci_up[i]) && is.finite(ci_up[i]) && ci_up[i] > 1) {
          icc_bound_problems <- c(icc_bound_problems,
                                  paste0("CI upper = ", round(ci_up[i], 3)))
        }
        if (!is.na(ci_lo[i]) && is.finite(ci_lo[i]) && ci_lo[i] < -1) {
          icc_bound_problems <- c(icc_bound_problems,
                                  paste0("CI lower = ", round(ci_lo[i], 3)))
        }
        if (length(icc_bound_problems) > 0) {
          flags[[i]] <- c(flags[[i]],
            paste0("[UNUSUAL] Raw-scale ICC CI bound outside [-1, 1]: ",
                   paste(icc_bound_problems, collapse = ", "),
                   " - the symmetric Wald interval escapes the parameter space, consider icc_to_es = 'bonett'", msuf))
        }
      }
    }
  }
  return(flags)
}


#########

#' Check for statistically unusual but not impossible values
#' @param es numeric vector of effect sizes
#' @param se numeric vector of standard errors
#' @param measure character, the effect size measure
#' @param opts list of flag options (thresholds)
#' @param n_sample numeric vector of sample sizes (can contain NA)
#' @param info_used character vector of method names (optional, for traceability)
#' @param n_exp numeric vector of experimental group sizes (optional, for C10)
#' @param n_nexp numeric vector of control group sizes (optional, for C10)
#' @return list of character vectors (one per row)
#' @noRd
.flag_plausibility <- function(es, se, measure, opts, n_sample,
                                info_used = NULL,
                                n_exp = NULL, n_nexp = NULL,
                                exp = FALSE,
                                prop_to_es = "raw") {
  n <- length(es)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  for (i in seq_len(n)) {
    msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""

    # C1: SMD magnitude too large
    if (measure %in% c("d", "g", "dw", "gw") &&
        !is.na(es[i]) && is.finite(es[i]) &&
        abs(es[i]) > opts$smd_max) {
      flags[[i]] <- c(flags[[i]],
        paste0("[UNUSUAL] Large SMD: |", measure, "| = ",
               round(abs(es[i]), 3), " (threshold: ", opts$smd_max, ")", msuf))
    }

    # C2: SE equals zero
    if (!is.na(se[i]) && is.finite(se[i]) && se[i] == 0) {
      flags[[i]] <- c(flags[[i]], paste0("[INVALID] SE is zero", msuf))
    }

    # C4: Correlation magnitude very high
    if (measure == "r" &&
        !is.na(es[i]) && is.finite(es[i]) &&
        abs(es[i]) > opts$r_max && abs(es[i]) <= 1) {
      flags[[i]] <- c(flags[[i]],
        paste0("[UNUSUAL] High correlation: |r| = ",
               round(abs(es[i]), 3), " (threshold: ", opts$r_max, ")", msuf))
    }

    # C5: Log OR/RR/IRR/HR magnitude very large
    # threshold is on the log scale, transform back when exp = TRUE
    if (measure %in% c("logor", "logrr", "logirr", "loghr") &&
        !is.na(es[i]) && is.finite(es[i])) {
      log_val <- if (exp) {
        if (es[i] > 0) log(es[i]) else NA_real_  # es <= 0 is handled by B2
      } else {
        es[i]
      }
      if (!is.na(log_val) && abs(log_val) > opts$log_or_max) {
        flags[[i]] <- c(flags[[i]],
          paste0("[UNUSUAL] Large log ", toupper(sub("log", "", measure)),
                 ": |value| = ", round(abs(log_val), 3),
                 " (threshold: ", opts$log_or_max, ")", msuf))
      }
    }

    # C6: Very small sample size (informational - gated by enable_informational)
    if (isTRUE(opts$enable_informational) &&
        !is.na(n_sample[i]) && is.finite(n_sample[i]) &&
        n_sample[i] < opts$n_min) {
      flags[[i]] <- c(flags[[i]],
        paste0("[INFO] Small sample size: n = ", n_sample[i],
               " (threshold: ", opts$n_min, ")"))
    }

    # C7: Near-perfect Cronbach's alpha (raw scale only: ES in (0, 1])
    if (measure == "alpha" &&
        !is.na(es[i]) && is.finite(es[i]) &&
        es[i] > 0 && es[i] <= 1 && es[i] > opts$alpha_max) {
      flags[[i]] <- c(flags[[i]],
        paste0("[UNUSUAL] Near-perfect alpha: ", round(es[i], 3),
               " (threshold: ", opts$alpha_max, ")", msuf))
    }

    # C8: Near-perfect ICC (raw scale only: ES in (0, 1])
    if (measure == "icc" &&
        !is.na(es[i]) && is.finite(es[i]) &&
        es[i] > 0 && es[i] <= 1 && es[i] > opts$icc_max) {
      flags[[i]] <- c(flags[[i]],
        paste0("[UNUSUAL] Near-perfect ICC: ", round(es[i], 3),
               " (threshold: ", opts$icc_max, ")", msuf))
    }

    # C9: Extreme transformed proportion (logit or Freeman-Tukey scale)
    # On raw scale, B4 catches out-of-bounds. This flag targets transformed scales only.
    if (measure == "prop" && prop_to_es != "raw" &&
        !is.na(es[i]) && is.finite(es[i])) {
      is_extreme <- FALSE
      # Logit scale: |logit| > 6 corresponds to p < 0.0025 or p > 0.9975
      if (prop_to_es == "logit" && abs(es[i]) > 6) is_extreme <- TRUE
      # Freeman-Tukey scale: bounded to [0, pi/2 ~ 1.5708].
      # Values near 0 or near pi/2 indicate extreme proportions.
      # es < 0.05 ~ p < 0.002, es > 1.52 ~ p > 0.998
      if (prop_to_es == "freeman_tukey" &&
          (es[i] < 0.05 || es[i] > 1.52)) is_extreme <- TRUE
      if (is_extreme) {
        flags[[i]] <- c(flags[[i]],
          paste0("[UNUSUAL] Extreme transformed proportion: value = ",
                 round(es[i], 3), msuf))
      }
    }

    # C10: Extreme group size imbalance (two-group measures only)
    if (measure %in% c("d", "g", "dw", "gw", "logor", "logrr", "nnt", "rd") &&
        !is.null(n_exp) && !is.null(n_nexp) &&
        length(n_exp) >= i && length(n_nexp) >= i &&
        !is.na(n_exp[i]) && !is.na(n_nexp[i]) &&
        is.finite(n_exp[i]) && is.finite(n_nexp[i]) &&
        n_exp[i] > 0 && n_nexp[i] > 0) {
      grp_ratio <- max(n_exp[i], n_nexp[i]) / min(n_exp[i], n_nexp[i])
      if (grp_ratio > opts$group_ratio_max) {
        flags[[i]] <- c(flags[[i]],
          sprintf("[UNUSUAL] Extreme group size imbalance: %d vs %d (ratio: %.1f, threshold: %g)",
                  as.integer(n_exp[i]), as.integer(n_nexp[i]),
                  grp_ratio, opts$group_ratio_max))
      }
    }
  }
  return(flags)
}


#########

#' Detect cross-row outliers using IQR method
#' @param es numeric vector of effect sizes
#' @param se numeric vector of standard errors
#' @param opts list of flag options (thresholds)
#' @param info_used character vector of method names (optional, for traceability)
#' @param measure character, the effect size measure (optional, for value labels)
#' @return list of character vectors (one per row)
#' @noRd
.flag_cross_row_outliers <- function(es, se, opts, info_used = NULL,
                                      measure = NULL, n_total = NULL,
                                      n_exp = NULL, n_nexp = NULL,
                                      sd_exp = NULL, sd_nexp = NULL,
                                      suffix = "", exp = FALSE) {
  n <- length(es)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  if (!isTRUE(opts$enable_cross_row)) return(flags)

  # outlier decision on the log scale when exp = TRUE, messages keep the natural value
  log_measures <- c("logor", "logrr", "logirr", "loghr")
  on_exp_scale <- isTRUE(exp) && !is.null(measure) && measure %in% log_measures
  es_work <- es
  if (on_exp_scale) {
    es_work <- ifelse(!is.na(es) & is.finite(es) & es > 0, log(es), NA_real_)
  }

  # Label for the ES value (e.g. "d" or "ES"; natural-scale name when exp=TRUE)
  es_label <- if (!is.null(measure) && nzchar(measure)) measure else "ES"
  if (on_exp_scale) es_label <- sub("^log", "", es_label)

  # minimum deviation from median, avoids flagging homogeneous pools
  if (!is.null(opts$outlier_min_deviation)) {
    min_dev <- opts$outlier_min_deviation
  } else if (is.null(measure) || !nzchar(measure)) {
    min_dev <- 0
  } else {
    min_dev <- switch(measure,
      "d" = , "g" = 1.0, # 1 SD
      "logor" = , "logrr" = , "logirr" = , "loghr" = 1.0, # log scale
      "r" = 0.3,
      "z" = 0.5,
      "rd" = 0.3,
      0 # no floor for md, nnt, prop, alpha, icc
    )
  }

  # D1: ES outlier (decision on es_work - log scale for exp-displayed ratios)
  es_valid <- es_work[!is.na(es_work) & is.finite(es_work)]
  if (length(es_valid) >= 4) {
    es_med <- stats::median(es_valid)
    es_iqr <- stats::IQR(es_valid)
    if (es_iqr > 0) {
      for (i in seq_len(n)) {
        if (!is.na(es_work[i]) && is.finite(es_work[i]) &&
            abs(es_work[i] - es_med) > opts$iqr_mult * es_iqr &&
            abs(es_work[i] - es_med) > min_dev) {
          msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""
          flags[[i]] <- c(flags[[i]],
            paste0("[UNUSUAL] ES outlier, IQR method: ", es_label, " = ",
                   round(es[i], 3), msuf))
        }
      }
    }
  }

  # D2: SE outlier. normalise SE before the IQR: d/g by expected SE,
  # logOR/logRR by 2x2 reconstruction, else sqrt(N)
  se_for_iqr <- se
  if (measure %in% c("d", "g") && !is.null(n_exp) && !is.null(n_nexp)) {
    has_both <- !is.na(n_exp) & !is.na(n_nexp) & is.finite(n_exp) & is.finite(n_nexp) &
                n_exp > 0 & n_nexp > 0 & !is.na(es) & is.finite(es) & !is.na(se) & is.finite(se)
    if (any(has_both)) {
      expected_se <- sqrt(1/n_exp[has_both] + 1/n_nexp[has_both] +
                          es[has_both]^2 / (2 * (n_exp[has_both] + n_nexp[has_both])))
      expected_se[expected_se == 0] <- NA_real_
      se_for_iqr[has_both] <- se[has_both] / expected_se
    }
    # Balanced-arms 2/sqrt(N) fallback for rows WITHOUT both arm sizes. This must run
    # whether or not ANY row had arm sizes: a pure pearson_r/fisher_z/spearman_r d or g
    # pool carries only n_sample (n_exp/n_nexp stay NA in raw_data), so nesting this inside
    # if(any(has_both)) left every row with its raw, N-dependent SE and falsely flagged a
    # legitimately small-N study as an SE outlier. Mirrors the branch-level no_data
    # fallback in the logOR/logRR branch below.
    not_norm <- !has_both & !is.na(se) & is.finite(se)
    if (any(not_norm)) {
      if (!is.null(n_total)) {
        approx <- not_norm & !is.na(n_total) & is.finite(n_total) & n_total > 0
        se_for_iqr[approx] <- se[approx] / (2 / sqrt(n_total[approx]))
        # only null out un-normalisable rows once SOME row has been normalised, to avoid
        # collapsing the whole pool to NA (matches the logOR branch's any()-guard)
        if (any(has_both) || any(approx)) se_for_iqr[not_norm & !approx] <- NA_real_
      } else if (any(has_both)) {
        se_for_iqr[not_norm] <- NA_real_
      }
    }
  } else if (measure %in% c("logor", "logrr") &&
             !is.null(n_exp) && !is.null(n_nexp) &&
             suffix != "_adjusted") {
    # 2x2 reconstruction accounts for event-rate variation across studies
    # skip adjusted estimates: multivariable SEs are larger than 2x2-implied SEs
    # 4/sqrt(N) ~ SE_logOR at a 50% event rate (puts fallback rows on the same ~1 scale)
    balanced_se_const <- 4
    has_data <- !is.na(n_exp) & !is.na(n_nexp) & is.finite(n_exp) & is.finite(n_nexp) &
                n_exp > 1 & n_nexp > 1 & !is.na(es) & is.finite(es) &
                !is.na(se) & is.finite(se) & se > 0
    if (any(has_data)) {
      for (j in which(has_data)) {
        # es is already the natural-scale OR when exp=TRUE
        or_j <- if (on_exp_scale) es[j] else exp(es[j])
        var_j <- se[j]^2
        n_exp_j <- as.integer(round(n_exp[j]))
        n_nexp_j <- as.integer(round(n_nexp[j]))
        # Skip if group sizes are too large (simulation is O(n_exp))
        if (n_exp_j > 5000 || n_nexp_j > 5000) {
          # Fallback: sqrt(N) normalization, scaled to the reconstruction (~1) pool
          se_for_iqr[j] <- se[j] * sqrt(n_exp_j + n_nexp_j) / balanced_se_const
          next
        }
        recon <- tryCatch(
          .estimate_n_from_or_and_n_exp(or_j, var_j, n_exp_j, n_nexp_j),
          error = function(e) NULL
        )
        if (!is.null(recon) && nrow(recon) > 0 &&
            all(is.finite(c(recon$n_cases_exp[1], recon$n_cases_nexp[1],
                            recon$n_controls_exp[1], recon$n_controls_nexp[1])))) {
          a <- recon$n_cases_exp[1]; b <- recon$n_cases_nexp[1]
          c_val <- recon$n_controls_exp[1]; d <- recon$n_controls_nexp[1]
          if (measure == "logor") {
            # logOR SE = sqrt(1/a + 1/b + 1/c + 1/d)
            cells <- c(a, b, c_val, d)
            if (any(cells == 0)) cells <- cells + 0.5  # continuity correction
            exp_se <- sqrt(sum(1 / cells))
          } else {
            # logRR SE = sqrt(1/a - 1/(a+c) + 1/b - 1/(b+d))
            if (a > 0 && b > 0 && (a + c_val) > 0 && (b + d) > 0) {
              exp_se <- sqrt(1/a - 1/(a + c_val) + 1/b - 1/(b + d))
            } else {
              exp_se <- NA_real_
            }
          }
          if (!is.na(exp_se) && exp_se > 0) {
            se_for_iqr[j] <- se[j] / exp_se
          } else {
            # Fallback (scaled to the reconstruction ~1 pool)
            se_for_iqr[j] <- se[j] * sqrt(n_exp_j + n_nexp_j) / balanced_se_const
          }
        } else {
          # Fallback (scaled to the reconstruction ~1 pool)
          se_for_iqr[j] <- se[j] * sqrt(n_exp_j + n_nexp_j) / balanced_se_const
        }
      }
    }
    # rows without arm counts: sqrt(N) fallback; no N at all -> excluded
    no_data <- !has_data & !is.na(se) & is.finite(se)
    if (any(no_data)) {
      if (!is.null(n_total)) {
        has_n <- no_data & !is.na(n_total) & is.finite(n_total) & n_total > 0
        se_for_iqr[has_n] <- se[has_n] * sqrt(n_total[has_n]) / balanced_se_const
        if (any(has_data) || any(has_n)) se_for_iqr[no_data & !has_n] <- NA_real_
      } else if (any(has_data)) {
        se_for_iqr[no_data] <- NA_real_
      }
    }
  } else if (!is.null(n_total)) {
    has_n <- !is.na(n_total) & is.finite(n_total) & n_total > 0
    if (any(has_n)) {
      se_for_iqr[has_n] <- se[has_n] * sqrt(n_total[has_n])
      se_for_iqr[!has_n] <- NA_real_
    }
  }
  # SE <= 0 is invalid (C2/A3 flag it) and must not anchor the outlier pool:
  # a zero would sit in the median/IQR as a spurious extreme-low observation
  se_for_iqr[!is.na(se) & is.finite(se) & se <= 0] <- NA_real_
  se_iqr_valid <- se_for_iqr[!is.na(se_for_iqr) & is.finite(se_for_iqr)]
  se_ratio_min <- opts$se_outlier_min_ratio
  if (length(se_iqr_valid) >= 4) {
    se_med <- stats::median(se_iqr_valid)
    se_iqr <- stats::IQR(se_iqr_valid)
    if (se_iqr > 0 && se_med > 0) {
      for (i in seq_len(n)) {
        if (!is.na(se_for_iqr[i]) && is.finite(se_for_iqr[i]) &&
            abs(se_for_iqr[i] - se_med) > opts$iqr_mult * se_iqr) {
          # require at least se_outlier_min_ratio x from the median
          ratio <- se_for_iqr[i] / se_med
          if (ratio > 1 / se_ratio_min && ratio < se_ratio_min) next
          msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""
          flags[[i]] <- c(flags[[i]],
            paste0("[UNUSUAL] SE outlier, IQR method: SE = ",
                   round(se[i], 3), msuf))
        }
      }
    }
  }

  # D2b: ratio-to-median fallback, catches scale mixing the IQR misses
  if (length(se_iqr_valid) >= 4 && se_med > 0) {
    se_ratio_extreme <- opts$se_ratio_extreme
    # decision on the normalised SE, fold-difference displayed on the raw SE scale
    se_raw_med <- stats::median(se[!is.na(se) & is.finite(se) & se > 0])
    for (i in seq_len(n)) {
      if (any(grepl("SE outlier", flags[[i]]))) next
      if (!is.na(se_for_iqr[i]) && is.finite(se_for_iqr[i])) {
        ratio <- se_for_iqr[i] / se_med
        if (ratio < 1 / se_ratio_extreme || ratio > se_ratio_extreme) {
          msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""
          if (!is.na(se_raw_med) && se_raw_med > 0 && se[i] < se_raw_med) {
            fold <- se_raw_med / se[i]
            cmp <- if (is.finite(fold)) {
              paste0(round(fold, 1), "x tighter than the cohort median SE (",
                     round(se_raw_med, 3), ")")
            } else {
              paste0("far tighter than the cohort median SE (",
                     round(se_raw_med, 3), ")")
            }
          } else if (!is.na(se_raw_med) && se_raw_med > 0) {
            cmp <- paste0(round(se[i] / se_raw_med, 1),
                          "x the cohort median SE (", round(se_raw_med, 3), ")")
          } else {
            cmp <- "far from the cohort SE distribution"
          }
          flags[[i]] <- c(flags[[i]],
            paste0("[UNUSUAL] SE outlier, ratio method: SE = ",
                   round(se[i], 3), " (", cmp, ")", msuf))
        }
      }
    }
  }

  # D3 (V27): spread = pooled arm SD or se*sqrt(N); flags SE entered as SD
  # md/mdw by default, d/g/dw/gw only via sd_outlier_smd (mixed instrument scales)
  sd_default <- c("md", "mdw")
  sd_optin   <- c("d", "g", "dw", "gw")
  sd_apply <- !is.null(measure) &&
    (measure %in% sd_default ||
     (isTRUE(opts$sd_outlier_smd) && measure %in% sd_optin))
  if (sd_apply) {
    K_sd <- if (!is.null(opts$sd_outlier_ratio)) opts$sd_outlier_ratio else 5
    # effective N: n_total when available, else n_exp (within-subject: n_nexp NA)
    n_eff <- if (!is.null(n_total)) n_total else rep(NA_real_, n)
    if (!is.null(n_exp)) {
      miss <- (is.na(n_eff) | !is.finite(n_eff)) & !is.na(n_exp) & is.finite(n_exp)
      n_eff[miss] <- n_exp[miss]
    }
    pooled_sd <- if (!is.null(sd_exp) && !is.null(sd_nexp)) {
      sqrt((sd_exp^2 + sd_nexp^2) / 2)
    } else rep(NA_real_, n)
    implied_sd <- if (!is.null(n_eff)) se * sqrt(n_eff) else rep(NA_real_, n)
    spread <- ifelse(!is.na(pooled_sd) & is.finite(pooled_sd) & pooled_sd > 0,
                     pooled_sd, implied_sd)
    sp_ok <- !is.na(spread) & is.finite(spread) & spread > 0
    if (sum(sp_ok) >= 4) {
      sp_med <- stats::median(spread[sp_ok])
      if (sp_med > 0) {
        for (i in seq_len(n)) {
          if (any(grepl("SE outlier", flags[[i]]))) next  # D2 already named this row
          if (sp_ok[i] && spread[i] < sp_med / K_sd) {
            msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""
            flags[[i]] <- c(flags[[i]],
              paste0("[UNUSUAL] SD outlier: spread = ", round(spread[i], 3),
                     " is ", round(sp_med / spread[i], 1),
                     "x below the pool median (", round(sp_med, 3),
                     ") - an SD this small for the scale can mean a standard error was entered as a standard deviation",
                     msuf))
          }
        }
      }
    }
  }

  return(flags)
}


#########

#' Flag rows when studies in the pool have conflicting significant directions
#'
#' Often a group labeling error rather than genuine heterogeneity; all rows on
#' both sides are flagged. A side qualifies by the count floor
#' (>= direction_conflict_min rows and >= direction_conflict_pct % of the pool)
#' or by a lone significant row that is >= direction_conflict_lone_pct % of the
#' pool.
#'
#' @param es numeric vector of effect sizes (length n)
#' @param ci_lo numeric vector of CI lower bounds (length n)
#' @param ci_up numeric vector of CI upper bounds (length n)
#' @param opts flag options list (uses enable_cross_row, direction_conflict_min,
#'   direction_conflict_pct, and direction_conflict_lone_pct)
#' @param info_used character vector of method names for traceability (optional)
#' @param measure character, effect size measure type (optional)
#' @param exp logical, whether values are on exponentiated scale
#' @return list of character vectors (one per row)
#' @noRd
.flag_cross_row_direction_conflict <- function(es, ci_lo, ci_up, opts,
                                                info_used = NULL,
                                                measure = NULL,
                                                exp = FALSE) {
  n <- length(es)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  if (!isTRUE(opts$enable_cross_row)) return(flags)

  # signed measures only (md/mdw included)
  signed_measures <- c("d", "g", "dw", "gw", "md", "mdw", "r", "z",
                       "logor", "logrr", "logirr", "loghr")
  if (is.null(measure) || !measure %in% signed_measures) return(flags)

  # Null reference: 1 on exponentiated scale for log-scale measures, 0 otherwise
  log_measures <- c("logor", "logrr", "logirr", "loghr")
  null_val <- if (isTRUE(exp) && measure %in% log_measures) 1 else 0

  es_label <- if (!is.null(measure) && nzchar(measure)) measure else "ES"
  if (isTRUE(exp) && measure %in% log_measures) es_label <- sub("^log", "", es_label)

  sig_pos <- which(!is.na(ci_lo) & is.finite(ci_lo) & ci_lo > null_val)
  sig_neg <- which(!is.na(ci_up) & is.finite(ci_up) & ci_up < null_val)

  n_pos <- length(sig_pos)
  n_neg <- length(sig_neg)

  min_threshold <- if (!is.null(opts$direction_conflict_min)) {
    opts$direction_conflict_min
  } else {
    2
  }

  pct_threshold <- if (!is.null(opts$direction_conflict_pct)) {
    opts$direction_conflict_pct
  } else {
    10
  }

  lone_pct_threshold <- if (!is.null(opts$direction_conflict_lone_pct)) {
    opts$direction_conflict_lone_pct
  } else {
    20
  }

  n_valid <- sum(!is.na(es) & is.finite(es))
  if (n_valid == 0) return(flags)

  pct_pos <- 100 * n_pos / n_valid
  pct_neg <- 100 * n_neg / n_valid

  side_qualifies <- function(n_side, pct_side) {
    (n_side >= min_threshold && pct_side >= pct_threshold) ||
      (n_side >= 1 && pct_side >= lone_pct_threshold)
  }
  if (!side_qualifies(n_pos, pct_pos) || !side_qualifies(n_neg, pct_neg)) {
    return(flags)
  }

  # minority side flagged [UNUSUAL], majority side [INFO]; tie: all [UNUSUAL]
  emit_info <- !isFALSE(opts$direction_conflict_info)
  majority <- if (n_pos > n_neg) "pos" else if (n_neg > n_pos) "neg" else "tie"

  build_msg <- function(i, dir) {
    msuf     <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""
    dir_word <- if (dir == "pos") "positive" else "negative"
    n_opp    <- if (dir == "pos") n_neg else n_pos
    core <- paste0("this row is significantly ", dir_word, " (",
                   es_label, " = ", round(es[i], 3),
                   ", CI [", round(ci_lo[i], 3), ", ", round(ci_up[i], 3), "]", msuf, ")")
    is_minority <- (majority == "tie") || (majority != dir)
    if (is_minority) {
      paste0("[UNUSUAL] Direction outlier: ", core,
             ", opposite to the pool majority (", n_opp,
             " significant row(s) favour the other direction) - verify the arm/group labels and the extracted value against the primary (likely extraction error or a genuine outlier)")
    } else {
      paste0("[INFO] Direction conflict in pool: ", core,
             ", agreeing with the majority direction while ", n_opp,
             " significant row(s) point the opposite way (see the UNUSUAL outlier flag)")
    }
  }

  for (dir in c("pos", "neg")) {
    rows <- if (dir == "pos") sig_pos else sig_neg
    is_minority_side <- (majority == "tie") || (majority != dir)
    if (!is_minority_side && !emit_info) next  # majority side suppressed
    for (i in rows) flags[[i]] <- c(flags[[i]], build_msg(i, dir))
  }

  return(flags)
}


#########

#' Flag rows that share a non-NA study_id with another row
#'
#' Same trial entered more than once (full publication + subanalysis, abstract
#' + full paper). Informational: data untouched, the reviewer aggregates or
#' drops.
#'
#' @param study_id character vector (or NULL) of per-row study IDs
#' @param opts flag options list (uses enable_cross_row)
#' @param info_used character vector of method names for traceability (optional)
#' @return list of character vectors (one per row)
#' @noRd
.flag_cross_row_duplicates <- function(study_id, opts, info_used = NULL,
                                       group_key = NULL) {
  n <- length(study_id)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  if (!isTRUE(opts$enable_cross_row)) return(flags)
  if (is.null(study_id) || n == 0)    return(flags)

  sid   <- as.character(study_id)
  valid <- !is.na(sid) & nzchar(sid)
  if (sum(valid) < 2) return(flags)

  # Scope duplication WITHIN group: a repeated study_id counts as a duplicate only
  # when the rows also share the same group. Multivariate / multi-outcome data
  # legitimately repeats a study_id across outcome groups, so keying on study_id
  # alone would flag it wholesale; keying on (group, study_id) does not. A NULL
  # group_key collapses to a single group -> the previous whole-dataset behaviour.
  grp <- if (is.null(group_key)) rep("__all__", n) else as.character(group_key)
  multi_grp <- length(unique(grp[valid])) > 1L
  dup_key <- paste(grp, sid, sep = "\r")

  tab      <- table(dup_key[valid])
  dup_keys <- names(tab[tab >= 2])
  if (length(dup_keys) == 0) return(flags)

  # no method suffix, the duplicate is a row-level fact
  for (dk in dup_keys) {
    rows   <- which(dup_key == dk & valid)
    dup_id <- sid[rows[1]]
    grp_txt <- if (multi_grp) paste0(" within group '", .group_label(grp[rows[1]]), "'") else ""
    for (i in rows) {
      others <- setdiff(rows, i)
      flags[[i]] <- c(flags[[i]],
        paste0("[INFO] Duplicate study_id '", dup_id, "'", grp_txt,
               ": row shares study_id with ",
               paste(.row_ref(others, sid), collapse = ", "),
               " - verify these are independent observations. ",
               "Use aggregate_df() if rows should be pooled, ",
               "or drop one row if they describe the same trial."))
    }
  }
  flags
}


#########

#' Check if SE is plausible given the total sample size
#'
#' For standardized metrics the SE has a known relationship with sample size,
#' which catches variance or SD entered as SE and model-mismatch errors
#' (per-unit-slope HR pooled with binary HRs). F1: SE below m / sqrt(N), with m
#' measure-specific. F2: SE above a measure-aware ceiling. F3 (loghr/logirr
#' only): SE below sqrt(1/n_exp + 1/n_nexp), the minimum of
#' sqrt(1/d_exp + 1/d_nexp) at d = n events; not applied to logor/logrr, whose
#' SE depends on event and non-event counts
#' (SE_logOR = sqrt(1/a + 1/b + 1/c + 1/d)) and approaches 0 as events near N,
#' so an N-only floor would false-alarm on routine 2x2 data. Skipped for
#' unstandardized measures (md, nnt, rd, alpha, icc, prop).
#'
#' @param es numeric vector of effect sizes
#' @param se numeric vector of standard errors
#' @param measure character, the effect size measure
#' @param n_total numeric vector of total sample sizes
#' @param info_used character vector of method names (optional, for traceability)
#' @param n_exp numeric vector of exposed-arm sample sizes (optional; required for F3)
#' @param n_nexp numeric vector of non-exposed-arm sample sizes (optional; required for F3)
#' @return list of character vectors (one per row)
#' @noRd
.flag_se_sample_size <- function(es, se, measure, n_total, info_used = NULL,
                                  n_exp = NULL, n_nexp = NULL) {
  n <- length(es)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  # Only apply to standardized metrics where SE scales predictably with N
  standardized <- c("d", "g", "dw", "gw", "r", "z", "logor", "logrr", "loghr", "logirr")
  if (!measure %in% standardized) return(flags)

  ratio_event_measures <- c("loghr", "logirr")

  for (i in seq_len(n)) {
    if (is.na(se[i]) || !is.finite(se[i]) || se[i] <= 0) next
    msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""

    # F1 / F2 require a usable n_total
    if (!is.na(n_total[i]) && is.finite(n_total[i]) && n_total[i] >= 4) {
      sqrt_n <- sqrt(n_total[i])

      # F1: implausibly small SE (floors: 0.5/sqrt(N) d/g, 0.3/sqrt(N) r/z, 0.1 ratios)
      se_floor_mult <- switch(measure,
        "d" = , "g" = , "dw" = , "gw" = 0.5,
        "r" = , "z" = 0.3,
        0.1  # logOR, logRR, logHR, logIRR: conservative (event-rate dependent)
      )
      se_floor <- se_floor_mult / sqrt_n
      # r: true SE ~ (1 - r^2)/sqrt(n - 1), scale the floor by (1 - r^2)
      if (measure == "r" && !is.na(es[i]) && is.finite(es[i]) && abs(es[i]) <= 1) {
        se_floor <- se_floor * (1 - es[i]^2)
      }
      if (se[i] < se_floor) {
        flags[[i]] <- c(flags[[i]],
          sprintf("[UNUSUAL] Implausibly small SE for N=%d: SE=%.4f (floor ~%.4f)",
                  as.integer(n_total[i]), se[i], se_floor))
      }

      # F2: implausibly large SE
      # ceiling: tight for d/g and r/z, permissive for event-rate dependent ratios
      se_ceiling <- switch(measure,
        "d" = , "g" = , "dw" = , "gw" = max(1.0, 8.0 / sqrt_n),
        "r" = , "z" = max(0.6, 4.0 / sqrt_n),
        max(2.0, 8.0 / sqrt_n)  # logor/logrr/loghr/logirr: event-rate dependent
      )
      if (se[i] > se_ceiling) {
        flags[[i]] <- c(flags[[i]],
          sprintf("[UNUSUAL] Implausibly large SE for N=%d: SE=%.3f (ceiling ~%.2f)%s",
                  as.integer(n_total[i]), se[i], se_ceiling, msuf))
      }
    }

    # F3: theoretical Cox/Wald floor from per-arm sizes (ratio-of-event-rate
    # measures only; requires both n_exp and n_nexp). Skip if F1 already fired
    # on this row to avoid duplicate "small SE" messages.
    if (measure %in% ratio_event_measures &&
        !is.null(n_exp) && !is.null(n_nexp) &&
        !is.na(n_exp[i]) && !is.na(n_nexp[i]) &&
        is.finite(n_exp[i]) && is.finite(n_nexp[i]) &&
        n_exp[i] > 0 && n_nexp[i] > 0 &&
        !any(grepl("Implausibly small SE", flags[[i]]))) {
      se_floor_f3 <- sqrt(1 / n_exp[i] + 1 / n_nexp[i])
      if (se[i] < se_floor_f3) {
        flags[[i]] <- c(flags[[i]],
          sprintf("[UNUSUAL] SE below theoretical minimum given group sizes: SE=%.4f, floor=%.4f (n_exp=%d, n_nexp=%d) - likely per-unit-slope model, wrong N, or wrong SE%s",
                  se[i], se_floor_f3,
                  as.integer(n_exp[i]), as.integer(n_nexp[i]), msuf))
      }
    }
  }

  return(flags)
}


#########

#' Check internal consistency across estimation methods
#' @param n_estimations numeric vector
#' @param dispersion numeric vector (SD of ES across methods)
#' @param overlap numeric vector (CI overlap between min/max)
#' @param diff_mm numeric vector (diff between min and max ES)
#' @param opts list of flag options
#' @param measure character, the effect size measure
#' @param min_info character vector of method names for min ES (optional)
#' @param max_info character vector of method names for max ES (optional)
#' @return list of character vectors (one per row)
#' @noRd
.flag_internal_consistency <- function(n_estimations, dispersion, overlap,
                                        diff_mm, opts, measure,
                                        min_info = NULL, max_info = NULL,
                                        min_es = NULL, max_es = NULL,
                                        exp = FALSE, es = NULL,
                                        r_defaulted = NULL) {
  n <- length(n_estimations)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  diff_max <- if (measure %in% c("d", "g", "dw", "gw", "md", "mdw")) {
    opts$diff_max_smd
  } else if (measure == "r") {
    opts$diff_max_r
  } else if (measure %in% c("logor", "logrr", "logirr", "loghr")) {
    opts$diff_max_logor
  } else {
    opts$diff_max_smd  # fallback
  }

  dispersion_max <- if (measure %in% c("d", "g", "dw", "gw", "md", "mdw")) {
    opts$dispersion_max_smd
  } else if (measure == "r") {
    opts$dispersion_max_r
  } else if (measure %in% c("logor", "logrr", "logirr", "loghr")) {
    opts$dispersion_max_logor
  } else {
    opts$dispersion_max_smd  # fallback
  }

  # recompute diff on the log scale when exp = TRUE; dispersion approximated
  # via delta method (MaxAD(log X) ~ MaxAD(X) / X)
  log_measures <- c("logor", "logrr", "logirr", "loghr")
  on_exp_scale <- isTRUE(exp) && measure %in% log_measures

  # all rows min ~ -max: direction convention, skip E flags
  all_reversed <- isTRUE(opts$ignore_sign_reversal) &&
    .rows_all_reversed(n_estimations, min_es, max_es, on_exp_scale)

  for (i in seq_len(n)) {
    if (is.na(n_estimations[i]) || n_estimations[i] < 2) next
    if (all_reversed) next

    # min/max method suffix for cross-method flags
    mmsuf <- .minmax_suffix(
      if (!is.null(min_info)) min_info[i] else NULL,
      if (!is.null(max_info)) max_info[i] else NULL
    )
    # r-sensitivity note for the SE/CI-driven checks (E2/E2b only; E1/E3 are
    # r-invariant point-estimate checks under the default Bonett pre-post SMD)
    rsuf <- if (!is.null(r_defaulted) && length(r_defaulted) >= i) {
      .r_sensitive_note(r_defaulted[i])
    } else {
      ""
    }

    # Working values on the threshold's (log) scale
    disp_work <- dispersion[i]
    diff_work <- diff_mm[i]
    scale_note <- ""
    if (on_exp_scale) {
      scale_note <- " (log scale)"
      disp_work <- if (!is.na(dispersion[i]) && is.finite(dispersion[i]) &&
                       !is.null(es) && !is.na(es[i]) && is.finite(es[i]) &&
                       es[i] > 0) {
        dispersion[i] / es[i]
      } else {
        NA_real_
      }
      diff_work <- if (!is.null(min_es) && !is.null(max_es) &&
                       !is.na(min_es[i]) && !is.na(max_es[i]) &&
                       is.finite(min_es[i]) && is.finite(max_es[i]) &&
                       min_es[i] > 0 && max_es[i] > 0) {
        log(max_es[i]) - log(min_es[i])
      } else {
        NA_real_
      }
    }

    # E1: High dispersion across ES estimates (measure-specific threshold)
    if (!is.na(disp_work) && is.finite(disp_work) &&
        disp_work > dispersion_max) {
      flags[[i]] <- c(flags[[i]],
        paste0("[DISCORDANT] High dispersion across estimates: max|ES - median| = ",
               round(disp_work, 3), scale_note,
               " (threshold: ", dispersion_max, ")"))
    }

    # E2: Zero overlap between min/max CIs
    if (!is.na(overlap[i]) && is.finite(overlap[i]) && overlap[i] == 0) {
      flags[[i]] <- c(flags[[i]],
        paste0("[DISCORDANT] Zero CI overlap between min and max estimates", mmsuf, rsuf))
    }

    # E2b: Low overlap between min/max CIs
    if (!is.na(overlap[i]) && is.finite(overlap[i]) &&
        overlap[i] > 0 && overlap[i] < opts$overlap_min) {
      flags[[i]] <- c(flags[[i]],
        paste0("[DISCORDANT] Low CI overlap between min/max estimates: ",
               round(overlap[i] * 100, 1), "% (threshold: ",
               opts$overlap_min * 100, "%)", mmsuf, rsuf))
    }

    # E3: Large difference between min and max ES
    if (!is.na(diff_work) && is.finite(diff_work) &&
        diff_work > diff_max) {
      flags[[i]] <- c(flags[[i]],
        paste0("[DISCORDANT] Large min-max difference: ",
               round(diff_work, 3), scale_note,
               " (threshold: ", diff_max, ")", mmsuf))
    }
  }
  return(flags)
}


#########

#' Determine if a V-flag is relevant to a given suffix ("_crude" or "_adjusted").
#'
#' Extracts the column name from the flag string: adjusted columns (ancova_*,
#' *_adj, cov_outcome_*, n_cov_*) go to _adjusted only, sample sizes to both,
#' everything else to _crude only.
#'
#' @param flag_str single V-flag string, e.g. "Asymmetric CI for 'md' ..."
#' @param suffix "_crude" or "_adjusted"
#' @return logical: TRUE if the flag should appear in this scope
#' @noRd
.v_flag_matches_scope <- function(flag_str, suffix) {
  # Extract column name from the quoted part, e.g. 'md'
  m <- regmatches(flag_str, regexpr("'[^']+'", flag_str))
  if (length(m) == 0) return(TRUE)  # Can't determine -> show everywhere
  col <- gsub("'", "", m)

  # Adjusted-only columns
  if (grepl("ancova|_adj$|^cov_outcome|^n_cov", col)) {
    return(suffix == "_adjusted")
  }

  # Shared columns (sample sizes): appear in both
  shared <- c("n_exp", "n_nexp", "n_sample", "n_cases", "n_controls",
              "n_cases_exp", "n_cases_nexp", "n_controls_exp", "n_controls_nexp")
  if (col %in% shared) return(TRUE)

  # Everything else is crude-only
  return(suffix == "_crude")
}

#########

#' Apply all quality/plausibility flags to the summary dataframe
#' @param res data.frame, the assembled summary result
#' @param measure character, the effect size measure
#' @param exp logical, whether the measure is on exp() scale
#' @param suffix character, "" or "_crude" or "_adjusted"
#' @param raw_data data.frame, the original input data
#' @param opts list of flag options (thresholds)
#' @param input_validation character vector (one per raw_data row), or NULL.
#'   Contains semicolon-separated validation flags from .validate_input_data().
#' @return data.frame with added flags column
#' @noRd
.flag_es_quality <- function(res, measure, exp, suffix, raw_data, opts,
                             input_validation = NULL,
                             alpha_to_es = "bonett",
                             icc_to_es = "bonett",
                             prop_to_es = "raw",
                             pre_post_to_smd = "bonett",
                             pool_sd = FALSE,
                             r_defaulted = NULL,
                             smd_denom = NULL,
                             es_order = NULL) {
  n <- nrow(res)
  flag_col <- paste0("flags", suffix)

  es <- suppressWarnings(as.numeric(as.character(res[[paste0("es", suffix)]])))
  se <- suppressWarnings(as.numeric(as.character(res[[paste0("se", suffix)]])))
  ci_lo <- suppressWarnings(as.numeric(as.character(res[[paste0("es_ci_lo", suffix)]])))
  ci_up <- suppressWarnings(as.numeric(as.character(res[[paste0("es_ci_up", suffix)]])))

  row_idx_map <- match(res$row_id, raw_data$row_id)
  # Per-row grouping key for the cross-row checks (D1/D2/D3, G, H, E4/E6/E7).
  # NULL flag_group -> a single "__all__" pool (unchanged whole-dataset behaviour).
  group_key <- .build_group_key(raw_data, opts$flag_group, idx = row_idx_map, n = n)

  # ---- Comparison representatives (route view) --------------------------------
  # Every CROSS-ROW check below assumes one row = one comparison: it compares a
  # study against its peers. Under main_es = FALSE a comparison occupies one row
  # per estimation route, so that invariant breaks and the checks silently change
  # their unit of analysis -- a study with k routes enters the IQR pool k times
  # (masking its own outlier status), counts k times toward the direction-conflict
  # floor, and matches its own study_id k times in the duplication check.
  #
  # Fix: decide every cross-row property on ONE representative row per
  # (row_id, scope), then broadcast the verdict back to that comparison's rows.
  # When main_es = TRUE every row is its own representative, cmp_map is the
  # identity and this is a no-op -- the default path is bit-for-bit unchanged.
  cmp_key <- paste(res$row_id,
                   if ("adjusted_input" %in% colnames(res)) res$adjusted_input else "",
                   sep = "\r")
  # The representative is the route the hierarchy would actually SELECT, not simply the
  # first row of the group. Under main_es = FALSE the rows of a comparison appear in the
  # internal list-construction order, which is not the hierarchy order, so taking the
  # first row picked an arbitrary route: every cross-row verdict was then computed from
  # an estimate the user never sees, and an outlier carried by the selected route could
  # go unreported (main_es = TRUE flagged it, main_es = FALSE did not).
  # es_order is the vector of info_used values ranked by the hierarchy, supplied by
  # summary(). Rows whose method is absent from it rank last, and when it is unavailable
  # the previous first-row behaviour is kept. Under main_es = TRUE each comparison has a
  # single row, so the choice is vacuous and the default path stays bit-for-bit identical.
  cmp_rank <- if (!is.null(es_order) && "info_used" %in% colnames(res)) {
    rk <- match(as.character(res$info_used), as.character(es_order))
    ifelse(is.na(rk), .Machine$integer.max, rk)
  } else {
    seq_len(n)
  }
  cmp_rep_i <- if (n > 0L) {
    sort(unname(vapply(split(seq_len(n), cmp_key),
                       function(ix) ix[which.min(cmp_rank[ix])], integer(1))))
  } else {
    integer(0)
  }
  cmp_map <- match(cmp_key, cmp_key[cmp_rep_i])     # every row -> its representative slot
  cmp_rep <- seq_len(n) %in% cmp_rep_i              # logical mask over all rows
  # Expand a length(cmp_rep_i) list of flag tokens back to all n rows.
  .cmp_broadcast <- function(sub) {
    out <- vector("list", n)
    for (i in seq_len(n)) out[[i]] <- sub[[cmp_map[i]]]
    out
  }
  # Per-row r_pre_post-default flag (aligned to raw_data rows, like input_validation)
  r_def_row <- if (!is.null(r_defaulted)) {
    as.logical(r_defaulted[row_idx_map])
  } else {
    rep(FALSE, n)
  }
  # per-row endpoint-SMD standardizer (Glass rows use control-df CIs -- A6)
  smd_denom_row <- if (!is.null(smd_denom)) {
    as.character(smd_denom[row_idx_map])
  } else {
    NULL
  }
  n_sample <- if ("n_sample" %in% colnames(raw_data)) {
    raw_data$n_sample[row_idx_map]
  } else {
    rep(NA_real_, n)
  }
  n_exp_vec <- if ("n_exp" %in% colnames(raw_data)) {
    raw_data$n_exp[row_idx_map]
  } else {
    rep(NA_real_, n)
  }
  n_nexp_vec <- if ("n_nexp" %in% colnames(raw_data)) {
    raw_data$n_nexp[row_idx_map]
  } else {
    rep(NA_real_, n)
  }
  # Total N: prefer n_sample, fall back to n_exp + n_nexp
  n_total <- ifelse(!is.na(n_sample), n_sample,
                    ifelse(!is.na(n_exp_vec) & !is.na(n_nexp_vec),
                           n_exp_vec + n_nexp_vec, NA_real_))

  n_est_col <- paste0("n_estimations", suffix)
  disp_col <- paste0("dispersion_es", suffix)
  overlap_col <- paste0("overlap_min_max", suffix)
  diff_col <- paste0("diff_min_max", suffix)

  n_estimations <- suppressWarnings(as.numeric(as.character(res[[n_est_col]])))
  dispersion <- suppressWarnings(as.numeric(as.character(res[[disp_col]])))
  overlap <- suppressWarnings(as.numeric(as.character(res[[overlap_col]])))
  diff_mm <- suppressWarnings(as.numeric(as.character(res[[diff_col]])))

  info_col <- paste0("info_used", suffix)
  min_info_col <- paste0("min_info", suffix)
  max_info_col <- paste0("max_info", suffix)
  info_used <- if (info_col %in% colnames(res)) as.character(res[[info_col]]) else NULL
  min_info_vec <- if (min_info_col %in% colnames(res)) as.character(res[[min_info_col]]) else NULL
  max_info_vec <- if (max_info_col %in% colnames(res)) as.character(res[[max_info_col]]) else NULL

  # baseline columns for the min NNT flag (B5)
  baseline_risk_vec <- if ("baseline_risk" %in% colnames(raw_data)) {
    raw_data$baseline_risk[match(res$row_id, raw_data$row_id)]
  } else {
    NULL
  }
  baseline_rate_vec <- if ("baseline_rate" %in% colnames(raw_data)) {
    raw_data$baseline_rate[match(res$row_id, raw_data$row_id)]
  } else {
    NULL
  }

  # covariate count, for the A6 CI-width check: .es_from_d() builds the d/g interval
  # of an ADJUSTED row on qt(.975, N - 2 - q), so A6 must expect the same df
  n_cov_row <- if ("n_cov_ancova" %in% colnames(raw_data)) {
    suppressWarnings(as.numeric(raw_data$n_cov_ancova[match(res$row_id, raw_data$row_id)]))
  } else {
    NULL
  }

  f_a <- .flag_numeric_integrity(es, se, ci_lo, ci_up, info_used,
                                  measure = measure, exp = exp,
                                  n_total = n_total, n_exp = n_exp_vec,
                                  n_nexp = n_nexp_vec,
                                  enable_informational = isTRUE(opts$enable_informational),
                                  r_defaulted = r_def_row,
                                  prop_to_es = prop_to_es,
                                  smd_denom = smd_denom_row,
                                  n_cov_ancova = n_cov_row)
  f_b <- .flag_bounds_violations(es, se, ci_lo, ci_up, measure, exp, info_used,
                                  baseline_risk = baseline_risk_vec,
                                  baseline_rate = baseline_rate_vec,
                                  alpha_to_es = alpha_to_es,
                                  icc_to_es = icc_to_es,
                                  prop_to_es = prop_to_es)
  f_c <- .flag_plausibility(es, se, measure, opts, n_sample, info_used,
                            n_exp = n_exp_vec, n_nexp = n_nexp_vec,
                            exp = exp, prop_to_es = prop_to_es)
  sd_exp_vec  <- if ("mean_sd_exp"  %in% colnames(raw_data)) {
    suppressWarnings(as.numeric(raw_data$mean_sd_exp[match(res$row_id, raw_data$row_id)]))
  } else NULL
  sd_nexp_vec <- if ("mean_sd_nexp" %in% colnames(raw_data)) {
    suppressWarnings(as.numeric(raw_data$mean_sd_nexp[match(res$row_id, raw_data$row_id)]))
  } else NULL
  # D1/D2/D3 and G run WITHIN each group (NULL flag_group -> one "__all__" pool).
  # NULL vector args index to NULL, which the checks already tolerate. Their messages
  # are self-contained, so .by_group can tag each with its group label.
  # Run on comparison representatives only (identity when main_es = TRUE), then
  # broadcast, so route replicates cannot inflate the pool or mask their own study.
  nR <- length(cmp_rep_i)
  esR <- es[cmp_rep_i]; seR <- se[cmp_rep_i]
  ci_loR <- ci_lo[cmp_rep_i]; ci_upR <- ci_up[cmp_rep_i]
  info_usedR <- if (is.null(info_used)) NULL else info_used[cmp_rep_i]
  n_totalR <- n_total[cmp_rep_i]
  n_expR <- n_exp_vec[cmp_rep_i]; n_nexpR <- n_nexp_vec[cmp_rep_i]
  sd_expR <- if (is.null(sd_exp_vec)) NULL else sd_exp_vec[cmp_rep_i]
  sd_nexpR <- if (is.null(sd_nexp_vec)) NULL else sd_nexp_vec[cmp_rep_i]
  group_keyR <- group_key[cmp_rep_i]

  f_d <- .cmp_broadcast(.by_group(group_keyR, nR, function(idx)
    .flag_cross_row_outliers(esR[idx], seR[idx], opts, info_usedR[idx], measure,
                             n_total = n_totalR[idx],
                             n_exp = n_expR[idx], n_nexp = n_nexpR[idx],
                             sd_exp = sd_expR[idx], sd_nexp = sd_nexpR[idx],
                             suffix = suffix, exp = exp)))
  f_g <- .cmp_broadcast(.by_group(group_keyR, nR, function(idx)
    .flag_cross_row_direction_conflict(esR[idx], ci_loR[idx], ci_upR[idx], opts,
                                       info_usedR[idx], measure, exp)))

  # H: Cross-row study duplication (reads study_id from raw_data;
  # always pass a length-n vector so f_dup[[i]] is safe downstream). Scoped WITHIN
  # group_key: sharing a study_id across different groups (e.g. one trial reporting
  # several outcomes) is expected in multivariate data and must NOT be flagged.
  study_id_vec <- if ("study_id" %in% colnames(raw_data)) {
    as.character(raw_data$study_id[match(res$row_id, raw_data$row_id)])
  } else {
    rep(NA_character_, n)
  }
  # Representatives only: otherwise a single comparison's k route rows all share
  # its study_id and the check reports the comparison as a duplicate of itself,
  # recommending aggregate_df() or dropping rows -- both destructive here.
  f_dup <- .cmp_broadcast(
    .flag_cross_row_duplicates(study_id_vec[cmp_rep_i], opts, info_usedR,
                               group_key = group_keyR))

  f_f <- .flag_se_sample_size(es, se, measure, n_total, info_used,
                              n_exp = n_exp_vec, n_nexp = n_nexp_vec)
  min_es_col <- paste0("min_es_value", suffix)
  max_es_col <- paste0("max_es_value", suffix)
  min_es_vec <- if (min_es_col %in% colnames(res)) {
    suppressWarnings(as.numeric(as.character(res[[min_es_col]])))
  } else NULL
  max_es_vec <- if (max_es_col %in% colnames(res)) {
    suppressWarnings(as.numeric(as.character(res[[max_es_col]])))
  } else NULL

  f_e <- .flag_internal_consistency(n_estimations, dispersion, overlap,
                                     diff_mm, opts, measure,
                                     min_info_vec, max_info_vec,
                                     min_es = min_es_vec, max_es = max_es_vec,
                                     exp = exp, es = es,
                                     r_defaulted = r_def_row)

  # E4: Cross-row NNT type mixing (risk-based vs rate-based), scoped WITHIN group_key
  # (only rows that would actually be pooled together can be "mixed").
  f_nnt_mix <- vector("list", n)
  for (i in seq_len(n)) f_nnt_mix[[i]] <- character(0)
  if (measure == "nnt" && !is.null(info_used)) {
    rate_methods <- "cases_time"
    multi_grp <- length(unique(group_key)) > 1L
    for (g in unique(group_key)) {
      gi <- which(group_key == g & !is.na(info_used) & nchar(info_used) > 0 & !is.na(es))
      # Whether the POOL mixes NNT types is decided on the routes that would
      # actually be selected; a non-selected cases_time route on one comparison
      # must not make every other comparison report "Mixed NNT types".
      gi_rep <- gi[cmp_rep[gi]]
      if (length(gi_rep) < 2) next
      vg <- info_used[gi_rep]
      if (any(vg %in% rate_methods) && any(!vg %in% rate_methods)) {
        gtxt <- if (multi_grp) paste0(" (within group '", .group_label(g), "')") else ""
        for (i in gi) {
          is_rate <- info_used[i] %in% rate_methods
          this_type <- if (is_rate) "rate-based (person-time)" else "risk-based"
          other_type <- if (is_rate) "risk-based" else "rate-based (person-time)"
          msuf <- .method_suffix(info_used[i])
          f_nnt_mix[[i]] <- paste0(
            "[DISCORDANT] Mixed NNT types: this row uses ", this_type,
            " NNT", msuf, ", but other rows use ", other_type,
            " NNT - these have different units and should not be pooled together", gtxt)
        }
      }
    }
  }

  # E6: Cross-row SMD standardizer mixing (change-SD vs raw-score-SD metric)
  #
  # [INFO], not [DISCORDANT]: DISCORDANT is reserved for a single row whose
  # several input sources disagree with one another. This is a cross-row
  # property of the pool -- every row may be individually correct.
  #
  # An SMD standardized by the CHANGE SD (morris_dz on pre/post or mean-change
  # data) is not on the same scale as an SMD standardized by a raw-score SD:
  # under equal pre/post SDs, sd_change = sd_raw * sqrt(2(1-r)), so the two
  # differ by a factor 1/sqrt(2(1-r)) (they coincide only at r = 0.5). The
  # Cochrane Handbook (v6, section 10.5.2) advises against combining
  # change-score SMDs with post-intervention SMDs, "because the SDs used in the
  # standardization reflect different things". The package's own docs say the
  # same of d_rm vs d_z.
  #
  # Rows are classified by the standardizer their info_used + pre_post_to_smd
  # imply. Endpoint (means_sd etc.), baseline-SD (bonett) and raw-metric d_rm
  # rows all live on a raw-score SD, so they are pooled into one "raw-score-SD"
  # class and mixing them raises NOTHING -- d_rm is precisely the transformation
  # onto that metric. Only morris_dz rows sit on the change-SD metric.
  f_std_mix <- vector("list", n)
  for (i in seq_len(n)) f_std_mix[[i]] <- character(0)
  smd_measures <- c("d", "g", "dw", "gw")
  # Per-row effective standardizer method: a per-row pre_post_to_smd column (which the
  # pipeline applies per row) overrides the dataset-wide scalar. Reading only the scalar
  # would miss a mixed column -- the very configuration most likely to mix standardizers.
  eff_method <- rep(as.character(pre_post_to_smd), n)
  if (!is.null(raw_data) && "pre_post_to_smd" %in% colnames(raw_data)) {
    col_method <- as.character(raw_data$pre_post_to_smd[row_idx_map])
    eff_method <- ifelse(is.na(col_method) | col_method == "", eff_method, col_method)
  }
  if (isTRUE(opts$enable_cross_row) && measure %in% smd_measures &&
      !is.null(info_used) && any(eff_method == "morris_dz", na.rm = TRUE)) {
    # methods whose standardizer is the change SD when pre_post_to_smd is dz
    change_metric_methods <- c(
      "means_sd_pre_post", "means_se_pre_post", "means_ci_pre_post",
      "mean_change_sd", "mean_change_se", "mean_change_ci", "mean_change_pval",
      "paired_t", "paired_t_pval", "paired_f", "paired_f_pval",
      "means_sd_pre_post_single_group", "means_se_pre_post_single_group",
      "means_ci_pre_post_single_group",
      "mean_change_sd_single_group", "mean_change_se_single_group",
      "mean_change_ci_single_group", "mean_change_pval_single_group",
      "paired_t_single_group"
    )
    # A row is on the change-SD metric iff it used a change-eligible method AND its
    # effective standardizer is morris_dz (a means_sd_pre_post row on bonett/d_rm is
    # on the raw-score metric, not the change-SD metric).
    row_is_change_all <- (info_used %in% change_metric_methods) & (eff_method == "morris_dz")
    valid_idx <- which(!is.na(info_used) & nchar(info_used) > 0 & !is.na(es))
    # Mixing is only meaningful among rows that would actually be pooled: decide it
    # WITHIN each group_key so a change-SD outcome and a raw-SD outcome in a
    # multivariate sheet are not cross-flagged.
    multi_grp <- length(unique(group_key)) > 1L
    for (g in unique(group_key[valid_idx])) {
      gsel <- valid_idx[group_key[valid_idx] == g]
      # Trigger on the routes that would actually be selected: under the route
      # view a single comparison legitimately offers both an endpoint and a
      # pre/post route, which is a within-comparison estimand difference (E1/E3),
      # not the cross-row standardizer mixing this check is about.
      ic <- row_is_change_all[gsel[cmp_rep[gsel]]]
      if (!(any(ic) && any(!ic))) next
      gtxt <- if (multi_grp) paste0(" (within group '", .group_label(g), "')") else ""
      for (i in gsel) {
        row_is_change <- row_is_change_all[i]
        this_metric <- if (row_is_change) "the change-SD metric (d_z)" else "a raw-score-SD metric"
        other_metric <- if (row_is_change) "a raw-score-SD metric" else "the change-SD metric (d_z)"
        msuf <- .method_suffix(info_used[i])
        f_std_mix[[i]] <- paste0(
          "[INFO] Mixed SMD standardizers: this row is on ", this_metric, msuf,
          ", but other rows are on ", other_metric,
          " - these differ by a factor 1/sqrt(2(1-r)) and should not be pooled ",
          "(Cochrane Handbook 10.5.2). Use pre_post_to_smd = 'morris_drm' (or ",
          "'bonett') to put pre/post rows on the raw-score-SD metric of the ",
          "endpoint rows", gtxt)
      }
    }
  }

  # E7: paired t/F rows coexisting with pooled-standardizer rows
  #
  # A paired t (or F) statistic identifies each arm's mean_change / sd_change
  # ratio but NOT the two arms' SD ratio, so the pooled standardizing SD is not
  # recoverable: these routes necessarily standardize each arm by its own SD and
  # subtract. When the user has opted into a pooled standardizer (pool_sd = TRUE)
  # for the rows that CAN be pooled, the two constructions differ whenever a
  # study's arm SDs differ. Informational: the paired-t rows are not wrong, they
  # are simply the best obtainable from the reported statistic. Silent under the
  # default (pool_sd = FALSE), where every row uses the per-arm construction.
  f_paired_t_mix <- vector("list", n)
  for (i in seq_len(n)) f_paired_t_mix[[i]] <- character(0)
  if (isTRUE(opts$enable_cross_row) && measure %in% smd_measures &&
      !is.null(info_used) && isTRUE(pool_sd)) {
    per_arm_only_methods <- c("paired_t", "paired_t_pval", "paired_f", "paired_f_pval")
    pooled_capable_methods <- c(
      "means_sd_pre_post", "means_se_pre_post", "means_ci_pre_post",
      "mean_change_sd", "mean_change_se", "mean_change_ci", "mean_change_pval"
    )
    valid_idx <- which(!is.na(info_used) & nchar(info_used) > 0 & !is.na(es))
    # P25: under pool_sd = TRUE the per-arm fallback contradicts the user's
    # explicit pooling request even when NO pooled-capable row coexists (a pool
    # made only of paired t/F rows used to be silent); disclose it either way,
    # with a message that fits each situation. Scoped WITHIN group_key so the
    # coexistence is decided among rows that would actually be pooled.
    multi_grp <- length(unique(group_key)) > 1L
    for (g in unique(group_key[valid_idx])) {
      gsel <- valid_idx[group_key[valid_idx] == g]
      # Coexistence decided on selected routes (see E4/E6 above).
      gsel_rep <- gsel[cmp_rep[gsel]]
      has_per_arm <- any(info_used[gsel_rep] %in% per_arm_only_methods)
      has_pooled  <- any(info_used[gsel_rep] %in% pooled_capable_methods)
      if (!has_per_arm) next
      gtxt <- if (multi_grp) paste0(" (within group '", .group_label(g), "')") else ""
      for (i in gsel) {
        if (!info_used[i] %in% per_arm_only_methods) next
        msuf <- .method_suffix(info_used[i])
        f_paired_t_mix[[i]] <- if (has_pooled) {
          paste0(
            "[INFO] Per-arm standardizer: a paired t/F statistic does not identify ",
            "the two arms' SD ratio", msuf, ", so this row standardizes each arm by ",
            "its own SD, while other rows in this pool use an SD pooled across arms ",
            "(you set pool_sd = TRUE). The two constructions coincide only when a ",
            "study's arm SDs are equal", gtxt)
        } else {
          paste0(
            "[INFO] Per-arm standardizer: a paired t/F statistic does not identify ",
            "the two arms' SD ratio", msuf, ", so this row standardizes each arm by ",
            "its own SD although you set pool_sd = TRUE. The pooled construction ",
            "is not recoverable from a paired t/F statistic", gtxt)
        }
      }
    }
  }

  # E5: min and max ES on opposite sides of the null (0, or 1 when exp = TRUE)
  # suppressed with E1-E3 when all rows show min ~ -max
  f_direction <- vector("list", n)
  for (i in seq_len(n)) f_direction[[i]] <- character(0)
  signed_measures <- c("d", "g", "dw", "gw", "md", "mdw", "r", "z",
                       "logor", "logrr", "logirr", "loghr")
  log_measures_e5 <- c("logor", "logrr", "logirr", "loghr")
  on_exp_scale_e5 <- isTRUE(exp) && measure %in% log_measures_e5
  null_ref <- if (on_exp_scale_e5) 1 else 0
  all_reversed <- isTRUE(opts$ignore_sign_reversal) &&
    .rows_all_reversed(n_estimations, min_es_vec, max_es_vec, on_exp_scale_e5)
  # near-null guard: both sides must sit >= direction_disagreement_min from the null
  dd_min <- if (!is.null(opts$direction_disagreement_min)) {
    opts$direction_disagreement_min
  } else {
    0.1
  }
  null_dist <- function(v) {
    if (on_exp_scale_e5) {
      if (v <= 0) Inf else abs(log(v))
    } else {
      abs(v)
    }
  }
  if (!all_reversed && measure %in% signed_measures) {
    min_es_col <- paste0("min_es_value", suffix)
    max_es_col <- paste0("max_es_value", suffix)
    if (min_es_col %in% colnames(res) && max_es_col %in% colnames(res)) {
      min_es <- suppressWarnings(as.numeric(as.character(res[[min_es_col]])))
      max_es <- suppressWarnings(as.numeric(as.character(res[[max_es_col]])))
      for (i in seq_len(n)) {
        if (!is.na(min_es[i]) && !is.na(max_es[i]) &&
            is.finite(min_es[i]) && is.finite(max_es[i]) &&
            min_es[i] < null_ref && max_es[i] > null_ref &&
            null_dist(min_es[i]) >= dd_min && null_dist(max_es[i]) >= dd_min) {
          mmsuf <- .minmax_suffix(
            if (!is.null(min_info_vec)) min_info_vec[i] else NULL,
            if (!is.null(max_info_vec)) max_info_vec[i] else NULL
          )
          f_direction[[i]] <- paste0(
            "[DISCORDANT] Direction disagreement: min ES = ", round(min_es[i], 3),
            ", max ES = ", round(max_es[i], 3),
            " - check group labeling", mmsuf)
        }
      }
    }
  }

  # Merge all flags per row (V-flags first, then A-F)
  res[[flag_col]] <- vapply(seq_len(n), function(i) {
    # Prepend input validation flags (V1-V4) if available
    v_flags <- character(0)
    if (!is.null(input_validation)) {
      row_idx <- match(res$row_id[i], raw_data$row_id)
      if (!is.na(row_idx) && nchar(input_validation[row_idx]) > 0) {
        v_flags <- strsplit(input_validation[row_idx], "; ")[[1]]
        # Filter V-flags by crude/adjusted scope when split
        if (suffix %in% c("_crude", "_adjusted")) {
          v_flags <- v_flags[vapply(v_flags, function(f) {
            .v_flag_matches_scope(f, suffix)
          }, logical(1))]
        }
      }
    }
    combined <- c(v_flags, f_a[[i]], f_b[[i]], f_c[[i]], f_d[[i]], f_e[[i]], f_f[[i]], f_g[[i]], f_dup[[i]], f_nnt_mix[[i]], f_std_mix[[i]], f_paired_t_mix[[i]], f_direction[[i]])
    if (length(combined) == 0) return("")
    paste(combined, collapse = "; ")
  }, character(1))

  return(res)
}
