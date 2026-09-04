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
    # One or more raw-data column names to scope the cross-row checks by. It now
    # reaches BOTH tiers: the post-computation ones (D1/D2/D3, G, H, E4/E6/E7/E8) via
    # .by_group() in summary(), and the input-data ones (V35, V36, V37, V38, V39, V42,
    # V43) via the matching .build_group_key() call inside .validate_input_data().
    # V23 is the single exception and stays dataset-wide on purpose: it hunts for a
    # statistics block cloned from one paper into another, so confining that search to
    # one outcome or subgroup would hide the case it exists for.
    # NULL (the default) treats the whole dataset as one pool.
    # Set it (for example "outcome", or c("outcome", "subgroup")) for multivariate or
    # multi-outcome extraction sheets, where deviations are only meaningful within a
    # group of comparable rows. Per-row and cross-method checks are unaffected.
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
    # D2c: opt-in. Flags a row that carries an SE but no sample size, whose RAW SE is
    # more than se_ratio_extreme-fold tighter than the pool median. Off by default: it
    # is the only check comparing un-normalised SEs, and is sound only because it is
    # confined to rows the normalised checks already dropped.
    se_outlier_missing_n = FALSE, # D2c
    se_missing_n_ratio = 5, # D2c fold threshold, on the RAW SE scale
    sd_ratio_max = 10, # V12
    sd_ratio_bl_ep_min = 0.7, # V18
    baseline_imbalance_max = 0.30, # V21
    templated_min_match = 3, # V23
    paired_as_indep_tol = 0.07, # V26
    sd_outlier_ratio = 5, # V27/D3
    sd_outlier_smd = FALSE, # V27/D3, off: smd pools mix instrument scales
    # V29: the DISTRIBUTIONAL band on (max - min)/SD, in multiples of xi(n), emitted
    # as [UNUSUAL]. The ceiling is scaled with n inside the check (see V29), because
    # no FLAT multiplier is simultaneously sensitive below n = 30 and specific above
    # n = 500. These two multipliers do not govern the whole check: V29 tests an
    # assumption-free arithmetic pair FIRST and emits [INVALID] on it, and only a row
    # that clears those bounds is compared against this band.
    range_sd_lo_mult = 0.5, # V29
    range_sd_hi_mult = 2.5, # V29
    margin_range_min = 0.30, # V35
    # V44: floor_position bands. Measured in simulations/studies/11_reliability_se.R
    # (11e, 20 response distributions x C in {2,3,5,7}): the variance error c is
    # 1.56-4.28 below 0.20 and 0.92-1.03 above 0.28, and the bands hold across every
    # category count. Spearman rho(floor_position, c) = -0.96 against -0.16 for the
    # category count alone, so this is the gate and the category count is not.
    floor_position_severe = 0.20, # V44
    floor_position_mild = 0.28 # V44
  )
}


#' Warn about flag_options entries that will not do what the caller expects.
#'
#' Two failures pass unnoticed, and both leave the user believing they tuned
#' something:
#'
#'  1. A misspelled name. `flag_options` is a plain list merged by name, so
#'     `alpha_maxx = 0.5` is accepted everywhere and simply never read.
#'  2. A Tier-1 option passed to `summary()`. The input-data checks run inside
#'     `convert_df()`, before `summary()` is ever called, so an option only they
#'     consume arrives too late. Ten options are Tier-1 only and have no effect at
#'     all this way. `enable_cross_row`, `enable_informational` and `flag_group` are
#'     read by both tiers, so passing them to `summary()` half-works: it retunes the
#'     post-computation checks while the input-data ones have already run on the old
#'     value. That partial effect is the more confusing of the two, so it is named
#'     separately.
#'
#' Warns rather than errors, matching the package convention for user input: one
#' mistyped option must not abort an analysis.
#'
#' @param flag_options the user's list
#' @param context "convert_df" or "summary"
#' @noRd
.validate_flag_options <- function(flag_options, context = "summary") {
  if (is.null(flag_options) || length(flag_options) == 0) return(invisible(NULL))
  nms <- names(flag_options)
  if (is.null(nms)) nms <- rep("", length(flag_options))

  known <- names(.default_flag_options())
  unknown <- setdiff(nms[nzchar(nms)], known)
  if (length(unknown)) {
    warning(paste0(
      "Unrecognised flag_options name(s): '", paste(unknown, collapse = "', '"),
      "'. They are ignored. See ?convert_df for the available options."), call. = FALSE)
  }

  if (identical(context, "summary")) {
    # read by convert_df() only -- passing them to summary() does nothing
    tier1_only <- c("sd_ratio_max", "sd_ratio_bl_ep_min", "baseline_imbalance_max",
                    "templated_min_match", "paired_as_indep_tol", "margin_range_min",
                    # V44 and V29 are input-data checks too, and were missing here, so
                    # passing their thresholds to summary() failed silently.
                    "floor_position_severe", "floor_position_mild",
                    "range_sd_lo_mult", "range_sd_hi_mult")
    # read by BOTH tiers -- passing them to summary() retunes only the second
    # flag_group joined this list when convert_df() began passing it to
    # .validate_input_data(): it now scopes the Tier-1 cross-row checks (V35-V43) as
    # well as the Tier-2 ones, so handing it to summary() alone silently groups half of
    # them. That is exactly the partial effect this warning exists to name.
    both_tiers <- c("enable_cross_row", "enable_informational", "flag_group")

    late <- intersect(nms, tier1_only)
    if (length(late)) {
      warning(paste0(
        "flag_options passed to summary() cannot reach the input-data checks, which ",
        "already ran inside convert_df(): '", paste(late, collapse = "', '"),
        "' had no effect. Pass them to convert_df() instead."), call. = FALSE)
    }
    partial <- intersect(nms, both_tiers)
    if (length(partial)) {
      warning(paste0(
        "'", paste(partial, collapse = "', '"), "' is read by both the input-data and ",
        "the post-computation checks. Passing it to summary() retunes only the latter; ",
        "the input-data checks already ran inside convert_df(). Pass it to convert_df() ",
        "to change both."), call. = FALSE)
    }
  }
  # Names were checked above; the VALUES were not. A character threshold made every
  # numeric comparison against it NA, silently disabling the check the user was trying
  # to tighten, and NA/NULL/length > 1 aborted the run inside a comparison with an
  # error naming neither the option nor the value. Warn and fall back to the default,
  # which is what the rest of this validator does for a name it does not recognise.
  defaults <- .default_flag_options()
  bad <- character(0)
  nulls <- character(0)
  for (nm in intersect(nms[nzchar(nms)], known)) {
    v <- flag_options[[nm]]
    dflt <- defaults[[nm]]
    ok <- if (is.null(v)) {
      # An explicit NULL reads as "use the default", but modifyList() DELETES the key
      # rather than leaving the default in place, so the option arrives as NULL and the
      # first comparison against it returns logical(0) -- `if (logical(0))` is an error
      # with no mention of the option. Dropped here so the default survives. Silent:
      # the user asked for the default and gets it.
      nulls <- c(nulls, nm)
      TRUE
    } else if (is.logical(dflt)) {
      is.logical(v) && length(v) == 1L && !is.na(v)
    } else if (is.numeric(dflt)) {
      is.numeric(v) && length(v) == 1L && !is.na(v) && is.finite(v)
    } else if (identical(nm, "flag_group")) {
      # flag_group names one OR MORE grouping columns: c("outcome", "subgroup") is a
      # documented, tested use. Only the element type and NA-freeness are constrained.
      is.null(v) || (is.character(v) && length(v) >= 1L && !anyNA(v))
    } else if (is.character(dflt) || is.null(dflt)) {
      is.null(v) || (is.character(v) && length(v) == 1L && !is.na(v))
    } else TRUE
    if (!ok) bad <- c(bad, nm)
  }
  if (length(bad)) {
    warning(paste0(
      "flag_options value(s) of the wrong type or length: '", paste(bad, collapse = "', '"),
      "'. Each threshold must be a single finite number, each switch a single TRUE/FALSE. ",
      "The default is used instead."), call. = FALSE)
  }

  invisible(c(bad, nulls))
}


#' Build a per-row grouping key for the cross-row checks from opts$flag_group.
#'
#' Returns a length-n character vector used to scope the cross-row checks of BOTH
#' tiers: D1/D2/D3, G, H and E4/E6/E7/E8 in summary(), and V35, V36, V37, V38, V39,
#' V42 and V43 in .validate_input_data(), which calls this helper on the raw input
#' data and splits on the result exactly as .by_group() does. Every one of those
#' messages asserts something about "the rest of the pool", so leaving a Tier-1 check
#' unscoped would state a pool composition the user has already told the package is
#' wrong -- and for the majority-guarded V37/V42 the grouping IS the remedy their own
#' comments describe.
#' V23 is the one cross-row check deliberately left unscoped: it hunts for a
#' statistics block cloned from one paper into another, so restricting the search to
#' one outcome or subgroup would hide the case it exists for.
#' When flag_group is NULL, absent, or names no present column, every row gets the
#' same key ("__all__"), so the whole dataset forms one pool. Rows with an NA in a
#' grouping column fall into a shared "<NA>" level.
#'
#' Multi-column keys are joined with "\r" rather than with the " / " shown to the
#' user, because a grouping value can itself contain " / ". Joining on it would let
#' distinct rows collide into one pool: `(g1 = "x / y", g2 = "z")` and
#' `(g1 = "x", g2 = "y / z")` both produce "x / y / z", merging their D1/D2/D3, G, H,
#' E4 and E6/E7/E8 comparison sets. "\r" cannot occur in a value read from a
#' spreadsheet cell, so it is collision-proof, and it is the separator `cmp_key` and
#' `dup_key` already use. .group_label() renders it back as " / " for display, and
#' every message that shows a key routes through it.
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
  do.call(paste, c(parts, sep = "\r"))
}


#' Split rows by group_key, run a cross-row flag function within each group, and
#' reassemble a length-n list of flag vectors (global row order preserved).
#'
#' `fn` receives a vector of global row indices for one group and must return a
#' list of character vectors of the same length, in that index order. When more than
#' one group is present, each emitted flag gets a " (within group '<g>')" suffix so
#' the reviewer knows the comparison set; the suffix contains no "; ", so the
#' downstream flag merge is unaffected. Messages must be self-contained, with no
#' absolute row-number references, which holds for the D-family and G checks.
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


#' Make a group key or study_id safe to interpolate into a flag message.
#'
#' Two substitutions, in this order:
#'
#' 1. "\r" -> " / ". .build_group_key() joins multi-column keys with "\r" so that
#'    a value containing " / " cannot make two distinct rows share a pool; this
#'    renders the key back into the form the user is meant to read.
#' 2. "; " -> ", ". The merged flag string is split on "; ", so a value containing
#'    that sequence (an outcome named "HAM-D; total score", a study_id like
#'    "Huang; 2017") would create an untagged token downstream -- and for the
#'    re-split Tier-1 messages it does worse than that, truncating the message and
#'    mis-routing the fragment into both the crude and adjusted scopes, because a
#'    fragment carrying no quoted column name matches every scope.
#'
#' Applied to grouping values at every display site and to study_id inside
#' .row_ref(), which is what V23, V36 and Category H interpolate.
#'
#' @param g character vector, a raw group key or study_id
#' @noRd
.group_label <- function(g) {
  g <- gsub("\r", " / ", as.character(g), fixed = TRUE)
  gsub("; ", ", ", g, fixed = TRUE)
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
  # digits is pinned rather than left to getOption("digits"): the default 7 makes the
  # decimal count depend on the MAGNITUDE of the value (1234.5678 reads as 3 dp) and on
  # a global display option the caller never set for this purpose, so the same study
  # could be flagged or not according to the console's print settings. 15 is inside
  # binary64's ~15.95 significant decimal digits, so it does not resurrect
  # representation noise: .count_decimals(0.1 + 0.2) stays 1.
  s <- format(x, scientific = FALSE, drop0trailing = TRUE, digits = 15)
  if (!grepl("\\.", s)) return(1L)  # integers assumed at least 1dp precision
  dp <- nchar(sub(".*\\.", "", s))
  min(dp, 10L)
}


#' Reporting precision of a (value, lower, upper) triplet
#'
#' .count_decimals() reads the decimals R actually stores, and R does not store
#' trailing zeros: a bound printed as 0.90 arrives as 0.9 and is read as 1 dp. Pricing
#' each member's rounding pad from its own decimals therefore makes V2/V3 behave
#' differently depending on whether a reported digit happened to be a zero. On a
#' [0, 1] reliability column the pad becomes +/- 0.05, wider than most real CI
#' half-widths, so an arithmetically impossible interval survives validation while the
#' identical error one digit over is caught and the row removed.
#'
#' A study printing 0.85 and 0.86 is printing to 2 dp, so its 0.90 is 0.90 and not
#' 0.9. Take the precision of the most precisely printed member of the triplet as the
#' row's reporting precision, and build every pad from that one number.
#'
#' @param val,lo,up numeric scalars (the point estimate and its two bounds)
#' @return integer, the number of decimals the rounding error is priced at
#' @noRd
.triplet_precision <- function(val, lo, up) {
  max(.count_decimals(val), .count_decimals(lo), .count_decimals(up))
}


#' Input columns whose presence means a row feeds an r_pre_post-consuming route
#'
#' A row is "r-consuming" when it carries pre/post, mean-change or paired data, that
#' is, when it reaches one of the 19 exported routes that take \code{r_pre_post_exp}
#' or \code{r_pre_post_nexp}. Those routes impute \code{r_pre_post} when the user
#' leaves it blank, which is a substantive assumption: under the change-SD
#' standardizers (\code{morris_dz}, \code{cooper}/\code{morris_drm}) the assumed
#' correlation scales the point estimate, not merely the SE.
#'
#' This is the single source of truth for that test. \code{convert_df()} consumes it
#' twice -- for the \code{verbose} note and for flag V6, which also sets the
#' \code{r_defaulted} attribute driving the \code{(r-sensitive: ...)} annotations on
#' A6/E2/E2b in \code{summary()} -- through the arm-aware \code{.r_consuming_arms()}
#' closure defined in \code{convert_df()}, which reports r-consuming data per ARM so
#' that each mask can be ANDed with its own arm's defaulted-r flag. (A plain row-level
#' OR reported an imputation on single-group rows, and on two-group rows whose
#' pre/post data sits in one arm only, that had actually supplied their correlation.)
#'
#' \code{.r_consuming_arms()} is evaluated separately at each of those two sites
#' rather than computed once and reused: \code{x} is reassigned by
#' \code{.validate_input_data()} in between, and under \code{correct_inputs = TRUE}
#' that step may set an invalid pre/post value to NA. A row whose only pre/post datum
#' was just invalidated no longer feeds an r-consuming route, so V6 must see the
#' post-validation data while the verbose note sees the input as the user supplied it.
#'
#' Keeping one list matters because \code{intersect(., colnames(x))} drops any name
#' that does not exist, so a list naming a column that was never created fails quietly
#' rather than erroring.
#'
#' Single-group routes: there are no \code{*_single_group} columns. The single-group
#' entry points reuse the \code{_exp} columns (see the \code{es_from_*_single_group}
#' calls in \code{convert_df()}), so they are already covered here.
#'
#' Endpoint columns: \code{mean_exp} and \code{mean_sd_exp} are excluded on purpose.
#' They are populated on ordinary two-group rows that never touch an r-consuming
#' route, so keying on them would make the note and V6 fire on data with no pre/post
#' component at all. The discriminating signal is the baseline (\code{mean_pre_*}),
#' change (\code{mean_change_*}) and paired (\code{paired_*}) columns.
#'
#' @return character vector of column names
#' @noRd
.r_consuming_columns <- function() {
  c(
    # Baseline arm statistics -- es_from_means_{sd,se,ci}_pre_post and their
    # single-group twins
    "mean_pre_exp", "mean_pre_nexp",
    "mean_pre_sd_exp", "mean_pre_sd_nexp",
    "mean_pre_se_exp", "mean_pre_se_nexp",
    "mean_pre_ci_lo_exp", "mean_pre_ci_up_exp",
    "mean_pre_ci_lo_nexp", "mean_pre_ci_up_nexp",
    # Change scores -- es_from_mean_change_{sd,se,ci,pval} and single-group twins
    "mean_change_exp", "mean_change_nexp",
    "mean_change_sd_exp", "mean_change_sd_nexp",
    "mean_change_se_exp", "mean_change_se_nexp",
    "mean_change_ci_lo_exp", "mean_change_ci_up_exp",
    "mean_change_ci_lo_nexp", "mean_change_ci_up_nexp",
    "mean_change_pval_exp", "mean_change_pval_nexp",
    # Paired test statistics -- es_from_paired_{t,f}, their p-value variants and
    # es_from_paired_t_single_group
    "paired_t_exp", "paired_t_nexp",
    "paired_t_pval_exp", "paired_t_pval_nexp",
    "paired_f_exp", "paired_f_nexp",
    "paired_f_pval_exp", "paired_f_pval_nexp"
  )
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
    # Regression standardisers. A negative sd_dv sign-flips d/g through
    # es_from_beta_std() while d_se stays positive, which is the exact failure class
    # the input guards exist to close, and neither column was checked here.
    "sd_dv", "sd_iv",
    # SEs
    "mean_se_exp", "mean_se_nexp", "mean_pre_se_exp", "mean_pre_se_nexp",
    "mean_change_se_exp", "mean_change_se_nexp",
    "ancova_mean_se_exp", "ancova_mean_se_nexp",
    "md_se", "ancova_md_se", "logor_se", "logrr_se", "logirr_se",
    "rd_se", "omega_se", "icc_se", "cronbach_alpha_se", "linreg_b_se",
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
    # omega is the one measure whose standard error comes entirely from the reported
    # interval, since there is no (n, k) variance for it, so an unvalidated CI is an
    # unvalidated weight. Without this entry a transposed interval, a point estimate
    # outside its own interval, or bounds outside [0, 1] would each produce a
    # plausible omega with a wrong SE and an empty flag string. .ci_lower() and
    # .ci_upper() absorb the transposition, so this is the last chance to notice.
    #
    # keep_val: what V3 ("value outside its own CI") destroys. Everywhere else the
    # interval is the row's only dispersion source, so a value that contradicts it
    # leaves nothing usable and the whole triplet goes. The three reliability
    # coefficients are the exception: each has an interval-FREE standard error (the
    # (n, k) closed form for alpha and icc; se = NA, row visible but unpooled, for
    # omega), so refusing the interval alone keeps a study whose weight then comes
    # from n and k, which no mis-transcription of the interval can corrupt. V2 is not
    # given this treatment: an inverted interval cannot even be read bound-for-bound,
    # and a transposition that severe usually means the whole block was taken off the
    # wrong row, so the coefficient is not above suspicion either.
    list(val = "omega",     lo = "omega_ci_lo",     up = "omega_ci_up",     scale = "additive", keep_val = TRUE),
    list(val = "icc",       lo = "icc_ci_lo",       up = "icc_ci_up",       scale = "additive", keep_val = TRUE),
    list(val = "cronbach_alpha", lo = "cronbach_alpha_ci_lo", up = "cronbach_alpha_ci_up", scale = "additive", keep_val = TRUE),
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
                                  alpha_to_es = "bonett",
                                  icc_to_es = "bonett",
                                  omega_to_es = "bonett",
                                  templated_min_match = 3,
                                  measure = NULL,
                                  paired_as_indep_tol = 0.07,
                                  floor_position_severe = 0.20,
                                  floor_position_mild = 0.28,
                                  margin_range_min = 0.30,
                                  alpha_se_source = "closed_form",
                                  omega_se_source = "reported",
                                  range_sd_lo_mult = 0.5,
                                  range_sd_hi_mult = 2.5,
                                  flag_group = NULL,
                                  correct_inputs = TRUE) {
  n <- nrow(x)
  row_issues <- vector("list", n)
  for (i in seq_len(n)) row_issues[[i]] <- character(0)

  # Comparison sets for the Tier-1 CROSS-ROW checks (V35, V36, V37, V38, V39, V42,
  # V43), scoped by flag_group exactly as the Tier-2 ones are by .by_group(). Every
  # one of those messages asserts a fact about "the rest of the pool", so computing
  # them over the whole dataset states a pool composition the user has explicitly told
  # the package is wrong -- and for the majority-guarded checks (V37, V42) grouping is
  # the very remedy their own comments describe: two strata each with a constant k
  # have no minority to flag, while their union does.
  #
  # V23 stays deliberately dataset-wide: it looks for a block of numbers cloned from
  # one paper into another, and confining that search to one outcome or subgroup would
  # hide the case it exists for.
  v_group_key <- .build_group_key(x, flag_group, n = n)
  v_groups <- split(seq_len(n), v_group_key)
  v_multi <- length(v_groups) > 1L
  # The same suffix .by_group() appends, for the same reason. It contains no "; ", so
  # the downstream flag merge is unaffected, and it is appended AFTER the message's own
  # leading quoted column name, so .v_flag_matches_scope() still routes on that.
  v_grp_note <- function(gk) {
    if (v_multi) paste0(" (within group '", .group_label(names(v_groups)[gk]), "')") else ""
  }

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
                    # omega_se is the ONLY route by which an omega row gets a weight,
                    # so a zero there silently removes the study. The alpha route
                    # surfaces the equivalent through the Tier-2 "SE is zero" check,
                    # which omega cannot reach: .positive_or_na() has already NA'd the
                    # value before any effect size exists to check.
                    "omega_se",
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
    # omega = (sum lambda)^2 / ((sum lambda)^2 + sum theta) is a square divided by
    # itself plus a sum of variances, so unlike alpha it is structurally non-negative.
    # No measurement pattern produces a negative omega, and a 1.34M-run simulation
    # found no negative omega_total against 107 negative alphas. Its lower bound is
    # therefore 0, rather than the -Inf alpha gets, since a negative average inter-item
    # covariance really can drive alpha below zero (see V25). The upper bound stays 1:
    # an estimate above it is a Heywood case, which es_from_omega() preserves on the
    # raw scale.
    list(col = "omega",           lo = 0,    up = 1, label = "omega"),
    # icc bounded below by -1/(k-1) >= -1
    list(col = "icc",             lo = -1, up = 1, label = "ICC"),
    list(col = "prop",            lo = 0, up = 1, label = "proportion"),
    # Arm-level proportions consumed by es_from_2x2_prop()
    list(col = "prop_cases_exp",  lo = 0, up = 1, label = "proportion"),
    list(col = "prop_cases_nexp", lo = 0, up = 1, label = "proportion"),
    # Reliability CI bounds. The point estimate is bounded here; the bounds it is read
    # between were not, so a mistyped bound reached the SE cascade unchecked.
    list(col = "icc_ci_lo",       lo = -1, up = 1, label = "ICC CI bound"),
    list(col = "icc_ci_up",       lo = -1, up = 1, label = "ICC CI bound"),
    # Pre-post correlations: a correlation, and consumed as one by 19 routes
    list(col = "r_pre_post_exp",  lo = -1, up = 1, label = "correlation"),
    list(col = "r_pre_post_nexp", lo = -1, up = 1, label = "correlation"),
    # small_margin_prop is the SMALLER of the two 2x2 margin proportions, so its domain
    # is (0, 0.5]: at 0 the Bonett correction divides by zero and above 0.5 it is not
    # the smaller margin at all, which silently flips the sign of the correlation and
    # returns a negative standard error. open_lo because this list is otherwise
    # inclusive and 0 is not admissible here.
    list(col = "small_margin_prop", lo = 0, up = 0.5, open_lo = TRUE,
         label = "small margin proportion"),
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
    lo_bad <- if (isTRUE(bc$open_lo)) vals <= bc$lo else vals < bc$lo
    bad <- which(!is.na(vals) & (lo_bad | vals > bc$up))
    if (length(bad) > 0) {
      # Describe the valid range; one-sided when a bound is infinite
      range_txt <- if (is.infinite(bc$lo)) {
        sprintf("<= %g", bc$up)
      } else if (is.infinite(bc$up)) {
        sprintf(">= %g", bc$lo)
      } else if (isTRUE(bc$open_lo)) {
        sprintf("(%g, %g]", bc$lo, bc$up)
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

  # V32: eta-squared must lie in [0, 1). A value >= 1 makes the Cohen's d conversion
  # diverge, whether through d = 2*sqrt(eta2 / (1 - eta2)) for a raw eta-squared or
  # through the implied ANCOVA F = eta2 * df / (1 - eta2) for an adjusted one, so such
  # a row is set to missing (correct_inputs = TRUE) or flagged (FALSE). Negatives are
  # already caught by the non-negative-column check above. A large but valid
  # eta-squared that implies a huge SMD is surfaced downstream by the
  # post-computation Large-SMD check, so no separate plausibility tier is added here.
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

  # ICC facts every ICC-scoped check below needs, computed once and above the first
  # of them. icc_measure_ok and icc_type_eff used to be local to V31, which sits
  # below V40 and V41 -- and both of those need them: V40 to know whether the route
  # actually drops a 1.00, V41 to know that its single-measures floor does not apply
  # to an average-measures value.
  #
  # icc_measure_ok scopes the ICC-named checks to an ICC run. Tier-1 checks key on
  # which columns exist rather than on the requested measure, so a COSMIN-style sheet
  # holding alpha, omega and ICC side by side (one sheet, several properties, one pool
  # per property) would otherwise carry ICC pooling advice into every run. A NULL or
  # empty measure leaves direct callers unaffected.
  icc_measure_ok <- is.null(measure) || !nzchar(measure) || identical(measure, "icc")
  # .normalise_icc_type() rather than bare string equality: es_from_icc() resolves ten
  # spellings to "agreement" ("Agreement", "ICC(2,1)", "absolute agreement", "two-way
  # random", ...), and comparing the raw cell against the literal "agreement" matches
  # only two of them. convert_df() normalises the column as well, but this call is kept
  # so that direct callers of .validate_input_data() get the same answer from one
  # definition across three call sites. warn = FALSE because es_from_icc() owns the
  # user-facing warning.
  icc_type_eff <- if ("icc_type" %in% colnames(x)) {
    tt <- as.character(x[["icc_type"]])
    tt[is.na(tt)] <- "agreement"
    .normalise_icc_type(tt, warn = FALSE)
  } else {
    rep("agreement", nrow(x))
  }
  icc_k_col <- if ("n_measurements" %in% colnames(x)) {
    suppressWarnings(as.numeric(x[["n_measurements"]]))
  } else {
    rep(NA_real_, nrow(x))
  }
  icc_is_avg_col <- if ("icc" %in% colnames(x)) {
    !is.na(x[["icc"]]) & .icc_is_average(icc_type_eff)
  } else {
    rep(FALSE, nrow(x))
  }

  # V40: a reliability coefficient of exactly 1 has no transform, so the row is
  # dropped, and without this flag it is dropped without a message. V11 bounds alpha
  # at (-Inf, 1] and omega/icc at [0, 1] inclusive, so 1.00 passes validation; the
  # route then NAs it against a strict `< 1` gate, and es_guidance reports "No partial
  # input data found" on a fully populated row. Compare alpha = 1.02, one rounding
  # step away, which does get explained. alpha = 1.00 is a routine printed value, from
  # a short subscale or from a paper rounding .996.
  #
  # Warn-only and ROUTE-aware. Whether 1.00 survives is a property of the route, not
  # of the transform, and the three routes disagree:
  #
  #   es_from_cronbach_alpha()  keeps 1 on "raw" ONLY. Its se_defined gate refuses the
  #                             standard error on every non-raw scale at alpha = 1, so
  #                             hakstian_whalen returns es = 1 with se = NA and the row
  #                             leaves the pool just as the bonett one does -- it simply
  #                             leaves with an estimate attached.
  #                                                 (R/es_from_ALPHA.R, valid_es/se_defined)
  #   es_from_omega()           keeps 1 on "raw" only
  #                                                 (R/es_from_OMEGA.R, transformable)
  #   es_from_icc()             drops it on EVERY scale, its gate icc < 1 being
  #                             scale-independent   (R/es_from_ICC.R, valid_es)
  #
  # What a row needs in order to be POOLED is an estimate AND a standard error, so the
  # table is keyed on that, not on whether the transform exists. The
  # "log(0) is undefined" parenthetical is true only of the Bonett transform:
  # 1 - (1 - a)^(1/3) has no logarithm and is perfectly defined at a = 1, which is why
  # the hakstian_whalen row is described by its missing variance instead.
  #
  # An earlier version of this table asserted that hakstian_whalen "is in fact retained
  # at es = 1, se = 0". That was true before es_from_cronbach_alpha() gained se_defined
  # in the same release, and stopped being true when it did. V40 therefore stayed
  # silent on the HW boundary row -- the unexplained drop it exists to prevent -- while
  # its own bonett remedy string sent the reader to that very scale.
  #
  # If either route changes its boundary convention, this table has to change with it.
  for (rel1 in list(c("cronbach_alpha", "alpha", "alpha"),
                    c("omega", "omega", "omega"),
                    c("icc", "ICC", "icc"))) {
    cl <- rel1[1]; lab <- rel1[2]; which_scale <- rel1[3]
    if (!cl %in% colnames(x)) next
    if (identical(which_scale, "icc") && !icc_measure_ok) next
    sc <- as.character(switch(which_scale, "alpha" = alpha_to_es,
                              "omega" = omega_to_es, icc_to_es))
    drops <- switch(which_scale,
                    "alpha" = !identical(sc, "raw"),
                    "omega" = !identical(sc, "raw"),
                    TRUE)                    # icc: no scale keeps the boundary
    if (!drops) next
    why <- if (identical(sc, "bonett")) {
      "has no bonett transform (log(0) is undefined)"
    } else if (identical(which_scale, "alpha")) {
      sprintf(paste0("has no standard error on the %s scale: es_from_cronbach_alpha() ",
                     "refuses the variance at the boundary on every scale but raw"), sc)
    } else if (identical(which_scale, "omega")) {
      sprintf("is refused on the %s scale by es_from_omega(), which admits the boundary on raw only", sc)
    } else {
      "is refused by es_from_icc() on every analysis scale, its validity gate icc < 1 being scale-independent"
    }
    # Only name a rescuing scale where one exists.
    remedy <- switch(which_scale,
                     "alpha" = " - or set alpha_to_es to raw, the only scale on which 1 is a usable boundary",
                     "omega" = " - or set omega_to_es to raw, where 1 is a usable boundary",
                     "")
    v1 <- suppressWarnings(as.numeric(x[[cl]]))
    for (i in which(!is.na(v1) & v1 == 1)) {
      row_issues[[i]] <- c(row_issues[[i]], sprintf(
        paste0("[UNUSUAL] '%s' = 1 exactly. A perfect %s %s, so this row yields no ",
               "effect size and drops out of the pool. Verify the value - a reported ",
               "1.00 is usually a rounded 0.99x%s"),
        cl, lab, why, remedy))
    }
  }

  # V41: the ICC's lower bound depends on k. An ICC is bounded below by -1/(k - 1)
  # rather than by -1, because with k measurements a negative correlation has to be
  # shared among the k - 1 others, so k = 3 cannot go below -0.5 and k = 6 cannot go
  # below -0.2. V11's .bounded_columns entry uses the loose -1 because that list holds
  # static scalar bounds and this one varies per row, so V41 is a separate check, in
  # the same way that V32 and V34 handle values inside the V11 range that are still
  # non-identified.
  #
  # It matters because the SE carries (1 + (k-1)rho), which shrinks as rho goes
  # negative: icc = -0.8 at k = 3 gives se = 0.0495 against 0.2144 for an ordinary
  # icc = 0.8, which is 18.8x the meta-analytic weight for an arithmetically
  # impossible value. k = 2 is unaffected, since -1/(k-1) is exactly -1 there, so the
  # common test-retest case is unchanged.
  #
  # The floor is the bound on a SINGLE-measurement ICC and must not be applied to an
  # average-measures one. For an ICC(2,k)/ICC(3,k) the reported quantity is rho_k, and
  # the inverse Spearman-Brown map rho_1 = rho_k / (k - (k-1) rho_k) sends every rho_k
  # in [-1, 1] to a rho_1 at or above -1/(2k-1), which is above -1/(k-1); so no
  # in-range average-measures ICC can violate the single-measures floor. V41 was
  # declaring such values "[INVALID] Impossible ICC" and, under the default
  # correct_inputs = TRUE, destroying a number the route itself computes correctly --
  # the row losing at the same time the [INFO] step-down note it would otherwise have
  # carried, which is the contradictory-flag pair the V41/V25 interlock exists to
  # avoid. V41's own justification does not transfer either: the stepped-down
  # rho_1 >= -1/(2k-1) bounds the (1 + (k-1)rho) SE factor at >= k/(2k-1) >= 0.5, so
  # at most ~4x weight inflation rather than the 18.8x the single-measures case gives.
  icc_below_floor <- rep(FALSE, nrow(x))
  if ("icc" %in% colnames(x) && "n_measurements" %in% colnames(x)) {
    icc_v <- x[["icc"]]
    k_v   <- suppressWarnings(as.numeric(x[["n_measurements"]]))
    floor_v <- -1 / (k_v - 1)
    icc_below_floor <- !is.na(icc_v) & !is.na(k_v) & is.finite(k_v) & k_v >= 2 &
      icc_v >= -1 & icc_v < floor_v &        # >= -1: below that is V11's job
      !icc_is_avg_col                        # average-measures: a different bound
    for (i in which(icc_below_floor)) {
      msg <- sprintf(
        "[INVALID] Impossible ICC for %g measurements: 'icc' = %g is below the lower bound -1/(k-1) = %g",
        k_v[i], icc_v[i], floor_v[i])
      if (correct_inputs) {
        x[["icc"]][i] <- NA_real_
        msg <- paste0(msg, ", set to missing")
      }
      row_issues[[i]] <- c(row_issues[[i]], msg)
    }
  }

  # V25: negative alpha/icc kept, warn only
  for (rel in list(c("cronbach_alpha", "alpha"), c("icc", "ICC"))) {
    if (!rel[1] %in% colnames(x)) next
    neg <- which(!is.na(x[[rel[1]]]) & x[[rel[1]]] < 0)
    # An ICC below -1/(k-1) is not mathematically possible, so V25 must not describe
    # it as such; emitting both would put two contradictory flags on one row. V41 has
    # already reported it.
    if (identical(rel[1], "icc")) neg <- neg[!icc_below_floor[neg]]
    for (i in neg) {
      row_issues[[i]] <- c(row_issues[[i]],
        sprintf("[UNUSUAL] Negative %s: '%s' = %g. Mathematically possible but indicates serious measurement problems or an extraction error - verify against the primary report",
                rel[2], rel[1], x[[rel[1]]][i]))
    }
  }

  # V31: agreement-type ICC rows use a one-way variance approximation that assumes
  # negligible between-rater variance, so the SE is anti-conservative when raters
  # differ systematically. Informational and always active, like V18 and V21. An
  # icc_type of NA, or an absent column, resolves to the package default "agreement".
  #
  # Scoped to measure = "icc". The note is about the ICC SE, so it does not belong on
  # an alpha or omega run, and a COSMIN-style extraction sheet holding alpha, omega
  # and ICC columns side by side (one sheet, several properties, one pool per
  # property) would otherwise carry it into every run. Tier-1 checks are keyed on
  # which columns exist rather than on the requested measure, so the measure test has
  # to be explicit here, as it already is for V26. A NULL or empty measure leaves
  # direct callers unaffected.
  if ("icc" %in% colnames(x) && icc_measure_ok) {
    # icc_measure_ok, icc_type_eff, icc_k_col and icc_is_avg_col are computed once
    # above V40, which needs them too. Average-measures ICC is resolved first: whether
    # such a row survives at all decides whether the V31 SE note below is worth
    # emitting.
    k_col    <- icc_k_col
    icc_v    <- x[["icc"]]
    is_avg   <- icc_is_avg_col
    can_step <- is_avg & !is.na(k_col) & is.finite(k_col) & k_col >= 2 &
      abs(icc_v) <= 1
    dropped  <- is_avg & !can_step        # es_from_icc() sets these to NA

    for (i in which(can_step)) {
      rho1 <- .icc_step_down(icc_v[i], k_col[i])
      row_issues[[i]] <- c(row_issues[[i]], paste0(
        "[INFO] Average-measures ICC (", round(icc_v[i], 4), ", k = ", k_col[i],
        ") stepped down to the single-measurement ICC ", round(rho1, 4),
        " with Spearman-Brown, so it is on the same scale as the rest of the pool"))
    }
    for (i in which(dropped)) {
      # No "; " anywhere in this message: .flag_es_quality() splits the merged flag
      # string on that sequence, and a fragment carrying no quoted column name
      # matches every scope, so it would land in BOTH crude and adjusted (the V23
      # defect). The leading quoted 'n_measurements' is what routes this one.
      row_issues[[i]] <- c(row_issues[[i]], paste0(
        "[INVALID] Average-measures ICC needs 'n_measurements' to be stepped down ",
        "to a single-measurement ICC (Spearman-Brown) - set to NA rather than pooled ",
        "as if it were single-measures"))
    }

    # V31. Use .icc_is_agreement() rather than == "agreement": an ICC(2,k) row is
    # stepped down to single-measures and then computed with the same one-way
    # agreement SE, so it needs the same note. Rows dropped just above are excluded,
    # since there is no SE there to be anti-conservative about and the [INVALID] is
    # the only actionable message.
    #
    # Severity is [INFO], and the exact opening string "[INFO] ICC agreement-type SE"
    # is pinned by tests_save/checked/test-flags.R. The note is advisory because
    # icc_agreement_se defaults to "compute": the row keeps its computed SE and stays
    # in the pool unless the user opts into "drop". The coverage behind the note is
    # 0.82 at n = 20 falling to 0.14 at n = 1000, so the approximation degrades as
    # studies get bigger.
    agr <- which(!is.na(icc_v) & .icc_is_agreement(icc_type_eff) & !dropped)
    for (i in agr) {
      row_issues[[i]] <- c(row_issues[[i]],
        "[INFO] ICC agreement-type SE assumes negligible between-rater variance and can be anti-conservative when raters differ systematically, increasingly so at larger n (measured 95% coverage 0.82 at n = 20, 0.32 at n = 200, 0.14 at n = 1000). Supply the study's own 'icc_se' or 'icc_ci_lo'/'icc_ci_up', or set icc_agreement_se = 'drop' to leave such rows out of the pool (see ?es_from_icc)")
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
      # One precision for the whole triplet (see .triplet_precision): pricing each
      # bound from its own stored decimals turns a printed 0.90 into 1 dp and inflates
      # the pad tenfold, which is what let [0.90, 0.86] through while the identical
      # inversion at [0.91, 0.87] was caught. The 1e-12 keeps the boundary from being
      # decided by floating point.
      dp <- .triplet_precision(val[i], lo[i], up[i])
      pad <- 10^(-dp)
      lo[i] - up[i] > pad + 1e-12
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

    # V3: value outside its CI, re-read after the V2 cleanup, and rounding-aware. The
    # point value and the CI bounds may be reported to different decimal precisions,
    # so a value that merely rounds to a bound must not be destroyed: md = 0.1 at 1 dp
    # and lo = 0.15 at 2 dp are both consistent with a true md in [0.145, 0.155). Flag
    # the row, and NA it under correct_inputs, only when the value lies outside by
    # more than the combined rounding error of the value and the offending bound.
    val <- x[[tri$val]]; lo <- x[[tri$lo]]; up <- x[[tri$up]]
    cand <- which(!is.na(val) & !is.na(lo) & !is.na(up) & lo <= up &
                    (val < lo | val > up))
    outside <- cand[vapply(cand, function(i) {
      # Same one-precision-per-triplet rule as V2 above: half a unit in the last place
      # for the value plus half a unit for the offending bound, both priced at the
      # triplet's precision.
      dp <- .triplet_precision(val[i], lo[i], up[i])
      pad <- 10^(-dp)
      (val[i] < lo[i] - pad - 1e-12) || (val[i] > up[i] + pad + 1e-12)
    }, logical(1))]
    if (length(outside) > 0) {
      if (correct_inputs) {
        # keep_val triplets (see .ci_triplets): refuse the interval, keep the
        # coefficient, which has an interval-free SE to fall back on.
        keep_val <- isTRUE(tri$keep_val)
        if (!keep_val) x[[tri$val]][outside] <- NA_real_
        x[[tri$lo]][outside] <- NA_real_
        x[[tri$up]][outside] <- NA_real_
        for (i in outside) {
          row_issues[[i]] <- c(row_issues[[i]], if (keep_val) {
            sprintf(paste0("[INVALID] Value outside CI for '%s': the interval is set to ",
                           "missing and the coefficient kept, since its standard error ",
                           "does not have to come from the interval"), tri$val)
          } else {
            sprintf("[INVALID] Value outside CI for '%s' set to missing", tri$val)
          })
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

  # V33: reported proportion against its own counts. Rounding-aware and warn only,
  # with the data preserved, because two independently transcribed values disagree and
  # neither can be shown to be the wrong one. Both prop routes still compute, and the
  # generic cross-method discordance thresholds are far too coarse on the [0, 1] scale
  # to catch this.
  #
  # V33a: prop differs from n_cases/n_sample by more than the rounding error of the
  # reported prop. V33b, which runs only when n_cases is absent, is the count analogue
  # of the mean-granularity (GRIM) test: prop is not achievable as k/n_sample for any
  # integer k, again beyond the reported rounding. Messages carry the quoted 'prop' so
  # .v_flag_matches_scope routes them to the crude scope only.
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
      # The same 1e-12 slack V2/V3/V33 carry. Without it a value sitting EXACTLY on the
      # pad fires or stays silent according to binary representation rather than
      # anything about the data; with it the boundary case is uniformly silent, which
      # is the conservative direction for a warn-only check.
      if (mv[i] < lo[i] - pad_lo - 1e-12 || mv[i] > hi[i] + pad_hi + 1e-12) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INVALID] Mean outside its reported range for '%s': mean = %g not in [min = %g, max = %g]",
                  arm$sfx, mv[i], lo[i], hi[i]))
      }
    }
  }

  # V29: two nested gates on (max - min)/SD, in an if/else, so at most one of them
  # fires on a row:
  #
  #  (a) [INVALID], tested first -- the ratio is outside the bounds that hold for ANY
  #      n numbers whatsoever, whatever their distribution (derived at the test site
  #      below). Not implausible, arithmetically impossible.
  #  (b) [UNUSUAL], reached only when (a) passed -- the ratio is outside
  #      range_sd_lo_mult*xi .. hi_mult_n*xi, where xi(n) = 2*qnorm((n - .375)/(n + .25))
  #      is the Wan (2014) expectation for roughly normal data.
  #
  # Which gate is the binding one depends on n, and the two are not redundant because
  # they answer different questions. Measured over 2 <= n <= 1200 at the default
  # multipliers: the arithmetic CEILING sqrt(2(n-1)) is the tighter for n <= 71 -- at
  # n = 5 it is 2.83 against 2.5*xi = 5.90, so the band's ceiling is provably inert
  # there -- and again over n = 263-448, where the 2.5 -> 5 ramp below is
  # mid-interpolation. The arithmetic FLOOR is the tighter only for n <= 24: it tends
  # to 2 as n grows (1.99 at n = 500) while 0.5*xi keeps climbing past it (3.02 at
  # n = 500), so above n = 24 the distributional floor is what catches a collapsed
  # ratio.
  #
  # The CEILING of (b) has to grow with n. xi(n) grows like 2*qnorm(1 - 1/n) while the
  # studentised range of skewed or heavy-tailed data grows much faster, so a flat 2.5
  # is not distribution-free: measured over 1000 replicates of error-free data it fires
  # on 0.000 of every distribution tested at n <= 100, but on 6.0% / 50.2% of
  # lognormal(sigma = 2) and 2.7% / 16.5% of t(3) samples at n = 200 / 500. Simply
  # widening it is not the answer either -- a flat 2.5 -> 5 clears those false alarms
  # through n = 1000 but collapses detection of a genuine SD-as-SE from 0.997 to 0.000
  # at n = 10 and from 1.000 to 0.806 at n = 30. No flat multiplier is both.
  #
  # So: keep the validated 2.5 up to n = 100, and interpolate to 5 by n = 400, where
  # the heavy-tail false alarms measure 0.000. Detection is untouched, because the
  # errors this check names inflate the ratio by a factor of sqrt(n) -- about 22 x xi
  # at n = 500 against a ceiling of 5 x xi -- so the gap between the null and the error
  # widens far faster than the ceiling does.
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
      hi_mult_n <- range_sd_hi_mult * max(1, min(2, sqrt(nv[i] / 100)))
      # (max - min)/SD has EXACT bounds for any n values whatsoever, independent of
      # any distribution: it is maximised by one value at each extreme with the other
      # n-2 at the mean (one outlier against n-1 equal values reaches only sqrt(n)) and
      # minimised by a half-and-half two-point split (David, Hartley & Pearson 1954).
      #   upper  sqrt(2(n-1))          lower  sqrt(n(n-1) / (floor(n/2) ceiling(n/2)))
      # A ratio outside those is not implausible, it is arithmetically impossible, so
      # it is [INVALID] rather than [UNUSUAL]. This matters most at small n, where the
      # multiplier band below is provably inert: at n = 5 the upper gate 2.5*xi is 5.90
      # but the arithmetic maximum is only 2.83, so no value could ever trip it.
      hard_hi <- sqrt(2 * (nv[i] - 1))
      hard_lo <- sqrt(nv[i] * (nv[i] - 1) / (floor(nv[i] / 2) * ceiling(nv[i] / 2)))
      if (r_obs > hard_hi * (1 + 1e-8) || r_obs < hard_lo * (1 - 1e-8)) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[INVALID] Range/SD ratio arithmetically impossible for '%s': (max - min)/SD = %.3f, but for n = %d it must lie in [%.3f, %.3f] whatever the distribution. One of min, max, SD or n is wrong",
                  arm$sfx, r_obs, as.integer(nv[i]), hard_lo, hard_hi))
      } else if (r_obs < range_sd_lo_mult * xi || r_obs > hi_mult_n * xi) {
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[UNUSUAL] Range/SD ratio implausible for '%s': (max - min)/SD = %.2f (expected ~%.1f for n = %d, accepted band %.1f to %.1f). Possible SD-as-SE, variance-as-SD, or range/IQR mix-up",
                  arm$sfx, r_obs, xi, as.integer(nv[i]),
                  range_sd_lo_mult * xi, hi_mult_n * xi))
      }
    }
  }

  # V26: within-subject measure but reported CI matches the independent-groups SE
  # implied critical value compared to z and t. warn only
  #
  # The check is stated as needing no assumed pre-post correlation, and that is true
  # only above a band. For a CORRECTLY computed paired interval the implied critical
  # value is qt(.975, n-1) * sqrt(1 - r), which re-enters the z/t band whenever r is
  # low: measured, it false-fires for r below about 0.135 as n grows and up to r = 0.33
  # at n = 11, i.e. on exactly the outcomes (biomarkers, single-item VAS, symptom
  # counts) where a paired analysis is least redundant.
  #
  # Remedy taken: where the row supplies its own r_pre_post, reconstruct the paired SE
  # too and stand down when the reported half-width is ALSO consistent with it -- fire
  # only when the interval matches the independent-groups SE and does not match the
  # paired one. Rows that supply no r_pre_post keep the r-free behaviour, since there
  # is nothing to reconstruct with (defaulting r here would reintroduce the assumption
  # the check is advertised as not making).
  within_meas <- !is.null(measure) && tolower(measure) %in% c("mdw", "dw", "gw")
  if (within_meas) {
    z_crit <- stats::qnorm(0.975)
    r_pp_exp  <- if ("r_pre_post_exp" %in% colnames(x)) {
      suppressWarnings(as.numeric(x[["r_pre_post_exp"]]))
    } else rep(NA_real_, n)
    r_pp_nexp <- if ("r_pre_post_nexp" %in% colnames(x)) {
      suppressWarnings(as.numeric(x[["r_pre_post_nexp"]]))
    } else rep(NA_real_, n)
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
        # Suppress a match that the row's OWN reported correlation explains: with
        # sd_diff = sqrt(sd1^2 + sd2^2 - 2 r sd1 sd2) the paired SE is sd_diff/sqrt(n)
        # and a correct paired interval implies qt(.975, n-1) on it.
        if (isTRUE(matches)) {
          r_i <- if (!is.na(r_pp_exp[i])) r_pp_exp[i] else r_pp_nexp[i]
          if (!is.na(r_i) && is.finite(r_i) && abs(r_i) < 1) {
            var_d <- sde^2 + sdn^2 - 2 * r_i * sde * sdn
            if (is.finite(var_d) && var_d > 0) {
              paired_se <- sqrt(var_d) / sqrt(ne)
              crit_paired <- ((up - lo) / 2) / paired_se
              t_paired <- stats::qt(0.975, ne - 1)
              if ((abs(crit_paired - z_crit) / z_crit <= paired_as_indep_tol) ||
                  (abs(crit_paired - t_paired) / t_paired <= paired_as_indep_tol)) {
                matches <- FALSE
              }
            }
          }
        }
        if (isTRUE(matches)) {
          row_issues[[i]] <- c(row_issues[[i]],
            # The leading quoted column is load-bearing: .v_flag_matches_scope() routes
            # a Tier-1 message on the FIRST quoted token and returns TRUE for BOTH
            # scopes when it finds none, so without it a flag raised from the adjusted
            # user CI was copied verbatim into flags_crude, where no user CI exists.
            sprintf("[UNUSUAL] Within-subject design (%s) but the reported CI in '%s' matches an independent-groups SE (n_exp = n_nexp = %g, implied critical value %.2f). Paired pre/post data analysed as two independent groups overestimates the variance - verify the SE formula (a paired/within-subject SE should be used)",
                    tolower(measure), lo_col, ne, implied_crit))
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
        # Quoted: .v_flag_matches_scope() routes on the first quoted token, and with
        # the names bare it finds none and copies the flag into BOTH scopes -- so the
        # ancova_mean_sd_* pair landed in flags_crude, where no ANCOVA estimate exists,
        # and the mean_sd_* pair in flags_adjusted.
        row_issues[[i]] <- c(row_issues[[i]],
          sprintf("[UNUSUAL] Extreme SD ratio between arms: '%s' = %g, '%s' = %g (ratio: %.1f, threshold: %g). Possible variance/SD confusion or extraction error",
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
            # quoted for the scope-routing reason given at V12 above
            sprintf("[INFO] Identical SDs across arms: '%s' = '%s' = %g. Verify this is not the pooled or between-group SD",
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
            # quoted for the scope-routing reason given at V12 above
            sprintf("[INFO] Identical means across arms: '%s' = '%s' = %g. Verify this is not a copy-paste error",
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
          # The leading quoted 'mean_pre_sd_exp' is what routes this: it is a
          # crude-scope pre/post property and has nothing to say about the ANCOVA
          # estimate, but with no column named at all it was copied into both scopes.
          sprintf("[INFO] 'mean_pre_sd_exp'/'mean_pre_sd_nexp' over 'mean_sd_exp'/'mean_sd_nexp': SD_baseline / SD_endpoint = %.2f (below %.2f). Pre-post SMD formulas (Bonett, d_rm, d_av) target a baseline-SD-standardised estimand structurally different from the ANCOVA-on-endpoint-SD target at this variance regime - choice of estimand non-trivial",
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

  # V44: floor effect in the response distribution, which is the property that
  # governs how wrong the (n, k) closed-form reliability variance is.
  #
  # It is not the category count. The obvious gate would be "dichotomous items are
  # bad", but over 20 response distributions (simulations study 11e) the number of
  # categories has Spearman rho = -0.16 with the variance error c, because what
  # matters is where the item mass sits rather than how many bins it is cut into.
  # Three five-point formats alone span c = 0.92 to 1.69. Item skew is never
  # reported, but a floor effect is visible in the reported total-score mean, and its
  # standardised position
  #
  #     floor_position = (scale_mean / n_items - scale_min) / (n_response_categories - 1)
  #
  # has rho = -0.96 and bands identically across C = 2, 3, 5 and 7:
  #   < 0.20      c = 1.56-4.28, coverage .65-.87
  #   0.20-0.28   c = 1.08-1.44, coverage .90-.94
  #   > 0.28      c = 0.92-1.03, coverage .94-.97
  #
  # scale_min is required rather than defaulted. A 0-based instrument (PHQ-9, 0-3)
  # and a 1-based one (Likert, 1-5) differ by 1/(C-1), about 0.25, which spans two
  # bands, and the origin cannot be inferred from the mean: a moderate PHQ-9 sample
  # with an item mean of 1.67 reads as floor_position 0.22 under a 1-based assumption
  # when the truth is 0.56. Defaulting would false-alarm on one of the most common
  # instruments in the literature, so the check stays silent without scale_min rather
  # than guessing.
  #
  # Warn only, data preserved: a floor effect is a property of the construct in that
  # population, not an extraction error. What it changes is which outputs to trust,
  # and the pooled coefficient is not among the casualties. Study 11c found that
  # pooling with the closed-form variance lands within 0.0095 of oracle-variance
  # pooling, while I2 is inflated by roughly 1 - 1/c even at zero true heterogeneity.
  if (all(c("scale_mean", "n_items", "n_response_categories", "scale_min") %in% colnames(x))) {
    rel_present <- rep(FALSE, n)
    for (cf in c("cronbach_alpha", "omega")) {
      if (cf %in% colnames(x)) rel_present <- rel_present | !is.na(x[[cf]])
    }
    # Only the rows whose SE will actually COME from the closed form. A study that
    # reported its own SE or CI is unaffected by any of this, and flagging it would
    # be noise.
    se_reported <- rep(FALSE, n)
    for (cc in c("cronbach_alpha_se", "cronbach_alpha_ci_lo", "cronbach_alpha_ci_up",
                 "omega_se", "omega_ci_lo", "omega_ci_up")) {
      if (cc %in% colnames(x)) se_reported <- se_reported | !is.na(x[[cc]])
    }
    # ... and only where the SE-source switch for THAT coefficient asks for it. The
    # two switches point opposite ways by default -- alpha's closed form is on,
    # omega's is off -- so this cannot be one flag for the sheet: on the default omega
    # path there is no closed-form variance computed at all, the row is not even
    # pooled, and the message would describe a calculation that never happened. Per
    # coefficient rather than ORed across them, so a sheet holding both reports the
    # alpha row and stays quiet on the omega one.
    uses_closed_form <- rep(FALSE, n)
    for (cf in list(list(col = "cronbach_alpha", src = alpha_se_source,
                         own = c("cronbach_alpha_se", "cronbach_alpha_ci_lo",
                                 "cronbach_alpha_ci_up")),
                    list(col = "omega", src = omega_se_source,
                         own = c("omega_se", "omega_ci_lo", "omega_ci_up")))) {
      if (!cf$col %in% colnames(x)) next
      if (!identical(as.character(cf$src), "closed_form")) next
      own_se <- rep(FALSE, n)
      for (cc in cf$own) if (cc %in% colnames(x)) own_se <- own_se | !is.na(x[[cc]])
      uses_closed_form <- uses_closed_form | (!is.na(x[[cf$col]]) & !own_se)
    }
    sm  <- suppressWarnings(as.numeric(x$scale_mean))
    ki  <- suppressWarnings(as.numeric(x$n_items))
    ncg <- suppressWarnings(as.numeric(x$n_response_categories))
    smin <- suppressWarnings(as.numeric(x$scale_min))
    checkable <- which(rel_present & !se_reported &
                        !is.na(sm) & !is.na(ki) & !is.na(ncg) & !is.na(smin) &
                        ki > 0 & ncg > 1)
    for (i in checkable) {
      item_mean <- sm[i] / ki[i]
      smax <- smin[i] + ncg[i] - 1
      # An item mean outside the scale's own range is an unambiguous transcription
      # error, either an item mean entered where a total was asked for or the wrong
      # scale_min. Report that rather than banding a meaningless position.
      if (item_mean < smin[i] || item_mean > smax) {
        row_issues[[i]] <- c(row_issues[[i]], sprintf(
          "[INVALID] Scale mean outside its own range: 'scale_mean'/'n_items' = %.3f is not within [%s, %s] implied by 'scale_min' and 'n_response_categories' - check that 'scale_mean' is the TOTAL score mean and that 'scale_min' is the per-item minimum",
          item_mean, format(smin[i], scientific = FALSE), format(smax, scientific = FALSE)))
        next
      }
      fp <- (item_mean - smin[i]) / (ncg[i] - 1)
      if (!is.finite(fp) || fp >= floor_position_mild) next
      # The out-of-range check above is arithmetic and applies whatever the SE comes
      # from; the band below is a statement about the closed-form variance only.
      if (!uses_closed_form[i]) next
      severe <- fp < floor_position_severe
      # The message must contain no "; ", because .flag_es_quality splits the merged
      # flag string on "; " and an internal one would break the flag apart and
      # mis-route the pieces. The leading quoted 'scale_mean' is required too: the
      # scope router (.v_flag_matches_scope) reads the first quoted column name to
      # decide whether a flag belongs to the crude or the adjusted column under
      # split_adjusted, and a message quoting nothing matches every scope and is
      # duplicated into both. V23 documents the same rule.
      row_issues[[i]] <- c(row_issues[[i]], sprintf(
        "%s Floor effect in the response distribution ('scale_mean' implies a standardised mean position of %.2f, below %.2f). The closed-form (n, k) reliability standard error is anti-conservative here: measured c (true variance / reported variance) = %s, 95%% coverage %s. The POOLED coefficient and its Hartung-Knapp CI are essentially unaffected - what is not trustworthy is I-squared, tau-squared and the prediction interval, since I-squared is inflated by about 1 - 1/c even at zero true heterogeneity",
        if (severe) "[UNUSUAL]" else "[INFO]",
        fp, if (severe) floor_position_severe else floor_position_mild,
        if (severe) "1.6-4.3" else "1.1-1.4",
        if (severe) "0.65-0.87" else "0.90-0.94"))
    }
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
            # Leading quoted column for the reason given at V18: this is computed
            # from the baseline arm means, a crude-scope property, and its own advice
            # is to prefer the ANCOVA statistics -- so duplicating it into
            # flags_adjusted told the reader to switch to the estimate in front of them.
            sprintf("[INFO] 'mean_pre_exp' vs 'mean_pre_nexp': |standardised baseline imbalance| = %.2f (above %.2f). Endpoint-only SMD formulas inherit ANCOVA omitted-variable bias of magnitude r * imb_std - ANCOVA-adjusted statistics recommended at this imbalance level",
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
    # Which residual-SD column actually fired, quoted and leading: this is an
    # ANCOVA-scope property, and with no quoted column .v_flag_matches_scope() copied
    # it into flags_crude, where no back-transformation is being performed.
    which_sd <- "ancova residual SD"
    for (col_name in ancova_resid_cols) {
      if (col_name %in% colnames(x) && !is.na(x[[col_name]][i]) && x[[col_name]][i] > 0) {
        which_sd <- col_name
        break
      }
    }
    row_issues[[i]] <- c(row_issues[[i]], paste0(
      "[UNUSUAL] '", which_sd, "': ANCOVA residual SD provided but ", why,
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

  # V36: reliability induction, meaning the same reliability coefficient reported by
  # two different studies.
  #
  # "Reliability induction" (Vacha-Haase 1998; Vacha-Haase & Thompson 2011) is the
  # central criticism of reliability-generalization work: a primary study reports the
  # alpha printed in the test manual, or in the original validation paper, instead of
  # computing alpha in its own sample. Those rows are not independent estimates of
  # anything and must not be pooled, since the premise of RG is that reliability is a
  # property of the scores in a sample rather than a fixed property of the instrument.
  #
  # This is V23 applied to the reliability columns, and it needs its own entropy gate.
  # V23 requires at least templated_min_match non-integer values in the matched block,
  # which a reliability block can never supply: it has one continuous column (the
  # coefficient) and one integer column (n_items). Alpha is also usually printed to 2
  # decimals, where collisions are common by chance, since there are only about 30
  # plausible two-decimal values in [.70, .99] and a shared ".87" in a 30-study pool is
  # unremarkable. Hussey et al. (2025) make the same point from the other side: alphas
  # pile up at round thresholds.
  #
  # Measured false-positive rate, with alphas drawn independently from
  # runif(0.78, 0.93) so that no induction is present, over 400 pools per cell. The
  # figures are the share of pools containing at least one flag, then the share of rows
  # flagged:
  #
  #     k    2 dp          3 dp          4 dp
  #    10    23.8% / 5.0%  23.0% / 5.0%   4.8% / 1.0%
  #    30    92.0% /15.0%  95.3% /15.8%  26.5% / 1.9%
  #   100   100.0% /41.1% 100.0% /44.1%  95.8% / 6.4%
  #
  # These are reproducible rather than quoted: simulations/studies/10_reliability.R,
  # v36_false_positive_rate(), at a fixed seed. Re-run it before citing them.
  #
  # Two conclusions follow, and they decide the design.
  #
  # (1) Tightening the decimal gate does not work. Requiring 4 decimals still flags a
  #     100-study pool 95.8% of the time. The driver is pool size, not precision: with
  #     about 150 plausible 3-decimal values in [.78, .93], a 30-study pool has
  #     C(30,2)/150, roughly 3, expected collisions, so a match is the norm rather
  #     than the exception. The gate is therefore left as it is and the message is
  #     worded accordingly: V36 is a prompt to check the primary sources, not a
  #     finding, and its count must never be reported as a number of induced values.
  #
  # (2) The n_sample half of the gate earns its place, which was not obvious.
  #     Repeating the sweep with non-round sample sizes (runif over 50:500 rather than
  #     multiples of 50) drops the 2-decimal rate from 23.8 / 92.0 / 100.0% to
  #     1.8 / 5.3 / 49.3% at k = 10 / 30 / 100. Clause (b) is therefore not the driver
  #     of false positives; round sample sizes are, because they collide with each
  #     other as readily as 2-decimal alphas do.
  #
  # Two ways to clear the gate, both keyed on genuine improbability rather than on a
  # count of decimals alone:
  #   (a) the coefficient carries at least 3 decimals, since 0.8734 shared by two
  #       independent samples is a near-certain copy;
  #   (b) the coefficient carries at least 2 decimals and n_sample is identical too,
  #       because two different studies agreeing on both alpha and N is far stronger
  #       evidence than either alone, and is the signature of a manual-quoted value.
  if (isTRUE(enable_cross_row) && n >= 2) {
    rel_blocks <- list(
      alpha = c(coef = "cronbach_alpha", n = "n_sample"),
      omega = c(coef = "omega",          n = "n_sample"),
      icc   = c(coef = "icc",            n = "n_sample")
    )
    sid_r <- if ("study_id" %in% colnames(x)) as.character(x[["study_id"]]) else rep(NA_character_, n)
    for (bl in names(rel_blocks)) {
      cc <- unname(rel_blocks[[bl]]["coef"])
      nc <- unname(rel_blocks[[bl]]["n"])
      if (!cc %in% colnames(x)) next
      coef_v <- suppressWarnings(as.numeric(x[[cc]]))
      n_v <- if (nc %in% colnames(x)) suppressWarnings(as.numeric(x[[nc]])) else rep(NA_real_, n)
      for (gk in seq_along(v_groups)) {           # one pool per flag_group stratum
      ok <- v_groups[[gk]][is.finite(coef_v[v_groups[[gk]]])]
      if (length(ok) < 2) next
      for (key in unique(coef_v[ok][duplicated(coef_v[ok])])) {
        grp <- ok[coef_v[ok] == key]
        dec <- .count_decimals(key)
        for (i in grp) {
          others <- setdiff(grp, i)
          # a repeated value WITHIN one study is legitimate (subscales, timepoints)
          # and is Category H job; only cross-study repetition is induction
          others <- others[is.na(sid_r[i]) | is.na(sid_r[others]) | sid_r[others] != sid_r[i]]
          if (length(others) == 0) next
          same_n <- others[is.finite(n_v[i]) & is.finite(n_v[others]) & n_v[others] == n_v[i]]
          if (dec >= 3) {
            why <- sprintf("to %d decimals", dec)
            hit <- others
          } else if (dec >= 2 && length(same_n) > 0) {
            why <- sprintf("to %d decimals, and n_sample is identical too (%s)", dec,
                           format(n_v[i], scientific = FALSE))
            hit <- same_n
          } else next
          row_issues[[i]] <- c(row_issues[[i]], paste0(sprintf(
            # Lead with the quoted column name so .v_flag_matches_scope() routes
            # this to the crude scope only, following the rule V23 documents.
            # Unquoted, the router cannot identify the column and copies the flag into
            # flags_adjusted, a scope with no estimates at all for these measures.
            paste0("[INFO] Same '%s' (%s) reported by %s, %s. This is a PROMPT TO CHECK ",
                   "THE PRIMARY SOURCES, not evidence of anything: plausible coefficients ",
                   "occupy a narrow band, so in a pool of %d studies coincidental repeats are ",
                   "expected (measured with no induction present: at 3 decimals, 92%% of ",
                   "30-study pools contain at least one). What it is worth checking for is ",
                   "reliability induction - a study quoting a test manual or an earlier ",
                   "validation paper instead of computing the coefficient in its own sample. ",
                   "Induced values are not independent estimates and should not be pooled. ",
                   "Never report the number of these flags as a count of induced values"),
            cc, format(key, scientific = FALSE),
            paste(.row_ref(hit, sid_r), collapse = ", "), why, length(ok)),
            v_grp_note(gk)))
        }
      }
      }
    }
  }

  # V37: the item count is not constant across a pool that shares one instrument.
  #
  # A reliability-generalization review is by construction about one questionnaire,
  # and a questionnaire has a fixed number of items. A row whose n_items differs from
  # its peers is therefore a short form, a different version, or a transcription slip,
  # and in all three cases it should not be pooled without comment, because k enters
  # the sampling variance of every transform, Bonett and Hakstian-Whalen alike.
  #
  # Nothing else can see this. n_items appears only in .positive_columns(), which
  # checks its sign, and the SE consequence is far too small for the D2 SE-outlier
  # check: across all k the Bonett SE moves by at most sqrt(2), and an 8-versus-18
  # mix-up moves it 1.039x, well under D2's 3x gate.
  #
  # False-positive guard: a review may legitimately span several instruments, in which
  # case k varies by design and this check is meaningless. It therefore fires only when
  # there is a dominant item count, with more than half the rows agreeing, and then
  # only on the rows that disagree with it. A pool with no majority k is treated as
  # multi-instrument and left alone.
  if (isTRUE(enable_cross_row) && n >= 3 && "n_items" %in% colnames(x)) {
    k_v <- suppressWarnings(as.numeric(x[["n_items"]]))
    for (gk in seq_along(v_groups)) {             # one pool per flag_group stratum
      ok <- v_groups[[gk]][is.finite(k_v[v_groups[[gk]]]) & k_v[v_groups[[gk]]] > 0]
      if (length(ok) >= 3) {
      tab <- table(k_v[ok])
      modal_k <- as.numeric(names(tab)[which.max(tab)])
      if (max(tab) > length(ok) / 2 && length(tab) > 1) {
        for (i in ok[k_v[ok] != modal_k]) {
          row_issues[[i]] <- c(row_issues[[i]], paste0(sprintf(
            # quoted for the same scope-routing reason as V36 above
            paste0("[UNUSUAL] Item count differs from the rest of the pool: 'n_items' = %s ",
                   "where %d of %d studies report %s. In a reliability-generalization review the ",
                   "instrument is fixed, so this is a short form, a different version, or an ",
                   "extraction error - k enters the sampling variance of every transform, so ",
                   "verify it before pooling"),
            format(k_v[i], scientific = FALSE), max(tab), length(ok),
            format(modal_k, scientific = FALSE)), v_grp_note(gk)))
        }
      }
      }
    }
  }

  # V42: the measurement/rater count is not constant across an ICC pool.
  #
  # The ICC analogue of V37, with a weaker justification, which is why the majority
  # guard below is essential here rather than merely convenient. For alpha the
  # instrument fixes n_items, so any disagreement is a short form or a slip. For an ICC
  # the number of raters legitimately varies between studies, and for a single-measures
  # ICC that is not an estimand problem at all: k enters the sampling variance, which
  # is exactly where it belongs, so a varying k is correctly handled by the SE.
  #
  # What this catches is the transcription error, a k that disagrees with a pool which
  # otherwise agrees. In a test-retest review at k = 2 throughout, a stray k = 3 is
  # suspicious; in an inter-rater review with genuinely mixed panels there is no
  # majority and the check stays silent by construction.
  #
  # Nothing else can see it. n_measurements is otherwise only sign-checked, and the SE
  # consequence sits under D2's 3x gate while being far larger than the case V37 was
  # written for: at n = 100 the ICC SE moves 1.26x to 2.36x across k = 2..10 depending
  # on rho (1.394x at rho = 0.80, k = 2..6), against 1.039x for an 8-versus-18
  # item-count mix-up on alpha.
  #
  # Gated on measure for V31's reason: this is ICC advice, and a COSMIN sheet carrying
  # alpha and icc side by side must not collect it on the alpha run.
  if (isTRUE(enable_cross_row) && n >= 3 && icc_measure_ok &&
      "n_measurements" %in% colnames(x) && "icc" %in% colnames(x)) {
    km_v <- suppressWarnings(as.numeric(x[["n_measurements"]]))
    icc_p <- suppressWarnings(as.numeric(x[["icc"]]))
    for (gk in seq_along(v_groups)) {             # one pool per flag_group stratum
      gidx <- v_groups[[gk]]
      ok <- gidx[is.finite(km_v[gidx]) & km_v[gidx] >= 2 & is.finite(icc_p[gidx])]
      if (length(ok) >= 3) {
      tab <- table(km_v[ok])
      modal_k <- as.numeric(names(tab)[which.max(tab)])
      if (max(tab) > length(ok) / 2 && length(tab) > 1) {
        for (i in ok[km_v[ok] != modal_k]) {
          row_issues[[i]] <- c(row_issues[[i]], paste0(sprintf(
            paste0("[UNUSUAL] Measurement count differs from the rest of the pool: ",
                   "'n_measurements' = %s where %d of %d studies report %s. k enters the ",
                   "ICC sampling variance directly, so verify it before pooling - and if the ",
                   "designs really do differ, check that the reported ICCs are all ",
                   "single-measurement rather than a mix with average-measures values"),
            format(km_v[i], scientific = FALSE), max(tab), length(ok),
            format(modal_k, scientific = FALSE)), v_grp_note(gk)))
        }
      }
      }
    }
  }

  # V43: a pool mixing absolute-agreement and consistency ICCs.
  #
  # The ICC analogue of V38, and the strongest member of that family, because here the
  # arithmetic is byte-identical: icc_type has no computational effect on the estimate
  # or the standard error (asserted to 1e-15 in tests_save/checked/test-icc.R), so no
  # numeric check anywhere in the package can reveal the mix. Without this flag such a
  # pool is completely silent.
  #
  # ICC(2,1) and ICC(3,1) are different estimands. Absolute agreement charges
  # systematic rater differences against the reliability and consistency does not, so a
  # consistency ICC is at least as large as the agreement ICC on the same data, and
  # averaging them estimates neither.
  #
  # [INFO] and always active, for V38's reasons: nothing is mis-extracted, since each
  # value is correct for what its study reported, so there is nothing to verify and
  # [UNUSUAL] would misdescribe it. What is wrong is pooling them, which is an analyst
  # choice. And two different types in one pool are different estimands by definition,
  # so the check cannot false-fire.
  #
  # Only the model axis is compared. The "_average" suffix is stripped first, because
  # the single-versus-average distinction is resolved by the Spearman-Brown step-down
  # and reported by its own flag, so an ICC(2,1)/ICC(2,k) pool is not a mixed estimand.
  # An absent or NA icc_type resolves to "agreement", the documented default, so a pool
  # of blanks plus explicit "consistency" rows is a genuine mix and correctly fires.
  #
  # Gated on measure for V31's reason (see V42 just above).
  if (n >= 2 && icc_measure_ok && "icc" %in% colnames(x)) {
    icc_p <- suppressWarnings(as.numeric(x[["icc"]]))
    it <- if ("icc_type" %in% colnames(x)) {
      tt <- as.character(x[["icc_type"]])
      tt[is.na(tt)] <- "agreement"
      sub("_average$", "", .normalise_icc_type(tt, warn = FALSE))
    } else {
      rep("agreement", nrow(x))
    }
    it_known_all <- is.finite(icc_p) & !is.na(it)
    for (gk in seq_along(v_groups)) {             # one pool per flag_group stratum
    it_known <- rep(FALSE, n)
    it_known[v_groups[[gk]]] <- it_known_all[v_groups[[gk]]]
    present <- unique(it[it_known])
    if (length(present) > 1) {
      lbl <- c(agreement = "absolute agreement, ICC(2,1)",
               consistency = "consistency, ICC(3,1)")
      for (i in which(it_known)) {
        row_issues[[i]] <- c(row_issues[[i]], paste0(sprintf(
          paste0("[INFO] Pool mixes ICC estimands: this row is '%s' (%s) while the pool ",
                 "also contains '%s'. Consistency ignores systematic rater differences and ",
                 "absolute agreement charges them against the reliability, so a consistency ",
                 "ICC is the larger of the two on the same data and their average estimates ",
                 "neither. Note the arithmetic is identical for both, so nothing else in the ",
                 "output reveals this - split the analysis by icc_type, or enter it as a ",
                 "moderator and report the contrast"),
          it[i], unname(lbl[it[i]]),
          paste(setdiff(present, it[i]), collapse = "', '")), v_grp_note(gk)))
      }
    }
    }
  }

  # V38: a pool mixing different omegas.
  #
  # omega_type is an estimand, not a label. omega_total is the proportion of total
  # score variance due to all common factors; omega_hierarchical is the proportion due
  # to the general factor alone and is systematically smaller; omega_asymptotic and
  # subscale omegas are different quantities again. Averaging them produces a number
  # that estimates none of them.
  #
  # Severity is [INFO], matching E6 and E8: nothing here is mis-extracted, so there is
  # nothing for the reviewer to verify and [UNUSUAL] would misdescribe it. Both values
  # are correct; what is wrong is pooling them. That is an analyst choice, which is
  # what [INFO] marks in this package (compare E6, E8, V18 and V21).
  #
  # This is the omega analogue of E6 (mixed SMD standardizers) and E8 (mixed z
  # transforms), and it belongs at Tier 1 because it is a property of the input: the
  # estimand is fixed by what the primary study reported, not by any conversion the
  # package performs. It is always active, with no enable_cross_row gate, because
  # unlike the induction and item-count checks it cannot false-fire: two different
  # omega_type values in one pool are different estimands by definition.
  if (n >= 2 && "omega_type" %in% colnames(x) && "omega" %in% colnames(x)) {
    # Compare the estimand rather than its spelling: 'total', 'Total' and 'omega_t'
    # are one estimand, and flagging them as mixed would be a false positive on
    # capitalisation alone. Normalise without warning, since es_from_omega() already
    # warns about unrecognised values and this check must not warn twice.
    ot <- .normalise_omega_type(x[["omega_type"]], warn = FALSE)
    om <- suppressWarnings(as.numeric(x[["omega"]]))
    ot[is.na(ot) & is.finite(om)] <- "total"
    # "unspecified" is now the fallback for an unrecognised omega_type, so it must be
    # excluded here for the same reason V39 excludes it: a value we could not read is
    # not evidence of a mixed estimand, and counting it would fire on any typo.
    ot_known_all <- is.finite(om) & !is.na(ot) & ot != "unspecified"
    for (gk in seq_along(v_groups)) {             # one pool per flag_group stratum
    ot_known <- rep(FALSE, n)
    ot_known[v_groups[[gk]]] <- ot_known_all[v_groups[[gk]]]
    present <- unique(ot[ot_known])
    if (length(present) > 1) {
      for (i in which(ot_known)) {
        row_issues[[i]] <- c(row_issues[[i]], paste0(sprintf(
          paste0("[INFO] Pool mixes omega estimands: this row reports omega_type = ",
                 "'%s' while the pool also contains '%s'. omega_total (all common ",
                 "factors) and omega_hierarchical (general factor only) answer different ",
                 "questions and omega_h is systematically smaller, so their average ",
                 "estimates neither - split the analysis by omega_type, or enter it as a ",
                 "moderator and report the contrast"),
          ot[i], paste(setdiff(present, ot[i]), collapse = "', '")), v_grp_note(gk)))
      }
    }
    }
  }

  # V39: a pool mixing omega estimators.
  #
  # Omega is not one computation, and which one produced the number moves it far more
  # than most moderators do. Zinbarg et al. (2006) ran several estimators over
  # identical data: the first principal component overestimated omega_h by about +0.40
  # on average, against about +0.02 for a hierarchical (bifactor) CFA. Revelle &
  # Zinbarg (2009) record EFA returning omega_h = .04 where CFA returns exactly 0.0 on
  # the same data. A bias of that size is larger than any moderator effect a
  # reliability-generalization review is likely to report, so a pool mixing estimators
  # can manufacture a finding that is purely an artefact of which software each author
  # ran.
  #
  # The treatment of "unspecified" is the whole design. Most primary studies do not
  # name their estimator, so counting "unspecified" as a level would make this fire on
  # nearly every real dataset and be tuned out within a week. It fires only when two
  # different known estimators coexist, which is a fact about the data rather than
  # about what the reviewer failed to code. A pool of known plus unspecified stays
  # silent here; the unspecified rows are a reporting-quality observation for the
  # review to make, not a mixing error.
  #
  # Severity is [INFO], matching V38 and its Tier-2 cousins E6 and E8: nothing is
  # mis-extracted and there is nothing to verify, since every value is correct for the
  # estimator that produced it. What is wrong is pooling them, an analyst choice.
  if (n >= 2 && "omega_estimator" %in% colnames(x) && "omega" %in% colnames(x)) {
    oe <- .normalise_omega_estimator(x[["omega_estimator"]], warn = FALSE)
    om_v <- suppressWarnings(as.numeric(x[["omega"]]))
    known_all <- is.finite(om_v) & !is.na(oe) & oe != "unspecified"
    for (gk in seq_along(v_groups)) {             # one pool per flag_group stratum
    known <- rep(FALSE, n)
    known[v_groups[[gk]]] <- known_all[v_groups[[gk]]]
    present <- unique(oe[known])
    if (length(present) > 1) {
      for (i in which(known)) {
        row_issues[[i]] <- c(row_issues[[i]], paste0(sprintf(
          paste0("[INFO] Pool mixes omega estimators: this row was estimated by '%s' ",
                 "while the pool also contains '%s'. The estimator moves omega more than ",
                 "most moderators do - a first principal component overestimates omega_h by ",
                 "~0.40 where a bifactor CFA overestimates it by ~0.02 (Zinbarg et al. 2006) ",
                 "- so a mixed pool can produce a moderator effect that is an artefact of ",
                 "which software each author ran. Split by estimator, or enter it as a ",
                 "moderator and report the contrast"),
          oe[i], paste(setdiff(present, oe[i]), collapse = "', '")), v_grp_note(gk)))
      }
    }
    }
  }

  # V45: a pool mixing raw (covariance-matrix) and standardised (correlation-matrix)
  # alpha. The alpha analogue of V38/V43, and the third member of that family.
  #
  # They are different coefficients, not two names for one. Raw alpha is computed from
  # the covariance matrix and standardised alpha from the correlation matrix -- alpha
  # on standardised items -- and they coincide only when the item variances are equal,
  # separating as those variances spread. psych::alpha() prints both (raw_alpha and
  # std.alpha), so an extractor routinely has the pair in front of them with nothing
  # in the sheet recording which one was taken.
  #
  # Nothing numeric can reveal the mix: alpha_type does not enter the point estimate
  # or the (n, k) closed-form SE, so a pool holding both is silent everywhere else,
  # exactly as V43 is the only thing that can see an agreement/consistency ICC mix.
  #
  # "unspecified" is not a level, following V39 and NOT V43. V43 resolves a blank to
  # the package default because icc_type HAS one the user accepted by not overriding
  # it; there is no package default here, the coefficient comes from the study, and
  # most primary reports never say which they computed. Counting blanks as raw would
  # make this fire on nearly every real dataset, on a guess, and it would be tuned out.
  #
  # [INFO] for V38 and V39's reason: nothing is mis-extracted and there is nothing to
  # verify, since each value is correct for the matrix that produced it. What is wrong
  # is pooling them, which is an analyst choice. Always active -- two known types in
  # one pool are different estimands by definition, so it cannot false-fire.
  if (n >= 2 && "alpha_type" %in% colnames(x) && "cronbach_alpha" %in% colnames(x)) {
    at <- .normalise_alpha_type(x[["alpha_type"]], warn = FALSE)
    al_v <- suppressWarnings(as.numeric(x[["cronbach_alpha"]]))
    known_all <- is.finite(al_v) & !is.na(at) & at != "unspecified"
    for (gk in seq_along(v_groups)) {           # one pool per flag_group stratum
      known <- rep(FALSE, n)
      known[v_groups[[gk]]] <- known_all[v_groups[[gk]]]
      present <- unique(at[known])
      if (length(present) > 1) {
        for (i in which(known)) {
          row_issues[[i]] <- c(row_issues[[i]], paste0(sprintf(
            paste0("[INFO] Pool mixes alpha types: this row reports the %s-matrix ",
                   "alpha while the pool also contains %s-matrix alpha. Raw and ",
                   "standardised alpha are different coefficients, equal only when the ",
                   "item variances are, and alpha_type enters neither the estimate nor ",
                   "the standard error, so nothing else can see the mix. Split by type, ",
                   "or enter it as a moderator and report the contrast"),
            at[i], paste(setdiff(present, at[i]), collapse = " and ")), v_grp_note(gk)))
        }
      }
    }
  }

  # V35: 2x2 event rates spanning a wide range in a correlation pool.  #
  # This is a disclosure about the tetrachoric route's precision, not a data error and
  # not a suggestion to use a different estimand; there is no alternative, as
  # ?convert_df explains under table_2x2_to_cor.
  #
  # The tetrachoric is estimated by solving a bivariate-normal probability for rho, and
  # how sharply the data identify rho depends on the margins: dp11/drho at rho = 0.30
  # falls from 0.167 at a 50% event rate to 0.068 at 10% and 0.016 at 2%, a roughly
  # tenfold collapse in identification. Its standard error inflates correspondingly, so
  # a rare-event study contributes far less precision to an inverse-variance pool than
  # its sample size suggests. That is correct behaviour, since the information really
  # is not there, but in a pool mixing common and rare outcomes it quietly concentrates
  # the weight on the common-outcome studies. Worth stating; nothing to fix.
  #
  # This is not a phi gate. phi is not offered as an alternative, and margin
  # heterogeneity is precisely the condition under which phi would be least defensible.
  if (isTRUE(enable_cross_row) && n >= 2 &&
      !is.null(measure) && length(measure) == 1 && measure %in% c("r", "z")) {
    cell_cols <- c("n_cases_exp", "n_controls_exp", "n_cases_nexp", "n_controls_nexp")
    if (all(cell_cols %in% colnames(x))) {
      cm <- suppressWarnings(matrix(as.numeric(as.matrix(x[, cell_cols, drop = FALSE])),
                                    nrow = n, dimnames = list(NULL, cell_cols)))
      ok2x2_all <- rowSums(!is.finite(cm)) == 0 & rowSums(cm) > 0
      ev <- rep(NA_real_, n)
      ev[ok2x2_all] <- (cm[ok2x2_all, "n_cases_exp"] + cm[ok2x2_all, "n_cases_nexp"]) /
                       rowSums(cm[ok2x2_all, , drop = FALSE])
      for (gk in seq_along(v_groups)) {           # one pool per flag_group stratum
      ok2x2 <- rep(FALSE, n)
      ok2x2[v_groups[[gk]]] <- ok2x2_all[v_groups[[gk]]]
      if (sum(ok2x2) >= 2) {
        rng <- range(ev[ok2x2])
        # Fire only when the pool BOTH spans a wide range and reaches into the region
        # where identification is materially weaker; a pool sitting entirely at 0.3-0.6
        # has nothing to disclose.
        if ((rng[2] - rng[1]) >= margin_range_min && (rng[1] < 0.10 || rng[2] > 0.90)) {
          for (i in which(ok2x2)) {
            row_issues[[i]] <- c(row_issues[[i]], paste0(sprintf(
              paste0("[INFO] Event rates across the 2x2 rows span %.2f to %.2f. The ",
                     "tetrachoric correlation is weakly identified at extreme margins ",
                     "(identification falls ~10x from a 50%% to a 2%% event rate), so ",
                     "its standard error inflates there and rare-outcome rows carry ",
                     "less weight than their sample size implies. This row's event ",
                     "rate is %.2f. Correct behaviour, not an extraction error"),
              rng[1], rng[2], ev[i]), v_grp_note(gk)))
          }
        }
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

#' Which transform each route put in the `z` column
#'
#' `measure = "z"` is meant to hold one quantity, Fisher's z, so that the column can
#' be pooled. It does not. The correlation and binary families report
#' \eqn{z = atanh(r)}, while the SMD family under \code{smd_to_cor = "viechtbauer"}
#' (the default) reports a variance-stabilising transform, which is a different
#' function of the correlation. Measured on the same data at a point-biserial
#' \eqn{\rho = 0.75}, that route reports z = 1.0925 while \code{atanh()} of its own r
#' is 1.7468.
#'
#' This is not the estimand difference between \code{viechtbauer} (the biserial) and
#' \code{lipsey_cooper} (the point-biserial), which is a documented choice and is what
#' study 01 of the simulation programme measures. It is the separate fact that the z
#' column is not \code{atanh()} of the r column on every route, so a pool mixing
#' families mixes two transforms.
#'
#' The classification is derived from the output rather than from a list of route
#' names. Every method frame carries both \code{r} and \code{z}, so each route is
#' classified by asking whether its own z equals \code{atanh()} of its own r. A
#' hardcoded map of route names would misclassify any new route without saying so,
#' and would drift out of date as route names change.
#'
#' @param res the list of per-method data frames built by \code{convert_df()}
#' @return named character vector, names are \code{info_used} values, entries
#'   "fisher", "vst" or NA when a route offered no row to classify
#' @noRd
.z_transform_by_route <- function(res) {
  out <- character(0)
  for (k in seq_along(res)) {
    f <- res[[k]]
    if (!is.data.frame(f) || !all(c("r", "z", "info_used") %in% names(f))) next
    ok <- which(is.finite(f$r) & is.finite(f$z) & abs(f$r) < 1 & abs(f$r) > 1e-8)
    if (!length(ok)) next
    # Classify on ALL usable rows, not the first: a route that agreed with atanh()
    # on one row and not another would be a third thing, and must not be recorded
    # as either.
    is_fisher <- abs(f$z[ok] - atanh(f$r[ok])) < 1e-8
    lab <- if (all(is_fisher)) "fisher" else if (!any(is_fisher)) "vst" else "mixed"
    for (nm in unique(as.character(f$info_used[ok]))) {
      if (is.na(nm) || !nzchar(nm)) next
      # A route seen twice must not flip label; disagreement is recorded as "mixed"
      # rather than resolved silently.
      out[nm] <- if (!is.null(out[nm]) && !is.na(out[nm]) && out[nm] != lab) "mixed" else lab
    }
  }
  out
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
      # Use .group_label() here rather than the raw study_id: a study_id containing
      # "; ", such as "Huang; 2017", is split by the flag merge, which truncates the
      # message and mis-routes the leading fragment. That fragment carries no quoted
      # column name, so it matches every scope and lands in both crude and adjusted.
      # This covers V23 (:1379), V36 (:1455) and Category H (:2767) in one place.
      sprintf("%s (row %d)", .group_label(sid_i), as.integer(i))
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

#' Routes whose nnt / rd output is on a RATE scale rather than a risk scale
#'
#' \code{es_from_cases_time()} (Mayne et al. 2006) writes an incidence rate difference
#' into the \code{rd} columns and a person-time NNT -- a DURATION -- into the
#' \code{nnt} ones. Neither is a probability: an IRD is events per person-time and is
#' unbounded, and |NNT| < 1 is ordinary rather than impossible once the denominator is
#' time. Every risk-scale bound in \code{.flag_bounds_violations()} (B3, B6, B6b) has
#' to stand down on these rows, and E4 exists to warn that the two scales are in one
#' pool. One definition, so the three checks and the warning cannot drift apart.
#'
#' @return character vector of info_used values
#' @noRd
.rate_scale_methods <- function() "cases_time"


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
                                     n_cov_ancova = NULL,
                                     n_covariates = NULL) {
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
    #
    # Skipped for nnt, whose CI is genuinely discontinuous (Altman 1998). NOT skipped
    # for the exp-displayed ratios: that the SE is on the log scale while the CI is on
    # the natural one is not a reason to abandon the check, it is the transform to
    # apply -- log(ci_up) - log(ci_lo) against 2 z se, exactly the back-transform C5,
    # D1, E3 and E5 already perform in this file. Skipping it left the single most
    # likely user-input error for OR/RR/IRR/HR, a standard error given on the natural
    # scale instead of the log scale, with no Tier-2 detector at all, while the
    # identical mistake on the additive scale was caught.
    log_scale_a6 <- isTRUE(exp) && !is.null(measure) &&
      measure %in% c("logor", "logrr", "logirr", "loghr")
    skip_a6 <- identical(measure, "nnt") || (isTRUE(exp) && !log_scale_a6)
    if (log_scale_a6 &&
        (is.na(ci_lo[i]) || is.na(ci_up[i]) ||
         !is.finite(ci_lo[i]) || !is.finite(ci_up[i]) ||
         ci_lo[i] <= 0 || ci_up[i] <= 0)) {
      skip_a6 <- TRUE          # no log-scale width to compare
    }
    if (!skip_a6 &&
        !is.na(es[i]) && !is.na(se[i]) && !is.na(ci_lo[i]) && !is.na(ci_up[i]) &&
        is.finite(es[i]) && is.finite(se[i]) && is.finite(ci_lo[i]) && is.finite(ci_up[i]) &&
        se[i] > 0 && ci_lo[i] <= ci_up[i]) {
      ci_width <- if (log_scale_a6) {
        log(ci_up[i]) - log(ci_lo[i])
      } else {
        ci_up[i] - ci_lo[i]
      }
      # package CIs: qt for d/g/md/r and dw/gw/mdw, qnorm otherwise
      # user CIs: checked against z width, t width also accepted
      # "rp" builds its CI on the t distribution like the others, but on the
      # regression residual df rather than N - 2: es_from_linreg_t() uses
      # df = n_sample - n_covariates - 2 (R/es_from_REGRESSION.R:249) and
      # rp +/- qt(.975, df) * rp_se (:268-269). Omitting it from qt_measures makes A6
      # compare every rp CI against the Wald-z width and fire a false [DISCORDANT] at
      # small n, up to N = 23 at 1 covariate, N = 27 at 3 and N = 37 at 10. The
      # crossover scales with n_covariates, so adding "rp" without subtracting the
      # covariates would only move the boundary. "zp" is excluded because it is built
      # with qnorm (:276-277).
      qt_measures <- c("d", "g", "dw", "gw", "md", "mdw", "r", "rp")
      n_i <- if (!is.null(n_total)) n_total[i] else NA_real_
      # paired measures: df = n_pairs - 1 (n_total double counts the subjects)
      within_measures <- c("dw", "gw", "mdw")
      n_exp_i <- if (!is.null(n_exp)) n_exp[i] else NA_real_
      n_nexp_i <- if (!is.null(n_nexp)) n_nexp[i] else NA_real_
      n_cov_i <- if (!is.null(n_covariates) && length(n_covariates) >= i) {
        n_covariates[i]
      } else {
        NA_real_
      }
      if (!is.null(measure) && measure %in% within_measures &&
          !is.na(n_exp_i) && is.finite(n_exp_i) && n_exp_i > 2) {
        t_df <- n_exp_i - 1
        n_eff <- n_exp_i
      } else if (!is.null(measure) && identical(measure, "rp") &&
                 !is.na(n_i) && !is.na(n_cov_i) && is.finite(n_cov_i) &&
                 n_i - n_cov_i - 2 > 0) {
        t_df <- n_i - n_cov_i - 2
        n_eff <- n_i
      } else if (!is.na(n_i) && n_i > 4) {
        t_df <- n_i - 2
        n_eff <- n_i
      } else {
        t_df <- NA_real_
        n_eff <- n_i
      }
      is_user_es <- !is.null(info_used) && !is.na(info_used[i]) &&
                    grepl("^user_(es|input)", info_used[i])
      # P7: Glass rows (smd_denom = control/control_robust) build their d/g CIs on
      # the control-arm df, following the SMD1 convention, rather than on the pooled
      # N - 2. The Glass df applies only when the selected estimate actually came from
      # the endpoint-means family, the only methods that honour smd_denom. When the
      # hierarchy picked a non-means method (cohen_d, etasq, t/F, ANCOVA, medians and
      # so on), the CI is built on the pooled N - 2 regardless of the smd_denom input
      # column, so keying the Glass df off the raw input alone would false-fire the A6
      # CI-width flag on valid rows.
      glass_row <- !is_user_es && !is.null(smd_denom) &&
        length(smd_denom) >= i && !is.na(smd_denom[i]) &&
        smd_denom[i] %in% c("control", "control_robust") &&
        !is.null(measure) && measure %in% c("d", "g") &&
        !is.null(info_used) && !is.na(info_used[i]) &&
        info_used[i] %in% c("means_sd", "means_se", "means_ci") &&
        !is.na(n_nexp_i) && is.finite(n_nexp_i) && n_nexp_i > 2
      if (glass_row) t_df <- n_nexp_i - 1
      # P9: an adjusted row's d/g interval is built by .es_from_d() on
      # qt(.975, n_exp + n_nexp - 2 - n_cov_ancova), so A6 must expect that df too.
      # Otherwise it compares the row's own CI against a wider-df expectation and fires
      # [DISCORDANT] on a perfectly consistent interval. The gap clears the tolerance
      # only at small n (17% against a 10% tolerance at N = 9, q = 3), but it is a pure
      # false positive when it does. The adjusted routes are exactly those whose
      # info_used starts with "ancova" or ends in "_adj" (see the es_from_ancova_*,
      # es_from_cohen_d_adj and es_from_etasq_adj entry points). User-entered rows are
      # excluded, because the source's own CI construction is unknown and is checked
      # against the Wald-z width instead.
      if (!is_user_es && !is.null(info_used) && !is.na(info_used[i]) &&
          grepl("^ancova|_adj$", info_used[i]) &&
          !is.null(n_cov_ancova) && length(n_cov_ancova) >= i &&
          !is.na(n_cov_ancova[i]) && is.finite(n_cov_ancova[i]) &&
          n_cov_ancova[i] > 0 && !is.na(t_df) && t_df - n_cov_ancova[i] > 0) {
        t_df <- t_df - n_cov_ancova[i]
      }
      if (is_user_es || log_scale_a6) {
        # log_scale_a6: the package and every source build a ratio CI as
        # exp(log(es) +/- z se), so the expectation is the Wald-z width on the log
        # scale and none of the qt_measures df logic applies.
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
      # P10: EVERY package-computed `r` CI is a back-transformed z-scale interval,
      # asymmetric about r and narrower than 2 z se_r, never r +/- qt se:
      #   es_from_pearson_r() / es_from_fisher_z() / es_from_spearman_rho() and the
      #     user r/z branches -- Fisher's z, back-transformed as tanh(); pearson/fisher
      #     build z_se as 1/sqrt(n - 3) (metafor ZCOR), spearman by the delta map;
      #   .contingency_to_cor() / .tet_r(), .or_to_cor() -- Fisher's z, tanh();
      #   .smd_to_cor(), viechtbauer branch (the default) -- the variance-stabilising
      #     z with a = sqrt(dnorm(qnorm(p))) / (p(1-p))^(1/4), r = tanh(z/a)/a;
      #   .smd_to_cor(), lipsey_cooper branch -- Fisher's z, tanh() (it was the one
      #     route still returning r +/- qt se on the r scale, and escaped [-1, 1] at
      #     small n like the former pearson_r interval; both were switched).
      # The primary qt expectation formed above is therefore never the construction of
      # a package r row. The gap exceeds A6's own tolerance at small-to-moderate n, so
      # A6 was reporting a [DISCORDANT] extraction error on 19 of the 170 rows of the
      # package's own df.haza at measure = "r", all of them correct. Accept any of the
      # back-transformed widths as an alternative expectation: the Fisher width on the
      # delta-mapped z SE (exact for the tetrachoric, OR, spearman and lipsey_cooper
      # routes), the Fisher width on 1/sqrt(n - 3) (pearson_r / fisher_z: at n = 5 the
      # delta map gives 1/sqrt(n - 1) = 0.50 against 0.71, which fired falsely), and
      # the variance-stabilising width (0.1% from its delta reconstruction at n = 26).
      # This block does not receive smd_to_cor or the route, so all alternatives are
      # accepted on every package-computed r row. qnorm(.975) rather than the block's
      # own z_crit, because every back-transformed interval is built on the normal
      # quantile whatever df A6 would otherwise have chosen. "rp" is excluded:
      # es_from_linreg_t() really does build rp +/- qt(.975, n - q - 2) se, so it keeps
      # the t expectation.
      if (fires_a6 && !is_user_es && identical(measure, "r") &&
          is.finite(es[i]) && abs(es[i]) < 1 && se[i] > 0) {
        zc_bt <- stats::qnorm(0.975)
        alt_widths <- numeric(0)
        # Fisher / tetrachoric / spearman / lipsey_cooper: z = atanh(r), z_se = r_se / (1 - r^2)
        zs_f <- se[i] / (1 - es[i]^2)
        if (is.finite(zs_f)) {
          alt_widths <- c(alt_widths,
                          tanh(atanh(es[i]) + zc_bt * zs_f) -
                          tanh(atanh(es[i]) - zc_bt * zs_f))
        }
        # pearson_r / fisher_z: z_se = 1/sqrt(n - 3)
        if (!is.na(n_i) && is.finite(n_i) && n_i > 3) {
          zs_n <- 1 / sqrt(n_i - 3)
          alt_widths <- c(alt_widths,
                          tanh(atanh(es[i]) + zc_bt * zs_n) -
                          tanh(atanh(es[i]) - zc_bt * zs_n))
        }
        # Variance-stabilising (viechtbauer): z = a atanh(a r), r = tanh(z/a)/a
        p_bt <- if (!is.na(n_exp_i) && !is.na(n_nexp_i) &&
                    is.finite(n_exp_i) && is.finite(n_nexp_i) &&
                    n_exp_i + n_nexp_i > 0) {
          n_exp_i / (n_exp_i + n_nexp_i)
        } else {
          0.5
        }
        if (p_bt > 0 && p_bt < 1) {
          a_bt <- sqrt(stats::dnorm(stats::qnorm(p_bt))) / (p_bt * (1 - p_bt))^(1 / 4)
          ar <- a_bt * es[i]
          if (is.finite(a_bt) && a_bt > 0 && abs(ar) < 1) {
            z_bt  <- a_bt * atanh(ar)
            zs_bt <- (a_bt^2 / (1 - ar^2)) * se[i]
            if (is.finite(zs_bt)) {
              alt_widths <- c(alt_widths,
                              (tanh((z_bt + zc_bt * zs_bt) / a_bt) -
                               tanh((z_bt - zc_bt * zs_bt) / a_bt)) / a_bt)
            }
          }
        }
        alt_widths <- alt_widths[is.finite(alt_widths) & alt_widths > 0]
        if (length(alt_widths) &&
            any(abs(ci_width - alt_widths) / alt_widths <= tol)) {
          fires_a6 <- FALSE
        }
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

  # Rows whose estimate is on the rate scale, not the risk scale (see
  # .rate_scale_methods): every [-1, 1] / |NNT| >= 1 bound below is a statement about
  # probabilities and does not apply to them. With info_used absent (direct callers,
  # tests_save/checked/test-flags.R) nothing is treated as rate-scale, so the risk
  # bounds keep their pre-existing reach.
  rate_row <- rep(FALSE, n)
  if (!is.null(info_used) && length(info_used) > 0) {
    iu <- rep(as.character(info_used), length.out = n)
    rate_row <- !is.na(iu) & iu %in% .rate_scale_methods()
  }

  for (i in seq_len(n)) {
    if (is.na(es[i])) next
    msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""

    # B1: Correlation outside [-1, 1]
    #
    # "rp" is included alongside "r": a partial correlation is bounded identically
    # and its interval is built identically (es +/- t*se on the r scale), so every
    # correlation check below applies to it unchanged. The gates used to read
    # `measure == "r"`, which left every rp pool silently unchecked -- e.g.
    # linreg_t = 2.5, n_sample = 8, n_covariates = 1 gives rp = 0.745 with
    # CI [0.234, 1.256] and no bound flag at all.
    cor_scale <- measure %in% c("r", "rp")
    cor_lab   <- if (identical(measure, "rp")) "partial r" else "r"
    if (cor_scale && is.finite(es[i]) && abs(es[i]) > 1) {
      flags[[i]] <- c(flags[[i]],
        paste0("[INVALID] ", cor_lab, " outside [-1, 1]: ", cor_lab, " = ",
               round(es[i], 3), msuf))
    }

    # B1b: raw-scale Wald CI escaping [-1, 1] while the point estimate is valid.
    # Same geometry as B7b/B8b for alpha and ICC: r is bounded identically and its
    # interval is built identically (es +/- z*se), so a symmetric Wald interval can
    # leave the parameter space near |r| = 1 even when r itself is fine. B1 above tests
    # only the point estimate, so without this the bound would pass through unnoticed.
    #
    # [INFO] rather than [UNUSUAL], on purpose. The data and the point estimate are
    # sound, so there is nothing to verify and nothing extracted wrongly. This is a
    # property of the interval metaConvert constructs, which makes the flag a
    # disclosure rather than a request to check the extraction.
    #
    # The tetrachoric route does not reach here: its r-scale interval is the tanh
    # back-transform of the Fisher-z interval (.tet_r in internal_multiple_formulas.R),
    # which cannot leave (-1, 1). That construction took the escape rate from 45.8% of
    # tables to 0.0% with coverage unchanged (96.1% -> 96.2%). Every package-computed
    # r interval is now built that way (pearson_r, fisher_z, spearman, the OR
    # conversions and both smd_to_cor branches), so B1b is reachable only through a
    # USER-supplied r interval, which is handed back verbatim: it stays as the
    # backstop for that case.
    if (cor_scale && is.finite(es[i]) && abs(es[i]) <= 1) {
      lo_out <- !is.na(ci_lo[i]) && is.finite(ci_lo[i]) && ci_lo[i] < -1
      up_out <- !is.na(ci_up[i]) && is.finite(ci_up[i]) && ci_up[i] > 1
      if (lo_out || up_out) {
        # Report BOTH bounds when both escape. This used to print only the upper
        # one, which understates a doubly-invalid interval.
        side <- if (lo_out && up_out) {
          paste0("lower bound is below -1 (", round(ci_lo[i], 3),
                 ") and upper bound exceeds 1 (", round(ci_up[i], 3), ")")
        } else if (up_out) {
          paste0("upper bound exceeds 1 (", round(ci_up[i], 3), ")")
        } else {
          paste0("lower bound is below -1 (", round(ci_lo[i], 3), ")")
        }
        flags[[i]] <- c(flags[[i]],
          paste0("[INFO] Correlation CI ", side,
                 " while ", cor_lab, " is valid - the symmetric Wald interval escapes the parameter space, interpret the interval with caution", msuf))
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
    #
    # True of a RISK-based NNT, where |RD| <= 1 forces |NNT| >= 1, and false of the
    # person-time NNT es_from_cases_time() computes: there the NNT is a duration, and
    # any incidence rate difference above one event per unit of time gives |NNT| < 1.
    if (measure == "nnt" && !rate_row[i] && is.finite(es[i]) &&
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

      # The bound is SIDED. With the baseline risk BR fixed and the treated risk in
      # [0, 1], RD = BR - p_t ranges over [BR - 1, BR]: the benefit side is capped at
      # BR, giving |NNT| >= 1/BR, but the harm side is capped at 1 - BR, giving
      # |NNH| >= 1/(1 - BR), a completely different number at small BR. Applying the
      # benefit-side floor to both signs stamped [INVALID] on every NNH between
      # 1/(1-BR) and 1/BR -- at BR = 0.05 a false window from 1.05 to 20, i.e. on
      # essentially every appreciable harm signal, including tables the package itself
      # had just produced.
      # ...and it is SCALE-SPECIFIC. The risk floor is a statement about probabilities
      # (RD = BR - p_t, both in [0, 1]); a person-time NNT is 1/IRD, a DURATION in
      # person-time units with no probability interpretation and its own floor below.
      # Neither bounds the other, so each branch is gated on rate_row -- the same gate
      # B3, B6 and B6b already carry. Without it a legitimate person-time NNT of 1.0
      # person-year was stamped [INVALID] against a risk floor of 3.3, and (because
      # baseline_rate is now derived from n_cases_nexp / time_nexp for EVERY row) a
      # risk-based 2x2 NNT was stamped against the person-time floor, in a message that
      # contradicted itself by printing "person-time NNT ... (from 2x2)".
      if (!rate_row[i] && !is.na(br) && br > 0 && br < 1) {
        min_nnt <- if (es[i] > 0) 1 / br else 1 / (1 - br)
        side_txt <- if (es[i] > 0) "baseline_risk" else "1 - baseline_risk"
        if (abs(es[i]) < min_nnt) {
          flags[[i]] <- c(flags[[i]],
            paste0("[INVALID] Minimum possible NNT = ", round(min_nnt, 1),
                   " (", side_txt, " = ",
                   round(if (es[i] > 0) br * 100 else (1 - br) * 100, 1), "%)",
                   " - reported |NNT| = ", round(abs(es[i]), 3), msuf))
        }
      }
      # Benefit side only. On the harm side the incidence rate ratio is unbounded
      # above, so |IRD| is unbounded and a rate-based NNH has no positive lower bound
      # at all.
      if (rate_row[i] && !is.na(brate) && brate > 0 && es[i] > 0) {
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
    # Not applied to a rate-scale row: es_from_cases_time() writes an incidence rate
    # DIFFERENCE into the rd columns, and 3 events per person-year is arithmetically
    # valid rather than impossible.
    if (measure == "rd" && !rate_row[i] && is.finite(es[i]) && abs(es[i]) > 1) {
      flags[[i]] <- c(flags[[i]],
        paste0("[INVALID] RD outside [-1, 1]: rd = ", round(es[i], 3), msuf))
    }

    # B6b: RD CI bound outside [-1, 1]
    if (measure == "rd" && !rate_row[i]) {
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
                                prop_to_es = "raw",
                                alpha_to_es = "bonett",
                                icc_to_es = "bonett") {
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
    # "rp" included for the same reason as B1/B1b: a partial correlation is on the
    # same bounded scale, so r_max means the same thing for it.
    if (measure %in% c("r", "rp") &&
        !is.na(es[i]) && is.finite(es[i]) &&
        abs(es[i]) > opts$r_max && abs(es[i]) <= 1) {
      flags[[i]] <- c(flags[[i]],
        paste0("[UNUSUAL] High correlation: |",
               if (identical(measure, "rp")) "partial r" else "r", "| = ",
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

    # C7: near-perfect Cronbach's alpha, compared on the coefficient scale rather
    # than the analysis scale. An `es > 0 && es <= 1` guard is unreachable under the
    # Bonett default, since ln(1 - alpha) is negative for every valid alpha, so the
    # check would be dead exactly where it is needed. Worse, the Bonett values that do
    # satisfy it, (0.99, 1], are strongly negative alphas, which it would then announce
    # as "near-perfect".
    if (measure == "alpha" && !is.na(es[i]) && is.finite(es[i])) {
      a_raw <- .to_reliability_scale(es[i], alpha_to_es)
      if (is.finite(a_raw) && a_raw > 0 && a_raw <= 1 && a_raw > opts$alpha_max) {
        flags[[i]] <- c(flags[[i]],
          paste0("[UNUSUAL] Near-perfect alpha: ", round(a_raw, 3),
                 " (threshold: ", opts$alpha_max, ")", msuf))
      }
    }

    # C8: Near-perfect ICC. Same coefficient-scale comparison as C7.
    if (measure == "icc" && !is.na(es[i]) && is.finite(es[i])) {
      i_raw <- .to_reliability_scale(es[i], icc_to_es)
      if (is.finite(i_raw) && i_raw > 0 && i_raw <= 1 && i_raw > opts$icc_max) {
        flags[[i]] <- c(flags[[i]],
          paste0("[UNUSUAL] Near-perfect ICC: ", round(i_raw, 3),
                 " (threshold: ", opts$icc_max, ")", msuf))
      }
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
                                      n_cases = NULL, n_controls = NULL,
                                      baseline_risk = NULL,
                                      suffix = "", exp = FALSE,
                                      rel_scale = "bonett") {
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
      # dw/gw are standardised SMDs in the same SD units as d/g, and rp/zp are on the
      # same bounded and Fisher scales as r/z, so they inherit the same floors. Falling
      # through to the 0 default (which the comment below reserves for md/nnt/prop,
      # whose units are arbitrary) made D1 flag homogeneous dw/gw/rp/zp pools at
      # deviations the identical code calls non-diagnostic for d/g/r/z. Every sibling
      # gate -- C1, F1/F2, E1/E3's thresholds, Category G's signed_measures -- already
      # groups them this way.
      "d" = , "g" = , "dw" = , "gw" = 1.0, # 1 SD
      "logor" = , "logrr" = , "logirr" = , "loghr" = 1.0, # log scale
      "r" = , "rp" = 0.3,
      "z" = , "zp" = 0.5,
      "rd" = 0.3,
      # Reliability coefficients live on a universal bounded scale, so a floor is
      # meaningful here in a way it is not for md or nnt, whose units are arbitrary.
      # Without one, a 20-study pool whose reported alphas span .91-.93, which is
      # textbook homogeneity, raises 6 "ES outlier" flags. Values are in analysis-scale
      # units: on the Bonett ln(1-x) scale 0.5 is about a .05 swing in the coefficient
      # near .9, and the Hakstian-Whalen scale is roughly 6x tighter, so it gets its
      # own floor. Override per analysis with flag_options$outlier_min_deviation.
      "alpha" = , "icc" = , "omega" = if (identical(rel_scale, "hakstian_whalen")) 0.08 else if (identical(rel_scale, "raw")) 0.05 else 0.5,
      0 # no floor for md, nnt, prop
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
    # Balanced-arms 2/sqrt(N) fallback for rows without both arm sizes. This must run
    # whether or not any row had arm sizes: a pure pearson_r/fisher_z/spearman_r d or g
    # pool carries only n_sample, leaving n_exp and n_nexp NA in raw_data, so nesting
    # this inside if(any(has_both)) would leave every row with its raw, N-dependent SE
    # and falsely flag a legitimately small-N study as an SE outlier. This mirrors the
    # branch-level no_data fallback in the logOR/logRR branch below.
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
    #
    # The balanced-event constant is the SE this measure has at a 50% event rate with
    # equal arms, which is what puts the fallback rows on the reconstruction's ~1
    # scale: 4/sqrt(N) for logOR (sqrt(4 x 4/N)) and 2/sqrt(N) for logRR
    # (sqrt(4/N - 2/N + 4/N - 2/N)). Reusing 4 for RR would leave a whole RR pool at
    # half scale.
    balanced_se_const <- if (identical(measure, "logrr")) 2 else 4
    # The reconstruction is OR-only. .estimate_n_from_or_and_n_exp() solves for a table
    # whose ODDS ratio is or_j and whose Var(logOR) is var_j; handing it a risk ratio
    # and a Var(logRR) returns cells that are not the study's, and the logRR SE
    # evaluated on them is not the study's expected SE -- measured normalised SEs of
    # 0.68 to 14.7 instead of ~1, which both false-flags ordinary low-event-rate trials
    # and inflates the IQR enough to mask real errors. A 2x2 with arm sizes n1, n2 and
    # a given RR has OR = RR (n2 - b)/(n1 - a), which is not identified from
    # (RR, Var(logRR), n1, n2) alone, so there is no correcting the input either; logRR
    # goes to the balanced-event sqrt(N) fallback, which at least keeps one scale.
    use_recon <- !identical(measure, "logrr")
    has_data <- !is.na(n_exp) & !is.na(n_nexp) & is.finite(n_exp) & is.finite(n_nexp) &
                n_exp > 1 & n_nexp > 1 & !is.na(es) & is.finite(es) &
                !is.na(se) & is.finite(se) & se > 0
    if (!use_recon && any(has_data)) {
      # A risk-ratio table is pinned by (RR, n_exp, n_nexp, n_cases) in closed form and
      # WITHOUT the reported variance, so the expected SE below is independent of the
      # SE being tested -- the same property that makes the odds-ratio solve usable:
      #   a/n1 = RR * c/n2 and a + c = n_cases  =>  a = RR*n1*n_cases / (n2 + RR*n1)
      # Expected SE is then Katz: sqrt(1/a - 1/n1 + 1/c - 1/n2).
      #
      # Without a case margin the only option left is the balanced-event sqrt(N)
      # constant, which is blind to the event rate: se*sqrt(N)/2 equals sqrt((1-p)/p)
      # at balanced arms, i.e. 1 at p = 0.5 but 4.36 at p = 0.05, so a pool spanning
      # ordinary event rates spreads over 5x with every table arithmetically exact.
      # Measured on 20 exact tables at RR = 0.8, rates .03 to .50, that false-flagged
      # the 3% row while the identical data at measure = "or" stayed clean.
      for (j in which(has_data)) {
        rr_j <- if (on_exp_scale) es[j] else exp(es[j])
        nm_j <- if (!is.null(n_cases) && length(n_cases) >= j) n_cases[j] else NA_real_
        exp_se <- NA_real_
        if (!is.na(nm_j) && is.finite(nm_j) && nm_j > 0 &&
            is.finite(rr_j) && rr_j > 0) {
          a_j <- rr_j * n_exp[j] * nm_j / (n_nexp[j] + rr_j * n_exp[j])
          c_j <- nm_j - a_j
          if (is.finite(a_j) && is.finite(c_j) && a_j > 0 && c_j > 0 &&
              a_j < n_exp[j] && c_j < n_nexp[j]) {
            v <- 1 / a_j - 1 / n_exp[j] + 1 / c_j - 1 / n_nexp[j]
            if (is.finite(v) && v > 0) exp_se <- sqrt(v)
          }
        }
        se_for_iqr[j] <- if (!is.na(exp_se) && exp_se > 0) {
          se[j] / exp_se
        } else {
          se[j] * sqrt(n_exp[j] + n_nexp[j]) / balanced_se_const
        }
      }
    }
    if (use_recon && any(has_data)) {
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
        # The case margin is what makes this check possible at all. Given only
        # (OR, n_exp, n_nexp) the table is NOT identified, so .estimate_n_from_or_and_n_exp()
        # falls back to searching for the table whose Var(logOR) best matches var_j --
        # and the expected SE read off that table is then, by construction, the reported
        # SE. Measured: every normalised value in a pure-2x2 logOR pool came out 1.0 to
        # ~16 significant digits, so IQR(se_for_iqr) was exactly 0 and the IQR gate below
        # never opened; multiplying one row's SE by 8 moved its normalised value from
        # 1.000 to 1.021, and 0 of 300 corrupted SEs were caught at x2, x3, x5 or x8.
        # You cannot test a number against a quantity derived from that same number.
        #
        # Supplying a second margin lets .solve_2x2_from_or() pin the table from the OR
        # alone -- the function documents that solving "ignores var entirely" -- so the
        # expected SE becomes independent of the reported one and the comparison means
        # something. Only a SOLVED table is used; a searched one is refused in favour of
        # the balanced-event sqrt(N) fallback, which is blind to the event rate but at
        # least is not circular.
        n_cases_j <- if (!is.null(n_cases) && length(n_cases) >= j) n_cases[j] else NA_real_
        n_controls_j <- if (!is.null(n_controls) && length(n_controls) >= j) n_controls[j] else NA_real_
        br_j <- if (!is.null(baseline_risk) && length(baseline_risk) >= j) baseline_risk[j] else NA_real_
        recon <- tryCatch(
          .estimate_n_from_or_and_n_exp(or_j, var_j, n_exp_j, n_nexp_j,
                                        n_cases = n_cases_j, n_controls = n_controls_j,
                                        baseline_risk = br_j),
          error = function(e) NULL
        )
        if (!is.null(recon) && nrow(recon) > 0 && isTRUE(attr(recon, "solved")) &&
            all(is.finite(c(recon$n_cases_exp[1], recon$n_cases_nexp[1],
                            recon$n_controls_exp[1], recon$n_controls_nexp[1])))) {
          a <- recon$n_cases_exp[1]; b <- recon$n_cases_nexp[1]
          c_val <- recon$n_controls_exp[1]; d <- recon$n_controls_nexp[1]
          # logOR SE = sqrt(1/a + 1/b + 1/c + 1/d). Only logOR reaches here: the
          # solver matches an ODDS ratio and a Var(logOR), so its cells mean nothing
          # for a risk-ratio pool (see use_recon above).
          cells <- c(a, b, c_val, d)
          if (any(cells == 0)) cells <- cells + 0.5  # continuity correction
          exp_se <- sqrt(sum(1 / cells))
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

  # D2c: the row D2 and D2b never saw. OPT-IN (se_outlier_missing_n), default FALSE.
  #
  # WHY IT EXISTS. D2/D2b compare an SE only after normalising it by sample size, so a
  # row carrying no usable N is set to NA in se_for_iqr above and leaves the pool
  # SILENTLY -- no flag, no note, nothing separating "checked and fine" from "never
  # checked". MEASURED on the one real odds-ratio pool in the tree
  # (papers/errors_framework/validation_set/42125682-ERROR-OR, Shan et al. 2026,
  # PNI and pCR, 6 studies): its forest row 5 (Wang 2025) is a per-1-unit CONTINUOUS
  # logistic coefficient copied into a pool of high-vs-low binary contrasts, it takes
  # the LARGEST weight in the meta-analysis (22.8%) by virtue of an SE 9x tighter than
  # any sibling, and dropping it moves the pooled OR from 1.90 to 2.21. Every default
  # check is silent on it, and NO THRESHOLD REACHES IT: at se_ratio_extreme = 2 -- a
  # gate loose enough to flag half of any pool -- the row is still silent, because it
  # is not in the pool at all. Supply its N and it becomes an ordinary 1.8% near-miss
  # (normalised ratio 0.1018 against the 0.1000 D2b requires).
  #
  # THE MISSING N IS PART OF THE ERROR SIGNATURE, which is what makes this worth a
  # check rather than an accepted blind spot: that row has no sample size precisely
  # because it is a per-unit coefficient with no group counts -- the same fact that
  # makes it the wrong estimand. The check that should catch it is disabled by the
  # thing that makes it wrong.
  #
  # WHY IT IS SCOPED TO UNNORMALISABLE ROWS rather than being a general raw-SE
  # comparison. Raw SEs are comparable only at a common N; a pool mixing n = 50 with
  # n = 5000 carries a 10x raw spread with nothing wrong in it, which is the whole
  # reason D2 normalises. Confining this to rows D2 already dropped bounds the exposure
  # to rows that are unchecked by construction, so it can ADD a flag but can never
  # change one D2/D2b already raised.
  #
  # DEFAULT OFF, because the trade is a judgement call -- though the obvious
  # alternative lever is measurably worse: lowering se_ratio_extreme from 10 to 9 over
  # the 22 validation-set pool-passes adds two flags, one of them in a SAFE pool
  # (41977094 forest row 8), and still does not reach this row.
  if (isTRUE(opts$se_outlier_missing_n) && length(se_iqr_valid) >= 4) {
    se_raw_pool <- se[!is.na(se_for_iqr) & is.finite(se_for_iqr) &
                      !is.na(se) & is.finite(se) & se > 0]
    if (length(se_raw_pool) >= 4) {
      raw_med <- stats::median(se_raw_pool)
      # D2c carries its OWN fold threshold rather than reusing se_ratio_extreme.
      # se_ratio_extreme governs a NORMALISED ratio; this compares RAW SEs, so the
      # two are not the same quantity and a shared number would be a coincidence.
      # Measured: the motivating row's raw fold is 9.4x, so reusing 10 misses it by
      # 6% -- the same near-miss as D2b, one scale along. Calibration below.
      X <- if (!is.null(opts$se_missing_n_ratio)) opts$se_missing_n_ratio else 5
      if (!is.na(raw_med) && raw_med > 0) {
        for (i in seq_len(n)) {
          # only rows the normalised checks could not see, and only with a usable SE
          if (!is.na(se_for_iqr[i])) next
          if (is.na(se[i]) || !is.finite(se[i]) || se[i] <= 0) next
          if (any(grepl("SE outlier", flags[[i]]))) next
          fold <- raw_med / se[i]
          if (is.finite(fold) && fold > X) {
            msuf <- if (!is.null(info_used)) .method_suffix(info_used[i]) else ""
            flags[[i]] <- c(flags[[i]], paste0(
              "[UNUSUAL] SE outlier, no sample size: SE = ", round(se[i], 3), " is ",
              round(fold, 1), "x tighter than the cohort median (", round(raw_med, 3),
              "), and with no N the normalised SE checks skipped this row", msuf))
          }
        }
      }
    }
  }

  # D3 (V27): spread = pooled arm SD, or the SD implied by the SE; flags SE-as-SD.
  # md/mdw by default, d/g/dw/gw only via sd_outlier_smd (mixed instrument scales)
  #
  # sd_outlier_ratio (K), measured on the validation-set extractions with the branch
  # selection below. EVERY pool has two extractions and they do NOT agree, so the pass
  # has to be named: "forest" is the statistics the meta-analysis published, "team" is
  # the re-extraction from the primary sources. A check meant to catch a REVIEWER's
  # transcription error is calibrated on the forest pass, because that is the artefact
  # the error is in.
  #
  #   pool                 measure  k   forest      team
  #   41065428-ERROR-MD    mdw      12  row 12      (silent)
  #   41946661-ERROR-MDw   mdw      11  (silent)    (silent)
  #   41571219-ERROR-MD    md        5  (silent)    (silent)
  #   41977094-SAFE-MD     md       11  row 8       (silent)
  #   41313891-ERROR-SMD   g        10  rows 4, 6   (silent)   (opt-in)
  #
  # 41065428 forest is the motivating case and it fires: spread 2.4 against a pool
  # median of 14.198, a ratio of 5.9. An earlier version of this comment reported that
  # case as 1.57 against 7.80 (ratio 4.98, under the gate) and concluded that "the
  # default-on branch currently fires on nothing in that pool" -- those are the TEAM
  # numbers, and the conclusion is false for the pass the flag is aimed at. The
  # published calibration of 2.4 / 14.2 was right all along.
  #
  # 41313891 forest still fires under the opt-in after the standardised-scope fix
  # below, because both of its rows transcribe arm SDs; only rows WITHOUT arm SDs are
  # excluded there. 41977094 is a SAFE pool and its forest pass flags one row: that
  # verdict is unchanged by any of this (a crude md row with both arm sizes takes the
  # same branch it always did) and is the specificity cost of K = 5 on a 5.9-ratio
  # gate, not a regression.
  #
  # Reproduce with papers/errors_framework/validation_set/<pool>/extraction_*_{forest,
  # team}.xlsx through convert_df(measure = <ma_measure>, flags = TRUE).
  sd_default <- c("md", "mdw")
  sd_optin   <- c("d", "g", "dw", "gw")
  sd_apply <- !is.null(measure) &&
    (measure %in% sd_default ||
     (isTRUE(opts$sd_outlier_smd) && measure %in% sd_optin))
  if (sd_apply) {
    K_sd <- if (!is.null(opts$sd_outlier_ratio)) opts$sd_outlier_ratio else 5
    is_within <- measure %in% c("mdw", "dw", "gw")
    is_std    <- measure %in% sd_optin
    pooled_sd <- if (!is.null(sd_exp) && !is.null(sd_nexp)) {
      sqrt((sd_exp^2 + sd_nexp^2) / 2)
    } else rep(NA_real_, n)
    # The SE-implied SD has to come out on the SAME scale as the pooled arm SD, or a
    # pool where some rows transcribed their arm SDs and some did not is comparing two
    # quantities. Which inversion is right depends on the DESIGN, so it is selected by
    # the measure rather than by which columns a row happens to carry:
    #
    #   standardised (d/g/dw/gw)  -> there is NO SE-implied SD. SE_g =
    #     sqrt(1/n1 + 1/n2 + g^2/(2N)), so se / sqrt(1/n1+1/n2) reduces to
    #     sqrt(1 + g^2/4): a dimensionless 1.00-1.46 carrying no outcome-scale
    #     information at all. Comparing that against branch (i)'s raw instrument SD
    #     flagged EVERY SE-only row in any pool whose instrument SD exceeded ~5.
    #     Such rows are therefore excluded from the pool, exactly as D2 excludes rows
    #     it cannot normalise the same way as their peers.
    #   within-subject (mdw)      -> se * sqrt(n)          (SE = SD_change/sqrt(n))
    #   two independent groups    -> se / sqrt(1/n1 + 1/n2)
    #   n_total only, two groups  -> se * sqrt(n_total)/2  (balanced arms)
    #
    # Branch (i)'s source is chosen to match, at the call site: the change-score SD
    # for a within-subject measure, the ANCOVA residual SD in the adjusted scope.
    implied_sd <- rep(NA_real_, n)
    if (!is_std) {
      if (is_within) {
        # Paired design: SE = SD_change / sqrt(n), so the inverse is se * sqrt(n).
        # n_nexp is often populated on a within-subject row from a two-arm trial, so
        # the branch is chosen by the MEASURE, not by which columns happen to be
        # filled in -- selecting on availability applied the independent-groups
        # inversion and came out sqrt(2) too small on exactly those rows.
        if (!is.null(n_exp)) {
          ok_one <- !is.na(se) & is.finite(se) &
            !is.na(n_exp) & is.finite(n_exp) & n_exp > 0
          implied_sd[ok_one] <- se[ok_one] * sqrt(n_exp[ok_one])
        }
        if (!is.null(n_total)) {
          ok_tot <- is.na(implied_sd) & !is.na(se) & is.finite(se) &
            !is.na(n_total) & is.finite(n_total) & n_total > 0
          implied_sd[ok_tot] <- se[ok_tot] * sqrt(n_total[ok_tot])
        }
      } else {
        # Two independent groups: SE_MD = SD sqrt(1/n1 + 1/n2).
        if (!is.null(n_exp) && !is.null(n_nexp)) {
          ok_two <- !is.na(se) & is.finite(se) &
            !is.na(n_exp) & is.finite(n_exp) & n_exp > 0 &
            !is.na(n_nexp) & is.finite(n_nexp) & n_nexp > 0
          implied_sd[ok_two] <- se[ok_two] / sqrt(1 / n_exp[ok_two] + 1 / n_nexp[ok_two])
        }
        if (!is.null(n_total)) {
          # balanced-arms fallback, mirroring D2's own 2/sqrt(N) constant
          ok_tot <- is.na(implied_sd) & !is.na(se) & is.finite(se) &
            !is.na(n_total) & is.finite(n_total) & n_total > 0
          implied_sd[ok_tot] <- se[ok_tot] * sqrt(n_total[ok_tot]) / 2
        }
      }
    }
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
  signed_measures <- c("d", "g", "dw", "gw", "md", "mdw", "r", "z", "rp", "zp",
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
    if (majority == "tie") {
      # No majority exists, so neither side can be described as opposing one.
      # Calling both rows of a 1-vs-1 split the outlier would put two contradictory
      # flags on the same pool. Both rows still warrant [UNUSUAL], since an evenly
      # split pool is just as suspicious and there is no basis for preferring
      # either side.
      paste0("[UNUSUAL] Direction conflict: ", core,
             ", and the pool is evenly split (", n_opp,
             " significant row(s) point each way, with no majority direction) - verify the arm/group labels and the extracted value against the primary (likely extraction error or a genuine outlier)")
    } else if (majority != dir) {
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
    dup_id <- .group_label(sid[rows[1]])   # see .row_ref(): "; " in a study_id
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

  # Only apply to standardized metrics where SE scales predictably with N.
  #
  # "rp" and "zp" belong here for the same reason "r" and "z" do: SE(rp) is
  # (1 - rp^2)/sqrt(n - q - 3) and SE(zp) is 1/sqrt(n - q - 3), both as predictable
  # in N as their zero-covariate counterparts. Omitting them makes the whole F family
  # return early on a partial correlation, so F1 and F2 never run on those rows at
  # all. At the fall-through defaults the switches below would otherwise give them
  # (floor 0.1/sqrt(N), ceiling max(2, 8/sqrt(N))):
  #   * F2 misses an SD-entered-as-SE error outright, since an SE of 0.85-0.97 is
  #     caught for "r" (ceiling 0.6) and passes for "rp" (ceiling 2.0);
  #   * F1's flat floor does not shrink with (1 - rp^2), so at N = 200, q = 2 and
  #     rp = 0.95 the floor (0.00707) sits above the true SE (0.00696), and would
  #     false-fire on a valid row if the family ran.
  # Giving rp and zp the r/z treatment rather than the ratio default avoids both.
  standardized <- c("d", "g", "dw", "gw", "r", "z", "rp", "zp",
                    "logor", "logrr", "loghr", "logirr")
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
        "r" = , "z" = , "rp" = , "zp" = 0.3,
        0.1  # logOR, logRR, logHR, logIRR: conservative (event-rate dependent)
      )
      se_floor <- se_floor_mult / sqrt_n
      # r/rp: true SE ~ (1 - r^2)/sqrt(n - 1), scale the floor by (1 - r^2). The
      # partial correlation is on the same bounded scale and its true SE is
      # (1 - rp^2)/sqrt(n - q - 3), i.e. slightly LARGER than r's for the same
      # value, so this floor is conservative for it. z/zp are on the unbounded
      # Fisher scale and get no such scaling.
      if (measure %in% c("r", "rp") && !is.na(es[i]) && is.finite(es[i]) && abs(es[i]) <= 1) {
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
        "r" = , "z" = , "rp" = , "zp" = max(0.6, 4.0 / sqrt_n),
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
                                        r_defaulted = NULL, dispersion_log = NULL) {
  n <- length(n_estimations)
  flags <- vector("list", n)
  for (i in seq_len(n)) flags[[i]] <- character(0)

  diff_max <- if (measure %in% c("d", "g", "dw", "gw", "md", "mdw")) {
    opts$diff_max_smd
  } else if (measure %in% c("r", "rp")) {
    # rp is on the same bounded [-1, 1] scale as r, so the r threshold is the right
    # one; the SMD fallback it used before is 3.33x looser (1.00 vs 0.30), which made
    # E3 effectively unreachable for a partial-correlation pool. "z"/"zp" deliberately
    # stay on the fallback -- z is already there, and moving one of the Fisher-scale
    # pair without the other would be worse than leaving both.
    opts$diff_max_r
  } else if (measure %in% c("logor", "logrr", "logirr", "loghr")) {
    opts$diff_max_logor
  } else {
    opts$diff_max_smd  # fallback
  }

  dispersion_max <- if (measure %in% c("d", "g", "dw", "gw", "md", "mdw")) {
    opts$dispersion_max_smd
  } else if (measure %in% c("r", "rp")) {
    opts$dispersion_max_r
  } else if (measure %in% c("logor", "logrr", "logirr", "loghr")) {
    opts$dispersion_max_logor
  } else {
    opts$dispersion_max_smd  # fallback
  }

  # recompute diff and dispersion on the log scale when exp = TRUE
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
      diff_work <- if (!is.null(min_es) && !is.null(max_es) &&
                       !is.na(min_es[i]) && !is.na(max_es[i]) &&
                       is.finite(min_es[i]) && is.finite(max_es[i]) &&
                       min_es[i] > 0 && max_es[i] > 0) {
        log(max_es[i]) - log(min_es[i])
      } else {
        NA_real_
      }
      # `dispersion` arrives on the NATURAL scale, and the old delta approximation
      # divided it by the hierarchy-SELECTED estimate rather than by the median it was
      # measured around. That is not the quantity the message names, the error factor
      # median/es_selected is unbounded (measured 1.5x on a realistic 2x2 + means_sd
      # row, 3.9x on v = c(3, 30)), it runs in both directions -- a large selected
      # estimate hides a genuine discordance -- and it is not invariant under swapping
      # the exposure and reference codings, which is the one thing the log scale is
      # there to give. Two replacements, in order of exactness:
      #   k = 2: max|v - median(v)| of two values IS half their range, so the log-scale
      #          dispersion is exactly diff_work/2. This also restores the E1 == E3
      #          identity at k = 2 that this file documents.
      #   k > 2: the per-method values are not carried here, so keep a delta form but
      #          centre it on sqrt(min max), the geometric mean, which is
      #          reversal-invariant and independent of which estimate was selected.
      #          (Exactness for k > 2 needs .dispersion_stat(log(v)) computed where the
      #          per-method values live, in .generate_df().)
      disp_work <- if (!is.null(dispersion_log) && length(dispersion_log) >= i &&
                       !is.na(dispersion_log[i]) && is.finite(dispersion_log[i])) {
        # exact: .dispersion_stat(log(v)), computed in .generate_df() where the
        # per-method values live
        dispersion_log[i]
      } else if (!is.na(n_estimations[i]) && n_estimations[i] == 2 &&
                 !is.na(diff_work) && is.finite(diff_work)) {
        # k = 2: max|v - median(v)| of two values IS half their range, so the
        # log-scale dispersion is exactly diff_work/2 without the per-method values.
        diff_work / 2
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
                             omega_to_es = "bonett",
                             prop_to_es = "raw",
                             pre_post_to_smd = "bonett",
                             pool_sd = FALSE,
                             r_defaulted = NULL,
                             smd_denom = NULL,
                             es_order = NULL,
                             z_transform = NULL) {
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
  # Every cross-row check below assumes one row is one comparison, since it compares
  # a study against its peers. Under main_es = FALSE a comparison occupies one row
  # per estimation route, which breaks that invariant and changes the checks' unit of
  # analysis without saying so: a study with k routes enters the IQR pool k times,
  # masking its own outlier status, counts k times toward the direction-conflict
  # floor, and matches its own study_id k times in the duplication check.
  #
  # Every cross-row property is therefore decided on one representative row per
  # (row_id, scope), and the verdict is broadcast back to that comparison's rows.
  # When main_es = TRUE every row is its own representative, cmp_map is the identity,
  # and this step does nothing, so the default path is unchanged.
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
  disp_log_col <- paste0("dispersion_es_log", suffix)
  dispersion_log <- if (disp_log_col %in% colnames(res)) {
    suppressWarnings(as.numeric(as.character(res[[disp_log_col]])))
  } else NULL
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
  # B5's person-time floor needs the rate es_from_cases_time() ACTUALLY used, which is
  # the baseline_rate column when supplied and n_cases_nexp / time_nexp otherwise
  # (R/es_from_stand_IRR.R). Reading only the column left every auto-computed row with
  # no lower bound at all once B3 stopped applying the risk-scale one to it.
  baseline_rate_vec <- if ("baseline_rate" %in% colnames(raw_data)) {
    raw_data$baseline_rate[match(res$row_id, raw_data$row_id)]
  } else {
    NULL
  }
  if (all(c("n_cases_nexp", "time_nexp") %in% colnames(raw_data))) {
    idx_br <- match(res$row_id, raw_data$row_id)
    cases_n <- suppressWarnings(as.numeric(raw_data$n_cases_nexp[idx_br]))
    time_n  <- suppressWarnings(as.numeric(raw_data$time_nexp[idx_br]))
    derived <- ifelse(!is.na(cases_n) & !is.na(time_n) & is.finite(time_n) & time_n > 0,
                      cases_n / time_n, NA_real_)
    if (is.null(baseline_rate_vec)) {
      baseline_rate_vec <- derived
    } else {
      miss_br <- is.na(baseline_rate_vec)
      baseline_rate_vec[miss_br] <- derived[miss_br]
    }
  }

  # covariate count, for the A6 CI-width check: .es_from_d() builds the d/g interval
  # of an ADJUSTED row on qt(.975, N - 2 - q), so A6 must expect the same df
  n_cov_row <- if ("n_cov_ancova" %in% colnames(raw_data)) {
    suppressWarnings(as.numeric(raw_data$n_cov_ancova[match(res$row_id, raw_data$row_id)]))
  } else {
    NULL
  }

  # Same, for the regression covariate count: a partial-correlation ("rp") row has
  # its interval built on qt(.975, n_sample - n_covariates - 2), which is a different
  # column from n_cov_ancova and a different df, so A6 needs it separately.
  n_covariates_row <- if ("n_covariates" %in% colnames(raw_data)) {
    suppressWarnings(as.numeric(raw_data$n_covariates[match(res$row_id, raw_data$row_id)]))
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
                                  n_cov_ancova = n_cov_row,
                                  n_covariates = n_covariates_row)
  f_b <- .flag_bounds_violations(es, se, ci_lo, ci_up, measure, exp, info_used,
                                  baseline_risk = baseline_risk_vec,
                                  baseline_rate = baseline_rate_vec,
                                  alpha_to_es = alpha_to_es,
                                  icc_to_es = icc_to_es,
                                  prop_to_es = prop_to_es)
  f_c <- .flag_plausibility(es, se, measure, opts, n_sample, info_used,
                            n_exp = n_exp_vec, n_nexp = n_nexp_vec,
                            exp = exp, prop_to_es = prop_to_es,
                            alpha_to_es = alpha_to_es, icc_to_es = icc_to_es)
  # D3's branch (i) must be on the same scale as its SE-implied branch (ii), and that
  # scale depends on the measure and on the scope:
  #   within-subject (mdw/dw/gw) -> the CHANGE-score SD, which is what the paired SE
  #     inverts to. The endpoint SD is a different quantity entirely
  #     (SD_change = SD sqrt(2 - 2r), i.e. 0.63x at the package's default r = 0.8), and
  #     feeding it in split one pool across two estimands and false-fired at the
  #     default gate on error-free data.
  #   adjusted scope -> the ANCOVA residual SD, which is what the adjusted SE inverts
  #     to; the crude endpoint SD is a marginal quantity and is larger by 1/sqrt(1-R^2).
  .d3_sd_col <- function(nm) if (nm %in% colnames(raw_data)) {
    suppressWarnings(as.numeric(raw_data[[nm]][match(res$row_id, raw_data$row_id)]))
  } else NULL
  .d3_sd_pick <- function(arm) {
    cands <- if (measure %in% c("mdw", "dw", "gw")) {
      paste0("mean_change_sd_", arm)
    } else if (identical(suffix, "_adjusted")) {
      c(paste0("ancova_mean_sd_", arm), paste0("mean_sd_", arm))
    } else {
      paste0("mean_sd_", arm)
    }
    out <- NULL
    for (nm in cands) {
      v <- .d3_sd_col(nm)
      if (is.null(v)) next
      if (is.null(out)) out <- v else out[is.na(out)] <- v[is.na(out)]
    }
    out
  }
  sd_exp_vec  <- .d3_sd_pick("exp")
  sd_nexp_vec <- .d3_sd_pick("nexp")
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
  # D2's logOR reconstruction needs a margin the reported SE did not supply; see the
  # circularity note in .flag_cross_row_outliers(). n_cases/n_controls are taken from
  # the raw input when present, and derived from the 2x2 cells when they are not.
  .d2_margin <- function(nm, cells) {
    v <- if (nm %in% colnames(raw_data)) {
      suppressWarnings(as.numeric(raw_data[[nm]][match(res$row_id, raw_data$row_id)]))
    } else rep(NA_real_, nrow(res))
    for (cc in cells) {
      if (!cc %in% colnames(raw_data)) return(v)
    }
    parts <- lapply(cells, function(cc)
      suppressWarnings(as.numeric(raw_data[[cc]][match(res$row_id, raw_data$row_id)])))
    derived <- Reduce(`+`, parts)
    v[is.na(v)] <- derived[is.na(v)]
    v
  }
  n_cases_vec <- .d2_margin("n_cases", c("n_cases_exp", "n_cases_nexp"))
  n_controls_vec <- .d2_margin("n_controls", c("n_controls_exp", "n_controls_nexp"))
  baseline_risk_vec <- if ("baseline_risk" %in% colnames(raw_data)) {
    suppressWarnings(as.numeric(raw_data$baseline_risk[match(res$row_id, raw_data$row_id)]))
  } else rep(NA_real_, nrow(res))
  n_casesR <- n_cases_vec[cmp_rep_i]
  n_controlsR <- n_controls_vec[cmp_rep_i]
  baseline_riskR <- baseline_risk_vec[cmp_rep_i]
  group_keyR <- group_key[cmp_rep_i]

  f_d <- .cmp_broadcast(.by_group(group_keyR, nR, function(idx)
    .flag_cross_row_outliers(esR[idx], seR[idx], opts, info_usedR[idx], measure,
                             n_total = n_totalR[idx],
                             n_exp = n_expR[idx], n_nexp = n_nexpR[idx],
                             sd_exp = sd_expR[idx], sd_nexp = sd_nexpR[idx],
                             n_cases = n_casesR[idx], n_controls = n_controlsR[idx],
                             baseline_risk = baseline_riskR[idx],
                             suffix = suffix, exp = exp,
                             rel_scale = switch(measure, "alpha" = alpha_to_es,
                                               "icc" = icc_to_es,
                                               "omega" = omega_to_es, "bonett"))))
  f_g <- .cmp_broadcast(.by_group(group_keyR, nR, function(idx)
    .flag_cross_row_direction_conflict(esR[idx], ci_loR[idx], ci_upR[idx], opts,
                                       info_usedR[idx], measure, exp)))

  # H: cross-row study duplication, reading study_id from raw_data. Always pass a
  # length-n vector so f_dup[[i]] is safe downstream. Scoped within group_key,
  # because sharing a study_id across different groups, as when one trial reports
  # several outcomes, is expected in multivariate data and must not be flagged.
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
                                     r_defaulted = r_def_row,
                                     dispersion_log = dispersion_log)

  # E4: Cross-row NNT type mixing (risk-based vs rate-based), scoped WITHIN group_key
  # (only rows that would actually be pooled together can be "mixed").
  f_nnt_mix <- vector("list", n)
  for (i in seq_len(n)) f_nnt_mix[[i]] <- character(0)
  if (measure %in% c("nnt","rd") && !is.null(info_used)) {
    rate_methods <- .rate_scale_methods()
    # The gate covers both measures; the WORDING has to follow the measure, or an rd
    # pool is told about a mix of "NNT types" it never asked for and never hears the
    # risk-difference / incidence-rate-difference distinction that is the real hazard.
    rd_pool <- identical(measure, "rd")
    head_lbl  <- if (rd_pool) "Mixed difference types" else "Mixed NNT types"
    rate_sing <- if (rd_pool) "an incidence rate difference (events per person-time)" else "a rate-based (person-time) NNT"
    risk_sing <- if (rd_pool) "a risk difference (a probability)" else "a risk-based NNT"
    rate_plur <- if (rd_pool) "incidence rate differences" else "rate-based (person-time) NNTs"
    risk_plur <- if (rd_pool) "risk differences" else "risk-based NNTs"
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
          this_type  <- if (is_rate) rate_sing else risk_sing
          other_type <- if (is_rate) risk_plur else rate_plur
          msuf <- .method_suffix(info_used[i])
          f_nnt_mix[[i]] <- paste0(
            "[DISCORDANT] ", head_lbl, ": this row is ", this_type, msuf,
            ", but other rows are ", other_type,
            " - these have different units and should not be pooled together", gtxt)
        }
      }
    }
  }

  # E6: Cross-row SMD standardizer mixing (change-SD vs raw-score-SD metric)
  #
  # [INFO] rather than [DISCORDANT], since DISCORDANT is reserved for a single row
  # whose several input sources disagree with one another. This is a cross-row
  # property of the pool, and every row may be individually correct.
  #
  # An SMD standardized by the change SD (morris_dz on pre/post or mean-change data)
  # is not on the same scale as an SMD standardized by a raw-score SD: under equal
  # pre/post SDs, sd_change = sd_raw * sqrt(2(1-r)), so the two differ by a factor
  # 1/sqrt(2(1-r)) and coincide only at r = 0.5. The Cochrane Handbook (v6, section
  # 10.5.2) advises against combining change-score SMDs with post-intervention SMDs,
  # "because the SDs used in the standardization reflect different things". The
  # package's own docs say the same of d_rm against d_z.
  #
  # Rows are classified by the standardizer their info_used and pre_post_to_smd
  # imply. Endpoint (means_sd and similar), baseline-SD (bonett) and raw-metric d_rm
  # rows all live on a raw-score SD, so they are pooled into one raw-score-SD class
  # and mixing them raises nothing; d_rm is precisely the transformation onto that
  # metric. Only morris_dz rows sit on the change-SD metric.
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
  # A paired t (or F) statistic identifies each arm's mean_change / sd_change ratio
  # but not the two arms' SD ratio, so the pooled standardizing SD is not
  # recoverable: these routes have to standardize each arm by its own SD and
  # subtract. When the user has opted into a pooled standardizer (pool_sd = TRUE) for
  # the rows that can be pooled, the two constructions differ whenever a study's arm
  # SDs differ. Informational, because the paired-t rows are not wrong; they are
  # simply the best obtainable from the reported statistic. Silent under the default
  # (pool_sd = FALSE), where every row uses the per-arm construction.
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
    # P25: under pool_sd = TRUE the per-arm fallback contradicts the user's explicit
    # pooling request even when no pooled-capable row coexists, so a pool made only
    # of paired t/F rows is disclosed too, with a message that fits each situation.
    # Scoped within group_key, so the coexistence is decided among rows that would
    # actually be pooled.
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

  # E8: Cross-row z-transform mixing (Fisher's z vs a variance-stabilising z)
  #
  # [INFO] rather than [DISCORDANT], and scoped within group_key, on the same
  # reasoning as E6 above: every row may be individually correct, and this is a
  # property of the pool. E6 is the direct precedent and this check mirrors it.
  #
  # measure = "z" is supposed to hold one quantity so the column can be pooled. It
  # does not. The correlation and binary families report z = atanh(r), while the SMD
  # family under smd_to_cor = "viechtbauer", the default, reports a
  # variance-stabilising transform, a different function of the correlation. Measured
  # on identical data at a point-biserial rho = 0.75, that route reports z = 1.0925
  # while atanh() of its own r is 1.7468. A review holding both SMD studies and
  # correlation studies, which is the ordinary case for a z meta-analysis, mixes the
  # two by default with no user choice involved.
  #
  # This is not the estimand difference between viechtbauer (biserial) and
  # lipsey_cooper (point-biserial), which is a documented choice about which
  # correlation to estimate. It is about the column not being atanh() of the r column
  # on every route.
  #
  # Rows are classified from `z_transform`, built in convert_df() by asking each route
  # whether its own z equals atanh() of its own r, rather than from a hardcoded list
  # of route names.
  f_ztrans_mix <- vector("list", n)
  for (i in seq_len(n)) f_ztrans_mix[[i]] <- character(0)
  if (isTRUE(opts$enable_cross_row) && measure == "z" &&
      !is.null(info_used) && !is.null(z_transform) && length(z_transform)) {
    row_trans <- unname(z_transform[as.character(info_used)])
    valid_idx <- which(!is.na(info_used) & nchar(info_used) > 0 & !is.na(es) &
                         !is.na(row_trans))
    multi_grp <- length(unique(group_key)) > 1L
    for (g in unique(group_key[valid_idx])) {
      gsel <- valid_idx[group_key[valid_idx] == g]
      # Decided on the routes that would actually be SELECTED, as E4/E6/E7 do: under
      # the route view one comparison legitimately offers several routes, which is a
      # within-comparison difference (E1/E3), not cross-row mixing.
      tr <- row_trans[gsel[cmp_rep[gsel]]]
      if (length(unique(tr[!is.na(tr)])) < 2L) next
      gtxt <- if (multi_grp) paste0(" (within group '", .group_label(g), "')") else ""
      for (i in gsel) {
        this_t <- row_trans[i]
        if (is.na(this_t)) next
        this_lab <- switch(this_t,
                           fisher = "Fisher's z, atanh(r)",
                           vst = "a variance-stabilising transform, not atanh(r)",
                           "an inconsistent transform")
        other_lab <- if (identical(this_t, "fisher"))
          "a variance-stabilising transform" else "Fisher's z"
        msuf <- .method_suffix(info_used[i])
        f_ztrans_mix[[i]] <- paste0(
          "[INFO] Mixed z transforms: this row is on ", this_lab, msuf,
          ", but other rows are on ", other_lab,
          " - these are different functions of the correlation and should not be ",
          "pooled. Convert through a single family, or set smd_to_cor = ",
          "'lipsey_cooper' to put the SMD rows on Fisher's z", gtxt)
      }
    }
  }

  # E5: min and max ES on opposite sides of the null (0, or 1 when exp = TRUE)
  # suppressed with E1-E3 when all rows show min ~ -max
  f_direction <- vector("list", n)
  for (i in seq_len(n)) f_direction[[i]] <- character(0)
  signed_measures <- c("d", "g", "dw", "gw", "md", "mdw", "r", "z", "rp", "zp",
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
    combined <- c(v_flags, f_a[[i]], f_b[[i]], f_c[[i]], f_d[[i]], f_e[[i]], f_f[[i]], f_g[[i]], f_dup[[i]], f_nnt_mix[[i]], f_std_mix[[i]], f_paired_t_mix[[i]], f_ztrans_mix[[i]], f_direction[[i]])
    if (length(combined) == 0) return("")
    paste(combined, collapse = "; ")
  }, character(1))

  return(res)
}
