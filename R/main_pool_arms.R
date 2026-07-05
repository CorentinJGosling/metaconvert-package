#' Pool or split arms in multi-arm trials
#'
#' @param x a data.frame in metaConvert wide format (as used by \code{\link{convert_df}}).
#' @param study_id character string: column name identifying which rows belong
#'   to the same multi-arm trial.
#' @param pool_side character string: column name whose values indicate which
#'   side to pool (\code{"exp"} or \code{"nexp"}) or \code{NA} for single-comparison rows
#'   that should be left unchanged.
#' @param method either \code{"pool"} (Cochrane-recommended pooling of summary
#'   statistics) or \code{"split"} (divide the shared arm's sample size).
#'   See details.
#' @param verbose logical: whether to print warnings about columns set to NA
#'   and non-identical shared-side values. Default is \code{TRUE}.
#'
#' @details
#' Multi-arm trials create dependency when multiple experimental arms share a
#' common control group (or vice versa). This function resolves the dependency
#' before effect size computation with \code{\link{convert_df}}.
#'
#' With \code{method = "pool"} (recommended), the arms to be pooled are combined
#' into a single composite arm using Cochrane Handbook formulas (Table 23.3.a):
#' \itemize{
#'   \item Sample sizes and event counts: summed
#'   \item Means: weighted average by sample size
#'   \item SDs: Cochrane pooling formula
#'   \item SEs: recomputed from pooled SD and pooled n
#'   \item CIs, medians, paired statistics: set to NA (cannot be pooled)
#' }
#' Returns one row per multi-arm group.
#'
#' With \code{method = "split"}, the shared arm's sample size (and event counts)
#' is divided equally across comparisons. Each comparison retains its own arm
#' data unchanged, so the output has the same number of rows as the input but
#' with adjusted sample sizes on the shared side. Simpler but slightly
#' conservative.
#'
#' @return A data.frame suitable for \code{\link{convert_df}}. The
#'   \code{pool_side} column is removed from the output.
#'
#' @references
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA
#' (editors). Cochrane Handbook for Systematic Reviews of Interventions.
#' Chapter 23, Table 23.3.a.
#'
#' @export pool_arms
#'
#' @md
#'
#' @examples
#' # Two experimental arms vs one shared control
#' dat <- data.frame(
#'   study_id = c("Study1", "Study1", "Study2"),
#'   pool_side = c("exp", "exp", NA),
#'   n_exp = c(30, 25, 40),
#'   n_nexp = c(28, 28, 35),
#'   mean_exp = c(12.5, 14.2, 10.0),
#'   mean_sd_exp = c(3.1, 2.8, 4.0),
#'   mean_nexp = c(10.1, 10.1, 9.5),
#'   mean_sd_nexp = c(3.0, 3.0, 3.8)
#' )
#'
#' # Pool experimental arms into one composite arm
#' pooled <- pool_arms(dat, study_id = "study_id",
#'                     pool_side = "pool_side", method = "pool")
#'
#' # Or split the shared control's sample size
#' split_dat <- pool_arms(dat, study_id = "study_id",
#'                        pool_side = "pool_side", method = "split")
#'
pool_arms <- function(x,
                      study_id,
                      pool_side,
                      method = c("pool", "split"),
                      verbose = TRUE) {

  method <- match.arg(method)

  if (!is.data.frame(x)) {
    stop("'x' must be a data.frame.")
  }
  if (!study_id %in% colnames(x)) {
    stop(paste0("Column '", study_id, "' not found in 'x'."))
  }
  if (!pool_side %in% colnames(x)) {
    stop(paste0("Column '", pool_side, "' not found in 'x'."))
  }

  ps_vals <- x[[pool_side]]
  valid_vals <- ps_vals %in% c("exp", "nexp") | is.na(ps_vals)
  if (!all(valid_vals)) {
    bad <- unique(ps_vals[!valid_vals])
    stop(paste0("Column '", pool_side, "' contains invalid values: ",
                paste(bad, collapse = ", "),
                ". Allowed values are 'exp', 'nexp', or NA."))
  }

  # Check: within each study_id group, all non-NA pool_side values must be the same
  groups <- split(ps_vals, x[[study_id]])
  for (grp_name in names(groups)) {
    vals <- na.omit(groups[[grp_name]])
    if (length(unique(vals)) > 1) {
      stop(paste0("Study '", grp_name,
                  "' has mixed pool_side values ('exp' and 'nexp'). ",
                  "All rows of a multi-arm trial must pool the same side."))
    }
  }

  x$.orig_row <- seq_len(nrow(x))

  is_passthrough <- is.na(x[[pool_side]])
  x_pass <- x[is_passthrough, , drop = FALSE]
  x_multi <- x[!is_passthrough, , drop = FALSE]

  if (nrow(x_multi) == 0) {
    res <- x_pass
    res[[pool_side]] <- NULL
    res$.orig_row <- NULL
    return(res)
  }

  x_split <- split(x_multi, x_multi[[study_id]])

  if (method == "pool") {
    processed <- lapply(x_split, function(grp) {
      side <- na.omit(unique(grp[[pool_side]]))[1]
      if (nrow(grp) == 1) {
        if (verbose) {
          warning(paste0("Study '", grp[[study_id]][1],
                         "' has only 1 row with pool_side = '", side,
                         "'. Passing through unchanged."))
        }
        return(grp)
      }
      .pool_one_group(grp, side = side, pool_side_col = pool_side,
                       verbose = verbose)
    })
  } else {
    processed <- lapply(x_split, function(grp) {
      side <- na.omit(unique(grp[[pool_side]]))[1]
      if (nrow(grp) == 1) {
        if (verbose) {
          warning(paste0("Study '", grp[[study_id]][1],
                         "' has only 1 row with pool_side = '", side,
                         "'. Passing through unchanged."))
        }
        return(grp)
      }
      .split_one_group(grp, side = side, pool_side_col = pool_side,
                        verbose = verbose)
    })
  }

  x_processed <- do.call(rbind, processed)

  all_cols <- union(colnames(x_pass), colnames(x_processed))
  for (col in setdiff(all_cols, colnames(x_pass))) {
    x_pass[[col]] <- NA
  }
  for (col in setdiff(all_cols, colnames(x_processed))) {
    x_processed[[col]] <- NA
  }

  res <- rbind(x_pass[, all_cols, drop = FALSE],
               x_processed[, all_cols, drop = FALSE])

  res <- res[order(res$.orig_row), , drop = FALSE]

  res[[pool_side]] <- NULL
  res$.orig_row <- NULL
  rownames(res) <- seq_len(nrow(res))

  return(res)
}

#########

# Column families to pool for a given side ("exp" or "nexp")
.pool_cols <- function(side) {
  s <- side  # "exp" or "nexp"

  list(
    n = paste0("n_", s),

    means = list(
      list(mean = paste0("mean_", s),        n = paste0("n_", s)),
      list(mean = paste0("mean_pre_", s),     n = paste0("n_", s)),
      list(mean = paste0("mean_change_", s),  n = paste0("n_", s)),
      list(mean = paste0("plot_mean_", s),    n = paste0("n_", s)),
      list(mean = paste0("ancova_mean_", s),  n = paste0("n_", s)),
      list(mean = paste0("plot_ancova_mean_", s), n = paste0("n_", s))
    ),

    sds = list(
      list(sd = paste0("mean_sd_", s),
           mean = paste0("mean_", s),
           n = paste0("n_", s)),
      list(sd = paste0("mean_pre_sd_", s),
           mean = paste0("mean_pre_", s),
           n = paste0("n_", s)),
      list(sd = paste0("mean_change_sd_", s),
           mean = paste0("mean_change_", s),
           n = paste0("n_", s))
    ),

    ses = list(
      list(se = paste0("mean_se_", s),
           sd = paste0("mean_sd_", s),
           n = paste0("n_", s)),
      list(se = paste0("mean_pre_se_", s),
           sd = paste0("mean_pre_sd_", s),
           n = paste0("n_", s)),
      list(se = paste0("mean_change_se_", s),
           sd = paste0("mean_change_sd_", s),
           n = paste0("n_", s)),
      list(se = paste0("ancova_mean_se_", s),
           sd = paste0("ancova_mean_sd_", s),    # note: may not exist directly
           n = paste0("n_", s))
    ),

    cis = c(
      paste0("mean_ci_lo_", s), paste0("mean_ci_up_", s),
      paste0("mean_pre_ci_lo_", s), paste0("mean_pre_ci_up_", s),
      paste0("mean_change_ci_lo_", s), paste0("mean_change_ci_up_", s),
      paste0("ancova_mean_ci_lo_", s), paste0("ancova_mean_ci_up_", s),
      paste0("plot_mean_ci_lo_", s), paste0("plot_mean_ci_up_", s),
      paste0("plot_ancova_mean_ci_lo_", s), paste0("plot_ancova_mean_ci_up_", s)
    ),

    counts = c(
      paste0("n_cases_", s),
      paste0("n_controls_", s)
    ),

    proportions = list(
      list(prop = paste0("prop_cases_", s),
           cases = paste0("n_cases_", s),
           n = paste0("n_", s))
    ),

    time = paste0("time_", s),

    na_cols = c(
      paste0("med_", s), paste0("q1_", s), paste0("q3_", s),
      paste0("min_", s), paste0("max_", s),
      paste0("paired_t_", s), paste0("paired_t_pval_", s),
      paste0("paired_f_", s), paste0("paired_f_pval_", s),
      paste0("r_pre_post_", s),
      paste0("mean_change_pval_", s)
    ),

    # bars: SD derived from bounds, pooled, then reconstructed around the mean
    plot_sds = list(
      list(lo = paste0("plot_mean_sd_lo_", s),
           up = paste0("plot_mean_sd_up_", s),
           mean = paste0("plot_mean_", s),
           n = paste0("n_", s)),
      list(lo = paste0("plot_ancova_mean_sd_lo_", s),
           up = paste0("plot_ancova_mean_sd_up_", s),
           mean = paste0("plot_ancova_mean_", s),
           n = paste0("n_", s))
    ),

    plot_ses = list(
      list(lo = paste0("plot_mean_se_lo_", s),
           up = paste0("plot_mean_se_up_", s),
           mean = paste0("plot_mean_", s),
           sd_lo = paste0("plot_mean_sd_lo_", s),
           sd_up = paste0("plot_mean_sd_up_", s),
           n = paste0("n_", s)),
      list(lo = paste0("plot_ancova_mean_se_lo_", s),
           up = paste0("plot_ancova_mean_se_up_", s),
           mean = paste0("plot_ancova_mean_", s),
           sd_lo = paste0("plot_ancova_mean_sd_lo_", s),
           sd_up = paste0("plot_ancova_mean_sd_up_", s),
           n = paste0("n_", s))
    ),

    # ancova_mean_sd_pooled is shared across arms, handled apart
    ancova_sds = list(
      list(sd = paste0("ancova_mean_sd_", s),
           mean = paste0("ancova_mean_", s),
           n = paste0("n_", s))
    )
  )
}

# Cochrane SD pooling for two arms
.pool_sd_cochrane <- function(n1, sd1, m1, n2, sd2, m2) {
  if (any(is.na(c(n1, sd1, m1, n2, sd2, m2)))) return(NA_real_)
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2 +
          n1 * n2 / (n1 + n2) * (m1^2 + m2^2 - 2 * m1 * m2)) /
         (n1 + n2 - 1))
}

# Iterative Cochrane SD pooling for k arms
.pool_sd_cochrane_multi <- function(ns, sds, means) {
  k <- length(ns)
  if (k == 0) return(NA_real_)
  if (k == 1) return(sds[1])

  n_cum <- ns[1]
  sd_cum <- sds[1]
  m_cum <- means[1]

  for (i in 2:k) {
    sd_new <- .pool_sd_cochrane(n_cum, sd_cum, m_cum, ns[i], sds[i], means[i])
    m_new <- (n_cum * m_cum + ns[i] * means[i]) / (n_cum + ns[i])
    n_cum <- n_cum + ns[i]
    sd_cum <- sd_new
    m_cum <- m_new
  }

  return(sd_cum)
}

# Get column value from group, handling missing columns
.get_col <- function(grp, col) {
  if (col %in% colnames(grp)) grp[[col]] else rep(NA_real_, nrow(grp))
}

# Pool one multi-arm group (method = "pool")
.pool_one_group <- function(grp, side, pool_side_col, verbose = TRUE) {

  other_side <- ifelse(side == "exp", "nexp", "exp")
  cols_info <- .pool_cols(side)

  result <- grp[1, , drop = FALSE]
  na_warned <- character(0)

  n_col <- cols_info$n
  if (n_col %in% colnames(grp)) {
    ns <- grp[[n_col]]
    result[[n_col]] <- sum(ns, na.rm = TRUE)
    if (all(is.na(ns))) result[[n_col]] <- NA_real_
  }

  n_pooled <- if (n_col %in% colnames(result)) result[[n_col]] else NA_real_

  for (m_info in cols_info$means) {
    mc <- m_info$mean
    nc <- m_info$n
    if (!mc %in% colnames(grp)) next
    vals <- grp[[mc]]
    weights <- .get_col(grp, nc)
    if (all(is.na(vals))) {
      result[[mc]] <- NA_real_
    } else if (all(is.na(weights))) {
      result[[mc]] <- mean(vals, na.rm = TRUE)
      if (verbose) warning(paste0("No sample sizes for weighting '", mc,
                                   "'. Using unweighted mean."))
    } else {
      ok <- !is.na(vals) & !is.na(weights)
      result[[mc]] <- sum(vals[ok] * weights[ok]) / sum(weights[ok])
    }
  }

  for (sd_info in c(cols_info$sds, cols_info$ancova_sds)) {
    sdc <- sd_info$sd
    mc <- sd_info$mean
    nc <- sd_info$n
    if (!sdc %in% colnames(grp)) next
    sd_vals <- grp[[sdc]]
    m_vals <- .get_col(grp, mc)
    n_vals <- .get_col(grp, nc)
    ok <- !is.na(sd_vals) & !is.na(m_vals) & !is.na(n_vals)
    if (sum(ok) < 2) {
      result[[sdc]] <- if (sum(ok) == 1) sd_vals[ok] else NA_real_
    } else {
      result[[sdc]] <- .pool_sd_cochrane_multi(n_vals[ok], sd_vals[ok], m_vals[ok])
    }
  }

  for (se_info in cols_info$ses) {
    sec <- se_info$se
    sdc <- se_info$sd
    if (!sec %in% colnames(grp)) next
    sd_pooled <- if (sdc %in% colnames(result)) result[[sdc]] else NA_real_
    if (!is.na(sd_pooled) && !is.na(n_pooled) && n_pooled > 0) {
      result[[sec]] <- sd_pooled / sqrt(n_pooled)
    } else {
      result[[sec]] <- NA_real_
    }
  }

  for (ci_col in cols_info$cis) {
    if (ci_col %in% colnames(grp)) {
      if (verbose && any(!is.na(grp[[ci_col]])) && !ci_col %in% na_warned) {
        na_warned <- c(na_warned, ci_col)
      }
      result[[ci_col]] <- NA_real_
    }
  }

  for (cnt_col in cols_info$counts) {
    if (cnt_col %in% colnames(grp)) {
      vals <- grp[[cnt_col]]
      result[[cnt_col]] <- if (all(is.na(vals))) NA_real_ else sum(vals, na.rm = TRUE)
    }
  }

  for (p_info in cols_info$proportions) {
    pc <- p_info$prop
    cc <- p_info$cases
    nc <- p_info$n
    if (!pc %in% colnames(grp)) next
    cases_pooled <- if (cc %in% colnames(result)) result[[cc]] else NA_real_
    if (!is.na(cases_pooled) && !is.na(n_pooled) && n_pooled > 0) {
      result[[pc]] <- cases_pooled / n_pooled
    } else {
      result[[pc]] <- NA_real_
    }
  }

  tc <- cols_info$time
  if (tc %in% colnames(grp)) {
    vals <- grp[[tc]]
    result[[tc]] <- if (all(is.na(vals))) NA_real_ else sum(vals, na.rm = TRUE)
  }

  # plot SD bars: derive SD, pool, reconstruct bounds
  for (ps_info in cols_info$plot_sds) {
    lo_col <- ps_info$lo
    up_col <- ps_info$up
    m_col <- ps_info$mean
    n_c <- ps_info$n
    if (!lo_col %in% colnames(grp) && !up_col %in% colnames(grp)) next

    # Derive SDs from plot bars: SD = up - mean (or mean - lo)
    m_vals <- .get_col(grp, m_col)
    up_vals <- .get_col(grp, up_col)
    lo_vals <- .get_col(grp, lo_col)
    n_vals <- .get_col(grp, n_c)

    derived_sds <- ifelse(!is.na(up_vals) & !is.na(m_vals),
                          up_vals - m_vals,
                          ifelse(!is.na(m_vals) & !is.na(lo_vals),
                                 m_vals - lo_vals, NA_real_))

    ok <- !is.na(derived_sds) & !is.na(m_vals) & !is.na(n_vals)
    if (sum(ok) >= 2) {
      sd_pooled <- .pool_sd_cochrane_multi(n_vals[ok], derived_sds[ok], m_vals[ok])
      m_pooled <- result[[m_col]]
      if (lo_col %in% colnames(grp)) result[[lo_col]] <- m_pooled - sd_pooled
      if (up_col %in% colnames(grp)) result[[up_col]] <- m_pooled + sd_pooled
    } else {
      if (lo_col %in% colnames(grp)) result[[lo_col]] <- NA_real_
      if (up_col %in% colnames(grp)) result[[up_col]] <- NA_real_
    }
  }

  # plot SE bars: recompute from pooled SD
  for (pse_info in cols_info$plot_ses) {
    lo_col <- pse_info$lo
    up_col <- pse_info$up
    m_col <- pse_info$mean
    sd_lo_col <- pse_info$sd_lo
    sd_up_col <- pse_info$sd_up
    n_c <- pse_info$n
    if (!lo_col %in% colnames(grp) && !up_col %in% colnames(grp)) next

    # SE = SD / sqrt(n), use pooled SD from plot_sd bars if available
    m_pooled <- if (m_col %in% colnames(result)) result[[m_col]] else NA_real_
    sd_up_val <- if (sd_up_col %in% colnames(result)) result[[sd_up_col]] else NA_real_
    sd_pooled <- if (!is.na(sd_up_val) && !is.na(m_pooled)) sd_up_val - m_pooled else NA_real_

    if (!is.na(sd_pooled) && !is.na(n_pooled) && n_pooled > 0) {
      se_pooled <- sd_pooled / sqrt(n_pooled)
      if (lo_col %in% colnames(grp)) result[[lo_col]] <- m_pooled - se_pooled
      if (up_col %in% colnames(grp)) result[[up_col]] <- m_pooled + se_pooled
    } else {
      if (lo_col %in% colnames(grp)) result[[lo_col]] <- NA_real_
      if (up_col %in% colnames(grp)) result[[up_col]] <- NA_real_
    }
  }

  for (na_col in cols_info$na_cols) {
    if (na_col %in% colnames(grp)) {
      if (verbose && any(!is.na(grp[[na_col]])) && !na_col %in% na_warned) {
        na_warned <- c(na_warned, na_col)
      }
      result[[na_col]] <- NA_real_
    }
  }

  if (verbose && length(na_warned) > 0) {
    warning(paste0("Study '", grp[[colnames(grp)[1]]][1],
                   "': the following columns cannot be pooled and were set to NA: ",
                   paste(na_warned, collapse = ", ")))
  }

  # shared side taken from first row; warn if values differ across arms
  other_cols <- grep(paste0("_", other_side, "$"), colnames(grp), value = TRUE)
  for (oc in other_cols) {
    vals <- grp[[oc]]
    if (length(unique(na.omit(vals))) > 1 && verbose) {
      warning(paste0("Shared-side column '", oc,
                     "' differs across arms of study '",
                     grp[[colnames(grp)[1]]][1],
                     "'. Using first row's value."))
    }
  }

  # n_sample follows the pooled n; pooled-SD columns are no longer valid
  n_other <- if (paste0("n_", other_side) %in% colnames(result)) {
    result[[paste0("n_", other_side)]]
  } else {
    NA_real_
  }
  if ("n_sample" %in% colnames(result) && !is.na(n_pooled) && !is.na(n_other)) {
    result[["n_sample"]] <- n_pooled + n_other
  }
  if ("mean_sd_pooled" %in% colnames(result)) {
    result[["mean_sd_pooled"]] <- NA_real_
  }
  if ("ancova_mean_sd_pooled" %in% colnames(result)) {
    result[["ancova_mean_sd_pooled"]] <- NA_real_
  }

  return(result)
}

# Split one multi-arm group (method = "split")
.split_one_group <- function(grp, side, pool_side_col, verbose = TRUE) {

  other_side <- ifelse(side == "exp", "nexp", "exp")
  k <- nrow(grp)

  # Divide the shared (other) side's sample size by k
  n_col <- paste0("n_", other_side)
  if (n_col %in% colnames(grp) && !all(is.na(grp[[n_col]]))) {
    n_shared <- grp[[n_col]][1]
    if (!is.na(n_shared)) {
      n_each <- floor(n_shared / k)
      remainder <- n_shared - n_each * k
      grp[[n_col]] <- n_each
      # Distribute remainder to first rows
      if (remainder > 0) {
        grp[[n_col]][seq_len(remainder)] <- n_each + 1
      }
    }
  }

  # Divide shared side 2x2 counts by k
  for (cnt_col in c(paste0("n_cases_", other_side),
                     paste0("n_controls_", other_side))) {
    if (cnt_col %in% colnames(grp) && !all(is.na(grp[[cnt_col]]))) {
      cnt_shared <- grp[[cnt_col]][1]
      if (!is.na(cnt_shared)) {
        cnt_each <- floor(cnt_shared / k)
        remainder <- cnt_shared - cnt_each * k
        grp[[cnt_col]] <- cnt_each
        if (remainder > 0) {
          grp[[cnt_col]][seq_len(remainder)] <- cnt_each + 1
        }
      }
    }
  }

  if ("n_sample" %in% colnames(grp)) {
    n_exp_vals <- if ("n_exp" %in% colnames(grp)) grp[["n_exp"]] else NA_real_
    n_nexp_vals <- if ("n_nexp" %in% colnames(grp)) grp[["n_nexp"]] else NA_real_
    grp[["n_sample"]] <- ifelse(!is.na(n_exp_vals) & !is.na(n_nexp_vals),
                                 n_exp_vals + n_nexp_vals, NA_real_)
  }

  return(grp)
}
