## -----------------------------------------------------------------------------
## ADEMP performance measures (Morris, White & Crowther 2019, Stat Med 38:2074).
##
## Every measure carries its Monte Carlo standard error. The old pipeline reported
## three numbers per cell with no MCSE and no n_sim, which made it impossible to
## tell a real difference from sampling noise in the simulation itself.
##
## THE TWO TARGETS
## ---------------
## Each replication records two reference values:
##
##   theta_pop     the population parameter the data were generated from.
##   theta_sample  the same statistic recomputed on that replication's own
##                 simulated sample, i.e. what the original authors would have
##                 reported had they analysed the other metric on the same people.
##
## Bias against theta_sample answers "how faithful is the conversion?" and is the
## meta-analyst's question. Bias against theta_pop answers "does the converted
## value recover the truth?" and is the only target against which COVERAGE and
## VARIANCE CALIBRATION mean anything -- a nominal 95% interval is not built to
## cover a random quantity, so coverage of theta_sample has no 0.95 reference
## point and must never be ranked on |x - 0.95|.
##
## Both are reported for every estimator. Where they disagree, the gap IS the
## estimand mismatch and should be reported as such rather than called bias.
##
## RMSE AND EmpSE DO NOT COMPOSE ON A RANDOM TARGET
## ------------------------------------------------
## `emp_se` below is sd(est) -- the dispersion of the ESTIMATOR, which is the
## ADEMP definition. `rmse` is sqrt(mean((est - target)^2)). The familiar
##
##     RMSE^2 = bias^2 + EmpSE^2
##
## holds only when `target` is a fixed constant. Against theta_sample it expands
## to bias^2 + var(est) + var(target) - 2*cov(est, target), and because
## theta_sample is recomputed on the SAME simulated participants as the estimate,
## that covariance is large: the deviation varies much less than the estimate.
## RMSE can therefore sit well BELOW EmpSE with nothing wrong. Measured on study
## 01a (n = 25, rho = 0, p_exp = 0.5): cor(est, target) = 0.80, sd(est) = 0.2547,
## sd(est - target) = 0.1522, and rmse = 0.1522. Against the fixed theta_pop in
## the same cell the identity is restored (0.2550 vs 0.2551).
##
## So: compare rmse only with other rmse on the SAME target, and do not
## reconstruct it from bias and emp_se unless the target is theta_pop.
## -----------------------------------------------------------------------------

#' Performance measures for one (condition x method) cell
#'
#' @param est,se,ci_lo,ci_up  per-replication estimates and interval
#' @param target              the reference value (length 1 or length(est))
#' @return one-row data.frame
performance <- function(est, se, ci_lo, ci_up, target) {
  n_total <- length(est)
  if (length(target) == 1L) target <- rep(target, n_total)

  ok <- is.finite(est) & is.finite(target)
  n_valid <- sum(ok)

  ## Non-estimability is a first-class result, not something to na.rm away.
  ## Grant's OR->RR transform, for instance, is undefined whenever
  ## rr * br_guess >= 1 -- exactly in the misspecification these studies exist
  ## to measure. Silently dropping those replications evaluates a method
  ## conditional on it not having failed.
  nonest <- 1 - n_valid / n_total

  if (n_valid < 2L) {
    return(data.frame(
      n_sim = n_total, n_valid = n_valid, nonest_rate = nonest,
      bias = NA_real_, bias_mcse = NA_real_,
      emp_se = NA_real_, emp_se_mcse = NA_real_,
      mod_se = NA_real_, se_ratio = NA_real_, se_ratio_mcse = NA_real_,
      rmse = NA_real_, coverage = NA_real_, coverage_mcse = NA_real_,
      ci_width = NA_real_
    ))
  }

  e <- est[ok]; t <- target[ok]
  dev <- e - t

  bias      <- mean(dev)
  emp_se    <- stats::sd(e)
  bias_mcse <- emp_se / sqrt(n_valid)

  ## MCSE of the empirical SE (Morris et al. eq. 8)
  emp_se_mcse <- emp_se / sqrt(2 * (n_valid - 1))

  ## Model-based SE: what the estimator claims about itself.
  ok_se  <- ok & is.finite(se) & se >= 0
  mod_se <- if (sum(ok_se) > 1L) sqrt(mean(se[ok_se]^2)) else NA_real_
  se_ratio <- if (is.finite(mod_se) && emp_se > 0) mod_se / emp_se else NA_real_
  ## Relative-error MCSE for the ratio, dominated by the EmpSE term.
  se_ratio_mcse <- if (is.finite(se_ratio)) se_ratio / sqrt(2 * (n_valid - 1)) else NA_real_

  rmse <- sqrt(mean(dev^2))

  ok_ci <- ok & is.finite(ci_lo) & is.finite(ci_up)
  if (sum(ok_ci) > 1L) {
    cov_i    <- ci_lo[ok_ci] <= target[ok_ci] & ci_up[ok_ci] >= target[ok_ci]
    coverage <- mean(cov_i)
    coverage_mcse <- sqrt(coverage * (1 - coverage) / sum(ok_ci))
    ci_width <- mean(ci_up[ok_ci] - ci_lo[ok_ci])
  } else {
    coverage <- coverage_mcse <- ci_width <- NA_real_
  }

  data.frame(
    n_sim = n_total, n_valid = n_valid, nonest_rate = nonest,
    bias = bias, bias_mcse = bias_mcse,
    emp_se = emp_se, emp_se_mcse = emp_se_mcse,
    mod_se = mod_se, se_ratio = se_ratio, se_ratio_mcse = se_ratio_mcse,
    rmse = rmse, coverage = coverage, coverage_mcse = coverage_mcse,
    ci_width = ci_width
  )
}

#' Summarise a raw per-replication frame into one row per condition x method x target
#'
#' @param raw   long frame with columns: method, est, se, ci_lo, ci_up,
#'              plus the target columns and the condition columns
#' @param cond_cols  names of the condition (factor) columns
#' @param targets    named vector of target columns, e.g.
#'                   c(population = "theta_pop", sample = "theta_sample")
summarise_raw <- function(raw, cond_cols, targets) {
  out <- list()
  grp <- do.call(interaction,
                 c(unname(as.list(raw[cond_cols])), list(raw$method),
                   list(drop = TRUE, sep = "\r")))
  for (tgt_name in names(targets)) {
    tgt_col <- targets[[tgt_name]]
    if (!tgt_col %in% names(raw)) next
    pieces <- lapply(split(raw, grp), function(d) {
      cbind(
        d[1, cond_cols, drop = FALSE],
        method = d$method[1],
        target = tgt_name,
        performance(d$est, d$se, d$ci_lo, d$ci_up, d[[tgt_col]]),
        row.names = NULL
      )
    })
    out[[tgt_name]] <- do.call(rbind, pieces)
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' Flag cells whose difference between two methods is within Monte Carlo noise.
#' Guards against reading a "winner" off a difference the simulation cannot resolve.
mc_resolvable <- function(bias_a, mcse_a, bias_b, mcse_b, k = 2) {
  abs(bias_a - bias_b) > k * sqrt(mcse_a^2 + mcse_b^2)
}
