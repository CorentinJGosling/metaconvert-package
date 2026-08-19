## =============================================================================
## Sparsity diagnostic for the binary studies (04, 05).
##
## WHY THIS EXISTS. Both binary studies show `se_ratio` > 1 for EVERY method --
## 1.17 to 1.66 on the full grid, including `transpose`, which merely copies
## logor_se through. A pattern that uniform across unrelated methods points at
## something structural, and the README long carried the +0.5 continuity
## correction as the suspect, undiagnosed.
##
## THE CORRECTION IS INNOCENT, AND THE CAUSAL DIRECTION WAS BACKWARDS. Split the
## grid by how likely a zero cell actually is and the pattern resolves:
##
##   region                              conditions   se_ratio
##   dense  P(any zero cell) < 1e-4         105       0.989 - 1.000
##   sparse everything else                 255       1.298 - 1.410
##
## Where zero cells essentially never occur, every route is calibrated. The
## inflation appears exactly where the tables get sparse, so the driver is
## SPARSITY ITSELF (Woolf's variance bias in log OR / log RR), and the continuity
## correction -- which only ever applies in the sparse region -- SUPPRESSES part of
## it. Removing the correction makes the ratio worse (1.20-1.71 against 1.08-1.32),
## which is the opposite of what "the correction inflates the SE" predicts.
##
## SECOND USE: SEPARATING MISCALIBRATION FROM SELECTION. `metafor_conv2x2` sits at
## 1.62 on the full grid and looks like the worst-calibrated route. It is at 0.992
## in the dense region. Its full-grid figure is a selection artifact: it fails to
## reconstruct a table in up to 82.6% of a cell (mean non-estimability 0.152,
## against 0.000 for every other route) and is scored on the survivors only.
## Any method with a non-trivial `nonest_rate` must be read this way before it is
## called miscalibrated.
##
## The diagnostic is here rather than in a one-off script so the next person can
## re-run it instead of taking the paragraph on trust.
## =============================================================================

#' Probability that at least one cell of the 2x2 is empty
#'
#' Union BOUND (the sum of the four marginal probabilities), not the exact union.
#' The four events are near-disjoint at the probabilities where this matters, and
#' over-stating the probability is the conservative direction: it can only move a
#' condition OUT of the dense region, never into it.
#'
#' @param n total sample size; p_exp proportion allocated to the exposed arm
#' @param br baseline (non-exposed) risk; rr the true risk ratio
#' @return numeric vector, same length as the inputs
zero_cell_prob <- function(n, p_exp, br, rr) {
  n_exp <- round(n * p_exp)
  n_nexp <- n - n_exp
  p1 <- rr * br            # risk in the exposed arm
  p0 <- br
  bad <- !is.finite(p1) | p1 <= 0 | p1 >= 1 | !is.finite(p0) | p0 <= 0 | p0 >= 1
  out <- (1 - p1)^n_exp + p1^n_exp + (1 - p0)^n_nexp + p0^n_nexp
  out[bad] <- NA_real_     # rr*br >= 1 is not a probability; build_grid_04 drops these
  pmin(out, 1)
}

#' Split an aggregate into the dense and sparse regions and summarise a statistic
#'
#' @param agg     an aggregate data.frame (studies 04 / 05 shape: n, p_exp, br, rr)
#' @param stat    the column to summarise, default "se_ratio"
#' @param target  which target rows to use, default "sample"
#' @param threshold  P(any zero cell) below which a condition counts as dense
#' @return one row per method: dense/sparse means and the counts behind them
#'
#' A method whose `nonest_rate` is non-trivial is reported with it, because its
#' summary is conditional on the replications that survived -- see the file header.
sparsity_split <- function(agg, stat = "se_ratio", target = "sample",
                           threshold = 1e-4) {
  need <- c("n", "p_exp", "br", "rr", "method", stat)
  missing_cols <- setdiff(need, names(agg))
  if (length(missing_cols))
    stop("sparsity_split(): this aggregate has no ", paste(missing_cols, collapse = "/"),
         " column, so it is not one of the binary studies (04 / 05).", call. = FALSE)
  if ("target" %in% names(agg) && !is.null(target)) agg <- agg[agg$target == target, ]
  if (!nrow(agg)) stop("sparsity_split(): no rows for target '", target, "'", call. = FALSE)

  pz <- zero_cell_prob(agg$n, agg$p_exp, agg$br, agg$rr)
  dense <- !is.na(pz) & pz < threshold

  out <- lapply(sort(unique(agg$method)), function(m) {
    i <- agg$method == m
    d <- agg[[stat]][i & dense]
    s <- agg[[stat]][i & !dense]
    data.frame(
      method      = m,
      dense_mean  = mean(d, na.rm = TRUE),
      dense_n     = sum(!is.na(d)),
      sparse_mean = mean(s, na.rm = TRUE),
      sparse_n    = sum(!is.na(s)),
      nonest_mean = if ("nonest_rate" %in% names(agg)) mean(agg$nonest_rate[i], na.rm = TRUE) else NA_real_,
      nonest_max  = if ("nonest_rate" %in% names(agg)) max(agg$nonest_rate[i], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "stat") <- stat
  attr(res, "threshold") <- threshold
  attr(res, "n_dense_conditions") <- sum(dense)
  attr(res, "n_sparse_conditions") <- sum(!dense)
  res
}
