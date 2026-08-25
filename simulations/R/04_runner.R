## -----------------------------------------------------------------------------
## The run engine.
##
## Contract for a study: supply a grid of conditions and two functions.
##
##   generate(cond, nrep) -> data.frame with nrep rows of SUMMARY STATISTICS
##                           (the aggregate data a meta-analyst would extract),
##                           plus the target columns theta_pop / theta_sample.
##
##   estimate(dat, cond)  -> long data.frame with columns
##                           method, est, se, ci_lo, ci_up (nrep rows per method),
##                           produced by calling metaConvert's es_from_*().
##
## The split matters for speed. generate() is a per-replication loop because each
## replication is a distinct dataset. estimate() must NOT be: metaConvert's
## es_from_*() functions are vectorised over rows, so one call on an nrep-row
## frame replaces nrep calls on one-row frames. Measured on es_from_means_sd,
## 10,000 replications:
##
##     per-replication loop      9.00 s
##     one vectorised call       0.15 s      (60x)
##     the DGP itself            1.09 s      <- the real cost
##
## The ES layer is therefore free relative to data generation, which removes the
## only reason the original scripts had for re-typing the package formulas inline.
## -----------------------------------------------------------------------------

#' Run one simulation study
#'
#' @param name       study name, used for seeds and output filenames
#' @param grid       data.frame of conditions, one row each
#' @param generate   function(cond, nrep) -> summary-statistic frame
#' @param estimate   function(dat, cond)  -> long frame of method estimates
#' @param nrep       replications per condition
#' @param cores      parallel workers (1 = sequential, easier to debug)
#' @param save_raw   keep the per-replication frame on disk. Aggregates cannot be
#'                   re-derived without it, and the old pipeline's raw output lived
#'                   on a drive that is not in the repo, which is why none of its
#'                   published numbers can be regenerated.
#' @param targets named vector mapping a target label to the column holding it.
#'   Defaults to the two every study must record. Studies where a third estimand
#'   is in play (e.g. SMD->r, where Viechtbauer targets the biserial and
#'   Lipsey-Cooper the point-biserial) pass it here so the estimand gap can be
#'   separated from bias instead of being reported as bias.
run_study <- function(name, grid, generate, estimate,
                      nrep = SIM_DEFAULTS$nrep,
                      cores = SIM_DEFAULTS$cores,
                      save_raw = TRUE,
                      verbose = TRUE,
                      targets = c(population = "theta_pop", sample = "theta_sample")) {

  stopifnot(is.data.frame(grid), nrow(grid) > 0)
  cond_cols <- names(grid)
  t0 <- Sys.time()

  run_one <- function(i) {
    cond <- grid[i, , drop = FALSE]
    ## Deterministic per-condition stream: re-running condition i alone
    ## reproduces exactly, sequentially or in parallel.
    set.seed(condition_seed(name, i), kind = "L'Ecuyer-CMRG")
    dat <- generate(cond, nrep)
    est <- estimate(dat, cond)
    ## carry the targets and the condition down to every method row
    n_meth <- nrow(est) / nrow(dat)
    for (tcol in targets) {
      ## A target column already supplied by estimate() is a PER-METHOD target
      ## ("this method's own estimand") and must not be overwritten. This is the
      ## matched-scale benchmark of papers/pre_post_ipd: it isolates computational
      ## error from estimand mismatch, and it is required whenever two routes
      ## return values on different transforms -- e.g. smd_to_cor's z output is
      ## Fisher's z under "lipsey_cooper" but a variance-stabilising transform
      ## under "viechtbauer", so no single column can score both.
      if (tcol %in% names(est)) next
      if (!tcol %in% names(dat))
        stop(sprintf("generate() for '%s' did not return the target column '%s'",
                     name, tcol))
      est[[tcol]] <- rep(dat[[tcol]], times = n_meth)
    }
    cbind(cond[rep(1, nrow(est)), , drop = FALSE], est, row.names = NULL)
  }

  idx <- seq_len(nrow(grid))
  if (cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    ## Exporting only the closures is not enough: a closure serialises its parent
    ## environment by REFERENCE, so a worker receives its own empty .GlobalEnv and
    ## cannot resolve as_method(), performance(), or any study-level helper. The
    ## whole global environment has to go across. all.names = TRUE is required --
    ## helpers named with a leading dot are otherwise skipped.
    parallel::clusterExport(cl, ls(envir = .GlobalEnv, all.names = TRUE),
                            envir = .GlobalEnv)
    parallel::clusterExport(cl, c("grid", "generate", "estimate", "nrep", "name"),
                            envir = environment())
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(
        pkgload::load_all(Sys.getenv("METACONVERT_ROOT"), quiet = TRUE))
    })
    raw_list <- parallel::parLapply(cl, idx, run_one)
  } else {
    raw_list <- lapply(idx, function(i) {
      if (verbose && (i %% 25 == 0 || i == 1))
        message(sprintf("  [%s] condition %d/%d", name, i, nrow(grid)))
      run_one(i)
    })
  }

  raw <- do.call(rbind, raw_list)

  if (save_raw) {
    dir.create(dir_raw(), recursive = TRUE, showWarnings = FALSE)
    f <- dir_raw(sprintf("%s_raw_nrep%d.rds", name, nrep))
    saveRDS(raw, f, compress = "xz")
    if (verbose) message("  raw -> ", f)
  }

  agg <- summarise_raw(raw, cond_cols, targets)
  agg$study <- name
  agg$nrep  <- nrep

  dir.create(dir_agg(), recursive = TRUE, showWarnings = FALSE)
  fa <- dir_agg(sprintf("%s_nrep%d.csv", name, nrep))
  utils::write.csv(agg, fa, row.names = FALSE)

  # Record WHICH CODE produced this file (see R/07_provenance.R). This is what makes a
  # re-run able to clear the freshness flag: the record changes even when the CSV comes
  # back byte-identical, so there is something to commit. Under the old commit-date
  # scheme an inert code change left the test red with no way to clear it.
  write_provenance(basename(fa))

  if (verbose) {
    message(sprintf("  [%s] %d conditions x %d reps in %.1f min -> %s",
                    name, nrow(grid), nrep,
                    as.numeric(difftime(Sys.time(), t0, units = "mins")), fa))
    nz <- agg$nonest_rate > 0
    if (any(nz, na.rm = TRUE))
      message(sprintf("  NOTE: %d cells have non-estimable replications (max %.1f%%)",
                      sum(nz, na.rm = TRUE), 100 * max(agg$nonest_rate, na.rm = TRUE)))
  }

  invisible(list(raw = raw, agg = agg))
}

#' Helper: turn a metaConvert es_from_* result into the long method frame.
#'
#' @param res   data.frame returned by an es_from_*() call (one row per replication)
#' @param measure  which output family to extract, e.g. "r", "z", "g", "logor"
#' @param method   label for this route
as_method <- function(res, measure, method) {
  data.frame(
    method = method,
    est    = res[[measure]],
    se     = res[[paste0(measure, "_se")]],
    ci_lo  = res[[paste0(measure, "_ci_lo")]],
    ci_up  = res[[paste0(measure, "_ci_up")]],
    stringsAsFactors = FALSE
  )
}
