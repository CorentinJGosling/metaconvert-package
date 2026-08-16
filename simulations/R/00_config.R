## -----------------------------------------------------------------------------
## Global configuration. No absolute paths anywhere: everything resolves from the
## simulations/ root, which is located by walking up from this file.
## -----------------------------------------------------------------------------

sim_root <- function() {
  ## works whether sourced from simulations/, simulations/R/ or studies/
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  while (!file.exists(file.path(d, "README.md")) ||
         !dir.exists(file.path(d, "studies"))) {
    parent <- dirname(d)
    if (identical(parent, d)) stop("could not locate the simulations/ root")
    d <- parent
  }
  d
}

sim_path <- function(...) file.path(sim_root(), ...)

## Where the package lives. The simulations MUST run against the development
## source, never a stale installed copy -- the whole point is to test what ships.
PKG_ROOT <- normalizePath(file.path(sim_root(), ".."), winslash = "/", mustWork = FALSE)

load_metaconvert <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(PKG_ROOT, quiet = TRUE, export_all = FALSE)
  } else {
    stop("pkgload is required so the simulations run against R/ source, ",
         "not an installed build. install.packages('pkgload')")
  }
  invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Seeding policy
##
## One base seed for the project. Each (study, condition) gets a deterministic
## stream derived from it, so any single condition can be re-run in isolation and
## reproduce bit-for-bit, with or without parallelism. L'Ecuyer-CMRG is the only
## generator in base R with independent substreams.
## -----------------------------------------------------------------------------
SIM_BASE_SEED <- 20260729L

condition_seed <- function(study, condition_index) {
  ## a stable 31-bit integer from the study name + condition index
  h <- sum(utf8ToInt(study) * seq_along(utf8ToInt(study))) * 7919L
  as.integer((SIM_BASE_SEED + h + condition_index * 104729L) %% .Machine$integer.max)
}

## -----------------------------------------------------------------------------
## Defaults
## -----------------------------------------------------------------------------
SIM_DEFAULTS <- list(
  nrep  = 10000L,
  cores = max(1L, parallel::detectCores() - 2L),
  ## Nominal CI level used by every estimator and every coverage computation.
  level = 0.95,
  ## Minimum valid replications below which performance() withholds a statistic
  ## rather than publishing it (roadmap 3.3). Applies PER FAMILY: the SE columns
  ## are keyed on n_valid_se, the interval columns on n_valid_ci. Worst-case
  ## coverage MCSE at 30 is sqrt(0.25/30) = 0.091, so a value computed from fewer
  ## is not a measurement. Measured against the shipped aggregates this withholds
  ## 30 of 19,129 cells (0.16%), all in study 05's `grant` and study 03a's `phi`.
  min_valid = 30L
)

dir_raw <- function(...) sim_path("data", "raw", ...)
dir_agg <- function(...) sim_path("data", "aggregated", ...)
