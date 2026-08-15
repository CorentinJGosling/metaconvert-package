## =============================================================================
## Entry point. Run from the simulations/ root:
##
##   source("run_all.R")          # loads the harness, defines run_01 ... run_08
##   run_03(nrep = 1000, cores = 1)   # one study, quick
##   run_everything()                 # the full set at production nrep
##
## Nothing here writes outside simulations/data/. No absolute paths.
## =============================================================================

.here <- normalizePath(dirname(sys.frame(1)$ofile %||% getwd()),
                       winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(.here, "run_all.R"))) .here <- normalizePath(getwd(), winslash = "/")
setwd(.here)

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
for (f in sort(list.files("studies", pattern = "[.]R$", full.names = TRUE))) source(f)

## Workers need to know where the package source lives.
Sys.setenv(METACONVERT_ROOT = PKG_ROOT)

run_everything <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  runners <- ls(envir = .GlobalEnv, pattern = "^run_[0-9]{2}$")
  message(sprintf("running %d studies at nrep = %d on %d core(s)",
                  length(runners), nrep, cores))
  for (r in sort(runners)) {
    message("\n=== ", r, " ===")
    do.call(r, list(nrep = nrep, cores = cores))
  }
}

message("harness loaded. studies available: ",
        paste(sort(ls(envir = .GlobalEnv, pattern = "^run_[0-9]{2}$")), collapse = ", "))
