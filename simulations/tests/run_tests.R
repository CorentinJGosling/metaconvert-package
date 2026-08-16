## =============================================================================
## Test runner for the simulation programme.
##
##   cd simulations
##   Rscript tests/run_tests.R
##
## simulations/ is not an R package -- it is .Rbuildignore'd and never shipped --
## so its code cannot be reached by tests/testthat/ in the package root without
## making package tests depend on non-package files. It gets its own tree instead.
## Uses testthat if available (for the reporter and expect_* messages) and falls
## back to stopifnot() so the suite runs on a bare R.
##
## THIS IS A THIRD SUITE. The package's two (tests/testthat/ and
## tests_save/checked/) do not execute it. Run it explicitly after touching
## anything under simulations/.
## =============================================================================

.sim_root <- normalizePath(
  if (basename(getwd()) == "tests") ".." else ".", winslash = "/", mustWork = TRUE)

if (!file.exists(file.path(.sim_root, "run_all.R")))
  stop("run from simulations/ or simulations/tests/ -- run_all.R not found")

## Source the study and harness files WITHOUT running anything: every file here
## only defines functions at top level.
for (f in sort(list.files(file.path(.sim_root, "R"), pattern = "[.]R$", full.names = TRUE)))
  source(f)
for (f in sort(list.files(file.path(.sim_root, "studies"), pattern = "[.]R$", full.names = TRUE)))
  source(f)

has_testthat <- requireNamespace("testthat", quietly = TRUE)

test_files <- sort(list.files(file.path(.sim_root, "tests"),
                              pattern = "^test-.*[.]R$", full.names = TRUE))
if (!length(test_files)) stop("no test files found in simulations/tests/")

if (has_testthat) {
  res <- testthat::test_dir(file.path(.sim_root, "tests"), reporter = "summary",
                            env = environment(), load_helpers = FALSE, stop_on_failure = FALSE)
  df <- as.data.frame(res)
  cat("\nSIMULATION TESTS: PASS", sum(df$passed), " FAIL", sum(df$failed),
      " ERROR", sum(df$error), " SKIP", sum(df$skipped), "\n")
  if (sum(df$failed) + sum(df$error) > 0) quit(status = 1L)
} else {
  message("testthat not installed -- running the stopifnot() fallbacks only")
  for (f in test_files) source(f)
  cat("\nSIMULATION TESTS: all stopifnot() invariants held\n")
}
