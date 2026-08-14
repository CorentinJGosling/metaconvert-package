#!/usr/bin/env Rscript
# Reverse-dependency check WITHOUT revdepcheck.
#
# Why this exists: revdepcheck builds a private library under revdep/library/ and
# installs into it with the usual "unpack to a temp name, then rename into place"
# dance. On Windows that rename fails intermittently with
#
#   file.rename(tmpInstPath, instPath) : ... 'Accès refusé' / 'Access denied'
#
# because a scanner, an indexer or another R session holds a handle on the
# directory for a moment. revdepcheck runs with warnings-as-errors, so a single
# flake aborts the whole run. It also refuses to install metaConvert at all if the
# calling session has metaConvert loaded ("package is in use").
#
# This script does the same job the plain way: install the development build into a
# throwaway library under the system temp directory (outside Documents, so away from
# whatever is holding handles), then R CMD check each reverse dependency against it.
#
# USAGE -- from the R console, in a session that has NOT loaded metaConvert
# (restart R first; devtools::load_all() counts as loading it):
#
#   source("revdep/check-revdeps.R")
#   check_revdeps()                                  # builds the tarball if needed
#   check_revdeps("metaConvert_2.0.1.tar.gz")        # reuse an existing tarball
#
# Or from a shell, using the full path to Rscript (it is usually not on PATH on
# Windows):
#
#   & "C:\Program Files\R\R-4.5.1\bin\x64\Rscript.exe" revdep/check-revdeps.R
#
# Set RSTUDIO_PANDOC first if you are outside RStudio and want vignettes checked:
#   Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

check_revdeps <- function(tarball = NULL, pkg_root = ".") {

if ("metaConvert" %in% loadedNamespaces()) {
  stop("metaConvert is loaded in this session. Restart R and run this again ",
       "without calling devtools::load_all() or library(metaConvert) first.")
}

repos <- getOption("repos")
if (is.null(repos[["CRAN"]]) || repos[["CRAN"]] == "@CRAN@") {
  repos <- c(CRAN = "https://cloud.r-project.org")
}

pkg_root <- normalizePath(pkg_root, winslash = "/")
stopifnot(file.exists(file.path(pkg_root, "DESCRIPTION")))
desc <- read.dcf(file.path(pkg_root, "DESCRIPTION"))
pkg <- desc[1, "Package"]
ver <- desc[1, "Version"]

tarball <- if (!is.null(tarball)) normalizePath(tarball, winslash = "/") else {
  cand <- file.path(pkg_root, sprintf("%s_%s.tar.gz", pkg, ver))
  if (!file.exists(cand)) {
    message("Building ", basename(cand), " ...")
    system2(file.path(R.home("bin"), "R"), c("CMD", "build", shQuote(pkg_root)))
  }
  normalizePath(cand, winslash = "/")
}
message("Development build: ", tarball)

# Throwaway library OUTSIDE the project tree.
lib <- file.path(tempdir(), "revdep-lib")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
src <- file.path(tempdir(), "revdep-src")
dir.create(src, showWarnings = FALSE, recursive = TRUE)

message("Installing ", pkg, " ", ver, " into ", lib)
install.packages(tarball, lib = lib, repos = NULL, type = "source", quiet = TRUE)
if (!nzchar(system.file(package = pkg, lib.loc = lib))) {
  stop("Could not install the development build into the throwaway library.")
}

ap <- available.packages(repos = repos)
revdeps <- tools::package_dependencies(pkg, db = ap, reverse = TRUE,
                                       which = c("Depends", "Imports", "LinkingTo", "Suggests"))[[pkg]]
revdeps <- sort(unique(revdeps))
if (!length(revdeps)) {
  message("No reverse dependencies on CRAN.")
  return(invisible(data.frame()))
}
message("Reverse dependencies (", length(revdeps), "): ", paste(revdeps, collapse = ", "))

# The revdep's own dependencies come from the user library; the throwaway library is
# put FIRST so the development metaConvert wins over any installed copy.
all_libs <- unique(c(lib, .libPaths(), strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep)[[1]]))
all_libs <- all_libs[nzchar(all_libs)]
old_libs <- Sys.getenv("R_LIBS")
old_force <- Sys.getenv("_R_CHECK_FORCE_SUGGESTS_", unset = NA)
Sys.setenv(R_LIBS = paste(all_libs, collapse = .Platform$path.sep))
# A revdep's *suggested* packages missing from this machine is a local gap, not a
# regression in metaConvert. Without this, R CMD check turns it into a hard ERROR.
Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = "false")
on.exit({
  Sys.setenv(R_LIBS = old_libs)
  if (is.na(old_force)) Sys.unsetenv("_R_CHECK_FORCE_SUGGESTS_") else
    Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = old_force)
}, add = TRUE)
message("Library path for the checks:\n  ", paste(all_libs, collapse = "\n  "))

# Preflight: can a CHILD process actually resolve the revdep's hard dependencies from
# that library path?  If not, R CMD check reports "Package required but not available"
# and the run looks like the revdep is broken -- which is exactly how the stale entry
# in revdep/problems.md was produced. Diagnose it here instead of misattributing it.
hard_deps <- function(p) {
  d <- tools::package_dependencies(p, db = ap, which = c("Depends", "Imports", "LinkingTo"))[[p]]
  setdiff(d, c(pkg, "R", rownames(installed.packages(lib.loc = .Library))))
}
resolvable <- function(deps) {
  if (!length(deps)) return(character(0))
  f <- tempfile(fileext = ".R")
  writeLines(sprintf(
    'cat(paste(Filter(function(p) !nzchar(system.file(package = p)), c(%s)), collapse = " "))',
    paste(sprintf('"%s"', deps), collapse = ", ")), f)
  out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(f)), stdout = TRUE)
  trimws(paste(out, collapse = " "))
}

results <- data.frame(package = character(), version = character(),
                      status = character(), stringsAsFactors = FALSE)
for (rd in revdeps) {
  message("\n=== ", rd, " ===")
  got <- try(download.packages(rd, destdir = src, type = "source", repos = repos), silent = TRUE)
  if (inherits(got, "try-error") || !nrow(got)) {
    results <- rbind(results, data.frame(package = rd, version = NA, status = "DOWNLOAD FAILED"))
    next
  }
  tgz <- got[1, 2]
  missing <- resolvable(hard_deps(rd))
  if (nzchar(missing)) {
    message(rd, ": SKIPPED - hard dependencies not installed: ", missing)
    results <- rbind(results, data.frame(
      package = rd, version = sub(".*_(.*)\\.tar\\.gz$", "\\1", basename(tgz)),
      status = paste("DEPS MISSING:", missing)))
    next
  }
  outdir <- file.path(src, paste0(rd, "-check"))
  dir.create(outdir, showWarnings = FALSE)
  system2(file.path(R.home("bin"), "R"),
          c("CMD", "check", "--no-manual", paste0("--library=", shQuote(lib)),
            paste0("--output=", shQuote(outdir)), shQuote(tgz)),
          stdout = NULL, stderr = NULL)
  log <- file.path(outdir, paste0(rd, ".Rcheck"), "00check.log")
  st <- if (file.exists(log)) {
    ln <- grep("^Status:", readLines(log, warn = FALSE), value = TRUE)
    if (length(ln)) sub("^Status: ", "", tail(ln, 1)) else "NO STATUS LINE"
  } else "CHECK DID NOT RUN"
  message(rd, ": ", st)
  results <- rbind(results, data.frame(
    package = rd, version = sub(".*_(.*)\\.tar\\.gz$", "\\1", basename(tgz)), status = st))
}

cat("\n\n#### Reverse dependency results for ", pkg, " ", ver, " ####\n", sep = "")
print(results, row.names = FALSE)
bad <- results[grepl("ERROR|WARNING|FAILED|DID NOT", results$status), ]
skipped <- results[grepl("DEPS MISSING", results$status), ]
cat("\n", if (nrow(bad)) "NEEDS ATTENTION:" else "All clear (notes are the revdeps' own).", "\n", sep = "")
if (nrow(bad)) print(bad, row.names = FALSE)
if (nrow(skipped)) {
  cat("\nSkipped because THIS MACHINE lacks their dependencies -- a local gap, NOT a\n",
      "regression in ", pkg, ". install.packages() the names below and re-run.\n", sep = "")
  print(skipped, row.names = FALSE)
}
cat("\nCheck logs under: ", src, "\n", sep = "")

invisible(results)
}   # end check_revdeps()

# Run automatically when invoked as `Rscript revdep/check-revdeps.R [tarball]`,
# but do nothing on source() so the console user can call check_revdeps() themselves.
if (!interactive() && sys.nframe() == 0L) {
  .args <- commandArgs(trailingOnly = TRUE)
  check_revdeps(if (length(.args)) .args[1] else NULL)
}
