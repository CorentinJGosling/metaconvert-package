## =============================================================================
## PROVENANCE OF THE SHIPPED AGGREGATES
##
## THE DEFECT THIS REPLACES. tests/test-aggregates-fresh.R used to ask "was the code
## changed after the result was last regenerated?" by comparing GIT COMMIT DATES. The
## question is the right one; the proxy breaks in the one case that actually happens:
## a code change that does not change any number.
##
##   1. R/es_from_stand_OR.R is edited (a message() string) -> 4 aggregates go stale.
##   2. Regenerate them -> the CSVs come back BYTE-IDENTICAL (measured: study 06,
##      md5 2c7ba28dc5623ad51a64aeb992b6a497 before and after).
##   3. Byte-identical means git has nothing to commit, so the aggregate's commit date
##      never advances -> the test stays red FOREVER.
##
## The honest remedy left no trace, so the test could not be cleared by doing the work.
## Its own header warns against exactly that outcome ("permanently red and would be
## ignored, which is worse than no test").
##
## WHAT THIS DOES INSTEAD. Every run records the CONTENT HASH of each file that fed it,
## in data/aggregated/PROVENANCE.csv. The test compares recorded hashes with current
## ones. Regenerating rewrites the record even when the CSV itself does not change, so
## there IS something to commit and the flag clears.
##
## It does NOT try to guess whether a change is inert -- that is undecidable without
## re-running. An inert edit still raises the flag; what changes is that the remedy now
## works. The price of clearing a false alarm becomes "re-run that study", which is the
## right price and is minutes for every study except 04 and 05.
##
## Hashes come from tools::md5sum(), which is base R. This is change detection, not
## security, so md5 is the right tool and adds no dependency.
## =============================================================================

## The provenance file itself, and the module you are reading, are EXCLUDED from every
## dependency set. This is the one declared exception in the scheme and it is
## structural, not a convenience: this file only RECORDS what ran, so it cannot change
## a number, and letting it into the dependency set would invalidate all 19 manifests
## every time a comment in it was edited -- reintroducing the exact failure mode the
## module exists to remove. Everything else is derived, never listed.
PROVENANCE_EXCLUDED <- c("07_provenance.R")

PROVENANCE_FILE <- "PROVENANCE.csv"


#' Hash one file's contents
#'
#' @param path file to hash
#' @return an unnamed md5 string, or NA_character_ when the file is missing
.prov_hash <- function(path) {
  if (!length(path) || is.na(path) || !file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}


#' Which simulations/R/ files can actually affect a run
#'
#' A file under simulations/R/ is a dependency only if something OUTSIDE it -- a study,
#' or another harness file -- calls a function it defines. A pure analysis helper
#' (06_sparsity.R, used only by the write-up and by tests) cannot change any aggregate,
#' and treating it as a dependency would mark every aggregate stale the moment it was
#' added. That is roadmap item 4.1's rule, kept verbatim from the test it moved out of.
#'
#' @return character vector of absolute paths
.live_harness_files <- function() {
  hf <- list.files(sim_path("R"), pattern = "[.]R$", full.names = TRUE)
  hf <- hf[!basename(hf) %in% PROVENANCE_EXCLUDED]
  others <- c(list.files(sim_path("studies"), pattern = "[.]R$", full.names = TRUE), hf)
  keep <- vapply(hf, function(f) {
    l <- readLines(f, warn = FALSE)
    defs <- unlist(regmatches(l, regexpr("^[.]?[A-Za-z][A-Za-z0-9_.]*(?= *<- *function)",
                                         l, perl = TRUE)))
    defs <- unique(defs)
    if (!length(defs)) return(TRUE)          # no functions: cannot rule it out
    for (g in setdiff(others, f)) {
      txt <- paste(readLines(g, warn = FALSE), collapse = "\n")
      if (any(vapply(defs, function(d) grepl(d, txt, fixed = TRUE), logical(1)))) return(TRUE)
    }
    FALSE
  }, logical(1))
  unname(hf[keep])
}


#' The study file behind an aggregate basename
#'
#' Study 01 writes 01a/01b, study 03 writes 03a/03b, so the two-digit prefix is the key.
#'
#' @param agg_basename e.g. "04_or_to_rr_nrep1000.csv"
#' @return absolute path, or NA_character_ when no single study matches
.study_file_for <- function(agg_basename) {
  num <- substr(agg_basename, 1, 2)
  hits <- list.files(sim_path("studies"), pattern = paste0("^", num, "_.*[.]R$"),
                     full.names = TRUE)
  if (length(hits) == 1L) hits else NA_character_
}


#' The package files a study exercises, derived from the es_from_*() calls in it
#'
#' Derived, never declared: a hand-written map would rot exactly as the column list in
#' roadmap 1.2 did.
#'
#' @param study_file absolute path to a studies/*.R file
#' @return character vector of absolute paths under the package's R/
.pkg_files_for <- function(study_file) {
  if (!length(study_file) || is.na(study_file)) return(character(0))
  txt <- readLines(study_file, warn = FALSE)
  calls <- unique(unlist(regmatches(txt, gregexpr("es_from_[A-Za-z0-9_]+", txt))))
  if (!length(calls)) return(character(0))
  pkg_root <- normalizePath(sim_path(".."), winslash = "/", mustWork = FALSE)
  r_files <- list.files(file.path(pkg_root, "R"), pattern = "[.]R$", full.names = TRUE)
  defs <- lapply(r_files, function(f) {
    l <- readLines(f, warn = FALSE)
    unlist(regmatches(l, regexpr("^es_from_[A-Za-z0-9_]+(?= *<- *function)", l, perl = TRUE)))
  })
  names(defs) <- r_files
  unname(names(defs)[vapply(defs, function(d) any(calls %in% d), logical(1))])
}


#' Every file an aggregate depends on
#'
#' @param agg_basename e.g. "04_or_to_rr_nrep1000.csv"
#' @return character vector of absolute paths (study + live harness + package files)
prov_dependencies <- function(agg_basename) {
  study <- .study_file_for(agg_basename)
  if (is.na(study)) return(character(0))
  unique(c(study, .live_harness_files(), .pkg_files_for(study)))
}


#' Path a dependency is recorded under
#'
#' Recorded RELATIVE to the package root so the record survives a clone, a different
#' checkout directory and a different machine -- the property the commit-date scheme
#' was reaching for and got by a route that broke.
.prov_relpath <- function(path) {
  pkg_root <- normalizePath(sim_path(".."), winslash = "/", mustWork = FALSE)
  sub(paste0("^", pkg_root, "/"), "",
      normalizePath(path, winslash = "/", mustWork = FALSE))
}


#' Record the provenance of one aggregate
#'
#' Called by run_study() immediately after the CSV is written, and usable directly for
#' the aggregates that are not produced by run_study().
#'
#' @param agg_basename basename of the CSV just written
#' @return invisibly, the rows written for this aggregate
write_provenance <- function(agg_basename) {
  deps <- prov_dependencies(agg_basename)
  if (!length(deps)) {
    warning("no dependencies derived for ", agg_basename, "; provenance not recorded")
    return(invisible(NULL))
  }
  rows <- data.frame(
    aggregate = agg_basename,
    dependency = vapply(deps, .prov_relpath, character(1)),
    md5 = vapply(deps, .prov_hash, character(1)),
    stringsAsFactors = FALSE)
  rows <- rows[order(rows$dependency), ]

  f <- dir_agg(PROVENANCE_FILE)
  if (file.exists(f)) {
    old <- utils::read.csv(f, stringsAsFactors = FALSE)
    old <- old[old$aggregate != agg_basename, , drop = FALSE]
    # Drop rows for aggregates that no longer exist. Without this the record silently
    # accumulates entries for deleted files -- a probe run at a throwaway nrep leaves
    # its rows behind forever, and check_provenance() cannot see them because it
    # iterates over the files present, not over the record. Found by counting: the
    # record claimed 21 aggregates for 19 files.
    present <- basename(list.files(dir_agg(), pattern = "[.]csv$"))
    old <- old[old$aggregate %in% present, , drop = FALSE]
    rows <- rbind(old, rows)
  }
  rows <- rows[order(rows$aggregate, rows$dependency), ]
  utils::write.csv(rows, f, row.names = FALSE)
  invisible(rows[rows$aggregate == agg_basename, , drop = FALSE])
}


#' Compare recorded provenance with the current tree
#'
#' @return a data.frame of mismatches, one row per (aggregate, dependency). Empty when
#'   every shipped aggregate matches the code now in the tree.
check_provenance <- function() {
  f <- dir_agg(PROVENANCE_FILE)
  aggs <- basename(list.files(dir_agg(), pattern = "[.]csv$"))
  aggs <- setdiff(aggs, PROVENANCE_FILE)
  # Aggregates with no study file are legacy .txt-era leftovers and are not tracked.
  aggs <- aggs[!vapply(aggs, function(a) is.na(.study_file_for(a)), logical(1))]

  if (!file.exists(f)) {
    return(data.frame(aggregate = aggs, dependency = NA_character_,
                      issue = "no PROVENANCE.csv at all", stringsAsFactors = FALSE))
  }
  rec <- utils::read.csv(f, stringsAsFactors = FALSE)
  out <- list()

  # Rows describing a file that is no longer here. Reported rather than ignored: an
  # orphan means the record is describing something the tree does not contain, which
  # is exactly the kind of quiet drift this module exists to surface.
  for (a in setdiff(unique(rec$aggregate), c(aggs, PROVENANCE_FILE)))
    out[[length(out) + 1]] <- data.frame(aggregate = a, dependency = NA_character_,
                                         issue = "recorded but the aggregate is gone",
                                         stringsAsFactors = FALSE)

  for (a in aggs) {
    r <- rec[rec$aggregate == a, , drop = FALSE]
    if (!nrow(r)) {
      out[[length(out) + 1]] <- data.frame(aggregate = a, dependency = NA_character_,
                                           issue = "not recorded", stringsAsFactors = FALSE)
      next
    }
    now <- prov_dependencies(a)
    now_rel <- vapply(now, .prov_relpath, character(1))

    # A dependency that appeared since the run: the study now exercises code it did not.
    for (d in setdiff(now_rel, r$dependency))
      out[[length(out) + 1]] <- data.frame(aggregate = a, dependency = d,
                                           issue = "new dependency", stringsAsFactors = FALSE)
    # One that disappeared: the recorded run used code the study no longer touches.
    for (d in setdiff(r$dependency, now_rel))
      out[[length(out) + 1]] <- data.frame(aggregate = a, dependency = d,
                                           issue = "dependency removed", stringsAsFactors = FALSE)
    # And the ones present in both, whose CONTENT moved.
    for (d in intersect(r$dependency, now_rel)) {
      cur <- .prov_hash(file.path(normalizePath(sim_path(".."), winslash = "/",
                                                mustWork = FALSE), d))
      if (!identical(cur, r$md5[r$dependency == d][1]))
        out[[length(out) + 1]] <- data.frame(aggregate = a, dependency = d,
                                             issue = "content changed since the run",
                                             stringsAsFactors = FALSE)
    }
  }
  if (!length(out))
    return(data.frame(aggregate = character(0), dependency = character(0),
                      issue = character(0), stringsAsFactors = FALSE))
  do.call(rbind, out)
}
