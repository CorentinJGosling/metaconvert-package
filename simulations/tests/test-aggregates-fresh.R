## =============================================================================
## Roadmap item 3.5 -- a shipped aggregate must not be older than the code that
## produced it.
##
## WHY THIS IS KEYED ON GIT AND NOT ON FILE MTIMES. The roadmap's wording is "assert
## each shipped CSV is newer than the study file and than the package files it
## exercises", and the obvious reading -- file.mtime() -- cannot do that job in a git
## repository: a clone or a checkout stamps every file with the checkout time, so the
## comparison is either trivially true or arbitrarily false depending on the order the
## files happened to be written. It would have been a test that passes everywhere and
## detects nothing.
##
## Last-commit time is stable across clones and is the question actually being asked:
## was the code changed after the result was last regenerated? Uncommitted edits are
## handled explicitly -- a modified dependency makes the committed aggregate stale, and
## a modified aggregate counts as just-regenerated.
##
## THE PACKAGE DEPENDENCIES ARE DERIVED, NOT DECLARED. Each study file is scanned for
## the es_from_*() routes it calls, and each route is traced to the R/ file that
## defines it. A hand-written map would rot exactly as the column list in roadmap 1.2
## did. Depending on the whole of R/ instead was rejected: it would mark every
## aggregate stale on any package commit, so the test would be permanently red and
## would be ignored, which is worse than no test.
## =============================================================================

.pkg_root <- normalizePath(file.path(.sim_root, ".."), winslash = "/", mustWork = FALSE)

.git_ok <- function() {
  if (nzchar(Sys.which("git")) == FALSE) return(FALSE)
  out <- suppressWarnings(system2("git", c("-C", shQuote(.pkg_root), "rev-parse", "--git-dir"),
                                  stdout = TRUE, stderr = FALSE))
  length(out) > 0 && !is.na(out[1])
}

## Last commit epoch touching a path (0 when the path is untracked / never committed).
.last_commit <- function(path) {
  rel <- sub(paste0("^", .pkg_root, "/"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
  out <- suppressWarnings(system2("git", c("-C", shQuote(.pkg_root), "log", "-1",
                                           "--format=%ct", "--", shQuote(rel)),
                                  stdout = TRUE, stderr = FALSE))
  if (!length(out) || is.na(out[1]) || !nzchar(out[1])) return(0)
  as.numeric(out[1])
}

## Does the working tree hold uncommitted changes to this path?
.is_modified <- function(path) {
  rel <- sub(paste0("^", .pkg_root, "/"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
  out <- suppressWarnings(system2("git", c("-C", shQuote(.pkg_root), "status", "--porcelain",
                                           "--", shQuote(rel)),
                                  stdout = TRUE, stderr = FALSE))
  length(out) > 0 && any(nzchar(out))
}

## Which simulations/R/ files can actually affect a run: those defining at least one
## function that is referenced outside the file itself, in studies/ or elsewhere in R/.
.live_harness_files <- function() {
  hf <- list.files(file.path(.sim_root, "R"), pattern = "[.]R$", full.names = TRUE)
  others <- c(list.files(file.path(.sim_root, "studies"), pattern = "[.]R$", full.names = TRUE), hf)
  keep <- vapply(hf, function(f) {
    l <- readLines(f, warn = FALSE)
    defs <- unlist(regmatches(l, regexpr("^[.]?[A-Za-z][A-Za-z0-9_.]*(?= *<- *function)",
                                         l, perl = TRUE)))
    defs <- unique(defs)
    if (!length(defs)) return(TRUE)          # no functions: cannot rule it out
    for (g in setdiff(others, f)) {
      txt <- paste(readLines(g, warn = FALSE), collapse = "
")
      if (any(vapply(defs, function(d) grepl(d, txt, fixed = TRUE), logical(1)))) return(TRUE)
    }
    FALSE
  }, logical(1))
  hf[keep]
}

## aggregate stem -> study file. Study 01 writes 01a/01b, study 03 writes 03a/03b, etc,
## so the numeric prefix is the key.
.study_file_for <- function(agg_basename) {
  num <- substr(agg_basename, 1, 2)   # "04_or_to_rr_nrep1000.csv" -> "04"
  hits <- list.files(file.path(.sim_root, "studies"),
                     pattern = paste0("^", num, "_.*[.]R$"), full.names = TRUE)
  if (length(hits) == 1L) hits else NA_character_
}

## The package files a study exercises, derived from the es_from_*() calls in it.
.pkg_files_for <- function(study_file) {
  txt <- readLines(study_file, warn = FALSE)
  calls <- unique(unlist(regmatches(txt, gregexpr("es_from_[A-Za-z0-9_]+", txt))))
  if (!length(calls)) return(character(0))
  r_files <- list.files(file.path(.pkg_root, "R"), pattern = "[.]R$", full.names = TRUE)
  defs <- lapply(r_files, function(f) {
    l <- readLines(f, warn = FALSE)
    m <- regmatches(l, regexpr("^es_from_[A-Za-z0-9_]+(?= *<- *function)", l, perl = TRUE))
    unlist(m)
  })
  names(defs) <- r_files
  keep <- vapply(defs, function(d) any(calls %in% d), logical(1))
  names(defs)[keep]
}

test_that("every shipped aggregate is at least as new as the code that produced it", {
  skip_if_not(.git_ok(), "not a git checkout")

  aggs <- list.files(file.path(.sim_root, "data", "aggregated"),
                     pattern = "[.]csv$", full.names = TRUE)
  skip_if(length(aggs) == 0, "no aggregates present")

  # The harness files that actually feed a run. DERIVED, not listed: a file under
  # simulations/R/ is a dependency only if something outside it -- a study, or
  # another harness file -- calls a function it defines. A pure analysis helper
  # (06_sparsity.R, used only by the write-up and by tests) cannot change any
  # aggregate, and treating it as a dependency would mark all 17 stale the moment
  # it was added, which is how a freshness test becomes noise and gets ignored.
  harness <- .live_harness_files()

  stale <- character(0)
  for (a in aggs) {
    base <- basename(a)
    study <- .study_file_for(base)
    if (is.na(study)) next                       # legacy .txt-era files, no study
    deps <- c(study, harness, .pkg_files_for(study))

    # A regenerated-but-uncommitted aggregate is fresh by definition.
    if (.is_modified(a)) next

    a_t <- .last_commit(a)
    for (d in deps) {
      d_t <- .last_commit(d)
      # An uncommitted dependency change makes the COMMITTED aggregate stale.
      if (.is_modified(d)) d_t <- as.numeric(Sys.time())
      if (d_t > a_t)
        stale <- c(stale, sprintf("%s is older than %s", base,
                                  sub(paste0("^", .pkg_root, "/"), "", d)))
    }
  }
  expect_equal(length(stale), 0L,
               info = paste0("stale aggregates:\n  ", paste(unique(stale), collapse = "\n  ")))
})

test_that("the dependency derivation actually finds something", {
  # A silent empty dependency set would make the test above vacuous.
  s4 <- .study_file_for("04_or_to_rr_nrep1000.csv")
  expect_true(!is.na(s4))
  expect_match(basename(s4), "^04_")
  pk <- .pkg_files_for(s4)
  expect_gt(length(pk), 0L)
  expect_true(any(grepl("es_from_stand_OR[.]R$", pk)))   # study 04 calls es_from_or_se()

  s5 <- .study_file_for("05_rr_to_or_nrep1000.csv")
  expect_true(any(grepl("es_from_stand_RR[.]R$", .pkg_files_for(s5))))
})

test_that("every aggregate maps to exactly one study file", {
  aggs <- list.files(file.path(.sim_root, "data", "aggregated"),
                     pattern = "_nrep[0-9]+[.]csv$", full.names = TRUE)
  skip_if(length(aggs) == 0, "no aggregates present")
  for (a in aggs)
    expect_false(is.na(.study_file_for(basename(a))),
                 info = paste("no study file found for", basename(a)))
})

test_that("the shipped aggregates carry the columns the harness now emits", {
  # Guards the 3.5 regeneration itself: a file written before item 3.3 lacks the two
  # replication counts, and would silently be read as if its coverage rested on
  # n_valid.
  aggs <- list.files(file.path(.sim_root, "data", "aggregated"),
                     pattern = "_nrep[0-9]+[.]csv$", full.names = TRUE)
  skip_if(length(aggs) == 0, "no aggregates present")
  for (a in aggs) {
    nm <- names(utils::read.csv(a, nrows = 1))
    expect_true(all(c("n_valid", "n_valid_se", "n_valid_ci") %in% nm),
                info = paste(basename(a), "is missing the replication-count columns"))
  }
})
