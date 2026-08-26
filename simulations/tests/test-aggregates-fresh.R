## =============================================================================
## Roadmap item 3.5 -- a shipped aggregate must not be older than the code that
## produced it. REWRITTEN: the question is unchanged, the instrument is not.
##
## WHY THE COMMIT-DATE VERSION HAD TO GO. It compared the last-commit time of each
## aggregate with that of its dependencies. Sound in principle, and it deliberately
## avoided file.mtime() (a checkout stamps every file, so mtimes answer nothing in a
## repository). But it broke in the case that actually occurs -- a code change that
## changes no number:
##
##   1. R/es_from_stand_OR.R is edited (only a message() string) -> 4 aggregates stale.
##   2. Regenerate them -> the CSVs come back BYTE-IDENTICAL. Measured, not assumed:
##      study 06 returned md5 2c7ba28dc5623ad51a64aeb992b6a497 before and after.
##   3. Byte-identical means git has nothing to commit, so the aggregate's commit date
##      never advances -> the test can never go green again.
##
## The only honest remedy left no trace. A test whose failure cannot be cleared by
## doing the work is a test that gets ignored -- which its own header warned against.
##
## THE INSTRUMENT NOW. Every run records the md5 of each file that fed it, in
## data/aggregated/PROVENANCE.csv (see R/07_provenance.R). This test compares recorded
## hashes with current ones. Regeneration rewrites the record even when the CSV does
## not move, so the flag clears.
##
## What this does NOT do is guess whether a change was inert -- that is undecidable
## without re-running. An inert edit still raises the flag. What changed is that the
## remedy works.
## =============================================================================

.pkg_root <- normalizePath(file.path(.sim_root, ".."), winslash = "/", mustWork = FALSE)


test_that("every shipped aggregate records the code that produced it, and that code has not moved", {
  aggs <- list.files(file.path(.sim_root, "data", "aggregated"),
                     pattern = "_nrep[0-9]+[.]csv$|_reps[0-9]+[.]csv$")
  skip_if(length(aggs) == 0, "no aggregates present")

  bad <- check_provenance()
  msg <- if (nrow(bad)) paste0(
    "\n  ", paste(sprintf("%-42s %-38s %s", bad$aggregate, bad$dependency, bad$issue),
                  collapse = "\n  "),
    "\n\n  Re-run the affected study to refresh both the aggregate and its record:",
    "\n    cd simulations && Rscript -e \"source('run_all.R'); run_04(nrep = 1000, cores = 26)\"") else ""
  expect_equal(nrow(bad), 0L, info = paste0("provenance mismatches:", msg))
})


test_that("the provenance record is not vacuous", {
  # A record that listed nothing would make the test above pass on any tree at all.
  f <- file.path(.sim_root, "data", "aggregated", PROVENANCE_FILE)
  skip_if(!file.exists(f), "no PROVENANCE.csv yet")
  rec <- utils::read.csv(f, stringsAsFactors = FALSE)
  expect_true(all(c("aggregate", "dependency", "md5") %in% names(rec)))
  expect_gt(nrow(rec), 0L)

  # every aggregate must record its own study file, the harness, and package code
  for (a in unique(rec$aggregate)) {
    d <- rec$dependency[rec$aggregate == a]
    expect_true(any(grepl("^simulations/studies/", d)), info = paste(a, "study file"))
    expect_true(any(grepl("^simulations/R/", d)),       info = paste(a, "harness"))
    expect_true(any(grepl("^R/", d)),                   info = paste(a, "package code"))
  }
  # hashes must be real md5s, not NA from a missing file
  expect_true(all(grepl("^[0-9a-f]{32}$", rec$md5)))
})


test_that("a changed dependency is actually detected", {
  # The teeth of the scheme. Without this, a check_provenance() that silently returned
  # zero rows would look identical to a clean tree -- exactly how the old test failed.
  f <- file.path(.sim_root, "data", "aggregated", PROVENANCE_FILE)
  skip_if(!file.exists(f), "no PROVENANCE.csv yet")
  rec <- utils::read.csv(f, stringsAsFactors = FALSE)
  skip_if(nrow(rec) == 0, "empty record")

  # Corrupt one recorded hash in a COPY of the record, point the module at it, and
  # confirm it is reported. The real file is restored on exit whatever happens.
  # raw bytes: writeLines() re-encodes line endings (see the note in
  # test-published-numbers.R -- a test must not mutate the tree it inspects)
  backup <- readBin(f, "raw", file.info(f)$size)
  on.exit({ writeBin(backup, f) }, add = TRUE)
  rec$md5[1] <- paste0(rep("0", 32), collapse = "")
  utils::write.csv(rec, f, row.names = FALSE)

  bad <- check_provenance()
  expect_gt(nrow(bad), 0L)
  expect_true(any(bad$issue == "content changed since the run"))
})


test_that("a record describing a file that is gone is reported, not ignored", {
  # FOUND BY COUNTING, not by design: after some throwaway probe runs the record held
  # 21 aggregates for 19 files. check_provenance() walks the files PRESENT, so orphan
  # rows were invisible to it -- the record could drift from the tree and still report
  # a clean bill. write_provenance() now prunes them and check_provenance() names them.
  f <- file.path(.sim_root, "data", "aggregated", PROVENANCE_FILE)
  skip_if(!file.exists(f), "no PROVENANCE.csv yet")
  # raw bytes: writeLines() re-encodes line endings (see the note in
  # test-published-numbers.R -- a test must not mutate the tree it inspects)
  backup <- readBin(f, "raw", file.info(f)$size)
  on.exit({ writeBin(backup, f) }, add = TRUE)

  rec <- utils::read.csv(f, stringsAsFactors = FALSE)
  ghost <- rec[1, , drop = FALSE]
  ghost$aggregate <- "99z_a_study_that_does_not_exist_nrep1.csv"
  utils::write.csv(rbind(rec, ghost), f, row.names = FALSE)

  bad <- check_provenance()
  expect_true(any(bad$issue == "recorded but the aggregate is gone"))
  expect_true(any(bad$aggregate == "99z_a_study_that_does_not_exist_nrep1.csv"))
})


test_that("the record describes exactly the aggregates that are present", {
  f <- file.path(.sim_root, "data", "aggregated", PROVENANCE_FILE)
  skip_if(!file.exists(f), "no PROVENANCE.csv yet")
  rec <- utils::read.csv(f, stringsAsFactors = FALSE)
  present <- setdiff(basename(list.files(file.path(.sim_root, "data", "aggregated"),
                                         pattern = "[.]csv$")), PROVENANCE_FILE)
  present <- present[!vapply(present, function(a) is.na(.study_file_for(a)), logical(1))]
  expect_setequal(unique(rec$aggregate), present)
})


test_that("the provenance-writing module is excluded from every dependency set", {
  # Structural, and the one declared exception in the scheme. 07_provenance.R only
  # RECORDS what ran, so it cannot change a number; letting it in would invalidate
  # every record whenever a comment in it was edited -- reintroducing the exact
  # failure mode the module exists to remove.
  expect_true("07_provenance.R" %in% PROVENANCE_EXCLUDED)
  expect_false(any(grepl("07_provenance[.]R$", .live_harness_files())))
  expect_false(any(grepl("07_provenance[.]R$", prov_dependencies("04_or_to_rr_nrep1000.csv"))))
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
