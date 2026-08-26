## =============================================================================
## Numbers quoted in README.md must still be the numbers in the CSV.
##
## PROVENANCE.csv ties each aggregate to the code that produced it. It says nothing
## about figures COPIED OUT of that aggregate into prose -- and that is exactly where
## the rot was found: README's study-09 table published bonett at 0.0185 / 0.2571 /
## 1.051 against the shipped 0.0197 / 0.2868 / 1.053, and 2x2_tetrachoric's coverage as
## 0.928 where the CSV said 0.961. The second figure carried a claim reused further
## down -- "the best method on bias is the WORST on coverage" -- which was false on the
## data as shipped. Both had been wrong since the item 3.5 regeneration, unnoticed,
## because nothing connected the paragraph to the file.
##
## NEITHER SIDE OF THIS COMPARISON IS A TRANSCRIPTION. The expected values are
## recomputed from the aggregate by R/08_published_numbers.R; the actual values are
## parsed out of the markdown table. A test that stored "should be 0.0195" would be the
## same copied number in a new place, and would rot the same way.
## =============================================================================

.needs_pkg <- function() if (!METACONVERT_AVAILABLE) skip("metaConvert not loadable")


test_that("every pinned README table matches the aggregate it was computed from", {
  skip_if(!file.exists(file.path(.sim_root, "README.md")), "no README")

  ids <- published_table_ids()          # named: id -> document it lives in
  for (i in seq_along(ids)) {
    id <- names(ids)[i]
    published <- read_pinned_table(id, sim_path(ids[[i]]))
    recomputed <- published_table(id)

    expect_equal(nrow(published), nrow(recomputed),
                 info = paste(id, "- row count"))
    expect_equal(ncol(published), ncol(recomputed),
                 info = paste(id, "- column count"))

    for (r in seq_len(nrow(published))) {
      pub_key <- .pub_clean(published[[1]][r])
      rec_key <- .pub_clean(as.character(recomputed[[1]][r]))
      # the documents label rows "2x2_tetrachoric (full table)", "(reference)",
      # "(default)" -- all prose, none of it part of the method name
      pub_key <- trimws(sub("\\((full table|reference|default)\\)", "", pub_key))
      expect_equal(pub_key, rec_key,
                   info = sprintf("%s - row %d key", id, r))

      for (cc in seq(2, ncol(published))) {
        pub <- .pub_clean(published[[cc]][r])
        rec <- .pub_clean(as.character(recomputed[[cc]][r]))
        # numeric cells compared as numbers so trailing-zero differences do not fail;
        # anything else (ranks, "0.961 (worst 0.939)") compared as text
        pn <- suppressWarnings(as.numeric(pub))
        rn <- suppressWarnings(as.numeric(rec))
        if (!is.na(pn) && !is.na(rn)) {
          expect_equal(pn, rn, tolerance = 1e-9,
                       info = sprintf("%s - row %d (%s), column %d (%s)",
                                      id, r, rec_key, cc, names(published)[cc]))
        } else {
          expect_equal(pub, rec,
                       info = sprintf("%s - row %d (%s), column %d (%s)",
                                      id, r, rec_key, cc, names(published)[cc]))
        }
      }
    }
  }
})


test_that("a drifted README figure is actually detected", {
  # The teeth. Without this, a recipe that silently returned the parsed table would
  # pass on any README at all -- which is how the figures rotted in the first place.
  skip_if(!file.exists(file.path(.sim_root, "README.md")), "no README")
  p <- file.path(.sim_root, "README.md")
  backup <- readLines(p, warn = FALSE)
  on.exit({ writeLines(backup, p, useBytes = TRUE) }, add = TRUE)

  # Find the row by its KEY and perturb whatever number is in it. Locating it by the
  # current value would have put a hardcoded 0.0195 in this file -- the same copied
  # number this whole mechanism exists to eliminate, just moved somewhere new.
  l <- backup
  i <- grep("^\\|\\s*\\*\\*`bonett`\\*\\*\\s*\\|", l)
  expect_equal(length(i), 1L)
  l[i] <- sub("([0-9])\\.([0-9]{3})([0-9])", "\\1.\\2\\3\\3", l[i])   # append a digit
  writeLines(l, p, useBytes = TRUE)

  published <- read_pinned_table("study09a-result")
  recomputed <- published_table("study09a-result")
  row <- which(.pub_clean(published[[1]]) == "bonett")
  expect_equal(length(row), 1L)
  expect_false(isTRUE(all.equal(
    as.numeric(.pub_clean(published[[2]][row])),
    as.numeric(recomputed[[2]][row]))))
})


test_that("the marker parser finds exactly one table per id, and reads its shape", {
  skip_if(!file.exists(file.path(.sim_root, "README.md")), "no README")
  ids <- published_table_ids()
  for (i in seq_along(ids)) {
    t <- read_pinned_table(names(ids)[i], sim_path(ids[[i]]))
    expect_gt(nrow(t), 0L)
    expect_gt(ncol(t), 1L)
    # the ---- separator row must not survive as data
    expect_false(any(grepl("^-+$", .pub_clean(t[[1]]))))
  }
  expect_error(read_pinned_table("no-such-table-id"), "found 0 times")
})


test_that("the recipes read the aggregate, not a stored copy of the numbers", {
  .needs_pkg()
  # If a recipe returned hardcoded values it would pass the comparison above while
  # verifying nothing. Point it at a perturbed copy of the CSV and it must move.
  f <- dir_agg("09a_or_to_cor_CONT_nrep1000.csv")
  backup <- readLines(f, warn = FALSE)
  on.exit({ writeLines(backup, f, useBytes = TRUE) }, add = TRUE)

  before <- published_table("study09a-result")
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d$bias <- d$bias * 2
  utils::write.csv(d, f, row.names = FALSE)
  after <- published_table("study09a-result")

  expect_false(isTRUE(all.equal(before$mean_abs, after$mean_abs)))
})
