# =============================================================================
# Roadmap item 0.1 -- simulations/ is versioned in git but must never ship.
#
# Un-ignoring simulations/ in .gitignore (so the Monte Carlo programme is under
# version control) is a DIFFERENT mechanism from .Rbuildignore (which keeps it
# out of the built package). This file pins that separation, so a future edit to
# one cannot silently change the other: without ^simulations$ in .Rbuildignore
# the tarball would gain ~16 MB of aggregated CSVs and fail CRAN's size limit.
#
# The same guard is applied to the other large development-only directories that
# are versioned but non-shipping, since they share the failure mode.
# =============================================================================

# Locate the package root from the test's working directory (tests/testthat/).
.pkg_root <- function() {
  p <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(p, "DESCRIPTION"))) return(p)
  # devtools::test() can run from the package root directly
  p2 <- normalizePath(".", winslash = "/", mustWork = FALSE)
  while (!file.exists(file.path(p2, "DESCRIPTION")) && dirname(p2) != p2) {
    p2 <- dirname(p2)
  }
  p2
}

test_that(".Rbuildignore excludes every versioned development-only directory", {
  root <- .pkg_root()
  skip_if_not(file.exists(file.path(root, ".Rbuildignore")),
              ".Rbuildignore not present (not a source checkout)")

  patterns <- readLines(file.path(root, ".Rbuildignore"), warn = FALSE)
  patterns <- trimws(patterns)
  patterns <- patterns[nzchar(patterns) & !startsWith(patterns, "#")]

  # Directories that are (or may be) tracked by git but must not reach the tarball.
  non_shipping <- c("simulations", "tests_save", "archive", "web", "data-raw", "papers")

  for (d in non_shipping) {
    # R CMD build applies .Rbuildignore patterns as Perl regexes to the path
    # relative to the package root.
    matched <- any(vapply(patterns,
                          function(p) grepl(p, d, perl = TRUE),
                          logical(1)))
    expect_true(
      matched,
      info = paste0("'", d, "' is not excluded by any .Rbuildignore pattern; ",
                    "it would be included in the built package.")
    )
  }
})

test_that("simulations/ is NOT git-ignored, so the validation programme is versioned", {
  root <- .pkg_root()
  skip_if_not(file.exists(file.path(root, ".gitignore")),
              ".gitignore not present")

  gi <- trimws(readLines(file.path(root, ".gitignore"), warn = FALSE))
  gi <- gi[nzchar(gi) & !startsWith(gi, "#")]

  # A bare `simulations` or `simulations/` entry would take the whole programme
  # back out of version control -- the exact state roadmap item 0.1 fixed.
  expect_false(
    any(gi %in% c("simulations", "simulations/", "/simulations", "/simulations/")),
    info = paste("The root .gitignore excludes simulations/ again.",
                 "The Monte Carlo programme would be untracked and unrecoverable.")
  )
})

test_that("simulations/.gitignore keeps raw output out but aggregated results in", {
  root <- .pkg_root()
  f <- file.path(root, "simulations", ".gitignore")
  skip_if_not(file.exists(f), "simulations/ not present in this checkout")

  gi <- trimws(readLines(f, warn = FALSE))
  gi <- gi[nzchar(gi) & !startsWith(gi, "#")]

  # ~315 MB of per-replication .rds, fully regenerable from the seeds.
  expect_true("data/raw/" %in% gi || "data/raw" %in% gi,
              info = "simulations/data/raw/ must stay untracked (~315 MB, regenerable).")

  # The aggregated CSVs ARE the shipped result and must be versioned.
  expect_false(any(grepl("^!?data/aggregated", gi)),
               info = "simulations/data/aggregated/ must remain tracked.")
})
