## =============================================================================
## Roadmap item 3.4 -- run_everything() aborted on run_99.
##
## run_everything() (run_all.R) discovers runners by NAME:
##
##   runners <- ls(envir = .GlobalEnv, pattern = "^run_[0-9]{2}$")
##   for (r in sort(runners)) do.call(r, list(nrep = nrep, cores = cores))
##
## so membership of the registry is decided by a regex, while the calling
## convention is decided by each function's formals -- and nothing connected the
## two. run_99 (a deterministic wiring check, not a Monte Carlo study) took no
## arguments, matched the pattern, and died with "unused arguments (nrep = ...,
## cores = ...)". It sorts last, so no results were lost and the failure was easy
## to keep ignoring -- which is the real cost: a full run that always ends in an
## error teaches whoever runs it that the error at the end is normal.
##
## The fix gives run_99 the ignored signature rather than excluding it from the
## pattern, so the registry keeps ONE contract. This file asserts that contract
## for every member, so a new runner with a mismatched signature fails here in
## seconds instead of at the end of a multi-hour run_everything().
## =============================================================================

## The discovery expression, kept character-identical to run_all.R's so a change
## there without a change here is visible.
.registry_pattern <- "^run_[0-9]{2}$"

.registry <- function() sort(ls(envir = .GlobalEnv, pattern = .registry_pattern))

test_that("run_everything()'s own discovery expression is the one tested here", {
  # If run_all.R's pattern is edited, this test must be edited too -- otherwise
  # the file below would silently check a different set of functions.
  src <- readLines(file.path(.sim_root, "run_all.R"), warn = FALSE)
  hit <- grep('ls\\(envir = \\.GlobalEnv, pattern = "\\^run_\\[0-9\\]\\{2\\}\\$"\\)', src)
  expect_gt(length(hit), 0)
})

test_that("the registry is non-empty and contains the nine studies plus run_99", {
  reg <- .registry()
  expect_gte(length(reg), 10L)
  expect_true(all(sprintf("run_%02d", 1:9) %in% reg))
  expect_true("run_99" %in% reg)
})

test_that("EVERY registered runner accepts (nrep, cores)", {
  # The contract run_everything() relies on. run_99 may ignore them -- it is not a
  # Monte Carlo study -- but it must accept them.
  for (r in .registry()) {
    f <- get(r, envir = .GlobalEnv)
    expect_true(is.function(f), info = r)
    fmls <- names(formals(f))
    expect_true("nrep" %in% fmls, info = paste(r, "must accept nrep; has:",
                                               paste(fmls, collapse = ", ")))
    expect_true("cores" %in% fmls, info = paste(r, "must accept cores; has:",
                                                paste(fmls, collapse = ", ")))
  }
})

test_that("every registered runner survives the exact call run_everything() makes", {
  # Arity only -- do NOT execute the studies (hours of Monte Carlo). match.call()
  # against the function raises the same "unused arguments" error a real call
  # would, without running the body.
  for (r in .registry()) {
    f <- get(r, envir = .GlobalEnv)
    err <- tryCatch({
      match.call(f, as.call(list(as.name(r), nrep = 10, cores = 1)))
      NA_character_
    }, error = function(e) conditionMessage(e))
    expect_true(is.na(err),
                info = paste0("do.call(", r, ", list(nrep=, cores=)) would fail: ", err))
  }
})

test_that("run_99 really does ignore what it is handed", {
  # The signature must be inert, not quietly load-bearing: a runner that accepted
  # nrep and then used it as a replication count would be a different bug.
  fmls <- formals(run_99)
  expect_true(all(c("nrep", "cores") %in% names(fmls)))
  expect_null(fmls$nrep)
  expect_null(fmls$cores)
  body_txt <- paste(deparse(body(run_99)), collapse = " ")
  expect_false(grepl("\\bnrep\\b", body_txt))
  expect_false(grepl("\\bcores\\b", body_txt))
})

test_that("run_99 executes end to end through the run_everything() calling convention", {
  # It is deterministic and fast, so unlike the studies it CAN be run here -- which
  # is the only way to show the abort is really gone rather than merely re-typed.
  skip_if_not(isTRUE(get0("METACONVERT_AVAILABLE", ifnotfound = FALSE)),
              "metaConvert not loadable")
  out <- tryCatch({
    utils::capture.output(suppressMessages(suppressWarnings(
      do.call("run_99", list(nrep = 1000, cores = 1)))))
    NA_character_
  }, error = function(e) conditionMessage(e))
  expect_true(is.na(out), info = paste("run_99 aborted:", out))
})

## --- stopifnot() fallback for a bare R (no testthat) ------------------------
if (!requireNamespace("testthat", quietly = TRUE)) {
  for (.r in sort(ls(envir = .GlobalEnv, pattern = "^run_[0-9]{2}$"))) {
    .f <- get(.r, envir = .GlobalEnv)
    stopifnot(all(c("nrep", "cores") %in% names(formals(.f))))
  }
}
