## =============================================================================
## Roadmap item 4.4 -- README.md's pre/post kernel table said the wrong thing, and
## said it in prose, which is why nobody noticed it go stale.
##
## The claim was "study 08 exercises one of three kernels", with
## .single_group_pre_post_to_smd marked "never run". It is run twice on every row:
## the two-group entry point is a DISPATCHER, and under pool_sd = FALSE (its
## default, and study 08's setting) it calls the single-group kernel once per arm
## and combines. A second thing the table missed entirely: the paired-t/F family
## called NEITHER kernel -- it reimplemented the arithmetic inline, a fourth pre/post
## code path.
##
## THAT FOURTH PATH NO LONGER EXISTS. The five paired-t/F routes now delegate through
## .paired_t_to_smd(), so there are three kernels again and the duplication is gone.
## The change surfaced HERE FIRST: the block below asserted "calls neither kernel", it
## failed, and that is precisely the job item 4.4 gave it. Making the claim executable
## meant restructuring the dispatch could not silently invalidate the paragraph
## describing it -- which is how the original claim had rotted.
##
## These tests install call counters on the two kernels and assert which ones each
## route actually reaches.
##
## They assert DISPATCH, not arithmetic. The numbers each path produces are pinned
## by the package suites (tests_save/checked/test-pooled-variance-calibration.R for
## the pooled kernel, test-paired-cross-validation.R and
## test-INTERNAL-FULL-LIFECYCLE.R for the paired routes, plus
## tests/testthat/test-paired-t-delegation.R for what the delegation must preserve).
## =============================================================================

.needs_pkg <- function() if (!METACONVERT_AVAILABLE) skip("metaConvert not loadable")

## Counting wrapper. trace() is the obvious tool and is deliberately NOT used: it
## edits the function in place and leaves the namespace altered if a test errors
## before untrace(). Rebinding through a local wrapper and restoring on.exit()
## cannot leak, and it counts calls the same way.
.count_kernel_calls <- function(expr) {
  ns <- asNamespace("metaConvert")
  nm <- c(single_group = ".single_group_pre_post_to_smd",
          pooled       = ".pooled_pre_post_to_smd")
  orig <- lapply(nm, get, envir = ns)
  ## pkgload::load_all() may leave namespace bindings unlocked, so restore whatever
  ## state was actually there rather than assuming they were locked.
  was_locked <- vapply(nm, bindingIsLocked, logical(1), env = ns)
  n <- c(single_group = 0L, pooled = 0L)

  on.exit({
    for (k in names(nm)) {
      if (bindingIsLocked(nm[[k]], ns)) unlockBinding(nm[[k]], ns)
      assign(nm[[k]], orig[[k]], envir = ns)
      if (was_locked[[k]]) lockBinding(nm[[k]], ns)
    }
  }, add = TRUE)

  for (k in names(nm)) {
    if (was_locked[[k]]) unlockBinding(nm[[k]], ns)
    local({
      key <- k
      f <- orig[[key]]
      assign(nm[[key]],
             function(...) { n[[key]] <<- n[[key]] + 1L; f(...) },
             envir = ns)
    })
  }

  force(expr)
  n
}

## Study 08's own inputs, at the row counts the assertions need. Values are
## arbitrary but COMPLETE: a row with a missing cell is skipped by the route and
## would silently make every count zero, i.e. reproduce the wrong claim.
.pp_rows <- function(k) {
  pick <- function(v) rep_len(v, k)
  list(n_exp = pick(c(40, 55)), n_nexp = pick(c(42, 51)),
       mean_pre_exp = pick(c(10, 12)),        mean_exp = pick(c(13, 16)),
       mean_pre_sd_exp = pick(c(2, 2.5)),     mean_sd_exp = pick(c(2.2, 2.7)),
       mean_pre_nexp = pick(c(10.1, 12.2)),   mean_nexp = pick(c(11, 13)),
       mean_pre_sd_nexp = pick(c(2.1, 2.4)),  mean_sd_nexp = pick(c(2.3, 2.6)),
       r_pre_post_exp = pick(0.8), r_pre_post_nexp = pick(0.8))
}


test_that("study 08's route runs the single-group kernel TWICE PER ROW, not never", {
  .needs_pkg()
  ## The exact call study 08 makes (studies/08_pre_post_to_smd.R:
  ## es_from_means_sd_pre_post with pool_sd left at its default). If this returned
  ## single_group = 0 the README's old "never run" would have been right.
  for (k in c(1L, 2L, 3L)) {
    n <- .count_kernel_calls(
      do.call(es_from_means_sd_pre_post,
              c(.pp_rows(k), list(pre_post_to_smd = "bonett"))))
    expect_equal(unname(n[["single_group"]]), 2L * k)
    expect_equal(unname(n[["pooled"]]), 0L)
  }
  ## Not an artefact of the default: naming it explicitly gives the same counts.
  n <- .count_kernel_calls(
    do.call(es_from_means_sd_pre_post,
            c(.pp_rows(2L), list(pre_post_to_smd = "bonett", pool_sd = FALSE))))
  expect_equal(unname(n[["single_group"]]), 4L)
  expect_equal(unname(n[["pooled"]]), 0L)
})


test_that("pool_sd = TRUE swaps kernels entirely -- one pooled call, no per-arm calls", {
  .needs_pkg()
  for (k in c(1L, 2L)) {
    n <- .count_kernel_calls(
      do.call(es_from_means_sd_pre_post,
              c(.pp_rows(k), list(pre_post_to_smd = "bonett", pool_sd = TRUE))))
    expect_equal(unname(n[["pooled"]]), 1L * k)
    expect_equal(unname(n[["single_group"]]), 0L)
  }
  ## ... and no study sets it. This is the half of the README claim that survived.
  ## grep is the right instrument: a study that started passing pool_sd would make
  ## the "unsimulated" row wrong in exactly the way item 3 went wrong.
  src <- unlist(lapply(
    list.files(file.path(.sim_root, "studies"), pattern = "[.]R$", full.names = TRUE),
    readLines, warn = FALSE))
  expect_false(any(grepl("pool_sd", src[!grepl("^\\s*#", src)], fixed = TRUE)))
})


test_that("the single-group entry points reach the kernel directly, once per row", {
  .needs_pkg()
  ## 7 of the 8 *_single_group routes. These are what "the single-group family is
  ## unsimulated" means: the kernel runs, but only ever as the two-group route's
  ## per-arm engine, never through the entry points below.
  sg <- list(n_exp = 40, r_pre_post_exp = 0.8)
  calls <- list(
    es_from_means_sd_pre_post_single_group = c(sg, list(
      mean_pre_exp = 10, mean_exp = 13, mean_pre_sd_exp = 2, mean_sd_exp = 2.2)),
    es_from_means_se_pre_post_single_group = c(sg, list(
      mean_pre_exp = 10, mean_exp = 13, mean_pre_se_exp = .3, mean_se_exp = .35)),
    es_from_means_ci_pre_post_single_group = c(sg, list(
      mean_pre_exp = 10, mean_exp = 13,
      mean_pre_ci_lo_exp = 9.4, mean_pre_ci_up_exp = 10.6,
      mean_ci_lo_exp = 12.3, mean_ci_up_exp = 13.7)),
    es_from_mean_change_sd_single_group = c(sg, list(
      mean_change_exp = 3, mean_change_sd_exp = 1.4)),
    es_from_mean_change_se_single_group = c(sg, list(
      mean_change_exp = 3, mean_change_se_exp = .22)),
    es_from_mean_change_ci_single_group = c(sg, list(
      mean_change_exp = 3, mean_change_ci_lo_exp = 2.5, mean_change_ci_up_exp = 3.5)),
    es_from_mean_change_pval_single_group = c(sg, list(
      mean_change_exp = 3, mean_change_pval_exp = .001)))

  for (f in names(calls)) {
    n <- .count_kernel_calls(do.call(f, calls[[f]]))
    expect_equal(unname(n[["single_group"]]), 1L, info = f)
    expect_equal(unname(n[["pooled"]]), 0L, info = f)
  }
})


test_that("the paired-t/F family now DELEGATES to the single-group kernel", {
  .needs_pkg()
  ## THIS BLOCK PREVIOUSLY ASSERTED THE OPPOSITE, and the change is the point.
  ##
  ## When item 4.4 wrote this file, these five routes reimplemented the single-group
  ## arithmetic in place -- a fourth pre/post code path that no three-row table could
  ## represent, and that a change to either kernel would leave untouched. That
  ## duplication was removed: they now route through .paired_t_to_smd(), which hands the
  ## kernel an equivalent synthetic single-group problem.
  ##
  ## So the counts flip from 0 to the ordinary ones -- 2 per row for the two-group
  ## routes (one per arm), 1 for the single-group route. Worth recording HOW the change
  ## surfaced: this assertion failed, which is exactly the job 4.4 gave it. The claim
  ## had been made executable, so restructuring the dispatch could not silently
  ## invalidate the paragraph describing it.
  tg <- list(n_exp = 40, n_nexp = 42, r_pre_post_exp = 0.8, r_pre_post_nexp = 0.8)
  calls <- list(
    es_from_paired_t      = c(tg, list(paired_t_exp = 4.2, paired_t_nexp = 1.8)),
    es_from_paired_t_pval = c(tg, list(paired_t_pval_exp = .001, paired_t_pval_nexp = .04)),
    es_from_paired_f      = c(tg, list(paired_f_exp = 17.6, paired_f_nexp = 3.2)),
    es_from_paired_f_pval = c(tg, list(paired_f_pval_exp = .001, paired_f_pval_nexp = .04)),
    es_from_paired_t_single_group = list(paired_t_exp = 4.2, n_exp = 40,
                                         r_pre_post_exp = 0.8))
  expected <- c(es_from_paired_t = 2L, es_from_paired_t_pval = 2L,
                es_from_paired_f = 2L, es_from_paired_f_pval = 2L,
                es_from_paired_t_single_group = 1L)
  for (f in names(calls)) {
    n <- .count_kernel_calls(do.call(f, calls[[f]]))
    expect_equal(unname(n[["single_group"]]), expected[[f]], info = f)
    ## still never the pooled kernel: a paired t does not identify the two arms' SD
    ## ratio, so the pooled standardizer is unrecoverable from it and these routes
    ## expose no pool_sd argument
    expect_equal(unname(n[["pooled"]]), 0L, info = f)
  }
})


test_that("the 19 routes taking pre_post_to_smd split 7 / 8 / 5 across the paths", {
  .needs_pkg()
  ns <- asNamespace("metaConvert")
  ex <- sort(grep("^es_from_", getNamespaceExports("metaConvert"), value = TRUE))
  pp <- ex[vapply(ex, function(f)
    "pre_post_to_smd" %in% names(formals(get(f, envir = ns))), logical(1))]
  ## The README quotes 19. Asserting it here rather than trusting the prose is the
  ## whole point: adding a pre/post route now edits a test, not just a paragraph.
  expect_equal(length(pp), 19L)

  ## Only the two-group means / mean-change routes can pool. The single-group and
  ## paired families expose no pool_sd -- a paired t does not identify the two arms'
  ## SD ratio, so the pooled standardizer is unrecoverable from it.
  poolable <- pp[vapply(pp, function(f)
    "pool_sd" %in% names(formals(get(f, envir = ns))), logical(1))]
  expect_equal(length(poolable), 7L)
  expect_false(any(grepl("_single_group$|^es_from_paired_", poolable)))

  ## 8 single-group-named routes, but only 7 are on the kernel path: the eighth,
  ## es_from_paired_t_single_group, is on the inline paired path above.
  expect_equal(length(grep("_single_group$", pp)), 8L)
  expect_equal(length(grep("^es_from_paired_", pp)), 5L)
  expect_true("es_from_paired_t_single_group" %in% grep("_single_group$", pp, value = TRUE))
})


test_that("study 99 is the only simulation touching the inline paired path, and covers neither the single-group family nor pool_sd", {
  .needs_pkg()
  ## The README's "how it is covered" column claims exactly this split. Reading the
  ## study sources is the only way to keep it honest: study 99 is a deterministic
  ## wiring check that writes no aggregate to inspect.
  read_study <- function(f) {
    x <- readLines(file.path(.sim_root, "studies", f), warn = FALSE)
    x[!grepl("^\\s*#", x)]
  }
  s99 <- read_study("99_route_equivalence.R")
  expect_true(any(grepl("es_from_paired_t(", s99, fixed = TRUE)))
  expect_true(any(grepl("es_from_paired_f(", s99, fixed = TRUE)))
  expect_false(any(grepl("_single_group", s99, fixed = TRUE)))

  others <- setdiff(list.files(file.path(.sim_root, "studies"), pattern = "[.]R$"),
                    "99_route_equivalence.R")
  rest <- unlist(lapply(others, read_study))
  expect_false(any(grepl("es_from_paired_", rest, fixed = TRUE)))
  expect_false(any(grepl("_single_group", rest, fixed = TRUE)))
})
