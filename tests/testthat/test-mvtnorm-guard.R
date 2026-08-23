# The tetrachoric correlation from a 2x2 table is solved by metafor::escalc(measure =
# "RTET"), which needs 'mvtnorm'. mvtnorm is a Suggests of metafor, NOT an Imports, so
# installing metafor does not bring it in: an ordinary installation can lack it. Before
# the guard, escalc() raised an error, .tet_r()'s tryCatch swallowed it, and every
# 2x2-derived r and Fisher's z came back NA with nothing said -- indistinguishable from
# data that genuinely cannot support them.
#
# The notice is a message rather than a warning because the call sites wrap this route in
# suppressWarnings() (a numerical failure is expected on some tables and must not spam).

clear_notice <- function() {
  if (exists("mvtnorm_missing", envir = metaConvert:::.mcv_notices)) {
    rm("mvtnorm_missing", envir = metaConvert:::.mcv_notices)
  }
}


test_that("with mvtnorm available the tetrachoric route is silent and numeric", {
  skip_if_not_installed("mvtnorm")

  expect_silent(res <- es_from_2x2(143, 52, 41, 164))
  expect_false(is.na(res$r))
  expect_false(is.na(res$z))
})


test_that("without mvtnorm the user is told why r and z are NA", {
  local_mocked_bindings(.has_mvtnorm = function() FALSE, .package = "metaConvert")
  clear_notice()

  expect_message(res <- es_from_2x2(143, 52, 41, 164), "mvtnorm", fixed = TRUE)
  expect_true(is.na(res$r))
  expect_true(is.na(res$z))

  clear_notice()
  expect_message(es_from_2x2(143, 52, 41, 164), "returned as NA", fixed = TRUE)
})


test_that("without mvtnorm the other effect size measures still compute", {
  local_mocked_bindings(.has_mvtnorm = function() FALSE, .package = "metaConvert")
  clear_notice()

  res <- suppressMessages(es_from_2x2(143, 52, 41, 164))
  # the odds ratio does not go through the tetrachoric route and must be unaffected
  expect_false(is.na(res$logor))
  expect_false(is.na(res$logor_se))
})


test_that("the notice is emitted once per session, not once per row", {
  local_mocked_bindings(.has_mvtnorm = function() FALSE, .package = "metaConvert")
  clear_notice()

  n_msg <- 0
  withCallingHandlers(
    es_from_2x2(c(143, 120, 90), c(52, 60, 70), c(41, 50, 60), c(164, 150, 140)),
    message = function(m) {
      if (grepl("mvtnorm", conditionMessage(m), fixed = TRUE)) n_msg <<- n_msg + 1
      invokeRestart("muffleMessage")
    })
  expect_equal(n_msg, 1)
})


test_that("the phi and chi-squared routes report it too (they go through the 2x2)", {
  local_mocked_bindings(.has_mvtnorm = function() FALSE, .package = "metaConvert")

  clear_notice()
  expect_message(es_from_phi(phi = 0.35, n_cases = 90, n_exp = 100, n_sample = 200),
                 "mvtnorm", fixed = TRUE)

  clear_notice()
  expect_message(es_from_chisq(chisq = 24.5, n_sample = 200, n_cases = 90, n_exp = 100),
                 "mvtnorm", fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Regression pin for CRAN's --no-suggests flavour (2.0.1 audit finding #3).
#
# This rotted TWICE. The files the audit named were guarded in 4955631; commit
# 2b6f141 (roadmap 1.3) then added tests/testthat/test-table-2x2-to-cor.R with 12
# unguarded blocks. Measured with mvtnorm genuinely removed from the library:
# PASS 20 / FAIL 6 before, PASS 13 / FAIL 0 / SKIP 5 after. Not an error -- a
# silent FAIL, which is still an R CMD check ERROR on that flavour.
#
# Scope of this pin, stated honestly: it only re-checks the five blocks measured
# to need mvtnorm. It CANNOT catch a newly added unguarded file, because whether
# a block survives mvtnorm's absence is a runtime property of its assertions --
# a static scan was tried and was unsound both ways (it flagged 4 blocks that
# pass fine and missed the es_from_phi block whose assertion reads a message).
# The general case is a platform property and belongs to the platform check:
# .github/workflows/rhub.yaml now carries noSuggests in its default config, and
# the local reproduction is _R_CHECK_DEPENDS_ONLY_=true R CMD check (see CLAUDE.md).
# ---------------------------------------------------------------------------
test_that("the five mvtnorm-dependent blocks of test-table-2x2-to-cor.R stay guarded", {
  f <- "test-table-2x2-to-cor.R"
  skip_if(!file.exists(f), paste(f, 'not visible from the test working directory'))
  src <- readLines(f, warn = FALSE)

  # Block titles measured to fail when mvtnorm is absent (r/z come back NA).
  needs_mvtnorm <- c(
    "the tetrachoric r CI stays inside",
    "the r CI is the tanh back-transform of the z CI",
    "reverse_2x2 reflects the CI rather than inverting it",
    "es_from_phi() warns loudly when ONE call mixes both estimands",
    "the two es_from_phi regimes really do differ")

  starts <- grep("^test_that[(]", src)
  ends   <- c(starts[-1] - 1L, length(src))
  unguarded <- character(0)
  for (title in needs_mvtnorm) {
    i <- which(vapply(starts, function(s) grepl(title, src[s], fixed = TRUE), logical(1)))
    if (length(i) != 1L) { unguarded <- c(unguarded, paste0(title, ' <block not found>')); next }
    blk <- src[starts[i]:ends[i]]
    if (!any(grepl('skip_if_not_installed[(]["]mvtnorm', blk)))
      unguarded <- c(unguarded, title)
  }
  expect_equal(unguarded, character(0))
})
