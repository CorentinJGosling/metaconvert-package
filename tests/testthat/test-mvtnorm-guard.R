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
