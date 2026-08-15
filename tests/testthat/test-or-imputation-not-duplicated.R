# =============================================================================
# convert_df() must not run the OR-alone SE imputation on a row that reports its
# own uncertainty.
#
# es_from_or() has no reported uncertainty to work from, so it IMPUTES var(logOR)
# by enumerating every 2x2 compatible with the case/control margins and averaging
# (.se_from_or). That average runs ~1.4x wide.
#
# On a row that ALSO reports an SE or a CI, the route is a strictly dominated
# duplicate of es_from_or_se()/es_from_or_ci(): identical point estimate, worse
# variance. Before this fix convert_df() ran both, so a single study with
# `or + logor_se` produced two effect sizes -- one with the reported SE and one with
# an imputed SE ~1.6x larger. The hierarchy picked the right one, but the dominated
# row still inflated n_estimations, was eligible for selection under
# es_selected = "minimum"/"maximum", and entered the Category-E cross-method
# dispersion and CI-overlap checks as a spurious second "method".
#
# p-values are deliberately NOT suppressors: the hierarchy ranks es_odds_ratio ABOVE
# es_odds_ratio_pval (or_list_L2 in R/main_convert_df.R), so on a pval-only row the
# marginal reconstruction is the better source and must still run.
# =============================================================================

.or_row <- function(...) {
  base <- list(study_id = "s1", or = 2, n_cases = 60, n_controls = 140,
               n_exp = 100, n_nexp = 100)
  do.call(data.frame, c(base, list(...)))
}

# Does the OR-alone route contribute an estimate for this input?
.imputation_ran <- function(d) {
  res <- suppressMessages(convert_df(d, measure = "logor", verbose = FALSE))
  z <- res[["es_odds_ratio"]]
  !is.null(z) && any(!is.na(z$logor_se))
}

test_that("a reported logor_se suppresses the OR-alone imputation", {
  expect_false(.imputation_ran(.or_row(logor_se = 0.30)))
})

test_that("a reported OR confidence interval suppresses the imputation", {
  expect_false(.imputation_ran(.or_row(or_ci_lo = 1.11, or_ci_up = 3.60)))
})

test_that("a reported log-scale confidence interval suppresses the imputation", {
  expect_false(.imputation_ran(.or_row(logor_ci_lo = log(1.11), logor_ci_up = log(3.60))))
})

test_that("an OR with NO reported uncertainty still gets the imputation", {
  # The route must not be disabled outright -- imputing is the whole point when
  # nothing else is available.
  expect_true(.imputation_ran(.or_row()))
})

test_that("a p-value does NOT suppress the imputation (it ranks below it)", {
  expect_true(.imputation_ran(.or_row(or_pval = 0.02)))
})

test_that("a half-specified CI does not suppress the imputation", {
  # One bound alone cannot produce an SE, so it must not disable the fallback.
  expect_true(.imputation_ran(.or_row(or_ci_lo = 1.11)))
  expect_true(.imputation_ran(.or_row(or_ci_up = 3.60)))
})

test_that("the surviving estimate carries the REPORTED standard error", {
  # Default split_adjusted = TRUE / format = "wide" suffixes the output columns.
  d <- .or_row(logor_se = 0.30)
  res <- suppressMessages(convert_df(d, measure = "logor", verbose = FALSE))
  out <- suppressMessages(summary(res))
  expect_equal(out$se_crude[1], 0.30, tolerance = 1e-8,
               info = "summary() must report the study's own SE, not an imputed one.")
})

test_that("suppression is per-row, not per-dataset", {
  # A dataset mixing reported and unreported uncertainty must impute for exactly the
  # rows that need it.
  d <- data.frame(
    study_id  = c("reported", "unreported"),
    or        = c(2, 2),
    logor_se  = c(0.30, NA),
    n_cases   = c(60, 60), n_controls = c(140, 140),
    n_exp     = c(100, 100), n_nexp = c(100, 100)
  )
  res <- suppressMessages(convert_df(d, measure = "logor", verbose = FALSE))
  z <- res[["es_odds_ratio"]]
  expect_true(is.na(z$logor_se[1]), info = "row 1 reports an SE; must not impute")
  expect_false(is.na(z$logor_se[2]), info = "row 2 reports none; must impute")
})

test_that("suppressing the duplicate does not change the selected effect size", {
  # The hierarchy already ranked es_odds_ratio_se first, so the ESTIMATE must be
  # untouched by this change -- only the spurious extra row disappears.
  d <- .or_row(logor_se = 0.30)
  out <- suppressMessages(summary(suppressMessages(
    convert_df(d, measure = "logor", verbose = FALSE))))
  # summary() rounds its display columns, hence the loose tolerance.
  expect_equal(out$es_crude[1], log(2), tolerance = 1e-3)
  expect_equal(as.character(out$info_used_crude[1]), "or_se")
})
