# Cross-row checks under main_es = FALSE must judge a comparison on the route the
# HIERARCHY selects, not on whichever route happens to come first in the internal list.
#
# .flag_es_quality() reduces each comparison to one representative row before running the
# cross-row checks (outliers, direction conflict, duplication). That representative used
# to be `which(!duplicated(cmp_key))`, i.e. the first row of the group. Under
# main_es = FALSE the rows of a comparison appear in the fixed internal list order, which
# does NOT track the hierarchy, so the verdict was computed from an estimate the user
# never sees: an outlier carried by the selected route went unreported in the route view
# while main_es = TRUE flagged it.

make_dat <- function() {
  set.seed(42)
  n <- 12
  data.frame(
    study_id    = paste0("S", seq_len(n)),
    n_exp       = 40,
    n_nexp      = 40,
    mean_exp    = 12 + round(rnorm(n, 0, 0.2), 2),
    mean_nexp   = 10,
    mean_sd_exp = 4,
    mean_sd_nexp = 4,
    # S12's reported t is inconsistent with its own means and SDs
    student_t   = c(round(rnorm(n - 1, 2.2, 0.35), 2), 13.0)
  )
}

outlier_flags <- function(hierarchy, main_es) {
  s <- summary(
    convert_df(make_dat(), measure = "g", es_selected = "hierarchy",
               hierarchy = hierarchy, main_es = main_es,
               split_adjusted = FALSE, verbose = FALSE),
    verbose = FALSE, include_raw = FALSE)
  hit <- grepl("outlier", s$flags, ignore.case = TRUE)
  list(n = sum(hit), studies = unique(as.character(s$study_id[hit])))
}


test_that("the route view flags the outlier the hierarchy actually selects", {
  # student_t is selected, and S12's t is the inconsistent one
  a <- outlier_flags("student_t > means_sd", main_es = TRUE)
  b <- outlier_flags("student_t > means_sd", main_es = FALSE)

  expect_equal(a$studies, "S12")
  # the route view broadcasts the verdict to that comparison's rows, so the count may
  # differ, but it must name the same study -- and must not be silent
  expect_gt(b$n, 0)
  expect_equal(b$studies, "S12")
})


test_that("the route view stays silent when the selected route is sound", {
  # means_sd is selected; it is consistent for every study, so nothing should fire even
  # though the inconsistent student_t route is still present in the frame
  a <- outlier_flags("means_sd > student_t", main_es = TRUE)
  b <- outlier_flags("means_sd > student_t", main_es = FALSE)

  expect_equal(a$n, 0)
  expect_equal(b$n, 0)
})


test_that("the two views agree on which studies are flagged, under both hierarchies", {
  for (h in c("student_t > means_sd", "means_sd > student_t")) {
    a <- outlier_flags(h, main_es = TRUE)
    b <- outlier_flags(h, main_es = FALSE)
    expect_equal(sort(b$studies), sort(a$studies),
                 info = paste("hierarchy:", h))
  }
})


test_that("the default main_es = TRUE path is untouched by the representative choice", {
  # every comparison is a single row there, so the representative is that row whatever
  # rule is used; this guards the no-op claim against future edits
  d <- make_dat()
  s <- summary(convert_df(d, measure = "g", split_adjusted = FALSE, verbose = FALSE),
               verbose = FALSE, include_raw = FALSE)
  expect_equal(nrow(s), nrow(d))
  expect_equal(as.character(s$study_id), as.character(d$study_id))
})
