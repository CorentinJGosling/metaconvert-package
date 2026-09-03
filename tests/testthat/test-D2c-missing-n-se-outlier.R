# D2c -- the SE-outlier check for rows the sample-size-normalised checks never saw.
#
# D2/D2b normalise an SE by N before comparing it, so a row carrying no usable N is set
# to NA and leaves the pool SILENTLY: no flag, and nothing separating "checked and fine"
# from "never checked". The motivating case is real -- validation_set/save/
# 42125682-ERROR-OR forest row 5 (Wang 2025), a per-1-unit continuous logistic
# coefficient pooled with high-vs-low binary contrasts, holding the largest weight in
# the meta-analysis by virtue of an SE 9.4x tighter than any sibling. NO THRESHOLD
# reaches it: even se_ratio_extreme = 2 is silent, because the row is not in the pool.

or_pool <- function(tight_n = NA) {
  data.frame(
    study_id = paste0("S", 1:6),
    or       = c(2.89, 2.16, 5.49, 2.22, 1.10, 1.06),
    or_ci_lo = c(1.40, 1.26, 2.11, 1.22, 1.03, 0.65),
    or_ci_up = c(5.97, 3.73, 14.28, 4.05, 1.17, 1.73),
    n_sample = c(174, 110, 237, 1170, tight_n, 624),
    stringsAsFactors = FALSE)
}
flags_of <- function(d, ...) {
  s <- as.data.frame(suppressMessages(suppressWarnings(
    summary(convert_df(d, measure = "or", verbose = FALSE), flags = TRUE,
            flag_options = list(...)))))
  if ("flags_crude" %in% names(s)) s$flags_crude else s$flags
}

test_that("D2c is off by default and changes nothing", {
  f <- flags_of(or_pool())
  expect_false(any(grepl("unchecked row", f)))
  # ROW 5 is invisible to the normalised checks at ANY threshold -- the property that
  # makes D2c necessary rather than a threshold tweak. Asserted on row 5 alone: at a
  # gate as loose as 2 the normalised checks legitimately flag OTHER rows (a pool
  # spanning n = 110 to 1170 has a wide normalised spread), which is not the point.
  for (X in c(10, 5, 2)) {
    expect_false(grepl("SE outlier", flags_of(or_pool(), se_ratio_extreme = X)[5]),
                 info = paste("se_ratio_extreme =", X))
  }
})

test_that("D2c flags the N-less row whose raw SE is far tighter than the pool", {
  f <- flags_of(or_pool(), se_outlier_missing_n = TRUE)
  expect_true(grepl("unchecked row", f[5]))
  expect_match(f[5], "9.4x tighter", fixed = TRUE)
  expect_match(f[5], "carries no sample size", fixed = TRUE)
  expect_match(f[5], "[UNUSUAL]", fixed = TRUE)
  # nothing else in the pool
  expect_false(any(grepl("unchecked row", f[-5])))
})

test_that("D2c never fires on a row the normalised checks could see", {
  # same pool, but the row now carries an N: it is normalisable, so D2c must stay out
  # of it and leave the verdict to D2/D2b. This is what bounds D2c to rows that are
  # unchecked by construction -- it can add a flag, never change an existing one.
  f <- flags_of(or_pool(tight_n = 303), se_outlier_missing_n = TRUE)
  expect_false(any(grepl("unchecked row", f)))
})

test_that("se_missing_n_ratio governs the fold threshold", {
  # measured raw fold for this row is 9.4x
  expect_true(grepl("unchecked row",
    flags_of(or_pool(), se_outlier_missing_n = TRUE, se_missing_n_ratio = 5)[5]))
  expect_false(any(grepl("unchecked row",
    flags_of(or_pool(), se_outlier_missing_n = TRUE, se_missing_n_ratio = 10))))
})
