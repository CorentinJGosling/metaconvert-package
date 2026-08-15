# =============================================================================
# Roadmap item 1.2 -- V6 / r_defaulted must fire on EVERY r-consuming data shape.
#
# convert_df() imputes r_pre_post = 0.8 when the user leaves it blank. That is a
# substantive assumption -- under the change-SD standardizers (morris_dz, cooper /
# morris_drm) the assumed correlation scales the POINT estimate, not just the SE --
# so the package raises flag V6, sets attr(res, "r_defaulted"), and uses that
# attribute to annotate A6/E2/E2b as "(r-sensitive: ...)" in summary().
#
# The bug this file pins: the V6 block carried its own hardcoded column list in
# which 7 of 11 names were not columns at all --
#   mean_post_exp, mean_post_nexp   (endpoint columns are mean_exp / mean_nexp)
#   mean_pre_single_group, mean_post_single_group, mean_change_single_group
#                                   (no *_single_group columns exist; the
#                                    single-group routes reuse the _exp columns)
#   paired_t, paired_f              (real names are paired_t_exp / paired_t_nexp,
#                                    paired_f_exp / paired_f_nexp)
# They were filtered out by intersect(., colnames(x)), so the whole check went
# silent on paired-t/F rows -- precisely the rows where r moves the point estimate.
#
# Both call sites now share .r_consuming_columns() / .rows_with_r_consuming_data().
# =============================================================================

# --- helpers ---------------------------------------------------------------

# Minimal row skeleton: every column convert_df needs, all NA except what a test
# fills in. Built from the real extraction sheet so a column name that does not
# exist cannot be introduced by this test file either.
.v6_row <- function(...) {
  sheet <- data_extraction_sheet()
  x <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = ncol(sheet)))
  names(x) <- names(sheet)
  vals <- list(...)
  for (nm in names(vals)) {
    if (!nm %in% names(x)) {
      stop("test bug: '", nm, "' is not a real input column", call. = FALSE)
    }
    x[[nm]] <- vals[[nm]]
  }
  x$study_id <- "s1"
  x
}

.r_defaulted_of <- function(x, ...) {
  res <- suppressMessages(convert_df(x, measure = "g", verbose = FALSE, ...))
  attr(res, "r_defaulted")
}

.v6_fired <- function(x, ...) {
  res <- suppressMessages(convert_df(x, measure = "g", verbose = FALSE, ...))
  s <- suppressMessages(summary(res, flags = TRUE))
  fl <- unlist(s[, grep("^flags", names(s)), drop = FALSE])
  any(grepl("Default r_pre_post", fl, fixed = TRUE))
}

# Every shape below reaches a route that takes r_pre_post_exp/_nexp.
.shapes <- list(
  "pre/post means + SDs (two-group)" = list(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 10, mean_pre_sd_exp = 3, mean_exp = 13, mean_sd_exp = 3,
    mean_pre_nexp = 10, mean_pre_sd_nexp = 3, mean_nexp = 11, mean_sd_nexp = 3
  ),
  "pre/post means + SEs" = list(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 10, mean_pre_se_exp = 0.42, mean_exp = 13, mean_se_exp = 0.42,
    mean_pre_nexp = 10, mean_pre_se_nexp = 0.42, mean_nexp = 11, mean_se_nexp = 0.42
  ),
  "mean change + SD" = list(
    n_exp = 50, n_nexp = 50,
    mean_change_exp = 3, mean_change_sd_exp = 4,
    mean_change_nexp = 1, mean_change_sd_nexp = 4
  ),
  "mean change + SE (no mean_change_sd)" = list(
    n_exp = 50, n_nexp = 50,
    mean_change_exp = 3, mean_change_se_exp = 0.57,
    mean_change_nexp = 1, mean_change_se_nexp = 0.57
  ),
  "mean change + CI" = list(
    n_exp = 50, n_nexp = 50,
    mean_change_exp = 3, mean_change_ci_lo_exp = 1.9, mean_change_ci_up_exp = 4.1,
    mean_change_nexp = 1, mean_change_ci_lo_nexp = -0.1, mean_change_ci_up_nexp = 2.1
  ),
  "mean change + p-value" = list(
    n_exp = 50, n_nexp = 50,
    mean_change_exp = 3, mean_change_pval_exp = 0.001,
    mean_change_nexp = 1, mean_change_pval_nexp = 0.20
  ),
  "paired t (the shape the bug lost)" = list(
    n_exp = 50, n_nexp = 50, paired_t_exp = 5.3, paired_t_nexp = 1.8
  ),
  "paired t p-value" = list(
    n_exp = 50, n_nexp = 50,
    paired_t_pval_exp = 0.001, paired_t_pval_nexp = 0.08
  ),
  "paired F" = list(
    n_exp = 50, n_nexp = 50, paired_f_exp = 28.1, paired_f_nexp = 3.2
  ),
  "paired F p-value" = list(
    n_exp = 50, n_nexp = 50,
    paired_f_pval_exp = 0.001, paired_f_pval_nexp = 0.08
  ),
  "single-group pre/post (reuses _exp columns)" = list(
    n_exp = 50,
    mean_pre_exp = 10, mean_pre_sd_exp = 3, mean_exp = 13, mean_sd_exp = 3
  ),
  "single-group paired t" = list(
    n_exp = 50, paired_t_exp = 5.3
  )
)

# --- the regression this file exists for ------------------------------------

test_that("the r-consuming column list contains only real input columns", {
  # This is the assertion that would have caught the original bug directly.
  real <- names(data_extraction_sheet())
  phantoms <- setdiff(metaConvert:::.r_consuming_columns(), real)
  expect_equal(
    phantoms, character(0),
    info = paste("These names are not input columns, so intersect() silently drops",
                 "them and the check never fires on those rows:",
                 paste(phantoms, collapse = ", "))
  )
})

test_that("r_defaulted fires for every r-consuming data shape", {
  for (nm in names(.shapes)) {
    x <- do.call(.v6_row, .shapes[[nm]])
    rd <- .r_defaulted_of(x)
    expect_true(
      isTRUE(any(rd)),
      info = paste0("r_defaulted is FALSE for '", nm, "'. The assumed r_pre_post ",
                    "is unflagged, so summary() will not annotate A6/E2/E2b as ",
                    "r-sensitive on this row.")
    )
  }
})

test_that("V6 flag text reaches summary() for every r-consuming data shape", {
  for (nm in names(.shapes)) {
    x <- do.call(.v6_row, .shapes[[nm]])
    expect_true(
      .v6_fired(x),
      info = paste0("No 'Default r_pre_post' flag for '", nm, "'.")
    )
  }
})

test_that("paired-t and paired-F rows specifically are flagged (the lost shapes)", {
  # Called out separately because these are the rows where the assumed r scales the
  # POINT estimate, not merely the SE -- the worst rows to lose.
  for (nm in c("paired t (the shape the bug lost)", "paired F",
               "paired t p-value", "paired F p-value",
               "single-group paired t")) {
    x <- do.call(.v6_row, .shapes[[nm]])
    expect_true(isTRUE(any(.r_defaulted_of(x))), info = nm)
  }
})

# --- negative controls: the flag must stay quiet when it should -------------

test_that("a row with no pre/post, change or paired data is NOT flagged", {
  x <- .v6_row(
    n_exp = 50, n_nexp = 50,
    mean_exp = 13, mean_sd_exp = 3, mean_nexp = 11, mean_sd_nexp = 3
  )
  expect_false(isTRUE(any(.r_defaulted_of(x))),
               info = "Endpoint-only two-group data must not raise V6.")
  expect_false(.v6_fired(x))
})

test_that("a row supplying r_pre_post explicitly is NOT flagged", {
  x <- .v6_row(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 10, mean_pre_sd_exp = 3, mean_exp = 13, mean_sd_exp = 3,
    mean_pre_nexp = 10, mean_pre_sd_nexp = 3, mean_nexp = 11, mean_sd_nexp = 3,
    r_pre_post_exp = 0.65, r_pre_post_nexp = 0.65
  )
  expect_false(isTRUE(any(.r_defaulted_of(x))),
               info = "A user-supplied r_pre_post must not be reported as defaulted.")
  expect_false(.v6_fired(x))
})

test_that("endpoint columns are excluded from the r-consuming list by design", {
  # mean_exp / mean_sd_exp are populated on ordinary two-group rows that never touch
  # an r-consuming route; including them would make V6 fire on non-paired data.
  cols <- metaConvert:::.r_consuming_columns()
  expect_false(any(c("mean_exp", "mean_nexp", "mean_sd_exp", "mean_sd_nexp") %in% cols))
})

test_that("the verbose note and flag V6 agree on which rows are r-consuming", {
  # The two call sites drifted apart before this fix; assert they now agree.
  x_paired <- do.call(.v6_row, .shapes[["paired t (the shape the bug lost)"]])
  expect_message(
    convert_df(x_paired, measure = "g", verbose = TRUE),
    "r_pre_post was not reported"
  )
  expect_true(isTRUE(any(.r_defaulted_of(x_paired))))

  x_plain <- .v6_row(n_exp = 50, n_nexp = 50, mean_exp = 13, mean_sd_exp = 3,
                     mean_nexp = 11, mean_sd_nexp = 3)
  expect_no_message(
    convert_df(x_plain, measure = "g", verbose = TRUE),
    message = "r_pre_post was not reported"
  )
  expect_false(isTRUE(any(.r_defaulted_of(x_plain))))
})
