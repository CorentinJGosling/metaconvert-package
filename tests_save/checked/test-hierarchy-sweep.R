# Long-running file (> 25 s): executed locally and in CI (devtools::test and
# testthat set NOT_CRAN=true); skipped wholesale on CRAN to respect check-time
# limits. The fast pre/post files still run on CRAN.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
# ==============================================================================
# Regression sweep: every measure must survive convert_df() through both the
# `es_selected = "auto"` branch and the `es_selected = "hierarchy"` branch.
# ------------------------------------------------------------------------------
# The `length(res) == expected_length` invariant in main_convert_df.R is easy to
# break silently when new methods are added to one branch but not the other.
# Example: logvr / logcvr was missing RD_stand after RD was introduced, so
# `convert_df(df.haza, measure = "logvr")` stopped with
#   "Expected 76, got 73"
# This sweep catches any future occurrence of the same class of bug.
# ==============================================================================

library(testthat)
library(metaConvert)

test_that("every measure is reachable through both es_selected branches", {
  measures_76 <- c("d", "g", "md", "or", "logor", "rr", "logrr",
                   "irr", "logirr", "nnt", "rd", "r", "z", "logvr", "logcvr")
  measures_other <- list(
    "dw" = 9, "gw" = 9, "mdw" = 9,
    "prop" = 3, "alpha" = 2, "icc" = 2,
    "rp" = 5, "zp" = 5
  )

  for (m in measures_76) {
    for (sel in c("auto", "hierarchy")) {
      res <- expect_no_error(
        convert_df(df.haza, measure = m, es_selected = sel, verbose = FALSE),
        message = sprintf("convert_df failed for measure='%s', es_selected='%s'", m, sel)
      )
      expect_equal(length(res), 76L,
                   info = sprintf("measure='%s', es_selected='%s'", m, sel))
      expect_s3_class(res, "metaConvert")
    }
  }

  for (m in names(measures_other)) {
    for (sel in c("auto", "hierarchy")) {
      res <- expect_no_error(
        convert_df(df.haza, measure = m, es_selected = sel, verbose = FALSE),
        message = sprintf("convert_df failed for measure='%s', es_selected='%s'", m, sel)
      )
      expect_equal(length(res), measures_other[[m]],
                   info = sprintf("measure='%s', es_selected='%s'", m, sel))
    }
  }
})

} # end NOT_CRAN gate
