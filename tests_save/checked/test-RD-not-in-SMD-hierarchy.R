# ==============================================================================
# Regression: RD is NOT a source for the SMD (d/g) or correlation (r/z)
# families in ANY es_selected mode.
#
# Rationale (R/es_from_stand_RD.R, R/main_convert_df.R): OR and the raw 2x2
# table are the sole gateway into the standardized families via the Cox
# transform d = log(OR) * sqrt(3)/pi. RD -> SMD/r/z is not identified without an
# assumed baseline_risk, and because baseline_risk is treated as a fixed known
# constant the resulting SEs are anti-conservative. RD_stand is therefore
# excluded from the d/g/md and r/z hierarchies.
#
# This exclusion previously held ONLY under es_selected = "auto"; the
# "hierarchy" / "minimum" / "maximum" modes leaked RD into d/r via the generic
# fallthrough branch, making the same measure = "r" call mode-dependent. This
# test locks the fix: RD must be excluded from SMD/r/z in every mode, while its
# legitimate ratio conversions (OR / RR / NNT / RD) are preserved.
# ==============================================================================

library(testthat)
library(metaConvert)

# RD is the ONLY binary source: no OR, no 2x2 cells, no means, no correlation.
rd_only <- data.frame(
  rd = c(-0.15, 0.20),
  rd_se = c(0.05, 0.06),
  baseline_risk = c(0.30, 0.25),
  n_exp = c(100, 120),
  n_nexp = c(100, 110)
)

test_that("RD does not feed the SMD/correlation families in any selection mode", {
  for (sel in c("auto", "minimum", "maximum")) {
    for (m in c("d", "g", "r", "z")) {
      res <- summary(convert_df(rd_only, measure = m, es_selected = sel,
                                verbose = FALSE))
      expect_true(
        all(is.na(res$es_crude)),
        info = sprintf("measure='%s', es_selected='%s' must be NA (RD excluded from SMD/r/z)", m, sel)
      )
    }
  }
})

test_that("RD still feeds its legitimate ratio conversions (OR/RR/NNT/RD) in every mode", {
  for (sel in c("auto", "minimum", "maximum")) {
    for (m in c("or", "rr", "nnt", "rd")) {
      res <- summary(convert_df(rd_only, measure = m, es_selected = sel,
                                verbose = FALSE))
      expect_false(
        any(is.na(res$es_crude)),
        info = sprintf("measure='%s', es_selected='%s' must be RD-derived (non-NA)", m, sel)
      )
      expect_true(
        all(res$info_used_crude %in% c("rd_se", "rd_ci", "rd_pval")),
        info = sprintf("measure='%s', es_selected='%s' must be sourced from RD", m, sel)
      )
    }
  }
})

test_that("es_from_rd_se() the FUNCTION does not produce d/g/r/z (removed everywhere)", {
  # RD -> SMD/correlation is removed at the function level too, so the function
  # and convert_df() agree: RD yields the ratio family (OR/RR/NNT/RD) only.
  res <- es_from_rd_se(rd = -0.15, rd_se = 0.05, baseline_risk = 0.30,
                       n_exp = 100, n_nexp = 100)
  # SMD / correlation columns must be absent.
  expect_false("d" %in% names(res))
  expect_false("g" %in% names(res))
  expect_false("r" %in% names(res))
  expect_false("z" %in% names(res))
  # Ratio-family conversions must remain.
  expect_false(is.na(res$logor))
  expect_false(is.na(res$logrr))
  expect_false(is.na(res$rd))
  expect_false(is.na(res$nnt))
})
