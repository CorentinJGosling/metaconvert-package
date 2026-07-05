test_that("es_from_chisq accepts vector yates_chisq", {
  # Two rows: row 1 uncorrected Pearson, row 2 Yates-corrected
  chisq_vec <- c(0.2224, 0.2680)
  n        <- c(60, 60)
  n_exp    <- c(30, 30)
  n_cases  <- c(56, 56)

  res_mixed   <- es_from_chisq(chisq = chisq_vec, n_sample = n, n_cases = n_cases,
                               n_exp = n_exp, yates_chisq = c(FALSE, TRUE))
  res_pearson <- es_from_chisq(chisq = chisq_vec, n_sample = n, n_cases = n_cases,
                               n_exp = n_exp, yates_chisq = c(FALSE, FALSE))
  res_yates   <- es_from_chisq(chisq = chisq_vec, n_sample = n, n_cases = n_cases,
                               n_exp = n_exp, yates_chisq = c(TRUE, TRUE))

  # Row 1 in the mixed call should match the all-Pearson value
  expect_equal(res_mixed$logor[1], res_pearson$logor[1], tolerance = 1e-6)
  # Row 2 in the mixed call should match the all-Yates value
  expect_equal(res_mixed$logor[2], res_yates$logor[2], tolerance = 1e-6)
})

test_that("es_from_chisq_pval accepts vector yates_chisq", {
  pvals    <- c(0.605, 0.605)
  n        <- c(60, 60)
  n_exp    <- c(30, 30)
  n_cases  <- c(56, 56)

  res_mixed   <- es_from_chisq_pval(chisq_pval = pvals, n_sample = n,
                                    n_cases = n_cases, n_exp = n_exp,
                                    yates_chisq = c(FALSE, TRUE))
  res_pearson <- es_from_chisq_pval(chisq_pval = pvals, n_sample = n,
                                    n_cases = n_cases, n_exp = n_exp,
                                    yates_chisq = FALSE)
  res_yates   <- es_from_chisq_pval(chisq_pval = pvals, n_sample = n,
                                    n_cases = n_cases, n_exp = n_exp,
                                    yates_chisq = TRUE)

  expect_equal(res_mixed$logor[1], res_pearson$logor[1], tolerance = 1e-6)
  expect_equal(res_mixed$logor[2], res_yates$logor[2], tolerance = 1e-6)
  # Pearson and Yates back-derive different cells for this p-value
  expect_false(isTRUE(all.equal(res_pearson$logor[1], res_yates$logor[1])))
})

test_that("es_from_chisq errors on yates_chisq length mismatch", {
  expect_error(
    es_from_chisq(chisq = c(0.2, 0.3), n_sample = c(60, 60),
                  n_cases = c(56, 56), n_exp = c(30, 30),
                  yates_chisq = c(TRUE, FALSE, TRUE)),
    "yates_chisq"
  )
})

test_that("convert_df honors per-row yates_chisq column and falls back to scalar for NAs", {
  df <- data.frame(
    study_id = c("A", "B"),
    author   = c("Pearson", "Yates"),
    year     = c(2020, 2020),
    n_sample = c(60, 60),
    n_exp    = c(30, 30),
    n_cases  = c(56, 56),
    chisq_pval  = c(0.605, 0.605),
    yates_chisq = c(FALSE, TRUE)
  )

  mc <- convert_df(df, measure = "or", es_selected = "hierarchy", verbose = FALSE)
  out_per_row <- mc[["es_chisq_pval"]]

  # Same data without the per-row column, with scalar yates_chisq = TRUE,
  # should reproduce row 2 of the per-row run (which used TRUE).
  df_scalar <- df[, setdiff(colnames(df), "yates_chisq")]
  mc_yates <- convert_df(df_scalar, measure = "or", es_selected = "hierarchy",
                         yates_chisq = TRUE, verbose = FALSE)
  out_yates <- mc_yates[["es_chisq_pval"]]

  # Per-row row 1 (FALSE) should NOT equal the all-Yates row 1
  expect_false(isTRUE(all.equal(out_per_row$logor[1], out_yates$logor[1])))
  # Per-row row 2 (TRUE) SHOULD equal the all-Yates row 2
  expect_equal(out_per_row$logor[2], out_yates$logor[2], tolerance = 1e-6)

  # NA fallback: c(FALSE, NA) + scalar TRUE -> row 2 picks up Yates from scalar
  df_na <- df
  df_na$yates_chisq <- c(FALSE, NA)
  mc_na <- convert_df(df_na, measure = "or", es_selected = "hierarchy",
                      yates_chisq = TRUE, verbose = FALSE)
  out_na <- mc_na[["es_chisq_pval"]]
  expect_equal(out_na$logor[1], out_per_row$logor[1], tolerance = 1e-6)
  expect_equal(out_na$logor[2], out_per_row$logor[2], tolerance = 1e-6)
})

test_that("data_extraction_sheet includes yates_chisq column (not chisq_yates)", {
  sheet <- data_extraction_sheet(measure = "or", verbose = FALSE)
  expect_true("yates_chisq" %in% colnames(sheet))
  expect_false("chisq_yates" %in% colnames(sheet))
})

test_that("convert_df does not silently consume an old chisq_yates column", {
  # Guard against the previous naming. If a user (or a stale xlsx) still has
  # `chisq_yates`, it must NOT be silently treated as `yates_chisq`.
  df <- data.frame(
    study_id = c("A", "B"),
    author   = c("Pearson", "Yates"),
    year     = c(2020, 2020),
    n_sample = c(60, 60),
    n_exp    = c(30, 30),
    n_cases  = c(56, 56),
    chisq_pval  = c(0.605, 0.605),
    chisq_yates = c(FALSE, TRUE)  # old column name; should be ignored
  )

  # Setting the scalar yates_chisq = TRUE; with the old column being ignored,
  # both rows should use Yates and give identical logor.
  mc <- convert_df(df, measure = "or", es_selected = "hierarchy",
                   yates_chisq = TRUE, verbose = FALSE)
  out <- mc[["es_chisq_pval"]]
  expect_equal(out$logor[1], out$logor[2], tolerance = 1e-6)
})
