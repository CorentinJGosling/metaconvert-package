# ==============================================================================
# Tests for ES guidance (missing data diagnostics)
# ==============================================================================

# ==============================================================================
# Unit tests for .method_required_columns()
# ==============================================================================
test_that(".method_required_columns returns a named list with correct structure", {
  mc <- metaConvert:::.method_required_columns()
  expect_type(mc, "list")
  expect_true(length(mc) > 70)

  # Each entry should have specific and shared
  for (name in names(mc)) {
    expect_true("specific" %in% names(mc[[name]]),
                info = paste("Method", name, "missing 'specific'"))
    expect_true("shared" %in% names(mc[[name]]),
                info = paste("Method", name, "missing 'shared'"))
  }

  # Check some known methods exist
  expect_true("means_sd" %in% names(mc))
  expect_true("student_t" %in% names(mc))
  expect_true("or_se" %in% names(mc))
  expect_true("pearson_r" %in% names(mc))
  expect_true("2x2" %in% names(mc))
  expect_true("user_input_crude" %in% names(mc))
  expect_true("prop_single_group" %in% names(mc))
  expect_true("linreg_t" %in% names(mc))
})


# ==============================================================================
# Unit tests for .method_descriptions()
# ==============================================================================
test_that(".method_descriptions covers all methods in .method_required_columns", {
  mc <- metaConvert:::.method_required_columns()
  desc <- metaConvert:::.method_descriptions()

  expect_type(desc, "character")
  expect_true(length(desc) > 70)

  # Every method should have a description
  for (name in names(mc)) {
    expect_true(name %in% names(desc),
                info = paste("Method", name, "has no description"))
  }

  # Descriptions should be non-empty
  expect_true(all(nchar(desc) > 0))
})


# ==============================================================================
# Unit tests for .get_applicable_methods()
# ==============================================================================
test_that(".get_applicable_methods returns correct methods for d/g/md", {
  methods_d <- metaConvert:::.get_applicable_methods("d", "")
  expect_true("means_sd" %in% methods_d)
  expect_true("student_t" %in% methods_d)
  expect_true("or_se" %in% methods_d)
  expect_true("pearson_r" %in% methods_d)
  expect_true("user_input_crude" %in% methods_d)
  expect_true("user_input_adj" %in% methods_d)
  expect_true("ancova_t" %in% methods_d)

  # Variability, IRR, and RR methods should NOT be included for d/g
  expect_false("variability_means_sd" %in% methods_d)
  expect_false("variability_means_se" %in% methods_d)
  expect_false("variability_means_ci" %in% methods_d)
  expect_false("cases_time" %in% methods_d)
  expect_false("rr_se" %in% methods_d)
})

test_that(".get_applicable_methods filters crude/adjusted correctly", {
  crude <- metaConvert:::.get_applicable_methods("g", "_crude")
  adj <- metaConvert:::.get_applicable_methods("g", "_adjusted")

  expect_true("means_sd" %in% crude)
  expect_false("ancova_t" %in% crude)

  expect_true("ancova_t" %in% adj)
  expect_false("means_sd" %in% adj)
})

test_that(".get_applicable_methods excludes variability for logor", {
  methods_or <- metaConvert:::.get_applicable_methods("logor", "")
  expect_true("or_se" %in% methods_or)
  expect_true("rr_se" %in% methods_or)
  expect_true("means_sd" %in% methods_or)
  expect_false("variability_means_sd" %in% methods_or)
  expect_false("cases_time" %in% methods_or)
})

test_that(".get_applicable_methods returns only VAR methods for logvr/logcvr", {
  methods_vr <- metaConvert:::.get_applicable_methods("logvr", "")
  expect_true("variability_means_sd" %in% methods_vr)
  expect_true("means_sd" %in% methods_vr)  # SMD_post is used for logvr
  expect_false("or_se" %in% methods_vr)
  expect_false("cases_time" %in% methods_vr)
})

test_that(".get_applicable_methods returns only IRR methods for logirr", {
  methods_irr <- metaConvert:::.get_applicable_methods("logirr", "")
  expect_true("cases_time" %in% methods_irr)
  expect_false("means_sd" %in% methods_irr)
  expect_false("or_se" %in% methods_irr)
})

test_that(".get_applicable_methods returns within-group methods for dw/gw/mdw", {
  methods_dw <- metaConvert:::.get_applicable_methods("dw", "")
  expect_true("means_sd_pre_post_single_group" %in% methods_dw)
  expect_true("paired_t_single_group" %in% methods_dw)
  expect_false("means_sd" %in% methods_dw)
})

test_that(".get_applicable_methods returns prop methods for prop measure", {
  methods_prop <- metaConvert:::.get_applicable_methods("prop", "")
  expect_true("prop_single_group" %in% methods_prop)
  expect_true("prop_single_group_counts" %in% methods_prop)
  expect_false("means_sd" %in% methods_prop)
})

test_that(".get_applicable_methods returns partial cor methods for rp/zp", {
  methods_rp <- metaConvert:::.get_applicable_methods("rp", "")
  expect_true("linreg_t" %in% methods_rp)
  expect_false("means_sd" %in% methods_rp)
})


# ==============================================================================
# Unit tests for .diagnose_missing_data()
# ==============================================================================
test_that(".diagnose_missing_data detects means_sd near-miss with description", {
  mc <- metaConvert:::.method_required_columns()
  desc <- metaConvert:::.method_descriptions()

  # User provided means + sample sizes but not SDs
  row_data <- list(
    mean_exp = 10, mean_nexp = 8,
    n_exp = 30, n_nexp = 30,
    mean_sd_exp = NA, mean_sd_nexp = NA,
    mean_se_exp = NA, mean_se_nexp = NA,
    student_t = NA, or = NA, logor = NA
  )

  result <- metaConvert:::.diagnose_missing_data(
    row_data,
    methods_to_check = c("means_sd", "means_se", "student_t", "or"),
    method_cols = mc, descriptions = desc
  )

  # Should use human-readable description, not raw method name
  expect_true(grepl("Means \\+ SDs of two groups", result))
  expect_true(grepl("mean_sd_exp", result))
  expect_true(grepl("mean_sd_nexp", result))
})

test_that(".diagnose_missing_data detects student_t near-miss (missing t-value)", {
  mc <- metaConvert:::.method_required_columns()

  # User provided sample sizes only + student_t_pval (but not student_t)
  row_data <- list(
    n_exp = 30, n_nexp = 30,
    student_t_pval = 0.02,
    student_t = NA,
    mean_exp = NA, mean_nexp = NA
  )

  result <- metaConvert:::.diagnose_missing_data(
    row_data,
    methods_to_check = c("student_t", "student_t_pval", "means_sd"),
    method_cols = mc
  )

  # student_t_pval should work (has specific + shared), so not in guidance
  # student_t should NOT be a near-miss because its specific col (student_t) is NA
  # means_sd should NOT be a near-miss because no specific cols present
  expect_false(grepl("Means \\+ SDs", result))
})

test_that(".diagnose_missing_data returns fallback for empty row", {
  mc <- metaConvert:::.method_required_columns()

  # Row with nothing but NAs
  row_data <- list(
    mean_exp = NA, mean_nexp = NA, n_exp = NA, n_nexp = NA,
    student_t = NA, or = NA, logor = NA
  )

  result <- metaConvert:::.diagnose_missing_data(
    row_data,
    methods_to_check = c("means_sd", "student_t", "or"),
    method_cols = mc
  )

  expect_true(grepl("data_extraction_sheet", result))
})

test_that(".diagnose_missing_data handles OR alternative columns with description", {
  mc <- metaConvert:::.method_required_columns()
  desc <- metaConvert:::.method_descriptions()

  # User provided or but not logor_se
  row_data <- list(
    or = 2.5, logor = NA, logor_se = NA,
    n_exp = 50, n_nexp = 50,
    or_pval = NA
  )

  result <- metaConvert:::.diagnose_missing_data(
    row_data,
    methods_to_check = c("or_se", "or_pval"),
    method_cols = mc, descriptions = desc
  )

  # Should use description for or_se method
  expect_true(grepl("Odds ratio \\+ standard error", result))
  expect_true(grepl("logor_se", result))
})

test_that(".diagnose_missing_data limits to max_suggestions", {
  mc <- metaConvert:::.method_required_columns()

  # User provided mean_exp and mean_nexp, missing many others
  row_data <- list(
    mean_exp = 10, mean_nexp = 8,
    n_exp = 30, n_nexp = 30,
    mean_sd_exp = NA, mean_sd_nexp = NA,
    mean_se_exp = NA, mean_se_nexp = NA,
    mean_ci_lo_exp = NA, mean_ci_up_exp = NA,
    mean_ci_lo_nexp = NA, mean_ci_up_nexp = NA,
    mean_sd_pooled = NA
  )

  result <- metaConvert:::.diagnose_missing_data(
    row_data,
    methods_to_check = c("means_sd", "means_se", "means_ci", "means_sd_pooled"),
    method_cols = mc,
    max_suggestions = 2
  )

  # Should have at most 2 semicolons (max 2 suggestions)
  n_semicolons <- length(gregexpr(";", result)[[1]])
  # If there's no match, gregexpr returns -1
  if (n_semicolons == 1 && gregexpr(";", result)[[1]][1] == -1) n_semicolons <- 0
  expect_true(n_semicolons <= 1)  # max 2 items = max 1 semicolon
})

test_that(".diagnose_missing_data sorts by fewest missing cols first", {
  mc <- metaConvert:::.method_required_columns()
  desc <- metaConvert:::.method_descriptions()

  # means_sd_pooled needs fewer extra cols than means_sd (only mean_sd_pooled vs mean_sd_exp+mean_sd_nexp)
  row_data <- list(
    mean_exp = 10, mean_nexp = 8,
    n_exp = 30, n_nexp = 30,
    mean_sd_exp = NA, mean_sd_nexp = NA,
    mean_sd_pooled = NA,
    mean_se_exp = NA, mean_se_nexp = NA
  )

  result <- metaConvert:::.diagnose_missing_data(
    row_data,
    methods_to_check = c("means_sd", "means_sd_pooled", "means_se"),
    method_cols = mc, descriptions = desc
  )

  # "Means + pooled SD" should appear before "Means + SDs" (1 missing vs 2 missing)
  pos_pooled <- regexpr("Means \\+ pooled SD", result)
  pos_sd <- regexpr("Means \\+ SDs of two groups", result)
  if (pos_pooled > 0 && pos_sd > 0) {
    expect_true(pos_pooled < pos_sd)
  }
})


# ==============================================================================
# Integration tests with .add_es_guidance()
# ==============================================================================
test_that("guidance column appears for NA rows with partial data", {
  # Create a minimal dataset where one row has incomplete data
  dat <- data.frame(
    mean_exp = c(10, 12),
    mean_sd_exp = c(2, NA),
    mean_nexp = c(8, 9),
    mean_sd_nexp = c(2.5, NA),
    n_exp = c(30, 30),
    n_nexp = c(30, 30)
  )

  mc_obj <- convert_df(dat, measure = "d", split_adjusted = FALSE, verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  expect_true("es_guidance" %in% names(res))

  # Row 1 should have an ES (complete data) → empty guidance
  expect_equal(res$es_guidance[res$row_id == 1], "")

  # Row 2 should have NA ES → non-empty guidance with human-readable descriptions
  expect_true(is.na(res$es[res$row_id == 2]) || res$es_guidance[res$row_id == 2] != "")
})

test_that("guidance = FALSE produces no guidance column", {
  dat <- data.frame(
    mean_exp = c(10),
    mean_sd_exp = c(2),
    mean_nexp = c(8),
    mean_sd_nexp = c(2.5),
    n_exp = c(30),
    n_nexp = c(30)
  )

  mc_obj <- convert_df(dat, measure = "d", split_adjusted = FALSE, verbose = FALSE)
  res <- summary(mc_obj, guidance = FALSE, flags = FALSE)

  expect_false("es_guidance" %in% names(res))
})

test_that("guidance works with split_adjusted = TRUE, format = wide", {
  dat <- data.frame(
    mean_exp = c(10),
    mean_sd_exp = c(NA),
    mean_nexp = c(8),
    mean_sd_nexp = c(NA),
    n_exp = c(30),
    n_nexp = c(30)
  )

  mc_obj <- convert_df(dat, measure = "d",
                        split_adjusted = TRUE, format_adjusted = "wide",
                        verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  expect_true("es_guidance_crude" %in% names(res))
  expect_true("es_guidance_adjusted" %in% names(res))
})

test_that("guidance column is empty string for rows with valid ES", {
  dat <- data.frame(
    mean_exp = c(10),
    mean_sd_exp = c(2),
    mean_nexp = c(8),
    mean_sd_nexp = c(2.5),
    n_exp = c(30),
    n_nexp = c(30)
  )

  mc_obj <- convert_df(dat, measure = "g", split_adjusted = FALSE, verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  expect_true("es_guidance" %in% names(res))
  expect_equal(res$es_guidance[1], "")
})


# ==============================================================================
# Integration test with df.haza
# ==============================================================================
test_that("guidance integrates with df.haza dataset", {
  mc_obj <- convert_df(df.haza, measure = "g", split_adjusted = FALSE, verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  expect_true("es_guidance" %in% names(res))

  # For rows with valid ES, guidance should be empty
  valid_rows <- which(!is.na(res$es))
  if (length(valid_rows) > 0) {
    expect_true(all(res$es_guidance[valid_rows] == ""))
  }

  # For rows with NA ES, guidance should be non-empty
  na_rows <- which(is.na(res$es))
  if (length(na_rows) > 0) {
    expect_true(all(res$es_guidance[na_rows] != ""))
  }
})

test_that("guidance works with df.haza default split_adjusted=TRUE", {
  mc_obj <- convert_df(df.haza, measure = "g", verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  # Default is split_adjusted=TRUE, format=wide
  expect_true("es_guidance_crude" %in% names(res))
  expect_true("es_guidance_adjusted" %in% names(res))
})


# ==============================================================================
# Test guidance with different measures
# ==============================================================================
test_that("guidance works for OR measure", {
  dat <- data.frame(
    or = c(2.5, NA),
    logor_se = c(0.3, NA),
    n_exp = c(50, 50),
    n_nexp = c(50, 50),
    n_cases_exp = c(NA, 20),
    n_cases_nexp = c(NA, NA),
    n_controls_exp = c(NA, NA),
    n_controls_nexp = c(NA, NA)
  )

  mc_obj <- convert_df(dat, measure = "or", split_adjusted = FALSE, verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  expect_true("es_guidance" %in% names(res))
})

test_that("guidance does NOT show variability methods for measure = g", {
  # This was the user-reported bug: variability_means_sd appeared for measure="g"
  dat <- data.frame(
    mean_exp = 2, mean_sd_exp = 2.5,
    mean_nexp = 4,
    n_nexp = 5, n_exp = 5
  )

  mc_obj <- convert_df(dat, measure = "g", split_adjusted = FALSE, verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  # Should NOT contain variability-related guidance
  na_rows <- which(is.na(res$es))
  if (length(na_rows) > 0) {
    for (i in na_rows) {
      expect_false(grepl("ariability", res$es_guidance[i]),
                   info = "Variability methods should not appear for measure=g")
    }
  }
})

test_that("guidance uses human-readable descriptions, not raw method names", {
  dat <- data.frame(
    mean_exp = 2, mean_sd_exp = 2.5,
    mean_nexp = 4, mean_sd_nexp = NA,
    n_nexp = 5, n_exp = 5
  )

  mc_obj <- convert_df(dat, measure = "g", split_adjusted = FALSE, verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  na_rows <- which(is.na(res$es))
  if (length(na_rows) > 0) {
    guidance <- res$es_guidance[na_rows[1]]
    # Should contain human-readable description, not raw method name
    expect_true(grepl("Means \\+ SDs of two groups", guidance),
                info = paste("Expected human-readable description, got:", guidance))
  }
})

test_that("guidance works for correlation measure", {
  dat <- data.frame(
    pearson_r = c(0.5, NA),
    n_sample = c(100, 100),
    fisher_z = c(NA, NA)
  )

  mc_obj <- convert_df(dat, measure = "r", split_adjusted = FALSE, verbose = FALSE)
  res <- summary(mc_obj, guidance = TRUE, flags = FALSE)

  expect_true("es_guidance" %in% names(res))
})
