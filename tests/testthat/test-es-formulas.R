## Tests for es_formulas(): sensitivity to the CONVERSION FORMULA (or_to_rr, cor_to_smd,
## pre_post_to_smd, ...), as distinct from the choice of source statistics covered by
## test-main-es-routes.R.

test_that("es_formulas() reports every or_to_rr formula for one OR", {
  dat <- data.frame(study_id = c("A", "B"),
                    or = c(2.5, 1.8), or_ci_lo = c(1.6, 1.2),
                    or_ci_up = c(3.9, 2.7),
                    n_exp = 100, n_nexp = 100, baseline_risk = 0.2)
  fx <- es_formulas(convert_df(dat, measure = "rr", verbose = FALSE), verbose = FALSE)

  expect_s3_class(fx, "metaConvert_formulas")
  expect_true("or_to_rr" %in% fx$parameter)
  # all five documented formulas are represented
  expect_setequal(fx$formula[fx$parameter == "or_to_rr"],
                  c("metaumbrella_cases", "metaumbrella_exp", "transpose",
                    "grant", "dipietrantonj"))
  # exactly one active formula per (comparison, parameter)
  n_active <- tapply(fx$is_active, paste(fx$row_id, fx$parameter), sum)
  expect_true(all(n_active == 1))
  # transpose returns RR = OR by construction
  tr <- fx[fx$parameter == "or_to_rr" & fx$formula == "transpose" & fx$row_id == 1, ]
  expect_equal(tr$es, 2.5, tolerance = 1e-6)
})

test_that("es_formulas() records a real spread and keeps it out of the flags", {
  dat <- data.frame(study_id = c("C1", "C2"),
                    pearson_r = c(.30, .45), n_sample = c(100, 120))
  obj <- convert_df(dat, measure = "g", verbose = FALSE)
  fx <- es_formulas(obj, verbose = FALSE)

  expect_true("cor_to_smd" %in% fx$parameter)
  smry <- attr(fx, "summary")
  expect_true(all(smry$spread[!is.na(smry$spread)] > 0))
  # the label is [FORMULA], never [DISCORDANT]
  expect_false(any(grepl("DISCORDANT", fx$formula, fixed = TRUE)))
})

test_that("summary(formulas = TRUE) adds a [FORMULA] message and attaches the table", {
  dat <- data.frame(study_id = c("C1", "C2"),
                    pearson_r = c(.30, .45), n_sample = c(100, 120))
  obj <- convert_df(dat, measure = "g", verbose = FALSE)

  s_off <- suppressMessages(summary(obj))
  s_on  <- suppressMessages(summary(obj, formulas = TRUE))

  # opt-in: default output untouched
  expect_false(any(grepl("[FORMULA]", s_off$flags_crude, fixed = TRUE)))
  expect_null(attr(s_off, "formulas"))

  expect_true(any(grepl("[FORMULA]", s_on$flags_crude, fixed = TRUE)))
  expect_s3_class(attr(s_on, "formulas"), "metaConvert_formulas")
  # effect sizes themselves are unchanged by the disclosure
  expect_equal(as.numeric(s_off$es_crude), as.numeric(s_on$es_crude))
})

test_that("the [FORMULA] message never contains the flag separator", {
  dat <- data.frame(study_id = c("C1", "C2"),
                    pearson_r = c(.30, .45), n_sample = c(100, 120))
  fx <- es_formulas(convert_df(dat, measure = "g", verbose = FALSE), verbose = FALSE)
  tok <- metaConvert:::.formula_tokens(fx)
  expect_true(length(tok) > 0)
  # '; ' is how .flag_es_quality separates merged messages; an embedded one would
  # create an unlabelled fragment
  expect_false(any(grepl("; ", tok, fixed = TRUE)))
  expect_true(all(grepl("^\\[FORMULA\\]", tok)))
})

test_that("scale-changing parameters report values but no range", {
  dat <- data.frame(study_id = c("P1", "P2"),
                    cronbach_alpha = c(0.85, 0.90),
                    n_sample = 200, n_items = 10)
  fx <- es_formulas(convert_df(dat, measure = "alpha", verbose = FALSE),
                    verbose = FALSE)
  expect_true(all(fx$scale_change))
  # bonett returns ln(1 - alpha) and raw returns alpha: direct comparison is invalid
  expect_true(all(is.na(fx$deviation)))
  expect_true(all(is.na(attr(fx, "summary")$spread)))
  # the values themselves are still reported
  expect_equal(sort(unique(fx$formula)), c("bonett", "raw"))
  expect_true(all(is.finite(fx$es)))
})

test_that("es_formulas() validates its inputs", {
  dat <- data.frame(pearson_r = .3, n_sample = 100)
  obj <- convert_df(dat, measure = "g", verbose = FALSE)
  expect_error(es_formulas(dat), "must be an object of class 'metaConvert'")
  expect_error(es_formulas(obj, parameters = "not_a_parameter"), "Unknown conversion parameter")
})

test_that("parameters that change nothing are excluded", {
  # means and SDs only: no conversion formula contributes to g here
  dat <- data.frame(study_id = c("M1", "M2"), n_exp = 50, n_nexp = 50,
                    mean_exp = c(10, 11), mean_nexp = 8,
                    mean_sd_exp = 3, mean_sd_nexp = 3)
  fx <- es_formulas(convert_df(dat, measure = "g", verbose = FALSE), verbose = FALSE)
  expect_false("cor_to_smd" %in% fx$parameter)
  expect_false("or_to_rr" %in% fx$parameter)
})

test_that("es_formulas() runs whatever the layout of the object", {
  # the re-evaluations are fed the checked data, which already carries the
  # result columns .check_data() appends; a comparison that cannot be estimated
  # used to abort them all, and every parameter was then silently skipped
  dat <- data.frame(study_id = c("S1", "S2"), pearson_r = c(.3, NA),
                    n_sample = c(100, NA))
  ref <- es_formulas(convert_df(dat, measure = "g", verbose = FALSE), verbose = FALSE)
  expect_gt(nrow(ref), 0)

  for (obj in list(convert_df(dat, measure = "g", split_adjusted = FALSE, verbose = FALSE),
                   convert_df(dat, measure = "g", format_adjusted = "long", verbose = FALSE),
                   convert_df(dat, measure = "g", main_es = FALSE, verbose = FALSE))) {
    fx <- es_formulas(obj, verbose = FALSE)
    expect_equal(nrow(fx), nrow(ref))
    expect_equal(fx$es, ref$es)
  }
})

test_that("each documented alias of a formula is recognised as the active one", {
  dat <- data.frame(study_id = "A", n_exp = 30, n_nexp = 32,
                    mean_exp = 10, mean_sd_exp = 2,
                    mean_nexp = 8.5, mean_sd_nexp = 2.6)
  for (v in c("pooled", "glass", "control", "glass_robust", "control_robust")) {
    obj <- suppressMessages(convert_df(dat, measure = "g", smd_denom = v, verbose = FALSE))
    fx <- es_formulas(obj, parameters = "smd_denom", verbose = FALSE)
    expect_equal(sum(fx$is_active), 1L, info = v)
    expect_false(is.na(attr(fx, "summary")$active[1]), info = v)
    expect_false(any(grepl("applied: 'NA'", metaConvert:::.formula_tokens(fx), fixed = TRUE)),
                 info = v)
    # the reference of `deviation` is the formula that was applied
    expect_equal(fx$deviation[fx$is_active], 0)
  }
})

test_that("a per-row conversion column is honoured comparison by comparison", {
  dat <- data.frame(study_id = c("A", "B"), or = 2.5, or_ci_lo = 1.6, or_ci_up = 3.9,
                    n_exp = 100, n_nexp = 100, n_cases = 120, n_controls = 80,
                    baseline_risk = 0.2, or_to_rr = c("transpose", "grant"))
  obj <- convert_df(dat, measure = "rr", split_adjusted = FALSE, verbose = FALSE)
  es <- as.numeric(suppressMessages(summary(obj))$es)
  fx <- es_formulas(obj, verbose = FALSE)

  expect_equal(fx$formula[fx$is_active], c("transpose", "grant"))
  # the active formula reproduces the estimate convert_df() returned for that row
  expect_equal(fx$es[fx$is_active], round(es, 3))
  expect_equal(fx$deviation[fx$is_active], c(0, 0))
  tok <- metaConvert:::.formula_tokens(fx)
  expect_true(grepl("applied: 'transpose'", tok[["1"]], fixed = TRUE))
  expect_true(grepl("applied: 'grant'", tok[["2"]], fixed = TRUE))
})

test_that("a scale-changing parameter is reported only where it applied", {
  # df.psychom mixes alpha, icc, correlation and proportion rows: icc_to_es must
  # not be reported for the comparisons that carry no ICC
  obj <- convert_df(df.psychom, measure = "icc", verbose = FALSE)
  s <- suppressMessages(summary(obj, formulas = TRUE))
  fx <- es_formulas(obj, verbose = FALSE)
  estimated <- !is.na(suppressWarnings(as.numeric(as.character(s$es_crude))))

  expect_true(all(is.finite(fx$es)))
  expect_setequal(unique(fx$row_id), s$row_id[estimated])
  # no comparison without an effect size carries a [FORMULA] message
  expect_false(any(grepl("[FORMULA]", s$flags_crude[!estimated], fixed = TRUE)))

  # and a parameter that cannot touch the requested measure is not reported
  dat <- data.frame(n_exp = 30, n_nexp = 32, mean_exp = 10, mean_sd_exp = 2,
                    mean_nexp = 8.5, mean_sd_nexp = 2.6)
  fx2 <- es_formulas(convert_df(dat, measure = "g", verbose = FALSE),
                     parameters = "alpha_to_es", verbose = FALSE)
  expect_equal(nrow(fx2), 0)
})

test_that("aliased pre-post standardisers are not evaluated twice", {
  # 'cooper' is rewritten to 'morris_drm' before any computation
  expect_false("cooper" %in% metaConvert:::.formula_parameters()$pre_post_to_smd)
  dat <- data.frame(study_id = "P", n_exp = 36, n_nexp = 35,
                    mean_pre_exp = 10, mean_pre_sd_exp = 3,
                    mean_exp = 12, mean_sd_exp = 3.2,
                    mean_pre_nexp = 10.2, mean_pre_sd_nexp = 3.1,
                    mean_nexp = 11, mean_sd_nexp = 3.3,
                    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5)
  obj <- convert_df(dat, measure = "g", es_selected = "hierarchy",
                    hierarchy = "means_sd_pre_post", verbose = FALSE)
  fx <- es_formulas(obj, parameters = "pre_post_to_smd", verbose = FALSE)
  expect_setequal(fx$formula,
                  c("bonett", "morris_dz", "morris_drm", "morris_dav"))
  expect_equal(attr(fx, "summary")$n_formulas, 4L)
  # a user who set the alias still gets one active formula
  obj2 <- convert_df(dat, measure = "g", pre_post_to_smd = "cooper",
                     es_selected = "hierarchy", hierarchy = "means_sd_pre_post",
                     verbose = FALSE)
  fx2 <- es_formulas(obj2, parameters = "pre_post_to_smd", verbose = FALSE)
  expect_equal(fx2$formula[fx2$is_active], "morris_drm")
})

test_that("a formula that alters only the standard error is still reported", {
  # smd_var selects the sampling-variance formula: the point estimate is common
  # to both options, the weight the study receives is not
  dat <- data.frame(mean_exp = 10, mean_sd_exp = 2, n_exp = 8,
                    mean_nexp = 8.5, mean_sd_nexp = 2.6, n_nexp = 9)
  fx <- es_formulas(convert_df(dat, measure = "g", verbose = FALSE),
                    parameters = "smd_var", verbose = FALSE)
  expect_setequal(fx$formula, c("borenstein", "hedges_olkin"))
  expect_equal(diff(range(fx$es)), 0)
  expect_gt(attr(fx, "summary")$se_spread, 0)
  expect_true(grepl("standard error", metaConvert:::.formula_tokens(fx)[[1]], fixed = TRUE))
})

test_that("the number of formulas that produced an estimate is reported", {
  # 'metaumbrella_cases' and 'dipietrantonj' need the case/control margins
  dat <- data.frame(or = 2.5, or_ci_lo = 1.6, or_ci_up = 3.9,
                    n_exp = 100, n_nexp = 100, baseline_risk = 0.2)
  fx <- es_formulas(convert_df(dat, measure = "rr", verbose = FALSE), verbose = FALSE)
  smry <- attr(fx, "summary")
  expect_equal(smry$n_formulas, 5L)
  expect_equal(smry$n_estimated, 3L)
  expect_false(smry$active_estimated)
  tok <- metaConvert:::.formula_tokens(fx)
  expect_true(grepl("across 3 of 5 formulas", tok[[1]], fixed = TRUE))
  expect_true(grepl("yields no estimate here", tok[[1]], fixed = TRUE))
})
