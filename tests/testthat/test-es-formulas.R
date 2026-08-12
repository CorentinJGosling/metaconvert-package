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
