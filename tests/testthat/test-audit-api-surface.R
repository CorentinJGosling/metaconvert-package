# Red half of red-green for the AUDIT-2026-08-28 "API surface" findings.
# Every test in this file asserts the CORRECT behaviour and therefore FAILS
# against the working tree as of the audit. Do not weaken an assertion to make
# it pass; fix the code it names.

# Values enumerated by man/data_extraction_sheet.Rd's \item{measure}.
# Kept here as a fallback for the (installed) case where man/ is not on disk.
.audit_documented_des_measures <- c(
  "d", "g", "md", "dw", "gw", "mdw", "r", "z", "rp", "zp", "or", "rr",
  "irr", "hr", "logor", "logrr", "logirr", "loghr", "logvr", "logcvr",
  "nnt", "rd", "prop", "alpha", "omega", "icc", "all"
)

# Locate the package source tree (present when running from the repo, absent
# under R CMD check on the built tarball).
.audit_src_file <- function(...) {
  p <- file.path("..", "..", ...)
  if (file.exists(p)) normalizePath(p) else NA_character_
}

# Pull the quoted tokens out of one \item{<arg>}{...} block of an Rd file.
.audit_rd_item_tokens <- function(rd_path, arg) {
  txt <- paste(readLines(rd_path, warn = FALSE), collapse = " ")
  txt <- gsub("[[:space:]]+", " ", txt)
  pat <- paste0("\\\\item\\{", arg, "\\}\\{[^}]*\\}")
  item <- regmatches(txt, regexpr(pat, txt))
  if (length(item) == 0L) return(character(0))
  unique(gsub('"', "", regmatches(item, gregexpr('"[A-Za-z]+"', item))[[1]]))
}


# AUDIT-2026-08-28 lines 621-630 and 654-663.
# data_extraction_sheet() has no `else` branch in its measure dispatch, so the
# four log-scale measure names its own Rd documents leave `dat` unassigned and
# the function dies at `dup <- duplicated(colnames(dat))` with the internal
# error "object 'dat' not found". The preferred fix (audit line 630, "Option (a)
# is preferable since the columns genuinely exist") aliases them to the
# identical-column non-log measures, so that is what is asserted here.
test_that("AUDIT-blocker: data_extraction_sheet() serves the documented log-scale measures", {
  aliases <- c(logor = "or", logrr = "rr", logirr = "irr", loghr = "hr")

  for (log_m in names(aliases)) {
    plain_m <- aliases[[log_m]]

    plain <- data_extraction_sheet(
      measure = plain_m, extension = "data.frame", verbose = FALSE
    )

    got <- tryCatch(
      data_extraction_sheet(
        measure = log_m, extension = "data.frame", verbose = FALSE
      ),
      error = function(e) e
    )

    # 1. It must not die with the internal "object 'dat' not found"
    #    (locale-independent check: the condition must not be an error at all).
    expect_false(
      inherits(got, "error"),
      info = paste0(
        "data_extraction_sheet(measure = \"", log_m, "\") errored: ",
        if (inherits(got, "error")) conditionMessage(got) else ""
      )
    )

    # 2. logor/logrr/logirr/loghr take exactly the same input columns as
    #    or/rr/irr/hr, so the sheet must be the same sheet.
    if (!inherits(got, "error")) {
      expect_identical(colnames(got), colnames(plain))
    }
  }
})


# AUDIT-2026-08-28 line 626.
# The \usage default vector for `measure` prints 16 values while the same help
# page's \arguments enumerates 27; the default vector must cover every value the
# page documents as valid.
test_that("AUDIT-doc: data_extraction_sheet()'s measure default vector covers its documented set", {
  rd <- .audit_src_file("man", "data_extraction_sheet.Rd")
  documented <- if (!is.na(rd)) {
    .audit_rd_item_tokens(rd, "measure")
  } else {
    .audit_documented_des_measures
  }
  expect_gt(length(documented), 0L)

  default_vec <- eval(formals(data_extraction_sheet)$measure)

  expect_setequal(default_vec, documented)
})


# AUDIT-2026-08-28 lines 1139-1150 (incl. the Verifier correction).
# man/metaConvert-package.Rd points users at see_input_data() to see the inputs
# for "the 22 effect size measures", but the by-name guard accepts only 11 of
# them; and the guard's own error message enumerates a stale exclusion set
# ("rd", "alpha", "icc", "prop", "hr") that omits omega, dw, gw, mdw, rp, zp.
test_that("AUDIT-doc: see_input_data()'s coverage claim and guard message match what it accepts", {
  supported <- eval(formals(see_input_data)$measure)
  n_supported_measures <- length(setdiff(supported, "all"))

  # The 22 pipeline measures the landing page counts.
  documented <- c(
    "d", "g", "md", "dw", "gw", "mdw", "r", "z", "rp", "zp", "or", "rr",
    "irr", "hr", "logvr", "logcvr", "nnt", "rd", "prop", "alpha", "omega", "icc"
  )
  uncovered <- setdiff(documented, supported)

  # (a) The guard's message must not present an exclusion list that leaves out
  #     measures which are in fact excluded. Ask for one that the current
  #     message does not name at all.
  msg <- tryCatch(see_input_data(measure = "dw"), error = conditionMessage)
  expect_type(msg, "character")

  quoted <- gsub('"', "", regmatches(msg, gregexpr('"[A-Za-z]+"', msg))[[1]])
  enumerated_exclusions <- setdiff(intersect(quoted, documented), c("dw", supported))

  if (length(enumerated_exclusions) > 0L) {
    # It enumerates; the enumeration must be complete.
    expect_setequal(enumerated_exclusions, setdiff(uncovered, "dw"))
  } else {
    succeed()
  }

  # (b) The package landing page must not promise more measures than
  #     see_input_data() actually serves.
  rd <- .audit_src_file("man", "metaConvert-package.Rd")
  if (!is.na(rd)) {
    txt <- gsub("[[:space:]]+", " ", paste(readLines(rd, warn = FALSE), collapse = " "))
    claim <- regmatches(
      txt, regexpr("[0-9]+ effect size measures[^.]*see_input_data", txt)
    )
    if (length(claim) == 1L) {
      claimed_n <- as.numeric(regmatches(claim, regexpr("[0-9]+", claim)))
      expect_lte(claimed_n, n_supported_measures)
    }
  }
})


# AUDIT-2026-08-28 lines 478-487.
# .parameters_for_measure() omits conversion parameters that demonstrably move
# the estimate -- cor_to_smd and pre_post_to_smd for the `ratio` class,
# smd_denom for the `cor_like` class -- so es_formulas() returns 0 rows and
# prints the explicit all-clear "No conversion formula altered any effect size".
# The expected numbers below were derived independently of es_formulas(), by
# running convert_df() directly once per scalar setting.
test_that("AUDIT-blocker: es_formulas() evaluates every conversion parameter that moves the estimate", {
  # returns NA_real_ for a formula the table does not contain, so a missing
  # parameter shows up as a clean failure rather than a subscript error
  pull <- function(tab, param, formula) {
    tab <- as.data.frame(tab)
    tab <- tab[tab$parameter == param & tab$formula == formula, , drop = FALSE]
    if (nrow(tab) != 1L) return(NA_real_)
    tab$es
  }

  # (a) the parameter lists themselves. es_formulas(object) with no `parameters`
  #     argument iterates exactly .parameters_for_measure(measure), so these
  #     three memberships ARE the defect.
  ratio_params <- metaConvert:::.parameters_for_measure("logor")
  expect_true("cor_to_smd" %in% ratio_params)
  expect_true("pre_post_to_smd" %in% ratio_params)
  expect_true("smd_denom" %in% metaConvert:::.parameters_for_measure("r"))

  # (b) case 1, end to end on the DEFAULT parameter selection: an OR pool built
  #     from a correlation. cor_to_smd moves the OR 4.026 -> 6.221 (+55%), yet
  #     es_formulas() currently reports 0 rows.
  #     The viechtbauer value is exp(transf.rtod(r, n1i, n2i) / J(N - 2) * pi/sqrt(3)):
  #     this row SUPPLIES n_exp = n_nexp = 60, so the r -> d map is inverted at the
  #     finite-sample h = m/n1 + m/n2 the forward d -> r map uses (see .rtod_delta()),
  #     not at the balanced-arms constant 4 that gave the former pin of 4.0734.
  d1 <- data.frame(study_id = "A", pearson_r = 0.45, n_sample = 120,
                   n_exp = 60, n_nexp = 60)
  o1 <- suppressMessages(convert_df(d1, measure = "or", verbose = FALSE,
                                    split_adjusted = FALSE))
  f1 <- as.data.frame(es_formulas(o1, verbose = FALSE, digits = 12))
  expect_gt(nrow(f1), 0L)
  expect_true("cor_to_smd" %in% f1$parameter)
  expect_equal(pull(f1, "cor_to_smd", "viechtbauer"), 4.02584061437363,
               tolerance = 1e-8)
  expect_equal(pull(f1, "cor_to_smd", "viechtbauer"),
               exp(metafor::transf.rtod(0.45, 60, 60) / metaConvert:::.d_j(118) *
                     pi / sqrt(3)),
               tolerance = 1e-10)
  expect_equal(pull(f1, "cor_to_smd", "cooper"), 6.22117518089525,
               tolerance = 1e-8)

  # (c) case 2: OR built from a paired t. pre_post_to_smd moves the OR
  #     1.999 -> 2.169 (+8.5%). Named explicitly here to keep the file fast --
  #     assertion (a) is what pins it into the default selection.
  d2 <- data.frame(study_id = "A", paired_t_exp = 3.5, paired_t_nexp = 1.2,
                   n_exp = 30, n_nexp = 32,
                   r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6)
  o2 <- suppressMessages(convert_df(d2, measure = "or", verbose = FALSE,
                                    split_adjusted = FALSE))
  f2 <- as.data.frame(es_formulas(o2, parameters = "pre_post_to_smd",
                                  verbose = FALSE, digits = 12))
  expect_equal(pull(f2, "pre_post_to_smd", "morris_drm"), 1.99876292762459,
               tolerance = 1e-8)
  expect_equal(pull(f2, "pre_post_to_smd", "morris_dz"), 2.16900892932206,
               tolerance = 1e-8)

  # (d) case 3: r built from means/SDs. smd_denom moves r 0.327 -> 0.444 (+36%),
  #     by far the largest influence on that pool, and it is the one parameter
  #     es_formulas() leaves out of its 4-row table for this dataset.
  d3 <- data.frame(study_id = "A", mean_exp = 25, mean_sd_exp = 7,
                   mean_nexp = 22, mean_sd_nexp = 4, n_exp = 40, n_nexp = 45)
  o3 <- suppressMessages(convert_df(d3, measure = "r", verbose = FALSE,
                                    split_adjusted = FALSE))
  f3 <- as.data.frame(es_formulas(o3, parameters = "smd_denom",
                                  verbose = FALSE, digits = 12))
  expect_equal(pull(f3, "smd_denom", "pooled"), 0.326983425572447,
               tolerance = 1e-8)
  expect_equal(pull(f3, "smd_denom", "glass"), 0.444444657087177,
               tolerance = 1e-8)
})


# AUDIT-2026-08-28 lines 489-498 (and the second half of 1087-1098).
# omega_to_es is a convert_df() argument with three options and is already
# recorded in attr(res, "conversion_args"), but es_formulas() has no entry for
# it, so an omega object gets an explicit "no conversion formula can influence
# the 'omega' measure" all-clear and `parameters = "omega_to_es"` errors.
test_that("AUDIT-blocker: es_formulas() sees omega_to_es", {
  d <- data.frame(study_id = "s1", omega = 0.86, omega_ci_lo = 0.80,
                  omega_ci_up = 0.91, n_sample = 200, n_items = 10)
  o <- suppressMessages(convert_df(d, measure = "omega", verbose = FALSE,
                                   split_adjusted = FALSE))

  # (a) omega_to_es must be a parameter es_formulas() knows about
  expect_true("omega_to_es" %in% names(metaConvert:::.formula_parameters()))
  expect_true("omega_to_es" %in% metaConvert:::.parameters_for_measure("omega"))

  # (b) ln(1 - w), w and 1 - (1 - w)^(1/3) are three different scales, so
  #     omega_to_es must be registered as scale-changing (deviation/spread NA)
  expect_true("omega_to_es" %in% metaConvert:::.scale_changing_parameters())

  # (c) naming it explicitly must not be rejected as unknown
  expect_error(es_formulas(o, parameters = "omega_to_es", verbose = FALSE), NA)

  # (d) the default run must report the three formulas, on their own scales.
  #     Expected values derived by hand from omega = 0.86:
  #       bonett          = log(1 - 0.86)            = -1.9661128563728
  #       raw             = 0.86
  #       hakstian_whalen = 1 - (1 - 0.86)^(1/3)     =  0.4807505898149
  f <- as.data.frame(es_formulas(o, verbose = FALSE, digits = 12))
  expect_true("omega_to_es" %in% f$parameter)

  rows <- f[f$parameter == "omega_to_es", , drop = FALSE]
  grab <- function(formula) {
    r <- rows[rows$formula == formula, , drop = FALSE]
    if (nrow(r) != 1L) NA_real_ else r$es
  }
  expect_setequal(rows$formula, c("bonett", "raw", "hakstian_whalen"))
  expect_equal(grab("bonett"), log(1 - 0.86), tolerance = 1e-8)
  expect_equal(grab("raw"), 0.86, tolerance = 1e-8)
  expect_equal(grab("hakstian_whalen"), 1 - (1 - 0.86)^(1 / 3), tolerance = 1e-8)
  expect_true(all(rows$scale_change))
  expect_true(all(is.na(rows$deviation)))
})


# AUDIT-2026-08-28 lines 500-509 and 1087-1098 (incl. the Verifier correction,
# which confirms the ICC half).
# attr(res, "conversion_args") omits alpha_se_source / omega_se_source /
# icc_agreement_se, so .run_with_parameter() silently reverts them to their
# defaults and es_formulas() prints standard errors the analysis refused.
test_that("AUDIT-blocker: es_formulas() honours the reliability SE-source switches", {
  three <- c("alpha_se_source", "omega_se_source", "icc_agreement_se")

  # --- alpha: alpha_se_source = "reported" with no reported SE/CI -> es is NA
  da <- data.frame(study_id = c("P1", "P2"), cronbach_alpha = c(0.85, 0.90),
                   n_sample = 200, n_items = 10)
  oa <- suppressMessages(convert_df(da, measure = "alpha",
                                    alpha_se_source = "reported",
                                    verbose = FALSE, split_adjusted = FALSE))

  expect_true(all(three %in% names(attr(oa, "conversion_args"))))
  expect_identical(attr(oa, "conversion_args")$alpha_se_source, "reported")

  suppressMessages(invisible(utils::capture.output(
    sa <- summary(oa, flags = FALSE, guidance = FALSE)
  )))
  expect_true(all(is.na(sa$es)))          # the analysis the user actually ran

  fa <- as.data.frame(es_formulas(oa, verbose = FALSE, digits = 12))
  # the sensitivity table must describe THAT analysis: no reported SE anywhere,
  # so either there is nothing to report or every SE is NA.
  expect_true(nrow(fa) == 0L || all(is.na(fa$se)))

  # --- icc: icc_agreement_se = "drop" -> the computed agreement SE is refused
  di <- data.frame(study_id = c("I1", "I2"), icc = c(0.80, 0.75),
                   n_sample = 100, n_measurements = 3)
  oi <- suppressMessages(convert_df(di, measure = "icc",
                                    icc_agreement_se = "drop",
                                    verbose = FALSE, split_adjusted = FALSE))

  expect_identical(attr(oi, "conversion_args")$icc_agreement_se, "drop")

  suppressMessages(invisible(utils::capture.output(
    si <- summary(oi, flags = FALSE, guidance = FALSE)
  )))
  expect_true(all(is.na(si$es)))

  fi <- as.data.frame(es_formulas(oi, verbose = FALSE, digits = 12))
  expect_true(nrow(fi) == 0L || all(is.na(fi$se)))
})
