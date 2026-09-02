# Red tests for the reliability findings of AUDIT-2026-08-28-findings.md.
#
# Every test in this file asserts the CORRECT behaviour, so it fails against the
# state of the tree at the time the audit was written and passes once the finding
# is repaired. Where the audit's "Verifier correction" field disagrees with the
# "Fix" field, the verifier's version is the one encoded here.
#
# Expected numbers are derived from the closed forms, not copied from the package:
#   inverse Spearman-Brown   rho_1 = rho_k / (k - (k - 1) rho_k)
#   ICC variance (Bonett)    se_raw = (1 - rho) * sqrt(2 (1 + (k-1) rho)^2 / (k (k-1) (n-1)))
#   alpha/omega (Bonett)     se     = sqrt(2 k / ((k - 1) (n - 2)))
#   Wan (2014) range/SD      xi(n)  = 2 * qnorm((n - 0.375) / (n + 0.25))

# Does `flags` carry an [INVALID] token that names `col`? V31's message mentions
# 'icc_ci_lo'/'icc_ci_up' in passing, so a bare grepl() on the column name is not
# enough to show the input bound was rejected.
.has_invalid_naming <- function(flags, col) {
  toks <- unlist(strsplit(paste(flags, collapse = "; "), "; ", fixed = TRUE))
  any(grepl("[INVALID]", toks, fixed = TRUE) &
        grepl(paste0("'", col, "'"), toks, fixed = TRUE))
}

# Locate a file in the package SOURCE tree (not the installed package). Returns
# NA_character_ when the sources are not beside the tests, e.g. R CMD check of a
# built tarball.
.audit_src <- function(rel) {
  cands <- c(
    tryCatch(testthat::test_path("..", "..", rel), error = function(e) NA_character_),
    file.path("..", "..", rel),
    rel
  )
  cands <- cands[!is.na(cands)]
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[[1]] else NA_character_
}


# AUDIT 970-981 (verifier correction at line 981). An ICC confidence-interval bound
# is never range-checked: .bounded_columns holds `icc` but not `icc_ci_lo`/`icc_ci_up`,
# and the CI triplet tests ordering/containment only. So an impossible bound (> 1)
# reaches the SE construction and pre-empts the closed-form (n, k) route with a
# finite, plausible-looking SE. On an average-measures row it is worse still: the
# Mobius map .icc_step_down() has a pole at k/(k-1), so the bound is sent negative
# and the SE is built from a sign-flipped interval, after the transposition guard
# has already been spent.
test_that("AUDIT-icc-ci-bounds: an out-of-range ICC CI bound is rejected, not turned into an SE", {
  # (1) the step-down itself must refuse a value past its pole (k - (k-1)x <= 0)
  expect_true(is.na(metaConvert:::.icc_step_down(1.20, 10)))
  # in-range bounds are untouched
  expect_equal(metaConvert:::.icc_step_down(0.80, 10), 0.8 / (10 - 9 * 0.8),
               tolerance = 1e-12)

  # (2) direct route (exported calculator; never sees .validate_input_data()).
  # The point estimate steps down as usual; the impossible interval yields no SE.
  direct <- es_from_icc(icc = 0.90, n_measurements = 10, icc_type = "ICC(2,k)",
                        icc_ci_lo = 0.80, icc_ci_up = 1.20, icc_to_es = "raw")
  expect_equal(direct$icc, 0.9 / (10 - 9 * 0.9), tolerance = 1e-9)
  expect_true(is.na(direct$icc_se))     # today: 0.4555477, from |(-1.5) - 0.2857| / (2 z)
  expect_true(is.na(direct$icc_ci_lo))
  expect_true(is.na(direct$icc_ci_up))

  # (3) pipeline, average-measures row. With the impossible bound neutralised the
  # SE must fall through to the closed form on the STEPPED-DOWN rho.
  d_avg <- data.frame(study_id = "A", icc = 0.90, n_measurements = 10,
                      icc_type = "ICC(2,k)", icc_ci_lo = 0.80, icc_ci_up = 1.20,
                      n_sample = 100)
  s_avg <- summary(convert_df(d_avg, measure = "icc", icc_to_es = "raw",
                              verbose = FALSE, split_adjusted = FALSE), flags = TRUE)
  rho1 <- 0.9 / (10 - 9 * 0.9)
  se_avg <- (1 - rho1) * sqrt(2 * (1 + 9 * rho1)^2 / (10 * 9 * (100 - 1)))
  expect_equal(s_avg$se[1], se_avg, tolerance = 1e-7)   # today: 0.4555477
  expect_true(.has_invalid_naming(s_avg$flags[1], "icc_ci_up"))

  # (4) the same wrong SE arises with no step-down and no sign flip at all
  d_one <- data.frame(study_id = "A", icc = 0.90, icc_ci_lo = 0.80, icc_ci_up = 1.20,
                      n_sample = 100, n_measurements = 2)
  s_one <- summary(convert_df(d_one, measure = "icc", icc_to_es = "raw",
                              verbose = FALSE, split_adjusted = FALSE), flags = TRUE)
  se_one <- (1 - 0.9) * sqrt(2 * (1 + 1 * 0.9)^2 / (2 * 1 * (100 - 1)))
  expect_equal(s_one$se[1], se_one, tolerance = 1e-7)   # today: 0.1020427
  expect_true(.has_invalid_naming(s_one$flags[1], "icc_ci_up"))
})


# AUDIT 983-994 (verifier correction at line 994). V38 is documented as always-on,
# but convert_df() rewrites omega_type through .normalise_omega_type() (NA ->
# "unspecified") before validation, and V38 excludes "unspecified". The rescue line
# `ot[is.na(ot) & is.finite(om)] <- "total"` is therefore dead, and a pool of blank
# rows + explicit 'hierarchical' rows is pooled silently. The ICC sibling V43 fires
# on the identical shape because convert_df() preserves icc_type's NA on purpose.
test_that("AUDIT-V38-blank-omega-type: a blank omega_type is the documented default 'total' and mixes with 'hierarchical'", {
  d <- data.frame(study_id = c("A", "B", "C"), omega = c(.86, .72, .88),
                  omega_se = 0.02, omega_type = c(NA, "hierarchical", NA),
                  n_sample = 300, n_items = 12)
  s <- summary(convert_df(d, measure = "omega", verbose = FALSE,
                          split_adjusted = FALSE), flags = TRUE)
  expect_true(all(grepl("mixes omega estimands", s$flags)))

  # control: with the default written out explicitly the flag already fires, so the
  # gap is exactly the blank cell, not the check itself
  d2 <- d
  d2$omega_type <- c("total", "hierarchical", "total")
  s2 <- summary(convert_df(d2, measure = "omega", verbose = FALSE,
                           split_adjusted = FALSE), flags = TRUE)
  expect_true(all(grepl("mixes omega estimands", s2$flags)))

  # and the ICC analogue, which handles blanks correctly, must keep doing so
  d3 <- data.frame(study_id = c("A", "B"), icc = c(.80, .75), icc_se = c(.05, .06),
                   icc_type = c(NA, "consistency"), n_sample = 100, n_measurements = 2)
  s3 <- summary(convert_df(d3, measure = "icc", verbose = FALSE,
                           split_adjusted = FALSE), flags = TRUE)
  expect_true(any(grepl("agreement", s3$flags) & grepl("consistency", s3$flags)))
})


# AUDIT 996-1007 (verifier correction at line 1007). Under omega_se_source =
# "closed_form" the SE is computed from n_sample and n_items, yet the missing-data
# guidance (.method_required_columns()'s omega entry) names only omega_se, so a row
# lacking n_sample is dropped with advice to hunt for a statistic no primary study
# reports. The guidance must be conditional: omega defaults to "reported", where
# naming n_sample would advise a route that cannot work.
test_that("AUDIT-omega-closed-form-guidance: guidance names n_sample when omega_se_source = 'closed_form'", {
  d <- data.frame(study_id = c("A", "B"), omega = c(.86, .90), n_items = c(12, 12))

  s_cf <- summary(convert_df(d, measure = "omega", verbose = FALSE,
                             split_adjusted = FALSE,
                             omega_se_source = "closed_form"), guidance = TRUE)
  expect_true(all(is.na(s_cf$es)))                       # the row IS dropped today
  expect_true(all(grepl("n_sample", s_cf$es_guidance)))  # ... with no mention of why

  # it must stay conditional: on the default "reported" route n_sample cannot help
  s_rep <- summary(convert_df(d, measure = "omega", verbose = FALSE,
                              split_adjusted = FALSE), guidance = TRUE)
  expect_false(any(grepl("n_sample", s_rep$es_guidance)))

  # the advice would be actionable: n_sample alone yields Bonett's (n, k) SE
  se_expected <- sqrt(2 * 12 / ((12 - 1) * (300 - 2)))
  got <- es_from_omega(0.86, n_sample = 300, n_items = 12,
                       omega_to_es = "bonett", omega_se_source = "closed_form")
  expect_equal(got$omega_se, se_expected, tolerance = 1e-9)
})


# AUDIT 1320-1331 (verifier correction at line 1331). ?es_from_omega and the V11
# comment both describe a raw-scale Heywood escape hatch (omega >= 1 preserved so it
# stays inspectable) that the default convert_df() path prevents: V11 NAs the value
# before the route sees it. The verifier rejects the filed fix (skipping the V11
# check would let omega = 1.05 into the pool at full weight) and endorses EITHER
# amending the documentation to name convert_df(correct_inputs = FALSE), OR a
# scale-aware variant that keeps the value on the raw scale under a warn-only flag.
# The assertion is the disjunction, so it goes green under either resolution and is
# red only while the stated behaviour and the actual behaviour disagree.
test_that("AUDIT-omega-heywood: the raw-scale Heywood escape hatch is either reachable or documented as not being reachable", {
  # the route itself does what it says, on a direct call
  direct <- es_from_omega(omega = 1.05, omega_se = 0.02, omega_to_es = "raw")
  expect_equal(direct$omega, 1.05, tolerance = 1e-12)

  d <- data.frame(study_id = c("A", "B"), omega = c(1.05, 0.88),
                  omega_se = c(0.02, 0.02), n_sample = 300, n_items = 10)
  s <- summary(convert_df(d, measure = "omega", omega_to_es = "raw",
                          verbose = FALSE, split_adjusted = FALSE))
  behaviour_ok <- isTRUE(all.equal(s$es[s$study_id == "A"], 1.05, tolerance = 1e-12))

  src <- .audit_src(file.path("R", "es_from_OMEGA.R"))
  skip_if(is.na(src), "package sources not available beside the tests")
  txt <- readLines(src, warn = FALSE)
  start <- grep("Heywood cases", txt, fixed = TRUE)
  skip_if(length(start) == 0, "Heywood paragraph not found in R/es_from_OMEGA.R")
  para <- paste(txt[start[[1]]:min(start[[1]] + 8L, length(txt))], collapse = " ")
  doc_ok <- grepl("correct_inputs", para, fixed = TRUE)

  expect_true(doc_ok || behaviour_ok)
})


# AUDIT 1217-1226. V44's gate is `rel_present & !se_reported`; alpha_se_source and
# omega_se_source are not passed into .validate_input_data() at all. So on the
# DEFAULT omega path (omega_se_source = "reported") the flag describes "the
# closed-form (n, k) reliability standard error" and an inflated I-squared for rows
# that have no SE, no weight and are not in the analysis. The two sources must be
# tracked per coefficient, since rel_present ORs alpha and omega together.
test_that("AUDIT-V44-se-source: V44 fires only where the SE really comes from the closed form", {
  d_om <- data.frame(study_id = c("A", "B"), omega = c(.88, .90), n_items = 10,
                     n_sample = c(200, 220), scale_mean = c(12, 13), scale_min = 1,
                     n_response_categories = 5, omega_type = "total")

  # default omega_se_source = "reported": no closed-form SE exists, no ES, no pool
  s_rep <- summary(convert_df(d_om, measure = "omega", verbose = FALSE,
                              split_adjusted = FALSE), flags = TRUE)
  expect_true(all(is.na(s_rep$es)))
  expect_false(any(grepl("Floor effect", s_rep$flags)))

  # opt in, and the flag is exactly right: it must keep firing
  s_cf <- summary(convert_df(d_om, measure = "omega", verbose = FALSE,
                             split_adjusted = FALSE,
                             omega_se_source = "closed_form"), flags = TRUE)
  expect_true(all(!is.na(s_cf$es)))
  expect_true(all(grepl("Floor effect", s_cf$flags)))

  d_al <- data.frame(study_id = c("A", "B"), cronbach_alpha = c(.88, .90), n_items = 10,
                     n_sample = c(200, 220), scale_mean = c(12, 13), scale_min = 1,
                     n_response_categories = 5)

  # alpha's default IS the closed form, so the flag stands there
  s_al <- summary(convert_df(d_al, measure = "alpha", verbose = FALSE,
                             split_adjusted = FALSE), flags = TRUE)
  expect_true(all(grepl("Floor effect", s_al$flags)))

  # ... and must go quiet when the user refuses it
  s_al_rep <- summary(convert_df(d_al, measure = "alpha", verbose = FALSE,
                                 split_adjusted = FALSE,
                                 alpha_se_source = "reported"), flags = TRUE)
  expect_true(all(is.na(s_al_rep$se)))
  expect_false(any(grepl("Floor effect", s_al_rep$flags)))

  # per coefficient, not ORed: alpha row (closed form, default) fires, omega row
  # (reported, default) does not, in the same sheet
  d_mix <- data.frame(study_id = c("A", "B"), cronbach_alpha = c(.88, NA),
                      omega = c(NA, .90), n_items = 10, n_sample = c(200, 220),
                      scale_mean = c(12, 13), scale_min = 1, n_response_categories = 5)
  s_mix <- summary(convert_df(d_mix, measure = "alpha", verbose = FALSE,
                              split_adjusted = FALSE), flags = TRUE)
  expect_true(grepl("Floor effect", s_mix$flags[s_mix$study_id == "A"]))
  expect_false(grepl("Floor effect", s_mix$flags[s_mix$study_id == "B"]))
})


# AUDIT 1228-1239 (verifier correction at line 1239). V29's upper gate is a FLAT
# 2.5 * xi(n), which is not distribution-free: on error-free skewed / heavy-tailed
# data it fires at 48% (lognormal s=2) and 18% (t3) at n = 500. The verifier rules
# out both "suppress the high side above n = 200" (V29 detects a genuine SD-as-SE
# with probability 1.000 there) and a merely larger flat multiplier (hi_mult = 5
# destroys detection at n < 30), leaving exactly one sound fix: an n-AWARE ceiling.
# The three rows below pin that fix from both sides. The null row is an error-free
# lognormal(0, 2) draw at n = 500 (ratio 17.69, the 75th percentile of that null,
# whose 99.9th percentile is 22.3); a true SD-as-SE inflates the ratio by sqrt(n),
# to 133 at n = 500 and 16.7 at n = 20, so a ceiling exists that separates them.
test_that("AUDIT-V29-n-aware: the range/SD ceiling spares skewed data at large n without losing SD-as-SE", {
  xi <- function(n) 2 * stats::qnorm((n - 0.375) / (n + 0.25))

  d <- data.frame(
    study_id     = c("skewed_n500", "sd_as_se_n500", "sd_as_se_n20"),
    n_exp        = c(500, 500, 20),
    n_nexp       = c(500, 500, 20),
    mean_exp     = c(655, 50, 50),
    mean_nexp    = c(650, 49, 49),
    mean_sd_exp  = c(2438, 0.30, 2.39),
    mean_sd_nexp = c(2438, 0.30, 2.39),
    min_exp      = c(0.4, 30, 30),
    max_exp      = c(43135, 70, 70)
  )
  ratio <- (d$max_exp - d$min_exp) / d$mean_sd_exp
  # the null row is above today's flat ceiling, which is why it false-fires
  expect_gt(ratio[1], 2.5 * xi(500))
  # the two error rows carry the sqrt(n) signature V29 exists to catch
  expect_equal(ratio[2] / ratio[1], 7.5, tolerance = 0.15)
  expect_gt(ratio[3], 2.5 * xi(20))

  s <- summary(convert_df(d, measure = "md", verbose = FALSE,
                          split_adjusted = FALSE), flags = TRUE)
  # Match EITHER V29 message. The two error rows carry the sqrt(n) signature, which
  # pushes the ratio past the arithmetic maximum sqrt(2(n-1)) -- 114 against 31.6 at
  # n = 500, and above 9.3 against 6.16 at n = 20 -- so they are now reported as
  # arithmetically impossible rather than merely implausible. That is a strictly
  # stronger statement about the same rows; asserted separately below.
  fired <- vapply(d$study_id, function(id) {
    grepl("Range/SD ratio", s$flags[s$study_id == id], fixed = TRUE)
  }, logical(1))

  expect_false(unname(fired[["skewed_n500"]]))   # error-free skewed cost outcome
  expect_true(unname(fired[["sd_as_se_n500"]]))  # must not be lost to a wider gate
  expect_true(unname(fired[["sd_as_se_n20"]]))   # must not be lost at small n either

  # ...and both are impossible for ANY distribution, not just improbable under normality
  impossible <- vapply(d$study_id, function(id) {
    grepl("arithmetically impossible", s$flags[s$study_id == id], fixed = TRUE)
  }, logical(1))
  expect_true(unname(impossible[["sd_as_se_n500"]]))
  expect_true(unname(impossible[["sd_as_se_n20"]]))
  expect_false(unname(impossible[["skewed_n500"]]))
})
