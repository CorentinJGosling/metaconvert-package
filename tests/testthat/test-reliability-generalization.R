# Tests for the reliability-generalization (RG) feature set:
#   * alpha_to_es = "hakstian_whalen"   (the transform the RG method literature uses)
#   * C7 / C8 near-perfect alpha/ICC, made scale-aware
#   * V36 reliability induction, V37 item-count constancy
#   * measure = "omega" (McDonald's omega) and V38 mixed-estimand guard

# ---------------------------------------------------------------------------
# Hakstian-Whalen
# ---------------------------------------------------------------------------
test_that("hakstian_whalen matches metafor's AHW bit for bit", {
  skip_if_not_installed("metafor")
  for (g in list(c(0.85, 10, 200), c(0.60, 5, 40), c(0.95, 20, 250))) {
    a <- g[1]; k <- g[2]; n <- g[3]
    pkg <- es_from_cronbach_alpha(a, n, k, alpha_to_es = "hakstian_whalen")
    ora <- metafor::escalc(measure = "AHW", ai = a, mi = k, ni = n)
    expect_equal(pkg$alpha, as.numeric(ora$yi), tolerance = 1e-12)
    expect_equal(pkg$alpha_se^2, as.numeric(ora$vi), tolerance = 1e-12)
  }
})

test_that("hakstian_whalen uses metafor's INCREASING orientation", {
  lo <- es_from_cronbach_alpha(0.70, 200, 10, alpha_to_es = "hakstian_whalen")$alpha
  hi <- es_from_cronbach_alpha(0.95, 200, 10, alpha_to_es = "hakstian_whalen")$alpha
  # Rodriguez & Maeda write T = (1-a)^(1/3), which DECREASES; storing that would
  # make transf.iahw() wrong on our output and flip every moderator coefficient
  expect_lt(lo, hi)
  expect_equal(hi, 1 - (1 - 0.95)^(1 / 3), tolerance = 1e-12)
})

test_that("hakstian_whalen weights differ by alpha where bonett cannot", {
  wr <- function(m) {
    x <- es_from_cronbach_alpha(c(0.85, 0.98), 200, 10, alpha_to_es = m)
    (x$alpha_se[1] / x$alpha_se[2])^2
  }
  # Bonett's SE does not involve alpha at all, so equal (n, k) => equal weight
  expect_equal(wr("bonett"), 1, tolerance = 1e-12)
  expect_gt(wr("hakstian_whalen"), 3)
})

test_that("reliability_backtransform inverts hakstian_whalen without swapping bounds", {
  skip_if_not_installed("metafor")
  a <- c(0.78, 0.85, 0.91)
  e <- es_from_cronbach_alpha(a, 200, 10, alpha_to_es = "hakstian_whalen")
  bt <- reliability_backtransform(e$alpha, e$alpha_ci_lo, e$alpha_ci_up,
                                  method = "hakstian_whalen")
  expect_equal(bt$reliability, a, tolerance = 1e-12)
  expect_equal(bt$reliability, as.numeric(metafor::transf.iahw(e$alpha)), tolerance = 1e-12)
  # increasing map: lower bound stays the lower bound
  expect_true(all(bt$reliability_ci_lo < bt$reliability))
  expect_true(all(bt$reliability_ci_up > bt$reliability))
})

test_that("hakstian_whalen back-transforms the prediction interval too", {
  skip_if_not_installed("metafor")
  e <- es_from_cronbach_alpha(c(.78, .85, .91, .88), 200, 10, alpha_to_es = "hakstian_whalen")
  m <- metafor::rma(yi = e$alpha, sei = e$alpha_se, method = "REML")
  bt <- reliability_backtransform(m, method = "hakstian_whalen")
  # regression: the PI used to pass through untransformed, landing BELOW the estimate
  expect_lt(bt$reliability_pi_lo, bt$reliability)
  expect_gt(bt$reliability_pi_up, bt$reliability)
  expect_true(bt$reliability_pi_lo > 0 && bt$reliability_pi_up < 1)
})

test_that("hakstian_whalen runs through convert_df()", {
  d <- data.frame(study_id = paste0("S", 1:3), cronbach_alpha = c(.82, .88, .91),
                  n_sample = c(200, 300, 250), n_items = 12)
  s <- summary(convert_df(d, measure = "alpha", alpha_to_es = "hakstian_whalen",
                          verbose = FALSE, split_adjusted = FALSE), digits = 15)
  expect_equal(nrow(s), 3L)
  expect_false(anyNA(s$es))
  expect_equal(reliability_backtransform(s$es, method = "hakstian_whalen")$reliability,
               c(.82, .88, .91), tolerance = 1e-10)
})

test_that("an unknown alpha_to_es is still rejected", {
  expect_error(es_from_cronbach_alpha(0.8, 200, 10, alpha_to_es = "fisher_z"),
               "alpha_to_es")
})

# ---------------------------------------------------------------------------
# C7 / C8 scale-awareness
# ---------------------------------------------------------------------------
alpha_flags <- function(a, scale = "bonett") {
  d <- data.frame(study_id = "x", cronbach_alpha = a, n_sample = 300, n_items = 10)
  suppressWarnings(summary(convert_df(d, measure = "alpha", alpha_to_es = scale,
                                      verbose = FALSE, split_adjusted = FALSE),
                           flags = TRUE, digits = 15))$flags
}

test_that("C7 fires for a near-perfect alpha on EVERY scale", {
  # it was guarded by `es > 0 && es <= 1`, unreachable on the always-negative
  # Bonett scale -- i.e. dead under the shipped default
  for (sc in c("bonett", "raw", "hakstian_whalen")) {
    expect_match(alpha_flags(0.995, sc), "Near-perfect alpha", info = sc)
    expect_match(alpha_flags(0.995, sc), "0.995", info = sc)
  }
})

test_that("C7 stays silent for an ordinary alpha", {
  for (sc in c("bonett", "raw", "hakstian_whalen")) {
    expect_false(grepl("Near-perfect", alpha_flags(0.85, sc)), info = sc)
  }
})

test_that("C7 no longer reports a strongly NEGATIVE alpha as near-perfect", {
  # ln(1 - alpha) for alpha in [-1.718, -1.691) lands in (0.99, 1], which passed
  # the old guard and was announced as "Near-perfect alpha: 0.99"
  f <- alpha_flags(-1.70, "bonett")
  expect_false(grepl("Near-perfect", f))
  expect_match(f, "[Nn]egative")
})

test_that("C8 fires for a near-perfect ICC on the Bonett scale", {
  # icc_agreement_se is passed explicitly because this block depends on it: under
  # "drop" an absolute-agreement row has no SE, summary() therefore reports no
  # effect size for it, and every Tier-2 check that reads the ES -- C8 included --
  # goes silent. "compute" is the 2.1.0 default, so this states the dependency
  # rather than overriding anything.
  d <- data.frame(study_id = "x", icc = 0.995, n_sample = 100, n_measurements = 3)
  s <- suppressWarnings(summary(convert_df(d, measure = "icc", verbose = FALSE,
                                           icc_agreement_se = "compute",
                                           split_adjusted = FALSE),
                                flags = TRUE, digits = 15))
  expect_match(s$flags, "Near-perfect ICC")
})

# ---------------------------------------------------------------------------
# V36 reliability induction
# ---------------------------------------------------------------------------
rel_flags <- function(d) {
  suppressWarnings(summary(convert_df(d, measure = "alpha", verbose = FALSE,
                                      split_adjusted = FALSE),
                           flags = TRUE, digits = 15))$flags
}

test_that("V36 flags an identical multi-decimal alpha across DIFFERENT studies", {
  f <- rel_flags(data.frame(study_id = c("Smith2019", "Jones2021", "Lee2020"),
                            cronbach_alpha = c(0.8734, 0.8734, 0.7912),
                            n_sample = c(210, 355, 180), n_items = 20))
  expect_match(f[1], "reliability induction")
  expect_match(f[2], "reliability induction")
  expect_false(grepl("reliability induction", f[3]))
})

test_that("V36 ignores a 2-decimal collision when n differs", {
  # ~30 plausible 2dp values in [.70,.99]: a shared .87 is unremarkable
  f <- rel_flags(data.frame(study_id = c("Smith2019", "Jones2021", "Lee2020"),
                            cronbach_alpha = c(0.87, 0.87, 0.79),
                            n_sample = c(210, 355, 180), n_items = 20))
  expect_false(any(grepl("reliability induction", f)))
})

test_that("V36 fires on a 2-decimal collision when n_sample is identical too", {
  f <- rel_flags(data.frame(study_id = c("Smith2019", "Jones2021", "Lee2020"),
                            cronbach_alpha = c(0.87, 0.87, 0.79),
                            n_sample = c(210, 210, 180), n_items = 20))
  expect_match(f[1], "reliability induction")
  expect_match(f[1], "n_sample is identical")
})

test_that("V36 leaves repeated values WITHIN one study to Category H", {
  f <- rel_flags(data.frame(study_id = c("Smith2019", "Smith2019", "Lee2020"),
                            cronbach_alpha = c(0.8734, 0.8734, 0.7912),
                            n_sample = c(210, 210, 180), n_items = 20))
  expect_false(any(grepl("reliability induction", f)))
  expect_match(f[1], "Duplicate study_id")
})

# ---------------------------------------------------------------------------
# V37 item-count constancy
# ---------------------------------------------------------------------------
test_that("V37 flags the row whose item count differs from a dominant k", {
  f <- rel_flags(data.frame(study_id = paste0("S", 1:5),
                            cronbach_alpha = c(.88, .85, .90, .86, .87),
                            n_sample = c(200, 210, 190, 205, 195),
                            n_items = c(20, 20, 8, 20, 20)))
  expect_match(f[3], "Item count differs")
  expect_match(f[3], "4 of 5")
  expect_false(any(grepl("Item count differs", f[-3])))
})

test_that("V37 stays silent when the pool has no dominant k (multi-instrument)", {
  f <- rel_flags(data.frame(study_id = paste0("S", 1:5),
                            cronbach_alpha = c(.88, .85, .90, .86, .87),
                            n_sample = c(200, 210, 190, 205, 195),
                            n_items = c(10, 12, 20, 8, 33)))
  expect_false(any(grepl("Item count differs", f)))
})

test_that("V37 stays silent when k is constant", {
  f <- rel_flags(data.frame(study_id = paste0("S", 1:4),
                            cronbach_alpha = c(.88, .85, .90, .86),
                            n_sample = c(200, 210, 190, 205), n_items = 20))
  expect_false(any(grepl("Item count differs", f)))
})

# ---------------------------------------------------------------------------
# McDonald's omega
# ---------------------------------------------------------------------------
test_that("es_from_omega derives the SE from a reported natural-scale SE", {
  o <- es_from_omega(omega = 0.86, omega_se = 0.021)
  expect_equal(o$omega, log(1 - 0.86), tolerance = 1e-12)
  # delta method on ln(1 - w): SE_T = SE_w / (1 - w)
  expect_equal(o$omega_se, 0.021 / (1 - 0.86), tolerance = 1e-12)
})

test_that("es_from_omega derives the SE from a reported CI by transforming the BOUNDS", {
  o <- es_from_omega(omega = 0.86, omega_ci_lo = 0.81, omega_ci_up = 0.90)
  expect_equal(o$omega_se,
               (log(1 - 0.81) - log(1 - 0.90)) / (2 * qnorm(0.975)),
               tolerance = 1e-12)
})

test_that("a bare omega keeps its point estimate but gets no SE", {
  o <- es_from_omega(omega = 0.86)
  expect_equal(o$omega, log(1 - 0.86), tolerance = 1e-12)
  expect_true(is.na(o$omega_se))
})

test_that("es_from_omega round-trips on all three scales", {
  for (m in c("bonett", "raw", "hakstian_whalen")) {
    o <- es_from_omega(omega = 0.86, omega_se = 0.021, omega_to_es = m)
    expect_equal(reliability_backtransform(o$omega, method = m)$reliability,
                 0.86, tolerance = 1e-12, info = m)
  }
})

test_that("a Heywood omega >= 1 is NA on transformed scales and kept on raw", {
  expect_true(is.na(es_from_omega(omega = 1.03, omega_se = 0.02)$omega))
  expect_equal(es_from_omega(omega = 1.03, omega_se = 0.02, omega_to_es = "raw")$omega, 1.03)
})

test_that("es_from_omega recycles scalars and validates arguments", {
  o <- es_from_omega(omega = c(.78, .86, .91), omega_se = 0.02)
  expect_equal(nrow(o), 3L)
  expect_false(anyNA(o$omega_se))
  expect_error(es_from_omega(c(.7, .8, .9), omega_se = c(0.02, 0.03)), "omega_se")
  # omega_type comes from a data column, so it warns and falls back per row
  # rather than aborting the run; only the *_to_es API arguments hard-stop
  expect_warning(es_from_omega(0.8, omega_type = "bogus"), "Unrecognised omega_type")
  expect_error(es_from_omega(0.8, omega_to_es = "fisher_z"), "omega_to_es")
})

test_that('measure = "omega" runs end to end and supports meta-regression', {
  skip_if_not_installed("metafor")
  dat <- data.frame(
    study_id = paste0("S", 1:6),
    omega = c(0.86, 0.91, 0.78, 0.88, 0.83, 0.90),
    omega_se = c(0.021, NA, 0.035, NA, 0.028, NA),
    omega_ci_lo = c(NA, 0.87, NA, 0.84, NA, 0.86),
    omega_ci_up = c(NA, 0.94, NA, 0.92, NA, 0.94),
    omega_type = "total", n_sample = c(300, 450, 180, 260, 320, 400),
    n_items = 12, language = c("en", "fr", "en", "fr", "en", "fr"))
  s <- summary(convert_df(dat, measure = "omega", verbose = FALSE, split_adjusted = FALSE),
               digits = 15)
  expect_equal(nrow(s), 6L)
  expect_false(anyNA(s$es))
  expect_false(anyNA(s$se))
  expect_true(all(s$info_used == "omega"))
  # arbitrary moderator columns survive, so meta-regression is possible
  expect_true("language" %in% names(s))
  mm <- metafor::rma(yi = s$es, sei = s$se, mods = ~ language, data = s, method = "REML")
  expect_equal(length(mm$beta), 2L)
  bt <- reliability_backtransform(metafor::predict.rma(mm))
  expect_equal(nrow(bt), 6L)
  expect_true(all(bt$reliability > 0 & bt$reliability < 1))
})

test_that("the omega extraction sheet carries the omega columns", {
  sh <- data_extraction_sheet(measure = "omega", extension = "data.frame")
  expect_true(all(c("omega", "omega_se", "omega_ci_lo", "omega_ci_up", "omega_type")
                  %in% names(sh)))
})

test_that("guidance names what a bare omega row is missing", {
  d <- data.frame(study_id = c("A", "B"), omega = c(0.86, 0.90),
                  omega_se = c(0.02, NA), n_sample = 300, n_items = 12)
  s <- summary(convert_df(d, measure = "omega", verbose = FALSE, split_adjusted = FALSE),
               digits = 15, guidance = TRUE)
  expect_match(s$es_guidance[2], "omega_se")
})

test_that("V38 flags a pool that mixes omega estimands", {
  d <- data.frame(study_id = c("A", "B", "C"), omega = c(.86, .72, .88),
                  omega_se = 0.02, omega_type = c("total", "hierarchical", "total"),
                  n_sample = 300, n_items = 12)
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE),
                                flags = TRUE, digits = 15))
  expect_true(all(grepl("mixes omega estimands", s$flags)))
  # severity matches its stated analogues E6/E8: nothing is mis-extracted, so this
  # is a disclosure about pooling, not a suspicion about the data
  expect_true(all(grepl("[INFO] Pool mixes omega", s$flags, fixed = TRUE)))
  expect_false(any(grepl("[UNUSUAL] Pool mixes omega", s$flags, fixed = TRUE)))
  d2 <- d; d2$omega_type <- "total"
  s2 <- suppressWarnings(summary(convert_df(d2, measure = "omega", verbose = FALSE,
                                            split_adjusted = FALSE),
                                 flags = TRUE, digits = 15))
  expect_false(any(grepl("mixes omega estimands", s2$flags)))
})

test_that("a raw omega cannot pass through user_es_* into a transformed pool", {
  d <- data.frame(study_id = c("native", "via_user"),
                  omega = c(0.88, NA), omega_se = c(0.02, NA), n_sample = 300, n_items = 12,
                  user_es_crude = c(NA, 0.91), user_se_crude = c(NA, 0.02),
                  user_es_original_measure_crude = c(NA, "omega"),
                  user_es_target_measure_crude = c(NA, "omega"))
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE), digits = 15))
  expect_true(is.na(s$es[s$study_id == "via_user"]))
  expect_equal(s$es[s$study_id == "native"], log(1 - 0.88), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Defects found by reading the supplied omega literature (Beland, Flora,
# Malkewitz, Revelle & Zinbarg, Zinbarg et al.) against the implementation.
# ---------------------------------------------------------------------------
test_that("a raw-scale CI whose upper bound is exactly 1.00 still yields an SE", {
  # "omega = .88 [.79, 1.00]" is a routine bootstrap print-out at high reliability.
  # The ci < 1 guard is needed only where log(0) / (0)^(1/3) would be taken.
  o <- es_from_omega(omega = 0.88, omega_ci_lo = 0.79, omega_ci_up = 1.00,
                     omega_to_es = "raw")
  expect_true(is.finite(o$omega_se))
  expect_equal(o$omega_se, (1.00 - 0.79) / (2 * qnorm(0.975)), tolerance = 1e-12)
  # on the transformed scales the bound genuinely has no image, so NA is correct
  expect_true(is.na(es_from_omega(omega = 0.88, omega_ci_lo = 0.79,
                                  omega_ci_up = 1.00)$omega_se))
})

test_that("a mis-cased or synonymous omega_type does not abort convert_df()", {
  d <- data.frame(study_id = c("a", "b", "c"), omega = c(.86, .88, .90),
                  omega_se = 0.02, omega_type = c("total", "Total", "omega_t"),
                  n_sample = 300, n_items = 12)
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE), digits = 15))
  expect_equal(nrow(s), 3L)
  expect_false(anyNA(s$es))
  # all three normalise to the same estimand, so V38 must NOT fire
  expect_false(any(grepl("mixes omega estimands", s$flags)))
  # an unrecognised value warns and falls back rather than stopping
  expect_warning(es_from_omega(0.8, omega_se = 0.02, omega_type = "bogus"),
                 "Unrecognised omega_type")
  # the fallback is deliberately "unspecified", NOT the default estimand: mapping an
  # unreadable value onto "total" would silently relabel a hierarchical omega and
  # disarm V38, the flag that exists for exactly that error
  expect_equal(suppressWarnings(
    es_from_omega(0.8, omega_se = 0.02, omega_type = "bogus")$omega_type), "unspecified")
})

test_that("a mis-cased icc_type does not abort convert_df() either", {
  expect_equal(suppressWarnings(
    es_from_icc(0.8, 50, 3, icc_type = "Agreement"))$icc_type, "agreement")
  expect_equal(suppressWarnings(
    es_from_icc(0.8, 50, 3, icc_type = "Consistency"))$icc_type, "consistency")
  expect_warning(es_from_icc(0.8, 50, 3, icc_type = "bogus"), "Unrecognised icc_type")
})

test_that("a negative omega is impossible and is now caught", {
  # omega = (sum lambda)^2 / ((sum lambda)^2 + sum theta) is structurally >= 0,
  # unlike alpha. It used to pass every check: es = +0.113, se = 0.045, no flag.
  d <- data.frame(study_id = "x", omega = -0.12, omega_se = 0.05,
                  n_sample = 300, n_items = 12)
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE),
                                flags = TRUE, digits = 15))
  expect_true(is.na(s$es))
  expect_match(s$flags, "INVALID")
  expect_match(s$flags, "omega")
  # a valid omega is untouched
  d2 <- d; d2$omega <- 0.86
  s2 <- suppressWarnings(summary(convert_df(d2, measure = "omega", verbose = FALSE,
                                            split_adjusted = FALSE), digits = 15))
  expect_equal(s2$es, log(1 - 0.86), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# omega_estimator provenance + V39 mixed-estimator flag
# ---------------------------------------------------------------------------
omega_est_flags <- function(est) {
  d <- data.frame(study_id = paste0("S", seq_along(est)),
                  omega = 0.85 + seq_along(est) / 100, omega_se = 0.02,
                  omega_estimator = est, n_sample = 300, n_items = 12)
  suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                      split_adjusted = FALSE),
                           flags = TRUE, digits = 15))$flags
}

test_that("V39 fires when two KNOWN omega estimators are pooled", {
  f <- omega_est_flags(c("cfa_bifactor", "first_pc"))
  expect_true(all(grepl("mixes omega estimators", f)))
  # severity matches V38 / E6 / E8: nothing is mis-extracted
  expect_true(all(grepl("[INFO] Pool mixes omega estimators", f, fixed = TRUE)))
})

test_that("V39 ignores 'unspecified' -- the design point", {
  # most primary studies do not name their estimator; counting "unspecified" as a
  # level would make this fire on nearly every real dataset
  expect_false(any(grepl("mixes omega estimators",
                         omega_est_flags(c("cfa_bifactor", "unspecified", "unspecified")))))
  expect_false(any(grepl("mixes omega estimators",
                         omega_est_flags(c("unspecified", "unspecified")))))
})

test_that("V39 stays silent on one estimator, however it is spelled", {
  expect_false(any(grepl("mixes omega estimators",
                         omega_est_flags(c("cfa_bifactor", "cfa_bifactor")))))
  expect_false(any(grepl("mixes omega estimators",
                         omega_est_flags(c("CFA_Bifactor", "bifactor", "Hi-CF")))))
})

test_that("omega_estimator normalises the synonyms a reviewer would actually type", {
  n <- function(x) suppressWarnings(.normalise_omega_estimator(x))
  expect_equal(n("bifactor"), "cfa_bifactor")
  expect_equal(n("Schmid-Leiman"), "efa_schmid_leiman")
  expect_equal(n("PCA"), "first_pc")
  expect_equal(n("first PF"), "first_pf")
  expect_equal(n("CFA"), "cfa_1factor")
  expect_equal(n("not reported"), "unspecified")
  # psych::omega() defaults to EFA + Schmid-Leiman, so the package name determines it
  expect_equal(n("psych_omega"), "efa_schmid_leiman")
  expect_equal(n("psych"), "efa_schmid_leiman")
  # but semTools / MBESS are CFA-based WITHOUT fixing the model, so they must NOT be
  # guessed into a level the mixing flag would then act on
  expect_equal(n("semTools"), "unspecified")
  expect_equal(n("MBESS"), "unspecified")
  # unrecognised values warn and fall back rather than aborting the run
  expect_warning(es_from_omega(0.8, omega_se = 0.02, omega_estimator = "bogus"),
                 "Unrecognised omega_estimator")
})

test_that("omega_estimator reaches the extraction sheet and survives convert_df()", {
  sh <- data_extraction_sheet(measure = "omega", extension = "data.frame")
  expect_true("omega_estimator" %in% names(sh))
  d <- data.frame(study_id = c("a", "b"), omega = c(.86, .88), omega_se = 0.02,
                  omega_estimator = c("bifactor", "PCA"), n_sample = 300, n_items = 12)
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE), digits = 15))
  expect_equal(nrow(s), 2L)
  expect_false(anyNA(s$es))
})

# ---------------------------------------------------------------------------
# Defects found by a hostile end-to-end dry run of a real RG analysis.
# Every one of these produced a wrong weight or a silent study loss.
# ---------------------------------------------------------------------------
test_that("omega confidence intervals are validated (its SE comes ENTIRELY from them)", {
  d <- data.frame(study_id = c("ok", "transposed", "outside"), omega = c(.86, .86, .95),
                  omega_ci_lo = c(.81, .90, .70), omega_ci_up = c(.90, .81, .80),
                  n_sample = 300, n_items = 12)
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE), flags = TRUE, digits = 15))
  expect_false(is.na(s$es[1]))
  expect_true(is.na(s$es[2])); expect_match(s$flags[2], "Inverted CI for 'omega'")
  expect_true(is.na(s$es[3])); expect_match(s$flags[3], "outside CI for 'omega'")
})

test_that("a non-positive omega_se is flagged rather than silently dropping the study", {
  d <- data.frame(study_id = c("neg", "zero", "ok"), omega = c(.86, .80, .90),
                  omega_se = c(-0.03, 0, 0.02), n_sample = 300, n_items = 12)
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE), flags = TRUE, digits = 15))
  expect_match(s$flags[1], "Negative input: 'omega_se'")
  expect_match(s$flags[2], "Zero SE: 'omega_se'")
  expect_false(is.na(s$es[3]))
})

test_that("an unrecognised omega_type falls back to 'unspecified', not to an estimand", {
  # mapping it onto "total" silently relabels a hierarchical omega AND disarms V38
  n <- function(x) suppressWarnings(.normalise_omega_type(x))
  expect_equal(n("omega_hierarchical"), "hierarchical")
  expect_equal(n("general factor"), "hierarchical")
  expect_equal(n("Omega Hierarchical"), "hierarchical")
  expect_equal(n("bogus"), "unspecified")
  # and V38 ignores unspecified, as V39 does
  d <- data.frame(study_id = c("a", "b"), omega = c(.86, .72), omega_se = 0.02,
                  omega_type = c("total", "bogus"), n_sample = 300, n_items = 12)
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE), flags = TRUE, digits = 15))
  expect_false(any(grepl("mixes omega estimands", s$flags)))
})

test_that("a coefficient of exactly 1 is explained instead of vanishing (V40)", {
  d <- data.frame(study_id = c("a", "b"), cronbach_alpha = c(1.00, .87),
                  n_sample = 300, n_items = 18)
  s <- suppressWarnings(summary(convert_df(d, measure = "alpha", verbose = FALSE,
                                           split_adjusted = FALSE), flags = TRUE, digits = 15))
  expect_true(is.na(s$es[1]))
  expect_match(s$flags[1], "= 1 exactly")
  expect_match(s$flags[1], "no bonett transform")
  # on the raw scale 1 is a usable boundary, so V40 must stay silent there
  s_raw <- suppressWarnings(summary(convert_df(d, measure = "alpha", alpha_to_es = "raw",
                                               verbose = FALSE, split_adjusted = FALSE),
                                    flags = TRUE, digits = 15))
  expect_equal(s_raw$es[1], 1)
  expect_false(grepl("= 1 exactly", s_raw$flags[1]))
})

test_that("summary() returns NORMALISED omega_estimator so meta-regression works", {
  skip_if_not_installed("metafor")
  d <- data.frame(study_id = paste0("S", 1:6), omega = c(.86, .88, .90, .84, .91, .87),
                  omega_se = 0.02,
                  omega_estimator = c("psych", "bifactor CFA", "PCA",
                                      "psych_omega", "cfa_bifactor", "first PC"),
                  n_sample = 300, n_items = 12)
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", verbose = FALSE,
                                           split_adjusted = FALSE), digits = 15))
  expect_equal(sort(unique(s$omega_estimator)),
               c("cfa_bifactor", "efa_schmid_leiman", "first_pc"))
  m <- metafor::rma(yi = s$es, sei = s$se, mods = ~ s$omega_estimator, method = "REML")
  expect_equal(length(m$beta), 3L)
})

test_that("D1 does not flag a homogeneous reliability pool, but still catches a real outlier", {
  set.seed(1)
  a <- round(runif(20, .91, .93), 3)
  d <- data.frame(study_id = paste0("S", 1:20), cronbach_alpha = a,
                  n_sample = sample(150:400, 20), n_items = 18)
  for (sc in c("bonett", "hakstian_whalen", "raw")) {
    s <- suppressWarnings(summary(convert_df(d, measure = "alpha", alpha_to_es = sc,
                                             verbose = FALSE, split_adjusted = FALSE),
                                  flags = TRUE, digits = 15))
    expect_equal(sum(grepl("ES outlier", s$flags)), 0L, info = sc)
  }
  d2 <- d; d2$cronbach_alpha[7] <- 0.30
  s2 <- suppressWarnings(summary(convert_df(d2, measure = "alpha", verbose = FALSE,
                                            split_adjusted = FALSE), flags = TRUE, digits = 15))
  expect_true(grepl("ES outlier", s2$flags[7]))
})

test_that("V36 and V37 route to the crude scope only", {
  d <- data.frame(study_id = paste0("S", 1:5), cronbach_alpha = c(.88, .85, .90, .86, .87),
                  n_sample = c(200, 210, 190, 205, 195), n_items = c(20, 20, 8, 20, 20))
  s <- suppressWarnings(summary(convert_df(d, measure = "alpha", verbose = FALSE),
                                flags = TRUE, digits = 15))
  expect_gt(sum(nzchar(s$flags_crude)), 0)
  # measure = "alpha" has no adjusted hierarchy, so nothing may land there
  expect_equal(sum(nzchar(s$flags_adjusted)), 0L)
})

test_that("the omega passthrough refusal names omega_to_es, not prop_to_es", {
  d <- data.frame(study_id = c("n", "u"), omega = c(0.88, NA), omega_se = c(0.02, NA),
                  n_sample = 300, n_items = 12,
                  user_es_crude = c(NA, 0.91), user_se_crude = c(NA, 0.02),
                  user_es_original_measure_crude = c(NA, "omega"),
                  user_es_target_measure_crude = c(NA, "omega"))
  w <- character(0)
  withCallingHandlers({convert_df(d, measure = "omega", verbose = FALSE,
                                  split_adjusted = FALSE); NULL},
                      warning = function(x) { w <<- c(w, conditionMessage(x))
                                              invokeRestart("muffleWarning") })
  expect_true(any(grepl("omega_to_es", w)))
  expect_false(any(grepl("prop_to_es", w)))
  # and omega_to_es = "raw" really does unlock it, as the message says
  s <- suppressWarnings(summary(convert_df(d, measure = "omega", omega_to_es = "raw",
                                           verbose = FALSE, split_adjusted = FALSE), digits = 15))
  expect_equal(s$es[2], 0.91)
})

test_that("V31 is scoped to measure = 'icc' and does not leak into other runs", {
  # A COSMIN-style sheet holds alpha, omega and ICC columns side by side and is
  # run once per property. Tier-1 checks key on which COLUMNS exist, not on the
  # requested measure, so without an explicit test the ICC standard-error note
  # appeared on every alpha and omega run off the same sheet.
  d <- data.frame(study_id = c("a", "b"),
                  n_sample = c(300, 120),
                  cronbach_alpha = c(0.87, NA), n_items = c(9, NA),
                  omega = c(0.89, NA), omega_se = c(0.02, NA),
                  icc = c(NA, 0.94), n_measurements = c(NA, 2),
                  icc_type = c(NA, "agreement"))
  f <- function(m) suppressWarnings(summary(
    convert_df(d, measure = m, split_adjusted = FALSE, verbose = FALSE),
    flags = TRUE, digits = 15))$flags
  expect_false(any(grepl("ICC agreement-type SE", f("alpha"))))
  expect_false(any(grepl("ICC agreement-type SE", f("omega"))))
  expect_true(any(grepl("ICC agreement-type SE", f("icc"))))
})


# ---------------------------------------------------------------------------
# Reliability roadmap 1.3: one definition of what icc_type means.
#
# es_from_icc() normalised its icc_type (case-fold + alias map + fallback) while
# the V31 check used bare string equality against the literal 'agreement'. So of
# the ten spellings that resolve to agreement, only the two already spelled
# exactly 'agreement' fired V31 -- i.e. the rows that GET the anti-conservative
# agreement SE were largely the rows that were not warned about it.
# ---------------------------------------------------------------------------

agreement_spellings <- c('agreement', 'Agreement', 'AGREEMENT', '  agreement  ',
                         'absolute agreement', 'Absolute Agreement',
                         'ICC(2,1)', 'icc21', 'two-way random', 'Two-Way Random')
consistency_spellings <- c('consistency', 'Consistency', 'ICC(3,1)', 'icc31',
                           'two-way mixed')

test_that('.normalise_icc_type resolves every documented spelling', {
  expect_equal(.normalise_icc_type(agreement_spellings),
               rep('agreement', length(agreement_spellings)))
  expect_equal(.normalise_icc_type(consistency_spellings),
               rep('consistency', length(consistency_spellings)))
})

test_that('.normalise_icc_type keeps digits, so ICC(3,1) is NOT relabelled agreement', {
  # The old inline map stripped non-letters, which made its own 'icc21'/'icc31'
  # entries unreachable: every digit-bearing spelling collapsed to the key 'icc',
  # missed the map and fell through to the agreement default -- so a CONSISTENCY
  # ICC entered as 'ICC(3,1)' was silently relabelled, with a warning calling it
  # unrecognised. This is the assertion that fails against the pre-fix tree.
  expect_equal(.normalise_icc_type('ICC(3,1)'), 'consistency')
  expect_equal(.normalise_icc_type('icc31'), 'consistency')
  expect_silent(.normalise_icc_type('ICC(3,1)'))
})

test_that('.normalise_icc_type is idempotent and falls back with a warning', {
  # convert_df() normalises the column and es_from_icc() normalises again, so
  # every output must also be a valid input.
  once <- .normalise_icc_type(c(agreement_spellings, consistency_spellings))
  expect_equal(.normalise_icc_type(once), once)
  expect_warning(out <- .normalise_icc_type('banana'), 'Unrecognised icc_type')
  expect_equal(out, 'agreement')
  expect_silent(.normalise_icc_type('banana', warn = FALSE))
})

test_that('V31 and es_from_icc() agree on every spelling, end to end', {
  fires_v31 <- function(s) {
    d <- data.frame(icc = 0.8, n_sample = 50, n_measurements = 2, icc_type = s)
    f <- suppressWarnings(summary(convert_df(d, measure = 'icc', verbose = FALSE),
                                  flags = TRUE))$flags_crude
    any(grepl('ICC agreement-type SE', f, fixed = TRUE))
  }
  route <- function(s) {
    suppressWarnings(es_from_icc(icc = 0.8, n_sample = 50, n_measurements = 2,
                                 icc_type = s))$icc_type[1]
  }
  for (s in c(agreement_spellings, consistency_spellings, 'banana')) {
    expect_equal(fires_v31(s), route(s) == 'agreement', info = paste('icc_type =', s))
  }
})

test_that('convert_df stores the normalised icc_type but leaves NA as NA', {
  d <- data.frame(icc = c(0.8, 0.7, 0.6),
                  n_sample = 50, n_measurements = 2,
                  icc_type = c('ICC(3,1)', 'Absolute Agreement', NA))
  raw <- attr(suppressWarnings(convert_df(d, measure = 'icc', verbose = FALSE)), 'raw_data')
  expect_equal(raw$icc_type, c('consistency', 'agreement', NA))
})


# ---------------------------------------------------------------------------
# Reliability roadmap 1.2: an average-measures ICC is a different estimand.
#
# 'average' / 'ICC(2,k)' / 'icc2k' used to be unrecognised: they warned, fell back
# to 'agreement', and were computed as single-measures. The SE ratio stays near 1
# (1.05-1.39x) so nothing downstream looks wrong, while the point estimate is off
# by up to 1.9 log units. Now stepped down with Spearman-Brown where k is known,
# and NA'd with [INVALID] where it is not.
# ---------------------------------------------------------------------------

test_that('.normalise_icc_type recognises the average-measures forms as their own level', {
  expect_equal(.normalise_icc_type(c('average', 'Average', 'ICC(2,k)', 'icc2k',
                                     'average measures', 'avg')),
               rep('agreement_average', 6))
  expect_equal(.normalise_icc_type(c('ICC(3,k)', 'icc3k', 'two-way mixed average')),
               rep('consistency_average', 3))
  # still idempotent with the new levels
  once <- .normalise_icc_type(c('ICC(2,k)', 'icc3k', 'ICC(2,1)', 'consistency'))
  expect_equal(.normalise_icc_type(once), once)
  # and they no longer warn
  expect_silent(.normalise_icc_type('ICC(2,k)'))
})

test_that('the Spearman-Brown step-down reproduces the roadmap 1.2 cost table', {
  # k | ICC_avg | true single-measure ICC_1 | correct es = log(1 - ICC_1)
  tab <- list(c(2, 0.90, 0.8181818), c(5, 0.90, 0.6428571), c(10, 0.95, 0.6551724))
  for (r in tab) {
    k <- r[1]; a <- r[2]; rho1 <- r[3]
    expect_equal(.icc_step_down(a, k), rho1, tolerance = 1e-6)
    got <- es_from_icc(icc = a, n_sample = 50, n_measurements = k, icc_type = 'average')
    expect_equal(got$icc, log(1 - rho1), tolerance = 1e-6,
                 info = paste('k =', k, 'ICC_avg =', a))
  }
  # the un-stepped value is what the pre-fix code returned; assert we moved off it
  expect_false(isTRUE(all.equal(
    es_from_icc(icc = 0.95, n_sample = 50, n_measurements = 10, icc_type = 'average')$icc,
    es_from_icc(icc = 0.95, n_sample = 50, n_measurements = 10, icc_type = 'agreement')$icc)))
})

test_that('the step-down inverts Spearman-Brown exactly', {
  # round trip: build an average-of-k ICC from a known single-measure one
  for (k in c(2, 3, 5, 10)) {
    for (rho1 in c(-0.2, 0.05, 0.4, 0.8, 0.99)) {
      rho_k <- k * rho1 / (1 + (k - 1) * rho1)
      expect_equal(.icc_step_down(rho_k, k), rho1, tolerance = 1e-10,
                   info = paste('k =', k, 'rho1 =', rho1))
    }
  }
})

test_that('a stepped-down row reports the single-measures estimand it computed', {
  for (t in c('average', 'ICC(2,k)', 'icc2k')) {
    expect_equal(es_from_icc(0.9, 50, 5, icc_type = t)$icc_type, 'agreement')
  }
  for (t in c('ICC(3,k)', 'icc3k')) {
    expect_equal(es_from_icc(0.9, 50, 5, icc_type = t)$icc_type, 'consistency')
  }
})

test_that('without n_measurements an average-measures ICC is NA, not single-measures', {
  got <- es_from_icc(icc = 0.9, n_sample = 50, n_measurements = NA, icc_type = 'average')
  expect_true(is.na(got$icc))
  expect_true(is.na(got$icc_se))
})

test_that('convert_df flags the step-down and the drop, and V31 skips the dropped row', {
  d <- data.frame(icc = c(0.90, 0.95, 0.90, 0.80), n_sample = 50,
                  n_measurements = c(5, 10, NA, 2),
                  icc_type = c('ICC(2,k)', 'icc2k', 'average', 'agreement'))
  # icc_agreement_se stated explicitly: this block tests the STEP-DOWN arithmetic,
  # which needs the row to survive into summary(). Under 'drop' an agreement row
  # carries no SE and no effect size is reported. 'compute' is the 2.1.0 default.
  s <- summary(convert_df(d, measure = 'icc', icc_agreement_se = 'compute',
                          verbose = FALSE), flags = TRUE)
  f <- s$flags_crude

  expect_equal(s$es_crude[1], log(1 - .icc_step_down(0.90, 5)), tolerance = 1e-4)
  expect_equal(s$es_crude[2], log(1 - .icc_step_down(0.95, 10)), tolerance = 1e-4)
  expect_true(is.na(s$es_crude[3]))
  expect_equal(s$es_crude[4], log(1 - 0.80), tolerance = 1e-4)

  expect_true(grepl('stepped down to the single-measurement ICC 0.6429', f[1], fixed = TRUE))
  expect_true(grepl('stepped down to the single-measurement ICC 0.6552', f[2], fixed = TRUE))
  expect_true(grepl("[INVALID] Average-measures ICC needs 'n_measurements'", f[3], fixed = TRUE))
  # the dropped row must NOT also carry the [INFO] step-down report
  expect_false(grepl('[INFO] Average-measures ICC (', f[3], fixed = TRUE))
  # Every '; '-separated token must be a TAGGED flag. '; ' is the merge separator,
  # so a message that contains it internally is split into an untagged fragment --
  # and a fragment carrying no quoted column name matches every scope, landing in
  # BOTH crude and adjusted. That is the V23 defect (audit #41) in general form.
  for (j in seq_along(f)) {
    if (!nzchar(f[j])) next
    toks <- strsplit(f[j], '; ', fixed = TRUE)[[1]]
    expect_true(all(startsWith(toks, '[')),
                info = paste('row', j, '- untagged token:',
                             paste(toks[!startsWith(toks, '[')], collapse = ' | ')))
  }

  # V31 follows the step-down onto the agreement rows...
  expect_true(grepl('ICC agreement-type SE', f[1], fixed = TRUE))
  expect_true(grepl('ICC agreement-type SE', f[2], fixed = TRUE))
  # ...but not onto the row that produced no estimate at all
  expect_false(grepl('ICC agreement-type SE', f[3], fixed = TRUE))
})

test_that('a consistency average-measures row is not given the agreement SE note', {
  d <- data.frame(icc = 0.9, n_sample = 50, n_measurements = 5, icc_type = 'ICC(3,k)')
  f <- summary(convert_df(d, measure = 'icc', verbose = FALSE), flags = TRUE)$flags_crude
  expect_true(grepl('stepped down', f[1], fixed = TRUE))
  expect_false(grepl('agreement-type SE', f[1], fixed = TRUE))
})


# ---------------------------------------------------------------------------
# Reliability roadmap 1.4: a reported ICC standard error or interval is usable.
#
# icc_se / icc_ci_lo / icc_ci_up were SILENTLY IGNORED. The ICC literature reports
# intervals routinely, and this matters twice over: it is the only uncertainty a
# study without a sample size supplies at all, and for an agreement-type ICC the
# package's own (n, k) SE is the one whose coverage runs 0.82 at n = 20 down to
# 0.14 at n = 1000, so a reported interval is the better source when both exist.
# Precedence mirrors es_from_omega(): reported SE > reported CI > computed (n, k).
# ---------------------------------------------------------------------------

test_that('the SE precedence is reported SE, then reported CI, then computed (n, k)', {
  nk <- es_from_icc(0.80, 50, 2)
  se <- es_from_icc(0.80, 50, 2, icc_se = 0.05)
  ci <- es_from_icc(0.80, 50, 2, icc_ci_lo = 0.70, icc_ci_up = 0.87)
  both <- es_from_icc(0.80, 50, 2, icc_se = 0.05, icc_ci_lo = 0.70, icc_ci_up = 0.87)

  # the point estimate is the same throughout; only the variance source changes
  expect_equal(nk$icc, log(1 - 0.80)); expect_equal(se$icc, log(1 - 0.80))
  # a raw-scale reported SE is delta-mapped onto the Bonett scale: se / (1 - rho)
  expect_equal(se$icc_se, 0.05 / (1 - 0.80), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(se$icc_se, nk$icc_se)))   # it really overrode (n, k)
  expect_false(isTRUE(all.equal(ci$icc_se, nk$icc_se)))
  expect_equal(both$icc_se, se$icc_se, tolerance = 1e-12) # SE beats CI
})

test_that('an ICC reported with a CI and NO sample size is usable', {
  # the motivating quote from the corpus: 'ICC 0.94 (95% CI 0.86-0.98)'
  got <- es_from_icc(icc = 0.94, icc_ci_lo = 0.86, icc_ci_up = 0.98)
  expect_equal(got$icc, log(1 - 0.94), tolerance = 1e-10)
  expect_false(is.na(got$icc_se))
  # pre-fix this row was dropped whole: n_sample/n_measurements gated the ESTIMATE
  expect_true(is.finite(got$icc_se))
})

test_that('the CI is transformed at the BOUNDS, not symmetrised first', {
  got <- es_from_icc(0.94, icc_ci_lo = 0.86, icc_ci_up = 0.98)
  expect_equal(got$icc_se,
               abs(log(1 - 0.98) - log(1 - 0.86)) / (2 * qnorm(0.975)),
               tolerance = 1e-12)
  # a transposed interval is the same interval (.ci_lower/.ci_upper guards)
  expect_equal(es_from_icc(0.94, icc_ci_lo = 0.98, icc_ci_up = 0.86)$icc_se,
               got$icc_se, tolerance = 1e-12)
})

test_that('a bare ICC keeps its estimate and gets se = NA (the omega precedent)', {
  got <- es_from_icc(icc = 0.80)
  expect_equal(got$icc, log(1 - 0.80), tolerance = 1e-10)
  expect_true(is.na(got$icc_se))
})

test_that('a reported CI travels with the Spearman-Brown step-down', {
  # otherwise the row ends up with an ICC(1) point estimate and an ICC(k) interval
  sd1 <- function(x, k) x / (k - (k - 1) * x)
  got <- es_from_icc(icc = 0.90, n_measurements = 5, icc_type = 'ICC(2,k)',
                     icc_ci_lo = 0.85, icc_ci_up = 0.94)
  expect_equal(got$icc, log(1 - sd1(0.90, 5)), tolerance = 1e-10)
  expect_equal(got$icc_se,
               abs(log(1 - sd1(0.94, 5)) - log(1 - sd1(0.85, 5))) / (2 * qnorm(0.975)),
               tolerance = 1e-10)
})

test_that('a reported SE travels with the step-down by the SB derivative', {
  # d(rho_1)/d(rho_k) = k / (k - (k-1) rho_k)^2
  sd1 <- function(x, k) x / (k - (k - 1) * x)
  k <- 5; rho_k <- 0.90; raw_se <- 0.02
  dfac <- k / (k - (k - 1) * rho_k)^2
  got <- es_from_icc(icc = rho_k, n_measurements = k, icc_type = 'icc2k',
                     icc_se = raw_se)
  expect_equal(got$icc_se, raw_se * dfac / (1 - sd1(rho_k, k)), tolerance = 1e-10)
})

test_that('on the raw scale a reported SE passes through unchanged', {
  got <- es_from_icc(0.80, 50, 2, icc_se = 0.05, icc_to_es = 'raw')
  expect_equal(got$icc, 0.80)
  expect_equal(got$icc_se, 0.05, tolerance = 1e-12)
})

test_that('a non-positive reported SE and an out-of-range CI bound degrade safely', {
  # .positive_or_na(): a negative SE is not a standard error
  expect_true(is.na(es_from_icc(0.80, icc_se = -0.05)$icc_se))
  expect_true(is.na(es_from_icc(0.80, icc_se = 0)$icc_se))
  # a bound of exactly 1 has no log(1 - .) on the Bonett scale -> fall through
  expect_true(is.na(es_from_icc(0.90, icc_ci_lo = 0.80, icc_ci_up = 1.00)$icc_se))
  # ...but it is fine on the raw scale, where fwd() is the identity
  expect_false(is.na(es_from_icc(0.90, icc_ci_lo = 0.80, icc_ci_up = 1.00,
                                 icc_to_es = 'raw')$icc_se))
})

test_that('the new arguments recycle length-1 values and reject a bad length', {
  v <- es_from_icc(icc = c(0.7, 0.8, 0.9), n_sample = 50, n_measurements = 2,
                   icc_se = 0.05)
  expect_equal(length(v$icc_se), 3L)
  expect_true(all(!is.na(v$icc_se)))
  expect_error(es_from_icc(c(.7, .8, .9), 50, 2, icc_se = c(0.05, 0.06)), 'icc_se')
})

test_that('convert_df threads icc_se / icc_ci through the pipeline', {
  # row 3 is the closed-form (n, k) route. It is marked 'consistency' because the
  # agreement default withholds exactly that SE (roadmap 1.1) -- which would test
  # the drop, not the threading.
  d <- data.frame(icc = c(0.80, 0.94, 0.80),
                  n_sample = c(50, NA, 50), n_measurements = c(2, NA, 2),
                  icc_type = c('agreement', 'agreement', 'consistency'),
                  icc_se = c(0.05, NA, NA),
                  icc_ci_lo = c(NA, 0.86, NA), icc_ci_up = c(NA, 0.98, NA))
  s <- summary(convert_df(d, measure = 'icc', verbose = FALSE))
  expect_equal(s$se_crude[1], 0.05 / (1 - 0.80), tolerance = 1e-4)   # reported SE
  expect_equal(s$se_crude[2],
               abs(log(1 - 0.98) - log(1 - 0.86)) / (2 * qnorm(0.975)),
               tolerance = 1e-4)                                     # reported CI
  expect_equal(s$se_crude[3], sqrt(2 * (1 + 0.80)^2 / (2 * 1 * 49)), tolerance = 1e-4)
  # the columns survive .check_data() rather than being dropped
  raw <- attr(convert_df(d, measure = 'icc', verbose = FALSE), 'raw_data')
  expect_true(all(c('icc_se', 'icc_ci_lo', 'icc_ci_up') %in% colnames(raw)))
})

test_that('the extraction sheet and the guidance name the new columns', {
  sheet <- data_extraction_sheet(measure = 'icc', extension = 'data.frame',
                                 verbose = FALSE)
  expect_true(all(c('icc_se', 'icc_ci_lo', 'icc_ci_up') %in% colnames(sheet)))
  g <- summary(convert_df(data.frame(icc = 0.8), measure = 'icc', verbose = FALSE),
               guidance = TRUE)$es_guidance_crude
  expect_true(any(grepl('icc_se|icc_ci', g)))
})


# ---------------------------------------------------------------------------
# Reliability roadmap 1.1: the agreement-type SE is not merely 'approximate'.
#
# Measured (two-way DGP, ICC(2,1) = 0.80, k = 2, rater variance 50% of the
# non-subject budget, 5000 reps): 95% coverage 0.820 at n = 20, 0.667 at n = 50,
# 0.321 at n = 200, 0.144 at n = 1000 -- it DEGRADES with n, because ICC(2,1)
# inherits MSC's k-1 df while the reported SE shrinks like 1/sqrt(n). The old
# documented 'around 0.74-0.76' is not a range the estimator occupies at any n.
#
# Two remedies the roadmap offered were rejected by measurement:
#   * flipping the default to 'consistency' is NUMERICALLY INERT (the two share
#     one SE formula -- asserted below, and pinned at 1e-15 in tests_save);
#   * dropping the SE inside es_from_icc() would overturn a deliberate archived
#     convention and change every direct caller.
# So the route keeps computing (agreement_se), convert_df() can decline to POOL it
# (icc_agreement_se), and BOTH default to 'compute' in 2.1.0. That default is not a
# preference: tests_save/checked pins the computed value both at the route
# (test-icc.R:48, deliberately, with an honesty note) and through the pipeline
# (test-icc.R convert_df blocks), plus V31's exact [INFO] string (test-flags.R V25).
# Defaulting to 'drop' failed 6 archived assertions, which by this project's rule is
# a signal about the change, not about the tests. 'drop' therefore ships opt-in.
# ---------------------------------------------------------------------------

test_that('flipping icc_type would have been numerically inert', {
  for (p in list(c(0.8, 50, 2), c(0.5, 120, 3), c(0.95, 30, 5), c(-0.2, 80, 2))) {
    a <- es_from_icc(p[1], p[2], p[3], icc_type = 'agreement')
    cc <- es_from_icc(p[1], p[2], p[3], icc_type = 'consistency')
    expect_equal(a$icc, cc$icc, tolerance = 1e-15)
    expect_equal(a$icc_se, cc$icc_se, tolerance = 1e-15)
  }
})

test_that("es_from_icc() still computes by default: the archived convention holds", {
  rho <- 0.80; k <- 2; n <- 50
  got <- es_from_icc(icc = rho, n_sample = n, n_measurements = k, icc_type = 'agreement')
  expect_equal(got$icc_se, sqrt(2 * (1 + (k - 1) * rho)^2 / (k * (k - 1) * (n - 1))),
               tolerance = 1e-10)
})

test_that("agreement_se = 'drop' withholds the COMPUTED SE and nothing else", {
  d <- es_from_icc(0.80, 50, 2, icc_type = 'agreement', agreement_se = 'drop')
  expect_true(is.na(d$icc_se))
  expect_equal(d$icc, log(1 - 0.80), tolerance = 1e-10)   # estimate survives
  # consistency rows are untouched: for them the formula is exact at leading order
  expect_false(is.na(es_from_icc(0.80, 50, 2, icc_type = 'consistency',
                                 agreement_se = 'drop')$icc_se))
  # a stepped-down ICC(2,k) is computed with the same approximation, so it drops too
  expect_true(is.na(es_from_icc(0.90, 50, 5, icc_type = 'ICC(2,k)',
                                agreement_se = 'drop')$icc_se))
  expect_false(is.na(es_from_icc(0.90, 50, 5, icc_type = 'ICC(3,k)',
                                 agreement_se = 'drop')$icc_se))
  expect_error(es_from_icc(0.8, 50, 2, agreement_se = 'banana'), 'agreement_se')
})

test_that("a REPORTED SE or CI survives agreement_se = 'drop'", {
  # this is the whole point of doing 1.4 first: the drop targets the approximation,
  # not the row
  expect_equal(es_from_icc(0.80, 50, 2, icc_type = 'agreement', icc_se = 0.05,
                           agreement_se = 'drop')$icc_se,
               0.05 / (1 - 0.80), tolerance = 1e-10)
  expect_equal(es_from_icc(0.94, icc_ci_lo = 0.86, icc_ci_up = 0.98,
                           icc_type = 'agreement', agreement_se = 'drop')$icc_se,
               abs(log(1 - 0.98) - log(1 - 0.86)) / (2 * qnorm(0.975)),
               tolerance = 1e-10)
})

test_that('convert_df defaults to drop and offers a documented escape hatch', {
  d <- data.frame(icc = c(0.85, 0.72), n_sample = c(145, 67),
                  n_measurements = 2, icc_type = c('agreement', 'consistency'))
  drop <- summary(convert_df(d, measure = 'icc', icc_agreement_se = 'drop',
                             verbose = FALSE))
  comp <- summary(convert_df(d, measure = 'icc', icc_agreement_se = 'compute',
                             verbose = FALSE))
  expect_true(is.na(drop$se_crude[1]))        # agreement withheld
  expect_false(is.na(drop$se_crude[2]))       # consistency kept
  expect_false(is.na(comp$se_crude[1]))       # escape hatch restores it
  expect_equal(comp$se_crude[2], drop$se_crude[2], tolerance = 1e-10)
  # NOT "the point estimate is preserved": summary() reports no effect size for a
  # row it cannot weight, so the agreement row's es_crude is NA under 'drop' too.
  # That is the same thing a bare omega does (es_from_omega gives an estimate with
  # se = NA, and summary() still shows es = NA), so this is the package's existing
  # convention rather than something introduced here -- but it does mean the
  # Tier-1 V31 flag is the only trace such a row leaves.
  expect_true(is.na(drop$es_crude[1]))
  expect_false(is.na(comp$es_crude[1]))
  expect_equal(drop$es_crude[2], comp$es_crude[2], tolerance = 1e-10)
  # and the choice is recorded on the object
  expect_equal(attr(convert_df(d, measure = 'icc', verbose = FALSE),
                    'icc_agreement_se'), 'compute')          # the 2.1.0 default
  expect_equal(attr(convert_df(d, measure = 'icc', icc_agreement_se = 'drop',
                               verbose = FALSE), 'icc_agreement_se'), 'drop')
})

test_that('V31 is [UNUSUAL] and names the two ways out', {
  f <- summary(convert_df(data.frame(icc = 0.8, n_sample = 50, n_measurements = 2),
                          measure = 'icc', verbose = FALSE), flags = TRUE)$flags_crude
  # V31 keeps the [INFO] severity and the exact opening string tests_save pins
  # (test-flags.R V25), but now carries the MEASURED coverage and names both ways
  # out. Escalating it to [UNUSUAL] broke that pin, so the wording is load-bearing.
  expect_true(grepl('[INFO] ICC agreement-type SE', f[1], fixed = TRUE))
  expect_true(grepl('0.14 at n = 1000', f[1], fixed = TRUE))
  expect_true(grepl('icc_se', f[1], fixed = TRUE))
  expect_true(grepl("icc_agreement_se = 'drop'", f[1], fixed = TRUE))
  # a consistency row gets no such note
  f2 <- summary(convert_df(data.frame(icc = 0.8, n_sample = 50, n_measurements = 2,
                                      icc_type = 'consistency'),
                           measure = 'icc', verbose = FALSE), flags = TRUE)$flags_crude
  expect_false(grepl('ICC agreement-type SE', f2[1], fixed = TRUE))
})

test_that('the documented coverage claim is the measured one, not 0.74-0.76', {
  rd <- readLines('../../man/es_from_icc.Rd', warn = FALSE)
  skip_if(length(rd) == 0, 'Rd not readable from the test working directory')
  # the measured numbers are present...
  expect_true(any(grepl('0.820', rd, fixed = TRUE)))
  expect_true(any(grepl('0.144', rd, fixed = TRUE)))
  expect_true(any(grepl('degrades as n grows', rd, fixed = TRUE)))
  # ...and the old figure survives ONLY as a retraction, never as the claim
  expect_false(any(grepl('coverage around 0.74-0.76', rd, fixed = TRUE)))
  hit <- grep('0.74-0.76', rd, fixed = TRUE, value = TRUE)
  expect_true(all(grepl('Earlier versions', hit, fixed = TRUE)))
  # the \details block must exist at all: roxygen drops it silently on a brace or
  # quote error and check_man() still reports 'No issues detected'
  expect_true(any(grepl('Where the standard error comes from', rd, fixed = TRUE)))
  expect_true(any(grepl('Scale note', rd, fixed = TRUE)))
})
