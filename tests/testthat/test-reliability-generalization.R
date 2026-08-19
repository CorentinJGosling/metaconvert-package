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
  d <- data.frame(study_id = "x", icc = 0.995, n_sample = 100, n_measurements = 3)
  s <- suppressWarnings(summary(convert_df(d, measure = "icc", verbose = FALSE,
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
  expect_error(es_from_omega(0.8, omega_type = "bogus"), "omega_type")
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
