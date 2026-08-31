## alpha_se_source / omega_se_source -- the pair of switches over the closed-form
## (n, k) reliability standard error.
##
## THE ASYMMETRY THESE EXIST TO FIX. alpha computes an SE from Bonett's
## 2k/((k-1)(n-2)) whenever the study reported none; omega refuses to and returns
## se = NA. The stated reason (?es_from_omega) is that Bonett's variance assumes
## essential tau-equivalence, which omega drops. Measured
## (simulations/studies/11_reliability_se.R, 11b) that rationale predicts nothing:
## the omega arm does not degrade as the loadings spread, and the tau-equivalent
## cell -- where the rationale says the formula IS valid -- is marginally the
## least well calibrated of the three. And 11a finds the two coefficients fail and
## succeed TOGETHER (|delta| <= 0.03 at k >= 8, against per-coefficient levels
## spanning 0.47 to 1.10).
##
## So the switches are symmetric: omega gains the ability to HAVE the closed form,
## alpha the ability to SKIP it. Neither default changes existing behaviour --
## that is the property most of this file pins.

test_that("neither default changes existing behaviour", {
  # omega: bare coefficient stays unpooled (se = NA), the fail-loud convention
  expect_true(is.na(es_from_omega(0.88, n_sample = 300, n_items = 8)$omega_se))
  # alpha: bare coefficient still gets the closed form
  expect_equal(es_from_cronbach_alpha(0.88, 300, 8)$alpha_se,
               sqrt(2 * 8 / ((8 - 1) * (300 - 2))), tolerance = 1e-12)
})


test_that("omega_se_source = 'closed_form' reproduces alpha's formula exactly", {
  # The whole claim is that this is the SAME formula, which is what makes the
  # published RG literature reproducible: metafor's measure = "ABT" with omega in
  # the alpha slot is what those papers ran.
  for (scale in c("bonett", "raw", "hakstian_whalen")) {
    w <- es_from_omega(0.88, n_sample = 300, n_items = 8,
                       omega_to_es = scale, omega_se_source = "closed_form")
    a <- es_from_cronbach_alpha(0.88, 300, 8, alpha_to_es = scale)
    expect_equal(w$omega_se, a$alpha_se, tolerance = 1e-12,
                 info = paste("scale:", scale))
    # and the point estimate was already identical -- only the SE route differed
    expect_equal(w$omega, a$alpha, tolerance = 1e-12, info = paste("scale:", scale))
  }
})


test_that("a study's own SE or CI still wins over the closed form", {
  # The switches govern the FALLBACK only. A reported SE must survive either way,
  # or turning the option on would silently discard study-reported precision.
  w <- es_from_omega(0.88, omega_se = 0.02, n_sample = 300, n_items = 8,
                     omega_se_source = "closed_form")
  w_ref <- es_from_omega(0.88, omega_se = 0.02, n_sample = 300, n_items = 8)
  expect_equal(w$omega_se, w_ref$omega_se, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(w$omega_se, sqrt(2 * 8 / (7 * 298)))))

  a <- es_from_cronbach_alpha(0.88, 300, 8, cronbach_alpha_se = 0.02,
                              alpha_se_source = "reported")
  expect_true(is.finite(a$alpha_se))

  # a reported CI likewise
  a_ci <- es_from_cronbach_alpha(0.88, 300, 8, cronbach_alpha_ci_lo = 0.84,
                                 cronbach_alpha_ci_up = 0.92,
                                 alpha_se_source = "reported")
  expect_true(is.finite(a_ci$alpha_se))
})


test_that("alpha_se_source = 'reported' refuses the closed form", {
  # The consequential half: alpha's closed form is ON by default and inherits
  # Bonett's normality assumption, which measured on dichotomous items is ~37%
  # too small (coverage .78, not improving with n). A dichotomous-instrument
  # review needs a way to decline it.
  a <- es_from_cronbach_alpha(0.88, 300, 8, alpha_se_source = "reported")
  expect_true(is.na(a$alpha_se))
  # the POINT estimate survives, so the row stays visible and countable
  expect_true(is.finite(a$alpha))
  expect_equal(a$alpha, log(1 - 0.88), tolerance = 1e-12)
})


test_that("both switches are honoured through convert_df()", {
  se_of <- function(d, m, ...) {
    s <- suppressMessages(summary(convert_df(d, measure = m, verbose = FALSE, ...)))
    if ("se_crude" %in% names(s)) s$se_crude else s$se
  }
  dw <- data.frame(study_id = c("A", "B"), omega = c(.90, .88),
                   n_sample = c(300, 250), n_items = 8)
  da <- data.frame(study_id = c("A", "B"), cronbach_alpha = c(.90, .88),
                   n_sample = c(300, 250), n_items = 8)

  expect_true(all(is.na(se_of(dw, "omega"))))
  expect_true(all(is.finite(se_of(dw, "omega", omega_se_source = "closed_form"))))
  expect_true(all(is.finite(se_of(da, "alpha"))))
  expect_true(all(is.na(se_of(da, "alpha", alpha_se_source = "reported"))))

  # and the two agree, which is the point of offering the omega route at all
  expect_equal(se_of(dw, "omega", omega_se_source = "closed_form"),
               se_of(da, "alpha"), tolerance = 1e-12)
})


test_that("an unknown level is rejected rather than silently ignored", {
  # A typo must not fall through to the default: that would silently reinstate
  # exactly the behaviour the caller asked to change.
  expect_error(es_from_omega(0.88, n_sample = 300, n_items = 8,
                             omega_se_source = "bonett_nk"))
  expect_error(es_from_cronbach_alpha(0.88, 300, 8, alpha_se_source = "none"))
})


test_that("the closed-form route keeps its own (n, k) guards under the switch", {
  # n <= 2 makes (n - 2) non-positive and k < 2 leaves the coefficient undefined.
  # Turning the option ON must not bypass either.
  expect_true(is.na(es_from_omega(0.88, n_sample = 2, n_items = 8,
                                  omega_se_source = "closed_form")$omega_se))
  expect_true(is.na(es_from_omega(0.88, n_sample = 300, n_items = 1,
                                  omega_se_source = "closed_form")$omega_se))
  # and with no n at all it simply stays NA rather than erroring
  expect_true(is.na(es_from_omega(0.88, n_items = 8,
                                  omega_se_source = "closed_form")$omega_se))
})


test_that("the closed form is vectorised and row-wise under the switch", {
  w <- es_from_omega(c(0.90, 0.88, 0.85),
                     omega_se = c(NA, 0.02, NA),
                     n_sample = c(300, 250, 400), n_items = c(8, 8, 12),
                     omega_se_source = "closed_form")
  expect_length(w$omega_se, 3L)
  # rows 1 and 3 from the closed form, row 2 from its own reported SE
  expect_equal(w$omega_se[1], sqrt(2 * 8 / (7 * 298)), tolerance = 1e-12)
  expect_equal(w$omega_se[3], sqrt(2 * 12 / (11 * 398)), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(w$omega_se[2], sqrt(2 * 8 / (7 * 248)))))
})
