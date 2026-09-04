# Red half of red-green for four findings of AUDIT-2026-08-28-findings.md, subsystem
# `2x2-core`. Every block here is expected to FAIL until the corresponding defect is
# fixed; none of them touches R/ or man/.
#
#   audit 379-388  es_from_chisq() / es_from_chisq_pval() switch correlation estimand
#                  silently -- the 2.1.0 disclosure went into es_from_phi() only.
#   audit 390-399  the +0.5 continuity correction reaches the tetrachoric solve.
#   audit 401-410  the mvtnorm-absent fallback overwrites correct d/g and misstates the
#                  cause.
#   audit 736-747  es_from_2x2()/_sum()/_prop() have no guard on impossible cells.

# .notify_once() is session-scoped, so a block that asserts on a one-time notice must
# start from a clean slate or it silently tests nothing.
clear_notices <- function() {
  keys <- ls(envir = metaConvert:::.mcv_notices, all.names = TRUE)
  if (length(keys)) rm(list = keys, envir = metaConvert:::.mcv_notices)
  invisible(NULL)
}


# --------------------------------------------------------------------------------
# audit 379-388 -- es_from_chisq() carries a byte-for-byte copy of the es_from_phi()
# fallback (R/es_from_PHI_CHISQ.R:243) with NO disclosure, so one call returns the
# TETRACHORIC correlation for the rows that have n_cases/n_exp and the phi coefficient
# itself for the rows that do not, both labelled info_used = "chisq", in silence.
# --------------------------------------------------------------------------------
test_that("AUDIT-major (379): es_from_chisq() discloses the estimand switch, as es_from_phi() does", {
  skip_if_not_installed("mvtnorm")

  # chisq = n * phi^2 = 200 * 0.2182^2 = 9.52216, i.e. the same data as the phi call
  # below; row 1 has the margins, row 2 does not.
  chisq_mixed <- function() {
    es_from_chisq(chisq = c(9.52216, 9.52216), n_sample = c(200, 200),
                  n_cases = c(60, NA), n_exp = c(100, NA))
  }

  # Reference: the phi entry point on the same data DOES warn. (Passes today.)
  clear_notices()
  expect_message(
    es_from_phi(phi = c(0.2182, 0.2182), n_sample = c(200, 200),
                n_cases = c(60, NA), n_exp = c(100, NA)),
    "TETRACHORIC", fixed = TRUE)

  # The two chisq rows really are two different quantities: row 1 is the tetrachoric
  # (metafor RTET on the reconstructed table), row 2 is phi = sqrt(chisq/n) = 0.2182.
  # 0.3548452 / 0.2182 = 1.626. (Passes today -- it is the reason the notice is owed.)
  clear_notices()
  res <- suppressMessages(chisq_mixed())
  expect_equal(res$r[2], sqrt(9.52216 / 200), tolerance = 1e-6)
  expect_gt(abs(res$r[1] / res$r[2]), 1.4)

  # THE DEFECT: no message at all. Fails today (0 messages emitted).
  clear_notices()
  expect_message(chisq_mixed(), "(?i)tetrachoric", perl = TRUE)

  # es_from_chisq_pval() delegates to es_from_chisq() and inherits the silence.
  # pchisq(9.52216, 1, lower.tail = FALSE) = 0.002030055
  clear_notices()
  expect_message(
    es_from_chisq_pval(chisq_pval = c(0.002030055, 0.002030055), n_sample = c(200, 200),
                       n_cases = c(60, NA), n_exp = c(100, NA)),
    "(?i)tetrachoric", perl = TRUE)

  # A call in which NO row can be reconstructed is equally undisclosed: every r is phi,
  # not the tetrachoric the documentation promises.
  clear_notices()
  expect_message(es_from_chisq(chisq = 9.52216, n_sample = 200), "(?i)tetrachoric",
                 perl = TRUE)
})


# --------------------------------------------------------------------------------
# audit 390-399 -- es_from_2x2() adds 0.5 to all four cells IN PLACE (R/es_from_2x2.R:
# 120-123) and then builds dat2x2 from the corrected cells (:156-161), so the
# tetrachoric solve sees the shrunk table. man/es_from_2x2.Rd:47-48 says the adjustment
# is "used for the OR/RR only", and the @details claims metafor parity for the
# correlation -- metafor deliberately does NOT apply add/to to RTET.
# --------------------------------------------------------------------------------
test_that("AUDIT-major (390): the +0.5 correction does not reach the tetrachoric route", {
  skip_if_not_installed("mvtnorm")

  # Independent expectation: metafor::escalc(measure = "RTET") on the RAW table, which
  # is the contract tests_save/checked/test-2x2.R pins at 1e-10.
  #   a/b/c/d = 0/20/10/10 -> yi = -1.0000000, vi = 22209.727125, sqrt(vi) = 149.029283
  raw1 <- metafor::escalc(measure = "RTET", ai = 0, bi = 20, ci = 10, di = 10)
  res1 <- es_from_2x2(n_cases_exp = 0, n_controls_exp = 20,
                      n_cases_nexp = 10, n_controls_nexp = 10)

  # Today: r = -0.8564096, r_se = 0.1134748 -- exactly escalc() on (0.5, 20.5, 10.5,
  # 10.5), i.e. the corrected table. The inverse-variance weight is 77.7 instead of
  # 4.5e-05, a factor of 1313.
  expect_equal(res1$r, as.numeric(raw1$yi), tolerance = 1e-8)
  expect_equal(res1$r_se, sqrt(as.numeric(raw1$vi)), tolerance = 1e-6)

  #   a/b/c/d = 0/50/8/42 -> yi = -1, sqrt(vi) = 4653.163 (mcv today: 0.17063)
  raw2 <- metafor::escalc(measure = "RTET", ai = 0, bi = 50, ci = 8, di = 42)
  res2 <- es_from_2x2(n_cases_exp = 0, n_controls_exp = 50,
                      n_cases_nexp = 8, n_controls_nexp = 42)
  expect_equal(res2$r, as.numeric(raw2$yi), tolerance = 1e-8)
  expect_equal(res2$r_se, sqrt(as.numeric(raw2$vi)), tolerance = 1e-6)

  # The correction must still apply to the OR arm, which is what the Rd documents.
  # log((0.5 * 10.5) / (20.5 * 10.5)) = log(0.5 / 20.5) = -3.7135721. Passes today and
  # must keep passing: the fix is to stop feeding the corrected cells to the solve, not
  # to drop the correction. The RD arm takes the same correction as the OR/RR arms
  # (metafor::escalc(measure = "RD") default), so it is pinned to metafor, not raw.
  expect_equal(res1$logor, log(0.5 / 20.5), tolerance = 1e-8)
  rd_m <- metafor::escalc(measure = "RD", ai = 0, bi = 20, ci = 10, di = 10)
  expect_equal(res1$rd, -as.numeric(rd_m$yi), tolerance = 1e-12)
  expect_equal(res1$rd_se, sqrt(as.numeric(rd_m$vi)), tolerance = 1e-12)
})


# --------------------------------------------------------------------------------
# audit 401-410 -- the fallback gate keys on the SYMPTOM (is.na(es$r)) rather than the
# CAUSE (the table could not be rebuilt), so with mvtnorm absent -- CRAN's own
# noSuggests flavour -- a row whose 2x2 WAS reconstructed has its correct d/g (built
# from the log OR, which never needed mvtnorm) overwritten by the phi-as-Pearson
# fallback, and the notice then states a cause that is false.
# tests/testthat/test-mvtnorm-guard.R exercises this path but never checks d.
# --------------------------------------------------------------------------------
test_that("AUDIT-major (401): a missing mvtnorm leaves d/g untouched and states the true cause", {
  local_mocked_bindings(.has_mvtnorm = function() FALSE, .package = "metaConvert")
  clear_notices()

  # Independent expectation, from the reconstructed table only (no mvtnorm anywhere in
  # this derivation):
  #   conv.2x2(ri = c(.3,.5), ni = c(120,200), n1i = c(40,100), n2i = c(20,60))
  #     -> a/b/c/d = 13/27/7/73 and 53/47/7/93
  #   logOR   = log(ad/bc)                       = 1.6136617836, 2.7068336559
  #   d       = logOR * sqrt(3)/pi               = 0.8896583687, 1.4923556096
  #   d_se    = sqrt((1/a+1/b+1/c+1/d) * 3/pi^2) = 0.2867524535, 0.2426813974
  #   g       = d * J(df = N - 2), exact gamma J = 0.8839897457, 1.4866944070
  ct <- suppressWarnings(metafor::conv.2x2(ri = c(0.3, 0.5), ni = c(120, 200),
                                           n1i = c(40, 100), n2i = c(20, 60)))
  a <- ct$ai; b <- ct$bi; cc <- ct$ci; dd <- ct$di
  logor_exp <- log(a * dd / (b * cc))
  d_exp     <- logor_exp * sqrt(3) / pi
  d_se_exp  <- sqrt((1 / a + 1 / b + 1 / cc + 1 / dd) * 3 / pi^2)
  df_t      <- (a + b + cc + dd) - 2
  J         <- exp(lgamma(df_t / 2) - lgamma((df_t - 1) / 2)) * sqrt(2 / df_t)

  msgs <- character(0)
  res <- withCallingHandlers(
    es_from_phi(phi = c(0.3, 0.5), n_sample = c(120, 200),
                n_cases = c(20, 60), n_exp = c(40, 100)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })

  # The table WAS reconstructed -- logor proves it and is unaffected by the mock.
  expect_equal(res$logor, logor_exp, tolerance = 1e-8)

  # THE DEFECT: d falls 0.8896584 -> 0.4962261 (-44%) and g 0.8839897 -> 0.4930643,
  # purely because a correlation the row never asked for could not be solved.
  expect_equal(res$d, d_exp, tolerance = 1e-8)
  expect_equal(res$g, d_exp * J, tolerance = 1e-8)
  expect_equal(res$d_se, d_se_exp, tolerance = 1e-8)

  # r/z are the only outputs mvtnorm is needed for: they must be NA, not silently
  # replaced by phi itself (0.3 / 0.5), which is a different estimand.
  expect_true(all(is.na(res$r)))
  expect_true(all(is.na(res$z)))

  # The honest notice (mvtnorm absent) must stand; the fabricated cause must not fire --
  # n_cases and n_exp were both supplied for every row.
  expect_true(any(grepl("mvtnorm", msgs, fixed = TRUE)))
  expect_false(any(grepl("could not be reconstructed", msgs, fixed = TRUE)))

  # Same gate, same consequence, at R/es_from_PHI_CHISQ.R:243.
  clear_notices()
  chi <- suppressMessages(es_from_chisq(chisq = 9.52216, n_sample = 200,
                                        n_cases = 60, n_exp = 100))
  expect_equal(chi$d, chi$logor * sqrt(3) / pi, tolerance = 1e-8)
  expect_true(is.na(chi$r))
})


# --------------------------------------------------------------------------------
# audit 736-747 (incl. the Verifier correction, which widens the scope to es_from_2x2()
# itself and to the DEFAULT convert_df() path) -- no guard on impossible cells. The RR
# arm uses only the row margins, so a table that cannot exist still yields a finite,
# correctly-signed logRR with a finite SE and CI.
# --------------------------------------------------------------------------------
test_that("AUDIT-major (736): impossible 2x2 cells are neutralised, not converted", {
  # (1) es_from_2x2_sum(): n_controls_exp = 40 - 50 = -10.
  # Today: logrr = log((50/40)/(10/40)) = 1.609438, logrr_se = 0.2645751, rd = -1,
  # nnt = -1 -- an "exposed risk" of 1.25.
  s <- suppressWarnings(es_from_2x2_sum(n_cases_exp = 50, n_exp = 40,
                                        n_cases_nexp = 10, n_nexp = 40))
  expect_true(is.na(s$logrr))
  expect_true(is.na(s$logrr_se))
  expect_true(is.na(s$logor_se))
  expect_true(is.na(s$rd))
  expect_true(is.na(s$nnt))

  # (2) es_from_2x2() itself is exported and equally unguarded (Verifier correction).
  # Today: logrr = 0.7178398 (se 0.1561738), rd = -0.525 (se 0.07489576),
  # nnt = -1.904762 -- and |rd| < 1, so flag B6 stays silent too.
  z <- suppressWarnings(es_from_2x2(n_cases_exp = 41, n_controls_exp = -1,
                                    n_cases_nexp = 20, n_controls_nexp = 20))
  expect_true(is.na(z$logrr))
  expect_true(is.na(z$logrr_se))
  expect_true(is.na(z$rd))
  expect_true(is.na(z$rd_se))
  expect_true(is.na(z$nnt))

  # (3) es_from_2x2_prop(): a proportion entered as a percentage. round(1.2 * 100) = 120
  # cases out of 100, so n_controls_exp = -20.
  # Today: logrr = 3.178054, logrr_se = 0.4339739, rd = -1.15, nnt = -0.8695652.
  p <- suppressWarnings(es_from_2x2_prop(prop_cases_exp = 1.2, prop_cases_nexp = 0.05,
                                         n_exp = 100, n_nexp = 100))
  expect_true(is.na(p$logrr))
  expect_true(is.na(p$logrr_se))
  expect_true(is.na(p$rd))
  expect_true(is.na(p$nnt))

  # (4) The prop route is exposed on the DEFAULT convert_df() path: the impossible count
  # is manufactured inside the route by round(prop * n), so no column-keyed Tier-1 check
  # sees it. Today: es_crude = 24, se_crude = 0.4339739, info_used_crude = "2x2_prop",
  # flags_crude = "" -- selected and pooled as an ordinary study.
  dat <- data.frame(prop_cases_exp = 1.2, prop_cases_nexp = 0.05,
                    n_exp = 100, n_nexp = 100)
  out <- suppressWarnings(suppressMessages(
    summary(convert_df(dat, measure = "rr", verbose = FALSE), flags = TRUE)))
  # the guard contains ...
  expect_true(is.na(out$es_crude[1]))
  # ... and the flag reports (prop_cases_exp/_nexp belong in .bounded_columns, beside
  # baseline_risk, whose entry exists for this exact percentage mix-up).
  flag_txt <- paste(unlist(out[1, grepl("^flags", names(out)), drop = FALSE]),
                    collapse = " ")
  expect_match(flag_txt, "prop_cases_exp", fixed = TRUE)
})
