# Red half of red-green for the Tier-1 input-validation findings of
# AUDIT-2026-08-28-findings.md (subsystems: flags-tier1, alpha, icc).
#
# Every test_that() below asserts the CORRECT behaviour, so it FAILS on the
# tree as audited and PASSES once the defect is repaired. Nothing here is a
# "differs from the buggy value" test: where a number is asserted it is
# recomputed inline from the published closed form.
#
# Convention: convert_df(split_adjusted = FALSE) so the summary columns are the
# unsuffixed es / se / flags / es_guidance.

# ---------------------------------------------------------------------------
.tier1_summary <- function(d, ...) {
  suppressMessages(
    summary(convert_df(d, verbose = FALSE, split_adjusted = FALSE, ...),
            flags = TRUE, guidance = FALSE)
  )
}

.tier1_flag <- function(s, id) s$flags[match(id, s$study_id)]


# AUDIT line 676 -- .count_decimals() reads a bound printed as 0.90 as 1 decimal,
# so V2's rounding pad becomes 0.10 and an arithmetically INVERTED reliability CI
# ([0.90, 0.86]) survives validation and builds the row's inverse-variance weight,
# while the same 0.04 inversion one digit over ([0.91, 0.87]) is caught and NA'd.
test_that("AUDIT-676: V2 catches an inverted icc CI whose bounds end in a zero", {
  d <- data.frame(study_id = c("A", "B", "C"),
                  icc = 0.85, n_sample = 120, n_measurements = 2,
                  icc_ci_lo = c(0.90, 0.91, NA),
                  icc_ci_up = c(0.86, 0.87, NA))
  s <- .tier1_summary(d, measure = "icc")

  # Control: the non-round twin IS caught today, and must stay caught.
  expect_true(grepl("Inverted CI for 'icc'", .tier1_flag(s, "B"), fixed = TRUE))
  expect_true(is.na(s$se[match("B", s$study_id)]))

  # The defect: row A carries the identical 0.04 inversion and is silent.
  expect_true(grepl("Inverted CI for 'icc'", .tier1_flag(s, "A"), fixed = TRUE),
              info = "V2 must flag icc_ci = [0.90, 0.86] exactly as it flags [0.91, 0.87]")
  expect_true(is.na(s$es[match("A", s$study_id)]),
              info = "correct_inputs = TRUE must NA the inverted-CI row, as it does for B")
  expect_true(is.na(s$se[match("A", s$study_id)]),
              info = "an inverted CI must not be allowed to set the meta-analytic weight")

  # Row C (no CI) pins the closed form the surviving bad interval displaces:
  # Bonett-scale icc SE = sqrt(2 (1 + (k-1) rho)^2 / (k (k-1) (n-1))).
  k <- 2; n <- 120; rho <- 0.85
  se_closed <- sqrt(2 * (1 + (k - 1) * rho)^2 / (k * (k - 1) * (n - 1)))
  expect_equal(s$se[match("C", s$study_id)], se_closed, tolerance = 1e-8)
})


# AUDIT line 511 -- same .count_decimals() hole on V3 ("Value outside CI"): with
# alpha = 0.85 and a reported CI of [0.90, 0.95] the pad is 0.5e-2 + 0.5e-1 = 0.055,
# so the point estimate can sit 0.05 outside its own interval unflagged -- and
# es_from_ALPHA.R now prefers the CI over the closed form, so that self-contradictory
# interval fixes the study's weight (0.1768 vs the closed-form 0.1059).
test_that("AUDIT-511: V3 catches an alpha lying outside a CI whose bound ends in a zero", {
  d <- data.frame(study_id = "s2", cronbach_alpha = 0.85, n_sample = 200, n_items = 10,
                  cronbach_alpha_ci_lo = 0.90, cronbach_alpha_ci_up = 0.95)
  s <- .tier1_summary(d, measure = "alpha")

  expect_true(grepl("Value outside CI for 'cronbach_alpha'", s$flags, fixed = TRUE),
              info = "alpha = 0.85 is 0.05 below its own reported lower bound of 0.90")

  # Once the contradictory interval is refused, the SE must be Bonett's closed form
  # sqrt(2k / ((k-1)(n-2))) rather than the CI-derived 0.1768265 (a 0.36x weight).
  k <- 10; n <- 200
  se_closed <- sqrt(2 * k / ((k - 1) * (n - 2)))
  expect_equal(s$se, se_closed, tolerance = 1e-8)
  expect_equal(s$es, log(1 - 0.85), tolerance = 1e-8)
})


# AUDIT line 918 (+ Verifier correction) -- V40 skips only the "raw" scale, so under
# alpha_to_es = "hakstian_whalen" a cronbach_alpha of exactly 1 is RETAINED
# (es = 1, se = 0) while the flag tells the user it "yields no effect size and drops
# out of the pool", and blames "log(0)" for a transform that has no logarithm.
test_that("AUDIT-918: V40 does not claim a retained hakstian_whalen alpha = 1 was dropped", {
  d <- data.frame(study_id = c("s1", "s2", "s3"),
                  cronbach_alpha = c(0.85, 1, 0.90), n_sample = 200, n_items = 10)
  s <- .tier1_summary(d, measure = "alpha", alpha_to_es = "hakstian_whalen")

  flag <- .tier1_flag(s, "s2")
  es   <- s$es[match("s2", s$study_id)]

  # es_from_cronbach_alpha()'s se_defined gate refuses the variance at the boundary on
  # every scale but raw, so hakstian_whalen alpha = 1 DOES leave the pool -- it just
  # leaves with an estimate attached. V40 must therefore say so on this scale too.
  #
  # This was written as a disjunction, "either the row is dropped or the flag must not
  # say it was". Once the route started dropping the row that passed on the first
  # branch and stopped testing the flag at all, which is how V40 came to be silent here
  # without anything failing. Both halves are now asserted separately.
  expect_true(is.na(es),
              info = "hakstian_whalen alpha = 1 has no SE, so the row is not poolable")
  expect_true(grepl("drops out of the pool", flag, fixed = TRUE),
              info = "V40 must explain the drop on hakstian_whalen, not only on bonett")
  expect_true(grepl("raw", flag, fixed = TRUE),
              info = "the remedy must name raw, the only scale that keeps the boundary")
  expect_false(grepl("hakstian_whalen, where 1 is a usable boundary", flag, fixed = TRUE),
               info = "the old remedy sent the reader to the scale that swallows the row")

  # 1 - (1 - a)^(1/3) is defined at a = 1; the parenthetical is false on this scale.
  expect_false(grepl("log(0) is undefined", flag, fixed = TRUE),
               info = "the hakstian_whalen transform contains no logarithm")
})


# AUDIT line 931 (+ Verifier correction) -- valid_es admits alpha == 1 on every
# non-bonett scale, and the reported-SE delta map divides by 3 (1 - a)^(2/3) = 0, so
# the route emits alpha_se = Inf and a [-Inf, Inf] CI, breaking its own documented
# "degrade to NA rather than emit Inf or NaN" contract. (The verifier's repair: on
# non-raw scales a == 1 must give se = NA from ALL THREE SE routes -- never 0, which
# would be an infinite weight.)
test_that("AUDIT-931: alpha = 1 with a reported SE gives NA, not Inf, on hakstian_whalen", {
  res <- es_from_cronbach_alpha(1, 200, 10, cronbach_alpha_se = 0.02,
                                alpha_to_es = "hakstian_whalen")

  expect_true(is.na(res$alpha_se),
              info = "delta map 1/(3 (1-a)^(2/3)) is singular at a = 1; must degrade to NA")
  expect_true(is.na(res$alpha_ci_lo))
  expect_true(is.na(res$alpha_ci_up))

  # Control: the identity map on the raw scale is unaffected and must stay so.
  raw <- es_from_cronbach_alpha(1, 200, 10, cronbach_alpha_se = 0.02, alpha_to_es = "raw")
  expect_equal(raw$alpha_se, 0.02, tolerance = 1e-12)
})


# AUDIT line 944 (+ Verifier correction) -- V41 applies the SINGLE-measures floor
# -1/(k-1) to the raw icc column without consulting icc_type. For an average-measures
# ICC the inverse Spearman-Brown step-down sends every rho_k in [-1, 1] to a rho_1 at
# or above -1/(2k-1), so no in-range ICC(2,k) can violate the single-measures floor --
# yet the row is declared "[INVALID] Impossible ICC" and NA'd, destroying a value the
# route itself computes correctly.
test_that("AUDIT-944: V41 does not null a legal average-measures ICC", {
  d <- data.frame(study_id = "a", icc = -0.6, n_sample = 60, n_measurements = 3,
                  icc_type = "ICC(2,k)")
  s <- .tier1_summary(d, measure = "icc", icc_to_es = "raw")

  expect_false(grepl("Impossible ICC", s$flags, fixed = TRUE),
               info = "ICC(2,k) = -0.6 at k = 3 steps down to -0.1429, well above the -0.5 floor")

  # Oracle, recomputed from the two published closed forms the route composes:
  #   inverse Spearman-Brown  rho_1 = rho_k / (k - (k-1) rho_k)
  #   raw-scale SE            (1 - rho) sqrt(2 (1 + (k-1) rho)^2 / (k (k-1) (n-1)))
  k <- 3; n <- 60; rho_k <- -0.6
  rho_1 <- rho_k / (k - (k - 1) * rho_k)
  se_oracle <- (1 - rho_1) * sqrt(2 * (1 + (k - 1) * rho_1)^2 / (k * (k - 1) * (n - 1)))

  expect_equal(s$es, rho_1, tolerance = 1e-8)          # -0.1428571
  expect_equal(s$se, se_oracle, tolerance = 1e-8)      #  0.06135886

  # And the route, called directly, already agrees -- pinning that the pipeline,
  # not the arithmetic, is what has to change.
  direct <- es_from_icc(icc = -0.6, n_sample = 60, n_measurements = 3,
                        icc_type = "ICC(2,k)", icc_to_es = "raw")
  expect_equal(direct$icc, rho_1, tolerance = 1e-8)
  expect_equal(direct$icc_se, se_oracle, tolerance = 1e-8)
})


# AUDIT line 957 (+ Verifier correction) -- V40's raw-scale short-circuit rests on
# "on the raw scale 1.00 is a usable boundary (es = 1, se = 0)", true for alpha/omega
# but false for es_from_icc(), whose gate icc < 1 is scale-independent. So an ICC of
# exactly 1 under icc_to_es = "raw" vanishes from the pool with no flag; and on the
# DEFAULT bonett scale V40 fires but routes the user into that same silent failure by
# recommending the raw scale.
test_that("AUDIT-957: icc = 1 is explained on the raw scale and not misdirected on bonett", {
  d <- data.frame(study_id = c("a", "b"), icc = c(1, 0.8),
                  n_sample = 60, n_measurements = 2)

  s_raw <- .tier1_summary(d, measure = "icc", icc_to_es = "raw")
  expect_true(is.na(s_raw$es[match("a", s_raw$study_id)]))   # the route drops it (correct)
  expect_true(grepl("= 1 exactly", .tier1_flag(s_raw, "a"), fixed = TRUE),
              info = "the raw-scale drop of icc = 1 is currently unexplained (no V40 flag)")

  s_bon <- .tier1_summary(d, measure = "icc")
  expect_false(grepl("set the analysis scale to 'raw'", .tier1_flag(s_bon, "a"), fixed = TRUE),
               info = "the raw scale does not rescue icc = 1, so V40 must not recommend it here")
})


# AUDIT line 709 -- V31 is gated with icc_measure_ok, but V42 and V43, added in the
# same release, key only on "icc" %in% colnames(x). A COSMIN sheet carrying alpha and
# icc side by side, run at measure = "alpha", is therefore decorated with ICC pooling
# advice about a measure the run does not estimate.
test_that("AUDIT-709: an alpha run on a COSMIN sheet emits no ICC-named Tier-1 flag", {
  d <- data.frame(study_id = paste0("S", 1:4),
                  cronbach_alpha = c(.88, .91, .85, .90), n_items = 20,
                  n_sample = c(100, 120, 140, 160),
                  icc = c(.80, .85, .78, .82), n_measurements = c(2, 2, 2, 3),
                  icc_type = c("agreement", "agreement", "consistency", "agreement"))

  s_alpha <- .tier1_summary(d, measure = "alpha")
  expect_false(any(grepl("Pool mixes ICC estimands", s_alpha$flags, fixed = TRUE)),
               info = "V43 must be gated on measure the way V31 is")
  expect_false(any(grepl("Measurement count differs", s_alpha$flags, fixed = TRUE)),
               info = "V42 must be gated on measure the way V31 is")

  # Control: on the ICC run both checks are in scope and must keep firing.
  s_icc <- .tier1_summary(d, measure = "icc")
  expect_true(any(grepl("Pool mixes ICC estimands", s_icc$flags, fixed = TRUE)))
  expect_true(any(grepl("Measurement count differs", s_icc$flags, fixed = TRUE)))
})


# AUDIT line 698 -- .validate_input_data() takes no flag_group and convert_df() passes
# none, so V23/V35/V36/V37/V38/V39/V42/V43 are always computed over the whole dataset
# while the documentation exempts only V23. Here the two design strata each hold a
# constant k, so V42's own majority guard would keep it silent within group.
test_that("AUDIT-698: flag_group scopes the Tier-1 cross-row measurement-count check", {
  d <- data.frame(study_id = paste0("T", 1:6),
                  design = c(rep("test-retest", 4), rep("inter-rater", 2)),
                  icc = c(.80, .82, .84, .86, .75, .77),
                  n_measurements = c(2, 2, 2, 2, 4, 4),
                  n_sample = c(50, 60, 70, 80, 90, 100))

  # Control: ungrouped, k = 4 is the minority of the whole pool, so V42 fires.
  s_all <- .tier1_summary(d, measure = "icc")
  expect_true(any(grepl("Measurement count differs", s_all$flags, fixed = TRUE)))

  # Grouped by design, k is constant inside each stratum -> V42 must be silent.
  s_grp <- .tier1_summary(d, measure = "icc", flag_options = list(flag_group = "design"))
  expect_false(any(grepl("Measurement count differs", s_grp$flags, fixed = TRUE)),
               info = "V42 is computed over the whole dataset regardless of flag_group")
})


# AUDIT line 687 -- V26's implied critical value for a CORRECT paired CI is
# qt(.975, n-1) sqrt(1 - r), which re-enters the z/t band whenever the pre-post
# correlation is low, so the check false-fires on correctly analysed within-subject
# rows -- refuting its documented "needs no assumed r_pre_post" robustness.
test_that("AUDIT-687: V26 does not fire on a correctly-computed paired CI at low r", {
  n <- 10; r <- 0.30; sd_pre <- 10; sd_post <- 10

  # A textbook paired CI: half-width = qt(.975, n-1) * sd_diff / sqrt(n).
  se_paired <- sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post) / sqrt(n)
  half <- qt(.975, n - 1) * se_paired
  d <- data.frame(n_exp = n, n_nexp = n, mean_sd_exp = sd_pre, mean_sd_nexp = sd_post,
                  r_pre_post_exp = r, r_pre_post_nexp = r,
                  user_ci_lo_crude = 5 - half, user_ci_up_crude = 5 + half,
                  user_ci_lo_adj = NA_real_, user_ci_up_adj = NA_real_)
  iss <- metaConvert:::.validate_input_data(d, verbose = FALSE, measure = "mdw")$issues

  expect_false(any(grepl("Within-subject design", iss, fixed = TRUE)),
               info = "this CI IS the paired one; V26 must not call it an independent-groups SE")

  # True-positive control: a CI genuinely built with the independent-groups SE must
  # still fire, so the fix cannot simply be to weaken V26 out of existence.
  se_indep <- sqrt(sd_pre^2 / n + sd_post^2 / n)
  half_i <- qt(.975, 2 * n - 2) * se_indep
  d_bad <- d
  d_bad$user_ci_lo_crude <- 5 - half_i
  d_bad$user_ci_up_crude <- 5 + half_i
  iss_bad <- metaConvert:::.validate_input_data(d_bad, verbose = FALSE, measure = "mdw")$issues
  expect_true(any(grepl("Within-subject design", iss_bad, fixed = TRUE)))
})
