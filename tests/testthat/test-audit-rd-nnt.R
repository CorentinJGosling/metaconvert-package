## Red half of red-green for two findings of AUDIT-2026-08-28-findings.md.
##
##  * lines 434-443 -- measure = "rd" pools incidence rate differences with risk
##    differences. Fix (1) (generalising E4's gate from `measure == "nnt"` to
##    `measure %in% c("nnt", "rd")`) is ALREADY on disk at R/internal_flags.R:3753,
##    so only the two residuals are exercised here: B6/B6b still stamp `[INVALID]
##    RD outside [-1, 1]` on a `cases_time` row (an IRD is events per person-time,
##    unbounded, and [-1, 1] is not its parameter space), and the message E4 emits
##    inside an `rd` pool still calls the quantity an "NNT".
##
##  * lines 1241-1251 -- Tier-1 messages that interpolate their column names
##    UNQUOTED. `.v_flag_matches_scope()` (R/internal_flags.R:3480) returns TRUE for
##    BOTH suffixes when `regexpr("'[^']+'", flag_str)` finds nothing, so every such
##    message is duplicated into `flags_crude` AND `flags_adjusted` under the default
##    `split_adjusted = TRUE`, and each lands in the one scope it does not belong to.
##    Only the sites whose routing actually MOVES once the name is quoted are tested:
##    V12, V14, V17, V18, V21, V22 and the `r_pre_post` note built at
##    R/main_convert_df.R:588. The audit also lists V7/V15/V16 (and the double-zero
##    message), but those name columns in `.v_flag_matches_scope()`'s `shared` list
##    (n_exp, n_sample, n_cases_exp, ...), which route to BOTH scopes anyway -- so
##    quoting them is provably a no-op and there is nothing to assert.
##
## Every assertion below states the CORRECT behaviour and therefore fails today.

# ---------------------------------------------------------------------------
# helpers: the flags column is a single "; "-delimited string per row
# ---------------------------------------------------------------------------

.arn_tokens <- function(s) {
  if (length(s) != 1L || is.na(s) || !nzchar(s)) return(character(0))
  strsplit(s, "; ", fixed = TRUE)[[1]]
}

# tokens matching every pattern in `include` and none in `exclude`
.arn_pick <- function(s, include, exclude = character(0)) {
  tk <- .arn_tokens(s)
  if (!length(tk)) return(character(0))
  keep <- rep(TRUE, length(tk))
  for (p in include)  keep <- keep &  grepl(p, tk)
  for (p in exclude)  keep <- keep & !grepl(p, tk)
  tk[keep]
}

.arn_has <- function(s, include, exclude = character(0)) {
  length(.arn_pick(s, include, exclude)) > 0L
}

# ---------------------------------------------------------------------------
# fixture 1 -- an `rd` pool mixing a genuine 2x2 risk difference (row A) with a
# person-time incidence rate difference (row B), plus a directly entered rd > 1
# (row C) that must KEEP its B6 flag so a fix cannot be over-broad.
#
# Independently derived, from the counts alone:
#   A: RD  = 55/100 - 40/100 = 0.15
#      SE  = sqrt(.40*.60/100 + .55*.45/100) = sqrt(0.004875) = 0.0698212
#   B: IRD = 60/10 - 30/10 = 3 events per person-time unit
#      SE  = sqrt(n_cases_exp/time_exp^2 + n_cases_nexp/time_nexp^2)
#          = sqrt(30/100 + 60/100) = sqrt(0.9) = 0.9486833
# ---------------------------------------------------------------------------

.arn_rd_df <- data.frame(
  study_id        = c("A", "B", "C"),
  n_cases_exp     = c(40, 30, NA),
  n_controls_exp  = c(60, NA, NA),
  n_cases_nexp    = c(55, 60, NA),
  n_controls_nexp = c(45, NA, NA),
  time_exp        = c(NA, 10, NA),
  time_nexp       = c(NA, 10, NA),
  rd              = c(NA, NA, 1.5),
  rd_se           = c(NA, NA, 0.2),
  stringsAsFactors = FALSE
)

.arn_rd <- as.data.frame(
  summary(convert_df(.arn_rd_df, measure = "rd", verbose = FALSE), flags = TRUE)
)

# AUDIT-2026-08-28-findings.md:437,439,443 (fix 2)
# B6/B6b stamp "[INVALID] RD outside [-1, 1]" on a cases_time row, but that row
# carries an incidence rate DIFFERENCE (events per person-time), which is unbounded;
# an IRD of 3 events/person-year is arithmetically valid, not invalid.
test_that("AUDIT-rd-B6: B6/B6b do not fire on a cases_time (person-time) row", {
  i_rate <- which(.arn_rd$info_used_crude == "cases_time")
  i_rd   <- which(.arn_rd$info_used_crude == "rd_se")
  expect_length(i_rate, 1L)
  expect_length(i_rd, 1L)

  # anchor: the value B6 is objecting to is the correct IRD and its Poisson SE
  expect_equal(.arn_rd$es_crude[i_rate], 3, tolerance = 1e-8)
  expect_equal(.arn_rd$se_crude[i_rate], sqrt(0.9), tolerance = 1e-7)

  # the defect: an IRD is not a probability, so [-1, 1] is not its parameter space
  expect_false(.arn_has(.arn_rd$flags_crude[i_rate], "RD outside \\[-1, 1\\]"))
  expect_false(.arn_has(.arn_rd$flags_crude[i_rate], "RD CI bound outside \\[-1, 1\\]"))

  # guard against an over-broad fix: a genuine risk difference of 1.5 IS impossible
  # and must still be flagged
  expect_true(.arn_has(.arn_rd$flags_crude[i_rd], "RD outside \\[-1, 1\\]"))
})

# AUDIT-2026-08-28-findings.md:437,443 (fix 1, residual)
# E4's gate now reads `measure %in% c("nnt", "rd")`, but the message it emits was
# not adapted: inside a `measure = "rd"` pool it still calls both quantities "NNT"
# ("Mixed NNT types ... risk-based NNT ... rate-based (person-time) NNT"), naming a
# measure the user did not ask for and never mentioning the risk-difference /
# incidence-rate-difference distinction that is the actual hazard.
test_that("AUDIT-rd-E4: the mixing message in an rd pool talks about rate differences, not NNT", {
  i_rate <- which(.arn_rd$info_used_crude == "cases_time")
  mix <- .arn_pick(.arn_rd$flags_crude[i_rate], "should not be pooled together")
  expect_length(mix, 1L)

  expect_false(grepl("NNT", mix, fixed = TRUE))
  expect_true(grepl("rate difference", mix, ignore.case = TRUE))

  # guard: a measure = "nnt" pool must keep the NNT wording
  nnt <- as.data.frame(
    summary(convert_df(.arn_rd_df[1:2, ], measure = "nnt", verbose = FALSE),
            flags = TRUE)
  )
  j <- which(nnt$info_used_crude == "cases_time")
  nnt_mix <- .arn_pick(nnt$flags_crude[j], "should not be pooled together")
  expect_length(nnt_mix, 1L)
  expect_true(grepl("NNT", nnt_mix, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# fixture 2 -- one convert_df() call carrying all seven unquoted-name Tier-1 sites.
#
#   row A: pre/post + endpoint + an ANCOVA residual SD with no cov_outcome_r
#          -> V18 (SD_pre/SD_post = 2/5 = 0.40 < 0.70), V21 (|10-8|/2 = 1.00 > 0.30),
#             V22 (ancova_mean_sd_pooled without cov_outcome_r), r_pre_post note
#   row B: an extreme SD ratio on BOTH the crude pair (20x) and the ANCOVA pair (30x)
#          -> V12 twice, distinguishable by the printed ratio / the "ancova" prefix
#   row C: identical arm SDs and identical arm means on both pairs
#          -> V14 twice, V17 twice
#
# enable_informational = TRUE is required for V14/V17 only.
# ---------------------------------------------------------------------------

.arn_scope_df <- data.frame(
  study_id              = c("A", "B", "C"),
  n_exp                 = c(30, 30, 30),
  n_nexp                = c(30, 30, 30),
  mean_pre_exp          = c(10, NA, NA),
  mean_pre_nexp         = c(8,  NA, NA),
  mean_pre_sd_exp       = c(2,  NA, NA),
  mean_pre_sd_nexp      = c(2,  NA, NA),
  mean_exp              = c(12, 12, 7),
  mean_nexp             = c(9,  9,  7),
  mean_sd_exp           = c(5,  1,  3),
  mean_sd_nexp          = c(5,  20, 3),
  ancova_mean_exp       = c(12, 12, 11),
  ancova_mean_nexp      = c(9,  9,  11),
  ancova_mean_sd_pooled = c(4,  NA, NA),
  ancova_mean_sd_exp    = c(NA, 2,  4),
  ancova_mean_sd_nexp   = c(NA, 60, 4),
  cov_outcome_r         = c(NA, 0.5, 0.5),
  stringsAsFactors = FALSE
)

.arn_scope <- as.data.frame(
  summary(
    convert_df(.arn_scope_df, measure = "g", verbose = FALSE,
               flag_options = list(enable_informational = TRUE)),
    flags = TRUE, split_adjusted = TRUE
  )
)

.arn_crude <- function(i) .arn_scope$flags_crude[i]
.arn_adj   <- function(i) .arn_scope$flags_adjusted[i]

# AUDIT-2026-08-28-findings.md:1244,1246 (V12, R/internal_flags.R:1339)
# "Extreme SD ratio between arms: %s = %g, %s = %g" interpolates the pair's column
# names unquoted, so the crude-pair flag and the ANCOVA-pair flag are BOTH emitted
# into both scopes. The ANCOVA one belongs in flags_adjusted only, the mean_sd one
# in flags_crude only.
test_that("AUDIT-scope-V12: extreme-SD-ratio routes on the pair it names", {
  i <- 2L
  # crude pair (mean_sd_exp / mean_sd_nexp, ratio 20.0) -> crude only
  expect_true(.arn_has(.arn_crude(i), "Extreme SD ratio", exclude = "ancova"))
  expect_false(.arn_has(.arn_adj(i),  "Extreme SD ratio", exclude = "ancova"))
  # ANCOVA pair (ancova_mean_sd_exp / ancova_mean_sd_nexp, ratio 30.0) -> adjusted only
  expect_true(.arn_has(.arn_adj(i),   c("Extreme SD ratio", "ancova")))
  expect_false(.arn_has(.arn_crude(i), c("Extreme SD ratio", "ancova")))
})

# AUDIT-2026-08-28-findings.md:1244 (V14, R/internal_flags.R:1349-1364)
# Same unquoted-%s construction as V12: "Identical SDs across arms: %s = %s = %g".
test_that("AUDIT-scope-V14: identical-SDs routes on the pair it names", {
  i <- 3L
  expect_true(.arn_has(.arn_crude(i), "Identical SDs across arms", exclude = "ancova"))
  expect_false(.arn_has(.arn_adj(i),  "Identical SDs across arms", exclude = "ancova"))
  expect_true(.arn_has(.arn_adj(i),   c("Identical SDs across arms", "ancova")))
  expect_false(.arn_has(.arn_crude(i), c("Identical SDs across arms", "ancova")))
})

# AUDIT-2026-08-28-findings.md:1244 (V17, R/internal_flags.R:1366-1387)
# "Identical means across arms: %s = %s = %g", column names unquoted.
test_that("AUDIT-scope-V17: identical-means routes on the pair it names", {
  i <- 3L
  expect_true(.arn_has(.arn_crude(i), "Identical means across arms", exclude = "ancova"))
  expect_false(.arn_has(.arn_adj(i),  "Identical means across arms", exclude = "ancova"))
  expect_true(.arn_has(.arn_adj(i),   c("Identical means across arms", "ancova")))
  expect_false(.arn_has(.arn_crude(i), c("Identical means across arms", "ancova")))
})

# AUDIT-2026-08-28-findings.md:1244,1246 (V18, R/internal_flags.R:1406)
# "SD_baseline / SD_endpoint = %.2f" names no column at all, so the flag is copied
# into flags_adjusted. It is a crude-scope pre/post property (mean_pre_sd_* over
# mean_sd_*) and has nothing to say about the ANCOVA estimate.
test_that("AUDIT-scope-V18: baseline/endpoint SD-ratio note is crude-scope only", {
  i <- 1L
  pat <- "SD_baseline|Pre-post SMD formulas"
  expect_true(.arn_has(.arn_crude(i), pat))
  expect_false(.arn_has(.arn_adj(i), pat))
})

# AUDIT-2026-08-28-findings.md:1244,1246 (V21, R/internal_flags.R:1558)
# "|standardised baseline imbalance| = %.2f" names no column; it is computed from
# mean_pre_* / mean_pre_sd_*, i.e. crude scope. Its own advice is to prefer the
# ANCOVA-adjusted statistics, so duplicating it into flags_adjusted tells the
# reader to switch to the estimate they are already looking at.
test_that("AUDIT-scope-V21: baseline-imbalance note is crude-scope only", {
  i <- 1L
  pat <- "standardised baseline imbalance"
  expect_true(.arn_has(.arn_crude(i), pat))
  expect_false(.arn_has(.arn_adj(i), pat))
})

# AUDIT-2026-08-28-findings.md:1244,1246 (V22, R/internal_flags.R:1594)
# "ANCOVA residual SD provided but cov_outcome_r is missing" names both columns
# unquoted, so this ANCOVA-scope property is also emitted into flags_crude, where
# no ANCOVA back-transformation is being performed.
test_that("AUDIT-scope-V22: ANCOVA-residual-SD note is adjusted-scope only", {
  i <- 1L
  pat <- "ANCOVA residual SD provided but"
  expect_true(.arn_has(.arn_adj(i), pat))
  expect_false(.arn_has(.arn_crude(i), pat))
})

# AUDIT-2026-08-28-findings.md:1241-1250, site at R/main_convert_df.R:588
# "[INFO] Default r_pre_post = %s used (not provided by user)" is appended straight
# into validation$issues with the column name unquoted, so it too lands in both
# scopes. r_pre_post is consumed only by the pre/post, mean-change and paired
# routes (.rows_with_r_consuming_data) -- all crude; the ANCOVA scope's covariate
# correlation is cov_outcome_r, a different column entirely.
test_that("AUDIT-scope-rprepost: the defaulted-r_pre_post note is crude-scope only", {
  i <- 1L
  # the [DISCORDANT] CI-overlap token also mentions "the default r_pre_post"; the
  # note itself is the one that OPENS with [INFO] Default
  pat <- c("^\\[INFO\\] Default", "r_pre_post")
  expect_true(.arn_has(.arn_crude(i), pat))
  expect_false(.arn_has(.arn_adj(i), pat))
})
