# Flag-message delimiter protocol: the per-row V-flags are joined with "; " in
# .validate_input_data() and split back on "; " in .flag_es_quality() (for the
# crude/adjusted scope routing). A literal "; " INSIDE a message therefore
# shatters it into fragments, and any fragment without a quoted 'column' is
# routed to BOTH scopes (e.g. the V23 tail used to leak into the adjusted
# column under split_adjusted = TRUE). Invariant tested here: every token
# obtained by splitting an emitted flag string on "; " is a complete flag,
# i.e. starts with an "[INVALID]"/"[UNUSUAL]"/"[INFO]"/"[DISCORDANT]" tag.

expect_clean_tokens <- function(s) {
  toks <- strsplit(s, "; ", fixed = TRUE)[[1]]
  ok <- grepl("^\\[(INVALID|UNUSUAL|INFO|DISCORDANT)\\]", toks)
  expect_true(all(ok), info = paste0(
    "message fragment without a leading [TAG] (internal '; ' in a flag): ",
    paste(toks[!ok], collapse = " || ")))
}

test_that("V15/V18/V21/V22/V23/V26 messages survive the '; ' split intact", {
  dat <- data.frame(
    study_id      = c("s1", "s2", "s3", "s4", "s5", "s6"),
    # row 1: V15 (n_exp + n_nexp != n_sample)
    n_exp         = c(10, 20, NA, NA, NA, 30),
    n_nexp        = c(10, 20, NA, NA, NA, 30),
    n_sample      = c(25, NA, NA, NA, NA, NA),
    # row 2: V18 (SD_baseline/SD_endpoint = 0.1) + V21 (baseline imbalance 2.5)
    mean_pre_exp     = c(NA, 10, NA, NA, NA, NA),
    mean_pre_nexp    = c(NA, 5, NA, NA, NA, NA),
    mean_pre_sd_exp  = c(NA, 2, NA, NA, NA, NA),
    mean_pre_sd_nexp = c(NA, 2, NA, NA, NA, NA),
    # rows 3+4: V23 (byte-identical endpoint block, different study_id)
    mean_exp      = c(NA, NA, 12.34, 12.34, NA, NA),
    mean_nexp     = c(NA, NA, 11.22, 11.22, NA, NA),
    mean_sd_exp   = c(NA, 20, 5.67, 5.67, NA, 10),
    mean_sd_nexp  = c(NA, 20, 4.56, 4.56, NA, 10),
    # row 5: V22 (ANCOVA residual SD, no cov_outcome_r column in the data)
    ancova_md_sd  = c(NA, NA, NA, NA, 5, NA),
    # row 6: V26 (mdw row whose CI half-width matches the independent-groups SE)
    user_ci_lo_crude = c(NA, NA, NA, NA, NA, -5.06),
    user_ci_up_crude = c(NA, NA, NA, NA, NA, 5.06)
  )
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE, measure = "mdw")
  issues <- res$issues

  # each targeted check actually fired
  expect_match(issues[1], "Sample size mismatch", fixed = TRUE)
  expect_match(issues[2], "SD_baseline / SD_endpoint", fixed = TRUE)
  expect_match(issues[2], "standardised baseline imbalance")
  expect_match(issues[3], "Identical endpoint summary statistics")
  expect_match(issues[5], "ANCOVA residual SD", fixed = TRUE)
  expect_match(issues[6], "Within-subject design", fixed = TRUE)

  # the invariant: no message shatters on the '; ' merge/split round-trip
  for (s in issues[nzchar(issues)]) expect_clean_tokens(s)
})

test_that("Tier-2 messages (H, B5, B7, D3) contain no internal '; '", {
  opts <- metaConvert:::.default_flag_options()

  # H: duplicate study_id
  f_h <- metaConvert:::.flag_cross_row_duplicates(c("a", "a"), opts)
  for (s in unlist(f_h)) expect_clean_tokens(s)

  # B5: NNT below the baseline-risk floor (1/0.2 = 5 > |2|)
  f_b <- metaConvert:::.flag_bounds_violations(
    es = 2, se = NA_real_, ci_lo = NA_real_, ci_up = NA_real_,
    measure = "nnt", exp = FALSE, baseline_risk = 0.2)
  expect_match(unlist(f_b), "Minimum possible NNT")
  for (s in unlist(f_b)) expect_clean_tokens(s)

  # B7: bonett-scale alpha > 0 (implied negative alpha)
  f_a <- metaConvert:::.flag_bounds_violations(
    es = 0.5, se = NA_real_, ci_lo = NA_real_, ci_up = NA_real_,
    measure = "alpha", exp = FALSE, alpha_to_es = "bonett")
  expect_match(unlist(f_a), "negative Cronbach's alpha")
  for (s in unlist(f_a)) expect_clean_tokens(s)

  # D3: SD spread outlier (md pool, one arm-SD pair ~20x below the median;
  # SEs kept unremarkable so D2 does not claim the row first)
  f_d <- metaConvert:::.flag_cross_row_outliers(
    es = c(1, 1.2, 0.9, 1.1, 1.0), se = c(2, 2.2, 1.9, 2.1, 2.0),
    opts = opts, measure = "md",
    n_total = rep(50, 5), n_exp = rep(25, 5), n_nexp = rep(25, 5),
    sd_exp = c(14, 15, 13, 14, 0.7), sd_nexp = c(14, 15, 13, 14, 0.7))
  expect_true(any(grepl("SD outlier", unlist(f_d))))
  for (s in unlist(f_d)) expect_clean_tokens(s)
})

test_that("V23 routes to the crude scope only (no unrouted fragment leaks to adjusted)", {
  dat <- data.frame(
    study_id     = c("s1", "s2"),
    mean_exp     = c(12.34, 12.34), mean_nexp    = c(11.22, 11.22),
    mean_sd_exp  = c(5.67, 5.67),   mean_sd_nexp = c(4.56, 4.56)
  )
  res <- metaConvert:::.validate_input_data(dat, verbose = FALSE)
  v23 <- strsplit(res$issues[1], "; ", fixed = TRUE)[[1]]
  v23 <- v23[grepl("Identical endpoint", v23)]
  expect_length(v23, 1)
  expect_true(metaConvert:::.v_flag_matches_scope(v23, "_crude"))
  expect_false(metaConvert:::.v_flag_matches_scope(v23, "_adjusted"))
})
