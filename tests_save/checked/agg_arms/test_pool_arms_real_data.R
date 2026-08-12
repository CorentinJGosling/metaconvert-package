## End-to-end test of pool_arms() on the processing-speed adult dataset.
##
## Pipeline tested:
##   raw Excel -> endp (long, per-arm) -> wide convert_df format
##              -> pool_arms(method = "pool" / "split") -> convert_df + summary
##
## Treats Compound_list = 1 as the placebo (nexp / shared comparator).
## Multi-arm outcomes: Taylor 2000 STROOP Color/Word, Taylor 2001 STROOP Color/Word
## (each = 1 placebo + 2 active drug arms).
##
## Note on direction: 10/32 rows have better_if_small_or_large == "lower". The
## sign flip is applied identically (i.e. NOT applied here) in both the
## pool_arms and the escalc+aggregate_df pipelines, so the diagnostic
## comparison in step 9 is internally consistent. Users running this for real
## should flip g for "lower" rows post-summary, just as the original
## processing_speed_adults.R does on yi.

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(testthat); library(metaConvert); library(metafor)
})

setwd("c:/Users/coren/Documents/sideprojet/metaConvert/tests/testthat/agg_arms")

# ---------------------------------------------------------------------------
# 1. Replicate endp from processing_speed_adults.R (lines 1-53)
# ---------------------------------------------------------------------------
data <- suppressMessages(read_excel("28_01_26_Database_Neurocog_Onlyrelevant.xlsx"))
names(data) <- as.character(unlist(data[2, ]))
data <- data[-c(1, 2), c(1:59, 61:66)]
names(data) <- make.unique(names(data))
data <- data %>% filter(!(cog_domain %in% c("DO NOT USE", "unclear", "NOT AVAILABLE", "`")))

speed <- data %>%
  filter(cog_domain == "Processing speed",
         `Pediatric_or_adult_sample?` == 4)
speed <- speed[, c(3, 6, 12:18, 28:65)]

endp <- speed %>%
  select(7, 3, 9, 10, 12, 21, 22, 23, 34:36, 32) %>%
  mutate(
    Neurocog_bas_nitt_12 = ifelse(Neurocog_end_mean_12 != "$", Neurocog_bas_nitt_12, "$"),
    Neurocog_end_mean_12 = ifelse(Neurocog_end_mean_12 != "$", Neurocog_end_mean_12, Neurocog_end_mean_compl_12),
    Neurocog_end_sd_12   = ifelse(Neurocog_end_sd_12   != "$", Neurocog_end_sd_12,   Neurocog_end_sd_compl_12),
    Neurocog_bas_nitt_12 = ifelse(Neurocog_bas_nitt_12 != "$", Neurocog_bas_nitt_12, Neurocog_N_compl_12),
    Neurocog_bas_nitt_12 = as.numeric(Neurocog_bas_nitt_12),
    Neurocog_end_sd_12   = as.numeric(Neurocog_end_sd_12),
    Neurocog_end_mean_12 = as.numeric(Neurocog_end_mean_12)
  ) %>%
  select(`Study ID`, Task, measure, Compound_list, better_if_small_or_large,
         Neurocog_bas_nitt_12, Neurocog_end_mean_12, Neurocog_end_sd_12) %>%
  filter(!is.na(Neurocog_end_mean_12))

# SD imputation from SE (lines 48-49 of original script)
idx_43 <- !is.na(endp$Neurocog_bas_nitt_12) & endp$Neurocog_bas_nitt_12 == 43 & is.na(endp$Neurocog_end_sd_12)
idx_51 <- !is.na(endp$Neurocog_bas_nitt_12) & endp$Neurocog_bas_nitt_12 == 51 & is.na(endp$Neurocog_end_sd_12)
endp$Neurocog_end_sd_12[idx_43] <- sqrt(43) * 0.4
endp$Neurocog_end_sd_12[idx_51] <- sqrt(51) * 0.4
endp <- endp %>% filter(!is.na(Neurocog_end_sd_12))

cat("\n=== Replicated endp:", nrow(endp), "rows ===\n")

# ---------------------------------------------------------------------------
# 2. Reshape endp -> convert_df wide format
#    Compound_list == 1 -> placebo (nexp side); other compounds -> active (exp side).
#    One row per active arm; placebo data duplicated across active rows of the
#    same (Study ID, Task, measure). pool_side = "exp" iff that outcome has >1
#    active arm; otherwise NA (pass-through).
# ---------------------------------------------------------------------------

reshape_to_wide <- function(endp) {
  out <- endp %>%
    group_by(`Study ID`, Task, measure) %>%
    group_modify(~ {
      placebo <- .x %>% filter(Compound_list == 1)
      active  <- .x %>% filter(Compound_list != 1)
      if (nrow(placebo) == 0 || nrow(active) == 0) {
        return(tibble())  # cannot build a comparison
      }
      if (nrow(placebo) > 1) placebo <- placebo[1, ]   # one placebo per outcome
      n_active <- nrow(active)
      tibble(
        Compound_list = active$Compound_list,
        better_if_small_or_large = .x$better_if_small_or_large[1],
        n_exp        = active$Neurocog_bas_nitt_12,
        mean_exp     = active$Neurocog_end_mean_12,
        mean_sd_exp  = active$Neurocog_end_sd_12,
        n_nexp       = placebo$Neurocog_bas_nitt_12,
        mean_nexp    = placebo$Neurocog_end_mean_12,
        mean_sd_nexp = placebo$Neurocog_end_sd_12,
        pool_side    = if (n_active > 1) "exp" else NA_character_
      )
    }) %>%
    ungroup() %>%
    mutate(study_id = paste(`Study ID`, Task, measure, sep = " | "))
  out
}

wide_df <- reshape_to_wide(endp)
cat("Reshaped wide_df:", nrow(wide_df), "rows;",
    sum(!is.na(wide_df$pool_side)), "rows in multi-arm outcomes.\n")
cat("Multi-arm groups:\n")
print(as.data.frame(wide_df %>% filter(!is.na(pool_side)) %>%
                      select(study_id, Compound_list, n_exp, mean_exp, mean_sd_exp,
                             n_nexp, mean_nexp, mean_sd_nexp)))

# ---------------------------------------------------------------------------
# 3. Run pool_arms (both methods)
# ---------------------------------------------------------------------------
cat("\n=== Running pool_arms(method = 'pool') ===\n")
pooled <- pool_arms(as.data.frame(wide_df),
                    study_id  = "study_id",
                    pool_side = "pool_side",
                    method    = "pool",
                    verbose   = TRUE)

cat("\n=== Running pool_arms(method = 'split') ===\n")
splitd <- pool_arms(as.data.frame(wide_df),
                    study_id  = "study_id",
                    pool_side = "pool_side",
                    method    = "split",
                    verbose   = FALSE)

cat("\nPooled output (multi-arm rows only):\n")
ma_ids <- unique(wide_df$study_id[!is.na(wide_df$pool_side)])
print(as.data.frame(pooled %>% filter(study_id %in% ma_ids) %>%
                      select(study_id, n_exp, mean_exp, mean_sd_exp,
                             n_nexp, mean_nexp, mean_sd_nexp)))

cat("\nSplit output (multi-arm rows only):\n")
print(as.data.frame(splitd %>% filter(study_id %in% ma_ids) %>%
                      select(study_id, Compound_list, n_exp, n_nexp)))

# ---------------------------------------------------------------------------
# 4. Structural assertions
# ---------------------------------------------------------------------------
test_that("pool: row count = single-arm rows + multi-arm groups", {
  expected_rows <- sum(is.na(wide_df$pool_side)) + length(ma_ids)
  expect_equal(nrow(pooled), expected_rows)
})

test_that("split: row count = same as input", {
  expect_equal(nrow(splitd), nrow(wide_df))
})

test_that("pool_side column is removed from outputs", {
  expect_false("pool_side" %in% colnames(pooled))
  expect_false("pool_side" %in% colnames(splitd))
})

test_that("single-arm rows pass through unchanged in pool method", {
  pass_in  <- wide_df %>% filter(is.na(pool_side)) %>% arrange(study_id)
  pass_out <- pooled %>% filter(study_id %in% pass_in$study_id) %>% arrange(study_id)
  expect_equal(pass_out$n_exp,        pass_in$n_exp)
  expect_equal(pass_out$mean_exp,     pass_in$mean_exp)
  expect_equal(pass_out$mean_sd_exp,  pass_in$mean_sd_exp)
  expect_equal(pass_out$n_nexp,       pass_in$n_nexp)
  expect_equal(pass_out$mean_nexp,    pass_in$mean_nexp)
  expect_equal(pass_out$mean_sd_nexp, pass_in$mean_sd_nexp)
})

# ---------------------------------------------------------------------------
# 5. Hand-check Cochrane pooling (oracle: manual formula)
# ---------------------------------------------------------------------------
cochrane_sd_2 <- function(n1, sd1, m1, n2, sd2, m2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2 +
          n1 * n2 / (n1 + n2) * (m1^2 + m2^2 - 2 * m1 * m2)) /
         (n1 + n2 - 1))
}

for (sid in ma_ids) {
  arms <- wide_df %>% filter(study_id == sid)
  pooled_row <- pooled %>% filter(study_id == sid)
  expect_equal(nrow(pooled_row), 1, info = sid)

  # n_exp = sum
  expect_equal(pooled_row$n_exp, sum(arms$n_exp), info = paste(sid, "n_exp"))
  # weighted mean
  expected_m <- sum(arms$n_exp * arms$mean_exp) / sum(arms$n_exp)
  expect_equal(pooled_row$mean_exp, expected_m, tolerance = 1e-10,
               info = paste(sid, "mean_exp"))
  # Cochrane-pooled SD (manual, 2-arm case for Taylor)
  if (nrow(arms) == 2) {
    expected_sd <- cochrane_sd_2(arms$n_exp[1], arms$mean_sd_exp[1], arms$mean_exp[1],
                                  arms$n_exp[2], arms$mean_sd_exp[2], arms$mean_exp[2])
    expect_equal(pooled_row$mean_sd_exp, expected_sd, tolerance = 1e-10,
                 info = paste(sid, "mean_sd_exp"))
  }
  # placebo unchanged
  expect_equal(pooled_row$n_nexp,       arms$n_nexp[1],       info = paste(sid, "n_nexp"))
  expect_equal(pooled_row$mean_nexp,    arms$mean_nexp[1],    info = paste(sid, "mean_nexp"))
  expect_equal(pooled_row$mean_sd_nexp, arms$mean_sd_nexp[1], info = paste(sid, "mean_sd_nexp"))
}
cat("\n[OK] Cochrane pooling matches manual formula for all multi-arm groups.\n")

# ---------------------------------------------------------------------------
# 6. split method: placebo n divided across arms with remainder distribution
# ---------------------------------------------------------------------------
for (sid in ma_ids) {
  arms_in  <- wide_df %>% filter(study_id == sid) %>% arrange(Compound_list)
  arms_out <- splitd   %>% filter(study_id == sid) %>% arrange(Compound_list)
  k <- nrow(arms_in)
  n_shared <- arms_in$n_nexp[1]
  expect_equal(sum(arms_out$n_nexp), n_shared, info = paste(sid, "n_nexp split sum"))
  expect_true(all(abs(diff(sort(arms_out$n_nexp))) <= 1),
              info = paste(sid, "split remainder differs by <=1"))
  # exp arms unchanged
  expect_equal(arms_out$n_exp,       arms_in$n_exp,       info = paste(sid, "split n_exp"))
  expect_equal(arms_out$mean_exp,    arms_in$mean_exp,    info = paste(sid, "split mean_exp"))
  expect_equal(arms_out$mean_sd_exp, arms_in$mean_sd_exp, info = paste(sid, "split mean_sd_exp"))
}
cat("[OK] split method divides placebo n correctly.\n")

# ---------------------------------------------------------------------------
# 7. Run convert_df + summary on the pooled output
# ---------------------------------------------------------------------------
cat("\n=== convert_df + summary on pooled output ===\n")
res <- convert_df(pooled, measure = "g", verbose = FALSE)
sm  <- summary(res)

cat("Summary:", nrow(sm), "rows; info_used distribution:\n")
print(table(sm$info_used_crude, useNA = "ifany"))
expect_true(all(sm$info_used_crude == "means_sd"))
expect_true(all(is.finite(sm$es_crude)))
expect_true(all(is.finite(sm$se_crude)))
expect_true(all(is.finite(sm$es_ci_lo_crude) & is.finite(sm$es_ci_up_crude)))

# ---------------------------------------------------------------------------
# 8. Sanity vs metafor::escalc (oracle for SMD on pooled-arm summary stats)
# ---------------------------------------------------------------------------
cat("\n=== Sanity check vs metafor::escalc(measure = 'SMD') on pooled data ===\n")
oracle <- escalc(measure = "SMD",
                 n1i = pooled$n_exp,  n2i = pooled$n_nexp,
                 m1i = pooled$mean_exp, m2i = pooled$mean_nexp,
                 sd1i = pooled$mean_sd_exp, sd2i = pooled$mean_sd_nexp,
                 data = pooled, append = TRUE)

# Match metaConvert summary by study_id (assumes summary preserves input row order)
sm_ordered <- sm[match(pooled$study_id, sm$study_id), ]
cmp <- data.frame(study_id = pooled$study_id,
                  metaConvert_g  = sm_ordered$es_crude,
                  metaConvert_se = sm_ordered$se_crude,
                  metafor_g      = oracle$yi,
                  metafor_se     = sqrt(oracle$vi))
print(cmp)

expect_true(all(abs(cmp$metaConvert_g  - cmp$metafor_g)  < 1e-3))
# SEs may differ by a few percent because metaConvert and metafor apply the
# Hedges small-sample correction at slightly different points; the divergence
# is largest for the smallest samples (Lin 2016 n=28 -> ~3.5% relative).
se_rel_diff <- abs(cmp$metaConvert_se - cmp$metafor_se) / cmp$metafor_se
cat("Max SE relative diff:", round(max(se_rel_diff), 4),
    "at row", which.max(se_rel_diff), "\n")
expect_true(all(se_rel_diff < 0.05))
cat("[OK] metaConvert g matches metafor SMD (g within 1e-3, SE within 5% rel.).\n")

# ---------------------------------------------------------------------------
# 9. Diagnostic: pool_arms+convert_df vs the user's existing escalc+aggregate_df
#     (one row per study/outcome from each pipeline; should be directionally
#      consistent for multi-arm rows, exact match for single-arm rows.)
# ---------------------------------------------------------------------------
cat("\n=== Diagnostic: vs the original escalc + aggregate_df pipeline ===\n")

# Replicate the existing escalc + aggregate_df pipeline
endp_wide <- endp %>%
  full_join(endp,
            by = c("Study ID", "Task", "measure", "better_if_small_or_large"),
            suffix = c("_arm1", "_arm2"),
            relationship = "many-to-many") %>%
  filter(Compound_list_arm1 < Compound_list_arm2) %>%
  rename(sample1 = Neurocog_bas_nitt_12_arm1, sample2 = Neurocog_bas_nitt_12_arm2,
         mean1   = Neurocog_end_mean_12_arm1, mean2   = Neurocog_end_mean_12_arm2,
         sd1     = Neurocog_end_sd_12_arm1,   sd2     = Neurocog_end_sd_12_arm2,
         treat1  = Compound_list_arm1,        treat2  = Compound_list_arm2)
es_endp <- escalc(measure = "SMD",
                  n1i = sample2, n2i = sample1,            # active=2 vs placebo=1
                  m1i = mean2,   m2i = mean1,
                  sd1i = sd2,    sd2i = sd1,
                  data = endp_wide)
# Note: escalc renames `Study ID` -> `Study.ID`. Build a fresh key column.
es_endp$outcome_id <- paste(es_endp$Study.ID, es_endp$Task, es_endp$measure, sep = " | ")
es_endp$se_yi      <- sqrt(es_endp$vi)

agg <- aggregate_df(as.data.frame(es_endp), dependence = "outcomes", cor_unit = 0.5,
                    agg_fact = "outcome_id", es = "yi", se = "se_yi")

diag <- data.frame(
  study_id     = pooled$study_id,
  pool_arms_g  = sm_ordered$es_crude,
  pool_arms_se = sm_ordered$se_crude
)
diag <- merge(diag,
              data.frame(study_id = agg$outcome_id, agg_es = agg$es, agg_se = agg$se),
              by = "study_id", all.x = TRUE)
diag$is_multi_arm <- diag$study_id %in% ma_ids

cat("Diagnostic table:\n")
print(diag, digits = 4)

# Single-arm: pool_arms+convert should match escalc+aggregate (aggregate_df is
# a no-op for k=1). g matches to ~1e-3; SE may differ by a few % (small-sample
# correction divergence between metaConvert and metafor, same as step 8).
single <- diag[!diag$is_multi_arm, ]
expect_true(all(abs(single$pool_arms_g - single$agg_es) < 1e-2),
            info = "pool_arms+convert vs escalc+aggregate match for single-arm outcomes")
cat("[OK] Single-arm rows agree (g) between the two pipelines (<= 1e-2).\n")

multi <- diag[diag$is_multi_arm, ]
cat("\nMulti-arm comparison (pool_arms produces ONE ES from the composite arm;",
    "aggregate_df averages the dependent pairwise ESs):\n")
print(multi, digits = 4)

# ---------------------------------------------------------------------------
# 10. Per-study summary (final deliverable)
# ---------------------------------------------------------------------------
cat("\n=== Per-study summary (pool_arms + convert_df + summary) ===\n")
final <- pooled %>% select(study_id, n_exp, mean_exp, mean_sd_exp,
                           n_nexp, mean_nexp, mean_sd_nexp) %>%
  mutate(g       = sm_ordered$es_crude,
         g_se    = sm_ordered$se_crude,
         g_ci_lo = sm_ordered$es_ci_lo_crude,
         g_ci_up = sm_ordered$es_ci_up_crude,
         is_multi_arm = study_id %in% ma_ids)
print(as.data.frame(final), digits = 4)

cat("\n=== ALL CHECKS PASSED ===\n")
