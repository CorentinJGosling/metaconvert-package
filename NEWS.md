# metaConvert 2.1.0

## Breaking changes (pre/post and change-score effect sizes)

- **`pool_sd` now defaults to `TRUE`.** For the two-group pre/post and mean-change
  converters (`es_from_means_sd/se/ci_pre_post()`, `es_from_mean_change_sd/se/ci/pval()`)
  and for `convert_df()`, the between-group SMD is now the difference in mean change
  divided by a **single standardizing SD pooled across arms** (Morris, 2008).
  Previously each arm's change was standardized by that arm's **own** SD and the two
  within-group values were subtracted. That difference is a valid between-group SMD
  only when the two arms' SDs are equal; when they differ it is biased (in a simulated
  case with arm change-SDs 6.5 vs 10.7 it overstated the effect by ~60%). Set
  `pool_sd = FALSE` to restore the previous behaviour.
- **`es_from_paired_t()` / `es_from_paired_f()` cannot pool** and are unchanged: a paired
  t identifies each arm's `mean_change / sd_change` ratio but not the ratio of the two
  arms' SDs, so the pooled standardizer is not recoverable from the reported statistic.
  These routes keep the per-arm construction; `summary(..., flags = TRUE)` now raises an
  informational flag when such rows share a pool with pooled-standardizer rows.

## Bug fixes

- **Pooled `morris_dz` variance.** The pooled change-SD-standardized SMD
  (`pre_post_to_smd = "morris_dz"`, `pool_sd = TRUE`) reused the `morris_drm` variance,
  whose `2(1 - r)` factor belongs only to estimators whose point estimate carries the
  `sqrt(2(1 - r))` raw-score rescaling. At r = 0.6 this understated the standard error by
  ~11% (95% CI coverage 92-93% instead of 95%); below r = 0.5 it overstated it. The
  variance is now the two-sample Hedges (1981) form and matches
  `metafor::escalc(measure = "SMD")` on change scores.
- **Paired-t `morris_dz` variance.** `es_from_paired_t()` and
  `es_from_paired_t_single_group()` built the variance on the uncorrected d, so
  `var(g) = J^2/n + g^2/(2n)`, disagreeing with the package's own pre/post kernel and with
  `metafor::escalc(measure = "SMCC")`. The same study entered as a paired t or as a mean
  change now returns the same standard error.
- **Per-row `pre_post_to_smd`.** A vector mixing `"morris_dz"` and `"morris_drm"` rows was
  applied wholesale as `morris_drm` in the paired-t routes; the method is now selected
  per row.
- **Degenerate inputs.** A standardizing SD of zero (or a non-finite one) and a pre-post
  correlation with `|r| >= 1` now yield `NA` instead of `Inf`/`NaN`.

## New

- **Estimand-mixing flags** in `summary(..., flags = TRUE)` for `measure = "d"/"g"/"dw"/"gw"`,
  both `[INFO]` and both gated by `enable_cross_row` (they are properties of the pool, not
  of any one row -- every row may be individually correct, so they are not `[DISCORDANT]`):
  one when change-SD-standardized rows (`morris_dz`) share a pool with raw-score-SD rows,
  which the Cochrane Handbook (v6, §10.5.2) advises against combining; and one for the
  per-arm-standardized paired-t rows described above. Note that `morris_drm` rows mixed
  with endpoint rows raise **nothing**: `d_rm` is precisely the transformation onto the
  raw-score metric, so that combination is the legitimate one.
- **`convert_df()` now reports its silent substitutions** (when `verbose = TRUE`): the
  coercion of `pre_post_to_smd = "bonett"/"morris_dav"` to `"cooper"` for mean-change and
  paired data, and the imputation of an unreported `r_pre_post` for rows that consume it
  (under the default `"cooper"` standardizer the assumed correlation scales the point
  estimate, not just the standard error).

## Documentation

- Corrected the `@details` of the four single-group mean-change converters, which
  documented a sign-flipped argument mapping (`mean_pre_exp = mean_change_exp`,
  `mean_post = 0`) that was the reverse of the code.

## Tests

- The pre/post and change-score test suite (21 files) is now wired into `tests/testthat/`
  and runs under `R CMD check`; previously only `test-ANCOVA-MD.R` ran. Long-running files
  are gated on `NOT_CRAN`.

# metaConvert 2.0.0
- Implemented the metaDETECT framework of automated checks: a quality-flag system in summary() ('flags' argument) that detects inconsistent input data and implausible effect sizes
- Added new effect size measures
- Added psychometric and regression effect size conversions
- Added a guidance system in summary() ('guidance' argument) that indicates which columns are missing when an effect size cannot be estimated
- A few other improvements in various functions

# metaConvert 1.1.0
- Added single-group pre/post converters and the within-group measures dw, gw and mdw
- Added the `pool_sd` argument to the two-group pre/post converters (SD pooled across arms)
- Improved the handling of paired/pre-post designs

# metaConvert 1.0.3
- Updated the citation
- Improved the checkings for the 95% CI asymmetry

# metaConvert 1.0.2
- Fixed a bug for the online app
- Updated the vignette
- Corrected typos in the documentation

# metaConvert 1.0.1
- Added an 'auto' argument in the hierarchy of convert_df()
- Improved the way Chi-sq and Phi are converted to other measures
- Improved the compareDF function to handle different rows ordering
- Improved the aggregate_df() function to handle different time-points
- Corrected typos in the documentation

# metaConvert 1.0.0
- First version released on CRAN
