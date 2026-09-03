# metaConvert 2.0.1

## Results change: re-run analyses produced with 2.0.0

This release fixes conversion formulas, so a dataset re-run under 2.0.1 can return
different effect sizes, different standard errors, or a different estimation route than
it did under 2.0.0. **If you have results from 2.0.0 that are not yet published, re-run
them.** Nothing here requires a change to your code or your data.

What moves, and roughly by how much:

- **Risk ratio from an odds ratio** (`or_to_rr = "metaumbrella_cases"`/`"_exp"`, the
  default): the reported log-OR standard error was ignored whenever all four margins were
  known, so the RR carried the precision of the reconstructed table instead of the
  study's own. A reported SE of 0.05 and one of 1.5 returned the same RR standard error.
  It is now delta-propagated. Non-significant odds ratios could previously be returned as
  significant risk ratios.
- **Odds ratio route ranking**: an OR reported with a p-value now uses the p-value
  (an exact standard error) rather than one imputed from the margins, which was about
  2.2x too wide.
- **Adjusted SMD route ranking**: an adjusted mean difference with its residual SD now
  outranks adjusted means with standard errors, which attenuate (measured ratio 0.94).
- **Correlations from a 2x2 table with a zero cell**: the +0.5 continuity correction is no
  longer applied to the tetrachoric route, which is documented as OR/RR only. Such a row
  now matches `metafor`'s boundary estimate and its large variance instead of an interior
  estimate whose standard error was 17x to 27,000x too small.
- **Fisher's z from a large standardised mean difference** (`smd_to_cor = "viechtbauer"`,
  the default): above roughly |d| = 2.6 the biserial correlation leaves [-1, 1] and z used
  to saturate at a constant that depends only on the group-size split, so two unrelated
  studies could report the same value. Such rows now return `NA` with a message naming
  `smd_to_cor = "lipsey_cooper"`, which is bounded and keeps every study.
- **`cor_to_smd = "mathur"`**: the documented default `unit_type = "raw_scale"` was never
  applied, so omitting the argument silently returned NA effect sizes. Seven entry points.
- **Correlations and z from phi or chi-square when `mvtnorm` is absent**: d and g were
  overwritten with values from the phi fallback even though they had been computed
  correctly, and the message named a false cause. d moved by up to 44%.
- **Pre-post `morris_dz`**: returned NA when `|r_pre_post| >= 1` although d_z does not
  depend on r. `es_from_mean_change_sd()` with r = 1 now also returns a finite d_z.
- **Per-row `pre_post_to_smd`**: the column was ignored by the mean-change and paired t/F
  routes, which silently used a different standardiser.
- **Standard errors and confidence bounds** are now `NA` rather than a wrong finite number
  in several degenerate cases: a negative `sd_dv` in the standardised-regression route, a
  non-positive `or`/`rr`, `baseline_risk = 0`, an impossible 2x2 cell derived from
  proportions, and a risk-difference confidence limit of exactly 0 (which returned Inf).
- **Risk-ratio standard error converted from an odds ratio** (`or_to_rr =
  "metaumbrella_cases"`/`"_exp"`, the default): the reported log-OR standard error was
  carried onto the log RR by differentiating through the table reconstruction with all
  four margins held *fixed*. That is the conditional standard error, and the case margin
  is an observed statistic rather than a design constant -- unlike for the odds ratio, it
  is not ancillary for the risk ratio, so conditioning on it discards real variability.
  The multiplier is now the ratio of the two *marginal* standard errors at the
  reconstructed table (Katz over Woolf), which is the product-binomial model
  `es_from_2x2()`, `metafor::escalc(measure = "RR")` and the Cochrane Handbook all use.
  On a 90/10 versus 60/40 table the old form returned 0.0711 against 0.0882: a 1.54x
  inflation of the study's inverse-variance weight, and the same study disagreeing with
  itself when entered as a 2x2 table instead. The two forms agree to within a fraction of
  a percent at ordinary 2-25% event rates and separate as events become common.
- **Hedges' g standard error from an odds ratio reported with one standard error**:
  `es_from_or_se()`/`es_from_or_ci()` pass a single standard error for a whole vector of
  d, and the variance selection is an `ifelse()`, which returns a result the length of
  its *test*. `g_se` was silently `NA` for every row after the first.
- **Cronbach's alpha and McDonald's omega with a self-contradictory interval**: a
  confidence interval that does not bracket its own point estimate no longer yields a
  standard error. `es_from_icc()` already refused it, so the same input previously
  produced a standard error from two routes and `NA` from the third.
- **A ratio effect size entered through `user_es_*` with a non-positive confidence
  bound**: `log(0)` is not finite, so no log-scale standard error exists -- but the
  natural-scale half-width was being exported in its place, and nothing downstream told
  the two apart. `or = 2` with a CI of `[0, 4.5]` was returned as `se = 1.148` with a CI
  of `[0.211, 18.98]`, and no quality flag fired on it. The standard error is now `NA`,
  leaving the row visible and unpooled.
- **Within-group `gw` through `user_es_*` with a sample size supplied and `df <= 1`**:
  Hedges' J is undefined there and the route warns that the row yields `NA`; it was
  returning an *uncorrected* `gw` instead.
- **Person-time risk difference and NNT at `baseline_rate = 0`**: the derived standard
  error collapses to exactly 0 -- an infinite inverse-variance weight. Both are now `NA`.
  A genuinely null rate difference with a usable standard error is unaffected.
- **`r_pre_post` rejected by input validation** is re-filled with the value passed to
  `convert_df()` rather than falling through to each route's own default of 0.8, and is
  treated as defaulted by the V6 note and the r-sensitivity annotations. Measured on one
  row: `g = 0.0727` under the route default against `0.2887` for the requested `r = 0.3`.
- **`or_pval` of exactly 1** joins 0 as unusable: `qnorm(1/2)` is 0, so the recovered
  standard error is infinite and the interval unbounded.

Diagnostic checks that now behave differently (metaDETECT):

- A6 (CI width versus standard error) now runs for `or`/`rr`/`irr`/`hr`, so new flags will
  appear where none did; it no longer false-fires on package-computed correlations.
- B3/B5/B6 no longer fire on person-time NNT rows; B5 uses the correct harm-side bound.
- D1 gained minimum-deviation floors for `dw`/`gw`/`rp`/`zp`; D2's risk-ratio
  normalisation and D3's spread scale were re-derived; V29's ceiling is now sample-size
  aware.
- V26 is suppressed where `r_pre_post` is supplied. Its previous claim to need no assumed
  pre-post correlation was incorrect and is retracted: it false-fires on correctly
  computed paired intervals when that correlation is below roughly 0.15 to 0.35.
- E1's dispersion statistic for the exponentiated ratio measures is computed exactly on
  the log scale instead of through a natural-to-log approximation. The approximation is
  attained at a different element from the exact statistic in about a fifth of rows, so
  no scalar divisor repairs it; it flipped roughly 11% of E1 decisions.
- V29 reports the sample-size-inflated ratios as `[INVALID]` "arithmetically impossible"
  rather than `[UNUSUAL]` "implausible" when they fall outside the assumption-free
  bounds -- a strictly stronger statement about the same rows.
- `flag_group` and `enable_cross_row` now reach the Tier-1 cross-row checks (V35-V43) as
  well as the Tier-2 ones. Passing either to `summary()` alone therefore retunes only
  half of them, and `summary()` now says so.
- Invalid `flag_options` values fall back to the documented default. A character
  threshold previously made every comparison against it `NA`, silently disabling the
  check the user was trying to tighten, while `NA`, `NULL` or a longer vector aborted the
  run inside a comparison naming neither the option nor the value.
- The six per-arm direction columns of the paired t and F routes
  (`reverse_paired_t_pval_exp`, `reverse_paired_f_exp` and their siblings) are registered
  inputs and reachable from `convert_df()`. Those routes' own documentation named them as
  the remedy for two arms moving in opposite directions, but they could previously be
  reached only by calling the route directly.
- New opt-in cross-row check `se_outlier_missing_n` (D2c), off by default. The
  sample-size-normalised SE-outlier checks drop a row that reports a standard error but
  no sample size, silently -- nothing separated "checked and found fine" from "never
  checked", and no threshold reached such a row. When enabled, a dropped row whose raw
  standard error is more than `se_missing_n_ratio`-fold (default 5) tighter than the
  pool median is flagged. Confined to rows the normalised checks could not see, so it
  can add a flag but never change one they raised.
- `es_guidance` no longer names columns that cannot produce an estimate: `es_from_or()`
  reconstructs from the case margin rather than the arm sizes, an NNT or risk difference
  from a ratio also needs `baseline_risk`, and `dw`/`gw`/`mdw`/`prop`/`alpha`/`omega`/
  `icc` have no adjusted scope at all.

- Added the `flag_group` option to group across multiple rows the checks
performed by metaDETECT
- Added and refined metaDETECT quality flags for reliability pools, mixed transforms, and 2x2 table margins
- Added `es_formulas()` to evaluate impact of choice of formula on effect size estimates 
- Improved speed performance of the convert_df() 
- Improved the reliability paths (alpha/ICC) 
- Corrected formulas and conversion explanations in the documentation
- `compare_df()` validates its `output` argument, and refuses a non-unique
  `ordering_columns` key rather than silently dropping the duplicated rows behind a
  base-R replacement warning.
- `summary()` resets row names, so the positions cited in flag messages ("row 52") match
  the rows as printed.
- Documentation corrections where the text described code the package no longer runs --
  chiefly the scope of the +0.5 continuity correction in `es_from_2x2()`, the risk-ratio
  standard error in `es_from_or_se()`, and the magnitude of the bias the paired t/F
  routes carry when the two arms move in opposite directions (biased *by* `2|d_nexp|` is
  not the same as *halved*: the estimate collapses to exactly zero when the two arms'
  magnitudes are equal).

# metaConvert 2.0.0
- Implemented the metaDETECT framework of automated checks
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
