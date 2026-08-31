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

- Added the `flag_group` option to group across multiple rows the checks
performed by metaDETECT
- Added and refined metaDETECT quality flags for reliability pools, mixed transforms, and 2x2 table margins
- Added `es_formulas()` to evaluate impact of choice of formula on effect size estimates 
- Improved speed performance of the convert_df() 
- Improved the reliability paths (alpha/ICC) 
- Corrected formulas and conversion explanations in the documentation

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
