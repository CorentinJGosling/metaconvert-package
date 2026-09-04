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
- **A correlation entered through `user_es_*` with a confidence interval**: the
  precision is now read on the Fisher-z scale, not the r scale. A correlation CI is
  almost always `tanh(z +/- 1.96/sqrt(n - 3))` -- asymmetric on the r scale -- and taking
  its r-scale half-width as the SE and delta-mapping that to z did not return
  `1/sqrt(n - 3)`: at n = 30 the Fisher's z standard error came out 0.966x (r = .3) to
  1.067x (r = .9) of the correct value, inverse-variance weights 7% too large to 12%
  too small, and the z interval was shifted. Transforming the bounds first reproduces
  `1/sqrt(n - 3)` exactly; the r-scale SE is its delta-method image, and a CI-only
  point estimate is the centre on the z scale rather than the arithmetic midpoint of
  the r bounds. A user-supplied SE is untouched.
- **Incidence rate ratio with zero events in either arm**: `es_from_cases_time()`
  returned `logirr = -Inf` (or `+Inf`) with `se = Inf`, so the row reached `summary()`
  as an estimate of 0 with an infinite standard error and was lost. Both counts now
  receive +0.5 when either is zero, the `metafor::escalc(measure = "IRR")` default and
  the rate analogue of the correction the 2x2 route already applies: 0 versus 14 events
  over equal person-time now returns `-3.367 (se 1.438)`, matching metafor. The rate
  difference and the person-time NNT stay on the raw counts.
- **`es_from_pearson_r()` / `es_from_fisher_z()` / `es_from_spearman_rho()` with arm
  sizes**: under `cor_to_smd = "viechtbauer"` (the default) the r -> d map inverted the
  biserial correlation at a 50/50 split even when `n_exp` and `n_nexp` were supplied
  (they reached only the standard error and Hedges' J). A d of 0.5 converted to r at
  20/80 and back came out 0.467 (-6.6%), and 0.421 (-15.9%) at 10/90. Supplied arm sizes
  are now handed to `metafor::transf.rtod()`, whose inversion is then exact; rows with
  `n_sample` alone are unchanged. The same applies to a `user_es_*` correlation.
- **Pre/post and paired standard errors now follow `smd_var` on every branch.** The
  `morris_drm` branch of the single-group kernel (the default for mean-change and paired
  t/F data) applied Hedges' J-squared to its whole variance (Borenstein) while the
  `morris_dz`, `bonett` and `morris_dav` branches used metafor's form, so at
  `r_pre_post = 0.5` the two raw-score branches returned the same estimate with standard
  errors differing by 19% at n = 5 and 3% at n = 25, and the documented identity
  `Var(d_rm) = 2(1 - r) Var(d_z)` did not hold. All four branches, single-group, two-group
  and pooled, and all 17 pre/post and paired routes now read `smd_var`, as the two-group
  routes already did: under the default `"borenstein"` every branch is
  `J^2 (L + d^2 C)`, so the `morris_dz`/`bonett`/`morris_dav` standard errors shrink by J
  (about 3% at n = 10, under 1% at n = 50); under `"hedges_olkin"` every branch is
  bit-exact with `metafor::escalc(measure = "SMCC"/"SMCR"/"SMCRH"/"SMCRPH")`.
- **Confidence intervals of correlations** are the back-transformed Fisher-z interval
  on every route (`tanh(z +/- 1.96 z_se)`, the construction of `cor.test()`, metafor
  and CMA). `es_from_pearson_r()`, `es_from_fisher_z()`, `es_from_spearman_rho()`, a
  `user_es_*` correlation and the `smd_to_cor = "lipsey_cooper"` branch returned
  `r +/- t se` on the r scale, which leaves [-1, 1] at small n: `r = -0.226, n = 5`
  gave `[-1.74, 1.28]`, `r = 0.85, n = 8` gave `[0.59, 1.11]`. The other correlation
  routes already back-transformed. The `[INFO]` "Wald interval escapes the parameter
  space" flag now fires only on a user-supplied interval, which is handed back as
  entered.
- **Odds ratio to risk ratio: the reconstructed 2x2 table is no longer rounded to whole
  counts unconditionally.** A reported OR of 2.00 on margins 40/60 and 30/70 admits no
  integer table; rounding replaced it with 16/24/14/46, whose odds ratio is 2.19, and
  the risk ratio inherited the error (log RR 0.539 against 0.477). A whole-count table
  is now taken only when its odds ratio rounds back to the reported one at the reported
  precision -- so a table printed to 2 or 3 decimals is still recovered exactly (97% and
  100% of the time), a continuity-corrected table still comes back with its
  half-integer cells, and an adjusted or otherwise non-integer OR keeps the exact root.
- **Risk difference from a 2x2 table with an empty cell** takes the same +0.5 correction
  as the odds ratio and the risk ratio, matching `metafor::escalc(measure = "RD")`
  (`add = 1/2, to = "only0"`). On the raw counts the empty arm contributed exactly 0 to
  the variance (`0/59` vs `8/48` returned se 0.0468 against metafor's 0.0490), and a
  double-zero table returned `NA` where metafor returns a finite estimate. Tables
  without an empty cell are unchanged.
- **Degenerate p-values on every p-value route**, not only `or_pval`: `p <= 0` (the
  numeric transcription of "p < .001") inverted to an infinite statistic -- `d = Inf`
  with `se = Inf` on the t/F/point-biserial/paired-t/chi-square routes -- and `p >= 1`
  to an infinite standard error on the `rr`, `rd`, `md`, `ancova_md`, `mean_change`
  and `linreg_b` routes. Both now return `NA`. A `p` of exactly 1 is kept where it is
  the ordinary boundary `t = 0` / `chi-square = 0` (a null estimate with a finite
  standard error). The guard lives in the routes, so direct callers get the same
  behaviour as `convert_df()`.

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
- New input column `alpha_type`, recording whether a reported Cronbach's alpha is the
  raw (covariance-matrix) or standardised (correlation-matrix) coefficient, with a new
  `[INFO]` flag (V45) when a pool mixes the two. They are different coefficients, equal
  only when the item variances are, and `alpha_type` enters neither the estimate nor the
  (n, k) standard error -- so nothing numeric could reveal the mix. Accepts `"raw"` /
  `"std.alpha"` and similar as synonyms for `"covariance"` / `"correlation"`; the stored
  values name the computation because `alpha_to_es = "raw"` already means the
  untransformed Bonett scale, an unrelated sense of the same word. Provenance only: it
  changes no existing result, and a blank is not treated as a level.
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
