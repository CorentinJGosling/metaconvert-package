# metaConvert 2.1.0

## Behaviour of the two-group pre/post SMD (no default change)

- **`pool_sd` still defaults to `FALSE`.** The default two-group construction is
  unchanged: each arm's change is standardized by that arm's **own** SD and the two
  within-group values are subtracted, their sampling variances adding because the arms
  are independent. This is Morris's (2008) **d_ppc1** (from Becker, 1988), and it is what
  `metafor` users do for this design (see the metafor-project Morris 2008 page: compute
  `escalc(measure = "SMCR")` per arm, then subtract the estimates and add the variances).
  Viechtbauer describes it as the more broadly applicable of the two options because it
  does not assume the arms' true standardizing SDs are equal.
- **`pool_sd = TRUE` is now correct and fully supported** as an opt-in. It divides the
  difference in mean change by a **single SD pooled across arms** — Morris's (2008)
  **d_ppc2** (his eq. 8-9), which he recommends as more efficient, at the cost of
  assuming the arms' true SDs are equal. The two options target the same estimand when
  that assumption holds (as randomization implies at baseline) and different estimands
  when it does not. **This is a deliberate analytic choice, not a technical detail**, so
  the package does not make it for you.
- Both options are now pinned to the published values: `pool_sd = FALSE` reproduces
  Morris (2008) Table 5 column *d_ppc1*, and `pool_sd = TRUE` reproduces column *d_ppc2*
  (`tests/testthat/test-EXTERNAL-MORRIS-TABLE5.R`).
- `es_from_paired_t()` / `es_from_paired_f()` **cannot** pool: a paired t identifies each
  arm's `mean_change / sd_change` ratio but not the ratio of the two arms' SDs, so the
  pooled standardizer is not recoverable from the reported statistic. These routes always
  use the per-arm construction. If you set `pool_sd = TRUE` and such rows share a pool
  with poolable rows, `summary(..., flags = TRUE)` raises an informational flag.

## Bug fixes

- **All four pooled variances rebuilt on one rule.** The pooled two-group variances now
  follow the rule metafor's own `escalc()` encodes for every pre/post measure
  (`SMD`, `SMCC`, `SMCR`, `SMCRH`, `SMCRP`):

      Var(g) = Var(numerator)/SD_standardizer^2 + g^2/(2*nu_std),   J = J(nu_std)

  with **no leading `J^2`** and `Var(d)` obtained by substituting `d` for `g` (not by
  dividing `Var(g)` by `J^2`). Three defects are fixed by this:
  - *Spurious leading `J^2`* on `morris_dz`/`morris_drm`/`morris_dav` (the `LS2`
    convention). It understated the sampling variance by `1 - J^2` — 9% at n = 10/arm,
    3% at n = 30/arm — which over-weighted small studies in an inverse-variance
    meta-analysis. Pooled `morris_dz` is now **bit-exact** with
    `metafor::escalc(measure = "SMD")` on change scores.
  - *Homoscedasticity-fragile numerator* on `bonett` and `morris_dav`. Both used the
    identity `Var(change) = 2*sigma^2*(1 - r)`, valid only when `SD_pre = SD_post`. The
    numerator variance is now taken from the **empirical pooled change SD** (metafor's
    heteroscedasticity-robust `SMCRH`/`SMCRPH` treatment, and already what this
    package's single-group `bonett` did). This is the largest of the three: at
    `SD_pre/SD_post = 0.64`, `bonett`'s variance was **43% too small** and its 95% CI
    covered only **86%** of the time; it now covers at nominal rate in every regime
    tested (Monte Carlo, `tests/testthat/test-pooled-variance-calibration.R`).
  - *Wrong degrees of freedom for `morris_dav`.* Its standardizer is the quadratic mean
    of two correlated SDs, so its effective df is `nu = 2m/(1 + r^2)` (Cousineau, 2020,
    eq. 2), not `m`. `J` and the CI are now evaluated at `nu` — as metafor's `SMCRP` and
    this package's own single-group `dav` branch already did. The `(1 + r^2)/4`
    coefficient in the variance is *correct and sourced*; it is exactly `1/(2*nu)`.

  No published sampling variance exists for the pooled two-group form of `d_z`, `d_rm`
  or `d_av`; these are delta-method generalizations of the single-group results, and are
  Monte-Carlo calibrated in the test suite.

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
