# metaConvert 2.0.0

## New: SMD sampling-variance convention (`smd_var`) — opt-in, non-breaking

- **New argument `smd_var`** on `convert_df()` and on every `es_from_*()` that builds a
  Cohen's *d* / Hedges' *g* from the large-sample formula (means/SE/CI, pooled SD, Cohen's
  *d*, Hedges' *g*, Student *t*, ANOVA *F*, eta-squared, unstandardized/standardized
  regression coefficient, point-biserial *r*, median/range→SD, digitized plots, raw mean
  difference + SD/SE/CI/p-value (`es_from_md_*`), and the ANCOVA-adjusted analogues). It selects the sampling-variance convention for the
  standardized mean difference (values are author names, matching the package's other
  method arguments):
  - `"borenstein"` (**default, unchanged behaviour**) — `v_g = cm^2 (1/n1 + 1/n2 + d^2/2N)`,
    bit-exact with `metafor::escalc(measure = "SMD", vtype = "LS2")` (Borenstein et al.,
    2009).
  - `"hedges_olkin"` (alias `"viechtbauer"`) — `v_g = 1/n1 + 1/n2 + g^2/2N`, bit-exact with
    metafor's own default `vtype = "LS"` (Hedges & Olkin, 1985; Viechtbauer, 2007). This
    unifies the endpoint/between-group family with the pre/post family (already on this
    convention) and with metafor's default.
  - The two differ by a `cm^2` (`J^2`) factor: the Hedges-Olkin variance is ~2–4% larger at
    *n* ≈ 30/arm and ~8% larger at *n* ≈ 10/arm. The default is **not** changed, so existing
    results and downstream numbers are unaffected; opt in explicitly.
  - Like `smd_to_cor`, `smd_var` can be set globally or per row via a `smd_var` column.

## New: control-SD standardizer / Glass's delta (`smd_denom`) — opt-in, non-breaking

- **New argument `smd_denom`** on the endpoint means family (`es_from_means_sd()`,
  `es_from_means_se()`, `es_from_means_ci()`) and `convert_df()`, letting the endpoint
  mean difference be standardized by a denominator other than the pooled SD:
  - `"pooled"` (**default, unchanged**) — pooled endpoint SD (Cohen's *d* / Hedges' *g*).
  - `"glass"` (alias `"control"`) — the control (non-experimental) endpoint SD (**Glass's
    delta**), with the Hedges correction and every variance term on the control degrees of
    freedom (`n_nexp - 1`). Bit-exact with `metafor::escalc(measure = "SMD1")` when the
    `vtype` is matched (`smd_var = "borenstein"` = `vtype = "LS2"`;
    `smd_var = "hedges_olkin"` = `vtype = "LS"`, which is metafor's *default* `vtype` —
    so plain `escalc(measure = "SMD1")` output matches only the `"hedges_olkin"` setting);
    honours `smd_var` for its homoscedastic variance.
  - `"glass_robust"` (alias `"control_robust"`) — Glass's delta with the
    heteroscedasticity-consistent variance, bit-exact with `metafor::escalc(measure = "SMD1H")`.
  - Like `smd_var` / `smd_to_cor`, `smd_denom` can be set globally or per row via a
    `smd_denom` column.
  - Motivation: holding the numerator (endpoint mean difference) fixed while varying the
    denominator lets a reviewer separate **estimand mismatch** (different standardizer)
    from **estimation error** when benchmarking the endpoint estimator against the
    pre/post estimators.
- Every `smd_var` × function and every `smd_denom` × variance combination is validated
  bit-exactly against metafor in `tests/testthat/test-endpoint-smd-variance.R`.

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
- Both options are pinned to the published **point estimates**: `pool_sd = FALSE`
  reproduces Morris (2008) Table 5 column *d_ppc1*, and `pool_sd = TRUE` reproduces column
  *d_ppc2* (`tests/testthat/test-EXTERNAL-MORRIS-TABLE5.R`). The **sampling variances** use
  the heteroscedasticity-robust `SMCRH`/`SMCRPH` (Bonett 2008) numerator rather than the
  homoscedastic form on the metafor-project Morris-2008 page, so they agree with that page
  exactly when `SD_pre = SD_post` within each arm and depart from it — by up to ~2.4x on
  Morris's own five studies — when they differ (this is deliberate; see Bug fixes).
- `es_from_paired_t()` / `es_from_paired_f()` **cannot** pool: a paired t identifies each
  arm's `mean_change / sd_change` ratio but not the ratio of the two arms' SDs, so the
  pooled standardizer is not recoverable from the reported statistic. These routes always
  use the per-arm construction. If you set `pool_sd = TRUE` and such rows share a pool
  with poolable rows, `summary(..., flags = TRUE)` raises an informational flag.

## Behaviour change: raw mean-difference CI from the endpoint means family

- **`es_from_means_sd()` (and its SE/CI/median-range/plot delegates) now build the
  `md_ci_lo`/`md_ci_up` bounds with the Welch–Satterthwaite degrees of freedom** instead of
  the pooled `n1 + n2 - 2`, matching `t.test(var.equal = FALSE)` bit-exactly. Only the raw
  MD confidence bounds change (d/g and every SE are untouched), and only when the arm
  variances/sizes are unbalanced enough for the Welch df to matter. The `es_from_md_*()`
  and ANCOVA-md routes keep pooled/adjusted-df CIs (they carry no per-arm SDs); the
  `summary()` A6 CI-width check accepts the whole pooled-to-Welch df band for
  package-computed md rows.

## Formula-audit fixes (July 2026)

A pre-release formula audit (multi-agent review with adversarial verification; tracker in
`AUDIT-FIX-TRACKER.md`) produced the following fixes, each paired with a test that fails on
the pre-fix code:

- **Spearman → SMD under `cor_to_smd = "cooper"`**: the d/g/logOR standard errors now carry
  the same Bonett–Wright Spearman delta correction as the r/z outputs in the same row
  (Cooper's `v_d = 4 v_r/(1-r^2)^3` is exactly linear in `v_r`); they were 5–15%
  anticonservative and internally contradictory with the corrected `r_se`. The `"mathur"`
  path is genuinely `r_se`-free and is unchanged. `es_from_spearman_rho()` also gains the
  `n_sample <- n_exp + n_nexp` fallback on direct calls.
- **Proportions**: the two surviving hard `stop()`s (prop outside [0,1]; `n_cases >
  n_sample`) are now per-row NA + warning, so one bad cell no longer aborts an entire
  `convert_df()` run (any measure) under `correct_inputs = FALSE`. Boundary (p = 0/1)
  raw/logit SEs now use the continuity-corrected denominator `n + 1`, matching
  `metafor`'s `PR`/`PLO` convention exactly (previously `sqrt((n+1)/n)` too large).
- **`es_disattenuate()`**: scalar `n_sample`/reliabilities now recycle across a vector `r`
  (previously rows 2+ silently got NA SEs); rows with an extreme corrected r
  (|r_c| > 0.999) now set the Fisher-z outputs to NA instead of emitting clamp artifacts
  (`z = atanh(0.9999)` regardless of input); the extreme-row warning text describes what
  actually happens; |r| > 1 inputs warn.
- **`compute_sem()`**: the same-sample branch now returns the exact chi-square CI
  (`df = (n-1)(k-1)`; the Wald interval covered 89.7% at n = 10, k = 2);
  `n_measurements = 1` and `icc > 1` degrade to NA with a warning (previously Inf / a
  NaN-with-SE-0 pair); a warning fires when `icc_se > 1 - icc`, the signature of
  `es_from_icc()`'s default Bonett (ln(1-ICC)-scale) SE being passed where the raw-scale
  SE is required (documented in both Rd files). `reliability_change_score()` warns on a
  negative (population-impossible) result and documents the equal-variances AND
  equal-reliabilities assumption.
- **ICC "agreement" SE scope documented + flagged**: the shared leading-order SE is exact
  for consistency ICC(3,1) but is a one-way approximation for absolute-agreement ICC(2,1)
  that assumes negligible between-rater variance — with real rater variance it is
  anti-conservative (simulated coverage ~0.74–0.76, worsening with n) and no exact fix is
  possible from summary data. The Rd/vignette now say so, and a new informational flag
  (V31) marks agreement-type rows in `summary(flags = TRUE)`.
- **Direct-call edge guards for `es_from_cronbach_alpha()` / `es_from_icc()`**: impossible
  or degenerate inputs (alpha > 1, alpha = 1 under Bonett, |ICC| ≥ 1, n too small for the
  SE denominator, fewer than 2 items/raters) now return NA instead of Inf/NaN.
- **Quality flags no longer contradict the package's own output**: the A6 CI-width check
  accepts the Welch-df band for package-computed md rows, expects `qt(.975, n_nexp - 1)`
  for Glass rows (SMD1 convention), and compares clamped raw-proportion CIs against the
  [0,1]-clipped expected width; new [UNUSUAL] flags fire when a raw-scale alpha/ICC Wald CI
  bound escapes the parameter space (with a pointer to the Bonett scale); E7 now also
  discloses the per-arm fallback when a pool contains *only* paired t/F rows under
  `pool_sd = TRUE`.
- **`convert_df(smd_denom = "glass")` scoping made explicit**: only the endpoint means
  family honours `smd_denom`; a one-time message now lists the methods that keep the
  pooled-SD standardizer, and the `@param` text states it.
- **Smaller consistency fixes**: `es_from_pt_bis_r()` now honours `smd_to_cor`; a per-row
  `smd_var` NA falls back to the default (mirroring `smd_denom`); the deprecated `measure`
  alias of `es_from_user_crude()`/`es_from_user_adj()` warns and no longer overrides an
  explicitly supplied `user_es_target_measure_*`.
- **Test infrastructure**: the psychometric suites (alpha, ICC, disattenuation, SEM/SDC,
  Spearman, proportions) previously lived only in the build-ignored `tests_save/` archive —
  R CMD check ran none of them. They now run from `tests/testthat/`, with the circular
  expectations replaced by external anchors (`metafor` `ABT`/`ARAW`/`PR`/`PLO`/`PFT`,
  hand-derived F-route ICC variance, `psychmeta` disattenuation comparisons), boundary-SE
  assertions added against `metafor`, the pooled-variance Monte Carlo graded against the
  `qt` interval the package actually emits over a q-grid (0.64/1/1.56), and a corrected
  Morris/Bonett citation split.

## Bug fixes

- **`es_from_cohen_d_adj()` now actually applies the covariate adjustment.** The exported
  function passed `n_cov_ancova` / `cov_outcome_r` to the internal effect-size builder but
  omitted `adjusted = TRUE`, so `cov_outcome_r` was silently ignored and the sampling
  variance stayed on the *crude* form (its documented formula — Cooper Table 12.3,
  `v_d = (n1+n2)/(n1 n2) (1 - r^2) + d^2/2N` with `df = n1 + n2 - 2 - n_cov_ancova` — was
  never used). At `cov_outcome_r = 0.7` the reported SE was ~40% too wide. The
  `es_from_ancova_*()` family was already correct; this aligns the standalone adjusted-*d*
  input with it. Results change **only** for rows supplying `cohen_d_adj` together with a
  non-zero `cov_outcome_r`.

- **`es_from_etasq_adj()` point estimate moved to the marginal SD scale (and its variance
  made coherent with it).** A partial eta-squared from an ANCOVA is natively defined on the
  *residual* (covariate-adjusted) SD scale; the old `d = 2*sqrt(etasq_adj/(1-etasq_adj))`
  stayed on that scale, so pairing it with the marginal-scale Cooper 12.3 variance would
  understate the SE of its own estimator by ~45% at `cov_outcome_r = 0.7` (simulated 95% CI
  coverage 0.71), and the function disagreed with `es_from_ancova_f()` fed the algebraically
  equivalent statistic by `1/sqrt(1-r^2) * sqrt(N/df)`. The adjusted eta-squared is now
  converted to the ANCOVA *F* it implies (`F = etasq_adj * df / (1 - etasq_adj)`,
  `df = n1 + n2 - 2 - n_cov_ancova`) and routed through `es_from_ancova_f()`, making the
  point estimate, variance, and CIs exactly consistent with the rest of the
  `es_from_ancova_*` family. When `cov_outcome_r` is missing the output is now NA (the
  marginal scale is unidentified) instead of a silent residual-scale value. The crude
  `es_from_etasq()` is unchanged.

- **`smd_denom` is honoured per row.** A per-row `smd_denom` column (or a length-*n* vector)
  previously collapsed to the first row's value (`match.arg(...[1])`) and applied it to every
  row; each row is now standardised by its own choice, and an unrecognised value raises a
  clear error naming the argument. `smd_denom` also joins `smd_var` / `smd_to_cor` as a
  recognised per-row column in `convert_df()`. Note: under `smd_denom = "control"` (Glass's
  delta) the derived `r`/`z`/OR inherit the control-SD standardiser, so they coincide with the
  pooled-SD conversion only when the arm SDs are equal.

- **Restored backward-compatible argument names in `es_from_user_crude()` /
  `es_from_user_adj()`.** The pre-2.0 argument names `measure` (target measure) and
  `user_es_measure_crude` / `user_es_measure_adj` (entered measure) are re-accepted as
  deprecated aliases for `user_es_target_measure_*` and `user_es_original_measure_*`. This
  fixes a reverse-dependency failure in *metaumbrella*, which still calls the previous
  names.

- **All four pooled variances rebuilt.** Each pooled two-group variance follows the LS
  *pattern* of the corresponding `metafor` pre/post measure — an empirical leading term
  plus a `g^2` term over `2N` — with `J = J(nu)` at the standardizer df and `Var(d)`
  obtained as `Var(g)/J^2` (an exact identity, since `g = J*d` with `J` a deterministic
  constant). Note this is a pattern, **not** a single closed-form rule: metafor puts `n`
  (not the standardizer df) in the `g^2` denominator of its non-heteroscedastic LS forms
  while evaluating `J` at the standardizer df — the two df deliberately differ; and its
  heteroscedastic variants (`SMCRH`/`SMCRPH`) use `n - 1` in both terms. Three defects
  are fixed:
  - *Spurious leading `J^2`* on `morris_dz`/`morris_drm`/`morris_dav` (the `LS2`
    convention). It understated the sampling variance by `1 - J^2` — 9% at n = 10/arm,
    3% at n = 30/arm — which over-weighted small studies in an inverse-variance
    meta-analysis. Pooled `morris_dz` is now **bit-exact** with
    `metafor::escalc(measure = "SMD")` on change scores.
  - *Homoscedasticity-fragile numerator* on `bonett` and `morris_dav`. Both used the
    identity `Var(change) = 2*sigma^2*(1 - r)`, valid only when `SD_pre = SD_post`. The
    numerator variance is now taken from the **empirical pooled change SD** (metafor's
    heteroscedasticity-robust `SMCRH`/`SMCRPH` treatment, = Bonett 2008 eq. 10/19). At
    `SD_pre/SD_post = 0.64`, `bonett`'s variance was **43% too small** and its 95% CI
    covered only **86%** of the time; with the df-weighted `r_avg` it now reduces
    *exactly* to Viechtbauer's published two-group variance
    `2(1 - r)(1/nT + 1/nC) + g^2/(2N)` under homoscedasticity (even when `n1 != n2`) and
    covers at the nominal rate in every regime tested **given a correctly specified
    `r_pre_post`** (Monte Carlo, `tests/testthat/test-pooled-variance-calibration.R`;
    when the default `r_pre_post = 0.8` is imputed and the true correlation differs,
    coverage degrades regardless of the variance formula — see `?convert_df`).
  - *Homoscedastic `morris_dav` variance.* Its `g^2` term used the coefficient
    `(1 + r^2)/(4N)` — metafor's homoscedastic `SMCRP`, the `SD_pre = SD_post` special
    case. It now uses the two-group fourth-moment form of **Bonett (2008) eq. 19**
    (metafor `SMCRPH` carried across two independent arms), so both its terms are
    heteroscedasticity-robust. `nu = 2m/(1 + r^2)` (Cousineau, 2020, eq. 2) is retained
    but used **only** for the Hedges bias correction `J`, not (as the previous NEWS
    claimed) as `1/(2*nu)` in the `g^2` coefficient — that identity was false
    (`1/(2*nu) = (1 + r^2)/(4m)`, not `/(4N)`).

  A directly published two-group variance exists for pooled `d_z` (it *is*
  `metafor::escalc(measure = "SMD", vtype = "LS")` on the change scores) and for `d_av`
  (Bonett 2008 eq. 19); `d_rm` follows from `d_z` by the exact rescaling
  `d_rm = d_z*sqrt(2(1-r))` (Caldwell & Vigotsky 2020). The robust-`bonett` departure from
  Viechtbauer's homoscedastic form under heteroscedasticity is Monte-Carlo calibrated in
  the test suite.

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
  paired data, and the imputation of an unreported `r_pre_post` for rows that consume it.
  The imputation note is method-aware: for change-SD standardizers (`morris_dz`, `cooper`/
  `morris_drm`) the assumed correlation scales the **point estimate**; for the default
  `bonett` and for `morris_dav` the point estimate is r-free and only the SE/CI move. **The
  SE and CI are unreliable when `r_pre_post` is defaulted rather than reported** — at a true
  `r` of 0.3 with the default 0.8, the reported pre/post variance is ~68% too small (95% CI
  coverage ~0.75), a larger error than the heteroscedasticity the robust numerator fixes.
  Supply `r_pre_post_exp`/`r_pre_post_nexp`, or run a sensitivity analysis over plausible
  values.

## Documentation

- Corrected the `@details` of the four single-group mean-change converters, which
  documented a sign-flipped argument mapping (`mean_pre_exp = mean_change_exp`,
  `mean_post = 0`) that was the reverse of the code.

## Tests

- The pre/post and change-score test suite (21 files) is now wired into `tests/testthat/`
  and runs under `R CMD check`; previously only `test-ANCOVA-MD.R` ran. Long-running files
  are gated on `NOT_CRAN`.

## Summary of major additions

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
