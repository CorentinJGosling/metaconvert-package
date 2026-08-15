# metaConvert (development version)

## Bug fixes (silent checks that never fired)

Both bugs below share one mechanism: a hardcoded list of column names containing
names that are not columns. Because every consumer filters with
`intersect(., colnames(x))` or `%in% colnames(x)`, a wrong name is dropped with no
error and no warning — the check simply never runs. A package-wide sweep of the ten
files carrying such lists found exactly these two; the remaining suspects
(`.positive_columns()`'s `logirr_se`, `.ci_triplets()`'s `logirr`/`loghr` triplets)
name output-only quantities with no input analogue and break nothing.

- **`data_extraction_sheet()` emitted `all_info_expected` where the pipeline reads
  `info_expected`.** The string appeared at exactly one line in the package and was
  read by nothing. A user who filled in the sheet exactly as generated had that
  value silently discarded: `.check_data()` manufactured `info_expected` as all-NA
  and `summary()` dropped it, so the documented `info_expected` output column never
  populated for anyone using the generated sheet. **If you keep extraction files
  created with an earlier version, rename that column to `info_expected`** — the old
  name is not accepted, and as before it will be ignored without a message.

- **Flag V6 (`default r_pre_post used`) never fired on paired-*t* / paired-*F*
  rows.** The eligibility test used its own column list in which 7 of 11 names did
  not exist (`mean_post_exp`, `mean_post_nexp`, three `*_single_group` names, and
  unsuffixed `paired_t` / `paired_f`). Pre/post and mean-change rows were unaffected
  — the list's four real names still caught them — but rows whose pre/post
  information is a paired test statistic raised no V6 flag, did not appear in
  `attr(res, "r_defaulted")`, and did not receive the `(r-sensitive: ...)`
  annotations on A6/E2/E2b in `summary()`. These are the rows where the assumption
  bites hardest: paired-*t*/*F* routes standardise by the change SD, so the assumed
  `r_pre_post` scales the **point estimate**, not merely the standard error.
  Affected results are unchanged in value; what was missing is the warning that they
  depend on an imputed correlation. The verbose note used a separate, correct-but-
  incomplete list (it omitted `mean_pre_se_*`, `mean_pre_ci_*` and
  `mean_change_ci_*`); both call sites now share one list.

## Internal

- New internal helpers `.r_consuming_columns()` / `.rows_with_r_consuming_data()`
  give `convert_df()` a single source of truth for "does this row feed a route that
  consumes `r_pre_post`?".
- New regression guards: `tests/testthat/test-extraction-sheet-column-names.R`
  asserts every name `data_extraction_sheet()` emits, for every `measure`, is a name
  `.check_data()` recognises; `tests/testthat/test-V6-r-defaulted-columns.R` covers
  12 r-consuming data shapes plus negative controls. Both assert the column lists
  contain no phantom names, which is the check that would have caught either bug.


# metaConvert 2.0.1

## Results change: re-run analyses produced with 2.0.0

This release fixes bugs in the conversion formulas themselves, so a dataset re-run
under 2.0.1 can return different effect sizes, different standard errors, or a
different estimation route than it did under 2.0.0. **If you have results from 2.0.0
that are not yet published, re-run them.** Nothing here requires a change to your
code or your data.

Six families of results move. In rough order of how much:

| what moves | how much | who is affected |
|---|---|---|
| Pearson *r* and Fisher's *z* obtained from an **odds ratio** | `\|r\|` always larger, by ~10-50% | anyone converting ORs to correlations |
| `r_se` / `z_se` on any design that is not a simple two-group comparison | ANCOVA up to -49%, regression -22%, pre-post -8%, Glass +1.5% | correlation meta-analyses mixing designs |
| *r* and *z* **point estimates** on covariate-adjusted rows | up to ~10% | ANCOVA rows with covariates |
| *d*, *g*, OR, *r*, *z* for any study entered as an **eta-squared** | -11% to +67%: slightly *down* with equal arms, sharply *up* the more unequal they are (+66% at a 9:1 split) | anyone using the `etasq` column |
| which route `summary()` picks for **adjusted** rows | changes the estimate where the study reports several | ANCOVA rows reporting both an SD and a test statistic |
| `logrr_se` when `or_to_rr = "dipietrantonj"` was used on any row | standard errors were taken from the wrong row | anyone who set that option |
| *r* / *z* confidence bounds on `reverse_2x2 = TRUE` rows | the interval came back **inverted** and is now the right way round | anyone flipping the direction of a 2x2 table |

Separately, inputs that are arithmetically impossible (a negative standard deviation,
a transposed confidence interval, `\|cov_outcome_r\| = 1`) previously produced a
confidently wrong number -- sometimes with the **sign reversed**. They now return `NA`
and are flagged. Valid inputs are unaffected.

Effect sizes that do **not** change: mean differences, and *d*/*g*/OR/RR/NNT computed
from means and SDs, from a *t* or *F* statistic, or from a 2x2 table.

Details for each item are in the sections below.

## Bug fixes (d -> r / z conversion)

Both fixes live in the internal `.smd_to_cor()` and affect **every** route whose
sampling variance differs from the crude two-group design, not only the ANCOVA
family. Rows whose `d_se` is the default large-sample formula are bit-identical
to before, and still bit-identical to `metafor::escalc(measure = "RBIS")`; rows
that carry their own standard error are not. The `r_se` / `z_se` movers, measured
against the released 2.0.0:

| route | `r_se` / `z_se` |
|---|---|
| every ANCOVA route | -2% at `cov_outcome_r = 0.3`, -9% at 0.5, -21% at 0.7, **-49% at 0.9** |
| `es_from_linreg_b_se()` / `_ci()` / `_pval()` / `_t()` | **-22%** |
| the paired / pre-post family (`es_from_paired_t()`, `es_from_mean_change_sd()`, `es_from_means_sd_pre_post()`, ...) | -8% to -2% |
| `smd_denom = "control"` (Glass's delta) | +1.5% |
| user-supplied standard errors | depends on the value supplied |

Covered by `tests/testthat/test-smd-to-cor-variance.R`.

- **`r_se` and `z_se` ignored the precision of the effect size they came from.** Under
  the default `smd_to_cor = "viechtbauer"` the correlation variance is Soper's closed
  form, a function of `n`, `p` and `r` only, so it reported the crude study's precision
  regardless of the actual design -- `z_se` was literally frozen at `1/sqrt(N - 1)`. An
  ANCOVA row whose `d_se` fell 54% as `cov_outcome_r` rose from 0 to 0.9 returned an
  unchanged `r_se`, so the same study received a covariate-aware weight in a Hedges'
  *g* meta-analysis and a covariate-blind one in a correlation meta-analysis. The
  closed form is now rescaled by `vd / vd_crude`, which is exactly 1 whenever `d_se`
  is the default large-sample formula. Rescaling the closed form, rather than
  replacing it with a plain delta propagation, is what makes an adjusted row
  *design-transparent*: at the same marginal *d* it now behaves exactly like the crude
  row it is pooled with, and agreement with `metafor` is preserved. Monte Carlo 95% CI
  coverage for an ANCOVA correlation improves from **0.993 to 0.948** (SE/empirical-SD
  from 1.38 to 1.00), and at `cov_outcome_r = 0.89` from 1.000 to 0.945 (2.07 to 0.98);
  Fisher's *z* behaves likewise. `smd_to_cor = "lipsey_cooper"` already propagated `vd`
  and is unchanged.

  Not changed by this release, but worth stating so it is not misattributed to it:
  Soper's closed form is itself optimistic at large `|r|`. At a marginal *d* of about
  2 the 95% interval covers ~0.90 rather than 0.95, in the crude two-group case as
  much as in the adjusted one. That is inherited from the estimator `metafor` uses for
  `measure = "RBIS"` and is unaffected by the rescaling.
- **The correlation point estimate depended on the source study's covariate count.**
  `h` was built from the ANCOVA residual df `N - 2 - n_cov_ancova`, but the *d* reaching
  the conversion is always on the marginal (unadjusted) SD scale, per Cooper's
  eq. 12.23/12.24 convention. The result was neither the marginal point-biserial *r*
  (which needs `N - 2`) nor the partial one (which would additionally need
  `h * (1 - cov_outcome_r^2)`): an identical *d* returned `r = 0.4236` at
  `n_cov_ancova = 0` but `r = 0.4506` at `n_cov_ancova = 5`, a 6.4% swing driven purely
  by how many covariates the source happened to adjust for. `h` now uses `N - 2`, so an
  adjusted and an unadjusted study reporting the same marginal *d* return the same
  correlation. `n_cov_ancova` now enters only the confidence-interval degrees of
  freedom, where it belongs -- including the `lipsey_cooper` *r* interval, which
  previously ignored it while the *d*/*g* interval on the same row did not.

## Bug fixes (ANCOVA / eta-squared family)

These change returned values. All were found by auditing the family against Cooper,
Hedges & Valentine (2019) ch. 12 (eqs. 12.23-12.26, tables 12.1 and 12.3) and
Lai & Kelley (2012), and are covered by `tests/testthat/test-ancova-route-equivalence.R`.

- **`es_from_etasq()` was wrong for unequal group sizes.** It used the closed form
  `d = 2 * sqrt(etasq / (1 - etasq))`, which is the equal-*n*, large-*n* limit: the
  constant `2` stands in for `sqrt((n_exp + n_nexp - 2) * (1/n_exp + 1/n_nexp))`, so
  the returned *d* did not depend on how the sample split between arms. It now inverts
  the exact identity `etasq = F / (F + N - 2)` and applies Cooper table 12.1, and so
  reproduces the Cohen's *d* computed from the raw data to machine precision at any
  split. **The correction goes in both directions.** With unequal arms the old value
  was too small -- 39% too small at `n_exp/n_nexp = 10/90`, 12% at 20/60 -- but with
  equal arms it was too *large* by `sqrt(n/(n-1))`, so a balanced study's *d* now moves
  slightly **down**: 1.0% at 50 per arm, 0.25% at 200, 12% at 5. Since most designs are
  balanced, expect a small downward shift far more often than a large upward one.
  `es_from_etasq()` now agrees exactly with `es_from_anova_f()` and
  `es_from_pt_bis_r()` given the same study; `es_from_etasq_adj()` already went
  through *F* and is unchanged. This route feeds *d*, *g*, OR, *r* and *z* alike, so a
  study estimated from an eta-squared changes on every one of them.
- **The adjusted hierarchy now prefers `es_from_ancova_md_sd()` to the ANCOVA test
  statistics.** Given an adjusted mean difference with its residual SD, that route's
  point estimate is unaffected by covariate imbalance, whereas `ancova_t` / `ancova_f`
  / `ancova_t_pval` / `ancova_f_pval` / `etasq_adj` all attenuate by
  `1/sqrt(1 + D/(1/n_exp + 1/n_nexp))` (Lai & Kelley eq. 5). It was previously ranked
  below all of them, so `summary()` selected the attenuated estimate when both were
  available. `es_from_ancova_md_se()` / `_md_ci()` / `_md_pval()` deliberately stay
  below the test statistics: they must invert a reported SE, CI or *p*-value and are
  in the same attenuation tier.

  The adjusted **means** hierarchy had the same defect and is reordered the same way:
  `es_from_ancova_means_sd_pooled_adj()` is handed the adjusted pooled SD directly, so
  its point estimate is imbalance-immune, but it was ranked *below*
  `es_from_ancova_means_se()` / `_means_ci()`, which must recover the standardizer from
  a reported SE or CI and attenuate. A study reporting both now returns the
  unattenuated estimate (on one imbalanced example, *g* 0.542 to 0.572).
- **An omitted `cov_outcome_r` is no longer silently replaced by 0.5.**
  `es_from_cohen_d_adj()` is the only adjusted route where `cov_outcome_r` enters the
  variance but not the point estimate, so R never forced the argument and an internal
  fallback of `R = 0.5` applied, understating the standard error by ~12% and inflating
  the study's inverse-variance weight by ~29%. Omitting it now blanks every
  variance-bearing output (`d_se`, `g`, `r`, `z`, `logor` and their CIs), leaving only
  the point estimate `d`, which does not depend on `R`. This affects **direct calls
  only**: `convert_df()` already supplied an all-NA column and returned `NA`, so the
  pipeline is unchanged. Crude (non-adjusted) routes are unaffected, and an omitted
  `n_cov_ancova` still defaults to `q = 0` (it only costs degrees of freedom).
- **A transposed confidence interval or a non-positive dispersion no longer reverses the
  sign of the effect size, nor produces a negative standard error.** An interval width
  taken as `ci_up - ci_lo` is negative when the bounds are transposed, and a standard
  deviation or standard error can be entered negative. Where the route *divides* by that
  quantity the result was a **sign-flipped** effect size while the mean difference kept
  its original sign, leaving the two pointing opposite ways (`md = +4` returning
  `d = -2`); where it only *propagates* it, a negative standard error reached the output
  and became a negative sampling variance, i.e. a wrong-signed inverse-variance weight.
  Widths are now taken as absolute, so a transposed interval yields exactly the same
  result as the correctly ordered one, and a non-positive dispersion returns `NA`, since
  it admits no sensible recovery. Seventeen routes were affected across four families:
  the mean-difference routes (`es_from_md_sd()` / `_se()` / `_ci()` and their
  `es_from_ancova_md_*()` counterparts), the pooled-SD routes
  (`es_from_means_sd_pooled()`, `es_from_ancova_means_sd_pooled_adj()` /
  `_crude()`), the regression routes (`es_from_linreg_b_se()` / `_ci()`), and the
  ratio and user-input routes (`es_from_or_se()` / `_ci()`, `es_from_rr_se()` / `_ci()`,
  `es_from_rd_se()` / `_ci()`, `es_from_user_crude()` / `es_from_user_adj()`, the last
  of which could also return a transposed interval). The guards now live in one internal
  helper rather than inline in each route, and
  `tests/testthat/test-input-guards.R` sweeps every exported route that takes a
  dispersion or an interval, so a new route cannot acquire the defect unnoticed. Valid
  inputs are unchanged -- verified identical on `df.haza` across all nine measures --
  and a genuinely negative mean difference still yields a negative effect size.

  Inside `convert_df()` the default `correct_inputs = TRUE` already blanked these
  values, so the pipeline default was never affected. The exposure was direct calls to
  the exported `es_from_*()` functions and `convert_df(correct_inputs = FALSE)`, where
  raw values are deliberately preserved: that combination previously returned
  `g = -1.985` for a study whose *g* is `+1.985`, with a healthy positive standard error
  and an `[INVALID]` flag that named the column but did not stop the number being
  emitted.
- **The CI-width quality check (A6) mis-expected the degrees of freedom on adjusted
  rows.** `.es_from_d()` builds an adjusted *d*/*g* interval on
  `qt(.975, n_exp + n_nexp - 2 - n_cov_ancova)`, but the check compared it against the
  crude `N - 2` width, so a perfectly self-consistent interval could be flagged
  `[DISCORDANT] CI width inconsistent with SE`. The gap only clears the check's own
  tolerance at small samples with several covariates -- 17% against a 10% tolerance at
  `N = 9, q = 3` -- but it is a pure false positive when it does. User-entered rows are
  unaffected: their source's CI construction is unknown and is still checked against
  the Wald-z width.
- **`etasq_adj` gains its own `reverse_etasq_adj` column**, matching
  `reverse_ancova_means` / `reverse_ancova_md`. Previously the adjusted eta-squared route
  shared `reverse_etasq` with the crude one, so a study reporting the two in opposite
  directions could not be encoded. Like every other `reverse_*` column it defaults to
  `FALSE` when absent or `NA`. **This changes existing behaviour**: a dataset that set
  `reverse_etasq = TRUE` and relied on it flipping the adjusted eta-squared route as well
  must now also set `reverse_etasq_adj = TRUE`.
- **`cov_outcome_r` is now range-checked.** Values outside `[-1, 1]` are flagged
  `[INVALID]` and set to missing like any other correlation column. `|R| = 1` is
  additionally caught (new check **V34**): it is in range but non-identified, since
  eq. 12.24 divides by `sqrt(1 - R^2)` and eq. 12.26 multiplies by `(1 - R^2)`.
  Of the 16 routes that accept `cov_outcome_r`, 14 returned an effect size of 0 with a
  standard error of 0 (infinite meta-analytic weight) at `|R| = 1`, and the two that
  build the SMD from an already marginal SD
  (`es_from_ancova_means_sd_pooled_crude()`, `es_from_cohen_d_adj()`) returned a
  plausible estimate with a silently deflated SE and no flag. Note that V34 is a
  `convert_df()` input check: a **direct** call to one of these functions with
  `cov_outcome_r = 1` still returns `0`, as it does for any other non-identified input.

## Documentation

- `cov_outcome_r` is now documented as the pooled **within-group** correlation (the R
  for which `MSE_ancova = MSW (1 - R^2)`) across all 13 ANCOVA entry points, with an
  explicit warning against supplying the total-sample correlation or the square root of
  the whole model's R-squared. Both mistakes bias the effect size *and* its standard
  error by the same factor, leaving the *p*-value unchanged, so no quality flag can
  detect them.
- `es_from_cohen_d_adj()`'s `cohen_d_adj` now states that the value must be
  standardized on the **marginal** (unadjusted) within-group SD, and how to rescale a
  residual-scale *d* if that is what the source reported.

## New features

- **`es_formulas()`** reports the effect size obtained under every alternative
  conversion formula, with one row per (comparison, conversion parameter, formula).
  Whereas `convert_df(main_es = FALSE)` evaluates variability in the **source
  statistics** from which an effect size is derived, `es_formulas()` evaluates the
  structural uncertainty of the **analytical transformation** applied to those
  statistics: the five methods for converting an odds ratio into a risk ratio
  (`or_to_rr`), the four distinct pre-post standardisation methods (`pre_post_to_smd`,
  whose `"cooper"` and `"morris_drm"` values name the same formula), and the
  `rr_to_or` / `or_to_cor` / `cor_to_smd` / `smd_to_cor` / `smd_denom` / `smd_var` /
  `prop_to_es` / `alpha_to_es` / `icc_to_es` alternatives. The resulting range can be
  substantial and was previously unaccounted for: a single `or = 2.5` with its
  confidence interval yields risk ratios from 1.58 to 2.50 depending solely on
  `or_to_rr`. Because alternative formulae encode different statistical assumptions
  rather than algebraic equivalencies, discrepancies among them constitute
  methodological uncertainty rather than extraction errors; results are therefore
  labelled **`[FORMULA]`** and never `[DISCORDANT]`, and must not be interpreted as a
  metaDETECT quality signal. Parameters whose options alter the *analysis scale* rather
  than the formula applied on a fixed scale (`alpha_to_es`, `icc_to_es`, `prop_to_es`)
  report their values side by side but return `deviation` and `spread` as `NA`, since
  direct mathematical comparison between `ln(1 - alpha)` and `alpha` is invalid.
  **The rows returned must never be pooled**: they are structurally dependent
  re-evaluations of identical primary study data, not independent estimates.
- `summary()` gains a **`formulas`** argument (default `FALSE`). When `TRUE`, a concise
  `[FORMULA]` message is appended to the flags column for each affected comparison, and
  the complete table of alternative estimates is attached as the `"formulas"` attribute.
  Disabled by default, so existing output is unchanged.

## Fixes that change numeric output

- `or_to_rr`/`or_to_cor` could return a standard error belonging to a **different row**.
  The conversion helpers return a 1x4 matrix, except the `or_to_rr = "dipietrantonj"`
  branch, which returned a data.frame on **every** path, whether or not its 2x2
  reconstruction succeeded; `mapply()` then produced a list-matrix instead of a numeric
  matrix, and the positional extraction in `es_from_or_se()` recycled values across
  rows. With `or_to_rr` supplied as a per-row column, one `dipietrantonj` row anywhere
  but the first made every row inherit the *first* row's `logrr_se`, leaving each
  standard error inconsistent with its own confidence interval. On four studies whose
  correct standard errors are 0.235, 0.250, 0.751 and 0.825, all four were reported as
  0.235. The point estimates were unaffected. Order-dependent (a leading
  `dipietrantonj` row left the others intact) and invisible to the check that compares a
  confidence interval's width against its standard error, which is skipped for
  exp-scale measures. Any analysis in which `or_to_rr` resolved to
  `dipietrantonj` for at least one comparison should have its weights re-checked. Fixed
  at both ends: the failing branch now returns a matrix, and extraction goes through a
  length-safe helper.

## Fixes to `convert_df(main_es = FALSE)` (the per-estimation-route view)

This branch had no test, vignette or example coverage, so `R CMD check` had never
executed it in any release.

- `es_selected = "minimum"`/`"maximum"` combined with `main_es = FALSE` overwrote
  **every** route row with the same comparison-level value, leaving the row count
  inflated k-fold while destroying all per-route information. Pooling the result shrank
  the standard error by roughly `sqrt(k)` and, because route counts differ between
  studies, moved the pooled point estimate. The combination is now rejected with an
  explanatory error.
- `summary()` no longer deletes comparisons that no estimation route could reach. They
  were silently dropped, so the row count did not match the input, the console banner
  reported a vacuous `ES estimated: n/n (100%)`, and `es_guidance` -- which only writes
  where the effect size is `NA` -- became unreachable, so no guidance was ever produced
  in this view. On `df.haza` with `measure = "g"` one comparison is reachable by no
  estimation route at all: it was dropped entirely, and is now returned with its
  guidance message (row count 288 to 289). This also affected the documented
  `split_adjusted = TRUE, format_adjusted = "long"` mode.
- `summary()` no longer errors when no comparison is estimable at all (a sheet holding
  only sample sizes, i.e. the first thing a new user tries).
- Cross-row checks are now decided on **one representative row per comparison** -- the
  row carrying the route the hierarchy actually selects -- and the verdict broadcast,
  instead of treating each estimation route as an independent study. Taking the first
  row of the comparison instead, as an earlier draft of this fix did, is not equivalent:
  the rows of a comparison appear in a fixed internal order that does not track the
  hierarchy, so the verdict was computed from an estimate the user never sees. On twelve
  studies whose twelfth reports a *t* statistic inconsistent with its own means and
  standard deviations, `main_es = TRUE` flagged it as an outlier and the route view
  reported nothing; both now flag it under `hierarchy = "student_t > means_sd"`, and
  both correctly stay silent under `"means_sd > student_t"`, where the selected route is
  sound.
  Previously a comparison matched its own `study_id` k times and was reported as a
  duplicate of itself, recommending `aggregate_df()` or dropping rows -- both
  destructive (12 studies with 12 distinct identifiers produced 48 duplicate flags).
  More seriously, the cross-row outlier pools for the effect size, its
  standard error and its spread were inflated by a study's
  own replicates, which could **mask** a genuine outlier, and the direction-conflict
  check counted routes toward its study floor. The checks for mixed
  NNT types and for mixed pre-post standardizers now decide pool-level
  properties from the routes that would actually be selected, so comparisons with no
  pre/post data are no longer told to change `pre_post_to_smd`.
- `split_adjusted = FALSE` is now honoured: crude and adjusted routes are compared in a
  single pool rather than always being split. (`format_adjusted` remains inert in this
  view -- the output is one row per route by construction.)

## Other fixes

- The correlation and Fisher's *z* obtained from a 2x2 table came back as `NA`, with
  nothing said, when the `mvtnorm` package was absent. The tetrachoric correlation is
  solved by `metafor::escalc(measure = "RTET")`, which needs `mvtnorm`; because
  `mvtnorm` is a *suggested* rather than an imported dependency of `metafor`, installing
  `metafor` does not bring it in and an ordinary installation can lack it. The failure
  was indistinguishable from data that genuinely cannot support a tetrachoric
  correlation. metaConvert now reports it once per session and names the remedy. This
  affects every route that reaches the tetrachoric conversion -- `es_from_2x2()`,
  `es_from_2x2_sum()`, `es_from_2x2_prop()`, `es_from_phi()`, `es_from_chisq()` and
  `es_from_chisq_pval()`. All other effect size measures were, and remain, unaffected.
- A user column named `blank` or `dat_long` crashed `summary()`. Both are internal
  merge-carrier names and `merge()` inferred them into the join key, producing a 0-row
  result. `blank` broke the **default** `main_es = TRUE` path. Both merges now pass
  `by = "row_id"` explicitly, and `.generate_df()` asserts its `row_id` contract.
- `es_consistency` labelled `dispersion_es` as `SD:`; it is the maximum absolute
  deviation from the median, not a standard deviation, and is now labelled `max dev:`.
- In the route view the console banner reported route rows as "studies" and announced
  "Methods selected" although `main_es = FALSE` selects nothing; it now reports
  comparisons and routes separately.
- `or_to_cor = "bonett"` (the `convert_df()` default) now actually runs. Bonett's
  coefficient needs `small_margin_prop`, a column that had no auto-derivation, so when
  it was left blank -- the normal case -- the row silently kept the `"lipsey_cooper"`
  correlation computed earlier, with no message and no flag. `small_margin_prop` is now
  derived as `min(n_exp, n_nexp, n_cases, n_controls) / n_sample` when the margins are
  available and consistent (Bonett & Price 2005 define it as the smallest marginal
  proportion); a user-supplied value still takes precedence. The two margin pairs
  (`n_exp`/`n_nexp` and `n_cases`/`n_controls`) each sum to `n_sample`, so supplying
  `n_sample` plus *one member of each pair* is enough -- all four combinations work and
  the missing margins are back-filled. **This changes the Pearson
  `r` and Fisher's `z` obtained from an odds ratio.** Example (`or = 2`,
  `logor_se = 0.2`, `n_exp = n_nexp = 50`, `n_cases = 40`, `n_controls = 60`):
  `r` 0.187680633693 -> 0.258600822080, `r_se` 0.0522456910529 -> 0.0715514798191.
  The shift is always an increase in `|r|`, and its size depends on the margins: over a
  grid of odds ratios from 1.2 to 6 crossed with exposed fractions of 0.3 to 0.7 and case
  fractions of 0.1 to 0.5 at `n_sample = 200` (175 cells), it ranges from +8% to +49%,
  with a median of +36%. Cohen's *d* and Hedges' *g* derived from an OR are unchanged (they use the
  logistic transform, not `or_to_cor`), and the 2x2 route is unaffected. The three
  worked examples of Bonett & Price (2005) now reproduce to every printed digit.
- A row that requests `or_to_cor = "bonett"` but cannot use it (margins absent,
  non-positive, or not summing to `n_sample`) now keeps the `lipsey_cooper` value for
  `r`/`z` and reports the substitution, in every case. Previously such a row returned
  the `lipsey_cooper` value or `NA` depending on its position in the data frame and on
  whether any other row in the same call qualified, because the blanking step ran only
  when at least one row did. Deriving `small_margin_prop` (above) would have made that
  blanking effective for the first time, so a margin-poor row that release 2.0.0
  estimated would have come back as `NA` and the study would have been dropped from the
  analysis. Nothing is now lost: the correlation is still returned, obtained by the
  documented fall-back method, and a message names the affected rows and the inputs
  `"bonett"` would need. Supply `small_margin_prop`, or `n_sample` with one of
  `n_exp`/`n_nexp` and one of `n_cases`/`n_controls`, to use `"bonett"` on those rows.
- `reverse_2x2 = TRUE` reflected the tetrachoric confidence interval incorrectly: the
  bounds were negated without being swapped, so the interval came back inverted
  (`lo > up`) and the point estimate fell outside it.
  `es_from_2x2(143, 52, 41, 164, reverse_2x2 = TRUE)` returned `r = -0.745874586075`
  with CI `(-0.659120734630, -0.832628437520)`; it now returns
  CI `(-0.832628437520, -0.659120734630)`. Point estimates and standard errors are
  unchanged, and only the correlation and Fisher's *z* intervals were affected.
  This concerned every route that reaches the tetrachoric conversion through the 2x2
  table, not `reverse_2x2` alone: `es_from_2x2()` and `es_from_2x2_sum()`
  (`reverse_2x2`), `es_from_2x2_prop()` (`reverse_prop`), `es_from_phi()`
  (`reverse_phi`), `es_from_chisq()` and `es_from_chisq_pval()` (`reverse_chisq`,
  `reverse_chisq_pval`), and the corresponding `convert_df()` columns. For example
  `es_from_phi(phi = 0.35, n_cases = 90, n_exp = 100, n_sample = 200,
  reverse_phi = TRUE)` returned `r = -0.513079984` with CI
  `(-0.336972867, -0.689187101)` and now returns `(-0.689187101, -0.336972867)`.
- `or_to_rr = "dipietrantonj"` no longer aborts the whole `convert_df()` run when the
  2x2 reconstruction has no real solution for a single row (which is routine when the
  outcome is not rare). That row now yields `NA` and a warning explaining why.
- Standard errors that are zero or negative no longer enter the IQR pool used by the
  cross-row standard-error outlier check, where they acted as spurious extreme-low
  observations and could mask genuine outliers.
- Several quality-flag messages contained an internal `"; "`, the separator used to join
  flags, which split them into untagged fragments in the `flags` column. The messages
  were reworded; no check changed its firing conditions.

## New diagnostic checks and options

- `flag_options$flag_group`: one or more input-column names used to scope the
  cross-row quality checks (ES/SE/SD outliers, direction conflict, study duplication,
  NNT-type and standardizer mixing). With the default `NULL` the whole dataset remains
  one pool. Intended for multivariate / multi-outcome datasets, where a repeated
  `study_id` or a deviating estimate is only meaningful within a group of comparable
  rows.
- New input-validation check: a reported `prop` inconsistent with its own
  `n_cases`/`n_sample`, or (when `n_cases` is absent) not achievable as a fraction of an
  integer `n_sample` -- the count analogue of the GRIM test. Rounding aware, warn only,
  data preserved.
- New post-computation check: a raw-scale correlation confidence-interval bound escaping
  `[-1, 1]` while the point estimate is valid.
- `or_to_rr` / `rr_to_or = "dipietrantonj"` now warn when several candidate 2x2 tables
  are compatible with the reported effect and `baseline_risk` is missing, instead of
  silently taking the first one.

## Performance

- The tetrachoric reconstruction is memoised and the correlation-to-SMD conversion is
  genuinely vectorised instead of being reached through a per-row `mapply()`.
  `es_from_2x2()` on 2000 rows drawn from 60 distinct tables: 44.8 s -> 1.3 s;
  `es_from_pearson_r()` on 10000 rows: 22.7 s -> 0.6 s. Input with no repeated rows is
  not slowed down. Results are unchanged.

## Documentation

- Corrected the Rd formula for Bonett's coefficient in `es_from_or_se()` (the `1 -` sat
  inside the `/5` numerator, and raw counts were shown where the code uses proportions),
  and a missing square in the Digby Fisher-*z* confidence interval.
- `or_to_rr = "grant"` / `rr_to_or = "grant"`: the estimator is credited to Zhang & Yu
  (1998, *JAMA* 280:1690-1691), now cited. Grant (2014) is retained for the
  communication framing. The argument value `"grant"` is unchanged.
- `es_from_or()`: the imputed OR standard error is the square root of the mean
  *variance* over the compatible tables, not the mean standard error.
- `es_from_ancova_t()`, `es_from_ancova_f()`, `es_from_ancova_t_pval()`,
  `es_from_ancova_f_pval()`, `es_from_ancova_md_se()`, `es_from_ancova_md_ci()`,
  `es_from_ancova_md_pval()`, `es_from_ancova_means_se()` and `es_from_etasq_adj()`: the
  note claiming the point estimate remains unbiased under covariate imbalance was wrong
  on these routes, which recover the standardizer with the covariate-leverage term set
  to zero. Both the effect size and its standard error are attenuated by the same factor,
  so the p-value is unaffected but the magnitude is understated. The size of the
  attenuation depends on how much covariate leverage the supplied input carries, which
  splits the routes into two tiers. The seven that invert the *combined* mean-difference
  standard error -- `es_from_ancova_t()`, `es_from_ancova_f()`,
  `es_from_ancova_t_pval()`, `es_from_ancova_f_pval()`, `es_from_ancova_md_se()`,
  `es_from_ancova_md_ci()`, `es_from_ancova_md_pval()` -- together with
  `es_from_etasq_adj()`, are attenuated by `1/sqrt(1 + D/(1/n_exp + 1/n_nexp))`: about 3%
  at a standardised covariate imbalance of 0.5 and 11% at 1.0. `es_from_ancova_means_se()`
  sees only each arm's own leverage and is attenuated about half as much, by
  `1/sqrt(1 + n_j D/4)` -- 5.8% at an imbalance of 1.0. All figures are negligible in
  randomised designs. Prefer `es_from_ancova_means_sd()` or `es_from_ancova_md_sd()` when
  adjusted means, or an adjusted mean difference with the residual SD, are reported:
  those routes are unattenuated.
- The `flags` column was documented under the name `es_flags`, which does not exist.
- `es_from_or_se()`: the documented closed form for the number needed to treat omitted
  the `(1 - or)` factor in its denominator. Because that factor is negative for
  `or > 1`, the printed expression did not merely misscale the result, it returned the
  wrong sign: at `or = 2`, `baseline_risk = 0.2` the package gives `nnt = -7.5` (a
  number needed to harm) where the documented shortcut gave `+7.5`. The computed values
  were always correct; only the documentation was wrong.
- `es_from_means_sd_pre_post_single_group()` and the other pre-post converters: the
  documented sampling variance for `pre_post_to_smd = "morris_dav"` was the
  homoscedastic form that the implementation deliberately rejects. The variance in use
  is Bonett (2008, eq. 10), equal to metafor's `SMCRPH` rather than `SMCRP`; the
  documented expression is anti-conservative whenever the pre and post standard
  deviations differ. The Rd now states the implemented expression and says why the
  homoscedastic form is not used.
- `es_from_rr_se()`: the four `rr_to_or` branches listed the wrong required inputs. The
  default branch was labelled `or_to_rr = "metaumbrella_exp"` (neither the argument nor
  a valid value; it is `rr_to_or = "metaumbrella"`) and was said to need
  `n_exp`/`n_nexp` where it needs `n_cases`/`n_controls`; the `dipietrantonj` branch
  listed `logrr_ci_lo` twice and omitted `n_exp`/`n_nexp`; the `grant` branch was said
  to produce a risk ratio rather than an odds ratio. The worked example followed the
  incorrect requirement and so returned `NA` for the odds ratio, risk difference and
  number needed to treat that its own value table promises; it now supplies the inputs
  the default branch uses.
- Nine help pages linked to `metaconvert.org/html/input.html`, which does not exist
  (the other 73 use `metaconvert.org/input.html`). The dead address also appeared in
  the table returned by `see_input_data()`. All 82 occurrences were corrected.

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
