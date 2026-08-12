# metaConvert 2.0.1

## New features

- **`es_formulas()`** reports the effect size obtained under every alternative
  conversion formula, with one row per (comparison, conversion parameter, formula).
  Whereas `convert_df(main_es = FALSE)` evaluates variability in the **source
  statistics** from which an effect size is derived, `es_formulas()` evaluates the
  structural uncertainty of the **analytical transformation** applied to those
  statistics: the five methods for converting an odds ratio into a risk ratio
  (`or_to_rr`), the five pre-post standardisation methods (`pre_post_to_smd`), and the
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
  `dipietrantonj` row left the others intact) and invisible to the A6 CI-width check,
  which is skipped for exp-scale measures. Any analysis in which `or_to_rr` resolved to
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
  where the effect size is `NA` -- became unreachable (51 non-empty guidance messages on
  `df.haza` fell to 0). This also affected the documented
  `split_adjusted = TRUE, format_adjusted = "long"` mode.
- `summary()` no longer errors when no comparison is estimable at all (a sheet holding
  only sample sizes, i.e. the first thing a new user tries).
- Cross-row checks are now decided on **one representative row per comparison** and the
  verdict broadcast, instead of treating each estimation route as an independent study.
  Previously a comparison matched its own `study_id` k times and was reported as a
  duplicate of itself, recommending `aggregate_df()` or dropping rows -- both
  destructive (12 studies with 12 distinct identifiers produced 48 duplicate flags).
  More seriously, the ES/SE/spread outlier pools (D1/D2/D3) were inflated by a study's
  own replicates, which could **mask** a genuine outlier, and the direction-conflict
  check (G) counted routes toward its study floor. E4/E6/E7 now decide pool-level
  properties from the routes that would actually be selected, so comparisons with no
  pre/post data are no longer told to change `pre_post_to_smd`.
- `split_adjusted = FALSE` is now honoured: crude and adjusted routes are compared in a
  single pool rather than always being split. (`format_adjusted` remains inert in this
  view -- the output is one row per route by construction.)

## Other fixes

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
  Over a grid of realistic odds ratios and margins the shift in `|r|` ranges from +8%
  to +58%. Cohen's *d* and Hedges' *g* derived from an OR are unchanged (they use the
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
  unchanged. The same correction was applied to the internal phi and chi-squared
  helpers.
- `or_to_rr = "dipietrantonj"` no longer aborts the whole `convert_df()` run when the
  2x2 reconstruction has no real solution for a single row (which is routine when the
  outcome is not rare). That row now yields `NA` and a warning explaining why.
- Standard errors that are zero or negative no longer enter the IQR pool used by the
  cross-row SE-outlier check (D2), where they acted as spurious extreme-low
  observations and could mask genuine outliers.
- Several quality-flag messages contained an internal `"; "`, the separator used to join
  flags, which split them into untagged fragments in the `flags` column. The messages
  were reworded; no check changed its firing conditions.

## New features

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
  to zero. Both the effect size and its standard error are attenuated by the same factor
  (about 3% at a standardised covariate imbalance of 0.5, 11% at 1.0), so the p-value is
  unaffected but the magnitude is understated. Negligible in randomised designs. Prefer
  `es_from_ancova_means_sd()` or `es_from_ancova_md_sd()` when adjusted means, or an
  adjusted mean difference with the residual SD, are reported.
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
