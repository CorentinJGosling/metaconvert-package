# Audit remediation roadmap

Consolidates two independent audits of the `simulations/` programme and the
`metaConvert` package it validates. Items are ordered for execution: each is
self-contained, but Phase 0 must land first and Phase 3.5 depends on 3.1.

Every item carries a **test unit**. Package items get a `testthat` file under
`tests/testthat/`. Simulation items get either an in-script `stopifnot()`
invariant or a file under `simulations/tests/`, since that tree has no suite yet.

Status legend: ☐ not started · ◐ in progress · ☑ done

---

## Phase 0 — Safety net (blocking)

### 0.1 ☑ `simulations/` is untracked — DONE (commit `ef38361`, branch `audit-remediation`)

> Landed with two additions found during execution: the root `.gitignore`'s bare
> `archive/` and `app` patterns match at **any depth**, so both
> `simulations/archive/` and the whole of `simulations/app/` needed explicit
> re-inclusion or they would have been silently dropped. 90 files, 17 MB tracked;
> `papers/`, `students/` and `data/raw/` excluded. Not pushed.

**Symptom.** `.gitignore:30` contains `simulations`; `git ls-files simulations`
returns 0 files. The entire Monte Carlo programme — nine studies, the harness,
the app, the aggregates — is one `git clean` from gone, and no fix below is
recoverable until this changes.

**Lay version.** The work exists only on this disk. Nothing is backed up by git.

**Strategy.** Remove `simulations` from `.gitignore`, keep `^simulations$` in
`.Rbuildignore` (so it never ships to CRAN), add a `simulations/.gitignore` that
excludes `data/raw/*.rds` (large, regenerable, seeded) but **keeps**
`data/aggregated/*.csv` (the shipped result). Commit before touching anything else.

**Test unit.** `tests/testthat/test-buildignore.R` — assert `simulations/` is
excluded from the built tarball (`pkgbuild` or a `.Rbuildignore` regex check), so
un-ignoring it in git cannot leak it into the package.

---

## Phase 1 — Package bugs (code changes)

### 1.1 ☐ `or_to_rr = "metaumbrella_exp"` picks the wrong 2×2 branch

**Symptom.** Coverage 0.497 in study 04 — the worst of any shipped option.

**Lay version.** To turn an odds ratio into a risk ratio, this method has to guess
what the original 2×2 table was. It searches all tables compatible with the
reported OR and its standard error. When the two arms are the same size, **two
different tables fit equally well** — one is the other rotated 180°
(`(a,b,c,d) → (d,c,b,a)`), which preserves both the OR and its variance exactly.
The code breaks the tie with `order(...)[1]`, i.e. *whichever came first in the
enumeration*. That is a coin flip decided by loop order, and because `order()` is
stable and low event rates push the rotated candidate to a lower index, the wrong
one wins **~61–65%** of the time (81% in the worst cell). The wrong branch returns
`RR = OR / RR_true`.

**Scope, which the earlier write-up got wrong.**
- The `metaumbrella_cases` sibling has the **mirror** defect (non-identified when
  `n_cases == n_controls`). It scores 0.975 only because study 04's low baseline
  risks keep those margins unbalanced.
- Study 04's grid is `p_exp = 0.5` **only**, so equal arms occur in 100% of
  replicates. At `p_exp = 0.40` coverage is 0.955. The 0.497 headline is
  conditional on a grid constant — see 4.2.
- "The disambiguating `n_cases`/`n_controls` are discarded" is not the fix: that
  is the documented contract of the `_exp` variant, and `_cases` fails the same
  way with the roles swapped.

**Strategy.** Detect non-identifiability rather than silently resolving it.
Compare the top candidates; when two are within tolerance on `(var_sim - var)^2`
and are rotations of each other, either (a) return `NA` with a warning, or
(b) pick by an explicit, documented rule. The package already warns in exactly
this situation in the `dipietrantonj` branch
([internal_multiple_formulas.R:218-227](../R/internal_multiple_formulas.R#L218-L227)),
so the silence here is an internal inconsistency, not a considered contract.
Decide (a) vs (b) before coding.

**Test unit.** `tests/testthat/test-or-to-rr-identifiability.R` — construct a
table with `n_exp == n_nexp`, assert the rotation is bit-exactly tied, and assert
the new behaviour (NA + warning, or the documented pick). Mirror it for
`metaumbrella_cases` at `n_cases == n_controls`.

---

### 1.2 ☑ V6's column list is 64% phantom — DONE (commit `6904aa0`)

> Landed. Precise impact, narrower than first stated: the 4 surviving real names
> (`mean_pre_*`, `mean_change_*`) still caught pre/post and mean-change rows, so the
> live regression class is rows whose pre/post datum is a **paired test statistic**
> (`paired_t_*`, `paired_f_*`, p-value variants) — still the worst class, since the
> assumed r scales the point estimate there. Both call sites now share
> `.r_consuming_columns()` (30 names, 0 phantom); the verbose list was itself
> incomplete (missing `mean_pre_se_*`, `mean_pre_ci_*`, `mean_change_ci_*`).
> 39 new assertions; full suite 1478 pass / 0 fail.

### 1.5 ☐ NEW — `data_extraction_sheet()` hands users a column name nothing reads

Found by the phantom-column sweep run for 1.2 (15 agents over the 10 files that
carry hardcoded column lists). **Exactly one further instance of the bug class
exists**, and it is user-facing.

**Symptom.** [data_extraction.R:57](../R/data_extraction.R#L57) puts
`all_info_expected` in `cols_req`, spliced into the header of **every** measure
branch of the extraction sheet. That string appears at that one line and nowhere
else in `R/`. Every consumer reads `info_expected`
([internal_check_data.R:47](../R/internal_check_data.R#L47), `:195`, `:333`, `:341`;
[functions_summary.R:598](../R/functions_summary.R#L598), `:635`;
[metaConvert.R:249](../R/metaConvert.R#L249)).

**Verified empirically.** A user filling the sheet exactly as generated
(`all_info_expected = "means_sd"`) gets `info_expected` **absent from the summary
output entirely** — `.check_data()` manufactures it all-NA and `functions_summary.R`
drops it. The documented QC feature (`man/summary.metaConvert.Rd:71`) has never
worked for anyone who used the generated sheet.

**Strategy — decision needed.** Emit `info_expected` from the sheet **and** have
`.check_data()` accept `all_info_expected` as a deprecated alias (renaming it with a
one-time message), so extraction files already saved in the wrong name start working
rather than silently continuing to fail.

**Test unit.** `tests/testthat/test-extraction-sheet-column-names.R` — assert every
name `data_extraction_sheet()` emits, for every `measure`, is a name `.check_data()`
recognises; assert the round trip populates `info_expected`; assert the legacy alias
is accepted. This generalises the bug class into a permanent guard.

**Cleared by the same sweep (no action):** `.positive_columns()`'s `logirr_se` and
`.ci_triplets()`'s `logirr`/`loghr` triplets name output-only quantities — IRR and HR
have no raw-data input columns — so they are dead weight behind an existing
`%in% colnames(x)` guard, not a coverage hole. `internal_guidance.R`'s ~80-method
map, `internal_check_data.R`, `functions_summary.R`, `internal_generate_df.R`,
`main_aggregate_df.R` and `main_compare_df.R` are clean.

**Symptom.** [main_convert_df.R:478-483](../R/main_convert_df.R#L478-L483) lists 11
columns; **7 do not exist**: `mean_post_exp`, `mean_post_nexp`,
`mean_pre_single_group`, `mean_post_single_group`, `mean_change_single_group`,
`paired_t`, `paired_f`. (Real names are `paired_t_exp`/`paired_t_nexp` etc.)

**Lay version.** When a user does not report the pre-post correlation, the package
substitutes 0.8 and is supposed to flag that. The flag decides "does this row have
pre/post data?" by looking for column names — and most of the names it looks for
were never column names. So on a row whose only pre/post data is a paired *t* or
*F*, the check finds nothing, `r_defaulted` stays FALSE, and **the V6 flag, the
`r_defaulted` attribute, and the downstream `(r-sensitive: ...)` annotations on
A6/E2/E2b all fail together**. Those are the worst rows to miss: on paired-t/F
routes the assumed correlation scales the **point estimate**, not just the SE.

**Strategy.** One line. The correct list already exists 130 lines earlier as
`.r_consuming_cols` ([main_convert_df.R:350-359](../R/main_convert_df.R#L350-L359)),
used by the verbose message. Hoist it and reuse it in both places. Also add the
`*_sd_*`/`*_se_*`/`*_pval_*` variants so a row carrying only `mean_change_sd_exp`
is caught.

**Test unit.** `tests/testthat/test-V6-r-defaulted-columns.R` — one row per
pre/post data shape (pre/post means, mean-change, mean-change SD only, paired t,
paired F, single-group), assert V6 fires and `attr(res, "r_defaulted")` is TRUE for
each; assert a row with no pre/post data does **not** fire.

---

### 1.3 ☐ `table_2x2_to_cor` is inert in `convert_df()` and silently accepts anything

**Symptom.** `convert_df(d, measure = "r", table_2x2_to_cor = "banana")` runs
without error and returns byte-identical output. The argument is commented out at
all three call sites ([main_convert_df.R:1118, 1123, 1128](../R/main_convert_df.R#L1118))
and excluded from the per-row column loop (`:394`). The guard at
[es_from_2x2.R:95-102](../R/es_from_2x2.R#L95-L102) fires only on direct calls.

**Lay version.** The manual advertises a choice the main workflow cannot make, and
a typo in it is silently ignored rather than reported.

**Compounding.** `man/convert_df.Rd:60` says "For now only 'tetrachoric' is
available", which reads as a temporary restriction rather than a dead argument.
The `stop()` message at
[internal_multiple_formulas.R:668-670](../R/internal_multiple_formulas.R#L668-L670)
advertises `cooper_delta`/`cooper_std`/`lipsey`, none of which any reachable code
implements. And there is no workaround: `es_from_phi(phi = 0.3015)` returns the
tetrachoric 0.4576, not phi.

**Strategy — decision needed.** Either **(a)** validate-and-document: reject any
value but `"tetrachoric"` in `convert_df()`, fix the Rd and the stale `stop()`
string, and state the latent-continuous assumption plainly; or **(b)** wire the
argument through and re-open a correct `phi` option (the disabled `lipsey` block
at `:611-642` has a documented margin error *and* no `return()`, so it would need
rewriting, not un-commenting). (a) is a bug fix; (b) is a feature with a
methodological argument attached — see 5.4.

**Test unit.** `tests/testthat/test-table-2x2-to-cor-validation.R` — assert an
invalid value errors through `convert_df()` (not only through `es_from_2x2()`),
and assert the tolerated set matches the documented set.

---

### 1.4 ☐ `rr_to_or = "grant"` returns a finite estimate with a NaN standard error, silently

**Symptom.** `es_from_rr_ci(rr = 1.5, rr_ci_lo = 0.9, rr_ci_up = 2.5,
baseline_risk = 0.5, rr_to_or = "grant")` returns `logor = 1.099` (a plausible
OR = 3.0), `logor_se = NaN`, and a half-open CI.

**Lay version.** Grant's formula needs `rr × baseline_risk < 1`. The *upper CI
bound* is by construction the largest of the three numbers, so it breaks first —
and the failure is wrapped in `suppressWarnings()`, so the user gets a
respectable-looking effect size with no usable uncertainty and no message.

**Note on scope.** The "SE fails ~9× more often than the point estimate" figure is
a property of one simulated draw distribution, not a general rate. The defect is
the **silence**, not the frequency; frame it that way.

**Strategy.** Emit a warning naming the condition when the CI-side transform is
undefined, and return NA consistently across the triplet (point, SE, both bounds)
rather than a mixed finite/NaN row.

**Test unit.** `tests/testthat/test-rr-to-or-grant-domain.R` — assert the
warning fires, assert no row is returned with a finite estimate and a non-finite
SE, and assert the in-domain case is unchanged (regression guard).

---

## Phase 2 — Package documentation defects

### 2.1 ☐ ANCOVA `@param cov_outcome_r` states the opposite of what the route does

**Symptom.** [es_from_ANCOVA_means.R:390-395](../R/es_from_ANCOVA_means.R#L390-L395)
tells users a mis-specified `cov_outcome_r` biases the ES and the SE "by the same
factor, so the p-value is unchanged and no quality flag can detect the error" —
on `es_from_ancova_means_sd_pooled_crude`, where the measured ratios are
**d 1.000, SE 0.463, t 2.16**.

**Lay version.** The warning is correct for the 13 routes that take a residual SD.
It was copy-pasted onto the two routes where it inverts — the ones that receive an
already-marginal SD, so the guess enters the **variance only** and nothing cancels.
On those, a wrong guess is *maximally* visible in the p-value, not invisible.

**Strategy.** Split the `@param` text: keep the existing wording on the residual-SD
routes; write a corrected paragraph for `es_from_ancova_means_sd_pooled_crude` and
`es_from_cohen_d_adj`, stating that ES is untouched, the SE is deflated (0.46 at a
0.9 guess against a true 0, coverage 0.633), and the p-value **does** move. Cross-
reference flag V22 and note it fires only when the value is *missing*, not wrong.

**Test unit.** `tests/testthat/test-ancova-cov-r-misspecification.R` — assert on
the pooled-crude route that ES is invariant to `cov_outcome_r` while SE is not
(the fact the doc denies), and the converse on `es_from_ancova_md_sd`.

---

### 2.2 ☐ `man/convert_df.Rd` on `table_2x2_to_cor` — folded into 1.3.

### 2.3 ☐ `or_to_cor` default diverges: `es_from_or_se()` alone ships `"pearson"`

**Symptom.** `convert_df()` and four of five OR entry points default to
`"bonett"`; [es_from_stand_OR.R:198](../R/es_from_stand_OR.R#L198) defaults to
`"pearson"`. `es_from_user_crude`/`_adj` also use `"pearson"`.

**Lay version.** The same table converts differently depending on which door you
came in. The pattern — one sibling out of five — reads as an oversight from the
2.0.1 default change, not a design choice.

**Scope correction.** The "18.7×" figure is a ratio of mean |bias| in one
simulation cell with each method on its own estimand, and its denominator is 1.67
MCSE from zero. The discrepancy a direct-call user actually sees between the two
defaults on the same rare-outcome table is **1.30× (0.068 r-units)**. Use that
number. Note also this is unreachable through `convert_df()`, which passes
`or_to_cor` explicitly.

**Strategy.** Align `es_from_or_se`, `es_from_user_crude`, `es_from_user_adj` to
`"bonett"`; add a NEWS entry flagging the behaviour change.

**Test unit.** `tests/testthat/test-method-defaults-consistency.R` — a formals
scan asserting every exported function's default for each shared method argument
matches `convert_df()`'s. Generic, catches future drift on all 19 arguments.

---

### 2.4 ☐ `unit_type`: documented values exclude the shipped default; no validation

**Symptom.** `?convert_df` says "must be either 'sd' or 'value'"; the default is
`"raw_scale"`; the only read is `ifelse(unit_type == "sd", ...)`. A typo silently
selects raw units.

**Strategy.** Make the default one of the documented values and validate the
argument. Arithmetic is correct — this is docs + a guard.

**Test unit.** Covered by `test-method-defaults-consistency.R` (2.3) plus an
`expect_error` on an invalid value.

---

### 2.5 ☐ `measure = "z"` returns two different transforms depending on route

**Symptom.** `smd_to_cor = "lipsey_cooper"` returns Fisher's `atanh(r)`;
`"viechtbauer"` returns the Jacobs & Viechtbauer variance-stabilising transform.
At ρ = 0.75, p = 0.5 they differ by 0.249. Undocumented in `?convert_df`.

**Lay version.** A review whose rows resolve to different routes pools two
quantities that are not on the same scale.

**Strategy.** Document it on `convert_df` and `es_from_means_sd`; consider a
metaDETECT check that fires when a `measure = "z"` pool mixes routes.

**Test unit.** `tests/testthat/test-z-scale-provenance.R` — assert the two routes
return different z for the same data and that whatever provenance marker we add
(attribute or flag) is present.

---

## Phase 3 — Simulation code bugs

### 3.1 ☐ Study 05's `vanderweele_rr_squared` computes 2·log(OR), not 2·log(RR)

**Symptom.** `gen_2x2_or()` overwrites `theta_sample` with `log(or)` at
[05_rr_to_or.R:36](studies/05_rr_to_or.R#L36); `estimate_rr_to_or()` then reads it
as the log RR at [:67](studies/05_rr_to_or.R#L67). Intent is pinned by the comment
("OR ≈ RR^2") and by `se <- 2 * dat$logrr_se` at `:68`, which is only right for the
RR reading. No defensible alternative interpretation.

**Lay version.** The method is meant to square the risk ratio. It squares the odds
ratio instead, because a variable was reused for something else two functions
earlier.

**Consequence, and why it is worse than "one bad row".** Because study 05 scores
against `theta_sample = log(OR)`, the bias is `2·logOR − logOR` = exactly
`log(OR)`, which is why mean |bias| 0.687 ≈ mean |logOR|. Correcting it does
**not** clean the row (0.687 → ~0.443, still worst in study 05) — it **inverts the
baseline-risk gradient**. The shipped file makes the method look worst at br = 0.50
(coverage 0.513), exactly where VanderWeele's minimax argument says it should be
best (corrected: 0.915). That artifact is already published to students, with an
explanation invented for it, in
`students/tasks/task-05-rr-to-or.md:470-484`.

**Strategy.** Read `2 * log(dat$rr)` (or stop overwriting `theta_sample`; prefer a
distinct `theta_sample_rr` so the two scales cannot collide again). Cross-check the
estimator against the source PDFs in
`data-raw/RR to OR and vice versa/` (`VanderWeele-Optimalapproximateconversions-2020.pdf`)
before regenerating. Then regenerate 05, and rewrite the student task section.

**Test unit.** `simulations/tests/test-study-05-invariants.R` — assert that on a
known table the candidate equals `2*log(RR)` and not `2*log(OR)`; plus a
`stopifnot()` invariant in `gen_2x2_or()` that `theta_sample` is the logOR and a
separate named column holds the logRR.

---

### 3.2 ☐ Study 03's `theta_own` is scale-matched, not estimand-matched

**Symptom.** [03_2x2_to_cor.R:153](studies/03_2x2_to_cor.R#L153) is
`theta_own <- if (meas == "r") theta_pop else atanh(theta_pop)`, so for every
`(r)` route `own` **is** `population` (bias columns byte-identical).

**Lay version.** The `own` column exists to remove estimand mismatch so that what
is left is computational error. In study 03 it does not do that: under CAT the
tetrachoric routes are still scored against phi. So `own` is correct for the
home-ground route in each file and mislabelled for the off-mechanism one.

**Strategy.** Record a genuine per-method estimand (phi for the phi routes, the
latent ρ for the tetrachoric routes, on both scales), as studies 01/09 already do.
Note study 03 uses plain `atanh()` where a VST route would need `vs_transform()` —
check before generalising.

**Test unit.** `simulations/tests/test-own-target-is-estimand.R` — for every study,
assert that `own` differs from `population` for at least one method wherever the
study registers methods with different estimands (a generic anti-regression).

---

### 3.3 ☐ `performance()` publishes cells built from 3–10 replications while advertising `n_valid` 380–996

**Symptom.** `n_valid` and `nonest_rate` are computed from the **point estimate**
only ([03_performance.R:56-64](R/03_performance.R#L56-L64)); `sum(ok_se)` and
`sum(ok_ci)` are computed, used as denominators, then discarded. Five cells report
`se_ratio` 1.22–1.68 and coverage from 3, 4, 5, 7 and 10 replications.

**Correction to the earlier write-up.** "No column says so" is false — the 26
`grant` rows have `mod_se`/`se_ratio`/`coverage`/`ci_width` all NA. The defect is
the *partial* cells, which look fully populated.

**Strategy.** Add `n_valid_se` and `n_valid_ci` to the returned frame; suppress (or
mark) cells below a minimum replication count.

**Test unit.** `simulations/tests/test-performance-counts.R` — feed a vector with
known non-finite SEs and assert the new counts are right and the point-estimate
count is unaffected.

---

### 3.4 ☐ `run_everything()` aborts on `run_99`

**Symptom.** `run_99` takes no arguments; `run_everything()` calls every
`run_[0-9]{2}` with `nrep`/`cores`. It sorts last, so nothing is lost — but the
run ends in an error.

**Strategy.** Give `run_99` the ignored signature, or exclude it from the pattern.

**Test unit.** `simulations/tests/test-runner-registry.R` — assert every
registered runner accepts `(nrep, cores)`.

---

### 3.5 ☐ Regenerate stale aggregates

**Depends on 3.1, 3.2.** Study 03's `nrep1000` files predate the `theta_own` fix
(confirmed by target-column diff, not just mtime). Study 09's predate two commits
to the package file it evaluates. Study 05 must be regenerated after 3.1.

**Test unit.** `simulations/tests/test-aggregates-fresh.R` — assert each shipped
CSV is newer than the study file and than the package files it exercises.

---

## Phase 4 — Reporting corrections (before anything is published)

### 4.1 ☐ The Haldane–Anscombe diagnosis is backwards

The claim that the +0.5 continuity correction *causes* the `se_ratio > 1` pattern
in studies 04/05 is wrong in causal direction. Removing it makes the ratio
**worse** (1.20–1.71 vs shipped 1.08–1.32). The driver is **sparsity itself**
(Woolf's variance bias); the correction partly suppresses it. The dense-region
result stands (105/360 conditions with P(any zero cell) < 1e-4 → every route
0.99–1.00) but attributes to sparsity, not the correction.

Also: **`metafor_conv2x2` is not miscalibrated.** It sits at 0.992 in the dense
region. Its full-grid 1.62 is a selection artifact — 43% of replications fail to
reconstruct and it is scored on the surviving 57%. Of the "two methods that stay
miscalibrated", only `metaumbrella_exp` does.

**Test unit.** N/A (analysis). Record the dense-region restriction as a reusable
diagnostic in `simulations/R/`.

### 4.2 ☐ `metaumbrella_exp`'s 0.497 is conditional on `p_exp = 0.5`

Studies 04 and 05 hold `p_exp = 0.5`, so equal arms occur in 100% of replicates —
precisely the non-identified case (1.1). Coverage is 0.955 at `p_exp = 0.40`. Add
`p_exp` to the grid and state the conditionality in `README.md`.

### 4.3 ☐ App: study 01b's default target is the one wrong target

The app's target-ordering heuristic defaults 01b to `fisherz_biserial`, and the
sidebar advises preferring shared targets over `own` — exactly backwards for a VST
route. (01b is not the app's default *study*, only the default target within it.)
Fix the ordering for VST-scaled studies, or suppress the mismatched shared cell.

### 4.4 ☐ README kernel table: `.single_group_pre_post_to_smd` **is** run

`README.md:83` says "never run". It is called twice per row by study 08 (once per
arm, `pool_sd = FALSE` path, `internal_multiple_formulas.R:1054`/`:1061`). Only
`.pooled_pre_post_to_smd` is unsimulated — and it is pinned by
`tests_save/checked/test-pooled-variance-calibration.R`, so "unsimulated" is the
accurate word, not "uncovered".

### 4.5 ☐ Strike two stale items

- `README.md:971-974` open item 0c ("the app is broken, 16 reads point at
  `./data_agg/`") — the fix landed; the item is stale.
- `studies/06_or_se_imputation.R:171-173` accuses study 04 of a `conv.2x2` margin
  error. The call at `04_or_to_rr.R:92-94` is correct and verified in a comment
  immediately above it. Delete the note.

### 4.6 ☐ Downgrade five findings to documented design consequences

`r_pre_post = 0.8`, `.se_from_or` mean-over-tables, tetrachoric-only, the ANCOVA
variance limitation, and the `or_to_cor` divergence are **already analysed in the
repo** (`man/convert_df.Rd:76` publishes the 68%/0.75 figure; `README.md:544-547`
and `06_or_se_imputation.R:113-118` already say the median variant is "a region
map, not a one-line replacement"). Presenting them as external discoveries
undersells the repo's own analysis and will not survive a referee who reads it.

Specific numbers to retire: "median gives 1.006" (a cancellation of +1.17 against
0.56–0.79; typical cell error ~11%, worst-case *worse* than the package's);
"1.19→1.70" (a br ≥ 0.3 subgrid; grid-wide 1.17→1.47); "18.7×" (use 1.30×, see
2.3); "all five pre/post routes deflate the SE" (`morris_dz` does not — its
`se_ratio < 1` is `emp_se` doubling).

---

## Phase 5 — Coverage gaps (new work, scoped separately)

### 5.1 ☐ `pool_sd = TRUE` kernel is unsimulated
`.pooled_pre_post_to_smd` has genuinely distinct arithmetic (own pooled
standardizer, own variance forms, own df). A "study 08b" reusing study 08's DGP
and grid takes pre/post from 1 kernel to 2.

### 5.2 ☐ Ten method arguments have zero simulation coverage
`smd_var`, `smd_denom`, `pool_sd`, `yates_chisq`, `prop_to_es`, `alpha_to_es`,
`icc_to_es`, `icc_type`, `es_selected`, `selection_auto`/`hierarchy`. The
psychometric family is entirely unsimulated, including the ICC SE the package
itself flags as anti-conservative (V31).

### 5.3 ☐ Entry-point coverage is 1 route deep
`convert_df()` is never called in `simulations/`. Monte Carlo reaches 12 of 89
exported `es_from_*`; 27 have no coverage anywhere (sims + tests). Per argument:
`smd_to_cor` 1 of 60, `pre_post_to_smd` 1 of 19, `cor_to_smd` 1 of 9,
`or_to_cor` 1 of 7.

### 5.4 ☐ Should VanderWeele's conversions become package options?
`data-raw/RR to OR and vice versa/` holds the sources. VanderWeele's √OR needs no
baseline risk — the situation most meta-analysts converting an OR are actually in —
and matched the baseline-risk-requiring methods for common outcomes in study 04.
Currently a simulation candidate only, not a package option. Same question for
`metafor::conv.2x2` as an `or_to_rr` route.

---

## Execution order

0.1 → 1.2 → 2.1 → 1.1 → 3.1 → 1.3 → 1.4 → 3.4 → 3.3 → 3.2 → 2.3 → 2.4 → 2.5 →
3.5 → 4.x → 5.x

Rationale: safety net first; then the two one-line fixes with real user impact
(1.2, 2.1) to bank early wins; then the one unambiguous bug (1.1) and the one
corrupting published numbers (3.1); then the remaining code and doc work; then
regeneration, which must come after every study fix; then reporting; then new
studies.
