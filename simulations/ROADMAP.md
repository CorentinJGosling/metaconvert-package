# Audit remediation roadmap

Consolidates two independent audits of the `simulations/` programme and the
`metaConvert` package it validates. Items are ordered for execution: each is
self-contained, but Phase 0 must land first and Phase 3.5 depends on 3.1.

Every item carries a **test unit**. Package items get a `testthat` file under
`tests/testthat/`. Simulation items get either an in-script `stopifnot()`
invariant or a file under `simulations/tests/`, since that tree has no suite yet.

Status legend: ☐ not started · ◐ in progress · ☑ done

**Verification rule (learned the hard way, twice).** Before declaring any item done:
run `tests_save/checked/` as well as `tests/testthat/` — the main suite does not cover
the cross-package agreement cases and missed a 19-assertion regression in 1.1. And
check testthat's **`error`** column, not just `failed`: an erroring `test_that` block
reports `failed = 0` while silently skipping every assertion after the error. In 1.3
that hid a whole block of 24. Compare the assertion TOTAL against the previous run;
a count that falls with no failures means assertions stopped executing.

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

### 1.1 ☑ `or_to_rr = "metaumbrella_exp"` picks the wrong 2×2 branch — DONE (commits `935767a` + `b675c09`)

> **Verified:** `tests/testthat` 2228 pass / 0 fail; `tests_save/checked` 7191 pass /
> 0 fail (exactly the pre-existing baseline).
>
> ⚠️ **Process note — a regression shipped in `935767a` and was caught only by the
> archived suite.** The `es_from_or` suppression dropped the route even when the user
> had explicitly named `"or"` in `hierarchy`, breaking 19 assertions across
> `test-ES-NNT.R` and `test-ES-OR.R` while `tests/testthat` stayed green. Fixed in
> `b675c09` by suppressing only when `"or"` is not a token of the user's hierarchy.
> **Run `tests_save/checked/` BEFORE declaring an item done, not after** — the main
> suite does not cover the cross-package agreement cases.
>
> General rule this established: "don't compute a dominated estimate" and "don't remove
> what the user explicitly asked for" conflict, and the second wins.

> **Reframed during execution: don't tie-break, SOLVE.** With a second margin all four
> margins are known, so the OR determines the table by a quadratic and the reported
> variance is not used at all. Branch recovery 1.000 at every event rate 0.03–0.97;
> mean |logRR| error 0.2295 → 0.0002. Both variants fixed; `_cases` had the exact
> mirror defect (34.7% wrong when `n_cases == n_controls`).
>
> Two things this overturned in my own earlier report: (a) **unequal arms were never
> safe** — the 0.2% failure rate held only for a bit-exact variance; with a rounded OR
> the current rule's MAE is 0.10–0.19 there too; (b) **problems A and B were coupled** —
> `es_from_or()` fed the enumeration an imputed variance ~1.4× wide, so the search
> missed the table before any tie arose. Solving ignores the variance and removes both.
>
> **No prior shipped.** The adversarial pass killed the minority-event rule: 3.65× worse
> on common outcomes, +44% pooled-RR bias end-to-end, and a bet on outcome labelling.
> Rows with no second margin keep today's behaviour bit-for-bit (verified).
>
> ⚠️ **Consequence for the simulation programme:** studies 04 and 05 pass `n_cases`/
> `n_controls` to `es_from_or_se()`, so both metaumbrella routes now solve exactly —
> measured bias **0.00000** on study 04's shape, against the shipped README's 0.431
> mean |bias| / 0.496 coverage. **`data/aggregated/04_*` and `05_*` are now stale** and
> must be regenerated (folded into item 3.5). The "conditional on `p_exp = 0.5`" caveat
> (item 4.2) still applies to any *pre-fix* number quoted in the write-up.

### 1.4 ☑ `rr_to_or = "grant"` returns a finite estimate with a NaN standard error — DONE

> Grant's transform needs `RR × baseline_risk < 1` for the point estimate **and both CI
> limits**; the upper limit is the largest of the three, so it fails first, and every
> `log()` was inside `suppressWarnings()`. A row could return `logor = 1.099` beside
> `logor_se = NaN` and a half-open interval. Now NA throughout, with a warning naming
> the offending value and pointing at `metaumbrella`/`transpose`. In-domain output is
> bit-identical; swept 400 random draws for a mixed finite/NA quartet — **0 found**.
>
> **The mirror is safe and was deliberately left alone**, checked rather than assumed:
> `or_to_rr = "grant"` computes `or / (1 − BR + BR·or)`, whose denominator is positive
> for every `or > 0`, `BR ∈ (0,1)`. Verified over a 30-point grid and pinned by a test,
> so a later "fix for symmetry" cannot add a guard that can never fire.
>
> ⚠️ **Consequence for study 05:** 13 of its 360 `grant` cells reported a point estimate
> for 534–997 replications while ≥99% of their standard errors were non-finite. Those
> rows now go NA, so `nonest_rate` reports the failure instead of the cell looking
> half-populated. Folds into the item 3.5 regeneration alongside 04/05.

### 1.7 ☑ NEW — `baseline_risk` unguarded on the exported routes — DONE

Found by the literature verification (below). `baseline_risk` outside `[0, 1)` made
both `(1 - BR)` and `(1 - RR*BR)` change sign together, so the grant conversions stayed
finite and positive while returning a **negative SE and a transposed CI**, silently, in
both directions. A finiteness test cannot catch it — which is why the 1.4 guard, which
tests `all(is.finite(...))`, did not. New `.baseline_risk_or_na()` in
`internal_guards.R`, applied at all 13 entry points. Also: RD/NNT now NA when the
implied exposed risk `rr*br > 1` (was emitting `rd = -0.900`, `nnt = -1.111` from a
non-existent risk pair, with B6 silent since `|rd| < 1`); stale `stop()` advertising
`grant_2x2`/`grant_CI` corrected; `.or_to_rr` given the terminal `else` it lacked.
78 assertions.

### 1.8 ☑ NEW — `metaumbrella_exp`'s documented gate does not require what it needs — DONE

> **Decision taken: gate only in the tied cell** — a second margin is required when
> `n_exp == n_nexp` (mirror: `n_cases == n_controls`), not universally. This **reverses
> the earlier "pick a branch, don't return NA"** call, and the reversal is deliberate:
> that call was made when the alternative was a *rule that guesses*; here the alternative
> is a number that is wrong about two times in three.
>
> **The cliff is sharp, which is what makes the narrow gate defensible.** Measured hit
> rate by `|n_exp - n_nexp|`: **0.384** (0) → 0.995 (1) → 0.987 (2) → 0.995 (5) → 0.982
> (10); with the OR rounded to 2 dp, 0.379 → 0.955 → 0.967 → 0.950 → 0.951. The `_cases`
> mirror by `|n_cases - n_controls|`: **0.745** (0) → 0.995 → 0.995 → 0.993 → 0.993. One
> participant of imbalance makes the rotation inadmissible, so exact equality is the
> whole failure set. **The three tied-cell rates are not comparable to each other**
> (0.384 / 0.745 / the 0.324 first reported all come from different draw distributions);
> what reproduces is the size of the drop at zero.
>
> **Mechanism confirmed rather than assumed** (782 draws at 50/50, no second margin):
> the helper returns the true table 37.0%, its **rotation 66.6%**, neither 0.8%. So the
> failure is the rotation tie, not general imprecision. Both recovery paths are exact:
> `+ n_cases` → 1.000, `+ baseline_risk` → 1.000.
>
> **The gate keys on whether the solve fired, not on whether a margin was supplied.**
> The helpers now mark their result `solved = TRUE/FALSE`. That catches the case a
> re-derived condition would miss: a margin that cannot describe this table (a multi-arm
> case margin, or a root rounding onto an empty cell) falls through to the search and is
> just as tied.
>
> **Left in `internal_flags.R` deliberately.** D2's SE-outlier check calls the same helper
> but uses only the reconstructed *variance*. The rotation permutes the four cells, so
> `var(logOR)` is invariant — measured equal to within 2.2e-16 over 5000 tables, and the
> first measurement's 0.703 "bit-identical" rate was summation-order rounding, not a real
> gap. Gating there would push good rows onto a cruder fallback for nothing. **Noted, not
> changed:** `var(logRR)` is *not* rotation-invariant (measured 0/2000 equal), so the
> `measure = "logrr"` arm of that same check is mildly exposed — a pre-existing
> flag-normalisation question, out of scope here.
>
> **Blast radius: zero existing rows.** No shipped dataset (`df.haza`, `df.short`,
> `df.psychom`, `df.compare1/2`) carries an `or`/`logor` input column at all. Across the
> whole of `tests_save/checked` the gate fired **0 times** (traced), so the archived suite
> never covered this configuration — the new tests are its only coverage. `tests/testthat`
> 2817 pass / 0 fail / 0 error / 0 skip (2427 baseline + 390 new);
> `tests_save/checked` 7193 / 0 / 0 / 0, exactly the baseline.
>
> ⚠️ **Process note.** The first blast-radius run reported `PASS 5608` against a 7193
> baseline with 0 failures — the roadmap's own "count that falls with no failures"
> signature. Cause: I had wrapped `test_dir()` in `suppressWarnings(suppressMessages(...))`
> to keep the trace output readable. Re-run without the wrapper: 7193, and the trace
> count unaffected. **Do not wrap the suite to tidy its output.**
>
> Also fixed en route: the roxygen `\%` escape. Written as `\%` in a roxygen comment it
> reaches the Rd as `\\%` — a literal backslash followed by an Rd comment, which silently
> truncates the rest of the line. The convention in this package is a **bare** `%`.

### 1.8-old (the original write-up, kept for the record)

The 1.1 solve fixes the reconstruction **only when a second margin is supplied**. The
documented minimal input ([es_from_stand_OR.R:47](../R/es_from_stand_OR.R#L47)) is
`or + logor_se + n_exp + n_nexp`, which does not include one. Measured on that input at
**equal arm sizes** — 1:1 randomisation, the modal RCT design — **69.2%** of 400 tables
come back wrong by >0.01 in log RR (quantiles 0.181 / 0.652 / 1.258 / 1.855 / 2.575 at
the 50/75/90/95/99th, max 3.019). Worked case: true table 1/49/25/25, true RR **0.040**,
returned **0.510** — a 12.8-fold error, silently and unflagged. Supplying either
`n_cases` or `baseline_risk` recovers it exactly. `metaumbrella_cases` is 35.2% wrong at
balanced case margins.

**Strategy — decision needed.** Either tighten the gate so the method only fires when it
can be identified (a second margin or `baseline_risk`), or let it fire and warn loudly
when the arms are equal and nothing disambiguates. The first changes which rows produce
an RR; the second keeps them but marks them.

**Test unit.** Extend `test-or-to-rr-identifiability.R` with the equal-arms /
minimal-input case, asserting whichever contract is chosen.

### 1.9 ☑ NEW — study 04's `vanderweele_sqrt_or` measures a straw man — DONE

> **Decision taken: fix the CI, keep the arm.** Verified against the PDF rather than
> from memory (`pdftools::pdf_text`, p.747–749), and the check changed what I would
> otherwise have written:
>
> - **The point estimate was never the straw man.** Corollary 1 (p.747) proves √OR *is*
>   the optimal bias-ratio minimax conversion whenever the outcome probabilities lie in
>   an interval symmetric about 0.5. Only the interval was wrong.
> - **1.25 is not a constant of the method.** It is `1/sqrt(1 - 4v^2)` for probabilities
>   in `[0.5-v, 0.5+v]` (p.748 + Appendix), which is *exactly* 1.25 at the paper's
>   `[0.2, 0.8]`. Implemented as `vw_bias_ratio(w, u)` so the provenance is at the call
>   site, not a magic number.
> - **The paper's own worked example is now a test.** p.749: OR 2.3 (1.5–3.4) → RR 1.5,
>   conservative CI (1.0–2.3). Reproduced to the published decimal.
>
> **Measured effect of the fix** (n = 150/arm, 20 000 reps, coverage of the true log RR):
>
> | region | naive (forbidden) CI | fixed CI | mean width |
> |---|---|---|---|
> | inside `[0.2, 0.8]` (guarantee applies) | 0.907 | **0.999** | 0.49 → 0.93 |
> | outside the band (premise fails) | 0.739 | 0.904 | 1.00 → 1.45 |
>
> Every in-band cell lands ≥ 0.996, matching the paper's "at least 95%, in general
> conservative". The shipped aggregate has this arm at **0.931 in-band / 0.850 outside /
> 0.870 overall**, so 3.5's regeneration should move in-band to ~1.00. The width roughly
> doubles — that is the price of the conservative interval and belongs in the write-up.
>
> ⚠️ **Scope limit the fix does NOT remove, and the README must respect.** The guarantee
> holds only where **both** outcome probabilities are in `[0.2, 0.8]`. In this grid
> `p0 = br`, `p1 = rr*br`, so the guaranteed region is just `br = 0.30` × `rr ∈
> {0.75, 1, 2}` and `br = 0.50` × `rr ∈ {0.5, 0.75, 1}` — **6 of 18 (br, rr) combinations,
> 90 of 360 rows**. Every `br ≤ 0.15` cell fails outright. Both columns are in the
> aggregate, so the split is recoverable post-hoc with no code change; the coverage column
> for this arm must be read in two regions, exactly like item 4.1's dense-region
> restriction. Outside the band even the widened interval under-covers badly (0.557 at
> `br = 0.30, rr = 0.25`) because the *point estimate* is biased there (log RR bias 0.54–0.61),
> which no interval rule can repair.
>
> **`logrr_se` deliberately left as the delta-method SE** (`SE(logOR)/2`). The paper gives
> no variance, so the interval is not `es ± z·se` — safe because `performance()` takes
> coverage/width from the CI columns and `se_ratio` from the SE column independently
> ([03_performance.R:88-101](R/03_performance.R#L88-L101)), verified rather than assumed. A
> test pins the inconsistency so a later "tidy-up" that rebuilds the CI from the SE fails.
>
> **New third suite: `simulations/tests/`** (owner's call), with `run_tests.R` and
> `test-study-04-vanderweele-ci.R` — 38 assertions, all passing. The package's two suites
> do **not** execute it; run `Rscript tests/run_tests.R` from `simulations/` after touching
> anything there. Phase 3 items 3.1/3.3/3.4 already name files in this directory.
>
> **Aggregates for 04 are now stale for a second reason** (the first being 1.1). Folds into 3.5.

### 1.9-old (the original write-up, kept for the record)

`simulations/studies/04_or_to_rr.R:79-86` builds `logrr <- log(or)/2`,
`se <- logor_se/2` and a plain symmetric interval. VanderWeele (2020) p.748 explicitly
rules that out — *"it cannot be applied directly to the confidence interval of the odds
ratio to obtain 95% coverage over repeated samples of the true risk ratio"* — and
prescribes dividing/multiplying the bounds by 1.25 for a conservative interval. So the
candidate's **coverage** results describe an operation the paper forbids, not the method
it proposes. The point-estimate results are unaffected. Fix the candidate and
regenerate 04 (folds into 3.5), or drop the coverage column for that arm.

### 1.10 ☑ NEW — the 1.1 solve destroyed continuity-corrected tables — DONE

**Found by the 3.5 regeneration, not by any test.** Study 04's `metaumbrella_cases`
mean |bias| went **0.006 → 0.031** after the fixes, concentrated entirely at
`br = 0.01` (50 of 360 cells grew by >0.01; `n_valid` changed in **0** cells, so it
was not a selection effect). I had predicted this route would not move.

**Cause.** Item 1.1's `.solve_2x2_from_or()` ended with `a <- round(a)`. A table that
has received a **+0.5 continuity correction** — which `es_from_2x2()` itself emits for
a zero cell, and which study 04 supplies — has **half-integer cells**, and nothing in
the inputs reveals it: adding 0.5 to all four cells adds exactly 1 to every margin.
Worked case: the corrected table `0.5/50.5/2.5/48.5` came back as `1/50/2/49`, giving
log RR **−0.693** against the true **−1.609** — an error of **0.916**. At `br = 0.01`,
**75.8%** of study 04's replications carry the correction.

**The obvious fix was measured and rejected.** Rounding to the nearest 0.5 fixes the
corrected case but costs the ordinary one. Over 900 tables of each kind:

| rule | integer tables (OR at 2 dp) | corrected tables |
|---|---|---|
| `round()` to 1 (was) | 96.2% exact | **0% exact, mean err 0.869, only 423/900 solved** |
| round to nearest 0.5 | 91.8% | 100% exact |
| no rounding | 5.0% | 100% exact |
| **prefer integer, take the half only when the exact root sits on one** | **96.1%** | **100% exact** |

The chosen rule costs **0.1 percentage points** on the ordinary case and removes the
failure entirely on the corrected one. The zero-cell guard is relaxed in step: a
half-integer solution is a corrected table, where a cell of 0.5 is the corrected value
of a legitimate zero and its variance is finite, so requiring `>= 1` there would
discard exactly the tables the branch exists to reconstruct. A genuine zero margin
still declines to the enumeration.

**Verified:** `test-or-to-rr-identifiability.R` **995 pass / 0 fail**; archived OR/RR
suites (`test-ES-OR.R`, `test-ES-NNT.R`, `test-rr-to-or-reconstruction.R`,
`test-reverse-ratio-ci.R`) **315 pass / 0 fail / 0 error**. Studies 04 and 05 are
re-run against this fix under 3.5.

⚠️ **This is why 3.5 matters as a check and not only as a chore**: a package defect
that four test suites and two audits missed was exposed by regenerating a simulation
and noticing a number move in the direction nobody predicted.

### 1.6 ☑ NEW — delete the shadowed duplicate reconstruction helpers — DONE (commit `e1b7f59`)

`R/estimate_n_from_es.R:41` and `:121` define `.estimate_n_from_or_and_n_cases` and
`.estimate_n_from_or_and_n_exp` a second time. DESCRIPTION has no `Collate` field, so
alphabetical sourcing lets `internal_multiple_formulas.R` ("i") overwrite
`estimate_n_from_es.R` ("e") — verified by formals. The stale copies now carry a
prominent `DEAD CODE -- DO NOT EDIT THIS COPY` block, but they should be deleted.
Deferred because removing shared-name code deserves its own reviewed change, and
because the same file also holds `.estimate_n_from_irr` and `.estimate_n_from_rr`,
which are **not** duplicated and **are** live — so the file itself must stay.

**Test unit.** Extend `test-or-to-rr-identifiability.R` with an assertion that each
helper is defined exactly once across `R/`.

### 1.1-old (superseded, kept for the record)

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

### 1.5 ☑ NEW — `data_extraction_sheet()` hands users a column name nothing reads — DONE (commit `e951cfc`)

> Landed as the correct name only, no deprecated alias (owner's call): files already
> saved under `all_info_expected` keep being ignored exactly as today, so no existing
> result changes. Flagged in NEWS.md so users can rename. 121 assertions; guard
> confirmed to have teeth (1 orphan pre-fix, 0 post-fix). The two `app/` copies are
> stale — `app.R` calls `metaConvert::` from the installed package and inherits the fix.

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

### 1.3 ☑ `table_2x2_to_cor` is inert in `convert_df()` and silently accepts anything — DONE (commit `2b6f141`)

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

### 2.1 ☑ ANCOVA `@param cov_outcome_r` states the opposite of what the route does — DONE (commit `554c2b3`)

> Split confirmed empirically across all 16 routes taking `cov_outcome_r` (r = 0 vs 0.9,
> n = 60/60): **14** give d ratio 0.4359 / SE ratio 0.4359 / **t ratio 1.0000** (wording
> correct), **2** give d ratio 1.0000 / SE ratio 0.4823 / **t ratio 2.0733** (wording
> inverted). Of the 15 copies of the shared wording, exactly **one** sat on a marginal-SD
> route (`es_from_ancova_means_sd_pooled_crude`); `es_from_cohen_d_adj` has separate text
> that was silent rather than wrong and gained the same warning. 40 assertions pinning
> behaviour, not prose. Full suite 1639 pass / 0 fail.

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

### 2.6 ☐ NEW — `or_to_cor`: the consistency question raised by 1.3

The phi decision (1.3) rests on margin dependence. That principle applies to three of
the four `or_to_cor` options too, so the position needs stating rather than leaving
implicit. Measured from `09a_or_to_cor_CONT_nrep1000.csv`:

| option | mean \|bias\| vs **own** | vs **population tetrachoric** | worst coverage | margin drift* |
|---|---|---|---|---|
| `2x2_tetrachoric` (reference) | 0.0163 | 0.0163 | 0.849 | **0.0109** |
| **`bonett`** (default) | 0.0185 | 0.0185 | **0.940** | **0.0212** |
| `digby` | 0.0201 | 0.0201 | 0.906 | 0.1205 |
| `pearson` | 0.0223 | 0.0223 | 0.902 | 0.1245 |
| `lipsey_cooper` | 0.0231 | **0.0870** | **0.000** | 0.1127 |

\* range of bias across the margin grid at ρ = 0.5, n = 300 — i.e. how much the answer
moves when only the margins change, at a fixed true correlation.

**Three conclusions.**

1. **The default is right and the menu is defensible.** All five compute their own
   estimand correctly (|bias| 0.016–0.023). `bonett` is the least margin-dependent
   option and has the best worst-case coverage, so the menu is "best by default,
   alternatives available" — unlike a phi option, which would have had no dominating
   sibling. This is what distinguishes 1.3's refusal from keeping these.
2. **`pearson` and `digby` are ~6× more margin-dependent than the default** (0.12 vs
   0.02). That is the same pathology phi was refused for, differing in degree rather
   than in kind, and it is currently undocumented. They need the warning phi is being
   denied for.
3. **`lipsey_cooper` is different in kind and needs more than a warning.** It targets
   the point-biserial via the Cox transform, not the tetrachoric: against the
   population tetrachoric its mean |bias| is 4.7× the default's and its worst-case
   coverage is **0.000**. A review whose rows resolve to a mix of `lipsey_cooper` and
   the other three pools two different estimands in one column — the same defect as
   the `es_from_phi` swap fixed in 1.3, at review scale.

**Strategy.** Document the drift on `pearson`/`digby`; document `lipsey_cooper` as a
different estimand and consider a cross-row flag when a pool mixes it with the
tetrachoric-targeting options. Do **not** remove them — they are legacy-compatible and
correctly implemented. Folds together with 2.3 below, which aligns the entry-point
defaults.

**Test unit.** Extend `test-method-defaults-consistency.R` (2.3) with an assertion that
every `or_to_cor` option is reachable and documented, plus a flag test for the mixed
`lipsey_cooper` pool if that check is added.

### 2.3 ☑ `or_to_cor` default diverges: `es_from_or_se()` alone ships `"pearson"` — DONE

> **The scan found 16 divergences, not 1.** Besides the three `or_to_cor` = `"pearson"`
> functions (`es_from_or_se`, `es_from_user_crude`, `es_from_user_adj`), **13**
> `es_from_mean_change_*` / `es_from_paired_*` functions default to
> `pre_post_to_smd = "cooper"` against `convert_df()`'s `"bonett"`.
>
> **Those 13 are legitimate and must not be "fixed".** Verified by calling one:
> `es_from_mean_change_sd(..., pre_post_to_smd = "bonett")` **errors** — bonett needs a
> baseline SD that change-score and paired test-statistic data do not carry. A blanket
> "every default must match `convert_df()`" test would have demanded 13 changes that
> break every direct call to those routes.
>
> ⚠️ **THE ITEM'S PREMISE WAS WRONG, and measuring caught it.** The write-up assumed both
> defaults produce a number differing by 1.30×. On the **minimal documented input**
> `(or, logor_se)` they do not:
>
> | `or_to_cor` | r |
> |---|---|
> | `pearson` (old default) | 0.3463 |
> | `digby` | 0.3307 |
> | `bonett` (new default) | **NA** |
> | `lipsey_cooper` | NA |
>
> So aligning turns `es_from_or_se(or, logor_se)` from 0.3463 into **NA**. That is a real
> capability loss for direct callers — **and it is exactly what `convert_df()` already
> returns for the identical row**, which is the whole point of the item. Landed as-is
> (owner's call), pinned by a test, and documented in NEWS with the escape hatch
> (`or_to_cor = "pearson"` explicitly). On a margin-complete table the change is
> 0.3463 → 0.3205.
>
> **The test derives its allowlist instead of hardcoding one.** A divergence on
> `pre_post_to_smd` is justified exactly when that function's own
> `.validate_pre_post_to_smd(allowed_methods = ...)` excludes `convert_df()`'s default.
> Self-maintaining, and it cannot hide real drift. **Delegation had to be followed
> transitively**: `es_from_paired_f`, `es_from_paired_t_pval` and `es_from_paired_f_pval`
> declare no `allowed_methods` at all — they hand off (`_f_pval` → `_t_pval` →
> `es_from_paired_t`), and a one-level search called three legitimate defaults "drift".
> Their refusal is also pinned behaviourally.
>
> **Verified:** `tests/testthat` **2921 pass / 0 fail / 0 error / 0 skip** (2817 + 104
> new — no existing assertion changed outcome) and `tests_save/checked` **7193 / 0 / 0 /
> 0**, exactly the baseline. Guard has teeth: reverting the alignment gives **16
> failures**.

### 2.7 ☐ NEW — the `or_to_cor` fallback picks the estimand-mismatched method

Found while measuring 2.3. When a row lacks bonett's inputs, `es_from_or_se()` falls back
to **`lipsey_cooper`** ([es_from_stand_OR.R:473-485](../R/es_from_stand_OR.R#L473-L485)) —
which by item 2.6's own table targets the **point-biserial via the Cox transform, not the
tetrachoric**, with worst-case coverage **0.000** against the population tetrachoric.
Meanwhile `pearson` and `digby` *do* target the tetrachoric and need only
`(or, logor_se)` — they were available and were not chosen.

Measured: with `n_sample` present but no margins, bonett→fallback gives r = 0.2449 where
`pearson` gives 0.3463 and true-bonett (margins supplied) gives 0.3205. So the fallback
moves the answer **further** from the method actually requested than the estimand-matched
alternatives would.

This is **pre-existing and affects `convert_df()` too**, not something 2.3 introduced —
which is why it was not folded in: changing the fallback to `pearson` would move results
for every margin-poor row in the package and needs its own measured decision.

**Test unit.** Extend `test-method-defaults-consistency.R` with the fallback chain's
target, once the chain is decided.

### 2.8 ☐ NEW — the `or_to_cor` fallback message tells users to supply something that does not work

The message says *"For 'bonett', supply 'small_margin_prop', **or** 'n_sample' together
with one of ('n_exp','n_nexp') and one of ('n_cases','n_controls')"*
([es_from_stand_OR.R:483-484](../R/es_from_stand_OR.R#L483-L484)). The eligibility test at
[:454-456](../R/es_from_stand_OR.R#L454) is a **conjunction**: `small_margin_prop` **and**
`n_sample` **and** `n_exp` **and** `n_cases` must all be non-NA.

Measured: `small_margin_prop = 0.3` alone → NA; `small_margin_prop + n_sample` → **still
falls back**; only `n_sample + n_cases + n_exp` reaches bonett (which then *derives*
`small_margin_prop` itself at [:436-439](../R/es_from_stand_OR.R#L436)). So the one input
the message names first is the one that cannot work on its own.

Either the message or the condition is wrong; deciding which is the item.

### 2.9 ☐ NEW — a delegating route reports a failure in the function it delegated to

`es_from_paired_f(pre_post_to_smd = "bonett")` errors with *"Invalid 'pre_post_to_smd'
argument in **es_from_paired_t()**"* — naming a function the user never called. Same for
`es_from_paired_t_pval()` and `es_from_paired_f_pval()`. Cosmetic, but it sends a user
looking in the wrong place. `.validate_pre_post_to_smd()` already takes `func_name`; the
delegating routes just pass the callee's.

### 2.4 ☑ `unit_type`: documented values exclude the shipped default; no validation — DONE

> Folded into 2.3 (owner's call): same defect class, shared test file.
>
> **Confirmed before fixing:** `convert_df(x, measure = "g", unit_type = "banana")` ran
> without error. The argument is read in exactly one place, as `unit_type == "sd"`, so
> every other value — a typo included — silently selected raw units. The docs were
> inconsistent with each other as well as with the code: **six** `@param` blocks said
> `"sd"` or `"raw_scale"`, **two** said `"sd"` or `"value"`, and the shipped default
> `"raw_scale"` appeared in neither of those two.
>
> **Landed:** tolerated set `c("sd", "raw_scale", "value")` — `"value"` kept as a
> documented synonym of `"raw_scale"` because the manual promised it and it already
> behaved that way (everything non-`"sd"` means raw units). All **8** `@param` blocks
> rewritten to one description that names the default. Validation added to `convert_df()`
> and all **7** exported takers.
>
> ⚠️ **A test caught a design error.** Warning on a bad **column** while leaving the value
> in place was not enough: `convert_df()` hands the column straight to the exported
> `es_from_*()` routes, whose own argument check then **stops**. So the run aborted
> anyway, two calls later, with a message naming the wrong layer. Bad column cells are now
> **neutralised to NA** after the warning — NA is already the internal "not supplied"
> value and means raw units, exactly what the warning says those rows will get.
>
> ⚠️⚠️ **THE ARCHIVED SUITE CAUGHT A REGRESSION THAT `tests/testthat` DID NOT — the third
> time in this project.** The main suite was **2919 / 0 / 0 / 0** while
> `tests_save/checked` was **7188 pass / 2 fail / 1 error** — and note the assertion
> total had *fallen* by 5 from the 7193 baseline, the exact signature this file's header
> warns about: the single **error** in `test-ES-COR.R` abandoned the rest of its block.
>
> Cause: the tolerated set had been built **from the documentation**, and the archived
> tests use `unit_type = "raw_data"` — a **fourth** spelling that appears in no `@param`
> block and worked only because everything that is not `"sd"` means raw units. Rebuilt
> from a census of the whole repository: `raw_scale` (67 uses), `sd` (29), `value` (8),
> `raw_data` (4). All four are accepted; a new test re-runs that census so a fifth
> spelling entering the sources fails there rather than in the archived suite.
>
> The guard keeps its value regardless: what it catches is a value **intended** as `"sd"`
> but not spelled `"sd"` (e.g. `"SD"`), which silently produced raw units and is the only
> direction in which this argument can go wrong unnoticed.
>
> After the fix: `tests/testthat` **2921 / 0 / 0 / 0**, `tests_save/checked`
> **7193 / 0 / 0 / 0** — both totals exact against baseline.

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

### 2.4-old (the original write-up, kept for the record — see the ☑ entry above)

**Symptom.** `?convert_df` says "must be either 'sd' or 'value'"; the default is
`"raw_scale"`; the only read is `ifelse(unit_type == "sd", ...)`. A typo silently
selects raw units.

**Strategy.** Make the default one of the documented values and validate the
argument. Arithmetic is correct — this is docs + a guard.

> Landed differently from this plan in one respect: the default was **kept** at
> `"raw_scale"` and the *documentation* corrected to match it, rather than the default
> being changed to `"value"`. `unit_type` is also a data **column**, so changing the
> shipped default would have silently altered the meaning of existing extraction files
> for no behavioural gain — `"value"` and `"raw_scale"` are the same thing to the one
> line that reads the argument.

**Test unit.** Covered by `test-method-defaults-consistency.R` (2.3) plus an
`expect_error` on an invalid value.

---

### 2.5 ☑ `measure = "z"` returns two different transforms depending on route — DONE

> ⚠️ **The item understated the scope, and the owner's challenge is what exposed it.**
> The write-up frames this as a consequence of *choosing* `smd_to_cor = "lipsey_cooper"`
> over `"viechtbauer"`. Measured, the split is **by input family, at package defaults,
> with no user choice involved**:
>
> | family | what lands in `z` |
> |---|---|
> | `pearson_r`, `fisher_z`, `or_se`, `2x2` | Fisher's z = `atanh(r)` |
> | `means_sd`, `cohen_d`, `student_t` | a **variance-stabilising** transform |
>
> At a point-biserial ρ = 0.75 the SMD route returns z = **1.0925** where `atanh()` of
> its own r is **1.7468**. One ordinary `convert_df(measure = "z")` call over a mixed
> sheet produced 0.564 (VST) beside 0.549 / 0.332 / 0.199 (Fisher) — the first row sits
> 0.174 from where the others' scale puts it.
>
> **"Why only z — isn't this the same for pre/post, mean-change, ANCOVA, paired F?"**
> (owner). Answered by measurement, and the answer changed the framing:
>
> | family | can two scales share one column? | already flagged? |
> |---|---|---|
> | mean-change vs endpoint | yes | **E6 fires** ✔ measured |
> | paired-F vs endpoint | yes | **E6 fires** ✔ measured |
> | pre/post means vs endpoint | yes | E6's method list covers it |
> | **ANCOVA vs endpoint** | **no** — ANCOVA lands in the *adjusted* columns, endpoint in the *crude* | nothing to flag |
> | **z** | yes | **nothing** ✔ measured |
>
> So it was never "only z": the package already implements this principle for the SMD
> family (E6/E7) and separates ANCOVA structurally. **z was the one gap**, and filling it
> restores consistency rather than singling z out.
>
> **Landed as E8**, deliberately mirroring E6 — same `[INFO]` severity (a pool property;
> every row is individually correct), same `cmp_rep` selected-route logic, same
> `group_key` scoping, same `enable_cross_row` gate, same message shape.
>
> **Classification is derived from the output, never from a route list.**
> `.z_transform_by_route()` asks each method frame whether its own `z` equals `atanh()`
> of its own `r`, and records the verdict per `info_used` in `attr(res, "z_transform")`
> (threaded like `r_defaulted` / `smd_denom_used`). A hardcoded map is exactly what
> rotted in item 1.2. It also means the marker tracks `smd_to_cor` automatically —
> verified: `means_sd` is `"vst"` under the default and `"fisher"` under
> `lipsey_cooper`.
>
> **Kept distinct from the estimand question**, in code comments, docs and tests:
> `viechtbauer` estimates the *biserial* and `lipsey_cooper` the *point-biserial*, which
> is a documented choice study 01 already measures. E8 is only about the `z` column not
> being `atanh()` of the `r` column on every route.
>
> Documented on `convert_df`'s `@param smd_to_cor` and in the `es_from_means_*` details;
> E8 added to the flag catalogue in `CLAUDE.md`; NEWS entry.
>
> **Verified:** `tests/testthat` **2964 pass / 0 fail / 0 error / 0 skip** (2921 + 43 new)
> and `tests_save/checked` **7193 / 0 / 0 / 0**, exactly the baseline. The remedy the flag
> recommends is itself tested: `smd_to_cor = "lipsey_cooper"` silences E8 and puts every
> route on Fisher's z.

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

### 3.1 ☑ Study 05's `vanderweele_rr_squared` computes 2·log(OR), not 2·log(RR) — DONE

> **Bug confirmed by running it, not by reading it.** At the point the candidate reads
> `theta_sample`, `all.equal(theta_sample, log(or))` is TRUE and `all.equal(theta_sample,
> log(rr))` is FALSE. The signature bias reproduces exactly: at (br = 0.30, rr = 2,
> n = 300, 20 000 reps) mean bias **1.2601** against mean |log OR| **1.2601** — the same
> number twice, because bias = 2·logOR − logOR = logOR. Corrected, that cell's bias is
> 0.1346.
>
> **Fix.** The sample log RR now has its own column, `logrr_sample`, set in
> `gen_2x2_or()` *before* the target overwrite and never reused; the candidate reads it.
> A `stopifnot()` beside the overwrite pins both scales. **Deliberately NOT asserted:
> that the two columns differ** — a replication with `a/n_exp == c/n_nexp` gives
> log OR = log RR = 0 legitimately, so that check could abort a valid run at small `nrep`.
>
> **Scope extension (owner's call): the mirrored CI, item 1.9's counterpart.** The
> candidate also built the symmetric interval 1.9 just removed from study 04. It now uses
> the transformed RR limits widened by the bias-ratio bound. **The factor is derived, not
> cited** — VanderWeele writes only the OR→RR direction — **but it is exact**, and the
> derivation is one line: `RR²/OR = p₁(1−p₁)/(p₀(1−p₀))`, so over `[0.5−v, 0.5+v]` the
> bound is `1/(1−4v²)`, exactly the **square** of the OR→RR factor `1/sqrt(1−4v²)`.
> At [0.2, 0.8] that is **1.25² = 1.5625**. Verified numerically over a 0.0005 grid:
> measured max 1.562499, attained at p₀ = 0.8, p₁ = 0.5; bands [0.1,0.9] / [0.3,0.7] /
> [0.4,0.6] all match `1/(1−4v²)` to 6 dp. Measured coverage of the true log OR, in-band
> cells: **naive 0.874–0.958 → widened 0.994–1.000**, width roughly doubling.
>
> **Student task marked, not rewritten.** `students/tasks/task-05-rr-to-or.md` had an
> explanation invented for the artifact ("built for outcomes happening about half the
> time — the true baseline risk here is 0.15"). That is not the cause, and the paragraph
> plus the `vanderweele_rr_squared` table row now carry a dated correction saying so.
> The full rewrite waits for 3.5, since the numbers themselves change on regeneration.
> ⚠️ `students/` is **gitignored**, so that edit is on disk only and is not in any commit.
>
> **Verified:** `simulations/tests` 78 pass / 0 fail / 0 error / 0 skip (38 from 1.9 + 40
> new). No package file touched, so the two package suites stand from the 1.8 commit.
>
> ⚠️ **Process note.** The first run of the new end-to-end test **errored** (`could not
> find function "es_from_rr_se"` — the runner sources the study files but had not loaded
> the package), which is precisely the "erroring block reports failed = 0 while
> abandoning its remaining assertions" trap in this file's header. `run_tests.R` now
> calls `load_metaconvert()` and exposes `METACONVERT_AVAILABLE` so such tests **skip**
> rather than error if it is unavailable.

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

### 3.2 ☑ Study 03's `theta_own` is scale-matched, not estimand-matched — DONE

> **Confirmed in the shipped files, not inferred:** for every `(r)` route in both
> `03a_CAT` and `03b_CONT`, the `own` bias column is **bit-identical** to `population`
> — 200 cells each, max |difference| exactly **0**. The `(z)` routes differ by up to
> 0.223, which is purely the `atanh` scale correction. So the column that exists to
> remove estimand mismatch was removing nothing for half the rows.
>
> **Size of the mis-scoring it left in place**, over study 03's own grid:
>
> | direction | mean \|gap\| | max |
> |---|---|---|
> | CONT: phi scored against ρ | 0.0925 | 0.2429 |
> | CAT: tetrachoric scored against phi | 0.1341 | 0.4186 |
>
> Worst condition (CAT, ρ = 0.5, p_exp = 0.3, p_case = 0.1): the tetrachoric's true
> estimand is **0.9186** while it was scored against **0.50**.
>
> **Fixed** by computing the genuine per-method estimand, as study 09 already did. In
> each mechanism exactly one method's estimand equals ρ — CONT: tetrachoric = ρ, phi =
> the dichotomised table's phi; CAT: phi = ρ, tetrachoric = the latent ρ implied — which
> is precisely why one shared `theta_pop` column cannot serve both.
>
> **The mechanism is stamped on the data (`dat$dgm`), not captured in a closure**, so it
> survives serialisation to a parallel worker regardless of how `run_study()` exports its
> arguments; `estimate_2x2()` **errors** if the stamp is missing rather than silently
> reverting to scale-only matching.
>
> **The roadmap's `atanh`-vs-VST warning checked and cleared:** `es_from_2x2()`'s `z` is
> Fisher's z exactly (`max |z − atanh(r)| = 0` over 20 random tables) and
> `phi_candidate()` builds `z <- atanh(r)`, so plain `atanh` is right *here*. A test pins
> it, so if the package ever changes that transform this fails instead of the targets
> going quietly wrong.
>
> **Duplicate helpers removed (the simulation-side 1.6).** Study 03's `phi_to_p11` /
> `phi_attainable` were functionally identical to study 09's `.phi_to_p11` /
> `.phi_attainable` — verified on 500 random inputs each. All five latent-2×2 helpers now
> live once in **`simulations/R/05_latent_2x2.R`**, bodies moved **byte-identical**
> (`deparse()`-compared against `git show HEAD`), so study 09's numbers cannot move.
>
> **Verified:** `simulations/tests` **367 pass / 0 fail / 0 error / 0 skip** (204 + 163
> new). Guard has teeth: reverting study 03 alone gives **11 failures**. The generic
> anti-regression the item asks for runs over all seven studies exposing a `theta_own`
> and reports which separate `own` from `population` — currently **all seven**. No
> package file touched.
>
> ⚠️ Study 03's aggregates are now stale for a second reason (the first being the
> `theta_own` change this item fixes). Folds into 3.5, which also needs study 09
> re-run only if one wants the two-column schema from 3.3.

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

### 3.3 ☑ `performance()` publishes cells built from 3–10 replications while advertising `n_valid` 380–996 — DONE

> **Census run over the raw files before fixing** (19,129 cells, grouping keys taken
> from each study's own aggregate). The earlier write-up's five cells were the tip:
>
> | statistic published from | cells | share |
> |---|---|---|
> | < 100 replications | 54 | 0.28% |
> | < 30 | 30 | 0.16% |
> | < 10 | 4 | 0.02% |
> | < 5 | 2 | 0.01% |
>
> Worst cell: `n_valid` 440, coverage from **3**. 124 cells have a CI count below half
> their `n_valid`, 44 below a tenth. **Not confined to `grant`** (40 cells): study 03a's
> `phi (r)` / `phi (z)` account for 14 more, which the earlier write-up missed.
> 798 cells already reported `coverage = NA`, so the *fully* empty case was handled — it
> is the partial ones that looked populated.
>
> **Landed:** `n_valid_se` and `n_valid_ci` returned beside `n_valid` (computed before
> the early return, so they appear on the give-up path too — otherwise that row could not
> be `rbind`ed beside a normal one). Statistics withheld below `SIM_DEFAULTS$min_valid`
> (**30**, owner's call), **per family, keyed on that family's own count**: a cell with
> 1000 usable point estimates and 3 usable intervals keeps its bias and loses its
> coverage. Counts always reported, so suppression is explainable and a reader can apply
> their own floor.
>
> **No floor on the point-estimate family, and that is measured rather than assumed:**
> `n_valid` has always been printed beside `bias`, and **0 of 57,026** shipped cells
> report a bias from fewer than 100 valid estimates. A rule there would govern nothing.
>
> **Consumers checked, not assumed safe.** `summarise_raw()` is the only caller. The app
> reads one CSV at a time (no cross-file `rbind`, so mixed schemas during the 3.5
> rollout are harmless) and derives condition columns as "everything before `method`" —
> the new columns land after it, so they cannot be mistaken for grid factors.
> README's schema paragraph updated.
>
> ⚠️ **A test caught a real design slip**: `min_valid = Inf` disabled the floor instead of
> withholding everything, because the guard was `is.finite(min_valid) && > 0`. Now
> `!is.na(min_valid) && > 0`. Only NA or a non-positive value disables it.
>
> **Verified:** `simulations/tests` **204 pass / 0 fail / 0 error / 0 skip** (128 + 76
> new). No package file touched.

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

### 3.4 ☑ `run_everything()` aborts on `run_99` — DONE

> **Reproduced before fixing**: the registry `ls(pattern = "^run_[0-9]{2}$")` returns
> `run_01 … run_09, run_99`; `run_99`'s formals were empty, and
> `do.call("run_99", list(nrep = 1000, cores = 1))` raised *"arguments inutilisés
> (nrep = 1000, cores = 1)"*.
>
> **Fixed by giving `run_99` the ignored signature, not by excluding it from the
> pattern.** Excluding it would create a second class of runner that a future `run_98`
> would have to know about; the ignored signature keeps **one** contract — every
> registered runner accepts `(nrep, cores)` — which is now asserted. `run_99` already
> calls `load_metaconvert()` itself, so it is safe in sequence.
>
> **The real cost was never lost results** (it sorts last): it is that a full run always
> ended in an error, which teaches whoever runs it that the error at the end is normal.
>
> `simulations/tests/test-runner-registry.R` (50 assertions) checks every registry
> member's arity via `match.call()` — arity only, never executing the studies — plus
> that `run_99`'s new arguments are **inert** (absent from its body, so a later edit
> cannot quietly make them load-bearing), that `run_all.R` still uses the pattern this
> file tests, and that `run_99` survives the exact `do.call()` end to end, since it is
> deterministic and fast enough to actually run.
>
> **Guard confirmed to have teeth**: with the fix reverted, 123 pass / **5 fail**; with
> it, **128 pass / 0 fail / 0 error / 0 skip**. No package file touched.

**Symptom.** `run_99` takes no arguments; `run_everything()` calls every
`run_[0-9]{2}` with `nrep`/`cores`. It sorts last, so nothing is lost — but the
run ends in an error.

**Strategy.** Give `run_99` the ignored signature, or exclude it from the pattern.

**Test unit.** `simulations/tests/test-runner-registry.R` — assert every
registered runner accepts `(nrep, cores)`.

---

### 3.5 ☑ Regenerate stale aggregates — DONE

All 17 aggregates regenerated at the shipped **nrep = 1000** (owner's call: keeping the
replication count fixed means every before/after difference is a change of CODE, not of
Monte Carlo error), plus the two study-03 `nrep = 300` variants so they were not left
orphaned. Full pass 120.8 min on 26 cores, no errors; studies 04 and 05 then re-run a
second time (83.8 min) against the item-1.10 rounding fix, in an isolated git worktree
so concurrent editing of the main checkout could not reach the 26 parallel workers
mid-run.

**Every predicted change landed, and one prediction of mine was wrong.**

*Study 05 (item 3.1) — the artifact is gone.* The descending baseline-risk gradient the
student task had invented an explanation for **was the bug**:

| br | coverage before | after |
|---|---|---|
| 0.01 | 1.000 | 1.000 |
| 0.15 | 0.870 | 0.988 |
| 0.30 | 0.658 | 0.985 |
| **0.50** | **0.509** | **0.987** |

Mean \|bias\| 0.687 → 0.440, against the ~0.443 predicted. `grant` moved 0.222 → 0.476,
which is item 1.4 working: it now returns NA instead of an estimate with a NaN variance,
so the surviving cells are a different set.

*Study 04 (items 1.1 / 1.8 / 1.9 / 1.10) — the headline reverses.*

| method | shipped | after 1.1/1.8 | after 1.10 |
|---|---|---|---|
| `metaumbrella_exp` | 0.36493 | 0.05785 | **0.00000** |
| `metaumbrella_cases` | 0.00596 | 0.03110 | **0.00000** |
| every other method | — | unchanged | **bit-identical** |

Coverage for `metaumbrella_exp`: **0.505 → 1.000**. The two variants are now *identical
to each other*, because both solve the same table instead of searching for it — so the
README's "two `metaumbrella` variants with very different reliability" no longer
describes the package. VanderWeele's √OR keeps its point estimates exactly (0.20998
throughout) with in-band coverage **0.931 → 0.996** (width 0.853 → 1.299), which is
item 1.9 doing precisely what it was scoped to do.

*Study 03 (item 3.2).* `own` is no longer a copy of `population` in any of the four
files (max \|diff\| 0.4186 CAT, 0.2429 CONT — the exact estimand gaps measured when the
item was built).

*The item-3.3 floor withheld **64** cells*, not the ~30 predicted — 52 in study 05, 12
in study 03a. Explainable: item 1.4 lowered `n_valid_ci` in the grant cells, pushing
more of them below 30. (A first count of 1,394 was wrong: it included cells whose
coverage was *always* NA, which is what study 06 is about.)

⚠️ **My "study 09 is a control that must not move" prediction was wrong, and this file
had already said so** — item 3.5's own text records that study 09's aggregates predate
package commits. Cause pinned: commit `4955631` (13 Aug) versus aggregates from 12 Aug,
which changed rows lacking bonett's inputs from *blanked to NA* to *keep the
lipsey_cooper value*. Effect: `n_valid` **964 → 1000** in 138 cells, none lost, all at
low event rates. Verified rather than inferred: with a degenerate margin,
`or_to_cor = "bonett"` returns r = 0.30533, bit-identical to `lipsey_cooper`, and fires
the fallback message. **So study 09's `bonett` row is now a mixture of two estimands**
(~3.6% of its replications), which is finding **2.7** appearing in the results, and is
why its mean \|bias\| rose 0.0498 → 0.0546 — the old figure was conditional on not having
failed.

**Test unit landed**, and keyed on **git commit times, not file mtimes**: the roadmap's
literal wording (`newer than the study file`) cannot work in a git repository, because a
clone stamps every file with the checkout time, so an mtime test would pass everywhere
and detect nothing. Package dependencies are **derived** by tracing each study's
`es_from_*` calls to their defining files, rather than hand-listed — depending on the
whole of `R/` instead would mark every aggregate stale on any package commit and the
test would be ignored. Teeth confirmed: dirtying `studies/07_ancova_to_smd.R` makes it
name `07a`/`07b` as stale; reverting restores **407 pass / 0 fail / 0 error**.

**Also updated:** the README's headline figures (the "coverage 0.497, worst shipped
option" claim, the no-go map table, and the "study 03's files are stale" note), and the
student task's `vanderweele_rr_squared` narrative — rewritten to teach the two lessons
the episode actually carries (a number that looks odd deserves a check, not a story; an
artifact often tracks the design rather than the method) rather than quietly deleted.
⚠️ `students/` is **gitignored**, so that rewrite is on disk only and appears in no
commit.

⚠️ **This item earned its keep as a check, not a chore**: regenerating the simulations
exposed a package defect (item 1.10) that four test suites and two audits had missed.

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
