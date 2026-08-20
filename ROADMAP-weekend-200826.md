# Reliability workstream — roadmap, 2026-08-20

Scope: the **reliability family** — Cronbach's alpha, McDonald's omega, ICC. Sister
document to `simulations/ROADMAP-WEEKEND-200826.md`, which covers the OR/RR
conversion + simulations stream. The two share a branch (`audit-remediation`) and
nothing else; no item here touches a file that one touches.

Status legend: ☐ not started · ◐ in progress · ☑ done

**Verification rule** (inherited from `simulations/ROADMAP.md`, and it earned its keep
again this week): run `tests_save/checked/` as well as `tests/testthat/`, and check
testthat's **`error`** column, not just `failed`. Compare the assertion TOTAL against
the previous run — a count that falls with no failures means assertions stopped
executing.

**Second verification rule, learned this week.** `devtools::check_man()` does **not**
tell you a help page rendered. It validates `\usage` against the formals. It reported
"No issues detected" for `man/es_from_omega.Rd` while that file's `\details{}` block
was **completely empty** — ~130 lines silently dropped. Verify documentation by
reading the generated `.Rd`, and read roxygen's own output: its error line
(`✖ ... has mismatched braces or quotes`) contains neither "Writing" nor "warn", so a
grep for those misses it.

## Baselines as of this document

| suite | command | count |
|:--|:--|--:|
| package, main | `devtools::test()` | **3210** pass / 0 fail / 0 error |
| package, archived | `tests_save/checked/` (~25 min) | **5613** pass / 0 fail / 0 error |
| `devtools::check_man()` | | clean |

Reliability commits this stream: `4f9eb4d` · `515ce68` · `cc945e4` · `6638146` ·
`1a521c2` · `17c5813` · `7e385fb`. Nothing pushed.

---

# Part 1 — ICC audit against `data-raw/icc`

Eight papers read in full: Bhat & Beretvas (2026, *Res Synth Methods*) — the direct
methods paper — plus Gnambs (2014), Polit (2014), Scharfen et al. (2018), Xie et al.
(2026), Elliott et al. (2020), the pediatric-HRV retest paper, and Child Abuse &
Neglect (2024). 23 findings confirmed by reproduction, 3 refuted.

## 1.0 ☑ The good news, established first

**metaConvert's ICC sampling variance is the one Bhat & Beretvas (2026) recommend.**
Their Eq (2.9)–(2.10) Fisher-TF variance, delta-mapped onto the `ln(1 − ICC)` scale,
is *algebraically identical* to `R/es_from_ICC.R:152-154`:

```
V(T) = [4(1+(k−1)ρ)²/k²] · 0.5·k/((k−1)(n−1)) = 2(1+(k−1)ρ)²/(k(k−1)(n−1))
n=50, k=2, ICC=0.80 → delta-of-FisherTF 0.257142857143
                      metaConvert       0.257142857143   diff −5.6e−17
```

And the formula is **nominal for the form it claims**: ICC(3,1) coverage matched
across all 8 tested cells at 20,000 reps (e.g. n=50, k=2, ρ=0.889 → coverage 0.946).
A claim that it was wrong was raised and **refuted**.

So the arithmetic is right. Everything below is about *which estimand* it is right
for, and about what the package lets a user do by accident.

## 1.1 ☐ BLOCKER — the default `icc_type = "agreement"` is far worse than documented, and degrades with n

`agreement` = ICC(2,1). The shared SE is a one-way approximation that assumes
negligible between-rater variance. `R/es_from_ICC.R:35-36` says coverage is "around
0.74–0.76". Measured (two-way DGP, ICC(2,1)=.80, k=2, rater variance 50% of the error
budget, 15,000 reps/cell, Shrout & Fleiss estimator):

| n | empirical SD | reported SE | SE understated | weight inflated | **coverage** |
|--:|--:|--:|--:|--:|--:|
| 20 | 0.5593 | 0.4130 | 1.35× | 1.8× | 0.832 |
| 50 | 0.4624 | 0.2584 | 1.79× | 3.2× | 0.664 |
| 200 | 0.4179 | 0.1283 | 3.26× | 10.6× | 0.323 |
| 1000 | 0.3997 | 0.0574 | 6.97× | **48.6×** | **0.138** |

The empirical SD barely shrinks (0.559 → 0.400) because ICC(2,1) inherits MSC, which
carries **k−1 = 1 df regardless of n**, while the reported SE shrinks like 1/√n. So
coverage gets *worse* as studies get bigger — the opposite of what a reader assumes
from a bounded "0.74–0.76". Dose–response at n=50, k=2: rater variance 0% → coverage
0.948; 20% → 0.899; 50% → 0.662; 80% → 0.371.

Fisher-TF and `atanh` fail identically here (0.317 / 0.318 at n=200), which confirms
this is an **estimand/df problem, not a transform problem**.

**Do.** (a) Replace the "0.74–0.76" claim with the measured range and state that
coverage degrades with n. (b) Decide between two options, both of which need no new
inputs: flip the default to `consistency` (the form the SE is exact for), or return
`se = NA` for agreement rows — the treatment `es_from_omega()` already gives a bare
omega, leaving the row visible and countable but out of the pool. **Recommend (b)**:
it is honest, and it makes V31 actionable instead of advisory.

**Test unit.** `tests/testthat/test-icc-agreement-se.R` — pin the documented range
against a seeded Monte Carlo, and assert the chosen behaviour for agreement rows.

## 1.2 ☐ BLOCKER — an average-measures ICC is silently computed as single-measures

I previously reported that `icc_type = "average"` errors. **It does not.** It warns,
falls back to `"agreement"` at `R/es_from_ICC.R:112`, and returns a wrong estimand:

```
icc_type='average'           -> ACCEPTED, stored 'agreement', es=-2.3026 se=0.2589
icc_type='ICC(2,k)'          -> ACCEPTED, stored 'agreement', es=-2.3026 se=0.2589
icc_type='icc2k'             -> ACCEPTED, stored 'agreement', es=-2.3026 se=0.2589
```

Cost of entering an ICC(2,k) as if single-measures, n=50
(Spearman–Brown: ρ₁ = ρ_k / (k − (k−1)ρ_k)):

| k | reported ICC_avg | true ICC_single | es wrong | es right | **error** |
|--:|--:|--:|--:|--:|--:|
| 2 | 0.90 | 0.8182 | −2.3026 | −1.7047 | −0.598 |
| 5 | 0.90 | 0.6429 | −2.3026 | −1.0296 | −1.273 |
| 10 | 0.95 | 0.6552 | −2.9957 | −1.0647 | **−1.931** |

The SE ratio stays near 1 (1.05–1.39×), so nothing downstream looks wrong. The
**point estimate** is off by up to 1.9 log units.

This is why the COLIVE extraction sheet carries a non-metaConvert `icc_unit` column:
the package cannot represent the distinction, so it must be captured outside it.

**Do.** Recognise `average` / `average_measures` / `icc2k` / `icc3k` as levels.
Minimum viable: `es = NA` + an `[INVALID]` flag naming Spearman–Brown. Better: step
the value down with ρ₁ = ρ_k/(k−(k−1)ρ_k) using the `n_measurements` already
collected, record it in `info_used` (`icc_avg_sb`) the way `omega_estimator` records
provenance. Either way the fallback at `:112` must stop swallowing unknown levels
into `agreement`.

**Test unit.** `tests/testthat/test-icc-average-measures.R`.

## 1.3 ☐ MAJOR — V31 misses 8 of 10 agreement spellings

`es_from_icc()` normalises `icc_type` (case-fold + alias map + fallback).
`R/internal_flags.R:608-615` uses **bare string equality**. Measured through the full
pipeline:

```
'agreement'          -> resolves 'agreement' | V31 fires: TRUE
'Agreement'          -> resolves 'agreement' | V31 fires: FALSE   <-- silent miss
'ICC(2,1)'           -> resolves 'agreement' | V31 fires: FALSE   <-- silent miss
'absolute agreement' -> resolves 'agreement' | V31 fires: FALSE   <-- silent miss
'two-way random'     -> resolves 'agreement' | V31 fires: FALSE   <-- silent miss
```

So the rows that get the anti-conservative SE are largely the ones that *don't* get
warned about it.

**Do.** Extract the normaliser into `.normalise_icc_type()` and call it from both
sites — the pattern already used for `.normalise_omega_type()` /
`.normalise_omega_estimator()`. Two places independently deciding what "agreement"
means is the same rot documented for the hardcoded route list in roadmap 1.2.

**Test unit.** Assert the two agree over the full alias set including case and
whitespace. Fold into `test-reliability-generalization.R`.

## 1.4 ☐ MAJOR — no way to enter a reported ICC SE or CI

Omega has `omega_se` / `omega_ci_lo` / `omega_ci_up`. ICC has nothing: supplied
`icc_se` / `icc_ci_*` columns are **silently ignored**. Yet the ICC literature reports
intervals routinely — every ICC quote in the COLIVE corpus that carries uncertainty
carries a CI, not an SE, e.g. *"ICC 0.94 (95% CI 0.86–0.98)"*.

This also matters because of 1.1: when the package's own SE is untrustworthy for
agreement-type ICCs, the study's reported interval is the better source.

**Do.** Mirror the omega design exactly: `icc_se` (natural scale, delta-mapped),
else CI read at the **bounds**, else the computed `(n, k)` SE. Register in
`internal_check_data.R`, `data_extraction.R`, `.positive_columns()`, `.ci_triplets()`
and `internal_guidance.R`.

**Test unit.** `tests/testthat/test-icc-reported-uncertainty.R`.

## 1.5 ☐ MAJOR — Bonett-scale weighting is biased downward at low ICC

15 studies, n = 20..120, k = 2, ICC(1,1) drawn from its one-way ANOVA sampling
distribution; the *same* simulated estimates feed all three pools:

| true ρ | metaConvert (Bonett IV) | Fisher-TF (recommended) | unweighted mean |
|--:|--:|--:|--:|
| 0.10 | 0.0768 (**−0.0232**) | 0.0987 (−0.0013) | 0.0965 |
| 0.20 | 0.1800 (−0.0200) | 0.1985 (−0.0015) | 0.1945 |
| 0.40 | 0.3860 (−0.0140) | 0.3983 (−0.0017) | 0.3911 |
| 0.85 | 0.8486 (−0.0014) | 0.8497 (−0.0003) | 0.8448 |

Relative parameter bias against Bhat & Beretvas's own yardstick (|RPB| > 0.05 is
"substantial"): Bonett **−0.216 at ρ=0.10** and **−0.092 at ρ=0.20**; Fisher-TF
+0.000 and +0.002.

The two variances are algebraically identical (§1.0) — the difference is *where the
weight is applied*. Bonett weights on `ln(1−ρ)`, which is strongly ρ-dependent at low
ρ, so imprecise low estimates dominate. High-reliability pools (the PROM case,
ρ ≳ 0.6) are barely affected; **low-ICC pools are the neuroimaging / task-fMRI
regime** (Elliott et al. 2020) and are materially biased.

**Do.** Add `icc_to_es = "fisher_tf"` implementing Eq (2.9)/(2.10) directly, with
`reliability_backtransform(method = "fisher_tf")` doing
ρ = (e^{2z} − 1)/(e^{2z} + m̄₀ − 1). The transform is study-specific (z depends on
that study's own k), so **m₀ must be stored per row** to back-transform — a genuine
design constraint, not a detail. Note the transform is monotone *increasing*, so no
bound swap.

**Test unit.** `tests/testthat/test-icc-fisher-tf.R` — reproduce Eq (2.9)/(2.10) and
the RPB table at a fixed seed.

## 1.6 ☐ MINOR — V11's ICC lower bound ignores k

`R/internal_flags.R:439-440` — the comment states the correct bound and the code uses
the loose one:

```r
# icc bounded below by -1/(k-1) >= -1
list(col = "icc", lo = -1, up = 1, label = "ICC"),
```

An ICC of −0.80 with k = 3 is arithmetically impossible (minimum −0.5), passes
validation, and is flagged "Mathematically possible" — which is false for that row.
The weight consequence is large because the SE carries (1+(k−1)ρ):

```
icc=-0.8, k=3 (IMPOSSIBLE) : se = 0.0451  weight = 492
icc= 0.8, k=3 (ordinary)   : se = 0.1954  weight =  26     ratio 18.8x
```

**Do.** Make the bound k-aware where `n_measurements` is present.
**Test unit.** fold into `test-reliability-generalization.R`.

## 1.7 ☐ MINOR — no cross-row checks for ICC form or measurement count

- No **V37 analogue for `n_measurements`**, though k moves the ICC SE by **1.435×**
  across a realistic mix — versus the 1.039× case V37 was written for.
- No **V38 analogue for `icc_type`**: a pool mixing agreement and consistency ICCs is
  completely silent, even though the arithmetic is byte-identical so nothing else can
  reveal it. Four of the eight supplied papers are cited for exactly this anti-pattern.

**Test unit.** fold into `test-reliability-generalization.R`.

## 1.8 ☐ MINOR — vignette defects specific to ICC

Three, all in `vignettes/Psychometrics.Rmd`:

1. The `I²` and funnel-plot cautions are scoped to **alpha only**; the ICC pool
   behaves identically and gets no warning.
2. The ICC section **hand-rolls the back-transform** that
   `reliability_backtransform()` exists to make safe — while the alpha and omega
   sections use the helper.
3. It **endorses pooling agreement and consistency ICCs together**, against the
   supplied literature.

## 1.9 ☐ MINOR — the retest interval appears nowhere

Not in the extraction sheet, not in `df.psychom`, not in the vignette, not as a flag.
Gnambs (2014) and Scharfen et al. (2018) both model it as the primary moderator, and
Polit (2014) makes interval reporting a core recommendation. A `retest_interval_days`
column plus a vignette paragraph would close it. *(No package code required — it
passes through as an arbitrary moderator today.)*

## 1.10 — Documented limitation, not a defect

**Balanced designs only.** `n_measurements` is a single scalar; there is no
per-cluster size vector, so Bhat & Beretvas's m₀ (Eq 2.3) cannot be computed from
unequal cluster sizes. For test-retest and inter-rater work with a fixed k this is
almost always fine. State it rather than fix it.

---

# Part 2 — alpha / omega queue

## 2.1 ☐ MAJOR — V36 (reliability induction) has a high false-positive rate

Confirmed on **randomly generated** alphas with no induction anywhere: at 3 decimals
with varied n, 99.3% of pools flagged and 27.3% of rows; at 2 decimals with round
sample sizes, 100% of pools and 45% of rows. Plausible alphas occupy a ~0.15-wide
band, so the birthday paradox beats the entropy gate.

It reproduced live in the final dry run: two pairs flagged from `runif(.78, .93)`.

**Do.** Either tighten (require ≥ 4 decimals, or ≥ 3 decimals **and** identical n),
or re-word the message from a finding to a prompt. Document the FP rate either way —
V36 must never be reported as a count of induced values.

## 2.2 ☐ MAJOR — no reported-CI entry route for alpha or ICC

Same gap as §1.4, for alpha. A study reporting *"α = .88, 95% CI [.85, .91]"* with no
item count cannot be entered at all, and `user_es_*` is (correctly) refused.

## 2.3 ☐ MINOR — `es_summary_crude` is scale-blind

On the Bonett scale it prints `alpha = -1.204 [-1.385, -1.023]` — an impossible alpha
— and labels all three scales identically as "alpha". Never paste that string into a
supplement. Either label the scale or back-transform for display.

## 2.4 ☐ MINOR — `aggregate_df()` returns 4 columns

`row_index`, `study_id`, `es`, `se`. Every moderator must be re-attached via
`col_mean` / `col_fact`, or merged back by `study_id`. The arithmetic is right; the
loss is silent. Document it in the vignette.

## 2.5 ☐ MINOR — `flag_options` is silently ignored by `summary()` for Tier-1

Tier-1 options reach the checks only from `convert_df()`. Passed to `summary()` they
are dropped with no warning, as is a misspelled option name (`alpha_maxx`).

## 2.6 ☐ MINOR — `data_extraction_sheet(measure = "omega")` omits `n_items`

So V37's item-count check is unreachable from that template. Use `measure = "all"` or
add the column by hand.

## 2.7 ☐ INFO — `DESCRIPTION` does not mention omega

---

# Part 3 — the structural gap

## 3.1 ☐ Zero simulation coverage of the psychometric family

`simulations/studies/` holds nine studies, all conversion routes; grep for
`alpha|icc|omega` returns nothing. Already logged as `simulations/ROADMAP.md` §5.2.

This week produced several numbers that **live only in commit messages and code
comments** — the ICC(2,1) coverage table (§1.1), the low-ICC RPB table (§1.5), the
omega-vs-Bonett variance band `[0.899, 1.003]`, the V36 false-positive rates (§2.1).
Every one is a simulation that was run once and thrown away.

**Do.** A `studies/10_reliability.R` reproducing §1.1, §1.5 and §2.1 at a fixed seed.
Until then, cite those numbers as *claims from this audit*, not as reproducible
results — and do **not** put them in a paper without regenerating them.

---

# Part 4 — limits to state in the paper, not to fix

1. **`I²` is uninterpretable** on a Bonett reliability pool. The variance is
   coefficient-free, so `I²` tracks sample size: one fixed set of α ∈ [.80, .86] gives
   0% at n = 50 and 97% at n = 5000. Report the prediction interval.
2. **No publication-bias machinery**, correctly — α's expectation does not depend on
   n, so funnel / Egger / trim-and-fill are not meaningful. The RG analogue is
   selective *reporting* of reliability, which no software can supply.
3. **No multilevel / RVE support inside the package.** `rma.mv(~1|study_id/row_id)` +
   `robust()` works and agrees with `aggregate_df()` to ~0.0002, but it is your code.
4. **Omega has no aggregate-data variance.** Not a gap to close — the alpha variances
   assume tau-equivalence, which omega exists to drop. The SE must come from the study.
5. **Tier-1 checks key on which columns exist, not on the requested measure.** V31 was
   scoped this week (`7e385fb`); the general property remains.

---

# Suggested order

**1.3** (30 min, removes a silent miss) → **1.2** (the wrong-estimand blocker) →
**1.1** (the documentation is actively misleading) → **1.4** (unlocks the ICC data
that actually exists) → **1.6**, **1.7** → **2.1** → **1.5** (largest, and only
matters for low-ICC pools) → **3.1**.

**1.1–1.4 are the ones that would change a published number.** Everything in Part 2
is quality-of-life; everything in Part 4 is a sentence in the methods section.
