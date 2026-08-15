# Source verification of metaConvert's conversion formulas

Every conversion method reachable through a `convert_df()` method argument, checked
term by term against its source paper in `../data-raw/`. Verification means one of:
the paper's own worked example reproduced through the package, or the paper's
formula hand-coded in R and compared to package output at ≥8 significant digits.

Method: 5 independent verification passes (one per area), then **every claimed
mismatch attacked by an independent skeptic** instructed to default to REFUTED.
That guard exists because earlier in this work three "the package is wrong" claims
were escalated and all three were false. It earned its keep: **3 of the 7 mismatch
claims were refuted outright, 1 partially.**

Two refutations and the automated synthesis did not run (session limit), so two
claims below are marked **UNVERIFIED** — do not act on them without checking.

---

## 1. Verdict table

| method / quantity | source | status |
|---|---|---|
| `smd_to_cor = "viechtbauer"` | Jacobs & Viechtbauer 2017 *RSM* 8:161-180, eq 17/19 | **verified** — worked example exact (a=0.893, z=0.39, CI (−0.08, 0.82)) |
| `smd_to_cor = "lipsey_cooper"` | Cooper/Borenstein (source not in repo) | z-scale variance is a documented metaConvert original |
| `cor_to_smd = "cooper"` | Mathur & VanderWeele 2020 *Epidemiology* 31(2):e16, eq 1.1 | **verified** exact |
| `cor_to_smd = "mathur"` | same, eq 1.2 (estimate **and** SE) | **verified** exact, 222 assertions |
| `cor_to_smd = "viechtbauer"` | `metafor::transf.rtod` (not J&V 2017 — different direction) | **verified** |
| `or_to_cor = "bonett"` | Bonett & Price 2005 *JEBS* 30(2):213, eq 3/4/9 | **verified** — all 3 worked examples exact; *more* faithful than the paper's own prose (see §4) |
| `or_to_cor = "pearson"` | Bonett 2007 *Am Psychol* 62(3):254, `cos{π/(1+OR^½)}` | **verified** exact |
| `or_to_cor = "digby"` | Digby 1983 *Biometrics* 39:753, `H = (OR^0.75−1)/(OR^0.75+1)` | **verified** exact incl. delta-method SE and his 0.375 constant; CI upgraded to B&P 2007 eq 11 |
| `or_to_cor = "lipsey_cooper"` | Wilson 2001 for the OR→d half; z-variance a metaConvert original | partial |
| `table_2x2_to_cor = "tetrachoric"` | `metafor:::.rtet` ML (verified to ~1e-5 against independent ML fit) | **verified**; CI is contested (§3) |
| `or_to_rr = "grant"` | Grant 2014 *BMJ* 348:f7450 | **verified** — all 351 cells of his Table 1 within the 3-dp print bound (max diff 5.0e-04), both sensitivity ranges digit-for-digit |
| — is `grant` = Zhang & Yu 1998? | **YES** — confirmed | one estimator, not two. Docs mis-attribute (§5) |
| `or_to_rr = "dipietrantonj"` | Di Pietrantonj 2006 *Stat Med* 25:2299, eq 18/19 | **verified** — his worked examples to 5-6 sig digits |
| `rr_to_or = "dipietrantonj"` | same, eq 10/11 | **verified** |
| `or_to_rr`/`rr_to_or` `metaumbrella_*`, `transpose` | **no methods paper** — software-only (Gosling 2023) | code says so honestly |
| VanderWeele √OR | VanderWeele 2020 *Biometrics* 76(3):746 | **minimax result confirmed**; genuinely **absent** from metaConvert |
| `.se_from_or()` (OR-SE imputation) | **no source** — metaConvert/metaumbrella invention | see §6, the strongest finding |
| ANCOVA marginal back-transform `σ_res/√(1−r²)` | Lai & Kelley 2012 *BJMSP* 65(2):350, eq 3 | **verified** |
| ANCOVA `n_cov_ancova` df (N−J−q) | Lai & Kelley Table 1 | **verified** — reproduces `lm()$df.residual` exactly at q=1 |
| `pre_post_to_smd = "bonett"` | Bonett 2008 *Psych Methods* 13(2):99 | **verified** — bit-exact with `metafor SMCRH` (max diff 8.9e-16 over a 144-cell grid) |
| `pre_post_to_smd = "morris_dz"` | Viechtbauer 2007 Table 1 eq 31 / `metafor SMCC` | **verified** — diff 0 |
| `pre_post_to_smd = "morris_drm"` (= `"cooper"`) | Caldwell & Vigotsky 2020 eq 13; var = Morris 2000 eq 13 | **verified** — diff 0 |
| `pre_post_to_smd = "morris_dav"` | Bonett 2008 eq 10 / `metafor SMCRPH` | **verified** — diff 2.2e-16 |
| `nu = 2(n−1)/(1+r²)` for d_av's J | Cousineau 2020 *TQMP* 16(4):418, eq 2 | **verified** |
| two-group `pool_sd = FALSE / TRUE` | Morris 2008 d_ppc1 / d_ppc2 / d_ppc3 | **verified** — all 5 studies of his Table 1 |

**No arithmetic error was found in any point estimate or sampling variance that
has a source.** Everything below is either a CI-construction choice, a
reachability/validation gap, a documentation error, or a sourceless method.

> **Superseded (NEWS.md 2.0.1).** That sentence over-reached: it summarised only the
> routes in the table above, and `es_from_etasq()` was never among them. Its crude
> point estimate `d = 2·sqrt(η²/(1−η²))` was the equal-*n*, large-*n* limit and was
> **39% too small at n = 10/90** — an arithmetic error in a point estimate. It now
> inverts `η² = F/(F + N − 2)` exactly and agrees with `es_from_anova_f()` and
> `es_from_pt_bis_r()` to 9 digits. §7 came one step from catching this: it flags that
> the formula's cited source ("Cohen 1988") is really Aaron, Kromrey & Ferron (1998),
> but prescribed "rename" rather than re-derive. Two variance defects in
> `.smd_to_cor()` were also fixed (the `viechtbauer` branch discarded `vd`; its `h`
> used the residual df) — neither is contradicted by §4, whose five disproved claims
> touch none of this.

---

## 2. Mismatches that survived refutation

### 2a. `reverse_2x2` inverts the tetrachoric CI instead of reflecting it
`R/internal_multiple_formulas.R:545-550` and `:742-747`. Bounds are negated in
place with no swap, so `lo > up` and the point estimate falls outside its own
interval:

```
es_from_2x2(143,52,41,164)                    r =  0.7458746  CI ( 0.6591207,  0.8326284)
es_from_2x2(143,52,41,164, reverse_2x2=TRUE)  r = -0.7458746  CI (-0.6591207, -0.8326284)  <- inverted
correct                                                        CI (-0.8326284, -0.6591207)
```
The identical bug on the OR path was already found and fixed, with a comment
explaining exactly this (`R/es_from_stand_OR.R:351-353`: *"Swapping alone left
inverted, wrong-signed bounds."*). Apply the same negate-and-swap. Not silent —
`summary()` reports `[INVALID] Inverted CI` — but every reversed 2×2 row is affected.
*Skeptic could not refute; called it "airtight on the 2×2 path".*

### 2b. `or_to_cor = "bonett"` silently degrades to `lipsey_cooper`
`R/es_from_stand_OR.R:328-339`. Bonett's `c` needs `small_margin_prop`, a
user-entered column with **no auto-derivation anywhere**. When it is blank — the
normal case — the row falls out of the `nn_miss` gate and the `lipsey_cooper` value
computed earlier is left in place. `convert_df()`'s default is `or_to_cor = "bonett"`
(`R/main_convert_df.R:279`), so the default method usually does not run:

```
or=2, logor_se=0.2, n_exp=50, n_nexp=50, n_cases=40, n_controls=60, measure="r"
  returns r = 0.188  (= lipsey_cooper)     Bonett with small_margin_prop=0.4: r = 0.259
```
38% understatement, `info_used` still says only `"or_se"`, no flag.
*Skeptic refuted the framing but not the defect: "the framing is refuted, the defect is not."*
Fix: either derive `small_margin_prop` from the margins when available, or fall back
explicitly and record it in `info_used`.

### 2c. Docblock asserts a numeric check the code does not produce
`R/internal_multiple_formulas.R:1233-1234` claims single-group `morris_dav`
reproduces Bonett 2008 Example 2's `var = 0.0148` "to the printed digit". Package
gives `0.0146674131`. Bonett's own arithmetic gives `0.0147449109`, so **0.0148 is
his rounding slip, not a package error** — and recomputing eq 10 with the package's
bias-corrected `g` gives `0.0146674131`, identical to 12 digits. The formula is
right; the comment is wrong. *Skeptic could not refute the narrow point.*
Fix the comment, don't touch the code.

---

## 3. RESOLVED — the tetrachoric CI, but the diagnosis was mis-attributed

Refutation **failed on the substance, succeeded on the framing.** The defect is real
but it is not the one originally described, and the originally prescribed fix is wrong.

**Confirmed:** the r-scale interval escapes [-1,1] (Example 2 `f=(4,6,1,89)` gives
(0.650, 1.078)), at 11.8% of random tables / 27.7% under dichotomised
bivariate-normal sampling, and it is **not flagged** — B1 (`internal_flags.R:1467`)
tests the point estimate only.

**Mis-attributed:** B&P 2005 is *not* the source for this route — `man/es_from_2x2.Rd`
cites Cooper/Cochrane/Lipsey/etc. and the roxygen says it follows **metafor**. B&P 2005
is cited for `or_to_cor = "bonett"`, and *that* route already implements B&P eqs 8a/8b
faithfully. `metafor::escalc(measure="RTET")` returns a bit-identical CI, and
`tests_save/checked/test-2x2.R` asserts that parity to 1e-10 — a deliberate contract.
B&P Example 1 also reproduces exactly (`.335/.048/(.240,.429)`, the paper's printed SAS
values). And switching to eqs 8a/8b would be wrong: they belong to a different
estimator (ρ̂* = .835 vs ML ρ̂ = .864) and would decouple the reported `r_se` from the CI.

**Why it is still a defect — the argument is internal, not external.** B7b/B8b
(`internal_flags.R:1576-1618`) already flag precisely this pattern for alpha and ICC:
a raw-scale Wald bound escaping the parameter space while the point estimate is valid.
`r` is bounded identically and has no analogue, so "metafor does it too" is no defence.
The indefensible part is the silence.

A coverage simulation (12 cells × 400 reps) does back the statistical criticism —
current Wald-r mean 0.9236 / worst 0.8700, the worst of four candidates, independently
reproducing B&P's .921/.821; `tanh(z_ci)` is best at 0.9606/0.9325. But switching would
break metafor parity and make A6 fire on every tetrachoric row (halfwidth/se = 2.4714
vs 1.9600), so it is an author decision, not a drive-by fix.

**Fix:** add the `r` analogue of B7b/B8b — purely additive, breaks nothing. Full
instructions in `HANDOFF-package-fixes.md` FIX 6, which also records why Option 2
(replacing the interval) is deliberately deferred.

## 3-old. Superseded framing (kept for the record)

`R/internal_multiple_formulas.R:574-575` builds the r-scale interval as
`r ± 1.96·√vr`. That is **Bonett & Price 2005 eq (2)** — the Wald interval those
authors present *as the method to be improved on* (their reported coverage: avg
.921, worst .821). It is unbounded, so on the paper's own Example 2 `f=(4,6,1,89)`
metaConvert returns `CI = (0.650, 1.078)` — upper bound > 1 — where B&P eqs 8a/8b
give (0.488, 0.956). A 382-table sweep put out-of-range bounds at 14.4%. The
in-bounds `tanh(z_ci)` is already computed two lines away.

Not flagged: Category B tests `|r| > 1` on the point estimate only, so
`summary(flags = TRUE)` returns an empty flag string for that row.

**This is the one finding I would most want a second opinion on before acting**,
because the claim is that a shipped CI is the wrong one of two intervals in the same
paper, and its refutation never ran.

### 3b. RESOLVED — the ANCOVA unbiasedness note IS false, on 8 routes

Refutation **failed**; confirmed by two independent Monte Carlo runs (40,000 reps,
and a 4,000-rep replication with a different DGP and seed) plus a 15-digit algebraic
identity. `R/es_from_ANCOVA_statistics.R:61` back-derives
`d <- ancova_t * sqrt(1/n_exp + 1/n_nexp) * sqrt(1 - cov_outcome_r^2)`, which recovers
`MD_adj` only if the leverage term `D` is zero, so **d is attenuated by
`1/sqrt(1 + δ_x²/4)`**:

| δ_x | adjusted-means route | t route | ratio | predicted |
|---|---|---|---|---|
| 0.0 | +0.4% | +0.1% | 0.9975 | 1.0000 |
| 0.5 | +1.2% | **−2.1%** | 0.9671 | 0.9701 |
| 1.0 | +0.6% | **−10.2%** | 0.8925 | 0.8944 |

Affected: `ancova_t`, `ancova_f`, `ancova_t_pval`, `etasq_adj`, `ancova_md_se`,
`ancova_md_ci`, `ancova_md_pval`, and `ancova_means_se` (the last attenuated LESS -- 0.942 vs 0.893 at delta_x = 1; see 99_route_equivalence.R). Clean:
`ancova_means_sd`, `..._pooled_adj`, `ancova_md_sd` — so the `@note` at
`R/es_from_ANCOVA_means.R:62` is correct where it sits.

**metaConvert is not deviating from the literature standard** — `compute.es`'s
`a.tes`/`a.fes`/`a.pes` implement the same Cooper table 12.3 formula and share the
limitation, which is why the archived cross-package test passes. Only the
unbiasedness *claim* is false. d and its SE shrink together, so the p-value and
"CI excludes zero" stay exactly right; the magnitude is ~10% understated at δ_x = 1
and the study takes ~+26% excess inverse-variance weight, and it does not wash out in
pooling (fixed-effect pooled bias −10.8% vs −0.1%).

Wording fix only — `D` needs the covariate group means and within-group SS_x, which
the wide format cannot carry. `CLAUDE.md:87` also needs its "and the point estimate
stays unbiased" clause scoped to the means/MD-with-SD routes. Full instructions in
`HANDOFF-package-fixes.md` FIX 5.

**Sourcing caveat:** `data-raw/ANCOVA to SMD/Lai_and_Kelley_...pdf` has **no text
layer** — `pdftools::pdf_text()` hangs on it. An earlier pass claimed to have "read
all 21 pages of extracted text"; that claim is not credible and should be discounted.
The eq (5) structure was instead confirmed against `lm()`'s exact ANCOVA SE to a ratio
of 1.000000000000000, which is stronger evidence than the prose would have been. If
the paper must be quoted verbatim, it needs OCR.

---

## 4. Refuted claims — do not re-raise

| claim | why it is wrong |
|---|---|
| `.se_from_or` docs say "mean of the standard error", code takes mean of the variance | The claimed quantity, direction and magnitude were all arithmetically wrong. Only an unquantified one-word copy-edit survives. |
| DPJ root selection ignores the paper's cut-off (eq 21) | Mis-attribution. Those helpers are metaumbrella's **option B**, documented at `R/es_from_stand_OR.R:42-54` and asserted by existing tests — a different, deliberate method, not metaConvert's DPJ implementation. |
| `small_margin_prop` needs bounds checking | Three of its four load-bearing assertions false; both proposed remedies defective, one actively harmful. |
| `pre-post-testing.md` staleness | The narrow fact survives; everything escalation-worthy does not — and it points at the copy that reaches nobody while missing the one that ships. |

And, from earlier in this work: `viechtbauer` d/g slots are **not** swapped (the
shipped convention is closer to unbiased — simulation disproved the "fix");
`unit_type = "raw_scale"` does **not** halve the effect (it correctly selects raw
units, which is Mathur's Δ); forwarding `n1i`/`n2i` to fix unequal groups is
**contested** — Pustejovsky argues against it for experimental designs.

---

## 5. Methods with no identifiable source

The real finding for the paper: metaConvert offers these and no publication backs them.

1. **`or_to_rr = "metaumbrella_cases"` / `"metaumbrella_exp"` / `"transpose"`**, and
   `rr_to_or = "metaumbrella"` / `"transpose"` — cited only to the metaumbrella
   software announcement. The grid-search reconstructions have no published accuracy
   characterisation anywhere. This is study 04/05's contribution.
2. **`.se_from_or()`** — see §6.
3. **`grant`'s standard error.** Grant 2014 gives **no** variance, SE or CI anywhere.
   metaConvert's is a Greenland back-out from the transformed CI bounds; measured at
   0.7-5% of the exact delta-method SE across p₀ ∈ [0.05,0.7] × OR ∈ [0.25,8], so
   defensible, and the Rd discloses it. But it is not Grant's.
4. **`smd_to_cor = "lipsey_cooper"`'s z-scale variance** and **`or_to_cor =
   "lipsey_cooper"`'s** — self-documented metaConvert originals.
5. **VanderWeele's √OR is absent**, and the gap bites exactly where he aimed it:
   an adjusted OR with no reported incidence, where metaConvert's only surviving
   route is `transpose`.

---

## 6. `.se_from_or()` — the strongest single result

Imputing `SE(logOR)` from an OR plus the case/control margins is a **non-identified**
problem, and no paper in the folder attempts it: Di Pietrantonj and Veroniki both
*require* `SE(logOR)` as an input; Bonett 2007 imputes cells but needs all four
marginal proportions, which makes the table uniquely identified. `.se_from_or()` is a
byte-for-byte copy of `metaumbrella:::.estimate_se_from_or()`.

It takes the arithmetic **mean of the logOR variance** over the one-parameter family
of compatible tables. Both papers document the geometry that dooms that choice:
`Var(logOR)` is strictly convex in the free cell with a unique closed-form minimiser
(DPJ eq 20/21 + Fig 3; Veroniki eq 3 + Fig 1), so a uniform mean is dominated by the
two divergent tails.

- 1.72× inflation reproduced exactly (`1.7179266088` at br=0.5/OR=1/N=400)
- **unbounded in N**: 2.0226 at N=4000; 2.2078 at N=4000/br=0.5/OR=5
- median instead of mean gives `1.1547005384` = √(4/3)
- **decisive**: the "mean" is not a well-defined functional. Reweighting the *same*
  family uniformly over the implied control-group risk rather than over the free cell
  changes the imputed variance from `0.20663518` to `0.75993193`.

That last point is the one to publish: the estimator's value depends on an arbitrary
choice of parameterisation, which no amount of tuning fixes.

---

## 7. Documentation fixes (code is right)

| location | problem |
|---|---|
| `R/es_from_stand_OR.R:113` | Rd LaTeX for Bonett's `c` puts `1 −` inside the `/5` numerator and uses raw **counts** where paper and code use **proportions**. Code (`internal_multiple_formulas.R:602`) is eq (3) exactly. |
| `R/es_from_stand_OR.R:103-104` | Rd for the digby z CI omits the square: shows `√((c²/4)·logor_se)`, code has `logor_se²`. |
| `or_to_rr = "grant"` docs | Should cite **Zhang & Yu 1998** as the estimator, Grant 2014 as the communication paper. |
| `R/internal_multiple_formulas.R:1233` | Bonett Example 2 comment (see §2c). |
| `.pooled_pre_post_to_smd` docblock | Attributes the whole pooled d_av variance to Bonett eq 19; only the `g²` fourth-moment coefficient is. The leading term is metafor's n-based pooled form, departing from eq 19 by up to **38% at strongly unequal arm sizes** — a regime the MC calibration test only exercises at n = 12/11. |
| `data-raw/PRE-POST-SMD/pre-post-testing.md` | No OUTDATED banner; maps `morris_dav` to metafor `SMCRP` when the code implements `SMCRPH` (4.2% apart). `pre-post-strategy.md` *does* carry a banner and self-corrects. |
| `data-raw/Cohen 1988.pdf` | **Not Cohen 1988.** It is Aaron, Kromrey & Ferron (1998), "Equating r-based and d-based Effect Size Indices". Zero occurrences of "covarian". Rename. |
| Cooper eq 12.26 | The source for the ANCOVA `d` variance is **not in the repo**; could not be checked. |

Confirmed correct: **CLAUDE.md's statement of the ANCOVA limitation is right.** The
omitted term is the leverage term it names, verified by reconstructing Lai & Kelley
eq (5) and reproducing `lm()`'s exact ANCOVA SE to a ratio of 1.000000000000.
Magnitude: 0.5% SE understatement at standardised covariate imbalance 0.2, 3.1% at
0.5, 11.8% at 1.0.

---

## 8. What to do

**Package fixes**
1. `reverse_2x2` / `reverse_phi`: negate **and swap** the tetrachoric CI bounds (§2a).
2. `or_to_cor = "bonett"`: stop silently falling through to `lipsey_cooper` (§2b) —
   derive `small_margin_prop` where possible, else fall back explicitly and say so in
   `info_used`.
3. Decide the tetrachoric CI question (§3) — **get a second opinion first**.
4. Consider implementing the DPJ negative-discriminant fallback (impose `D = 0`);
   Veroniki puts it at ~15% of studies reported to 2 dp, where metaConvert now returns NaN.
5. Consider adding VanderWeele's √OR as an `or_to_rr` value.
6. Blank `n_cov_ancova` silently yields `g = NA`; blank `baseline_risk` in the DPJ
   root choice silently takes the rare-event root. Both deserve a warning.

**Documentation** — everything in §7. Cheap, and a referee will spot-check them.

**What the simulation should measure** (the literature does not settle these)
- the accuracy of the metaumbrella reconstructions (no published characterisation)
- `.se_from_or()`'s inflation, and the parameterisation-dependence in §6
- `grant`'s back-out SE vs the exact delta-method SE
- VanderWeele √OR vs the baseline-risk-requiring routes across event rates
- the pooled d_av leading-term departure at unequal arm sizes
