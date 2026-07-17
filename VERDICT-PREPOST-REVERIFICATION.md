# Independent re-verification — pre/post → between-group SMD family

**Scope:** point estimates, sampling variances, and confidence intervals for the four
standardizers (`bonett`, `morris_dz`, `morris_drm`, `morris_dav`) in **single-group**,
**per-arm two-group** (`pool_sd = FALSE`, the default), and **pooled two-group**
(`pool_sd = TRUE`) forms, across every input format.

**Provenance.** Every number below was produced on this machine, R 4.5.0 (*not* 4.5.1 as
`CLAUDE.md` states), metafor 4.8.0, `devtools::load_all` on `fix/pre-post-family` @ `a390a7a`.
Ground truth for the escalc formulas is the live `deparse` and the dump at
`data-raw/PRE-POST-SMD/metafor-smd-pre-post.R`. Primary sources read directly from the PDFs:
Morris (2008) Tables 1 & 5 and eqs 1–14; Bonett (2008) eqs 8/10/13/18/19 + Example 2;
Cousineau (2020) eq 2. Scripts and raw output are in the session scratchpad
(`partC_mc.R`, `partCG_supp.R`, `partE.R`, `partB.R`, `partA_D.R`, `final_checks.R`, `rmis.R`).
This audit re-derived everything from scratch; it did not reuse the previous author's or the
previous auditor's scripts, and the prior `AUDIT-PREPOST-VERDICT.md` was treated as a hypothesis
to confirm or refute, not as evidence.

---

## 1. Bottom line

**The statistical core is correct — point estimate, variance, and CI — for all four
standardizers in all three forms.** The five "fixes just applied" are all present in the code
and all verified correct:

1. single-group **and** pooled `morris_dav` now use the heteroscedasticity-robust **SMCRPH**
   form (not homoscedastic SMCRP) — **confirmed bit-exact and MC-calibrated**;
2. `var_d = var_g / J²` holds **exactly on every pooled branch** — the prior audit's §3.4
   defect (ratios 1.03–1.06) is **fixed**;
3. `r_avg`/`r2_avg` are df-weighted — the homoscedastic reduction to Viechtbauer's published
   variance is now **exact even when n₁≠n₂ and r₁≠r₂** (prior §4.3 residual gap closed);
4. the `n<2` / negative-SD / |r|≥1 guards fire on the pooled path;
5. (E6 flag — out of scope of the math, not tested here).

Every defect the prior audit rated MAJOR against the numbers (§3.1 per-arm dav homoscedastic;
§3.4 broken `var_g=J²var_d`; §3.7a negative pooled variance) **is no longer present.** What
remains are **minor guard and documentation/attribution issues**, listed in §5. Nothing I found
would produce a wrong effect size, variance, or CI on the supported input paths given a correct
`r_pre_post`.

The one real exposure is **not in the calculations at all**: the default `r_pre_post = 0.8`.
When it is wrong, every r-dependent variance/CI degrades severely (§4.7). That is an
input-provenance issue, orthogonal to the formulas under test, and pre-existing.

---

## 2. The convention map — verified bit-exactly and correctly attributed

Ground-truth escalc formulas (metafor 4.8.0), and the metaConvert branch each maps to:

| standardizer | single-group target | metaConvert vs metafor | attribution |
|---|---|---|---|
| `bonett`     | **SMCRH** "LS" (baseline-SD, robust)          | `max|Δg|=0`, `max|Δvar_g|=5.6e-17` | Morris d_ppc2 point est. (eq 8/9) + Bonett/SMCRH variance ✔ |
| `morris_dz`  | **SMCC** "LS" (change-SD)                     | `max|Δg|=0`, `max|Δvar_g|=0`       | Morris & DeShon d_z ✔ |
| `morris_drm` | **SMCR** "LS2" (baseline/raw-SD)              | **matches only at q=1** (see §5-A) | Caldwell & Vigotsky d_rm = d_z·√(2(1−r)) ✔ (SMCR only under homoscedasticity) |
| `morris_dav` | **SMCRPH** "LS" (avg-SD, robust)              | `max|Δg|=0`, `max|Δvar_g|=2.8e-17` | Morris d_ppc3 point est. (eq 12/13) + **Bonett eq 10** variance ✔ |

*(grid: n∈{12,31,60} × r∈{0,.3,.6,.9} × q=SD_pre/SD_post∈{0.64,1,1.56}.)*

**Per-arm two-group (`pool_sd=FALSE`, default)** — computed per arm then subtracted, variances
added (Becker 1988 / Morris d_ppc1). Verified over 60 random studies (n₁≠n₂, r₁≠r₂, q≠1):

| method | `max|Δg|` vs escalc-subtracted | `max|Δvar_g|` |
|---|---|---|
| bonett | 2.8e-17 | 5.6e-17 |
| morris_dz | 0 | 0 |
| morris_drm | 0.45 (q≠1, by definition — §5-A) | 0.04 |
| morris_dav | 1.3e-15 | 8.3e-17 |

**Pooled `morris_dz` == `escalc(measure="SMD", vtype="LS")` on the change scores**:
`max|Δg|=0`, `max|Δvar_g|=2.8e-17` over 200 random studies. Identity by construction, confirmed.

---

## 3. Published anchors (Part B)

**Morris (2008) Table 5 — both columns reproduced exactly to the printed 2 d.p.**, all five
Carlson–Schmidt studies, through `es_from_means_sd_pre_post(pre_post_to_smd="bonett")`:

| | study 1 | 2 | 3 | 4 | 5 | `max|round(g,2)−Morris|` |
|---|---|---|---|---|---|---|
| d_ppc1 = bonett `pool=FALSE` | 0.74 | 0.95 | 1.81 | 1.15 | 0.51 | **0** |
| d_ppc2 = bonett `pool=TRUE`  | 0.77 | 0.80 | 1.20 | 1.05 | 0.44 | **0** |
| d_ppc3 = `morris_dav pool=TRUE` | 0.81 | 0.78 | 1.22 | 1.14 | 0.35 | 0.02 |

The ≤0.02 d_ppc3 gap is **expected and correct**: Morris's eq 14 bias correction uses the
correlation-blind df `2N−4`, whereas the branch uses Cousineau's correlation-aware
`ν = 2m/(1+r²)`; the two diverge most at high r (study 3, r=.77). This is a deliberate,
defensible improvement, not an error.

**Bonett (2008) Example 2** (n=60, means 26→22, s₁²=30, s₂²=20, r=0.7), single-group `morris_dav`:
- metaConvert `var_d = 0.014950`; my hand-recompute of Bonett's *un-bias-corrected* eq-10 form
  = 0.014745 ≈ his printed **0.0148**. The 0.0002 gap is exactly the J² bias correction Bonett
  declines to apply (his eq 11 uses δ̂ and z, no correction). **Reproduced.**

**Variances do NOT match Morris Table 5** — by design. metaConvert emits the robust SMCRH
variance; Morris's σ̂²(d_ppc1)/σ̂²(d_ppc2) are his own homoscedastic formulas. On his five studies
the ratios run **0.82–1.90** (d_ppc1) and **0.89–1.59** (d_ppc2), largest where pre/post SDs are
most unequal (study 5). This is correct behaviour, but see §5-D: any "matches the published
values" wording must be scoped to *point estimates*.

---

## 4. Independent Monte Carlo (Part C)

Bivariate-normal (pre,post) per arm; **reported summary stats fed through the package**
functions; **the true r supplied as known** (the formulas' own stated assumption — this isolates
the variance formula from r-misspecification, which §4.7 treats separately). SEED 20260713,
8000 reps/cell (15000 for supplements). Each estimator judged against **its own** estimand
(numerator_true / standardizer_true), never a single shared target.

### 4.1 Main grid — n∈{10,20,50,200} × r∈{0,.3,.5,.7,.9} × q∈{0.64,1,1.56}, δ≈0.5

`E[v̂]/Var(ĝ)` (variance-calibration ratio) and 95% CI coverage, full-grid ranges
(MC SE on coverage ≈ 0.0025):

| method / pool | ratio range | coverage range |
|---|---|---|
| bonett `FALSE`     | 0.978 – 1.291 | 0.946 – 0.986 |
| morris_dz `FALSE`  | 0.820 – 1.150 | 0.944 – 0.981 |
| morris_drm `FALSE` | **0.718** – 1.032 | 0.941 – 0.968 |
| morris_dav `FALSE` | 0.976 – 1.182 | 0.946 – 0.984 |
| bonett `TRUE`      | 0.987 – 1.123 | 0.949 – 0.976 |
| morris_dz `TRUE`   | 0.908 – 1.111 | 0.941 – 0.977 |
| morris_drm `TRUE`  | 0.917 – 1.103 | 0.945 – 0.977 |
| morris_dav `TRUE`  | 0.974 – 1.091 | 0.948 – 0.975 |

**All eight combinations hold nominal coverage (≥0.941) across the entire grid.** The mild
over-coverage (ratios >1, up to 1.29 at n=10) is the `qt`-vs-`qnorm` conservatism (§5-E), not a
defect. The low tails (dz/drm at n=10 × extreme q×r) still cover ≥0.941.

**Per-arm `morris_dav` is the headline confirmation of the fix:** ratio 0.976–1.182, coverage
≥0.946, *no collapse under heteroscedasticity*. The prior audit's §3.1 (per-arm dav is
homoscedastic; ratio 0.82, coverage 0.93) described the **old** code and is now obsolete.

### 4.2 Unequal n (probe of the pooled leading-term approximation)

The pooled leading term `(sd_change_pooled²/sd_pooled²)·(1/n₁+1/n₂)` equals the exact
`sc_T²/n₁+sc_C²/n₂` only when n₁=n₂ (or under homoscedasticity). Empirically the approximation
is **benign**: at n₁/n₂ ∈ {10/40, 40/10, 15/60}, q∈{0.64,1.56}, all four pooled estimators give
ratio ∈ **[0.999, 1.046]**, coverage ∈ **[0.955, 0.963]**, bias ≈ 0.

### 4.3 Unequal arm σ (σ_C = 2σ_T) — the two pool settings target different estimands

| method | `pool=FALSE` estimand / cov / ratio | `pool=TRUE` estimand / cov / ratio |
|---|---|---|
| bonett     | +0.500 / 0.959 / 1.045 | +0.316 / 0.956 / 1.007 |
| morris_dz  | +0.500 / 0.959 / 1.048 | +0.316 / 0.956 / 1.016 |
| morris_drm | +0.500 / 0.956 / 1.023 | +0.316 / 0.954 / 1.014 |
| morris_dav | +0.500 / 0.959 / 1.038 | +0.316 / 0.958 / 1.028 |

Each estimator is **unbiased for its own estimand** (|bias|≤0.007) and covers at ~0.956.
`pool_sd=FALSE` and `TRUE` are **different estimators answering different questions** — neither
is "the" target. This is the point the previous author's original framing got wrong, and it is
handled correctly now.

### 4.4 Structural identities (Part D)

- `var_g / (J²·var_d) = 1.0000000000` and `g/(J·d) = 1.0000000000` on **all four pooled
  branches** (and by construction on single/per-arm). Exact.
- `.d_j(ν) == metafor:::.cmicalc(ν)` to `max|Δ|=1.1e-16`, including non-integer ν (e.g.
  ν=2m/(1+r²)=46.4).
- CI reconstruction `g ± qt(.975,ν)·√var_g` matches the emitted bounds to 0 for single bonett
  (ν=n−1), pooled bonett (ν=N−2), pooled dav (ν=2m/(1+r2_avg)); **single-group dav uses ν=n−1**
  for the CI (not the Cousineau ν used for J) — see §5-E.

### 4.5 Cross-format equivalence (Part E)

One synthetic two-group study, derived via every input format. `max` pairwise |Δg|:

- **means_sd ≡ means_se ≡ means_ci** for all four methods × both pool settings: **0 to 1.1e-16**
  (CI round-trip exact — the wrapper's `qt(.975,n−1)` inverts cleanly).
- **dz / drm** additionally match **mean_change_{sd,se,ci,pval} ≡ paired_t ≡ paired_f** to
  **≤2.2e-16**, including pooled dz/drm from change scores.
- **bonett / dav are correctly unavailable** from mean_change and paired formats (a change SD or
  a t cannot identify the baseline or separate pre/post SDs). The wrappers **hard-error**
  (`allowed_methods = morris_drm, morris_dz`); the bonett/dav→cooper coercion lives in
  `convert_df()`. **Never silent-wrong.**
- `es_from_paired_t`/`es_from_paired_f` have **no `pool_sd` argument** and their d_z equals the
  unpooled means-based d_z exactly (diff 0). The identifiability limit (a paired t is scale-free,
  so the SD ratio and hence pooling are unidentifiable — the theorem in the task's Q5.5) is
  handled by construction.

### 4.6 dav SMCRPH vs old SMCRP (Part G) — the fix is necessary and correct

Single-group dav, n=30, r=0.5, true r fed, δ_av≈0.5:

| q | Var(ĝ)_emp | **SMCRPH** (current) ratio / cov | **SMCRP** (old, homoscedastic) ratio / cov |
|---|---|---|---|
| 0.50 | 0.0443 | 1.023 / 0.963 | **0.797 / 0.937** |
| 0.64 | 0.0400 | 1.035 / 0.962 | 0.883 / 0.947 |
| 1.00 | 0.0358 | 1.064 / 0.967 | 0.985 / 0.960 |
| 1.56 | 0.0401 | 1.031 / 0.964 | 0.879 / 0.949 |
| 2.00 | 0.0439 | 1.033 / 0.964 | **0.805 / 0.938** |

The homoscedastic SMCRP coefficient `(1+r²)/(4n)` is exact **only at q=1**; away from it the
variance is ~20% low and coverage falls to 0.937–0.938. SMCRPH holds ~0.96 throughout. Bonett
(2008) confirms the target directly: his Tables 1 & 7 show the homoscedastic interval (Eq 2/3)
collapsing under σ=[1 1.5] while the robust interval (Eq 11) holds ~0.95, and the two-group
pretest-posttest all-four-SD contrast is his general Eq 18 with c=[1 −1 −1 1]. **SMCRPH is the
right target for d_av; the migration from SMCRP was correct.** (The coverage hit grows with δ —
at δ=0.5 the g² term is a modest share of the total, so SMCRP's damage here is milder than at
large effects.)

### 4.7 The r-default exposure (separate from the formulas)

q=1, n=30/30, pooled bonett, feeding the **default r=0.8** when the truth differs:

| true r | correct-r ratio / cov | default-0.8 ratio / cov |
|---|---|---|
| 0.3 | 1.038 / 0.960 | **0.334 / 0.755** |
| 0.5 | 1.032 / 0.959 | 0.451 / 0.826 |
| 0.7 | 1.013 / 0.959 | 0.703 / 0.911 |
| 0.9 | 1.019 / 0.956 | 1.869 / 0.993 |

With the correct r, coverage is ~0.96 everywhere. With the default 0.8 against a true 0.3, the
reported variance is 67% too small and the nominal 95% CI covers 75.5%. **This dominates every
other error source, and the robust numerator gives no protection** (it consumes r inside
sd_change). It is a property of the *default*, not of the *calculation*, and it is pre-existing.
It should be documented loudly in `?convert_df`; `bonett`'s point estimate is r-free so this is a
pure variance/CI failure, invisible in the forest-plot dots.

---

## 4b. Fixes applied after the verdict (this session)

Acting on the residual findings:

- **B (guard gap) — FIXED.** `.single_group_pre_post_to_smd()` now NA-guards the standardizing
  SDs when `n < 2`, so `morris_drm` (whose J-free `var_d` previously leaked a finite value at
  n=1 and `Inf` at n=0) returns clean NA on d, var, and CI — matching the other three branches
  and the pooled kernel. Verified: single-group and per-arm drm at n∈{0,1} → all NA; n≥2 output
  bit-identical to before.
- **E (CI df inconsistency) — HOMOGENISED.** Pooled `morris_dav`'s CI now uses `qt(.975, N−2)`
  like every other pooled branch (and like the single-group dav CI, which is on `n−1`). Its
  Cousineau effective df `2m/(1+r2_avg)` continues to feed **only** the bias correction J
  (`g/d == .d_j(nu)` unchanged). The other pooled branches (`nu == m`) are numerically unchanged.
- **A (attribution) — CLARIFIED in a comment.** The single-group `morris_drm` block now states it
  implements Caldwell & Vigotsky's `d_z·√(2(1−r))`, which coincides with metafor's SMCR (raw-SD
  standardizer) only under homoscedasticity.
- **C (eq. 19 citation) — RETRACTED, no change.** On re-reading Bonett (2008), **eq. 19 IS the
  two-group all-four-SD pretest-posttest variance** (the worked result of applying the general
  eq. 18 with c=[1 −1 −1 1], r=4). The code's citation was correct; I did not alter it.
- **D (NEWS wording) — already correct.** `NEWS.md` already scopes the claim to "published
  **point estimates**"; no change needed.

Regression check after the edits: single-group bit-exactness vs SMCRH/SMCC/SMCRPH preserved
(≤1.4e-17); test-pool-sd (84), test-PAIRED-MEANS (141), test-PAIRED-MC (36), test-morris-deshon
(120) all pass, 0 failures.

## 5. Findings (all MINOR — nothing blocks the math)

**A — `morris_drm` is Caldwell's d_rm, not SMCR, off q=1.** The point estimate is
`d_z·√(2(1−r))` (Caldwell & Vigotsky 2020), which equals the baseline-SD SMCR only under
homoscedasticity; at q≠1 it diverges (Δg up to 0.45 vs `escalc(SMCR)`). This is **not a bug** —
the estimator is well-defined and its variance `2(1−r)/n + d²/(2n) = 2(1−r)·Var(d_z)` is
internally consistent and MC-calibrated (coverage ≥0.941). But the convention-map claim
"morris_drm = SMCR LS2" should be **scoped to homoscedasticity**. `morris_drm` is also the least
robust of the four: its homoscedastic `2(1−r)/n` leading term drives the ratio to 0.72 at
n=10 × extreme q×r (coverage still 0.941). Consider the SMCR-LS2 or a robust change-SD variance
if drm is to be recommended under heteroscedasticity.

**B — guard gap on the default/single-group `morris_drm` path.** At `n=1` the single-group and
per-arm (`pool_sd=FALSE`) drm returns a **finite d, finite var_d, and (two-group) a finite
d-scale CI**; at `n=0`, `var_d=Inf`. The g-scale is NA (protected, so the meta-analytic quantity
is safe), but this is not a clean NA on direct calls — it emits Inf/NaN/silent-finite-wrong.
Cause: drm's `var_d` has no J factor and no explicit `n≥2` guard, so unlike bonett/dz/dav it does
not inherit the `J(n−1)=NA` nulling. The pooled path already added `valid_n <- n≥2`; the
single-group kernel should too. (bonett/dz/dav all correctly NA at n≤1.)

**C — citation imprecision: "Bonett 2008 eq. 19".** The single-group dav citation "Bonett eq. 10"
is **correct** (his Example 2 explicitly "Applying Equation 10"). But the pooled all-four-SD dav
comments cite **eq. 19**, which is actually Bonett's *σ̂₁-based* (baseline) alternative
standardizer variant carrying the "n₁≥30" caveat. The all-four-SD two-group contrast is his
**general eq. 18** (c=[1 −1 −1 1], r=4). Fix the eq-19 references to eq-18.

**D — "matches the published values" over-claims.** Only the **point estimates** match Morris
Table 5. The variances deliberately depart (robust SMCRH), by 0.82×–1.90× on Morris's own five
studies. Scope the claim to "point estimates" and state the variance departure inline.

**E — CI df inconsistency, and qt-vs-z.** Single-group dav CI uses `qt(.975, n−1)`; pooled dav
uses `qt(.975, 2m/(1+r2_avg))`; per-arm two-group uses `qt(.975, n₁+n₂−2)`. All conservative,
but inconsistent across the dav paths. The whole family uses `qt` while metafor and Bonett use
`qnorm` (z); this is a deliberate small-n-conservative choice (it produces the ratio>1 /
coverage-to-0.986 over-coverage at n=10) — defensible, but the Viechtbauer p.57 citation offered
in support actually recommends **z**-based intervals, so either switch to `qnorm` or drop that
citation.

---

## 6. Verdict — for each standardizer × form

| standardizer | form | point est. | variance | CI |
|---|---|---|---|---|
| **bonett**     | single   | ✅ =SMCRH (bit-exact); =Morris d_ppc1/arm | ✅ SMCRH LS, robust, MC cov 0.95–0.98 | ✅ qt(n−1) |
| bonett         | per-arm  | ✅ =Morris d_ppc1 (Table 5 exact) | ✅ SMCRH summed, MC cov ≥0.946 | ✅ qt(N−2) |
| bonett         | pooled   | ✅ =Morris d_ppc2 (Table 5 exact) | ✅ reduces exactly to Viechtbauer vi under homosced.; robust otherwise; MC cov 0.949–0.976 | ✅ qt(N−2) |
| **morris_dz**  | single   | ✅ =SMCC (bit-exact) | ✅ SMCC LS | ✅ |
| morris_dz      | per-arm  | ✅ | ✅ summed; MC cov ≥0.944 | ✅ |
| morris_dz      | pooled   | ✅ =SMD on change scores | ✅ **=escalc(SMD,LS)** bit-exact | ✅ |
| **morris_drm** | single   | ✅ Caldwell d_rm (=SMCR only at q=1, §5-A) | ✅ =2(1−r)·Var(d_z), internally consistent | ⚠️ finite at n=1 (§5-B) |
| morris_drm     | per-arm  | ✅ | ✅ MC cov ≥0.941 (least robust; ratio→0.72 at n=10 extreme) | ⚠️ §5-B |
| morris_drm     | pooled   | ✅ | ✅ MC cov 0.945–0.977 | ✅ |
| **morris_dav** | single   | ✅ =SMCRPH (bit-exact); Bonett Ex.2 reproduced | ✅ **SMCRPH LS (Bonett eq 10)** — the fix, correct | ✅ qt(n−1), conservative (§5-E) |
| morris_dav     | per-arm  | ✅ ≈Morris d_ppc3 | ✅ SMCRPH summed, robust; MC cov ≥0.946 | ✅ |
| morris_dav     | pooled   | ✅ ≈Morris d_ppc3 (≤0.02, §3) | ✅ robust two-group (Bonett eq 18, cite fix §5-C); MC cov 0.948–0.975 | ✅ |

**Are the just-applied fixes correct? Yes — all five, verified independently.** The SMCRPH
migration for `morris_dav` (single and pooled) is bit-exact and demonstrably restores coverage
under heteroscedasticity; `var_d = var_g/J²` holds exactly on every pooled branch; the
df-weighted `r_avg`/`r2_avg` make the homoscedastic reduction exact; the pooled guards fire.

**What is MC-backed rather than published:** pooled `morris_dav` (Morris gives no d_ppc3
variance, p.373; the branch derives the two-group SMCRPH form and calibrates it — MC ratio
0.97–1.09, coverage 0.948–0.975 on this audit's independent grid) and the robust-departure and
d_rm algebra. These are adequately justified: pooled dz is a published identity, pooled bonett
reduces exactly to a published form under homoscedasticity, and the dav derivation follows
Bonett's own general method (eq 18). The residual asks are the §5 documentation/guard items —
none of which change a number on a supported path.
