# Independent adversarial audit — `fix/pre-post-family`

**Scope:** the pre/post (paired) → between-group SMD family.
**Provenance:** every number below was re-derived on this machine (R 4.5.0 — *not* 4.5.1 as `CLAUDE.md` says; metafor 4.8.0; `devtools::load_all` on `fix/pre-post-family` @ `a390a7a`). Monte Carlo grids are the auditor's own, written from scratch; the author's MC scripts were not used and are not in the repo. Primary sources were read directly (Morris 2008 PDF; Cousineau 2020 PDF; Viechtbauer 2007 PDF; Harrer 2025 PDF; Caldwell & Vigotsky 2020 PDF; the metafor-project Morris-2008 page; `deparse(metafor::escalc)`). Where a source could not be obtained, that is stated.

---

## 1. Bottom line

**MERGE WITH CHANGES.**

The statistical core of the rewrite is **right**, and it is right for the right reason. The four pooled variances are calibrated across a grid the author never ran (E[v̂]/Var(ĝ) ∈ [0.99, 1.12], coverage ∈ [0.95, 0.98] over q ∈ {0.64, 1, 1.56} × r ∈ {0, .3, .5, .7, .9} × n ∈ {10, 20, 50, 200}), while the published comparator it replaces collapses to a variance ratio of **0.39** with **79%** coverage at q = 0.64, r = 0.9. The point estimates reproduce **both** columns of Morris (2008) Table 5 exactly. The revert of `pool_sd = FALSE` to default was correct. The test suite is green (0 failures, 0 errors).

It cannot ship as written. There is **one real numerical defect on the default path** (§3.1), and the **documentation asserts as fact four things that are false** (§3.3) — including an arithmetic identity that is off by a factor of N/(N−2). The author's pattern is not that his code is wrong; it is that **his stated justifications keep outrunning his evidence.** That is what must be fixed, because it is what produced the two previous retractions.

And there is a **larger problem than the one the branch set out to solve** (§3.2), which the branch neither caused nor fixed.

---

## 2. What is correct

1. **The point estimates are exactly right, in both modes.** Morris's five Carlson–Schmidt studies, through `es_from_means_sd_pre_post(pre_post_to_smd = "bonett")`:

   | | `pool_sd=FALSE` | Morris *d*<sub>ppc1</sub> | `pool_sd=TRUE` | Morris *d*<sub>ppc2</sub> |
   |---|---|---|---|---|
   | | 0.7421, 0.9524, 1.8054, 1.1514, 0.5050 | 0.74, 0.95, 1.81, 1.15, 0.51 | 0.7684, 0.8010, 1.2045, 1.0476, 0.4389 | 0.77, 0.80, 1.20, 1.05, 0.44 |

   `max |round(g,2) − published| = 0` for both columns.

2. **The `pool_sd = FALSE` default is the metafor-conformant one, and the retracted "the default is a bug / 61% biased" position was flatly wrong.** Viechtbauer's own page computes `escalc("SMCR")` per arm, subtracts, and *adds* the variances — that **is** `pool_sd = FALSE`. He presents it first, meta-analyses only it, and says the pooled route *"assumes that the true pretest SDs are equal … The approach given above does not make that assumption and therefore is more broadly applicable."* Pooling buys **efficiency at the cost of an assumption**; it does not remove a bias. The revert stands.

3. **The heteroscedasticity-robust numerator is necessary, not a preference.** See §5.3.

4. **Pooled `morris_dz` is bit-exact with `escalc(measure="SMD", vtype="LS")`** on change scores (max |Δvi| = 2.8e−17 over 200 randomized studies with n₁≠n₂, r₁≠r₂, q≠1). This is an identity by construction, not luck.

5. **Three of the four single-group branches are bit-exact with metafor's LS** (`bonett`↔SMCRH, `morris_dz`↔SMCC, `morris_dav`↔SMCRP; Δvi = 0 to machine precision).

6. **`.d_j()` is the exact Gamma J**, agreeing with `metafor:::.cmicalc` to ≤ 5.6e−17, including at non-integer df.

7. **The E6/E7 flags fire on the intended data and stay silent on clean data** — verified: a two-group `means_sd` dataset with no pre/post columns emits **zero** messages and **zero** flags.

8. **Per-row `pre_post_to_smd` genuinely works per row**; the `bonett`/`morris_dav` → `cooper` coercion for mean-change and paired-t data is correctly implemented and messaged.

---

## 3. What is wrong

### 3.1 MAJOR — per-arm `morris_dav` is homoscedastic. The rewrite fixed only the opt-in path.

`R/internal_multiple_formulas.R:1090`, in the kernel of the **`pool_sd = FALSE` default path**:

```r
var_g <- 2 * (1 - r_pre_post) / n + g^2 * (1 + r_pre_post^2) / (4 * n)
```

That is metafor's **SMCRP "LS"** verbatim — whose leading term `2(1−r)/n` is the **homoscedastic** one. `SMCRPH` is the robust variant. The rewrite made all four *pooled* branches robust and left this one alone.

**My MC** (20k reps/cell, seed 20260712, true *r* supplied — so this is pure heteroscedasticity, not r-misspecification):

| q | r | n/arm | E[v̂]/Var(ĝ) | coverage |
|---|---|---|---|---|
| 0.64 | 0.7 | 10 | 0.817 | 0.949 |
| 0.64 | 0.7 | 20 | 0.814 | **0.936** |
| 0.64 | 0.7 | 50 | 0.816 | **0.929** |
| 1.56 | 0.7 | 50 | 0.835 | **0.932** |

**Coverage degrades as n grows** — the signature of a misspecified variance, not a small-sample artifact. The 0.82 is predicted exactly: at q = 0.64, r = 0.7 the assumed `2(1−r) = 0.600` versus the true `σ²_change/σ²_p = 0.5136/0.7048 = 0.729`; ratio **0.823**. An independent MC run reproduced this at 0.819 / 0.926.

Two comments above it are also false (`:1077–1078`): *"Recommended by Morris (2008) as general-purpose measure"* — Morris p.384 says verbatim **"the use of d_ppc3 is not recommended"**; and *"Robust to variance heterogeneity"* — SMCRP is precisely the **non**-robust variant.

**Fix:** use metafor's SMCRPH leading term. (MC-verified to restore ratio 0.99–1.02, coverage 0.954–0.958.)

### 3.2 MAJOR — the imputed `r_pre_post = 0.8` dominates every other error source, and the docblock's coverage claim is false at the package's own default.

The docblock (`:838–839`) asserts *"the robust form holds coverage at 0.95 in every regime tested."* **It does not.** `convert_df` defaults to `r_pre_post = 0.8`, and studies almost never report *r*. My MC at **q = 1** — zero heteroscedasticity, so the entire rewrite is irrelevant — with the default 0.8 imputed:

| true *r* | `pool_sd=TRUE` ratio / cov | `pool_sd=FALSE` ratio / cov |
|---|---|---|
| **0.3** | **0.324 / 0.748** | **0.352 / 0.772** |
| 0.5 | 0.447 / 0.819 | 0.473 / 0.838 |
| 0.7 | 0.731 / 0.919 | 0.752 / 0.918 |
| 0.8 | 1.004 / 0.958 | 1.044 / 0.961 |
| 0.9 | 1.922 / 0.995 | 1.932 / 0.994 |

At a true *r* of 0.3 the reported variance is **68% too small** and the nominal 95% CI covers **75%** of the time. Analytically: at q = 1, `sd_change² = 2s²(1−r)`, so imputing 0.8 for a true 0.3 shrinks the variance by 0.2/0.7 = 0.286. This is **larger than the ~43% heteroscedasticity effect the whole rewrite is built around**, and the robust numerator gives **no protection** — it *consumes* the imputed *r* inside `sd_change`.

`bonett`'s point estimate is *r*-free, so this is a pure variance/CI failure — arguably worse, since it is invisible in the forest plot's dots and lives entirely in its whiskers.

**This is pre-existing and not caused by the branch.** It is not a reason to block. But shipping a docblock advertising nominal coverage "in every regime tested" while the package's own default imputation destroys coverage at q = 1 is not defensible. **Scope the claim, and document the exposure in `?convert_df`.**

### 3.3 MAJOR (documentation, but load-bearing — it *is* the branch's stated justification)

The `.pooled_pre_post_to_smd` docblock (`:818–854`, echoed in `NEWS.md`) states four false things.

**(a) The "single rule" no branch implements.** It claims all four branches follow `Var(g) = Var(num)/SD_std² + g²/(2·ν_std)` with `J = J(ν_std)`, *"the single rule that metafor's own escalc() encodes for every pre/post measure."* **Every branch divides the g² term by 2N, not 2·ν_std.** Verified (n₁=n₂=20, r=0.5): pooled bonett returns `var_g = 0.11628525 = T1 + g²/80`, not `0.11652318 = T1 + g²/76`. **The code is right and the comment is wrong** — metafor's non-H LS forms put *n* in the g² denominator while evaluating *J* at *n−1*; the two df **deliberately differ**. And metafor has **no single rule**: the H variants use (n−1) in *both* terms, and SMCR/SMCRP's leading terms are homoscedastic — which is *why* SMCRH/SMCRPH exist.

**(b) A false arithmetic identity** (`:940`): *"the g^2 term is divided by 2*nu — which is exactly the (1 + r^2)/(4N) coefficient below."* With ν = 2m/(1+r²), `1/(2ν) = (1+r²)/(4m)`, **not** `(1+r²)/(4N)`. At n₁=30, n₂=40, r=0.5: 0.00459559 vs the coded 0.00446429 — **2.94% apart** (11.1% at N=20). `NEWS.md` repeats it verbatim.

**(c) "No published sampling variance exists for the pooled two-group d_z"** — self-contradicted eight lines below, where the `morris_dz` comment correctly identifies it as `escalc("SMD", vtype="LS")`. True only for pooled *d_rm* and *d_av*.

**(d) The Viechtbauer p.57 citation argues against the code beneath it.** It is cited to justify LS forms and `qt()` intervals. That page recommends **dL1z, gL1z, dL2z** — all **z**-based — and p.51 states normal bounds are *more* accurate than t-based ones. metafor uses z. The code uses `qt`. Either switch to `qnorm` or drop the citation.

### 3.4 MODERATE — `Var(g) = J²·Var(d)` is violated on the pooled path

`g = J·d` with `J` a deterministic constant, so `Var(g) = J²·Var(d)` is exact algebra. Measured (n=15/arm, r=0.5, SDs 8/10), `var_g / (J²·var_d)`:

| | bonett | morris_dz | morris_drm | morris_dav |
|---|---|---|---|---|
| `pool_sd=FALSE` | 1.000000 ✓ | 1.000000 ✓ | 1.000000 ✓ | 1.000000 ✓ |
| `pool_sd=TRUE` | **1.055104** ✗ | **1.055104** ✗ | **1.055104** ✗ | **1.034066** ✗ |

The docstring endorses this (`:827–828`): *"Var(d) obtained by substituting d for g (NOT by dividing Var(g) by J^2 — that is only an identity under the LS2 convention)."* The **identity** holds under every convention; what is convention-dependent is whether the *plug-in estimators* respect it. metafor never faces this because `escalc` emits one `vi` for one `yi`. metaConvert emits both, so they must cohere. Magnitude ≈ 3/(2ν): 5.6% at n=15/arm. Confined to the *d* scale (`g`, `g_se`, `g_ci` are unaffected, and those are what get meta-analysed).

**Fix:** keep `var_g` (it is the MC-calibrated one); derive `var_d = var_g / J²`. Delete the false parenthetical.

### 3.5 MODERATE — pooled `morris_dav` is an SMCRP/SMCRPH chimera

```r
T1    <- (sd_change_pooled^2 / sd_pooled^2) * T_change   # SMCRPH's ROBUST first term
var_g <- T1 + g^2 * (1 + r_avg^2) / (4 * N)              # SMCRP's HOMOSCEDASTIC second term
```

SMCRPH's second term is `g²(sd1⁴+sd2⁴+2r²sd1²sd2²)/(8·sdp⁴·(n−1))`, which collapses to `(1+r²)/4` **only when SD_pre = SD_post**. The coefficient is **9.5% too small** at q = 0.64 and q = 1.56 (r=0.5). **Do not overclaim it:** the g² term is a small share of `vi`, so the effect on the total is 0.03% (g=0.2) to 1.4% (g=1.5). It is an internal-consistency defect, not a numerical emergency — but the branch is robust in its numerator and quietly reverts to homoscedasticity in the standardizer-uncertainty term, *the very term Bonett (2008) exists to fix*.

### 3.6 MODERATE — `NEWS.md` overclaims "pinned to the published values"

Only the **point estimates** are pinned. On Morris's own five studies the variances depart by up to **2.39×** (`pool_sd=FALSE`) and **1.75×** (`pool_sd=TRUE`) from the metafor page. Study 3 matches to 1e−8 because it is the **only** study with SD_pre = SD_post in both arms. The departure is deliberate and correct — but a user benchmarking against the metafor page will see a 2.4× discrepancy with no warning. Reword to *"pinned to the published **point estimates**"* and state the variance departure inline.

### 3.7 MINOR

| # | Defect | Location |
|---|---|---|
| a | The `m > 0` guard does not protect `T_change = N/(n_exp·n_nexp)`. With `n_exp = 0`: finite g, `se = Inf`. Worse, `(n_exp − 1) = −1` gives a **negative pooling weight** → negative pooled variance → NaN. **Guard `n_exp >= 2 & n_nexp >= 2`.** | `:872, 886–897, 943` |
| b | `.guard_standardizer` guards only the *assembled* denominator; a negative per-arm SD survives on the **direct-call** path (which the paper's own scripts use). | `:749` |
| c | The r-imputation message says *"Under the default 'cooper'/'morris_drm' standardizer the assumed correlation scales the POINT estimate"* — but `convert_df`'s default is **`bonett`**, whose point estimate is **r-free**. Doubly wrong; same class of defect `96b6da3` fixed for the coercion message. | `main_convert_df.R:344–353` |
| d | **Direct self-contradiction introduced by this branch:** the new guidance message says *"morris_dav … Morris recommends **AGAINST** it (2008, p.384)"* while `es_from_PAIRED_MEANS.R:60` still says *"morris_dav, **recommended**"*. The message is right; the roxygen is stale. | `es_from_PAIRED_MEANS.R:60` |
| e | CI df conventions differ across the three `dav` paths (`N−2` / `n−1` / `2m/(1+r²)`). All conservative. Switching to `qnorm` dissolves the class. | |
| f | `test-formula-verification.R:15–17` attributes to **"Bonett (2008)"** a citation whose title, journal, volume, pages and DOI all belong to **Morris**. | |
| g | Stale test comments assert the **opposite of the current default** (`test-PAIRED-MEANS.R:10–21`: *"convert_df() now defaults to pool_sd = TRUE"*). | |
| h | `test-pooled-variance-calibration.R:81` grades coverage with **`qnorm`** while the package emits **`qt`** intervals — it does not validate the interval the package ships. | |

---

## 4. What is unverifiable

**4.1 Bonett (2008), *Psych Methods* 13(2):99–109 — the PDF could not be obtained.** Paywalled; no mirror. **No page or equation number from that article appears anywhere in this audit, and none was invented.**

The attribution chain nevertheless closes, via two artifacts *Bonett himself wrote*: his textbook (*Statistical Methods for Psychologists* Vol. 1) and his R package **`statpsych`**, whose help files cite `Bonett2008` **by DOI** for exactly these functions. His code reproduces metafor's SMCRH/SMCRPH standard errors to 7 s.f. on his own worked example. So:

- **metafor's SMCRH/SMCRPH *are* Bonett (2008).** The branch's Bonett citation for the robust numerator is **correct**.
- **Bonett *does* pool degrees of freedom across arms** (`df <- sum(n) − a`). So **ν = N−2 for bonett/dz/drm is directly supported** — it is a citation, not a derivation.
- **But Bonett NEVER uses ν = 2m/(1+r²). Anywhere.** Not in the textbook, not in one line of `statpsych`. He absorbs the correlation into the SE's **fourth-moment term**, never into a df. **ν = 2m/(1+r²) cannot be back-attributed to Bonett** — it remains the author's own composition of Cousineau's per-arm df under Bonett's additivity principle.
- **And Bonett's g² term for a pooled standardizer carries the fourth-moment weights** that `morris_dav` drops (§3.5). *`morris_dav` omits the exact correction Bonett wrote the paper to supply.*

**4.2 The pooled two-group variances.** Precisely:
- **Pooled *d_z*: a published variance DOES exist** — `escalc("SMD","LS")` on change scores. The code is bit-exact with it. (Harrer 2025 eq. 13 also publishes one, but it carries a spurious leading `2(1−r)` that belongs to raw-score standardizers, not a change-SD standardizer. **The rewrite made the right change here — for a reason the author never stated.** Deleting the Harrer citation was correct.)
- **Pooled *d_rm*: follows from *d_z* by exact algebra** (Caldwell & Vigotsky eq. 13 is a *definition*). A derivation, not an invention.
- **Pooled *d_av*: genuinely unsourced.** Morris gives **no** variance for d_ppc3 and says so (p.373). Cousineau gives none. Bonett gives a *different* standardizer and a *different* leading term. **This one rests entirely on the author's own derivation plus his own Monte Carlo.**

**4.3 Is MC calibration adequate evidence for shipping? — Yes, conditionally.** It is adequate when (i) the formula **reduces analytically to a published form where one exists** — it does: at SD_pre = SD_post *and* r_exp = r_nexp it collapses **exactly** to Viechtbauer's `2(1−r)(1/n_T + 1/n_C) + g²/(2N)`, so it is a **generalization, not a replacement**; (ii) the calibration is confirmed on an **independent grid** — it is (this audit's); (iii) the acceptance criteria are **data-independent** — they are.

**But three conditions currently fail:** the derivation is **asserted, not written down** (the docblock states a rule, which is false, instead of deriving the result — a reviewer cannot check an assertion); the MC test **grades a CI the package does not emit**; and the q = 0.64 test point is **justified in the code by reference to the paper's own median q**. Choosing your test point from your result is a bad habit even when harmless — and it is harmless *only because* a q-sweep shows the choice is not load-bearing. **Put the sweep in the test.**

**Caveat on "reduces exactly":** it does **not**, in general. The reduction requires `r₁ = r₂` **or** `n₁ = n₂`. With both unequal (n = 10/90, r = .3/.8) the ratio is **0.9837** — a 1.6% gap. The culprit is `r_avg`, an **n-weighted** mean of raw r. **Switching `r_avg` to a df-weighted mean `((n_i−1)r_i)/(N−2)` makes the reduction exact in all homoscedastic cases** — the very consistency property the author claimed but did not achieve. That is a principled reason to change it, and it also settles the Fisher-z discrepancy (§7).

---

## 5. The six math questions

**5.1 Is "no leading J²" what metafor does?** **Yes for the LS forms, and only for them.** No LS form carries a leading `cm²` — the author's narrow point is right. But *"the single rule metafor encodes for every pre/post measure"* is false three ways: the g² divisor is *n*, not ν_std; the H variants use (n−1) in **both** terms; and SMCR/SMCRP's leading terms are **homoscedastic**, which is the whole reason SMCRH/SMCRPH exist. **And it is emphatically not true of Morris:** eq. 25 **does** carry a leading `c_P²` (p.373: *"The variance of d_ppc2 is c_P² times the variance of g_ppc2"*). "No leading J²" is a **metafor-LS fact, not a universal one.**

**5.2 Is `(1+r²)/(4N)` really `1/(2ν)`? Is Cousineau eq. 2 the source? Do single-arm dfs add?**
- **The identity is FALSE** (§3.3b). The code is right; the docblock is wrong.
- **The FORM is genuinely Cousineau's; the APPLICATION is not.** Eq. 2 verbatim: *"ν = 2(n − 1)/(1 + ρ²) are the degrees of freedom."* Attribution of the *single-arm* df is accurate. **But** Cousineau is a strictly single-group paper; he gives **no sampling variance at all**; and **his derivation assumes homoscedasticity, explicitly** (*"The demonstration assumes that the variances are homogeneous"*). His ν is the σ_pre = σ_post **special case** of the general Satterthwaite df. Using it inside a branch **sold as heteroscedasticity-robust** is incoherent — and **metafor does not make this mistake** (SMCRP uses Cousineau's ν; SMCRPH uses the general one).
- **The dfs DO legitimately add** — and not merely because `2(n₁−1)/(1+r²) + 2(n₂−1)/(1+r²) = 2m/(1+r²)`. A proper Satterthwaite calculation on the code's actual (n_i−1)-weighted pooled standardizer returns the same number, **because the pooling weights happen to be ∝ ν_i**. And **Bonett independently sanctions df-addition across arms.** So ν = N−2 for bonett/dz/drm is citable.
- **Verdict: HALF TRUE. Keep the derivation; label it a derivation, not a citation.**

**5.3 Is replacing Viechtbauer's published `vi` with the robust form justified?** **USE ROBUST. Do not offer both.** (Note the premise is shaky: that formula is *not* in Viechtbauer 2007 — that paper contains no two-group pre/post design. The comparator is the *metafor-project page*, which he himself calls non-canonical.) At n = 20/arm, δ = 0.5, **E[v̂]/Var(ĝ) / coverage**:

| | r=0 | r=.3 | r=.5 | r=.7 | r=.9 |
|---|---|---|---|---|---|
| q=0.64 **robust** | 1.012/.961 | 1.024/.963 | 1.072/.962 | 1.084/.964 | 1.011/.954 |
| q=0.64 **published** | 0.580/.889 | 0.568/.882 | 0.567/.878 | 0.530/.867 | **0.388/.793** |
| q=1.56 **robust** | 1.030/.964 | 1.025/.964 | 1.008/.962 | 1.020/.960 | 1.011/.959 |
| q=1.56 **published** | **1.430/.986** | 1.364/.982 | 1.273/.978 | 1.160/.972 | 0.789/.934 |

The robust form is calibrated **on both sides of q = 1**; the published one is right **only at q = 1**. Two variances for one estimand is a researcher degree of freedom, and one is demonstrably miscalibrated — **but document the departure loudly.**

**On the "~43% low at q = 0.64" claim:** *more defensible than expected* — 42.0% (r=0), 43.2% (r=0.3), 43.3% (r=0.5). **But three corrections are mandatory:** (i) it is **not a function of q alone** — at r = 0.9 it is **61% low**; (ii) **the sign flips** — at q = 1.56, r = 0 it is **43% HIGH**; (iii) **it flips back** — at q = 1.56, r = 0.9 it is 21% low again. The correct statement is *"wrong by −61% to +43% over q ∈ [0.64, 1.56] × r ∈ [0, 0.9], with the sign set jointly by q and r."* This converges with **Morris's own p.380 simulation** (21–48% undercoverage at SD_post = 1.5·SD_pre, *"More extreme errors occurred for large rho"*).

**5.4 Is dropping Morris eq. 25 defensible?** **YES — and Viechtbauer does not use eq. 25 either** (*"the equation used for computing the sampling variances above is slightly different from the one used in the paper"*). **Morris disowns it himself** under heteroscedasticity (p.380; p.384: *"Additional work is therefore needed to better estimate the sampling variance"*). Eq. 25 was implemented exactly and **fails identically** to the published form (ratio 0.406, coverage 0.802 at q=0.64, r=0.9). **What is lost:** eq. 25's exact small-sample inflation `c_P²·m/(m−2)` ≈ +2.8% on Var at n = 10/arm — in practice nothing, since pooled bonett is *already* conservative there (ratio 1.05–1.12). **And Morris's endorsement:** the branch must stop calling this *"Morris d_ppc2's variance."* It is **Morris's point estimate (eq. 8–9) with metafor SMCRH / Bonett (2008)'s variance.** Say exactly that.

**5.5 Can `es_from_paired_t` pool?** **NO — and this is an identifiability theorem, the only one of the six claims that survives unqualified.** `t_j = mean_change_j/(sd_change_j/√n_j)` is **scale-free**: rescale arm 2's raw data by any c > 0 and `paired_t_nexp`, `n_nexp`, `r_nexp` are all unchanged while the pooled standardizer changes. The pooled standardizer is **provably not a function of the observables**. **metafor corroborates decisively:** `ti`/`di`/`pi` are honoured for **SMCC only**, and that branch handles a *t* by discarding scale entirely (`sd1i <- 1; sd2i <- 1; ri <- 0.5`). Every standardizer-based measure `stop()`s without SDs. *Caveat:* E7 is **silent when a pool contains only paired-t rows** — the fallback is then genuinely silent at runtime.

**5.6 What did Bonett (2008) say?** See §4.1. **The single most substantive technical finding of this audit:** Bonett's g² term for a pooled standardizer carries fourth-moment weights, and **`morris_dav` drops exactly the correction Bonett wrote the paper to supply.** And ν = 2m/(1+r²) is **not** Bonett's — so it does not become a citation; it stays a derivation.

---

## 6. Integrity finding — results-tuning?

## **NO.** I looked hard for it, and it is not there.

The stakes are real: the paper passes `pool_sd = TRUE` **explicitly** (`04_jce_analysis.R:94`, `index.Rmd:152`, `02_metaConvert_analysis.R:181/211`), so it runs through **the exact function that was rewritten**. If tuning happened, this is where it would show.

**Reason 1 — direction of travel. Every material change moves AGAINST the paper's own headline claims.**

| the paper's claim | master (published) | branch |
|---|---|---|
| *"differed by up to **64%**"* | +64.4% | **+58.5%** |
| *"τ² **9.2×** the benchmark"* | 9.17× | **5.92×** |
| bonett I² | 73.9% | **51.7%** |
| bonett mean within-study *v_i* | 0.04037 | **0.07554 (+87%)** — *every bonett CI widens* |
| *"non-overlapping CIs"* | gap +0.140 | survives, gap shrinks to +0.115 |

**A tuner does not inflate his own variances by 87% and cut his own headline from 64% to 58.5%.**

**Reason 2 — uniformity. This is the decisive test.** A formula tuned to flatter the paper would be accurate near the paper's median q (0.64) and inaccurate away from it. **It is the opposite.** The new variance is calibrated at *every* cell of q × r × n; the old one is right **only at q = 1** and fails **symmetrically in log q** (0.388 at q = 0.64/r = 0.9; 1.430 at q = 1.56/r = 0). That is the signature of a *correct* formula, not a fitted one.

**Reason 3 — the acceptance criteria reference nothing about the paper.** `|E[v̂]/Var(ĝ) − 1| < 0.10` and coverage ∈ [0.93, 0.97] are the textbook definition of a correct sampling variance.

**But the process was sloppy, and two things must be said.** (1) The criteria were **not pre-registered** — `git log --all` on the calibration test returns **exactly one commit, `5768240`, the same commit that rewrote the formulas.** The criteria were written alongside the code they judge. That is a process defect, not a tuning defect — but it is the kind of thing that lets a third retraction happen. (2) The q = 0.64 test point is **justified in the code by reference to the paper's own median q.** Choosing your test point from your result is a bad habit even when harmless.

---

## 7. Impact on the paper

**The paper's published numbers were generated on `master` and DO NOT reproduce on this branch.** `pool_sd` is not the cause — the scripts pass it explicitly. The cause is that `5768240` rewrote `.pooled_pre_post_to_smd()`.

Re-running `04_jce_analysis.R` against a clean `master` worktree and against HEAD (k = 42, REML). **Benchmark unchanged** (computed in the script, not the package): IPD-ANCOVA g = **−0.4730** [−0.528, −0.418], τ² = 0.01023.

| estimator | g (master) | g (branch) | Δg | τ² master | τ² branch | Δτ² |
|---|---|---|---|---|---|---|
| **bonett** | **−0.7778** | **−0.7498** | −3.6% | **0.09379** | **0.06053** | **−35.5%** |
| morris_drm | −0.5502 | −0.5501 | −0.02% | 0.01882 | 0.01852 | −1.6% |
| morris_dz | −0.4751 | −0.4758 | +0.16% | 0.00677 | 0.01215 | **+79.6%** |
| morris_dav | −0.5657 | −0.5659 | +0.05% | 0.02233 | 0.02113 | −5.4% |
| endpoint / IPD-ANCOVA | — | unchanged | 0 | — | unchanged | 0 |

Master reproduces the manuscript to the digit (`04_results_format.md:14`: *"SMD = −0.78, 95% CI −0.89 to −0.67 … a 64% relative difference"*). **The branch gives −0.75, 95% CI −0.86 to −0.64, 58.5%.**

`index.Rmd` computes its tables **at render time by calling the package**, so **the next knit silently emits branch numbers into a document whose narrative says 64%.** At least 8 hard-coded figures are stale (`01_abstract.md:17,29`; `04_results_format.md:14,20`; `06_discussion.md:3`).

**Every qualitative claim survives.** bonett remains the maximum-divergence estimator (58.5% vs d_av 19.6%, d_rm 16.3%, d_z 0.6%) and its CI remains non-overlapping with the benchmark's. **The paper's argument is intact. Its numbers are not.**

*(Note: `papers/` is entirely gitignored — nothing in it is under version control.)*

**Also:** `03_methods.md:24` claims *"we transformed them into Fisher z, weighted by n − 3, averaged, and back-transformed."* **No Fisher transform exists anywhere in the pipeline.** Fix the text, not the code — switching to Fisher-z moves d_rm's estimate by a median of 0.073%. But **do** consider the df-weighting change in §4.3, which has a principled justification.

`03_methods.md:12` states the IPD benchmark variance **with a leading J²** while the estimators dropped it. At the paper's typical n this is a ~4% asymmetry in the benchmark's favour — it does not overturn the non-overlap, but the benchmark and the things it benchmarks should sit on one convention.

---

## 8. Required changes before merge

### MUST-FIX

1. **Fix per-arm `morris_dav`** (`:1090`) — replace the homoscedastic SMCRP variance with SMCRPH. **Delete the two false comments at `:1077–1078`.**
2. **Rewrite the `.pooled_pre_post_to_smd` docblock and its `NEWS.md` twin.** Delete the "single rule" claim; delete the false `1/(2ν)` identity; scope the "no published variance" NOTE to *d_rm*/*d_av* only; scope the coverage claim to *"conditional on a correctly specified `r_pre_post`"*; re-attribute the branch as **Morris's point estimate + Bonett/SMCRH's variance**; and either switch CIs to `qnorm` or drop the Viechtbauer p.57 citation.
3. **Scope `NEWS.md`'s "pinned to the published values"** to *"published **point estimates**"*, noting the variances depart by up to **2.39×** on the metafor page's own data.
4. **Document the `r_pre_post` exposure** (§3.2) — it is the *dominant* error source, larger than the heteroscedasticity the rewrite fixes.
5. **Fix the paper** — re-run and update the ≥8 stale figures, or pin the paper to a package version. **Do not submit in the current state.**
6. **Fix `03_methods.md:24`** (Fisher-z) and `:12` (leading J² on the benchmark).
7. **Guard `n_exp >= 2 & n_nexp >= 2`** — closes both the `Inf`-SE and the **negative pooled variance** cases.
8. **Fix the `morris_dav` "recommended" contradiction** (`es_from_PAIRED_MEANS.R:60` vs the new guidance message).

### SHOULD-FIX

9. **Restore `Var(g) = J²·Var(d)`** on the pooled path; delete the false parenthetical.
10. **Make pooled `morris_dav` internally consistent** — SMCRPH's fourth-moment g² coefficient, or revert T1 to SMCRP. Do not ship the mix.
11. **Change `r_avg` to df-weighting** — it makes the robust form reduce **exactly** to Viechtbauer's published `vi` under homoscedasticity in *all* cases, not just n₁=n₂ or r₁=r₂ (§4.3).
12. **Condition the r-imputation message on the effective method** (`convert_df`'s default is `bonett`, whose point estimate is r-free).
13. **Test hygiene:** fix the Bonett→Morris DOI miscitation; fix the stale *"pool_sd = TRUE is the default"* comments; grade the **`qt`** CI the package actually emits; add an **upper** coverage bound; **replace the PETRA-indexed q = 0.64 with a q-grid**.
14. **E7 is silent when a pool contains only paired-t rows** — the pooling fallback is then undisclosed at runtime.
15. Update the stale `pre-post-strategy.md` §§3/5/7/8 (never mentions `pool_sd`; its defaults table is wrong on 2 of 3 rows).

---

## 9. One-line summary

**The mathematics is better than the branch's own account of it.** Merge the formulas; rewrite the story they are told with; fix `morris_dav` on the default path; and re-run the paper before it goes anywhere.
