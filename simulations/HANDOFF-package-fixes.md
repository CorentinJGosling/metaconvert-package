# Handoff prompt — metaConvert package fixes

Paste everything below the line into a fresh Claude Code session in
`c:\Users\coren\Documents\sideprojet\metaConvert`. It is self-contained.

When it finishes, paste its **REPORT** section back into the original session for
cross-checking. The report format is specified so the claims can be checked against
`git diff` and re-run rather than taken on trust.

---

You are making a small number of targeted fixes to the metaConvert R package. Each
fix has already been verified against its source paper by a prior session, and the
findings are written up in `simulations/SOURCE-VERIFICATION.md` — **read that file
first**, especially §2 (confirmed defects), §4 (claims already refuted — do not
re-raise them), and §7 (documentation fixes).

## Non-negotiable ground rules

1. **The prior on any "the package is wrong" hunch is that the package is RIGHT.**
   In the previous session five such hunches were investigated and *four were false*.
   One of them (a d/g slot "fix") changed shipped numbers for the worse and had to be
   reverted after simulation disproved it. If you form a new suspicion beyond the
   listed fixes, do **not** change code — write it in the report as a question.
2. **No behaviour change without evidence.** Fixes 1 and 2 change user-facing numbers
   and that is intended. Everything else must be provably output-identical.
3. **Verify numerically, not by eye.** Print ≥8 significant digits, before and after.
4. **Do not edit `tests_save/checked/`.** Those archived tests encode deliberate
   conventions and have already caught one bad change. If one fails, that is a signal
   your fix is wrong — investigate, do not "fix" the test.
5. R is at `/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe`. Load the dev source with
   `devtools::load_all(".")`, never an installed build. Internals: `metaConvert:::.name`.
6. PDF text layers in `data-raw/` **drop minus signs** and mangle fractions. Never
   conclude anything from a sign that could be an extraction artefact.

## Baseline to record before touching anything

Run and save the output — you will need it for before/after comparison:

```r
suppressPackageStartupMessages({library(testthat); devtools::load_all(".", quiet=TRUE)})
setwd("tests_save/checked"); if (file.exists("setup.R")) source("setup.R")
for (f in c("test-2x2.R","test-ES-OR.R","test-ES-RR.R","test-ES-COR.R","test-CORPB.R",
            "test-PHI-CHISQ.R","test-ES-SMD.R","test-ES-USER-CONVERSION.R",
            "test-ES-NNT.R","test-ES-NNT-RD.R","test-ANCOVA-MEANS-F-T-pval.R")) {
  r <- try(as.data.frame(test_file(f, reporter="silent")), silent=TRUE)
  if (inherits(r,"try-error")) { cat(sprintf("%-32s ERROR\n", f)); next }
  cat(sprintf("%-32s pass=%-6d fail=%d\n", f, sum(r$passed), sum(r$failed)+sum(r$error)))
}
```
Also run the current `tests/testthat/` suite (11 files) and record pass/fail.

---

## FIX 1 — `reverse_2x2` inverts the tetrachoric CI instead of reflecting it

**Locations:** `R/internal_multiple_formulas.R:545-550` (`.contingency_to_cor`) and
`:742-747` (`.phi_to_cor`, tetrachoric branch).

**Defect:** CI bounds are negated in place with no swap, so `lo > up` and the point
estimate falls outside its own interval.

```
es_from_2x2(143,52,41,164)                     r =  0.7458746  CI ( 0.6591207,  0.8326284)
es_from_2x2(143,52,41,164, reverse_2x2 = TRUE)  r = -0.7458746  CI (-0.6591207, -0.8326284)
                                                        correct CI (-0.8326284, -0.6591207)
```

**The fix already exists elsewhere in the package** — the identical bug on the OR path
was fixed with a comment explaining exactly this at `R/es_from_stand_OR.R:351-353`
("Swapping alone left inverted, wrong-signed bounds."). Apply the same negate-and-swap
to both the r-scale and z-scale bounds.

**Verify:** for a grid of tables and both `reverse_2x2` values, assert
`ci_lo <= est & est <= ci_up` and that reversing maps `(lo, up) -> (-up, -lo)` exactly.
Do the same for `phi_to_cor = "tetrachoric"` via whichever exported function reaches it.
Add a test in `tests/testthat/test-reverse-ci-reflection.R`.

## FIX 2 — `or_to_cor = "bonett"` silently degrades to `lipsey_cooper`

**Locations:** `R/es_from_stand_OR.R:328-339` (the `nn_miss` gate); the value left in
place comes from the `.es_from_d` call at `:247-251`.

**Defect:** Bonett's `c` needs `small_margin_prop`, a user-entered column with **no
auto-derivation anywhere** (grep confirms: only `data_extraction.R:400`,
`internal_check_data.R:147`, and pass-through). When it is blank — the normal case —
the row falls out of the gate and the earlier `lipsey_cooper` r survives.
`convert_df()` defaults to `or_to_cor = "bonett"` (`R/main_convert_df.R:279`), so the
**default method usually does not run**:

```
or=2, logor_se=0.2, n_exp=50, n_nexp=50, n_cases=40, n_controls=60, measure="r"
  currently returns r = 0.188   (= lipsey_cooper, 0.1876806337)
  Bonett with small_margin_prop = 0.4: r = 0.259  (0.25860082)
```
38% understatement, `info_used` still reports only `"or_se"`, no flag.

**Decide between two fixes and say which you chose and why:**
- (a) **derive** `small_margin_prop` when the margins are available. Bonett & Price
  2005 p.216 define `pmin` as *the smallest marginal proportion*, so with
  `n_exp`, `n_nexp`, `n_cases`, `n_controls`, `n_sample` present it is computable as
  `min(n_exp, n_nexp, n_cases, n_controls)/n_sample`. Verify against the paper's three
  worked examples before trusting this.
- (b) **fall back explicitly**: keep `lipsey_cooper` but record it (e.g. in
  `info_used`) and raise a flag, so the user knows the requested method did not run.

(a) is better if and only if the derivation is exactly what the paper means. Check it.
Do **not** do both. Do not change the default.

**Verify:** reproduce all three Bonett & Price (2005) worked examples (Examples 1-3,
pp. 220-222) through the package. Note the prior session's finding: metaConvert
computes `c` from **uncorrected** marginal proportions and thereby matches the paper's
*printed* numbers (Example 2: rho*=0.83115717, se=0.10765128, CI (0.48766686,
0.95605983) vs printed .831/.108/(.488,.956)), whereas the paper's own prose ("after
0.5 has been added to each cell") would give .835/.107/(.492,.958). **Preserve that
behaviour** — do not "correct" it toward the prose.

## FIX 3 — documentation only, no behaviour change

All verified as code-correct / docs-wrong. Change only comments and Rd/roxygen.

| location | fix |
|---|---|
| `R/es_from_stand_OR.R:113` | Rd LaTeX for Bonett's `c` puts `1 −` inside the `/5` numerator and uses raw **counts** where paper and code use **proportions**. Code at `internal_multiple_formulas.R:602` is B&P eq (3) exactly: `(1 - abs(n_exp/n_sample - n_cases/n_sample)/5 - (1/2 - small_margin_prop)^2)/2`. Make the Rd match the code. |
| `R/es_from_stand_OR.R:103-104` | Rd for the digby z CI omits a square: shows `sqrt((c^2/4) * logor_se)`, code has `logor_se^2`. |
| `or_to_rr = "grant"` docs/`@references` | The estimator is **Zhang & Yu (1998)**, *JAMA* 280:1690-1691. Confirmed: Grant 2014 *BMJ* 348:f7450 prints the same expression and claims no originality (he cites Prasad 2008 / Shrier 2006). Cite Zhang & Yu for the formula, Grant for the communication framing. Do **not** add a separate `zhang_yu` method — it is one estimator. |
| `R/es_from_stand_OR.R` `@details` for the imputed OR SE | Says "the variance of the OR is equal to the **mean of the standard error** of all possible situations"; the code takes the mean of the **variance** then square-roots (`internal_multiple_formulas.R:322`). One-word correction: "mean of the variance". No behaviour change. |
| `R/internal_multiple_formulas.R:1233-1234` | Comment claims single-group `morris_dav` reproduces Bonett 2008 Example 2's `var = 0.0148` "to the printed digit". Package gives `0.0146674131`. Bonett's own arithmetic gives `0.0147449109`, so **0.0148 is his rounding slip** — the formula identity is exact (0 diff vs `metafor SMCRPH`). Reword to state the formula identity, not the printed number. |
| `.pooled_pre_post_to_smd` docblock | Attributes the whole pooled `d_av` variance to Bonett 2008 eq 19; only the `g^2` fourth-moment coefficient is eq 19's. The leading term is metafor's n-based pooled form and departs from eq 19 by up to **38% at strongly unequal arm sizes**. State this. |
| `data-raw/PRE-POST-SMD/pre-post-testing.md` | No OUTDATED banner and stale: maps `morris_dav` to metafor `SMCRP` when the code implements `SMCRPH` (4.2% apart). Add a banner or correct it. (`pre-post-strategy.md` already has a banner.) |
| `data-raw/Cohen 1988.pdf` | **Not Cohen 1988** — it is Aaron, Kromrey & Ferron (1998), "Equating r-based and d-based Effect Size Indices" (ERIC ED433353). Zero occurrences of "covarian". Rename the file. |

## FIX 4 — silent-NA warnings (only if cheap and clearly safe)

Both are missing-input paths that currently return NA with no explanation.
- blank `n_cov_ancova` silently yields `g = NA`
- blank `baseline_risk` in the Di Pietrantonj root choice silently takes the smaller
  (rare-event) root

A one-line warning each. **If adding a warning breaks any archived test** (several
capture warnings), stop and report instead — do not suppress the test.

## FIX 5 — the ANCOVA unbiasedness note is false on 8 routes (wording only)

**Confirmed by simulation, twice, independently.** No code fix is possible or wanted.

The reported ANCOVA `t` embeds the *exact* SE
`s_res * sqrt(1/n_exp + 1/n_nexp + D)`, `D = (x̄_exp − x̄_nexp)²/SS_x`
(Lai & Kelley 2012 eq 5). `R/es_from_ANCOVA_statistics.R:61` back-derives

```r
d <- ancova_t * sqrt(1/n_exp + 1/n_nexp) * sqrt(1 - cov_outcome_r^2)
```

which recovers `MD_adj` only if `D = 0`. So **d is attenuated by
`1/sqrt(1 + D/(1/n_exp + 1/n_nexp))` ≈ `1/sqrt(1 + δ_x²/4)`** for equal arms. Verified
as an algebraic identity to 15 digits, and by Monte Carlo (two independent runs, one at
40,000 reps, one at 4,000 with a different DGP and seed):

| δ_x | adjusted-means route | t route | ratio | predicted |
|---|---|---|---|---|
| 0.0 | +0.4% | +0.1% | 0.9975 | 1.0000 |
| 0.5 | +1.2% | **−2.1%** | 0.9671 | 0.9701 |
| 1.0 | +0.6% | **−10.2%** | 0.8925 | 0.8944 |

The means route's small residual is d's ordinary upward small-sample bias, flat in δ_x.
The t route's is monotone in covariate imbalance. It also persists at δ_x = 0 for small
n (E[D] > 0 even under exact balance) and is unchanged whether the true or sample `r`
is supplied (`sqrt(1-r²)` cancels).

**Affected routes (8), verified from identical `lm()` fits at δ_x = 1:**
`ancova_t`, `ancova_f`, `ancova_t_pval`, `etasq_adj` (partial η² = t²/(t²+df) exactly),
`ancova_md_se`, `ancova_md_ci`, `ancova_md_pval` (mirror-image: `ancova_md_sd <-
ancova_md_se / sqrt(1/n_exp + 1/n_nexp)` *inflates the standardizer* by the same
factor), and `ancova_means_se`, which is attenuated LESS than the others (0.942 vs 0.893 at delta_x = 1) because it sees only per-arm leverage, never the combined SE.

**Clean routes:** `es_from_ancova_means_sd`, `..._pooled_adj`, `ancova_md_sd`. The
`@note` at `R/es_from_ANCOVA_means.R:62` is therefore **correct where it sits** — do
not change it.

**Important context for the wording:** metaConvert is **not** deviating from the
literature standard. `compute.es::a.tes/a.fes/a.pes` implement the same Cooper table
12.3 formula and share the limitation, which is why
`tests_save/checked/test-ANCOVA-MEANS-F-T-pval.R` passes — it is a cross-package
agreement test on the formula, not a validity test against raw data. Only the
*unbiasedness claim* is false.

Also note **what is not broken**: d and its SE shrink by the *same* factor, so the
implied z, the p-value, and "does the CI exclude zero" are all exactly right (verified:
z bias +0.000% at every δ_x). What is wrong is the magnitude (~10% understated at
δ_x = 1) and the inverse-variance weight (~+26% excess), and it does not wash out —
fixed-effect pooled bias −10.8% on the t route vs −0.1% on the means route.

**What to do:** amend the `@note` on `R/es_from_ANCOVA_statistics.R:36-43` and mirror
it onto `es_from_etasq_adj` (`R/es_from_ETASQ.R:134`), and add it to the three
contaminated `ancova_md_*` variants and `ancova_means_se`/`_ci`, which currently carry
no note at all. State: the point estimate is *not* unbiased under covariate imbalance
on these routes; give the attenuation factor and the ≈3% / ≈11% figures at δ_x = 0.5 / 1.0;
note it is negligible in randomised designs; note the p-value is unaffected but the
magnitude and the pooling weight are not; and point users to
`es_from_ancova_means_sd()` / `es_from_ancova_md_sd()` when adjusted means or the MD
with residual SD are available.

Also fix **`CLAUDE.md` line 87**, which says the omitted term is unrecoverable "and the
point estimate stays unbiased" — scope that clause to the means/MD-with-SD routes.

Leave the hierarchy alone: `R/main_convert_df.R:1372` already ranks
`means_adjusted_list_L19` above `ancova_adjusted_list_L18`, so the attenuated routes are
selected only when adjusted means are unavailable. That is the right behaviour.

> **Both halves of this section have since been superseded (NEWS.md 2.0.1).**
>
> *The `@note`*: it was later shortened, then restored in condensed form. `?es_from_ancova_t`
> is again the canonical statement and still carries the `1/sqrt(1 + imbalance²/4)` factor and
> the "about 3% / about 11%" anchors that `students/tasks/task-07` sends RAs to find; the
> sibling routes now point at it rather than repeating it.
>
> *The hierarchy*: this verdict was wrong, and this document's own line 177 shows why — it
> lists `ancova_md_sd` among the unattenuated routes, yet `md_adjusted_list_L20` sat **below**
> `ancova_adjusted_list_L18`. A study reporting an adjusted MD with its residual SD *and* an
> ANCOVA F therefore got the attenuated estimate (g = 0.996 selected over an available 1.084
> at δ_x = 1). `L20` is now split: `es_ancova_md_sd` alone moves above `L18`, while
> `md_se`/`md_ci`/`md_pval` stay below it, since those must invert a reported SE/CI/p-value
> and are in the same attenuation tier as the test statistics.

## FIX 6 — tetrachoric CI bounds escape [-1,1] silently (do OPTION 1 only)

**Read this whole section before touching anything. The original diagnosis was
mis-attributed and its prescribed fix is wrong. A skeptic pass established what is
actually defective and what the correct minimal fix is.**

### What is confirmed
`internal_multiple_formulas.R:574-575` builds `r_lo/r_up <- r ± qnorm(.975)*sqrt(vr)`.
On Bonett & Price 2005's own Example 2, `f = (4,6,1,89)`:
`r = 0.8641990826`, `r_ci = (0.6503674747, 1.07803069)` — **upper bound > 1**.
Out-of-bounds frequency: 11.8% of random tables, 27.7% under dichotomised
bivariate-normal sampling, 100% of a systematic n=50 scan with small off-diagonal
cells. And it is **not flagged**: `internal_flags.R:1467` (B1) tests `abs(es) > 1` on
the **point estimate only**, so `summary(..., flags = TRUE)` prints
`r = 0.864 [0.650, 1.078]` with an empty flag string.

### What was wrong in the original diagnosis — do not repeat it
- **B&P 2005 is not the source of this route.** `man/es_from_2x2.Rd` cites
  Cooper/Cochrane/Lipsey/Sedgwick/Altman/Wen and the roxygen says it "relied on the
  implementation of the formulas of the 'metafor' package". B&P 2005 *is* cited at
  `R/es_from_stand_OR.R:151` for `or_to_cor = "bonett"` — and **that** route already
  implements B&P eqs 8a/8b faithfully (`r_lo <- cos(pi/(1 + or_ci_lo^c))`). The package
  already ships B&P's recommended interval where B&P's estimator is used.
- **metafor does the identical thing.** `summary(escalc(measure = "RTET"))` returns
  `ci.lb = 0.6504, ci.ub = 1.0780` — bit-identical.
- **It is a tested contract.** `tests_save/checked/test-2x2.R` asserts
  `expect_equal(es_ci_lo_crude, comp_res_r$ci.lb, tolerance = 1e-10)` and the same for
  `ci.ub`. Deliberate metafor parity, not a slip.
- **B&P Example 1 reproduces exactly**: metaConvert gives
  `0.3350968379 / 0.0481943695 / (0.2406376094, 0.4295560664)` = the paper's printed
  SAS `.335 / .048 / (.240, .429)`.
- **Switching to eqs 8a/8b would be wrong**: they belong to a *different estimator*
  (ρ̂* = .835 vs the ML ρ̂ = .864), and pairing them with the ML point estimate would
  decouple the reported `r_se` from the CI. B&P also state (p.215) that
  "normalizing and variance-stabilizing transformations have not been developed for
  the tetrachoric correlation", which undercuts the transformation-based fix.

### Why it is still a genuine defect — the argument is INTERNAL
metaConvert already flags exactly this pattern for alpha and ICC: `internal_flags.R:1576-1585`
(B7b) and `:1600-1618` (B8b) fire `[UNUSUAL]` when a raw-scale Wald bound escapes the
parameter space while the point estimate is valid, with the message noting the
symmetric Wald interval cannot be trusted there. `r` is bounded identically, built
identically, and has **no analogue**. So "metafor does it too" is not a defence — it is
an inconsistency against metaConvert's own published quality standard, and the
indefensible part is the *silent* pass-through.

(For context, a coverage simulation over 12 cells × 400 reps found the current interval
the worst of four candidates — mean 0.9236, worst 0.8700 — independently reproducing
B&P's reported .921/.821. So the statistical criticism is real. But see below.)

### DO THIS — Option 1 only
Add the `r` analogue of B7b/B8b in `R/internal_flags.R`: flag `[UNUSUAL]` when a
correlation CI bound leaves [-1,1] while `|es| <= 1`, with a message pointing out that
the symmetric Wald interval on the r scale can escape the parameter space near |r| = 1.
Purely additive. Breaks no test. Preserves metafor parity. Makes the flag framework
self-consistent. **This alone fixes the indefensible part.**

Verify: `f = (4,6,1,89)` must now raise the flag; a well-behaved table must not; run the
`test-flag-*` files plus `test-2x2.R` and confirm counts unchanged apart from your new
assertions. Add `tests/testthat/test-r-ci-out-of-bounds.R`.

### DO NOT DO Option 2 in this session
Replacing the r CI with `tanh(z_ci)` (already computed at `:572-573`) has the best
coverage of the four candidates (0.9606 / worst 0.9325), but it: breaks the two metafor
parity assertions in `tests_save/checked/test-2x2.R`; is grossly asymmetric
(0.4300 below vs 0.1093 above the estimate) with halfwidth/se = 2.4714 instead of 1.9600,
so **A6 would start firing on every tetrachoric row** unless a construction-aware
exception is added alongside the existing Glass/md/clamped-proportion ones; and needs a
documentation note that the r CI is no longer metafor's. That is a design decision for
the author, not a drive-by fix. **Report it as a recommendation with these numbers and
stop.**

## DO NOT TOUCH
- `viechtbauer` d/g slot assignment — investigated, the shipped convention is correct,
  and a "fix" was reverted after simulation. See `SOURCE-VERIFICATION.md` §4.
- `unit_type = "raw_scale"` arithmetic — correct; docs-only issue already recorded.
- `cor_to_smd` unequal-groups generalisation — contested in the literature
  (Pustejovsky argues against it for experimental designs). Not a bug.
- The `simulations/` folder. Out of scope for this session.

---

## REPORT (paste this back, filled in)

```
### Baseline
tests_save/checked (11 suites):   pass=____  fail=____
tests/testthat (11 files):        pass=____  fail=____

### FIX 1 reverse_2x2 CI reflection
files+lines changed:
before: es_from_2x2(143,52,41,164, reverse_2x2=TRUE) CI = (________, ________)
after :                                              CI = (________, ________)
invariant ci_lo <= est <= ci_up across grid:  PASS/FAIL  (n tables = ____)
reversal maps (lo,up) -> (-up,-lo) exactly:   PASS/FAIL
phi_to_cor="tetrachoric" also fixed:          YES/NO   (how verified: ______)
new test file + assertion count:              ____

### FIX 2 or_to_cor="bonett"
option chosen: (a) derive small_margin_prop / (b) explicit fallback
why:
B&P 2005 worked examples reproduced through the package:
  Ex1  rho*=________  se=________  CI=(________, ________)   paper: ____________
  Ex2  rho*=________  se=________  CI=(________, ________)   paper: .831/.108/(.488,.956)
  Ex3  rho*=________  se=________  CI=(________, ________)   paper: ____________
uncorrected-margins behaviour preserved (Ex2 gives .831 not .835):  YES/NO
the motivating case now returns r = ________  (was 0.188, Bonett = 0.259)
info_used / flag when the method cannot run:  ____________

### FIX 3 documentation
one line per item: location -> what changed.  Confirm zero behaviour change:
  method to prove it: ____________   result: ____________

### FIX 4 warnings
added / skipped, and why:

### FIX 6 tetrachoric out-of-bounds CI flag (Option 1 ONLY)
flag added at internal_flags.R:____
f=(4,6,1,89) now flagged:                     YES/NO   message: ____________
well-behaved table NOT flagged:               YES/NO
test-2x2.R metafor parity assertions intact:  YES/NO   (must be YES)
new test file + assertion count:              ____
Option 2 attempted:                           MUST BE NO

### FIX 5 ANCOVA note wording
files+lines amended (expect ~5 files):
CLAUDE.md line 87 scoped:                     YES/NO
es_from_ANCOVA_means.R:62 left unchanged:     YES/NO  (it is correct — must stay)
confirm zero behaviour change:                method: ________  result: ________

### Regression after all fixes
tests_save/checked (11 suites):   pass=____  fail=____   (baseline: ____/____)
tests/testthat:                   pass=____  fail=____   (baseline: ____/____)
any suite whose count CHANGED, and the explanation:

### Questions / suspicions NOT acted on
(anything you thought was wrong but did not change — with the evidence)

### git diff --stat
(paste)
```

Report honestly. If a fix did not work or you are unsure, say so — an accurate
"blocked" is far more useful than a confident wrong claim. Do not report a test as
passing unless you saw it pass.
