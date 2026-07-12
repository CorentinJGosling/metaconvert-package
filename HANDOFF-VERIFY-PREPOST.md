# Verification prompt — pre/post SMD rework in metaConvert

Copy everything below the line into a **fresh session**, run from the repo root
(`c:/Users/coren/Documents/maison_juillet/metaConvert`). It is deliberately adversarial:
the goal is to find what is *wrong*, not to confirm what was done.

---

You are a senior statistician performing an **independent, adversarial audit** of a change
set in the R package **metaConvert**, before it ships to production and before an
accompanying paper is submitted. Assume the previous author was competent but
**over-confident, and demonstrably wrong more than once** — he twice made confident claims he
later had to retract (details in §2.3). **Your job is to try to break the work.** Verify
everything from primary sources. Do **not** trust any claim in this prompt, in the code
comments, in the commit messages, or in `NEWS.md` — the author you are auditing wrote all of
them.

## 0. What the package is

metaConvert converts study-level summary statistics into effect sizes for meta-analysis.
`convert_df()` runs every applicable `es_from_*()` converter on each row; `summary()` then
picks the best estimate per row from a hierarchy and can attach quality flags. The area under
audit is the **pre/post (paired) family**: turning pre/post means, mean-change scores, or
paired t/F statistics into a **between-group** SMD.

## 1. The central design question (understand this before reading code)

For a pretest–posttest–control design there are two defensible ways to build a between-group
SMD, and metaConvert supports both via the `pool_sd` argument:

- **`pool_sd = FALSE` — the DEFAULT, and the package's long-standing behaviour.** Compute a
  standardized mean change *within each arm*, subtract them, and **add their sampling
  variances** (valid because the arms are independent). This is **Morris (2008) d_ppc1**, from
  **Becker (1988)**. It is exactly what metafor users do for this design — see
  <https://www.metafor-project.org/doku.php/analyses:morris2008>, where Viechtbauer computes
  `escalc(measure="SMCR")` per arm and then `yi = yT - yC; vi = vT + vC`. He calls it **"more
  broadly applicable"** because it does *not* assume the two arms' true standardizing SDs are
  equal.
- **`pool_sd = TRUE` — opt-in.** Divide the difference in mean change by **one SD pooled across
  arms**: **Morris (2008) d_ppc2** (his eq. 8–9). Morris **recommends** it (more efficient), but
  it *assumes* the arms' true standardizing SDs are equal.

The two coincide when that assumption holds (as randomization implies at baseline) and target
**different estimands** when it fails. So this is a genuine analytic choice, and the package
must not make it silently. **A key thing to check: is that framing right?**

## 2. What was changed (branch `fix/pre-post-family`; `master` untouched)

Diff it yourself: `git diff master..HEAD -- R/` and `git log --oneline master..HEAD`.

### 2.1 The pooled variance formulas were rebuilt (the main substantive change)

All four pooled standardizers — `bonett` (pooled baseline SD), `morris_dz` (pooled change SD),
`morris_drm` (change SD rescaled by `sqrt(2(1-r))`), `morris_dav` (pooled average SD) — now
follow one rule the author claims to have read off metafor's own `escalc()` source:

    Var(g) = Var(numerator)/SD_std^2 + g^2/(2*nu_std),   J = J(nu_std)

with **no leading `J^2`**, and `Var(d)` obtained by substituting `d` for `g` (not `Var(g)/J^2`).
Concretely:
- removed a spurious `2(1-r)` factor from pooled `morris_dz`;
- removed a leading `J^2` from `dz`/`drm`/`dav` (the "LS2" convention; claimed to understate `vi`
  by ~9% at n = 10/arm);
- **`bonett` and `dav` now take the numerator variance from the *empirical pooled change SD*** instead
  of the identity `Var(change) = 2*sigma^2*(1-r)` (which assumes `SD_pre = SD_post` *within* arm).
  The author claims this is metafor's heteroscedasticity-robust `SMCRH`/`SMCRPH` treatment, that at
  `q = SD_pre/SD_post = 0.64` the old form understated `Var(g_bonett)` by ~43% with CI coverage 0.86,
  and that the new form **reduces exactly** to Viechtbauer's published `vi` when `SD_pre = SD_post`;
- `dav`'s `J` and CI now use `nu = 2m/(1+r^2)` (claimed: **Cousineau 2020 eq. 2**) instead of `m`;
- **`bonett` no longer uses Morris eq. 25.** This is the most consequential judgement call in the set.
  **Scrutinise it hardest.**

### 2.2 Other changes
- `es_from_paired_t{,_single_group}()` `morris_dz` variance moved to the metafor `SMCC` convention;
- per-row `pre_post_to_smd` vectors now applied per row (a mixed vector previously fell wholesale into `drm`);
- guards: a zero/non-finite standardizing SD and `|r| >= 1` now yield `NA` instead of `Inf`/`NaN`;
- `verbose` messages for the silent `bonett`→`cooper` coercion and the `r_pre_post = 0.8` imputation;
- two new `[INFO]` cross-row flags (E6 = change-SD rows pooled with raw-score-SD rows; E7 = per-arm
  paired-t rows sharing a pool with pooled rows);
- **the pre/post test suite was moved into `tests/testthat/`** — it previously never ran under
  `R CMD check` (only `test-ANCOVA-MD.R` did).

### 2.3 Claims the author already had to retract — assume more remain

1. He first called the per-arm path (`pool_sd = FALSE`) **a bug**, "not a valid between-group SMD",
   "~61% biased", and **flipped the default to `TRUE`**. That was wrong: it is Becker (1988) /
   Morris d_ppc1, the approach Viechtbauer presents *first*. The "61% bias" figure was **circular** —
   he simulated arms with genuinely unequal SDs, declared the pooled estimator to be the target, and
   measured deviation from it. **The default has since been reverted to `FALSE`.** Check the revert is
   complete and that nothing still assumes pooling.
2. He claimed **"no published sampling variance exists for the pooled two-group form"**. False:
   Viechtbauer publishes one on the metafor-project page,
   `vi = 2*(1-r)*(1/nT + 1/nC) + yi^2/(2*N)`. The author's "robust" formula **departs** from it under
   heteroscedasticity. Decide whether that departure is justified or whether the package should simply
   match the published formula.

## 3. Repository map

**Source (`R/`)**
- `R/internal_multiple_formulas.R` — **the heart of this audit.**
  - `.pre_post_to_smd()` — two-group dispatcher (`pool_sd` switch).
  - `.pooled_pre_post_to_smd()` — **the rewritten function.** Four standardizers.
  - `.single_group_pre_post_to_smd()` — within-arm kernel (4 methods). **Not** rewritten.
  - `.guard_standardizer()`, `.guard_r_pre_post()`, `.d_j()` (Hedges' J).
- `R/es_from_PAIRED_MEANS.R` — `es_from_means_{sd,se,ci}_pre_post()`.
- `R/es_from_PAIRED_MC.R` — `es_from_mean_change_{sd,se,ci,pval}()`.
- `R/es_from_PAIRED_SINGLE_GROUP.R` — the 8 single-group converters.
- `R/es_from_PAIRED_STATISTICS.R` — `es_from_paired_{t,f}()` + `_pval` variants.
- `R/main_convert_df.R` — pipeline, defaults, `verbose` messages.
- `R/internal_flags.R` — quality flags (V18 = baseline/endpoint SD ratio; new E6/E7).
- `R/functions_summary.R` — `summary.metaConvert()`.

**Tests (`tests/testthat/`)** — ~27 files, ~2,300 assertions.
- `test-EXTERNAL-MORRIS-TABLE5.R` — **the strongest anchor.** Reproduces **both** columns of
  Morris (2008) Table 5 (`d_ppc1` and `d_ppc2`) and pins the pooled variance to Viechtbauer's
  published `vi` on the one homoscedastic study.
- `test-pooled-variance-calibration.R` — Monte-Carlo calibration of the four pooled variances +
  a claimed bit-exactness check against metafor. The author claims this file **is** the
  specification for `d_z`/`d_rm`/`d_av`, since no published two-group variance exists for them.
- `test-pool-sd.R`, `test-prepost-fixes.R`, `test-EXTERNAL-PAIRED-{MEANS,MC,SINGLE-GROUP}.R`,
  `test-EXTERNAL-TOSTER.R`, `test-morris-deshon.R`, `test-formula-verification.R`.
- `tests_save/checked/` — 48 further files, **not** wired into `R CMD check`.
- Run: `Sys.setenv(NOT_CRAN="true"); devtools::load_all("."); testthat::test_dir("tests/testthat")`
  (long files gated on `NOT_CRAN`). R is at `/c/Program Files/R/R-4.5.0/bin/x64/Rscript.exe`
  (**4.5.0**, not 4.5.1 as `CLAUDE.md` claims).

**Literature — read the PDFs; do not rely on memory.**
- `data-raw/PRE-POST-SMD/morris2007.pdf` — **Morris (2008)**, ORM 11(2):364-386. Eq. 8-10 = d_ppc2
  point estimate; **eq. 25** = its variance; eq. 12-14 = d_ppc3 (`morris_dav`); **Table 5** = the
  worked example; **p.380 Fig. 3** = behaviour under unequal pre/post SDs; **p.384** = his
  recommendation *against* d_ppc3.
- `data-raw/PRE-POST-SMD/Caldwell_2020_PeerJ_SMD.pdf` — Caldwell & Vigotsky (2020); eq. 13 is
  `d_rm = d_z * sqrt(2(1-r))`.
- `data-raw/Pre post SMD/2020 - QUANT.pdf` — **Cousineau (2020)**, TQMP 16(4):418-421; **eq. 2** is
  the claimed source of the `(1+r^2)` effective df.
- `papers/pre_post_ipd/Literature/2007 VIECHT.pdf` — **Viechtbauer (2007)**, JEBS 32(1):39-60. The
  key source on SMD variance conventions (large-sample vs exact vs unbiased) and CI coverage.
- `papers/pre_post_ipd/Literature/{Shieh,zhang,AHRQ,gibbons1993}.pdf`.
- `data-raw/PRE-POST-SMD/metafor-smd-pre-post.R` — a verbatim dump of metafor's `escalc()` pre/post
  branch (SMCC/SMCR/SMCRH/SMCRP/SMCRPH). **Ground truth for what metafor does.** Also read the live
  source: `deparse(metafor::escalc)`.
- **`https://www.metafor-project.org/doku.php/analyses:morris2008`** — Viechtbauer's worked example.
  **Read this before anything else.** It defines what "the metafor way" is for this design and it is
  what the users of this package will compare against.
- `data-raw/PRE-POST-SMD/{pre-post-strategy.md,pre-post-testing.md}` — the team's own design and
  validation notes. **PARTLY STALE** (they still describe the old architecture); treat as evidence of
  intent, not of current behaviour. §3 and §8 of `pre-post-strategy.md` are wrong under
  `pool_sd = TRUE` and need rewriting.
- `data-raw/internal_multiple_formulas2.R` (~lines 874-886) — carries a comment marking a variance
  formula "WRONG … corrected formula by W Viechtbauer".
- **`C:/Users/coren/Downloads/harrer.pdf`** — Harrer et al. (2025), PLOS Ment Health 2(7):e0000347.
  The code *used to* cite "Harrer eq. 13/14"; those comments were deleted as unverifiable. **You now
  have the PDF: check eqs. 12-14 yourself** and decide whether the deletion was right and whether
  anything should cite it. (Their SMD_CS/CS, SMD_CS/BL, SMD_CS/EP variants map onto `morris_dz`,
  `bonett`, and an unimplemented hybrid.)
- **`C:/Users/coren/Downloads/ostinelli.pdf`** — Ostinelli et al. (2024), Res Synth Methods
  15:758-768, on combining endpoint and change data.

**The team's own paper (`papers/pre_post_ipd/`)** — an IPD-ANCOVA benchmark of these very estimators
(target: JCE). See `index.Rmd`, `03_methods.md`, `05_supplement.md`, `scripts/` (esp.
`04_jce_analysis.R`, `01_extract_functions.R`), and `REVIEW-critical-appraisal.md` (an existing
internal peer review). **`index.Rmd` pools with `rma(yi = g, sei = res$g_se)`, so every pooled
estimate, CI and tau^2 in the paper depends on the variance formulas under audit.** Headline claims
include: bonett ~64% higher than the IPD benchmark "with non-overlapping CIs"; bonett tau^2 = 0.094
(9.2× the benchmark); median baseline/endpoint SD ratio q = 0.64.

## 4. Your tasks

**A. Verify the mathematics from primary sources.** For each pooled branch, derive the sampling
variance yourself and check the code. Confirm or refute, with page/equation citations:
- Is the "no leading `J^2`" rule actually what metafor does? (Read the `escalc` source.)
- Is `(1+r^2)/(4N)` really `1/(2*nu)` with `nu = 2m/(1+r^2)`? Is Cousineau (2020) eq. 2 really its
  source, and does his single-arm df legitimately add across two arms the way the code assumes?
- **Is replacing Viechtbauer's published `vi` (`2(1-r)(1/nT+1/nC) + yi^2/(2N)`) with the
  "robust" empirical-change-SD form justified?** Verify the author's two claims: (i) that the robust
  form reduces *exactly* to Viechtbauer's when `SD_pre = SD_post`; (ii) that Viechtbauer's form is
  ~43% low at `q = 0.64`. Then decide: should the package match the published formula, use the robust
  one, or offer both? Note Morris himself (p.380, Fig. 3) reports his variances degrade under unequal
  pre/post SDs.
- Is dropping **Morris eq. 25** for `bonett` defensible? What exactly is lost?
- **`es_from_paired_t` cannot pool** (a paired t identifies each arm's `mean_change/sd_change` but not
  the ratio of the arms' SDs). Is that reasoning correct?
- **Read Bonett (2008)**, *"Confidence intervals for standardized linear contrasts of means"*, Psych
  Methods 13(2):99-109 — **not in the repo; obtain it.** "Linear contrasts" is multi-group language.
  If Bonett already pools effective df across arms, `nu = 2m/(1+r^2)` stops being the author's
  derivation and becomes a citation. **Highest-value open lookup.**

**B. Re-run the Monte Carlo INDEPENDENTLY.** Do **not** reuse the author's scripts (they are not in
the repo). Write your own: draw bivariate-normal (pre, post) data for two arms, compute each
estimator and its variance **through the actual package functions**, and check (i) `E[vi]` against the
true sampling variance of `g` and (ii) 95% CI coverage. Sweep `n`, `r`, effect size, **and** the SD
ratio `q = SD_pre/SD_post` (include `q = 1`, `0.64`, `1.56`). **Do this for `pool_sd = FALSE` too** —
it is the default, and its calibration matters more than the pooled path's.

**C. Integrity check — look for results-tuning.** The author had access to `papers/pre_post_ipd`,
whose headline findings depend on these variances. Determine whether any formula choice was made to
flatter that dataset:
- The heteroscedastic MC was run at `q = 0.64` — PETRA's median. Is the chosen formula better *only*
  near `q = 0.64`, or uniformly across `q`? (Uniform ⇒ not tuned.)
- Which direction do the changes push the paper's headline claims (bonett's 64% gap "with
  non-overlapping CIs"; bonett tau^2 = 9.2× benchmark)? A change that *weakens* the paper's claims is
  evidence *against* tuning; one that strengthens them deserves suspicion.
- Were the criteria (unbiased `E[vi]`, nominal coverage) data-independent, or reverse-engineered?
Report your finding plainly either way.

**D. Quantify the impact on the paper.** Re-run `papers/pre_post_ipd/scripts/` (at minimum
`04_jce_analysis.R`) against the new code; report how the pooled estimates, CIs and tau^2 move. Check
whether the scripts pass `pool_sd` explicitly. Also check `03_methods.md:12`, which states the IPD
benchmark variance **with a leading `J^2`** — if the estimators dropped it, the benchmark now sits on
a different convention from the things it benchmarks.

**E. Audit the test suite for circularity.** Many tests "verify" a formula by re-typing the source
expression; those verify nothing. Identify every remaining transcription and say which tests would
actually catch a regression. Verify the claimed **bit-exact** agreement between pooled `morris_dz`
and `metafor::escalc(measure = "SMD", vtype = "LS")` on change scores. **Also: derive the same effect
size from every input format (pre/post means, mean change, paired t, SE, CI, p-value) and confirm they
all agree** — a cross-format equivalence failure would be a serious defect.

**F. Should the single-group branches change too?** They were left alone. Single-group `morris_drm`
still uses the LS2 pattern (`var_g = J^2 * var_d`) while the other three use LS. The team documented
that as deliberate (Morris & DeShon specify `var(d)`). Is it defensible now that the pooled branches
all use LS? Note the team's own list of gaps says *"cooper g_se never externally validated"*.

## 5. Known-open questions — do not let these slide

1. **No published sampling variance exists for the pooled two-group `d_z`, `d_rm` or `d_av`** (only
   `bonett`/d_ppc2 has one). The implemented formulas are delta-method derivations by the author,
   backed only by his own MC. Is the derivation sound? Is MC calibration adequate evidence?
2. Bonett (2008) has not been read (see A).
3. `r_avg` is an n-weighted mean of **raw** correlations, not Fisher-z, not df-weighted. No source
   prescribes a rule. It now enters `J` (via `nu`) as well as the variance.
   `papers/.../03_methods.md:24` claims Fisher-z pooling the code does not do — one must change.
4. A homoscedasticity assumption still sits in `morris_drm`'s `2(1-r)` rescaling and in `dav`'s
   `(1+r^2)` df, even after the numerator was made robust. How much does that matter at `q = 0.64`?
   metafor offers `SMCRPH` for the heteroscedastic d_av case.
5. `pre-post-strategy.md` §3/§8 and `pre-post-testing.md` are stale and partly false under
   `pool_sd = TRUE`; the package's "validated against metafor" claim does not cover the pooled path
   (metafor has **no** pooled two-group pre/post measure).
6. The E6/E7 flags and the new `verbose` messages are the author's own inventions — check they fire
   when they should and stay silent otherwise (especially: no spurious messages on datasets with no
   pre/post data).

## 6. Deliverable

A written verdict: what is correct, what is wrong, what is unverifiable. State clearly whether this
branch is safe to merge and list every change you would require first. **Show your working** —
formulas, page citations, and your own Monte-Carlo numbers.
