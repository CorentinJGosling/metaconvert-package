# Audit fix tracker — 2026-07-16 formula audit

Source: 8-reviewer multi-agent audit + adversarial verification (findings digest:
`c:/tmp/audit_digest.md`, keyed `[N]`). Every fix of a MAJOR/MODERATE problem must land
together with a test that FAILS on the pre-fix behaviour.

Status: `[ ]` open · `[x]` fixed+tested · `[~]` partially addressed (see note) · `[u]` user decision

## A — Wrong numbers (blockers)

- [x] **P1** [8] MAJOR (default path): `es_from_etasq_adj` applies the marginal-scale ANCOVA
  variance (`adjusted=TRUE`) to a residual-scale point estimate `d = 2*sqrt(etasq_adj/(1-etasq_adj))`.
  Emitted SE understates its own estimator's sampling SD by ~45% at r=0.7 (sim coverage 0.712).
  Cross-route contradiction with `es_from_ancova_f` reproduced exactly (d ratio 1.2006 =
  `1/sqrt(1-r^2)*sqrt(N/df_err)`). Fix: convert η²ₐ → implied F → route through the ancova_f
  pathway (marginal-scale d, coherent variance). File: `R/es_from_ETASQ.R`.
- [x] **P2** [14] MAJOR (default path): ICC(2,1) `icc_type="agreement"` SE understates sampling
  variance whenever between-rater variance > 0 (MC coverage 0.74–0.76, worsens with n; ICC(2,1)
  is not √n-consistent for fixed k while the SE shrinks 1/√n). Not fully fixable from summary
  data. Fix: correct the roxygen claim ("same at leading order" is false with rater variance),
  fix vignette line ("SE is computed correctly for each study based on its icc_type" — icc_type
  has zero computational effect), add an informational flag on agreement rows. Files:
  `R/es_from_ICC.R`, `vignettes/Psychometrics.Rmd`, `R/internal_flags.R`.
- [x] **P3** [20] MAJOR (opt-in `cor_to_smd="cooper"`): Spearman → d/g/logOR SEs omit the
  Spearman delta correction (5–15% understated) while r/z in the same row carry it; the
  justifying comment ("cooper … does not route through r_se") is factually false — the cooper
  branch is `sqrt(4*r_se^2/(1-r^2)^3)`, linear in r_se. Fix: include "cooper" in the delta-scale
  override; correct comment. File: `R/es_from_spearman_COR.R`.
- [x] **P4** [0] MAJOR (opt-in `smd_var="hedges_olkin"`): the raw-MD family
  (`es_from_md_sd/md_se/md_ci/md_pval`) hard-codes the LS2 d_se and accepts no `smd_var`, so a
  mixed dataset silently mixes variance conventions (contradicts NEWS "every es_from_*() that
  builds a d/g from the large-sample formula"). Fix: add + thread `smd_var`. Files:
  `R/es_from_stand_MD.R`, `R/main_convert_df.R`.
- [x] **P5** [23] MINOR (numeric, conservative direction): boundary (p=0/1) proportion SEs use
  raw n instead of the continuity-corrected n+1; deviates from the metafor PR/PLO convention the
  roxygen claims to match by exactly `sqrt((n+1)/n)`. Fix: corrected denominator. File:
  `R/es_from_PROP.R`.

## B — Package contradicting itself (false flags, crash, silent mixing)

- [x] **P6** [2][9] MODERATE (default path): A6 CI-width check flags the package's own new
  Welch-df md CI as `[DISCORDANT]` (expects `qt(.975, N-2)`). Fix: accept implied critical
  values in `[qt(.975,N-2), qt(.975,min(n1,n2)-1)]` for package-computed md rows. File:
  `R/internal_flags.R`.
- [x] **P7** [3] MODERATE (opt-in glass): A6 flags Glass rows (CIs correctly on
  `qt(.975, n_nexp-1)` per SMD1). Fix: store per-row `smd_denom` as attribute in `convert_df`,
  thread into flags, use control df for glass rows. Files: `R/main_convert_df.R`,
  `R/functions_summary.R`, `R/internal_flags.R`.
- [x] **P8** [21] MODERATE (default path): A6 flags every clamped raw-proportion CI (bound at
  0/1 — all boundary rows at any n, and p ≲ 4/n). Fix: compare against the [0,1]-clipped
  expected width when the clamp signature is present. Files: `R/internal_flags.R`,
  `R/functions_summary.R`.
- [x] **P9** [22] MODERATE: two surviving hard `stop()`s in `es_from_PROP.R` (prop outside [0,1];
  n_cases > n_sample) abort the ENTIRE `convert_df()` under `correct_inputs=FALSE` — for any
  measure, killing all valid rows. Fix: per-row NA + warning (centralised-validation convention).
  File: `R/es_from_PROP.R`.
- [x] **P10** [10] MODERATE: `convert_df(smd_denom="glass")` silently mixes Glass and pooled-SD
  estimands across methods (only means_sd/se/ci honour it). Fix: one-time scoping message +
  `@param` doc. File: `R/main_convert_df.R`.
- [x] **P11** [1] MODERATE (doc): the default-path md CI change (pooled df → Welch-Satterthwaite)
  is deliberate and tested but absent from NEWS. Fix: NEWS behaviour-change entry. File: `NEWS.md`.

## C — Psychometric API robustness

- [x] **P12** [28] MODERATE: `es_disattenuate` scalar `n_sample` + vector `r` silently NAs SEs
  for rows 2..n (missing recycling). File: `R/es_disattenuate.R`.
- [x] **P13** [29] MODERATE: clamped rows (|r_c|>0.999) NA the r-CI but emit fabricated z outputs
  (z = atanh(0.9999) regardless of input; z_se exploded). Fix: NA the z columns too. File:
  `R/es_disattenuate.R`.
- [x] **P14** [30] MODERATE: `compute_sem(icc_se=)` expects RAW-scale SE but `es_from_icc()`
  default (bonett) returns ln(1-ICC)-scale under the same name → sem_se inflated 5.33×.
  Fix: document scale + cross-reference + heuristic warning. Files: `R/psychometric_utils.R`,
  `R/es_from_ICC.R` (cross-ref).
- [x] **P15** [24] MAJOR-direct/MINOR-pipeline: `es_from_spearman_rho(n_exp=, n_nexp=)` without
  `n_sample` returns uncorrected d/g SEs + NA r/z SEs (missing `n_sample <- n_exp+n_nexp`
  fallback). File: `R/es_from_spearman_COR.R`.
- [x] **P16** [17] MINOR: alpha/ICC direct-call edges emit Inf/NaN instead of NA (n_sample≤2,
  n_items<2, alpha>1, icc outside [-1,1] / icc=1). Files: `R/es_from_ALPHA.R`, `R/es_from_ICC.R`.
- [x] **P17** [18] MINOR: raw-method Wald CIs can exceed 1 for alpha/ICC; B7/B8 check only the
  point estimate. Fix: CI-bound `[UNUSUAL]` check. File: `R/internal_flags.R`.
- [x] **P18** [32] MINOR: `compute_sem` k=1 → SE=Inf silently; icc>1 → sem=NaN with se=0.
  Fix: guards → NA + warning. File: `R/psychometric_utils.R`.
- [x] **P19** [33] MINOR: `compute_sem` Wald CI undercovers at small df (89.7% at n=10,k=2).
  Fix: exact chi-square CI in the same-sample branch (df known). File: `R/psychometric_utils.R`.
- [x] **P20** [34] MINOR: disattenuate warning text wrong in the 0.999<|r_c|<0.9999 band; no
  |r|>1 input guard. File: `R/es_disattenuate.R`.
- [x] **P21** [35] MINOR: `reliability_change_score` silent on negative output (NaNs downstream
  in es_disattenuate); Rd omits the equal-reliabilities assumption. File: `R/psychometric_utils.R`.

## D — Endpoint feature consistency

- [x] **P22** [5] MINOR: deprecated `measure` alias silently overrides an explicit
  `user_es_target_measure_*`; no deprecation warning. File: `R/es_from_USER.R`.
- [x] **P23** [11] MINOR: `es_from_pt_bis_r` accepts but ignores `smd_to_cor`. File:
  `R/es_from_point_biserial_COR.R`.
- [x] **P24** [12] MINOR: per-row `smd_var` NA errors while `smd_denom` NA defaults — align
  (NA → default). File: `R/internal_es_from_d.R`.
- [x] **P25** [40] MINOR: E7 silent when a pool contains ONLY paired-t rows under `pool_sd=TRUE`.
  File: `R/internal_flags.R`.

## E — Tests

- [x] **P26** [16][25][31] MAJOR (infrastructure): ALL psychometric tests (alpha, icc,
  disattenuate, psychometric-utils, spearman, prop) live only in `tests_save/checked/`, which
  `.Rbuildignore` excludes — zero executed coverage; and test-alpha/icc/spearman are circular
  (expectations recompute the source formulas). Fix: promote to `tests/testthat/` with external
  anchors (metafor ABT/ARAW, PR/PLO/PFT, psychmeta) and `skip_if_not_installed` guards.
- [x] **P27** [25] MINOR: `test-EXTERNAL-PROP.R` omits SE-vs-metafor assertions at the p=0/1
  boundaries (exactly where the package deviated). Add after P5.
- [x] **P28** [39] MINOR: `test-pooled-variance-calibration.R` grades a `qnorm` CI the package
  does not emit; single results-derived q=0.64 point; heteroscedastic block missing upper
  coverage bound.
- [x] **P29** [38] MINOR: `test-formula-verification.R:15-17` attributes Morris (2008) metadata
  to "Bonett (2008)".

## F — Documentation

- [x] **P30** [4] MINOR: NEWS "bit-exact with escalc(measure='SMD1')" needs the vtype qualifier.
- [x] **P31** [37] MINOR: NEWS "covers at the nominal rate in every regime tested" must be scoped
  "given a correctly specified r_pre_post".
- [x] **P32** [15][26][36] MINOR: CLAUDE.md stale in three places (ICC SEs wrong-scale/wrong-form;
  Spearman SE missing the (1+r²/2) factor; compute_sem documents the removed icc_se default).
- [x] **P33** [41] MINOR: `data-raw/PRE-POST-SMD/pre-post-strategy.md` documents old formulas and
  wrong defaults.
- [x] **P34** [48] TRIVIAL: docblock "The three pooled forms … " names two.
  File: `R/internal_multiple_formulas.R:852-854`.
- [x] **P35** [13] TRIVIAL: document that `glass_robust` has a single variance form (smd_var
  ignored). File: `R/internal_es_from_d.R` roxygen/@param.
- [x] **P36** [27] TRIVIAL: dangling Miller (1978) inverse-FT citation (no back-transform exists,
  by design). File: `R/es_from_PROP.R` references.
- [x] **P37** RESOLVED (user decision 2026-07-17): keep version **2.0.0** everywhere (DESCRIPTION + NEWS already read 2.0.0). NOTE for release day: the existing git tag `v2.0.0` points at pre-feature content and must be moved (or deleted and re-created) at the actual release commit before any CRAN/GitHub release references it.
  Original issue: version renumbered 2.1.0 → 2.0.0 while git tag `v2.0.0` exists [6] USER DECISION: version renumbered 2.1.0 → 2.0.0 while git tag `v2.0.0` exists
  at different content. Re-tag at the release commit, or restore 2.1.0.

## Explicitly verified NON-problems (no action)

Endpoint core bit-exact vs metafor SMD LS2/LS, SMD1/SMD1H (arm mapping correct, control-df J,
Var(d)=Var(g)/J² both conventions); default endpoint path bit-identical to master except the
deliberate Welch md CI; pre/post family bit-exactness intact at current tree; alpha Bonett/raw
bit-exact vs metafor ABT/ARAW (sign convention documented); ICC(3,1) consistency SE correct
(MC 0.94–0.95); Spearman default path well calibrated (Bonett-Wright factor present);
prop interior points bit-exact vs PR/PLO/PFT incl. exact double-arcsine; disattenuation core
formulas match psychmeta; compute_sem same-sample variance MC-validated (0.9%);
Lord change-score formula correct as the stated special case.

---

## Final status (2026-07-17)

All P1-P36 items FIXED + TESTED (P37 left as a user decision). Full active suite:
**3,501 passed / 0 failed / 0 warnings / 0 skipped** (pre-fix baseline: 2,536).
Every MAJOR/MODERATE fix landed with a test that failed on the pre-fix tree (TDD proof
captured in the fix-agent reports). New/changed test files: test-etasq-adj-scale.R (376
assertions, cross-route oracle), test-spearman-fixes.R, test-prop-fixes.R,
test-disattenuate-fixes.R, test-psychometric-utils-fixes.R, test-alpha-icc-anchored.R
(metafor ABT/ARAW + F-route oracles), test-flag-consistency-fixes.R, md-family additions
to test-endpoint-smd-variance.R, plus the 8 promoted (and de-circularized) psychometric
suites and the repaired pooled-variance calibration test (qt grading, q-grid, upper bound).
