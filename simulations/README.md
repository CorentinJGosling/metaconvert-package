# metaConvert simulations — validating the conversion formulas the package chooses for you

This folder holds the Monte Carlo programme behind `metaConvert`'s **conversion
method arguments**: the ~16 places where `convert_df()` must pick one of several
competing formulas for turning one piece of input data into one effect size, and
where the package currently ships a default backed by a citation rather than by
evidence we generated.

```r
convert_df(df, measure = "r",
           table_2x2_to_cor = "tetrachoric",   # vs lipsey, cooper_delta, cooper_std
           smd_to_cor       = "viechtbauer",   # vs lipsey_cooper
           or_to_cor        = "bonett",        # vs pearson, digby
           or_to_rr         = "metaumbrella_cases",
           pre_post_to_smd  = "bonett", ...)   # vs d_av, d_z, endpoint
```

Each of those defaults silently picks an **estimand** on the user's behalf. The
aim here is to say, for every one of them, how much the choice costs, in which
region of the design space, and on which of the three things a meta-analyst
actually pools: the point estimate, its standard error, and its interval.

---

---

## State of play — read this first if you are picking the project up

Last full run: **all 8 studies at `nrep = 1000`, 26 cores, no errors**, ~40k
aggregated rows in `data/aggregated/*_nrep1000.csv`.

### Companion documents

| file | what it holds |
|---|---|
| `SOURCE-VERIFICATION.md` | every conversion method checked term-by-term against its source paper, with the adversarial refutations. **Read §4 before proposing any package fix** — it lists five "the package is wrong" claims that were investigated and disproved |
| `HANDOFF-package-fixes.md` | the package-fix brief (6 fixes, all applied; FIX 5's `@note` was later shortened and FIX 6's "leave the hierarchy alone" verdict was reversed — see NEWS.md 2.0.1) and the report format |
| `students/` | onboarding + 8 task sheets for non-statistician RAs, plus `SUPERVISOR-NOTES.md` |

### What the nrep = 1000 run established

Scored against each method's **own** estimand, everything is healthy: bias ≤ 0.02
mean, `se_ratio` 0.96–1.05, coverage 0.94–0.97 on the continuous studies. No
non-finite values, no empty cells, no absurd magnitudes anywhere.

Low-coverage rows exist in every study and are, without exception, a method scored
against the *other* estimand. That is the documented behaviour, not error. Study 07's
`se_ratio < 0.5` occurs only where the covariate correlation is over-guessed (0.992
when correct); study 08's 0.65–1.62 tracks `r_guess − r_true` monotonically, centred
on 1.005. Both are the misspecification surfaces working as designed.

Confirmed at full grid, and **now fixed**: `or_to_rr = "metaumbrella_exp"` had
coverage **0.497** — the worst shipped option, worse than the 0.756 seen on the
smaller grid. The cause was a non-identified 2x2 reconstruction (roadmap 1.1/1.8) and
a rounding rule that destroyed continuity-corrected tables (roadmap 1.10). After both
fixes its mean |bias| against the sample log RR is **exactly 0.00000** and its full-grid
coverage is **0.976**; `metaumbrella_cases` is identical to it, because the two now
solve the same table rather than searching for it.
`metafor_conv2x2` still fails to reconstruct in up to 83% of a cell (negative
discriminant; full-grid non-estimability 0.152). `grant`'s non-estimability reaches
0.620 in study 05 (`rr × br_guess ≥ 1`), now counted rather than silently dropped.

### Three things to know before continuing

**1. Study 03's `nrep = 1000` files were stale; they have been regenerated
(roadmap 3.5), together with every other aggregate in `data/aggregated/`.** The
history is kept because the lesson is general. The bug was in the *simulation*, not
the package: the `(z)`
routes were scored against `theta_pop = rho`, an r-scale target, so `atanh(rho) − rho`
was being reported as bias. Proof: the correlation between `phi (z)`'s apparent bias
and the pure scale gap was **0.9929**. After the fix (a scale-matched `theta_own`),
`phi (z)` bias fell 0.0347 → 0.0056 and its coverage rose to 0.9618, identical to
`phi (r)`. `tetrachoric (z)`'s remaining ≈0.23 is the genuine CAT estimand mismatch.
*Lesson worth carrying: whenever a study reports both r and z, check the target scale.*

**2. The `se_ratio` > 1 pattern is diagnosed: it is sparsity, and the continuity
correction is innocent.** In both binary studies `se_ratio` is > 1 for every method
(1.17–1.66 on the full grid), including `transpose`, which merely copies `logor_se`
through. This README used to name the +0.5 continuity correction as the likely cause
and say it had not been chased. It has now been chased, and **the causal direction was
backwards**. Split the grid by how likely a zero cell actually is:

| region | conditions | `se_ratio`, study 04 |
|---|---|---|
| **dense** — P(any zero cell) < 1e-4 | 105 | **0.989 – 1.000**, every method |
| **sparse** — everything else | 255 | 1.298 – 1.410 |

Where zero cells essentially never occur, every route is calibrated. The inflation
appears exactly where the tables get sparse, so the driver is **sparsity itself**
(Woolf's variance bias in the log OR / log RR), and the +0.5 correction — which only
ever applies in the sparse region — *suppresses* part of it. Removing the correction
makes the ratio worse (1.20–1.71 against 1.08–1.32), which is the opposite of what
"the correction inflates the SE" predicts. `transpose` being affected is the tell: a
route that copies the standard error through cannot have a per-method bug.

Reproduce it with `sparsity_split()` in `R/06_sparsity.R`; the split is pinned by
`tests/test-sparsity-diagnostic.R`.

**2b. Two routes must be read with `nonest_rate` beside them, or they look
miscalibrated when they are not.** `metafor_conv2x2` sits at 1.62 full-grid and looks
like the worst performer; it is at **0.992** in the dense region. Its full-grid figure
is a **selection artifact** — it fails to reconstruct a table in up to **82.6%** of a
cell (mean non-estimability 0.152) and is scored on the survivors. `grant` in study 05
is the same case: dense-region 1.10, but non-estimability reaching **1.000** in whole
cells where `rr × br_guess ≥ 1`. Every other route reconstructs essentially always
(non-estimability < 0.01), which is why only these two need the caveat.

**2c. Every number in studies 04 and 05 is measured at 1:1 allocation, and 1:1 is
the favourable end for 2b and for the split above.** Both studies hold `p_exp = 0.5`,
so the two arms are always equal. That is not a neutral choice. P(any zero cell)
depends on the two arms *separately* and is driven by the smaller one, so balance
maximises the dense region — moving participants out of either arm can only shrink it:

| `p_exp` | dense conditions | of 360 |
|---|---|---|
| 0.10 | 25 | 7% |
| 0.25 | 65 | 18% |
| **0.50 — the shipped grid** | **105** | **29%** |
| 0.75 | 85 | 24% |
| 0.90 | 45 | 13% |

The same mechanism drives 2b. `metafor::conv.2x2` reconstructs **every** integer
table and fails on roughly **half** the continuity-corrected ones (enumerating the
whole table space at `n = 50`: 0.00 integer, 0.45–0.68 corrected, at every
allocation), so its non-estimability is the *corrected-table share* times that
per-corrected rate — and the share is smallest at balance too, for the same reason the
dense region is largest there: 0.148 at 1:1 against 0.362 at 1:9 or 9:1.
Its shipped 0.152 is therefore a best case, as is `grant`'s in study 05.

**The no-go map itself is not conditional.** The point estimates in the table
further down are allocation-free: with both margins supplied the odds ratio determines
the 2×2 by a quadratic, and the solve recovers the table exactly (max error `0` over
400 draws per allocation, `p_exp` = 0.10 to 0.90). The original form of this item —
"`metaumbrella_exp`'s 0.497 is conditional on `p_exp = 0.5`" — was about the
non-identified tie the old *search* faced at equal arms; items 1.1/1.8/1.10 replaced
the search with the solve, so that reading is retired rather than corrected.

Reproduce both halves with `dense_region_by_allocation()` and
`corrected_table_share()` in `R/06_sparsity.R`; pinned by
`tests/test-allocation-scope.R`, whose first assertion fails on purpose if anyone adds
a second `p_exp` level — at which point this paragraph should be deleted, not edited.

**3. Pre-post coverage is thinner than the eight-study structure suggests, but the
old form of this note named the wrong kernels, and named them in prose, which is how
it went stale unnoticed.** Study 08 does exercise all five `pre_post_to_smd` values;
what it under-exercises is the code beneath them — and not the code the old table named.
The two-group entry point is a dispatcher rather than a kernel.
Under `pool_sd = FALSE` — its default, and study 08's setting — it calls the
single-group kernel once per arm and combines, so study 08 runs that kernel **twice on
every row**. Measured by rebinding both kernels to counting wrappers and calling study
08's own route, `es_from_means_sd_pre_post`:

| route call | `.single_group_pre_post_to_smd` | `.pooled_pre_post_to_smd` |
|---|---|---|
| `pool_sd = FALSE` (the default), 2 rows | **4** — two per row, one per arm | 0 |
| `pool_sd = TRUE`, 1 row | 0 | **1** |

Dispatch is at `R/internal_multiple_formulas.R:1389` (pooled) and `:1400`/`:1407` (per
arm); the per-row `mapply` that makes the count scale with rows is at
`R/es_from_PAIRED_MEANS.R:154`.

19 exported functions take `pre_post_to_smd`. They used to split across **four** code
paths: the fourth was the paired-t/F family, which called neither kernel and
reimplemented the arithmetic inline. **That path has since been removed** — those five
routes now delegate through `.paired_t_to_smd()`, so there are three kernels and one
implementation of each. Where they land now:

| code path | exported routes | run by the simulation? | otherwise covered by |
|---|---|---|---|
| `.pre_post_to_smd` — the two-group dispatcher plus the `pool_sd = FALSE` combination step (difference of per-arm d/g, summed variances, CI on `n_exp + n_nexp - 2`) | the 7 two-group means / mean-change routes | **yes** — study 08, every row | study 99's deterministic cross-route identity block |
| `.single_group_pre_post_to_smd` — the per-arm engine; all five methods | the same 7 indirectly, plus 7 of the 8 `*_single_group` routes directly | **yes, but only as the per-arm engine.** Its own entry points are unsimulated | `tests_save/checked/test-EXTERNAL-PAIRED-SINGLE-GROUP.R`, `test-PAIRED-SINGLE-GROUP-2.R`, `test-paired-cross-validation.R` |
| `.pooled_pre_post_to_smd` — `pool_sd = TRUE`: own pooled standardizer, own variance forms, own df | the 7 two-group routes that accept `pool_sd` | **no — genuinely never run** | `tests_save/checked/test-pooled-variance-calibration.R`: 6 blocks, **78 assertions**, including bit-exact agreement with `metafor::escalc` and Monte-Carlo calibration of the one branch with no closed-form reference. Plus `test-pool-sd.R` |
| the paired-t/F family — **delegates** to the single-group kernel since the duplication was removed | `es_from_paired_t`, `_t_pval`, `_f`, `_f_pval`, `_t_single_group` | not by Monte Carlo; the four two-group ones are in study 99's equivalence block | `tests/testthat/test-paired-t-delegation.R` (what the delegation must preserve), plus the cross-route pins in `tests_save/checked/test-paired-cross-validation.R` and `test-INTERNAL-FULL-LIFECYCLE.R` |

**"Unsimulated" is the accurate word; "uncovered" is not.** Every path above has package
test coverage, and the one the old table called "never run" has the most of any of them.
What the simulation cannot currently say anything about is *calibration under a
data-generating mechanism* for paths 3 and 4, and for path 2 reached through its own
entry points.

The fourth path **was** a duplicate implementation: it reproduced the single-group
kernel's arithmetic in place rather than calling it. The two agreed, with a maximum
|difference| over d, SE, g and g's SE of `0` for `morris_dz` and `1.1e-16` for
`morris_drm`, but only for as long as both were kept in step by hand. The duplication
has been removed, and with it the risk of drift; the delegation moves no
number a user can see (max |difference| `3.55e-15` over a 576-row grid). Those routes
still offer only two of the five methods and still **error** on the other three rather
than silently coercing, because a paired t does not identify the separate pre/post SDs
the other three need.

The removal was **detected by the test above**, not planned alongside it: that test
asserted "calls neither kernel", the refactor made it false, and it failed. That is what
item 4.4 bought — the claim had been made executable, so restructuring the dispatch
could not silently invalidate the paragraph describing it.

A `pool_sd` study reusing study 08's DGP and grid with a new `estimate()` would take
pre/post from **2 kernels to 3**. It cannot be called "study 08b": study 08 already emits
`08a_pre_post_to_smd_d` and `08b_pre_post_to_smd_g`, both shipped in `data/aggregated/`.

One correction to an earlier claim. `pool_sd = TRUE` was described here as "the route carrying the
known docblock defect". That is stale: the misattribution of d_av's leading variance term
to Bonett (2008) eq. 19 was corrected in `4955631` (2026-08-13), *before* this paragraph
was first tracked. `R/internal_multiple_formulas.R:1465-1479` now states exactly which
part of eq. 19 the code does use (the fourth-moment g^2 coefficient, reproduced to 1e-14)
and quantifies the departure of the part it does not: ~2% at n = 50/50, ~3% at 30/30,
~10% at 12/11, ~32% at 100/4, and ~38% at 100/10 under heteroscedasticity. That is a
documented departure with a stated calibration range, not a defect.

Pinned by `tests/test-prepost-kernel-coverage.R`, which installs the counters above and
asserts the dispatch in both directions, the 7/8/5 route split, and which studies touch
which path. Encoding the earlier claim in it instead — `single_group = 0`, and no fourth
path — turns the file red. Prose is what rotted here; this is the same statement made
executable.

### Compute budget

`run_04` is the bottleneck at **41 min** for `nrep = 1000` on 26 cores (estimraw
reconstruction); `run_05` 26 min; everything else ≤ 8 min. A production run at
`nrep = 10000` is roughly **1.5–2 h** for all eight. Do not raise `nrep` to buy
precision — MCSE is already far below the contrasts of interest (see the note further
down); spend compute on missing design factors instead.

### Package fixes that landed during this work

All verified, all regression-clean (1304-assertion subset per-file identical; the
full archived suite is 7189 with `NOT_CRAN=true` — **the baseline command must set
it**, or 12 files / 1585 assertions are silently skipped):

- `.mapply_memo()` — memoises the tetrachoric solve. 10.8× on `es_from_2x2`, bit-exact
- `.cor_to_smd_vec()` — vectorises the `cor_to_smd` routes. 29×, estimates bit-exact
- two crash guards in the Di Pietrantonj reconstruction (it aborted a whole
  `convert_df()` run on common-outcome data)
- `reverse_2x2`/`reverse_phi`/`reverse_chisq` now **reflect** the tetrachoric CI
  instead of negating it in place
- `or_to_cor = "bonett"` derives `small_margin_prop` from the margins instead of
  silently falling back to `lipsey_cooper` (it was the `convert_df()` default and
  usually did not run)
- the ANCOVA `@note` and `CLAUDE.md` now state that the point estimate is **not**
  unbiased on the 8 statistic-inverting routes
- new `[INFO]` flag B1b when a correlation CI bound escapes [-1, 1]

---

## Why this is not Poom (2022) or Jacobs & Viechtbauer (2017)

Both are in `papers/` and both must be cited prominently — they are the closest
prior art and they overlap parts of this design.

| | Prior work | This programme |
|---|---|---|
| **Poom & af Wåhlberg 2022** (*RSM* 13:508-519) | Conversion accuracy of formulas *in the abstract*, same-sample benchmark, point estimates | — |
| **Jacobs & Viechtbauer 2017** (*RSM* 8:161-180) | Biserial *r* and its sampling variance: bias, variance, coverage | — |
| **Object of study** | a formula as written in a paper | a formula **as implemented in a package**, reachable through a documented argument |
| **Nuisance parameters** | assumed known | **crossed: true value × the value the analyst guesses** |
| **Coverage target** | population parameter | both population parameter **and** same-sample statistic, reported separately |
| **Output** | a table | a table + the package defaults it does or does not support |

The framing is therefore **software validation and decision support**, not formula
discovery. That makes the prior art an asset: reproducing Jacobs & Viechtbauer's
biserial results is positive evidence that `smd_to_cor = "viechtbauer"` is
implemented correctly, which is exactly what a package paper needs to show.

The genuinely unoccupied ground, in decreasing order of strength:

1. **Nuisance-parameter misspecification surfaces.** No published simulation
   crosses the *true* baseline risk (or pre-post correlation, or covariate-outcome
   correlation) against the value the analyst *guesses*. This is the situation
   every meta-analyst converting an OR to an RR is actually in.
2. **Reconstruction algorithms that exist only in software.** `metaumbrella`'s
   grid searches for recovering a 2×2 table from an OR have no closed form and no
   published accuracy characterisation anywhere.
3. **Calibrating metaDETECT's Discordant checks.** Category E asserts thresholds
   (`dispersion_max`, `direction_disagreement_min`, `outlier_min_deviation`, …)
   that cannot currently distinguish "flag fired because the extraction was wrong"
   from "flag fired because two correct formulas legitimately disagree that much".
   Running every applicable route on *error-free* simulated data gives the null
   distribution of cross-method disagreement — i.e. the false-positive rate of
   every one of those checks. This single analysis serves all three projects.

---

## The two design decisions that define this rebuild

### 1. Two targets, always reported separately

Every replication records **both**:

| column | meaning |
|---|---|
| `theta_pop` | the population parameter the data were generated from |
| `theta_sample` | the same statistic recomputed on that replication's own sample — what the trialist would have reported had they analysed the other metric on the same participants |

The previous pipeline used `theta_sample` only. That is a defensible choice for
**bias** — it is the meta-analyst's question, and it is the same device the
[PETRA IPD paper](../papers/pre_post_ipd/) uses deliberately. But it makes
**coverage and variance calibration meaningless**, because a nominal 95% interval
is not built to cover a random quantity. Coverage of `theta_sample` equals 0.95
only when the conversion discrepancy happens to equal the sampling SE, so ranking
on `|coverage - 0.95|` **rewards miscalibration**. In the old shipped files, 6,706
of 9,330 rows have that statistic at exactly 1.0000, and the Shiny ranking put
Cooper first precisely because its SE was ~2.7× too wide.

So: bias against `theta_sample` (argued, not assumed), coverage and SE
calibration against `theta_pop`. Where the two disagree, **the gap is estimand
mismatch and must be named as such, not called bias**. Lipsey–Cooper's apparent
bias in the old SMD→r file is entirely this: it targets the point-biserial while
the benchmark was the latent correlation, which is why its "bias" is flat in *n*
(−0.189, −0.189, −0.183, −0.185, −0.182 from n=25 to 300). Estimand gaps do not
shrink with sample size; that flatness is the diagnostic signature.

#### Consequence: on a `*_sample` target, `RMSE² ≠ bias² + EmpSE²`

This trips people up when reading `emp_se` and `rmse` side by side, so it is worth
stating explicitly. The familiar decomposition

```
RMSE² = bias² + EmpSE²
```

holds **only when the target is a fixed constant**. `EmpSE` is defined by ADEMP as
the dispersion of the *estimator*, `sd(est)`, which is what `performance()`
computes — correctly. But `RMSE` is `sqrt(mean((est − target)²))`, and against a
random target that expands to

```
bias² + var(est) + var(target) − 2·cov(est, target)
```

`theta_sample` is recomputed on the *same* simulated participants as the estimate,
so the two are strongly positively correlated and the covariance term is large.
The deviation is therefore much less variable than the estimate itself, and RMSE
can sit well **below** EmpSE without anything being wrong.

Measured on study 01a, `n = 25`, `rho = 0`, `p_exp = 0.5` (nrep = 1000):

| method | cor(est, target) | sd(est) | sd(est − target) | RMSE |
|---|---|---|---|---|
| `viechtbauer` | 0.802 | 0.2547 | 0.1522 | 0.1522 |
| `lipsey_cooper` | 0.801 | 0.1932 | 0.1261 | 0.1261 |

RMSE equals `sd(est − target)`, not `sd(est)`. Against the fixed `theta_pop` in
the same cell the identity is restored exactly: RMSE 0.2550 vs
`sqrt(bias² + sd(est)²)` 0.2551 for `viechtbauer`, and 0.1934 vs 0.1935 for
`lipsey_cooper`.

**Practical rule.** Compare `rmse` only with other `rmse` values on the same
target. Do not reconstruct it from `bias` and `emp_se` unless the target is
`theta_pop`, and do not read `rmse < emp_se` as a defect.

#### Worked example: why two methods can look equally good and not be

Scored against `biserial_sample`, `viechtbauer` and `lipsey_cooper` cross at
around `rho ≈ 0.35` at n = 25 and `rho ≈ 0.15` at n = 300, and near the crossing
their RMSEs are nearly equal. That similarity is an artefact, produced by two
opposing forces that happen to cancel:

1. **`lipsey_cooper` carries a growing estimand bias.** It targets the
   point-biserial, which at `p_exp = 0.5` is `0.7979 × biserial`, so against the
   biserial its bias is forced to be `−0.2021 × rho`. Observed at n = 300:
   −0.0196, −0.0509, −0.1041, −0.1523 at rho = 0.10, 0.25, 0.50, 0.75, against
   −0.0202, −0.0505, −0.1011, −0.1516 predicted. None of it is computational error.
2. **`viechtbauer` has ~1.25× the dispersion, purely because its estimand is
   ~1.25× larger.** The EmpSE ratio at n = 300 is 1.258, 1.258, 1.259, 1.258,
   1.255 across rho, against `1/0.7979 = 1.253`. Relative precision is identical.

The bias term is constant in *n* while the noise term shrinks, so the crossing
moves left as *n* grows and the estimand gap always wins eventually. Scored on
`own`, both are unbiased (|bias| ≤ 0.004 at every rho, n = 300) and their RMSEs
differ by exactly that 1.25 scale factor. **Equal RMSE against a foreign estimand
is not evidence of equal performance** — which is what the app's *Estimand check*
panel is for, and why it defaults to `own`.

### 2. Call the package, never re-type the formula

The old scripts contained zero calls to `metaConvert`. Every formula was retyped
inline, and the copy had already drifted:

| location | simulation computed | package ships |
|---|---|---|
| `archive/legacy_scripts/1_SIM_SMD_to_COR.R:56` | `vd / (vd + a)` | `vd / (d² + a)` — and `R/internal_multiple_formulas.R` documents that the first form is `esc`'s and is *intentionally rejected* as Fisher-inconsistent |
| `archive/legacy_scripts/8_COR_to_SMD.R:38` | `4 * sqrt(r_se)` | `4 * r_se²` |

A study whose thesis is "here is evidence for our package's defaults" cannot
compute something other than what the package does. Each study calls `es_from_*()`
directly, and the harness loads the package from **source**
(`pkgload::load_all`), never an installed build.

---

## The speed problem, and why it is solved

Calling `metaConvert` used to be slower than the inline reimplementation. It was,
but for two fixable reasons — neither of them intrinsic to the package.

**(a) The functions are vectorised; the old code called them one row at a time.**
Measured on `es_from_means_sd`, 10,000 replications:

| | elapsed |
|---|---|
| per-replication loop | 9.00 s |
| **one vectorised call** | **0.15 s** |
| the DGP itself (mvrnorm + summaries) | 1.09 s |

So the ES layer costs ~14% of data generation. The harness enforces this split:
`generate()` loops over replications, `estimate()` never does.

**(b) `es_from_2x2()` spent 99.8% of its time in the tetrachoric ML solve**
(`escalc` → `.rtet` → `optim` + `mvtnorm::pmvnorm`), ~20 ms per row, mapped row by
row. Count data repeats heavily — 251 distinct tables in 1,500 draws at n=50/arm —
so identical argument tuples are now solved once and expanded back, via
`.mapply_memo()` in `R/internal_multiple_formulas.R`.

| n per arm | distinct tables / 10k | speedup |
|---|---|---|
| 25 | 202 | ~50× |
| 50 | 363 | ~28× |
| 100 | 645 | ~17× |
| 300 | 1,516 | ~7× |

Measured end to end: **4,000 rows at n=50/arm went from ~88 s to 8.1 s, bit-exact.**
This is a package-level fix, so every `convert_df()` user gets it, not just this
folder. Covered by `tests/testthat/test-mapply-memo.R`; the five archived suites
touching these paths (`test-2x2`, `test-ES-OR`, `test-ES-COR`, `test-PHI-CHISQ`,
`test-ES-RR`) pass 587/587 unchanged.

---

## Layout

```
simulations/
├── README.md              this file
├── run_all.R              entry point: source it, then run_03(), run_everything()
├── R/                     the harness
│   ├── 00_config.R          paths (no absolutes), seed policy, defaults
│   ├── 03_performance.R     ADEMP measures + MCSE, dual targets
│   └── 04_runner.R          run_study(): grid × reps, seeded, parallel, saves raw
├── studies/               one file per conversion area (see table below)
├── data/
│   ├── raw/                 per-replication output (large; regenerable, seeded)
│   └── aggregated/          summary tables — the shipped result
├── app/                   Shiny explorer (app.R, rsconnect deployment record)
├── papers/                Jacobs & Viechtbauer 2017, Poom & af Wåhlberg 2022
└── archive/               everything superseded — kept, not deleted
    ├── legacy_scripts/      the 11 original MonteCarlo scripts + save/ + cancel/
    └── app_versions/        app2.R
```

`data/aggregated/` currently still holds the nine **legacy** `*_AGG*.txt` files.
They are kept only so the app keeps running and so the defects below can be
checked against them; they are superseded and must not be quoted.

### Seeding and reproducibility

`SIM_BASE_SEED = 20260729`. Each `(study, condition)` derives a deterministic
L'Ecuyer-CMRG stream from it, so any single condition re-runs bit-for-bit in
isolation, sequentially or in parallel. The old scripts had **no `set.seed`
anywhere**, wrote raw output to `D:/simulations/data/` which is not in the repo,
and read filenames they never wrote — none of the nine shipped aggregate files can
be regenerated by anyone, including their author.

---

## The eight conversion areas

| # | Area | Package argument(s) | Status |
|---|---|---|---|
All are implemented and run: `source("run_all.R")` registers `run_01` … `run_11`, plus
`run_99` (a deterministic wiring check, not a Monte Carlo study). Studies 01–09 cover the
conversion routes; studies 10-11 cover the reliability family, which had no coverage at all
until roadmap item 3.1.

| # | Area | Package argument(s) | File |
|---|---|---|---|
| 01 | SMD → *r* / *z* | `smd_to_cor` | `studies/01_smd_to_cor.R` |
| 02 | *r* → SMD | `cor_to_smd` | `studies/02_cor_to_smd.R` |
| 03 | 2×2 → *r* / *z* | `table_2x2_to_cor` | `studies/03_2x2_to_cor.R` |
| 04 | OR → RR | `or_to_rr` | `studies/04_or_to_rr.R` |
| 05 | RR → OR | `rr_to_or` | `studies/05_rr_to_or.R` |
| 06 | OR standard-error imputation | (`es_from_or`) | `studies/06_or_se_imputation.R` |
| 07 | ANCOVA means → SMD | `cov_outcome_r` guess | `studies/07_ancova_to_smd.R` |
| 08 | pre/post → SMD | `pre_post_to_smd` | `studies/08_pre_post_to_smd.R` |
| 09 | OR → *r* / *z* | `or_to_cor` | `studies/09_or_to_cor.R` |
| 10 | reliability (ICC, alpha) | `icc_to_es`, `icc_agreement_se` | `studies/10_reliability.R` |
| 11 | reliability **standard errors** (alpha vs omega) + the floor-effect detector | `alpha_to_es`; the omega SE that is *not* shipped | `studies/11_reliability_se.R` |

### Study 10: the reliability family had zero coverage, and it mattered

Nine studies covered the conversion routes. `grep -i "alpha|icc|omega"` over `studies/`
and `R/` returned **nothing**. That is worse than a missing row in a table, because the
reliability audit produced several numbers that lived only in commit messages and code
comments — each one a simulation run once and thrown away. A number recorded that way
cannot be re-checked when the code moves underneath it, which is exactly what happened
to the *"coverage around 0.74–0.76"* claim that shipped in `?es_from_icc` and was wrong
at every *n*.

Study 10 makes three of those numbers re-runnable at a fixed seed.

**10a — the ICC(2,1) coverage table.** A two-way random-effects DGP, which is what an
absolute-agreement ICC assumes and what the package's shared SE formula does *not*: the
formula is a one-way approximation that drops the between-rater term. At `nrep = 1000`,
`k = 2`, ICC = 0.80, rater variance 50% of the non-subject budget:

| n | empirical SD | 95% coverage | SE ratio |
|---:|---:|---:|---:|
| 20 | 0.573 | 0.825 | 0.721 |
| 50 | 0.459 | 0.654 | 0.563 |
| 200 | 0.421 | 0.295 | 0.305 |
| 1000 | 0.393 | 0.142 | 0.146 |

Coverage **degrades as n grows**, because ICC(2,1) inherits the between-rater mean
square (k−1 df) while the reported SE shrinks like 1/√n. The grid also contains
`rater_share = 0`, and that control is the point of the design: there the formula's own
assumption holds and coverage returns to nominal (0.939 / 0.953 / 0.944 / 0.959), so the
degradation is attributable to rater variance rather than to a broken DGP. The
comparison is internal and does not rest on trusting this study's data generation.

The same run also settles a design question by demonstration rather than assertion:
`agreement` and `consistency` differ by **exactly 0** in bias and coverage across all 32
cells, which is why flipping the default `icc_type` was rejected as a remedy — it
relabels rows and changes no number.

**`v36_false_positive_rate()`** and **`icc_se_by_k()`** are standalone deterministic
functions rather than `run_study()` calls, in the style of `99_route_equivalence.R`:
they measure a *flag's* behaviour, not an estimator's, so they have no target column and
do not fit the generate/estimate contract. The first refutes tightening V36's decimal
gate (a 100-study pool still flags 95.8% of the time at 4 decimals — the driver is pool
size, not precision); the second quantifies V42's premise.

**What is left out, and why.** Roadmap 3.1 also asks for the low-ICC
Bonett-vs-Fisher-TF pooling comparison (item 1.5). `icc_to_es = "fisher_tf"` **is not
implemented**, so simulating it would mean re-typing the candidate formula inline and
scoring the package against a private copy of a method it does not ship — the exact
practice this rebuild exists to stop. When 1.5 ships, add a `10c` that calls the real
route.

### Study 11: the standard error metaConvert refuses to compute

`es_from_cronbach_alpha()` computes an SE from a closed form in `(n, k)` alone — Bonett's
`2k/((k−1)(n−2))`. `es_from_omega()` **refuses to**, and returns `se = NA` unless the
primary study reported an SE or a CI. The rationale in `?es_from_omega` is that Bonett's
variance descends from the Feldt (1965) / Kristof (1963) *F* result, which assumes
essential tau-equivalence — the assumption omega exists to drop.

The applied literature ignores the distinction. Villacura-Herrera et al. (2025,
*Work & Stress* 39(2) 169–196) ran 13 omega reliability generalisations through
`escalc(measure = "ABT", ai = omega, mi = n_items, ni = n)` — Bonett's *alpha* variance,
applied to omega. No primary study reported an omega SE. So metaConvert cannot reproduce
any published omega RG meta-analysis, and the question is whether the refusal is earned.

**No package change was needed to find out.** Feeding `omega_hat` to
`es_from_cronbach_alpha()` is bit-identical to metafor's `ABT` up to metaConvert's
deliberate sign flip, so the candidate is scored as *shipped code* — honouring the same
rule that keeps study 10's `fisher_tf` comparison out of the tree.

**11a — the head-to-head.** Both coefficients from **one covariance matrix per
replication** (the fairness rule; see below), 42 cells, `nrep = 1000`. The primary measure
is *centred* coverage — coverage about the estimator's own mean across replications, which
isolates SE calibration and needs no plim, so it is immune to the four-way ordinal estimand
confound. Δ = omega − alpha, paired over replication index:

| k | max &#124;Δ se_ratio&#124; | max &#124;Δ coverage&#124; | mean paired *r* |
|---:|---:|---:|---:|
| 3 | 0.105 | 0.083 | 0.831 |
| 8 | 0.026 | 0.021 | 0.990 |
| 20 | 0.030 | 0.024 | 0.997 |

At k ≥ 8 the two coefficients **fail and succeed together**, against per-coefficient
`se_ratio` levels spanning **0.468 to 1.095**. The asymmetric treatment in the package is
not supported at the item counts real instruments use. At k = 3 the differential is real
(5.5 MCSE) and runs in *both* directions across cells; the pairing benefit also collapses
there (*r* = 0.83 vs 0.99), so k = 3 is reported as a caveat, not as a gate.

**Response format is the axis, not the coefficient.** Exact marginal moments, so these are
arithmetic rather than measurement:

All rows below are k = 8, n = 200, from `11a_delta_paired_reps1000.csv`:

| format | excess kurtosis | `se_ratio` (α / ω) | centred coverage (α / ω) |
|---|---:|---:|---:|
| continuous normal *(bound)* | 0.00 | 1.064 / 1.053 | .962 / .958 |
| 5-pt symmetric | −0.50 | 1.045 / 1.042 | .961 / .958 |
| 7-pt floor-skew | −0.42 | 0.980 / 0.986 | .943 / .944 |
| **5-pt floor-skew (BAT-like)** | **−0.32** | **0.948 / 0.941** | **.941 / .936** |
| 5-pt severe floor (65% bottom) | +2.10 | 0.783 / 0.768 | .889 / .882 |
| binary, 10% endorsement | +5.11 | 0.632 / 0.614 | .784 / .777 |
| continuous χ²(1.5) *(bound)* | +8.00 | 0.520 / 0.494 | .686 / .665 |

Five-point categorisation **caps** the reachable excess kurtosis near 2, so the continuous
cell — where the formula is worst — has no empirical referent in RG data and is labelled a
bound, not a condition. Response format is also the only candidate gate a data extractor
can actually read off a primary study; kurtosis is not reported by anyone.

**The benchmark gate.** `bin_10` is Maydeu-Olivares, Coffman & Hartmann (2007) Table 5,
where normal-theory alpha coverage is ≈ 0.79 and **flat in n**. Measured **0.784 / 0.783**
at n = 200 / 1000. `run_11()` asserts this before any omega number is quoted and warns
loudly if it fails: the decision turns on a differential of 0.03 against levels spanning
0.47 to 1.10, so a generator bug would be invisible in the differential and decisive in
the level.

**11b — the shipped rationale, tested directly.** If tau-equivalence were the operative
mechanism, the omega arm would degrade as the loadings spread. It does not:

| loadings | `se_ratio` ω | coverage ω |
|---|---:|---:|
| tau-equivalent (all .70) | 0.984 | .952 |
| congeneric (.40–.90) | 1.039 | .961 |
| extreme (.25–.95) | 1.008 | .959 |

The tau-equivalent cell — where the rationale says the formula *is* valid — is marginally
the least well calibrated of the three. **The argument shipped in `?es_from_omega` predicts
nothing that happens.**

**11c — propagation into a pooled RG, at zero new model fits.** It resamples 11a's raw
frame, so the sampling distributions it pools are the measured, heavy-tailed ones. Writing
`c` for the factor by which the true variance exceeds the reported one, four results are
algebraic and are used here as *validity gates*, not findings: μ̂ is unbiased for any
weights (c never enters `w = 1/(τ̂² + v)`, a function of `(n, k)` alone); the REML pooled SE
is self-correcting; Hartung–Knapp is exactly invariant to a uniform rescaling of the
weights; and the PI is too *wide*. All four hold. What is contingent:

- **Pooled inference survives.** Coverage with the Bonett variance vs coverage with
  **oracle** variances differs by at most **0.0095** across the unconfounded cells — inside
  Monte Carlo error. Only fixed-effect pooling takes the damage (coverage 0.082–0.935).
- **The heterogeneity statistics do not.** I² inflation tracks the predicted `1 − 1/c`
  almost exactly: for binary items, predicted 62.3%, **measured 55.5% (K=10) and 62.7%
  (K=30) against a true I² of exactly zero.** A homogeneous pool reports substantial
  heterogeneity out of nothing. For the BAT-like format, predicted 11.5%, measured 12.9 /
  10.1.
- **The PI over-covers** (≥ .95 in 33 of 36 cells), which partly *masks* the well-known
  small-K PI shortfall — the opposite of the intuitive prediction.
- **The one real bias route is Spearman–Brown weighting.** `v ∝ k/(k−1)` decreases in k, so
  longer instruments get more weight *and* have higher alpha — a confound the formula
  creates itself. Measured |bias| up to **0.032 log units**. Confounding the SE *error*
  with the effect cannot bias μ̂ at all, because the error never enters the weights.

So the headline is neither "the formula is fine" nor "the formula is broken": **the pooled
coefficient survives, the heterogeneity statistics do not** — which exonerates the headline
numbers of published RG work while indicting the heterogeneity narratives built on them.

**The fairness rule, and why it has its own test.** A pilot comparing a Pearson-covariance
alpha against a WLSMV/polychoric omega measured Δ = −0.107 and would have concluded the
formula is materially worse for omega. On the *same* covariance matrix the same contrast is
−0.010. The estimator swap manufactures the entire effect the study exists to detect, so
`.rel_coefficients()` computes both from one `S` and
`tests/test-study-11-reliability-se.R` asserts it.

**Why `se_ratio` is secondary and coverage is primary.** The sample SD of a transformed
reliability is a fourth-moment quantity, so its Monte Carlo error is
`sqrt((κ+2)/(4R))`, not the `sqrt(1/(2(R−1)))` every normal-theory formula quotes.
`rel_stability()` records κ per cell and blanks `se_ratio` where it exceeds 2, and reports
the robust IQR/1.349 spread beside the SD — where they disagree, the SD is the broken one.

**11d — deterministic audit of the claims in `?es_from_omega`.** Seconds, no Monte Carlo,
and it settles shipped documentation outright. `sqrt(2k/((k−1)(n−2)))` at k = 3, n = 200 is
**0.12309** — finite, and the **largest** value in the k-sweep, decreasing monotonically to
`sqrt(2/(n−2))` = 0.10050. **The shipped sentence "it diverges at k = 3" is false as
written.** The block also reproduces the bifactor 14–24% figure as pure arithmetic and
flags its unstated denominator mismatch (`n` vs `n−2`), and pins that Bonett's variance is
coefficient-*free* (identical `(n, k)` ⇒ identical weight whether α = .70 or .98) where
Hakstian–Whalen's spreads the weight by more than 3×.

**11e — can a *reported* number identify the regime?** 11a leaves the practical question
open: response format drives the variance error, but the **number of categories does not
identify it** — the three five-point formats span c = 0.92 to 1.69, because what matters is
where the item mass sits, not how many bins it is cut into. Item skew is never reported, so
that looked like the end of the road.

It is not. A floor effect is visible in the **total-score mean**, which RG extractors
routinely have. On a *C*-point scale of *k* items, rescale it to

```
floor_position = (scale_mean / k − 1) / (C − 1)      ∈ [0, 1]
```

— 0 at a hard floor, 0.5 when symmetric, and comparable across scales of different length
*and* different *C*. 11e sweeps 20 distributions (C ∈ {2,3,5,7} × five degrees of floor) and
asks which reported quantity actually predicts c:

| predictor | Spearman ρ with c | R² (log c) |
|---|---:|---:|
| number of response categories | −0.15 | 0.05 |
| **floor_position** | **−0.96** | **0.77** |
| + total-score SD (Bhatia–Davis normalised) | — | 0.77 (**+0.002**) |
| + category count on top of floor_position | — | 0.79 (+0.02) |

**The mean carries the signal; the SD adds nothing measurable.** That is a useful negative
result — the extraction sheet needs *one* new number per study, not two. And the rule holds
across every category count, which is what makes it shippable:

| floor_position | c | centred coverage | category counts in band |
|---|---:|---:|---|
| < 0.20 | 1.56 – 4.28 | .653 – .874 | 2, 3, 5, 7 |
| 0.20 – 0.28 | 1.08 – 1.44 | .895 – .941 | 2, 3, 5, 7 |
| 0.28 – 0.40 | 0.92 – 1.03 | .939 – .954 | 2, 3, 5, 7 |
| > 0.40 | 0.84 – 1.01 | .952 – .966 | 2, 3, 5, 7 |

Each band contains all four scale types. The only residual *C* effect is at the extreme
floor, where a dichotomous item is worse than a 7-point one at the same position
(c = 4.28 vs 2.21 at floor_position ≈ 0.05–0.10) — which is why the category count is worth
collecting as well, but as an *ingredient* of floor_position rather than as the gate.

**One fitter caveat worth knowing.** The production fitter is `factanal(covmat = )`, pinned
against lavaan to 1e-5 and 30× faster — which is what makes 11c free. But `factanal` floors
uniquenesses at `lower = 0.005`, so it cannot produce a Heywood case; an unregularised
lavaan fit at k = 3 does, with a much heavier tail. The k = 3 numbers here are therefore
those of a *regularised* estimator, and that is a property of the fitter, not of the data.

### Study 09: why `or_to_cor` needed a study of its own

`or_to_cor` was the last conversion argument with no simulation evidence, and it
is the one whose behaviour changed most recently: before 2.0.1 the documented
default `"bonett"` never actually ran — it requires `small_margin_prop`, that
column had no auto-derivation, and a blank value fell through to the
`lipsey_cooper` result. The fix moves |*r*| by +8% to +58% (see the package
`NEWS.md`), a larger shift than several arguments that already had studies.

Its four options do not share an estimand, which is the whole design problem:

| option | estimand |
|---|---|
| `pearson`, `digby`, `bonett` | the **tetrachoric** correlation — the correlation of the two latent continuous normals that were dichotomised |
| `lipsey_cooper` | the **point-biserial** correlation between the binary grouping variable and a latent continuous *outcome*, reached through the Cox logit transform |

Neither equals **phi**, the correlation of the two observed binary variables. At
ρ = 0.5, `p_exp` = 0.3, `p_case` = 0.1 the three quantities are 0.500, 0.379 and
0.257. Scoring all four options against one reference therefore reports estimand
mismatch as bias — the error study 01 documents for `smd_to_cor` — so the study
records a per-method scale- and estimand-matched target (`own`) alongside
`population` and `sample`. The gap between `own` and `population` **is** the
estimand mismatch and must be reported as such.

Two data-generating mechanisms, mirroring study 03: `CONT` (bivariate normal,
both variables dichotomised — the tetrachoric methods are on home ground) and
`CAT` (a genuine 4-cell multinomial — nothing was dichotomised and phi is the
estimand). `es_from_2x2(table_2x2_to_cor = "tetrachoric")` is carried alongside
as a reference route, so the practical question "if I have the full 2×2, should I
convert the OR at all?" is answerable from the same output.

#### Study 09 result: the 2.0.1 default change is justified, and it is a margin effect

`CONT` mechanism, nrep = 1000, 120 conditions, scored on each method's **own**
estimand (so the numbers below are computational error, not estimand mismatch):

<!-- pinned: study09a-result -->
| method | mean \|bias\| | worst \|bias\| | coverage | SE ratio | not identified |
|---|---|---|---|---|---|
| **`bonett`** | **0.0195** | 0.2848 | 0.964 | 1.061 | 0/120 |
| `2x2_tetrachoric` (full table) | 0.0197 | **0.1900** | **0.951** | **0.990** | **67/120** |
| `digby` | 0.0201 | 0.2103 | 0.965 | 1.054 | 1/120 |
| `pearson` | 0.0223 | 0.1865 | 0.965 | 1.054 | 1/120 |
| `lipsey_cooper` | 0.0231 | 0.1937 | 0.945 | 1.051 | 1/120 |

Coverage and the SE ratio are taken over the conditions in which the method's estimate
is **identified**; the last column counts the conditions in which it is not. Only the
full-2×2 route has any, and it has them by design — see the note below.

`bonett` is the best of the four OR-only options, which supports its promotion to
the default. Its 0.0195 against the full-2×2 route's 0.0197 is **not** an ordering:
comparing the two condition by condition over the 120 cells gives a mean paired
difference of −0.000186 with a standard error of 0.00263 (t = −0.07), against a
per-cell Monte Carlo error of about 0.006. The two are tied on this metric, and the
column is quoted to a precision the data does not support. Where the full-2×2 route
separates is elsewhere: it has the **smallest worst-case bias** of the five (0.1900),
and it is the only one whose interval is honestly calibrated where it answers at all
(SE ratio 0.990 against 1.05–1.06; coverage 0.951 against a nominal 0.95). So the
advice is unchanged — **convert the table, not the odds ratio, whenever the table is
reported** — with one caveat the older versions of this table hid, stated below. The
margin over `digby` is thin — 0.0195 against 0.0201 — so the case for the default
rests on the margin-dependence result below, not on this average.

> **Two of these cells moved, for two different reasons, and both are worth stating.**
>
> The `bonett` row changed with the `or_to_cor` fallback (roadmap 2.7). Its
> replications are not all bonett: where the sampled 2×2 is degenerate, bonett cannot
> run and the row takes a stand-in — which used to be `lipsey_cooper` and is now
> `digby`. So this row has always been a mixture, and it is now a closer one: mean
> \|bias\| 0.0197 → 0.0195 here, 0.0183 → 0.0181 in 09b. Every other method is
> bit-identical across that change (delta exactly 0 on bias and coverage), which is the
> control showing the change stayed inside the fallback.
>
> The `2x2_tetrachoric` coverage was simply out of date: this table published 0.928 while
> the shipped aggregate said 0.961, and had done since the item 3.5 regeneration. The
> `bonett` figures were stale too (0.0185 / 0.2571 / 1.051 against 0.0197 / 0.2868 /
> 1.053). Nothing detected it, because nothing tied the prose to the CSV — which is
> what `data/aggregated/PROVENANCE.csv` now does for the code, and does not yet do for
> numbers quoted in the text.

> **The `2x2_tetrachoric` row is not comparable with the one this table carried before
> the tetrachoric fix, and the reason is worth reading before using the route.**
>
> The tetrachoric solve used to be handed the **+0.5-corrected** table. That correction
> is a ratio-measure device — it exists because a zero cell makes an odds ratio
> undefined — and `metafor` pointedly refuses to apply it to `measure = "RTET"`, not
> even under `to = "all"`. On a table with a zero cell the tetrachoric maximum
> likelihood lands on the **boundary**, r = ±1, with an enormous variance, and that
> pair *is* the estimator reporting that the correlation is not identified by these
> counts. Correcting the table first replaced that refusal with an interior estimate
> carrying an ordinary-looking standard error that took no account of the shrinkage:
> on a/b/c/d = 0/20/10/10 it gave r = −0.856 with SE = 0.113 against `metafor`'s
> r = −1 with SE = 149.03 — an inverse-variance weight 1.7 million times too large, on
> a study whose correlation the data cannot pin down. The route is now bit-exact with
> `metafor` on every table, zero cell or not.
>
> So the route did not get worse; it stopped hiding something, and this table has to
> report the something. In **67 of the 120 conditions** the boundary is reached often
> enough that an averaged model standard error is meaningless: the mean over all 120
> is 5344.150, and the median across conditions is 861, so no trimmed or median
> summary recovers it either. Those conditions are exactly the ones in which a
> reconstructed 2×2 hits a zero cell — **30/30 at n = 25**, 22/30 at n = 50, 3/30 at
> n = 300; **33/40 at a 10% case rate**. Restricted to the conditions where the
> estimate is identified, the route is the best-calibrated of the five. In the others
> it declines to answer, which is the correct behaviour and means an OR-only method is
> what you need there.
>
> The threshold that separates the two regimes (`se_ratio <= 2`, in
> `R/08_published_numbers.R`) is not doing the comparison's work: it removes 67
> conditions from `2x2_tetrachoric` and 0, 1, 1 and 1 from `bonett`, `digby`,
> `pearson` and `lipsey_cooper`, moving those four means by at most 0.008.

The averages understate the case, because Bonett's correction is a *margin*
correction and averaging over balanced and unbalanced conditions dilutes it. At
ρ = 0.5, n = 300 (MCSE ≈ 0.003 throughout):

| `p_exp` | `p_case` | `bonett` | `digby` | `pearson` |
|---|---|---|---|---|
| 0.5 | 0.5 | −0.005 | −0.025 | −0.003 |
| 0.3 | 0.3 | −0.007 | −0.009 | +0.014 |
| 0.3 | 0.1 | −0.012 | +0.056 | +0.081 |
| 0.5 | 0.1 | +0.007 | +0.096 | **+0.122** |

At perfectly balanced margins `bonett` and `pearson` agree to within Monte Carlo
noise — as they must, since Bonett's exponent *c* reduces to 1/2 there and the two
formulas become identical. As the outcome becomes rare the gap opens to a factor
of roughly 17. **The cost of the pre-2.0.1 behaviour was therefore concentrated
exactly where meta-analyses of binary outcomes usually sit: rare events.**

For `lipsey_cooper` the two targets separate cleanly. At ρ = 0.5, n = 300 its bias
against its **own** estimand ranges −0.044 to +0.055, while its bias against the
**population** tetrachoric ranges −0.047 to −0.159. Almost all of the latter is
estimand mismatch, not computational error, and reporting it as bias — which a
single-reference benchmark would do — would be wrong.

### Study 08 result: the simulation and the IPD paper agree, independently

`q = SD_baseline / SD_endpoint`. True SMD on the endpoint-SD scale = 0.5, correct
`r_pre_post`, no baseline imbalance, n = 100, nrep = 1200, *g* scale, vs the
population parameter:

| q | **bonett** (default) | cooper / morris_drm | morris_dav | endpoint |
|---|---|---|---|---|
| 0.5 | **+0.501 (cov 0.52)** | +0.076 (0.92) | +0.133 (0.88) | −0.001 (0.95) |
| 0.7 | **+0.215 (0.82)** | +0.060 (0.93) | +0.079 (0.92) | −0.001 (0.95) |
| 1.0 | +0.000 (0.95) | −0.002 (0.95) | −0.000 (0.94) | −0.001 (0.95) |
| 1.4 | **−0.143 (0.81)** | −0.101 (0.90) | −0.089 (0.91) | −0.001 (0.95) |

Two things. **At q = 1 every method is unbiased with nominal coverage** — that was
the legacy DGP's only cell, so the old design was structurally incapable of
finding anything. And **`pre_post_to_smd = "bonett"`, the shipped default, is the
worst performer at both ends**: at q = 0.5 it doubles a true SMD of 0.5.

The PETRA IPD paper finds median q = 0.64 in real ADHD trials and indicts Bonett
as the most divergent estimator. That is the same conclusion from a completely
independent evidence stream — which is exactly the division of labour argued
above: the IPD establishes that q ≈ 0.64 exists in the world, the simulation
establishes what it costs. Together they make the case for changing the default
that neither makes alone. See open item 5.

### Study 02: `cor_to_smd`, verified against its source papers

The three options are three **different estimands**, each correctly implemented.
Sources in `../data-raw/R to SMD/`.

**Mathur MB & VanderWeele TJ (2020)**, *Epidemiology* 31(2):e16–e18 — its two
equations *are* two of the shipped routes:

| paper | code |
|---|---|
| Eq (1.1) `d = 2r/√(1−r²)` | `cor_to_smd = "cooper"` |
| Eq (1.2) `d = rΔ/(s_x√(1−r²))`, `SE = abs(d)·√(1/(r²(N−3)) + 1/(2(N−1)))` | `cor_to_smd = "mathur"` — term for term, SE included |

The paper states plainly that Eq (1.1) "was derived for a Pearson correlation
computed between a **binary** exposure X and a continuous outcome Y, also called a
point-biserial correlation". So `cooper` is scoped to genuine two-group data, and
`viechtbauer` (`metafor::transf.rtod`, the *biserial* map) is the route for a
dichotomised continuous exposure. Neither approximates the other.

It also supplies an exact consistency check worth asserting: applying Eq (1.1) to a
continuous-X correlation "coincides with the effect size associated with an increase
in X of two standard deviations", so `mathur` at `unit_type = "sd"`,
`unit_increase_iv = 2` must reproduce `cooper`'s point estimate exactly
(`r·2s_x/(s_x√(1−r²)) = 2r/√(1−r²)`) — while the two SEs "will, in general, still
not coincide".

**Correction to the note further down this file:** `unit_type = "raw_scale"` is
*not* a silent arithmetic error. With `unit_type != "sd"` the code sets Δ to raw
units of X, which is exactly the paper's Δ. The defect is **documentation only** —
`?convert_df` says the argument "must be either 'sd' or 'value'" while the shipped
default is a third spelling and nothing validates it, so a typo silently selects
raw units. Fix the docs or the default; the arithmetic is right.

**Group balance — a scope limit, not a bug.** Both routes are written at p = 0.5.
At the population value (no sampling noise), true two-group d = 0.8 with r the
matching point-biserial:

| split | `cooper` | general inverse `r/√((1−r²)p(1−p))` |
|---|---|---|
| p = 0.50 | **0.800** | 0.800 |
| p = 0.15 | **0.571** | 0.800 |

For `cooper` that is faithful to Eq (1.1) as published.

**Do not "fix" this by exposing p without reading Pustejovsky first.** My initial
recommendation — forward `n1i`/`n2i` to `transf.rtod`, since
`es_from_pearson_r()` already takes `n_exp`/`n_nexp` — is a *design-dependent*
methodological choice, not a bug fix, and Pustejovsky disputes it.
`data-raw/Converting-from-d-to-r-to-z PUSTEJOVSKY.pdf` identifies the balanced case
exactly: "The minimum possible value of the factor `a = (n1+n2)²/(n1n2)` is achieved
when the two groups are of equal size. In this case, **a = 4**, which is equivalent
to assuming that the absolute value of the treatment-control differential is w = 2."
(`a = 4` ⇔ `d = 2r/√(1−r²)` ⇔ `cooper`.) But he then argues against the
`a`-based generalisation for experiments: the Borenstein and Hunter–Schmidt formulas
"replace w² with the factor `a`… **However, in controlled experiments, these
proportions are arbitrary (and often equal) and do not provide any information about
the magnitude of the treatment-control differential** on the continuous variable X."

So the right split is by design, and it is a genuine open question, not a defect:
- **dichotomisation / extreme-groups** (groups formed by cutting a continuous X):
  the proportions carry real information, and using p is appropriate.
- **controlled experiment** (allocation ratio is an arbitrary design choice): using
  p imports information that is not there; Pustejovsky parameterises by `w`, the
  treatment–control differential on X, instead.

Note metaConvert already handles unequal groups in the **reverse** direction —
`smd_to_cor = "lipsey_cooper"` computes `a <- ((n_exp+n_nexp)^2)/(n_exp*n_nexp)` and
uses `r = d/√(d²+a)`. So the asymmetry is r→d only. Study 02 should measure the cost
of `a = 4` under each mechanism and report it as a design-conditional
recommendation, with Pustejovsky's `w` parameterisation as a third candidate — not
declare a winner.

(The author was already aware of this: `data-raw/R to SMD/R to D.docx` bookmarks
"convert correlation r to cohens d unequal groups of known size".)

### Study 06 result: the shipped OR-SE imputation is systematically too wide

When a paper reports an OR and the case/control margins but no CI, `es_from_or()`
enumerates every compatible 2×2 and averages the logOR **variance**. The variance
is convex in the cell counts, so the mean is dragged up by near-degenerate tables.
Measured against the SE the realised table would have given (br = 0.5, OR = 1,
n = 400):

| estimator | median imputed/true SE |
|---|---|
| shipped mean-of-variance | **1.72×** |
| median-of-variance (one-word change) | **1.15×** |

An SE 1.7× too wide under-weights that study by ~3× in an inverse-variance
meta-analysis. The ordering reverses in very thin cells, so this is a region map
rather than a drop-in replacement — but the default is clearly not the right
choice in the data-rich regime. Study 06 also shows the `4/√N` constant that
**D2's cross-row SE-outlier fallback divides by** is off by ~4× at br = 0.02, so
rare-outcome pools remain bimodal under the current normalisation.

### Study 01 result: both `smd_to_cor` routes are correct — the difference is estimand

`viechtbauer` targets the **biserial** correlation (the latent continuous variable);
`lipsey_cooper` targets the **point-biserial** (the binary grouping variable).
Scored against each route's own estimand (nrep = 400, ρ = 0.5, n = 100, p = 0.5):

| scale | method | bias | MCSE | se_ratio | coverage |
|---|---|---|---|---|---|
| r | lipsey_cooper | −0.003 | 0.004 | 0.99 | 0.94 |
| r | viechtbauer | +0.002 | 0.005 | 0.98 | 0.94 |
| z | lipsey_cooper | −0.001 | 0.005 | 0.99 | 0.96 |
| z | viechtbauer | +0.005 | 0.005 | 0.97 | 0.96 |

Both are unbiased with nominal coverage. Scored against the *other* estimand,
each picks up ~0.10 of apparent bias — which is what the legacy pipeline reported
as Lipsey-Cooper "bias". Consistent with Jacobs & Viechtbauer (2017), so this
doubles as an implementation check.

⚠️ **The two `z` outputs are on different transforms and must not be pooled.**
`lipsey_cooper` returns Fisher's z, `atanh(r)`. `viechtbauer` returns Jacobs &
Viechtbauer's *variance-stabilising* transform
`z = (a/2)·log((1+a·r)/(1−a·r))`, `a = √dnorm(qnorm(p))/(p(1−p))^¼`, with
`vz = 1/(N−1)` rescaled by `vd/vd_crude` (`R/internal_multiple_formulas.R:844-851`; the rescaling is exactly 1 on crude two-group rows, so this study is unaffected — see NEWS.md 2.0.1). So under
`measure = "z"` the value returned depends on which route fired, and a review
whose rows resolve to different routes pools incommensurable quantities. This is
not documented in `?convert_df` and is a candidate metaDETECT check.

### Study 04/05 result: a no-go map for OR ⇄ RR

Scored against the population log RR (nrep = 1000, n = 300, correct baseline risk),
**after the 1.1 / 1.8 / 1.9 / 1.10 fixes**:

| | br = 0.50 (common) | br = 0.15 |
|---|---|---|
| metaumbrella_cases | **−0.011** (cov 0.95) | **−0.018** (cov 0.96) |
| metaumbrella_exp | **−0.011** (0.95) | **−0.018** (0.96) |
| metafor conv.2x2 | −0.011 (0.95) | −0.016 (0.96) |
| dipietrantonj | −0.012 (0.95) | −0.018 (0.96) |
| grant / Zhang-Yu | −0.018 (0.95) | −0.025 (0.96) |
| **vanderweele √OR** | **+0.141 (0.94)** | +0.153 (0.90) |
| transpose (OR as RR) | **−0.311 (0.73)** | −0.030 (0.95) |

Three headlines, and the first two are changes from the pre-fix table.

**`metaumbrella_exp` is no longer the worst performer; it is now identical to
`metaumbrella_cases`.** Both previously searched for the 2x2 table; both now solve it,
so they return the same table and the same answer. `metaumbrella_exp`'s old coverage of
0.22 at br = 0.15 was the non-identified rotation (roadmap 1.1/1.8), and its residual
bias after that fix was continuity-corrected tables being rounded to integers
(roadmap 1.10). The sentence this README used to carry — "the package offers two
`metaumbrella` variants with very different reliability and says nothing about it" —
no longer describes the package.

**VanderWeele's √OR still needs no baseline risk, but its interval is now the
conservative one the paper actually prescribes** (roadmap 1.9): the previous arm built
a symmetric interval, which VanderWeele (2020, p.748) explicitly rules out. Its point
estimates are unchanged. Within the band the method is defined for — both outcome
probabilities in [0.2, 0.8] — coverage is **0.996** at width 1.299; outside it,
**0.941**. Read its coverage in those two regions separately, never pooled.

**`transpose` remains the one clearly disqualifying choice for a common outcome**
(bias −0.311, coverage 0.73 at br = 0.50).

On misspecification (true br = 0.15, analyst guesses 0.50), **only Grant moves**
(bias −0.033 → −0.293); `metaumbrella_*` and `dipietrantonj` are unchanged because
they never read `baseline_risk`. The legacy design had the same property but
presented all lines as if they were being stress-tested.

⚠️ **Read this table with 2c above.** The bias column is allocation-free — the 2×2
solve is exact at every `p_exp` — but the coverage figures beside it are not: both
studies run at 1:1 only, which is the allocation where zero cells are rarest.

### Reference result from study 03 (nrep = 2,000; ρ = φ = 0.5, n = 100, event rate 0.3)

Each method is near-unbiased under its own data-generating mechanism and badly
wrong under the other. Bias is against the population parameter; MCSE ≈ 0.002–0.005
throughout, so every contrast below is resolved by two orders of magnitude.

| mechanism | method | route | bias | se_ratio | coverage | RMSE |
|---|---|---|---|---|---|---|
| **CAT** (φ is the estimand) | phi (r) | candidate | **−0.002** | 1.21 | **0.983** | 0.079 |
| | tetrachoric (r) | package | **+0.260** | 0.99 | **0.262** | 0.277 |
| **CONT** (latent ρ is the estimand) | tetrachoric (r) | package | **−0.002** | 0.97 | **0.918** | 0.140 |
| | phi (r) | candidate | **−0.188** | 1.02 | **0.448** | 0.211 |

The separation the dual-target design was built to expose is visible here:
`se_ratio ≈ 0.97–0.99` for tetrachoric under **both** mechanisms — its standard
error is well calibrated — while coverage collapses to 0.26. The failure is
purely estimand, not variance estimation. A single-target design reports this as
"bias" and cannot tell the two apart.

**The actionable consequence:** metaConvert offers only `"tetrachoric"`, and gives
the user no way to signal which mechanism applies. Either re-open
`table_2x2_to_cor` with a correct φ implementation, or document plainly that the
route assumes a continuous latent variable and should not be used for genuinely
dichotomous outcomes.

⚠️ **If you un-comment the `# TO DO` Lipsey block** in
`R/internal_multiple_formulas.R:484-508`, note that its φ denominator has the same
error as the legacy simulation: the last factor is `(n_controls_exp + n_cases_nexp)`
where φ requires the controls margin `(n_controls_exp + n_controls_nexp)`.
`phi_candidate()` in `studies/03_2x2_to_cor.R` is the corrected version.

### Package bugs found by building these studies

Both were found by running real data through the package rather than a
reimplementation of it — the argument for the "call the package" rule.

1. **`dipietrantonj` aborted the entire `convert_df()` run** on common-outcome
   data (fixed, `R/internal_multiple_formulas.R`). When `estimraw::estim_raw()`
   cannot reconstruct a 2×2 (negative discriminant → NaN cells, routine at
   br ≈ 0.5), `which.min()` returned `integer(0)` and `estim[[integer(0)]]` raised
   *"attempt to select less than one element"*. One unreconstructable row killed
   the whole dataset. Now returns NA for that row. Present in **both** the
   `or_to_rr` and `rr_to_or` branches; both guarded.
2. **The `table_2x2_to_cor = "lipsey"` `# TO DO` block carries a φ denominator
   error** (comment corrected in place, code still disabled). Its last factor was
   `(n_controls_exp + n_cases_nexp)` — a diagonal — where φ requires the controls
   margin `(n_controls_exp + n_controls_nexp)`. Not live, so no user is affected,
   but it would have shipped the moment someone un-commented it. The same error is
   in the legacy simulation and is what produced the "Lipsey bias" in the archived
   aggregates.

### `cor_to_smd`: performance fixed, and one slot-assignment change awaiting a decision

**Performance (done, no behaviour change).** `.cor_to_smd()` was reached through a
per-row `mapply()` at four call sites, and its `"viechtbauer"` branch — the
**default** — built a fresh `metafor::conv.delta()` call per row although
`conv.delta` is vectorised. Memoisation does not help here (correlations are
continuous, so rows are all distinct), so `.cor_to_smd_vec()` groups rows by
method and computes each group in one call: **3.22 s → 0.11 s per 4,000 rows
(29×)**. Estimates bit-identical; SEs agree to ~1e-15 (`conv.delta` differentiates
numerically, so a vectorised call lands ~1 ULP away).

**Slot assignment: investigated, and the shipped convention is correct.** It looks
inconsistent — `cooper` puts its map value in `d` and derives `g = J·d`, while
`viechtbauer` puts `transf.rtod` in `g` and backs out `d = g/J`. The tempting
argument is that `transf.rtod` takes no `n`, so it is a population map and its
output "must" be an uncorrected `d`.

That argument is wrong, and simulation settles it. Which slot a value belongs in
is a question about **estimator bias**, not about the algebra of the transform.
Bivariate normal, ρ = 0.5, median split, true δ = 0.870126, nrep = 2×10⁵ — bias of
the reported *g* against δ:

| n | Hedges *g* from raw data (truth) | `g = transf.rtod` (shipped) | `g = J·transf.rtod` |
|---|---|---|---|
| 10 | −0.004 | **+0.033** | −0.055 |
| 15 | −0.002 | **+0.019** | −0.034 |
| 25 | −0.000 | **+0.012** | −0.017 |
| 50 | −0.001 | **+0.005** | −0.008 |
| 200 | −0.000 | **+0.002** | −0.002 |

The shipped convention is closer to unbiased at every n. The reason: this route's
small-sample bias is **not** the Hedges bias. `E[r̂]` is biased downward, which
partly cancels the upward bias the pooled-SD denominator induces in a directly
computed *d*, so applying J on top over-corrects. Neither convention is exactly
unbiased — the residual is the uncancelled remainder, and it is small (< 0.02 for
n ≥ 25).

Two consequences. `tests_save/checked/test-ES-COR.R` (VIECHT block), which
asserts `g == conv.delta(transf.rtod)`, is **right** and was worth trusting. And
the `d`/`g` pair is internally consistent (`g = J·d`) on every route, so
`measure = "d"` and `measure = "g"` both return what they should. Recorded in
`R/internal_multiple_formulas.R` and locked by
`tests/testthat/test-cor-to-smd-d-g.R` (38 assertions) so the same "obvious" fix
is not attempted again.

### Method-set gaps a referee will find immediately

- **OR → RR omits VanderWeele's √OR and Zhang–Yu.** Grep the legacy scripts for
  `vanderweele|zhang|yu`: nothing. VanderWeele's is designed precisely for readers
  who lack the baseline risk — i.e. the novelty claim of area 04.
- **No `metafor::conv.2x2`**, the CRAN-standard reconstruction, even though
  metaConvert already depends on metafor. Adding it turns a weakness into the
  selling point: *metaumbrella*'s bespoke grid search vs the established
  implementation is a question every user of both packages has.
- In the legacy OR→RR design `br_guess` perturbs **only** Grant; every other line
  is pure replication being read as robustness.

---

## Defects inventory (legacy code, before the rebuild)

Recorded so the rebuild does not silently reproduce them, and so no number from
`archive/` is quoted without knowing what is wrong with it.

**Manufacture the conclusions — no legacy 2×2/OR→COR number is usable:**
- Lipsey phi uses the margin `(b+c)` where phi needs `(b+d)`
  (`6_SIM_2x2_to_COR_CAT.R:99-101`, `_CONT.R:62-64`). In the CAT script the
  benchmark *is* phi, so a correct Lipsey would be near-unbiased — the bug creates
  the reported bias.
- Cooper's `vd = sqrt(or_se² · 3/π²)` (`CAT:83`) takes the square root of a
  variance. This *is* the `bias_var` of 3.36–7.39, and hence Cooper's spurious win.
- z-coverage compares the lower bound to `z_raw_data` and the upper to
  `r_raw_data` (`CAT:292-293`, `CONT:248-249`).
- `psych::phi()` defaults to `digits = 2`, quantising the CAT benchmark to ±0.005
  against an r-grid of step 0.1.
- The `uniroot` escape hatch returning `9e9` (`CAT:45-53`) makes `rbinom` draw with
  prob > 1 and return NA for **100%** of replications in the (r=0.5, br=0.7) block —
  the 30 rows missing from the shipped CAT file. Study 03 replaces it with a
  closed-form attainability check that *reports* dropped conditions.

**Structural:**
- `8_COR_to_SMD.R` is non-runnable *and* dangerous: its export lines still write
  `MEANS_ADJ_to_SMD_sim2500d.txt` and `MEANS_ADJ_to_SMD_AGG10000.txt` — the file the
  app serves. It also correlates `grp_exp` against `grp_nexp`, two independent
  samples, so its benchmark is noise.
- `7_SIM_MEANS_PRE_POST_to_SMD.R` never ran (nrep = 100, no output). Its benchmark
  divides by the residual SD with the marginal back-transform commented out
  (`:45`), so it is on a different scale from every competitor. Prefer the
  `_SAVE` variant, which is despite its name the **newer, corrected** file.
- That DGP also draws both timepoints with unit variance in both arms, so
  **q = SD_baseline/SD_endpoint ≡ 1 by construction** — and Bonett's and d_av's
  deviation factors are both exactly 1 at q = 1. It asks the IPD paper's question
  in the one cell where the answer vanishes. Any rebuild must vary q.
- Silent method-specific missingness: Grant's transform is undefined whenever
  `rr · br_guess ≥ 1`, dropped by a bare `na.rm = TRUE`. Non-estimability is now a
  reported performance measure (`nonest_rate`), not a hidden NA.
- Variance denominators pooled over `p`, which is absent from the `group_by`
  (`2_SIM:396-397` and four others), inflating the `bias_var` denominator.

**App (`app/app.R`) — replaced.** The defects below were in the legacy app that
had been copied into `app/` and deployed at `ebiact.shinyapps.io/simulations`.
That file read none of the rebuilt aggregates and has been rewritten (see
"The results viewer" below); the list is kept as the record of what was wrong.
- Ranked on `abs(bias_ci - 0.95)`, which rewards miscalibration (see above).
- Axis label said "converted − generated"; the generated parameter is the one
  quantity the code never used.
- `:445` read `MEANS_PRE_POST_to_SMD_AGG2500a.txt`, which does not exist.
- Four `(z)` branches omitted `res$es` while offering `es` as a grouper — a crash path.
- `www/` was empty, so the hero image 404'd on every load.
- No LICENSE, CITATION, DOI, or link back to metaConvert; a paper could not cite it.

### The results viewer (`app/app.R`)

Rewritten against the rebuilt aggregates. Run it from the `simulations/` root:

```r
shiny::runApp("app")
```

It is **generic**: every file matching `data/aggregated/*_nrepN.csv` is discovered
at start-up, and each study's condition columns, methods, targets and available
replication counts are read from the file itself. A new study appears in the app
as soon as it has been run — nothing in the app is hard-coded to a study. (The
one exception is the table of human-readable study titles, and a study missing
from it still works; it simply shows its file stem.)

The page reads top-down as **study → answer → evidence**.

A **study header** names the conversion, the package argument it fixes, the
question it answers, and the size of the grid. Below it, a **headline row** of
three cards states the result before any chart is read:

| card | what it says |
|---|---|
| **Leads on this measure** | the method with the smallest mean deviation from the no-error value, over the conditions currently shown, with its average |
| **Its least favourable condition** | that method's single worst cell — what a review should plan for, rather than the average |
| **Spread across methods** | best-to-worst gap, set against the resolution limit (2 × median MCSE), and labelled **resolvable** or **within simulation noise** |

That third card is the one that stops a difference smaller than the
simulation's own noise being quoted as a finding. It is also where the
programme's central result becomes visible in one glance: in study 01a the
spread against the population parameter is 0.176 and resolvable, while against
each method's own estimand it is 0.001 — *within simulation noise*. The methods
are equally good at their own jobs and far apart on the quantity the review
actually wants, which is the estimand-mismatch argument in two numbers.

The headline is suppressed, and replaced by a reason, wherever it cannot be
read: against a `*_sample` target for coverage or the SE ratio, and on a study
that did not record the selected measure (study 06 imputes a standard error, so
it has no interval and no coverage — the note names the measures it *does*
carry). On an `own` target the lead card is relabelled "most accurate on its own
estimand", because "leads" would be a claim about methods held to different
standards.

Then four panels:

| panel | what it answers |
|---|---|
| **Plot** | the chosen performance measure against any condition, faceted by any other, one line per method, with ±1 Monte Carlo SE error bars and a dashed reference line at the no-error value (0 for bias, 0.95 for coverage, 1 for the SE ratio) |
| **Ranking** | a scorecard: bias, absolute bias, RMSE, SE ratio, coverage and worst-case coverage for every method side by side |
| **Estimand check** | one row per method with two markers — mean absolute bias against its **own** estimand, and against the **population** parameter — joined by a bar. The bar length *is* the estimand gap |
| **Data** | the filtered rows, searchable, with a CSV download |

Three things carry across the whole page. Method **colours are fixed by the
study, not by the current selection** — they come from the full method list, so
unticking one does not recolour the others and two screenshots of the same study
are comparable. The same swatch appears in the sidebar checkboxes, in the key
above the plot, and as the left edge of each ranking row, so the ggplot legend is
switched off and its space returned to the panels. And a **provenance footer**
names the source CSV, its aggregation date, and the fact that every number came
from calling the package rather than from a re-implemented formula.

The long reading guide — what the targets mean, why coverage is the default, how
the resolution band is built — sits behind the **How to read this** button rather
than permanently in the sidebar, where it cost about 200px of scroll on every
visit and was read once.

#### Why coverage is the default measure, and why the ranking shows everything

Bias alone is a poor summary, because a conversion can be nearly unbiased in the
point estimate and badly miscalibrated in its standard error — which matters more
in meta-analysis than almost anywhere else, since the standard error sets the
inverse-variance weight. **Coverage is the only ADEMP measure that degrades under
either failure**, so it is the app's default.

The two orderings genuinely disagree in these data. In study 09a, scored on
`own`, the *worst* method on bias is the best on coverage:

<!-- pinned: study09a-ranking -->
| method | rank by \|bias\| | coverage | rank by coverage |
|---|---|---|---|
| `2x2_tetrachoric` | 2nd | 0.951 (worst 0.939) | 1st |
| `bonett` | 1st | 0.964 | 3rd |
| `lipsey_cooper` | 5th | 0.945 | 2nd |

> This table previously read "the best method on bias is the *worst* on coverage",
> with `2x2_tetrachoric` at 0.928 (worst 0.849) and ranked 5th. Those figures matched
> no slice of the shipped aggregate — not `own`, not `population`, not `sample` — and
> the claim they supported was false on the data as shipped. The disagreement between
> the two orderings is real and is why the panel prints both; it simply runs the other
> way round. Ranking is by \|coverage − 0.95\|.
>
> Coverage here is the **identified-conditions** figure of the table above, for the
> reason given in its note: averaging in the conditions where `2x2_tetrachoric`'s
> interval is declaring "not identified" would rank the method on the width of a
> refusal. On that basis it is first on coverage (0.951 against a nominal 0.95) and
> tied-first on bias, and `bonett` is the reverse — so the two orderings still
> disagree, which is the point the panel exists to make.

Across the shipped studies, "best on bias" and "best on coverage" pick the same
method in only **3 of 10**, and the rank correlation between the two orderings runs
from −0.89 to +0.89 — it is not reliably positive, and in two studies it is strongly
negative:

<!-- pinned: rank-disagreement -->
| study | methods | Spearman \|rho\| | best on bias = best on coverage? |
|---|---|---|---|
| `02a_cor_to_smd_GROUPS` | 9 | 0.27 | no |
| `02b_cor_to_smd_CONT` | 9 | 0.29 | no |
| `03a_2x2_to_cor_CAT` | 4 | 0.89 | yes |
| `03b_2x2_to_cor_CONT` | 4 | -0.89 | no |
| `07a_ancova_to_smd_d` | 4 | 0.60 | no |
| `07b_ancova_to_smd_g` | 4 | 0.80 | yes |
| `08a_pre_post_to_smd_d` | 6 | -0.59 | no |
| `08b_pre_post_to_smd_g` | 6 | 0.24 | yes |
| `09a_or_to_cor_CONT` | 10 | 0.23 | no |
| `09b_or_to_cor_CAT` | 10 | 0.21 | no |

Spearman ρ between each method's mean |bias| and its **mean per-condition**
|coverage − 0.95|, on the `own` target, over the studies carrying at least three
methods — with two, ρ is ±1 by construction and says nothing. So no single number is
trusted: the Ranking panel prints all of them and the selected measure only decides
the sort order.

> **This paragraph used to quote three numbers in prose, and all three had rotted.**
> It read "08a: −0.59, 09a: −0.23, 09b: −0.29 ... only 3 of 12". Only the first was
> right. `09a` and `09b` are **positive** (+0.23, +0.21) — the sign appears to have
> been taken from the sentence's own claim that the correlation is "sometimes
> negative" rather than from the data, which inverts the point being made, since a
> positive ρ means the two rankings *agree* in those studies. `09b`'s 0.29 was the
> magnitude of the pre-tetrachoric-fix aggregate; `09a`'s 0.23 matched neither era
> when it was written. And "3 of 12" had the right numerator against the wrong
> denominator. Nothing detected any of it, because prose is not recomputed — which
> is why the figures are now a pinned table with a recipe in
> `R/08_published_numbers.R`, like every other published table here.

**The sort column is shown explicitly**, because it is not the deviation of the
mean beside it. It is the mean deviation *per condition* — a method covering 0.90
in half the conditions and 1.00 in the rest averages to a flawless 0.95 while
being miscalibrated in every one of them. Averaging the per-condition deviation
catches that; deviating the average does not.

Coverage and the SE ratio are meaningless against a `*_sample` target (a nominal
95% interval is not built to cover a random quantity), so the panel says so rather
than printing a number that cannot be read.

Three defaults are deliberate, because each prevents a misreading:

- **The Monte Carlo resolution band.** Every panel carries a shaded band of
  ±2 × the median MCSE of its own cells — the `k = 2` convention of
  `mc_resolvable()` in `R/03_performance.R`. A difference that stays inside the
  band is smaller than the simulation's own noise and is not a result. The band
  is computed **per panel, not once per plot**: MCSE shrinks with sample size, so
  a single global band would be far too narrow at small *n* and too wide at large
  *n*, misstating resolvability exactly where the question matters. Watching it
  contract from *n* = 25 to *n* = 300 is the point.
- **A shared vertical scale.** Free scales let each panel pick its own range, so
  the axis labels collide and panels cannot be compared. There is a checkbox for
  the rare case where free scales are what you want.
- **One reported scale at a time.** Several studies encode the output scale in
  the method name (`bonett (r)` and `bonett (z)` are one formula on two scales).
  Plotting both against one axis compares quantities that are not on the same
  scale, and the larger *z* values flatten the *r* values. Where such suffixes
  exist the app offers them as a switch and shows one group at a time.

The target selector is populated from the file rather than assumed, because the
target names are study-specific: most studies record `own` / `population` /
`sample`, but study 01 records `biserial_*` and `pointbiserial_*` separately and
study 08 adds `own_true_r`.

#### `own` is a diagnostic, not a ranking basis

The targets are two different kinds of thing, and listing them as one flat set of
options invites a real mistake:

- a **shared** target (`population`, `sample`) is *one* quantity applied to every
  method, which is how methods are compared against each other;
- **`own`** is a *different* quantity per method. Scoring on it answers "does this
  method compute its own value correctly?" and never "is that the value I want?".
  Two methods scored on `own` are being held to two different standards.

So the app defaults to the population parameter, groups `own` separately, and
labels it "Each method's own estimand". Ranking on `own` is not wrong, it just
answers a different question, and the panel says so when it is selected.

The difference is not cosmetic. In study 09a, `lipsey_cooper` has a mean absolute
bias of **0.021 against `own`** and **0.086 against the population parameter**,
with worst-case coverage falling from 0.897 to 0.482 — because it is computing a
point-biserial correctly while the analysis asked for a tetrachoric. Both numbers
are true and they answer different questions; only the second one tells you not to
use it here.

Estimand gaps are the norm rather than the exception, which is why `own` earns its
place: 10/10 methods in `09b_or_to_cor_CAT`, 6/6 in `08a`/`08b_pre_post_to_smd`,
8/9 in `02a_cor_to_smd_GROUPS`, 6/10 in `09a_or_to_cor_CONT`.

Visual identity follows metaconvert.org: rose accent, indigo for data, warm
neutrals, Work Sans with Cascadia Code for every number. The brand contains no
green, so the categorical scale for methods is green-free by construction and is
ordered so that neighbouring slots are far apart in hue (methods are taken
alphabetically, and study 09's `2x2_tetrachoric` and `lipsey_cooper` otherwise
landed on two near-identical reds).

Requires `shiny`, `bslib`, `ggplot2` and `DT`. No `plotly`, no `tidyverse`.

Dashboard headline numbers do not reconcile: the "41,000,000 datasets" card traces
to `archive/legacy_scripts/summary_files.xlsx`, whose two component errors are each
exactly 6,750,000 and cancel — the total is right by accident. "10 simulations" and
"33 formulas" match neither the registry (7 studies with output) nor the files
(48 distinct method labels).

---

## Running it

```r
setwd("simulations")
source("run_all.R")          # loads harness + studies, prints what is available

run_03(nrep = 300, cores = 1)              # quick check of one study
run_03()                                   # production: nrep = 10000, all cores
run_everything(nrep = 10000)               # everything
```

Outputs: `data/raw/<study>_raw_nrep<N>.rds` (per replication) and
`data/aggregated/<study>_nrep<N>.csv` (one row per condition × method × target,
with `bias`, `bias_mcse`, `emp_se`, `mod_se`, `se_ratio`, `rmse`, `coverage`,
`coverage_mcse`, `nonest_rate`, `n_valid`, `n_valid_se`, `n_valid_ci`).

**Read `n_valid_se` and `n_valid_ci`, not just `n_valid`.** A row's columns rest on
three different subsets of the replications: `bias`/`emp_se`/`rmse` on those whose
point estimate is finite, `mod_se`/`se_ratio` on those whose standard error is, and
`coverage`/`ci_width` on those whose interval is. A method can return a perfectly
good point estimate with a non-finite variance for most of a cell — `grant` does,
whenever `rr × br_guess ≥ 1` — so these counts differ, and only the first used to be
reported. Statistics resting on fewer than `SIM_DEFAULTS$min_valid` (30) valid
replications are now withheld as `NA`, per family, keyed on that family's own count;
the counts are always reported, so a withheld cell is explainable and a reader can
apply a different floor. Note also that the surviving replications are in general a
*selected* subset, so a statistic computed over a handful of them is conditional on
the estimator not having failed.

**Monte Carlo error is not the binding constraint.** MCSE on coverage at p=0.95 is
0.0044 at nrep=2500 and 0.0022 at 10000; the contrasts of interest are far larger.
Do not spend compute raising nrep — spend it on the design factors that are
currently missing (unequal arm SDs, q ≠ 1, non-normal outcomes).

**Requirements:** R ≥ 4.5, `pkgload`, `MASS`, `parallel`, and metaConvert's own
dependencies. The archived scripts additionally need `MonteCarlo`, which was
archived from CRAN in 2019 and is not installable — another reason the rebuild
does not use it.

---

## Open items before this becomes a paper

**Ordered, with the blocking one first.**

0. **The reporting plan (≤6 display objects) is still unwritten, and it blocks the
   production run.** 40k rows at `nrep = 1000` is already untabulatable; a production
   run without a plan produces the same pile at 10× the cost. Decide the displays
   *first*. The genre standards are nested-loop plots for bias and zip plots for
   coverage.
0b. Regenerate study 03 (stale, see STATE OF PLAY), then diagnose the binary
   `se_ratio > 1` pattern, then build a `pool_sd` study for the ONE unsimulated
   kernel (not two — see THREE THINGS #3; and not "study 08b", an id study 08
   already uses).
0c. The app is **broken**: all 16 of its reads point at `./data_agg/`, which the
   restructure moved to `data/aggregated/`, and it still ranks on
   `abs(bias_ci - 0.95)` at two sites and reads a schema the new studies do not
   produce. Retire it or rebuild against the new CSV after the production run.
0d. `simulations/` is still in `.gitignore:30` — **this work is unversioned**.
   Un-ignore it; keep `^simulations$` in `.Rbuildignore`.


1. **Reporting plan.** Nothing yet specifies the ≤6 display objects. 9,330 legacy
   rows are untabulatable; decide on nested-loop plots / zip plots for coverage
   **before** the production run, or it produces the same unpublishable pile.
2. **Pre-specify** the primary comparison per area and split confirmatory from
   exploratory. The IPD paper's reviewer raised exactly this at ~200× smaller scale.
3. **Declare the COI**: this evaluates the authors' own packages (metaConvert,
   metaumbrella). That is legitimate and interesting, but must be stated.
4. **Verify every citation against the PDF.** The sibling paper's review found a
   fabricated reference; the same referee pool applies here.
5. **Two shipped defaults are in tension with the evidence** and need resolving
   either way: `pre_post_to_smd = "bonett"` (which the IPD paper indicts as the most
   divergent estimator in a realistic evidence base) and `table_2x2_to_cor`
   accepting only `"tetrachoric"` (whose assumption the CAT/CONT split shows is
   consequential, with no argument exposed to the user to act on it).
6. **Doc inconsistencies to fix regardless**:
   - `prop_to_es = "raw"` vs `vignettes/Psychometrics.Rmd` calling `freeman_tukey`
     "recommended".
   - `or_to_cor` defaults to `"bonett"` in `convert_df()` but `"pearson"` in the
     exported `es_from_or_se()`, so the same data converts differently depending on
     the entry point.
   - `unit_type` (verified): `?convert_df` says it *"Must be either 'sd' or
     'value'"*, the default is `"raw_scale"` — neither — and the only use anywhere
     is `ifelse(unit_type == "sd", ...)` (`R/internal_multiple_formulas.R:1259`)
     with **no validation**. So the default happens to behave as `"value"`, and any
     typo silently does too. Either validate the argument or make the default one
     of the two documented values.
