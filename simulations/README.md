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

## STATE OF PLAY — read this first if you are picking the project up

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

Confirmed at full grid: **`or_to_rr = "metaumbrella_exp"` has coverage 0.497** — the
worst shipped option, and worse than the 0.756 seen on the smaller grid.
`metafor_conv2x2` fails to reconstruct in up to 83% of a cell (negative discriminant).
`grant`'s non-estimability reaches 0.620 in study 05 (`rr × br_guess ≥ 1`), now
counted rather than silently dropped.

### THREE THINGS THE NEXT PERSON SHOULD KNOW

**1. Study 03's `nrep = 1000` files are stale.** They were produced before a fix and
must be regenerated. The bug was in the *simulation*, not the package: the `(z)`
routes were scored against `theta_pop = rho`, an r-scale target, so `atanh(rho) − rho`
was being reported as bias. Proof: the correlation between `phi (z)`'s apparent bias
and the pure scale gap was **0.9929**. After the fix (a scale-matched `theta_own`),
`phi (z)` bias fell 0.0347 → 0.0056 and its coverage rose to 0.9618, identical to
`phi (r)`. `tetrachoric (z)`'s remaining ≈0.23 is the genuine CAT estimand mismatch.
*Lesson worth carrying: whenever a study reports both r and z, check the target scale.*

**2. An unexplained pattern, not yet diagnosed.** In **both** binary studies,
`se_ratio` is systematically **> 1 for every method** (1.13–1.62), including
`transpose`, which merely copies `logor_se` through. A pattern that uniform across
unrelated methods usually points at something structural — the +0.5 continuity
correction, or an OR-SE-vs-RR-spread comparison in `transpose`'s case — rather than a
per-method bug. **This has not been chased. Do not assume it is benign.**

**3. Pre-post coverage is much thinner than the eight-study structure suggests.**
Study 08 exercises all five `pre_post_to_smd` values but only **one of three kernels**:

| kernel | status |
|---|---|
| `.pre_post_to_smd` (two-group, `pool_sd = FALSE`) | tested |
| `.pooled_pre_post_to_smd` (`pool_sd = TRUE`) | **never run** |
| `.single_group_pre_post_to_smd` | **never run** |

19 exported functions take `pre_post_to_smd`; study 08 calls one
(`es_from_means_sd_pre_post`). Untested: the whole single-group family, all eight
`es_from_mean_change_*`, the five `es_from_paired_t/f*`, and the `_se`/`_ci` variants.
Worse, `pool_sd = TRUE` is the route carrying the known docblock defect (leading
variance term misattributed to Bonett eq 19, departing up to 38% at unequal arm
sizes) — and it is the one never exercised. A "study 08b" reusing the same DGP and
grid with a new `estimate()` would take this from 1 kernel to 3.

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
All eight are implemented and run: `source("run_all.R")` registers `run_01` … `run_08`.

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

| method | mean \|bias\| | worst \|bias\| | coverage | SE ratio |
|---|---|---|---|---|
| `2x2_tetrachoric` (full table) | 0.0163 | 0.2701 | 0.928 | 1.042 |
| **`bonett`** | **0.0185** | 0.2571 | 0.964 | 1.051 |
| `digby` | 0.0201 | 0.2103 | 0.966 | 1.062 |
| `pearson` | 0.0223 | 0.1865 | 0.965 | 1.062 |
| `lipsey_cooper` | 0.0231 | 0.1937 | 0.946 | 1.060 |

`bonett` is the best of the four OR-only options, which supports its promotion to
the default. The full-2×2 route still wins, so the advice is unchanged: **convert
the table, not the odds ratio, whenever the table is reported.**

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

Scored against the population log RR (nrep = 600, n = 300, correct baseline risk):

| | br = 0.50 (common) | br = 0.15 |
|---|---|---|
| metaumbrella_cases | **−0.009** (cov 0.93) | **−0.028** (cov 0.96) |
| dipietrantonj | −0.010 (0.93) | −0.028 (0.96) |
| grant / Zhang-Yu | −0.017 (0.94) | −0.033 (0.96) |
| **vanderweele √OR** | **+0.024 (0.93)** | +0.290 (0.64) |
| metafor conv.2x2 | +0.022 (0.94) | +0.290 (0.37) |
| metaumbrella_exp | +0.045 (0.90) | +0.471 (0.22) |
| transpose (OR as RR) | **−0.239 (0.81)** | −0.114 (0.95) |

Two headlines. **VanderWeele's √OR matches the baseline-risk-requiring methods for
a common outcome while needing no baseline risk at all** — exactly its minimax
claim — and degrades outside that regime. And **`metaumbrella_exp`, a shipped
option, is the worst performer in both columns** (coverage 0.22–0.90), consistent
with the legacy files' mean coverage of 0.756; the package offers two
`metaumbrella` variants with very different reliability and says nothing about it.

On misspecification (true br = 0.15, analyst guesses 0.50), **only Grant moves**
(bias −0.033 → −0.293); `metaumbrella_*` and `dipietrantonj` are unchanged because
they never read `baseline_risk`. The legacy design had the same property but
presented all lines as if they were being stress-tested.

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

Two things follow. `tests_save/checked/test-ES-COR.R` (VIECHT block), which
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

**App (`app/app.R`) — REPLACED.** The defects below were in the legacy app that
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
`own`, the best method on bias is the *worst* on coverage:

| method | rank by \|bias\| | coverage | rank by coverage |
|---|---|---|---|
| `2x2_tetrachoric` | 1st | 0.928 (worst 0.849) | 5th |
| `bonett` | 2nd | 0.964 | 2nd |
| `lipsey_cooper` | 5th | 0.946 | 1st |

Across the shipped studies the Spearman correlation between mean |bias| and
|coverage − 0.95| is weak and sometimes negative (08a: −0.59, 09a: −0.23,
09b: −0.29), and "best on bias" coincides with "best on coverage" in only 3 of 12.
So no single number is trusted: the Ranking panel prints all of them and the
selected measure only decides the sort order.

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
  method, which is what lets you compare methods against each other;
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
`coverage_mcse`, `nonest_rate`, `n_valid`).

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
   `se_ratio > 1` pattern, then build study 08b for the two untested kernels.
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
