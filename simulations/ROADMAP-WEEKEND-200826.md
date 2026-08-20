# Weekend handoff — 2026-08-20

Companion to `ROADMAP.md`, which stays the source of truth for the items themselves.
This file is what you need to **pick the work up on another machine**: what is where,
what is blocking, what was measured this session, and what the next item needs.

---

## 0. READ FIRST — the work is not on GitHub

```
branch    audit-remediation
upstream  NONE   <- the durable fact: git status -sb prints "## audit-remediation"
                    with no "...origin/…", i.e. this branch exists nowhere else
remote    origin  https://github.com/CorentinJGosling/metaconvert-package.git
ahead of  master by 35 commits, tip 11ab2f8   <- drifts; check with
                    git log --oneline master..HEAD | wc -l
```

**Nothing has been pushed.** On the other laptop the branch does not exist. Before you
travel, one of these has to happen — **(d) is the chosen route**; (a)–(c) are recorded
because they are what you would need if the folder ever has to travel over a network.

**(a) Push it.**

```bash
git push -u origin audit-remediation
```

⚠️ **That repository is PUBLIC.** The push would publish, among 102 `simulations/`
files: `simulations/README.md` and `simulations/ROADMAP.md`, which are candid write-ups
of defects in the shipped package (wrong 2×2 branch, NaN standard errors, a dead flag,
an app defaulting to a meaningless benchmark). That may be exactly what you want —
it is good, defensible work and the repo is the natural home — but it is a publication
decision, not a mechanical step, so it is yours to make rather than mine.

Checked, so you do not have to: **no `papers/` or `students/` file is tracked anywhere**
(`git ls-files | grep -E 'papers/|students/'` is empty), so the 2026-07 `ready/fix`
incident does not repeat. `.Rbuildignore` **is** tracked and line 7 is `^simulations$`,
so none of this reaches a CRAN tarball. Total pack size is trivial.

**(b) Keep it private — bundle it.** One file, copy it by any means, no server:

```bash
git bundle create ../metaconvert-audit-200826.bundle master..audit-remediation
# on the other laptop, inside a clone:
git fetch ../metaconvert-audit-200826.bundle audit-remediation:audit-remediation
```

**(c) Push to a private remote** you add as a second origin.

**(d) Copy the whole folder — CHOSEN.** Better than (a)–(c), not a fallback: it is the
only route that carries `.git/` *and* the gitignored files, and it needs no publication
decision. Verified before recommending it — the working tree is clean, and the only
absolute path in the codebase is inside a dead comment at `R/internal_check_data.R:193`,
so nothing breaks on a different machine.

- **Copy the folder itself, not its contents.** `.git` is hidden. Opening `metaConvert/`,
  Ctrl+A with hidden items not shown, and copying *that* gives you the files with no
  history and no branch — the one way this goes wrong. Right-click `metaConvert` → Copy.
- **786 MB** total: `.git` 163 MB, `simulations/data/raw/` 313 MB. Dropping `data/raw/`
  saves 313 MB and loses nothing you cannot regenerate (seeded, bit-for-bit; the
  aggregated CSVs are tracked).
- **Let any running R process finish first**, or the copy picks up a temp file.
- **The R library does not travel.** The other laptop needs `devtools`, `testthat`,
  `pkgload`, `metafor`, `estimraw`, `compareDF`, `rio`, `esc`, `MASS`, `mvtnorm`, `psych`,
  `MBESS`, `semTools` for the package and suites, plus `shiny`, `bslib`, `ggplot2`, `DT`
  for the app.

#### Copy-back protocol — one authoritative copy at a time

The plan is: work on the other laptop, then delete this folder and drop the returning one
in its place. That is clean and needs no merge, but it has exactly one failure mode:

> ⚠️ **The moment the copy leaves, THIS machine is frozen.** Any commit made here
> afterwards is destroyed by the copy-back, silently — the returning `.git` simply does
> not contain it, and there is no conflict to warn you.

So: no work here while the copy is away. And on the way back, **rename rather than
delete** — `metaConvert` → `metaConvert-OLD-200826` — confirm the returning copy runs
(`git log --oneline -1`, then `cd simulations && Rscript tests/run_tests.R`), and only
then remove the old one. Deleting first turns a bad copy into a lost repository.

### What (d) makes moot — but a push or bundle would leave behind

These are gitignored, so they live only on this disk and contain edits from this session.
A folder copy carries all of them; (a)–(c) do not:

| path | why it matters |
|---|---|
| `simulations/students/` | I added a "not a bug" entry to `task-04-or-to-rr.md` and `task-05-rr-to-or.md` about the single-valued `p_exp` column (roadmap 4.2). |
| `CLAUDE.md` (repo root) | gitignored; whatever it says about the package does not travel |
| `simulations/data/raw/*.rds` | ~315 MB, regenerable — every `(study, condition)` re-runs bit-for-bit from `SIM_BASE_SEED`. The **aggregated CSVs are tracked**, so nothing you need for analysis is missing. |
| `simulations/papers/` | copyrighted PDFs, deliberately excluded |

### Working tree at the moment of writing

**Clean** — `git status --porcelain` empty. Earlier in the session two files were open
(`R/internal_flags.R`, `tests/testthat/test-reliability-generalization.R`, the
reliability/omega workstream); they landed in `7e385fb`. Nothing of mine is left
unstaged. **Re-run `git status --porcelain` immediately before copying** — this drifted
several times while the file was being written.

---

## 1. Environment

```
Rscript   /c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe
R version 4.5.1        platform win32, Git Bash + PowerShell both available
Pandoc    RSTUDIO_PANDOC="C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
```

**Always set `NOT_CRAN=true`.** Without it 12 files / ~1585 assertions silently skip.

### The three suites

```bash
# 1. package, main            (~4 min)
Rscript -e "Sys.setenv(NOT_CRAN='true'); r <- as.data.frame(devtools::test(reporter='silent')); \
            cat('PASS',sum(r$passed),'FAIL',sum(r$failed),'ERROR',sum(r$error),'\n')"

# 2. package, archived reference   (~25 min — the slow one)
Rscript -e "Sys.setenv(NOT_CRAN='true'); pkgload::load_all('.', quiet=TRUE); library(testthat); \
            r <- as.data.frame(test_dir('tests_save/checked', reporter='silent', package='metaConvert')); \
            cat('PASS',sum(r$passed),'FAIL',sum(r$failed),'ERROR',sum(r$error),'\n')"

# 3. simulations — a THIRD suite neither of the above runs   (~1 min)
cd simulations && Rscript tests/run_tests.R
```

### Baselines as of this session

| suite | count | note |
|---|---|---|
| `tests/testthat/` | **3210** pass / 0 fail / 0 error | moved 3178 → 3207 → 3210 *during* the session as your reliability work landed. **Re-measure before trusting it**; it is not a fixed number while two people are in the tree. |
| `tests_save/checked/` | **7198** pass / 0 fail / 0 error | confirmed twice, second run finished after `7e385fb` landed. Was 7193 before this session; the +5 came from the reliability changes in `R/`, not from test edits — `tests_save/` itself is unmodified in git. Takes ~25 min. |
| `simulations/tests/` | **548** pass / 0 fail / 0 error / 0 skip | 438 at session start, +39 (item 4.2) +71 (item 4.3) |

**Check the `error` column, not just `failed`.** An erroring `test_that` block reports
`failed = 0` while abandoning every assertion after the error. Compare the assertion
TOTAL against the previous run: a count that falls with no failures means assertions
stopped executing.

---

## 2. Working agreement (this process is working — keep it)

1. One roadmap item at a time, in the roadmap's execution order.
2. **Before writing code**: state the problem in plain, non-statistical terms plus the
   strategy, and wait for approval. Do not batch items.
3. Every fix gets a test. Package → `tests/testthat/`. Simulation → `simulations/tests/`.
4. Before declaring an item done, run **both** package suites plus the simulation suite.
5. **Verify claims by running code, not by reasoning.** This is not a slogan — see §3.
6. Commit messages record WHY, the measured numbers, and anything the change overturns.
7. **Never `git add -A` / `git commit -a`.** Stage explicit paths, then confirm with
   `git diff --cached --name-only` before committing.

### Gotchas that have cost real time

- **Heredocs collapse `\\` to `\`.** This corrupted a roxygen block (`\rho` → literal CR)
  and a regex (`"\\1"` → octal escape). Write files via Python with explicit `newline=`,
  or use the editor tools.
- **Line endings differ per file.** `simulations/ROADMAP.md` is **CRLF**;
  `simulations/README.md`, `app/app.R` and the R files are **LF**. Normalise, edit, write
  back with the original ending — otherwise the diff shows the whole file as changed.
- `flags = TRUE` belongs to `summary()`, not `convert_df()`.
- `summary()` rounds `es_crude` to 2 dp — compare at that precision.
- **Never wrap a suite in `suppressWarnings`/`suppressMessages`** to tidy output; it
  silently drops assertions (once produced a false 5608-vs-7193 scare).
- Long simulation runs: use an isolated git worktree pinned to a commit, since both of
  us edit the same tree and 26 parallel workers each `pkgload::load_all()`.

### ⚠️ Parallel-work hazard, which actually fired

HEAD moved **twice while this handoff was being written** — `17c5813`, then `7e385fb`
("Scope V31 to measure = icc so it does not leak into alpha/omega runs"). Both are yours.
Treat every commit count and SHA in this file as *true at 2026-08-20, tip `7e385fb`*, and
re-read `git log` before relying on one.

The concrete casualty: commit `17c5813` swept in **four files belonging to roadmap item
4.2** alongside the reliability work. The content is intact and correct; only the commit message is shared. History was
**not** rewritten — the item's full reasoning and every number live in `ROADMAP.md` §4.2
and `README.md` §2c, which is where a reader would look anyway. If you want it split out
later, say so and it can be done non-destructively.

---

## 3. Done this session

Both items had their premise **dissolve or shift under measurement**. That keeps
happening, which is why rule 5 above is the important one.

### 4.2 ☑ Studies 04/05 run at 1:1 allocation only — *documented, not regenerated*

The item said "`metaumbrella_exp`'s 0.497 coverage is conditional on `p_exp = 0.5`; add
`p_exp` to the grid". After items 1.1/1.8/1.10 **0.497 no longer exists** — the shipped
aggregate has mean coverage 0.976 and `metaumbrella_exp` is bit-identical to
`metaumbrella_cases` (max |difference| **0** on bias, coverage and `se_ratio` across all
360 conditions). There was no reporting correction left to make.

What replaced it, all measured without a simulation run:

- **The point-estimate side is allocation-free.** Both routes recover the 2×2 exactly at
  `p_exp` 0.10 → 0.90: 400 draws each, exact-recovery **1.000**, max |log RR error| **0**.
  So the no-go map needs no caveat.
- **`p_exp = 0.5` is the most favourable allocation for item 4.1's sparsity diagnosis** —
  the part that really is conditional. Dense region over the same 360 cells:
  **25 / 65 / 105 / 85 / 45** at `p_exp` 0.10 / 0.25 / 0.50 / 0.75 / 0.90. The shipped
  grid sits at the maximum, and the two directions are not interchangeable.
- **Same mechanism drives 4.1's second caveat.** Enumerating the *whole* 2×2 table space
  at n = 50: `metafor::conv.2x2` reconstructs **every** integer table (failure 0.000 at
  every allocation) and fails on **0.45–0.68** of the continuity-corrected ones. Its
  non-estimability is therefore the corrected-table share times that rate — and the share
  is minimised at balance by a closed-form argument: **0.148** at 1:1 vs **0.362** at 1:9.

**Costing kept, in case you want the run later.** `p_exp = c(0.5, 0.25, 0.75)` would be
strictly additive: `expand.grid` varies the last factor slowest, so the current 360
conditions keep positions 1–360, keep their `condition_seed()` streams and reproduce
bit-for-bit. ~3.4 h wall (`run_04` 41 min + `run_05` 26 min, ×3, on 26 cores).

Landed: `R/06_sparsity.R` gains `dense_region_by_allocation()` and
`corrected_table_share()`; README §2c; `tests/test-allocation-scope.R` (39 assertions),
whose first assertion **fails on purpose** if anyone adds a `p_exp` level.

### 4.3 ☑ App: study 01b's default target — commit `4732667`

Confirmed as written. 01b records no target named "population", so the ordering rule fell
through to `fisherz_biserial` while the sidebar recommended shared benchmarks over `own`.
That target is `atanh(biserial)`, and **neither route produces it**: at ρ = 0.75, p = 0.5
the three quantities are **0.973** (target) / **0.724** (viechtbauer's variance-stabilising
transform) / **0.691** (lipsey_cooper's `atanh` of the point-biserial). Coverage on the
old default was **0.736 / 0.773**, against 0.944 / 0.934 on `fisherz_pointbiserial` —
which is now what the app opens on.

Two corrections from measurement:

- **Studies 03 and 09 look like the same defect and are not**, so they were left alone.
  Their `population` target is attained *exactly* by the route correct for the mechanism
  (03a `phi (r)`: bias 0.0035, coverage 0.962, identical to its `own`), and where nothing
  attains it — 09b, where φ is the estimand and the package ships no φ route — **that gap
  is the study's result**. An earlier draft of the fix would have deleted a finding.
- **The fix cannot be a data-derived rule.** "Every method misses this shared target by a
  margin large and flat in n" is the intended finding in 07/08/09b and the defect in 01b,
  and nothing in the aggregates separates them. Four candidate predicates were tried; all
  four misfired. Hence a declared map, `STUDY_RANK_ON`, beside the `STUDY_LABELS` map the
  app already carries — **one entry**.

Also landed: `unattained_targets()` marks a benchmark *no route estimates this* — derived,
not declared (a route whose own estimand is a named target produces a bit-identical bias
column against both, so exact vector equality decides it with no threshold). Neutral by
design: in 08 and 09b it fires on `population` and states those studies' results. Plus
`.scale_mismatch()`, warning when 03/09 are switched to the `(z)` routes with a shared
target selected.

⚠️ **A regression caught only by running it**: the first `unattained_targets()` labelled
every `*sample*` target — all 12 studies. `own` is a population quantity everywhere, so a
same-sample target can never equal it. Excluded now, with a test pinning it.

**Verified in the live app** via `shiny::testServer()`, not just unit tests: 01b renders
`fisherz_pointbiserial` with `checked="checked"`, the marker on `fisherz_biserial` alone;
01a still opens on `biserial_population`, 09b on `population` (marked); the scale warning
absent at `(r)`, present at `(z)`, absent again on `own`.

---

## 4. NEXT ITEM — 4.4, already measured, approval pending

**Everything below is measured. Do not re-measure it; go straight to implementation.**

### The problem in plain terms

`simulations/README.md:146-150` carries a table of which of three internal calculation
engines the simulation exercises, and marks two as **"never run"**. One of those two *is*
run — it is the engine the two-group route calls once per arm, so study 08 executes it
**twice on every row**.

Instrumented with `trace()` through the exact route study 08 uses
(`es_from_means_sd_pre_post`, `pool_sd` at its default):

| route | `.single_group_pre_post_to_smd` | `.pooled_pre_post_to_smd` |
|---|---|---|
| `pool_sd = FALSE`, 2 rows | **4 calls** (2 per row, one per arm) | 0 |
| `pool_sd = TRUE`, 1 row | 0 | **1 call** |

Dispatch: `R/internal_multiple_formulas.R:1400` and `:1407`.

**The second half of the item also holds.** `.pooled_pre_post_to_smd` genuinely is never
run by the simulation — but it is not untested: `tests_save/checked/test-pooled-variance-
calibration.R` puts **78 assertions** through it across 6 blocks (measured by running that
file alone), including bit-exact agreement with `metafor::escalc` and Monte-Carlo
calibration of the one branch with no closed-form reference. So **"unsimulated" is the
accurate word; "uncovered" is not.**

**A nuance the roadmap does not mention, worth writing in.** The single-group kernel is
executed *only* as the per-arm engine of the two-group route. Its own user-facing entry
points — the `es_from_PAIRED_SINGLE_GROUP` family, which calls it directly at
`R/es_from_PAIRED_SINGLE_GROUP.R:124` — remain unsimulated. The section's *conclusion*
(pre-post coverage is thin) survives; the mechanism is what is misstated. This matters
because it changes what a "study 08b" would buy.

### Plan

1. Rewrite the kernel table with accurate statuses and a "how it is covered" column.
2. Fix two knock-on sentences: `README.md:158` *"would take this from 1 kernel to 3"*, and
   `README.md:1069` *"build study 08b for the two untested kernels"* — it is one.
3. **Test unit**: a `simulations/tests/` file that instruments the kernels with `trace()`
   and asserts which ones study 08's route actually executes, both directions. That turns
   the README claim from prose into something executable, so it cannot rot again — which
   is exactly how it rotted.

---

## 5. The rest of the queue

Execution order from `ROADMAP.md`: **4.4 → 4.5 → 4.6 → 2.6–2.9 → 5.1–5.4**.

### 4.5 — strike stale items · *scope is bigger than the item says*

The item names two: `README.md` open item 0c ("the app is broken, 16 reads point at
`./data_agg/`" — fixed) and the note at `studies/06_or_se_imputation.R:171-173` accusing
study 04 of a `conv.2x2` margin error (the call is correct and verified in a comment
above it). **Three more are stale**, found while reading for 4.4, all in the same
"Ordered, with the blocking one first" list around `README.md:1060-1075`:

- 0b "Regenerate study 03 (stale)" — landed with item 3.5
- 0b "then diagnose the binary `se_ratio > 1` pattern" — landed with item 4.1
- 0d "`simulations/` is still in `.gitignore:30` — this work is unversioned" — landed
  with item 0.1

Sweeping all five together is the sensible move.

### 4.6 — downgrade five findings to documented design consequences

`r_pre_post = 0.8`, `.se_from_or` mean-over-tables, tetrachoric-only, the ANCOVA variance
limitation and the `or_to_cor` divergence are **already analysed in the repo**, so
presenting them as external discoveries undersells it and will not survive a referee.

⚠️ **Re-measure every number before quoting it.** The roadmap lists figures to retire
("median gives 1.006", "1.19→1.70", "18.7×", "all five pre/post routes deflate the SE" —
`morris_dz` does not), and **several have moved again** since that was written, because
items 1.1/1.8/1.10/3.1/3.5 changed the aggregates underneath them.

### 2.6–2.9 — the `or_to_cor` cluster

- **2.6** state the margin-dependence position: `bonett` (default) is the least
  margin-dependent at drift **0.0212**; `pearson`/`digby` are ~6× worse (**0.1245** /
  **0.1205**) — the same pathology φ was refused for, differing in degree not kind, and
  currently undocumented.
- **2.7** the fallback picks the **estimand-mismatched** method. With `n_sample` but no
  margins, bonett→fallback gives r = **0.2449**, where `pearson` gives 0.3463 and
  true-bonett 0.3205 — i.e. the fallback moves *further* from the requested method than
  the estimand-matched alternatives. Pre-existing, affects `convert_df()` too.
  ⚠️ **Now visible in the results**: study 09's `bonett` row is a mixture of two
  estimands — ~3.6% of replications are silently `lipsey_cooper` fallbacks. Verified: at a
  degenerate margin `or_to_cor = "bonett"` returns r = 0.30533, bit-identical to
  `lipsey_cooper`.
- **2.8** the fallback message names `small_margin_prop` first, and that input **cannot
  work on its own** — the eligibility test is a conjunction. Either the message or the
  condition is wrong; deciding which *is* the item.
- **2.9** cosmetic: `es_from_paired_f()` reports an error naming `es_from_paired_t()`.
  `.validate_pre_post_to_smd()` already takes `func_name`; the delegating routes just need
  to pass the callee's.

### Phase 5 — new work, scoped separately

5.1 simulate the `pool_sd = TRUE` kernel (study 08b) · 5.2 ten method arguments with zero
simulation coverage, the psychometric family entirely · 5.3 `convert_df()` is never called
in `simulations/`; Monte Carlo reaches 12 of 89 exported `es_from_*` · 5.4 should
VanderWeele's conversions become package options.

### Two ☐ that are actually done

`1.4` (line 472) and `2.2` (line 532) show unticked — they are **stale duplicate headings
kept for the record**. The ☑ entry above each is the live one. Check before acting.

---

## 6. One-minute restart on the other laptop

After a folder copy the branch is already checked out; these four lines confirm the copy
arrived intact and put the next item in front of you.

```bash
git status -sb                  # expect "## audit-remediation", clean. If it errors,
                                # .git did not come across -- copy the FOLDER, not its contents
git log --oneline -3            # 22e7fb0 / 7e385fb / 4732667 at the time of writing
cd simulations && Rscript tests/run_tests.R     # expect 548 / 0 / 0 / 0
sed -n '/^### 4.4/,/^### 4.5/p' ROADMAP.md      # the next item
```

If the suite errors on a missing package, install from the list in §0(d) — the R library
does not travel with the folder.

Then re-read §4 above and implement — the measurement is done.
