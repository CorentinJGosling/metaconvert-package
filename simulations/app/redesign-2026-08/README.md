# Redesign proposal — simulation results viewer (2026-08-25)

Three static mockups of a redesigned `simulations/app/app.R`, answering three
complaints about the current design:

1. the headline cards are unreadable,
2. the sidebar mixes scientific choices with plotting controls,
3. nobody understands the *Estimand check* tab.

Open **`variants.html`** in a browser. It carries all three directions behind a
switcher, plus the brief and a side-by-side comparison.

## Directions

| | | |
|---|---|---|
| **A** | Verdict first | No rail. The question is a sentence across the top, the plotting controls sit beside the figure. Signature: the **boundary grid** — 20 measured coverage cells tinted by shortfall, so *where the advice stops holding* is readable at a glance. |
| **B** | Two questions | The rail survives but only draws. The reference quantity becomes a tab strip inside the document. Signature: the **ledger** — routes × reference quantities, whose mirror is the whole result. |
| **C** | Results listing | No rail, no tabs, no cards. The control state is one query line; findings hang off numbered marks in a gutter. |

## Every number is measured

All figures come from `data/aggregated/01a_smd_to_cor_r_nrep1000.csv`, target
`biserial_population`, `p_exp = 0.30`, 20 conditions, 1,000 replications each —
the view the app opens on. Nothing was invented to fill a layout.

Writing them surfaced a content problem in the current app, and all three
directions are built so it cannot recur: **the headline names a leader that the
`Scored against` dropdown decides.** Mean coverage is 0.940 / 0.764 against the
biserial parameter and 0.803 / 0.940 against the point-biserial one. Scored
against its own estimand each route is near-exact (mean |bias| 0.004 and 0.005).
The two routes estimate different correlations; which one is wanted is a property
of the review, not of the package.

## Build

```sh
Rscript app/redesign-2026-08/make-figures.R    # run from simulations/
sh     app/redesign-2026-08/build.sh
```

`build.sh` concatenates the `_p*.html` parts with the SVG figures inlined.
Figure 2 has its glyph ids rewritten (`glyph-` → `gB-`) because both figures
carry `<use xlink:href="#glyph-N-M">` references and, inlined into one document,
the second figure would otherwise resolve every glyph against the first.
