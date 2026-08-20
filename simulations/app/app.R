## =============================================================================
## metaConvert simulation results
##
## Browses the ADEMP aggregates written by run_study() into
## simulations/data/aggregated/*.csv. Entirely generic: every study is discovered
## from disk and its condition columns, methods, targets and replication counts
## are read from the file, so a new study appears here as soon as it has been run.
##
## Visual identity follows metaconvert.org: rose accent, indigo for data, warm
## neutrals, Work Sans + Cascadia Code. The brand contains no green.
##
## Run from the simulations/ root:   shiny::runApp("app")
## =============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(DT)

## ---- brand -------------------------------------------------------------------
MC <- list(
  paper  = "#f7f7f3",
  white  = "#ffffff",
  ink    = "#211f38",
  ink2   = "#2c2a47",
  soft   = "#6b6a7d",
  faint  = "#8f8ea0",
  rule   = "#e6e6e0",
  rose   = "#d66268",
  rose_d = "#a23c42",
  amber  = "#8a6712",
  red    = "#9c2b2b"
)

## Categorical scale for methods. Anchored on the brand's rose and indigo and
## extended with amber, plum, blue and brown; deliberately green-free.
##
## Ordered so that neighbouring positions are far apart in hue. Methods are taken
## alphabetically, so the two rose tones must not land on adjacent slots -- with
## the deep rose in position 4, study 09's "2x2_tetrachoric" and "lipsey_cooper"
## came out as two near-identical reds.
MC_METHODS <- c("#d66268", "#2c2a47", "#8a6712", "#6b4370",
                "#4a6fa5", "#8a5a3c", "#a23c42", "#c98a5b",
                "#3f4a6b", "#b07d4f", "#7d5e8c", "#c4525f")

## ---- locate the aggregates ---------------------------------------------------
.find_agg_dir <- function() {
  cands <- c(file.path("..", "data", "aggregated"),
             file.path("data", "aggregated"),
             file.path(dirname(getwd()), "data", "aggregated"))
  for (p in cands) if (dir.exists(p)) return(normalizePath(p, winslash = "/"))
  stop("could not locate simulations/data/aggregated")
}
AGG_DIR <- .find_agg_dir()

## Human-readable titles. A study missing from this list still works; it simply
## shows its file stem.
STUDY_LABELS <- c(
  "01a_smd_to_cor_r"      = "SMD to correlation (r)",
  "01b_smd_to_cor_z"      = "SMD to correlation (z)",
  "02a_cor_to_smd_GROUPS" = "Correlation to SMD (grouped data)",
  "02b_cor_to_smd_CONT"   = "Correlation to SMD (continuous data)",
  "03a_2x2_to_cor_CAT"    = "2x2 table to correlation (categorical latent)",
  "03b_2x2_to_cor_CONT"   = "2x2 table to correlation (continuous latent)",
  "04_or_to_rr"           = "Odds ratio to risk ratio",
  "05_rr_to_or"           = "Risk ratio to odds ratio",
  "06_or_se_imputation"   = "Odds ratio standard error imputation",
  "07a_ancova_to_smd_d"   = "ANCOVA to SMD (Cohen's d)",
  "07b_ancova_to_smd_g"   = "ANCOVA to SMD (Hedges' g)",
  "08a_pre_post_to_smd_d" = "Pre-post to SMD (Cohen's d)",
  "08b_pre_post_to_smd_g" = "Pre-post to SMD (Hedges' g)",
  "09a_or_to_cor_CONT"    = "Odds ratio to correlation (continuous latent)",
  "09b_or_to_cor_CAT"     = "Odds ratio to correlation (categorical latent)"
)

## Which target this study should OPEN on, where the generic ordering below gets it
## wrong. A study missing from this list falls back to the ordering rule in
## target_ui(); this map only overrides, it never has to be complete.
##
## WHY IT HAS TO BE DECLARED RATHER THAN DERIVED (roadmap 4.3). The generic rule
## prefers a shared benchmark over `own`, which is right almost everywhere: the gap to
## the population parameter is usually the FINDING (study 07's covariate
## misspecification, study 08's standardiser choice, study 03a's tetrachoric-under-a-
## categorical-mechanism). But the same shape also occurs when the benchmark is simply
## on the wrong scale for the routes being scored, and no statistic in the file
## separates the two: both produce a bias that is large and flat in n. The judgement
## lives in the study file's header prose, so it is written down here.
##
## 01b IS THE CASE THAT NEEDS IT. Its targets were named symmetrically with 01a's, but
## the two scales are not symmetric: on the r scale `biserial_population` IS what
## viechtbauer estimates, whereas on the z scale viechtbauer returns a
## VARIANCE-STABILISING transform, not atanh() of anything. So `fisherz_biserial` --
## which the ordering rule picks, there being no target named "population" -- is a
## quantity NEITHER route estimates. At rho = 0.75, p_exp = 0.5 the three values are
## atanh(biserial) 0.973, viechtbauer's transform 0.724, atanh(point-biserial) 0.691.
## Measured consequence on the shipped aggregate, mean over the grid:
##
##   target                  lipsey_cooper          viechtbauer
##   fisherz_biserial        bias -0.108, cov 0.736 bias -0.092, cov 0.773
##   fisherz_pointbiserial   bias -0.000, cov 0.944 bias +0.015, cov 0.934
##
## i.e. the app opened on the one benchmark that makes both routes look broken while
## the sidebar told the reader to prefer it over `own`. `fisherz_pointbiserial` is a
## genuine shared benchmark -- lipsey_cooper estimates it exactly, and viechtbauer's
## gap to it is a real, interpretable estimand difference -- so that is the opening
## view. `fisherz_biserial` stays available and is labelled for what it is.
##
## STUDIES 03 AND 09 WERE CHECKED AND ARE NOT LISTED, deliberately. They look like the
## same defect and are not: their `population` target is attained exactly by the route
## that is correct for the mechanism (03a phi (r) 0.0035/0.962 = its own; 03b and 09a
## tetrachoric (r) likewise), and where no route attains it -- 09b, where phi is the
## estimand and metaConvert ships no phi route -- that gap IS the study's result. Both
## also carry a scale switch that opens on (r), which is the scale their shared targets
## are on. Changing their default would have deleted the finding.
STUDY_RANK_ON <- c(
  "01b_smd_to_cor_z" = "fisherz_pointbiserial"
)

## The scale the SHARED targets are expressed on, for the studies whose method names
## carry an "(r)" / "(z)" suffix. Their population/sample columns stay on the r scale
## whichever route is being scored -- studies/09_or_to_cor.R says so in terms
## ("`population`/`sample` stay on the r scale and are interpretable only for the (r)
## routes") -- so selecting the (z) methods against one of them compares two different
## transforms. The app opens on (r), so this only fires after the reader switches.
STUDY_SHARED_SCALE <- c(
  "03a_2x2_to_cor_CAT" = "r", "03b_2x2_to_cor_CONT" = "r",
  "09a_or_to_cor_CONT" = "r", "09b_or_to_cor_CAT"   = "r"
)

## The package argument each study evaluates. Shown beside the title so the app
## answers "which setting am I looking at?" without opening the study file.
STUDY_ARG <- c(
  "01a_smd_to_cor_r" = "smd_to_cor", "01b_smd_to_cor_z" = "smd_to_cor",
  "02a_cor_to_smd_GROUPS" = "cor_to_smd", "02b_cor_to_smd_CONT" = "cor_to_smd",
  "03a_2x2_to_cor_CAT" = "table_2x2_to_cor", "03b_2x2_to_cor_CONT" = "table_2x2_to_cor",
  "04_or_to_rr" = "or_to_rr", "05_rr_to_or" = "rr_to_or",
  "06_or_se_imputation" = "es_from_or", "07a_ancova_to_smd_d" = "cov_outcome_r",
  "07b_ancova_to_smd_g" = "cov_outcome_r", "08a_pre_post_to_smd_d" = "pre_post_to_smd",
  "08b_pre_post_to_smd_g" = "pre_post_to_smd",
  "09a_or_to_cor_CONT" = "or_to_cor", "09b_or_to_cor_CAT" = "or_to_cor"
)

## One line saying what the study is for. A reader who has never opened the
## programme should be able to tell, from the header alone, which decision the
## numbers below inform.
STUDY_Q <- c(
  "01a_smd_to_cor_r"      = "A study reports a standardised mean difference and the review pools correlations. Which conversion route should convert_df() take, and does the answer change on the r scale?",
  "01b_smd_to_cor_z"      = "The same conversion reported on the Fisher z scale, where the variance-stabilising transformation changes what a fixed interval width means.",
  "02a_cor_to_smd_GROUPS" = "A correlation has to become an SMD, and the underlying data were two groups. Which route recovers the grouped estimand?",
  "02b_cor_to_smd_CONT"   = "The same conversion when the predictor was genuinely continuous, so no grouping ever existed to recover.",
  "03a_2x2_to_cor_CAT"    = "A 2x2 table becomes a correlation when the latent variable really is categorical, so the tetrachoric assumption is violated by construction.",
  "03b_2x2_to_cor_CONT"   = "The same conversion when the dichotomy hides a genuinely continuous latent variable, which is the case tetrachoric was built for.",
  "04_or_to_rr"           = "Turning an odds ratio into a risk ratio requires a baseline risk that the report rarely gives. Which route survives the guess, and where does it fail outright?",
  "05_rr_to_or"           = "The reverse conversion, and the region where the reconstruction becomes non-estimable rather than merely inaccurate.",
  "06_or_se_imputation"   = "An odds ratio arrives without a standard error. How well does the shipped imputation reproduce the one the study did not report?",
  "07a_ancova_to_smd_d"   = "ANCOVA means become an SMD only through a covariate-outcome correlation the analyst has to guess. What does guessing wrong cost?",
  "07b_ancova_to_smd_g"   = "The same misspecification surface, carried through to the small-sample-corrected Hedges g.",
  "08a_pre_post_to_smd_d" = "Five ways to standardise a pre-post design, which do not estimate the same quantity. What is the price of picking one?",
  "08b_pre_post_to_smd_g" = "The same five routes after the small-sample correction, where the standardiser and the correction interact.",
  "09a_or_to_cor_CONT"    = "An odds ratio becomes a correlation, with a genuinely continuous latent variable behind the dichotomy.",
  "09b_or_to_cor_CAT"     = "The same conversion when the latent variable is categorical, so every route is estimating something it was not designed for."
)

## ---- performance measures ----------------------------------------------------
METRICS <- list(
  bias        = list(label = "Bias",                        mcse = "bias_mcse",     ref = 0,    better = "closer to 0"),
  rmse        = list(label = "Root mean squared error",     mcse = NA,              ref = 0,    better = "smaller"),
  coverage    = list(label = "Confidence interval coverage", mcse = "coverage_mcse", ref = 0.95, better = "closer to 0.95"),
  se_ratio    = list(label = "SE ratio (model / empirical)", mcse = "se_ratio_mcse", ref = 1,    better = "closer to 1"),
  emp_se      = list(label = "Empirical standard error",    mcse = "emp_se_mcse",   ref = NA,   better = "smaller"),
  ci_width    = list(label = "Mean interval width",         mcse = NA,              ref = NA,   better = "smaller"),
  nonest_rate = list(label = "Non-estimable rate",          mcse = NA,              ref = 0,    better = "smaller")
)

## Coverage and the SE ratio are calibration measures: a nominal 95% interval is
## built to cover a FIXED quantity. Against a *_sample target -- the statistic
## recomputed on that replication's own sample -- they are not interpretable, and
## the app says so rather than printing a number that cannot be read.
##
## `own` is NOT such a case, and must not be lumped in with it. Each method's own
## estimand is a fixed population quantity; coverage against it is the meaningful
## diagnostic "does this interval cover the value this method is estimating?".
## What is wrong with `own` is that it is a DIFFERENT quantity per method, so it
## cannot rank methods against each other -- a separate warning, raised
## separately.
## Is the selected shared target on a different scale from the selected methods?
##
## Studies whose method names carry an "(r)" / "(z)" suffix keep their population and
## sample columns on ONE scale whichever route is scored -- studies/09_or_to_cor.R:
## "`population`/`sample` stay on the r scale and are interpretable only for the (r)
## routes". Selecting the (z) methods against one of them compares two different
## transforms, and the resulting gap is arithmetic, not performance. The app opens on
## (r), so this only fires once the reader has switched the scale.
.scale_mismatch <- function(study, scale, target) {
  if (is.null(study) || is.null(scale) || is.null(target)) return(FALSE)
  if (!study %in% names(STUDY_SHARED_SCALE)) return(FALSE)
  if (grepl("^own", target)) return(FALSE)          # own is scale-matched by construction
  !identical(scale, unname(STUDY_SHARED_SCALE[[study]]))
}

.calibration_ok <- function(metric, target) {
  if (is.null(metric) || is.null(target)) return(TRUE)
  if (!metric %in% c("coverage", "se_ratio")) return(TRUE)
  !grepl("sample", target)
}

## ---- discovery ---------------------------------------------------------------
scan_files <- function() {
  f <- list.files(AGG_DIR, pattern = "_nrep[0-9]+[.]csv$", full.names = TRUE)
  if (!length(f)) return(NULL)
  stem <- sub("_nrep[0-9]+[.]csv$", "", basename(f))
  nrep <- as.integer(sub(".*_nrep([0-9]+)[.]csv$", "\\1", basename(f)))
  d <- data.frame(path = f, stem = stem, nrep = nrep, stringsAsFactors = FALSE)
  d <- d[order(d$stem, -d$nrep), ]
  d$label <- ifelse(d$stem %in% names(STUDY_LABELS), STUDY_LABELS[d$stem], d$stem)
  ## the leading number of the study file, used as the index chip
  d$idx <- toupper(sub("^([0-9]+[a-z]?)_.*$", "\\1", d$stem))
  d
}
FILES <- scan_files()

condition_cols <- function(df) {
  i <- match("method", names(df))
  if (is.na(i) || i < 2) return(character(0))
  names(df)[seq_len(i - 1)]
}

fmt_n <- function(x) format(as.integer(x), big.mark = ",")
fmt_v <- function(x, d = 3) if (!is.finite(x)) "--" else formatC(x, format = "f", digits = d)

## Target names are column values, so they arrive snake_cased. Prose reads them
## back as English; the places that refer to the value in the file keep the raw
## form inside a <code> span.
pretty_target <- function(z) {
  if (is.null(z) || !length(z)) return("")
  if (z == "own") "own estimand"
  else if (z == "own_true_r") "own estimand at the true r"
  else gsub("_", " ", z)
}

## Order the targets for the radio list, shared benchmarks first.
##
## Targets are two different kinds of thing, and presenting them as one flat list
## invites the mistake of ranking methods against `own`. A SHARED target is one
## benchmark applied to every method, and is what you need to choose between them;
## `own` is a DIFFERENT benchmark per method, a diagnostic rather than a basis for
## comparison. Shared first, population-scale ones ahead of the rest, `own` last.
target_order <- function(targets) {
  own_like <- grep("^own", targets, value = TRUE)
  shared <- setdiff(targets, own_like)
  pop <- grep("population", shared, value = TRUE)
  c(sort(pop), sort(setdiff(shared, pop)), sort(own_like))
}

## Which target the app opens on. A declared entry in STUDY_RANK_ON wins, but only if
## the study actually records it -- an entry naming a target that is not in the file
## would silently do nothing, which is how a hardcoded list rots. Otherwise: the first
## population-scale target, else the first in the ordering above.
default_target <- function(targets, study = NULL) {
  if (!length(targets)) return(NULL)
  if (!is.null(study) && study %in% names(STUDY_RANK_ON)) {
    want <- STUDY_RANK_ON[[study]]
    if (want %in% targets) return(want)
  }
  pop <- grep("population", setdiff(targets, grep("^own", targets, value = TRUE)),
              value = TRUE)
  if (length(pop)) sort(pop)[1] else target_order(targets)[1]
}

## Shared targets that NO method in this study estimates.
##
## Exact and threshold-free: run_study() writes one row per (condition, method,
## target), so a method whose `own` estimand IS a named shared target produces the
## identical bias column against both. Comparing the two ordered vectors therefore
## answers "does any route actually target this?" with no tolerance to choose.
##
## The answer is worth showing either way, and the label is deliberately neutral. In
## 01b it says `fisherz_biserial` is an artefact of the target naming; in 08 and 09b it
## says no route recovers the quantity the review wants, which is those studies' result.
##
## SAME-SAMPLE TARGETS ARE EXCLUDED, and must stay excluded. `own` is a population
## quantity in every study, so a `*sample*` target -- the same statistic recomputed on
## the replication's own draw -- can never equal it, and the test would fire on all 12
## studies while saying nothing: a method DOES estimate that statistic, it is simply
## random rather than fixed. The app already carries a dedicated warning for them.
unattained_targets <- function(agg) {
  tg <- unique(agg$target)
  if (!"own" %in% tg) return(character(0))
  shared <- setdiff(tg, grep("^own", tg, value = TRUE))
  shared <- shared[!grepl("sample", shared)]
  if (!length(shared) || !"bias" %in% names(agg)) return(character(0))
  k <- condition_cols(agg)
  if (!length(k)) return(character(0))
  vec <- function(m, t) {
    d <- agg[agg$method == m & agg$target == t, , drop = FALSE]
    d[do.call(order, d[k]), "bias"]
  }
  ms <- unique(agg$method)
  shared[!vapply(shared, function(t)
    any(vapply(ms, function(m) isTRUE(all.equal(vec(m, t), vec(m, "own"))), logical(1))),
    logical(1))]
}

## ---- plot theme --------------------------------------------------------------
## The y-axis title is deliberately absent: the measure is named in the key row
## directly above every plot, and dropping the rotated title returns ~30px of
## width to the panels, which is where the information is.
theme_mc <- function() {
  theme_minimal(base_size = 13, base_family = "Work Sans") +
    theme(
      text             = element_text(colour = MC$ink),
      axis.title.x     = element_text(size = 11.5, colour = MC$soft, margin = margin(t = 9)),
      axis.title.y     = element_blank(),
      axis.text        = element_text(size = 10, colour = MC$soft, family = "Cascadia Code"),
      panel.grid.major = element_line(colour = MC$rule, linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.spacing.x  = unit(1.5, "lines"),
      panel.spacing.y  = unit(1.7, "lines"),
      ## a boxed grey strip is the ggplot default, not a design decision. A bold
      ## left-aligned label sits in the panel's own reading order instead.
      strip.background = element_blank(),
      strip.text       = element_text(size = 11, colour = MC$ink, face = "bold",
                                      hjust = 0, margin = margin(0, 0, 6, 0)),
      ## the legend lives in the HTML key above the plot, where it is shared with
      ## the sidebar checkboxes and the ranking table
      legend.position  = "none",
      plot.margin      = margin(4, 12, 2, 2)
    )
}

## ---- UI ----------------------------------------------------------------------
## The stylesheet is a PLAIN string, never a sprintf() format. CSS is full of
## literal per-cent signs (width:100%, background-size:98% 62%, @media rules),
## and inside a format string each one is read as a conversion specification and
## aborts the app at start-up with "too few arguments". The custom properties are
## built from the MC palette by pasting, so the two still cannot drift.
##
## It is also a SINGLE-QUOTED R string, so no apostrophe may appear anywhere in
## it -- including inside a CSS comment -- or the string closes early and the
## file no longer parses.
CSS_VARS <- list(paper = MC$paper, ink = MC$ink, ink2 = MC$ink2, soft = MC$soft,
                 faint = MC$faint, rule = MC$rule, rose = MC$rose,
                 "rose-d" = MC$rose_d, amber = MC$amber, red = MC$red)

app_css <- paste0(
":root{",
paste0("--", names(CSS_VARS), ":", unlist(CSS_VARS), ";", collapse = ""),
'--rule2:#f0f0ea;--rose-l:#fdf3f3;--amber-l:#fbf6e9;--red-l:#fdf0f0;
--mono:"Cascadia Code","Cascadia Mono",ui-monospace,SFMono-Regular,Menlo,monospace;}

body{background:var(--paper);color:var(--ink);font-size:14px;}
::selection{background:rgba(214,98,104,.18);}

/* =========================================================== masthead ===== */
.navbar{
  background:#fff !important;border-bottom:1px solid var(--rule);
  box-shadow:none;min-height:56px;padding:0 22px;
}
.navbar>.container-fluid{padding:0;}
.navbar .navbar-brand{flex:1 1 auto;width:100%;margin:0;padding:0;}
.mast{display:flex;align-items:center;gap:13px;width:100%;min-height:56px;}
.mast .wm{font-weight:700;font-size:17px;letter-spacing:-.015em;color:var(--ink);}
.mast .bar{width:1px;height:17px;background:var(--rule);}
.mast .sub{font-size:13px;color:var(--soft);font-weight:400;}
.mast .right{margin-left:auto;display:flex;align-items:center;gap:14px;}
.mast .prov{font-family:var(--mono);font-size:11.5px;color:var(--faint);letter-spacing:.01em;}
.mast .prov b{color:var(--soft);font-weight:600;}
.helpbtn{
  border:1px solid var(--rule);background:#fff;border-radius:8px;
  padding:6px 12px;font-size:12.5px;font-weight:600;color:var(--ink);
  cursor:pointer;line-height:1.2;
}
.helpbtn:hover{border-color:var(--rose);color:var(--rose-d);background:var(--rose-l);}
.helpbtn:focus-visible{outline:none;box-shadow:0 0 0 3px rgba(214,98,104,.2);}

/* ========================================================== sidebar ======= */
.bslib-sidebar-layout>.sidebar{background:#fff;border-right:1px solid var(--rule);}
.bslib-sidebar-layout>.sidebar>.sidebar-content{padding:16px 20px 48px;gap:0;}
.bslib-sidebar-layout>.main{padding:22px 26px 8px;}

.sec{display:flex;align-items:center;gap:9px;margin:24px 0 11px;}
/* the sidebar collapse chevron floats over the top-right corner, so the first
   rule has to stop short of it */
.sec:first-child{margin-top:2px;padding-right:26px;}
.sec .t{
  font-size:10.5px;font-weight:700;letter-spacing:.11em;text-transform:uppercase;
  color:var(--soft);white-space:nowrap;
}
.sec .r{flex:1 1 auto;height:1px;background:var(--rule);}

.sidebar label,.sidebar .form-label,.sidebar .control-label{
  font-size:12.5px;font-weight:600;color:var(--ink);margin-bottom:5px;
}
.form-select,.form-control{
  border-color:var(--rule);border-radius:8px;font-size:13.5px;color:var(--ink);
  padding:7px 11px;background-color:#fff;
}
.form-select:focus,.form-control:focus{
  border-color:var(--rose);box-shadow:0 0 0 3px rgba(214,98,104,.16);
}
.shiny-input-container{margin-bottom:11px;width:100%;}
.shiny-input-container:last-child{margin-bottom:0;}

/* the study picker is the primary axis of the whole app, so it is the one
   control that gets weight */
.pick .form-select{font-weight:600;font-size:14px;padding:9px 11px;}
.nrep-chip{
  display:inline-flex;align-items:baseline;gap:7px;font-family:var(--mono);
  font-size:12.5px;color:var(--ink2);background:var(--paper);
  border:1px solid var(--rule);border-radius:7px;padding:6px 10px;
}
.nrep-chip span{font-size:10.5px;letter-spacing:.07em;text-transform:uppercase;color:var(--faint);}

/* radio + checkbox groups, given room to be clicked (44px targets are not
   achievable in a dense control rail, but 34px with generous hit padding is) */
.opts .form-check{
  padding:7px 10px 7px 31px;border-radius:8px;margin:0 0 1px;min-height:34px;
  display:flex;align-items:center;
}
.opts .form-check:hover{background:var(--paper);}
.opts .form-check-input{margin-left:-23px;margin-top:0;float:none;flex:none;}
.opts .form-check-label{font-size:13px;font-weight:500;color:var(--ink);cursor:pointer;}
.form-check-input{border-color:#c9c9c2;cursor:pointer;}
.form-check-input:checked{background-color:var(--rose);border-color:var(--rose);}
.form-check-input:focus{box-shadow:0 0 0 3px rgba(214,98,104,.18);border-color:var(--rose);}

/* method list doubles as the persistent legend */
.opts.methods .form-check{padding:5px 8px 5px 29px;min-height:30px;}
.mk{display:inline-flex;align-items:center;gap:8px;}
.mk i{width:9px;height:9px;border-radius:2px;flex:none;}
.mk span{font-family:var(--mono);font-size:12.5px;color:var(--ink2);}
.qa{display:flex;align-items:center;gap:7px;margin:0 0 6px 2px;}
.qa a{
  font-size:11.5px;font-weight:600;color:var(--rose-d);cursor:pointer;
  text-decoration:none;border-bottom:1px dashed rgba(162,60,66,.35);
}
.qa a:hover{color:var(--rose);}
.qa em{color:var(--rule);font-style:normal;}

.hint{font-size:12px;line-height:1.5;color:var(--soft);margin:2px 0 4px;}
.hint b{color:var(--ink2);font-weight:600;}
.hint a{color:var(--rose-d);cursor:pointer;font-weight:600;text-decoration:none;
  border-bottom:1px dashed rgba(162,60,66,.35);}

/* filters subset the data rather than mapping it, so they read differently */
.held{background:var(--paper);border:1px solid var(--rule);border-radius:10px;padding:11px 12px 3px;}
.held .form-select{background-color:#fff;font-family:var(--mono);font-size:12.5px;}
.held label{font-family:var(--mono);font-size:11.5px;font-weight:600;color:var(--soft);}

/* ==================================================== study header ======== */
.hdr{display:flex;align-items:flex-start;gap:14px;flex-wrap:wrap;margin:0 0 18px;}
.hdr .txt{flex:1 1 380px;min-width:0;}
.hdr .line{display:flex;align-items:center;gap:11px;flex-wrap:wrap;}
.ix{
  font-family:var(--mono);font-size:11.5px;font-weight:700;color:var(--rose-d);
  background:var(--rose-l);border:1px solid rgba(214,98,104,.24);border-radius:6px;
  padding:4px 8px;letter-spacing:.05em;flex:none;
}
.hdr h1{
  font-size:26px;font-weight:700;letter-spacing:-.022em;margin:0;
  color:var(--ink);line-height:1.14;
}
.arg{
  font-family:var(--mono);font-size:12.5px;color:var(--ink2);background:var(--paper);
  border:1px solid var(--rule);border-radius:6px;padding:3px 8px;
}
.hdr .q{font-size:13.5px;color:var(--soft);margin:8px 0 0;max-width:76ch;line-height:1.55;}

.facts{display:flex;flex:none;}
.facts .f{padding:2px 0 2px 18px;margin-left:18px;border-left:1px solid var(--rule);}
.facts .f:first-child{border-left:0;margin-left:0;padding-left:0;}
.facts .v{font-family:var(--mono);font-size:17px;font-weight:600;color:var(--ink);line-height:1.1;}
.facts .k{font-size:10px;text-transform:uppercase;letter-spacing:.09em;color:var(--faint);margin-top:4px;}

/* ========================================================= headline ======= */
.kpis{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin:0 0 18px;}
.kpi{
  background:#fff;border:1px solid var(--rule);border-radius:12px;
  padding:14px 16px 15px;position:relative;overflow:hidden;
}
.kpi.lead{border-color:rgba(214,98,104,.32);}
.kpi.lead:before{content:"";position:absolute;left:0;top:0;bottom:0;width:3px;background:var(--rose);}
.kpi .k{font-size:10px;font-weight:700;letter-spacing:.1em;text-transform:uppercase;color:var(--faint);}
.kpi .v{
  font-family:var(--mono);font-size:23px;font-weight:600;color:var(--ink);
  margin-top:7px;line-height:1.15;letter-spacing:-.01em;word-break:break-word;
}
.kpi.lead .v{color:var(--rose-d);}
.kpi .v u{font-size:13px;color:var(--soft);font-weight:400;text-decoration:none;}
.kpi .n{font-size:12.5px;color:var(--soft);margin-top:6px;line-height:1.45;}
.kpi .n b{color:var(--ink2);font-weight:600;font-family:var(--mono);}
.chip{
  display:inline-flex;align-items:center;gap:5px;font-size:11px;font-weight:700;
  border-radius:20px;padding:3px 9px;border:1px solid;letter-spacing:.02em;
  text-transform:uppercase;margin-top:8px;
}
.chip.solid{color:var(--ink2);background:rgba(44,42,71,.055);border-color:rgba(44,42,71,.16);}
.chip.warn{color:var(--amber);background:var(--amber-l);border-color:rgba(138,103,18,.22);}
.chip.bad{color:var(--red);background:var(--red-l);border-color:rgba(156,43,43,.2);}

/* ============================================================= tabs ======= */
.card{
  border:1px solid var(--rule);border-radius:14px;background:#fff;
  box-shadow:0 1px 2px rgba(33,31,56,.04);
}
.card>.card-header{background:#fff;border-bottom:1px solid var(--rule);padding:0 8px;}
.nav-tabs{border-bottom:0;gap:2px;}
.nav-tabs .nav-link{
  border:0;border-radius:0;color:var(--soft);font-weight:600;font-size:13.5px;
  padding:13px 14px;border-bottom:2px solid transparent;
}
.nav-tabs .nav-link:hover{color:var(--ink);background:transparent;}
.nav-tabs .nav-link.active{
  color:var(--ink) !important;background:transparent !important;
  border-bottom:2px solid var(--rose) !important;
}
.nav-tabs .nav-link:focus-visible{outline:none;box-shadow:inset 0 0 0 2px rgba(214,98,104,.3);}
.card-body{padding:20px 22px 22px;}

/* ============================================== key row above the plot ==== */
.key{
  display:flex;flex-wrap:wrap;gap:8px 18px;align-items:center;
  padding:0 0 13px;margin:0 0 16px;border-bottom:1px solid var(--rule);
}
.key .lbl{
  font-size:10px;font-weight:700;letter-spacing:.1em;text-transform:uppercase;
  color:var(--faint);
}
.key .m{display:inline-flex;align-items:center;gap:7px;font-family:var(--mono);
  font-size:12.5px;color:var(--ink2);}
.key .m i{width:10px;height:10px;border-radius:3px;flex:none;}
.key .band{margin-left:auto;display:inline-flex;align-items:center;gap:8px;
  font-size:12px;color:var(--soft);cursor:help;}
.key .band i{
  width:26px;height:12px;border-radius:3px;flex:none;
  background:rgba(44,42,71,.09);
  border-top:1px dashed rgba(44,42,71,.42);border-bottom:1px dashed rgba(44,42,71,.42);
}
.measline{font-size:13px;color:var(--soft);margin:-6px 0 14px;}
.measline b{color:var(--ink);font-weight:600;}

/* ============================================================ prose ======= */
.lede{font-size:13.5px;line-height:1.62;color:var(--soft);max-width:78ch;margin:0 0 16px;}
.lede b{color:var(--ink);font-weight:600;}
.lede i{color:var(--ink2);}
.lede code{font-family:var(--mono);font-size:12.5px;color:var(--ink2);
  background:var(--paper);padding:1px 5px;border-radius:4px;}
.note{
  display:flex;gap:11px;border-radius:10px;padding:12px 14px;font-size:13px;
  line-height:1.55;margin:0 0 16px;border:1px solid;max-width:82ch;
}
.note .ic{font-family:var(--mono);font-weight:700;flex:none;line-height:1.4;}
.note.warn{background:var(--amber-l);border-color:rgba(138,103,18,.22);color:#6d520e;}
.note.warn b{color:#54400a;font-weight:700;}
.note.info{background:var(--paper);border-color:var(--rule);color:var(--soft);}
.note.info b{color:var(--ink);font-weight:600;}

/* =========================================================== tables ======= */
table.dataTable{font-size:13px;}
table.dataTable thead th{
  font-size:10.5px;text-transform:uppercase;letter-spacing:.06em;font-weight:700;
  color:var(--faint);border-bottom:1px solid var(--rule) !important;padding:11px 10px;
}
table.dataTable tbody td{
  font-family:var(--mono);color:var(--ink2);padding:10px;
  border-top:1px solid var(--rule2);
}
table.dataTable tbody tr:hover td{background:var(--paper);}
/* a wrapped method name doubles the row height, and the abs-bias colour bar is
   drawn to the cell box, so it turned into a tall pink slab on those rows */
table.dataTable td.mth{font-family:var(--mono);font-weight:600;color:var(--ink);white-space:nowrap;}
table.dataTable td.rk{font-family:var(--mono);color:var(--faint);width:26px;}
.dataTables_wrapper{overflow-x:auto;}
.dataTables_wrapper .dataTables_filter input{
  border:1px solid var(--rule);border-radius:7px;padding:5px 9px;font-size:13px;
}
.dataTables_wrapper .dataTables_info,.dataTables_wrapper .dataTables_length,
.dataTables_wrapper .dataTables_paginate{font-size:12.5px;color:var(--soft);}

.btn{border-radius:8px;font-weight:600;font-size:13px;padding:7px 13px;}
.btn-outline-secondary{border-color:var(--rule);color:var(--ink);background:#fff;}
.btn-outline-secondary:hover{border-color:var(--rose);color:var(--rose-d);background:var(--rose-l);}
/* the card body is a flex column, so a plain width:auto anchor still stretches
   edge to edge; align-self is what actually shrink-wraps it */
#dl{display:inline-flex;align-items:center;gap:7px;width:auto;align-self:flex-start;
  margin:0 0 16px;}

/* =========================================================== footer ======= */
.foot{
  margin:20px 2px 10px;padding-top:14px;border-top:1px solid var(--rule);
  display:flex;flex-wrap:wrap;gap:6px 20px;align-items:baseline;
  font-size:11.5px;color:var(--faint);line-height:1.5;
}
.foot .mono{font-family:var(--mono);}
.foot b{color:var(--soft);font-weight:600;}

/* ============================================================ modal ======= */
.modal-content{border:1px solid var(--rule);border-radius:14px;}
.modal-header{border-bottom:1px solid var(--rule);padding:16px 22px;}
.modal-title{font-size:17px;font-weight:700;letter-spacing:-.01em;}
.modal-body{padding:18px 22px 22px;}
.modal-footer{border-top:1px solid var(--rule);padding:12px 22px;}
.doc h4{font-size:11px;font-weight:700;letter-spacing:.1em;text-transform:uppercase;
  color:var(--rose-d);margin:20px 0 8px;}
.doc h4:first-child{margin-top:0;}
.doc p{font-size:13.5px;line-height:1.62;color:var(--soft);margin:0 0 10px;max-width:74ch;}
.doc p b{color:var(--ink);font-weight:600;}
.doc code{font-family:var(--mono);font-size:12.5px;color:var(--ink2);
  background:var(--paper);padding:1px 5px;border-radius:4px;}
.doc ul{margin:0 0 10px;padding-left:18px;}
.doc li{font-size:13.5px;line-height:1.6;color:var(--soft);margin-bottom:5px;}
.doc li b{color:var(--ink);font-weight:600;}

/* ======================================================= responsive ======= */
@media (max-width:1250px){
  .kpis{grid-template-columns:1fr;}
  .facts{width:100%;}
  .facts .f:first-child{padding-left:0;}
}
@media (max-width:760px){
  .hdr h1{font-size:22px;}
  .facts{flex-wrap:wrap;gap:10px 0;}
  .mast .prov{display:none;}
  .bslib-sidebar-layout>.main{padding:16px 14px 8px;}
}
@media (prefers-reduced-motion:no-preference){
  .btn,.helpbtn,.nav-tabs .nav-link,.form-select,.form-control,.kpi,.opts .form-check{
    transition:background-color .16s ease,border-color .16s ease,color .16s ease,box-shadow .16s ease;
  }
}
')

## The full reading guide. It used to sit permanently in the sidebar as three
## paragraphs of small grey prose, where it cost ~200px of scroll on every visit
## and was read once. It is reference material, so it lives behind a button.
HOW_TO_READ <- HTML(
  '<div class="doc">
   <h4>What this is</h4>
   <p>A Monte Carlo evaluation of the conversion routes <code>convert_df()</code>
   chooses between. Each study fixes one package argument, generates data from a
   known truth, applies every competing formula to the same replication, and
   scores the results with the ADEMP measures of Morris, White &amp; Crowther
   (2019). Nothing here is a re-implementation: every number comes from calling
   the package itself.</p>

   <h4>Scored against</h4>
   <p>Targets are two different kinds of thing, and reading them as one flat list
   invites a real mistake.</p>
   <ul>
   <li><b>Population parameter</b> &mdash; the value the data were generated from.
   This is the quantity you actually want, so it is the default, and coverage and
   standard-error calibration are meaningful only here.</li>
   <li><b>Same-sample statistic</b> &mdash; the same statistic recomputed on that
   replication&rsquo;s own sample: how faithful is the conversion? Being random, it
   has no 0.95 reference, so coverage cannot be read against it.</li>
   <li><b>Each method&rsquo;s own estimand</b> &mdash; every method scored against a
   <i>different</i> quantity, the one it is designed to estimate. A diagnostic, not
   a ranking: it answers &ldquo;does this method compute its own value
   correctly?&rdquo;, never &ldquo;is that the value I want?&rdquo;. Two methods
   scored here are being held to two different standards.</li>
   </ul>
   <p>The gap between the last two is the estimand mismatch, and it is not bias
   &mdash; the method is computing its own quantity correctly, but that quantity is
   not the one being asked for. The <b>Estimand check</b> tab shows both at once,
   which is the only place the comparison is safe.</p>

   <p>A shared benchmark is sometimes marked <i>no route estimates this</i>. That is
   read off the data, not declared: a route whose own estimand IS a named target
   produces the identical bias column against both, so where no route does, nothing
   in the study is trying to produce that number. It means two different things and
   the tab you are on tells you which. In study 08 and in <b>OR to correlation
   (categorical latent)</b> it is the RESULT &mdash; the quantity a review wants and
   no conversion delivers. In <b>SMD to correlation (z)</b> it is an artefact of how
   the targets were named, which is why that study opens elsewhere.</p>

   <h4>Why coverage is the default measure</h4>
   <p>A conversion can be nearly unbiased in the point estimate and badly
   miscalibrated in its standard error &mdash; which matters more in meta-analysis
   than almost anywhere else, since the standard error sets the inverse-variance
   weight. <b>Coverage is the only ADEMP measure that degrades under either
   failure.</b> The two orderings genuinely disagree in these data, so the
   <b>Ranking</b> tab prints every measure and the selected one decides only the
   sort order.</p>

   <h4>The resolution band</h4>
   <p>Every panel carries a shaded band of &plusmn;2 &times; the median Monte Carlo
   standard error of its own cells. A difference that stays inside the band is
   smaller than the simulation&rsquo;s own noise and is <b>not a result</b>. The band
   is computed per panel, not once per plot: Monte Carlo error shrinks with sample
   size, so one global band would be far too narrow at small <i>n</i> and too wide
   at large <i>n</i>, misstating resolvability exactly where the question matters.
   Watching it contract across panels is the point.</p>

   <h4>Two deliberate defaults</h4>
   <ul>
   <li><b>A shared vertical scale.</b> Free scales let each panel pick its own
   range, so panels cannot be compared. There is a switch for the rare case where
   free scales are what you want.</li>
   <li><b>One reported scale at a time.</b> Several studies encode the output scale
   in the method name (<code>bonett (r)</code> and <code>bonett (z)</code> are one
   formula on two scales). Plotting both against one axis compares quantities that
   are not on the same scale, and the larger <i>z</i> values flatten the <i>r</i>
   values.</li>
   </ul>
   </div>')

ui <- page_sidebar(
  title = div(
    class = "mast",
    span(class = "wm", "metaConvert"),
    span(class = "bar"),
    span(class = "sub", "Simulation programme"),
    div(class = "right",
        span(class = "prov",
             if (is.null(FILES)) "no aggregates found"
             else HTML(paste0(fmt_n(length(unique(FILES$stem))), " studies &middot; <b>ADEMP</b>"))),
        tags$button(id = "howto", class = "helpbtn action-button", type = "button",
                    "How to read this"))
  ),
  window_title = "metaConvert simulation results",
  ## NOT fillable. In fill mode bslib stretches the tab card to the viewport and
  ## resizes plotOutput to the container, which (a) leaves a table-sized panel
  ## floating in 200px of white space and (b) silently discards the per-study
  ## plot height computed from the panel count. The page scrolls instead.
  fillable = FALSE,
  theme = bs_theme(
    version = 5,
    bg = MC$white, fg = MC$ink, primary = MC$rose,
    base_font = font_google("Work Sans"),
    code_font = font_collection("Cascadia Code", "Cascadia Mono", font_google("IBM Plex Mono")),
    "border-color" = MC$rule
  ),
  tags$head(tags$style(HTML(app_css))),

  sidebar = sidebar(
    width = 340, bg = MC$white,
    if (is.null(FILES)) {
      div(class = "note warn", span(class = "ic", "!"),
          div(HTML("No aggregate files found in <code>data/aggregated</code>. Run the
                    simulation programme first.")))
    } else {
      tagList(
        div(class = "sec", span(class = "t", "Study"), span(class = "r")),
        div(class = "pick",
            selectInput("study", NULL,
                        choices = stats::setNames(unique(FILES$stem),
                                                  sprintf("%s — %s",
                                                          FILES$idx[!duplicated(FILES$stem)],
                                                          FILES$label[!duplicated(FILES$stem)])))),
        uiOutput("nrep_ui"),

        div(class = "sec", span(class = "t", "Scored against"), span(class = "r")),
        uiOutput("target_ui"),

        div(class = "sec", span(class = "t", "Measure"), span(class = "r")),
        selectInput("metric", NULL,
                    choices = stats::setNames(names(METRICS),
                                              vapply(METRICS, `[[`, "", "label")),
                    selected = "coverage"),

        div(class = "sec", span(class = "t", "Panels"), span(class = "r")),
        uiOutput("xvar_ui"),
        uiOutput("facet_ui"),
        checkboxInput("free_y", "Free vertical scale per panel", FALSE),
        div(class = "hint", "A shared scale is the default so panels can be compared."),

        uiOutput("scale_ui"),

        div(class = "sec", span(class = "t", "Methods"), span(class = "r")),
        div(class = "qa",
            actionLink("all_m", "all"), tags$em("/"), actionLink("no_m", "none")),
        uiOutput("methods_ui"),

        uiOutput("filters_ui")
      )
    }
  ),

  uiOutput("study_header"),
  uiOutput("kpi_row"),

  ## Every card_body is fillable = FALSE. Left at the default, bslib treats
  ## DTOutput and plotOutput as fill items and stretches them to the card, so a
  ## two-row ranking table sat in 250px of white space and the data table grew an
  ## internal scrollbar of its own. Flow layout sizes each panel to its content.
  navset_card_tab(
    id = "tabs",
    nav_panel("Plot",
              card_body(fillable = FALSE,
                        uiOutput("key_row"), uiOutput("plot_slot"))),
    nav_panel("Ranking",
              card_body(fillable = FALSE,
                        uiOutput("rank_lede"), DTOutput("rank_tbl"))),
    nav_panel("Estimand check",
              card_body(fillable = FALSE,
                        uiOutput("estimand_lede"), uiOutput("estimand_slot"))),
    nav_panel("Data",
              card_body(fillable = FALSE,
                        uiOutput("data_lede"),
                        downloadButton("dl", "Download these rows",
                                       class = "btn-outline-secondary btn-sm"),
                        DTOutput("data_tbl")))
  ),

  uiOutput("footer")
)

## ---- server ------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$howto, {
    showModal(modalDialog(
      title = "How to read this", HOW_TO_READ, size = "l", easyClose = TRUE,
      footer = modalButton("Close")))
  })

  output$nrep_ui <- renderUI({
    req(input$study)
    n <- sort(unique(FILES$nrep[FILES$stem == input$study]), decreasing = TRUE)
    ## a select with one option is noise; show the fact instead
    if (length(n) == 1)
      tagList(div(class = "nrep-chip", span("reps"), fmt_n(n)),
              tags$div(style = "display:none;",
                       selectInput("nrep", NULL, choices = n, selected = n)))
    else
      selectInput("nrep", "Replications per condition", choices = n, selected = n[1])
  })

  raw <- reactive({
    req(input$study, input$nrep)
    p <- FILES$path[FILES$stem == input$study & FILES$nrep == as.integer(input$nrep)]
    req(length(p) == 1)
    utils::read.csv(p, stringsAsFactors = FALSE)
  })

  cond_cols <- reactive(condition_cols(raw()))

  ## Target names are study-specific: most studies record own/population/sample,
  ## but study 01 records biserial_* and pointbiserial_* separately and study 08
  ## adds own_true_r. Read them from the file rather than assuming.
  ## The ordering, the opening selection and the "no route estimates this" labelling
  ## live in target_order() / default_target() / unattained_targets() at the top of
  ## this file, so they can be exercised without starting Shiny (roadmap 4.3).
  output$target_ui <- renderUI({
    ag <- raw()
    tg <- unique(ag$target)
    own_like <- grep("^own", tg, value = TRUE)
    ord <- target_order(tg)
    unattained <- unattained_targets(ag)
    pretty <- vapply(ord, function(z) {
      base <- if (z == "own") "Each method's own estimand"
      else if (z == "own_true_r") "Own estimand, at the true r"
      else if (z == "population") "Population parameter"
      else if (z == "sample") "Same-sample statistic"
      else gsub("_", " ", z)
      if (z %in% unattained) paste0(base, " — no route estimates this") else base
    }, character(1))
    tagList(
      div(class = "opts",
          radioButtons("target", NULL,
                       choices = stats::setNames(ord, pretty),
                       selected = default_target(tg, input$study))),
      div(class = "hint",
          if (length(own_like))
            HTML(paste0("Shared benchmarks come first — use those to compare methods, ",
                        "except where one is marked <i>no route estimates this</i>. ",
                        "<b>Own estimand</b> scores each method against a different quantity. ",
                        "<a id=\"howto2\" class=\"action-button\">What that means</a>"))
          else NULL)
    )
  })

  observeEvent(input$howto2, {
    showModal(modalDialog(
      title = "How to read this", HOW_TO_READ, size = "l", easyClose = TRUE,
      footer = modalButton("Close")))
  })

  ## Several studies encode the output scale in the method name -- "bonett (r)"
  ## and "bonett (z)" are the same formula reported on the correlation and on the
  ## Fisher's z scale. Plotting both against one vertical axis compares numbers
  ## that are not on the same scale, and the z values, being larger, flatten the
  ## r values. Where such suffixes exist, offer them as a scale switch and show
  ## one scale at a time.
  method_scales <- reactive({
    m <- unique(raw()$method)
    sfx <- ifelse(grepl("\\s\\([^()]+\\)$", m), sub(".*\\s\\(([^()]+)\\)$", "\\1", m), NA)
    if (length(unique(stats::na.omit(sfx))) < 2) return(NULL)
    stats::setNames(sfx, m)
  })

  output$scale_ui <- renderUI({
    s <- method_scales()
    if (is.null(s)) return(NULL)
    lv <- sort(unique(stats::na.omit(s)))
    tagList(
      div(class = "sec", span(class = "t", "Reported scale"), span(class = "r")),
      div(class = "opts",
          radioButtons("scale", NULL, choices = lv, selected = lv[1], inline = TRUE)),
      div(class = "hint", "The same formulas on different scales, so one is shown at a time.")
    )
  })

  ## The method colours are fixed by the study, NOT by the current selection.
  ## Deriving them from the plotted subset (which is what scale_colour_manual
  ## does by default) means unticking one method silently recolours every other
  ## one, and a reader comparing two screenshots is comparing two palettes.
  all_methods <- reactive({
    s <- method_scales()
    m <- sort(unique(raw()$method))
    if (!is.null(s) && !is.null(input$scale)) m <- sort(names(s)[!is.na(s) & s == input$scale])
    m
  })

  method_pal <- reactive({
    m <- all_methods()
    stats::setNames(rep(MC_METHODS, length.out = max(1, length(m)))[seq_along(m)], m)
  })

  output$methods_ui <- renderUI({
    pal <- method_pal(); m <- names(pal)
    div(class = "opts methods",
        checkboxGroupInput(
          "methods", NULL,
          choiceNames = lapply(m, function(k)
            tags$span(class = "mk",
                      tags$i(style = paste0("background:", pal[[k]], ";")),
                      tags$span(k))),
          choiceValues = as.list(m),
          selected = m))
  })

  observeEvent(input$all_m, updateCheckboxGroupInput(session, "methods",
                                                     selected = all_methods()))
  observeEvent(input$no_m,  updateCheckboxGroupInput(session, "methods",
                                                     selected = character(0)))

  output$xvar_ui <- renderUI({
    cc <- cond_cols(); req(length(cc) > 0)
    nlev <- vapply(cc, function(k) length(unique(raw()[[k]])), integer(1))
    selectInput("xvar", "Horizontal axis", choices = cc, selected = cc[which.max(nlev)])
  })

  output$facet_ui <- renderUI({
    cc <- cond_cols(); req(length(cc) > 0)
    tagList(
      selectInput("facet_row", "Split panels by", choices = c("(none)", cc),
                  selected = if (length(cc) >= 2) cc[2] else "(none)"),
      selectInput("facet_col", "and by", choices = c("(none)", cc), selected = "(none)")
    )
  })

  output$filters_ui <- renderUI({
    cc <- cond_cols()
    rest <- setdiff(cc, c(input$xvar, input$facet_row, input$facet_col))
    if (!length(rest)) return(NULL)
    tagList(
      div(class = "sec", span(class = "t", "Held constant"), span(class = "r")),
      div(class = "held",
          lapply(rest, function(k) {
            v <- sort(unique(raw()[[k]]))
            selectInput(paste0("f_", k), k, choices = v, selected = v[1])
          }))
    )
  })

  ## input$methods is deliberately NOT req()d: unticking every method is a state
  ## the user can reach with one click, and req() would blank the whole page with
  ## no explanation. An empty selection returns zero rows and the views say so.
  dat <- reactive({
    d <- raw(); req(input$target)
    ms <- input$methods; if (is.null(ms)) ms <- character(0)
    d <- d[d$target == input$target & d$method %in% ms, , drop = FALSE]
    for (k in setdiff(cond_cols(), c(input$xvar, input$facet_row, input$facet_col))) {
      v <- input[[paste0("f_", k)]]
      if (!is.null(v)) d <- d[as.character(d[[k]]) == as.character(v), , drop = FALSE]
    }
    d
  })

  ## The condition count is a property of the design grid, not of how many
  ## methods happen to be ticked, so it is counted before the method filter --
  ## otherwise unticking a method would appear to shrink the simulation.
  n_conditions <- reactive({
    d <- raw(); req(input$target)
    d <- d[d$target == input$target, , drop = FALSE]
    for (k in setdiff(cond_cols(), c(input$xvar, input$facet_row, input$facet_col))) {
      v <- input[[paste0("f_", k)]]
      if (!is.null(v)) d <- d[as.character(d[[k]]) == as.character(v), , drop = FALSE]
    }
    cc <- cond_cols()
    if (!length(cc) || !nrow(d)) return(nrow(d))
    ## single-argument [ on a data.frame always returns a data.frame, so drop =
    ## FALSE is not only unnecessary here, it warns on every recalculation
    nrow(unique(d[cc]))
  })

  ## ---- study header ----------------------------------------------------------
  output$study_header <- renderUI({
    req(input$study)
    lab <- if (input$study %in% names(STUDY_LABELS)) STUDY_LABELS[[input$study]] else input$study
    arg <- STUDY_ARG[[input$study]]
    q   <- STUDY_Q[[input$study]]
    idx <- FILES$idx[match(input$study, FILES$stem)]
    nm  <- length(input$methods)
    ncond <- n_conditions()
    div(
      class = "hdr",
      div(class = "txt",
          div(class = "line",
              span(class = "ix", idx),
              h1(lab),
              if (!is.null(arg)) span(class = "arg", arg)),
          if (!is.null(q)) p(class = "q", q)),
      div(class = "facts",
          div(class = "f", div(class = "v", fmt_n(ncond)), div(class = "k", "Conditions")),
          div(class = "f", div(class = "v", fmt_n(input$nrep)), div(class = "k", "Replications")),
          div(class = "f", div(class = "v", fmt_n(nm)), div(class = "k", "Methods")))
    )
  })

  ## ---- headline ---------------------------------------------------------------
  ## The payoff of the whole page: which method leads on the selected measure,
  ## what its worst single condition looks like, and -- the question a reader
  ## should ask before quoting any of it -- whether the spread between methods is
  ## larger than the simulation's own noise. A difference inside the resolution
  ## limit is not a finding, and saying so here stops it being quoted as one.
  headline <- reactive({
    d <- dat(); m <- METRICS[[input$metric]]
    if (!nrow(d)) return(NULL)
    v <- suppressWarnings(as.numeric(d[[input$metric]]))
    if (!any(is.finite(v))) return(NULL)
    dev <- if (is.na(m$ref)) abs(v) else abs(v - m$ref)
    s <- split(seq_along(v), d$method)
    md <- vapply(s, function(i) mean(dev[i], na.rm = TRUE), numeric(1))
    md <- md[is.finite(md)]
    if (!length(md)) return(NULL)
    best <- names(md)[which.min(md)]
    i <- s[[best]]
    worst_i <- which.max(dev[i])
    res <- NA_real_
    if (!is.na(m$mcse) && m$mcse %in% names(d)) {
      e <- suppressWarnings(as.numeric(d[[m$mcse]]))
      res <- 2 * stats::median(e[is.finite(e)], na.rm = TRUE)
    }
    list(metric = m, best = best, best_dev = unname(md[best]),
         best_val = mean(v[i], na.rm = TRUE),
         worst_val = v[i][worst_i], worst_dev = dev[i][worst_i],
         spread = if (length(md) > 1) unname(max(md) - min(md)) else NA_real_,
         laggard = if (length(md) > 1) names(md)[which.max(md)] else NA_character_,
         nmethod = length(md), res = res)
  })

  ## Not every study records every measure -- study 06 imputes a standard error,
  ## so it has no interval and no coverage. With coverage as the default the
  ## whole headline silently disappeared on those studies, leaving the reason
  ## buried inside one tab. Name the measures the study does carry instead.
  available_metrics <- reactive({
    d <- dat(); if (!nrow(d)) return(character(0))
    Filter(function(k) k %in% names(d) &&
             any(is.finite(suppressWarnings(as.numeric(d[[k]])))),
           names(METRICS))
  })

  output$kpi_row <- renderUI({
    if (!length(input$methods) && !is.null(input$study))
      return(div(class = "note info", span(class = "ic", "i"),
                 div(HTML("No method is selected. Tick at least one in the
                           <b>Methods</b> list, or use <b>all</b>."))))
    h <- headline()
    if (is.null(h)) {
      av <- available_metrics()
      if (!length(av)) return(NULL)
      return(div(
        class = "note warn", span(class = "ic", "!"),
        div(HTML(paste0(
          "<b>", METRICS[[input$metric]]$label, " was not recorded for this study.</b> ",
          "Measures available here: ",
          paste(vapply(av, function(k) paste0("<b>", METRICS[[k]]$label, "</b>"),
                       character(1)), collapse = ", "), ".")))))
    }
    m <- h$metric
    ok <- .calibration_ok(input$metric, input$target)
    own <- grepl("^own", input$target)

    if (.scale_mismatch(input$study, input$scale, input$target)) {
      return(div(
        class = "note warn", span(class = "ic", "!"),
        div(HTML(paste0(
          "<b>The <code>", input$target, "</code> target is on the <b>",
          unname(STUDY_SHARED_SCALE[[input$study]]), "</b> scale, and you are viewing the <b>",
          input$scale, "</b> routes.</b> The two are different transforms of the same ",
          "quantity, so the gap below is arithmetic rather than a difference in ",
          "performance. Switch the reported scale back, or score against ",
          "<b>each method&rsquo;s own estimand</b>, which is scale-matched by construction.")))))
    }

    if (!ok) {
      return(div(
        class = "note warn", span(class = "ic", "!"),
        div(HTML(paste0(
          "<b>", m$label, " is not interpretable against the <code>", input$target,
          "</code> target.</b> A nominal 95% interval is not built to cover a quantity ",
          "that is itself recomputed on each replication&rsquo;s own sample. Switch to a ",
          "population parameter to read it, or choose bias or RMSE.")))))
    }

    resolvable <- is.finite(h$spread) && is.finite(h$res) && h$spread > h$res
    verdict <-
      if (!is.finite(h$spread)) span(class = "chip solid", "single method")
      else if (!is.finite(h$res)) span(class = "chip solid", "no Monte Carlo SE recorded")
      else if (resolvable) span(class = "chip solid", "resolvable")
      else span(class = "chip warn", "within simulation noise")

    tagList(
      if (own)
        div(class = "note info", span(class = "ic", "i"),
            div(HTML(paste0(
              "Every method below is scored against a <b>different</b> quantity, so this ",
              "says how accurately each computes its own value — not which to use. ",
              "See the <b>Estimand check</b> tab.")))),
      div(
        class = "kpis",
        div(class = "kpi lead",
            ## on `own` the methods are held to different standards, so "leads"
            ## would be a claim the data do not support
            div(class = "k", if (own) "Most accurate on its own estimand"
                             else "Leads on this measure"),
            div(class = "v", h$best),
            div(class = "n", HTML(paste0(
              m$label, " averages <b>", fmt_v(h$best_val),
              "</b>, a mean deviation of <b>", fmt_v(h$best_dev, 4), "</b> per condition.")))),
        div(class = "kpi",
            div(class = "k", "Its least favourable condition"),
            div(class = "v", fmt_v(h$worst_val)),
            div(class = "n", HTML(paste0(
              "The single worst cell for <b>", h$best, "</b> across the conditions shown. ",
              "This is what a review should plan for, not the average.")))),
        div(class = "kpi",
            div(class = "k", "Spread across methods"),
            div(class = "v", if (is.finite(h$spread)) fmt_v(h$spread, 4) else "--",
                if (is.finite(h$spread)) tags$u(paste0("  over ", h$nmethod, " methods")) else NULL),
            div(class = "n", HTML(paste0(
              if (is.finite(h$res))
                paste0("Resolution limit is <b>", fmt_v(h$res, 4), "</b> (2 &times; median Monte Carlo SE). ")
              else "",
              if (is.finite(h$spread) && !is.na(h$laggard))
                paste0("Widest gap: <b>", h$best, "</b> to <b>", h$laggard, "</b>.") else ""))),
            verdict)
      )
    )
  })

  ## ---- key row ---------------------------------------------------------------
  ## A persistent, HTML method key. ggplot draws its own legend in a font and at
  ## a size nothing else on the page uses; this one shares the swatch colours
  ## with the sidebar checkboxes and the ranking table, and survives being
  ## screenshotted with the plot.
  output$key_row <- renderUI({
    d <- dat()
    if (!nrow(d)) return(NULL)
    m <- METRICS[[input$metric]]
    pal <- method_pal()
    ms <- intersect(names(pal), unique(d$method))
    ## if the study never recorded this measure there is no reference line and no
    ## band to key, and claiming otherwise above an empty panel is a small lie
    have <- input$metric %in% available_metrics()
    tagList(
      div(class = "key",
          span(class = "lbl", "Methods"),
          lapply(ms, function(k)
            span(class = "m", tags$i(style = paste0("background:", pal[[k]], ";")), k)),
          if (have && !is.na(m$mcse))
            span(class = "band",
                 title = paste("Half-width is twice the median Monte Carlo SE of the cells",
                               "in that panel. Differences inside the band are smaller than",
                               "the simulation own noise and are not results."),
                 tags$i(), "resolution limit, per panel")
          else NULL),
      div(class = "measline",
          HTML(paste0("<b>", m$label, "</b> against the <b>", pretty_target(input$target),
                      "</b> target",
                      if (have && !is.na(m$ref))
                        paste0(", with the dashed line at ", m$ref, " marking no error")
                      else "", ".")))
    )
  })

  ## Panel count drives the height. A fixed 700px leaves a two-panel study
  ## floating in white space and squashes a nine-panel one.
  n_panels <- reactive({
    d <- dat(); if (!nrow(d)) return(1)
    fr <- if (identical(input$facet_row, "(none)")) NULL else input$facet_row
    fc <- if (identical(input$facet_col, "(none)")) NULL else input$facet_col
    keys <- c(fr, fc)
    if (!length(keys)) return(1)
    nrow(unique(d[keys]))
  })

  output$plot_slot <- renderUI({
    np <- n_panels()
    ncol <- if (np <= 1) 1 else 2
    rows <- ceiling(np / ncol)
    ## enough per row that the band and the error bars are legible, capped so a
    ## nine-condition facet does not become a two-screen scroll
    h <- min(1000, max(320, 260 * rows + 30))
    plotOutput("main_plot", height = paste0(h, "px"))
  })

  ## ---- main plot -------------------------------------------------------------
  output$main_plot <- renderPlot({
    d <- dat()
    validate(need(length(input$methods) > 0, "No method selected."))
    validate(need(nrow(d) > 0, "No rows match these filters."))
    m <- METRICS[[input$metric]]
    req(input$xvar)

    d$.x <- d[[input$xvar]]
    d$.y <- suppressWarnings(as.numeric(d[[input$metric]]))
    validate(need(any(is.finite(d$.y)),
                  paste0(m$label, " was not recorded for this study.")))

    fr <- if (identical(input$facet_row, "(none)")) NULL else input$facet_row
    fc <- if (identical(input$facet_col, "(none)")) NULL else input$facet_col
    keys <- c(fr, fc)

    p <- ggplot(d, aes(.x, .y, colour = method, group = method))

    ## SIGNATURE: the Monte Carlo resolution band. Morris, White & Crowther (2019)
    ## report every measure with its MCSE precisely so a difference smaller than
    ## the simulation's own noise is not mistaken for a finding; this draws that
    ## threshold instead of leaving the reader to compute it. Half-width is 2x the
    ## median MCSE, the k = 2 convention of mc_resolvable() in R/03_performance.R.
    ##
    ## Computed PER PANEL, not once for the whole plot: MCSE shrinks with the
    ## sample size and the replication count, so a single global band would be far
    ## too narrow in the small-n panels and too wide in the large-n ones -- i.e. it
    ## would misstate resolvability in exactly the panels where the question
    ## matters most. Per panel, the band visibly narrows as precision improves.
    if (!is.na(m$mcse) && m$mcse %in% names(d) && !is.na(m$ref)) {
      e <- suppressWarnings(as.numeric(d[[m$mcse]]))
      grp <- if (length(keys)) interaction(d[keys], drop = TRUE) else factor(rep(1, nrow(d)))
      band <- do.call(rbind, lapply(split(seq_len(nrow(d)), grp), function(i) {
        h <- 2 * stats::median(e[i][is.finite(e[i])], na.rm = TRUE)
        if (!is.finite(h) || h <= 0) return(NULL)
        cbind(d[i[1], keys, drop = FALSE],
              data.frame(ymin = m$ref - h, ymax = m$ref + h), row.names = NULL)
      }))
      if (!is.null(band) && nrow(band))
        p <- p + geom_rect(data = band, inherit.aes = FALSE,
                           aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
                           fill = MC$ink2, alpha = 0.07)
    }
    if (!is.na(m$ref))
      p <- p + geom_hline(yintercept = m$ref, linetype = "22",
                          colour = MC$ink2, linewidth = 0.5)

    if (!is.na(m$mcse) && m$mcse %in% names(d)) {
      e <- suppressWarnings(as.numeric(d[[m$mcse]]))
      p <- p + geom_linerange(aes(ymin = .y - e, ymax = .y + e),
                              alpha = 0.5, linewidth = 0.45)
    }

    p <- p +
      geom_line(alpha = 0.9, linewidth = 0.6) +
      geom_point(size = 2.2, stroke = 0) +
      scale_colour_manual(values = method_pal(), guide = "none") +
      scale_y_continuous(labels = function(x) format(x, trim = TRUE)) +
      labs(x = input$xvar, y = NULL) +
      theme_mc()

    if (length(keys)) {
      ## facet_wrap keeps the strips on TOP and horizontal. facet_grid puts row
      ## strips on the right rotated 90 degrees, which is unreadable at these
      ## label lengths.
      p <- p + facet_wrap(keys, labeller = label_both,
                          ncol = if (length(keys) > 1) NULL else 2,
                          scales = if (isTRUE(input$free_y)) "free_y" else "fixed")
    }
    p
  }, res = 108)

  ## ---- ranking ---------------------------------------------------------------
  ## A scorecard, not a single ranking. No one measure is sufficient: a
  ## conversion can be almost unbiased in the point estimate and badly
  ## miscalibrated in its standard error, and the two orderings genuinely
  ## disagree in these data. In study 09a the best method on bias
  ## (2x2_tetrachoric) is the WORST on coverage, and across the shipped studies
  ## the rank correlation between |bias| and |coverage - 0.95| is weak and
  ## sometimes negative. So all four measures are shown side by side; the
  ## selected measure only decides the sort order.
  output$rank_lede <- renderUI({
    m <- METRICS[[input$metric]]
    tagList(
      div(class = "lede",
          HTML(paste0(
            "Every measure for every method, averaged over the conditions currently shown, ",
            "sorted by <b>", m$label, "</b> (", m$better, "). Read them together: ",
            "<b>bias</b> is the point estimate only, <b>RMSE</b> adds the estimator&rsquo;s ",
            "variability, <b>SE ratio</b> asks whether the reported standard error matches ",
            "the real one, and <b>coverage</b> is the only measure that degrades under ",
            "<i>either</i> a wrong estimate or a wrong standard error — which is why it ",
            "is the default. The sort column is the mean deviation <i>per condition</i>, ",
            "which is not the deviation of the mean beside it: a method covering 0.90 half ",
            "the time and 1.00 the rest averages to a perfect 0.95 while being miscalibrated ",
            "throughout. <code>worst</code> is the least favourable single condition, which ",
            "is what a review should plan for."))),
      if (.scale_mismatch(input$study, input$scale, input$target))
        div(class = "note warn", span(class = "ic", "!"),
            div(HTML(paste0("The <b>", input$target, "</b> target is on the <b>",
                            unname(STUDY_SHARED_SCALE[[input$study]]), "</b> scale while these ",
                            "routes report <b>", input$scale, "</b>. This ranking is of two ",
                            "different transforms; switch the scale back or rank on ",
                            "<b>each method&rsquo;s own estimand</b>."))))
      else if (grepl("sample", input$target))
        div(class = "note warn", span(class = "ic", "!"),
            div(HTML(paste0("Coverage and the SE ratio are <b>not interpretable</b> against the ",
                            "<b>", input$target, "</b> target: a nominal 95% interval is not built ",
                            "to cover a quantity that is itself random. Switch to the population ",
                            "parameter to read them."))))
      else if (grepl("^own", input$target))
        div(class = "note info", span(class = "ic", "i"),
            div(HTML(paste0("Every method here is scored against a <b>different</b> quantity, so ",
                            "this table says how accurately each computes its own value — not ",
                            "which method to use. Switch to the population parameter to compare ",
                            "them, and see the Estimand check tab for the difference."))))
      else NULL
    )
  })

  rank_df <- reactive({
    d <- dat(); m <- METRICS[[input$metric]]
    validate(need(nrow(d) > 0, "No rows match these filters."))
    sortv <- suppressWarnings(as.numeric(d[[input$metric]]))
    sortdev <- if (is.na(m$ref)) abs(sortv) else abs(sortv - m$ref)
    num <- function(col) suppressWarnings(as.numeric(d[[col]]))
    s <- split(seq_along(sortv), d$method)
    out <- do.call(rbind, lapply(names(s), function(k) {
      i <- s[[k]]
      cov_i <- num("coverage")[i]
      data.frame(
        method       = k,
        conditions   = length(i),
        ## The sort key, shown because it is NOT the same as a deviation computed
        ## from the mean column below: a method covering 0.90 in half the
        ## conditions and 1.00 in the rest averages to a perfect 0.95 while being
        ## miscalibrated everywhere. Averaging the per-condition deviation catches
        ## that; deviating the average does not.
        sortkey      = round(mean(sortdev[i], na.rm = TRUE), 4),
        bias         = round(mean(num("bias")[i], na.rm = TRUE), 4),
        `abs bias`   = round(mean(abs(num("bias")[i]), na.rm = TRUE), 4),
        rmse         = round(mean(num("rmse")[i], na.rm = TRUE), 4),
        `SE ratio`   = round(mean(num("se_ratio")[i], na.rm = TRUE), 3),
        coverage     = round(mean(cov_i, na.rm = TRUE), 3),
        `worst coverage` = round(suppressWarnings(min(cov_i, na.rm = TRUE)), 3),
        `max non-estimable` = round(suppressWarnings(max(num("nonest_rate")[i], na.rm = TRUE)), 3),
        .sort        = mean(sortdev[i], na.rm = TRUE),
        check.names  = FALSE, stringsAsFactors = FALSE)
    }))
    out <- out[order(out$.sort), ]
    out$.sort <- NULL
    ## numeric columns only: is.finite() is FALSE for every element of a
    ## character vector, so applying this to `method` blanks the whole column.
    isnum <- vapply(out, is.numeric, logical(1))
    out[isnum] <- lapply(out[isnum], function(x) { x[!is.finite(x)] <- NA; x })
    names(out)[names(out) == "sortkey"] <-
      if (is.na(m$ref)) paste0("mean |", input$metric, "|")
      else paste0("mean |", input$metric, " - ", m$ref, "|")
    out <- cbind(`#` = seq_len(nrow(out)), out)
    rownames(out) <- NULL
    out
  })

  output$rank_tbl <- renderDT({
    r <- rank_df()
    pal <- method_pal()
    dt <- datatable(r, rownames = FALSE,
                    options = list(dom = "t", pageLength = 30, ordering = TRUE,
                                   columnDefs = list(
                                     list(className = "rk", targets = 0),
                                     list(className = "mth", targets = 1))))
    ## the same swatch as the plot key and the sidebar, carried as a left edge so
    ## the column still sorts as text
    known <- intersect(r$method, names(pal))
    if (length(known))
      dt <- formatStyle(dt, "method",
                        borderLeft = styleEqual(known,
                                                paste0("3px solid ", unname(pal[known])),
                                                default = "3px solid transparent"))
    ## bar on the absolute bias, and a colour cue on coverage: red as it falls
    ## away from nominal in either direction.
    dt <- formatStyle(dt, "abs bias",
                      background = styleColorBar(range(c(0, r[["abs bias"]]), na.rm = TRUE),
                                                 "rgba(214,98,104,.18)"),
                      backgroundSize = "98% 62%", backgroundRepeat = "no-repeat",
                      backgroundPosition = "center")
    dt <- formatStyle(dt, "coverage",
                      color = styleInterval(c(0.90, 0.93, 0.97),
                                            c("#9c2b2b", "#8a6712", "#211f38", "#8a6712")),
                      fontWeight = "bold")
    ## a column that prints 0.815 beside a bare 0 does not read as a table of
    ## comparable quantities. Fixed decimals per column, the sort key included
    ## (its name is built from the selected measure, so it is taken by position).
    dt <- formatRound(dt, intersect(names(r), c(names(r)[4], "bias", "abs bias", "rmse")), 4)
    dt <- formatRound(dt, intersect(names(r), c("SE ratio", "coverage",
                                                "worst coverage", "max non-estimable")), 3)
    dt
  })

  ## ---- estimand check --------------------------------------------------------
  output$estimand_lede <- renderUI({
    tagList(
      div(class = "lede",
          HTML(paste0(
            "Two questions have to be kept apart. <b>Is the method computing its own ",
            "quantity correctly?</b> — that is the dark marker. <b>Is that quantity the ",
            "one you asked for?</b> — the distance to the rose marker. A marker far from ",
            "zero with a short bar is a method with a real accuracy problem. A marker at ",
            "zero with a long bar is a method working perfectly on a different estimand, ",
            "and reporting that bar as bias is a mistake. Bars are the mean absolute bias ",
            "over the conditions currently shown."))))
  })

  estimand_df <- reactive({
    d0 <- raw()
    validate(need("own" %in% d0$target,
                  "This study did not record an 'own' target, so there is no estimand gap to show."))
    pop <- grep("population", unique(d0$target), value = TRUE)
    validate(need(length(pop) > 0, "This study records no population-scale target."))
    validate(need(length(input$methods) > 0, "No method selected."))
    a <- d0[d0$target == "own" & d0$method %in% input$methods, , drop = FALSE]
    b <- d0[d0$target == pop[1] & d0$method %in% input$methods, , drop = FALSE]
    key <- function(x) do.call(paste, c(x[c(cond_cols(), "method")], sep = "\r"))
    i <- match(key(a), key(b)); ok <- !is.na(i)
    validate(need(any(ok), "No matching own/population cells."))
    df <- data.frame(method = a$method[ok], bias_own = a$bias[ok], bias_pop = b$bias[i[ok]])

    s <- do.call(rbind, lapply(split(df, df$method), function(x) data.frame(
      method = x$method[1],
      own    = mean(abs(x$bias_own), na.rm = TRUE),
      pop    = mean(abs(x$bias_pop), na.rm = TRUE),
      stringsAsFactors = FALSE)))
    s <- s[is.finite(s$own) & is.finite(s$pop), , drop = FALSE]
    validate(need(nrow(s) > 0, "No comparable cells."))
    s$gap <- s$pop - s$own
    s <- s[order(s$pop), ]
    s$method <- factor(s$method, levels = s$method)
    s
  })

  output$estimand_slot <- renderUI({
    n <- tryCatch(nrow(estimand_df()), error = function(e) 3)
    ## each row carries two markers and two mono labels stacked vertically, so it
    ## needs more than a bar chart row would
    plotOutput("estimand_plot", height = paste0(max(260, 74 * n + 110), "px"))
  })

  output$estimand_plot <- renderPlot({
    s <- estimand_df()
    pal <- method_pal()

    ## One row per method: how far it misses its OWN estimand (computational
    ## error) and how far it misses the POPULATION parameter (the total). The
    ## distance between the two is the estimand gap. A dumbbell reads far more
    ## directly than the bias-vs-bias scatter this replaced, where the reader had
    ## to infer the gap from a point's distance off the identity line.
    lab <- function(v) formatC(v, format = "f", digits = 3)
    ggplot(s) +
      geom_segment(aes(y = method, yend = method, x = own, xend = pop),
                   colour = MC$rule, linewidth = 2.8, lineend = "round") +
      geom_point(aes(y = method, x = pop), colour = MC$rose, size = 3.8) +
      geom_point(aes(y = method, x = own), colour = MC$ink, size = 3.8) +
      ## own above the dot, population below: when a method has no estimand gap
      ## the two dots coincide, and labels on the same side would overlap into a
      ## dark dot carrying a rose number.
      geom_text(aes(y = method, x = own, label = lab(own)),
                colour = MC$ink, size = 3.1, vjust = -1.5, family = "Cascadia Code") +
      geom_text(aes(y = method, x = pop, label = lab(pop)),
                colour = MC$rose_d, size = 3.1, vjust = 2.3, family = "Cascadia Code") +
      ## start at zero so the dot positions and the gap length are both read
      ## against a true origin
      scale_x_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0.02, 0.12))) +
      labs(x = "Mean absolute bias", y = NULL,
           subtitle = "dark = against the method own estimand   |   rose = against the population parameter") +
      theme_mc() +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_text(family = "Cascadia Code", size = 11.5,
                                       colour = MC$ink),
            plot.subtitle = element_text(size = 10.5, colour = MC$soft,
                                         margin = margin(0, 0, 14, 0)))
  }, res = 108)

  ## ---- data ------------------------------------------------------------------
  output$data_lede <- renderUI({
    ## emp_se and rmse sit next to each other here, and on a random target they do
    ## not compose the way readers expect. Say so at the point of confusion.
    random_target <- !is.null(input$target) && !grepl("population", input$target)
    tagList(
      div(class = "lede",
          "The rows behind the current view, with every ADEMP measure and its Monte Carlo standard error."),
      if (random_target)
        div(class = "note info", span(class = "ic", "i"),
            div(HTML(paste0(
              "Note on <code>rmse</code> and <code>emp_se</code>: the <b>", input$target,
              "</b> target is recomputed on each replication&rsquo;s own sample, so it is not a ",
              "fixed constant and <code>RMSE&sup2; = bias&sup2; + EmpSE&sup2;</code> does not apply. ",
              "The estimate and the target are strongly correlated, so <code>rmse</code> can sit ",
              "below <code>emp_se</code> with nothing wrong. Compare <code>rmse</code> only with ",
              "other <code>rmse</code> values on the same target."))))
      else NULL
    )
  })

  output$data_tbl <- renderDT({
    d <- dat()
    num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], round, 4)
    datatable(d, rownames = FALSE, filter = "top",
              options = list(pageLength = 20, scrollX = TRUE))
  })

  output$dl <- downloadHandler(
    filename = function() paste0(input$study, "_nrep", input$nrep, "_", input$target, ".csv"),
    content = function(f) utils::write.csv(dat(), f, row.names = FALSE)
  )

  ## ---- provenance ------------------------------------------------------------
  ## A results viewer that accompanies a paper has to say where its numbers came
  ## from. This is the line a reviewer looks for.
  output$footer <- renderUI({
    req(input$study, input$nrep)
    p <- FILES$path[FILES$stem == input$study & FILES$nrep == as.integer(input$nrep)]
    if (!length(p)) return(NULL)
    div(class = "foot",
        span(HTML(paste0("<b>Source</b> <span class=\"mono\">", basename(p[1]), "</span>"))),
        span(HTML(paste0("<b>Aggregated</b> <span class=\"mono\">",
                         format(file.info(p[1])$mtime, "%Y-%m-%d"), "</span>"))),
        span(HTML(paste0("<b>Measures</b> ADEMP, after Morris, White &amp; Crowther (2019)"))),
        span(HTML(paste0("<b>Produced by</b> <span class=\"mono\">run_study()</span>, ",
                         "calling metaConvert directly — no formula is re-implemented here"))))
  })
}

shinyApp(ui, server)
