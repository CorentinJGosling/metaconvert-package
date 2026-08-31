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
  "09b_or_to_cor_CAT"     = "Odds ratio to correlation (categorical latent)",
  "10a_icc_agreement_coverage" = "ICC agreement-type interval coverage"
)

## Every lookup into the STUDY_* maps goes through this. `x[["absent"]]` on a
## named vector is an ERROR, not NULL, so a study that appears in
## data/aggregated before its prose is written used to abort the page with
## "subscript out of bounds" -- which is what 10a did the moment it shipped.
.lookup <- function(map, key) {
  if (is.null(key) || !length(key) || !key %in% names(map)) return(NULL)
  unname(map[[key]])
}

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
  "09a_or_to_cor_CONT" = "or_to_cor", "09b_or_to_cor_CAT" = "or_to_cor",
  "10a_icc_agreement_coverage" = "icc_type"
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
  "09b_or_to_cor_CAT"     = "The same conversion when the latent variable is categorical, so every route is estimating something it was not designed for.",
  "10a_icc_agreement_coverage" = "<code>icc_type</code> changes which ICC is being estimated but not the standard error metaConvert computes for it. What does the shared formula cost the agreement-type interval?"
)

## ---- what each route estimates -----------------------------------------------
## The block that replaced the two rhetorical questions of the previous design.
##
## IT IS DECLARED, NOT DERIVED, AND IT HAS TO BE. Nothing in an aggregate file
## says what a method targets: `method` is a bare string and `target` names the
## benchmark, not the estimand the route was built to return. The distinction is
## the whole subject of studies 01, 02, 03, 08 and 09 -- and it is settled in the
## QUESTION header of each studies/*.R, which is the authority these entries are
## transcribed from. Every line below can be checked against that header; nothing
## here is inferred from the numbers.
##
## `kind` separates what metaConvert SHIPS from what the study added as a
## comparator. A reader choosing an argument value needs to know that
## `metafor_rtod_np` is not one, and the aggregates do not record it (run_study()
## carries a `route` column that aggregate_df() drops).
##
## A study missing from ROUTES still works: the block degrades to the plain route
## selector, which is what the previous design had. Adding a study means adding a
## header line here, not touching any code.

.rt <- function(...) {
  a <- list(...)
  data.frame(method = names(a),
             what   = vapply(a, `[`, "", 1),
             ref    = vapply(a, `[`, "", 2),
             kind   = vapply(a, `[`, "", 3),
             row.names = NULL, stringsAsFactors = FALSE)
}

PKG <- "shipped"      ## a value convert_df() accepts for this argument
CND <- "comparator"   ## added by the study; NOT reachable from the package
OTH <- "other route"  ## shipped, but through a different argument

## The heading over the block. It states the relation between the routes, which
## is different in kind from study to study: 01/02/03/08/09 are estimand splits,
## 04/05/06 are competing approximations of one quantity, 07 is a choice of
## standardiser, 10 is a choice of ICC model.
ROUTE_HEAD <- c(
  "01a_smd_to_cor_r"      = "The two routes estimate different correlations",
  "01b_smd_to_cor_z"      = "The two routes estimate different correlations",
  "02a_cor_to_smd_GROUPS" = "Each route presupposes a different dichotomisation",
  "02b_cor_to_smd_CONT"   = "Each route presupposes a different dichotomisation",
  "03a_2x2_to_cor_CAT"    = "Two estimands, each reported on two scales",
  "03b_2x2_to_cor_CONT"   = "Two estimands, each reported on two scales",
  "04_or_to_rr"           = "Approximations of one risk ratio, differing in what they need",
  "05_rr_to_or"           = "Approximations of one odds ratio, differing in what they need",
  "06_or_se_imputation"   = "Five ways to recover the same standard error",
  "07a_ancova_to_smd_d"   = "Three ANCOVA denominators, and the fallback that needs no guess",
  "07b_ancova_to_smd_g"   = "Three ANCOVA denominators, and the fallback that needs no guess",
  "08a_pre_post_to_smd_d" = "One numerator, five standardisers",
  "08b_pre_post_to_smd_g" = "One numerator, five standardisers",
  "09a_or_to_cor_CONT"    = "Two estimands, each reported on two scales",
  "09b_or_to_cor_CAT"     = "Two estimands, each reported on two scales",
  "10a_icc_agreement_coverage" = "Two ICC models, one shared standard-error formula"
)

## studies/01_smd_to_cor.R, QUESTION block.
.r01 <- .rt(
  lipsey_cooper = c("point-biserial correlation",
                    "Cooper et al. (2019), eq. 12.40-12.42; Borenstein et al. (2021), eq. 54-56", PKG),
  viechtbauer   = c("biserial correlation",
                    "Jacobs & Viechtbauer (2017), eq. 5, 8, 13, 17-19", PKG))

## studies/02_cor_to_smd.R. The nine labels are one two-parameter family
## (k = correlation-scale multiplier, p = assumed group proportion); the three
## SHIPPED values of cor_to_smd are cooper, viechtbauer and mathur.
.r02 <- .rt(
  cooper           = c("d of a balanced split, treating r as point-biserial",
                       "Cooper et al. (2019), eq. 12.38-12.39; Mathur & VanderWeele (2020), eq. 1.1", PKG),
  cooper_g         = c("the same, Hedges-corrected", "as cooper, with J", PKG),
  viechtbauer      = c("d of a MEDIAN split, treating r as biserial",
                       "metafor::transf.rtod (conv.delta)", PKG),
  viechtbauer_g    = c("the same, Hedges-corrected", "as viechtbauer, with J", PKG),
  mathur_2sd       = c("d per 2 SD of a continuous exposure",
                       "Mathur & VanderWeele (2020), eq. 1.2", PKG),
  mathur_raw_1unit = c("d per 1 raw unit of the exposure",
                       "Mathur & VanderWeele (2020), eq. 1.2", PKG),
  pb_exact_p       = c("exact point-biserial inverse at the OBSERVED p",
                       "study comparator: cooper without the balanced-arm assumption", CND),
  biserial_split_p = c("exact biserial inverse at the observed p",
                       "study comparator: viechtbauer without the median-split assumption", CND),
  metafor_rtod_np  = c("biserial inverse with metafor arm sizes",
                       "study comparator: metafor::transf.rtod with n1i/n2i supplied", CND))

## studies/03_2x2_to_cor.R. The (r)/(z) suffix is the reported scale, not a
## different estimand; the app offers it as a scale switch.
.r03 <- .rt(
  "phi (r)"         = c("phi, the binary-binary correlation",
                        "study comparator: table_2x2_to_cor accepts only 'tetrachoric'", CND),
  "phi (z)"         = c("the same, on the Fisher z scale", "atanh(phi)", CND),
  "tetrachoric (r)" = c("tetrachoric correlation of the two latent normals",
                        "table_2x2_to_cor = 'tetrachoric' (mvtnorm)", PKG),
  "tetrachoric (z)" = c("the same, on the Fisher z scale", "atanh(tetrachoric)", PKG))

## studies/04_or_to_rr.R, "NOTE ON grant AND ZHANG-YU" and the input-need list.
.r04 <- .rt(
  metaumbrella_cases  = c("RR from the reconstructed table, cases margin",
                          "or_to_rr = 'metaumbrella_cases' (package default)", PKG),
  metaumbrella_exp    = c("RR from the reconstructed table, exposed margin",
                          "or_to_rr = 'metaumbrella_exp'", PKG),
  grant               = c("RR = OR / (1 - p0 + p0*OR); needs the baseline risk",
                          "Zhang & Yu (1998), restated by Grant (2014); or_to_rr = 'grant'", PKG),
  dipietrantonj       = c("the same transform applied through the interval",
                          "Di Pietrantonj (2006), Stat Med 25(13); or_to_rr = 'dipietrantonj'", PKG),
  transpose           = c("the OR used unchanged as an RR", "or_to_rr = 'transpose'", PKG),
  vanderweele_sqrt_or = c("RR ~ sqrt(OR), minimax over p in [0.2, 0.8]",
                          "study comparator: VanderWeele (2020), Biometrics 76(3):746-752", CND),
  metafor_conv2x2     = c("RR from metafor's table reconstruction",
                          "study comparator: metafor::conv.2x2", CND))

## studies/05_rr_to_or.R. Grant's inverse is undefined when RR * p0 >= 1; the
## rate of that is the nonest_rate column, and it is a result rather than a gap.
.r05 <- .rt(
  metaumbrella           = c("OR from the reconstructed table",
                             "rr_to_or = 'metaumbrella' (package default)", PKG),
  grant                  = c("OR = RR(1 - p0) / (1 - RR*p0); non-estimable when RR*p0 >= 1",
                             "Zhang & Yu (1998), restated by Grant (2014); rr_to_or = 'grant'", PKG),
  dipietrantonj          = c("the same transform applied through the interval",
                             "Di Pietrantonj (2006), Stat Med 25(13); rr_to_or = 'dipietrantonj'", PKG),
  transpose              = c("the RR used unchanged as an OR", "rr_to_or = 'transpose'", PKG),
  vanderweele_rr_squared = c("OR ~ RR^2, the inverse minimax conversion",
                             "study comparator: VanderWeele (2020), Biometrics 76(3):746-752", CND))

## studies/06_or_se_imputation.R. All five recover the same scalar,
## sqrt(1/a + 1/b + 1/c + 1/d) for the table behind the reported OR.
.r06 <- .rt(
  package_mean_var         = c("mean Woolf variance over every compatible 2x2",
                               ".se_from_or(); es_from_or(or, n_cases, n_controls)", PKG),
  median_var               = c("median of that same variance set",
                               "study comparator: one-line change to .se_from_or", CND),
  metafor_conv2x2_oracle   = c("reconstruction given the true exposed margin",
                               "study comparator: metafor::conv.2x2 (upper bound, not feasible)", CND),
  metafor_conv2x2_balanced = c("reconstruction assuming balanced arms",
                               "study comparator: metafor::conv.2x2 at n/2", CND),
  balanced_4_over_sqrtN    = c("4 / sqrt(N), the Woolf SE of an all-equal table",
                               "study comparator; the constant the D2 flag divides by", CND))

## studies/07_ancova_to_smd.R, "STANDARDISER CHOICES INCLUDED".
.r07 <- .rt(
  ancova_means_sd     = c("per-arm adjusted SDs, df-pooled, back-transformed with the guess",
                          "Cooper et al. (2019), eq. 12.24, table 12.3", PKG),
  ancova_md_sd        = c("ANCOVA root MSE as one pooled residual SD, back-transformed with the guess",
                          "Cooper et al. (2019), eq. 12.24, table 12.3", PKG),
  ancova_pooled_crude = c("adjusted means over the CRUDE pooled SD; the guess enters only the variance",
                          "Cooper et al. (2019), eq. 12.26", PKG),
  crude_means_sd      = c("the ANCOVA ignored; unadjusted endpoint means and SDs",
                          "es_from_means_sd(); no nuisance parameter", PKG))

## studies/08_pre_post_to_smd.R. Same numerator, five denominators.
.r08 <- .rt(
  bonett            = c("difference in mean change over the pooled BASELINE SD",
                        "Bonett (2008), Psych Methods 13(2); pre_post_to_smd = 'bonett' (default)", PKG),
  morris_dav        = c("over the quadratic mean of baseline and endpoint SD",
                        "Morris (2008), d_av; variance from Bonett (2008), eq. 10/19", PKG),
  morris_dz         = c("over the pooled CHANGE-SCORE SD", "Morris & DeShon (2002), d_z", PKG),
  morris_drm        = c("change SD rescaled by sqrt(2(1 - r))",
                        "Cooper et al. (2019); Morris & DeShon (2002), d_rm", PKG),
  cooper            = c("documented alias of morris_drm",
                        "same formula; ?es_from_means_sd_pre_post lists the two together", PKG),
  endpoint_means_sd = c("baseline discarded; endpoint means and SDs only",
                        "es_from_means_sd()", PKG))

## studies/09_or_to_cor.R, QUESTION block. Three approximations of one estimand
## plus one route targeting a different quantity, each on two scales.
.r09 <- .rt(
  "bonett (r)"           = c("tetrachoric, cosine exponent adjusted for the margins",
                             "Bonett & Price (2005); or_to_cor = 'bonett' (default)", PKG),
  "bonett (z)"           = c("the same, Fisher z scale", "atanh", PKG),
  "pearson (r)"          = c("tetrachoric, r = cos(pi / (1 + OR^0.5))", "or_to_cor = 'pearson'", PKG),
  "pearson (z)"          = c("the same, Fisher z scale", "atanh", PKG),
  "digby (r)"            = c("tetrachoric, r = (OR^0.75 - 1)/(OR^0.75 + 1)", "or_to_cor = 'digby'", PKG),
  "digby (z)"            = c("the same, Fisher z scale", "atanh", PKG),
  "lipsey_cooper (r)"    = c("point-biserial, via the Cox logit d = log(OR) sqrt(3)/pi",
                             "or_to_cor = 'lipsey_cooper'", PKG),
  "lipsey_cooper (z)"    = c("the same, Fisher z scale", "atanh", PKG),
  "2x2_tetrachoric (r)"  = c("tetrachoric from the FULL 2x2, never passing through the OR",
                             "es_from_2x2(table_2x2_to_cor = 'tetrachoric'); the reference route", OTH),
  "2x2_tetrachoric (z)"  = c("the same, Fisher z scale", "atanh", OTH))

## studies/10_reliability.R. icc_type changes the estimand and NOT the arithmetic:
## the two share one SE formula, which is the point 10a measures.
.r10 <- .rt(
  agreement   = c("ICC(2,1), absolute agreement; charges rater differences against reliability",
                  "icc_type = 'agreement' (package default)", PKG),
  consistency = c("ICC(3,1), consistency; does not", "icc_type = 'consistency'", PKG))

ROUTES <- list(
  "01a_smd_to_cor_r" = .r01, "01b_smd_to_cor_z" = .r01,
  "02a_cor_to_smd_GROUPS" = .r02, "02b_cor_to_smd_CONT" = .r02,
  "03a_2x2_to_cor_CAT" = .r03, "03b_2x2_to_cor_CONT" = .r03,
  "04_or_to_rr" = .r04, "05_rr_to_or" = .r05, "06_or_se_imputation" = .r06,
  "07a_ancova_to_smd_d" = .r07, "07b_ancova_to_smd_g" = .r07,
  "08a_pre_post_to_smd_d" = .r08, "08b_pre_post_to_smd_g" = .r08,
  "09a_or_to_cor_CONT" = .r09, "09b_or_to_cor_CAT" = .r09,
  "10a_icc_agreement_coverage" = .r10
)

## One row of the block, for a method that may or may not be declared. An
## undeclared method still gets its swatch and its tick, so a study added to
## data/aggregated before its header is transcribed here is usable immediately.
route_note <- function(study, method) {
  r <- ROUTES[[study]]
  if (is.null(r)) return(NULL)
  i <- match(method, r$method)
  if (is.na(i)) return(NULL)
  as.list(r[i, ])
}

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
## Which shared target, if any, each route is ALREADY estimating.
##
## Exact and threshold-free, and it is the same test the app has always used to
## label an unattained target: run_study() writes one row per (condition, method,
## target), so a route whose own estimand IS a named shared target produces the
## identical bias column against both. Comparing the two ordered vectors answers
## "is this column the quantity that route was built to return?" with no tolerance
## to choose.
##
## SAME-SAMPLE AND own-LIKE TARGETS ARE EXCLUDED. `own` is a population quantity in
## every study, so a `*sample*` target -- the same statistic recomputed on the
## replication draw -- can never equal it; and a match against another own-like
## column (study 08 own_true_r) says something about r-misspecification rather than
## about which shared quantity the route targets.
own_matches <- function(agg) {
  tg <- unique(agg$target)
  if (!"own" %in% tg || !"bias" %in% names(agg)) return(character(0))
  shared <- setdiff(tg, grep("^own", tg, value = TRUE))
  shared <- shared[!grepl("sample", shared)]
  k <- condition_cols(agg)
  if (!length(shared) || !length(k)) return(character(0))
  vec <- function(m, t) {
    d <- agg[agg$method == m & agg$target == t, , drop = FALSE]
    d[do.call(order, d[k]), "bias"]
  }
  ms <- unique(agg$method)
  stats::setNames(vapply(ms, function(m) {
    hit <- shared[vapply(shared, function(t)
      isTRUE(all.equal(vec(m, t), vec(m, "own"))), logical(1))]
    if (length(hit)) hit[1] else NA_character_
  }, character(1)), ms)
}

## Shared targets that NO route in this study estimates.
##
## The answer is worth showing either way, and the label is deliberately neutral. In
## 01b it says `fisherz_biserial` is an artefact of the target naming; in 08 and 09b it
## says no route recovers the quantity the review wants, which is those studies result.
unattained_targets <- function(agg) {
  tg <- unique(agg$target)
  if (!"own" %in% tg) return(character(0))
  shared <- setdiff(tg, grep("^own", tg, value = TRUE))
  shared <- shared[!grepl("sample", shared)]
  hit <- own_matches(agg)
  setdiff(shared, stats::na.omit(unname(hit)))
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
'--rule2:#eeeee9;--rose-l:#f6e7e7;--rose-m:#eec9ca;--amber-l:#fbf6e9;--red-l:#fdf0f0;
--mono:"Cascadia Code","Cascadia Mono",ui-monospace,SFMono-Regular,Menlo,monospace;}

body{background:var(--paper);color:var(--ink);font-size:14px;}
::selection{background:rgba(192,71,78,.16);}
.num,table td,table th,.mono,code{
  font-variant-numeric:tabular-nums lining-nums;
  font-feature-settings:"tnum" 1,"lnum" 1;}
code,.mono{font-family:var(--mono);}
:focus-visible{outline:2px solid var(--rose);outline-offset:2px;}

/* =========================================================== masthead ===== */
.mast{
  display:flex;align-items:center;gap:10px;background:#fff;
  border-bottom:1px solid var(--ink);padding:11px 26px;
}
.mast .wm{font-size:13.5px;font-weight:600;letter-spacing:-.01em;}
.mast .wm em{font-style:normal;color:var(--faint);font-weight:400;}
.mast .rt{margin-left:auto;font-family:var(--mono);font-size:11.5px;color:var(--faint);
  display:flex;align-items:center;gap:14px;}
.helpbtn{
  background:none;border:0;border-bottom:1px solid var(--rule);padding:0 0 1px;
  font-family:var(--mono);font-size:11.5px;color:var(--faint);cursor:pointer;
}
.helpbtn:hover{color:var(--ink);border-bottom-color:var(--ink);}

.wrap{max-width:1220px;margin:0 auto;padding:20px 26px 70px;background:#fff;
  min-height:calc(100vh - 46px);}

/* ====================================================== study picker ====== */
/* The only control on the page drawn as a FIELD. It is the one that changes
   every number below it, so it is the one that has to look changeable; every
   other control is a token, a tick or a tab. */
.pick{display:flex;align-items:center;gap:14px;flex-wrap:wrap;
  padding-bottom:15px;border-bottom:1px solid var(--rule);}
.pick .lb{font-size:11.5px;color:var(--faint);}
.pick .fld{min-width:430px;flex:0 1 auto;}
.pick .fld .form-group{margin:0;}
.pick .fld .form-select,.pick .fld select{
  border:1px solid var(--ink) !important;border-radius:0 !important;
  font-size:18px;font-weight:600;letter-spacing:-.015em;color:var(--ink);
  padding:7px 34px 7px 12px;height:auto;box-shadow:none !important;background-color:#fff;
}
.pick .fld .form-select:hover{border-color:var(--rose) !important;}
.pick .fld .form-select:focus{border-color:var(--rose) !important;}
.pick .ct{font-family:var(--mono);font-size:12px;color:var(--faint);}
.pick .arg{margin-left:auto;font-family:var(--mono);font-size:11.5px;color:var(--faint);}
.pick .reps{font-family:var(--mono);font-size:11.5px;color:var(--faint);}

.qline{font-size:15px;color:var(--soft);margin:15px 0 5px;max-width:88ch;line-height:1.55;}
.qline code{font-size:.9em;color:var(--ink2);}
.slice{font-family:var(--mono);font-size:11.5px;color:var(--faint);margin:0 0 6px;}
.slice b{color:var(--soft);font-weight:400;}

/* ======================================================== section rule ==== */
.sec-h{
  display:flex;align-items:baseline;gap:12px;margin:26px 0 0;padding-bottom:6px;
  border-bottom:1px solid var(--ink);font-size:12.5px;font-weight:600;color:var(--ink);
}
.sec-h .n{margin-left:auto;font-family:var(--mono);font-size:11px;color:var(--faint);
  font-weight:400;}
.sec-h .form-group{margin:0;}
.sec-h .form-select{
  border:0;border-bottom:1px solid var(--rule);border-radius:0;background-color:#fff;
  font-family:var(--mono);font-size:12.5px;font-weight:400;color:var(--ink);
  padding:2px 24px 3px 0;height:auto;box-shadow:none !important;
}
.sec-h .form-select:hover{border-bottom-color:var(--rose);color:var(--rose-d);}

/* ============================================ routes and their estimands == */
/* The block that replaced the two rhetorical questions, and it doubles as the
   route selector: the tick, the colour key and the definition are one row
   instead of three places that have to be kept in step. */
.routes .form-group{margin:0;}
.routes .control-label{display:none;}
/* Shiny gives every non-inline input container a fixed 300px width. Left alone
   it squeezed the whole routes table into a third of the page and the estimand
   overprinted its reference. */
.routes .shiny-input-container{width:100%;max-width:none;}
.routes .checkbox{margin:0;border-bottom:1px solid var(--rule2);}
.routes .checkbox:last-child{border-bottom:0;}
.routes .checkbox>label{
  display:grid;grid-template-columns:20px minmax(0,1fr);gap:0 8px;
  align-items:baseline;margin:0;padding:10px 0;width:100%;cursor:pointer;font-weight:400;
}
.routes .checkbox>label>input{margin:0;position:relative;top:2px;accent-color:var(--rose);
  width:14px;height:14px;}
/* Shiny renders choiceNames inside a <span> of its own. That span is the second
   grid cell, and being inline it collapsed the .rw grid inside it onto itself --
   which is what made the estimand and the reference overprint each other. */
.routes .checkbox>label>span{display:block;min-width:0;}
.routes .checkbox>label:hover .wt{color:var(--rose-d);}
.rw{display:grid;grid-template-columns:minmax(0,190px) minmax(0,1fr) minmax(0,30ch);
  gap:0 26px;align-items:baseline;}
.rw .nm{font-family:var(--mono);font-size:12.5px;color:var(--ink);white-space:nowrap;}
.rw .nm i{width:9px;height:9px;display:inline-block;margin-right:9px;
  vertical-align:middle;position:relative;top:-1px;}
.rw .wt{font-size:13px;color:var(--ink);font-weight:600;line-height:1.45;}
.rw .rf{font-size:11.5px;color:var(--faint);text-align:right;line-height:1.5;}
.rw .kd{display:block;color:var(--rose-d);}
@media(max-width:900px){
  .rw{grid-template-columns:1fr;gap:3px;}
  .rw .rf{text-align:left;}
}
.qa{font-size:11.5px;color:var(--faint);padding:8px 0 0;}
.qa a{color:var(--rose-d);text-decoration:none;border-bottom:1px solid var(--rose-m);
  cursor:pointer;}
.qa em{font-style:normal;padding:0 5px;}

/* ============================================== scored against, as tabs === */
/* The choice that decides the answer, so it navigates the page rather than
   sitting in a rail beside it. Shiny renders label>input+span, so the SPAN is
   the tab and the checked state is read off its sibling input. */
.estr{margin-top:22px;border-top:1px solid var(--ink);}
.estr .lb{display:flex;align-items:baseline;gap:14px;padding:7px 0 2px;
  font-size:11.5px;color:var(--faint);}
.estr .lb .as{margin-left:auto;}
.estr .form-group{margin:0;}
.estr .control-label{display:none;}
.estr .shiny-options-group{display:flex;flex-wrap:nowrap;overflow-x:auto;
  border-bottom:1px solid var(--rule);}
.estr .radio-inline{margin:0 !important;padding:0 !important;font-weight:400;
  display:block;flex:0 0 auto;}
.estr .radio-inline input{position:absolute;opacity:0;width:0;height:0;}
.estr .radio-inline span{
  display:block;font-family:var(--mono);font-size:12.5px;color:var(--soft);
  white-space:nowrap;padding:9px 16px 10px;border-bottom:2px solid transparent;
  margin-bottom:-1px;cursor:pointer;
}
.estr .radio-inline:first-child span{padding-left:0;}
.estr .radio-inline span:hover{color:var(--ink);}
.estr .radio-inline input:checked+span{color:var(--ink);border-bottom-color:var(--rose);}
.estr .radio-inline input:focus-visible+span{outline:2px solid var(--rose);outline-offset:-2px;}
.estr .unatt{color:var(--amber);}

/* ================================================================ ledger == */
/* Routes down, reference quantities across. Where a study records more than one
   shared benchmark this table IS the result, and it needs no sentence. */
table.ledger{border-collapse:collapse;font-size:12.5px;width:auto;max-width:100%;
  margin-top:20px;table-layout:auto;}
table.ledger th{font-weight:400;font-size:11.5px;color:var(--soft);
  padding:0 7px 8px 16px;text-align:right;white-space:nowrap;
  border-bottom:1px solid var(--ink);vertical-align:bottom;}
table.ledger th.l,table.ledger td.l{min-width:210px;padding-right:26px;}
table.ledger th.l{text-align:left;padding-left:0;color:var(--faint);}
table.ledger th.on{color:var(--rose-d);font-weight:600;background:var(--rose-l);
  padding-top:6px;}
table.ledger td{font-family:var(--mono);text-align:right;padding:9px 7px 9px 16px;
  color:var(--ink2);white-space:nowrap;border-bottom:1px solid var(--rule2);}
table.ledger td.l{text-align:left;padding-left:0;color:var(--ink);}
table.ledger td.l i{width:9px;height:9px;display:inline-block;margin-right:9px;
  vertical-align:middle;position:relative;top:-1px;}
table.ledger td.on{background:var(--rose-l);}
table.ledger td.ownc,.lednote .mk{box-shadow:inset 0 0 0 1.5px var(--ink2);}
.lednote .mk{display:inline-block;width:13px;height:13px;position:relative;
  top:2px;margin:0 5px 0 1px;}
table.ledger td.bad{color:var(--rose-d);}
table.ledger td .sm{display:block;font-size:11px;color:var(--faint);margin-top:2px;}
table.ledger tbody tr:last-child td{border-bottom:1.4px solid var(--ink);}
.lednote{font-size:12.5px;color:var(--soft);margin:11px 0 0;max-width:92ch;line-height:1.55;}
.lednote b{color:var(--ink);font-weight:600;}
.lednote .fn{color:var(--faint);}

/* ========================================================= plot section === */
/* The controls sit with the figure they redraw. In the previous design they
   were 340px of rail at the top of the page and the plot was two screens down,
   so a control whose object is off screen read as belonging to something else. */
.plot{display:grid;grid-template-columns:200px minmax(0,1fr);gap:0;margin-top:14px;}
.plot .rail{border-right:1px solid var(--rule);padding:14px 22px 16px 0;}
.plot .rail .form-group{margin:0 0 14px;}
.plot .rail label.control-label{font-size:11px;font-weight:600;color:var(--faint);
  margin:0 0 5px;}
.plot .rail .form-select,.plot .rail select{
  border:0;border-bottom:1px solid var(--rule);border-radius:0;background-color:#fff;
  font-family:var(--mono);font-size:12px;color:var(--ink2);padding:2px 20px 4px 0;
  height:auto;box-shadow:none !important;
}
.plot .rail .form-select:hover{border-bottom-color:var(--rose);color:var(--rose-d);}
.plot .rail .checkbox label{font-size:12px;color:var(--soft);font-weight:400;}
.plot .rail .checkbox input{accent-color:var(--rose);margin-right:6px;}
.plot .rail .radio-inline{font-size:12px;color:var(--soft);margin-right:12px;}
.plot .rail .radio-inline input{accent-color:var(--rose);margin-right:5px;}
.plot .rail .hint{font-size:11px;color:var(--faint);line-height:1.5;margin:-8px 0 14px;}
.plot .rail .btn{
  background:#fff;border:1px solid var(--rule);border-radius:0;color:var(--rose-d);
  font-family:var(--mono);font-size:11.5px;padding:5px 10px;align-self:flex-start;
  width:auto;display:inline-block;
}
.plot .rail .btn:hover{border-color:var(--rose);color:var(--rose-d);}
.plot .fig{padding:14px 0 0 26px;min-width:0;}
.key{font-family:var(--mono);font-size:11.5px;color:var(--soft);display:flex;
  flex-wrap:wrap;gap:6px 18px;margin:0 0 12px;align-items:center;}
.key .m i{width:9px;height:9px;display:inline-block;margin-right:7px;
  vertical-align:middle;position:relative;top:-1px;}
.key .band i{width:16px;height:9px;display:inline-block;margin-right:7px;
  background:rgba(44,42,71,.12);vertical-align:middle;position:relative;top:-1px;}
.key .band{color:var(--faint);}
.figcap{font-size:11.5px;color:var(--faint);margin-top:10px;line-height:1.6;
  max-width:96ch;}
.figcap b{color:var(--soft);font-weight:500;}
@media(max-width:820px){
  .plot{grid-template-columns:1fr;}
  .plot .rail{border-right:0;border-bottom:1px solid var(--rule);padding:14px 0;}
  .plot .fig{padding-left:0;}
  .pick .fld{min-width:0;flex:1 1 100%;}
}

/* ============================================================== notices === */
.note{display:flex;gap:11px;padding:11px 14px;margin:14px 0 0;font-size:13px;
  line-height:1.55;border:1px solid var(--rule);background:var(--paper);}
.note .ic{font-family:var(--mono);font-weight:600;flex:none;width:15px;text-align:center;}
.note code{font-size:.9em;}
.note.warn{border-color:var(--rose-m);background:var(--rose-l);color:var(--ink2);}
.note.warn .ic{color:var(--rose-d);}
.note.info .ic{color:var(--faint);}
.lede{font-size:13px;color:var(--soft);line-height:1.6;max-width:96ch;margin:12px 0 14px;}
.lede b{color:var(--ink);font-weight:600;}
.lede code{font-size:.9em;}

/* ============================================================== tables ==== */
table.dataTable{font-size:12.5px !important;border-collapse:collapse !important;}
table.dataTable thead th{
  border-bottom:1px solid var(--ink) !important;border-top:1.4px solid var(--ink) !important;
  font-weight:600;font-size:11.5px;color:var(--soft);padding:7px 10px !important;
}
table.dataTable tbody td{
  border-top:0 !important;border-bottom:1px solid var(--rule2) !important;
  padding:6px 10px !important;font-family:var(--mono);color:var(--ink2);
}
table.dataTable tbody tr.odd,table.dataTable tbody tr.even{background:#fff !important;}
table.dataTable tbody td.mth{font-family:var(--mono);color:var(--ink);}
table.dataTable tbody td.rk{color:var(--faint);}
.dataTables_wrapper .dataTables_filter input,
.dataTables_wrapper .dataTables_length select{
  border:1px solid var(--rule);border-radius:0;font-size:12px;
}
.dataTables_wrapper .dataTables_info,.dataTables_wrapper .dataTables_paginate{
  font-size:11.5px;color:var(--faint);
}

/* =============================================================== footer === */
.foot{
  margin-top:34px;padding-top:11px;border-top:1px solid var(--rule);
  display:flex;flex-wrap:wrap;gap:6px 26px;font-size:11.5px;color:var(--faint);
}
.foot b{color:var(--soft);font-weight:500;}

/* =============================================================== modal ==== */
.modal-content{border-radius:0;border:1px solid var(--ink);}
.modal-header,.modal-footer{border-color:var(--rule);}
.doc h4{font-size:13px;font-weight:600;margin:18px 0 6px;}
.doc h4:first-child{margin-top:0;}
.doc p,.doc li{font-size:13px;color:var(--soft);line-height:1.6;}
.doc ul{padding-left:18px;}
.doc code{font-size:.9em;color:var(--ink2);}

@media (prefers-reduced-motion:no-preference){
  .btn,.helpbtn,.form-select,.estr .radio-inline span,.rw .wt{
    transition:background-color .15s ease,border-color .15s ease,color .15s ease;
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
   not the one being asked for. The <b>ledger</b> under the estimand strip shows
   every route against every reference quantity at once, which is the only place
   the comparison is safe.</p>

   <p>A shared benchmark is sometimes marked <i>no route estimates this</i>. That is
   read off the data, not declared: a route whose own estimand IS a named target
   produces the identical bias column against both, so where no route does, nothing
   in the study is trying to produce that number. It means two different things. In study 08 and in <b>OR to correlation
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

## ---- UI ----------------------------------------------------------------------
## ONE COLUMN, and no sidebar. The previous design put every control in a 340px
## rail, which mixed two kinds of decision that must not be confused: choosing the
## reference quantity changes the RESULT, choosing the horizontal axis changes the
## DRAWING. They are now in different places, and the plot controls sit beside the
## plot rather than two screens above it.
##
## The study and the measure are STATIC selectInputs, not uiOutputs. A select
## whose renderUI reads its own input to set `selected` re-renders on every change,
## which drops focus mid-interaction and can loop; the two selects whose choices
## are known before the session starts are therefore built here, and only the
## things that genuinely depend on the chosen study are rendered.
ui <- page_fluid(
  title = "metaConvert simulation results",
  theme = bs_theme(
    version = 5,
    bg = MC$white, fg = MC$ink, primary = MC$rose,
    base_font = font_google("Work Sans"),
    code_font = font_collection("Cascadia Code", "Cascadia Mono", font_google("IBM Plex Mono")),
    "border-color" = MC$rule
  ),
  padding = 0, gap = 0,
  tags$head(tags$style(HTML(app_css))),

  div(class = "mast",
      span(class = "wm", "metaConvert ", tags$em("simulations")),
      div(class = "rt",
          span(if (is.null(FILES)) "no aggregates found"
               else HTML(paste0(fmt_n(length(unique(FILES$stem))), " studies &middot; ADEMP"))),
          tags$button(id = "howto", class = "helpbtn action-button", type = "button",
                      "how to read this"))),

  div(
    class = "wrap",
    if (is.null(FILES)) {
      div(class = "note warn", span(class = "ic", "!"),
          div(HTML("No aggregate files found in <code>data/aggregated</code>. Run the
                    simulation programme first.")))
    } else {
      tagList(
        ## ---- the question ------------------------------------------------
        ## The one control drawn as a field, because it is the one that changes
        ## every number below it.
        div(class = "pick",
            span(class = "lb", "Study"),
            div(class = "fld",
                selectInput("study", NULL, width = "100%",
                            choices = stats::setNames(
                              unique(FILES$stem),
                              sprintf("%s  %s",
                                      FILES$idx[!duplicated(FILES$stem)],
                                      FILES$label[!duplicated(FILES$stem)])))),
            uiOutput("study_pos", inline = TRUE),
            uiOutput("nrep_ui"),
            uiOutput("study_arg", inline = TRUE)),
        uiOutput("study_intro"),

        ## ---- the routes, and what each estimates -------------------------
        uiOutput("routes_head"),
        div(class = "routes", uiOutput("methods_ui")),
        div(class = "qa", actionLink("all_m", "all"), tags$em("/"),
            actionLink("no_m", "none")),

        ## ---- the reference quantity, as the page navigation --------------
        uiOutput("estimand_strip"),
        uiOutput("ledger_block"),

        ## ---- the evidence ------------------------------------------------
        div(class = "sec-h",
            div(style = "display:flex;align-items:baseline;gap:11px;",
                span("Across the grid:"),
                div(style = "width:272px;",
                    selectInput("metric", NULL, width = "100%",
                                choices = stats::setNames(
                                  names(METRICS), vapply(METRICS, `[[`, "", "label")),
                                selected = "coverage"))),
            uiOutput("metric_note", inline = TRUE)),
        div(class = "plot",
            div(class = "rail",
                uiOutput("xvar_ui"),
                uiOutput("facet_ui"),
                checkboxInput("free_y", "free vertical scale", FALSE),
                div(class = "hint", "A shared scale is the default so panels compare."),
                uiOutput("scale_ui"),
                uiOutput("filters_ui"),
                downloadButton("dl", "rows, CSV", class = "btn")),
            div(class = "fig",
                uiOutput("key_row"),
                uiOutput("plot_slot"),
                uiOutput("plot_cap"))),

        ## ---- the scorecard -----------------------------------------------
        div(class = "sec-h", span("Every measure, every route"),
            span(class = "n", "means over the conditions shown")),
        uiOutput("rank_lede"),
        DTOutput("rank_tbl"),

        ## ---- the rows behind it ------------------------------------------
        div(class = "sec-h", span("The rows behind this view")),
        uiOutput("data_lede"),
        DTOutput("data_tbl"),

        uiOutput("footer")
      )
    }
  )
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
    ## a select with one option is noise; show the fact instead. The hidden select
    ## still has to exist, because every reactive below reads input$nrep.
    if (length(n) == 1)
      tagList(span(class = "reps", paste0("nrep ", fmt_n(n))),
              tags$div(style = "display:none;",
                       selectInput("nrep", NULL, choices = n, selected = n)))
    else
      div(style = "width:120px;",
          selectInput("nrep", NULL, choices = n, selected = n[1], width = "100%"))
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
  ## target_ui is gone: the reference quantity is now the estimand strip, rendered
  ## inside the document by output$estimand_strip. target_order(), default_target()
  ## and unattained_targets() are unchanged and still live at the top of this file.

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
      radioButtons("scale", "Reported scale", choices = lv, selected = lv[1], inline = TRUE),
      div(class = "hint", "One scale at a time: the same formulas, different transforms.")
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
      lapply(rest, function(k) {
        v <- sort(unique(raw()[[k]]))
        selectInput(paste0("f_", k), paste("held at:", k), choices = v, selected = v[1])
      })
    )
  })

  ## input$methods is deliberately NOT req()d: unticking every method is a state
  ## the user can reach with one click, and req() would blank the whole page with
  ## no explanation. An empty selection returns zero rows and the views say so.
  dat <- reactive({
    d <- raw(); req(input$target)
    ms <- input$methods; if (is.null(ms)) ms <- character(0)
    d <- d[d$target == input$target & d$method %in% ms, , drop = FALSE]
    apply_held(d)
  })

  ## The condition count is a property of the design grid, not of how many
  ## methods happen to be ticked, so it is counted before the method filter --
  ## otherwise unticking a method would appear to shrink the simulation.
  n_conditions <- reactive({
    d <- raw(); req(input$target)
    d <- apply_held(d[d$target == input$target, , drop = FALSE])
    cc <- cond_cols()
    if (!length(cc) || !nrow(d)) return(nrow(d))
    ## single-argument [ on a data.frame always returns a data.frame, so drop =
    ## FALSE is not only unnecessary here, it warns on every recalculation
    nrow(unique(d[cc]))
  })

  ## Not every study records every measure -- study 06 imputes a standard error,
  ## so it has no interval and no coverage. With coverage as the default the plot
  ## would silently empty, so the measures the study DOES carry are named instead.
  available_metrics <- reactive({
    d <- dat(); if (!nrow(d)) return(character(0))
    Filter(function(k) k %in% names(d) &&
             any(is.finite(suppressWarnings(as.numeric(d[[k]])))),
           names(METRICS))
  })

  ## ---- the question ----------------------------------------------------------
  ## The study is the only control drawn as a field, because it is the only one
  ## that changes every number on the page. Prev/next and the position counter are
  ## there because the programme is read study by study, and a select alone gives
  ## no sense of where in it you are.
  ## The three facts that sit beside the study select. They are separate outputs
  ## rather than one re-rendered picker, so the select itself is static: a select
  ## whose renderUI reads its own value to set `selected` drops focus on every
  ## change and can loop.
  output$study_pos <- renderUI({
    req(input$study)
    stems <- unique(FILES$stem)
    span(class = "ct", sprintf("%d of %d", match(input$study, stems), length(stems)))
  })

  output$study_arg <- renderUI({
    req(input$study)
    a <- .lookup(STUDY_ARG, input$study)
    if (is.null(a)) return(NULL)
    span(class = "arg", a)
  })

  output$study_intro <- renderUI({
    req(input$study)
    q <- .lookup(STUDY_Q, input$study)
    cc <- cond_cols()
    req(input$target)
    ## The design grid, spelled out -- over the SAME rows the condition count is
    ## taken from, and before the method filter. Counting the levels over the whole
    ## file while counting the conditions over the held subset made study 06 read
    ## "2 conditions, 4 br x 17 rr". The method filter is deliberately not applied:
    ## unticking a route must not appear to shrink the simulation.
    d0 <- apply_held(raw()[raw()$target == input$target, , drop = FALSE])
    lv <- vapply(cc, function(k) length(unique(d0[[k]])), integer(1))
    shown <- setdiff(cc, names(held_values()))
    shown <- shown[lv[shown] > 1]
    grid <- if (length(shown))
      paste(sprintf("%d %s", lv[shown], shown), collapse = " × ") else ""
    held <- held_values()
    tagList(
      if (!is.null(q)) p(class = "qline", HTML(q)) else NULL,
      p(class = "slice",
        HTML(paste0(
          fmt_n(n_conditions()), " conditions",
          if (nzchar(grid)) paste0("  <b>", grid, "</b>") else "",
          if (length(held))
            paste0("  &middot;  held at ",
                   paste(sprintf("%s <b>%s</b>", names(held), unlist(held)), collapse = ", "))
          else "",
          "  &middot;  <b>", fmt_n(input$nrep), "</b> replications each")))
    )
  })

  ## The condition columns not on an axis and not a facet, with their current
  ## value. Used by dat(), by the slice line and by the ledger, so the three can
  ## never describe different subsets.
  held_values <- reactive({
    out <- list()
    for (k in setdiff(cond_cols(), c(input$xvar, input$facet_row, input$facet_col))) {
      v <- input[[paste0("f_", k)]]
      if (!is.null(v)) out[[k]] <- v
    }
    out
  })

  apply_held <- function(d) {
    for (k in names(held_values()))
      d <- d[as.character(d[[k]]) == as.character(held_values()[[k]]), , drop = FALSE]
    d
  }

  ## ---- routes and their estimands --------------------------------------------
  output$routes_head <- renderUI({
    req(input$study)
    h <- .lookup(ROUTE_HEAD, input$study); if (is.null(h)) h <- "Routes"
    n <- length(all_methods())
    div(class = "sec-h", span(h),
        span(class = "n", paste(n, if (n == 1) "route" else "routes")))
  })

  ## Each row carries the tick, the colour key and what the route estimates, so
  ## the three cannot drift apart. The estimand and the reference are DECLARED in
  ## ROUTES (transcribed from the QUESTION header of the matching studies/*.R);
  ## an undeclared route still gets its tick and its swatch.
  output$methods_ui <- renderUI({
    pal <- method_pal(); m <- names(pal)
    checkboxGroupInput(
      "methods", NULL,
      choiceNames = lapply(m, function(k) {
        nt <- route_note(input$study, k)
        tags$span(
          class = "rw",
          tags$span(class = "nm",
                    tags$i(style = paste0("background:", pal[[k]], ";")), k),
          tags$span(class = "wt", if (is.null(nt)) "" else nt$what),
          tags$span(class = "rf",
                    if (is.null(nt)) "" else tagList(
                      if (!identical(nt$kind, "shipped"))
                        tags$span(class = "kd", nt$kind) else NULL,
                      nt$ref)))
      }),
      choiceValues = as.list(m),
      selected = m)
  })

  ## Coverage is the default measure and study 06 does not record it. Rather than
  ## leaving the plot to answer with an error, the selection moves to the first
  ## measure the study DOES carry, and the section rule says which.
  ## Keyed on available_metrics() itself, not on (study, target): when the study
  ## changes, the estimand strip has not re-rendered yet, so dat() is momentarily
  ## empty and an observer keyed on the inputs sees no measures at all and does
  ## nothing. input$metric is isolated so writing to it cannot re-trigger this.
  observe({
    av <- available_metrics()
    m  <- isolate(input$metric)
    if (length(av) && !is.null(m) && !(m %in% av))
      updateSelectInput(session, "metric", selected = av[1])
  })

  ## ---- the reference quantity ------------------------------------------------
  ## What used to be a radio list in the rail, and a tab nobody opened. Rendered
  ## as a tab strip inside the document, because choosing it IS choosing the
  ## question -- and on studies 01, 02, 03, 08 and 09 it reverses the ranking.
  output$estimand_strip <- renderUI({
    ag <- raw()
    tg <- unique(ag$target)
    ord <- target_order(tg)
    unattained <- unattained_targets(ag)
    pretty <- vapply(ord, function(z) {
      base <- if (z == "own") "own estimand"
      else if (z == "own_true_r") "own estimand, true r"
      else if (z == "population") "population parameter"
      else if (z == "sample") "same-sample statistic"
      else sub("_sample$", " (samp.)", sub("_population$", " (pop.)", z))
      base
    }, character(1))
    names(ord) <- pretty
    div(
      class = "estr",
      p(class = "lb", "Scored against"),
      radioButtons("target", NULL, inline = TRUE,
                   choiceNames = lapply(seq_along(ord), function(i)
                     if (ord[[i]] %in% unattained)
                       tags$span(pretty[i], tags$span(class = "unatt", " · no route estimates this"))
                     else tags$span(pretty[i])),
                   choiceValues = as.list(unname(ord)),
                   selected = default_target(tg, input$study))
    )
  })

  ## ---- the ledger ------------------------------------------------------------
  ## Routes down, reference quantities across, over the conditions currently
  ## shown. On the estimand-split studies the mirror in it IS the result, and it
  ## replaces both the three headline cards and the Estimand check tab: the cards
  ## named a leader that this table shows to be a property of the selected column.
  ledger_df <- reactive({
    d <- apply_held(raw())
    ms <- input$methods; if (is.null(ms)) ms <- character(0)
    d <- d[d$method %in% ms, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    tg <- target_order(unique(d$target))
    ms <- intersect(names(method_pal()), unique(d$method))
    if (!length(ms) || !length(tg)) return(NULL)
    num <- function(col, i) if (col %in% names(d)) suppressWarnings(as.numeric(d[[col]][i])) else NA_real_
    cell <- function(m, t) {
      i <- which(d$method == m & d$target == t)
      if (!length(i)) return(list(lead = NA_real_, sub = NA_real_))
      list(lead = mean(num("coverage", i), na.rm = TRUE),
           sub  = mean(abs(num("bias", i)), na.rm = TRUE))
    }
    ## Which shared column each route is already estimating. Computed on the WHOLE
    ## study, not the held slice: it is a structural property of the route, and on a
    ## one- or two-cell slice two bias columns could agree by accident.
    hit <- own_matches(raw())
    ## The `own` column earns its place only where some route own estimand is NOT
    ## on the table. Where every route is marked -- study 01a -- the column is a
    ## duplicate of the marked cells, so it goes.
    if ("own" %in% tg && length(hit) && all(!is.na(hit[ms])) &&
        !identical(input$target, "own"))
      tg <- setdiff(tg, "own")
    grid <- lapply(ms, function(m) lapply(tg, function(t) cell(m, t)))
    have_cov <- any(vapply(unlist(grid, recursive = FALSE),
                           function(z) is.finite(z$lead), logical(1)))
    list(methods = ms, targets = tg, grid = grid, have_cov = have_cov, hit = hit)
  })

  output$ledger_block <- renderUI({
    L <- ledger_df()
    if (is.null(L))
      return(div(class = "note info", span(class = "ic", "i"),
                 div(HTML("No route is selected. Tick at least one above, or use <b>all</b>."))))
    if (length(L$targets) < 2 && !L$have_cov) return(NULL)
    pal <- method_pal()
    lab <- function(z) if (z == "own") HTML("its own<br>estimand")
      else HTML(gsub(" ", "<br>", pretty_target(z)))
    ## Coverage is the lead number because it is the only measure that degrades
    ## under EITHER a wrong estimate or a wrong standard error. Study 06 imputes a
    ## standard error and has no interval, so there the table leads on |bias|.
    head_lab <- if (L$have_cov) HTML("mean coverage <span style=\"font-weight:400\">/ mean |bias|</span>")
                else HTML("mean |bias|")
    body <- lapply(seq_along(L$methods), function(i) {
      m <- L$methods[i]
      tags$tr(
        tags$td(class = "l",
                tags$i(style = paste0("background:", pal[[m]], ";")), m),
        lapply(seq_along(L$targets), function(j) {
          z <- L$grid[[i]][[j]]
          on <- identical(L$targets[j], input$target)
          bad <- is.finite(z$lead) && abs(z$lead - 0.95) > 0.05
          ## The marker says "this column IS the quantity this route was built to
          ## return", which is the reading the note below asks for. On 09a it lands
          ## on every tetrachoric route and not on lipsey_cooper, which is that
          ## study result; on 09b it lands nowhere, which is 09b result.
          own <- identical(unname(L$hit[m]), L$targets[j])
          tags$td(class = paste(c("", if (on) "on", if (own) "ownc",
                                  if (L$have_cov && bad) "bad"), collapse = " "),
                  title = if (own) "the quantity this route was built to return" else NULL,
                  if (L$have_cov) fmt_v(z$lead) else fmt_v(z$sub, 4),
                  if (L$have_cov) tags$span(class = "sm", fmt_v(z$sub, 4)) else NULL)
        }))
    })
    marked <- length(L$hit) && any(!is.na(L$hit[L$methods]))
    tagList(
      tags$table(
        class = "ledger",
        tags$thead(tags$tr(
          tags$th(class = "l", head_lab),
          lapply(L$targets, function(t)
            tags$th(class = if (identical(t, input$target)) "on" else "", lab(t))))),
        tags$tbody(body)),
      p(class = "lednote", HTML(paste0(
        if (marked)
          paste0("<span class=\"mk\"></span> indicates the estimand the approach was ",
                 "built to capture. When a formula is near-exact on its target estimand ",
                 "but not on the selected one, the gap is an <b>estimand mismatch, ",
                 "not an estimation error</b>. ")
        else if ("own" %in% L$targets)
          paste0("No route here estimates any of the shared quantities, so its own ",
                 "estimand is shown as its own column. ")
        else ""
        # ,
        # if ("own" %in% L$targets && marked)
        #   paste0("The <b>own estimand</b> column is there for the routes with no mark: ",
        #          "the quantity they return is not one of the columns. ")
        # else "",
        # "<br>Select the column your review pools, then read the rows</b> ranking ",
        # "routes without fixing the column reports a setting as a finding.",
        # if (any(grepl("sample", L$targets)))
        #   paste0(" <span class=\"fn\">Sample columns score against a statistic recomputed ",
        #          "on each replication, so coverage there is not a calibration statement.</span>")
        # else ""
        ))
        )
    )
  })

  ## ---- the plot section ------------------------------------------------------
  output$metric_note <- renderUI({
    m <- METRICS[[input$metric]]
    span(class = "n",
         if (!is.na(m$ref)) paste0("reference ", m$ref) else m$better)
  })

  output$plot_cap <- renderUI({
    d <- dat(); m <- METRICS[[input$metric]]
    if (!nrow(d)) return(NULL)
    div(class = "figcap", HTML(paste0(
      "<b>", m$label, "</b> against <code>", input$xvar, "</code>",
      if (!identical(input$facet_row, "(none)"))
        paste0(", split by <code>", input$facet_row, "</code>") else "",
      ", scored against the <b>", pretty_target(input$target), "</b>. ",
      fmt_n(input$nrep), " replications per cell",
      if (!is.na(m$mcse) && !is.na(m$ref))
        "; the shaded band is twice the median Monte Carlo standard error of that panel, and a difference inside it is not a result"
      else "", ".")))
  })

  ## ---- key row ---------------------------------------------------------------
  ## A persistent, HTML method key. ggplot draws its own legend in a font and at
  ## a size nothing else on the page uses; this one shares the swatch colours
  ## with the sidebar checkboxes and the ranking table, and survives being
  ## screenshotted with the plot.
  output$key_row <- renderUI({
    d <- dat()
    m <- METRICS[[input$metric]]
    pal <- method_pal()
    ms <- intersect(names(pal), unique(d$method))
    ## if the study never recorded this measure there is no reference line and no
    ## band to key, and claiming otherwise above an empty panel is a small lie
    have <- input$metric %in% available_metrics()
    tagList(
      ## The three warnings that used to sit above the headline cards. They are
      ## about how the numbers may be read, so they belong with the figure rather
      ## than at the top of the page.
      if (.scale_mismatch(input$study, input$scale, input$target))
        div(class = "note warn", span(class = "ic", "!"),
            div(HTML(paste0(
              "The <code>", input$target, "</code> target is on the <b>",
              unname(STUDY_SHARED_SCALE[[input$study]]), "</b> scale and you are viewing the <b>",
              input$scale, "</b> routes. The two are different transforms of the same ",
              "quantity, so the gap below is arithmetic rather than a difference in ",
              "performance. Switch the reported scale back, or score against ",
              "<b>each route&rsquo;s own estimand</b>, which is scale-matched by construction."))))
      else if (!.calibration_ok(input$metric, input$target))
        div(class = "note warn", span(class = "ic", "!"),
            div(HTML(paste0(
              "<b>", m$label, "</b> is not interpretable against the <code>", input$target,
              "</code> target. A nominal 95% interval is not built to cover a quantity ",
              "that is itself recomputed on each replication&rsquo;s own sample. Switch to a ",
              "population parameter to read it, or choose bias or RMSE."))))
      else if (grepl("^own", input$target))
        div(class = "note info", span(class = "ic", "i"),
            div(HTML(paste0(
              "Every route is scored against a <b>different</b> quantity here, so this says ",
              "how accurately each computes its own value &mdash; not which to use. The ledger ",
              "above reads that column beside the shared ones."))))
      else if (!have)
        div(class = "note warn", span(class = "ic", "!"),
            div(HTML(paste0("<b>", m$label, "</b> was not recorded for this study. ",
                            "Measures available: ",
                            paste(vapply(available_metrics(),
                                         function(k) paste0("<b>", METRICS[[k]]$label, "</b>"),
                                         character(1)), collapse = ", "), "."))))
      else NULL,
      div(class = "key",
          lapply(ms, function(k)
            span(class = "m", tags$i(style = paste0("background:", pal[[k]], ";")), k)),
          if (have && !is.na(m$mcse) && !is.na(m$ref))
            span(class = "band", tags$i(), "resolution limit, per panel")
          else NULL)
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
    req(input$target)
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
                            "them; the ledger above shows both columns side by side."))))
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

  ## The Estimand check tab is gone. Its whole content -- how far each route
  ## misses its own estimand versus the population parameter -- is the ledger
  ## above, where it is read as a table rather than inferred from a dumbbell, and
  ## where it cannot be missed by a reader who never opened the third tab.

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
