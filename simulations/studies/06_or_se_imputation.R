## =============================================================================
## STUDY 06 -- imputing SE(log OR) from the OR and the case/control margins
##
## QUESTION
##   A paper reports an odds ratio and how many cases and controls it had. No
##   confidence interval, no p-value, no 2x2 table. The effect size is poolable
##   only if a standard error can be manufactured, and metaConvert manufactures
##   one: es_from_or(or, n_cases, n_controls) enumerates every 2x2 table
##   compatible with (OR, n_cases, n_controls) and takes the MEAN of the log-OR
##   VARIANCE over that set (.se_from_or, R/internal_multiple_formulas.R:309).
##
##   The legacy study ran that one estimator against the truth and called the
##   result a calibration check. With a single line there is nothing to choose
##   between, so it could not answer the question a user actually has: is this
##   the best thing to do with (OR, n_cases, n_controls), and if not, by how much
##   is the shipped default off? This rebuild adds three comparators (one of them
##   run in two variants, see candidate 2), the first of which is a one-line
##   change to the package.
##
## THE ESTIMAND IS A STANDARD ERROR, WHICH BREAKS THE USUAL ADEMP FRAMING
##   There is no interval and no sampling distribution attached to an imputed SE:
##   .se_from_or returns a single number derived deterministically from three
##   reported quantities. So ci_lo/ci_up are NA by construction, `se` (the
##   estimator's claim about its own variability) is NA for the same reason, and
##   performance() therefore reports coverage = NA and se_ratio = NA. Those cells
##   are not failures, they are inapplicable; the calibration question is settled
##   entirely by bias and rmse. Nothing is being hidden -- an imputation that came
##   with an interval would be a different (and better) method.
##
## WHY EVERYTHING IS SCORED ON THE LOG OF THE SE
##   The quantity being recovered spans an order of magnitude across this grid
##   (measured: SE = 0.13 at n = 1000, br = 0.5 up to SE = 1.21 at n = 100,
##   br = 0.10, p_exp = 0.25 -- a factor of nine between two live cells). On the
##   raw SE scale a single `bias` number mixes "wrong by 60%" at n = 50 with
##   "wrong by 60%" at n = 1000 as if they were different errors, and rmse
##   degenerates into a proxy for n. On the log scale, bias = mean log(SEhat/SEtrue)
##   is exactly the mean PROPORTIONAL error -- a bias of +0.47 says the imputed SE
##   is exp(0.47) = 1.6x too wide, whatever n is -- and over- and under-estimation
##   by the same factor are equidistant, which they are not on the raw scale. This
##   is also the scale on which the error propagates: an inverse-variance weight
##   is a function of log SE, so a constant log-scale bias is a constant weight
##   distortion. It is the same argument that puts studies 04/05 on log RR.
##
## ALL FIVE METHODS TARGET THE SAME SCALAR
##   Unlike study 01 (biserial vs point-biserial) there is no estimand split here:
##   every method is trying to recover sqrt(1/a + 1/b + 1/c + 1/d) for the one
##   table that actually produced the reported OR. They differ only in what
##   information they use to get there. So a single pair of targets is legitimate
##   and no per-method `theta_own` column is emitted.
## =============================================================================

## ---- targets -----------------------------------------------------------------
## log() that returns NA rather than -Inf/NaN/a warning for the non-positive and
## non-finite SEs the degenerate enumerations and zero-cell reconstructions
## produce, so those replications land in nonest_rate instead of poisoning a mean.
.log_pos <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- log(x[ok])
  out
}

## The binary DGP is gen_2x2_rr() from studies/04_or_to_rr.R -- same tables, same
## Haldane-Anscombe handling of zero cells, so the two studies share a
## data-generating mechanism and any difference between them is the estimator.
## Only the targets change: study 04 scores a log RR, this one scores a log SE.
gen_or_se <- function(cond, nrep) {
  dat <- gen_2x2_rr(cond, nrep)

  ## theta_sample: the SE the analyst would have computed had the 2x2 table been
  ## printed. gen_2x2_rr already returns exactly this as logor_se (Woolf's
  ## formula on that replication's own realised, zero-corrected table), so the
  ## imputation is scored against the number it exists to replace.
  dat$theta_sample <- .log_pos(dat$logor_se)

  ## theta_pop: the population SE implied by the true cell probabilities and n,
  ## i.e. Woolf's variance evaluated at the EXPECTED cell counts. This is the
  ## large-sample SE of the log OR estimator, and it is the only fixed reference
  ## in the cell -- theta_sample is itself a random quantity, so bias against it
  ## measures fidelity to the realised table while bias against theta_pop
  ## measures recovery of the sampling variability the meta-analysis is
  ## implicitly assuming when it weights this study.
  n_exp <- round(cond$n * cond$p_exp); n_nexp <- cond$n - n_exp
  p1 <- cond$rr * cond$br; p0 <- cond$br
  var_pop <- 1 / (n_exp * p1) + 1 / (n_exp * (1 - p1)) +
             1 / (n_nexp * p0) + 1 / (n_nexp * (1 - p0))
  dat$theta_pop <- .log_pos(sqrt(var_pop))
  dat
}

## ---- candidate 1: median instead of mean over the compatible tables ----------
##
## The enumeration below is .se_from_or's, transcribed unchanged
## (R/internal_multiple_formulas.R:309-330): split the CASES margin every possible
## way, solve the OR identity for controls_exp,
##   OR = cases_exp * controls_nexp / (cases_nexp * controls_exp)
##     => controls_exp = n_controls / (1 + cases_nexp * OR / cases_exp),
## discard splits whose implied controls_exp falls outside [1, n_controls - 1],
## and reduce the resulting vector of Woolf variances to one number. The package
## reduces with mean(); this reduces with median(). That is the ONLY difference,
## and passing stat = mean here reproduces es_from_or()'s logor_se bit-for-bit --
## verified against the package on the test condition before this file shipped.
##
## WHY THE MEDIAN. 1/a + 1/b + 1/c + 1/d is convex in the split, and it diverges
## at both ends of the enumeration (cases_exp -> 1 and cases_exp -> n_cases - 1
## are legal tables and are included). The mean of a convex function over the
## whole range is therefore dominated by a handful of near-degenerate tables that
## no journal would have published; the median is the obvious robust reduction
## and, as far as we can find, has never been compared against the mean.
##
## The prediction holds and the effect is large: where the true table is far from
## degenerate the package's mean inflates the SE by 1.28x-1.85x (worst at OR = 1,
## br = 0.5, n = 1000) and the median cuts that to 1.04x-1.19x. But the ordering
## REVERSES in thin cells -- at br = 0.10, n = 100, p_exp = 0.25 both are biased
## DOWN and the mean is closer (0.72x vs 0.66x), because there the enumeration's
## near-degenerate members are no longer implausible and the mean's upward pull
## partially offsets. So this is a region map, not a one-line replacement, and
## the grid has to span both regimes for the recommendation to be honest.
##
## The per-row loop is intrinsic: the enumeration length is n_cases - 1 and so
## differs by row. The package does the same thing (apply(dat_or_se, 1, ...)),
## so this is not the forbidden per-replication call into a vectorised es_from_*.
se_or_enumerate <- function(or, n_cases, n_controls, stat = stats::median) {
  vapply(seq_along(or), function(i) {
    o <- or[i]; nca <- n_cases[i]; nco <- n_controls[i]
    ## Same NA guard as the package. Degenerate margins (n_cases < 2) are NOT
    ## guarded, deliberately: they must fail identically in both routes so the
    ## mean/median contrast is not confounded with a different missingness rule.
    if (is.na(o) || is.na(nca) || is.na(nco)) return(NA_real_)
    cases_exp    <- 1:(nca - 1)
    cases_nexp   <- nca - cases_exp
    controls_exp <- round(nco / (1 + cases_nexp * o / cases_exp))
    controls_exp[which(controls_exp < 1 | controls_exp > nco - 1)] <- NA
    controls_nexp <- nco - controls_exp
    v <- 1 / cases_exp + 1 / cases_nexp + 1 / controls_exp + 1 / controls_nexp
    if (all(is.na(v))) return(NA_real_)
    sqrt(stat(v, na.rm = TRUE))
  }, numeric(1))
}

## ---- candidate 2: metafor::conv.2x2 reconstruction ---------------------------
##
## THIS METHOD ANSWERS A SLIGHTLY DIFFERENT QUESTION, and the reason is in
## conv.2x2's argument semantics. Verified empirically (a = 10, b = 20, c = 30,
## d = 40, OR = 2/3): conv.2x2(ori, ni = 100, n1i = 30, n2i = 40) returns that
## table exactly, so n1i is the ROW margin (the EXPOSED group total) and n2i is
## the COLUMN margin (the CASES total). One of each, not both exposure margins.
##
## Our scenario supplies n2i (= n_cases) and ni (= n_cases + n_controls) but NOT
## n1i: the exposure margin is precisely the thing the paper did not report and
## the thing .se_from_or enumerates over. conv.2x2 cannot be run without it, so
## it is run twice, bracketing what reconstruction can achieve:
##
##   _oracle    n1i = the TRUE exposed-group size. Not available to a real
##              analyst, and with both margins plus the OR the table is uniquely
##              determined, so this line is EXACT against theta_sample wherever
##              the integer rounding is exact (measured bias 0.000, rmse 0.000 at
##              n = 300, br = 0.30). It is therefore a CEILING, not a competitor:
##              its distance from zero is what rounding alone costs (visible only
##              in thin cells: bias -0.206 with 25% non-estimable at br = 0.10,
##              n = 100, p_exp = 0.25, where a reconstructed cell rounds to 0),
##              and its distance from _balanced is the price of the missing
##              exposure margin (0.093 on the log scale at p_exp = 0.25).
##   _balanced  n1i = round(N/2), the assumption an analyst who must run
##              conv.2x2 anyway is forced into. This is the feasible version, and
##              it is the one to quote against the package default. It coincides
##              with _oracle by construction at p_exp = 0.5; the two separate
##              only where the design is unbalanced, which is why p_exp is in the
##              grid.
##
## (Note in passing: studies/04_or_to_rr.R:90 calls conv.2x2(n1i = n_exp,
## n2i = n_nexp), i.e. it puts an exposure margin in the cases slot, and so
## reconstructs a table with the wrong column margin. Not fixed here -- flagged.)
metafor_conv_se <- function(or, n_exposed, n_cases, n_total) {
  na_out <- rep(NA_real_, length(or))
  rec <- try(suppressWarnings(metafor::conv.2x2(
    ori = or, ni = round(n_total), n1i = round(n_exposed), n2i = round(n_cases))),
    silent = TRUE)
  if (inherits(rec, "try-error") || is.null(rec)) return(na_out)
  ## Reconstruction returns integer cells; a cell of 0 gives an infinite Woolf
  ## variance, which .log_pos turns into NA and performance() reports as
  ## nonest_rate. No continuity correction is applied because the competing
  ## routes do not get one either at this stage.
  sqrt(1 / rec$ai + 1 / rec$bi + 1 / rec$ci + 1 / rec$di)
}

## ---- estimators ---------------------------------------------------------------
## as_method() wants a frame with <measure>, <measure>_se, <measure>_ci_lo,
## <measure>_ci_up. The estimand here is log SE, and it carries no interval and
## no self-reported precision, so the last three are NA by construction.
se_as_frame <- function(se_hat) {
  data.frame(logse = .log_pos(se_hat), logse_se = NA_real_,
             logse_ci_lo = NA_real_, logse_ci_up = NA_real_)
}

estimate_or_se <- function(dat, cond) {
  pieces <- list()
  n_total <- dat$n_cases + dat$n_controls

  ## --- package: mean of the Woolf variance over the compatible tables --------
  ## Only or / n_cases / n_controls are supplied, because only those three are
  ## available in the scenario and they are the only three .se_from_or reads.
  ##
  ## or_to_rr = "transpose" is a pure SPEED choice and cannot touch the estimand:
  ## logor_se is computed before any OR->RR conversion and is bit-identical under
  ## either setting (checked). The shipped default "metaumbrella_cases" runs a
  ## per-row grid search costing ~23 ms/row, which would put this grid at ~100
  ## CPU-hours at nrep = 10000 to produce a column this study never reads.
  res_pkg <- es_from_or(or = dat$or, n_cases = dat$n_cases,
                        n_controls = dat$n_controls, or_to_rr = "transpose")
  p <- as_method(se_as_frame(res_pkg$logor_se), "logse", "package_mean_var")
  p$route <- "package"
  pieces[["pkg"]] <- p

  ## --- candidate: median of the same variance set ---------------------------
  med <- se_or_enumerate(dat$or, dat$n_cases, dat$n_controls, stat = stats::median)
  p <- as_method(se_as_frame(med), "logse", "median_var")
  p$route <- "candidate"
  pieces[["med"]] <- p

  ## --- candidate: metafor reconstruction, oracle and feasible ---------------
  mo <- metafor_conv_se(dat$or, dat$n_exp, dat$n_cases, n_total)
  p <- as_method(se_as_frame(mo), "logse", "metafor_conv2x2_oracle")
  p$route <- "candidate"
  pieces[["mfo"]] <- p

  mb <- metafor_conv_se(dat$or, round(n_total / 2), dat$n_cases, n_total)
  p <- as_method(se_as_frame(mb), "logse", "metafor_conv2x2_balanced")
  p$route <- "candidate"
  pieces[["mfb"]] <- p

  ## --- candidate: the balanced-table heuristic 4/sqrt(N) --------------------
  ## SE(logOR) = 4/sqrt(N) is the exact Woolf SE of a table with all four cells
  ## equal to N/4. It is not offered as a conversion anywhere, but metaConvert
  ## RELIES on it: the Category D2 cross-row SE-outlier check divides the sqrt(N)
  ## fallback by this constant so that rows whose 2x2 cannot be reconstructed sit
  ## on the same ~1 scale as rows that can (CLAUDE.md, Category D). If the
  ## constant is systematically wrong away from br = 0.5, D2's pool is bimodal
  ## again and the flag is miscalibrated -- so its error surface here is a direct
  ## input to the flag's thresholds, not just a straw man. It uses only N, which
  ## the scenario always supplies.
  p <- as_method(se_as_frame(4 / sqrt(n_total)), "logse", "balanced_4_over_sqrtN")
  p$route <- "candidate"
  pieces[["bal"]] <- p

  do.call(rbind, pieces)
}

## ---- grid ---------------------------------------------------------------------
build_grid_06 <- function() {
  ## br must reach 0.5 (where the enumeration's compatible-table set is widest
  ## and 4/sqrt(N) is exact by construction) AND the rare end (where the true
  ## table is near-degenerate and the enumeration's near-degenerate members stop
  ## being implausible). p_exp is crossed because both the 4/sqrt(N) heuristic
  ## and conv.2x2_balanced assume a 50/50 exposure split, and a grid fixed at
  ## p_exp = 0.5 would test them only where their assumption is true.
  g <- expand.grid(
    or_true = c(0.25, 0.5, 1, 2, 4),
    br      = c(0.02, 0.10, 0.30, 0.50),
    n       = c(50, 100, 300, 1000),
    p_exp   = c(0.5, 0.25),
    KEEP.OUT.ATTRS = FALSE
  )

  ## gen_2x2_rr is parameterised by the RISK ratio, so translate the OR grid into
  ## the rr it expects. Going through the odds guarantees risk_exp is a
  ## probability for every (or_true, br), so unlike build_grid_04 no condition can
  ## be lost to rr * br >= 1.
  odds0 <- g$br / (1 - g$br)
  risk1 <- (g$or_true * odds0) / (1 + g$or_true * odds0)
  g$rr  <- risk1 / g$br

  ## Drop cells where the reported margins would be too thin for the study to
  ## have been published as an OR at all: fewer than 5 expected cases or 5
  ## expected controls. There the enumeration is a handful of tables, Woolf's
  ## formula is itself unreliable, and every method degenerates together, which
  ## tells us nothing about the contrast the study is for.
  n_exp <- round(g$n * g$p_exp); n_nexp <- g$n - n_exp
  e_cases    <- n_exp * risk1 + n_nexp * g$br
  e_controls <- n_exp * (1 - risk1) + n_nexp * (1 - g$br)
  ok <- pmin(e_cases, e_controls) >= 5
  if (any(!ok))
    message(sprintf("  [06_or_se_imputation] %d/%d conditions dropped: ", sum(!ok), nrow(g)),
            "fewer than 5 expected cases or controls (margins too thin for the ",
            "compatible-table enumeration to be meaningful)")
  g[ok, , drop = FALSE]
}

run_06 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  message("STUDY 06: imputing SE(log OR) from the OR and the case/control margins")
  run_study("06_or_se_imputation", build_grid_06(), gen_or_se, estimate_or_se,
            nrep, cores)
}
