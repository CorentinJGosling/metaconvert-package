## =============================================================================
## STUDY 04 -- odds ratio -> risk ratio   [`or_to_rr`]
##
## QUESTION
##   convert_df() converts an OR to an RR via `or_to_rr`, default
##   "metaumbrella_cases", alternatives "metaumbrella_exp", "transpose",
##   "grant", "dipietrantonj". All five are shipped and reachable.
##
##   The methods split into three families, and the split is the point:
##     needs the baseline risk    grant  (and dipietrantonj, via the CI)
##     needs the marginal counts  metaumbrella_cases, metaumbrella_exp
##     needs neither              transpose (OR used as RR), VanderWeele sqrt(OR)
##
##   A meta-analyst who has the baseline risk rarely has it EXACTLY -- it is read
##   off a different table, a different follow-up window, or another trial. So
##   the design crosses the TRUE baseline risk with the GUESSED one. Every
##   published comparison assumes the nuisance parameter is known; that is the
##   gap this study fills.
##
## NOTE ON "grant" AND ZHANG-YU
##   The package's "grant" is RR = OR / (1 - p0 + p0*OR)
##   (R/internal_multiple_formulas.R:56). That is the SAME estimator as
##   Zhang & Yu (1998, JAMA 280:1690-1691); Grant (2014, BMJ 348:f7450) presents
##   it for the same purpose. They are not two methods and must not be listed as
##   such -- cite both for the one estimator.
##
## CANDIDATES ADDED (not in metaConvert)
##   vanderweele  sqrt(OR). VanderWeele TJ (2020) "Optimal approximate
##                conversions of odds ratios and hazard ratios to risk ratios",
##                Biometrics 76(3):746-752 -- proves sqrt(OR) is the MINIMAX
##                conversion when the outcome probability lies in an interval
##                symmetric about 0.5. Requires no baseline risk, which is
##                exactly the situation this study is about, so omitting it
##                would be indefensible.
##                ITS INTERVAL IS NOT SYMMETRIC AND NOT DERIVED FROM ITS SE --
##                the paper forbids transforming the OR interval directly and
##                prescribes dividing/multiplying the transformed limits by the
##                bias ratio instead. See vanderweele_sqrt_or() below, and note
##                that its >= 95% guarantee covers only part of this grid.
##   metafor_conv metafor::conv.2x2() reconstructs the 2x2 from the OR and the
##                margins -- the CRAN-standard alternative to metaumbrella's
##                bespoke grid search. metaConvert already depends on metafor.
## =============================================================================

## ---- shared binary DGP (also used by study 05) ------------------------------
## risk_exp = rr * br must be a probability; conditions violating that are
## dropped by build_grid_04() and REPORTED, not silently NA'd as in the legacy
## code (where rbinom(prob > 1) returned NA for 100% of some cells).
gen_2x2_rr <- function(cond, nrep) {
  n <- cond$n; rr <- cond$rr; br <- cond$br; p <- cond$p_exp
  n_exp <- round(n * p); n_nexp <- n - n_exp

  a <- stats::rbinom(nrep, n_exp,  rr * br)   # cases, exposed
  c_ <- stats::rbinom(nrep, n_nexp, br)       # cases, unexposed
  b <- n_exp - a; d <- n_nexp - c_

  ## Haldane-Anscombe correction on tables with a zero cell. Applied to the
  ## ESTIMATORS' inputs only -- theta_sample is computed from the corrected table
  ## too, since that is the table an analyst would actually have worked from.
  zero <- a == 0 | b == 0 | c_ == 0 | d == 0
  ac <- a + 0.5 * zero; bc <- b + 0.5 * zero
  cc <- c_ + 0.5 * zero; dc <- d + 0.5 * zero

  or_hat     <- (ac * dc) / (bc * cc)
  logor_se   <- sqrt(1 / ac + 1 / bc + 1 / cc + 1 / dc)
  logrr_hat  <- log((ac / (ac + bc)) / (cc / (cc + dc)))

  data.frame(
    a = ac, b = bc, c = cc, d = dc,
    n_exp = ac + bc, n_nexp = cc + dc,
    n_cases = ac + cc, n_controls = bc + dc,
    or = or_hat, logor_se = logor_se,
    or_ci_lo = exp(log(or_hat) - stats::qnorm(.975) * logor_se),
    or_ci_up = exp(log(or_hat) + stats::qnorm(.975) * logor_se),
    corrected = zero,
    ## Targets are on the LOG scale: that is the analysis scale, the scale the
    ## SE lives on, and the scale a meta-analysis pools.
    theta_pop    = log(rr),
    theta_sample = logrr_hat
  )
}

## ---- candidate estimators (not in metaConvert) ------------------------------

## The bias-ratio bound behind VanderWeele's conservative interval. For outcome
## probabilities in [0.5 - v, 0.5 + v] the maximum bias ratio of the square-root
## conversion is 1/sqrt(1 - 4v^2) (2020, p.748 and Appendix). At the paper's
## worked band [w, u] = [0.2, 0.8] -- so v = 0.3 -- that is exactly 1.25. Written
## as the formula rather than the constant so the number's provenance, and its
## dependence on the assumed band, are visible at the call site.
vw_bias_ratio <- function(w = 0.2, u = 0.8) {
  stopifnot(length(w) == 1L, length(u) == 1L, w > 0, u < 1, w < u,
            isTRUE(all.equal(w + u, 1)))   # the bound is derived for a band symmetric about 0.5
  v <- (u - w) / 2
  1 / sqrt(1 - 4 * v^2)
}

## VanderWeele (2020) Biometrics 76(3):746-752.
##
## POINT ESTIMATE: sqrt(OR). Corollary 1 (p.747) proves this is the optimal
## bias-ratio minimax conversion whenever the outcome probabilities lie in an
## interval symmetric about 0.5, so the estimate itself needs no defending.
##
## INTERVAL: NOT a symmetric interval around the estimate. The paper rules that
## out in terms (p.748): "because the square-root transformation ... is an
## approximation, it cannot be applied directly to the confidence interval of the
## odds ratio to obtain 95% coverage over repeated samples of the true risk
## ratio." What it prescribes instead, same page, is to take the square root of
## each OR confidence limit and then widen: "if the square-root transformation of
## the lower limit of the odds ratio confidence interval is divided by 1.25 and
## the square-root transformation of the upper limit ... is multiplied by 1.25,
## then this resulting transformed confidence interval will have at least 95%
## coverage of the true risk ratio provided the outcome probabilities do indeed
## fall between w = 20% and u = 80% ... The coverage of the transformed interval
## will in general be conservative."
##
## Until 2026-08 this function built the forbidden symmetric interval, so the
## coverage column for this arm measured an operation the source paper rejects
## rather than the method it proposes. The point-estimate columns were unaffected.
##
## THE INTERVAL IS DELIBERATELY NOT `logrr +/- qnorm(.975) * logrr_se`. The paper
## gives no variance for the conversion, so `logrr_se` remains the delta-method SE
## of the point estimate (SE(log OR)/2) and is reported for what it is. That is
## safe here because performance() takes coverage and ci_width from the ci_lo/ci_up
## columns and se_ratio from the se column, independently (R/03_performance.R:88-101)
## -- but do NOT try to reconstruct one from the other.
##
## SCOPE, WHICH THE FIX DOES NOT REMOVE: the >= 95% guarantee holds only where
## BOTH outcome probabilities fall in [w, u]. In this study's grid p0 = br and
## p1 = rr * br, so the guaranteed region is only
##     br = 0.30 with rr in {0.75, 1, 2}      (p1 = 0.225 / 0.30 / 0.60)
##     br = 0.50 with rr in {0.5, 0.75, 1}    (p1 = 0.25 / 0.375 / 0.50)
## Every br <= 0.15 cell fails the premise outright (p0 < 0.2). Both br and rr are
## columns of the aggregate, so this split is recoverable post-hoc with no code
## change -- and the coverage column for this arm must be read in the two regions
## separately, exactly as item 4.1's dense-region restriction is.
vanderweele_sqrt_or <- function(or, logor_se, w = 0.2, u = 0.8) {
  logrr <- log(or) / 2
  se    <- logor_se / 2                      # delta method, for the ESTIMATE only

  ## sqrt() of each OR limit is the same thing as +/- z * se/2 on the log scale;
  ## the widening by the bias ratio is what the paper adds and this arm lacked.
  widen <- log(vw_bias_ratio(w, u))
  data.frame(logrr = logrr, logrr_se = se,
             logrr_ci_lo = logrr - stats::qnorm(.975) * se - widen,
             logrr_ci_up = logrr + stats::qnorm(.975) * se + widen)
}

## conv.2x2's margins are n1i = the ROW (exposure) margin and n2i = the COLUMN
## (cases) margin -- NOT the two exposure arms. Verified: for the table
## a/b/c/d = 10/20/30/40 (OR = 2/3), conv.2x2(ori, ni = 100, n1i = 30, n2i = 40)
## returns it exactly, whereas passing n2i = n_nexp = 70 returns 19/11/51/19.
metafor_conv_rr <- function(or, n_exp, n_cases, n_sample) {
  rec <- try(suppressWarnings(metafor::conv.2x2(
    ori = or, ni = n_sample, n1i = n_exp, n2i = n_cases)), silent = TRUE)
  if (inherits(rec, "try-error") || is.null(rec))
    return(data.frame(logrr = NA_real_, logrr_se = NA_real_,
                      logrr_ci_lo = NA_real_, logrr_ci_up = NA_real_)[rep(1, length(or)), ])
  a <- rec$ai; b <- rec$bi; c_ <- rec$ci; d <- rec$di
  logrr <- log((a / (a + b)) / (c_ / (c_ + d)))
  se <- sqrt(1 / a - 1 / (a + b) + 1 / c_ - 1 / (c_ + d))
  data.frame(logrr = logrr, logrr_se = se,
             logrr_ci_lo = logrr - stats::qnorm(.975) * se,
             logrr_ci_up = logrr + stats::qnorm(.975) * se)
}

## ---- estimators --------------------------------------------------------------
estimate_or_to_rr <- function(dat, cond) {
  br_guess <- cond$br_guess
  pieces <- list()

  ## --- package routes reachable from es_from_or_se() -------------------------
  for (m in c("metaumbrella_cases", "metaumbrella_exp", "transpose", "grant")) {
    res <- es_from_or_se(
      or = dat$or, logor_se = dat$logor_se,
      baseline_risk = rep(br_guess, nrow(dat)),      # <- the GUESS, not the truth
      n_exp = dat$n_exp, n_nexp = dat$n_nexp,
      n_cases = dat$n_cases, n_controls = dat$n_controls,
      n_sample = dat$n_exp + dat$n_nexp,
      or_to_rr = m
    )
    p <- as_method(res, "logrr", m); p$route <- "package"
    pieces[[m]] <- p
  }

  ## --- dipietrantonj needs the CI, so it comes through es_from_or_ci() -------
  res_dp <- es_from_or_ci(
    or = dat$or, or_ci_lo = dat$or_ci_lo, or_ci_up = dat$or_ci_up,
    baseline_risk = rep(br_guess, nrow(dat)),
    n_exp = dat$n_exp, n_nexp = dat$n_nexp,
    n_cases = dat$n_cases, n_controls = dat$n_controls,
    n_sample = dat$n_exp + dat$n_nexp,
    or_to_rr = "dipietrantonj"
  )
  p <- as_method(res_dp, "logrr", "dipietrantonj"); p$route <- "package"
  pieces[["dipietrantonj"]] <- p

  ## --- candidates ------------------------------------------------------------
  vw <- vanderweele_sqrt_or(dat$or, dat$logor_se)
  p <- as_method(vw, "logrr", "vanderweele_sqrt_or"); p$route <- "candidate"
  pieces[["vw"]] <- p

  mf <- metafor_conv_rr(dat$or, dat$n_exp, dat$n_cases, dat$n_exp + dat$n_nexp)
  p <- as_method(mf, "logrr", "metafor_conv2x2"); p$route <- "candidate"
  pieces[["mf"]] <- p

  do.call(rbind, pieces)
}

## ---- grid --------------------------------------------------------------------
build_grid_04 <- function() {
  ## The baseline-risk grid must reach 0.5. VanderWeele's sqrt(OR) is minimax for
  ## outcome probabilities in an interval SYMMETRIC ABOUT 0.5, so a grid topping
  ## out at 0.30 tests it only where it is not designed to work and would
  ## misrepresent it. Conversely the rare-outcome end (0.01) is where the naive
  ## "OR = RR" transpose is defensible. Both extremes must be present for the
  ## no-go map to mean anything.
  g <- expand.grid(
    rr       = c(0.25, 0.5, 0.75, 1, 2),
    n        = c(50, 100, 300),
    br       = c(0.01, 0.05, 0.15, 0.30, 0.50),
    br_guess = c(0.01, 0.05, 0.15, 0.30, 0.50),
    p_exp    = c(0.5),
    KEEP.OUT.ATTRS = FALSE
  )
  ok <- g$rr * g$br < 1
  if (any(!ok))
    message(sprintf("  [04_or_to_rr] %d/%d conditions dropped: rr*br >= 1 ",
                    sum(!ok), nrow(g)), "(exposed risk not a probability)")
  g[ok, , drop = FALSE]
}

run_04 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  message("STUDY 04: OR -> RR (or_to_rr), with baseline-risk misspecification")
  run_study("04_or_to_rr", build_grid_04(), gen_2x2_rr, estimate_or_to_rr, nrep, cores)
}
