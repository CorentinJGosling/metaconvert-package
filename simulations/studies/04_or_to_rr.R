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
vanderweele_sqrt_or <- function(or, logor_se) {
  ## RR ~= sqrt(OR); delta method on the log scale gives SE(log RR) = SE(log OR)/2
  logrr <- log(or) / 2
  se    <- logor_se / 2
  data.frame(logrr = logrr, logrr_se = se,
             logrr_ci_lo = logrr - stats::qnorm(.975) * se,
             logrr_ci_up = logrr + stats::qnorm(.975) * se)
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
