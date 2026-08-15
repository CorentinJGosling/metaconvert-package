## =============================================================================
## STUDY 07 -- ANCOVA-adjusted means -> SMD (d, g)   [`cov_outcome_r` guess]
##
## QUESTION
##   An ANCOVA reports adjusted means and a RESIDUAL standard deviation. The
##   residual SD is conditional on the covariate, so it is smaller than the
##   marginal SD by a factor sqrt(1 - R^2), and an SMD built on it is on a
##   different scale from every other study in the meta-analysis. metaConvert
##   puts it back on the marginal scale (Cooper eq. 12.24, table 12.3):
##
##       mean_sd = ancova_mean_sd / sqrt(1 - cov_outcome_r^2)
##
##   which needs the covariate-outcome multiple correlation. Trials essentially
##   never print it. The analyst therefore GUESSES, and the package's own V22
##   check exists precisely to say so ("ANCOVA residual SD without
##   cov_outcome_r ... metaConvert will silently underestimate the marginal SD,
##   typically 5-15%"). Leaving the column blank is not an escape: cov_outcome_r
##   propagates as NA, so the row yields no g at all.
##
##   So the design crosses the TRUE covariate-outcome correlation with the
##   GUESSED one, exactly as study 04 crosses `br` with `br_guess`. The truth is
##   used ONLY in the DGP; only the guess is ever handed to metaConvert. This is
##   the study with the least prior art: published ANCOVA-conversion work
##   (Cooper's table 12.3; the Cochrane Handbook 6.5.2.5, which does not flag the
##   back-transform at all) assumes the correlation known.
##
## SCALE -- and the legacy bug this replaces
##   theta_pop is on the MARGINAL scale: the population adjusted mean difference
##   divided by the population marginal SD of Y. The DGP gives Y unit within-arm
##   marginal variance, so theta_pop = delta exactly.
##   theta_sample is that replication's own adjusted mean difference divided by
##   that replication's CRUDE pooled SD -- the marginal-scale sample benchmark.
##
##   archive/legacy_scripts/7_SIM_MEANS_PRE_POST_to_SMD.R:45 divides its
##   benchmark by the RESIDUAL sigma with the marginal back-transform commented
##   out. That puts the benchmark on the conditional scale while every estimator
##   is on the marginal one, so the whole misspecification axis it was built to
##   measure is confounded with a constant factor 1/sqrt(1 - rho^2) -- 2.3x at
##   rho = 0.9. Nothing from that file is usable.
##
## STANDARDISER CHOICES INCLUDED (not input shapes)
##   The package exposes many ANCOVA entry points but only three distinct
##   DENOMINATORS, plus one non-ANCOVA baseline:
##
##   ancova_means_sd      per-arm adjusted SDs, df-pooled by the package, then
##                        back-transformed with the GUESS.
##   ancova_md_sd         the ANCOVA root MSE (df = N-2-n_cov) as a single
##                        pooled residual SD, back-transformed with the GUESS.
##                        Same target as above but a different sigma estimate:
##                        it carries no per-arm leverage inflation and one more
##                        df. The gap between these two is a pure input-shape
##                        effect and is worth measuring for that reason.
##   ancova_pooled_crude  adjusted means over the CRUDE pooled SD. The guess
##                        NEVER enters the point estimate -- only the variance,
##                        through Cooper's (1 - r^2) leading term. This is the
##                        route that isolates "wrong guess moves the estimate"
##                        from "wrong guess moves the interval".
##   crude_means_sd       ignore the ANCOVA entirely (es_from_means_sd on the
##                        unadjusted endpoint means/SDs). Under randomisation
##                        this targets the SAME estimand with no nuisance
##                        parameter at all, so it is the honest fallback the
##                        analyst always has. Its empirical SE should exceed the
##                        ANCOVA routes' by 1/sqrt(1 - rho^2) on the
##                        mean-difference scale; on the SMD scale the ratio is
##                        diluted by the denominator-variability term
##                        theta^2/(2 df), which is the same absolute size for
##                        both and therefore proportionally larger for the
##                        smaller ANCOVA variance. Measured (n = 100,
##                        delta = 0.5, correct guess, nrep = 2000):
##                        1.03 / 1.21 / 2.11 at rho = 0.3 / 0.6 / 0.9, against a
##                        naive 1.05 / 1.25 / 2.29. Getting that dilution right
##                        is the check that the apparatus is sound.
##
## WHAT MISSPECIFICATION COSTS, ANALYTICALLY
##   A route standardised on the back-transformed residual SD converges to
##   theta_pop * sqrt(1 - guess^2) / sqrt(1 - rho^2), so its bias against
##   theta_pop is delta * (sqrt(1-guess^2)/sqrt(1-rho^2) - 1) -- multiplicative,
##   FLAT IN n (the estimand-gap signature of the README), and identically zero
##   at delta = 0, which is why delta = 0 is in the grid: it is the cell where
##   the whole cost moves into the interval and none of it into the estimate.
##   Spot-checked at rho = 0.6, delta = 0.5, n = 100, nrep = 200: guess = 0 gives
##   bias +0.128 (predicted +0.125), guess = 0.9 gives -0.226 (predicted -0.228).
##
## VALIDATION ANCHORS (nrep = 60,000, delta = 0.5, n = 100, 1:1, correct guess)
##   These are the checks that must hold before any misspecification number in
##   this study is quoted; all were run before the file was committed.
##
##   * Against theta_own, g is unbiased for ancova_md_sd, ancova_pooled_crude and
##     crude_means_sd (+0.0001 to +0.0003, i.e. <= 0.34 MCSE). ancova_means_sd
##     sits at -0.0011 (-1.4 MCSE), consistently negative across cells: the
##     per-arm adjusted SDs it pools each carry the leverage inflation
##     sqrt(1 + n_j (xbar_j - xbar)^2 / Sxx), so its standardiser is marginally
##     the largest of the three. That is the input-shape effect, and it is the
##     reason ancova_means_sd and ancova_md_sd are both kept.
##   * d is biased UP by +0.0041 (+4.8 MCSE) where g is not, against the same
##     target. delta*(1/J - 1) = 0.00383 at df = 98. So the d/g split reproduces
##     the textbook small-sample correction and is not an artefact of the DGP.
##   * Coverage against theta_own is 0.951-0.952 and se_ratio 0.986-0.994 for
##     every route: the ANCOVA variance is well calibrated FOR ITS OWN ESTIMAND
##     however wrong the guess. Every coverage failure this study reports is
##     therefore estimand movement, never variance estimation -- the separation
##     the README's dual-target design exists to make. One known contributor to
##     the ~1% shortfall is Cooper's omitted leverage term: E[dx^2/Sxx] is about
##     (1/n_exp + 1/n_nexp)/(n-2), worth ~0.5% on the SE at n = 100.
##   * At delta = 0 with rho = 0.6 and guess = 0.3, every bias is 0.00000 while
##     ancova_pooled_crude reports se_ratio 1.176 and coverage 0.980. The cost of
##     the guess has moved entirely into the interval. That is why delta = 0 and
##     the crude-SD route are both in the design.
##
##   DELIBERATELY EXCLUDED as pure input shapes, all of which reduce
##   algebraically to one of the above before any ES is computed:
##     es_from_ancova_means_se / _ci   -> multiply by sqrt(n) / divide by 2qt,
##                                        then call es_from_ancova_means_sd
##     es_from_ancova_md_se / _ci / _pval -> recover ancova_md_sd, then
##                                        call es_from_ancova_md_sd
##     es_from_ancova_means_sd_pooled_adj -> identical to ancova_md_sd when fed
##                                        the root MSE, and identical to
##                                        ancova_means_sd when fed the df-pooled
##                                        per-arm SDs. It is a shape, not a
##                                        standardiser.
##   Including them would report the same three denominators seven times.
## =============================================================================

## ---- data-generating mechanism ----------------------------------------------
## Within each arm: X ~ N(0,1), Y = delta*T + rho*X + sqrt(1-rho^2)*e.
## So Var(Y | arm) = 1 (marginal), corr(X, Y | arm) = rho, and the residual SD is
## sqrt(1 - rho^2). Unit marginal variance is what makes theta_pop = delta.
##
## X is drawn identically in both arms, i.e. the covariate is balanced in
## expectation. That is the randomised-trial case Cooper's eq. 12.26 variance
## assumes; the package documents (@note on es_from_ancova_means_sd) that it
## omits the covariate-imbalance leverage term. Realised imbalance still varies
## from replication to replication, so the leverage term is present in the DGP's
## own standard errors and absent from the package's -- which is exactly the
## approximation being checked, not something to design away.
##
## The ANCOVA is fitted in closed form rather than by lm() per replication: the
## treatment coefficient, the residual sigma, the LS means at the grand covariate
## mean and their standard errors all agree with lm(y ~ trt + x) / predict(...,
## se.fit = TRUE) to 9+ digits (verified before this file was committed), and the
## closed form vectorises over all nrep replications at once, so generate() costs
## about what the rnorm draws cost. metaConvert is handed ONLY the printed
## summary statistics that come out of it -- adjusted means, their SEs, the root
## MSE, the crude SDs. Nothing sees rho, the raw data, or the fitted slope.
gen_ancova <- function(cond, nrep) {
  n     <- cond$n
  delta <- cond$delta
  rho   <- cond$r_cov          # TRUE covariate-outcome correlation; DGP only
  p     <- cond$p_exp

  n_exp <- round(n * p); n_nexp <- n - n_exp
  i1 <- seq_len(n_exp); i0 <- n_exp + seq_len(n_nexp)
  trt <- c(rep(1, n_exp), rep(0, n_nexp))

  x <- matrix(stats::rnorm(n * nrep), n, nrep)
  y <- delta * trt + rho * x +
    sqrt(1 - rho^2) * matrix(stats::rnorm(n * nrep), n, nrep)

  xb1 <- colMeans(x[i1, , drop = FALSE]); xb0 <- colMeans(x[i0, , drop = FALSE])
  yb1 <- colMeans(y[i1, , drop = FALSE]); yb0 <- colMeans(y[i0, , drop = FALSE])

  ## within-arm centring: everything below is a pooled WITHIN-group sum of
  ## squares, which is what the two-group ANCOVA estimates the slope from.
  xc <- x; yc <- y
  xc[i1, ] <- x[i1, , drop = FALSE] - rep(xb1, each = n_exp)
  xc[i0, ] <- x[i0, , drop = FALSE] - rep(xb0, each = n_nexp)
  yc[i1, ] <- y[i1, , drop = FALSE] - rep(yb1, each = n_exp)
  yc[i0, ] <- y[i0, , drop = FALSE] - rep(yb0, each = n_nexp)
  rm(x, y)

  Sxx <- colSums(xc^2); Sxy <- colSums(xc * yc); Syy <- colSums(yc^2)
  Syy1 <- colSums(yc[i1, , drop = FALSE]^2)
  Syy0 <- colSums(yc[i0, , drop = FALSE]^2)
  rm(xc, yc)

  b  <- Sxy / Sxx                       # pooled within-group slope
  dx <- xb1 - xb0; dy <- yb1 - yb0
  adj_md <- dy - b * dx                 # the ANCOVA treatment coefficient

  df_res    <- n - 3                    # 2 groups + 1 covariate
  sigma_res <- sqrt((Syy - b^2 * Sxx) / df_res)

  ## adjusted (LS) means, evaluated at the grand covariate mean, and their SEs --
  ## the "adjusted mean (SE)" line a trial paper actually prints.
  xbar   <- (n_exp * xb1 + n_nexp * xb0) / n
  adj_m1 <- yb1 - b * (xb1 - xbar)
  adj_m0 <- yb0 - b * (xb0 - xbar)
  se_m1  <- sigma_res * sqrt(1 / n_exp  + (xb1 - xbar)^2 / Sxx)
  se_m0  <- sigma_res * sqrt(1 / n_nexp + (xb0 - xbar)^2 / Sxx)

  data.frame(
    n_exp = n_exp, n_nexp = n_nexp,
    ## --- section 19 inputs: adjusted means + dispersion ---------------------
    ancova_mean_exp = adj_m1, ancova_mean_nexp = adj_m0,
    ## es_from_ancova_means_se converts SE -> SD by SE*sqrt(n) and then calls
    ## es_from_ancova_means_sd, so feeding SE*sqrt(n) here is bit-exact to the
    ## _se entry point -- which is how an extractor reaches this route in
    ## practice, since adjusted SDs are almost never printed and adjusted SEs are.
    ancova_mean_sd_exp = se_m1 * sqrt(n_exp),
    ancova_mean_sd_nexp = se_m0 * sqrt(n_nexp),
    ## --- section 20 inputs: adjusted MD + the root MSE ----------------------
    ancova_md = adj_md, ancova_md_sd = sigma_res,
    ## --- crude endpoint statistics, on the marginal scale -------------------
    mean_exp = yb1, mean_nexp = yb0,
    mean_sd_exp = sqrt(Syy1 / (n_exp - 1)),
    mean_sd_nexp = sqrt(Syy0 / (n_nexp - 1)),
    mean_sd_pooled = sqrt(Syy / (n - 2)),
    ## what the analyst could have computed had they had the raw data -- kept so
    ## the "correct guess" cells can be read as "correct POPULATION value", which
    ## is still not the sample value any single trial would have reported.
    r_cov_sample = Sxy / sqrt(Sxx * Syy),
    ## --- targets ------------------------------------------------------------
    ## MARGINAL scale in both cases. Var(Y | arm) = 1 by construction, so the
    ## population adjusted MD over the population marginal SD is just delta.
    ##
    ## CAVEAT, stated because it will otherwise be misread: on the d scale
    ## theta_sample is ALGEBRAICALLY IDENTICAL to what ancova_pooled_crude
    ## computes (same numerator, same denominator). That row's bias against the
    ## sample target is therefore exactly 0 and its coverage exactly 1 by
    ## construction, not by merit; only its population-target row carries
    ## information. On the g scale the two differ by the factor J.
    theta_pop    = delta,
    theta_sample = adj_md / sqrt(Syy / (n - 2))
  )
}

## ---- estimators: package only, one vectorised call per method ---------------
##
## d and g are on the same scale (both estimate the population SMD; g is the
## small-sample-corrected estimator of it), so a single set of targets serves
## both. They are still run as two separate studies, following 01_smd_to_cor.R,
## so the d-vs-g contrast is a row comparison inside one file rather than a
## measure column smuggled into `method`.
estimate_ancova_scale <- function(measure) {
  force(measure)
  function(dat, cond) {
    nr      <- nrow(dat)
    guess   <- rep(cond$r_cov_guess, nr)     # <- the GUESS, never cond$r_cov
    n_cov   <- rep(1, nr)                    # the DGP has exactly one covariate
    pieces  <- list()

    ## The estimand each back-transforming route converges to, on the marginal
    ## (theta_pop) scale: dividing by sigma_res/sqrt(1-guess^2) instead of by
    ## sigma_Y rescales the population SMD by sqrt(1-guess^2)/sqrt(1-rho^2).
    ## Emitting it as theta_own separates "the formula is implemented correctly"
    ## (bias vs own ~ 0 everywhere, whatever the guess) from "the guess is wrong"
    ## (bias vs pop). A non-zero bias against theta_own is a package bug, not a
    ## misspecification effect -- that is the whole reason for the third target.
    bt_estimand <- cond$delta *
      sqrt(1 - cond$r_cov_guess^2) / sqrt(1 - cond$r_cov^2)

    res <- es_from_ancova_means_sd(
      n_exp = dat$n_exp, n_nexp = dat$n_nexp,
      ancova_mean_exp = dat$ancova_mean_exp,
      ancova_mean_nexp = dat$ancova_mean_nexp,
      ancova_mean_sd_exp = dat$ancova_mean_sd_exp,
      ancova_mean_sd_nexp = dat$ancova_mean_sd_nexp,
      cov_outcome_r = guess, n_cov_ancova = n_cov
    )
    p <- as_method(res, measure, "ancova_means_sd")
    p$route <- "package"; p$theta_own <- bt_estimand
    pieces[["means_sd"]] <- p

    res <- es_from_ancova_md_sd(
      ancova_md = dat$ancova_md, ancova_md_sd = dat$ancova_md_sd,
      cov_outcome_r = guess, n_cov_ancova = n_cov,
      n_exp = dat$n_exp, n_nexp = dat$n_nexp
    )
    p <- as_method(res, measure, "ancova_md_sd")
    p$route <- "package"; p$theta_own <- bt_estimand
    pieces[["md_sd"]] <- p

    ## Crude pooled SD as the standardiser: the guess cannot touch the point
    ## estimate, only Cooper's (1 - r^2) variance term, so this row reads the
    ## misspecification cost off the INTERVAL alone. Its estimand is the marginal
    ## SMD whatever the guess, hence theta_own = theta_pop.
    res <- es_from_ancova_means_sd_pooled_crude(
      ancova_mean_exp = dat$ancova_mean_exp,
      ancova_mean_nexp = dat$ancova_mean_nexp,
      mean_sd_pooled = dat$mean_sd_pooled,
      cov_outcome_r = guess, n_cov_ancova = n_cov,
      n_exp = dat$n_exp, n_nexp = dat$n_nexp
    )
    p <- as_method(res, measure, "ancova_pooled_crude")
    p$route <- "package"; p$theta_own <- cond$delta
    pieces[["pooled_crude"]] <- p

    ## No covariate anywhere: unadjusted endpoint means and SDs. Under
    ## randomisation this is the same estimand, so it is the baseline against
    ## which "was adjusting worth needing a guess?" is answered. It should be
    ## unbiased regardless of rho and its empirical SE should be the ANCOVA
    ## routes' divided by sqrt(1 - rho^2).
    res <- es_from_means_sd(
      mean_exp = dat$mean_exp, mean_sd_exp = dat$mean_sd_exp,
      mean_nexp = dat$mean_nexp, mean_sd_nexp = dat$mean_sd_nexp,
      n_exp = dat$n_exp, n_nexp = dat$n_nexp
    )
    p <- as_method(res, measure, "crude_means_sd")
    p$route <- "package"; p$theta_own <- cond$delta
    pieces[["crude"]] <- p

    ## Nothing here is wrapped in try(). A guess of |r| >= 1 would make the
    ## back-transform undefined and must surface as nonest_rate, not disappear.
    do.call(rbind, pieces)
  }
}

## ---- grid --------------------------------------------------------------------
build_grid_07 <- function() {
  g <- expand.grid(
    delta       = c(0, 0.3, 0.8),
    n           = c(50, 100, 300),
    ## TRUE covariate-outcome correlation. r_cov = 0 is the internal check: the
    ## back-transform is a no-op, sigma_res = sigma_Y, and every route -- however
    ## it is standardised and whatever is guessed at guess = 0 -- must agree.
    r_cov       = c(0, 0.3, 0.5, 0.7),
    ## GUESSED value handed to cov_outcome_r. Spans both sides of every true
    ## value and contains each of them, so the diagonal r_cov_guess == r_cov is
    ## the correctly-specified reference inside the same table. 0 is the analyst
    ## who does no back-transform at all; 0.9 over-corrects by 2.3x on the SD.
    r_cov_guess = c(0, 0.3, 0.5, 0.7, 0.9),
    ## allocation ratio: 1:1, ~1:2, ~1:5.7
    p_exp       = c(0.5, 0.33, 0.15),
    KEEP.OUT.ATTRS = FALSE
  )
  ## The per-arm adjusted SD route pools two arm-level dispersion estimates, so
  ## an arm below ~10 makes the comparison between routes a small-sample artefact
  ## of one route rather than a standardiser difference.
  n_small <- pmin(round(g$n * g$p_exp), g$n - round(g$n * g$p_exp))
  ok <- n_small >= 10
  if (any(!ok))
    message(sprintf("  [07_ancova_to_smd] %d/%d conditions dropped: smaller arm < 10 ",
                    sum(!ok), nrow(g)),
            "(per-arm adjusted SD not comparable to a pooled residual SD)")
  g[ok, , drop = FALSE]
}

run_07 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  grid <- build_grid_07()

  targets <- c(own        = "theta_own",     # per-method, supplied by estimate()
               population = "theta_pop",     # marginal-scale truth: coverage/SE
               sample     = "theta_sample")  # same-sample marginal benchmark

  message("STUDY 07a: ANCOVA adjusted means -> Cohen's d, cov_outcome_r misspecified")
  a <- run_study("07a_ancova_to_smd_d", grid, gen_ancova,
                 estimate_ancova_scale("d"), nrep, cores, targets = targets)

  message("STUDY 07b: ANCOVA adjusted means -> Hedges' g, cov_outcome_r misspecified")
  b <- run_study("07b_ancova_to_smd_g", grid, gen_ancova,
                 estimate_ancova_scale("g"), nrep, cores, targets = targets)

  invisible(list(d = a, g = b))
}
