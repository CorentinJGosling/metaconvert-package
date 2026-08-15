## =============================================================================
## STUDY 08 -- pre/post means -> SMD   [`pre_post_to_smd`]
##
## QUESTION
##   es_from_means_sd_pre_post() (and convert_df() behind it) turns baseline and
##   endpoint means/SDs of two arms into a d/g through `pre_post_to_smd`, default
##   "bonett", alternatives "cooper", "morris_drm", "morris_dz", "morris_dav".
##   All five are shipped and reachable.
##
##   THEY DO NOT TARGET THE SAME QUANTITY. Every one of them divides the SAME
##   numerator (the difference in mean change) by a DIFFERENT standardiser:
##
##     bonett      pooled BASELINE SD                    S_b
##     morris_dav  quadratic mean of baseline & endpoint sqrt((S_b^2+S_p^2)/2)
##     morris_dz   pooled CHANGE-SCORE SD                S_c
##     morris_drm  change SD rescaled by sqrt(2(1-r))    S_c / sqrt(2(1-r))
##     cooper      documented alias of morris_drm
##
##   So a table that scores all five against one reference measures estimand
##   mismatch and prints it as bias. This study emits a PER-METHOD target so the
##   two separate, and a second per-method target so the r-misspecification cost
##   separates from both.
##
## THE FLAW THIS FILE EXISTS TO FIX
##   The legacy `7_SIM_MEANS_PRE_POST_to_SMD.R` never ran, and could not have
##   answered the question if it had: it drew both timepoints with unit variance
##   in both arms, so q = S_b / S_p was identically 1 by construction. Every
##   estimator's deviation factor is exactly 1 at q = 1 -- bonett's is 1/q,
##   d_av's is sqrt(2)/sqrt(q^2+1), d_z's is S_p/S_c, and at q = 1 with a correct
##   r all three collapse to unity. The study asked its question in the single
##   cell where the answer vanishes. Its benchmark was also divided by the ANCOVA
##   RESIDUAL SD with the marginal back-transform commented out, putting it on a
##   different scale from every competitor.
##
##   q IS THEREFORE THE PRIMARY DESIGN FACTOR HERE. papers/pre_post_ipd finds
##   median q = 0.64 across 42 ADHD trials and shows q drives its headline
##   result; the grid spans 0.5 to 1.4 so it covers the observed evidence base
##   AND the region q > 1 that IPD cannot supply -- endpoint SD smaller than
##   baseline SD, which is what a treatment that compresses variance produces.
##   That region is where a simulation earns its keep over a re-analysis.
##
## RELATION TO papers/pre_post_ipd
##   That paper is the empirical half of the same question: 42 real trials, one
##   observed q / r / imbalance each, IPD-ANCOVA benchmark. This is the
##   controlled half: q, r, the analyst's guess at r, the imbalance and n are
##   set by design and crossed, so the deviation surface can be read off
##   directly instead of regressed out of 42 points. Its deviation algebra
##   (scripts/verify_derivations.R) is what the theta_own columns below encode.
## =============================================================================

## ---- data-generating mechanism ----------------------------------------------
## An ANCOVA-parameterised two-arm trial. Endpoint SD is fixed at 1, so every
## quantity below is already on the endpoint-SD scale that theta_pop lives on.
##
##   pre_a  ~ N(mu_pre_a, (q * 1)^2)
##   post_a = delta * 1{a = exp} + beta2 * pre_a + eps,   eps ~ N(0, 1 - rho^2)
##   beta2  = rho / q
##
## This is deliberately NOT the change-score parameterisation used by
## verify_derivations.R (which centres pre at each arm's own mean, making the
## mean CHANGE the fixed parameter). Both give identical deviation algebra --
## they differ only in what "delta" names -- but this one makes the ANCOVA arm
## coefficient itself the grid parameter, so theta_pop = delta exactly rather
## than delta plus an imbalance-dependent offset.
##
## Checks that fall out of the parameterisation:
##   Var(post) = beta2^2 q^2 + (1 - rho^2) = 1          (endpoint SD = 1 in BOTH arms)
##   Cor(pre, post) = rho                                (in both arms)
##   population ANCOVA arm coefficient = delta           (common slope, so the
##                                                        arm term is the gap in
##                                                        arm-specific intercepts)
##   E[post_exp] - E[post_nexp] = delta + rho * imb_std  (endpoint-only MD carries
##                                                        the imbalance, which is
##                                                        the IPD paper's r*imb term)
##
## Both arms share q, rho and the endpoint SD. Unequal arm SDs are the design
## factor of the `pool_sd` study, not this one: under pool_sd = FALSE (the
## default kept here) each arm is standardised by its own SD, under TRUE by a
## cross-arm pooled SD, and the two estimands separate ONLY when the arms' true
## SDs differ. Crossing that factor in here would confound the two arguments.

## Per-arm sufficient statistics, drawn in replication blocks.
##
## Blocked rather than one n x nrep matrix per arm because run_study() executes
## conditions in parallel workers: at n = 200, nrep = 10000 the unblocked form
## holds ~100 MB live per worker, and the harness runs one worker per core.
## Blocking caps that at a few MB while leaving the RNG stream identical (the
## draws are consumed in the same order either way).
arm_sufficient_stats <- function(n, nrep, mu_pre, sigma_pre, beta2, sigma_resid,
                                 arm_effect, block = 2000L) {
  out <- matrix(NA_real_, nrep, 5L,
                dimnames = list(NULL, c("m_pre", "m_post", "ss_pre", "ss_post", "sp")))
  for (s in seq.int(1L, nrep, by = block)) {
    idx <- s:min(s + block - 1L, nrep)
    k   <- length(idx)
    pre  <- matrix(stats::rnorm(n * k, mu_pre, sigma_pre), n, k)
    post <- arm_effect + beta2 * pre +
            matrix(stats::rnorm(n * k, 0, sigma_resid), n, k)
    m_pre  <- colMeans(pre)
    m_post <- colMeans(post)
    pre_c  <- pre  - rep(m_pre,  each = n)
    post_c <- post - rep(m_post, each = n)
    out[idx, "m_pre"]   <- m_pre
    out[idx, "m_post"]  <- m_post
    out[idx, "ss_pre"]  <- colSums(pre_c^2)
    out[idx, "ss_post"] <- colSums(post_c^2)
    out[idx, "sp"]      <- colSums(pre_c * post_c)
  }
  out
}

gen_pre_post <- function(cond, nrep) {
  n     <- cond$n            # per arm
  q     <- cond$q            # S_baseline / S_endpoint, the primary factor
  rho   <- cond$r_true       # TRUE within-arm pre-post correlation
  imb   <- cond$imb_std      # standardised baseline imbalance, (mu_e - mu_n)/S_b
  delta <- cond$delta        # ANCOVA arm coefficient, endpoint-SD units

  sigma_pre   <- q                       # endpoint SD fixed at 1
  beta2       <- rho / q
  sigma_resid <- sqrt(1 - rho^2)

  e <- arm_sufficient_stats(n, nrep, mu_pre = imb * sigma_pre, sigma_pre = sigma_pre,
                            beta2 = beta2, sigma_resid = sigma_resid, arm_effect = delta)
  u <- arm_sufficient_stats(n, nrep, mu_pre = 0, sigma_pre = sigma_pre,
                            beta2 = beta2, sigma_resid = sigma_resid, arm_effect = 0)

  df_pooled <- 2 * n - 2

  ## theta_sample: the ANCOVA estimand recomputed on THIS replication's raw data.
  ##
  ## b_hat is the within-arm pooled slope and the arm coefficient is the
  ## endpoint MD adjusted by it -- algebraically identical to
  ## coef(lm(post ~ group + pre))["group"] by Frisch-Waugh-Lovell, because the
  ## group indicator sweeps out the arm means. Written in closed form because
  ## lm() would be 8.6M fits across the production design (432 conditions x
  ## 10,000 reps x the d and g runs). The identity was checked against
  ## coef(lm(post ~ grp + pre)) on 200 raw-IPD draws: max |difference| = 1.2e-15.
  b_hat     <- (e[, "sp"] + u[, "sp"]) / (e[, "ss_pre"] + u[, "ss_pre"])
  ancova_md <- (e[, "m_post"] - u[, "m_post"]) - b_hat * (e[, "m_pre"] - u[, "m_pre"])

  ## Divided by the CRUDE pooled endpoint SD, NOT by the ANCOVA residual sigma.
  ## The residual SD is conditional on the covariate and is smaller by
  ## sqrt(1 - rho^2) (here 0.60 at rho = 0.8), which is exactly the scale error
  ## that made the legacy benchmark incomparable to every estimator it scored.
  sd_p_hat <- sqrt((e[, "ss_post"] + u[, "ss_post"]) / df_pooled)

  ## Hedges' small-sample factor for the g-scale benchmark, at the pooled df, as
  ## in papers/pre_post_ipd. Written out rather than taken from the package: it
  ## is applied to the REFERENCE quantity, and a benchmark that inherits a bug
  ## from the code under test cannot detect it. (The estimators themselves are
  ## never re-implemented here -- that is the rule this file's comments defend.)
  J_pooled <- exp(lgamma(df_pooled / 2) - 0.5 * log(df_pooled / 2) -
                  lgamma((df_pooled - 1) / 2))

  data.frame(
    n_exp = n, n_nexp = n,
    mean_pre_exp  = e[, "m_pre"],  mean_exp  = e[, "m_post"],
    mean_pre_nexp = u[, "m_pre"],  mean_nexp = u[, "m_post"],
    mean_pre_sd_exp  = sqrt(e[, "ss_pre"]  / (n - 1)),
    mean_sd_exp      = sqrt(e[, "ss_post"] / (n - 1)),
    mean_pre_sd_nexp = sqrt(u[, "ss_pre"]  / (n - 1)),
    mean_sd_nexp     = sqrt(u[, "ss_post"] / (n - 1)),
    ## theta_pop: the ANCOVA estimand on the MARGINAL endpoint SD scale. delta is
    ## already in those units because the DGP fixes the endpoint SD at 1.
    theta_pop      = delta,
    theta_sample   = ancova_md / sd_p_hat,              # d scale (uncorrected)
    theta_sample_g = (ancova_md / sd_p_hat) * J_pooled  # g scale (Hedges LS2)
  )
}

## ---- per-method estimands ----------------------------------------------------
## Population probability limits of each shipped route, from the DGP parameters.
##
## With the endpoint SD fixed at 1 and both arms sharing q and rho:
##   population MD_change = delta + imb * (rho - q)      [the numerator all five share]
##   population MD_post   = delta + imb * rho            [the endpoint route's numerator]
##   S_b  = q                                            [pooled baseline SD]
##   S_av = sqrt((q^2 + 1) / 2)                          [quadratic mean]
##   S_c(r) = sqrt(q^2 + 1 - 2 r q)                      [change SD]
##
## S_c is the one that matters for misspecification. A meta-analyst does not
## observe the change SD; metaConvert RECONSTRUCTS it from the reported SDs and
## the supplied r (R/internal_multiple_formulas.R:1165, :1192). So d_z and d_rm
## inherit the guess in their POINT ESTIMATE, while bonett and d_av do not --
## for those two an r error is purely an SE/coverage error. That asymmetry is
## invisible in any design that fixes q = 1, because there
## sqrt(2(1-r))/S_c(r) = 1 for every r and even d_rm becomes r-free.
##
## `use_r` selects which correlation enters:
##   the GUESS  -> theta_own, the true probability limit of what the analyst ran
##   the TRUTH  -> theta_own_true_r, the same route had r been known
## The three-way gap is the whole point of the study:
##   est - theta_own                 finite-sample estimation error
##   theta_own - theta_own_true_r    cost of guessing r wrong
##   theta_own_true_r - theta_pop    estimand mismatch (the standardiser choice)
pre_post_estimands <- function(cond, use_r) {
  q <- cond$q; rho <- cond$r_true; imb <- cond$imb_std; delta <- cond$delta
  r <- if (use_r == "guess") cond$r_guess else cond$r_true

  md_change <- delta + imb * (rho - q)
  md_post   <- delta + imb * rho
  s_c       <- sqrt(q^2 + 1 - 2 * r * q)   # > 0 for all |r| < 1, q > 0
  s_av      <- sqrt((q^2 + 1) / 2)

  c(bonett            = md_change / q,
    morris_dav        = md_change / s_av,
    morris_dz         = md_change / s_c,
    morris_drm        = md_change * sqrt(2 * (1 - r)) / s_c,
    cooper            = md_change * sqrt(2 * (1 - r)) / s_c,
    endpoint_means_sd = md_post)
}

## ---- estimators: package only, one vectorised call per method ---------------
## Hedges' J is a finite-n bias correction, not an estimand shift (J -> 1), so
## theta_own is identical on the d and g scales; the J gap belongs in the bias
## column, which is precisely what the d-scale run measures.
estimate_pre_post_scale <- function(measure) {
  force(measure)
  function(dat, cond) {
    own      <- pre_post_estimands(cond, "guess")
    own_true <- pre_post_estimands(cond, "true")
    r_g      <- rep(cond$r_guess, nrow(dat))   # the GUESS, fed to both arms
    pieces   <- list()

    ## "cooper" is documented as an alias of "morris_drm" and is remapped at the
    ## first line of both kernels (.pooled_ and .single_group_pre_post_to_smd),
    ## so its row is bit-identical by construction. It is run anyway: the object
    ## of study is the ARGUMENT a user types, and a duplicated row is the
    ## regression test that catches the alias drifting.
    for (m in c("bonett", "cooper", "morris_drm", "morris_dz", "morris_dav")) {
      res <- es_from_means_sd_pre_post(
        mean_pre_exp = dat$mean_pre_exp, mean_exp = dat$mean_exp,
        mean_pre_sd_exp = dat$mean_pre_sd_exp, mean_sd_exp = dat$mean_sd_exp,
        mean_pre_nexp = dat$mean_pre_nexp, mean_nexp = dat$mean_nexp,
        mean_pre_sd_nexp = dat$mean_pre_sd_nexp, mean_sd_nexp = dat$mean_sd_nexp,
        n_exp = dat$n_exp, n_nexp = dat$n_nexp,
        r_pre_post_exp = r_g, r_pre_post_nexp = r_g,
        pre_post_to_smd = m
        ## pool_sd left at its default FALSE (Morris d_ppc1 / Becker 1988: each
        ## arm standardised by its own SD). It is a separate argument with its
        ## own documented rationale and its own estimand gap, which opens only
        ## when the arms' true SDs differ -- a factor this DGP holds equal on
        ## purpose. It deserves its own study, crossed with unequal arm SDs.
      )
      p <- as_method(res, measure, m)
      p$route            <- "package"
      p$scale            <- measure
      p$theta_own        <- own[[m]]
      p$theta_own_true_r <- own_true[[m]]
      pieces[[m]] <- p
    }

    ## REFERENCE ROUTE, not a `pre_post_to_smd` value: the endpoint-only SMD,
    ## which is the estimator whose estimand IS theta_pop's scale (numerator
    ## MD_post, denominator pooled endpoint SD). Included because without it the
    ## table has no line showing what "correctly specified" looks like: at
    ## imb_std = 0 its estimand equals theta_pop exactly, and under imbalance it
    ## picks up the rho * imb_std term the IPD paper derives for it. It also
    ## discards the baseline data entirely, which is the cost side of the
    ## comparison. Reached through a different function, so it is labelled
    ## distinctly and must not be read as a sixth value of the argument.
    res_e <- es_from_means_sd(
      mean_exp = dat$mean_exp, mean_sd_exp = dat$mean_sd_exp,
      mean_nexp = dat$mean_nexp, mean_sd_nexp = dat$mean_sd_nexp,
      n_exp = dat$n_exp, n_nexp = dat$n_nexp
    )
    p <- as_method(res_e, measure, "endpoint_means_sd")
    p$route            <- "package_reference"
    p$scale            <- measure
    p$theta_own        <- own[["endpoint_means_sd"]]
    p$theta_own_true_r <- own_true[["endpoint_means_sd"]]
    pieces[["endpoint"]] <- p

    ## Nothing is dropped or try()-wrapped. The standardisers cannot degenerate
    ## on this grid -- S_c(r)^2 = (q-1)^2 + 2(1-r)q > 0 for |r| < 1 -- so a
    ## non-zero nonest_rate in the output means a package-side guard fired and
    ## must be chased, not filtered.
    do.call(rbind, pieces)
  }
}

## ---- grid --------------------------------------------------------------------
build_grid_08 <- function() {
  g <- expand.grid(
    ## PRIMARY FACTOR. Spans both sides of 1. 0.7 is the IPD paper's subgroup cut
    ## point (its sample's upper-quartile q, and roughly where bonett's predicted
    ## deviation crosses |SMD| = 0.20); 0.5 sits below the observed median 0.64;
    ## 1.0 is the legacy design's single cell, kept so the collapse is visible in
    ## the same table; 1.4 is the variance-compressing region no IPD sample can
    ## supply.
    q        = c(0.5, 0.7, 1.0, 1.4),
    ## TRUE pre-post correlation, crossed with the analyst's GUESS. Guessing this
    ## is routine -- it is almost never reported -- and metaConvert's own default
    ## when the column is blank is a flat 0.8 (R/es_from_PAIRED_MEANS.R:128-131),
    ## the top of this grid. The diagonal r_true == r_guess gives the correctly
    ## specified reference within every other cell of the design.
    r_true   = c(0.2, 0.5, 0.8),
    r_guess  = c(0.2, 0.5, 0.8),
    ## delta = 0 isolates the imbalance term: the estimand-mismatch component is
    ## proportional to the truth, so at delta = 0 only the imbalance term of the
    ## deviation survives and the two contributions are separately identified.
    delta    = c(0, 0.5),
    n        = c(25, 50, 200),   # per arm
    ## Standardised baseline imbalance, the IPD paper's third deviation term.
    ## Built into the population, not left to sampling noise, so its coefficient
    ## is estimable rather than averaged away. 0.3 is beyond that paper's +-0.15
    ## subgroup cut, i.e. a bad-luck randomisation rather than a typical one.
    imb_std  = c(0, 0.3),
    KEEP.OUT.ATTRS = FALSE
  )

  ## Nothing is dropped, and that is a property of the design worth stating
  ## rather than leaving to inference: unlike study 04 (where rr * br >= 1 makes
  ## the exposed risk not a probability) every cell here is attainable, because
  ## the only feasibility constraints are |rho| < 1 and |r_guess| < 1 and both
  ## hold by construction. Any missingness in the output is therefore a package
  ## guard firing, not a design hole.
  n_correct <- sum(g$r_true == g$r_guess)
  message(sprintf(
    "  [08_pre_post_to_smd] %d conditions, 0 dropped (all cells attainable); %d (%.0f%%) have a correctly specified r_pre_post",
    nrow(g), n_correct, 100 * n_correct / nrow(g)))
  g
}

run_08 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  grid <- build_grid_08()

  ## SCALE MATCHING. d and g are not the same scale: g = J*d, and the five routes
  ## are not even corrected on a common df -- under the pool_sd = FALSE kernel
  ## used here, morris_dav's J alone runs on Cousineau's effective df
  ## nu = 2(n-1)/(1+r^2) (R/internal_multiple_formulas.R:1225) rather than n-1, so
  ## on the g scale the analyst's r guess moves d_av's POINT ESTIMATE through J
  ## even though the d-scale estimate is exactly r-invariant. Verified on a fixed
  ## row (n = 25/arm, q = 1, guess 0.2 -> 0.8): d unchanged to machine precision,
  ## g shifts by -3.2e-03 (~1%). Reporting d and g in
  ## one table would average that away, so they are run separately, each against
  ## its own sample-scale benchmark (theta_sample vs theta_sample_g). theta_pop
  ## is delta in both runs: g estimates it, d overstates it by the known factor
  ## 1/J, and the d run is where that cost is quantified for users who pool d.
  message("STUDY 08a: pre/post means -> SMD, d scale")
  a <- run_study("08a_pre_post_to_smd_d", grid, gen_pre_post,
                 estimate_pre_post_scale("d"), nrep, cores,
                 targets = c(own          = "theta_own",
                             own_true_r   = "theta_own_true_r",
                             population   = "theta_pop",
                             sample       = "theta_sample"))

  message("STUDY 08b: pre/post means -> SMD, g scale")
  b <- run_study("08b_pre_post_to_smd_g", grid, gen_pre_post,
                 estimate_pre_post_scale("g"), nrep, cores,
                 targets = c(own          = "theta_own",
                             own_true_r   = "theta_own_true_r",
                             population   = "theta_pop",
                             sample       = "theta_sample_g"))

  invisible(list(d = a, g = b))
}
