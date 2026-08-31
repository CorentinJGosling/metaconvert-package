## =============================================================================
## STUDY 11 -- the reliability standard error
##
##   11a  alpha vs omega, head to head, under six response formats
##   11b  loading shape (does the tau-equivalence rationale survive?)
##   11c  what a misstated per-study SE does to a POOLED reliability generalisation
##   11d  deterministic audit of the claims shipped in ?es_from_omega
##   11e  can a REPORTED total-score mean identify the variance regime?
##
## WHY THIS EXISTS
##   `es_from_cronbach_alpha()` computes a standard error from a closed form in
##   (n, k) alone -- Bonett's (2002) 2k/((k-1)(n-2)) on the ln(1-alpha) scale.
##   `es_from_omega()` REFUSES to, and returns se = NA unless the primary study
##   reported an SE or a CI. The rationale in ?es_from_omega @details is that
##   Bonett's variance descends from the Feldt (1965) / Kristof (1963) F result,
##   which assumes essential tau-equivalence -- the assumption omega exists to
##   drop.
##
##   The applied literature ignores the distinction wholesale. Villacura-Herrera
##   et al. (2025, Work & Stress 39(2) 169-196) ran 13 omega reliability
##   generalisations by feeding omega to metafor's
##   escalc(measure = "ABT", ai = omega, mi = n_items, ni = n) -- Bonett's ALPHA
##   variance, applied to omega. No primary study reported an omega SE; the
##   variance was manufactured from (n, k). metaConvert therefore cannot
##   reproduce any published omega RG meta-analysis.
##
##   The question is whether that refusal is supported. Three outcomes were
##   possible: (a) the closed form is equally wrong for both coefficients, so the
##   asymmetry is unjustified; (b) it is materially worse for omega, so the
##   refusal stands and gets numbers instead of an assumption; (c) it depends on
##   conditions a data extractor can actually observe, so the fallback ships
##   gated.
##
## HOW THE CANDIDATE IS SCORED WITHOUT TOUCHING THE PACKAGE
##   metaConvert ships no omega (n, k) standard error, and the programme's
##   cardinal rule (see the runner header, and study 10's refusal to simulate
##   `icc_to_es = "fisher_tf"`) is that a study must score SHIPPED code, never a
##   private copy of a formula. Both are honoured here by feeding omega_hat to
##   `es_from_cronbach_alpha()`. That is not a workaround: the call is
##   bit-identical to metafor's `measure = "ABT"` up to metaConvert's deliberate
##   sign flip, so it IS the candidate implementation, and it is literally the
##   computation the published RG papers ran. No package edit is required before
##   the study that decides whether to make one.
##
## THE FAIRNESS RULE, which decides the whole question
##   alpha_hat and omega_hat MUST come from the same covariance matrix in the
##   same replication, and be scored on the same surviving replications. A pilot
##   that compared a Pearson-covariance alpha against a WLSMV/polychoric omega
##   measured a gap of -0.107 and would have concluded (b); on the same S the
##   same contrast is -0.010. The estimator swap manufactures the entire effect
##   the study exists to detect. `.rel_coefficients()` computes both from one S
##   and is the single place this can regress.
##
## WHAT IS PRIMARY, AND WHY IT IS NOT se_ratio
##   se_ratio (mod_se / emp_se) is the natural measure and it is NOT ESTIMABLE in
##   the cell that matters most. At k = 3 the sampling distribution of
##   ln(1 - omega_hat) has excess kurtosis around 28, and its sample SD does not
##   converge: measured 0.1448 / 0.1661 / 0.1477 / 0.1616 / 0.1518 / 0.1495 at
##   R = 1k / 2k / 5k / 10k / 20k / 40k, still wandering at forty thousand
##   replications. The implied se_ratio bounces between 0.74 and 0.85 as a pure
##   function of nrep. A study reporting it there would ship an `n_items >= 4`
##   gate, and that gate would be false: coverage at k = 3 is 0.942 (omega) and
##   0.947 (alpha).
##
##   So COVERAGE is primary. It is bounded, always defined, converges at
##   1/sqrt(R) regardless of tail weight, and is the unit Kelley &
##   Pornprasertmanit (2016) and Maydeu-Olivares et al. (2007) report in, which
##   is what makes this study commensurable with the literature.
##
##   Two targets are recorded, and the distinction is the point:
##     `centred`     coverage about the estimator's OWN mean across replications.
##                   Isolates SE calibration exactly. Needs no plim, so it is
##                   immune to the four-way estimand confound below.
##     `population`  coverage about the coefficient's plim in that cell.
##                   SE calibration PLUS bias. Reported separately, never under
##                   the other's name.
##
##   se_ratio is still recorded, and `rel_stability()` reports the excess
##   kurtosis and the robust IQR/1.349 dispersion beside it so a reader can see
##   which cells' se_ratio is meaningless. Where the SD and the robust spread
##   disagree, the SD is the broken one.
##
## THE ESTIMAND CONFOUND
##   On ordinal data there are four different population values: latent alpha,
##   latent omega, the alpha of the categorised composite, and the plim of a
##   one-factor fit to the categorised items. Scoring a coefficient against the
##   wrong one turns bias into apparent under-coverage. `theta_own` is therefore
##   PER METHOD (the study-01 pattern) and is obtained by running the SAME
##   estimator on one large reference sample per condition -- the plim by
##   definition, for a consistent estimator, with no closed form to get wrong.
##
## THE VALIDATION GATE
##   Nothing here is interpretable unless the generator reproduces a published
##   result first. Binary items at 10% endorsement (excess kurtosis 5.11) are
##   Maydeu-Olivares et al. (2007) Table 5, where normal-theory alpha coverage is
##   0.79 and FLAT IN N. `rel_benchmark_gate()` asserts it and must pass before
##   any omega number is quoted. Given that the decision turns on a differential
##   of 0.03 against per-coefficient ratios spanning 0.49 to 1.07, an undetected
##   generator error would be invisible in the differential and decisive in the
##   level.
##
## WHAT 11c IS FOR, AND WHY IT IS THE NOVEL PART
##   The psychometric literature studies SINGLE-SAMPLE coverage. A metaConvert
##   user does not pool single samples; the SE becomes an inverse-variance
##   weight. Writing v_i for the reported variance and c*v_i for the true one,
##   the algebra (verified in 11c) says:
##
##     mu_hat        unbiased for ANY weights -- c never enters w_i = 1/(tau2+v_i),
##                   which is a function of (n_i, k_i) alone.
##     pooled SE     EXACTLY self-correcting under REML. The REML estimating
##                   equation for tau2 is a variance-calibration condition on
##                   precisely the weighted combination that appears in
##                   Var(mu_hat).
##     Hartung-Knapp EXACTLY invariant to a uniform rescaling of the weights.
##     tau2          inflated by (c-1)v. NOT self-correcting.
##     I2            inflated by exactly 1 - 1/c at tau2 = 0. NOT self-correcting.
##     PI            too WIDE, over-covers -- the opposite of the intuitive
##                   prediction, and it partly masks the small-K PI shortfall.
##     fixed effect  takes the full damage: se ratio exactly 1/sqrt(c).
##
##   Four of those are theorems and are used here as VALIDITY GATES, not
##   findings. What is contingent, and what 11c reports, is the MAGNITUDE of the
##   I2 inflation after REML absorption and whether Hartung-Knapp restores pooled
##   coverage on real (heavy-tailed, measured) sampling distributions rather than
##   on the normal errors the algebra assumes.
##
## COST
##   11a 42 cells, 11b 3 cells, at nrep = 1000: about 45k factanal fits at
##   ~1.1 ms including data generation. 11c does ZERO model fits -- it resamples
##   11a's raw frame. The whole study is minutes, not hours.
## =============================================================================


## ---- response formats -------------------------------------------------------
##
## The primary distribution factor is RESPONSE FORMAT, not continuous kurtosis,
## and that is a correction to the obvious design rather than a convenience.
## Real RG data is Likert. Five-point categorisation CAPS the reachable excess
## kurtosis near 2 (exact marginal moments below), so the continuous
## chi-square cell -- where the closed form is worst -- has no empirical referent
## in the data metaConvert users hold. Response format is also the only candidate
## outcome-(c) gate that a data extractor can actually read off a primary study;
## kurtosis is not reported by anyone.
##
## `cont_normal` and `cont_kurt8` are retained as explicitly labelled BOUNDS, not
## as conditions. `bin_10` is the published benchmark cell.
.response_formats <- function() {
  list(
    cont_normal = list(kind = "continuous", df = NA_real_),
    cont_kurt8  = list(kind = "continuous", df = 1.5),
    lik5_sym    = list(kind = "ordinal", p = c(.10, .20, .40, .20, .10)),
    lik5_skew   = list(kind = "ordinal", p = c(.45, .23, .16, .10, .06)),
    lik7_skew   = list(kind = "ordinal", p = c(.30, .22, .17, .12, .09, .06, .04)),
    lik5_sev    = list(kind = "ordinal", p = c(.65, .18, .09, .05, .03)),
    bin_10      = list(kind = "ordinal", p = c(.90, .10))
  )
}

#' Exact marginal skew and excess kurtosis of a categorical response
#'
#' Computed from the category probabilities, not simulated: these are the
#' numbers that justify calling the continuous cells unreachable, so they must
#' not carry Monte Carlo error of their own.
.rel_moments <- function(p) {
  x <- seq_along(p)
  m <- sum(p * x)
  v <- sum(p * (x - m)^2)
  c(skew = sum(p * (x - m)^3) / v^1.5,
    exkurt = sum(p * (x - m)^4) / v^2 - 3)
}

#' Loading vector for k items
#'
#' `congeneric` is the working default. `tau_equiv` and `extreme` exist only for
#' 11b, whose single job is to falsify the shipped tau-equivalence rationale.
.rel_loadings <- function(k, shape = "congeneric") {
  switch(shape,
    tau_equiv  = rep(0.70, k),
    congeneric = seq(0.40, 0.90, length.out = k),
    extreme    = seq(0.25, 0.95, length.out = k),
    stop("unknown loading shape: ", shape)
  )
}

#' One simulated dataset, returned as its covariance matrix
#'
#' Ordinal cells threshold the SAME latent normal deviates, so the latent model
#' is exactly the one-factor model in every cell and only the OBSERVED responses
#' are coarse. That is what makes response format a clean factor rather than a
#' confound with model misspecification.
.rel_draw <- function(n, lam, fmt) {
  k <- length(lam)
  th <- 1 - lam^2
  if (identical(fmt$kind, "continuous") && !is.na(fmt$df)) {
    ## standardised chi-square: skew sqrt(8/df), excess kurtosis 12/df
    rr <- function(m) (stats::rchisq(m, fmt$df) - fmt$df) / sqrt(2 * fmt$df)
  } else {
    rr <- function(m) stats::rnorm(m)
  }
  X <- outer(rr(n), lam) + matrix(rr(n * k), n, k) * rep(sqrt(th), each = n)
  if (identical(fmt$kind, "ordinal")) {
    thr <- stats::qnorm(cumsum(fmt$p)[-length(fmt$p)])
    X <- apply(X, 2, function(z) as.numeric(cut(z, c(-Inf, thr, Inf), labels = FALSE)))
  }
  S <- stats::cov(X)
  ## The item means ride along as an attribute rather than changing the return
  ## type, so every existing caller is untouched. Everything else 11e needs is
  ## already in S: the total score is a plain sum, so var(rowSums(X)) = sum(S)
  ## exactly -- there is no second pass over the data and no approximation.
  attr(S, "item_means") <- colMeans(X)
  S
}

## ---- the floor-effect sweep (11e) -------------------------------------------
##
## 11a establishes that RESPONSE FORMAT drives the variance error, and that the
## NUMBER OF CATEGORIES does not identify it: the three five-point formats span
## c = 0.92 to 1.69 because what matters is where the item mass sits, not how
## many bins it is cut into. Item skew is never reported, so that looked like the
## end of the road.
##
## It is not, and 11e is the test of the way round. A floor effect is visible in
## the TOTAL-SCORE MEAN, which reliability-generalisation extractors routinely
## have: on a C-point scale of k items the mean item response
## (mean_total / k) sits at 1 under a hard floor and at (1 + C)/2 under a
## symmetric distribution. Rescaled to
##
##     floor_position = (mean_total/k - 1) / (C - 1)   in [0, 1]
##
## it is comparable across scales of different length AND different numbers of
## categories -- so if it tracks c, one reported number identifies the regime
## that the category count cannot.
##
## The family below indexes distributions by (C, delta): C-1 cut points placed at
## equally spaced latent quantiles and shifted by delta, so delta = 0 is
## symmetric and delta > 0 pushes mass into the bottom category. This is a SWEEP,
## not a set of hand-picked cells: the relationship has to be estimated over a
## grid rather than interpolated from the five formats 11a happens to use.
.floor_formats <- function(C = c(2, 3, 5, 7), delta = c(0, 0.4, 0.8, 1.2, 1.6)) {
  g <- expand.grid(C = C, delta = delta, KEEP.OUT.ATTRS = FALSE)
  out <- list()
  for (i in seq_len(nrow(g))) {
    CC <- g$C[i]; dd <- g$delta[i]
    thr <- stats::qnorm(seq_len(CC - 1) / CC) + dd
    p <- diff(c(0, stats::pnorm(thr), 1))
    out[[sprintf("q%d_d%02d", CC, round(dd * 10))]] <- list(kind = "ordinal", p = p)
  }
  out
}

.all_formats <- function() c(.response_formats(), .floor_formats())

#' Both coefficients, from ONE covariance matrix
#'
#' THE FAIRNESS RULE LIVES HERE. Any change that computes omega from a different
#' moment matrix than alpha (a polychoric correlation, a WLSMV fit) invalidates
#' the head-to-head and must become a separate labelled method instead.
#'
#' The fitter is `factanal(covmat = )` on the correlation matrix with the
#' solution rescaled to the covariance scale. Pinned against lavaan to < 1e-5 in
#' tests/test-study-11-reliability-se.R; it agrees to ~1e-7 and is 30x faster,
#' which is what makes the bank affordable enough for 11c to be free.
.rel_coefficients <- function(S, n) {
  k <- ncol(S)
  alpha_hat <- (k / (k - 1)) * (1 - sum(diag(S)) / sum(S))

  sdv <- sqrt(diag(S))
  fa <- try(suppressWarnings(
    stats::factanal(covmat = stats::cov2cor(S), factors = 1, n.obs = n)),
    silent = TRUE)
  if (inherits(fa, "try-error")) {
    omega_hat <- NA_real_
  } else {
    lam <- as.numeric(fa$loadings[, 1]) * sdv
    th  <- as.numeric(fa$uniquenesses) * sdv^2
    s <- sum(lam)
    omega_hat <- s^2 / (s^2 + sum(th))
  }
  c(alpha = alpha_hat, omega = omega_hat)
}

#' The plim of each coefficient in one condition
#'
#' Obtained by running the SAME estimator on one large reference sample. For a
#' consistent estimator that IS the plim, and it sidesteps the analytic
#' categorised-covariance machinery whose errors would show up as bias rather
#' than as a crash. Cost is one draw per condition, not per replication.
.rel_plim <- function(lam, fmt, n_ref = 2e5) {
  S <- .rel_draw(n_ref, lam, fmt)
  .rel_coefficients(S, n_ref)
}


## ---- 11a / 11b: generate ----------------------------------------------------
gen_rel <- function(cond, nrep) {
  fmt <- .all_formats()[[as.character(cond$response)]]
  lam <- .rel_loadings(cond$k, as.character(cond$shape))
  n <- cond$n

  out <- vapply(seq_len(nrep), function(i) {
    S <- .rel_draw(n, lam, fmt)
    ## The two total-score moments an RG extractor would transcribe, recorded per
    ## replication so 11e can ask whether they identify the variance regime.
    ## var(rowSums) is exactly sum(S) for a plain sum score.
    c(.rel_coefficients(S, n),
      scale_mean = sum(attr(S, "item_means")),
      scale_sd = sqrt(sum(S)))
  }, numeric(4))

  dat <- data.frame(alpha_hat = out["alpha", ], omega_hat = out["omega", ],
                    scale_mean = out["scale_mean", ], scale_sd = out["scale_sd", ],
                    n_sample = n, n_items = cond$k)

  ## Per-method plims, carried as columns so estimate() can attach the right one
  ## to the right arm. Population alpha != population omega, and on ordinal data
  ## neither equals its latent counterpart.
  pl <- .rel_plim(lam, fmt)
  dat$plim_alpha <- pl["alpha"]
  dat$plim_omega <- pl["omega"]
  dat
}


## ---- 11a / 11b: estimate ----------------------------------------------------
##
## One vectorised call per arm, per the runner's contract. Both arms go through
## the SHIPPED alpha route; the omega arm differs only in which coefficient is
## handed to it, which is exactly the candidate implementation and exactly what
## Villacura-Herrera et al. ran.
##
## `hakstian_whalen` is deliberately absent from the Monte Carlo: its one
## decision-relevant property (a coefficient-DEPENDENT variance, so equal (n, k)
## no longer implies equal weight) is deterministic arithmetic and lives in 11d.
estimate_rel <- function(dat, cond) {
  arms <- list(alpha = dat$alpha_hat, omega = dat$omega_hat)
  pieces <- list()
  for (nm in names(arms)) {
    res <- es_from_cronbach_alpha(
      cronbach_alpha = arms[[nm]],
      n_sample = dat$n_sample, n_items = dat$n_items,
      alpha_to_es = "bonett"
    )
    pc <- as_method(res, "alpha", nm)

    ## PRIMARY target: the estimator's own mean across replications. Coverage
    ## about it isolates SE calibration and needs no plim at all.
    pc$theta_centred <- mean(pc$est[is.finite(pc$est)])
    ## SECONDARY target: this arm's own plim, on the analysis scale. Coverage
    ## about it is SE calibration PLUS bias, and is labelled as such.
    pc$theta_own <- log(1 - dat[[paste0("plim_", nm)]])

    ## carried so 11e can relate the variance error to what an extractor reports
    pc$scale_mean <- dat$scale_mean
    pc$scale_sd <- dat$scale_sd

    pieces[[nm]] <- pc
  }
  do.call(rbind, pieces)
}

## Everything in a raw frame that is NOT a condition column. run_study() binds
## the grid columns on the left and these on the right, so the groupers recover
## the factors by exclusion -- and a new per-replication column that is missing
## from this list would silently become a grouping FACTOR, splitting every cell
## into one group per replication. One list, used by both.
.REL_NON_COND <- c("method", "est", "se", "ci_lo", "ci_up",
                   "theta_centred", "theta_own", "scale_mean", "scale_sd")


build_grid_11a <- function() {
  g <- expand.grid(
    response = names(.response_formats()),
    n        = c(200, 1000),
    ## k = 3 is IN the grid, and its coverage is the evidence that the
    ## `n_items >= 4` gate everyone reaches for is unnecessary. Its se_ratio is
    ## not interpretable (see the header); rel_stability() marks it.
    k        = c(3, 8, 20),
    shape    = "congeneric",
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  g
}

#' 11e: the floor-effect sweep
#'
#' 20 response distributions spanning 2 to 7 categories and symmetric to hard
#' floor, at one k and two n. The point is to ESTIMATE the relationship between
#' what an extractor can report and the variance error, not to interpolate it
#' from the five formats 11a happens to use -- so the grid is a sweep and the
#' verdict is a comparison of two candidate predictors, not a threshold picked
#' to fit.
build_grid_11e <- function() {
  expand.grid(
    response = names(.floor_formats()),
    n        = c(200, 1000),
    k        = 8,
    shape    = "congeneric",
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
}

#' Does a reported total-score mean identify the variance regime?
#'
#' Returns one row per (cell x coefficient) with the variance error `c` beside
#' the two quantities an extractor could supply: the number of response
#' categories, and the standardised position of the mean item response
#'
#'     floor_position = (scale_mean / k - 1) / (C - 1)
#'
#' which is 0 at a hard floor and 0.5 for a symmetric distribution, and is
#' comparable across scales of different length and different C.
#'
#' `sd_over_max` is the second reported quantity, expressed as a fraction of the
#' largest total-score SD arithmetically possible at that mean. The bound is
#' Bhatia-Davis at the ITEM level, var <= (C - m)(m - 1), scaled up by k^2
#' because a floor compresses the items and the sum inherits it: near the floor
#' there is very little dispersion left to have, so a raw SD cannot be read
#' without its mean. Reporting the ratio rather than the SD is what makes the
#' column interpretable across instruments.
rel_floor_detector <- function(raw, min_valid = 30L) {
  fmts <- .all_formats()
  cond_cols <- setdiff(names(raw), .REL_NON_COND)
  grp <- interaction(raw[c(cond_cols, "method")], drop = TRUE, sep = "\r")
  parts <- lapply(split(raw, grp), function(d) {
    ok <- is.finite(d$est) & is.finite(d$se)
    if (sum(ok) < min_valid) return(NULL)
    d <- d[ok, ]
    p <- fmts[[as.character(d$response[1])]]$p
    C <- length(p); k <- d$k[1]

    se_ratio <- mean(d$se) / stats::sd(d$est)
    mean_item <- mean(d$scale_mean) / k
    max_item_var <- (C - mean_item) * (mean_item - 1)   # Bhatia-Davis, per item

    cbind(d[1, c(cond_cols, "method"), drop = FALSE], data.frame(
      n_valid = nrow(d),
      C = C,
      scale_mean = mean(d$scale_mean),
      scale_sd = mean(d$scale_sd),
      mean_item = mean_item,
      floor_position = (mean_item - 1) / (C - 1),
      sd_over_max = mean(d$scale_sd) / (k * sqrt(max(max_item_var, 1e-12))),
      se_ratio = se_ratio,
      c_variance_error = 1 / se_ratio^2,
      coverage_centred = mean(d$ci_lo <= mean(d$est) & d$ci_up >= mean(d$est)),
      row.names = NULL))
  })
  res <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
  rownames(res) <- NULL
  res[order(res$method, res$n, res$floor_position), ]
}

#' The verdict: which reported quantity actually identifies c?
#'
#' 11a already showed that the CATEGORY COUNT does not -- three five-point
#' formats spanned c = 0.92 to 1.69. This puts the candidates side by side on the
#' same sweep, so "collect the response scale" becomes a measured recommendation
#' rather than a plausible one. Spearman, because the relationship is expected to
#' be monotone but not linear, and a flag only ever needs an ordering.
rel_floor_verdict <- function(fd) {
  do.call(rbind, lapply(split(fd, list(fd$method, fd$n), drop = TRUE), function(d) {
    sp <- function(x) suppressWarnings(stats::cor(x, d$c_variance_error,
                                                  method = "spearman"))
    data.frame(method = d$method[1], n = d$n[1], n_cells = nrow(d),
               rho_categories = sp(d$C),
               rho_floor_position = sp(d$floor_position),
               rho_sd_over_max = sp(d$sd_over_max),
               c_min = min(d$c_variance_error), c_max = max(d$c_variance_error),
               row.names = NULL)
  }))
}

#' 11b: three cells whose only job is to falsify the shipped rationale
#'
#' If the tau-equivalence argument in ?es_from_omega were the operative one, the
#' omega arm would degrade as the loadings spread. One contrast settles it; a
#' factorial would be spending compute to restate the same comparison.
build_grid_11b <- function() {
  expand.grid(
    response = "cont_normal",
    n        = 500,
    k        = 8,
    shape    = c("tau_equiv", "congeneric", "extreme"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
}


## =============================================================================
## Diagnostics and gates -- deterministic given a raw frame, so they are plain
## functions rather than run_study() calls (the study-10 / run_99 pattern).
## =============================================================================

#' Is this cell's se_ratio worth reading?
#'
#' The sample SD of a transformed reliability is a fourth-moment quantity: its
#' Monte Carlo error is sqrt((kappa + 2)/(4R)), not the sqrt(1/(2(R-1))) every
#' normal-theory formula quotes. At kappa = 28 and R = 1000 the true relative
#' MCSE is 8.7%, not 2.24%. This records kappa per cell so a reader can see which
#' se_ratio values are noise, and reports the robust IQR/1.349 spread beside the
#' SD -- where they disagree, the SD is the broken one.
rel_stability <- function(raw, kurt_gate = 2) {
  cond_cols <- setdiff(names(raw), .REL_NON_COND)
  grp <- interaction(raw[c(cond_cols, "method")], drop = TRUE, sep = "\r")
  parts <- lapply(split(raw, grp), function(d) {
    e <- d$est[is.finite(d$est)]
    R <- length(e)
    if (R < 10) return(NULL)
    m <- mean(e); s <- stats::sd(e)
    kap <- mean((e - m)^4) / s^4 - 3
    rob <- stats::IQR(e) / 1.349
    cbind(d[1, c(cond_cols, "method"), drop = FALSE],
          data.frame(
            n_valid = R,
            emp_sd = s,
            robust_sd = rob,
            sd_over_robust = s / rob,
            ex_kurtosis = kap,
            ## the honest MCSE of an SD under heavy tails
            sd_mcse_rel = sqrt((kap + 2) / (4 * R)),
            se_ratio_readable = kap <= kurt_gate,
            row.names = NULL))
  })
  res <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
  rownames(res) <- NULL
  res
}

#' The gate that must pass before any omega number is quoted
#'
#' Maydeu-Olivares, Coffman & Hartmann (2007) Table 5: normal-theory coverage for
#' Cronbach's alpha on binary items at 10% endorsement (excess kurtosis 5.11) is
#' about 0.79 and does NOT improve with n. If our generator cannot reproduce
#' that, a DGP bug and a substantive finding are indistinguishable -- and since
#' the decision rests on a differential of 0.03 against levels spanning 0.49 to
#' 1.07, the bug would be invisible in the differential.
rel_benchmark_gate <- function(agg, lo = 0.76, hi = 0.83) {
  d <- agg[agg$response == "bin_10" & agg$k == 8 &
             agg$method == "alpha" & agg$target == "centred", ]
  d <- d[order(d$n), ]
  ## A MISSING benchmark cell must fail, not vanish. An earlier version built the
  ## result frame straight from `d` and so ERRORED on an empty selection -- which
  ## in a runner that only inspects `$pass` is one refactor away from a gate that
  ## is silently absent rather than red.
  if (nrow(d) < 2 || !all(is.finite(d$coverage))) {
    return(data.frame(n = NA_real_, coverage = NA_real_, coverage_mcse = NA_real_,
                      lo = lo, hi = hi, pass = FALSE, row.names = NULL))
  }
  pass <- all(d$coverage >= lo & d$coverage <= hi)
  data.frame(n = d$n, coverage = d$coverage, coverage_mcse = d$coverage_mcse,
             lo = lo, hi = hi, pass = pass, row.names = NULL)
}

#' The head-to-head differential, paired over replication index
#'
#' Pairing is what makes this cheap. alpha_hat and omega_hat are computed from
#' the SAME covariance matrix, so their transformed values correlate 0.988-0.997
#' at k = 8 and the paired Monte Carlo error on the difference is roughly a tenth
#' of the unpaired one. The correlation is REPORTED per cell rather than assumed,
#' because it collapses to about 0.60 at k = 3 where the pairing benefit
#' evaporates.
rel_delta <- function(raw) {
  cond_cols <- setdiff(names(raw), .REL_NON_COND)
  grp <- interaction(raw[cond_cols], drop = TRUE, sep = "\r")
  parts <- lapply(split(raw, grp), function(d) {
    a <- d[d$method == "alpha", ]
    w <- d[d$method == "omega", ]
    if (!nrow(a) || nrow(a) != nrow(w)) return(NULL)
    ok <- is.finite(a$est) & is.finite(w$est) & is.finite(a$se) & is.finite(w$se)
    if (sum(ok) < 30) return(NULL)
    a <- a[ok, ]; w <- w[ok, ]

    cov_c <- function(x) mean(x$ci_lo <= mean(x$est) & x$ci_up >= mean(x$est))
    ra <- mean(a$se) / stats::sd(a$est)
    rw <- mean(w$se) / stats::sd(w$est)
    ca <- cov_c(a); cw <- cov_c(w)
    R <- nrow(a)

    cbind(d[1, cond_cols, drop = FALSE], data.frame(
      n_valid = R,
      cor_paired = stats::cor(a$est, w$est),
      se_ratio_alpha = ra, se_ratio_omega = rw,
      delta_se_ratio = rw - ra,
      cov_centred_alpha = ca, cov_centred_omega = cw,
      delta_coverage = cw - ca,
      ## paired MCSE on the coverage difference: McNemar-style, from the
      ## discordant pairs rather than from two independent proportions
      delta_cov_mcse = {
        ia <- a$ci_lo <= mean(a$est) & a$ci_up >= mean(a$est)
        iw <- w$ci_lo <= mean(w$est) & w$ci_up >= mean(w$est)
        sqrt(max(sum(ia & !iw) + sum(iw & !ia), 1)) / R
      },
      row.names = NULL))
  })
  res <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
  rownames(res) <- NULL
  res
}


## =============================================================================
## 11c -- meta-analytic propagation. ZERO new model fits: it resamples 11a's raw
## frame, so the sampling distributions it pools are the measured, heavy-tailed
## ones rather than the normal errors the algebra assumes.
##
## CONSTRUCTION. One replication is one meta-analysis of K studies. Study j is
## assigned a cell (response, n, k), and one row is drawn from that cell's bank.
## The row contributes its SAMPLING DEVIATION d_j = est - plim, which has exactly
## the right distribution for that cell, and the study's value is
##
##     y_j = theta_j + d_j,     theta_j = mu + u_j,   u_j ~ N(0, tau2)
##
## so mu and tau2 are exact by construction and the k-dependence of the estimand
## is centred away rather than silently counted as heterogeneity. The reported
## variance v_j is the Bonett value the bank row carries; the ORACLE variance is
## that cell's measured emp_se^2.
##
## `confound` is where the only real bias routes live, and they are NOT the
## obvious one. Confounding the SE error with the effect cannot bias mu_hat,
## because the error never enters the weights: w_j = 1/(tau2 + v_j) is a function
## of (n_j, k_j) alone. What CAN bias it is a correlation between the effect and
## the REPORTED variance:
##   sb_tied  longer instruments have higher alpha by Spearman-Brown, AND get
##            more weight because v ~ k/(k-1) decreases in k. The formula creates
##            this confound itself.
##   n_tied   larger samples are more reliable -- an ordinary small-study effect,
##            to which a Bonett pool is maximally exposed because a
##            coefficient-free variance makes the pool an almost pure n-weighted
##            mean.
## =============================================================================

.rel_bank <- function(raw) {
  cond_cols <- intersect(c("response", "n", "k", "shape"), names(raw))
  b <- raw[raw$method %in% c("alpha", "omega") & is.finite(raw$est) & is.finite(raw$se), ]
  b$dev <- b$est - b$theta_own
  key <- do.call(paste, c(unname(as.list(b[cond_cols])), list(b$method), list(sep = "|")))
  split(b[c("dev", "se")], key)
}

.sb_alpha <- function(k, rbar = 0.35) k * rbar / (1 + (k - 1) * rbar)

#' @param M meta-analyses per condition
#' @param arm which coefficient's bank to pool ("omega" is the one at issue)
rel_meta_propagation <- function(raw, M = 2000L, arm = "omega",
                                 cores = 1L, seed_tag = "11c_meta") {
  if (!requireNamespace("metafor", quietly = TRUE))
    stop("metafor is required for 11c")

  bank <- .rel_bank(raw)
  ## the oracle variance of a cell is the measured dispersion of its deviations
  oracle_v <- vapply(bank, function(d) stats::var(d$dev), numeric(1))

  grid <- expand.grid(
    K            = c(10L, 30L),
    tau2         = c(0, 0.01, 0.06),
    response_mix = c("lik5_skew", "bin_10"),
    confound     = c("none", "sb_tied", "n_tied"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )

  one_cond <- function(i) {
    cond <- grid[i, ]
    set.seed(condition_seed(seed_tag, i), kind = "L'Ecuyer-CMRG")
    K <- cond$K; tau2 <- cond$tau2
    ks <- c(3, 8, 20); ns <- c(200, 1000)
    mu0 <- log(1 - 0.85)

    res <- vapply(seq_len(M), function(r) {
      k_j <- sample(ks, K, TRUE)
      n_j <- sample(ns, K, TRUE)
      ## theta_j: mu + u_j, with the confound written into the LOCATION so it
      ## correlates with the reported variance through (n, k)
      theta <- mu0 + stats::rnorm(K, 0, sqrt(tau2))
      if (cond$confound == "sb_tied")
        theta <- log(1 - .sb_alpha(k_j)) + stats::rnorm(K, 0, sqrt(tau2))
      if (cond$confound == "n_tied")
        theta <- mu0 - 0.15 * (n_j == 200) + stats::rnorm(K, 0, sqrt(tau2))
      ## mu_true is the UNWEIGHTED mean of the study-level effect law: the whole
      ## point of a confounded cell is that the weighted estimate need not equal it
      mu_true <- switch(cond$confound,
        none    = mu0,
        sb_tied = mean(log(1 - .sb_alpha(ks))),
        n_tied  = mu0 - 0.15 * 0.5)

      keys <- paste(cond$response_mix, n_j, k_j, "congeneric", arm, sep = "|")
      if (!all(keys %in% names(bank))) return(rep(NA_real_, 9))
      d <- vapply(keys, function(kk) {
        b <- bank[[kk]]; b$dev[sample.int(nrow(b), 1L)] }, numeric(1))
      v_rep <- vapply(keys, function(kk) {
        b <- bank[[kk]]; b$se[sample.int(nrow(b), 1L)]^2 }, numeric(1))
      v_orc <- oracle_v[keys]

      y <- theta + d

      ## NB the padding is load-bearing: an equal-effects fit does not carry the
      ## same slots as a REML one (tau2/I2 can come back NULL), and an unpadded
      ## c() would silently return a SHORT vector, shifting every downstream
      ## index by one rather than failing.
      fit <- function(vv, meth, tst) {
        f <- try(suppressWarnings(
          metafor::rma(yi = y, vi = vv, method = meth, test = tst)), silent = TRUE)
        if (inherits(f, "try-error")) return(rep(NA_real_, 5))
        num1 <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[1]
        c(num1(f$b), num1(f$ci.lb), num1(f$ci.ub), num1(f$tau2), num1(f$I2))
      }
      hk  <- fit(v_rep, "REML", "knha")
      z   <- fit(v_rep, "REML", "z")
      fe  <- fit(v_rep, "FE",   "z")
      orc <- fit(v_orc, "REML", "knha")

      ## true I2 from the oracle variances (Higgins-Thompson typical within-study
      ## variance), so the DEFECT is I2_rep - I2_true rather than I2_rep alone
      wt <- 1 / v_orc
      v_typ <- (K - 1) * sum(wt) / (sum(wt)^2 - sum(wt^2))
      i2_true <- 100 * tau2 / (tau2 + v_typ)

      ## prediction interval on the reported-variance REML+HK fit
      pi_cov <- NA_real_
      f <- try(metafor::rma(yi = y, vi = v_rep, method = "REML", test = "knha"),
               silent = TRUE)
      if (!inherits(f, "try-error")) {
        p <- try(stats::predict(f), silent = TRUE)
        if (!inherits(p, "try-error") && is.finite(p$pi.lb))
          pi_cov <- as.numeric(p$pi.lb <= mu_true && p$pi.ub >= mu_true)
      }

      c(cov_hk  = as.numeric(hk[2]  <= mu_true && hk[3]  >= mu_true),
        cov_z   = as.numeric(z[2]   <= mu_true && z[3]   >= mu_true),
        cov_fe  = as.numeric(fe[2]  <= mu_true && fe[3]  >= mu_true),
        cov_orc = as.numeric(orc[2] <= mu_true && orc[3] >= mu_true),
        bias_hk = hk[1] - mu_true,
        tau2_hk = hk[4],
        i2_rep  = hk[5],
        i2_true = i2_true,
        pi_cov  = pi_cov)
    }, numeric(9))

    m <- rowMeans(res, na.rm = TRUE)
    Mv <- sum(is.finite(res["cov_hk", ]))
    cbind(cond, data.frame(
      M_valid = Mv,
      cov_reml_hk = m["cov_hk"], cov_reml_z = m["cov_z"],
      cov_fe = m["cov_fe"], cov_oracle_hk = m["cov_orc"],
      cov_mcse = sqrt(m["cov_hk"] * (1 - m["cov_hk"]) / max(Mv, 1)),
      bias_mu = m["bias_hk"],
      tau2_est = m["tau2_hk"], tau2_true = tau2,
      i2_rep = m["i2_rep"], i2_true = m["i2_true"],
      i2_inflation = m["i2_rep"] - m["i2_true"],
      pi_coverage = m["pi_cov"],
      row.names = NULL))
  }

  idx <- seq_len(nrow(grid))
  parts <- if (cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, ls(envir = .GlobalEnv, all.names = TRUE),
                            envir = .GlobalEnv)
    parallel::clusterExport(cl, c("bank", "oracle_v", "grid", "M", "arm", "seed_tag"),
                            envir = environment())
    parallel::parLapply(cl, idx, one_cond)
  } else {
    lapply(idx, one_cond)
  }
  res <- do.call(rbind, parts)
  rownames(res) <- NULL
  res
}


## =============================================================================
## 11d -- deterministic audit of the claims shipped in ?es_from_omega @details.
## Seconds, no Monte Carlo. This is the cheapest decisive block in the study: it
## changes shipped documentation before anything else runs.
## =============================================================================

#' The k-sweep that refutes "it diverges at k = 3"
#'
#' sqrt(2k/((k-1)(n-2))) at k = 3 is finite and is the LARGEST value in the
#' sweep; the expression decreases monotonically to sqrt(2/(n-2)). The real k = 3
#' problem is the POINT estimate, not the SE: the one-factor model is
#' just-identified, roughly 18% of samples give an improper solution, and
#' ln(1 - omega) has excess kurtosis around 28. Its INTERVAL is fine.
bonett_se_by_k <- function(k = c(3, 4, 6, 8, 10, 20, 100), n = c(200, 1000)) {
  g <- expand.grid(k = k, n = n, KEEP.OUT.ATTRS = FALSE)
  g$bonett_se <- sqrt(2 * g$k / ((g$k - 1) * (g$n - 2)))
  g$limit_k_inf <- sqrt(2 / (g$n - 2))
  g <- g[order(g$n, g$k), ]
  g$is_max_for_n <- ave(g$bonett_se, g$n, FUN = function(x) x == max(x)) == 1
  rownames(g) <- NULL
  g
}

#' The bifactor identity, and the denominator mismatch nobody states
#'
#' ?es_from_omega asserts the Bonett form is 14-24% narrower than
#' sqrt(2m/((m-1)n)) for a bifactor omega_h. That is arithmetic, not a
#' measurement, and it is reproduced here so the claim is pinned to a CSV rather
#' than to a sentence. Note the two expressions do NOT differ only in k vs m:
#' one carries n - 2 and the other n, which is worth about 0.5% at n = 200 and
#' must not be presented as "the same formula with k replaced by m".
bifactor_se_identity <- function(m = 2:5, k = c(4, 8, 10, 20),
                                 n = c(100, 200, 500, 1500)) {
  g <- expand.grid(m = m, k = k, n = n, KEEP.OUT.ATTRS = FALSE)
  g$bonett_k <- sqrt(2 * g$k / ((g$k - 1) * (g$n - 2)))
  g$bifactor_m <- sqrt(2 * g$m / ((g$m - 1) * g$n))
  g$ratio <- g$bonett_k / g$bifactor_m
  g$pct_narrower <- 100 * (1 - g$ratio)
  ## the same comparison with the denominators matched, isolating the k-vs-m part
  g$bifactor_m_nm2 <- sqrt(2 * g$m / ((g$m - 1) * (g$n - 2)))
  g$denom_effect_pct <- 100 * (1 - g$bifactor_m / g$bifactor_m_nm2)
  rownames(g) <- NULL
  g
}

#' What the Bonett scale's coefficient-free variance costs, and Hakstian-Whalen
#'
#' The Bonett variance is alpha-FREE, so two studies with the same (n, k) get
#' identical weight whether alpha is .85 or .98. Hakstian-Whalen's is not. This
#' is the whole decision-relevant content of the HW scale and it needs no
#' simulation, which is why HW is absent from 11a.
hw_weight_ratio <- function(alpha = c(0.70, 0.85, 0.95, 0.98),
                            n = c(100, 200, 1000), k = c(5, 10, 20)) {
  load_metaconvert()
  g <- expand.grid(alpha = alpha, n = n, k = k, KEEP.OUT.ATTRS = FALSE)
  g$se_bonett <- vapply(seq_len(nrow(g)), function(i)
    es_from_cronbach_alpha(g$alpha[i], g$n[i], g$k[i], alpha_to_es = "bonett")$alpha_se,
    numeric(1))
  g$se_hw <- vapply(seq_len(nrow(g)), function(i)
    es_from_cronbach_alpha(g$alpha[i], g$n[i], g$k[i],
                           alpha_to_es = "hakstian_whalen")$alpha_se, numeric(1))
  ## within each (n, k), how far does the weight move across the alpha range?
  g$w_bonett <- 1 / g$se_bonett^2
  g$w_hw <- 1 / g$se_hw^2
  sp <- function(x) max(x) / min(x)
  g$weight_spread_bonett <- ave(g$w_bonett, g$n, g$k, FUN = sp)
  g$weight_spread_hw <- ave(g$w_hw, g$n, g$k, FUN = sp)
  rownames(g) <- NULL
  g
}

#' Does the tau-equivalence rationale survive its own premise?
#'
#' Under exact tau-equivalence alpha and omega are the same population quantity,
#' so the Bonett variance is the same variance for both by construction. The
#' shipped rationale says the formula stops applying to omega BECAUSE omega drops
#' tau-equivalence. If that were the operative mechanism, the omega arm would
#' degrade as the loadings spread -- which is what 11b tests. This function pins
#' the population half of the argument: the gap between population alpha and
#' population omega as a function of loading spread, so the estimand difference
#' the rationale rests on can be seen next to the variance behaviour it predicts.
tau_equivalence_premise <- function(k = c(4, 8, 20)) {
  shapes <- c("tau_equiv", "congeneric", "extreme")
  out <- lapply(k, function(kk) {
    do.call(rbind, lapply(shapes, function(sh) {
      lam <- .rel_loadings(kk, sh)
      th <- 1 - lam^2
      Sig <- outer(lam, lam) + diag(th)
      a <- (kk / (kk - 1)) * (1 - sum(diag(Sig)) / sum(Sig))
      w <- sum(lam)^2 / (sum(lam)^2 + sum(th))
      data.frame(k = kk, shape = sh,
                 loading_min = min(lam), loading_max = max(lam),
                 pop_alpha = a, pop_omega = w, omega_minus_alpha = w - a,
                 gap_in_bonett_se = (log(1 - a) - log(1 - w)) /
                   sqrt(2 * kk / ((kk - 1) * 198)))
    }))
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}


## =============================================================================
## Runner
## =============================================================================
run_11 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  dir.create(dir_agg(), recursive = TRUE, showWarnings = FALSE)

  ## ---- 11d first: seconds, and it settles shipped documentation claims -------
  message("STUDY 11d: deterministic audit of the ?es_from_omega @details claims")
  d_k  <- bonett_se_by_k()
  d_bf <- bifactor_se_identity()
  d_hw <- hw_weight_ratio()
  d_te <- tau_equivalence_premise()
  for (nmv in list(c("11d_bonett_se_by_k", "d_k"), c("11d_bifactor_identity", "d_bf"),
                   c("11d_hw_weight_ratio", "d_hw"), c("11d_tau_equivalence_premise", "d_te"))) {
    f <- dir_agg(sprintf("%s_reps1.csv", nmv[1]))
    utils::write.csv(get(nmv[2]), f, row.names = FALSE)
    write_provenance(basename(f))
  }
  message(sprintf("  k=3 Bonett SE at n=200: %.5f (largest in sweep: %s)",
                  d_k$bonett_se[d_k$k == 3 & d_k$n == 200][1],
                  d_k$is_max_for_n[d_k$k == 3 & d_k$n == 200][1]))

  ## ---- 11a: the head-to-head -------------------------------------------------
  message("STUDY 11a: alpha vs omega under six response formats")
  a <- run_study("11a_alpha_vs_omega", build_grid_11a(), gen_rel, estimate_rel,
                 nrep, cores,
                 targets = c(centred = "theta_centred", population = "theta_own"))

  ## THE GATE. Nothing downstream is interpretable if this fails.
  gate <- rel_benchmark_gate(a$agg)
  f <- dir_agg("11a_benchmark_gate_reps1.csv")
  utils::write.csv(gate, f, row.names = FALSE)
  write_provenance(basename(f))
  if (!isTRUE(gate$pass[1])) {
    warning("BENCHMARK GATE FAILED: alpha coverage on binary-10% items is ",
            paste(round(gate$coverage, 3), collapse = "/"),
            " outside [", gate$lo[1], ", ", gate$hi[1], "] ",
            "(Maydeu-Olivares et al. 2007 Table 5 reports ~0.79, flat in n). ",
            "The generator is suspect; do not interpret any omega number.")
  } else {
    message(sprintf("  benchmark gate PASSED: coverage %s (target ~0.79, flat in n)",
                    paste(round(gate$coverage, 3), collapse = " / ")))
  }

  ## paired differential + the se_ratio readability diagnostics
  dl <- rel_delta(a$raw)
  st <- rel_stability(a$raw)
  for (nmv in list(c("11a_delta_paired", "dl"), c("11a_stability", "st"))) {
    f <- dir_agg(sprintf("%s_reps%d.csv", nmv[1], nrep))
    utils::write.csv(get(nmv[2]), f, row.names = FALSE)
    write_provenance(basename(f))
  }
  message(sprintf("  |delta_se_ratio| max = %.3f ; |delta_coverage| max = %.3f",
                  max(abs(dl$delta_se_ratio), na.rm = TRUE),
                  max(abs(dl$delta_coverage), na.rm = TRUE)))

  ## ---- 11b: loading shape ----------------------------------------------------
  message("STUDY 11b: does the tau-equivalence rationale predict anything?")
  b <- run_study("11b_loading_shape", build_grid_11b(), gen_rel, estimate_rel,
                 nrep, cores,
                 targets = c(centred = "theta_centred", population = "theta_own"))

  ## ---- 11e: can a reported total-score mean identify the regime? -------------
  message("STUDY 11e: floor-effect sweep (does the scale mean identify c?)")
  e <- run_study("11e_floor_sweep", build_grid_11e(), gen_rel, estimate_rel,
                 nrep, cores,
                 targets = c(centred = "theta_centred", population = "theta_own"))
  fd <- rel_floor_detector(e$raw)
  vd <- rel_floor_verdict(fd)
  for (nmv in list(c("11e_floor_detector", "fd"), c("11e_floor_verdict", "vd"))) {
    f <- dir_agg(sprintf("%s_reps%d.csv", nmv[1], nrep))
    utils::write.csv(get(nmv[2]), f, row.names = FALSE)
    write_provenance(basename(f))
  }
  message(sprintf("  Spearman rho with c: categories %.2f | floor_position %.2f | sd_over_max %.2f",
                  mean(vd$rho_categories), mean(vd$rho_floor_position),
                  mean(vd$rho_sd_over_max)))

  ## ---- 11c: propagation into a pooled RG -------------------------------------
  message("STUDY 11c: meta-analytic propagation (zero new model fits)")
  cc <- rel_meta_propagation(a$raw, M = 2000L, arm = "omega", cores = cores)
  f <- dir_agg("11c_meta_propagation_reps2000.csv")
  utils::write.csv(cc, f, row.names = FALSE)
  write_provenance(basename(f))
  message(sprintf("  REML+HK coverage %.3f-%.3f | FE coverage %.3f-%.3f | I2 inflation up to %.1f pp",
                  min(cc$cov_reml_hk, na.rm = TRUE), max(cc$cov_reml_hk, na.rm = TRUE),
                  min(cc$cov_fe, na.rm = TRUE), max(cc$cov_fe, na.rm = TRUE),
                  max(cc$i2_inflation, na.rm = TRUE)))

  invisible(list(head_to_head = a, loading_shape = b, propagation = cc,
                 floor_sweep = e, floor_detector = fd, floor_verdict = vd,
                 delta = dl, stability = st, gate = gate,
                 claims = list(k_sweep = d_k, bifactor = d_bf,
                               hw = d_hw, tau_equiv = d_te)))
}
