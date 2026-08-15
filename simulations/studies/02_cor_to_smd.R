## =============================================================================
## STUDY 02 -- correlation -> SMD (d, g)   [`cor_to_smd`]
##
## QUESTION
##   convert_df(measure = "d"/"g") converts a reported correlation into a
##   standardised mean difference via `cor_to_smd`, default "viechtbauer",
##   alternatives "cooper" and "mathur" (R/es_from_REGRESSION.R:241, the shipped
##   set is exactly those three). All are reachable through es_from_pearson_r().
##
##   EVERY ONE OF THEM PRESUPPOSES A DICHOTOMISATION THAT WAS NEVER PERFORMED.
##   A correlation carries no group structure; a d is defined only once one
##   exists. So the argument does not choose an approximation to a common
##   quantity, it chooses WHICH grouping to pretend happened, and the routes
##   disagree about that:
##
##     cooper        d = 2r / sqrt(1 - r^2)   (R/internal_multiple_formulas.R:1278)
##                   the EXACT inverse of the point-biserial when r came from a
##                   genuine two-group comparison with EQUAL arms. Cooper eq.
##                   12.38-12.39.
##     viechtbauer   metafor::conv.delta(transf = metafor::transf.rtod)
##                   (R/internal_multiple_formulas.R:1267-1274). Verified against
##                   metafor 4.x: with n1i/n2i absent, transf.rtod computes
##                   r_pb = lambda*r, lambda = dnorm(0)/0.5 = 0.7978846, and
##                   returns 2*r_pb/sqrt(1 - r_pb^2). That is the BISERIAL
##                   inverse: it reads r as the correlation with a latent
##                   continuous variable and reports the d a MEDIAN SPLIT of that
##                   variable would produce. Not the same estimand as cooper --
##                   it is smaller by ~20-25% across the usable range.
##     mathur        d per `unit_increase_iv` units of the exposure (Mathur &
##                   VanderWeele 2020). A THIRD scale: it coincides with cooper
##                   only at unit_type = "sd", unit_increase_iv = 2, which the
##                   package's own docs state and its defaults do not supply.
##
## VERIFIED AGAINST THE SOURCE PAPERS (data-raw/R to SMD/)
##
##   Mathur MB & VanderWeele TJ (2020) "A Simple, Interpretable Conversion from
##   Pearson's Correlation to Cohen's d for Continuous Exposures", Epidemiology
##   31(2):e16-e18. Its two equations are exactly the two shipped routes:
##
##     Eq (1.1)  d = 2r/sqrt(1-r^2),  SE = sqrt(2/((N-1)(1-r^2)))
##               == cor_to_smd = "cooper". The paper states plainly that this
##               "was derived for a Pearson correlation computed between a BINARY
##               exposure X and a continuous outcome Y, also called a
##               point-biserial correlation".
##     Eq (1.2)  d = r*Delta/(s_x*sqrt(1-r^2))
##               SE = |d| * sqrt(1/(r^2*(N-3)) + 1/(2*(N-1)))
##               == cor_to_smd = "mathur", term for term, SE included:
##                 increase <- ifelse(unit_type == "sd", unit_increase_iv*sd_iv,
##                                    unit_increase_iv)
##                 d <- r * increase / (sd_iv * sqrt(1 - r^2))
##                 d_se <- abs(d) * sqrt(1/(r^2*(n_sample-3)) + 1/(2*(n_sample-1)))
##
##   The paper also supplies an exact internal consistency check this study should
##   assert: applying Eq (1.1) to a continuous-X correlation "coincides with the
##   effect size associated with an increase in X of two standard deviations",
##   i.e. mathur at unit_type = "sd", unit_increase_iv = 2 must reproduce cooper's
##   POINT ESTIMATE exactly -- algebraically r*2*s_x/(s_x*sqrt(1-r^2)) =
##   2r/sqrt(1-r^2) -- while the paper is equally explicit that the two SEs
##   "will, in general, still not coincide". Both halves are checkable.
##
##   So `unit_type = "raw_scale"` is NOT a silent error, contrary to what a first
##   pass suggests: with unit_type != "sd" the code sets Delta to raw units of X,
##   which IS the paper's Delta. The defect is documentation only -- ?convert_df
##   says unit_type "must be either 'sd' or 'value'" while the shipped default is
##   the third spelling "raw_scale", and nothing validates the argument, so a typo
##   silently selects raw units. Fix the docs or the default; the arithmetic is right.
##
##   Two DGPs make each estimand true in turn, as in studies 01 and 03: under
##   GROUPS the reported r is a genuine point-biserial, under CONT it is a
##   continuous-continuous Pearson r and the target is what a split would have
##   given. Scoring each route against its OWN estimand separates computational
##   error from the estimand gap; scoring it against theta_pop is what a user who
##   trusts the default actually gets.
##
## WHAT THIS DESIGN IS BUILT TO MEASURE, AND ONE THING IT IS NOT
##
##   1. Group balance. Both shipped routes are written at p = 0.5: cooper's `2` is
##      1/sqrt(p(1-p)) there, and viechtbauer's lambda is dnorm(qnorm(p))/
##      sqrt(p(1-p)) there. For cooper this is FAITHFUL TO ITS SOURCE -- Mathur &
##      VanderWeele Eq (1.1) is stated for the balanced point-biserial -- so it is
##      a scope limit, not a coding error. What is worth measuring is the COST:
##      at the population value, true two-group d = 0.8 with r the matching
##      point-biserial, cooper returns 0.800 at p = 0.5 but 0.571 at p = 0.15,
##      while the general inverse d = r/sqrt((1-r^2)*p*(1-p)) recovers 0.800.
##
##      BUT DO NOT REPORT "expose p" AS THE FIX. Pustejovsky
##      (data-raw/Converting-from-d-to-r-to-z PUSTEJOVSKY.pdf) pins the balanced
##      case -- "the minimum possible value of the factor a = (n1+n2)^2/(n1*n2) is
##      achieved when the two groups are of equal size. In this case, a = 4, which
##      is equivalent to assuming that the absolute value of the treatment-control
##      differential is w = 2" (a = 4 <=> d = 2r/sqrt(1-r^2) <=> cooper) -- and
##      then argues AGAINST generalising via a: Borenstein's and Hunter-Schmidt's
##      formulas "replace w^2 with the factor a ... However, in controlled
##      experiments, these proportions are arbitrary (and often equal) and do not
##      provide any information about the magnitude of the treatment-control
##      differential on the continuous variable X."
##
##      So the answer is design-conditional and this study must report it that way:
##        dichotomisation / extreme groups -> proportions are informative, use p
##        controlled experiment            -> allocation ratio is arbitrary; using
##                                            p imports absent information, and
##                                            Pustejovsky parameterises by w instead
##      Note metaConvert ALREADY handles unequal groups in the reverse direction:
##      smd_to_cor = "lipsey_cooper" computes a <- ((n_exp+n_nexp)^2)/(n_exp*n_nexp)
##      and uses r = d/sqrt(d^2 + a). The asymmetry is r->d only.
##      The `metafor_rtod_np` candidate below measures the cost of a = 4 under each
##      mechanism; it is evidence for a design-conditional recommendation, NOT a
##      winner. A w-parameterised candidate is the obvious third arm to add.
##
##   2. NOT a defect: the d/g slot assignment. It looks wrong -- cooper puts its
##      map value in d and derives g = J*d, while viechtbauer puts conv.delta's
##      output in g and backs out d = g/J -- and since transf.rtod carries no df
##      term it is tempting to call its output an uncorrected d and "fix" the
##      asymmetry. That was tried and is WRONG. Which slot a value belongs in is a
##      question about estimator bias, not about the algebra of the transform.
##      Simulation (bivariate normal, rho = 0.5, median split, true delta =
##      0.870126, nrep = 2e5), bias of the reported g against delta:
##
##        n      Hedges g from raw data   g = transf.rtod   g = J*transf.rtod
##        10           -0.004                 +0.033             -0.055
##        25           -0.000                 +0.012             -0.017
##       200           -0.000                 +0.002             -0.002
##
##      The shipped convention is closer to unbiased at every n, because E[r_hat]
##      is biased down and partly cancels the upward bias the pooled-SD
##      denominator gives a directly computed d -- so this route does not carry
##      the full Hedges bias and applying J on top over-corrects. g = J*d holds on
##      every route, so measure = "d" and measure = "g" both return what they
##      should. See tests/testthat/test-cor-to-smd-d-g.R and the note in
##      R/internal_multiple_formulas.R. Do not "fix" this again.
##
## PERFORMANCE NOTE (package, not simulation)
##   `cor_to_smd = "viechtbauer"` costs 18.1 s per 10,000 rows against 0.23 s for
##   "cooper" -- 75x -- because .cor_to_smd() is reached through a per-row
##   mapply() and each row builds its own conv.delta() call, while conv.delta is
##   itself fully vectorised (0.34 s for the same 10,000 rows called once). This
##   is the same defect .mapply_memo() was introduced for on the tetrachoric
##   path, and it is on the DEFAULT route.
## =============================================================================

## ---- the population maps every route is a special case of --------------------

## Multiplier taking a latent correlation to the point-biserial induced by
## dichotomising the other variable at its (1 - p) quantile:
##   r_pb = rho * dnorm(qnorm(p)) / sqrt(p(1-p))
## the same constant study 01 reads in the opposite direction. Symmetric in
## p <-> 1-p and maximised at p = 0.5, where it is lambda = 0.7978846.
.split_factor <- function(p) stats::dnorm(stats::qnorm(p)) / sqrt(p * (1 - p))

## Exact inverse of the point-biserial at group proportion p. From
##   r_pb = d*sqrt(p(1-p)) / sqrt(d^2*p(1-p) + 1)
## invert to
##   d = r_pb / (sqrt(p(1-p)) * sqrt(1 - r_pb^2)).
## At p = 0.5 this IS cooper's 2r/sqrt(1-r^2): the shipped formula is the
## balanced-arm special case, and nothing in ?es_from_pearson_r says so.
.pb_to_d <- function(r_pb, p) r_pb / (sqrt(p * (1 - p)) * sqrt(1 - r_pb^2))

## Composing the two gives every route in the study as one two-parameter family
## (k = correlation-scale multiplier, p = group proportion):
##   cooper      k = 1,               p = 0.5
##   viechtbauer k = lambda,          p = 0.5
##   pb_exact_p  k = 1,               p = observed
##   biserial_p  k = .split_factor(p) p = observed
.k_p_to_d <- function(r, k, p) .pb_to_d(k * r, p)

## Delta-method SE of that map. Taken from metafor::conv.delta rather than
## differentiated by hand: conv.delta reproduces the package's own cooper d_se
## to machine precision (checked), so the candidates differ from the shipped
## routes ONLY in the transform, never in the variance machinery.
.delta_frame <- function(r, r_se, k, p, df) {
  cd <- metafor::conv.delta(yi = r, vi = r_se^2,
                            transf = function(x) .k_p_to_d(x, k, p),
                            var.names = c("d", "d_var"))
  se <- sqrt(cd$d_var)
  ## qt(.975, n-2) matches .es_from_d()'s own CI construction
  ## (R/internal_es_from_d.R:180-181), so coverage differences are not an
  ## artefact of z-vs-t.
  data.frame(d = cd$d, d_se = se,
             d_ci_lo = cd$d - stats::qt(.975, df) * se,
             d_ci_up = cd$d + stats::qt(.975, df) * se)
}

## ---- DGP A: GROUPS -----------------------------------------------------------
## There really were two arms: N(delta, 1) vs N(0, 1). The meta-analyst reports
## the point-biserial r between arm and outcome, which is what a regression paper
## reporting "r" for a binary predictor reports. theta_pop = delta, the Cohen's d
## the trial actually estimates.
gen_cor_groups <- function(cond, nrep) {
  n <- cond$n; delta <- cond$delta
  n_exp <- round(n * cond$p_exp); n_nexp <- n - n_exp
  p <- n_exp / n                      # the REALISED proportion after rounding

  Y1 <- matrix(stats::rnorm(nrep * n_exp,  delta, 1), nrow = nrep)
  Y0 <- matrix(stats::rnorm(nrep * n_nexp, 0,     1), nrow = nrep)
  m1 <- rowMeans(Y1); m0 <- rowMeans(Y0)
  ybar <- (n_exp * m1 + n_nexp * m0) / n

  ## Row-wise sums of squares; the vector subtraction recycles down columns, so
  ## each row is centred on its own mean. Kept matrix-wide rather than looping
  ## cor()/sd() per replication -- generate() may loop, but it need not.
  SS_tot <- rowSums((Y1 - ybar)^2) + rowSums((Y0 - ybar)^2)
  SS_w   <- rowSums((Y1 - m1)^2)   + rowSums((Y0 - m0)^2)

  data.frame(
    ## the datum the meta-analyst extracts
    pearson_r = sqrt(n_exp * n_nexp / n) * (m1 - m0) / sqrt(SS_tot),
    n_sample = n, n_exp = n_exp, n_nexp = n_nexp,
    p_real = p,
    ## population value of the reported correlation: the argument every route's
    ## estimand map is evaluated at
    r_pop = delta * sqrt(p * (1 - p)) / sqrt(delta^2 * p * (1 - p) + 1),
    theta_pop = delta,
    ## the SMD the trialist would have reported from the same participants,
    ## uncorrected (so the g rows carry the Hedges factor as a visible offset)
    theta_sample = (m1 - m0) / sqrt(SS_w / (n - 2))
  )
}

## ---- DGP B: CONT -------------------------------------------------------------
## Two continuous variables, bivariate normal. The meta-analyst reports the
## Pearson r; the estimand is the d a split of X at its (1-p) quantile WOULD have
## produced -- the counterfactual dichotomisation the conversion silently invokes.
## At p = 0.5 the closed form reduces to the median-split identity
##   theta_pop = 2*lambda*rho / sqrt(1 - lambda^2*rho^2),  lambda = dnorm(0)/0.5,
## which is exactly .pb_to_d(lambda*rho, 0.5); the general-p version is the same
## expression with lambda replaced by .split_factor(p). Derivation: the pooled
## within-group variance of Y after splitting at the (1-p) quantile is
## 1 - r_pb^2 with r_pb = .split_factor(p)*rho, and the mean difference is
## rho*dnorm(qnorm(p))/(p(1-p)).
gen_cor_cont <- function(cond, nrep) {
  n <- cond$n; rho <- cond$rho
  n_up <- round(n * cond$p_split); n_lo <- n - n_up
  p <- n_up / n

  X <- matrix(stats::rnorm(nrep * n), nrow = nrep)
  Y <- rho * X + sqrt(1 - rho^2) * matrix(stats::rnorm(nrep * n), nrow = nrep)

  Xc <- X - rowMeans(X); Yc <- Y - rowMeans(Y)
  r  <- rowSums(Xc * Yc) / sqrt(rowSums(Xc^2) * rowSums(Yc^2))

  ## theta_sample: split THIS replication's own X at its n_lo-th order statistic
  ## and compute the d that split gives. Ties have probability zero for
  ## continuous X, so `X > cut` selects exactly n_up observations.
  cut  <- apply(X, 1, function(v) sort(v, partial = n_lo)[n_lo])
  up   <- X > cut
  s_up <- rowSums(Y * up)
  m_up <- s_up / n_up
  m_lo <- (rowSums(Y) - s_up) / n_lo
  SS_w <- rowSums((Y - (up * m_up + (!up) * m_lo))^2)

  data.frame(
    pearson_r = r,
    n_sample = n, n_exp = n_up, n_nexp = n_lo,
    p_real = p,
    r_pop = rho,
    theta_pop = .pb_to_d(.split_factor(p) * rho, p),
    theta_sample = (m_up - m_lo) / sqrt(SS_w / (n - 2))
  )
}

## ---- estimators --------------------------------------------------------------
## One shared function: every route's estimand is its own map applied to the
## population correlation the analyst reports, which is r_pop under BOTH
## mechanisms. The mechanism decides which map is correct, not which maps exist.
estimate_cor_to_smd <- function(dat, cond) {
  n_r  <- nrow(dat)
  p    <- dat$p_real[1]
  rho  <- dat$r_pop[1]
  df   <- dat$n_sample - 2
  ones <- rep(1, n_r)
  lam  <- .split_factor(0.5)   # what BOTH shipped routes hard-code
  k_p  <- .split_factor(p)

  pk <- function(res, measure, label, own) {
    x <- as_method(res, measure, label)
    x$route <- "package"; x$theta_own <- own; x
  }

  ## n_exp/n_nexp are supplied throughout. They do not reach the transform in any
  ## shipped route -- that is the point of finding 1 -- and since they always sum
  ## to n_sample the CI df is n-2 either way; passing them changes nothing and
  ## keeps one code path across the two mechanisms.
  call_pkg <- function(method, ...) es_from_pearson_r(
    pearson_r = dat$pearson_r, n_sample = dat$n_sample,
    n_exp = dat$n_exp, n_nexp = dat$n_nexp, cor_to_smd = method, ...)

  res_co <- call_pkg("cooper")
  res_vi <- call_pkg("viechtbauer")
  ## The ONE mathur setting that lands on the two-group d scale, per
  ## ?es_from_pearson_r: "coincides with mathur only in the special case
  ## unit_type = 'sd' and unit_increase_iv = 2". sd_iv = 1 matches both DGPs.
  res_m2 <- call_pkg("mathur", sd_iv = ones, unit_increase_iv = 2 * ones,
                     unit_type = rep("sd", n_r))
  ## What a user gets from the SIGNATURE DEFAULT unit_type = "raw_scale" with the
  ## natural unit_increase_iv = 1. "raw_scale" is in neither of the two values
  ## the argument's own documentation permits ("sd", "value"); .cor_to_smd() only
  ## tests `unit_type == "sd"`, so any other string silently behaves as "value"
  ## and the result is d per ONE unit of the exposure -- exactly half the
  ## two-group d. Kept as a route, with its own halved estimand, so the trap is
  ## measured rather than asserted.
  res_m1 <- call_pkg("mathur", sd_iv = ones, unit_increase_iv = ones,
                     unit_type = rep("raw_scale", n_r))

  own_pb   <- .k_p_to_d(rho, 1,   0.5)   # cooper / mathur_2sd estimand
  own_bis  <- .k_p_to_d(rho, lam, 0.5)   # viechtbauer estimand
  own_pb_p <- .k_p_to_d(rho, 1,   p)     # point-biserial inverse at observed p
  own_bis_p<- .k_p_to_d(rho, k_p, p)     # biserial inverse at observed split p

  pieces <- list(
    cooper           = pk(res_co, "d", "cooper",           own_pb),
    viechtbauer      = pk(res_vi, "d", "viechtbauer",      own_bis),
    mathur_2sd       = pk(res_m2, "d", "mathur_2sd",       own_pb),
    mathur_raw_1unit = pk(res_m1, "d", "mathur_raw_1unit", own_pb / 2),
    ## Same two routes read out of the g slot. Same estimand -- g is the
    ## small-sample-unbiased estimator of the same population SMD -- so this is
    ## scale-matched, and the d/g pair isolates which slot the Hedges factor
    ## belongs on (finding 2).
    cooper_g         = pk(res_co, "g", "cooper_g",         own_pb),
    viechtbauer_g    = pk(res_vi, "g", "viechtbauer_g",    own_bis)
  )

  ## ---- candidates: the two shipped maps freed of the balanced-arm assumption
  ## The package's own r_se is reused (res_co$r_se, i.e. Cooper eq. 12.27) so the
  ## only thing that differs is the transform.
  r_se <- res_co$r_se
  cnd <- function(k, label, own) {
    x <- as_method(.delta_frame(dat$pearson_r, r_se, k, p, df), "d", label)
    x$route <- "candidate"; x$theta_own <- own; x
  }
  pieces$pb_exact_p        <- cnd(1,   "pb_exact_p",       own_pb_p)
  pieces$biserial_split_p  <- cnd(k_p, "biserial_split_p", own_bis_p)

  ## metafor's transf.rtod used AS DOCUMENTED: n1i/n2i supplied, output taken as
  ## a d. conv.delta forwards ... to transf (verified), so this is literally the
  ## one-line change .cor_to_smd()'s viechtbauer branch would need. It differs
  ## from biserial_split_p only by metafor's finite-sample factor
  ## sqrt((N-2)/N) <= 1, which is small-sample bias and is therefore NOT folded
  ## into its target: same theta_own as biserial_split_p, so the factor shows up
  ## where it belongs.
  mf <- metafor::conv.delta(
    yi = dat$pearson_r, vi = r_se^2, transf = metafor::transf.rtod,
    n1i = dat$n_exp, n2i = dat$n_nexp, var.names = c("d", "d_var"))
  se_mf <- sqrt(mf$d_var)
  x <- as_method(data.frame(d = mf$d, d_se = se_mf,
                            d_ci_lo = mf$d - stats::qt(.975, df) * se_mf,
                            d_ci_up = mf$d + stats::qt(.975, df) * se_mf),
                 "d", "metafor_rtod_np")
  x$route <- "candidate"; x$theta_own <- own_bis_p
  pieces$metafor_rtod_np <- x

  ## Nothing is wrapped in try(). mathur's SE carries a 1/(r^2 (n-3)) term that
  ## diverges as r -> 0, and every .pb_to_d call is undefined at |k*r| >= 1
  ## (unreachable here: k <= lambda < 1). Whatever they produce -- Inf, NaN, NA --
  ## reaches performance() and surfaces as nonest_rate.
  do.call(rbind, pieces)
}

## ---- grids -------------------------------------------------------------------
## p and 1-p are exchangeable in both mechanisms: every quantity in play depends
## on the proportion only through p(1-p) and dnorm(qnorm(p)), both symmetric
## about 0.5. Gridding 0.7 as well as 0.3 (as study 01 does) would duplicate
## conditions, so the imbalance axis runs downward only, to a 85/15 split.
.check_arms <- function(g, pcol, label) {
  n_exp <- round(g$n * g[[pcol]]); n_nexp <- g$n - n_exp
  ok <- n_exp >= 2 & n_nexp >= 2
  if (any(!ok))
    message(sprintf("  [%s] %d/%d conditions dropped: an arm has < 2 observations",
                    label, sum(!ok), nrow(g)))
  else
    message(sprintf("  [%s] %d conditions, none dropped (all arms >= 2)",
                    label, nrow(g)))
  g[ok, , drop = FALSE]
}

build_grid_02_groups <- function() {
  ## delta = 0 is kept deliberately: it is where mathur's SE term 1/(r^2 (n-3))
  ## is closest to diverging, and where a conversion that presupposes a grouping
  ## should at least agree with every other route.
  .check_arms(expand.grid(
    delta = c(0, 0.2, 0.5, 0.8, 1.2),
    n     = c(25, 50, 100, 300),
    p_exp = c(0.5, 0.3, 0.15),
    KEEP.OUT.ATTRS = FALSE
  ), "p_exp", "02a_cor_to_smd_GROUPS")
}

build_grid_02_cont <- function() {
  .check_arms(expand.grid(
    rho     = c(0, 0.1, 0.25, 0.5, 0.75),
    n       = c(25, 50, 100, 300),
    p_split = c(0.5, 0.3, 0.15),
    KEEP.OUT.ATTRS = FALSE
  ), "p_split", "02b_cor_to_smd_CONT")
}

run_02 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()

  ## Both calls score on the d scale (g rows included as separate methods target
  ## the same population SMD, so no scale is mixed). Three targets throughout:
  ##   own        each route's own estimand -> isolates computational error
  ##   population the parameter -> the only valid target for coverage / se_ratio
  ##   sample     the same-sample SMD -> the meta-analyst's faithfulness question
  tg <- c(own = "theta_own", population = "theta_pop", sample = "theta_sample")

  message("STUDY 02a: r -> SMD, GROUPS mechanism (r is a genuine point-biserial)")
  a <- run_study("02a_cor_to_smd_GROUPS", build_grid_02_groups(),
                 gen_cor_groups, estimate_cor_to_smd, nrep, cores, targets = tg)

  message("STUDY 02b: r -> SMD, CONT mechanism (r is continuous-continuous; ",
          "the target is the d a split would have produced)")
  b <- run_study("02b_cor_to_smd_CONT", build_grid_02_cont(),
                 gen_cor_cont, estimate_cor_to_smd, nrep, cores, targets = tg)

  invisible(list(groups = a, cont = b))
}
