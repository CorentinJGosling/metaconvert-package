## =============================================================================
## STUDY 03 -- 2x2 contingency table -> correlation (r, z)
##
## QUESTION
##   convert_df(measure = "r") turns a 2x2 table into a correlation via
##   `table_2x2_to_cor`. The package currently accepts only "tetrachoric".
##   Tetrachoric assumes both binary variables are crude dichotomisations of
##   underlying CONTINUOUS normal variables. When the variables are genuinely
##   dichotomous (dead/alive, relapsed/not), that assumption is false and the phi
##   coefficient is the estimand.
##
##   So the choice is not "which formula is more accurate" -- it is "which
##   estimand did the trialist measure". This study quantifies the cost of
##   getting that wrong in each direction, which is the actionable result: it
##   tells a reviewer when the package default is safe and when it is not.
##
## DESIGN
##   Two data-generating mechanisms, same estimators applied to both:
##     CONT  a bivariate normal with correlation rho, both variables then
##           dichotomised at given marginals. Population target = rho.
##     CAT   a genuinely 4-cell multinomial with population phi = rho.
##           Population target = phi.
##
## TARGETS
##   theta_pop     the population correlation the data were generated from
##   theta_sample  the same correlation recomputed on that replication's own
##                 sample (Pearson r on the two 0/1 vectors = the sample phi)
## =============================================================================

## This file defines functions only. Load the harness first:
##   source("run_all.R")   from the simulations/ root
## then call run_03().

## ---- data-generating mechanisms ---------------------------------------------

## Solve the joint cell probability p11 that gives a target phi with fixed
## marginals p (row) and q (column). phi has a closed form here, so no root
## finding and no escape hatch: the original script used uniroot() with a `9e9`
## fallback that silently returned NA for 100% of replications in the
## (r = 0.5, br = 0.7) block, which is why 30 rows are missing from the shipped
## CAT file. Instead we check attainability up front and skip the condition.
phi_to_p11 <- function(phi, p, q) phi * sqrt(p * q * (1 - p) * (1 - q)) + p * q

phi_attainable <- function(phi, p, q) {
  p11 <- phi_to_p11(phi, p, q)
  ## all four cells must be non-negative probabilities
  p11 >= 0 && p11 <= min(p, q) && (p11 - p - q + 1) >= 0 &&
    (p - p11) >= 0 && (q - p11) >= 0
}

gen_cat <- function(cond, nrep) {
  n <- cond$n; rho <- cond$rho; p <- cond$p_exp; q <- cond$p_case
  p11 <- phi_to_p11(rho, p, q)
  prob <- c(p11, p - p11, q - p11, 1 - p - q + p11)   # exp&case, exp&ctrl, nexp&case, nexp&ctrl

  cells <- stats::rmultinom(nrep, size = n, prob = prob)
  a <- cells[1, ]; b <- cells[2, ]; c_ <- cells[3, ]; d <- cells[4, ]

  ## sample phi from the realised table (identical to Pearson r on the 0/1 data)
  theta_sample <- (a * d - b * c_) /
    sqrt(as.numeric(a + b) * (c_ + d) * (a + c_) * (b + d))

  data.frame(a = a, b = b, c = c_, d = d,
             theta_pop = rho, theta_sample = theta_sample)
}

gen_cont <- function(cond, nrep) {
  n <- cond$n; rho <- cond$rho; p <- cond$p_exp; q <- cond$p_case
  cut_x <- stats::qnorm(1 - p); cut_y <- stats::qnorm(1 - q)
  S <- rbind(c(1, rho), c(rho, 1))

  out <- vapply(seq_len(nrep), function(i) {
    z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = S)
    x <- z[, 1] >= cut_x; y <- z[, 2] >= cut_y
    c(a = sum(x & y), b = sum(x & !y), c = sum(!x & y), d = sum(!x & !y),
      ## For the CONT mechanism the meta-analyst's implicit target is the latent
      ## continuous correlation, so theta_sample is the sample correlation of the
      ## UNDERLYING variables -- what the trialist would have reported had they
      ## not dichotomised. This is the honest same-sample analogue.
      theta_sample = stats::cor(z[, 1], z[, 2]))
  }, numeric(5))

  dat <- as.data.frame(t(out))
  dat$theta_pop <- rho
  dat
}

## ---- estimators: metaConvert, called ONCE per condition ----------------------

## NOTE ON REACHABILITY -- read before interpreting any result here ------------
##
## As of metaConvert v2.0.0 there is exactly ONE 2x2 -> correlation method:
##   * es_from_2x2() hard-stops on anything but "tetrachoric"
##     (R/es_from_2x2.R:95; "cooper"/"lipsey" are commented out of the guard)
##   * the internal .contingency_to_cor() "lipsey" branch is a `# TO DO` stub,
##     entirely commented out, returning NULL (R/internal_multiple_formulas.R:484)
##   * "cooper" is absent from the internal altogether
##
## So the comparison a user would want cannot be made through the package at all.
## That is itself a headline finding, and it is why `route` is recorded on every
## row: the paper must be able to say which estimates came from shipped code and
## which from a candidate implementation.
##
## WARNING for whoever un-comments that TO DO block: its phi formula
## (internal_multiple_formulas.R:490-491) divides by
##   (n_cases_exp + n_controls_exp)(n_cases_nexp + n_controls_nexp) *
##   (n_cases_exp + n_cases_nexp)(n_controls_exp + n_cases_nexp)
## The last factor must be the CONTROLS margin (n_controls_exp + n_controls_nexp),
## not (n_controls_exp + n_cases_nexp). The same error is in the legacy simulation
## (archive/legacy_scripts/6_SIM_2x2_to_COR_CAT.R:99-101), which is what created
## the "Lipsey bias" reported in the old aggregates. `phi_candidate()` below is
## the corrected version.

## Candidate implementation, NOT currently in the package. phi = Pearson r on the
## two 0/1 indicators; SE via the Fisher-z route on the logOR variance, which is
## the construction the package's TO DO block intends.
phi_candidate <- function(a, b, c, d) {
  r <- (a * d - b * c) / sqrt(as.numeric(a + b) * (c + d) * (a + c) * (b + d))
  log_or <- log((a * d) / (c * b))
  v_log_or <- 1 / a + 1 / b + 1 / c + 1 / d
  z <- atanh(r)
  vz <- v_log_or * (z^2) / (log_or^2)
  z_lo <- z - stats::qnorm(.975) * sqrt(vz)
  z_up <- z + stats::qnorm(.975) * sqrt(vz)
  eff_n <- 1 / vz + 3
  vr <- (1 - r^2)^2 / (eff_n - 1)
  data.frame(r = r, r_se = sqrt(vr), r_ci_lo = tanh(z_lo), r_ci_up = tanh(z_up),
             z = z, z_se = sqrt(vz), z_ci_lo = z_lo, z_ci_up = z_up)
}

estimate_2x2 <- function(dat, cond) {
  routes <- list(
    tetrachoric = list(route = "package",
      f = function() es_from_2x2(n_cases_exp     = dat$a, n_controls_exp  = dat$b,
                                 n_cases_nexp    = dat$c, n_controls_nexp = dat$d,
                                 table_2x2_to_cor = "tetrachoric")),
    phi = list(route = "candidate",
      f = function() phi_candidate(dat$a, dat$b, dat$c, dat$d))
  )
  pieces <- list()
  for (m in names(routes)) {
    res <- routes[[m]]$f()
    for (meas in c("r", "z")) {
      p <- as_method(res, meas, paste0(m, " (", meas, ")"))
      p$route <- routes[[m]]$route
      # SCALE MATCHING. The (z) routes report a Fisher-z, so scoring them against
      # theta_pop -- which is on the r scale -- reports atanh(rho) - rho as if it
      # were bias. At rho = 0.75 that alone is +0.223, and it dominates: the
      # correlation between phi (z)'s apparent bias and the pure scale gap is 0.993.
      # theta_own carries the scale-matched target so the `own` column is honest;
      # `population` and `sample` remain on the r scale and are meaningful only for
      # the (r) routes.
      p$theta_own <- if (meas == "r") dat$theta_pop else atanh(dat$theta_pop)
      pieces[[paste(m, meas)]] <- p
    }
  }
  do.call(rbind, pieces)
}

## ---- grid --------------------------------------------------------------------
## p_exp  = proportion allocated to the exposed arm
## p_case = overall event rate (the "baseline risk" axis)
build_grid <- function() {
  g <- expand.grid(
    rho    = c(0, 0.1, 0.25, 0.5, 0.75),
    n      = c(25, 50, 100, 300),
    p_exp  = c(0.3, 0.5),
    p_case = c(0.1, 0.3, 0.5),
    KEEP.OUT.ATTRS = FALSE
  )
  ## Drop conditions where the requested phi is not attainable with these
  ## marginals, and RECORD how many were dropped rather than emitting silent NAs.
  ok <- mapply(phi_attainable, g$rho, g$p_exp, g$p_case)
  if (any(!ok))
    message(sprintf("  [03_2x2_to_cor] %d/%d conditions dropped: phi unattainable ",
                    sum(!ok), nrow(g)),
            "at those marginals (max |phi| is bounded by the margins)")
  g[ok, , drop = FALSE]
}

## ---- entry point -------------------------------------------------------------
run_03 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  grid <- build_grid()

  # `own` is the scale-matched target (r for the (r) routes, atanh(r) for the (z)
  # routes); `population`/`sample` stay on the r scale and are interpretable only
  # for the (r) routes. See the note in estimate_2x2().
  tg <- c(own = "theta_own", population = "theta_pop", sample = "theta_sample")

  message("STUDY 03a: 2x2 -> cor, CATEGORICAL latent (phi is the estimand)")
  a <- run_study("03a_2x2_to_cor_CAT", grid, gen_cat, estimate_2x2, nrep, cores,
                 targets = tg)

  message("STUDY 03b: 2x2 -> cor, CONTINUOUS latent (tetrachoric is the estimand)")
  b <- run_study("03b_2x2_to_cor_CONT", grid, gen_cont, estimate_2x2, nrep, cores,
                 targets = tg)

  invisible(list(cat = a, cont = b))
}

if (sys.nframe() == 0L) run_03()
