## =============================================================================
## STUDY 09 -- odds ratio -> correlation (r, z)   [`or_to_cor`]
##
## QUESTION
##   convert_df(measure = "r") converts an odds ratio into a correlation via
##   `or_to_cor`, which accepts four values, all shipped and all reachable:
##
##     "bonett"        (the package default since 2.0.1)
##     "pearson"
##     "digby"
##     "lipsey_cooper"
##
##   THEY DO NOT ALL TARGET THE SAME QUANTITY.
##
##     pearson / digby / bonett -> the TETRACHORIC correlation: the correlation
##         of the two underlying CONTINUOUS normal variables that were
##         dichotomised to produce the 2x2 table. All three are closed-form
##         approximations to the same estimand and differ only in accuracy:
##           pearson  r = cos(pi / (1 + OR^0.5))
##           digby    r = (OR^0.75 - 1) / (OR^0.75 + 1)
##           bonett   the pearson cosine with the exponent c adjusted for the
##                    table margins (Bonett & Price 2005), which is the whole
##                    point of the method -- it should pay off precisely where
##                    the margins are unbalanced and be a no-op where they are
##                    not (at p_exp = p_case = 0.5, bonett's c reduces to 1/2
##                    and it IS pearson).
##
##     lipsey_cooper -> a different estimand entirely. It routes through the Cox
##         logit transform d = log(OR) * sqrt(3)/pi (R/es_from_stand_OR.R:255)
##         and then applies the Lipsey-Cooper SMD -> r formula
##         (R/es_from_stand_OR.R:268). The Cox step treats the OUTCOME as a
##         latent continuous variable and returns a standardised mean difference
##         on it; the second step converts that to the POINT-BISERIAL correlation
##         between the binary GROUPING variable and that latent outcome:
##             r_pb = rho * dnorm(qnorm(1 - p_exp)) / sqrt(p_exp * (1 - p_exp))
##         The multiplier is below 1 for every p_exp (0.798 at its maximum,
##         p_exp = 0.5), so lipsey_cooper always targets a SMALLER quantity than
##         the tetrachoric -- and a different one again from phi, the
##         binary-binary correlation. At rho = 0.5, p_exp = 0.3, p_case = 0.1 the
##         three are 0.500 (tetrachoric), 0.379 (point-biserial) and 0.257 (phi).
##
##   Scoring all four against one reference therefore measures estimand mismatch
##   and reports it as bias -- the error study 01 documents for `smd_to_cor`.
##   This study records a per-method scale- AND estimand-matched target so the
##   two effects separate.
##
## WHY THIS STUDY EXISTS NOW
##   `or_to_cor` had no simulation study, yet it is the parameter whose default
##   behaviour changed in 2.0.1: before that release `or_to_cor = "bonett"` never
##   actually ran, because it requires `small_margin_prop` and that column had no
##   auto-derivation, so a blank value silently fell through to the
##   `lipsey_cooper` result. NEWS records the fix moving |r| by +8% to +58%.
##   That is a larger shift than several parameters that DO have studies, and it
##   was made without simulation evidence. This study supplies it.
##
## DESIGN
##   Two data-generating mechanisms, the same four estimators applied to both --
##   the structure of study 03, which asks the same estimand question of the
##   2x2 -> correlation route.
##
##     CONT  bivariate normal with correlation rho, both variables dichotomised
##           at the specified marginals. The tetrachoric estimand IS rho, so the
##           three cosine methods are on their home ground and lipsey_cooper is
##           targeting the (smaller) phi.
##     CAT   a genuinely 4-cell multinomial with population phi = rho. Nothing
##           was dichotomised, no latent normal exists, and phi is the estimand;
##           the tetrachoric methods are now the mismatched ones.
##
## TARGETS
##   theta_pop     the population correlation the data were generated from
##                 (rho = tetrachoric for CONT, rho = phi for CAT), r scale.
##   theta_sample  the same quantity recomputed on that replication's own sample.
##   theta_own     the per-method estimand, on that method's own scale:
##                   tetrachoric family -> population tetrachoric  (atanh() for z)
##                   lipsey_cooper      -> population phi          (atanh() for z)
##                 `own` isolates computational error; the gap between `own` and
##                 `population` IS the estimand mismatch and must be read as
##                 such, not as bias.
##
## CONTINUITY CORRECTION
##   The odds ratio is undefined when any cell is 0. As a meta-analyst would, we
##   add 0.5 to all four cells of an affected table (Haldane-Anscombe) before
##   computing OR and its SE. This is applied identically to every method, so it
##   cannot favour one, but it does attenuate all of them at small n and low
##   event rates -- visible as a shared negative bias in the n = 25 blocks.
## =============================================================================

## This file defines functions only. Load the harness first:
##   source("run_all.R")   from the simulations/ root
## then call run_09().

## ---- population quantities ---------------------------------------------------
##
## .biv_p11(), .p11_to_phi(), .phi_to_p11(), .phi_attainable() and
## .p11_to_tetrachoric() were defined here. They now live in R/05_latent_2x2.R,
## because study 03 needs the same quantities to state ITS per-method estimands
## (roadmap 3.2) and was carrying its own copies of two of them under different
## names. Bodies moved unchanged -- see that file. Nothing else in this study
## changes, and its numbers do not move.

## ---- data-generating mechanisms ---------------------------------------------

## CAT: genuinely 4-cell multinomial, population phi = rho.
gen_cat_or <- function(cond, nrep) {
  n <- cond$n; rho <- cond$rho; p <- cond$p_exp; q <- cond$p_case
  p11 <- .phi_to_p11(rho, p, q)
  prob <- c(p11, p - p11, q - p11, 1 - p - q + p11)  # a, b, c, d

  cells <- stats::rmultinom(nrep, size = n, prob = prob)
  a <- cells[1, ]; b <- cells[2, ]; c_ <- cells[3, ]; d <- cells[4, ]

  ## sample phi of the realised table (= Pearson r on the two 0/1 vectors)
  theta_sample <- (a * d - b * c_) /
    sqrt(as.numeric(a + b) * (c_ + d) * (a + c_) * (b + d))

  data.frame(a = a, b = b, c = c_, d = d,
             theta_pop = rho, theta_sample = theta_sample)
}

## CONT: bivariate normal dichotomised at the margins; tetrachoric estimand = rho.
gen_cont_or <- function(cond, nrep) {
  n <- cond$n; rho <- cond$rho; p <- cond$p_exp; q <- cond$p_case
  cut_x <- stats::qnorm(1 - p); cut_y <- stats::qnorm(1 - q)
  S <- rbind(c(1, rho), c(rho, 1))

  out <- vapply(seq_len(nrep), function(i) {
    z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = S)
    x <- z[, 1] >= cut_x; y <- z[, 2] >= cut_y
    c(a = sum(x & y), b = sum(x & !y), c = sum(!x & y), d = sum(!x & !y),
      ## the correlation the trialist would have reported without dichotomising
      theta_sample = stats::cor(z[, 1], z[, 2]))
  }, numeric(5))

  dat <- as.data.frame(t(out))
  dat$theta_pop <- rho
  dat
}

## ---- estimators --------------------------------------------------------------

## OR and its log-scale SE from a 2x2, with the Haldane-Anscombe correction
## applied to tables containing a zero cell.
.or_from_cells <- function(a, b, c, d) {
  zero <- (a == 0 | b == 0 | c == 0 | d == 0)
  a2 <- ifelse(zero, a + 0.5, a); b2 <- ifelse(zero, b + 0.5, b)
  c2 <- ifelse(zero, c + 0.5, c); d2 <- ifelse(zero, d + 0.5, d)
  list(or = (a2 * d2) / (b2 * c2),
       logor_se = sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2),
       corrected = zero)
}

estimate_or_to_cor <- function(dat, cond) {
  p <- cond$p_exp; q <- cond$p_case
  n_tot <- dat$a + dat$b + dat$c + dat$d
  oo <- .or_from_cells(dat$a, dat$b, dat$c, dat$d)

  ## observed margins of each realised table -- what a meta-analyst would report
  n_exp_v      <- dat$a + dat$b          # exposed
  n_cases_v    <- dat$a + dat$c          # events
  n_nexp_v     <- dat$c + dat$d
  n_controls_v <- dat$b + dat$d

  ## ---- per-method population estimands (constants within a condition) --------
  ## CONT: the latent correlation is rho by construction; phi follows from the
  ## bivariate-normal quadrant probability.
  ## CAT:  phi is rho by construction; the tetrachoric is solved numerically.
  if (identical(cond$dgm, "CONT")) {
    tetra_pop <- cond$rho
    phi_pop   <- .p11_to_phi(.biv_p11(cond$rho, p, q), p, q)
  } else {
    phi_pop   <- cond$rho
    tetra_pop <- .p11_to_tetrachoric(.phi_to_p11(cond$rho, p, q), p, q)
  }

  ## lipsey_cooper's estimand is NEITHER the tetrachoric NOR phi. The Cox
  ## transform d = log(OR)*sqrt(3)/pi treats the OUTCOME as a latent continuous
  ## variable and returns a standardised mean difference on it; the Lipsey-Cooper
  ## d -> r step then converts that to the POINT-BISERIAL correlation between the
  ## binary GROUPING variable and that latent outcome. For a normal latent
  ## dichotomised at the (1 - p_exp) quantile,
  ##     r_pb = rho_latent * dnorm(qnorm(1 - p_exp)) / sqrt(p_exp * (1 - p_exp))
  ## and the multiplier is < 1 for every p_exp (max 0.798 at p_exp = 0.5), so
  ## lipsey_cooper is always targeting something SMALLER than the tetrachoric.
  ## The residual against this target is the logistic-vs-normal latent mismatch
  ## the Cox transform carries, and isolating it is the point of scoring on
  ## `own`: benchmarked against the tetrachoric instead, that mismatch and the
  ## estimand gap are summed together and reported as one large "bias".
  pb_mult <- stats::dnorm(stats::qnorm(1 - p)) / sqrt(p * (1 - p))
  pb_pop  <- tetra_pop * pb_mult

  ## estimand each method targets, on the r scale
  own_r <- c(bonett = tetra_pop, pearson = tetra_pop, digby = tetra_pop,
             lipsey_cooper = pb_pop)

  pieces <- list()
  for (m in c("bonett", "pearson", "digby", "lipsey_cooper")) {
    res <- es_from_or_se(
      or = oo$or, logor_se = oo$logor_se,
      n_exp = n_exp_v, n_nexp = n_nexp_v,
      n_cases = n_cases_v, n_controls = n_controls_v,
      n_sample = n_tot,
      or_to_cor = m
    )
    for (meas in c("r", "z")) {
      pc <- as_method(res, meas, paste0(m, " (", meas, ")"))
      pc$route <- "package"
      ## SCALE + ESTIMAND MATCHING. The (z) routes report a Fisher's z, so
      ## scoring them against an r-scale reference reports atanh(rho) - rho as
      ## bias; and lipsey_cooper targets phi, not the tetrachoric, so scoring it
      ## against the tetrachoric reports the estimand gap as bias. theta_own
      ## carries the matched target for both reasons at once.
      pc$theta_own <- if (meas == "r") own_r[[m]] else atanh(own_r[[m]])
      pieces[[paste(m, meas)]] <- pc
    }
  }

  ## Reference route, for the actionable comparison "should I use the OR at all?"
  ## When the full 2x2 is available the package converts it directly, without
  ## passing through the odds ratio. Any accuracy the or_to_cor approximations
  ## give up relative to this is the price of having only the OR.
  res2 <- es_from_2x2(n_cases_exp = dat$a, n_controls_exp = dat$b,
                      n_cases_nexp = dat$c, n_controls_nexp = dat$d,
                      table_2x2_to_cor = "tetrachoric")
  for (meas in c("r", "z")) {
    pc <- as_method(res2, meas, paste0("2x2_tetrachoric (", meas, ")"))
    pc$route <- "package (full 2x2)"
    pc$theta_own <- if (meas == "r") tetra_pop else atanh(tetra_pop)
    pieces[[paste("2x2", meas)]] <- pc
  }

  do.call(rbind, pieces)
}

## ---- grid --------------------------------------------------------------------
## Deliberately the grid of study 03, so the two 2x2 -> correlation routes are
## directly comparable. p_exp = allocation to the exposed arm, p_case = overall
## event rate; together they set the margin imbalance that bonett corrects for.
build_grid_09 <- function(dgm) {
  g <- expand.grid(
    rho    = c(0, 0.1, 0.25, 0.5, 0.75),
    n      = c(25, 50, 100, 300),
    p_exp  = c(0.3, 0.5),
    p_case = c(0.1, 0.3, 0.5),
    KEEP.OUT.ATTRS = FALSE
  )
  g$dgm <- dgm
  ## CAT generates from phi directly, so the requested phi must be attainable
  ## with these margins. CONT generates from rho and is always attainable.
  if (identical(dgm, "CAT")) {
    ok <- mapply(.phi_attainable, g$rho, g$p_exp, g$p_case)
    if (any(!ok))
      message(sprintf("  [09_or_to_cor] %d/%d CAT conditions dropped: phi not ",
                      sum(!ok), nrow(g)),
              "attainable at those margins")
    g <- g[ok, , drop = FALSE]
  }
  g
}

## ---- entry point -------------------------------------------------------------
run_09 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  for (pkg in c("mvtnorm", "MASS")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("package '", pkg, "' is required by study 09")
  }

  ## `own` is the scale- and estimand-matched target and is the one to rank on.
  ## `population`/`sample` stay on the r scale and are interpretable only for the
  ## (r) routes; their gap from `own` is the estimand mismatch.
  tg <- c(own = "theta_own", population = "theta_pop", sample = "theta_sample")

  message("STUDY 09a: OR -> cor, CONTINUOUS latent (tetrachoric is the estimand)")
  a <- run_study("09a_or_to_cor_CONT", build_grid_09("CONT"),
                 gen_cont_or, estimate_or_to_cor, nrep, cores, targets = tg)

  message("STUDY 09b: OR -> cor, CATEGORICAL latent (phi is the estimand)")
  b <- run_study("09b_or_to_cor_CAT", build_grid_09("CAT"),
                 gen_cat_or, estimate_or_to_cor, nrep, cores, targets = tg)

  invisible(list(cont = a, cat = b))
}

if (sys.nframe() == 0L) run_09()
