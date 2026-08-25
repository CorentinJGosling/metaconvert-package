## =============================================================================
## STUDY 10 -- the reliability family (ICC, Cronbach's alpha)
##
## WHY THIS EXISTS
##   Until now `simulations/` had ZERO coverage of the psychometric family:
##   grep -i "alpha|icc|omega" over studies/ and R/ returned nothing, while nine
##   studies covered the conversion routes. That gap mattered more than a missing
##   row in a coverage table, because the reliability audit produced several
##   numbers that lived only in commit messages and code comments -- each one a
##   simulation run once and thrown away. Anything asserted that way cannot be
##   re-checked when the code moves underneath it, which is exactly what happened
##   to the "coverage around 0.74-0.76" claim this study now replaces.
##
##   Three quantities are made reproducible here, at a fixed seed:
##
##     10a  the ICC(2,1) agreement-type coverage table       (roadmap 1.1)
##     v36_false_positive_rate()  the reliability-induction FP rate (roadmap 2.1)
##     icc_se_by_k()              how far k moves the ICC SE       (V42's premise)
##
## WHAT IS *NOT* HERE, AND WHY
##   Roadmap 3.1 also asks for the low-ICC Bonett-vs-Fisher-TF pooling comparison
##   (roadmap 1.5). It is deliberately absent, because
##   `icc_to_es = "fisher_tf"` IS NOT IMPLEMENTED.
##   Simulating it would mean re-typing the candidate formula inline
##   and scoring the package against a private copy of a method it does not ship
##   -- the exact practice this simulation programme was rebuilt to stop (see the
##   runner's header). When 1.5 ships, add a `10c` that calls the real route.
##
##   Item 1.5 is also not merely unimplemented but currently unsafe as specified:
##   the bonett-vs-else branch at internal_flags.R:2144 fails OPEN, so a
##   fisher_tf pool would raise a false [INVALID] on every ICC > 0.762 at k = 2.
##   That has to be designed before it can be measured.
##
## STRUCTURE
##   10a uses the run_study() harness, because it scores an ESTIMATOR against a
##   known target and that is what the harness is for. The other two are
##   standalone deterministic functions in the style of 99_route_equivalence.R:
##   they measure a FLAG's behaviour, not an estimate's, so they have no target
##   column and do not fit the generate/estimate contract.
## =============================================================================


## ---- 10a: data-generating mechanism -----------------------------------------
##
## Two-way random-effects design, which is what an absolute-agreement ICC(2,1)
## assumes and what the package's SE formula does NOT:
##
##     y_ij = mu + a_i + b_j + e_ij
##       a_i ~ N(0, s2_subj)     subject
##       b_j ~ N(0, s2_rater)    rater      <- the term the SE formula omits
##       e_ij ~ N(0, s2_err)     residual
##
##     ICC(2,1) = s2_subj / (s2_subj + s2_rater + s2_err)
##
## `rater_share` is the fraction of the NON-SUBJECT variance carried by raters.
## At rater_share = 0 the design collapses to one-way and the package formula is
## correct; everything above 0 is the regime the formula approximates away.
##
## The meta-analyst never sees the variance components -- only the ICC estimate,
## n and k -- so that is all `generate()` hands on. The estimate is the
## Shrout & Fleiss (1979) ICC(2,1):
##
##     ICC = (MSR - MSE) / (MSR + (k-1) MSE + k (MSC - MSE) / n)
##
## which is what a primary study would have computed and printed.
gen_icc <- function(cond, nrep) {
  n <- cond$n; k <- cond$k
  rho <- cond$icc; rater_share <- cond$rater_share

  ## Solve the variance components for the requested ICC and rater share.
  ## Fix the total at 1: s2_subj = rho, and the remaining (1 - rho) splits into
  ## rater and residual by rater_share.
  s2_subj  <- rho
  s2_rater <- (1 - rho) * rater_share
  s2_err   <- (1 - rho) * (1 - rater_share)

  out <- vapply(seq_len(nrep), function(i) {
    a <- stats::rnorm(n, 0, sqrt(s2_subj))
    b <- stats::rnorm(k, 0, sqrt(s2_rater))
    y <- outer(a, b, "+") + matrix(stats::rnorm(n * k, 0, sqrt(s2_err)), n, k)

    gm <- mean(y); ri <- rowMeans(y); cj <- colMeans(y)
    MSR <- k * sum((ri - gm)^2) / (n - 1)                       # between subjects
    MSC <- n * sum((cj - gm)^2) / (k - 1)                       # between raters
    MSE <- sum((y - outer(ri, rep(1, k)) -
                  outer(rep(1, n), cj) + gm)^2) / ((n - 1) * (k - 1))

    icc_hat <- (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n)
    ## An ICC estimate can fall outside the parameter space in small samples.
    ## Clamp just inside it: es_from_icc() refuses |icc| >= 1 by design, and a
    ## row NA'd there would be scored as non-estimable rather than as the
    ## extreme value it is.
    c(icc = min(max(icc_hat, -0.999), 0.999), n_sample = n, n_measurements = k)
  }, numeric(3))

  dat <- as.data.frame(t(out))
  ## The estimand, on the analysis scale the package reports.
  dat$theta_pop    <- log(1 - rho)
  ## No meaningful same-sample target: the ICC is a variance-component ratio, so
  ## the "sample value" IS the estimate being scored. Record the population value
  ## in both slots rather than inventing a second quantity.
  dat$theta_sample <- log(1 - rho)
  dat
}


## ---- 10a: estimator ----------------------------------------------------------
## One vectorised call to the package route, per the runner's contract. Both
## icc_type levels are run: they are byte-identical by construction (the SE is
## one shared formula), and the point of including both is to DEMONSTRATE that
## identity in the output rather than assert it in a comment -- it is the reason
## flipping the default to "consistency" was rejected as a remedy for 1.1.
estimate_icc <- function(dat, cond) {
  pieces <- list()
  for (ty in c("agreement", "consistency")) {
    res <- es_from_icc(
      icc = dat$icc, n_sample = dat$n_sample,
      n_measurements = dat$n_measurements, icc_type = ty,
      ## the closed-form (n, k) SE is the object under study, so ask for it
      ## explicitly rather than depending on convert_df's default
      agreement_se = "compute"
    )
    pc <- as_method(res, "icc", ty)
    pc$route <- "package"
    pc$scale <- "bonett"
    pieces[[ty]] <- pc
  }
  do.call(rbind, pieces)
}


build_grid_10a <- function() {
  expand.grid(
    ## n spans the range where the problem is claimed to WORSEN, which is the
    ## whole point: a bounded "0.74-0.76" reads as a fixed penalty, and the
    ## defect is that coverage falls as n grows.
    n           = c(20, 50, 200, 1000),
    k           = c(2, 3),
    icc         = c(0.80),
    ## 0 is the null case the formula is exact for -- it is in the grid so the
    ## comparison is internal and does not rest on trusting the implementation.
    rater_share = c(0, 0.2, 0.5, 0.8),
    KEEP.OUT.ATTRS = FALSE
  )
}


run_10 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()

  message("STUDY 10a: ICC(2,1) agreement-type SE, coverage against rater variance")
  a <- run_study("10a_icc_agreement_coverage", build_grid_10a(),
                 gen_icc, estimate_icc, nrep, cores)

  # 10b SHIPPED WITHOUT A PRODUCER. The CSV was in data/aggregated/ but nothing in the
  # codebase wrote it -- so it could not be regenerated by run_everything(), and under
  # the provenance scheme it could not be recorded either. v36_false_positive_rate() is
  # seeded, and reproduces the shipped file byte-for-byte (md5
  # dcbf27425441079c85b66ad54bdbfbb8), so giving it the writer it always needed changes
  # no number. It is not a run_study() call -- it measures a FLAG, not an estimator, so
  # it has no target column -- hence the explicit write plus write_provenance().
  message("STUDY 10b: V36 reliability-induction false-positive rate")
  reps_10b <- 400L
  b <- v36_false_positive_rate(reps = reps_10b)
  dir.create(dir_agg(), recursive = TRUE, showWarnings = FALSE)
  f_b <- dir_agg(sprintf("10b_v36_false_positive_reps%d.csv", reps_10b))
  utils::write.csv(b, f_b, row.names = FALSE)
  write_provenance(basename(f_b))
  message("  -> ", f_b)

  invisible(list(icc_coverage = a, v36_fp = b))
}


## =============================================================================
## Deterministic checks -- not run_study(), because they measure a FLAG rather
## than an estimator and so have no target column. Same role as
## ancova_equivalence() in 99_route_equivalence.R.
## =============================================================================

#' Measured false-positive rate of V36 (reliability induction)
#'
#' V36 flags a reliability coefficient reported by two DIFFERENT studies, as a
#' possible sign of reliability induction. Every alpha here is drawn
#' independently, so there is NO induction anywhere and every flag raised is a
#' false positive.
#'
#' The result decides a design question. The roadmap proposed tightening the
#' entropy gate to require 4 decimals; that is refuted by the 4-decimal column,
#' because the driver is POOL SIZE, not precision -- with ~150 plausible
#' three-decimal values in [.78, .93] a 30-study pool has C(30,2)/150 ~ 3
#' expected collisions, so a match is the norm. No decimal gate fixes a birthday
#' problem, which is why the message was reworded instead.
#'
#' `n_round` matters more than it looks: round sample sizes (multiples of 50)
#' collide with each other as readily as two-decimal alphas do, and the
#' 2-decimal clause of the gate keys on an identical n_sample.
#'
#' @param k pool sizes to test
#' @param dp decimal places the coefficient is reported to
#' @param n_round whether sample sizes are round numbers
#' @param reps pools per cell
#' @return a data.frame: pool_fp = share of pools with at least one flag,
#'   row_fp = share of rows flagged
v36_false_positive_rate <- function(k = c(10, 30, 100), dp = c(2, 3, 4),
                                    n_round = c(TRUE, FALSE), reps = 400L) {
  load_metaconvert()
  grid <- expand.grid(k = k, dp = dp, n_round = n_round, KEEP.OUT.ATTRS = FALSE)
  out <- vector("list", nrow(grid))

  for (g in seq_len(nrow(grid))) {
    kk <- grid$k[g]; dd <- grid$dp[g]; nr <- grid$n_round[g]
    set.seed(condition_seed("10_v36_false_positive", g), kind = "L'Ecuyer-CMRG")
    pools <- 0L; rows <- 0L

    for (r in seq_len(reps)) {
      a <- round(stats::runif(kk, 0.78, 0.93), dd)
      n <- if (nr) sample(seq(50, 500, by = 50), kk, TRUE) else sample(50:500, kk, TRUE)
      d <- data.frame(study_id = paste0("S", seq_len(kk)),
                      cronbach_alpha = a, n_sample = n, n_items = 10)
      ## Tier 1 directly: V36 is an input-data check, so the whole effect-size
      ## machinery is irrelevant to it. Going through convert_df() + summary()
      ## instead makes the same sweep ~100x slower for an identical answer.
      iss <- suppressWarnings(
        metaConvert:::.validate_input_data(d, verbose = FALSE))$issues
      hit <- grepl("PROMPT TO CHECK THE PRIMARY SOURCES", iss, fixed = TRUE)
      if (any(hit)) pools <- pools + 1L
      rows <- rows + sum(hit)
    }
    out[[g]] <- data.frame(k = kk, dp = dd, n_round = nr,
                           pool_fp = pools / reps,
                           row_fp = rows / (reps * kk))
  }
  do.call(rbind, out)
}


#' How far the measurement count moves the ICC standard error
#'
#' The premise behind V42. V37 (its alpha counterpart) was written for an
#' 8-vs-18 item-count mix-up, which moves the alpha SE by 1.039x -- far under
#' D2's 3x outlier gate, so nothing else in the package can see it. This
#' quantifies the same thing for k, and the answer is rho-dependent, which is
#' why the roadmap's single "1.435x" could not be reproduced as a single number.
#'
#' Deterministic: no simulation, just the closed-form SE evaluated across k.
#'
#' @param rho ICC values to evaluate
#' @param k measurement counts
#' @param n sample size
#' @return a data.frame with the max/min SE ratio across k, per rho
icc_se_by_k <- function(rho = c(0.3, 0.5, 0.7, 0.8, 0.9),
                        k = 2:10, n = 100) {
  load_metaconvert()
  out <- lapply(rho, function(rr) {
    se <- vapply(k, function(kk) {
      es_from_icc(rr, n, kk, agreement_se = "compute")$icc_se
    }, numeric(1))
    data.frame(rho = rr, k_min = min(k), k_max = max(k), n = n,
               se_at_k_min = se[1], se_at_k_max = se[length(se)],
               ratio_max_min = max(se) / min(se))
  })
  do.call(rbind, out)
}


#' The identity that killed one of roadmap 1.1's two proposed remedies
#'
#' Flipping the default `icc_type` from "agreement" to "consistency" was offered
#' as a fix for the anti-conservative agreement SE. It cannot be one: the two
#' types share a single SE formula, so the change is numerically inert and only
#' relabels rows. Asserting that in a comment is weaker than showing it, and the
#' difference is exactly 0 rather than merely small.
#'
#' @return a data.frame of the paired estimates and their differences
icc_type_is_inert <- function() {
  load_metaconvert()
  cells <- expand.grid(icc = c(-0.2, 0.5, 0.8, 0.95), n = c(30, 120),
                       k = c(2, 3, 5), KEEP.OUT.ATTRS = FALSE)
  res <- lapply(seq_len(nrow(cells)), function(i) {
    p <- cells[i, ]
    A <- es_from_icc(p$icc, p$n, p$k, icc_type = "agreement",
                     agreement_se = "compute")
    C <- es_from_icc(p$icc, p$n, p$k, icc_type = "consistency",
                     agreement_se = "compute")
    data.frame(icc = p$icc, n = p$n, k = p$k,
               es_diff = abs(A$icc - C$icc),
               se_diff = abs(A$icc_se - C$icc_se))
  })
  do.call(rbind, res)
}
