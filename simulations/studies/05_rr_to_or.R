## =============================================================================
## STUDY 05 -- risk ratio -> odds ratio   [`rr_to_or`]
##
## QUESTION
##   The mirror of study 04. `rr_to_or` defaults to "metaumbrella", alternatives
##   "transpose", "grant", "dipietrantonj". Same three families, same
##   baseline-risk misspecification axis.
##
##   This direction matters more than it looks: OR is the metric most
##   meta-analytic software assumes, so RR -> OR is the conversion applied when
##   a review pools trials that reported RRs alongside trials that reported ORs.
##   An error here is applied to the studies that were reported CORRECTLY.
##
## NON-ESTIMABILITY IS A RESULT
##   Grant's inverse transform OR = RR*(1 - p0) / (1 - RR*p0) is undefined
##   whenever RR * p0 >= 1 -- i.e. precisely when the analyst over-guesses the
##   baseline risk for a protective-to-neutral RR. The legacy script dropped
##   those replications with a bare na.rm = TRUE, so Grant was evaluated
##   conditional on its own transform not blowing up. Here they are counted:
##   see `nonest_rate` in the aggregate output.
##
## CANDIDATE ADDED (not in metaConvert)
##   vanderweele_rr_squared  OR ~= RR^2, the inverse of the minimax conversion in
##   VanderWeele TJ (2020) Biometrics 76(3):746-752. Like its mirror in study 04,
##   its interval is the TRANSFORMED RR limits widened by the bias-ratio bound --
##   not a symmetric interval, which the paper rules out for a transformation that
##   is an approximation. The widening factor for THIS direction (1.5625 at the
##   [0.2, 0.8] band) is derived rather than cited; see vw_bias_ratio_sq() below.
##   Its guarantee, like study 04's, covers only the part of the grid where both
##   outcome probabilities lie in [0.2, 0.8]: br = 0.30 with rr in {0.75, 1, 2}
##   and br = 0.50 with rr in {0.5, 0.75, 1}. Read its coverage in two regions.
## =============================================================================

## Reuses gen_2x2_rr() from studies/04_or_to_rr.R, but the targets flip to the
## log odds ratio.
gen_2x2_or <- function(cond, nrep) {
  dat <- gen_2x2_rr(cond, nrep)

  ## gen_2x2_rr() returns theta_sample on the log RR scale. This function is about
  ## to overwrite it with the log OR, which is correct for THIS study's target --
  ## but any estimator needing the sample log RR must not read theta_sample after
  ## that point. It gets its own name, and the two are never the same column again.
  ## (Roadmap 3.1: vanderweele_rr_squared read theta_sample AFTER the overwrite and
  ## so squared the odds ratio instead of the risk ratio, for exactly this reason.)
  dat$logrr_sample <- dat$theta_sample

  dat$rr       <- exp(dat$logrr_sample)
  dat$logrr_se <- sqrt(1 / dat$a - 1 / (dat$a + dat$b) +
                       1 / dat$c - 1 / (dat$c + dat$d))
  dat$rr_ci_lo <- exp(dat$logrr_sample - stats::qnorm(.975) * dat$logrr_se)
  dat$rr_ci_up <- exp(dat$logrr_sample + stats::qnorm(.975) * dat$logrr_se)

  ## population log OR implied by (rr, br): risks are rr*br and br
  r1 <- cond$rr * cond$br; r0 <- cond$br
  dat$theta_pop    <- log((r1 / (1 - r1)) / (r0 / (1 - r0)))
  dat$theta_sample <- log(dat$or)          # sample logOR from the same table

  ## The invariant the 3.1 defect violated, asserted where it can be seen.
  ## NB deliberately NOT asserting that the two columns differ: a replication with
  ## a/n_exp == c/n_nexp gives log OR = log RR = 0 legitimately, so at small nrep
  ## that check could abort a valid run.
  stopifnot(
    isTRUE(all.equal(dat$theta_sample, log(dat$or))),    # target is the log OR
    isTRUE(all.equal(dat$logrr_sample, log(dat$rr)))     # and the log RR is kept apart
  )
  dat
}

## ---- candidate estimator (not in metaConvert) -------------------------------

## The bias-ratio bound for the INVERSE conversion, OR ~= RR^2.
##
## PROVENANCE, WHICH MATTERS HERE: VanderWeele (2020) states the 1.25 bound and
## the divide/multiply interval rule for the OR -> RR direction only (p.748). He
## does not write the RR -> OR direction. The factor below is therefore DERIVED,
## not cited -- but it is exact, not a heuristic, and the derivation is one line:
##
##     RR^2 / OR = (p1/p0)^2 * (p0(1 - p1)) / (p1(1 - p0)) = p1(1 - p1) / (p0(1 - p0))
##
## Over a band [w, u] = [0.5 - v, 0.5 + v] the numerator is at most 0.25 and the
## denominator at least w(1 - w), so the bias ratio is bounded by
## 0.25 / (w(1 - w)) = 1 / (1 - 4v^2) -- exactly the SQUARE of the OR -> RR factor
## 1 / sqrt(1 - 4v^2), which is what one expects of squaring a bounded ratio.
## At [0.2, 0.8] that is 1.25^2 = 1.5625. Verified numerically over a 0.0005 grid:
## measured max 1.5625, attained at p0 = 0.8, p1 = 0.5, and the closed form agrees
## everywhere; bands [0.1,0.9] / [0.3,0.7] / [0.4,0.6] match 1/(1-4v^2) to 6 dp.
vw_bias_ratio_sq <- function(w = 0.2, u = 0.8) vw_bias_ratio(w, u)^2   # from study 04

## OR ~= RR^2, with the conservative interval implied by the bound above.
##
## Same structure as study 04's vanderweele_sqrt_or(): the interval is the
## TRANSFORMED RR limits, then widened -- not `estimate +/- z * se`. The paper's
## objection to transforming an interval through an approximation applies in this
## direction too, so the mirror study must not use the symmetric interval that
## roadmap item 1.9 removed from study 04. `logor_se` remains the delta-method SE
## of the point estimate and is NOT the basis of the interval; performance() reads
## coverage from the CI columns and se_ratio from the SE column independently.
vanderweele_rr_squared <- function(logrr, logrr_se, w = 0.2, u = 0.8) {
  logor <- 2 * logrr
  se    <- 2 * logrr_se                     # delta method, for the ESTIMATE only
  widen <- log(vw_bias_ratio_sq(w, u))
  data.frame(logor = logor, logor_se = se,
             logor_ci_lo = logor - stats::qnorm(.975) * se - widen,
             logor_ci_up = logor + stats::qnorm(.975) * se + widen)
}

estimate_rr_to_or <- function(dat, cond) {
  br_guess <- cond$br_guess
  pieces <- list()

  for (m in c("metaumbrella", "transpose", "grant")) {
    res <- es_from_rr_se(
      rr = dat$rr, logrr_se = dat$logrr_se,
      baseline_risk = rep(br_guess, nrow(dat)),
      n_exp = dat$n_exp, n_nexp = dat$n_nexp,
      n_cases = dat$n_cases, n_controls = dat$n_controls,
      rr_to_or = m
    )
    p <- as_method(res, "logor", m); p$route <- "package"
    pieces[[m]] <- p
  }

  res_dp <- es_from_rr_ci(
    rr = dat$rr, rr_ci_lo = dat$rr_ci_lo, rr_ci_up = dat$rr_ci_up,
    baseline_risk = rep(br_guess, nrow(dat)),
    n_exp = dat$n_exp, n_nexp = dat$n_nexp,
    n_cases = dat$n_cases, n_controls = dat$n_controls,
    rr_to_or = "dipietrantonj"
  )
  p <- as_method(res_dp, "logor", "dipietrantonj"); p$route <- "package"
  pieces[["dipietrantonj"]] <- p

  ## Candidate: the inverse of VanderWeele's minimax conversion, OR ~= RR^2.
  ##
  ## READS logrr_sample, NOT theta_sample. theta_sample is the log OR in this study
  ## (it is the target), and until 2026-08 this line read it -- so the candidate
  ## computed 2*log(OR) while its own comment, and the SE on the next line, say
  ## log RR. Scored against theta_sample the bias was then 2*logOR - logOR = logOR
  ## exactly; measured mean bias 1.2601 against mean|logOR| 1.2601 at
  ## (br = 0.30, rr = 2, n = 300). See gen_2x2_or()'s stopifnot.
  vw <- vanderweele_rr_squared(dat$logrr_sample, dat$logrr_se)
  p <- as_method(vw, "logor", "vanderweele_rr_squared"); p$route <- "candidate"
  pieces[["vw"]] <- p

  do.call(rbind, pieces)
}

run_05 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  message("STUDY 05: RR -> OR (rr_to_or), with baseline-risk misspecification")
  run_study("05_rr_to_or", build_grid_04(), gen_2x2_or, estimate_rr_to_or, nrep, cores)
}
