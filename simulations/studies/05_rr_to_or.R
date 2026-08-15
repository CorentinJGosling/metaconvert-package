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
## =============================================================================

## Reuses gen_2x2_rr() from studies/04_or_to_rr.R, but the targets flip to the
## log odds ratio.
gen_2x2_or <- function(cond, nrep) {
  dat <- gen_2x2_rr(cond, nrep)
  dat$rr       <- exp(dat$theta_sample)
  dat$logrr_se <- sqrt(1 / dat$a - 1 / (dat$a + dat$b) +
                       1 / dat$c - 1 / (dat$c + dat$d))
  dat$rr_ci_lo <- exp(dat$theta_sample - stats::qnorm(.975) * dat$logrr_se)
  dat$rr_ci_up <- exp(dat$theta_sample + stats::qnorm(.975) * dat$logrr_se)

  ## population log OR implied by (rr, br): risks are rr*br and br
  r1 <- cond$rr * cond$br; r0 <- cond$br
  dat$theta_pop    <- log((r1 / (1 - r1)) / (r0 / (1 - r0)))
  dat$theta_sample <- log(dat$or)          # sample logOR from the same table
  dat
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
  logor <- 2 * dat$theta_sample
  se <- 2 * dat$logrr_se
  vw <- data.frame(logor = logor, logor_se = se,
                   logor_ci_lo = logor - stats::qnorm(.975) * se,
                   logor_ci_up = logor + stats::qnorm(.975) * se)
  p <- as_method(vw, "logor", "vanderweele_rr_squared"); p$route <- "candidate"
  pieces[["vw"]] <- p

  do.call(rbind, pieces)
}

run_05 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  message("STUDY 05: RR -> OR (rr_to_or), with baseline-risk misspecification")
  run_study("05_rr_to_or", build_grid_04(), gen_2x2_or, estimate_rr_to_or, nrep, cores)
}
