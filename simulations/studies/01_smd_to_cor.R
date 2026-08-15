## =============================================================================
## STUDY 01 -- SMD -> correlation (r, z)   [`smd_to_cor`]
##
## QUESTION
##   convert_df(measure = "r") converts a Cohen's d from two groups into a
##   correlation via `smd_to_cor`, default "viechtbauer", alternative
##   "lipsey_cooper". Both are shipped and both are reachable.
##
##   THEY DO NOT TARGET THE SAME QUANTITY:
##     lipsey_cooper -> the POINT-BISERIAL correlation, i.e. the correlation
##                      between the binary grouping variable and the outcome
##     viechtbauer   -> the BISERIAL correlation, i.e. the correlation between
##                      the underlying CONTINUOUS variable that was dichotomised
##                      to form the groups, and the outcome
##
##   Comparing them against a single reference therefore measures estimand
##   mismatch and calls it bias. The old pipeline did exactly that, benchmarked
##   on the latent correlation, and reported a Lipsey-Cooper "bias" that is flat
##   in n (-0.189 at n=25, -0.182 at n=300). Estimand gaps do not shrink with n;
##   that flatness is the signature.
##
##   This study records THREE targets so the two effects separate cleanly.
##
## VERIFIED AGAINST THE SOURCE PAPER
##   Jacobs & Viechtbauer (2017), Res Synth Methods 8:161-180
##   (papers/Jacobs_Viechtbauer_2017_RSM_biserial.pdf). The package implements it
##   faithfully -- checked term by term against internal_multiple_formulas.R:
##
##     Eq [17]  z = (a/2)*ln((1 + a*r_b)/(1 - a*r_b)),  a = sqrt(f(z_p))/(p*q)^(1/4)
##              -> a_viechtbauer <- sqrt(fzp)/(p*(1-p))^(1/4); same z
##     Var(z) = 1/(n-1)          -> vz_viechtbauer <- 1/(n_exp + n_nexp - 1)
##              (since NEWS.md 2.0.1 this is rescaled by vd/vd_crude, which is
##               exactly 1 for the crude two-group rows this study generates)
##     Eq [19]  inverse          -> same, applied to both bounds
##     truncate |r_b| > 1 to +-1 -> r_trunc <- ifelse(...)
##
##   The paper's worked example (r_b = 0.46, n = 20, p = 0.5) reproduces exactly:
##   a = 0.893, z = 0.39, CI (-0.08, 0.82). (The published PDF's text layer drops
##   minus signs, so its printed lower bounds read 0.06 / 0.08; 0.39 - 1.96*sqrt(1/19)
##   is negative.)
##
##   SCOPE: this paper is about estimating the BISERIAL CORRELATION from
##   dichotomised data, i.e. SMD -> r. It says nothing about r -> SMD. The
##   `cor_to_smd = "viechtbauer"` option (study 02) is named for
##   metafor::transf.rtod, Viechtbauer's function, NOT for this paper, and must
##   not be justified by citing it.
## =============================================================================

## ---- data-generating mechanism ----------------------------------------------
## Bivariate normal (X, Y) with correlation rho. X is dichotomised at the
## (1 - p) quantile to form the two groups; the meta-analyst sees only the
## per-group means, SDs and sample sizes.
gen_smd <- function(cond, nrep) {
  n <- cond$n; rho <- cond$rho; p <- cond$p_exp
  cut <- stats::qnorm(1 - p)
  S <- rbind(c(1, rho), c(rho, 1))

  ## Population point-biserial for a normal dichotomised at `cut`:
  ##   r_pb = rho * f / sqrt(p * (1 - p)),  f = dnorm(qnorm(p))
  ## which is exactly the Viechtbauer biserial transform read backwards.
  f <- stats::dnorm(cut)
  theta_pb_pop <- rho * f / sqrt(p * (1 - p))

  out <- vapply(seq_len(nrep), function(i) {
    z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = S)
    g <- z[, 1] >= cut
    y1 <- z[g, 2]; y0 <- z[!g, 2]
    c(n_exp = length(y1), n_nexp = length(y0),
      mean_exp = mean(y1), mean_nexp = mean(y0),
      sd_exp = stats::sd(y1), sd_nexp = stats::sd(y0),
      theta_sample = stats::cor(z[, 1], z[, 2]),          # latent, same sample
      theta_pb_sample = stats::cor(as.numeric(g), z[, 2]) # point-biserial, same sample
    )
  }, numeric(8))

  dat <- as.data.frame(t(out))
  dat$theta_pop    <- rho             # latent (biserial) population value
  dat$theta_pb_pop <- theta_pb_pop    # point-biserial population value

  ## Fisher-z counterparts, so the (z) methods are scored on their own scale.
  dat$theta_pop_z       <- atanh(dat$theta_pop)
  dat$theta_sample_z    <- atanh(dat$theta_sample)
  dat$theta_pb_pop_z    <- atanh(dat$theta_pb_pop)
  dat$theta_pb_sample_z <- atanh(dat$theta_pb_sample)
  ## Degenerate draws (an empty or singleton arm) leave NA SDs; they are kept so
  ## they surface as nonest_rate rather than vanishing.
  dat
}

## ---- estimators: package only, one vectorised call per method ---------------
##
## SCALE MATCHING. r and z are different scales, so they need different targets:
## comparing a Fisher-z estimate to an r-scale target reports atanh(rho) - rho as
## if it were bias (at rho = 0.5 that alone is +0.049). Each scale is therefore
## run separately against its own targets, rather than producing a grid of
## method x target cells half of which are meaningless.
## Viechtbauer's variance-stabilising transform, applied to a population value so
## the (z) route can be scored on its OWN scale rather than on Fisher's z.
vs_transform <- function(r, p) {
  a <- sqrt(stats::dnorm(stats::qnorm(p))) / (p * (1 - p))^(1 / 4)
  (a / 2) * log((1 + a * r) / (1 - a * r))
}

estimate_smd_scale <- function(measure) {
  force(measure)
  function(dat, cond) {
    p <- cond$p_exp; rho <- cond$rho
    pieces <- list()
    for (m in c("viechtbauer", "lipsey_cooper")) {
      res <- es_from_means_sd(
        mean_exp = dat$mean_exp, mean_sd_exp = dat$sd_exp,
        mean_nexp = dat$mean_nexp, mean_sd_nexp = dat$sd_nexp,
        n_exp = dat$n_exp, n_nexp = dat$n_nexp,
        smd_to_cor = m
      )
      pc <- as_method(res, measure, m)
      pc$route <- "package"
      pc$scale <- measure

      ## theta_own = the estimand THIS route actually targets, on THIS scale.
      pc$theta_own <- if (measure == "r") {
        if (m == "viechtbauer") dat$theta_pop else dat$theta_pb_pop
      } else {
        if (m == "viechtbauer") vs_transform(rho, p) else dat$theta_pb_pop_z
      }
      pieces[[m]] <- pc
    }
    do.call(rbind, pieces)
  }
}

build_grid_01 <- function() {
  expand.grid(
    rho   = c(0, 0.1, 0.25, 0.5, 0.75),
    n     = c(25, 50, 100, 300),
    p_exp = c(0.3, 0.5, 0.7),
    KEEP.OUT.ATTRS = FALSE
  )
}

run_01 <- function(nrep = SIM_DEFAULTS$nrep, cores = SIM_DEFAULTS$cores) {
  load_metaconvert()
  grid <- build_grid_01()

  ## biserial targets  = what viechtbauer   estimates
  ## point-biserial    = what lipsey_cooper estimates
  message("STUDY 01a: SMD -> correlation, r scale")
  a <- run_study("01a_smd_to_cor_r", grid, gen_smd, estimate_smd_scale("r"),
                 nrep, cores,
                 targets = c(own                      = "theta_own",
                             biserial_population      = "theta_pop",
                             biserial_sample          = "theta_sample",
                             pointbiserial_population = "theta_pb_pop",
                             pointbiserial_sample     = "theta_pb_sample"))

  message("STUDY 01b: SMD -> correlation, Fisher-z scale")
  b <- run_study("01b_smd_to_cor_z", grid, gen_smd, estimate_smd_scale("z"),
                 nrep, cores,
                 targets = c(own                      = "theta_own",
                             fisherz_biserial         = "theta_pop_z",
                             fisherz_pointbiserial    = "theta_pb_pop_z"))

  invisible(list(r = a, z = b))
}
