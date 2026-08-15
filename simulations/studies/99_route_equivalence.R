## =============================================================================
## STUDY 99 -- ROUTE EQUIVALENCE (package wiring check, not a simulation study)
##
## PURPOSE
##   metaConvert exposes the same underlying information through many different
##   input shapes. A paper may report adjusted means + SDs, or an adjusted mean
##   difference + its SE, or just an ANCOVA t -- and all three describe the SAME
##   analysis. Whenever the shapes are informationally equivalent, the package must
##   return the SAME effect size. If two routes disagree, that is a wiring bug: a
##   mis-derived intermediate, a wrong df, a dropped factor.
##
##   This is deliberately NOT a Monte Carlo study. It is a deterministic identity
##   check on ONE simulated dataset per scenario, so any disagreement is exact and
##   attributable, with no sampling noise to hide behind.
##
## WHY IT IS NEEDED
##   The eight simulation studies exercise one entry point per area. 19 exported
##   functions take `pre_post_to_smd`; study 08 calls one of them. The ANCOVA family
##   has 9 entry points; study 07 calls four. This file walks the rest.
##
## WHAT "EQUIVALENT" MEANS HERE -- read before interpreting a FAIL
##   Two routes are only required to agree when they carry the same information AND
##   target the same estimand. Three cases where disagreement is CORRECT:
##
##   1. ANCOVA statistic-inverting routes under covariate imbalance. `ancova_t`,
##      `_f`, `_pval`, `etasq_adj`, `ancova_md_se/_ci/_pval` and `ancova_means_se`
##      recover d THROUGH the D = 0 variance expression, so they are attenuated by
##      1/sqrt(1 + D/(1/n_exp + 1/n_nexp)) relative to the means/MD+SD routes.
##      Documented in R/es_from_ANCOVA_statistics.R. The scenarios below therefore
##      test equivalence at EXACT covariate balance (D = 0), where all routes must
##      agree, and separately CHECK the predicted attenuation under imbalance.
##
##   2. Pre-post: the mean-change and paired-t contexts do not carry the separate
##      pre/post SDs, so `bonett` and `morris_dav` are unavailable there (coerced by
##      .validate_pre_post_to_smd). Cross-family equivalence is therefore tested only
##      on `morris_drm` and `morris_dz`, which are defined from the change score.
##
##   3. Anything the package documents as a deliberate difference.
##
## HOW TO RUN
##   setwd("simulations"); source("run_all.R"); run_99()
## =============================================================================

## ---- helpers ----------------------------------------------------------------

# compare a list of route results on the columns that matter
eq_report <- function(lbl, results, cols = c("d", "d_se", "g", "g_se"), tol = 1e-8) {
  nm <- names(results)
  ref <- results[[1]]
  out <- data.frame(route = nm, stringsAsFactors = FALSE)
  for (cl in cols) {
    out[[cl]] <- vapply(results, function(r)
      if (is.null(r) || !cl %in% names(r)) NA_real_ else as.numeric(r[[cl]])[1], numeric(1))
  }
  # max abs deviation from the reference route, over the compared columns
  dev <- vapply(seq_along(results), function(i) {
    a <- unlist(out[i, cols]); b <- unlist(out[1, cols])
    ok <- is.finite(a) & is.finite(b)
    if (!any(ok)) NA_real_ else max(abs(a[ok] - b[ok]))
  }, numeric(1))
  out$max_dev_vs_ref <- dev
  out$verdict <- ifelse(is.na(dev), "NOT RUN",
                        ifelse(dev <= tol, "match", "**DIFFERS**"))
  cat("\n--- ", lbl, "  (reference route = ", nm[1], ", tolerance = ", tol, ")\n", sep = "")
  print(out, row.names = FALSE, digits = 10)
  invisible(out)
}

safe <- function(expr) tryCatch(expr, error = function(e) {
  structure(list(err = conditionMessage(e)), class = "route_error")
})

## =============================================================================
## PART 1 -- ANCOVA family
## =============================================================================
##
## One dataset. The covariate is centred WITHIN each arm, so the group means of x
## are exactly equal and the leverage term D is exactly 0. At D = 0 every ANCOVA
## route -- including the statistic-inverting ones -- must return the same d.

ancova_equivalence <- function(n1 = 60, n2 = 60, delta = 0.5, rho = 0.5,
                               balanced = TRUE, seed = 20260801) {
  set.seed(seed)
  x <- c(scale(rnorm(n1)), scale(rnorm(n2)))            # exact mean 0, sd 1 per arm
  if (!balanced) x[1:n1] <- x[1:n1] + 1.0               # standardised imbalance 1.0
  grp <- rep(c("exp", "nexp"), c(n1, n2))
  y <- delta * (grp == "exp") + rho * x + sqrt(1 - rho^2) * rnorm(n1 + n2)

  fit <- lm(y ~ grp + x)
  s_res <- summary(fit)$sigma
  co <- summary(fit)$coefficients
  # coefficient for nexp vs exp -> flip so md is exp - nexp
  md <- -co["grpnexp", "Estimate"]
  md_se <- co["grpnexp", "Std. Error"]
  df_res <- fit$df.residual

  # The package's D = 0 relation. Under exact balance this equals md_se above.
  md_se_D0 <- s_res * sqrt(1 / n1 + 1 / n2)

  # A real paper reports the t from its ACTUAL ANCOVA SE, which includes D. Feeding
  # that is what exposes the documented attenuation; feeding md/md_se_D0 instead
  # would make the statistic routes agree by construction and test nothing.
  t_stat <- md / md_se
  f_stat <- t_stat^2
  etasq <- f_stat / (f_stat + df_res)
  pval <- 2 * pt(-abs(t_stat), df_res)
  # likewise the CI and SE a paper prints come from the real ANCOVA SE
  ci_lo <- md - qt(.975, df_res) * md_se
  ci_up <- md + qt(.975, df_res) * md_se
  # Adjusted (least-squares) means and THEIR OWN standard errors, taken from the fit
  # at the overall covariate mean -- this is what a paper actually prints. Do NOT
  # derive the per-arm SEs from md_se: SE(adj_mean_j) carries the per-arm leverage
  # (x_bar_j - x_bar)^2 / SS_x, which is not the between-arm D, so md_se / sqrt(2)
  # would hand this route the full between-arm leverage and make it look like the
  # statistic routes. That mistake hides the fact that ancova_means_se is attenuated
  # LESS than ancova_t and friends: for equal arms the factor is
  # 1/sqrt(1 + delta_x^2/8) rather than 1/sqrt(1 + delta_x^2/4), i.e. half the D
  # CONTRIBUTION, not half the attenuation. Measured at delta_x = 1, n = 60/60:
  # 0.9419 vs 0.8929. (An earlier version of this comment said "half the magnitude",
  # which would predict 0.9697 and is wrong.)
  nd <- data.frame(grp = c("exp", "nexp"), x = mean(x))
  pr <- predict(fit, newdata = nd, se.fit = TRUE)
  adj_mean_exp  <- pr$fit[1];     adj_mean_nexp  <- pr$fit[2]
  adj_se_exp    <- pr$se.fit[1];  adj_se_nexp    <- pr$se.fit[2]

  cat(sprintf("\n[ANCOVA] n = %d/%d, delta = %.2f, rho = %.2f, %s\n",
              n1, n2, delta, rho,
              if (balanced) "EXACT covariate balance (D = 0)" else "imbalance = 1.0 SD"))
  cat(sprintf("         lm md = %.10f, lm SE = %.10f, D=0 SE = %.10f, ratio = %.10f\n",
              md, md_se, md_se_D0, md_se / md_se_D0))

  a <- list(
    md_sd      = safe(es_from_ancova_md_sd(ancova_md = md, ancova_md_sd = s_res,
                        cov_outcome_r = rho, n_cov_ancova = 1, n_exp = n1, n_nexp = n2)),
    means_sd   = safe(es_from_ancova_means_sd(n_exp = n1, n_nexp = n2,
                        ancova_mean_exp = adj_mean_exp, ancova_mean_nexp = adj_mean_nexp,
                        ancova_mean_sd_exp = s_res, ancova_mean_sd_nexp = s_res,
                        cov_outcome_r = rho, n_cov_ancova = 1)),
    md_se      = safe(es_from_ancova_md_se(ancova_md = md, ancova_md_se = md_se,
                        cov_outcome_r = rho, n_cov_ancova = 1, n_exp = n1, n_nexp = n2)),
    md_ci      = safe(es_from_ancova_md_ci(ancova_md = md, ancova_md_ci_lo = ci_lo,
                        ancova_md_ci_up = ci_up, cov_outcome_r = rho,
                        n_cov_ancova = 1, n_exp = n1, n_nexp = n2)),
    md_pval    = safe(es_from_ancova_md_pval(ancova_md = md, ancova_md_pval = pval,
                        cov_outcome_r = rho, n_cov_ancova = 1, n_exp = n1, n_nexp = n2)),
    ancova_t   = safe(es_from_ancova_t(ancova_t = t_stat, cov_outcome_r = rho,
                        n_cov_ancova = 1, n_exp = n1, n_nexp = n2)),
    ancova_f   = safe(es_from_ancova_f(ancova_f = f_stat, cov_outcome_r = rho,
                        n_cov_ancova = 1, n_exp = n1, n_nexp = n2)),
    etasq_adj  = safe(es_from_etasq_adj(etasq_adj = etasq, n_exp = n1, n_nexp = n2,
                        n_cov_ancova = 1, cov_outcome_r = rho)),
    means_se   = safe(es_from_ancova_means_se(n_exp = n1, n_nexp = n2,
                        ancova_mean_exp = adj_mean_exp, ancova_mean_nexp = adj_mean_nexp,
                        ancova_mean_se_exp = adj_se_exp,
                        ancova_mean_se_nexp = adj_se_nexp,
                        cov_outcome_r = rho, n_cov_ancova = 1))
  )
  a <- lapply(a, function(z) if (inherits(z, "route_error")) NULL else z)
  eq_report(if (balanced) "ANCOVA routes at EXACT balance (all must match)"
            else "ANCOVA routes under imbalance (statistic routes expected to differ)",
            a)
}

## =============================================================================
## PART 2 -- pre/post family, two-group
## =============================================================================
##
## One dataset. Every input shape is derived from the SAME per-arm summary
## statistics, so the routes carry identical information.
## Cross-family comparison uses morris_drm / morris_dz only (see header note 2).

prepost_equivalence <- function(n1 = 50, n2 = 50, r = 0.6,
                                method = "morris_drm", seed = 20260802) {
  set.seed(seed)
  gen <- function(n, d_shift) {
    S <- rbind(c(1, r), c(r, 1))
    z <- MASS::mvrnorm(n, c(0, d_shift), S)
    list(m_pre = mean(z[, 1]), m_post = mean(z[, 2]),
         sd_pre = sd(z[, 1]), sd_post = sd(z[, 2]), n = n)
  }
  A <- gen(n1, 0.6); B <- gen(n2, 0.2)

  # change-score summaries implied by the pre/post summaries and r
  sdc <- function(g) sqrt(g$sd_pre^2 + g$sd_post^2 - 2 * r * g$sd_pre * g$sd_post)
  mc  <- function(g) g$m_post - g$m_pre
  scA <- sdc(A); scB <- sdc(B); mcA <- mc(A); mcB <- mc(B)
  seA <- scA / sqrt(n1); seB <- scB / sqrt(n2)
  tA <- mcA / seA; tB <- mcB / seB
  ciloA <- mcA - qt(.975, n1 - 1) * seA; ciupA <- mcA + qt(.975, n1 - 1) * seA
  ciloB <- mcB - qt(.975, n2 - 1) * seB; ciupB <- mcB + qt(.975, n2 - 1) * seB
  pA <- 2 * pt(-abs(tA), n1 - 1); pB <- 2 * pt(-abs(tB), n2 - 1)

  cat(sprintf("\n[PRE-POST two-group] n = %d/%d, r = %.2f, pre_post_to_smd = '%s'\n",
              n1, n2, r, method))

  a <- list(
    means_sd    = safe(es_from_means_sd_pre_post(
                    mean_pre_exp = A$m_pre, mean_exp = A$m_post,
                    mean_pre_sd_exp = A$sd_pre, mean_sd_exp = A$sd_post,
                    mean_pre_nexp = B$m_pre, mean_nexp = B$m_post,
                    mean_pre_sd_nexp = B$sd_pre, mean_sd_nexp = B$sd_post,
                    n_exp = n1, n_nexp = n2, r_pre_post_exp = r, r_pre_post_nexp = r,
                    pre_post_to_smd = method)),
    means_se    = safe(es_from_means_se_pre_post(
                    mean_pre_exp = A$m_pre, mean_exp = A$m_post,
                    mean_pre_se_exp = A$sd_pre / sqrt(n1), mean_se_exp = A$sd_post / sqrt(n1),
                    mean_pre_nexp = B$m_pre, mean_nexp = B$m_post,
                    mean_pre_se_nexp = B$sd_pre / sqrt(n2), mean_se_nexp = B$sd_post / sqrt(n2),
                    n_exp = n1, n_nexp = n2, r_pre_post_exp = r, r_pre_post_nexp = r,
                    pre_post_to_smd = method)),
    change_sd   = safe(es_from_mean_change_sd(
                    mean_change_exp = mcA, mean_change_sd_exp = scA,
                    mean_change_nexp = mcB, mean_change_sd_nexp = scB,
                    r_pre_post_exp = r, r_pre_post_nexp = r,
                    n_exp = n1, n_nexp = n2, pre_post_to_smd = method)),
    change_se   = safe(es_from_mean_change_se(
                    mean_change_exp = mcA, mean_change_se_exp = seA,
                    mean_change_nexp = mcB, mean_change_se_nexp = seB,
                    r_pre_post_exp = r, r_pre_post_nexp = r,
                    n_exp = n1, n_nexp = n2, pre_post_to_smd = method)),
    change_ci   = safe(es_from_mean_change_ci(
                    mean_change_exp = mcA, mean_change_ci_lo_exp = ciloA,
                    mean_change_ci_up_exp = ciupA, mean_change_nexp = mcB,
                    mean_change_ci_lo_nexp = ciloB, mean_change_ci_up_nexp = ciupB,
                    r_pre_post_exp = r, r_pre_post_nexp = r,
                    n_exp = n1, n_nexp = n2, pre_post_to_smd = method)),
    change_pval = safe(es_from_mean_change_pval(
                    mean_change_exp = mcA, mean_change_pval_exp = pA,
                    mean_change_nexp = mcB, mean_change_pval_nexp = pB,
                    r_pre_post_exp = r, r_pre_post_nexp = r,
                    n_exp = n1, n_nexp = n2, pre_post_to_smd = method)),
    paired_t    = safe(es_from_paired_t(
                    paired_t_exp = tA, paired_t_nexp = tB, n_exp = n1, n_nexp = n2,
                    r_pre_post_exp = r, r_pre_post_nexp = r, pre_post_to_smd = method)),
    paired_tp   = safe(es_from_paired_t_pval(
                    paired_t_pval_exp = pA, paired_t_pval_nexp = pB,
                    n_exp = n1, n_nexp = n2,
                    r_pre_post_exp = r, r_pre_post_nexp = r, pre_post_to_smd = method)),
    paired_f    = safe(es_from_paired_f(
                    paired_f_exp = tA^2, paired_f_nexp = tB^2, n_exp = n1, n_nexp = n2,
                    r_pre_post_exp = r, r_pre_post_nexp = r, pre_post_to_smd = method))
  )
  errs <- vapply(a, function(z) inherits(z, "route_error"), logical(1))
  if (any(errs)) for (k in names(a)[errs])
    cat("   [error] ", k, ": ", a[[k]]$err, "\n", sep = "")
  a <- lapply(a, function(z) if (inherits(z, "route_error")) NULL else z)
  eq_report(paste0("pre/post two-group routes, method = ", method), a)
}

## =============================================================================
run_99 <- function() {
  load_metaconvert()
  cat("\n#########################################################\n")
  cat("## ROUTE EQUIVALENCE -- deterministic wiring check\n")
  cat("#########################################################\n")

  ancova_equivalence(balanced = TRUE)
  ancova_equivalence(balanced = FALSE)

  for (m in c("morris_drm", "morris_dz")) prepost_equivalence(method = m)

  cat("\n\nNOTE ON READING THIS: 'match' means the routes agree to 1e-8.\n")
  cat("'**DIFFERS**' is only a bug if the two routes carry the same information\n")
  cat("AND target the same estimand -- see the header of this file.\n")
  invisible(NULL)
}
