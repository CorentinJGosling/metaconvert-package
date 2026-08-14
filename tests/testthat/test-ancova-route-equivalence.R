## Fidelity of the ANCOVA / eta-squared family to Cooper ch. 12 (eqs. 12.23-12.26,
## table 12.1 row 3 and table 12.3) and Lai & Kelley (2012).
##
## The organising principle of these tests: several routes are fed DIFFERENT summary
## statistics computed from THE SAME raw dataset. Because they all describe the same
## study, they must return the same effect size and the same sampling variance. Any
## divergence is either a formula error or an undocumented estimand change.

# ---------------------------------------------------------------------------
# helper: build one dataset and extract every statistic the routes consume
# ---------------------------------------------------------------------------
.make_ancova_study <- function(n1, n2, tau = 0.8, beta = 1, s_e = 1,
                               dx = 0, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n1); x2 <- rnorm(n2)
  # centre each arm exactly, then impose the wanted covariate imbalance, so that
  # dx = 0 means EXACT balance (D = 0 in Lai & Kelley eq. 5) and not merely balance
  # in expectation. At D = 0 every route is algebraically identical.
  x1 <- x1 - mean(x1) + dx
  x2 <- x2 - mean(x2)
  x <- c(x1, x2)
  grp <- rep(c(1, 0), c(n1, n2))
  y <- tau * grp + beta * x + rnorm(n1 + n2, 0, s_e)

  fit  <- lm(y ~ grp + x)
  fit0 <- lm(y ~ grp)
  sm   <- summary(fit)$coefficients["grp", ]
  s_resid <- summary(fit)$sigma
  s_crude <- summary(fit0)$sigma
  b <- coef(fit); xbar <- mean(x)

  list(
    n1 = n1, n2 = n2, q = 1,
    md = unname(sm[1]), se = unname(sm[2]),
    t = unname(sm[3]), p = unname(sm[4]), f = unname(sm[3])^2,
    ci = confint(fit)["grp", ],
    s_resid = s_resid, s_crude = s_crude,
    # Cooper 12.24 / Lai & Kelley eq. 3: the R for which MSE_ancova = MSW (1 - R^2),
    # i.e. the pooled WITHIN-GROUP correlation -- not the total-sample one.
    r = sqrt(max(0, 1 - s_resid^2 / s_crude^2)),
    adj_exp  = unname(b[1] + b[2] + b[3] * xbar),
    adj_nexp = unname(b[1] + b[3] * xbar),
    etap2 = unname(sm[3])^2 / (unname(sm[3])^2 + n1 + n2 - 2 - 1)
  )
}

.all_ancova_routes <- function(s) {
  list(
    means_sd = es_from_ancova_means_sd(
      n_exp = s$n1, n_nexp = s$n2, ancova_mean_exp = s$adj_exp,
      ancova_mean_nexp = s$adj_nexp, ancova_mean_sd_exp = s$s_resid,
      ancova_mean_sd_nexp = s$s_resid, cov_outcome_r = s$r, n_cov_ancova = s$q),
    means_se = es_from_ancova_means_se(
      n_exp = s$n1, n_nexp = s$n2, ancova_mean_exp = s$adj_exp,
      ancova_mean_nexp = s$adj_nexp, ancova_mean_se_exp = s$s_resid / sqrt(s$n1),
      ancova_mean_se_nexp = s$s_resid / sqrt(s$n2), cov_outcome_r = s$r,
      n_cov_ancova = s$q),
    pooled_adj = es_from_ancova_means_sd_pooled_adj(
      ancova_mean_exp = s$adj_exp, ancova_mean_nexp = s$adj_nexp,
      ancova_mean_sd_pooled = s$s_resid, cov_outcome_r = s$r,
      n_cov_ancova = s$q, n_exp = s$n1, n_nexp = s$n2),
    pooled_crude = es_from_ancova_means_sd_pooled_crude(
      ancova_mean_exp = s$adj_exp, ancova_mean_nexp = s$adj_nexp,
      mean_sd_pooled = s$s_resid / sqrt(1 - s$r^2), cov_outcome_r = s$r,
      n_cov_ancova = s$q, n_exp = s$n1, n_nexp = s$n2),
    md_sd = es_from_ancova_md_sd(
      ancova_md = s$md, ancova_md_sd = s$s_resid, cov_outcome_r = s$r,
      n_cov_ancova = s$q, n_exp = s$n1, n_nexp = s$n2),
    md_se = es_from_ancova_md_se(
      ancova_md = s$md, ancova_md_se = s$se, cov_outcome_r = s$r,
      n_cov_ancova = s$q, n_exp = s$n1, n_nexp = s$n2),
    md_ci = es_from_ancova_md_ci(
      ancova_md = s$md, ancova_md_ci_lo = s$ci[1], ancova_md_ci_up = s$ci[2],
      cov_outcome_r = s$r, n_cov_ancova = s$q, n_exp = s$n1, n_nexp = s$n2),
    md_pval = es_from_ancova_md_pval(
      ancova_md = s$md, ancova_md_pval = s$p, cov_outcome_r = s$r,
      n_cov_ancova = s$q, n_exp = s$n1, n_nexp = s$n2),
    anc_t = es_from_ancova_t(
      ancova_t = s$t, cov_outcome_r = s$r, n_cov_ancova = s$q,
      n_exp = s$n1, n_nexp = s$n2),
    anc_f = es_from_ancova_f(
      ancova_f = s$f, cov_outcome_r = s$r, n_cov_ancova = s$q,
      n_exp = s$n1, n_nexp = s$n2),
    anc_t_pval = es_from_ancova_t_pval(
      ancova_t_pval = s$p, cov_outcome_r = s$r, n_cov_ancova = s$q,
      n_exp = s$n1, n_nexp = s$n2),
    anc_f_pval = es_from_ancova_f_pval(
      ancova_f_pval = s$p, cov_outcome_r = s$r, n_cov_ancova = s$q,
      n_exp = s$n1, n_nexp = s$n2),
    etasq_adj = es_from_etasq_adj(
      etasq_adj = s$etap2, n_exp = s$n1, n_nexp = s$n2, n_cov_ancova = s$q,
      cov_outcome_r = s$r)
  )
}

## ---------------------------------------------------------------------------
## 1. At EXACT covariate balance every ANCOVA route must return the same d, SE and g
## ---------------------------------------------------------------------------
test_that("all ANCOVA routes agree at exact covariate balance", {
  for (nn in list(c(60, 60), c(20, 60), c(15, 15))) {
    s <- .make_ancova_study(nn[1], nn[2], dx = 0, seed = 42)
    routes <- .all_ancova_routes(s)
    ref <- routes$md_sd

    for (nm in names(routes)) {
      expect_equal(routes[[nm]]$d,    ref$d,    tolerance = 1e-8,
                   info = paste("d mismatch on route", nm, "at n =", nn[1], "/", nn[2]))
      expect_equal(routes[[nm]]$d_se, ref$d_se, tolerance = 1e-8,
                   info = paste("d_se mismatch on route", nm))
      expect_equal(routes[[nm]]$g,    ref$g,    tolerance = 1e-8,
                   info = paste("g mismatch on route", nm))
      expect_equal(routes[[nm]]$g_se, ref$g_se, tolerance = 1e-8,
                   info = paste("g_se mismatch on route", nm))
      # the converted measures must agree too: r/z are deterministic functions of d
      # and the arm sizes, and their variances inherit vd (see
      # test-smd-to-cor-variance.R)
      expect_equal(routes[[nm]]$r,    ref$r,    tolerance = 1e-8,
                   info = paste("r mismatch on route", nm))
      expect_equal(routes[[nm]]$r_se, ref$r_se, tolerance = 1e-8,
                   info = paste("r_se mismatch on route", nm))
      expect_equal(routes[[nm]]$z_se, ref$z_se, tolerance = 1e-8,
                   info = paste("z_se mismatch on route", nm))
      expect_equal(routes[[nm]]$logor_se, ref$logor_se, tolerance = 1e-8,
                   info = paste("logor_se mismatch on route", nm))
    }
    # the MD-returning routes must also agree on the raw scale
    for (nm in c("means_sd", "means_se", "pooled_adj", "pooled_crude",
                 "md_sd", "md_se", "md_ci", "md_pval")) {
      expect_equal(routes[[nm]]$md,    ref$md,    tolerance = 1e-8, info = nm)
      expect_equal(routes[[nm]]$md_se, ref$md_se, tolerance = 1e-8, info = nm)
    }
  }
})

## ---------------------------------------------------------------------------
## 2. Cooper's own formulas, reproduced from the book
## ---------------------------------------------------------------------------
test_that("Cooper eq. 12.23/12.24 worked example (p. 230) is reproduced", {
  # Ybar1 = 103, Ybar2 = 100, S_Adjusted = 5.5, n1 = n2 = 50, R = 0.7, q = 1
  # Book: S_Within = 7.7015, d = 0.3895
  res <- es_from_ancova_md_sd(ancova_md = 3, ancova_md_sd = 5.5, cov_outcome_r = 0.7,
                              n_cov_ancova = 1, n_exp = 50, n_nexp = 50)
  expect_equal(5.5 / sqrt(1 - 0.7^2), 7.7015, tolerance = 1e-4)
  expect_equal(res$d, 0.3895, tolerance = 1e-4)
})

test_that("Cooper eq. 12.26 is the sampling variance used on every ANCOVA route", {
  s <- .make_ancova_study(60, 60, dx = 0, seed = 7)
  for (route in .all_ancova_routes(s)) {
    v_cooper <- (s$n1 + s$n2) * (1 - s$r^2) / (s$n1 * s$n2) +
      route$d^2 / (2 * (s$n1 + s$n2))
    expect_equal(route$d_se^2, v_cooper, tolerance = 1e-12)
  }
})

test_that("Hedges' J uses the ANCOVA error df (N - 2 - q)", {
  for (q in c(0, 1, 5)) {
    r <- es_from_ancova_md_sd(ancova_md = 2, ancova_md_sd = 3, cov_outcome_r = 0.4,
                              n_cov_ancova = q, n_exp = 30, n_nexp = 30)
    J <- metaConvert:::.d_j(30 + 30 - 2 - q)
    expect_equal(r$g, r$d * J, tolerance = 1e-12)
    expect_equal(r$g_se, r$d_se * J, tolerance = 1e-12)
  }
})

## ---------------------------------------------------------------------------
## 3. eta-squared -> d must round-trip against the raw data it came from
##    (regression test for the equal-n closed form that ignored both arm sizes)
## ---------------------------------------------------------------------------
test_that("es_from_etasq reproduces the raw-data Cohen's d for unequal groups", {
  set.seed(2024)
  for (nn in list(c(50, 50), c(20, 60), c(10, 90), c(40, 10), c(15, 85))) {
    n1 <- nn[1]; n2 <- nn[2]
    y1 <- rnorm(n1, 0.6, 1); y2 <- rnorm(n2, 0, 1)
    y <- c(y1, y2); g <- factor(rep(c("e", "c"), c(n1, n2)))

    av   <- anova(lm(y ~ g))
    eta2 <- av[["Sum Sq"]][1] / (av[["Sum Sq"]][1] + av[["Sum Sq"]][2])

    s_pooled <- sqrt(((n1 - 1) * var(y1) + (n2 - 1) * var(y2)) / (n1 + n2 - 2))
    d_raw <- (mean(y1) - mean(y2)) / s_pooled

    lab <- paste("n =", n1, "/", n2)
    expect_equal(es_from_etasq(eta2, n_exp = n1, n_nexp = n2)$d, d_raw,
                 tolerance = 1e-8, info = lab)
    # and the three routes carrying the same information must coincide
    expect_equal(es_from_etasq(eta2, n_exp = n1, n_nexp = n2)$d,
                 es_from_anova_f(av[["F value"]][1], n_exp = n1, n_nexp = n2)$d,
                 tolerance = 1e-10, info = lab)
    expect_equal(es_from_etasq(eta2, n_exp = n1, n_nexp = n2)$d,
                 es_from_pt_bis_r(sqrt(eta2), n_exp = n1, n_nexp = n2)$d,
                 tolerance = 1e-10, info = lab)
  }
})

test_that("es_from_etasq depends on the group-size split (it must not be constant)", {
  e2 <- 0.09
  d_bal   <- es_from_etasq(e2, n_exp = 50, n_nexp = 50)$d
  d_unbal <- es_from_etasq(e2, n_exp = 10, n_nexp = 90)$d
  expect_false(isTRUE(all.equal(d_bal, d_unbal)))
  expect_gt(d_unbal, d_bal)
})

test_that("es_from_etasq_adj still matches es_from_ancova_f exactly", {
  for (nn in list(c(20, 20), c(20, 60), c(10, 90))) {
    df_err <- nn[1] + nn[2] - 2 - 2
    e2 <- 0.2
    expect_equal(
      es_from_etasq_adj(etasq_adj = e2, n_exp = nn[1], n_nexp = nn[2],
                        n_cov_ancova = 2, cov_outcome_r = 0.5)$d,
      es_from_ancova_f(ancova_f = e2 * df_err / (1 - e2), cov_outcome_r = 0.5,
                       n_cov_ancova = 2, n_exp = nn[1], n_nexp = nn[2])$d,
      tolerance = 1e-12)
  }
})

## ---------------------------------------------------------------------------
## 4. Hierarchy: the imbalance-immune route must outrank the attenuating ones
## ---------------------------------------------------------------------------
test_that("convert_df prefers ancova_md_sd over the ANCOVA test statistics", {
  s <- .make_ancova_study(60, 60, dx = 1, seed = 3)   # deliberate covariate imbalance
  d <- data.frame(
    study_id = "s1", n_exp = s$n1, n_nexp = s$n2,
    n_cov_ancova = s$q, cov_outcome_r = s$r,
    ancova_md = s$md, ancova_md_sd = s$s_resid,
    ancova_t = s$t, ancova_f = s$f
  )
  out <- summary(convert_df(d, measure = "g", verbose = FALSE))
  expect_equal(as.character(out$info_used_adjusted), "ancova_md_sd")

  # the selected value must be the unattenuated one
  ref <- es_from_ancova_md_sd(ancova_md = s$md, ancova_md_sd = s$s_resid,
                              cov_outcome_r = s$r, n_cov_ancova = s$q,
                              n_exp = s$n1, n_nexp = s$n2)
  expect_equal(out$es_adjusted, round(ref$g, 3), tolerance = 1e-3)
})

test_that("ancova_md_se/_ci/_pval stay BELOW the test statistics (same attenuation tier)", {
  s <- .make_ancova_study(60, 60, dx = 1, seed = 3)
  d <- data.frame(
    study_id = "s1", n_exp = s$n1, n_nexp = s$n2,
    n_cov_ancova = s$q, cov_outcome_r = s$r,
    ancova_md = s$md, ancova_md_se = s$se, ancova_f = s$f
  )
  out <- summary(convert_df(d, measure = "g", verbose = FALSE))
  expect_equal(as.character(out$info_used_adjusted), "ancova_f")
})

test_that("the adjusted hierarchy still exposes exactly the expected method count", {
  # guards the L20a / L20b split against dropping or duplicating a route
  d <- data.frame(study_id = "s1", n_exp = 30, n_nexp = 30,
                  ancova_md = 2, ancova_md_sd = 3,
                  cov_outcome_r = 0.4, n_cov_ancova = 1)
  expect_no_error(convert_df(d, measure = "g", verbose = FALSE))
})

## ---------------------------------------------------------------------------
## 5. An omitted cov_outcome_r must not be silently invented
## ---------------------------------------------------------------------------
test_that("es_from_cohen_d_adj returns NA rather than assuming cov_outcome_r = 0.5", {
  omitted  <- es_from_cohen_d_adj(cohen_d_adj = 1, n_cov_ancova = 4,
                                  n_exp = 20, n_nexp = 20)
  supplied <- es_from_cohen_d_adj(cohen_d_adj = 1, n_cov_ancova = 4,
                                  cov_outcome_r = 0.5, n_exp = 20, n_nexp = 20)
  expect_true(is.na(omitted$d_se))
  expect_false(is.na(supplied$d_se))
  # the point estimate itself does not depend on R, so it is still returned
  expect_equal(omitted$d, 1)
  # ...and the blanking must be COHERENT: everything variance-bearing goes with it,
  # rather than leaving a hole in one measure while its siblings look fine
  for (nm in c("d_se", "g", "g_se", "r", "r_se", "z", "z_se", "logor_se"))
    expect_true(is.na(omitted[[nm]]), info = paste("expected NA:", nm))
})

test_that("an omitted n_cov_ancova defaults to q = 0 and returns a complete row", {
  # n_cov_ancova only costs degrees of freedom, so unlike cov_outcome_r it HAS a safe
  # default. Regression guard: it must not blank g while leaving d/d_se/r finite.
  res <- es_from_cohen_d_adj(cohen_d_adj = 0.5, n_exp = 50, n_nexp = 50,
                             cov_outcome_r = 0.5)
  for (nm in c("d", "d_se", "g", "g_se", "r", "r_se"))
    expect_false(is.na(res[[nm]]), info = paste("unexpected NA:", nm))
  expect_equal(res$g, res$d * metaConvert:::.d_j(50 + 50 - 2), tolerance = 1e-12)
})

test_that("the crude branch of .es_from_d is unaffected by the NA fallback", {
  crude <- es_from_cohen_d(cohen_d = 1, n_exp = 20, n_nexp = 20)
  expect_equal(crude$d_se, sqrt(40 / 400 + 1 / 80), tolerance = 1e-12)
  expect_false(is.na(crude$g_se))
})

## ---------------------------------------------------------------------------
## 6. cov_outcome_r bounds
## ---------------------------------------------------------------------------
test_that("out-of-range cov_outcome_r is flagged and set to missing", {
  d <- data.frame(study_id = c("a", "b", "c"), n_exp = 30, n_nexp = 30,
                  ancova_mean_exp = 10, ancova_mean_nexp = 8,
                  ancova_mean_sd_exp = 4, ancova_mean_sd_nexp = 4,
                  n_cov_ancova = 1, cov_outcome_r = c(0.5, 1, 1.2))
  out <- summary(convert_df(d, measure = "g", verbose = FALSE), flags = TRUE)

  expect_false(is.na(out$es_adjusted[1]))          # valid row survives
  expect_true(is.na(out$es_adjusted[2]))           # |R| = 1 removed
  expect_true(is.na(out$es_adjusted[3]))           # |R| > 1 removed
  expect_match(out$flags_adjusted[2], "Degenerate covariate-outcome correlation")
  expect_match(out$flags_adjusted[3], "Out-of-range covariate-outcome correlation")
  expect_false(grepl("Degenerate|Out-of-range covariate", out$flags_adjusted[1]))
})

test_that("cov_outcome_r = 1 no longer yields a zero-SE effect size", {
  d <- data.frame(study_id = "a", n_exp = 50, n_nexp = 50,
                  ancova_mean_exp = 10, ancova_mean_nexp = 8,
                  mean_sd_pooled = 4, n_cov_ancova = 1, cov_outcome_r = 1)
  out <- summary(convert_df(d, measure = "g", verbose = FALSE), flags = TRUE)
  # Previously this route returned a plausible ES with a silently deflated SE (2.13x)
  # and NO flag, so "se > 0" was already true then -- assert the actual contract:
  # V34 removes the row and names the offending input.
  expect_true(is.na(out$es_adjusted))
  expect_true(is.na(out$se_adjusted))
  expect_match(out$flags_adjusted, "Degenerate covariate-outcome correlation")
})

test_that("negative cov_outcome_r is valid (only R^2 is used)", {
  pos <- es_from_ancova_t(ancova_t = 2, cov_outcome_r = 0.7, n_cov_ancova = 1,
                          n_exp = 30, n_nexp = 30)
  neg <- es_from_ancova_t(ancova_t = 2, cov_outcome_r = -0.7, n_cov_ancova = 1,
                          n_exp = 30, n_nexp = 30)
  expect_equal(pos$d, neg$d)
  expect_equal(pos$d_se, neg$d_se)
})

## ---------------------------------------------------------------------------
## 7. Documented Lai & Kelley attenuation, as an exact closed form
## ---------------------------------------------------------------------------
test_that("covariate imbalance attenuates the test-statistic routes by the documented factor", {
  n1 <- n2 <- 60
  set.seed(7)
  x1 <- rnorm(n1); x2 <- rnorm(n2)
  x1 <- x1 - mean(x1) + 1; x2 <- x2 - mean(x2)   # standardised imbalance = 1
  x <- c(x1, x2); grp <- rep(c(1, 0), c(n1, n2))
  y <- 0.8 * grp + x + rnorm(120)

  fit <- lm(y ~ grp + x); fit0 <- lm(y ~ grp)
  s_resid <- summary(fit)$sigma
  r <- sqrt(1 - s_resid^2 / summary(fit0)$sigma^2)
  sm <- summary(fit)$coefficients["grp", ]
  SSx <- sum((x1 - mean(x1))^2) + sum((x2 - mean(x2))^2)
  D <- (mean(x1) - mean(x2))^2 / SSx              # Lai & Kelley eq. (5)

  ref <- es_from_ancova_md_sd(ancova_md = sm[1], ancova_md_sd = s_resid,
                              cov_outcome_r = r, n_cov_ancova = 1,
                              n_exp = n1, n_nexp = n2)$d
  d_t <- es_from_ancova_t(ancova_t = sm[3], cov_outcome_r = r, n_cov_ancova = 1,
                          n_exp = n1, n_nexp = n2)$d
  expect_equal(d_t / ref, 1 / sqrt(1 + D / (1 / n1 + 1 / n2)), tolerance = 1e-9)

  # the per-arm SE route sees only its own leverage -> the half-D tier
  xbar <- mean(x)
  se_j <- function(nj, xbj) s_resid * sqrt(1 / nj + (xbj - xbar)^2 / SSx)
  b <- coef(fit)
  d_mse <- es_from_ancova_means_se(
    n_exp = n1, n_nexp = n2,
    ancova_mean_exp  = unname(b[1] + b[2] + b[3] * xbar),
    ancova_mean_nexp = unname(b[1] + b[3] * xbar),
    ancova_mean_se_exp  = se_j(n1, mean(x1)),
    ancova_mean_se_nexp = se_j(n2, mean(x2)),
    cov_outcome_r = r, n_cov_ancova = 1)$d
  expect_equal(d_mse / ref, 1 / sqrt(1 + n1 * D / 4), tolerance = 1e-9)
})

test_that("ancova_md_sd stays unbiased under imbalance while ancova_t attenuates", {
  # Unbiasedness is a property of the EXPECTATION, not of a single sample: under
  # covariate imbalance grp and x become collinear, which inflates the sampling error
  # of the adjusted MD without biasing it. So this has to be a Monte Carlo check.
  # Deterministic via the fixed seed.
  skip_on_cran()
  reps <- 400
  delta_true <- 0.8 / sqrt(1^2 * 1 + 1^2)   # tau / sqrt(beta^2 var(x within) + sigma^2)

  mc <- function(dx) {
    set.seed(99)
    out <- vapply(seq_len(reps), function(i) {
      s <- .make_ancova_study(60, 60, dx = dx, seed = sample.int(1e6, 1))
      c(md_sd = es_from_ancova_md_sd(ancova_md = s$md, ancova_md_sd = s$s_resid,
                                     cov_outcome_r = s$r, n_cov_ancova = 1,
                                     n_exp = 60, n_nexp = 60)$d,
        anc_t = es_from_ancova_t(ancova_t = s$t, cov_outcome_r = s$r,
                                 n_cov_ancova = 1, n_exp = 60, n_nexp = 60)$d)
    }, numeric(2))
    rowMeans(out)
  }

  at_bal   <- mc(0)
  at_unbal <- mc(1)

  # the adjusted-MD + residual-SD route is imbalance-immune in expectation
  expect_lt(abs(at_bal[["md_sd"]]   / delta_true - 1), 0.05)
  expect_lt(abs(at_unbal[["md_sd"]] / delta_true - 1), 0.05)

  # the test-statistic route is NOT: it attenuates toward zero, as documented
  expect_lt(abs(at_bal[["anc_t"]] / delta_true - 1), 0.05)
  expect_lt(at_unbal[["anc_t"]] / at_unbal[["md_sd"]], 0.95)
})
