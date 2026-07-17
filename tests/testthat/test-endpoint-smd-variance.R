# Bit-exactness of the endpoint SMD variance convention (smd_var) and the
# control-SD standardiser (smd_denom, Glass's delta) against metafor::escalc().
#
# The package's poolable measure is Hedges' g, so every check compares the
# reported g variance (g_se^2) against metafor's vi. The point estimate g is
# compared against metafor's yi. These oracles are metafor's own formulas, so
# they are NOT circular re-statements of the package's arithmetic.

skip_if_not_installed("metafor")
library(metafor)

TOL <- 1e-10

# ---- scenarios (unbalanced + balanced, small + large) -----------------------
SCEN <- list(
  list(n1 = 20, n2 = 22, m1 = 2.3, m2 = 1.9, s1 = 1.2, s2 = 0.9),
  list(n1 = 10, n2 = 15, m1 = 5.1, m2 = 4.4, s1 = 2.0, s2 = 1.8),
  list(n1 = 40, n2 = 30, m1 = 0.8, m2 = 1.5, s1 = 3.1, s2 = 2.4),
  list(n1 = 30, n2 = 30, m1 = 10,  m2 = 8,   s1 = 5,   s2 = 6)
)

# =============================================================================
# Task A - pooled-SD SMD variance: LS2 (default) and LS
# =============================================================================
test_that("es_from_means_sd g variance is bit-exact with metafor SMD LS2 (default) and LS", {
  for (s in SCEN) {
    o_ls2 <- escalc(measure = "SMD", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
                    n1i = s$n1, n2i = s$n2, vtype = "LS2")
    o_ls  <- escalc(measure = "SMD", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
                    n1i = s$n1, n2i = s$n2, vtype = "LS")

    p_def <- es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1,
                              mean_nexp = s$m2, mean_sd_nexp = s$s2,
                              n_exp = s$n1, n_nexp = s$n2)                 # default = LS2
    p_ls  <- es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1,
                              mean_nexp = s$m2, mean_sd_nexp = s$s2,
                              n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin")

    expect_equal(p_def$g,       as.numeric(o_ls2$yi), tolerance = TOL)
    expect_equal(p_def$g_se^2,  as.numeric(o_ls2$vi), tolerance = TOL)   # default reproduces LS2
    expect_equal(p_ls$g_se^2,   as.numeric(o_ls$vi),  tolerance = TOL)   # opt-in reproduces LS
  }
})

test_that("means_sd_pooled, means_se and means_ci honour smd_var", {
  s <- SCEN[[1]]
  # pooled SD path: feed the same SD to both arms so sp == that SD
  sp <- sqrt(((s$n1 - 1) * s$s1^2 + (s$n2 - 1) * s$s2^2) / (s$n1 + s$n2 - 2))
  o_ls2 <- escalc(measure = "SMD", m1i = s$m1, m2i = s$m2, sd1i = sp, sd2i = sp,
                  n1i = s$n1, n2i = s$n2, vtype = "LS2")
  o_ls  <- escalc(measure = "SMD", m1i = s$m1, m2i = s$m2, sd1i = sp, sd2i = sp,
                  n1i = s$n1, n2i = s$n2, vtype = "LS")
  pp_ls2 <- es_from_means_sd_pooled(mean_exp = s$m1, mean_nexp = s$m2, mean_sd_pooled = sp,
                                    n_exp = s$n1, n_nexp = s$n2)
  pp_ls  <- es_from_means_sd_pooled(mean_exp = s$m1, mean_nexp = s$m2, mean_sd_pooled = sp,
                                    n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin")
  expect_equal(pp_ls2$g_se^2, as.numeric(o_ls2$vi), tolerance = TOL)
  expect_equal(pp_ls$g_se^2,  as.numeric(o_ls$vi),  tolerance = TOL)

  # means_se / means_ci: build SE and CI inputs equivalent to means_sd
  o2_ls2 <- escalc(measure = "SMD", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
                   n1i = s$n1, n2i = s$n2, vtype = "LS2")
  se <- es_from_means_se(mean_exp = s$m1, mean_se_exp = s$s1 / sqrt(s$n1),
                         mean_nexp = s$m2, mean_se_nexp = s$s2 / sqrt(s$n2),
                         n_exp = s$n1, n_nexp = s$n2)
  expect_equal(se$g_se^2, as.numeric(o2_ls2$vi), tolerance = TOL)
  se_ls <- es_from_means_se(mean_exp = s$m1, mean_se_exp = s$s1 / sqrt(s$n1),
                            mean_nexp = s$m2, mean_se_nexp = s$s2 / sqrt(s$n2),
                            n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin")
  o2_ls <- escalc(measure = "SMD", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
                  n1i = s$n1, n2i = s$n2, vtype = "LS")
  expect_equal(se_ls$g_se^2, as.numeric(o2_ls$vi), tolerance = TOL)
})

test_that("statistic-based d/g producers honour smd_var (oracle = escalc di = the d used)", {
  # For each producer, take the uncorrected d it computed and check g_se^2 equals
  # metafor's SMD vi built from that same di - this isolates the variance convention.
  s <- SCEN[[1]]
  producers <- list(
    cohen_d   = es_from_cohen_d(cohen_d = 0.5, n_exp = s$n1, n_nexp = s$n2),
    hedges_g  = es_from_hedges_g(hedges_g = 0.5, n_exp = s$n1, n_nexp = s$n2),
    student_t = es_from_student_t(student_t = 2.1, n_exp = s$n1, n_nexp = s$n2),
    anova_f   = es_from_anova_f(anova_f = 4.2, n_exp = s$n1, n_nexp = s$n2),
    etasq     = es_from_etasq(etasq = 0.12, n_exp = s$n1, n_nexp = s$n2),
    pt_bis_r  = es_from_pt_bis_r(pt_bis_r = 0.3, n_exp = s$n1, n_nexp = s$n2),
    beta      = es_from_beta_unstd(beta_unstd = 0.7, sd_dv = 1.4, n_exp = s$n1, n_nexp = s$n2)
  )
  producers_ls <- list(
    cohen_d   = es_from_cohen_d(cohen_d = 0.5, n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin"),
    hedges_g  = es_from_hedges_g(hedges_g = 0.5, n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin"),
    student_t = es_from_student_t(student_t = 2.1, n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin"),
    anova_f   = es_from_anova_f(anova_f = 4.2, n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin"),
    etasq     = es_from_etasq(etasq = 0.12, n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin"),
    pt_bis_r  = es_from_pt_bis_r(pt_bis_r = 0.3, n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin"),
    beta      = es_from_beta_unstd(beta_unstd = 0.7, sd_dv = 1.4, n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin")
  )
  for (nm in names(producers)) {
    p  <- producers[[nm]]
    pl <- producers_ls[[nm]]
    o_ls2 <- escalc(measure = "SMD", di = p$d, n1i = s$n1, n2i = s$n2, vtype = "LS2")
    o_ls  <- escalc(measure = "SMD", di = p$d, n1i = s$n1, n2i = s$n2, vtype = "LS")
    expect_equal(p$g_se^2,  as.numeric(o_ls2$vi), tolerance = TOL, info = paste(nm, "LS2"))
    expect_equal(pl$g_se^2, as.numeric(o_ls$vi),  tolerance = TOL, info = paste(nm, "LS"))
  }
})

test_that("smd_var aliases resolve and the metafor codes are no longer accepted", {
  s <- SCEN[[2]]
  base <- function(v) es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1,
                                       mean_nexp = s$m2, mean_sd_nexp = s$s2,
                                       n_exp = s$n1, n_nexp = s$n2, smd_var = v)$g_se
  expect_equal(base("viechtbauer"), base("hedges_olkin"))
  expect_error(base("LS2"), "smd_var")   # metafor vtype codes are no longer user-facing
  expect_error(base("LS"),  "smd_var")
  expect_error(base("nonsense"), "smd_var")
})

test_that("smd_denom accepts the 'glass' / 'glass_robust' aliases", {
  s <- SCEN[[1]]
  gse <- function(v) es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1,
                                      mean_nexp = s$m2, mean_sd_nexp = s$s2,
                                      n_exp = s$n1, n_nexp = s$n2, smd_denom = v)$g_se
  expect_equal(gse("glass"),        gse("control"))
  expect_equal(gse("glass_robust"), gse("control_robust"))
  expect_error(gse("nope"), "smd_denom")
})

test_that("adjusted (ANCOVA) SMD reduces to metafor crude SMD when there is no covariate", {
  # With cov_outcome_r = 0 and n_cov_ancova = 0 the adjusted leading term (1 - r^2)
  # collapses to the crude one, so es_from_cohen_d_adj must match crude metafor SMD.
  s <- SCEN[[3]]
  adj_ls2 <- es_from_cohen_d_adj(cohen_d_adj = 0.5, n_cov_ancova = 0, cov_outcome_r = 0,
                                 n_exp = s$n1, n_nexp = s$n2)
  adj_ls  <- es_from_cohen_d_adj(cohen_d_adj = 0.5, n_cov_ancova = 0, cov_outcome_r = 0,
                                 n_exp = s$n1, n_nexp = s$n2, smd_var = "hedges_olkin")
  o_ls2 <- escalc(measure = "SMD", di = 0.5, n1i = s$n1, n2i = s$n2, vtype = "LS2")
  o_ls  <- escalc(measure = "SMD", di = 0.5, n1i = s$n1, n2i = s$n2, vtype = "LS")
  expect_equal(adj_ls2$g_se^2, as.numeric(o_ls2$vi), tolerance = TOL)
  expect_equal(adj_ls$g_se^2,  as.numeric(o_ls$vi),  tolerance = TOL)
})

# =============================================================================
# Task B - Glass's delta (control-SD standardiser): SMD1 / SMD1H
# =============================================================================
test_that("smd_denom='control' is bit-exact with metafor SMD1 (LS2 default and LS)", {
  for (s in SCEN) {
    gc     <- es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1,
                               mean_nexp = s$m2, mean_sd_nexp = s$s2,
                               n_exp = s$n1, n_nexp = s$n2, smd_denom = "control")
    gc_ls  <- es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1,
                               mean_nexp = s$m2, mean_sd_nexp = s$s2,
                               n_exp = s$n1, n_nexp = s$n2, smd_denom = "control", smd_var = "hedges_olkin")
    o_ls2 <- escalc(measure = "SMD1", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
                    n1i = s$n1, n2i = s$n2, vtype = "LS2")
    o_ls  <- escalc(measure = "SMD1", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
                    n1i = s$n1, n2i = s$n2, vtype = "LS")
    expect_equal(gc$g,        as.numeric(o_ls2$yi), tolerance = TOL)
    expect_equal(gc$g_se^2,   as.numeric(o_ls2$vi), tolerance = TOL)
    expect_equal(gc_ls$g_se^2, as.numeric(o_ls$vi), tolerance = TOL)
  }
})

test_that("smd_denom='control' honoured through means_se and means_ci (SMD1)", {
  s <- SCEN[[3]]
  o <- escalc(measure = "SMD1", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
              n1i = s$n1, n2i = s$n2, vtype = "LS2")
  se <- es_from_means_se(mean_exp = s$m1, mean_se_exp = s$s1 / sqrt(s$n1),
                         mean_nexp = s$m2, mean_se_nexp = s$s2 / sqrt(s$n2),
                         n_exp = s$n1, n_nexp = s$n2, smd_denom = "control")
  ci_lo_e <- s$m1 - qt(.975, s$n1 - 1) * s$s1 / sqrt(s$n1)
  ci_up_e <- s$m1 + qt(.975, s$n1 - 1) * s$s1 / sqrt(s$n1)
  ci_lo_n <- s$m2 - qt(.975, s$n2 - 1) * s$s2 / sqrt(s$n2)
  ci_up_n <- s$m2 + qt(.975, s$n2 - 1) * s$s2 / sqrt(s$n2)
  ci <- es_from_means_ci(mean_exp = s$m1, mean_ci_lo_exp = ci_lo_e, mean_ci_up_exp = ci_up_e,
                         mean_nexp = s$m2, mean_ci_lo_nexp = ci_lo_n, mean_ci_up_nexp = ci_up_n,
                         n_exp = s$n1, n_nexp = s$n2, smd_denom = "control")
  expect_equal(se$g_se^2, as.numeric(o$vi), tolerance = TOL)
  expect_equal(ci$g_se^2, as.numeric(o$vi), tolerance = TOL)
})

test_that("adjusted SMD (es_from_cohen_d_adj / es_from_etasq_adj) applies the covariate adjustment", {
  # Regression guard: these functions previously dropped `adjusted = TRUE`, so cov_outcome_r
  # was ignored and the SE stayed crude (~40% too wide at r = 0.7).
  n1 <- 30; n2 <- 25; d <- 0.6; rr <- 0.7; ncov <- 3
  adj <- es_from_cohen_d_adj(cohen_d_adj = d, n_cov_ancova = ncov, cov_outcome_r = rr,
                             n_exp = n1, n_nexp = n2)
  lead <- (n1 + n2) / (n1 * n2) * (1 - rr^2)
  vd   <- lead + d^2 / (2 * (n1 + n2))
  df   <- n1 + n2 - 2 - ncov
  J    <- exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
  expect_equal(adj$d_se,     sqrt(vd),  tolerance = TOL)   # Cooper Table 12.3, adjusted
  expect_equal(adj$g_se^2,   J^2 * vd,  tolerance = TOL)   # LS2
  expect_equal(adj$g / adj$d, J,        tolerance = TOL)   # J on df = N - 2 - n_cov

  # LS (hedges_olkin) adjusted stays internally consistent
  adjL <- es_from_cohen_d_adj(cohen_d_adj = d, n_cov_ancova = ncov, cov_outcome_r = rr,
                              n_exp = n1, n_nexp = n2, smd_var = "hedges_olkin")
  expect_equal(adjL$g_se^2, lead + (J * d)^2 / (2 * (n1 + n2)), tolerance = TOL)

  # etasq_adj must now depend on the covariate
  ea0 <- es_from_etasq_adj(etasq_adj = 0.28, n_cov_ancova = 0, cov_outcome_r = 0,
                           n_exp = 20, n_nexp = 22)
  ea1 <- es_from_etasq_adj(etasq_adj = 0.28, n_cov_ancova = 3, cov_outcome_r = 0.6,
                           n_exp = 20, n_nexp = 22)
  expect_false(isTRUE(all.equal(ea0$d_se, ea1$d_se)))
})

test_that("smd_denom='control_robust' is bit-exact with metafor SMD1H", {
  for (s in SCEN) {
    gcr <- es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1,
                            mean_nexp = s$m2, mean_sd_nexp = s$s2,
                            n_exp = s$n1, n_nexp = s$n2, smd_denom = "control_robust")
    o_h <- escalc(measure = "SMD1H", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
                  n1i = s$n1, n2i = s$n2)
    expect_equal(gcr$g,      as.numeric(o_h$yi), tolerance = TOL)
    expect_equal(gcr$g_se^2, as.numeric(o_h$vi), tolerance = TOL)
  }
})

test_that("smd_denom is honoured PER ROW (a mixed column must not collapse to row 1)", {
  # Regression guard: previously a per-row smd_denom column silently took the first
  # row's value (match.arg(...[1])) and applied it to every row.
  d <- es_from_means_sd(
    mean_exp = c(2.3, 5.1), mean_sd_exp = c(1.2, 2.0),
    mean_nexp = c(1.9, 4.4), mean_sd_nexp = c(0.9, 1.8),
    n_exp = c(20, 40), n_nexp = c(22, 30),
    smd_denom = c("pooled", "control")
  )
  o1 <- escalc(measure = "SMD",  m1i = 2.3, m2i = 1.9, sd1i = 1.2, sd2i = 0.9,
               n1i = 20, n2i = 22, vtype = "LS2")
  o2 <- escalc(measure = "SMD1", m1i = 5.1, m2i = 4.4, sd1i = 2.0, sd2i = 1.8,
               n1i = 40, n2i = 30, vtype = "LS2")
  expect_equal(d$g[1],      as.numeric(o1$yi), tolerance = TOL)  # pooled
  expect_equal(d$g_se[1]^2, as.numeric(o1$vi), tolerance = TOL)
  expect_equal(d$g[2],      as.numeric(o2$yi), tolerance = TOL)  # control (Glass)
  expect_equal(d$g_se[2]^2, as.numeric(o2$vi), tolerance = TOL)

  # NA cell defaults to "pooled"; an unknown value errors on the argument name
  dn <- es_from_means_sd(mean_exp = c(2.3, 5.1), mean_sd_exp = c(1.2, 2.0),
                         mean_nexp = c(1.9, 4.4), mean_sd_nexp = c(0.9, 1.8),
                         n_exp = c(20, 40), n_nexp = c(22, 30), smd_denom = c(NA, "control"))
  expect_equal(dn$g[1], as.numeric(o1$yi), tolerance = TOL)
  expect_error(
    es_from_means_sd(mean_exp = 1, mean_sd_exp = 1, mean_nexp = 0, mean_sd_nexp = 1,
                     n_exp = 10, n_nexp = 10, smd_denom = "not_a_denom"),
    "smd_denom"
  )

  # and through the pipeline: a mixed smd_denom column drives each row independently
  df <- data.frame(n_exp = c(20, 40), n_nexp = c(22, 30),
                   mean_exp = c(2.3, 5.1), mean_sd_exp = c(1.2, 2.0),
                   mean_nexp = c(1.9, 4.4), mean_sd_nexp = c(0.9, 1.8),
                   smd_denom = c("control", "pooled"))
  g_mixed <- convert_df(df, measure = "g", verbose = FALSE)$es_means_sd_raw$g
  expect_equal(g_mixed[1], as.numeric(escalc(measure = "SMD1", m1i = 2.3, m2i = 1.9,
               sd1i = 1.2, sd2i = 0.9, n1i = 20, n2i = 22, vtype = "LS2")$yi), tolerance = TOL)
  expect_equal(g_mixed[2], as.numeric(o2 <- escalc(measure = "SMD", m1i = 5.1, m2i = 4.4,
               sd1i = 2.0, sd2i = 1.8, n1i = 40, n2i = 30, vtype = "LS2")$yi), tolerance = TOL)
})

test_that("mean-difference CI uses the Welch-Satterthwaite df (matches t.test) and SE matches metafor MD", {
  mk <- function(n, m, s) { x <- qnorm(ppoints(n)); (x - mean(x)) / sd(x) * s + m }
  for (s in SCEN) {
    p <- es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1,
                          mean_nexp = s$m2, mean_sd_nexp = s$s2,
                          n_exp = s$n1, n_nexp = s$n2)
    o  <- escalc(measure = "MD", m1i = s$m1, m2i = s$m2, sd1i = s$s1, sd2i = s$s2,
                 n1i = s$n1, n2i = s$n2)
    tt <- t.test(mk(s$n1, s$m1, s$s1), mk(s$n2, s$m2, s$s2), var.equal = FALSE)
    expect_equal(p$md_se^2, as.numeric(o$vi), tolerance = TOL)             # metafor MD variance
    expect_equal(c(p$md_ci_lo, p$md_ci_up), as.numeric(tt$conf.int),
                 tolerance = 1e-7)                                          # Welch CI (df + bounds)
  }
  # degenerate: both arm SDs zero -> CI collapses to the point estimate (no NaN from 0/0 df)
  z <- es_from_means_sd(mean_exp = 5, mean_sd_exp = 0, mean_nexp = 3, mean_sd_nexp = 0,
                        n_exp = 10, n_nexp = 10)
  expect_equal(z$md_ci_lo, z$md)
  expect_equal(z$md_ci_up, z$md)
})

test_that("Glass reverse_means flips the sign but keeps the control-SD variance", {
  s <- SCEN[[1]]
  g  <- es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1, mean_nexp = s$m2,
                         mean_sd_nexp = s$s2, n_exp = s$n1, n_nexp = s$n2, smd_denom = "control")
  gr <- es_from_means_sd(mean_exp = s$m1, mean_sd_exp = s$s1, mean_nexp = s$m2,
                         mean_sd_nexp = s$s2, n_exp = s$n1, n_nexp = s$n2, smd_denom = "control",
                         reverse_means = TRUE)
  expect_equal(gr$g, -g$g, tolerance = TOL)
  expect_equal(gr$g_se, g$g_se, tolerance = TOL)
})

# =============================================================================
# Backward compatibility - default pipeline output must not move
# =============================================================================
test_that("convert_df default (no smd args) still reproduces metafor SMD LS2", {
  df <- data.frame(study_id = c("a", "b"),
                   n_exp = c(20, 40), n_nexp = c(22, 30),
                   mean_exp = c(2.3, 5.1), mean_sd_exp = c(1.2, 2.0),
                   mean_nexp = c(1.9, 4.4), mean_sd_nexp = c(0.9, 1.8))
  # summary() rounds for display, so read the unrounded g/g_se from the object.
  obj <- convert_df(df, measure = "g", verbose = FALSE)
  ms  <- obj$es_means_sd_raw
  for (r in 1:2) {
    o <- escalc(measure = "SMD", m1i = df$mean_exp[r], m2i = df$mean_nexp[r],
                sd1i = df$mean_sd_exp[r], sd2i = df$mean_sd_nexp[r],
                n1i = df$n_exp[r], n2i = df$n_nexp[r], vtype = "LS2")
    expect_equal(ms$g[r],    as.numeric(o$yi), tolerance = TOL)
    expect_equal(ms$g_se[r], sqrt(as.numeric(o$vi)), tolerance = TOL)
  }
})

# =============================================================================
# P4 - raw-MD family (es_from_md_sd/md_se/md_ci/md_pval) honours smd_var
# =============================================================================
test_that("P4: es_from_md_sd honours smd_var (oracle = escalc SMD LS2/LS on the same d)", {
  for (s in SCEN) {
    md <- 0.42; md_sd <- 1.05
    d  <- md / md_sd
    # escalc inputs reconstructing the same uncorrected d: m1 = d, m2 = 0, sd = 1
    o_ls2 <- escalc(measure = "SMD", m1i = d, m2i = 0, sd1i = 1, sd2i = 1,
                    n1i = s$n1, n2i = s$n2, vtype = "LS2")
    o_ls  <- escalc(measure = "SMD", m1i = d, m2i = 0, sd1i = 1, sd2i = 1,
                    n1i = s$n1, n2i = s$n2, vtype = "LS")
    p_def <- es_from_md_sd(md = md, md_sd = md_sd, n_exp = s$n1, n_nexp = s$n2)
    p_ls  <- es_from_md_sd(md = md, md_sd = md_sd, n_exp = s$n1, n_nexp = s$n2,
                           smd_var = "hedges_olkin")
    expect_equal(p_def$g,      as.numeric(o_ls2$yi), tolerance = TOL)
    expect_equal(p_def$g_se^2, as.numeric(o_ls2$vi), tolerance = TOL)
    expect_equal(p_ls$g_se^2,  as.numeric(o_ls$vi),  tolerance = TOL)
    # regression guard: default must equal the old hard-coded d_se formula
    d_se_old <- sqrt((s$n1 + s$n2) / (s$n1 * s$n2) + d^2 / (2 * (s$n1 + s$n2)))
    expect_equal(p_def$d_se, d_se_old, tolerance = TOL)
  }
})

test_that("P4: md_se / md_ci / md_pval delegate smd_var consistently with md_sd", {
  s <- SCEN[[1]]
  md <- 0.42; md_sd <- 1.05
  md_se_val   <- md_sd * sqrt(1 / s$n1 + 1 / s$n2)
  md_ci_half  <- qt(.975, s$n1 + s$n2 - 2) * md_se_val
  t_val       <- md / md_se_val
  md_pval_val <- 2 * pt(abs(t_val), s$n1 + s$n2 - 2, lower.tail = FALSE)
  for (sv in c("borenstein", "hedges_olkin")) {
    ref <- es_from_md_sd(md = md, md_sd = md_sd, n_exp = s$n1, n_nexp = s$n2, smd_var = sv)
    v1  <- es_from_md_se(md = md, md_se = md_se_val, n_exp = s$n1, n_nexp = s$n2, smd_var = sv)
    v2  <- es_from_md_ci(md = md, md_ci_lo = md - md_ci_half, md_ci_up = md + md_ci_half,
                         n_exp = s$n1, n_nexp = s$n2, smd_var = sv)
    v3  <- es_from_md_pval(md = md, md_pval = md_pval_val, n_exp = s$n1, n_nexp = s$n2, smd_var = sv)
    for (v in list(v1, v2, v3)) {
      expect_equal(v$g,    ref$g,    tolerance = 1e-8)
      expect_equal(v$g_se, ref$g_se, tolerance = 1e-8)
    }
  }
})

test_that("P4: convert_df threads smd_var to the md_sd route", {
  df <- data.frame(study_id = "a", md = 0.42, md_sd = 1.05, n_exp = 8, n_nexp = 12)
  o_def <- convert_df(df, measure = "g", verbose = FALSE)
  o_ho  <- convert_df(df, measure = "g", verbose = FALSE, smd_var = "hedges_olkin")
  se_def <- o_def$es_md_sd$g_se[1]
  se_ho  <- o_ho$es_md_sd$g_se[1]
  expect_false(isTRUE(all.equal(se_def, se_ho)))  # must now change (P4 pre-fix: identical)
  d <- 0.42 / 1.05
  o_ls <- escalc(measure = "SMD", m1i = d, m2i = 0, sd1i = 1, sd2i = 1,
                 n1i = 8, n2i = 12, vtype = "LS")
  expect_equal(se_ho^2, as.numeric(o_ls$vi), tolerance = TOL)
})
