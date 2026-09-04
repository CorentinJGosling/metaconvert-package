# Third-round audit fixes (2.1.0). Every assertion here is pinned to an EXTERNAL
# reference -- metafor 5.0.1, a closed form with a citation, or an identity the
# package's own forward map must satisfy -- never to a value the code once produced.

.quiet <- function(expr) suppressMessages(suppressWarnings(expr))
.smry <- function(df, measure, ...) {
  as.data.frame(.quiet(summary(.quiet(convert_df(df, measure = measure, ...)),
                               flags = TRUE)))
}

# --- 1. user_es_* correlation with a CI: precision read on the Fisher-z scale ------

test_that("a user r with a Fisher CI recovers z_se = 1/sqrt(n - 3) and the atanh bounds", {
  for (n in c(30, 100)) for (r in c(0.3, 0.7, 0.9)) {
    lo <- tanh(atanh(r) - qnorm(.975) / sqrt(n - 3))
    up <- tanh(atanh(r) + qnorm(.975) / sqrt(n - 3))
    df <- data.frame(study_id = "A", user_es_original_measure_crude = "r",
                     user_es_crude = r, user_ci_lo_crude = lo, user_ci_up_crude = up,
                     n_sample = n)
    s <- .smry(df, "z")
    # metafor::escalc(measure = "ZCOR"): vi = 1/(n - 3)
    expect_equal(s$se_crude, 1 / sqrt(n - 3), tolerance = 1e-6)
    expect_equal(s$es_ci_lo_crude, atanh(lo), tolerance = 1e-6)
    expect_equal(s$es_ci_up_crude, atanh(up), tolerance = 1e-6)
    # the r-scale SE is the delta-method image of that z SE, so the pair is coherent
    sr <- .smry(df, "r")
    expect_equal(sr$se_crude, (1 - r^2) / sqrt(n - 3), tolerance = 1e-6)
    # the user's own r bounds are handed back verbatim
    expect_equal(c(sr$es_ci_lo_crude, sr$es_ci_up_crude), c(lo, up), tolerance = 1e-10)
  }
})

test_that("the adjusted user-r path applies the same rule, and a supplied se is never overwritten", {
  df <- data.frame(study_id = "A",
                   user_es_original_measure_adj = "r", user_es_adj = NA,
                   user_ci_lo_adj = 0.2, user_ci_up_adj = 0.8)
  s <- .smry(df, "z")
  expect_equal(s$se_adjusted, (atanh(0.8) - atanh(0.2)) / (2 * qnorm(.975)), tolerance = 1e-6)
  # CI-only point estimate: centre on the z scale, not the arithmetic r midpoint
  expect_equal(s$es_adjusted, (atanh(0.8) + atanh(0.2)) / 2, tolerance = 1e-6)
  sr <- .smry(df, "r")
  expect_equal(sr$es_adjusted, tanh((atanh(0.8) + atanh(0.2)) / 2), tolerance = 1e-6)

  df2 <- data.frame(study_id = "A", user_es_original_measure_crude = "r",
                    user_es_crude = 0.5, user_se_crude = 0.1,
                    user_ci_lo_crude = 0.2, user_ci_up_crude = 0.7)
  expect_equal(.smry(df2, "r")$se_crude, 0.1)
})

# --- 2. incidence rate ratio with a zero event count -------------------------------

test_that("a zero event count takes the +0.5 correction metafor applies, and stays finite", {
  x1 <- c(0, 12, 5); x2 <- c(14, 0, 6); t1 <- c(100, 90, 120); t2 <- c(100, 95, 110)
  r <- es_from_cases_time(n_cases_exp = x1, n_cases_nexp = x2, time_exp = t1, time_nexp = t2)
  m <- metafor::escalc(measure = "IRR", x1i = x1, x2i = x2, t1i = t1, t2i = t2)
  expect_equal(r$logirr, as.numeric(m$yi), tolerance = 1e-12)
  expect_equal(r$logirr_se, sqrt(as.numeric(m$vi)), tolerance = 1e-12)
  expect_true(all(is.finite(r$logirr)) && all(is.finite(r$logirr_se)))
  # the rate difference is on the RAW counts (finite at a zero count without the +0.5)
  expect_equal(r$rd, x2 / t2 - x1 / t1, tolerance = 1e-12)
  expect_equal(r$rd_se, sqrt(x1 / t1^2 + x2 / t2^2), tolerance = 1e-12)

  s <- .smry(data.frame(study_id = c("A", "B"), n_cases_exp = c(0, 5), n_cases_nexp = c(14, 6),
                        time_exp = c(100, 120), time_nexp = c(100, 110)), "irr")
  expect_equal(s$es_crude[1], exp(as.numeric(m$yi[1])), tolerance = 1e-10)
  expect_false(grepl("non-finite", s$flags_crude[1]))
})

# --- 3. r -> d inverts the biserial map at the supplied arm split ------------------

test_that("es_from_pearson_r() with arm sizes inverts .smd_to_cor()'s biserial exactly", {
  for (n1 in c(50, 30, 20, 10)) {
    n2 <- 100 - n1
    fwd <- .quiet(es_from_means_sd(mean_exp = 12, mean_sd_exp = 4, mean_nexp = 10,
                                   mean_sd_nexp = 4, n_exp = n1, n_nexp = n2))
    back <- .quiet(es_from_pearson_r(pearson_r = fwd$r, n_sample = 100,
                                     n_exp = n1, n_nexp = n2))
    # transf.rtod's output lands in the g slot (documented convention), and with
    # n1i/n2i it returns the uncorrected d of the forward map to machine precision
    expect_equal(back$g, fwd$d, tolerance = 1e-10)
    expect_equal(back$g, metafor::transf.rtod(fwd$r, n1i = n1, n2i = n2), tolerance = 1e-12)
  }
})

test_that("without arm sizes the balanced-arms inversion is unchanged", {
  r <- 0.2829377771
  a <- .quiet(es_from_pearson_r(pearson_r = r, n_sample = 100))
  expect_equal(a$g, metafor::transf.rtod(r), tolerance = 1e-12)
  # only n_exp -> the other arm is back-filled, so the conversion stays balanced
  # (compared on d: g additionally carries .es_from_d()'s J on the back-filled df)
  b <- .quiet(es_from_pearson_r(pearson_r = r, n_sample = 100, n_exp = 20))
  expect_equal(b$d, a$d, tolerance = 1e-12)
})

test_that("the user-input r branch also honours supplied arm sizes", {
  r <- 0.2829377771
  s <- .smry(data.frame(study_id = "A", user_es_original_measure_crude = "r",
                        user_es_crude = r, user_se_crude = 0.09, n_exp = 20, n_nexp = 80), "d")
  expect_equal(s$es_crude, metafor::transf.rtod(r, 20, 80) / metaConvert:::.d_j(98),
               tolerance = 1e-8)
})

# --- 4. degenerate p-values on every p-value route ---------------------------------

test_that("p <= 0 yields NA on every route; p >= 1 yields NA only where se = |es / z|", {
  se_ratio <- list(
    rr   = list(data.frame(study_id = "S", rr = 1.5, rr_pval = NA, n_exp = 50, n_nexp = 50), "rr", "rr_pval"),
    rd   = list(data.frame(study_id = "S", rd = 0.1, rd_pval = NA, n_exp = 50, n_nexp = 50), "rd", "rd_pval"),
    md   = list(data.frame(study_id = "S", md = 2, md_pval = NA, n_exp = 50, n_nexp = 50), "md", "md_pval"),
    or   = list(data.frame(study_id = "S", or = 1.5, or_pval = NA, n_exp = 50, n_nexp = 50), "or", "or_pval"),
    mc   = list(data.frame(study_id = "S", mean_change_exp = 3, mean_change_nexp = 1,
                           mean_change_pval_exp = NA, mean_change_pval_nexp = 0.3,
                           n_exp = 50, n_nexp = 50), "d", "mean_change_pval_exp")
  )
  for (nm in names(se_ratio)) {
    df <- se_ratio[[nm]][[1]]; measure <- se_ratio[[nm]][[2]]; col <- se_ratio[[nm]][[3]]
    for (p in c(0, 1)) {
      df[[col]] <- p
      s <- .smry(df, measure)
      # es_from_or() (imputed SE from margins) legitimately takes the OR row; the
      # p-value route itself must not contribute an Inf/0 se
      if (nm == "or") {
        expect_false(identical(s$info_used_crude, "or_pval"), info = paste(nm, p))
      } else {
        expect_true(is.na(s$es_crude) && is.na(s$se_crude), info = paste(nm, p))
      }
    }
    df[[col]] <- 0.5
    expect_true(is.finite(.smry(df, measure)$se_crude), info = nm)
  }

  stat <- list(
    t   = list(data.frame(study_id = "S", student_t_pval = NA, n_exp = 50, n_nexp = 50), "d", "student_t_pval"),
    f   = list(data.frame(study_id = "S", anova_f_pval = NA, n_exp = 50, n_nexp = 50), "d", "anova_f_pval"),
    pb  = list(data.frame(study_id = "S", pt_bis_r_pval = NA, n_sample = 100, n_exp = 50, n_nexp = 50), "d", "pt_bis_r_pval"),
    pt  = list(data.frame(study_id = "S", paired_t_pval_exp = NA, paired_t_pval_nexp = 0.3,
                          n_exp = 50, n_nexp = 50), "d", "paired_t_pval_exp"),
    chi = list(data.frame(study_id = "S", chisq_pval = NA, n_sample = 100, n_cases = 40, n_exp = 50), "or", "chisq_pval")
  )
  for (nm in names(stat)) {
    df <- stat[[nm]][[1]]; measure <- stat[[nm]][[2]]; col <- stat[[nm]][[3]]
    df[[col]] <- 0
    s <- .smry(df, measure)
    expect_true(is.na(s$es_crude) && is.na(s$se_crude), info = paste(nm, "p = 0"))
    # p = 1 is the ordinary t = 0 / chi-square = 0 boundary: a finite null estimate
    df[[col]] <- 1
    s <- .smry(df, measure)
    expect_true(is.finite(s$es_crude) && is.finite(s$se_crude), info = paste(nm, "p = 1"))
  }

  # direct routes (the metaumbrella path) behave the same as convert_df()
  expect_true(is.na(.quiet(es_from_or_pval(or = 2, or_pval = 1, n_exp = 40, n_nexp = 60))$logor_se))
  expect_true(is.na(.quiet(es_from_rr_pval(rr = 2, rr_pval = 0, n_exp = 40, n_nexp = 60))$logrr_se))
  expect_true(is.na(.quiet(es_from_student_t_pval(student_t_pval = 0, n_exp = 40, n_nexp = 60))$d))
})

# --- 9. risk difference takes the +0.5 correction on a zero-cell table -----------

test_that("es_from_2x2() rd / rd_se equal metafor's RD (add = 1/2, to = 'only0') on zero-cell tables", {
  set.seed(1)
  a <- sample(0:60, 300, TRUE); b <- sample(0:80, 300, TRUE)
  cc <- sample(0:60, 300, TRUE); d <- sample(0:80, 300, TRUE)
  keep <- (a + b) > 0 & (cc + d) > 0
  a <- a[keep]; b <- b[keep]; cc <- cc[keep]; d <- d[keep]
  r <- .quiet(es_from_2x2(n_cases_exp = a, n_controls_exp = b,
                          n_cases_nexp = cc, n_controls_nexp = d))
  m <- metafor::escalc(measure = "RD", ai = a, bi = b, ci = cc, di = d)
  expect_gt(sum(a == 0 | b == 0 | cc == 0 | d == 0), 10)
  # metaConvert reports control minus exposed
  expect_equal(r$rd, -as.numeric(m$yi), tolerance = 1e-12)
  expect_equal(r$rd_se, sqrt(as.numeric(m$vi)), tolerance = 1e-12)
})

# --- 4. smd_var reaches every pre/post branch ------------------------------------

test_that("hedges_olkin is bit-exact with metafor on all four single-group branches", {
  set.seed(3); K <- 60
  n <- sample(5:150, K, TRUE); m1 <- rnorm(K, 20, 5); m2 <- m1 + rnorm(K, 2, 3)
  s1 <- runif(K, 2, 8); s2 <- runif(K, 2, 8); r <- runif(K, -0.2, 0.95)
  run <- function(meth, sv, sd_post = s2) .quiet(es_from_means_sd_pre_post_single_group(
    mean_pre_exp = m1, mean_exp = m2, mean_pre_sd_exp = s1, mean_sd_exp = sd_post,
    n_exp = n, r_pre_post_exp = r, pre_post_to_smd = meth, smd_var = sv))
  ref <- list(
    bonett     = metafor::escalc(measure = "SMCRH",  m1i = m2, m2i = m1, sd1i = s1, sd2i = s2, ni = n, ri = r),
    morris_dz  = metafor::escalc(measure = "SMCC",   m1i = m2, m2i = m1, sd1i = s2, sd2i = s1, ni = n, ri = r),
    morris_dav = metafor::escalc(measure = "SMCRPH", m1i = m2, m2i = m1, sd1i = s2, sd2i = s1, ni = n, ri = r))
  J <- metaConvert:::.d_j(n - 1)
  for (meth in names(ref)) {
    ho <- run(meth, "hedges_olkin"); bo <- run(meth, "borenstein")
    expect_equal(ho$g_se, sqrt(as.numeric(ref[[meth]]$vi)), tolerance = 1e-12, info = meth)
    expect_equal(ho$g, bo$g, tolerance = 1e-12, info = meth)              # estimate untouched
    # borenstein = J^2 (L + d^2 C), i.e. the LS variance with g^2 -> d^2 and J^2 outside
    C <- if (meth == "morris_dz") 1 / (2 * n) else NA
    if (meth == "morris_dz") {
      expect_equal(bo$g_se, J * sqrt(1 / n + bo$d^2 * C), tolerance = 1e-12)
    }
  }
  # morris_drm under homoscedasticity is metafor SMCR
  ho <- run("morris_drm", "hedges_olkin", sd_post = s1)
  mR <- metafor::escalc(measure = "SMCR", m1i = m2, m2i = m1, sd1i = s1, ni = n, ri = r)
  expect_equal(ho$g_se, sqrt(as.numeric(mR$vi)), tolerance = 1e-12)
})

test_that("morris_drm and morris_dz agree exactly at r = 0.5 on both conventions", {
  for (sv in c("hedges_olkin", "borenstein")) {
    f <- function(m) .quiet(es_from_means_sd_pre_post_single_group(
      mean_pre_exp = 10, mean_exp = 12, mean_pre_sd_exp = 4, mean_sd_exp = 4, n_exp = 5,
      r_pre_post_exp = 0.5, pre_post_to_smd = m, smd_var = sv))
    expect_equal(f("morris_drm")$g, f("morris_dz")$g, tolerance = 1e-12)
    expect_equal(f("morris_drm")$g_se, f("morris_dz")$g_se, tolerance = 1e-12, info = sv)
  }
})

test_that("smd_var is threaded through convert_df() to the paired and mean-change routes", {
  df <- data.frame(study_id = "A", mean_change_exp = 3, mean_change_sd_exp = 4,
                   mean_change_nexp = 1, mean_change_sd_nexp = 4, n_exp = 10, n_nexp = 10)
  ho <- .smry(df, "g", smd_var = "hedges_olkin"); bo <- .smry(df, "g", smd_var = "borenstein")
  expect_equal(ho$es_crude, bo$es_crude, tolerance = 1e-12)
  # per-arm d_rm = (change / sd_change) * sqrt(2(1 - r)), r = 0.8 default, variances
  # added (pool_sd = FALSE): LS is metafor's SMCR form 2(1-r)/n + g^2/(2n) on each arm's
  # own g, LS2 is J^2 (2(1-r)/n + d^2/(2n))
  J <- metaConvert:::.d_j(9)
  d1 <- 3 / 4 * sqrt(2 * 0.2); d2 <- 1 / 4 * sqrt(2 * 0.2)
  expect_equal(ho$es_crude, J * (d1 - d2), tolerance = 1e-10)
  expect_equal(ho$se_crude, sqrt((2 * 0.2 / 10 + (J * d1)^2 / 20) + (2 * 0.2 / 10 + (J * d2)^2 / 20)), tolerance = 1e-10)
  expect_equal(bo$se_crude, J * sqrt((2 * 0.2 / 10 + d1^2 / 20) + (2 * 0.2 / 10 + d2^2 / 20)), tolerance = 1e-10)
  p <- .quiet(es_from_paired_t(paired_t_exp = 3, paired_t_nexp = 1, n_exp = 30, n_nexp = 25,
                               r_pre_post_exp = .5, r_pre_post_nexp = .5, smd_var = "hedges_olkin"))
  expect_equal(p$g_se, sqrt(as.numeric(metafor::escalc(measure = "SMCR", m1i = 3 / sqrt(30), m2i = 0, sd1i = 1, ni = 30, ri = .5)$vi) +
                            as.numeric(metafor::escalc(measure = "SMCR", m1i = 1 / sqrt(25), m2i = 0, sd1i = 1, ni = 25, ri = .5)$vi)),
               tolerance = 1e-10)
})

# --- 6. correlation intervals are back-transformed on every route -------------------

test_that("pearson_r / spearman / user r / lipsey_cooper intervals equal tanh of the z interval", {
  r <- c(-0.226, 0.6, 0.85); n <- c(5, 12, 8)
  e <- .quiet(es_from_pearson_r(pearson_r = r, n_sample = n))
  z <- metafor::escalc(measure = "ZCOR", ri = r, ni = n)
  bt <- summary(z, transf = metafor::transf.ztor)
  expect_equal(e$r_ci_lo, as.numeric(bt$ci.lb), tolerance = 1e-10)
  expect_equal(e$r_ci_up, as.numeric(bt$ci.ub), tolerance = 1e-10)
  expect_true(all(abs(c(e$r_ci_lo, e$r_ci_up)) <= 1))
  sp <- .quiet(es_from_spearman_rho(spearman_r = 0.9, n_sample = 8))
  expect_equal(c(sp$r_ci_lo, sp$r_ci_up), tanh(c(sp$z_ci_lo, sp$z_ci_up)), tolerance = 1e-12)
  lc <- .quiet(es_from_means_sd(mean_exp = 12, mean_sd_exp = 4, mean_nexp = 10, mean_sd_nexp = 4,
                                n_exp = 8, n_nexp = 8, smd_to_cor = "lipsey_cooper"))
  expect_equal(c(lc$r_ci_lo, lc$r_ci_up), tanh(c(lc$z_ci_lo, lc$z_ci_up)), tolerance = 1e-12)
  u <- .smry(data.frame(study_id = "A", user_es_original_measure_crude = "r",
                        user_es_crude = 0.9, user_se_crude = 0.08), "r")
  expect_equal(c(u$es_ci_lo_crude, u$es_ci_up_crude),
               tanh(atanh(0.9) + c(-1, 1) * qnorm(.975) * 0.08 / (1 - 0.81)), tolerance = 1e-10)
  # and A6 does not fire on the package's own small-n correlation rows
  s <- .smry(data.frame(study_id = c("A", "B", "C"), pearson_r = r, n_sample = n), "r")
  expect_false(any(grepl("CI width inconsistent", s$flags_crude)))
})
