## Guards against impossible dispersion / interval inputs (R/internal_guards.R).
##
## The point of this file is the SWEEP at the bottom: it walks every exported
## route that takes a standard deviation, a standard error or a confidence
## interval, corrupts one input at a time, and asserts that none of them can be
## made to emit a sign-flipped effect size, a negative standard error or a
## transposed output interval. The defect was originally fixed on six routes and
## missed on eleven more with exactly the same shape, so an enumerated list of
## known-bad cases is not enough -- a new route must be caught automatically.

test_that(".positive_or_na() and .ci_width() behave", {
  expect_equal(.positive_or_na(c(2, -2, 0, NA)), c(2, NA, NA, NA))
  expect_identical(.positive_or_na(NULL), NULL)
  expect_identical(.positive_or_na(numeric(0)), numeric(0))
  # a transposed interval describes the same interval, so the same width
  expect_equal(.ci_width(3, 5), .ci_width(5, 3))
  expect_equal(.ci_width(-1, 1), 2)
  expect_equal(.ci_lower(5, 3), 3)
  expect_equal(.ci_upper(5, 3), 5)
  expect_true(is.na(.ci_lower(NA, 3)))
})

test_that("a transposed CI gives the same result as the correctly ordered one", {
  a <- es_from_md_ci(md = 4, md_ci_lo = 3, md_ci_up = 5, n_exp = 50, n_nexp = 50)
  b <- es_from_md_ci(md = 4, md_ci_lo = 5, md_ci_up = 3, n_exp = 50, n_nexp = 50)
  expect_equal(a$d, b$d); expect_equal(a$d_se, b$d_se)

  a <- es_from_linreg_b_ci(linreg_b = 2, linreg_b_ci_lo = 1, linreg_b_ci_up = 3,
                           n_sample = 100, n_covariates = 1, sd_iv = 1.5)
  b <- es_from_linreg_b_ci(linreg_b = 2, linreg_b_ci_lo = 3, linreg_b_ci_up = 1,
                           n_sample = 100, n_covariates = 1, sd_iv = 1.5)
  expect_equal(a$d, b$d); expect_equal(a$d_se, b$d_se)
  expect_gt(a$d, 0)   # and it points the same way as the coefficient

  a <- suppressWarnings(es_from_or_ci(or = 2, or_ci_lo = 1.2, or_ci_up = 3.3, n_exp = 40, n_nexp = 60))
  b <- suppressWarnings(es_from_or_ci(or = 2, or_ci_lo = 3.3, or_ci_up = 1.2, n_exp = 40, n_nexp = 60))
  expect_equal(a$logor_se, b$logor_se)
  expect_gt(b$logor_se, 0)
})

test_that("a non-positive dispersion blanks the row instead of flipping the sign", {
  # the five routes that DIVIDE by the supplied dispersion
  expect_true(is.na(es_from_means_sd_pooled(mean_exp = 10, mean_nexp = 6,
    mean_sd_pooled = -2, n_exp = 50, n_nexp = 50)$d))
  expect_true(is.na(es_from_means_sd_pooled(mean_exp = 10, mean_nexp = 6,
    mean_sd_pooled = 0, n_exp = 50, n_nexp = 50)$d))
  expect_true(is.na(es_from_ancova_means_sd_pooled_adj(ancova_mean_exp = 10,
    ancova_mean_nexp = 6, ancova_mean_sd_pooled = -2, cov_outcome_r = .3,
    n_cov_ancova = 1, n_exp = 50, n_nexp = 50)$d))
  expect_true(is.na(es_from_ancova_means_sd_pooled_crude(ancova_mean_exp = 10,
    ancova_mean_nexp = 6, mean_sd_pooled = -2, cov_outcome_r = .3,
    n_cov_ancova = 1, n_exp = 50, n_nexp = 50)$d))
  expect_true(is.na(es_from_linreg_b_se(linreg_b = 4, linreg_b_se = -0.4,
    n_sample = 100, n_covariates = 2, sd_iv = 1)$d))
  expect_true(is.na(es_from_md_sd(md = 4, md_sd = -2, n_exp = 50, n_nexp = 50)$d))

  # ... and the three that only PROPAGATE it
  expect_true(is.na(suppressWarnings(es_from_or_se(or = 2, logor_se = -0.4,
    n_exp = 50, n_nexp = 50))$logor_se))
  expect_true(is.na(suppressWarnings(es_from_rr_se(rr = 1.6, logrr_se = -0.2,
    n_exp = 50, n_nexp = 50, baseline_risk = .25))$logrr_se))
  expect_true(is.na(suppressWarnings(es_from_rd_se(rd = .2, rd_se = -0.1,
    n_exp = 50, n_nexp = 50, baseline_risk = .3))$rd_se))
})

test_that("valid inputs are untouched by the guards", {
  a <- es_from_means_sd_pooled(mean_exp = 10, mean_nexp = 6, mean_sd_pooled = 2,
                               n_exp = 50, n_nexp = 50)
  expect_equal(a$d, 2)
  # a genuinely negative mean difference still gives a negative effect size
  b <- es_from_means_sd_pooled(mean_exp = 6, mean_nexp = 10, mean_sd_pooled = 2,
                               n_exp = 50, n_nexp = 50)
  expect_equal(b$d, -2)
  expect_gt(b$d_se, 0)
})

test_that("convert_df(correct_inputs = FALSE) no longer emits a sign-flipped ES", {
  d <- data.frame(study_id = c("ok", "bad"), mean_exp = c(10, 10),
                  mean_nexp = c(6, 6), mean_sd_pooled = c(2, -2),
                  n_exp = c(50, 50), n_nexp = c(50, 50))
  s <- suppressWarnings(suppressMessages(
    summary(convert_df(d, measure = "g", correct_inputs = FALSE, verbose = FALSE),
            flags = TRUE)))
  expect_gt(s$es_crude[1], 0)
  expect_true(is.na(s$es_crude[2]))
  # the Tier-1 flag still reports it: the guard contains, the flag informs
  expect_match(s$flags_crude[2], "mean_sd_pooled")
})

## ---------------------------------------------------------------------------
## The sweep.
## ---------------------------------------------------------------------------
test_that("no exported route can be made to emit a flipped ES, a negative SE or a transposed CI", {
  N1 <- 40; N2 <- 60
  B <- list(
    es_from_md_sd = list(md=4, md_sd=2, n_exp=N1, n_nexp=N2),
    es_from_md_se = list(md=4, md_se=0.4, n_exp=N1, n_nexp=N2),
    es_from_md_ci = list(md=4, md_ci_lo=3, md_ci_up=5, n_exp=N1, n_nexp=N2),
    es_from_ancova_md_sd = list(ancova_md=4, ancova_md_sd=2, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_ancova_md_se = list(ancova_md=4, ancova_md_se=0.4, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_ancova_md_ci = list(ancova_md=4, ancova_md_ci_lo=3, ancova_md_ci_up=5, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_means_sd = list(mean_exp=10, mean_nexp=6, mean_sd_exp=2, mean_sd_nexp=2.2, n_exp=N1, n_nexp=N2),
    es_from_means_sd_pooled = list(mean_exp=10, mean_nexp=6, mean_sd_pooled=2, n_exp=N1, n_nexp=N2),
    es_from_means_se = list(mean_exp=10, mean_nexp=6, mean_se_exp=.4, mean_se_nexp=.3, n_exp=N1, n_nexp=N2),
    es_from_means_ci = list(mean_exp=10, mean_ci_lo_exp=9, mean_ci_up_exp=11, mean_nexp=6, mean_ci_lo_nexp=5, mean_ci_up_nexp=7, n_exp=N1, n_nexp=N2),
    es_from_ancova_means_sd = list(ancova_mean_exp=10, ancova_mean_nexp=6, ancova_mean_sd_exp=2, ancova_mean_sd_nexp=2.2, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_ancova_means_sd_pooled_adj = list(ancova_mean_exp=10, ancova_mean_nexp=6, ancova_mean_sd_pooled=2, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_ancova_means_sd_pooled_crude = list(ancova_mean_exp=10, ancova_mean_nexp=6, mean_sd_pooled=2, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_ancova_means_se = list(ancova_mean_exp=10, ancova_mean_nexp=6, ancova_mean_se_exp=.4, ancova_mean_se_nexp=.3, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_ancova_means_ci = list(ancova_mean_exp=10, ancova_mean_ci_lo_exp=9, ancova_mean_ci_up_exp=11, ancova_mean_nexp=6, ancova_mean_ci_lo_nexp=5, ancova_mean_ci_up_nexp=7, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_linreg_b_se = list(linreg_b=2, linreg_b_se=.5, n_sample=100, n_covariates=1, sd_iv=1.5, n_exp=N1, n_nexp=N2),
    es_from_linreg_b_ci = list(linreg_b=2, linreg_b_ci_lo=1, linreg_b_ci_up=3, n_sample=100, n_covariates=1, sd_iv=1.5, n_exp=N1, n_nexp=N2),
    es_from_or_se = list(or=2, logor_se=.3, n_exp=N1, n_nexp=N2),
    es_from_or_ci = list(or=2, or_ci_lo=1.2, or_ci_up=3.3, n_exp=N1, n_nexp=N2),
    es_from_rr_se = list(rr=1.6, logrr_se=.2, n_exp=N1, n_nexp=N2, baseline_risk=.25),
    es_from_rr_ci = list(rr=1.6, rr_ci_lo=1.1, rr_ci_up=2.3, n_exp=N1, n_nexp=N2, baseline_risk=.25),
    es_from_rd_se = list(rd=.2, rd_se=.05, n_exp=N1, n_nexp=N2, baseline_risk=.3),
    es_from_rd_ci = list(rd=.2, rd_ci_lo=.1, rd_ci_up=.3, n_exp=N1, n_nexp=N2, baseline_risk=.3),
    es_from_mean_change_sd = list(mean_change_exp=2, mean_change_sd_exp=3, mean_change_nexp=.5, mean_change_sd_nexp=3, n_exp=N1, n_nexp=N2, r_pre_post_exp=.6, r_pre_post_nexp=.6),
    es_from_mean_change_se = list(mean_change_exp=2, mean_change_se_exp=.4, mean_change_nexp=.5, mean_change_se_nexp=.4, n_exp=N1, n_nexp=N2, r_pre_post_exp=.6, r_pre_post_nexp=.6),
    es_from_mean_change_ci = list(mean_change_exp=2, mean_change_ci_lo_exp=1, mean_change_ci_up_exp=3, mean_change_nexp=.5, mean_change_ci_lo_nexp=-.5, mean_change_ci_up_nexp=1.5, n_exp=N1, n_nexp=N2, r_pre_post_exp=.6, r_pre_post_nexp=.6),
    es_from_mean_change_sd_single_group = list(mean_change_exp=2, mean_change_sd_exp=3, n_exp=N1, r_pre_post_exp=.6),
    es_from_mean_change_se_single_group = list(mean_change_exp=2, mean_change_se_exp=.4, n_exp=N1, r_pre_post_exp=.6),
    es_from_mean_change_ci_single_group = list(mean_change_exp=2, mean_change_ci_lo_exp=1, mean_change_ci_up_exp=3, n_exp=N1, r_pre_post_exp=.6),
    es_from_means_sd_pre_post = list(mean_pre_exp=8, mean_exp=10, mean_pre_sd_exp=3, mean_sd_exp=3, mean_pre_nexp=8, mean_nexp=8.5, mean_pre_sd_nexp=3, mean_sd_nexp=3, n_exp=N1, n_nexp=N2, r_pre_post_exp=.6, r_pre_post_nexp=.6),
    es_from_means_se_pre_post = list(mean_pre_exp=8, mean_exp=10, mean_pre_se_exp=.4, mean_se_exp=.4, mean_pre_nexp=8, mean_nexp=8.5, mean_pre_se_nexp=.4, mean_se_nexp=.4, n_exp=N1, n_nexp=N2, r_pre_post_exp=.6, r_pre_post_nexp=.6),
    es_from_means_ci_pre_post = list(mean_pre_exp=8, mean_exp=10, mean_pre_ci_lo_exp=7, mean_pre_ci_up_exp=9, mean_ci_lo_exp=9, mean_ci_up_exp=11, mean_pre_nexp=8, mean_nexp=8.5, mean_pre_ci_lo_nexp=7, mean_pre_ci_up_nexp=9, mean_ci_lo_nexp=7.5, mean_ci_up_nexp=9.5, n_exp=N1, n_nexp=N2, r_pre_post_exp=.6, r_pre_post_nexp=.6),
    es_from_means_sd_pre_post_single_group = list(mean_pre_exp=8, mean_exp=10, mean_pre_sd_exp=3, mean_sd_exp=3, n_exp=N1, r_pre_post_exp=.6),
    es_from_means_se_pre_post_single_group = list(mean_pre_exp=8, mean_exp=10, mean_pre_se_exp=.4, mean_se_exp=.4, n_exp=N1, r_pre_post_exp=.6),
    es_from_means_ci_pre_post_single_group = list(mean_pre_exp=8, mean_exp=10, mean_pre_ci_lo_exp=7, mean_pre_ci_up_exp=9, mean_ci_lo_exp=9, mean_ci_up_exp=11, n_exp=N1, r_pre_post_exp=.6),
    es_from_plot_means = list(plot_mean_exp=10, plot_mean_sd_lo_exp=8, plot_mean_sd_up_exp=12, plot_mean_nexp=6, plot_mean_sd_lo_nexp=4, plot_mean_sd_up_nexp=8, n_exp=N1, n_nexp=N2),
    es_from_plot_ancova_means = list(plot_ancova_mean_exp=10, plot_ancova_mean_sd_lo_exp=8, plot_ancova_mean_sd_up_exp=12, plot_ancova_mean_nexp=6, plot_ancova_mean_sd_lo_nexp=4, plot_ancova_mean_sd_up_nexp=8, cov_outcome_r=.3, n_cov_ancova=1, n_exp=N1, n_nexp=N2),
    es_from_user_crude    = list(user_es_measure_crude="d", user_es_crude=.6, user_ci_lo_crude=.2, user_ci_up_crude=1.0, n_exp=N1, n_nexp=N2),
    es_from_user_crude_se = list(user_es_measure_crude="d", user_es_crude=.6, user_se_crude=.2, n_exp=N1, n_nexp=N2),
    es_from_user_crude_or = list(user_es_measure_crude="or", user_es_crude=2, user_ci_lo_crude=1.2, user_ci_up_crude=3.3, n_exp=N1, n_nexp=N2),
    es_from_user_adj      = list(user_es_measure_adj="d", user_es_adj=.6, user_ci_lo_adj=.2, user_ci_up_adj=1.0, n_exp=N1, n_nexp=N2),
    es_from_user_adj_se   = list(user_es_measure_adj="d", user_es_adj=.6, user_se_adj=.2, n_exp=N1, n_nexp=N2),
    es_from_user_adj_or   = list(user_es_measure_adj="or", user_es_adj=2, user_ci_lo_adj=1.2, user_ci_up_adj=3.3, n_exp=N1, n_nexp=N2)
  )
  ## es_from_user_* appears several times under disambiguating list names; strip
  ## the suffix for THOSE only -- es_from_rd_se and friends are real route names.
  call1 <- function(fn, args) {
    if (grepl("^es_from_user_(crude|adj)_(se|or)$", fn)) fn <- sub("_(se|or)$", "", fn)
    suppressWarnings(suppressMessages(try(do.call(fn, args), silent = TRUE)))
  }
  prim <- function(e) {
    for (k in c("d", "logor", "logrr", "rd", "md"))
      if (!is.null(e[[k]]) && !is.na(e[[k]][1])) return(e[[k]][1])
    NA_real_
  }
  neg_se <- function(e) {
    k <- grep("_se$", names(e), value = TRUE)
    v <- suppressWarnings(as.numeric(unlist(lapply(k, function(x) e[[x]][1]))))
    any(v[!is.na(v)] < 0)
  }
  inv_ci <- function(e) {
    any(vapply(grep("_ci_lo$", names(e), value = TRUE), function(lo) {
      up <- sub("_ci_lo$", "_ci_up", lo)
      if (!up %in% names(e)) return(FALSE)
      a <- suppressWarnings(as.numeric(e[[lo]][1])); b <- suppressWarnings(as.numeric(e[[up]][1]))
      !is.na(a) && !is.na(b) && a > b
    }, TRUE))
  }
  check <- function(fn, label, e, p0) {
    p1 <- prim(e)
    expect_false(!is.na(p0) && !is.na(p1) && sign(p1) != sign(p0),
                 label = paste0(fn, " [", label, "] flips the effect size"))
    expect_false(neg_se(e), label = paste0(fn, " [", label, "] returns a negative SE"))
    expect_false(inv_ci(e), label = paste0(fn, " [", label, "] returns a transposed CI"))
  }

  n_cases <- 0L
  for (fn in names(B)) {
    base <- B[[fn]]
    ok <- call1(fn, base)
    expect_false(inherits(ok, "try-error"), label = paste0(fn, " baseline call"))
    if (inherits(ok, "try-error")) next
    p0 <- prim(ok); nm <- names(base)

    disp <- setdiff(nm[grepl("_sd$|_se$|_sd_exp$|_sd_nexp$|_se_exp$|_se_nexp$|_sd_pooled$|^user_se_", nm)],
                    "pool_sd")
    for (a in disp) {
      b2 <- base; b2[[a]] <- -abs(b2[[a]])
      e <- call1(fn, b2); if (inherits(e, "try-error")) next
      check(fn, paste0(a, " < 0"), e, p0); n_cases <- n_cases + 1L
    }
    for (lo in nm[grepl("_ci_lo", nm)]) {
      up <- sub("_ci_lo", "_ci_up", lo); if (!up %in% nm) next
      b2 <- base; b2[[lo]] <- base[[up]]; b2[[up]] <- base[[lo]]
      e <- call1(fn, b2); if (inherits(e, "try-error")) next
      check(fn, paste0(lo, "/", up, " transposed"), e, p0); n_cases <- n_cases + 1L
    }
  }
  expect_gte(n_cases, 55L)   # guard against the sweep silently going empty
})
