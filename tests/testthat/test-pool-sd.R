test_that("pool_sd=FALSE (default) matches existing behavior", {
  res_default <- es_from_means_sd_pre_post(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 6,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 6,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "bonett"
  )

  res_explicit <- es_from_means_sd_pre_post(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 6,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 6,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "bonett",
    pool_sd = FALSE
  )

  expect_identical(res_default, res_explicit)
})

test_that("pool_sd=TRUE produces different results from pool_sd=FALSE", {
  res_unpooled <- es_from_means_sd_pre_post(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 6,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 6,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "bonett",
    pool_sd = FALSE
  )

  res_pooled <- es_from_means_sd_pre_post(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 6,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 6,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "bonett",
    pool_sd = TRUE
  )

  expect_false(identical(res_pooled$g, res_unpooled$g))
  expect_false(identical(res_pooled$g_se, res_unpooled$g_se))
})

test_that("pool_sd=TRUE with equal groups converges to pool_sd=FALSE for d", {
  # When both groups have identical SDs, n, and r, the pooled and unpooled
  # approaches should give the same d (but not necessarily the same variance)
  res_unpooled <- es_from_means_sd_pre_post(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 5,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 5,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "bonett",
    pool_sd = FALSE
  )

  res_pooled <- es_from_means_sd_pre_post(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 5,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 5,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "bonett",
    pool_sd = TRUE
  )

  # d values should be identical (same SD for both groups)
  # g values may differ slightly because J correction uses different df
  # (n-1 per group in unpooled vs N-2 pooled)
  expect_equal(res_pooled$d, res_unpooled$d, tolerance = 1e-6)
})

test_that("pool_sd bonett method: manual formula verification", {
  n1 <- 36; n2 <- 35
  mean_pre1 <- 20; mean_post1 <- 28
  sd_pre1 <- 5; sd_post1 <- 6
  mean_pre2 <- 20; mean_post2 <- 22
  sd_pre2 <- 4; sd_post2 <- 5
  r1 <- 0.6; r2 <- 0.7

  res <- es_from_means_sd_pre_post(
    n_exp = n1, n_nexp = n2,
    mean_pre_exp = mean_pre1, mean_exp = mean_post1,
    mean_pre_sd_exp = sd_pre1, mean_sd_exp = sd_post1,
    mean_pre_nexp = mean_pre2, mean_nexp = mean_post2,
    mean_pre_sd_nexp = sd_pre2, mean_sd_nexp = sd_post2,
    r_pre_post_exp = r1, r_pre_post_nexp = r2,
    pre_post_to_smd = "bonett",
    pool_sd = TRUE
  )

  # Manual calculation
  N <- n1 + n2
  m <- N - 2
  J <- metaConvert:::.d_j(m)
  mean_diff <- (mean_post1 - mean_pre1) - (mean_post2 - mean_pre2)
  sd_pooled <- sqrt(((n1 - 1) * sd_pre1^2 + (n2 - 1) * sd_pre2^2) / m)
  d_manual <- mean_diff / sd_pooled
  g_manual <- d_manual * J
  r_avg <- (n1 * r1 + n2 * r2) / N

  var_g_manual <- 2 * J^2 * (1 - r_avg) * (N / (n1 * n2)) * (m / (m - 2)) *
                  (1 + (n1 * n2 / N) * g_manual^2 / (2 * (1 - r_avg))) - g_manual^2

  expect_equal(res$d, d_manual, tolerance = 1e-8)
  expect_equal(res$g, g_manual, tolerance = 1e-8)
  expect_equal(res$g_se^2, var_g_manual, tolerance = 1e-8)
})

test_that("pool_sd morris_dz method: manual formula verification", {
  n1 <- 40; n2 <- 38
  mean_pre1 <- 15; mean_post1 <- 22
  sd_pre1 <- 4; sd_post1 <- 5
  mean_pre2 <- 15; mean_post2 <- 17
  sd_pre2 <- 3.5; sd_post2 <- 4.5
  r1 <- 0.5; r2 <- 0.6

  res <- es_from_means_sd_pre_post(
    n_exp = n1, n_nexp = n2,
    mean_pre_exp = mean_pre1, mean_exp = mean_post1,
    mean_pre_sd_exp = sd_pre1, mean_sd_exp = sd_post1,
    mean_pre_nexp = mean_pre2, mean_nexp = mean_post2,
    mean_pre_sd_nexp = sd_pre2, mean_sd_nexp = sd_post2,
    r_pre_post_exp = r1, r_pre_post_nexp = r2,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE
  )

  # Manual calculation
  N <- n1 + n2
  m <- N - 2
  J <- metaConvert:::.d_j(m)
  mean_diff <- (mean_post1 - mean_pre1) - (mean_post2 - mean_pre2)

  sd_change1 <- sqrt(sd_pre1^2 + sd_post1^2 - 2 * r1 * sd_pre1 * sd_post1)
  sd_change2 <- sqrt(sd_pre2^2 + sd_post2^2 - 2 * r2 * sd_pre2 * sd_post2)
  sd_pooled <- sqrt(((n1 - 1) * sd_change1^2 + (n2 - 1) * sd_change2^2) / m)

  d_manual <- mean_diff / sd_pooled
  g_manual <- d_manual * J

  # Two-sample SMD variance on the change-score scale (Hedges 1981; metafor
  # escalc measure = "SMD" on change scores). No 2(1-r) factor: the dz point
  # estimate carries no sqrt(2(1-r)) rescaling, so its variance cannot either
  # (the factor belongs to the raw-score-metric morris_drm branch only).
  var_g_manual <- J^2 * (N / (n1 * n2) + g_manual^2 / (2 * m))

  expect_equal(res$d, d_manual, tolerance = 1e-8)
  expect_equal(res$g, g_manual, tolerance = 1e-8)
  expect_equal(res$g_se^2, var_g_manual, tolerance = 1e-8)
})

test_that("pool_sd morris_dz variance matches metafor::escalc(measure = 'SMD') on change scores", {
  skip_if_not_installed("metafor")
  # External anchor for the pooled dz branch: entering the change scores as
  # two independent groups in escalc(measure = "SMD") targets exactly the
  # pooled-change-SD SMD. escalc uses vi = 1/n1 + 1/n2 + g^2/(2N) while the
  # package uses the J^2(N/(n1 n2) + g^2/(2m)) convention; both are equivalent
  # to first order, hence the 2% tolerance on the SE (point estimate exact).
  n1 <- 29; n2 <- 34
  mc1 <- -12.4; sdc1 <- 6.52
  mc2 <- -3.5; sdc2 <- 10.67

  res <- es_from_mean_change_sd(
    n_exp = n1, n_nexp = n2,
    mean_change_exp = mc1, mean_change_sd_exp = sdc1,
    mean_change_nexp = mc2, mean_change_sd_nexp = sdc2,
    r_pre_post_exp = 0.6, r_pre_post_nexp = 0.6,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE
  )

  ef <- metafor::escalc(
    measure = "SMD",
    m1i = mc1, sd1i = sdc1, n1i = n1,
    m2i = mc2, sd2i = sdc2, n2i = n2
  )

  expect_equal(res$g, as.numeric(ef$yi), tolerance = 1e-6)
  expect_equal(res$g_se, sqrt(as.numeric(ef$vi)), tolerance = 2e-2)
})

test_that("pool_sd morris_drm equals morris_dz rescaled by sqrt(2(1-r_avg))", {
  # Metric-transformation identity (Morris & DeShon 2002): the raw-score-metric
  # drm estimator is the change-SD-metric dz estimator times sqrt(2(1-r)), in
  # both the point estimate and the SE.
  n1 <- 29; n2 <- 34; r <- 0.6

  args <- list(
    n_exp = n1, n_nexp = n2,
    mean_change_exp = -12.4, mean_change_sd_exp = 6.52,
    mean_change_nexp = -3.5, mean_change_sd_nexp = 10.67,
    r_pre_post_exp = r, r_pre_post_nexp = r,
    pool_sd = TRUE
  )
  res_dz <- do.call(es_from_mean_change_sd, c(args, pre_post_to_smd = "morris_dz"))
  res_drm <- do.call(es_from_mean_change_sd, c(args, pre_post_to_smd = "morris_drm"))

  k <- sqrt(2 * (1 - r))
  expect_equal(res_drm$g, res_dz$g * k, tolerance = 1e-10)
  expect_equal(res_drm$g_se, res_dz$g_se * k, tolerance = 1e-10)
})

test_that("pool_sd morris_drm method: manual formula verification", {
  n1 <- 40; n2 <- 38
  mean_pre1 <- 15; mean_post1 <- 22
  sd_pre1 <- 4; sd_post1 <- 5
  mean_pre2 <- 15; mean_post2 <- 17
  sd_pre2 <- 3.5; sd_post2 <- 4.5
  r1 <- 0.5; r2 <- 0.6

  res <- es_from_means_sd_pre_post(
    n_exp = n1, n_nexp = n2,
    mean_pre_exp = mean_pre1, mean_exp = mean_post1,
    mean_pre_sd_exp = sd_pre1, mean_sd_exp = sd_post1,
    mean_pre_nexp = mean_pre2, mean_nexp = mean_post2,
    mean_pre_sd_nexp = sd_pre2, mean_sd_nexp = sd_post2,
    r_pre_post_exp = r1, r_pre_post_nexp = r2,
    pre_post_to_smd = "morris_drm",
    pool_sd = TRUE
  )

  # Manual calculation
  N <- n1 + n2
  m <- N - 2
  J <- metaConvert:::.d_j(m)
  mean_diff <- (mean_post1 - mean_pre1) - (mean_post2 - mean_pre2)
  r_avg <- (n1 * r1 + n2 * r2) / N

  sd_change1 <- sqrt(sd_pre1^2 + sd_post1^2 - 2 * r1 * sd_pre1 * sd_post1)
  sd_change2 <- sqrt(sd_pre2^2 + sd_post2^2 - 2 * r2 * sd_pre2 * sd_post2)
  sd_pooled <- sqrt(((n1 - 1) * sd_change1^2 + (n2 - 1) * sd_change2^2) / m)

  d_manual <- (mean_diff / sd_pooled) * sqrt(2 * (1 - r_avg))
  g_manual <- d_manual * J

  var_g_manual <- J^2 * (2 * (1 - r_avg) * N / (n1 * n2) + g_manual^2 / (2 * m))

  expect_equal(res$d, d_manual, tolerance = 1e-8)
  expect_equal(res$g, g_manual, tolerance = 1e-8)
  expect_equal(res$g_se^2, var_g_manual, tolerance = 1e-8)
})

test_that("pool_sd morris_dav method: manual formula verification", {
  n1 <- 40; n2 <- 38
  mean_pre1 <- 15; mean_post1 <- 22
  sd_pre1 <- 4; sd_post1 <- 5
  mean_pre2 <- 15; mean_post2 <- 17
  sd_pre2 <- 3.5; sd_post2 <- 4.5
  r1 <- 0.5; r2 <- 0.6

  res <- es_from_means_sd_pre_post(
    n_exp = n1, n_nexp = n2,
    mean_pre_exp = mean_pre1, mean_exp = mean_post1,
    mean_pre_sd_exp = sd_pre1, mean_sd_exp = sd_post1,
    mean_pre_nexp = mean_pre2, mean_nexp = mean_post2,
    mean_pre_sd_nexp = sd_pre2, mean_sd_nexp = sd_post2,
    r_pre_post_exp = r1, r_pre_post_nexp = r2,
    pre_post_to_smd = "morris_dav",
    pool_sd = TRUE
  )

  # Manual calculation
  N <- n1 + n2
  m <- N - 2
  J <- metaConvert:::.d_j(m)
  mean_diff <- (mean_post1 - mean_pre1) - (mean_post2 - mean_pre2)
  r_avg <- (n1 * r1 + n2 * r2) / N

  sd_av1 <- sqrt((sd_pre1^2 + sd_post1^2) / 2)
  sd_av2 <- sqrt((sd_pre2^2 + sd_post2^2) / 2)
  sd_pooled <- sqrt(((n1 - 1) * sd_av1^2 + (n2 - 1) * sd_av2^2) / m)

  d_manual <- mean_diff / sd_pooled
  g_manual <- d_manual * J

  var_g_manual <- J^2 * (2 * (1 - r_avg) * N / (n1 * n2) +
                         g_manual^2 * (1 + r_avg^2) / (4 * m))

  expect_equal(res$d, d_manual, tolerance = 1e-8)
  expect_equal(res$g, g_manual, tolerance = 1e-8)
  expect_equal(res$g_se^2, var_g_manual, tolerance = 1e-8)
})

test_that("pool_sd works through mean change delegation chain", {
  res_mc <- es_from_mean_change_sd(
    n_exp = 36, n_nexp = 35,
    mean_change_exp = 8.4, mean_change_sd_exp = 3.2,
    mean_change_nexp = 2.43, mean_change_sd_nexp = 2.8,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE
  )

  # Equivalent pre_post call with mean_pre=0
  res_pp <- es_from_means_sd_pre_post(
    n_exp = 36, n_nexp = 35,
    mean_pre_exp = 0, mean_exp = 8.4,
    mean_pre_sd_exp = 0, mean_sd_exp = 3.2,
    mean_pre_nexp = 0, mean_nexp = 2.43,
    mean_pre_sd_nexp = 0, mean_sd_nexp = 2.8,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE
  )

  expect_equal(res_mc$g, res_pp$g, tolerance = 1e-10)
  expect_equal(res_mc$g_se, res_pp$g_se, tolerance = 1e-10)
})

test_that("pool_sd flows through SE delegation chain", {
  res <- es_from_mean_change_se(
    n_exp = 50, n_nexp = 50,
    mean_change_exp = 5, mean_change_se_exp = 1,
    mean_change_nexp = 2, mean_change_se_nexp = 0.8,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE
  )

  expect_false(is.na(res$g))
  expect_false(is.na(res$g_se))
  expect_true(res$g_se > 0)
})

test_that("pool_sd flows through CI delegation chain", {
  res <- es_from_mean_change_ci(
    n_exp = 50, n_nexp = 50,
    mean_change_exp = 5,
    mean_change_ci_lo_exp = 3, mean_change_ci_up_exp = 7,
    mean_change_nexp = 2,
    mean_change_ci_lo_nexp = 0.5, mean_change_ci_up_nexp = 3.5,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE
  )

  expect_false(is.na(res$g))
  expect_false(is.na(res$g_se))
  expect_true(res$g_se > 0)
})

test_that("pool_sd flows through pval delegation chain", {
  res <- es_from_mean_change_pval(
    n_exp = 50, n_nexp = 50,
    mean_change_exp = 5, mean_change_pval_exp = 0.001,
    mean_change_nexp = 2, mean_change_pval_nexp = 0.05,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "morris_dz",
    pool_sd = TRUE
  )

  expect_false(is.na(res$g))
  expect_false(is.na(res$g_se))
  expect_true(res$g_se > 0)
})

test_that("pool_sd with very unequal sample sizes", {
  res <- es_from_means_sd_pre_post(
    n_exp = 100, n_nexp = 10,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 6,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 6,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "bonett",
    pool_sd = TRUE
  )

  expect_false(is.na(res$g))
  expect_false(is.na(res$g_se))
  expect_true(res$g_se > 0)
})

test_that("pool_sd with very unequal correlations", {
  res <- es_from_means_sd_pre_post(
    n_exp = 50, n_nexp = 50,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 6,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 6,
    r_pre_post_exp = 0.1, r_pre_post_nexp = 0.9,
    pre_post_to_smd = "bonett",
    pool_sd = TRUE
  )

  expect_false(is.na(res$g))
  expect_false(is.na(res$g_se))
  expect_true(res$g_se > 0)
  expect_true(res$g_ci_lo < res$g)
  expect_true(res$g_ci_up > res$g)
})

test_that("pool_sd bonett with small N (edge case m <= 2)", {
  # N = 4 gives m = 2, bonett variance has m/(m-2) which is undefined
  res <- es_from_means_sd_pre_post(
    n_exp = 2, n_nexp = 2,
    mean_pre_exp = 20, mean_exp = 25,
    mean_pre_sd_exp = 5, mean_sd_exp = 6,
    mean_pre_nexp = 20, mean_nexp = 21,
    mean_pre_sd_nexp = 5, mean_sd_nexp = 6,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.5,
    pre_post_to_smd = "bonett",
    pool_sd = TRUE
  )

  # d should still be computed, but variance should be NA
  expect_false(is.na(res$d))
  expect_true(is.na(res$g_se))
})

test_that("pool_sd works through convert_df pipeline", {
  dat <- data.frame(
    n_exp = c(36, 50),
    n_nexp = c(35, 50),
    mean_pre_exp = c(20, 30),
    mean_exp = c(28, 35),
    mean_pre_sd_exp = c(5, 7),
    mean_sd_exp = c(6, 8),
    mean_pre_nexp = c(20, 30),
    mean_nexp = c(22, 32),
    mean_pre_sd_nexp = c(4, 6),
    mean_sd_nexp = c(5, 7),
    r_pre_post_exp = c(0.6, 0.5),
    r_pre_post_nexp = c(0.7, 0.5)
  )

  res_unpooled <- convert_df(dat, measure = "g", verbose = FALSE,
                              pre_post_to_smd = "bonett", pool_sd = FALSE)
  res_pooled <- convert_df(dat, measure = "g", verbose = FALSE,
                            pre_post_to_smd = "bonett", pool_sd = TRUE)

  sum_unpooled <- summary(res_unpooled)
  sum_pooled <- summary(res_pooled)

  expect_true(nrow(sum_unpooled) > 0)
  expect_true(nrow(sum_pooled) > 0)
})

test_that("pool_sd works with all 4 methods in pre_post_to_smd", {
  methods <- c("bonett", "morris_dz", "morris_drm", "morris_dav")

  for (method in methods) {
    res <- es_from_means_sd_pre_post(
      n_exp = 40, n_nexp = 38,
      mean_pre_exp = 15, mean_exp = 22,
      mean_pre_sd_exp = 4, mean_sd_exp = 5,
      mean_pre_nexp = 15, mean_nexp = 17,
      mean_pre_sd_nexp = 3.5, mean_sd_nexp = 4.5,
      r_pre_post_exp = 0.5, r_pre_post_nexp = 0.6,
      pre_post_to_smd = method,
      pool_sd = TRUE
    )

    expect_false(is.na(res$g), info = paste("Method:", method))
    expect_false(is.na(res$g_se), info = paste("Method:", method))
    expect_true(res$g_se > 0, info = paste("Method:", method))
    expect_true(res$g_ci_lo < res$g, info = paste("Method:", method))
    expect_true(res$g_ci_up > res$g, info = paste("Method:", method))
  }
})

test_that("pool_sd cooper alias works", {
  res_cooper <- es_from_means_sd_pre_post(
    n_exp = 40, n_nexp = 38,
    mean_pre_exp = 15, mean_exp = 22,
    mean_pre_sd_exp = 4, mean_sd_exp = 5,
    mean_pre_nexp = 15, mean_nexp = 17,
    mean_pre_sd_nexp = 3.5, mean_sd_nexp = 4.5,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.6,
    pre_post_to_smd = "cooper",
    pool_sd = TRUE
  )

  res_drm <- es_from_means_sd_pre_post(
    n_exp = 40, n_nexp = 38,
    mean_pre_exp = 15, mean_exp = 22,
    mean_pre_sd_exp = 4, mean_sd_exp = 5,
    mean_pre_nexp = 15, mean_nexp = 17,
    mean_pre_sd_nexp = 3.5, mean_sd_nexp = 4.5,
    r_pre_post_exp = 0.5, r_pre_post_nexp = 0.6,
    pre_post_to_smd = "morris_drm",
    pool_sd = TRUE
  )

  expect_equal(res_cooper$g, res_drm$g, tolerance = 1e-10)
  expect_equal(res_cooper$g_se, res_drm$g_se, tolerance = 1e-10)
})
