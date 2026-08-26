# =============================================================================
# The paired-t / paired-F routes used to hold a SECOND COPY of the morris_dz and
# morris_drm arithmetic that .single_group_pre_post_to_smd already computes. The two
# agreed -- measured at 0 for morris_dz and 1.1e-16 for morris_drm -- but only because
# both were maintained in step. The old comment said so in as many words: "matching
# .single_group_pre_post_to_smd so the paired-t and mean-change routes return identical
# SEs on equivalent inputs". That is a promise, not a mechanism.
#
# They now delegate through .paired_t_to_smd(), so there is one implementation. These
# tests pin the things a delegation could quietly break -- none of which a comparison of
# ordinary numbers would have caught:
#
#   * the per-row VECTOR of methods (the kernel dispatches on a scalar `if`);
#   * the n < 2 guard;
#   * the warning stream, which is user-visible output;
#   * the closed forms the synthetic-input trick has to reproduce.
# =============================================================================

test_that("the delegated routes still produce the published closed forms", {
  # d_z = t / sqrt(n); d_rm = t * sqrt(2(1-r)/n). These are what the inline copy wrote
  # out, and what the kernel must return for the synthetic problem it is handed.
  for (n in c(2, 12, 40, 200)) {
    for (r in c(0.1, 0.5, 0.9)) {
      t1 <- 4.2; t2 <- 1.8; n2 <- 30
      dz <- es_from_paired_t(paired_t_exp = t1, paired_t_nexp = t2, n_exp = n,
                             n_nexp = n2, r_pre_post_exp = r, r_pre_post_nexp = r,
                             pre_post_to_smd = "morris_dz")
      expect_equal(dz$d, t1 / sqrt(n) - t2 / sqrt(n2), tolerance = 1e-12,
                   info = paste("dz", n, r))

      drm <- es_from_paired_t(paired_t_exp = t1, paired_t_nexp = t2, n_exp = n,
                              n_nexp = n2, r_pre_post_exp = r, r_pre_post_nexp = r,
                              pre_post_to_smd = "morris_drm")
      expect_equal(drm$d,
                   t1 * sqrt(2 * (1 - r) / n) - t2 * sqrt(2 * (1 - r) / n2),
                   tolerance = 1e-12, info = paste("drm", n, r))
    }
  }
})


test_that("a per-row VECTOR of methods still selects per row", {
  # The kernel dispatches on a SCALAR pre_post_to_smd, so a straight vectorised
  # delegation would have dropped every row into whichever branch the first element
  # picked -- silently, and only for users who mix methods within one call. This is the
  # single most likely way the refactor could have gone wrong.
  n1 <- c(20, 30, 40); n2 <- c(22, 33, 44)
  t1 <- c(2, 3, 4);    t2 <- c(1, 1.5, 2)
  r  <- c(.5, .6, .7)
  mixed <- es_from_paired_t(paired_t_exp = t1, paired_t_nexp = t2, n_exp = n1,
                            n_nexp = n2, r_pre_post_exp = r, r_pre_post_nexp = r,
                            pre_post_to_smd = c("morris_dz", "cooper", "morris_dz"))
  for (i in 1:3) {
    m <- c("morris_dz", "cooper", "morris_dz")[i]
    one <- es_from_paired_t(paired_t_exp = t1[i], paired_t_nexp = t2[i],
                            n_exp = n1[i], n_nexp = n2[i],
                            r_pre_post_exp = r[i], r_pre_post_nexp = r[i],
                            pre_post_to_smd = m)
    expect_equal(mixed$d[i], one$d, tolerance = 1e-12, info = paste(i, m))
    expect_equal(mixed$d_se[i], one$d_se, tolerance = 1e-12, info = paste(i, m))
  }
  # and the mixed call is genuinely mixed: row 2 differs from what dz would give
  all_dz <- es_from_paired_t(paired_t_exp = t1, paired_t_nexp = t2, n_exp = n1,
                             n_nexp = n2, r_pre_post_exp = r, r_pre_post_nexp = r,
                             pre_post_to_smd = "morris_dz")
  expect_false(isTRUE(all.equal(mixed$d[2], all_dz$d[2], tolerance = 1e-8)))
})


test_that("the n < 2 guard survives delegation, on both paired entry points", {
  # A within-subject arm with n < 2 has no estimable variance. morris_drm carries no
  # J(n-1) factor, so before the guard existed it stayed finite at n = 1 and leaked a
  # spurious d/g/logOR/r downstream. The guard now lives in the kernel; these assert it
  # still reaches both routes.
  for (m in c("morris_dz", "cooper")) {
    for (n in c(0, 1)) {
      two <- suppressWarnings(es_from_paired_t(
        paired_t_exp = 3, paired_t_nexp = 2, n_exp = n, n_nexp = 30,
        r_pre_post_exp = .5, r_pre_post_nexp = .5, pre_post_to_smd = m))
      expect_true(is.na(two$d), info = paste("two-group", m, n))
      expect_true(is.na(two$d_se), info = paste("two-group", m, n))

      one <- suppressWarnings(es_from_paired_t_single_group(
        paired_t_exp = 3, n_exp = n, r_pre_post_exp = .5, pre_post_to_smd = m))
      expect_true(is.na(one$d), info = paste("single-group", m, n))
      expect_true(is.na(one$d_se), info = paste("single-group", m, n))
    }
    # n = 2 is estimable and must NOT be caught by the guard
    ok <- es_from_paired_t(paired_t_exp = 3, paired_t_nexp = 2, n_exp = 2, n_nexp = 30,
                           r_pre_post_exp = .5, r_pre_post_nexp = .5,
                           pre_post_to_smd = m)
    expect_false(is.na(ok$d), info = paste("n = 2", m))
  }
})


test_that("the informative Hedges' J warning survives, and normal input stays silent", {
  # WARNING STREAM IS OUTPUT. Delegation moved where .d_j() is called, and skipping the
  # n < 2 rows meant the kernel never reached it for them -- which silently dropped the
  # "Hedges' J is undefined for df <= 1" warning that tells the caller g will be NA.
  # .paired_t_to_smd() calls .d_j() explicitly to keep it. Neither a value comparison
  # nor a NA-pattern comparison would have noticed.
  w <- NULL
  withCallingHandlers(
    invisible(es_from_paired_t(paired_t_exp = 3, paired_t_nexp = 2, n_exp = 1,
                               n_nexp = 30, pre_post_to_smd = "cooper")),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("Hedges' J correction is undefined", w, fixed = TRUE)))

  expect_silent(es_from_paired_t(paired_t_exp = 3, paired_t_nexp = 2, n_exp = 30,
                                 n_nexp = 30, pre_post_to_smd = "cooper"))
  expect_silent(es_from_paired_t_single_group(paired_t_exp = 3, n_exp = 30,
                                              pre_post_to_smd = "cooper"))
})


test_that("cooper stays an exact alias of morris_drm, and the other three are refused", {
  a <- es_from_paired_t(paired_t_exp = 3, paired_t_nexp = 2, n_exp = 30, n_nexp = 30,
                        pre_post_to_smd = "cooper")$d
  b <- es_from_paired_t(paired_t_exp = 3, paired_t_nexp = 2, n_exp = 30, n_nexp = 30,
                        pre_post_to_smd = "morris_drm")$d
  expect_identical(a, b)

  # The kernel supports five methods; these routes support two, because a paired t does
  # not identify the separate pre/post SDs the other three need. Delegation must not
  # quietly widen the menu.
  for (m in c("bonett", "morris_dav")) {
    expect_error(es_from_paired_t(paired_t_exp = 3, paired_t_nexp = 2, n_exp = 30,
                                  n_nexp = 30, pre_post_to_smd = m), "pre_post_to_smd")
    expect_error(es_from_paired_t_single_group(paired_t_exp = 3, n_exp = 30,
                                               pre_post_to_smd = m), "pre_post_to_smd")
  }
})


test_that("paired_t and mean_change agree on equivalent inputs, which is now structural", {
  # This equivalence was the reason the duplicated copy existed. It used to hold by
  # maintenance; it now holds because there is one implementation. Kept as the
  # end-to-end statement of what the refactor preserves.
  n1 <- 40; n2 <- 42; mc1 <- 3; sc1 <- 1.4; mc2 <- 0.9; sc2 <- 1.5; r <- .8
  t1 <- mc1 / (sc1 / sqrt(n1)); t2 <- mc2 / (sc2 / sqrt(n2))
  for (m in c("morris_dz", "morris_drm")) {
    a <- es_from_mean_change_sd(mean_change_exp = mc1, mean_change_sd_exp = sc1,
                                mean_change_nexp = mc2, mean_change_sd_nexp = sc2,
                                n_exp = n1, n_nexp = n2, r_pre_post_exp = r,
                                r_pre_post_nexp = r, pre_post_to_smd = m)
    b <- es_from_paired_t(paired_t_exp = t1, paired_t_nexp = t2, n_exp = n1,
                          n_nexp = n2, r_pre_post_exp = r, r_pre_post_nexp = r,
                          pre_post_to_smd = m)
    for (cl in c("d", "d_se", "g", "g_se")) {
      expect_equal(a[[cl]], b[[cl]], tolerance = 1e-12, info = paste(m, cl))
    }
  }
})
