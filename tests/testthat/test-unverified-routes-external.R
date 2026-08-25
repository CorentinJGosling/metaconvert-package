# =============================================================================
# Four exported routes had NO test naming them anywhere -- not in tests/testthat/,
# not in tests_save/checked/, not in simulations/, not in a vignette:
#
#   es_from_2x2_prop  es_from_pt_bis_r_pval  es_from_rr_pval  es_from_student_t_pval
#
# They were not literally unexercised: convert_df() reaches all 90 exported routes,
# and tests_save/checked/test-ES-UNTESTED-PATHS.R drives some of them through it. But
# that file says what it is in its own header -- "these are *not* formula-correctness
# tests" -- and the coverage it provides is that the route returns something finite.
#
# The gap that leaves is specific and worth naming: everything else pinning these
# routes compares metaConvert to metaConvert. Cross-route consistency cannot detect an
# error the routes SHARE, and three of these four delegate to a common callee, so a
# mistake in the shared step would be invisible to every existing assertion.
#
# So these tests deliberately use INDEPENDENT references only: base R's own qt/qnorm
# for the distributional inversions, and metafor / esc / compute.es for the effect
# sizes. Where metaConvert and a reference legitimately disagree, the test asserts the
# disagreement and its cause, rather than loosening a tolerance until it passes.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. es_from_student_t_pval -- p-value -> t -> SMD
# -----------------------------------------------------------------------------

test_that("es_from_student_t_pval inverts the p-value with base R's own qt()", {
  # The route's ONLY arithmetic is the inversion; everything after is delegated. So
  # the inversion is what has to be checked against something outside the package.
  for (p in c(0.5, 0.05, 0.013, 1e-4)) {
    for (nn in list(c(40, 55), c(12, 80), c(200, 200))) {
      n1 <- nn[1]; n2 <- nn[2]
      t_ref <- stats::qt(p / 2, df = n1 + n2 - 2, lower.tail = FALSE)
      got <- es_from_student_t_pval(student_t_pval = p, n_exp = n1, n_nexp = n2)
      ref <- es_from_student_t(student_t = t_ref, n_exp = n1, n_nexp = n2)
      expect_equal(got$d, ref$d, tolerance = 1e-14)
      expect_equal(got$d_se, ref$d_se, tolerance = 1e-14)
      expect_equal(got$info_used, "student_t_pval")
    }
  }
})


test_that("es_from_student_t_pval agrees with esc::esc_t on d and its SE", {
  skip_if_not_installed("esc")
  p <- 0.013; n1 <- 40; n2 <- 55
  t <- stats::qt(p / 2, df = n1 + n2 - 2, lower.tail = FALSE)
  got <- es_from_student_t_pval(student_t_pval = p, n_exp = n1, n_nexp = n2)
  ref <- esc::esc_t(t = t, grp1n = n1, grp2n = n2, es.type = "d")
  # Bit-exact: both compute d = t * sqrt(1/n1 + 1/n2) and the Borenstein variance.
  expect_equal(got$d, as.numeric(ref$es), tolerance = 1e-12)
  expect_equal(got$d_se, as.numeric(ref$se), tolerance = 1e-12)
})


test_that("es_from_student_t_pval agrees with compute.es::tes on d and Var(d)", {
  skip_if_not_installed("compute.es")
  p <- 0.013; n1 <- 40; n2 <- 55
  t <- stats::qt(p / 2, df = n1 + n2 - 2, lower.tail = FALSE)
  got <- es_from_student_t_pval(student_t_pval = p, n_exp = n1, n_nexp = n2)
  ref <- compute.es::tes(t = t, n.1 = n1, n.2 = n2, verbose = FALSE, dig = 12)
  expect_equal(got$d, ref$d, tolerance = 1e-10)
  expect_equal(got$d_se^2, ref$var.d, tolerance = 1e-10)
})


test_that("the g gap against esc is the EXACT Hedges J, not an error", {
  skip_if_not_installed("esc")
  # esc uses the approximation J ~= 1 - 3/(4*df - 1); metaConvert uses the exact
  # gamma-function form. They differ in the 6th decimal, and metaConvert is the more
  # accurate of the two. Asserting this keeps the difference from being "fixed" toward
  # the approximation by someone who sees only a failing comparison.
  p <- 0.013; n1 <- 40; n2 <- 55; df <- n1 + n2 - 2
  t <- stats::qt(p / 2, df = df, lower.tail = FALSE)
  d <- t * sqrt(1 / n1 + 1 / n2)
  J_exact  <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  J_approx <- 1 - 3 / (4 * df - 1)

  got <- es_from_student_t_pval(student_t_pval = p, n_exp = n1, n_nexp = n2)
  expect_equal(got$g, d * J_exact, tolerance = 1e-12)

  ref_g <- as.numeric(esc::esc_t(t = t, grp1n = n1, grp2n = n2, es.type = "g")$es)
  expect_equal(ref_g, d * J_approx, tolerance = 1e-12)
  # The whole discrepancy is J, and it is small and one-directional.
  expect_lt(abs(got$g - ref_g), 1e-5)
  expect_lt(J_exact, J_approx)
})


# -----------------------------------------------------------------------------
# 2. es_from_pt_bis_r_pval -- the same inversion, reached from a correlation's p
# -----------------------------------------------------------------------------

test_that("es_from_pt_bis_r_pval matches es_from_student_t_pval except for info_used", {
  # A point-biserial r and a two-sample t are the same test statistic on the same data,
  # so the two routes SHOULD coincide. Pinning it makes that a decision rather than an
  # accident, since the two bodies are separate copies of the same three lines.
  p <- 0.013; n1 <- 40; n2 <- 55
  a <- es_from_pt_bis_r_pval(pt_bis_r_pval = p, n_exp = n1, n_nexp = n2)
  b <- es_from_student_t_pval(student_t_pval = p, n_exp = n1, n_nexp = n2)
  shared <- setdiff(intersect(names(a), names(b)), "info_used")
  for (cl in shared) expect_equal(a[[cl]], b[[cl]], tolerance = 1e-14, info = cl)
  expect_equal(a$info_used, "pt_bis_r_pval")
})


test_that("smd_to_cor = 'lipsey_cooper' returns t/sqrt(t^2 + N), NOT the exact point-biserial", {
  # Cooper's published conversion is r = d / sqrt(d^2 + a) with a = (n1+n2)^2/(n1*n2),
  # which reduces algebraically to t / sqrt(t^2 + N). The EXACT point-biserial is
  # t / sqrt(t^2 + N - 2). metaConvert implements the published formula faithfully, so
  # the N-vs-df difference is inherited from the source, not introduced here.
  #
  # This is pinned because the difference is systematic, grows as N shrinks, and is
  # easy to mistake for a bug: a reader who checks the route against the textbook
  # point-biserial identity will find a mismatch and may "fix" a correct implementation.
  p <- 0.013
  for (nn in list(c(30, 30), c(50, 50), c(40, 55), c(12, 80))) {
    n1 <- nn[1]; n2 <- nn[2]; N <- n1 + n2; df <- N - 2
    t <- stats::qt(p / 2, df = df, lower.tail = FALSE)
    got <- es_from_pt_bis_r_pval(pt_bis_r_pval = p, n_exp = n1, n_nexp = n2,
                                 smd_to_cor = "lipsey_cooper")$r
    expect_equal(got, t / sqrt(t^2 + N), tolerance = 1e-14,
                 info = paste(n1, n2))
    # and it is NOT the df-based point-biserial -- a real, non-vanishing gap
    expect_gt(abs(got - t / sqrt(t^2 + df)), 1e-4)
  }
})


# -----------------------------------------------------------------------------
# 3. es_from_rr_pval -- recover SE(log RR) from RR + p
# -----------------------------------------------------------------------------

test_that("es_from_rr_pval round-trips metafor's own RR and Wald p-value exactly", {
  # The strongest available check: metafor produces log RR and its variance from a 2x2,
  # the Wald p-value follows, and feeding RR + that p back must return metafor's SE.
  # Nothing in this chain is metaConvert's own arithmetic except the inversion.
  tabs <- list(c(30, 70, 15, 85), c(12, 88, 5, 95), c(60, 40, 45, 55), c(5, 195, 2, 198))
  for (tb in tabs) {
    a <- tb[1]; b <- tb[2]; cc <- tb[3]; dd <- tb[4]
    ef <- metafor::escalc(measure = "RR", ai = a, bi = b, ci = cc, di = dd)
    logrr <- as.numeric(ef$yi[1]); se <- sqrt(as.numeric(ef$vi[1]))
    pv <- 2 * stats::pnorm(abs(logrr) / se, lower.tail = FALSE)

    got <- es_from_rr_pval(rr = exp(logrr), rr_pval = pv,
                           n_exp = a + b, n_nexp = cc + dd,
                           n_cases = a + cc, n_controls = b + dd)
    expect_equal(got$logrr, logrr, tolerance = 1e-12, info = paste(tb, collapse = "/"))
    expect_equal(got$logrr_se, se, tolerance = 1e-12, info = paste(tb, collapse = "/"))
  }
})


test_that("es_from_rr_pval returns NA -- not NaN, and not 0 -- for the SE at rr = 1", {
  # WRITTEN FROM MEASUREMENT, AGAINST THE CODE COMMENT. The comment above the formula
  # in R/es_from_stand_RR.R says this case "gives logrr_se = 0 (finite)". That is true
  # of the INTERMEDIATE -- |log(1)| / z is exactly 0 -- but not of the return value:
  # es_from_rr_se() then applies .positive_or_na(), and a zero SE is non-positive, so
  # the user sees NA.
  #
  # NA is the right answer. An SE of exactly 0 is an infinite inverse-variance weight,
  # which is precisely what that guard exists to prevent; and a p-value of 0.9 beside
  # rr = 1 carries no information about the standard error, so there is nothing to
  # report. The pair of assertions below distinguishes NA from NaN deliberately --
  # is.na() is TRUE for both, so only adding is.nan() pins which one is returned, and
  # NaN is what the old sign()-based denominator produced.
  got <- es_from_rr_pval(rr = 1, rr_pval = 0.9, n_exp = 100, n_nexp = 100)
  expect_equal(got$logrr, 0, tolerance = 1e-14)
  expect_true(is.na(got$logrr_se))
  expect_false(is.nan(got$logrr_se))

  # The comment also claims this mirrors es_from_or_pval. It does -- checked, not assumed.
  o <- es_from_or_pval(or = 1, or_pval = 0.9, n_exp = 100, n_nexp = 100)
  expect_true(is.na(o$logor_se))
  expect_false(is.nan(o$logor_se))
})


# -----------------------------------------------------------------------------
# 4. es_from_2x2_prop -- proportions -> counts -> 2x2
# -----------------------------------------------------------------------------

test_that("es_from_2x2_prop reproduces metafor on the cells its rounding implies", {
  # The route's own arithmetic is one rounding step; everything after is es_from_2x2.
  # So the check is: do the implied cells give metafor's log OR / log RR / RD exactly?
  # (r and z are deliberately not asserted here -- they route through the tetrachoric
  # solve, which needs mvtnorm, a Suggests.)
  grid <- expand.grid(pe = c(0.30, 0.55, 0.08), pn = c(0.15, 0.40),
                      N1 = c(100, 63), N2 = c(100, 87))
  for (i in seq_len(nrow(grid))) {
    pe <- grid$pe[i]; pn <- grid$pn[i]; N1 <- grid$N1[i]; N2 <- grid$N2[i]
    a <- round(pe * N1); cc <- round(pn * N2)
    b <- N1 - a; dd <- N2 - cc
    lbl <- sprintf("%.2f/%.2f %d/%d", pe, pn, N1, N2)

    got <- es_from_2x2_prop(prop_cases_exp = pe, prop_cases_nexp = pn,
                            n_exp = N1, n_nexp = N2)
    mo <- metafor::escalc(measure = "OR", ai = a, bi = b, ci = cc, di = dd)
    mr <- metafor::escalc(measure = "RR", ai = a, bi = b, ci = cc, di = dd)
    expect_equal(got$logor, as.numeric(mo$yi[1]), tolerance = 1e-12, info = lbl)
    expect_equal(got$logor_se, sqrt(as.numeric(mo$vi[1])), tolerance = 1e-12, info = lbl)
    expect_equal(got$logrr, as.numeric(mr$yi[1]), tolerance = 1e-12, info = lbl)
    expect_equal(got$logrr_se, sqrt(as.numeric(mr$vi[1])), tolerance = 1e-12, info = lbl)
  }
})


test_that("es_from_2x2_prop rounds the implied cells to integers, half away from zero", {
  # The rounding IS the route, and it is lossy: two different proportions can land on
  # the same table. Pinned so a change to the rule (floor, or keeping fractional cells)
  # is a visible decision -- fractional cells would change every downstream variance.
  got <- es_from_2x2_prop(prop_cases_exp = 1/3, prop_cases_nexp = 0.15,
                          n_exp = 30, n_nexp = 100)
  ref <- es_from_2x2(n_cases_exp = 10, n_controls_exp = 20,
                     n_cases_nexp = 15, n_controls_nexp = 85)
  expect_equal(got$logor, ref$logor, tolerance = 1e-14)   # 9.99 -> 10, not 9

  # a proportion that cannot be met exactly is silently snapped, and the route reports
  # the snapped table, not the requested proportion
  snapped <- es_from_2x2_prop(prop_cases_exp = 0.334, prop_cases_nexp = 0.15,
                              n_exp = 30, n_nexp = 100)
  expect_equal(snapped$logor, got$logor, tolerance = 1e-14)
  expect_equal(got$info_used, "2x2_prop")
})


# -----------------------------------------------------------------------------
# 5. es_from_cohen_d_adj -- the identity its @details names
# -----------------------------------------------------------------------------
# The @details block used to say the point estimate maps "via the point-biserial
# identity at n_exp + n_nexp - 2". Measured: it is the N-based Cooper form. The
# n_exp + n_nexp - 2 - n_cov_ancova df set the CONFIDENCE INTERVAL, not the point
# estimate. The block now says so, and quotes two numbers; these tests pin all of it,
# so the documentation cannot drift away from the code again.

test_that("es_from_cohen_d_adj's lipsey_cooper r uses the TOTAL N, not N - 2", {
  d <- 0.6
  for (nn in list(c(30, 30), c(40, 55), c(12, 80))) {
    n1 <- nn[1]; n2 <- nn[2]; N <- n1 + n2
    a_N <- N^2 / (n1 * n2)
    lbl <- paste(n1, n2)
    for (ncov in c(0, 3)) {
      got <- es_from_cohen_d_adj(cohen_d_adj = d, n_exp = n1, n_nexp = n2,
                                 n_cov_ancova = ncov, cov_outcome_r = 0.5,
                                 smd_to_cor = "lipsey_cooper")$r
      expect_equal(got, d / sqrt(d^2 + a_N), tolerance = 1e-13, info = lbl)
    }
    # and it is NOT the sample point-biserial identity, which the old text named
    t_stat <- d / sqrt(1 / n1 + 1 / n2)
    expect_gt(abs(d / sqrt(d^2 + a_N) - t_stat / sqrt(t_stat^2 + N - 2)), 1e-4)
  }
})


test_that("the correlation point estimate really is unchanged by n_cov_ancova", {
  # The half of the old sentence that WAS true -- kept, because it is the reason the
  # sentence exists and a future edit could easily lose it while fixing the other half.
  for (cor_opt in c("viechtbauer", "lipsey_cooper")) {
    r0 <- es_from_cohen_d_adj(cohen_d_adj = 0.6, n_exp = 40, n_nexp = 55,
                              n_cov_ancova = 0, cov_outcome_r = 0.5,
                              smd_to_cor = cor_opt)$r
    r3 <- es_from_cohen_d_adj(cohen_d_adj = 0.6, n_exp = 40, n_nexp = 55,
                              n_cov_ancova = 3, cov_outcome_r = 0.5,
                              smd_to_cor = cor_opt)$r
    expect_equal(r0, r3, tolerance = 1e-14, info = cor_opt)
  }
})


test_that("the two smd_to_cor options return the values ?es_from_cohen_d_adj quotes", {
  # The @details block states 0.3658 (viechtbauer, the biserial) and 0.2873
  # (lipsey_cooper) at n = 30/30, cohen_d_adj = 0.6. Numbers printed in documentation
  # are claims like any other, so they get an assertion.
  v <- es_from_cohen_d_adj(cohen_d_adj = 0.6, n_exp = 30, n_nexp = 30,
                           n_cov_ancova = 0, cov_outcome_r = 0.5)$r
  l <- es_from_cohen_d_adj(cohen_d_adj = 0.6, n_exp = 30, n_nexp = 30,
                           n_cov_ancova = 0, cov_outcome_r = 0.5,
                           smd_to_cor = "lipsey_cooper")$r
  expect_equal(round(v, 4), 0.3658)
  expect_equal(round(l, 4), 0.2873)
  expect_gt(abs(v - l), 0.05)   # the two options are genuinely different estimands
})
