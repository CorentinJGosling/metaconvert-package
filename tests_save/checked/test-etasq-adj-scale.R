# P1 (audit digest [8]): es_from_etasq_adj must emit an adjusted d on the
# MARGINAL (unadjusted-SD) scale, obtained by converting the partial
# eta-squared to the ANCOVA F it implies (eta_p^2 = F / (F + df_err)) and
# routing through the es_from_ancova_f computation path. Before the fix,
# the point estimate d = 2*sqrt(etasq_adj/(1-etasq_adj)) was on the RESIDUAL
# SD scale while the variance was the marginal-scale Cooper 12.3 form, so the
# SE understated its own estimator's sampling SD (~45% at r = 0.7) and the
# route contradicted es_from_ancova_f fed the algebraically equivalent
# statistic (d ratio = 1/sqrt(1-r^2) * sqrt(N/df_err)).

library(testthat)
library(metaConvert)

# Hedges' small-sample correction J, hand-derived (Hedges 1981), written
# independently of the package internals.
.J_hand <- function(df) exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))

# P1 test (1): cross-route identity with es_from_ancova_f -------------------
# Oracle: the ancova_f route, independently verified in the audit against
# Cooper (2019) table 12.3 and by simulation (95% CI coverage 0.945).
test_that("P1: es_from_etasq_adj(F/(F+df_err)) equals es_from_ancova_f(F) exactly", {
  ns <- list(c(20, 20), c(15, 25))
  cols <- c(
    "d", "d_se", "d_ci_lo", "d_ci_up",
    "g", "g_se", "g_ci_lo", "g_ci_up"
  )
  for (f in c(2, 4, 9)) {
    for (nn in ns) {
      for (r in c(0.3, 0.5, 0.7)) {
        for (q in c(1, 3)) {
          n1 <- nn[1]
          n2 <- nn[2]
          df_err <- n1 + n2 - 2 - q
          eta_p <- f / (f + df_err) # partial eta-squared implied by F

          res_eta <- es_from_etasq_adj(
            etasq_adj = eta_p, n_exp = n1, n_nexp = n2,
            n_cov_ancova = q, cov_outcome_r = r
          )
          res_f <- es_from_ancova_f(
            ancova_f = f, cov_outcome_r = r, n_cov_ancova = q,
            n_exp = n1, n_nexp = n2
          )

          lbl <- sprintf("F=%g, n=(%d,%d), r=%.1f, n_cov=%d", f, n1, n2, r, q)
          for (col in cols) {
            expect_equal(res_eta[[col]], res_f[[col]],
              tolerance = 1e-10,
              label = paste0("etasq_adj route ", col, " [", lbl, "]"),
              expected.label = paste0("ancova_f route ", col)
            )
          }
          expect_equal(res_eta$info_used, "etasq_adj")
        }
      }
    }
  }
})

# P1 test (2): the emitted d and variance are mutually scale-consistent -----
# Hand-derived formulas (Cooper 2019, table 12.3), not package recomputes:
#   d       = sqrt(F_implied) * sqrt(1/n1 + 1/n2) * sqrt(1 - r^2)
#   d_se^2  = (1 - r^2) * (1/n1 + 1/n2) + d^2 / (2N)   [default smd_var]
test_that("P1: emitted d_se^2 = (1-r^2)*(1/n1+1/n2) + d^2/(2N) with the emitted d", {
  n1 <- 20
  n2 <- 20
  q <- 1
  r <- 0.5
  f <- 4
  df_err <- n1 + n2 - 2 - q # 37
  eta_p <- f / (f + df_err) # 0.097561

  res <- es_from_etasq_adj(
    etasq_adj = eta_p, n_exp = n1, n_nexp = n2,
    n_cov_ancova = q, cov_outcome_r = r
  )

  # Marginal-scale point estimate, hand formula (digest [8]: 0.547723)
  d_hand <- sqrt(f) * sqrt(1 / n1 + 1 / n2) * sqrt(1 - r^2)
  expect_equal(res$d, d_hand, tolerance = 1e-10)

  # Variance identity with the emitted d (default smd_var = "borenstein")
  v_expected <- (1 - r^2) * (1 / n1 + 1 / n2) + res$d^2 / (2 * (n1 + n2))
  expect_equal(res$d_se^2, v_expected, tolerance = 1e-10)

  # Same identity at r = 0.7, where the pre-fix SE was ~45% anti-conservative
  r2 <- 0.7
  res2 <- es_from_etasq_adj(
    etasq_adj = eta_p, n_exp = n1, n_nexp = n2,
    n_cov_ancova = q, cov_outcome_r = r2
  )
  d_hand2 <- sqrt(f) * sqrt(1 / n1 + 1 / n2) * sqrt(1 - r2^2)
  expect_equal(res2$d, d_hand2, tolerance = 1e-10)
  v_expected2 <- (1 - r2^2) * (1 / n1 + 1 / n2) + res2$d^2 / (2 * (n1 + n2))
  expect_equal(res2$d_se^2, v_expected2, tolerance = 1e-10)
})

# P1 test (3): reverse_etasq flips d and swaps/negates CI bounds ------------
test_that("P1: reverse_etasq flips d/g and mirrors the CI bounds", {
  args <- list(
    etasq_adj = 0.15, n_exp = 18, n_nexp = 22,
    n_cov_ancova = 2, cov_outcome_r = 0.4
  )
  res_fwd <- do.call(es_from_etasq_adj, args)
  res_rev <- do.call(es_from_etasq_adj, c(args, list(reverse_etasq = TRUE)))

  expect_equal(res_rev$d, -res_fwd$d, tolerance = 1e-10)
  expect_equal(res_rev$d_se, res_fwd$d_se, tolerance = 1e-10)
  expect_equal(res_rev$d_ci_lo, -res_fwd$d_ci_up, tolerance = 1e-10)
  expect_equal(res_rev$d_ci_up, -res_fwd$d_ci_lo, tolerance = 1e-10)

  expect_equal(res_rev$g, -res_fwd$g, tolerance = 1e-10)
  expect_equal(res_rev$g_se, res_fwd$g_se, tolerance = 1e-10)
  expect_equal(res_rev$g_ci_lo, -res_fwd$g_ci_up, tolerance = 1e-10)
  expect_equal(res_rev$g_ci_up, -res_fwd$g_ci_lo, tolerance = 1e-10)
})

# P1 test (4): cov_outcome_r = NA -> the marginal scale is unidentified -----
# The function must degrade to NA, not silently emit a residual-scale value.
test_that("P1: cov_outcome_r = NA yields all-NA d/d_se/CI (no residual-scale leak)", {
  res_na <- es_from_etasq_adj(
    etasq_adj = 0.2, n_exp = 20, n_nexp = 20,
    n_cov_ancova = 1, cov_outcome_r = NA
  )
  expect_true(is.na(res_na$d))
  expect_true(is.na(res_na$d_se))
  expect_true(is.na(res_na$d_ci_lo))
  expect_true(is.na(res_na$d_ci_up))
  expect_true(is.na(res_na$g))
  expect_true(is.na(res_na$g_se))
  expect_true(is.na(res_na$g_ci_lo))
  expect_true(is.na(res_na$g_ci_up))
})

# P1 test (5): crude es_from_etasq must remain UNCHANGED --------------------
# Regression pin, hand-derived (Cohen 1988): d = 2*sqrt(eta/(1-eta)), crude
# variance v_d = 1/n1 + 1/n2 + d^2/(2N), df = n1 + n2 - 2.
test_that("P1: crude es_from_etasq regression pin (standard Cohen conversion)", {
  eta <- 0.28
  n1 <- 20
  n2 <- 22
  res <- es_from_etasq(etasq = eta, n_exp = n1, n_nexp = n2)

  d_hand <- 2 * sqrt(eta / (1 - eta))
  v_hand <- (1 / n1 + 1 / n2) + d_hand^2 / (2 * (n1 + n2))
  df <- n1 + n2 - 2
  J <- .J_hand(df)

  # literal numeric pin of the current crude value
  expect_equal(res$d, 1.247219129, tolerance = 1e-9)
  expect_equal(res$d, d_hand, tolerance = 1e-10)
  expect_equal(res$d_se, sqrt(v_hand), tolerance = 1e-10)
  expect_equal(res$d_ci_lo, d_hand - qt(.975, df) * sqrt(v_hand), tolerance = 1e-10)
  expect_equal(res$d_ci_up, d_hand + qt(.975, df) * sqrt(v_hand), tolerance = 1e-10)
  expect_equal(res$g, d_hand * J, tolerance = 1e-10)
  expect_equal(res$g_se, sqrt(v_hand) * J, tolerance = 1e-10)
  expect_equal(res$info_used, "etasq")
})

# P1 test (6): vectorization is preserved -----------------------------------
test_that("P1: vectorized es_from_etasq_adj matches row-wise scalar calls", {
  etasq_adj <- c(0.05, 0.10, 0.20)
  n_exp <- c(20, 15, 30)
  n_nexp <- c(20, 25, 28)
  n_cov <- c(1, 3, 2)
  r <- c(0.3, 0.5, 0.7)
  rev <- c(FALSE, TRUE, NA) # NA must be treated as FALSE

  res_vec <- es_from_etasq_adj(
    etasq_adj = etasq_adj, n_exp = n_exp, n_nexp = n_nexp,
    n_cov_ancova = n_cov, cov_outcome_r = r, reverse_etasq = rev
  )

  for (i in seq_along(etasq_adj)) {
    res_i <- es_from_etasq_adj(
      etasq_adj = etasq_adj[i], n_exp = n_exp[i], n_nexp = n_nexp[i],
      n_cov_ancova = n_cov[i], cov_outcome_r = r[i],
      reverse_etasq = isTRUE(rev[i])
    )
    for (col in c("d", "d_se", "d_ci_lo", "d_ci_up", "g", "g_se", "g_ci_lo", "g_ci_up")) {
      expect_equal(res_vec[[col]][i], res_i[[col]],
        tolerance = 1e-10,
        label = sprintf("row %d, column %s", i, col)
      )
    }
  }
})
