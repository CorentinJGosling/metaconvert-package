# =============================================================================
# Regression coverage for two audit fixes in es_from_spearman_rho():
#
#   P3  (digest [20]): with cor_to_smd = "cooper", the d/g/logOR SEs must carry
#        the same Bonett-Wright Spearman delta correction as the r/z SEs of the
#        same output row (the cooper d_se is exactly linear in r_se, so the
#        r_se_delta / r_se_pearson scale factor applies to it just as it does
#        on the default viechtbauer path). mathur stays excluded: its variance
#        genuinely has no r_se term.
#
#   P15 (digest [24]): direct calls supplying n_exp/n_nexp but no n_sample must
#        fall back to n_sample = n_exp + n_nexp (as es_from_pearson_r does
#        internally), instead of returning r_se/z_se = NA next to UNCORRECTED
#        d/g/logOR SEs.
#
# Design rule: every numeric oracle below is hand-derived IN this file from
# published formulas — never re-typed from the R/ source:
#   Rupinski & Dunlap (1996)   r_p = 2*sin(pi/6 * r_s)
#   Bonett & Wright (2000)     var(r_s) = (1 + r_s^2/2) * (1 - r_s^2)^2 / (n - 1)
#   delta-method derivative    d r_p / d r_s = (pi/3) * cos(pi/6 * r_s)
#   Cooper handbook / Borenstein  d = 2 r / sqrt(1 - r^2),
#                                 d_se = sqrt(4 * v_r / (1 - r^2)^3)
#   Mathur & VanderWeele (2020)   d_se = |d| * sqrt(1/(r^2 (n-3)) + 1/(2 (n-1)))
#   Hasselblad & Hedges (1995)    logor_se = d_se * pi / sqrt(3)
#   Hedges (1981) exact J         J(df) = gamma(df/2) / (sqrt(df/2) gamma((df-1)/2))
# =============================================================================

# ---- hand-derived oracle helpers (published formulas, independent of R/) ----

# Rupinski & Dunlap (1996) Spearman -> Pearson transform
.orc_rp <- function(rs) 2 * sin(pi / 6 * rs)

# delta-method SE of the converted Pearson r_p, with the Bonett & Wright (2000)
# Spearman sampling variance: SE = |d r_p/d r_s| * sqrt(var(r_s))
.orc_r_se_delta <- function(rs, n) {
  deriv <- (pi / 3) * cos(pi / 6 * rs)
  sqrt(deriv^2 * (1 + rs^2 / 2) * (1 - rs^2)^2 / (n - 1))
}

# Hedges (1981) exact small-sample correction factor at df degrees of freedom
.orc_hedges_j <- function(df) exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))

# -----------------------------------------------------------------------------
# P3: cooper d/g/logOR SEs must be delta-consistent with the corrected r_se
# -----------------------------------------------------------------------------
test_that("P3: cooper d_se/g_se/logor_se carry the Spearman delta correction", {
  n <- 100
  for (rs in c(0.2, 0.55, 0.85)) {
    res <- es_from_spearman_rho(
      spearman_r = rs, n_sample = n, cor_to_smd = "cooper"
    )

    rp <- .orc_rp(rs)
    r_se_delta <- .orc_r_se_delta(rs, n)

    # Cooper handbook: v_d = 4 * v_r / (1 - r^2)^3, fed the DELTA-corrected v_r
    d_se_expected <- sqrt(4 * r_se_delta^2 / (1 - rp^2)^3)

    expect_equal(res$d_se, d_se_expected, tolerance = 1e-10)

    # downstream measures inherit the same correction:
    # Hasselblad & Hedges (1995): logor_se = d_se * pi / sqrt(3)
    expect_equal(res$logor_se, d_se_expected * pi / sqrt(3), tolerance = 1e-10)
    # Hedges (1981): g_se = J(n - 2) * d_se
    expect_equal(res$g_se, .orc_hedges_j(n - 2) * d_se_expected, tolerance = 1e-10)

    # and the r_se in the SAME row is the delta-corrected SE (no contradiction)
    expect_equal(res$r_se, r_se_delta, tolerance = 1e-10)
  }
})

test_that("P3: within one cooper row, d_se/r_se matches es_from_pearson_r's cooper ratio", {
  # es_from_pearson_r is a different function, independently verified: for
  # cooper, its d_se is linear in its r_se, so the ratio d_se/r_se depends only
  # on r_p — the same ratio must hold inside the Spearman row once the delta
  # scale factor is applied to BOTH numerator and denominator.
  n <- 100
  for (rs in c(0.2, 0.55, 0.85)) {
    rp <- .orc_rp(rs)

    res_sp <- es_from_spearman_rho(
      spearman_r = rs, n_sample = n, cor_to_smd = "cooper"
    )
    res_pe <- es_from_pearson_r(
      pearson_r = rp, n_sample = n, cor_to_smd = "cooper"
    )

    expect_equal(res_sp$d_se / res_sp$r_se,
                 res_pe$d_se / res_pe$r_se,
                 tolerance = 1e-10)
    # the point estimates agree exactly (same converted r_p, same cooper d)
    expect_equal(res_sp$d, res_pe$d, tolerance = 1e-10)
  }
})

# -----------------------------------------------------------------------------
# mathur must stay excluded from the scale factor (its variance has no r_se term)
# -----------------------------------------------------------------------------
test_that("mathur SEs are unchanged by the cooper fix (published Mathur formula)", {
  n <- 100
  rs <- 0.55
  sd_iv <- 1.5
  unit_increase <- 1

  res <- es_from_spearman_rho(
    spearman_r = rs, n_sample = n, cor_to_smd = "mathur",
    sd_iv = sd_iv, unit_increase_iv = unit_increase, unit_type = "raw_scale"
  )

  # Mathur & VanderWeele: d = r * increase / (sd_iv * sqrt(1 - r^2)),
  # d_se = |d| * sqrt(1 / (r^2 (n - 3)) + 1 / (2 (n - 1))), applied to the
  # CONVERTED Pearson r_p; no r_se enters, so no delta scaling is expected.
  rp <- .orc_rp(rs)
  d_expected <- rp * unit_increase / (sd_iv * sqrt(1 - rp^2))
  d_se_expected <- abs(d_expected) *
    sqrt(1 / (rp^2 * (n - 3)) + 1 / (2 * (n - 1)))

  expect_equal(res$d, d_expected, tolerance = 1e-10)
  expect_equal(res$d_se, d_se_expected, tolerance = 1e-10)
  expect_equal(res$g_se, .orc_hedges_j(n - 2) * d_se_expected, tolerance = 1e-10)

  # pin the current numeric values as an extra regression anchor
  expect_equal(res$d, 0.4601260434951955, tolerance = 1e-10)
  expect_equal(res$d_se, 0.0885088301895168, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# default viechtbauer path: regression pins (verified-correct pre-fix values)
# -----------------------------------------------------------------------------
test_that("default viechtbauer path is unchanged (pins + metafor oracle)", {
  n <- 100
  rs <- 0.55
  res <- es_from_spearman_rho(spearman_r = rs, n_sample = n)

  # literal pins: current values, verified bit-exact against the
  # metafor::conv.delta oracle (audit digest [26])
  expect_equal(res$r_se, 0.0755226148886572, tolerance = 1e-10)
  expect_equal(res$z_se, 0.1114986391448512, tolerance = 1e-10)
  expect_equal(res$d_se, 0.1714666763684312, tolerance = 1e-10)
  expect_equal(res$g_se, 0.1701505131503400, tolerance = 1e-10)

  # in-test external oracle: Bonett-Wright delta variance fed to metafor's
  # r -> d delta transform; d_se recovered with the exact Hedges J at df = n - 2
  skip_if_not_installed("metafor")
  rp <- .orc_rp(rs)
  r_se_delta <- .orc_r_se_delta(rs, n)
  orc <- metafor::conv.delta(yi = rp, vi = r_se_delta^2,
                             transf = metafor::transf.rtod)
  J <- .orc_hedges_j(n - 2)
  expect_equal(res$r_se, r_se_delta, tolerance = 1e-10)
  expect_equal(res$z_se, r_se_delta / (1 - rp^2), tolerance = 1e-10)
  expect_equal(res$g_se, sqrt(orc$vi), tolerance = 1e-10)
  expect_equal(res$d_se, sqrt(orc$vi) / J, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# P15: n_exp/n_nexp without n_sample must behave exactly like n_sample = n1 + n2
# -----------------------------------------------------------------------------
test_that("P15: n_exp/n_nexp without n_sample falls back to n_exp + n_nexp (viechtbauer)", {
  res_arms <- es_from_spearman_rho(spearman_r = 0.55, n_exp = 50, n_nexp = 50)
  res_full <- es_from_spearman_rho(spearman_r = 0.55, n_sample = 100)

  expect_equal(res_arms$r_se, res_full$r_se, tolerance = 1e-10)
  expect_equal(res_arms$z_se, res_full$z_se, tolerance = 1e-10)
  # d / g and their SEs are NOT expected to coincide: supplied arm sizes are handed to
  # transf.rtod (n1i / n2i), which inverts the biserial at h = m/n1 + m/n2 rather than
  # at the balanced-arms constant 4, so the arms call is pinned to that oracle and the
  # n_sample-only call to the constant-4 one.
  rp_arms <- .orc_rp(0.55)
  expect_equal(res_arms$g, metafor::transf.rtod(rp_arms, 50, 50), tolerance = 1e-10)
  expect_equal(res_full$g, metafor::transf.rtod(rp_arms), tolerance = 1e-10)
  expect_equal(res_arms$g_se,
               sqrt(metafor::conv.delta(yi = rp_arms, vi = .orc_r_se_delta(0.55, 100)^2,
                                        transf = metafor::transf.rtod, n1i = 50, n2i = 50,
                                        var.names = c("g", "gv"))$gv),
               tolerance = 1e-10)
  expect_equal(res_arms$r_ci_lo, res_full$r_ci_lo, tolerance = 1e-10)
  expect_equal(res_arms$r_ci_up, res_full$r_ci_up, tolerance = 1e-10)
  expect_equal(res_arms$z_ci_lo, res_full$z_ci_lo, tolerance = 1e-10)
  expect_equal(res_arms$z_ci_up, res_full$z_ci_up, tolerance = 1e-10)

  # against the hand-derived oracle directly (not just internal equality)
  expect_equal(res_arms$r_se, .orc_r_se_delta(0.55, 100), tolerance = 1e-10)
})

test_that("P15 composes with P3: cooper path also corrected under the n fallback", {
  res_arms <- es_from_spearman_rho(spearman_r = 0.55, n_exp = 50, n_nexp = 50,
                                   cor_to_smd = "cooper")

  rp <- .orc_rp(0.55)
  r_se_delta <- .orc_r_se_delta(0.55, 100)
  d_se_expected <- sqrt(4 * r_se_delta^2 / (1 - rp^2)^3)

  expect_equal(res_arms$r_se, r_se_delta, tolerance = 1e-10)
  expect_equal(res_arms$d_se, d_se_expected, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# reverse_spearman_r still flips coherently after both fixes
# -----------------------------------------------------------------------------
test_that("reverse_spearman_r flips point estimates and preserves SEs (both paths)", {
  for (method in c("viechtbauer", "cooper")) {
    res_pos <- es_from_spearman_rho(spearman_r = 0.55, n_sample = 100,
                                    cor_to_smd = method)
    res_rev <- es_from_spearman_rho(spearman_r = 0.55, n_sample = 100,
                                    cor_to_smd = method,
                                    reverse_spearman_r = TRUE)
    res_neg <- es_from_spearman_rho(spearman_r = -0.55, n_sample = 100,
                                    cor_to_smd = method)

    # reversal == feeding the negated Spearman rho, on every column
    expect_equal(res_rev$r, res_neg$r, tolerance = 1e-10)
    expect_equal(res_rev$d, res_neg$d, tolerance = 1e-10)
    expect_equal(res_rev$z, res_neg$z, tolerance = 1e-10)
    expect_equal(res_rev$r_ci_lo, res_neg$r_ci_lo, tolerance = 1e-10)
    expect_equal(res_rev$r_ci_up, res_neg$r_ci_up, tolerance = 1e-10)

    # antisymmetry of the point estimates, symmetry of the SEs
    expect_equal(res_rev$r, -res_pos$r, tolerance = 1e-10)
    expect_equal(res_rev$d, -res_pos$d, tolerance = 1e-10)
    expect_equal(res_rev$z, -res_pos$z, tolerance = 1e-10)
    expect_equal(res_rev$r_se, res_pos$r_se, tolerance = 1e-10)
    expect_equal(res_rev$d_se, res_pos$d_se, tolerance = 1e-10)
    expect_equal(res_rev$g_se, res_pos$g_se, tolerance = 1e-10)
    expect_equal(res_rev$z_se, res_pos$z_se, tolerance = 1e-10)
  }

  # and the reversed cooper SE is still the delta-consistent one (P3 holds
  # under reversal: the correction is even in r_s)
  res_rev_c <- es_from_spearman_rho(spearman_r = 0.55, n_sample = 100,
                                    cor_to_smd = "cooper",
                                    reverse_spearman_r = TRUE)
  rp <- .orc_rp(0.55)
  r_se_delta <- .orc_r_se_delta(0.55, 100)
  expect_equal(res_rev_c$d_se, sqrt(4 * r_se_delta^2 / (1 - rp^2)^3),
               tolerance = 1e-10)
})
