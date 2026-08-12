## Reversing a tetrachoric correlation must REFLECT its confidence interval, not
## negate the bounds in place.
##
## The tetrachoric routes store (r, vr, r_ci_lo, r_ci_up, z, vz, z_ci_lo, z_ci_up) and
## used to apply `-` to each bound individually when `reverse_* = TRUE`. Negating
## without swapping leaves lo > up, so the interval no longer brackets its own point
## estimate:
##
##   es_from_2x2(143, 52, 41, 164)                      r =  0.7458746  CI ( 0.6591207,  0.8326284)
##   es_from_2x2(143, 52, 41, 164, reverse_2x2 = TRUE)   r = -0.7458746  CI (-0.6591207, -0.8326284)  <- inverted
##   correct                                                             CI (-0.8326284, -0.6591207)
##
## The same defect was found and fixed earlier on the OR path (.or_to_cor result block
## in es_from_stand_OR.R). Both tetrachoric intervals are symmetric Wald bounds around
## their own estimate, so reflecting is exactly equivalent to recomputing the interval
## around the negated estimate -- which is what these tests assert.

tetra_grid <- function() {
  g <- expand.grid(
    a = c(3, 12, 40, 143, 220),
    b = c(2, 9, 52, 130),
    c = c(1, 15, 41, 160),
    d = c(4, 20, 164, 300)
  )
  g[, c("a", "b", "c", "d")]
}

test_that("es_from_2x2 tetrachoric: the documented reference table reflects exactly", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  fwd <- es_from_2x2(143, 52, 41, 164, table_2x2_to_cor = "tetrachoric",
                     reverse_2x2 = FALSE)
  rev <- es_from_2x2(143, 52, 41, 164, table_2x2_to_cor = "tetrachoric",
                     reverse_2x2 = TRUE)

  # Forward values are unchanged by this fix (regression anchors, 8+ significant digits)
  expect_equal(fwd$r, 0.74587459, tolerance = 1e-7)
  expect_equal(fwd$r_ci_lo, 0.65912073, tolerance = 1e-7)
  expect_equal(fwd$r_ci_up, 0.83262844, tolerance = 1e-7)

  # Reversed: (lo, up) -> (-up, -lo)
  expect_equal(rev$r, -0.74587459, tolerance = 1e-7)
  expect_equal(rev$r_ci_lo, -0.83262844, tolerance = 1e-7)
  expect_equal(rev$r_ci_up, -0.65912073, tolerance = 1e-7)
  expect_equal(rev$z_ci_lo, -1.15912791, tolerance = 1e-7)
  expect_equal(rev$z_ci_up, -0.76805510, tolerance = 1e-7)

  # And the interval is the right way round
  expect_lt(rev$r_ci_lo, rev$r_ci_up)
  expect_lt(rev$z_ci_lo, rev$z_ci_up)
})

test_that("es_from_2x2 tetrachoric: CI brackets the estimate for both reverse settings", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  g <- tetra_grid()
  fwd <- suppressWarnings(es_from_2x2(g$a, g$b, g$c, g$d,
                                      table_2x2_to_cor = "tetrachoric",
                                      reverse_2x2 = rep(FALSE, nrow(g))))
  rev <- suppressWarnings(es_from_2x2(g$a, g$b, g$c, g$d,
                                      table_2x2_to_cor = "tetrachoric",
                                      reverse_2x2 = rep(TRUE, nrow(g))))
  ok <- !is.na(fwd$r) & !is.na(fwd$r_ci_lo) & !is.na(fwd$r_ci_up) &
        !is.na(rev$r) & !is.na(rev$r_ci_lo) & !is.na(rev$r_ci_up) &
        !is.na(fwd$z_ci_lo) & !is.na(rev$z_ci_lo)
  expect_gt(sum(ok), 200)  # the grid must actually exercise the path

  expect_true(all(fwd$r_ci_lo[ok] <= fwd$r[ok] & fwd$r[ok] <= fwd$r_ci_up[ok]))
  expect_true(all(rev$r_ci_lo[ok] <= rev$r[ok] & rev$r[ok] <= rev$r_ci_up[ok]))
  expect_true(all(fwd$z_ci_lo[ok] <= fwd$z[ok] & fwd$z[ok] <= fwd$z_ci_up[ok]))
  expect_true(all(rev$z_ci_lo[ok] <= rev$z[ok] & rev$z[ok] <= rev$z_ci_up[ok]))
})

test_that("es_from_2x2 tetrachoric: reversing maps (lo, up) -> (-up, -lo) exactly", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  g <- tetra_grid()
  fwd <- suppressWarnings(es_from_2x2(g$a, g$b, g$c, g$d,
                                      table_2x2_to_cor = "tetrachoric",
                                      reverse_2x2 = rep(FALSE, nrow(g))))
  rev <- suppressWarnings(es_from_2x2(g$a, g$b, g$c, g$d,
                                      table_2x2_to_cor = "tetrachoric",
                                      reverse_2x2 = rep(TRUE, nrow(g))))

  expect_equal(rev$r,       -fwd$r,       tolerance = 0)
  expect_equal(rev$r_ci_lo, -fwd$r_ci_up, tolerance = 0)
  expect_equal(rev$r_ci_up, -fwd$r_ci_lo, tolerance = 0)
  expect_equal(rev$z,       -fwd$z,       tolerance = 0)
  expect_equal(rev$z_ci_lo, -fwd$z_ci_up, tolerance = 0)
  expect_equal(rev$z_ci_up, -fwd$z_ci_lo, tolerance = 0)

  # SEs are direction-invariant
  expect_equal(rev$r_se, fwd$r_se, tolerance = 0)
  expect_equal(rev$z_se, fwd$z_se, tolerance = 0)
})

test_that("es_from_2x2 tetrachoric: reflected bounds equal the Wald interval rebuilt around -est", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  # An independent check that does not reference the forward interval at all: both
  # tetrachoric intervals are est +- qnorm(.975) * se, so the reversed row must satisfy
  # that identity about its own (negated) estimate.
  g <- tetra_grid()
  rev <- suppressWarnings(es_from_2x2(g$a, g$b, g$c, g$d,
                                      table_2x2_to_cor = "tetrachoric",
                                      reverse_2x2 = rep(TRUE, nrow(g))))
  ok <- !is.na(rev$r) & !is.na(rev$r_se) & !is.na(rev$z) & !is.na(rev$z_se)
  z975 <- stats::qnorm(.975)
  expect_equal(rev$r_ci_lo[ok], rev$r[ok] - z975 * rev$r_se[ok], tolerance = 1e-12)
  expect_equal(rev$r_ci_up[ok], rev$r[ok] + z975 * rev$r_se[ok], tolerance = 1e-12)
  expect_equal(rev$z_ci_lo[ok], rev$z[ok] - z975 * rev$z_se[ok], tolerance = 1e-12)
  expect_equal(rev$z_ci_up[ok], rev$z[ok] + z975 * rev$z_se[ok], tolerance = 1e-12)
})

test_that("es_from_2x2_sum tetrachoric also reflects (wrapper reaches the same route)", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  # 143/52/41/164 expressed as arm totals
  fwd <- es_from_2x2_sum(n_cases_exp = 143, n_exp = 195,
                         n_cases_nexp = 41, n_nexp = 205,
                         table_2x2_to_cor = "tetrachoric", reverse_2x2 = FALSE)
  rev <- es_from_2x2_sum(n_cases_exp = 143, n_exp = 195,
                         n_cases_nexp = 41, n_nexp = 205,
                         table_2x2_to_cor = "tetrachoric", reverse_2x2 = TRUE)
  expect_equal(rev$r_ci_lo, -fwd$r_ci_up, tolerance = 0)
  expect_equal(rev$r_ci_up, -fwd$r_ci_lo, tolerance = 0)
  expect_lt(rev$r_ci_lo, rev$r_ci_up)
  expect_true(rev$r_ci_lo <= rev$r && rev$r <= rev$r_ci_up)
})

test_that(".phi_to_cor tetrachoric reflects on reverse_phi", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  # .phi_to_cor / .chi_to_cor are currently unreachable from any exported function
  # (their call sites in es_from_PHI_CHISQ.R are commented out), so they are exercised
  # directly. They carried the identical negate-in-place defect and are fixed with it,
  # so that reviving those routes does not re-introduce an inverted CI.
  f <- metaConvert:::.phi_to_cor(
    phi = .5, n_sample = 400,
    n_cases_exp = 143, n_controls_exp = 52,
    n_cases_nexp = 41, n_controls_nexp = 164,
    phi_to_cor = "tetrachoric", reverse_phi = FALSE)
  r <- metaConvert:::.phi_to_cor(
    phi = .5, n_sample = 400,
    n_cases_exp = 143, n_controls_exp = 52,
    n_cases_nexp = 41, n_controls_nexp = 164,
    phi_to_cor = "tetrachoric", reverse_phi = TRUE)

  expect_equal(as.numeric(r[1]), -as.numeric(f[1]), tolerance = 0)
  expect_equal(as.numeric(r[3]), -as.numeric(f[4]), tolerance = 0)  # r_ci_lo <- -r_ci_up
  expect_equal(as.numeric(r[4]), -as.numeric(f[3]), tolerance = 0)
  expect_equal(as.numeric(r[7]), -as.numeric(f[8]), tolerance = 0)  # z_ci_lo <- -z_ci_up
  expect_equal(as.numeric(r[8]), -as.numeric(f[7]), tolerance = 0)
  expect_lt(as.numeric(r[3]), as.numeric(r[4]))
  expect_lt(as.numeric(r[7]), as.numeric(r[8]))
  expect_true(as.numeric(r[3]) <= as.numeric(r[1]) &&
              as.numeric(r[1]) <= as.numeric(r[4]))
  # vr / vz untouched by the reversal
  expect_equal(as.numeric(r[2]), as.numeric(f[2]), tolerance = 0)
  expect_equal(as.numeric(r[6]), as.numeric(f[6]), tolerance = 0)
})

test_that(".chi_to_cor tetrachoric reflects on reverse_chisq", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  f <- metaConvert:::.chi_to_cor(
    chisq = 100, n_sample = 400,
    n_cases_exp = 143, n_controls_exp = 52,
    n_cases_nexp = 41, n_controls_nexp = 164,
    chisq_to_cor = "tetrachoric", reverse_chisq = FALSE)
  r <- metaConvert:::.chi_to_cor(
    chisq = 100, n_sample = 400,
    n_cases_exp = 143, n_controls_exp = 52,
    n_cases_nexp = 41, n_controls_nexp = 164,
    chisq_to_cor = "tetrachoric", reverse_chisq = TRUE)

  expect_equal(as.numeric(r[1]), -as.numeric(f[1]), tolerance = 0)
  expect_equal(as.numeric(r[3]), -as.numeric(f[4]), tolerance = 0)
  expect_equal(as.numeric(r[4]), -as.numeric(f[3]), tolerance = 0)
  expect_equal(as.numeric(r[7]), -as.numeric(f[8]), tolerance = 0)
  expect_equal(as.numeric(r[8]), -as.numeric(f[7]), tolerance = 0)
  expect_lt(as.numeric(r[3]), as.numeric(r[4]))
  expect_lt(as.numeric(r[7]), as.numeric(r[8]))
  expect_true(as.numeric(r[3]) <= as.numeric(r[1]) &&
              as.numeric(r[1]) <= as.numeric(r[4]))
})

test_that("summary() no longer flags a reversed tetrachoric row as an inverted CI", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  dat <- data.frame(
    study_id = c("s1", "s2"),
    n_cases_exp = c(143, 120), n_controls_exp = c(52, 60),
    n_cases_nexp = c(41, 55), n_controls_nexp = c(164, 150),
    reverse_2x2 = c(TRUE, TRUE)
  )
  s <- summary(convert_df(dat, measure = "r", table_2x2_to_cor = "tetrachoric",
                          verbose = FALSE),
               flags = TRUE)
  flags <- paste(s$es_flags, collapse = " | ")
  expect_false(grepl("Inverted CI", flags, fixed = TRUE))
  expect_false(grepl("outside", flags, fixed = TRUE))
  expect_true(all(s$es_ci_lo <= s$es & s$es <= s$es_ci_up))
})
