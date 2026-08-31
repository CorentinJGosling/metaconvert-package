# ---------------------------------------------------------------------------
# Fisher's-z saturation on the smd_to_cor = "viechtbauer" route.
#
# NOT in AUDIT-2026-08-28-findings.md -- found on shipped data (df.short).
#
# .smd_to_cor() (R/internal_multiple_formulas.R, viechtbauer branch) computes the
# biserial r as sqrt(p*q)/f * r_pb, which is NOT bounded by 1: the factor
# sqrt(p*q)/f is 1.2533 at p = 0.5, so any point-biserial r_pb above 0.798 pushes
# it out of range. The branch then builds a CLAMPED copy, r_trunc, and uses it for
# the variance AND for the variance-stabilising z, while returning the UNCLAMPED
# r as the point estimate. Consequences:
#
#   * measure = "r"  -> r > 1 is returned and Category B fires
#                       "[INVALID] r outside [-1, 1]".  (working; pinned below)
#   * measure = "z"  -> z collapses onto the r_trunc = 1 saturation constant
#                       z_sat(p) = (a/2) * log((1+a)/(1-a)), a data-independent
#                       function of the arm-size split alone, with a plausible SE,
#                       a plausible CI, and NO flag of any kind: Category B has no
#                       z-scale counterpart, and z is unbounded so no bound check
#                       can catch it.
#
# On df.short this yields z = 1.2840605 for Lopez_2019 (row 6, route cohen_d) and
# z = 1.2842612 for Yu_2018 (row 28, route med_min_max) -- two different studies
# from two different routes pinned to the same constant. Both equal z_sat(p) for
# their own p to 7 decimals, so the number carries no information about the data.
#
# Correct behaviour: a z built from an out-of-range biserial r must be NA, or the
# row must carry an [INVALID] flag. It must never be silently reported.
# ---------------------------------------------------------------------------

# Saturation constant: the value the viechtbauer z takes once r is clamped to 1.
# Derived independently from the branch's own algebra (a is the variance-stabilising
# constant, z = (a/2) * log((1 + a*r)/(1 - a*r)) evaluated at r = 1).
.z_saturation <- function(p) {
  a <- sqrt(stats::dnorm(stats::qnorm(p))) / (p * (1 - p))^(1 / 4)
  (a / 2) * log((1 + a) / (1 - a))
}

# Independent closed form for the viechtbauer biserial r and its z, from
# Viechtbauer's r_pb -> biserial rescaling. Used as the "still correct" anchor.
.viechtbauer_ref <- function(d, n_exp, n_nexp) {
  df <- n_exp + n_nexp - 2
  h <- df / n_exp + df / n_nexp
  p <- n_exp / (n_exp + n_nexp)
  q <- 1 - p
  r_pb <- d / sqrt(d^2 + h)
  f <- stats::dnorm(stats::qnorm(p, lower.tail = FALSE))
  r <- sqrt(p * q) / f * r_pb
  a <- sqrt(stats::dnorm(stats::qnorm(p))) / (p * q)^(1 / 4)
  list(r = r, z = (a / 2) * log((1 + a * r) / (1 - a * r)))
}

.quiet <- function(expr) {
  val <- NULL
  invisible(utils::capture.output(
    val <- suppressWarnings(suppressMessages(force(expr)))
  ))
  val
}

# Two rows, same arm sizes, materially different effects. d = 3 gives the biserial
# r = 1.0482, d = 4 gives r = 1.1248 -- both out of range, both must not be
# reported as the same z.
.sat_df <- data.frame(
  study_id = c("S_d3", "S_d4"),
  n_exp    = c(30, 30),
  n_nexp   = c(30, 30),
  cohen_d  = c(3, 4),
  stringsAsFactors = FALSE
)


# ---------------------------------------------------------------------------
# CONTROL (expected to PASS today). The r-scale half of the bound check already
# works; this pins it so a fix to the z scale cannot be made by weakening it, and
# pins the valid-r arithmetic so a fix cannot be made by NA-ing everything.
# ---------------------------------------------------------------------------
test_that("AUDIT-z-saturation reference: measure='r' already flags the out-of-range biserial r", {
  s_r <- .quiet(summary(convert_df(.sat_df, measure = "r"), flags = TRUE))

  # Both rows return an impossible correlation ...
  expect_true(all(s_r$es_crude > 1))
  # ... and Category B says so, on both rows.
  expect_true(all(grepl("r outside [-1, 1]", s_r$flags_crude, fixed = TRUE)))

  # A legitimate d (biserial r well inside the parameter space) must keep its exact
  # Viechtbauer value -- an independently derived closed form, not a captured number.
  ok <- es_from_cohen_d(cohen_d = 1, n_exp = 30, n_nexp = 30,
                        smd_to_cor = "viechtbauer")
  ref <- .viechtbauer_ref(1, 30, 30)
  expect_equal(ok$r, ref$r, tolerance = 1e-9)   # 0.5681253106
  expect_equal(ok$z, ref$z, tolerance = 1e-9)   # 0.4996112009
})


# ---------------------------------------------------------------------------
# DEFECT (must FAIL today).
# On the z scale a saturated row is returned with no flag and no NA: the bound
# check that fires on the r scale has no z-scale counterpart.
# ---------------------------------------------------------------------------
test_that("AUDIT-z-saturation: a z built from an out-of-range biserial r must be NA or [INVALID]", {
  s_z <- .quiet(summary(convert_df(.sat_df, measure = "z"), flags = TRUE))

  for (i in seq_len(nrow(s_z))) {
    handled <- is.na(s_z$es_crude[i]) ||
      grepl("[INVALID]", s_z$flags_crude[i], fixed = TRUE)
    expect_true(
      handled,
      info = sprintf(
        "row %s: es = %s, flags = '%s'",
        s_z$study_id[i], format(s_z$es_crude[i]), s_z$flags_crude[i]
      )
    )
  }

  # The same two rows on df.short: Lopez_2019 and Yu_2018 both saturate on the shipped
  # example data (before the fix, via routes cohen_d and med_min_max respectively).
  #
  # Locate them by study_id ONLY, never by info_used. A correct fix makes the saturated
  # estimate unavailable, so the row either drops out of the z pool or falls through to
  # another route -- keying the scaffold on the pre-fix route would make this test
  # unsatisfiable together with its own subject, and would silently stop checking the
  # rows it exists to check.
  s_short <- .quiet(summary(convert_df(df.short, measure = "z"), flags = TRUE))
  hit <- s_short$study_id %in% c("Lopez_2019", "Yu_2018")
  hit[is.na(hit)] <- FALSE
  expect_gte(sum(hit), 2L)   # guards the row selection itself

  # These two study_ids also carry perfectly healthy rows from other routes (means_sd),
  # so the requirement is NOT "every row of this study is NA" -- it is that no row of it
  # reports the saturation constant. A row that came through a route whose biserial r
  # stayed in range is a correct estimate and must be left alone.
  for (i in which(hit)) {
    p <- s_short$n_exp[i] / (s_short$n_exp[i] + s_short$n_nexp[i])
    reports_sat <- !is.na(s_short$es_crude[i]) && !is.na(p) &&
      abs(s_short$es_crude[i] - .z_saturation(p)) < 1e-6
    handled <- !reports_sat ||
      grepl("[INVALID]", s_short$flags_crude[i], fixed = TRUE)
    expect_true(
      handled,
      info = sprintf(
        "df.short %s / %s: es = %s is the saturation constant and carries no [INVALID]; flags = '%s'",
        s_short$study_id[i], s_short$info_used_crude[i],
        format(s_short$es_crude[i]), s_short$flags_crude[i]
      )
    )
  }
})


# ---------------------------------------------------------------------------
# DEFECT (must FAIL today).
# The reported z is the r = 1 saturation constant, i.e. a function of the arm-size
# split only. Two rows with materially different data therefore return the same z.
# ---------------------------------------------------------------------------
test_that("AUDIT-z-saturation: saturated rows must not collapse onto the r=1 constant", {
  # Route level: identical n, d = 3 vs d = 4.
  out <- es_from_cohen_d(cohen_d = c(3, 4), n_exp = c(30, 30), n_nexp = c(30, 30),
                         smd_to_cor = "viechtbauer")

  z_sat_half <- .z_saturation(0.5)          # 1.2842611620
  expect_equal(round(z_sat_half, 7), 1.2842612)   # the constant is what I claim

  # Neither row may report the saturation constant as its effect size.
  for (i in 1:2) {
    is_sat <- !is.na(out$z[i]) && abs(out$z[i] - z_sat_half) < 1e-6
    expect_false(
      is_sat,
      info = sprintf("d = %s returned the r=1 saturation constant z = %.7f",
                     out$d[i], out$z[i])
    )
  }

  # d = 3 and d = 4 are materially different data; they must not share a z.
  both_present <- !is.na(out$z[1]) && !is.na(out$z[2])
  expect_false(
    both_present && abs(out$z[1] - out$z[2]) < 1e-3,
    info = sprintf("d = 3 -> z = %.7f, d = 4 -> z = %.7f", out$z[1], out$z[2])
  )

  # Shipped data: Lopez_2019 (cohen_d, p = 31/61) and Yu_2018 (med_min_max, p = 1/2)
  # are different studies from different routes; each currently returns exactly the
  # saturation constant for its own arm-size split.
  s_short <- .quiet(summary(convert_df(df.short, measure = "z"), flags = TRUE))
  # By study_id only -- see the note in the previous block. A fixed package may report
  # these rows through a different route, or not at all; what must never happen is that
  # one of them reports the saturation constant as its effect size.
  pick <- function(sid) s_short[which(s_short$study_id == sid), , drop = FALSE]
  lopez <- pick("Lopez_2019")
  yu    <- pick("Yu_2018")
  expect_gt(nrow(lopez), 0L)
  expect_gt(nrow(yu), 0L)

  for (df_row in list(lopez, yu)) {
    for (i in seq_len(nrow(df_row))) {
      row <- df_row[i, ]
      p <- row$n_exp / (row$n_exp + row$n_nexp)
      is_sat <- !is.na(row$es_crude) && !is.na(p) &&
        abs(row$es_crude - .z_saturation(p)) < 1e-6
      expect_false(
        is_sat,
        info = sprintf("%s / %s: z = %s is exactly z_sat(p = %s)",
                       row$study_id, row$info_used_crude,
                       format(row$es_crude), format(p))
      )
    }
  }

  # ... and the two consequently sit on top of each other despite unrelated data.
  lz <- lopez$es_crude[!is.na(lopez$es_crude)]
  yz <- yu$es_crude[!is.na(yu$es_crude)]
  collide <- length(lz) > 0 && length(yz) > 0 && any(abs(outer(lz, yz, "-")) < 1e-3)
  expect_false(
    collide,
    info = sprintf("Lopez_2019 z = %s, Yu_2018 z = %s",
                   paste(format(lz), collapse = "/"), paste(format(yz), collapse = "/"))
  )
})
