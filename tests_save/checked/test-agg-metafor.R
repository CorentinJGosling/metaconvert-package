# Comparison of aggregate_df() against metafor::aggregate.escalc() as the oracle.
#
# Mapping of metaConvert dependence types to metafor V-matrix structures:
#   subgroups <-> struct = "ID"  (independent, inverse-variance weighted)
#   outcomes  <-> struct = "CS"  (compound symmetry, unweighted, rho = cor_unit)
#   times     <-> struct = "CAR" (continuous-time AR1, weighted GLS,
#                                 phi = cor_unit, time = time_agg)
#
# These tests deliberately exercise the input shapes that the previous
# implementation mishandled (unsorted / string / factor cluster ids, and NA
# es/se) and that the original test-agg.R never covered (it used only
# pre-sorted single-column integer ids with no NA).

skip_if_not_installed("metafor")

## ---- helpers ---------------------------------------------------------------

# Run metafor::aggregate.escalc on a dependent-ES data.frame and return a named
# list keyed by cluster id (character), so comparisons are order-independent.
metafor_agg <- function(dat, struct, cor) {
  e <- metafor::escalc(yi = dat$es, sei = dat$se)
  e$agg      <- dat$agg
  e$time_agg <- dat$time_agg
  cl <- as.character(dat$agg)

  a <- switch(struct,
    "ID"  = metafor::aggregate.escalc(e, cluster = cl, struct = "ID",
                                      weighted = TRUE),
    "CS"  = metafor::aggregate.escalc(e, cluster = cl, struct = "CS",
                                      weighted = FALSE, rho = cor),
    "CAR" = metafor::aggregate.escalc(e, cluster = cl, struct = "CAR",
                                      phi = cor, time = e$time_agg)
  )
  key <- as.character(a$agg)
  list(yi = setNames(as.numeric(a$yi), key),
       vi = setNames(as.numeric(a$vi), key))
}

# Run aggregate_df and return the same keyed representation.
mc_agg <- function(dat, dependence, cor) {
  r <- aggregate_df(dat, dependence = dependence, agg_fact = "agg",
                    es = "es", se = "se", cor_unit = cor)
  key <- as.character(r$agg)
  list(yi = setNames(r$es, key),
       vi = setNames(r$se^2, key))
}

# Assert metaConvert == metafor for both the point estimate and the variance,
# aligning by cluster id.
expect_matches_metafor <- function(dat, struct, dependence, cor = 0.5,
                                   tol = 1e-8) {
  mf <- metafor_agg(dat, struct, cor)
  mc <- mc_agg(dat, dependence, cor)
  keys <- sort(names(mf$yi))
  expect_setequal(names(mc$yi), keys)
  expect_equal(mc$yi[keys], mf$yi[keys], tolerance = tol)
  expect_equal(mc$vi[keys], mf$vi[keys], tolerance = tol)
}

# A dependent-ES dataset with clusters of sizes 1, 3, 2, 2 (mixed) plus a
# per-row time index. No cor_unit column: the "times" path then falls back to
# the cor_unit argument (exercised below), matching metafor's phi. The
# column-precedence path is checked in its own test.
make_dat <- function(ids, seed = 123) {
  set.seed(seed)
  n <- length(ids)
  data.frame(
    agg      = ids,
    time_agg = ave(seq_len(n), as.character(ids), FUN = seq_along),
    es       = rnorm(n),
    se       = abs(rnorm(n)) + 0.3,
    stringsAsFactors = FALSE
  )
}

# cluster layout: sizes 1,3,2,2 -> 8 rows
LAYOUT <- c("k1", "k2", "k2", "k2", "k3", "k3", "k4", "k4")

## ---- struct = ID (subgroups) ----------------------------------------------

test_that("subgroups == metafor ID: sorted integer ids", {
  dat <- make_dat(rep(1:4, c(1, 3, 2, 2)))
  expect_matches_metafor(dat, "ID", "subgroups")
})

test_that("subgroups == metafor ID: UNSORTED integer ids", {
  dat <- make_dat(c(4, 2, 2, 2, 1, 1, 3, 3))
  expect_matches_metafor(dat, "ID", "subgroups")
})

test_that("subgroups == metafor ID: UNSORTED string ids", {
  dat <- make_dat(c("Zeta", "Alpha", "Alpha", "Alpha", "Mu", "Mu", "Beta", "Beta"))
  expect_matches_metafor(dat, "ID", "subgroups")
})

test_that("subgroups == metafor ID: FACTOR ids with unused levels", {
  dat <- make_dat(factor(rep(c("a", "b", "c", "d"), c(1, 3, 2, 2)),
                         levels = c("a", "b", "c", "d", "e", "f")))
  expect_matches_metafor(dat, "ID", "subgroups")
})

test_that("subgroups == metafor ID: NA es/se rows are dropped like metafor", {
  dat <- make_dat(rep(1:4, c(1, 3, 2, 2)))
  dat$es[3] <- NA   # inside cluster 2
  dat$se[7] <- NA   # inside cluster 4
  expect_matches_metafor(dat, "ID", "subgroups")
})

## ---- struct = CS (outcomes) -----------------------------------------------

test_that("outcomes == metafor CS: sorted integer ids", {
  dat <- make_dat(rep(1:4, c(1, 3, 2, 2)))
  expect_matches_metafor(dat, "CS", "outcomes", cor = 0.6)
})

test_that("outcomes == metafor CS: UNSORTED string ids", {
  dat <- make_dat(c("Zeta", "Alpha", "Alpha", "Alpha", "Mu", "Mu", "Beta", "Beta"))
  expect_matches_metafor(dat, "CS", "outcomes", cor = 0.6)
})

test_that("outcomes == metafor CS: FACTOR ids with unused levels", {
  dat <- make_dat(factor(rep(c("a", "b", "c", "d"), c(1, 3, 2, 2)),
                         levels = letters[1:6]))
  expect_matches_metafor(dat, "CS", "outcomes", cor = 0.3)
})

test_that("outcomes == metafor CS: NA es/se rows dropped like metafor", {
  dat <- make_dat(rep(1:4, c(1, 3, 2, 2)))
  dat$se[2] <- NA   # inside cluster 2
  expect_matches_metafor(dat, "CS", "outcomes", cor = 0.6)
})

## ---- struct = CAR (times) --------------------------------------------------

test_that("times == metafor CAR: sorted integer ids", {
  dat <- make_dat(rep(1:4, c(1, 3, 2, 2)))
  expect_matches_metafor(dat, "CAR", "times", cor = 0.5)
})

test_that("times == metafor CAR: UNSORTED string ids", {
  dat <- make_dat(c("Zeta", "Alpha", "Alpha", "Alpha", "Mu", "Mu", "Beta", "Beta"))
  expect_matches_metafor(dat, "CAR", "times", cor = 0.7)
})

test_that("times == metafor CAR: FACTOR ids with unused levels", {
  dat <- make_dat(factor(rep(c("a", "b", "c", "d"), c(1, 3, 2, 2)),
                         levels = letters[1:6]))
  expect_matches_metafor(dat, "CAR", "times", cor = 0.4)
})

test_that("times: cor_unit COLUMN takes precedence over the argument", {
  dat <- make_dat(rep(1:4, c(1, 3, 2, 2)))
  dat$cor_unit <- 0.55                       # per-row column, constant per cluster
  # metafor oracle built with phi = 0.55 (the column value)
  mf <- metafor_agg(dat, "CAR", cor = 0.55)
  # aggregate_df called with a DIFFERENT argument (0.1): the column must win
  r  <- aggregate_df(dat, dependence = "times", agg_fact = "agg",
                     es = "es", se = "se", cor_unit = 0.1)
  key <- as.character(r$agg)
  expect_equal(setNames(r$es, key)[names(mf$yi)], mf$yi, tolerance = 1e-8)
  expect_equal(setNames(r$se^2, key)[names(mf$vi)], mf$vi, tolerance = 1e-8)
})
