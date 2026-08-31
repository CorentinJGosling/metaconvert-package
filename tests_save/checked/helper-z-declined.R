# Helpers for the measure = "z" pins after the out-of-range biserial r fix (2.1.0).
#
# Two df.haza rows carry a Viechtbauer biserial correlation outside [-1, 1] (1.0527 and
# 1.0702). sqrt(p*q)/f is 1.2533 at p = 0.5, so r_viechtbauer is not bounded by 1 and any
# point-biserial above about 0.798 leaves the parameter space. Until 2.1.0 the internal
# r_trunc clamp bound at that point and the variance-stabilising z collapsed onto
# (a/2)*log((1+a)/(1-a)) -- a function of the arm-size split ALONE, carrying no
# information about the study, arriving with a plausible SE and a plausible CI. On the r
# scale Category B catches it ("[INVALID] r outside [-1, 1]"); z is unbounded, so nothing
# caught it on the metric that actually gets pooled.
#
# The package now declines the derived quantities for such a row (keeping the r-scale
# point estimate, which is what Category B flags). On a measure = "z" run those rows
# therefore either carry no estimate or fall through to the next route in the hierarchy.
#
# These pins previously asserted that EVERY row was estimated by the route under test.
# That is no longer true, and the honest re-pin is not to drop the rows but to say what
# is still guaranteed: the named route, optionally NA, and never some third route.

# `x` is a vector of info_used values. Asserts the named route, allowing NA for the
# declined rows, and forbidding any other route.
expect_route <- function(x, route) {
  testthat::expect_true(
    all(is.na(x) | x == route),
    label = paste0("every info_used is '", route, "' or NA (declined z)")
  )
  testthat::expect_true(
    any(!is.na(x) & x == route),
    label = paste0("at least one row is still estimated by '", route, "'")
  )
}

# Sign symmetry under a reverse_* flag, restricted to the rows the flag actually governs.
# A row whose z is declined falls through to a route the flag does not touch (df.haza row
# 18 lands on ancova_t when reverse_ancova_means is set), so it legitimately keeps its
# sign in both runs. Comparing it would assert that a reverse flag reverses a route it
# was never applied to.
expect_negates_on_route <- function(fwd, rev, route_prefix, col = "es_adjusted",
                                    info = "info_used_adjusted") {
  keep <- !is.na(fwd[[info]]) & startsWith(as.character(fwd[[info]]), route_prefix)
  testthat::expect_true(any(keep), label = paste0("some row still uses ", route_prefix))
  testthat::expect_equal(fwd[[col]][keep], -rev[[col]][keep], tolerance = 1e-10)
}
