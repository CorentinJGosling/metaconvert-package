# =============================================================================
# Roadmap item 1.7 -- baseline_risk is DIVIDED BY, so it must be guarded, not just
# flagged.
#
# baseline_risk is the proportion of cases in the non-exposed group, so it belongs to
# [0, 1). Grant's conversions are
#     OR = RR(1 - BR) / (1 - RR*BR)      RR = OR / (1 - BR + BR*OR)
# and outside that range (1 - BR) and (1 - RR*BR) change sign TOGETHER, so the ratio
# stays finite and positive. The route then returned a plausible-looking effect size
# with a NEGATIVE standard error and a transposed confidence interval, silently, in
# both directions:
#     es_from_rr_se(rr = 2, logrr_se = 0.2, baseline_risk = 20,  "grant") -> se = -0.00527
#     es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = 1.5, "grant") -> se = -0.04177
# A finiteness test cannot catch this because every value involved is finite -- which is
# why the domain guard added for the RR->OR NaN case (roadmap 1.4) did not close it.
#
# convert_df() does range-check the column (flag V11), but all 13 exported routes taking
# baseline_risk are callable directly and bypass it, and correct_inputs = FALSE
# deliberately preserves a bad value. Hence a guard in R/internal_guards.R, whose
# contract is that it is unconditional.
#
# UPPER BOUND EXCLUSIVE, unlike V11's [0, 1]: at BR = 1 the OR->RR branch returns
# logrr = 0 with logrr_se = 0 EXACTLY -- an infinite inverse-variance weight.
# LOWER BOUND INCLUSIVE: BR = 0 is the rare-disease limit where both conversions become
# the identity, which is correct rather than degenerate.
# =============================================================================

.se_cols <- function(x) {
  v <- unlist(x[grepl("_se$", names(x))])
  v[is.finite(v)]
}
.inverted_ci <- function(x) {
  los <- names(x)[grepl("_ci_lo$", names(x))]
  any(vapply(los, function(lo) {
    up <- sub("_ci_lo$", "_ci_up", lo)
    if (!up %in% names(x)) return(FALSE)
    a <- x[[lo]]; b <- x[[up]]
    isTRUE(is.finite(a) && is.finite(b) && a > b)
  }, logical(1)))
}

.BAD_BR <- c(20, 1.5, 1.0, -0.3, NaN, Inf, -Inf)

# --- the defect ------------------------------------------------------------

test_that("no out-of-range baseline_risk can produce a negative standard error", {
  for (br in .BAD_BR) {
    calls <- list(
      or_se = function() es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = br,
                                       or_to_rr = "grant"),
      rr_se = function() es_from_rr_se(rr = 2, logrr_se = 0.2, baseline_risk = br,
                                       rr_to_or = "grant"),
      or_ci = function() es_from_or_ci(or = 2, or_ci_lo = 1.3, or_ci_up = 3.1,
                                       baseline_risk = br, or_to_rr = "grant"),
      rr_ci = function() es_from_rr_ci(rr = 1.2, rr_ci_lo = 1.0, rr_ci_up = 1.4,
                                       baseline_risk = br, rr_to_or = "grant"),
      rd_se = function() es_from_rd_se(rd = 0.1, rd_se = 0.03, baseline_risk = br)
    )
    for (nm in names(calls)) {
      res <- suppressWarnings(calls[[nm]]())
      expect_true(all(.se_cols(res) >= 0),
                  info = paste0(nm, " returned a negative SE at baseline_risk = ", br))
    }
  }
})

test_that("no out-of-range baseline_risk can produce a transposed interval", {
  for (br in .BAD_BR) {
    for (f in list(
      function() es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = br, or_to_rr = "grant"),
      function() es_from_rr_se(rr = 2, logrr_se = 0.2, baseline_risk = br, rr_to_or = "grant"),
      function() es_from_rd_se(rd = 0.1, rd_se = 0.03, baseline_risk = br)
    )) {
      expect_false(.inverted_ci(suppressWarnings(f())),
                   info = paste0("inverted CI at baseline_risk = ", br))
    }
  }
})

test_that("the guard itself maps the domain correctly", {
  g <- metaConvert:::.baseline_risk_or_na
  expect_equal(g(c(0, 0.001, 0.5, 0.999)), c(0, 0.001, 0.5, 0.999))   # kept
  expect_true(all(is.na(g(c(1, 1.5, 20, -0.3, Inf, -Inf, NaN)))))     # neutralised
  expect_true(is.na(g(NA_real_)))                                     # NA in, NA out
  expect_length(g(numeric(0)), 0)
})

# --- the endpoints, which are deliberately treated differently -------------

test_that("baseline_risk = 0 is kept: it is the rare-disease limit", {
  o <- es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = 0, or_to_rr = "grant")
  expect_equal(o$logrr, log(2), tolerance = 1e-10)
  expect_equal(o$logrr_se, 0.2, tolerance = 1e-10)
  r <- es_from_rr_se(rr = 2, logrr_se = 0.2, baseline_risk = 0, rr_to_or = "grant")
  expect_equal(r$logor, log(2), tolerance = 1e-10)
})

test_that("baseline_risk = 1 is rejected: it yields SE = 0, an infinite weight", {
  o <- suppressWarnings(es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = 1,
                                      or_to_rr = "grant"))
  expect_true(is.na(o$logrr) || !isTRUE(o$logrr_se == 0))
})

# --- no regression in the valid range --------------------------------------

test_that("valid baseline_risk values are untouched", {
  o <- es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = 0.2, or_to_rr = "grant")
  expect_equal(o$logrr, log(2 / (1 - 0.2 + 0.2 * 2)), tolerance = 1e-10)
  expect_true(is.finite(o$logrr_se) && o$logrr_se > 0)

  r <- es_from_rr_se(rr = 1.2, logrr_se = 0.15, baseline_risk = 0.3,
                     rr_to_or = "grant")
  expect_equal(r$logor, log(1.2 * (1 - 0.3) / (1 - 1.2 * 0.3)), tolerance = 1e-10)
  expect_true(is.finite(r$logor_se) && r$logor_se > 0)
})

# --- the RD / NNT path -----------------------------------------------------

test_that("RD and NNT are NA when the implied exposed risk exceeds 1", {
  # rr * baseline_risk = 1.8: the risk pair does not exist. The OR path already
  # declined here; the RD path used to emit rd = -0.900 and nnt = -1.111, and B6 never
  # fired because |rd| < 1.
  r <- suppressWarnings(es_from_rr_se(rr = 2, logrr_se = 0.2, baseline_risk = 0.9))
  expect_true(is.na(r$rd))
  expect_true(is.na(r$nnt))
  expect_true(is.na(r$logor))   # unchanged: the grant guard already covered this
})

test_that("RD and NNT survive when the implied exposed risk is valid", {
  r <- es_from_rr_se(rr = 1.05, logrr_se = 0.2, baseline_risk = 0.5)  # implied 0.525
  expect_equal(r$rd, 0.5 * (1 - 1.05), tolerance = 1e-10)
  expect_true(is.finite(r$nnt))
  # and exactly at the boundary rr * br == 1 the risk pair is still admissible
  r2 <- es_from_rr_se(rr = 2, logrr_se = 0.2, baseline_risk = 0.5)
  expect_true(is.finite(r2$rd))
})

# --- terminal messages name only reachable methods -------------------------

test_that("the method-dispatch errors advertise only values that exist", {
  m <- tryCatch(
    metaConvert:::.rr_to_or(rr = 2, logrr_se = .2, rr_ci_lo = 1, rr_ci_up = 3,
                            n_cases = NA, n_controls = NA, n_exp = NA, n_nexp = NA,
                            baseline_risk = .2, rr_to_or = "nope"),
    error = function(e) conditionMessage(e))
  expect_match(m, "'grant'")
  expect_match(m, "'dipietrantonj'")
  expect_false(grepl("grant_2x2|grant_CI", m))
})

test_that(".or_to_rr errors on an unknown method instead of returning NULL", {
  expect_error(
    metaConvert:::.or_to_rr(or = 2, logor_se = .2, or_ci_lo = 1, or_ci_up = 3,
                            n_cases = NA, n_controls = NA, n_exp = NA, n_nexp = NA,
                            baseline_risk = .2, or_to_rr = "nope"),
    "not in tolerated values")
})
