##############################################################################
# Guards against arithmetically impossible dispersion / interval inputs.
#
# Two failure modes recur across the es_from_*() family. Both come from an input
# that is impossible rather than merely implausible, and neither admits a
# sensible recovery:
#
#  * a non-positive standard deviation or standard error. Where the route
#    DIVIDES by it, the effect size comes back SIGN-FLIPPED while the mean
#    difference it was built from keeps its own sign -- md = +4 returning
#    d = -2, with a perfectly healthy-looking positive d_se, so nothing
#    downstream can notice. Where the route only PROPAGATES it, the negative
#    standard error survives into the output and becomes a negative sampling
#    variance, i.e. a wrong-signed inverse-variance weight, with no other signal
#    at all.
#  * a transposed confidence interval (ci_lo > ci_up). Its width, taken as
#    ci_up - ci_lo, is negative and feeds both of the paths above.
#
# They are therefore neutralised on entry rather than patched route by route:
# .positive_or_na() blanks a non-positive dispersion, and .ci_width() takes the
# width as an absolute value -- a transposed interval describes the same
# interval, so its width is the same and it yields the same result as the
# correctly ordered one. .ordered_ci() does the matching job for interval bounds
# that are passed through to the output rather than recomputed from a standard
# error.
#
# Route new code through these helpers instead of writing the guard inline: the
# defect was originally fixed on six routes and missed on eleven others that
# have exactly the same shape, which is what a single audited entry point
# prevents. tests/testthat/test-input-guards.R sweeps every exported route that
# takes a dispersion or an interval and asserts none of them can be made to
# produce a sign flip, a negative standard error or a transposed output CI.
#
# NB these are guards against impossible input, not quality checks. Inside
# convert_df() the same values are also reported by the Tier-1 validation
# (.validate_input_data()), which flags them and -- under the default
# correct_inputs = TRUE -- sets them to NA before any route runs. The guards are
# what protects a DIRECT call to an exported es_from_*() function, and
# convert_df(correct_inputs = FALSE), where that validation deliberately
# preserves the raw values.
##############################################################################

# A standard deviation or standard error must be strictly positive; zero is
# rejected as well as negative, since dividing by it yields Inf rather than an
# effect size. NA in, NA out.
.positive_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  ifelse(!is.na(x) & x <= 0, NA_real_, x)
}

# A baseline risk is the proportion of cases in the NON-EXPOSED group, so it must lie
# in [0, 1). Anything else is neutralised to NA. NA in, NA out.
#
# WHY THE UPPER BOUND IS EXCLUSIVE, unlike the [0, 1] range check flag V11 applies to
# the same column: the value is not merely reported, it is DIVIDED BY. Grant's
# conversions are OR = RR(1 - BR)/(1 - RR*BR) and RR = OR/(1 - BR + BR*OR), and at
# BR = 1 the first has a zero numerator factor and the second collapses to OR/OR = 1 --
# measured, es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = 1, or_to_rr = "grant")
# returns logrr = 0 with logrr_se = 0 EXACTLY, i.e. an infinite inverse-variance weight.
# BR = 0 is kept: it is the rare-disease limit where both conversions become the
# identity (logrr = logor, SE preserved), which is correct rather than degenerate.
#
# WHY THIS IS A GUARD AND NOT LEFT TO THE FLAG. convert_df() range-checks the column
# (V11, internal_flags.R) but all 13 exported routes taking baseline_risk are callable
# directly and bypass it entirely, and correct_inputs = FALSE deliberately preserves a
# bad value so the user can inspect it. Unguarded, a baseline_risk of 1.5 or 20 makes
# BOTH (1 - BR) and (1 - RR*BR) negative, so their ratio stays positive and finite: the
# route then returns a plausible-looking effect size with a NEGATIVE standard error and
# a transposed confidence interval, silently. Measured before this guard:
#   es_from_rr_se(rr = 2, logrr_se = 0.2, baseline_risk = 20,  "grant") -> se = -0.00527
#   es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = 1.5, "grant") -> se = -0.04177
# both with lo > up. That is exactly the failure class this file exists to close, and
# a finiteness test cannot catch it because every value involved is finite.
.baseline_risk_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  ifelse(!is.na(x) & (!is.finite(x) | x < 0 | x >= 1), NA_real_, x)
}

# Full width of an interval, robust to transposed bounds.
.ci_width <- function(ci_lo, ci_up) abs(ci_up - ci_lo)

# Lower / upper bound of an interval, robust to transposed bounds. Used where
# the bounds are handed straight to the output instead of being rebuilt from a
# standard error.
.ci_lower <- function(ci_lo, ci_up) pmin(ci_lo, ci_up, na.rm = FALSE)
.ci_upper <- function(ci_lo, ci_up) pmax(ci_lo, ci_up, na.rm = FALSE)
