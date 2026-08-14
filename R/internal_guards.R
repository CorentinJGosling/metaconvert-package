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

# Full width of an interval, robust to transposed bounds.
.ci_width <- function(ci_lo, ci_up) abs(ci_up - ci_lo)

# Lower / upper bound of an interval, robust to transposed bounds. Used where
# the bounds are handed straight to the output instead of being rebuilt from a
# standard error.
.ci_lower <- function(ci_lo, ci_up) pmin(ci_lo, ci_up, na.rm = FALSE)
.ci_upper <- function(ci_lo, ci_up) pmax(ci_lo, ci_up, na.rm = FALSE)
