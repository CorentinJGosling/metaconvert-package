##############################################################################
# Guards against arithmetically impossible dispersion / interval inputs.
#
# Two failure modes recur across the es_from_*() family. Both come from an input
# that is impossible rather than merely implausible, and neither admits a
# sensible recovery:
#
#  * a non-positive standard deviation or standard error. Where the route
#    divides by it, the effect size comes back sign-flipped while the mean
#    difference it was built from keeps its own sign, so md = +4 returns
#    d = -2 with a perfectly healthy-looking positive d_se, and nothing
#    downstream can notice. Where the route only propagates it, the negative
#    standard error survives into the output and becomes a negative sampling
#    variance, that is, a wrong-signed inverse-variance weight, with no other
#    signal at all.
#  * a transposed confidence interval (ci_lo > ci_up). Its width, taken as
#    ci_up - ci_lo, is negative and feeds both of the paths above.
#
# They are therefore neutralised on entry rather than patched route by route.
# .positive_or_na() blanks a non-positive dispersion, .nonneg_or_na() blanks a
# negative one while keeping an exact zero (the paired-difference variances need
# that weaker predicate; see its comment), and .ci_width() takes the width as an
# absolute value: a transposed interval describes the same interval, so its width
# is the same and it yields the same result as the correctly ordered one.
# .ci_lower() / .ci_upper() (a pmin/pmax pair) do the matching job for interval
# bounds that are passed through to the output rather than recomputed from a
# standard error.
#
# Route new code through these helpers instead of writing the guard inline. The
# defect has the same shape on every route that takes a dispersion or an
# interval, and a single audited entry point is what keeps it from being fixed
# on some of them and missed on the rest. tests/testthat/test-input-guards.R
# sweeps all of those routes and asserts that none can be made to produce a sign
# flip, a negative standard error or a transposed output CI.
#
# These are guards against impossible input, not quality checks. Inside
# convert_df() the same values are also reported by the Tier-1 validation
# (.validate_input_data()), which flags them and, under the default
# correct_inputs = TRUE, sets them to NA before any route runs. The guards are
# what protects a direct call to an exported es_from_*() function, and
# convert_df(correct_inputs = FALSE), where that validation preserves the raw
# values on purpose.
##############################################################################

# A standard deviation or standard error must be strictly positive; zero is
# rejected as well as negative, since dividing by it yields Inf rather than an
# effect size. NA in, NA out.
.positive_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  ifelse(!is.na(x) & x <= 0, NA_real_, x)
}

# A standard deviation, unlike a standard error that is divided by, may legitimately be
# exactly zero, so the paired-difference variances need ">= 0, else NA" rather than the
# strict predicate above. Two callers depend on the difference: the mean-change wrappers
# zero the baseline slot by construction (es_from_mean_change_sd() passes both arms'
# mean_pre_sd = 0, es_from_mean_change_sd_single_group() passes mean_pre_sd_exp = 0), so
# .positive_or_na() there would blank the entire mean-change family.
#
# What must be neutralised is a NEGATIVE SD. It flips the sign of the cross term in
# var(post - pre) = (sd_pre^2 + sd_post^2 - 2*r*sd_pre*sd_post)/n, so the route returns a
# wrong-but-finite variance whose standard error is positive and plausible and whose CI
# is correctly ordered: none of the three properties tests/testthat/test-input-guards.R
# sweeps for can see it, and the reported weight is out by up to a factor of 3.
#
# This is the predicate .single_group_pre_post_to_smd() and .pooled_pre_post_to_smd()
# already write out inline for the standardizing SDs (R/internal_multiple_formulas.R);
# it is promoted here so the md/mdw arithmetic sitting beside them uses the same one
# rather than a second copy of it. A non-finite value is neutralised as well, exactly as
# the kernel does. NA in, NA out.
.nonneg_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  ifelse(is.finite(x) & x >= 0, x, NA_real_)
}

# Normalises one per-arm direction flag. An unsigned statistic -- an F, or a two-sided
# p-value -- carries no direction, so the routes that invert one assume both arms moved
# the same way unless the caller says otherwise. These flags are how the caller says
# otherwise, and they are optional: absent means "not reversed", which is the assumption
# the routes made before they existed, so an existing call is unaffected.
#
# `x` is NULL when the argument was missing at the call site. A scalar is recycled, as
# every other per-row argument in the package is, and an incompatible length is an error
# rather than a silent recycle: a direction flag that lands on the wrong row flips the
# sign of that row's effect size, which is exactly the failure these guards exist to stop.
# NA means "not stated", which is the same instruction as absent.
.arm_direction_flag <- function(x, n, arg_name) {
  if (is.null(x) || length(x) == 0) return(rep(FALSE, n))
  if (length(x) != n) {
    if (length(x) != 1) {
      stop("'", arg_name, "' has length ", length(x), " but ", n,
           " effect sizes were requested.", call. = FALSE)
    }
    x <- rep(x, n)
  }
  x <- as.logical(x)
  x[is.na(x)] <- FALSE
  x
}

# A baseline risk is the proportion of cases in the non-exposed group, so it must lie
# in [0, 1). Anything else is neutralised to NA. NA in, NA out.
#
# The upper bound is exclusive, unlike the [0, 1] range check flag V11 applies to the
# same column, because the value here is not merely reported but divided by. Grant's
# conversions are OR = RR(1 - BR)/(1 - RR*BR) and RR = OR/(1 - BR + BR*OR), and at
# BR = 1 the first has a zero numerator factor while the second collapses to
# OR/OR = 1: es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = 1,
# or_to_rr = "grant") returns logrr = 0 with logrr_se exactly 0, an infinite
# inverse-variance weight. BR = 0 is kept, being the rare-disease limit where both
# conversions become the identity (logrr = logor, SE preserved), which is correct
# rather than degenerate.
#
# This is a guard rather than a flag because convert_df() range-checks the column
# (V11, internal_flags.R) but all 13 exported routes taking baseline_risk are callable
# directly and bypass it entirely, and correct_inputs = FALSE preserves a bad value on
# purpose so the user can inspect it. Unguarded, a baseline_risk of 1.5 or 20 makes
# both (1 - BR) and (1 - RR*BR) negative, so their ratio stays positive and finite:
# the route then returns a plausible-looking effect size with a negative standard
# error and a transposed confidence interval, and says nothing. Without this guard:
#   es_from_rr_se(rr = 2, logrr_se = 0.2, baseline_risk = 20,  "grant") -> se = -0.00527
#   es_from_or_se(or = 2, logor_se = 0.2, baseline_risk = 1.5, "grant") -> se = -0.04177
# both with lo > up. That is the failure class this file exists to close, and a
# finiteness test cannot catch it because every value involved is finite.
.baseline_risk_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  ifelse(!is.na(x) & (!is.finite(x) | x < 0 | x >= 1), NA_real_, x)
}

# Full width of an interval, robust to transposed bounds.
.ci_width <- function(ci_lo, ci_up) abs(ci_up - ci_lo)

# A p-value that cannot be inverted.
#
# p <= 0 -- what "p < .001" becomes when it is transcribed as a number -- sends the
# inverted statistic to +Inf on every route: the t/F/point-biserial/paired-t routes
# then return d = Inf with se = Inf, and the ratio routes an se of exactly 0, i.e. an
# infinite inverse-variance weight. Neutralised on every p-value route, so that the
# exported routes (called directly by metaumbrella) behave as convert_df() does and
# do not depend on the caller's own screening.
#
# p >= 1 is degenerate only where the standard error is recovered as |estimate / z|:
# qnorm(1/2) = 0, so se = Inf and the interval is unbounded (or, rr, rd, md,
# ancova_md, mean_change, linreg_b). On the routes that invert p into a STATISTIC
# (student t, ANOVA F, ANCOVA t/F, paired t, point-biserial, chi-square) p = 1 is the
# ordinary boundary t = 0 / chi-square = 0 -- an effect size of exactly zero with a
# finite standard error -- and is kept; `se_from_ratio = FALSE` selects that case.
.pval_or_na <- function(p, se_from_ratio = TRUE) {
  if (is.null(p) || length(p) == 0) return(p)
  bad <- !is.na(p) & (p <= 0 | (se_from_ratio & p >= 1))
  ifelse(bad, NA_real_, p)
}

# Lower / upper bound of an interval, robust to transposed bounds. Used where
# the bounds are handed straight to the output instead of being rebuilt from a
# standard error.
.ci_lower <- function(ci_lo, ci_up) pmin(ci_lo, ci_up, na.rm = FALSE)
.ci_upper <- function(ci_lo, ci_up) pmax(ci_lo, ci_up, na.rm = FALSE)


# ---- unit_type -------------------------------------------------------------
#
# `unit_type` says whether `unit_increase_iv` is expressed in standard-deviation
# units or in the raw units of the independent variable. It is read in exactly one
# place, as `ifelse(unit_type == "sd", unit_increase_iv * sd_iv, unit_increase_iv)`
# in .cor_to_smd()/.cor_to_smd_vec(), and only on the `cor_to_smd = "mathur"` route.
# Any value other than "sd" therefore selects raw units, a typo included.
#
# The shipped default is "raw_scale". Six @param blocks document "sd" or "raw_scale"
# and two document "sd" or "value", so a user following ?convert_df could type a value
# that worked only by accident. All three names are accepted, with "value" and
# "raw_scale" as synonyms both meaning "already in raw units", and anything else is
# rejected rather than quietly behaving like "raw_scale".
#
# The argument stops; a column warns and is neutralised to NA. That asymmetry is the
# convention in this package: an argument is one decision by the analyst and can be
# fixed on the spot, whereas one bad cell must never abort a whole convert_df() run.
# The dipietrantonj and prop routes warn per row for the same reason.
#
# The column has to be neutralised rather than merely reported, because convert_df()
# hands it straight on to the exported es_from_*() routes, which validate their own
# argument and stop. Warning while leaving 'banana' in place would turn a per-row
# warning into a hard error two calls later: the run aborts anyway, with a message
# naming the wrong layer. NA is the right replacement, since it is already the
# internal "not supplied" value, and every non-"sd" value means raw units, which is
# exactly what the warning says those rows will get.
#
# The accepted set is built from what the repository actually uses, not from the
# documentation. A census found four spellings in real use: raw_scale (67 uses), sd
# (29), value (8) and raw_data (4). The last appears nowhere in the docs and works
# only because everything that is not "sd" means raw units, and
# tests_save/checked/test-ES-COR.R passes it, so all four are accepted.
#
# What the guard catches is a value intended as "sd" that is not spelled "sd", such
# as "SD". That produces raw units, and it is the only direction in which this
# argument can go wrong unnoticed.
.unit_type_values <- function() c("sd", "raw_scale", "value", "raw_data")

.validate_unit_type <- function(x, arg_name = "unit_type", column = FALSE) {
  if (is.null(x) || length(x) == 0) return(invisible(x))
  bad_i <- !is.na(x) & !(x %in% .unit_type_values())
  if (!any(bad_i)) return(invisible(x))
  bad <- unique(x[bad_i])
  msg <- paste0("'", paste(bad, collapse = "', '"), "' not in tolerated values for the '",
                arg_name, "' argument. Possible inputs are: ",
                paste0("'", .unit_type_values(), "'", collapse = ", "),
                " ('value' and 'raw_scale' are synonyms).")
  if (column) {
    warning(msg, " Those rows are treated as raw units, which is what every ",
            "non-'sd' value does; supply 'sd' if the increase is in SD units.",
            call. = FALSE)
    x[bad_i] <- NA
    return(invisible(x))
  }
  stop(msg, call. = FALSE)
}
