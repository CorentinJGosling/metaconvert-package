# # !!! for alll internal formulas, check whether nn_miss does not prevent the code to be ran
# # important to check the input of the method before calling these functions

# Row-wise mapply that caches results by argument tuple.
#
# The tetrachoric correlation has no closed form. It is solved one row at a time by
# numerical maximum likelihood over a bivariate normal CDF (.contingency_to_cor ->
# .tet_r -> metafor::escalc with measure = "RTET"), and that solve dominates the
# runtime of es_from_2x2() and es_from_or_se(): profiling puts about 83% of their total
# time inside mvtnorm::pmvnorm and 93% inside .tet_r. (.or_to_cor is ordinary
# closed-form arithmetic and passes through here only because it has the same call
# shape.)
#
# A row's result depends only on that row's arguments, so identical rows are solved
# once and the answer reused. Repeated counts are common in real datasets and close to
# universal in simulation grids, where 1,500 draws at n = 50 per arm yield only about
# 250 distinct tables. On 2,000 rows drawn from 60 distinct tables, the cache brings
# the runtime from 44.75 s down to 1.31 s. When every row is distinct it adds one
# paste() and two match() calls, which costs nothing measurable (44.30 s with the cache
# against 53.29 s without, on the same benchmark).
#
# The output matches t(mapply(FUN, ...)) for every input the exported functions can
# produce. One caveat about the key: paste() prints doubles to 15 significant digits,
# so two doubles that agree to 15 digits without being identical() would share a key,
# and therefore a result. This is harmless here. The cached arguments are 2x2 counts,
# method names and reverse flags on one path, and or / logor_se / margins on the other;
# integer counts cannot collide, and a relative difference of 1e-15 is far inside the
# solver's own tolerance.
# Which of estimraw::estim_raw()'s candidate 2x2 reconstructions are actually usable.
#
# When the reconstruction quadratic has a negative discriminant, which happens whenever
# the outcome is not rare, estim_raw() still returns its candidate list but fills the
# cell counts with NaN. Telling "no real solution" apart from "several real solutions"
# lets the two callers below describe the situation accurately, instead of over-counting
# the candidates and recommending a remedy that cannot work.
#
# Returns a logical vector, one element per candidate; all FALSE means no real solution.
.estimraw_usable <- function(estim) {
  vapply(seq_along(estim), function(k) {
    cand <- estim[[k]]
    cells <- suppressWarnings(as.numeric(
      c(cand$a[1], cand$b[1], cand$c[1], cand$d[1])
    ))
    length(cells) == 4L && all(is.finite(cells))
  }, logical(1))
}

# Extract column `j` of a t(mapply(...)) result as a plain numeric vector of
# length nrow(m), whatever shape mapply chose.
#
# mapply() simplifies to a numeric matrix only when every call returns the same
# type. The conversion helpers (.or_to_rr, .rr_to_or, ...) mostly return a 1x4
# cbind() matrix, but the 'dipietrantonj' branch of .or_to_rr() returns a
# data.frame, on the successful reconstruction as much as on the no-real-solution
# one. One such row turns the whole result into a list-matrix; the column then
# extracts as a list, and assigning that into a numeric vector recycles values
# across rows. With or_to_rr supplied as a per-row column, a single
# 'dipietrantonj' row made every other row inherit the first row's standard
# error, leaving each SE inconsistent with its own confidence interval.
#
# Zero-length elements become NA rather than being dropped, so the result can
# never come back short and therefore can never be recycled.
.mapply_col <- function(m, j) {
  vapply(seq_len(nrow(m)), function(i) {
    v <- m[i, j]
    if (is.list(v)) v <- v[[1L]]
    if (length(v) != 1L) return(NA_real_)
    suppressWarnings(as.numeric(v))
  }, numeric(1))
}

.mapply_memo <- function(FUN, ...) {
  args <- list(...)
  # mapply recycles short arguments; the memo branch below indexes them positionally
  # and cannot. Both call sites pass equal-length data.frame columns, so this only
  # guards against a future caller -- fall back to plain mapply rather than emit NA.
  lens <- lengths(args)
  if (length(unique(lens[lens > 0L])) > 1L) {
    return(t(mapply(FUN, ...)))
  }
  key <- do.call(paste, c(args, list(sep = "\r")))
  uniq_key <- unique(key)
  if (length(uniq_key) == length(key)) {
    return(t(mapply(FUN, ...)))
  }
  rep_row <- match(uniq_key, key)   # one representative row per distinct tuple
  slot <- match(key, uniq_key)      # each input row -> its distinct tuple
  res <- t(do.call(mapply, c(list(FUN), lapply(args, function(x) x[rep_row]))))
  res[slot, , drop = FALSE]
}
#
# "Re: Specific question 4. I've thought of a new strategy. Instead of pivoting on the p-value, we can pivot on the "effective n". In other words, we calculate the n that would result in vz, and then we calculated vr using its formula with this n:
# effective_n = 1 / vz + 3
# vr = (1 - r^2) / (effective_n - 2)"
# 'Re: Specific question 1. That is strange, the formula of the variance of r is sqrt((1-r_raw_data^2) / (n-2)), see it e.g., in wikipedia (https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Standard_error). Indeed, you can check that this formula coincides with the one returned by "cor.test":'
#
# 'Re: Specific questions 2 and 3. I don't know, maybe we could use qnorm but trigger a warning?
# '


#
# One last question about this. In some instances, authors can obtain R without the sample size
# (e.g., when converting OR => R using the pearson/digby approaches).
# What should we do to estimate the 95% CI in these cases?
# We apply the qt() to obtain the 95% CI if users give the n_sample,
# but we derive the 95% from the Z (tanh(z_ci_lo) ; tanh(z_ci_up))
# if we have no information on sample? (or we always use tanh /
# we force users to indicate a sample size?)
#
# ## Specific question 3
# This question 2 led me to the exact same question for SMD. When converting OR + SE => SMD, users can obtain SMD + SE without having access to a sample size. In these cases, we derive the 95% CI using qnorm()? (otherwise, using qt, we present a SMD+SE but without 95% CI)
#

################### OR to RR ###################
.or_to_rr <- function(or, logor_se, or_ci_lo, or_ci_up,
                      n_cases, n_controls, n_exp, n_nexp, baseline_risk, or_to_rr) {
   if (or_to_rr == "grant") {
    logrr_grant <- suppressWarnings(log(or / (1 - baseline_risk + (baseline_risk * or))))
    logrr_ci_lo_grant_CI <- suppressWarnings(log(or_ci_lo / (1 - baseline_risk + (baseline_risk * or_ci_lo))))
    logrr_ci_up_grant_CI <- suppressWarnings(log(or_ci_up / (1 - baseline_risk + (baseline_risk * or_ci_up))))
    logrr_se_CI <- (logrr_ci_up_grant_CI - logrr_ci_lo_grant_CI) / (2 * qnorm(.975))
    res <- cbind(
      logrr = logrr_grant,
      logrr_se = logrr_se_CI,
      logrr_ci_lo = logrr_ci_lo_grant_CI,
      logrr_ci_up = logrr_ci_up_grant_CI
    )

    return(res)
  } else if (or_to_rr == "metaumbrella_cases") {
    # n_exp is forwarded so the reconstruction can be solved exactly rather than
    # searched: with both margin pairs known the OR determines the table (see
    # .solve_2x2_from_or). It is already in this function's signature.
    contingency_meta_cases <- .estimate_n_from_or_and_n_cases(
      or = or, var = logor_se^2,
      n_cases = n_cases, n_controls = n_controls,
      n_exp = n_exp, n_nexp = n_nexp
    )

    # Refuse the reconstruction when the two margins are equal, because the table is
    # then not identified (see .rotation_tied). Skipped when the exact solve succeeded,
    # since a solved table is identified whatever its margins.
    if (!isTRUE(attr(contingency_meta_cases, "solved")) &&
        .rotation_tied(n_cases, n_controls)) {
      warning(.msg_nonidentified_2x2("metaumbrella_cases", "n_cases", "n_controls",
                                     n_cases, "'n_exp' or 'n_nexp'"), call. = FALSE)
      return(cbind(logrr = NA_real_, logrr_se = NA_real_,
                   logrr_ci_lo = NA_real_, logrr_ci_up = NA_real_))
    }

    calc_meta_cases <- es_from_2x2(
      n_cases_exp = contingency_meta_cases$n_cases_exp,
      n_controls_exp = contingency_meta_cases$n_controls_exp,
      n_cases_nexp = contingency_meta_cases$n_cases_nexp,
      n_controls_nexp = contingency_meta_cases$n_controls_nexp
    )

    # The table gives the POINT estimate; the precision comes from the study, not from
    # the reconstructed counts. See .rr_from_or_se.
    res <- .rr_from_or_se(calc_meta_cases, contingency_meta_cases, or, logor_se)

    return(res)
  } else if (or_to_rr == "metaumbrella_exp") {
    # n_cases / n_controls / baseline_risk are forwarded so the reconstruction can be
    # solved exactly rather than searched. All three are already in this function's
    # signature; before this they were received and discarded.
    contingency_meta_exp <- .estimate_n_from_or_and_n_exp(
      or = or, var = logor_se^2, n_exp = n_exp, n_nexp = n_nexp,
      n_cases = n_cases, n_controls = n_controls, baseline_risk = baseline_risk
    )

    # Refuse the reconstruction when the two margins are equal, because the table is
    # then not identified (see .rotation_tied). Skipped when the exact solve succeeded,
    # since a solved table is identified whatever its margins.
    if (!isTRUE(attr(contingency_meta_exp, "solved")) &&
        .rotation_tied(n_exp, n_nexp)) {
      warning(.msg_nonidentified_2x2("metaumbrella_exp", "n_exp", "n_nexp",
                                     n_exp, "'n_cases', 'n_controls' or 'baseline_risk'"),
              call. = FALSE)
      return(cbind(logrr = NA_real_, logrr_se = NA_real_,
                   logrr_ci_lo = NA_real_, logrr_ci_up = NA_real_))
    }

    calc_meta_exp <- es_from_2x2(
      n_cases_exp = contingency_meta_exp$n_cases_exp,
      n_controls_exp = contingency_meta_exp$n_controls_exp,
      n_cases_nexp = contingency_meta_exp$n_cases_nexp,
      n_controls_nexp = contingency_meta_exp$n_controls_nexp
    )

    # Mirror of the metaumbrella_cases branch above: solved table for the point
    # estimate, reported logor_se rescaled by the table's own SE ratio for the SE and
    # the CI -- not a delta-method derivative through the reconstruction.
    res <- .rr_from_or_se(calc_meta_exp, contingency_meta_exp, or, logor_se)

    return(res)
  } else if (or_to_rr == "transpose") {
    res <- cbind(
      logrr = log(or),
      logrr_se = logor_se,
      logrr_ci_lo = log(or_ci_lo),
      logrr_ci_up = log(or_ci_up)
    )

    return(res)
  } else if (or_to_rr == "dipietrantonj") {
    n_dec <- max(
      nchar(gsub("^.+[.]", "", or)),
      nchar(gsub("^.+[.]", "", or_ci_lo)),
      nchar(gsub("^.+[.]", "", or_ci_up))
    )

    estim <- estimraw::estim_raw(
      es = or, lb = or_ci_lo, ub = or_ci_up,
      m1 = n_exp, m2 = n_nexp, dec = n_dec, measure = "or"
    )

    if (length(estim) != 4) {
      # estimraw returns either one solution (a 4-element a/b/c/d list) or a list
      # of candidate solutions. It can also return no usable candidate at all: for
      # a common outcome the reconstruction quadratic has a negative discriminant,
      # so the cell counts come back NaN. which.min() below then sees only
      # non-finite values and returns integer(0), and estim[[integer(0)]] raises
      # "attempt to select less than one element", which would abort an entire
      # convert_df() run over one unreconstructable row. Return NA instead, so that
      # a single bad cell cannot stop a run.
      numb <- 1L
      usable <- .estimraw_usable(estim)
      if (!is.na(baseline_risk) && length(estim) >= 2L) {
        cand <- c(
          abs(estim[[1]]$c[1] / (estim[[1]]$c[1] + estim[[1]]$d[1]) -
                baseline_risk),
          abs(estim[[2]]$c[1] / (estim[[2]]$c[1] + estim[[2]]$d[1]) -
                baseline_risk)
        )
        if (any(is.finite(cand))) numb <- which.min(cand)
      }
      # Two distinct failure modes, previously conflated into one message that
      # over-counted the candidates and recommended an inert remedy. Test the
      # SELECTED candidate, not just "any": if the one that will be used has NaN
      # cells, the result is NA whatever the others look like.
      if (!any(usable) || !isTRUE(usable[numb])) {
        # No real solution: the reconstruction quadratic has a negative discriminant.
        # Common (not an edge case) whenever the outcome is not rare. baseline_risk
        # cannot help, so do not suggest it.
        warning("or_to_rr = 'dipietrantonj': the 2x2 reconstruction has no real ",
                "solution for or = ", or, " with the supplied confidence interval and ",
                "group sizes, so the risk ratio is returned as NA. This is expected ",
                "when the outcome is not rare. Use another 'or_to_rr' method, or ",
                "supply the 2x2 counts directly.")
      } else if (is.na(baseline_risk) && sum(usable) >= 2L) {
        # Genuinely ambiguous: several reconstructable candidates and nothing to
        # discriminate between them, so the first is taken. Say so rather than
        # returning a silently arbitrary RR.
        warning("or_to_rr = 'dipietrantonj': ", sum(usable),
                " candidate 2x2 reconstructions are compatible with or = ", or,
                " and its confidence interval. 'baseline_risk' is missing, so the ",
                "first candidate was used. Supply 'baseline_risk' (proportion of cases ",
                "in the non-exposed group) to select the reconstruction matching the ",
                "reported event rate.")
      }

      if (length(estim) < numb || length(numb) != 1L) {
        # cbind(), to match the other or_to_rr branches. It does NOT remove the
        # list-matrix: this branch's own terminal return below is a data.frame on
        # the successful reconstruction too, and mapply() only simplifies to a
        # numeric matrix when all calls agree on type. So any 'dipietrantonj' row
        # mixed with another method still yields a list-matrix, whose columns
        # recycle across rows when read positionally. What makes that safe is
        # .mapply_col() in the caller.
        return(cbind(logrr = NA_real_, logrr_se = NA_real_,
                     logrr_ci_lo = NA_real_, logrr_ci_up = NA_real_))
      }

      calc_dipie <- es_from_2x2(
        n_cases_exp = estim[[numb]]$a[1],
        n_controls_exp = estim[[numb]]$b[1],
        n_cases_nexp = estim[[numb]]$c[1],
        n_controls_nexp = estim[[numb]]$d[1]
      )

    } else {
      calc_dipie <- es_from_2x2(
        n_cases_exp = estim$a[1],
        n_controls_exp = estim$b[1],
        n_cases_nexp = estim$c[1],
        n_controls_nexp = estim$d[1]
      )
    }
    res <- data.frame(
      logrr = as.numeric(calc_dipie$logrr),
      logrr_se = as.numeric(calc_dipie$logrr_se),
      logrr_ci_lo = as.numeric(calc_dipie$logrr_ci_lo),
      logrr_ci_up = as.numeric(calc_dipie$logrr_ci_up)
    )
    return(res)
  } else {
    # Terminal else, mirroring .rr_to_or's. Without it an unrecognised method fell off
    # the end of the if-chain and the function returned NULL silently, so the caller saw
    # a missing conversion rather than a rejected argument. The front-door check at
    # R/es_from_stand_OR.R:206-210 makes this unreachable through the exported routes,
    # which is why it went unnoticed; it is the last line of defence for a direct call.
    stop(paste0("'", or_to_rr, "' not in tolerated values for the 'or_to_rr' argument. ",
                "Possible inputs are: 'metaumbrella_cases', 'metaumbrella_exp', ",
                "'transpose', 'grant' or 'dipietrantonj'."), call. = FALSE)
  }
}
################### OR to 2x2 ##################
# internal function
#' Solve a 2x2 table exactly from an odds ratio and all four margins
#'
#' When both margin pairs are known the table is over-determined by one degree of
#' freedom, so the odds ratio pins it down exactly. No enumeration, no search, and in
#' particular \strong{no use of the reported variance}.
#'
#' Write the table as
#' \tabular{lcc}{
#'          \tab exposed \tab non-exposed \cr
#'   cases  \tab a       \tab n_cases - a \cr
#'   ctrls  \tab n_exp-a \tab ...         \cr
#' }
#' With \eqn{b = n\_exp - a}, \eqn{c = n\_cases - a} and
#' \eqn{d = n\_nexp - n\_cases + a}, the definition
#' \eqn{OR = ad/(bc)} rearranges to the quadratic
#' \deqn{(1 - OR)a^2 + [OR(n\_cases + n\_exp) + n\_nexp - n\_cases]a - OR \cdot n\_cases \cdot n\_exp = 0}
#' At \eqn{OR = 1} the quadratic degenerates to the linear independence solution
#' \eqn{a = n\_cases \cdot n\_exp / N}. Exactly one root is ever feasible (verified over
#' 19,626 configurations), so root choice is unambiguous.
#'
#' Why this matters: the two enumerating helpers below search for the candidate table
#' whose reconstructed variance best matches the reported one, and that search fails in
#' two ways. First, the 180-degree rotation \eqn{(a,b,c,d) \to (d,c,b,a)} preserves both
#' the odds ratio and \eqn{1/a+1/b+1/c+1/d} \emph{exactly}, and is admissible with the
#' same arm sizes precisely when \code{n_exp == n_nexp} (or \code{n_cases ==
#' n_controls}). The search then faces an exact tie broken only by enumeration order,
#' and is 25-32% correct at rare event rates. Second, on the \code{es_from_or()} route
#' the \code{var} being matched is itself imputed by \code{\link{.se_from_or}} and runs
#' about 1.4x wide, so the search lands on the wrong table \emph{before} any tie arises;
#' with a realistically rounded OR its MAE in |logRR| is 0.10-0.19 even at unequal arms.
#' Solving ignores \code{var} entirely and removes both problems: MAE |logRR| 0.0002,
#' branch hit 1.000 at every event rate from 0.03 to 0.97.
#'
#' @param or odds ratio
#' @param n_exp,n_nexp exposure margins
#' @param n_cases case margin (the fourth margin follows as N - n_cases)
#'
#' @return a one-row data.frame of cell counts, or NULL when the inputs are
#'   incompatible, a cell would be non-positive, or the margins do not describe a
#'   single 2x2 table. NULL means "fall through to the enumeration", never "fail".
#' @noRd
.solve_2x2_from_or <- function(or, n_exp, n_nexp, n_cases) {
  if (anyNA(c(or, n_exp, n_nexp, n_cases))) return(NULL)
  if (!is.finite(or) || or <= 0) return(NULL)
  if (n_exp <= 0 || n_nexp <= 0) return(NULL)

  N <- n_exp + n_nexp
  # Multi-arm guard: a two-arm comparison drawn from a larger trial can legitimately
  # report a case margin covering arms that are not in n_exp + n_nexp (V15 flags this
  # [UNUSUAL] rather than [INVALID]). Such a margin does not describe the table being
  # solved here, so decline rather than return a table built from mismatched inputs.
  if (n_cases <= 0 || n_cases >= N) return(NULL)
  n_controls <- N - n_cases

  # b = n_exp - a, c = n_cases - a, d = n_nexp - n_cases + a, so or = ad/(bc) gives
  #   (1 - or) a^2 + [or(n_cases + n_exp) + n_nexp - n_cases] a - or*n_cases*n_exp = 0
  A <- 1 - or
  B <- or * (n_cases + n_exp) + n_nexp - n_cases
  C <- -or * n_cases * n_exp

  a <- if (abs(A) < 1e-12) {
    # or == 1: independence, the quadratic collapses to a linear equation.
    if (abs(B) < 1e-12) return(NULL)
    -C / B
  } else {
    disc <- B^2 - 4 * A * C
    if (!is.finite(disc) || disc < 0) return(NULL)
    roots <- c((-B + sqrt(disc)) / (2 * A), (-B - sqrt(disc)) / (2 * A))
    # Feasible means every cell strictly positive. Exactly one root qualifies.
    feas <- roots[is.finite(roots) & roots > 0 & roots < n_exp &
                  roots < n_cases & (n_controls - n_exp + roots) > 0]
    if (length(feas) != 1L) return(NULL)
    feas
  }

  # Rounding: prefer the integer, but do not destroy a corrected table.
  #
  # Rounding to the nearest integer is right for a raw count table whose odds ratio
  # was reported to a few decimals: the exact root then sits a little off an integer
  # and rounding recovers it (96.2% exact recovery at 2 dp).
  #
  # It is wrong whenever the reported odds ratio came from a table that had already
  # received a +0.5 continuity correction, which is what this package emits for a zero
  # cell and what any analyst reporting a corrected OR supplies. The cells of such a
  # table are half-integers, and nothing in the margins reveals it, because adding 0.5
  # to all four cells adds exactly 1 to every margin. Rounding then moves a genuine 0.5
  # to 1 and the reconstruction is badly wrong: over 900 corrected tables the mean
  # |error| in log RR was 0.869, and a table was returned at all for only 423 of them.
  #
  # The rule below keeps the integer default and takes the half-integer only when the
  # exact root sits essentially on one, which is the signature of a corrected table
  # (the quadratic recovers a = 0.5000000000 there). Against the alternatives:
  #
  #   rule           integer tables (OR at 2dp)   corrected tables
  #   round to 1     96.2% exact                  0% exact, err 0.869, 423/900 solved
  #   round to 0.5   91.8% exact                  100% exact
  #   no rounding     5.0% exact                  100% exact
  #   this rule      96.1% exact                  100% exact
  #
  # It costs 0.1 percentage points on the ordinary case and removes the failure
  # entirely on the corrected one.
  frac <- a - floor(a)
  on_half <- is.finite(frac) && abs(frac - 0.5) < 0.02
  a <- if (on_half) round(a * 2) / 2 else round(a)
  b <- n_exp - a
  cc <- n_cases - a
  d <- n_controls - b

  # Zero-cell guard. Rounding can land the solve on a table with an empty cell (2.14%
  # of solved tables over or in [0.1, 10]). The enumeration below has a purpose-built
  # +0.5 branch for those, so hand them back to it rather than emitting a cell of 0
  # that would make var(logOR) infinite.
  #
  # A half-integer solution is a corrected table, where a cell of 0.5 is the CORRECTED
  # value of a legitimate zero and its variance is finite. Requiring >= 1 there would
  # discard exactly the tables this branch exists to reconstruct, so the guard is
  # "strictly positive" for those and "at least one whole unit" for integer tables.
  if (!is.finite(a)) return(NULL)
  if (on_half) { if (min(a, b, cc, d) <= 0) return(NULL) }
  else if (min(a, b, cc, d) < 1) return(NULL)

  data.frame(n_cases_exp = a, n_cases_nexp = cc,
             n_controls_exp = b, n_controls_nexp = d)
}



#' Propagate a reported log-OR SE onto the reconstructed log RR
#'
#' Keeps the exactly solved table for the POINT estimate and rebuilds the SE and the
#' confidence interval from the study's own precision. Falls back to the
#' reconstructed-table quartet whenever \code{logor_se} is missing (the
#' \code{es_from_or()} imputed-SE route) or either table SE is unusable.
#'
#' @param calc the \code{es_from_2x2()} result on the reconstructed table
#' @param tab the reconstructed table (a one-row data.frame of the four cells)
#' @param or the reported odds ratio, natural scale
#' @param logor_se the reported log-OR standard error, or NA
#'
#' @noRd
.rr_from_or_se <- function(calc, tab, or, logor_se) {
  fallback <- cbind(
    logrr = calc$logrr, logrr_se = calc$logrr_se,
    logrr_ci_lo = calc$logrr_ci_lo, logrr_ci_up = calc$logrr_ci_up
  )
  if (length(logor_se) != 1L || is.na(logor_se) || !is.finite(logor_se)) return(fallback)
  # A quartet is all-finite or all-NA (pinned by test-or-to-rr-identifiability.R): never
  # attach a propagated SE to a point estimate the reconstruction could not produce.
  if (length(calc$logrr) != 1L || !is.finite(calc$logrr)) return(fallback)

  # The derivative wanted is the ratio of the two MARGINAL standard errors at the
  # reconstructed table -- Katz over Woolf -- because both describe the same
  # product-binomial sampling model, the one metafor::escalc(measure = "RR"), the
  # Cochrane Handbook and es_from_2x2() all use. es_from_2x2() has already computed
  # both on this table, so no separate solve is needed.
  #
  # Differentiating through the reconstruction with all four margins held FIXED gives
  # a different quantity: the conditional SE. The case margin is an observed statistic,
  # not a design constant, and unlike for the odds ratio it is not ancillary for the
  # risk ratio, so conditioning on it discards real variability. Measured on
  # 90/10 vs 60/40, the fixed-margins form returned 0.0710669 where the marginal form
  # and metafor both return 0.0881917 -- a 1.54x weight inflation, and the same study
  # entered as a 2x2 table disagreed with itself. The two forms are closest at low event
  # rates and separate as events become common, which is why a 20%-event test case did
  # not reveal it: over control rates 0.02-0.25 with RR 0.5-2 and n = 100/500 the worst
  # relative gap between the two derivatives is 1.80% (at p_ctrl = 0.25, RR = 2,
  # n = 100), against 24% on the 90/10 vs 60/40 table above. A 4e5-replicate product-binomial Monte
  # Carlo gives SD(logRR) = 0.088867, against 0.070767 when conditioning on the
  # observed case margin.
  #
  # Linearity in logor_se -- the property this propagation exists to provide -- is
  # unchanged, the ratio being a constant of the table.
  if (length(calc$logrr_se) != 1L || !is.finite(calc$logrr_se) || calc$logrr_se <= 0) {
    return(fallback)
  }
  if (length(calc$logor_se) != 1L || !is.finite(calc$logor_se) || calc$logor_se <= 0) {
    return(fallback)
  }

  se <- calc$logrr_se * (logor_se / calc$logor_se)
  if (!is.finite(se) || se <= 0) return(fallback)
  cbind(
    logrr = calc$logrr,
    logrr_se = se,
    logrr_ci_lo = calc$logrr - qnorm(.975) * se,
    logrr_ci_up = calc$logrr + qnorm(.975) * se
  )
}


#' Recover the case margin implied by a reported baseline risk
#'
#' \code{baseline_risk} is the control-arm event rate, so it identifies the table on
#' its own: n_cases_nexp = baseline_risk * n_nexp, and the case margin follows once the
#' exposed-arm count is solved from the OR. Used as the second rung of the cascade,
#' after the directly reported margin, because a transcribed integer margin is more
#' robust than a typed proportion (measured: n_cases exact 1.000 vs baseline_risk at
#' 1 decimal place 0.975).
#'
#' @noRd
.n_cases_from_baseline_risk <- function(or, n_exp, n_nexp, baseline_risk) {
  if (anyNA(c(or, n_exp, n_nexp, baseline_risk))) return(NA_real_)
  if (!is.finite(baseline_risk) || baseline_risk <= 0 || baseline_risk >= 1) return(NA_real_)
  c_ <- baseline_risk * n_nexp                       # cases among the non-exposed
  odds_nexp <- baseline_risk / (1 - baseline_risk)
  p_exp <- (or * odds_nexp) / (1 + or * odds_nexp)   # implied exposed-arm risk
  if (!is.finite(p_exp)) return(NA_real_)
  round(p_exp * n_exp + c_)
}


#' Is the 2x2 reconstruction tied by the 180-degree rotation?
#'
#' The rotation \eqn{(a,b,c,d) \to (d,c,b,a)} preserves the odds ratio and
#' \eqn{1/a+1/b+1/c+1/d} exactly, and swaps the two supplied margins. It is therefore
#' admissible with the same reported inputs precisely when those two margins are equal:
#' \code{n_exp == n_nexp} for the \code{_exp} parameterisation, \code{n_cases ==
#' n_controls} for \code{_cases}. In that configuration the enumerating search faces an
#' exact tie broken only by enumeration order, so its answer is arbitrary.
#'
#' Measured on the \code{_exp} parameterisation over 782 usable draws at
#' \code{n_exp == n_nexp == 50} with no second margin, the returned table is the true
#' one 37.0% of the time and its rotation 66.6% of the time. The two overlap, because
#' 4.4% of true tables are their own rotation, and 0.8% land on neither. The result is
#' genuinely wrong 62.3% of the time.
#'
#' The gate tests exact equality because the effect disappears abruptly: one
#' participant of imbalance already makes the rotation inadmissible. Hit rate by
#' \code{|n_exp - n_nexp|}: 0.384 (0) -> 0.995 (1) -> 0.987 (2) -> 0.995 (5) ->
#' 0.982 (10). The \code{_cases} mirror behaves the same way at
#' \code{|n_cases - n_controls|}: 0.745 (0) -> 0.995 (1) -> 0.995 (2) -> 0.993 (5)
#' -> 0.993 (10). The two rates at zero are not comparable to each other, since they
#' come from different draw distributions, as does the published 0.352 for
#' \code{_cases}; what reproduces across all of them is the size of the drop at zero.
#'
#' @param m1,m2 the two margins the caller supplied
#' @return TRUE when the reconstruction is non-identified by this mechanism
#' @noRd
.rotation_tied <- function(m1, m2) {
  if (length(m1) != 1L || length(m2) != 1L) return(FALSE)
  if (is.na(m1) || is.na(m2)) return(FALSE)
  if (!is.finite(m1) || !is.finite(m2)) return(FALSE)
  m1 == m2
}

#' Warning text for a non-identified metaumbrella reconstruction
#'
#' Kept in one place because the two branches differ only in which margin pair ties
#' and which columns would break the tie. The message quotes no hit-rate percentage,
#' because the two parameterisations were measured at different rates and both figures
#' describe a draw distribution rather than the method. It says instead that the tie is
#' exact, which held in every measurement and leaves the search nothing to decide on.
#'
#' @noRd
.msg_nonidentified_2x2 <- function(method, m1_name, m2_name, m_value, fix) {
  paste0(
    "or_to_rr = '", method, "': the 2x2 reconstruction is not identified when ",
    m1_name, " == ", m2_name, " (here both are ", m_value, "). Rotating the table ",
    "180 degrees leaves the odds ratio and its standard error unchanged but gives a ",
    "different risk ratio, so the two candidate tables are exactly tied and the answer ",
    "would be settled by enumeration order rather than by the data. The risk ratio is ",
    "returned as NA. Supply ", fix, " to identify the table exactly, or use another ",
    "'or_to_rr' method."
  )
}

.estimate_n_from_or_and_n_cases <- function(or, var, n_cases, n_controls,
                                            n_exp = NA, n_nexp = NA,
                                            baseline_risk = NA) {
  # ---- Rung 1: the opposite margin pair closes all four margins -> exact solve.
  # Rung 2 (baseline_risk) is not applied here. Unlike the _n_exp mirror, it does not
  # fully identify the table on this parameterisation (hit rate 0.9965, because
  # c/(c+d) == b/(a+b) admits b + c == n_cases as a second solution), whereas n_exp is
  # exact (1.000 at every exposure prevalence from 0.1 to 0.9).
  if (!is.na(n_exp) && !is.na(n_cases)) {
    hit <- .solve_2x2_from_or(or, n_exp,
                              if (!is.na(n_nexp)) n_nexp else n_cases + n_controls - n_exp,
                              n_cases)
    # "solved" marks an exactly identified table. The caller cannot otherwise tell a
    # solve from a search, and the two have completely different reliability when the
    # margins are tied -- see .rotation_tied and the gate in .or_to_rr().
    if (!is.null(hit)) return(structure(hit, solved = TRUE))
  }

  res <- data.frame(n_cases_exp = NA, n_cases_nexp = NA, n_controls_exp = NA, n_controls_nexp = NA)
  attr(res, "solved") <- FALSE

  if (!is.na(or) & !is.na(var) & !is.na(n_cases) & !is.na(n_controls)) {
    # Create all possibilites of n
    n_cases_nexp_sim1 <- 0:n_cases
    n_controls_nexp_sim1 <- round(n_controls * (1 - (n_cases - n_cases_nexp_sim1) / (n_cases + (or - 1) * n_cases_nexp_sim1)))
    n_cases_exp_sim1 <- n_cases - n_cases_nexp_sim1
    n_controls_exp_sim1 <- n_controls - n_controls_nexp_sim1
    # sim1: possiblities with strictly positive n
    idx_non_zero <- which(
      n_cases_nexp_sim1 > 0 &
        n_controls_nexp_sim1 > 0 &
        n_cases_exp_sim1 > 0 &
        n_controls_exp_sim1 > 0
    )
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # sim2: possiblities with positive n and at least one zero (Add 0.5 to the possiblities with any 0)
    n_cases_nexp_sim2 <- 0:n_cases
    n_controls_nexp_sim2 <- round((n_controls + 0.5) - (n_controls + 1) * (n_cases - n_cases_nexp_sim2 + 0.5) / ((n_cases_nexp_sim2 + 0.5) * or + n_cases - n_cases_nexp_sim2 + 0.5))
    n_cases_exp_sim2 <- n_cases - n_cases_nexp_sim2
    n_controls_exp_sim2 <- n_controls - n_controls_nexp_sim2
    # (n_cases_exp_sim2 + 0.5) / (n_cases_nexp_sim2 + 0.5) / (n_controls_exp_sim2 + 0.5) * (n_controls_nexp_sim2 + 0.5)
    # select the ones with some 0 but non-negative
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    # join both previous vectors
    n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)
    n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)
    n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)

    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0

    var_sim <- ifelse(some_zero,
      1 / ((n_cases + 1) - (n_cases_nexp_sim + 0.5)) + 1 / ((n_controls + 1) - (n_controls_nexp_sim + 0.5)) + 1 / (n_cases_nexp_sim + 0.5) + 1 / (n_controls_nexp_sim + 0.5),
      1 / (n_cases - n_cases_nexp_sim) + 1 / (n_controls - n_controls_nexp_sim) + 1 / n_cases_nexp_sim + 1 / n_controls_nexp_sim
    )

    # var_sim2 = 1 / ((n_cases+1) - (n_cases_nexp_sim+0.5)) + 1 / ((n_controls+1) - (n_controls_nexp_sim+0.5)) + 1 / (n_cases_nexp_sim+0.5) + 1 / (n_controls_nexp_sim+0.5)

    best <- order((var_sim - var)^2)[1]

    res$n_cases_nexp <- n_cases_nexp_sim[best]
    res$n_controls_nexp <- n_controls_nexp_sim[best]
    res$n_cases_exp <- n_cases - res$n_cases_nexp
    res$n_controls_exp <- n_controls - res$n_controls_nexp
  }

  return(res)
}

#' Estimate the n, using the variance, the number of exposed and non-exposed subjects
#'
#' @param or OR
#' @param var variance
#' @param n_exp number of exposed participants
#' @param n_nexp number of non exposed participants
#'
#' @noRd
.estimate_n_from_or_and_n_exp <- function(or, var, n_exp, n_nexp,
                                          n_cases = NA, n_controls = NA,
                                          baseline_risk = NA) {
  # ---- Rung 1: the case margin closes all four margins -> exact solve, no var used.
  n_cases_use <- if (!is.na(n_cases)) {
    n_cases
  } else if (!is.na(n_controls) && !is.na(n_exp) && !is.na(n_nexp)) {
    n_exp + n_nexp - n_controls
  } else {
    NA
  }
  if (!is.na(n_cases_use)) {
    hit <- .solve_2x2_from_or(or, n_exp, n_nexp, n_cases_use)
    # "solved" marks an exactly identified table -- see .estimate_n_from_or_and_n_cases.
    if (!is.null(hit)) return(structure(hit, solved = TRUE))
  }

  # ---- Rung 2: a reported baseline risk identifies the table just as well.
  if (!is.na(baseline_risk)) {
    hit <- .solve_2x2_from_or(or, n_exp, n_nexp,
                              .n_cases_from_baseline_risk(or, n_exp, n_nexp, baseline_risk))
    if (!is.null(hit)) return(structure(hit, solved = TRUE))
  }

  # ---- Otherwise: fall through to the enumeration below, with no prior applied.
  # Where neither rung fires the row is genuinely non-identified, and every candidate
  # rule tested amounted to a bet on outcome coding. "Assume events are the minority",
  # for instance, is 3.65x worse than doing nothing on common outcomes, wrong on 74% of
  # those rows, and produced a +44% pooled-RR bias end-to-end on a common-outcome
  # review. Recoding an outcome from "response" to "non-response" flips its answer on
  # identical data, so there is no principled default.
  #
  # The caller does refuse to publish a risk ratio built on this path when the two
  # supplied margins are equal, because the search is then not merely imprecise but
  # exactly tied. See .rotation_tied.
  res <- data.frame(n_cases_exp = NA, n_cases_nexp = NA, n_controls_exp = NA, n_controls_nexp = NA)
  attr(res, "solved") <- FALSE

  if (!is.na(or) & !is.na(var) & !is.na(n_exp) & !is.na(n_nexp)) {
    # first: uncorrected values with 0
    n_controls_exp_sim1 <- 0:n_exp
    n_controls_nexp_sim1 <- round(n_nexp / (1 + (n_exp - n_controls_exp_sim1) / (or * n_controls_exp_sim1)))
    n_cases_exp_sim1 <- n_exp - n_controls_exp_sim1
    n_cases_nexp_sim1 <- n_nexp - n_controls_nexp_sim1
    # we take the ones without 0 and non-negative
    idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # correcting by 0.5
    n_controls_exp_sim2 <- 0:n_exp
    n_controls_nexp_sim2 <- round((n_nexp + 0.5) - ((n_nexp + 1) * (n_exp - n_controls_exp_sim2 + 0.5)) / ((n_controls_exp_sim2 + 0.5) * or + n_exp - n_controls_exp_sim2 + 0.5))
    n_cases_exp_sim2 <- n_exp - n_controls_exp_sim2
    n_cases_nexp_sim2 <- n_nexp - n_controls_nexp_sim2

    # SELECT THE ONES THAT HAS SOME 0 BUT NO NEGATIVE ONES
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)
    n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)
    n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)


    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
    var_sim <- ifelse(some_zero,
      1 / ((n_exp + 1) - (n_controls_exp_sim + 0.5)) + 1 / (n_controls_exp_sim + 0.5) + 1 / ((n_nexp + 1) - (n_controls_nexp_sim + 0.5)) + 1 / (n_controls_nexp_sim + 0.5),
      1 / (n_exp - n_controls_exp_sim) + 1 / n_controls_exp_sim + 1 / (n_nexp - n_controls_nexp_sim) + 1 / n_controls_nexp_sim
    )

    # var_sim = 1 / (n_exp - n_controls_exp_sim) + 1 / n_controls_exp_sim + 1 / (n_nexp - n_controls_nexp_sim) + 1 / n_controls_nexp_sim
    best <- order((var_sim - var)^2)[1]
    res$n_controls_exp <- n_controls_exp_sim[best]
    res$n_controls_nexp <- n_controls_nexp_sim[best]
    res$n_cases_exp <- n_exp - res$n_controls_exp
    res$n_cases_nexp <- n_nexp - res$n_controls_nexp
  }
  return(res)
}

################### SE of OR ###################
.se_from_or <- function(x) {
  or <- as.numeric(x[1])
  n_cases <- as.numeric(x[2])
  n_controls <- as.numeric(x[3])

  res <- data.frame(value = NA, var = NA, se = NA)

  if (!is.na(or) & !is.na(n_cases) & !is.na(n_controls)) {
    cases_exp <- 1:(n_cases - 1)
    cases_nexp <- n_cases - cases_exp
    controls_exp <- round(n_controls / (1 + cases_nexp * or / cases_exp))
    controls_exp[which(controls_exp < 1 | controls_exp > n_controls - 1)] <- NA
    controls_nexp <- n_controls - controls_exp
    v_or_mean <- mean(1 / cases_exp + 1 / cases_nexp + 1 / controls_exp + 1 / controls_nexp, na.rm = TRUE)

    # An empty enumeration (every candidate NA'd by the feasibility filter above)
    # makes mean(numeric(0)) return NaN, not NA. .positive_or_na() cannot catch
    # that downstream (is.na(NaN) is TRUE, so its x <= 0 branch is skipped), so
    # the NaN would reach logor_se / logor_ci_lo / logor_ci_up beside a finite
    # logor. Decline the variance here instead; the OR itself is the user's own
    # input and is kept with se = NA (the bare-omega convention).
    if (!is.finite(v_or_mean)) return(res)

    res$value <- or
    res$var <- v_or_mean
    res$se <- sqrt(v_or_mean)
  }

  return(res)
}

################### RR to OR ###################
.rr_to_or <- function(rr, logrr_se, rr_ci_lo, rr_ci_up,
                      n_cases, n_controls, n_exp, n_nexp,
                      baseline_risk, rr_to_or) {
  if (rr_to_or == "grant") {
    # Grant's transform maps a risk ratio to an odds ratio as
    #   OR = RR (1 - BR) / (1 - RR * BR)
    # which is defined only while RR * BR < 1. It is applied to three values, the point
    # estimate and both interval bounds, and the upper bound is the largest of them, so
    # it leaves the domain first. The failure is therefore asymmetric and quiet: a row
    # can return a perfectly plausible log OR beside a NaN standard error and a
    # half-open interval, because every log() is wrapped in suppressWarnings(). A finite
    # estimate with no usable variance is worse than no estimate at all, since it looks
    # poolable and is not.
    #
    # The caller (es_from_stand_RR.R:153) routes rows here only when rr, baseline_risk
    # and both CI bounds are non-NA, so a non-finite result is always a domain violation
    # and never a missing input. Return the whole quartet as NA and say why.
    #
    # The mirror direction needs no such guard: the grant branch of .or_to_rr() computes
    # or / (1 - BR + BR*or), whose denominator is positive for every or > 0 and every BR
    # in (0, 1), so it has no domain boundary to cross.
    .grant_or <- function(x) {
      suppressWarnings(log(x * (1 - baseline_risk) / (1 - x * baseline_risk)))
    }
    logor_grant <- .grant_or(rr)
    logor_ci_lo_grant <- .grant_or(rr_ci_lo)
    logor_ci_up_grant <- .grant_or(rr_ci_up)
    logor_se <- (logor_ci_up_grant - logor_ci_lo_grant) / (2 * qnorm(.975))

    if (!all(is.finite(c(logor_grant, logor_se,
                         logor_ci_lo_grant, logor_ci_up_grant)))) {
      .offending <- c(rr = rr, rr_ci_lo = rr_ci_lo, rr_ci_up = rr_ci_up)
      .offending <- .offending[!is.finite(.grant_or(.offending))]
      warning("rr_to_or = 'grant': the transform is undefined because ",
              if (length(.offending)) {
                paste0(paste(sprintf("%s = %s", names(.offending),
                                     format(.offending)), collapse = " and "),
                       " times baseline_risk = ", format(baseline_risk),
                       " is not below 1")
              } else {
                paste0("baseline_risk = ", format(baseline_risk),
                       " does not admit a finite odds ratio for this interval")
              },
              ". Grant's conversion requires RR * baseline_risk < 1 for the point ",
              "estimate AND both confidence limits; the upper limit is the largest ",
              "of the three and so fails first. The odds ratio, its standard error ",
              "and both limits are returned as NA rather than a finite estimate with ",
              "a missing variance. Use rr_to_or = 'metaumbrella' or 'transpose', ",
              "which have no such restriction, or supply a smaller baseline_risk if ",
              "the one given is not the non-exposed event rate.",
              call. = FALSE)
      logor_grant <- NA_real_
      logor_se <- NA_real_
      logor_ci_lo_grant <- NA_real_
      logor_ci_up_grant <- NA_real_
    }

    res <- cbind(
      logor = logor_grant,
      logor_se = logor_se,
      logor_ci_lo = logor_ci_lo_grant,
      logor_ci_up = logor_ci_up_grant
    )
    return(res)
  } else if (rr_to_or == "metaumbrella") {
    raw_res <- .metaumbrella_rr_se_to_or(rr = rr, logrr_se = logrr_se, n_cases = n_cases, n_controls = n_controls)
    res <- cbind(
      logor = log(raw_res$value),
      logor_se = raw_res$se,
      logor_ci_lo = log(raw_res$value) - qnorm(.975) * raw_res$se,
      logor_ci_up = log(raw_res$value) + qnorm(.975) * raw_res$se
    )
    return(res)
  } else if (rr_to_or == "transpose") {
    res <- cbind(
      logor = log(rr),
      logor_se = logrr_se,
      logor_ci_lo = log(rr_ci_lo),
      logor_ci_up = log(rr_ci_up)
    )
    return(res)
  } else if (rr_to_or == "dipietrantonj") {
    n_dec <- max(
      nchar(gsub("^.+[.]", "", rr)),
      nchar(gsub("^.+[.]", "", rr_ci_lo)),
      nchar(gsub("^.+[.]", "", rr_ci_up))
    )

    estim <- estimraw::estim_raw(
      es = rr, lb = rr_ci_lo, ub = rr_ci_up,
      m1 = n_exp, m2 = n_nexp, dec = n_dec, measure = "rr"
    )
    if (length(estim) != 4) {
      # See the matching guard in the or_to_rr branch above: estimraw can return
      # no usable candidate, in which case which.min() yields integer(0) and
      # estim[[integer(0)]] aborts the whole convert_df() run. Return NA instead.
      numb <- 1L
      usable <- .estimraw_usable(estim)
      if (!is.na(baseline_risk) && length(estim) >= 2L) {
        cand <- c(
          abs(estim[[1]]$c[1] / (estim[[1]]$c[1] + estim[[1]]$d[1]) - baseline_risk),
          abs(estim[[2]]$c[1] / (estim[[2]]$c[1] + estim[[2]]$d[1]) - baseline_risk)
        )
        if (any(is.finite(cand))) numb <- which.min(cand)
      }
      # See the matching messages in the or_to_rr branch above.
      if (!any(usable) || !isTRUE(usable[numb])) {
        warning("rr_to_or = 'dipietrantonj': the 2x2 reconstruction has no real ",
                "solution for rr = ", rr, " with the supplied confidence interval and ",
                "group sizes, so the odds ratio is returned as NA. This is expected ",
                "when the outcome is not rare. Use another 'rr_to_or' method, or ",
                "supply the 2x2 counts directly.")
      } else if (is.na(baseline_risk) && sum(usable) >= 2L) {
        warning("rr_to_or = 'dipietrantonj': ", sum(usable),
                " candidate 2x2 reconstructions are compatible with rr = ", rr,
                " and its confidence interval. 'baseline_risk' is missing, so the ",
                "first candidate was used. Supply 'baseline_risk' (proportion of cases ",
                "in the non-exposed group) to select the reconstruction matching the ",
                "reported event rate.")
      }

      if (length(estim) < numb || length(numb) != 1L) {
        return(cbind(logor = NA_real_, logor_se = NA_real_,
                     logor_ci_lo = NA_real_, logor_ci_up = NA_real_))
      }

      calc_dipie <- es_from_2x2(
        n_cases_exp = estim[[numb]]$a[1],
        n_controls_exp = estim[[numb]]$b[1],
        n_cases_nexp = estim[[numb]]$c[1],
        n_controls_nexp = estim[[numb]]$d[1]
      )
    } else {
      calc_dipie <- es_from_2x2(
        n_cases_exp = estim$a[1],
        n_controls_exp = estim$b[1],
        n_cases_nexp = estim$c[1],
        n_controls_nexp = estim$d[1]
      )
    }
    res <- cbind(
      logor = calc_dipie$logor,
      logor_se = calc_dipie$logor_se,
      logor_ci_lo = calc_dipie$logor_ci_lo,
      logor_ci_up = calc_dipie$logor_ci_up
    )

    return(res)
  } else {
    # Keep this list in step with the branches above: the accepted values are 'grant'
    # and 'dipietrantonj'. There is no 'grant_2x2' or 'grant_CI'.
    stop(paste0("'", rr_to_or, "' not in tolerated values for the 'rr_to_or' argument. ",
                "Possible inputs are: 'metaumbrella', 'transpose', 'grant' or ",
                "'dipietrantonj'."), call. = FALSE)
  }
}

.metaumbrella_rr_se_to_or <- function(rr, logrr, logrr_se, n_cases, n_controls) {
  es <- data.frame(value = NA, se = NA)

  if (!is.na(rr) & !is.na(logrr_se) & !is.na(n_cases) & !is.na(n_controls)) {
    # uncorrected
    n_cases_nexp_sim1 <- 0:n_cases
    n_controls_nexp_sim1 <- round(n_cases_nexp_sim1 * ((rr * (n_cases + n_controls)) / (n_cases + (rr - 1) * n_cases_nexp_sim1) - 1))
    n_cases_exp_sim1 <- n_cases - n_cases_nexp_sim1
    n_controls_exp_sim1 <- n_controls - n_controls_nexp_sim1

    # we take only positives (no-zero)
    idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # corregint 0.5
    n_cases_nexp_sim2 <- 0:n_cases
    n_controls_nexp_sim2 <- ((n_cases + n_controls - n_cases_nexp_sim2 + 1) - (n_cases + n_controls + 2) * (n_cases - n_cases_nexp_sim2 + 0.5) / ((n_cases_nexp_sim2 + 0.5) * rr + n_cases - n_cases_nexp_sim2 + 0.5))
    n_cases_exp_sim2 <- n_cases - n_cases_nexp_sim2
    n_controls_exp_sim2 <- n_controls - n_controls_nexp_sim2

    # we take the ones with some 0 but non negative
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    #
    n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)
    n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)
    n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)

    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
    # var(log RR) = 1/a - 1/n1 + 1/c - 1/n2, where the two arm-total terms are
    # subtracted (delta method for the log of a binomial proportion; the same
    # expression as es_from_2x2()'s se_rr and the Cochrane Handbook). Note that this
    # differs from the log-OR pattern 1/a+1/b+1/c+1/d, in which all four terms are
    # added. Adding the arm-total terms here instead selects the wrong table along the
    # RR-constrained family and biases the reconstructed OR, so a round trip through a
    # known OR fails to recover it. The sweep variables map as
    # a = n_cases - n_cases_nexp_sim, c = n_cases_nexp_sim,
    # n1 (exposed total) = n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim),
    # n2 (non-exposed total) = n_cases_nexp_sim + n_controls_nexp_sim.
    var_sim <- ifelse(some_zero,
      1 / ((n_cases + 1) - (n_cases_nexp_sim + 0.5)) - 1 / ((n_cases + 1) + (n_controls + 1) - ((n_cases_nexp_sim + 0.5) + (n_controls_nexp_sim + 0.5))) +
        1 / (n_cases_nexp_sim + 0.5) - 1 / ((n_cases_nexp_sim + 0.5) + (n_controls_nexp_sim + 0.5)),
      1 / (n_cases - n_cases_nexp_sim) - 1 / (n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim)) +
        1 / n_cases_nexp_sim - 1 / (n_cases_nexp_sim + n_controls_nexp_sim)
    )

    best <- order((var_sim - logrr_se^2)^2)[1]
    n_cases_nexp <- n_cases_nexp_sim[best]
    n_controls_nexp <- n_controls_nexp_sim[best]
    n_cases_exp <- n_cases - n_cases_nexp
    n_controls_exp <- n_controls - n_controls_nexp

    # es$n_cases_nexp = n_cases_nexp
    # es$n_controls_nexp = n_controls_nexp
    # es$n_cases_exp = n_cases_exp
    # es$n_controls_exp = n_controls_exp

    cont_table <- es_from_2x2(
      n_cases_exp = n_cases_exp, n_controls_exp = n_controls_exp,
      n_cases_nexp = n_cases_nexp, n_controls_nexp = n_controls_nexp
    )

    es$value <- exp(cont_table$logor)
    es$se <- cont_table$logor_se
  }
  return(es)
}



################# 2x2 to R/Z ###################
.contingency_to_cor <- function(n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp,
                               table_2x2_to_cor, reverse_2x2) {

  # There is no "lipsey" (phi) branch here, and its absence is intentional.
  #
  # Phi is not a poolable estimand. Its attainable range is bounded by the margins, so
  # studies of the same association with different event rates report different phi
  # values, and pooling them manufactures heterogeneity that is pure margin artefact:
  # I^2 rises from 14% to 88% as the primary studies get larger, because tau^2 is pinned
  # by the margins while the within-study variance falls as 1/n.
  #
  # Nor could the choice be offered responsibly. A 2x2 table with fixed n has three free
  # parameters, and the dichotomised-bivariate-normal family also has three, so the
  # latent-normal model is saturated and no goodness-of-fit test can tell a user which
  # estimand applies to their data. See ?convert_df. Genuinely dichotomous variables
  # should use a binary measure (logor / rr / rd) rather than a correlation.
  if (table_2x2_to_cor == "tetrachoric") {
    res <- .tet_r(as.numeric(n_cases_exp),
                  as.numeric(n_controls_exp),
                  as.numeric(n_cases_nexp),
                  as.numeric(n_controls_nexp))
    res[res == "calculation failure"] <- NA

    # res[] <- lapply(res, function(x) as.numeric(as.character(x)))
    # On reverse, each CI must be negated and swapped (new_lo = -old_up,
    # new_up = -old_lo). Negating the bounds in place leaves lo > up, an inverted
    # interval that does not bracket the negated point estimate. The OR path handles the
    # same case in es_from_stand_OR.R, in the .or_to_cor result block. The z interval is
    # a symmetric Wald interval (z +- z*sqrt(vz)) and the r interval is its tanh
    # back-transform, so both are odd-symmetric about 0 because tanh(-x) = -tanh(x). The
    # reflected interval is therefore exactly the interval one would recompute around
    # the negated estimate. Old bounds are saved first, since res is overwritten.
    r_lo_raw <- res[3]; r_up_raw <- res[4]
    z_lo_raw <- res[7]; z_up_raw <- res[8]
    res[1] <- ifelse(reverse_2x2, -res[1], res[1])
    res[3] <- ifelse(reverse_2x2, -r_up_raw, r_lo_raw)
    res[4] <- ifelse(reverse_2x2, -r_lo_raw, r_up_raw)
    res[5] <- ifelse(reverse_2x2, -res[5], res[5])
    res[7] <- ifelse(reverse_2x2, -z_up_raw, z_lo_raw)
    res[8] <- ifelse(reverse_2x2, -z_lo_raw, z_up_raw)

    return(res)
  } else {
    # The previous message advertised 'cooper_delta', 'cooper_std' and 'lipsey'. None of
    # the three was implemented by any reachable code, so it named three methods the user
    # could not select.
    stop(paste0("'", table_2x2_to_cor, "' not in tolerated values for the ",
                "'table_2x2_to_cor' argument. The only possible input is 'tetrachoric'."),
         call. = FALSE)
  }
}

# Session-scoped record of which one-time notices have already been emitted, so an
# installation problem is reported once rather than once per row.
.mcv_notices <- new.env(parent = emptyenv())

# Wrapped rather than called inline so the missing-package branch is reachable in tests.
.has_mvtnorm <- function() requireNamespace("mvtnorm", quietly = TRUE)

# message(), not warning(): the tetrachoric call sites wrap this route in
# suppressWarnings() because a genuine numerical failure is expected on some tables and
# must not spam the console. That suppression would swallow a warning raised here, so the
# installation notice is signalled as a message, which survives it.
.notify_once <- function(key, ...) {
  if (isTRUE(get0(key, envir = .mcv_notices, ifnotfound = FALSE))) {
    return(invisible(NULL))
  }
  assign(key, TRUE, envir = .mcv_notices)
  message(...)
  invisible(NULL)
}

.tet_r <- function(n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp) {
  # metafor::escalc(measure = "RTET") solves the tetrachoric correlation by numerical ML
  # over a bivariate normal CDF, which needs 'mvtnorm'. mvtnorm is a Suggests of metafor
  # rather than an Imports, so installing metafor does not bring it in and a perfectly
  # ordinary installation can lack it. Without this guard escalc() raises an error, the
  # tryCatch below swallows it, and every 2x2-derived correlation and Fisher's z comes
  # back NA with nothing said, which looks exactly like data that genuinely cannot
  # support them.
  if (!.has_mvtnorm()) {
    .notify_once(
      "mvtnorm_missing",
      "The tetrachoric correlation estimated from a 2x2 table requires the 'mvtnorm' ",
      "package, which is not installed. The correlation and Fisher's z derived from a ",
      "2x2 table (and from phi, chi-squared and proportions, which are routed through ",
      "it) are therefore returned as NA. Run install.packages(\"mvtnorm\") to enable ",
      "them. All other effect size measures are unaffected.")
    return(cbind(NA, NA, NA, NA, NA, NA, NA, NA))
  }
  tryCatch(
    expr = {
      tet <- metafor::escalc(
        ai = n_cases_exp, bi = n_controls_exp,
        ci = n_cases_nexp, di = n_controls_nexp,
        measure = "RTET"
      )
      n_sample = n_cases_exp + n_controls_exp + n_cases_nexp + n_controls_nexp
      r <- tet$yi
      vr <- tet$vi
      z <- atanh(r)
      vz <- vr / ((1 - r^2)^2)

      z_lo <- z - qnorm(.975) * sqrt(vz)
      z_up <- z + qnorm(.975) * sqrt(vz)
      # The r-scale interval is the back-transformed z interval, not a symmetric Wald
      # interval on r. A symmetric interval r +- z*sqrt(vr) is unbounded and routinely
      # escapes the parameter space: over a realistic grid it left [-1, 1] on 45.8% of
      # tables, and flag B1b (the [INFO] disclosure that exists for this reason) fired
      # on 20.6%, rising to 67% at n = 50 with |r| >= 0.6. Since vz is already the
      # delta-method transform of vr on the line above, tanh() of the z bounds gives the
      # same interval on a scale where it cannot overshoot. Escape falls to 0.0%, and
      # coverage of the true correlation is unchanged (96.1% -> 96.2%).
      r_lo <- tanh(z_lo)
      r_up <- tanh(z_up)

      # BOUNDARY. metafor's ML returns the boundary estimate r = +-1 whenever a cell is
      # empty, together with an enormous vi -- that pair is how escalc() says the
      # correlation is not identified by this table, and both halves of it are honest
      # and must be kept. What is not honest is everything the transform produces from
      # it: z = atanh(+-1) = +-Inf, vz = vr / (1 - r^2)^2 = Inf, one CI bound is then
      # Inf - Inf = NaN, and tanh() carries the NaN onto the r scale. Measured on the
      # raw table (0, 20, 10, 10): r = -1, vr = 22209.7 (both wanted), z = -Inf,
      # vz = Inf, z_up = NaN, r_up = NaN. Those four are exported columns, and an Inf
      # standard error is not a large standard error -- it is a value no downstream
      # check, weight or plot can consume.
      #
      # So the point estimate and its variance pass through and the transformed
      # quantities are declined. Keyed on non-finite vz rather than on |r| >= 1 alone,
      # so a non-finite vi from any other cause is caught by the same line.
      degenerate <- !is.finite(r) | abs(r) >= 1 | !is.finite(vr) | !is.finite(vz)
      z[degenerate] <- NA_real_
      vz[degenerate] <- NA_real_
      z_lo[degenerate] <- NA_real_
      z_up[degenerate] <- NA_real_
      r_lo[degenerate] <- NA_real_
      r_up[degenerate] <- NA_real_
      # A non-finite r or vr is not an estimate at all, unlike the boundary pair above.
      #
      # vr == 0 is the third boundary shape and it breaks the premise of the paragraph
      # above. On a PERFECT-ASSOCIATION table -- both off-diagonal cells zero -- the
      # boundary r does NOT arrive with an enormous vi: escalc(measure = "RTET")
      # returns vi = 0 exactly. Enumerating every table with cells 0..6 and N >= 4
      # (2366 tables), vi == 0 occurs on exactly the 66 perfect-association tables and
      # nowhere else. Keeping that pair would export r = +-1 with r_se = 0, i.e. an
      # infinite inverse-variance weight, and rma() aborts on it -- which is precisely
      # what the risk-difference guard in es_from_2x2() already refuses for the same
      # tables. A zero variance is not a precise estimate, it is the same statement
      # "not identified" that the enormous-vi shape makes at the other extreme, so it
      # is declined the same way.
      # The point estimate is kept, matching the convention used for a reported omega
      # with no variance source: es present, se = NA, so the row stays visible and
      # countable while summary() leaves it out of the pool.
      r[!is.finite(r) | !is.finite(vr)] <- NA_real_
      vr[!is.finite(vr) | (!is.na(vr) & vr <= 0)] <- NA_real_

      dat <- cbind(r, vr, r_lo, r_up, z, vz, z_lo, z_up)
      return(dat)
    },
    error = function(e) {
      dat <- cbind(NA, NA, NA, NA, NA, NA, NA, NA)
      return(dat)
    },
    warning = function(w) {
      dat <- cbind(NA, NA, NA, NA, NA, NA, NA, NA)
      return(dat)
    },
    finally = {
    }
  )
}

################# OR to R/Z ###################
.or_to_cor <- function(or, logor_se,
                       n_cases,
                       n_exp,
                       small_margin_prop,
                       n_sample,
                       or_to_cor) {
  if (or_to_cor %in% c("bonett", "pearson")) {
    if (or_to_cor == "bonett") {
      # small_margin_prop is Bonett and Price's p_min, the SMALLEST of the four marginal
      # proportions of the underlying 2x2 table, so its domain is (0, 0.5]. Their
      # coefficient c is guaranteed to land in (0.275, 0.5] only on that domain; outside
      # it the guarantee is lost and, once c goes negative, r = cos(pi / (1 + or^c))
      # reverses its monotonicity in or while r_se picks up the sign of c. A percentage
      # entered instead of a proportion (30 for 0.30) gives c = -434.6, r = -1 for an
      # odds ratio of 3, and r_se = -1.4e-221 -- a NEGATIVE standard error, i.e. an
      # inverse-variance weight of ~5e442 that dominates any pool it enters. Entering
      # the LARGEST margin (0.9) instead of the smallest is the quiet case: r is ~12%
      # off with nothing at all to see.
      #
      # .bounded_columns() in R/internal_flags.R carries the same (0, 0.5] domain, so a
      # convert_df() run names the cell in a Tier-1 flag and, under correct_inputs =
      # TRUE, blanks it. This guard is the route-level half of the two-layer policy
      # written out at the top of R/internal_guards.R: es_from_or_se()/_ci()/_pval() are
      # exported, documented and called directly by metaumbrella, none of which passes
      # through the column-keyed Tier-1 validation, and under correct_inputs = FALSE the
      # value is preserved on purpose. The flag reports; the guard contains.
      small_margin_prop <- ifelse(
        !is.na(small_margin_prop) & (small_margin_prop <= 0 | small_margin_prop > 0.5),
        NA_real_, small_margin_prop
      )
      c <- (1 - abs(n_exp/n_sample - n_cases/n_sample) / 5 - (1 / 2 - small_margin_prop)^2) / 2
      # Belt and braces. With p_min inside its domain c cannot exceed 0.5, and it can
      # only reach 0 through margins so incoherent that n_exp/n_sample and
      # n_cases/n_sample differ by more than 3.75 -- which no bound on p_min can rule
      # out, since those two proportions are read from the margins directly. A
      # non-positive c is neutralised rather than clamped: a clamped c would return a
      # number that is not the Bonett estimate of anything.
      c <- .positive_or_na(c)
    } else {
      c <- 1 / 2
    }

    r <- cos(pi / (1 + or^c))
    r_se <- logor_se * (pi * c * or^c) * sin(pi / (1 + or^c)) / (1 + or^c)^2

    or_ci_lo <- exp(log(or) - qnorm(.975) * logor_se)
    or_ci_up <- exp(log(or) + qnorm(.975) * logor_se)
    r_lo <- cos(pi / (1 + or_ci_lo^c))
    r_up <- cos(pi / (1 + or_ci_up^c))

    z <- atanh(r)
    z_se <- sqrt(r_se^2 / ((1 - r^2)^2)) # delta method
    z_lo <- atanh(r_lo)
    z_up <- atanh(r_up)

    res <- cbind(r, r_se, r_lo, r_up, z, z_se, z_lo, z_up)
    return(res)
  } else if (or_to_cor == "digby") {
    c <- 3 / 4

    r <- (or^c - 1) / (or^c + 1)
    r_se <- sqrt((c^2 / 4) * (1 - r^2)^2 * logor_se^2)

    z <- atanh(r)
    z_se <- sqrt(r_se^2 / ((1 - r^2)^2)) # delta method
    z_lo <- z - qnorm(.975) * sqrt(c^2 / 4 * logor_se^2)
    z_up <- z + qnorm(.975) * sqrt(c^2 / 4 * logor_se^2)

    r_lo <- tanh(z_lo)
    r_up <- tanh(z_up)

    res <- cbind(r, r_se, r_lo, r_up, z, z_se, z_lo, z_up)
    return(res)
  }
}

################# SMD to R/Z ###################
.smd_to_cor <- function(d, vd, n_exp, n_nexp, smd_to_cor, n_cov_ancova) {
  # ---------------------------------------------------------------------------
  # Precision carried by the supplied d, relative to the crude two-group design.
  #
  # Both branches below convert d -> r through a deterministic map, so Var(r) must
  # inherit Var(d). The viechtbauer branch, however, uses Soper's (1914) large-sample
  # closed form for the biserial correlation, which depends only on n, p and r. That
  # closed form is not an independent estimator: to leading order it is the delta
  # propagation of the crude d variance, agreeing with it to about 0.3% at n = 400 and
  # exactly up to (n-1)/(n-2) at r = 0. Using it verbatim would therefore discard
  # whatever precision the d actually carries, whether that is the Cooper eq. 12.26
  # (1 - R^2) shrink on an ANCOVA row, the 2(1 - r_pre_post) factor on a pre-post row,
  # the control-group df on a Glass row, or a user-reported standard error. Rescaling by
  # vd / vd_crude restores it and leaves the crude case untouched, since the ratio is
  # then exactly 1 (vd_crude below is the same expression .es_from_d() uses). That
  # preserves agreement with metafor's measure = "RBIS" on the rows it applies to.
  #
  # HISTORICAL NOTE (fixed upstream, in .es_from_d()): vd_crude below is hard-coded to
  # the LS2 form, and .es_from_d() used to hand over whatever default variance smd_var
  # had selected. Under the opt-in smd_var = "hedges_olkin" (alias "viechtbauer",
  # metafor's LS) an ordinary crude row therefore arrived with
  # vd = (leading + (d*J)^2/(2N))/J^2 and got prec_ratio = 1.0257 at n = 30/30, d = 0.5,
  # inflating r_se^2 and z_se^2 by a pure estimator-convention factor and taking z_se^2
  # off the stabilised 1/(N - 1). .es_from_d() now distinguishes a SUPPLIED d_se (passed
  # through untouched, so a user-reported SE never depends on smd_var) from a DEFAULTED
  # one (handed over as vd_ls2, the LS2 expression) -- see the note above dat_r in
  # R/internal_es_from_d.R. So the "exactly 1" above now holds under either setting:
  # measured at n = 30/30, d = 0.5, r_se and z_se^2 are identical across smd_var,
  # z_se^2 = 1/(N - 1) = 0.01694915, and r_se^2 = 0.02273584 matches
  # metafor::escalc(measure = "RBIS") to 15 digits for both.
  vd_crude <- (n_exp + n_nexp) / (n_exp * n_nexp) + d^2 / (2 * (n_exp + n_nexp))
  prec_ratio <- ifelse(!is.na(vd) & is.finite(vd) & !is.na(vd_crude) & vd_crude > 0,
                       vd / vd_crude, 1)

  if (smd_to_cor == "viechtbauer") {
    # h encodes the finite-sample point-biserial identity r_pb = t / sqrt(t^2 + df),
    # which holds for the marginal two-group design at df = n_exp + n_nexp - 2. The d
    # arriving here is always on the marginal (unadjusted) SD scale, which is the point
    # of Cooper's eq. 12.23/12.24 convention: it keeps ANCOVA studies poolable with
    # unadjusted ones. Subtracting n_cov_ancova from df would shrink h and inflate r
    # into a quantity that is neither the marginal r_pb nor the partial one (the partial
    # would also require h * (1 - cov_outcome_r^2)), and would make the reported
    # correlation depend on how many covariates the source study happened to adjust for.
    # n_cov_ancova therefore enters only the confidence-interval degrees of freedom.
    df <- n_exp + n_nexp - 2
    h <- df / n_exp + df / n_nexp
    p <- n_exp / (n_exp + n_nexp)
    q <- n_nexp / (n_exp + n_nexp)
    r_pb <- d / sqrt(d^2 + h)

    f <- dnorm(qnorm(p, lower.tail = FALSE))
    r_viechtbauer <- sqrt(p * q) / f * r_pb
    r_trunc = ifelse(r_viechtbauer > 1, 1, ifelse(r_viechtbauer < -1, -1, r_viechtbauer))
    vr_viechtbauer <- 1 / (n_exp + n_nexp - 1) *
      (p * q / f^2 - (3 / 2 + (1 - p * qnorm(p, lower.tail = FALSE) / f) *
                              (1 + q * qnorm(p, lower.tail = FALSE) / f)) *
         r_trunc^2 + r_trunc^4)


    # ========= z ========= #
    fzp <- dnorm(qnorm(p))
    a_viechtbauer <- sqrt(fzp) / (p * (1 - p))^(1 / 4)
    z_viechtbauer <- (a_viechtbauer / 2) * log((1 + a_viechtbauer * r_trunc) /
                                               (1 - a_viechtbauer * r_trunc))
    vz_viechtbauer <- 1 / (n_exp + n_nexp - 1)

    # Carry the actual precision of d through to r and z (see the note at the top of
    # this function). Both are scaled by the same factor, so the variance-stabilising
    # relation between them is preserved exactly whatever the design: z is built so that
    # (dz/dr)^2 * vr = 1/(n-1) under the standard design.
    vr_viechtbauer <- vr_viechtbauer * prec_ratio
    vz_viechtbauer <- vz_viechtbauer * prec_ratio
    # ========= 95% CI ===== #


    z_lo_viechtbauer <- z_viechtbauer - qnorm(.975) * sqrt(vz_viechtbauer)
    z_up_viechtbauer <- z_viechtbauer + qnorm(.975) * sqrt(vz_viechtbauer)
    r_lo_viechtbauer <- (1 / a_viechtbauer) * ((exp(2 * z_lo_viechtbauer / a_viechtbauer) - 1) / (exp(2 * z_lo_viechtbauer / a_viechtbauer) + 1))
    r_up_viechtbauer <- (1 / a_viechtbauer) * ((exp(2 * z_up_viechtbauer / a_viechtbauer) - 1) / (exp(2 * z_up_viechtbauer / a_viechtbauer) + 1))

    # ---- Out-of-range biserial r: keep the point estimate, decline what is derived.
    #
    # sqrt(p*q)/f is 1.2533 at p = 0.5, so r_viechtbauer is NOT bounded by 1: any
    # point-biserial above 0.798 pushes it out of the parameter space. r_trunc exists so
    # that Soper's variance and the variance-stabilising z stay computable, but once it
    # binds, z collapses onto (a/2)*log((1+a)/(1-a)) -- a function of the arm-size split
    # ALONE, carrying no information about the data, and it arrives with a plausible SE
    # and a plausible CI. On the r scale Category B catches this ("[INVALID] r outside
    # [-1, 1]"); z is unbounded, so no bound check can see it there, which left the
    # POOLABLE metric as the one where the error was invisible. On the shipped df.short
    # that is Lopez_2019 (route cohen_d) and Yu_2018 (route med_min_max) -- two
    # unrelated studies from two different routes, both reported as the constant.
    #
    # Same policy as es_disattenuate() for the same situation (R/es_disattenuate.R:147):
    # keep the point estimate and its first-order SE for inspection, NA every derived
    # quantity, since a clamped value is an artifact identical whatever the input.
    saturated <- !is.na(r_viechtbauer) & abs(r_viechtbauer) > 1

    # Say so, once per session. Declining is right for this estimand -- a biserial r
    # outside [-1, 1] means the latent-normal model behind it does not fit -- but a row
    # that leaves a measure = "z" pool with no explanation is the part that could bias a
    # synthesis: the threshold is |d| ~ 2.6 at balanced arms (2.4 at 20/80), so what
    # drops out is always the upper tail, and dropping the upper tail pulls the pooled
    # estimate down. The remedy loses nothing, which is why it is named here rather than
    # left to the help page: smd_to_cor = "lipsey_cooper" uses d/sqrt(d^2 + a), which is
    # bounded by construction, so every study stays in and the affected rows come back
    # as distinct values rather than the shared constant the old clamp produced.
    if (any(saturated)) {
      .notify_once(
        "viechtbauer_biserial_out_of_range",
        # No count: .smd_to_cor() is called ONE ROW AT A TIME (via mapply in
        # .es_from_d(), and once per route in convert_df()), so sum(saturated) here is
        # the count within a single call -- always 1 -- not the number of affected
        # studies. It read "1 row(s)" on a dataset with 6 distinct saturated studies
        # across 9 routes. Counting honestly would mean accumulating across routes and
        # de-duplicating by study, which this scope cannot see; until then the notice
        # says only what it can support.
        "One or more rows have a biserial correlation that falls outside [-1, 1] ",
        "under smd_to_cor = \"viechtbauer\" (this happens above |d| of roughly 2.6). ",
        "Their Fisher's z is returned as NA rather than clamped, because a clamped value ",
        "is a function of the group-size split alone and is identical for any two such ",
        "studies. The correlation itself is kept and flagged. To keep these studies in a ",
        "z analysis, use smd_to_cor = \"lipsey_cooper\", whose correlation is bounded by ",
        "construction, so no row is dropped.")
    }

    r_lo_viechtbauer <- ifelse(saturated, NA_real_, r_lo_viechtbauer)
    r_up_viechtbauer <- ifelse(saturated, NA_real_, r_up_viechtbauer)
    z_viechtbauer    <- ifelse(saturated, NA_real_, z_viechtbauer)
    vz_viechtbauer   <- ifelse(saturated, NA_real_, vz_viechtbauer)
    z_lo_viechtbauer <- ifelse(saturated, NA_real_, z_lo_viechtbauer)
    z_up_viechtbauer <- ifelse(saturated, NA_real_, z_up_viechtbauer)

    res <- cbind(
      r_viechtbauer, vr_viechtbauer, r_lo_viechtbauer, r_up_viechtbauer,
      z_viechtbauer, vz_viechtbauer, z_lo_viechtbauer, z_up_viechtbauer
    )

    return(res)
  } else if (smd_to_cor == "lipsey_cooper") {
    a <- ((n_exp + n_nexp)^2) / (n_exp * n_nexp)
    p <- n_exp / (n_exp + n_nexp)
    r_lipsey <- d / sqrt(d^2 + 1 / (p * (1 - p)))
    vr_lipsey <- a^2 * vd / ((d^2 + a)^3)
    z_lipsey <- atanh(r_lipsey)
    # Delta method for z = atanh(r), r = d / sqrt(d^2 + a) with a = 1/(p*(1-p)):
    # dz/dd = 1/sqrt(d^2 + a), so Var(z) = vd / (d^2 + a). The additive denominator term
    # is the squared point estimate d^2, not the sampling variance vd. That makes
    # vz_lipsey the exact Fisher transform of vr_lipsey above (vz = vr/(1-r^2)^2).
    # esc::convert_d2r() returns vd/(vd + a) here instead, which is not Fisher-consistent
    # with its own r variance. metaConvert uses the consistent value, so the
    # lipsey_cooper z-SE differs from esc for large |d| and agrees with it as d -> 0.
    vz_lipsey <- vd / (d^2 + 1 / (p * (1 - p)))
    # Same error degrees of freedom as the d/g interval built for this row in
    # .es_from_d() (n_cov_ancova is 0 on every crude route, so this is a no-op there).
    df_ci_lipsey <- n_exp + n_nexp - 2 - n_cov_ancova
    r_lo_lipsey <- r_lipsey - qt(.975, df = df_ci_lipsey) * sqrt(vr_lipsey)
    r_up_lipsey <- r_lipsey + qt(.975, df = df_ci_lipsey) * sqrt(vr_lipsey)
    z_lo_lipsey <- z_lipsey - qnorm(.975) * sqrt(vz_lipsey)
    z_up_lipsey <- z_lipsey + qnorm(.975) * sqrt(vz_lipsey)

    res <- cbind(
      r_lipsey, vr_lipsey, r_lo_lipsey, r_up_lipsey,
      z_lipsey, vz_lipsey, z_lo_lipsey, z_up_lipsey
    )

    return(res)
  }
}

################# PHI to R/Z ###################
.phi_to_cor <- function(phi, n_sample,
                        n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp,
                        phi_to_cor, reverse_phi) {
  if (phi_to_cor == "lipsey") {
    r <- ifelse(reverse_phi, -phi, phi)
    z <- atanh(r)
    vz <- z^2 / (r^2 * n_sample)
    z_lo <- z - qnorm(.975) * sqrt(vz)
    z_up <- z + qnorm(.975) * sqrt(vz)
    effective_n = 1 / vz + 3
    vr = (1 - r^2)^2 / (effective_n - 1)
    # t <- qt(pnorm(z / sqrt(vz)), n_sample - 2)
    # vr <- (r / t)^2
    # r_lo = r - qnorm(.975)*sqrt(vr)
    # r_up = r + qnorm(.975)*sqrt(vr)
    r_lo <- tanh(z_lo)
    r_up <- tanh(z_up)
    res <- cbind(r, vr, r_lo, r_up, z, vz, z_lo, z_up)
    return(res)
  } else if (phi_to_cor == "tetrachoric") {
    # if (!requireNamespace("mvtnorm", quietly = TRUE) & !requireNamespace("metafor", quietly = TRUE)) {
    #   stop("Please install the 'mvtnorm' and 'metafor' packages to compute tetrachoric correlation.")
    # } else if (!requireNamespace("metafor", quietly = TRUE)) {
    #   stop("Please install the 'metafor' package to compute tetrachoric correlation.")
    # } else if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    #   stop("Please install the 'mvtnorm' package to compute tetrachoric correlation.")
    # }

    res <- .tet_r(
      as.numeric(n_cases_exp), as.numeric(n_controls_exp),
      as.numeric(n_cases_nexp), as.numeric(n_controls_nexp)
    )
    res[res == "calculation failure"] <- NA

    # Negate and swap each CI on reverse; see the matching block in
    # .contingency_to_cor() for why negating in place is wrong.
    r_lo_raw <- res[3]; r_up_raw <- res[4]
    z_lo_raw <- res[7]; z_up_raw <- res[8]
    res[1] <- ifelse(reverse_phi, -res[1], res[1])
    res[3] <- ifelse(reverse_phi, -r_up_raw, r_lo_raw)
    res[4] <- ifelse(reverse_phi, -r_lo_raw, r_up_raw)
    res[5] <- ifelse(reverse_phi, -res[5], res[5])
    res[7] <- ifelse(reverse_phi, -z_up_raw, z_lo_raw)
    res[8] <- ifelse(reverse_phi, -z_lo_raw, z_up_raw)

    return(res)
  } else {
    stop(paste0("'", phi_to_cor, "' not in tolerated values for the 'phi_to_cor' argument. Possible inputs are: 'tetrachoric', or 'lipsey'"))
  }
}

################# CHI to R/Z ###################
.chi_to_cor <- function(chisq, n_sample,
                        n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp,
                        chisq_to_cor, reverse_chisq) {
  if (chisq_to_cor == "lipsey") {
    r <- sqrt(chisq / n_sample)
    r <- ifelse(reverse_chisq, -r, r)
    z <- atanh(r)
    vz <- z^2 / (chisq)
    z_lo <- z - qnorm(.975) * sqrt(vz)
    z_up <- z + qnorm(.975) * sqrt(vz)
    # t <- qt(pnorm(z / sqrt(vz)), n_sample - 2)
    # vr <- (r / t)^2
    effective_n = 1 / vz + 3
    vr = (1 - r^2)^2 / (effective_n - 1)
    r_lo <- tanh(z_lo)
    r_up <- tanh(z_up)
    # r_lo = r - qt(.975, n_sample - 2)*sqrt(vr)
    # r_up = r + qt(.975, n_sample - 2)*sqrt(vr)
    # r_lo <- tanh(z_lo)
    # r_up <- tanh(z_up)
    res <- cbind(r, vr, r_lo, r_up, z, vz, z_lo, z_up)
    return(res)
  } else if (chisq_to_cor == "tetrachoric") {
    # if (!requireNamespace("mvtnorm", quietly = TRUE) & !requireNamespace("metafor", quietly = TRUE)) {
    #   stop("Please install the 'mvtnorm' and 'metafor' packages to compute tetrachoric correlation.")
    # } else if (!requireNamespace("metafor", quietly = TRUE)) {
    #   stop("Please install the 'metafor' package to compute tetrachoric correlation.")
    # } else if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    #   stop("Please install the 'mvtnorm' package to compute tetrachoric correlation.")
    # }

    res <- .tet_r(
      as.numeric(n_cases_exp), as.numeric(n_controls_exp),
      as.numeric(n_cases_nexp), as.numeric(n_controls_nexp)
    )
    res[res == "calculation failure"] <- NA

    # Negate and swap each CI on reverse; see the matching block in
    # .contingency_to_cor() for why negating in place is wrong.
    r_lo_raw <- res[3]; r_up_raw <- res[4]
    z_lo_raw <- res[7]; z_up_raw <- res[8]
    res[1] <- ifelse(reverse_chisq, -res[1], res[1])
    res[3] <- ifelse(reverse_chisq, -r_up_raw, r_lo_raw)
    res[4] <- ifelse(reverse_chisq, -r_lo_raw, r_up_raw)
    res[5] <- ifelse(reverse_chisq, -res[5], res[5])
    res[7] <- ifelse(reverse_chisq, -z_up_raw, z_lo_raw)
    res[8] <- ifelse(reverse_chisq, -z_lo_raw, z_up_raw)

    return(res)
  } else {
    stop(paste0("'", chisq_to_cor, "' not in tolerated values for the 'chisq_to_cor' argument. Possible inputs are: 'tetrachoric', or 'lipsey'"))
  }
}


################# PRE POST helpers ##############
#' Guard a standardizing SD before it is used as a denominator
#'
#' Returns NA for any standardizer that is zero, negative or non-finite, so the
#' resulting SMD is NA rather than an unnoticed Inf or NaN. The guard applies to
#' the denominator only: the mean-change wrappers legitimately pass
#' mean_pre_sd = 0, because the pre slot is zeroed by construction, so guarding
#' raw inputs would break them, whereas a zero standardizer is always degenerate.
#' This is reachable when a study reports SD = 0, or when r_pre_post = 1 with
#' equal pre/post SDs makes the change SD collapse to 0.
#'
#' @noRd
.guard_standardizer <- function(sd_value) {
  ifelse(is.finite(sd_value) & sd_value > 0, sd_value, NA_real_)
}

#' Guard a pre-post correlation
#'
#' |r| >= 1 is not a valid correlation and drives sd_change to 0 (or a negative
#' radicand) in every pre-post formula. Strictly at the bounds so that r = 0.99
#' still computes.
#'
#' @noRd
.guard_r_pre_post <- function(r) {
  ifelse(is.finite(r) & abs(r) < 1, r, NA_real_)
}

#' The pre-post cross term 2*r*SD_pre*SD_post, made zero-aware
#'
#' Var(post - pre) = SD_pre^2 + SD_post^2 - 2*r*SD_pre*SD_post. Several callers zero
#' one slot by construction: the mean-change wrappers pass mean_pre_sd = 0 (the change
#' SD arrives through mean_post_sd), and .paired_t_to_smd() passes mean_pre_sd = 0 with
#' mean_post_sd = 1. The cross term is then mathematically zero whatever r is -- with
#' Var(pre) = 0 the change SD is exactly SD_post -- but in R `0 * NA` is NA, so an r
#' that .guard_r_pre_post() has NA'd (|r| >= 1, e.g. a correlation typed as a
#' percentage) used to poison the whole radicand and null a d_z that carries no r at
#' all. Contribute an exact 0 on those rows instead.
#'
#' Only the cross term is neutralised: morris_drm still returns NA through its
#' sqrt(2(1-r)) factor and morris_dav through its (1 + r^2) effective df, which is
#' correct -- those estimands do use r.
#'
#' @noRd
.pre_post_cross <- function(r, sd_pre, sd_post) {
  ifelse(!is.na(sd_pre) & !is.na(sd_post) & (sd_pre == 0 | sd_post == 0),
         0, 2 * r * sd_pre * sd_post)
}


#' Convert a paired t statistic into a within-group SMD
#'
#' The paired-t and paired-F routes delegate to .single_group_pre_post_to_smd() rather
#' than keeping their own copy of the morris_dz and morris_drm arithmetic, so that all
#' pre-post routes share one implementation.
#'
#' A paired t carries no pre-test or post-test SD, so the kernel is given an equivalent
#' problem instead: a change score of t/sqrt(n) with a change SD of 1, and the pre-test
#' slot set to zero exactly as the mean-change wrappers do. The kernel's sd_diff then
#' reduces to 1, which yields d = t/sqrt(n) under morris_dz and d = t*sqrt(2(1-r)/n)
#' under morris_drm. The variances come from the same branch. That reduction holds for
#' ANY r, including an out-of-range one that .guard_r_pre_post() has NA'd: the kernel's
#' cross term is zero-aware (.pre_post_cross), so the r-free morris_dz row stays finite.
#'
#' The kernel is called row by row through mapply because it dispatches on a single
#' pre_post_to_smd value, whereas these routes accept one method per row. A plain
#' vectorised call would apply the first row's method to every row.
#'
#' @param paired_t vector of paired t values
#' @param n vector of per-arm sample sizes
#' @param r_pre_post vector of pre-post correlations
#' @param pre_post_to_smd scalar or vector of method names
#' @return a matrix with the kernel's eight columns, one row per input row
#' @noRd
.paired_t_to_smd <- function(paired_t, n, r_pre_post, pre_post_to_smd) {
  k <- max(length(paired_t), length(n), length(r_pre_post), length(pre_post_to_smd))
  rec <- function(x) rep_len(x, k)
  nn <- rec(n)

  # An arm with n < 2 has no estimable variance, and the kernel returns NA for it, but
  # it gets there through qt(.975, n - 1) and .d_j(n - 1), which emit "NaNs produced"
  # along the way. Such arms are skipped below, so the results are unchanged and the
  # warning stream stays clean. Wrapping the kernel in suppressWarnings() would instead
  # hide genuine warnings from the other rows.
  #
  # .d_j() is the source of the "Hedges' J correction is undefined for df <= 1"
  # warning, which is worth keeping: it tells the caller that g will be NA on those
  # rows. Because the skip means the kernel never reaches .d_j() for those rows, it is
  # called explicitly here on the per-arm df (n - 1). Apart from the warning the call
  # has no effect.
  invisible(.d_j(nn - 1))

  # A row with no paired statistic yields all-NA from every branch of the kernel, so
  # sending it through the per-row mapply buys nothing. convert_df() fills unsupplied
  # columns with NA and calls all five paired routes on every dataset, so a sheet
  # carrying no paired data at all used to pay nine full row-by-row loops purely to
  # produce NA: ~26% of the R time of a convert_df() run at 3000 rows. Skipping those
  # rows is output-identical (checked branch by branch) and leaves the .d_j() warning
  # stream above untouched, since that call sits outside the filter on purpose.
  ok <- is.finite(nn) & nn >= 2 & !is.na(rec(paired_t))
  out <- matrix(NA_real_, nrow = k, ncol = 8L)
  if (!any(ok)) {
    colnames(out) <- c("d", "var_d", "d_ci_lo", "d_ci_up",
                       "g", "var_g", "g_ci_lo", "g_ci_up")
    return(out)
  }
  res <- t(mapply(.single_group_pre_post_to_smd,
    mean_pre        = 0,
    mean_post       = (rec(paired_t) / sqrt(nn))[ok],
    mean_pre_sd     = 0,
    mean_post_sd    = 1,
    n               = nn[ok],
    r_pre_post      = rec(r_pre_post)[ok],
    pre_post_to_smd = rec(pre_post_to_smd)[ok]))
  if (is.null(dim(res))) res <- matrix(res, nrow = sum(ok), byrow = FALSE)
  out[ok, ] <- res
  storage.mode(out) <- "double"
  colnames(out) <- c("d", "var_d", "d_ci_lo", "d_ci_up",
                     "g", "var_g", "g_ci_lo", "g_ci_up")
  out
}

################# PRE POST to SMD ##############
.pre_post_to_smd <- function(mean_pre_exp, mean_pre_sd_exp,
                             mean_exp, mean_sd_exp,
                             mean_pre_nexp, mean_pre_sd_nexp,
                             mean_nexp, mean_sd_nexp,
                             n_exp, n_nexp,
                             r_pre_post_exp, r_pre_post_nexp,
                             pre_post_to_smd,
                             pool_sd = FALSE) {
  # pool_sd = TRUE: standardizing SD pooled across groups (Morris 2008, eq. 9/13)
  if (pool_sd) {
    return(.pooled_pre_post_to_smd(
      mean_pre_exp = mean_pre_exp, mean_pre_sd_exp = mean_pre_sd_exp,
      mean_exp = mean_exp, mean_sd_exp = mean_sd_exp,
      mean_pre_nexp = mean_pre_nexp, mean_pre_sd_nexp = mean_pre_sd_nexp,
      mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp,
      n_exp = n_exp, n_nexp = n_nexp,
      r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,
      pre_post_to_smd = pre_post_to_smd
    ))
  }

  res_exp <- .single_group_pre_post_to_smd(
    mean_pre = mean_pre_exp, mean_post = mean_exp,
    mean_pre_sd = mean_pre_sd_exp, mean_post_sd = mean_sd_exp,
    n = n_exp, r_pre_post = r_pre_post_exp,
    pre_post_to_smd = pre_post_to_smd
  )

  res_nexp <- .single_group_pre_post_to_smd(
    mean_pre = mean_pre_nexp, mean_post = mean_nexp,
    mean_pre_sd = mean_pre_sd_nexp, mean_post_sd = mean_sd_nexp,
    n = n_nexp, r_pre_post = r_pre_post_nexp,
    pre_post_to_smd = pre_post_to_smd
  )

  d_final <- res_exp[,"d"] - res_nexp[,"d"]
  g_final <- res_exp[,"g"] - res_nexp[,"g"]
  vd_final <- res_exp[,"var_d"] + res_nexp[,"var_d"]
  vg_final <- res_exp[,"var_g"] + res_nexp[,"var_g"]

  d_ci_lo <- d_final - sqrt(vd_final) * qt(.975, n_exp + n_nexp - 2)
  d_ci_up <- d_final + sqrt(vd_final) * qt(.975, n_exp + n_nexp - 2)
  g_ci_lo <- g_final - sqrt(vg_final) * qt(.975, n_exp + n_nexp - 2)
  g_ci_up <- g_final + sqrt(vg_final) * qt(.975, n_exp + n_nexp - 2)

  res <- cbind(
    d_final, vd_final, d_ci_lo, d_ci_up,
    g_final, vg_final, g_ci_lo, g_ci_up
  )

  return(res)
}

################# POOLED TWO-GROUP PRE POST to SMD ##############
#' Between-group pre/post SMD, standardizing SD pooled across arms
#'
#' Point estimates are Morris (2008): the numerator is the difference in mean change
#' between arms; the standardizer is pooled across arms (baseline SD for bonett /
#' d_ppc2 eq. 8-9; change SD for d_z/d_rm; average SD for d_av / d_ppc3 eq. 12-14).
#'
#' Each branch's variance follows metafor's LS *pattern* for the corresponding measure,
#' an empirical leading term plus a g^2 term. It is a pattern rather than a single
#' closed-form rule: metafor puts n (or n1+n2) in the g^2 denominator of its
#' non-heteroscedastic LS forms rather than the standardizer df, while evaluating J at
#' the standardizer df, and the two df differ by design. The g^2 term below is therefore
#' over 2N, and J is J(nu) with nu the standardizer df. Var(d) = Var(g)/J^2 exactly,
#' since g = J*d with J a deterministic constant.
#'
#' The leading term is heteroscedasticity-robust, being built from the empirical pooled
#' change SD rather than from Var(change) = 2*sigma^2*(1-r), which is valid only when
#' SD_pre = SD_post. For bonett it reduces exactly to Viechtbauer's published two-group
#' pre/post variance vi = 2(1-r)(1/nT + 1/nC) + g^2/(2N) when SD_pre = SD_post within
#' each arm, the df-weighted r_avg making that reduction exact even when n1 != n2. Under
#' heteroscedasticity it departs from the homoscedastic form, which understates the
#' variance, with coverage falling to about 0.86 at SD_pre/SD_post around 0.64. d_av
#' takes its fourth-moment g^2 coefficient from Bonett (2008) eq. 19 but not its leading
#' term; see the d_av entry below.
#'
#' Publication status of each pooled variance:
#'   - bonett: reduces to Viechtbauer's published vi (metafor-project Morris-2008 page)
#'             under homoscedasticity; the robust departure is metafor SMCRH / Bonett
#'             (2008) generalized across arms.
#'   - d_z:    exactly metafor::escalc(measure = "SMD", vtype = "LS") on the change
#'             scores (Hedges 1981), a published two-group variance.
#'   - d_rm:   d_rm = d_z * sqrt(2(1-r)) (Caldwell & Vigotsky 2020 eq. 13, a definition),
#'             so Var(d_rm) = 2(1-r)*Var(d_z) for known r.
#'   - d_av:   partly Bonett (2008) eq. 19 (two-group mixed design, all-four-SD
#'             standardizer). Only the fourth-moment g^2 coefficient comes from eq. 19,
#'             and it is reproduced exactly: eq. 19's first bracket
#'             [(s1^4+s2^4+2 r12^2 s1^2 s2^2)/df1 + (s3^4+s4^4+2 r34^2 s3^2 s4^2)/df2]
#'             /(32 s^4) is identical to g2_coef below (checked to 1e-14 over a grid of
#'             arm sizes). The leading term is not eq. 19's. Eq. 19 uses the df-based
#'             sum Sc1^2/(n1-1) + Sc2^2/(n2-1), whereas this code uses metafor's n-based
#'             pooled form (sd_change_pooled^2/sd_pooled^2) * N/(n1*n2), the two-sample
#'             SMD leading term. Under homoscedasticity the ratio of the two is
#'             (1/(n1-1) + 1/(n2-1)) / (1/n1 + 1/n2), so eq. 19 is larger by about 2% at
#'             n = 50/50, 3% at 30/30, 10% at 12/11 and 32% at 100/4. With
#'             heteroscedastic SDs the departure reaches about 38% at n = 100/10 and is
#'             unbounded as min(n1, n2) -> 2. The n-based form is the conventional and
#'             more conservative (smaller-variance) choice, but at strongly unequal arm
#'             sizes it is not eq. 19 and should not be described as such.
#' The pooled forms without a directly published two-group source (the robust bonett
#' departure, d_rm's algebra, and d_av's leading term) are Monte-Carlo calibrated in
#' tests_save/checked/test-pooled-variance-calibration.R, whose coverage grid is
#' n1, n2 in {10, 29, 90}. That covers imbalance up to 9:1, where the departure above is
#' about 10%. Ratios beyond that, and a small min(n) combined with a large imbalance,
#' are uncalibrated.
#'
#' CIs use qt(.975, m) with m = N - 2, the pooled standardizer df, on every branch. That
#' includes d_av, whose Cousineau effective df nu = 2m/(1 + r2_avg) feeds the bias
#' correction J alone and not the interval; the single-group kernel follows the same
#' split, putting its d_av CI on n - 1 rather than on its Cousineau df. qt is a
#' small-sample-conservative choice and slightly over-covers at n around 10, whereas
#' metafor and Bonett (2008) use a normal (z) critical value.
#'
#' References:
#'   Morris (2008) ORM 11(2):364-386 -- d_ppc1/d_ppc2/d_ppc3 point estimates (eq. 6-14)
#'   Hedges (1981); Viechtbauer (2007) JEBS 32(1):39-60 -- two-sample SMD "LS"
#'   Bonett (2008) Psych Methods 13(2):99-109 eq. 10/19 -- heteroscedasticity-robust var
#'   Caldwell & Vigotsky (2020) PeerJ 8:e10314 eq. 13 -- d_rm = d_z * sqrt(2(1-r))
#'   Cousineau (2020) TQMP 16(4):418-421 eq. 2 -- the (1+r^2) effective df (d_av J only)
#'
#' @noRd
.pooled_pre_post_to_smd <- function(mean_pre_exp, mean_pre_sd_exp,
                                     mean_exp, mean_sd_exp,
                                     mean_pre_nexp, mean_pre_sd_nexp,
                                     mean_nexp, mean_sd_nexp,
                                     n_exp, n_nexp,
                                     r_pre_post_exp, r_pre_post_nexp,
                                     pre_post_to_smd) {
  if (pre_post_to_smd == "cooper") {
    pre_post_to_smd <- "morris_drm"
  }

  r_pre_post_exp <- .guard_r_pre_post(r_pre_post_exp)
  r_pre_post_nexp <- .guard_r_pre_post(r_pre_post_nexp)

  # A negative SD is invalid input (see .single_group_pre_post_to_smd): NA it out so a
  # negative raw SD cannot flip the sign of the sd_change cross term and survive.
  mean_pre_sd_exp  <- ifelse(is.finite(mean_pre_sd_exp)  & mean_pre_sd_exp  >= 0, mean_pre_sd_exp,  NA_real_)
  mean_sd_exp      <- ifelse(is.finite(mean_sd_exp)      & mean_sd_exp      >= 0, mean_sd_exp,      NA_real_)
  mean_pre_sd_nexp <- ifelse(is.finite(mean_pre_sd_nexp) & mean_pre_sd_nexp >= 0, mean_pre_sd_nexp, NA_real_)
  mean_sd_nexp     <- ifelse(is.finite(mean_sd_nexp)     & mean_sd_nexp     >= 0, mean_sd_nexp,     NA_real_)

  N <- n_exp + n_nexp
  m <- N - 2 # pooled degrees of freedom
  # Each arm needs >= 2 observations: the (n-1) pooling weights must be positive and
  # the between-arm factor N/(n1*n2) must be finite. m > 0 alone does not ensure this
  # (e.g. n_exp = 0, n_nexp = 30 gives m = 28 but a negative weight and T_change = Inf).
  valid_n <- is.finite(n_exp) & is.finite(n_nexp) & n_exp >= 2 & n_nexp >= 2
  m <- ifelse(is.finite(m) & m > 0 & valid_n, m, NA_real_)

  change_exp <- mean_exp - mean_pre_exp
  change_nexp <- mean_nexp - mean_pre_nexp
  mean_diff <- change_exp - change_nexp

  # df-weighted mean of the per-arm correlations. df-weighting (not n-weighting) is
  # what makes the robust bonett/dz variances reduce exactly to Viechtbauer's published
  # vi = 2(1-r)(1/nT + 1/nC) + g^2/(2N) under homoscedasticity when the arms' SDs are
  # equal but n1 != n2 (an n-weighted mean leaves a residual there).
  r_avg <- ((n_exp - 1) * r_pre_post_exp + (n_nexp - 1) * r_pre_post_nexp) / m

  # Pooled change SD: the empirical scale of the numerator. Used by every branch,
  # so that no branch assumes SD_pre = SD_post.
  sd_change_exp <- sqrt(mean_pre_sd_exp^2 + mean_sd_exp^2 -
                        .pre_post_cross(r_pre_post_exp, mean_pre_sd_exp, mean_sd_exp))
  sd_change_nexp <- sqrt(mean_pre_sd_nexp^2 + mean_sd_nexp^2 -
                         .pre_post_cross(r_pre_post_nexp, mean_pre_sd_nexp, mean_sd_nexp))
  sd_change_pooled <- sqrt(((n_exp - 1) * sd_change_exp^2 +
                            (n_nexp - 1) * sd_change_nexp^2) / m)

  # Var(mean_diff) on the change scale, in units of sd_change_pooled^2
  T_change <- N / (n_exp * n_nexp)

  if (pre_post_to_smd == "bonett") {
    # Morris (2008) d_ppc2: numerator = difference in mean change, standardizer =
    # baseline SD pooled across arms (eq. 8-9). Variance is the SMCRH form, the
    # numerator's empirical variance expressed in standardizer units.
    sd_pooled <- .guard_standardizer(sqrt(((n_exp - 1) * mean_pre_sd_exp^2 +
                                           (n_nexp - 1) * mean_pre_sd_nexp^2) / m))
    nu <- m
    J <- .d_j(nu)
    d <- mean_diff / sd_pooled
    g <- d * J

    T1 <- (sd_change_pooled^2 / sd_pooled^2) * T_change
    var_g <- T1 + g^2 / (2 * N)
    # g = J*d with J a deterministic constant, so Var(d) = Var(g)/J^2 exactly.
    var_d <- var_g / J^2

  } else if (pre_post_to_smd == "morris_dz") {
    # Change-score metric. Once the change SD is pooled across arms this is exactly
    # an independent-groups Hedges g computed on the change scores, so its variance
    # is metafor::escalc(measure = "SMD", vtype = "LS"), with no r and no 2(1-r) term.
    sd_pooled <- .guard_standardizer(sd_change_pooled)
    nu <- m
    J <- .d_j(nu)
    d <- mean_diff / sd_pooled
    g <- d * J

    var_g <- T_change + g^2 / (2 * N)
    var_d <- var_g / J^2

  } else if (pre_post_to_smd == "morris_drm") {
    # Raw-score metric: d_rm = d_z * sqrt(2(1-r)) (Caldwell & Vigotsky 2020 eq. 13).
    # r is a known constant, so Var(d_rm) = 2(1-r) * Var(d_z): the 2(1-r) here is a
    # deterministic rescaling rather than a homoscedasticity assumption.
    sd_pooled <- .guard_standardizer(sd_change_pooled)
    nu <- m
    J <- .d_j(nu)
    k <- sqrt(2 * (1 - r_avg))
    d <- (mean_diff / sd_pooled) * k
    g <- d * J

    T1 <- 2 * (1 - r_avg) * T_change
    var_g <- T1 + g^2 / (2 * N)
    var_d <- var_g / J^2

  } else if (pre_post_to_smd == "morris_dav") {
    # Morris (2008) d_ppc3: standardizer = quadratic mean of the pre and post SDs,
    # pooled across arms (eq. 12-13). Variance is the two-group heteroscedasticity-
    # robust form Bonett (2008) eq. 19 (the mixed-design "difference-in-differences"
    # standardized by all four cell SDs), equivalently metafor SMCRPH carried across
    # two independent arms: a robust change-SD leading term plus a per-arm
    # fourth-moment g^2 term. Note that this is not the homoscedastic SMCRP
    # coefficient (1+r^2)/(4N), which is only its SD_pre = SD_post special case.
    sd_av_exp <- sqrt((mean_pre_sd_exp^2 + mean_sd_exp^2) / 2)
    sd_av_nexp <- sqrt((mean_pre_sd_nexp^2 + mean_sd_nexp^2) / 2)
    sd_pooled <- .guard_standardizer(sqrt(((n_exp - 1) * sd_av_exp^2 +
                                           (n_nexp - 1) * sd_av_nexp^2) / m))
    # nu is used for the bias correction J alone. The (1+r^2) effective df is
    # homoscedastic (Cousineau 2020 eq. 2), so we feed it the df-weighted mean of r^2
    # (its two-arm generalization); the variance itself uses the fourth moments below.
    r2_avg <- ((n_exp - 1) * r_pre_post_exp^2 + (n_nexp - 1) * r_pre_post_nexp^2) / m
    nu <- 2 * m / (1 + r2_avg)
    J <- .d_j(nu)
    d <- mean_diff / sd_pooled
    g <- d * J

    T1 <- (sd_change_pooled^2 / sd_pooled^2) * T_change
    # Bonett (2008) eq. 19 fourth-moment g^2 coefficient (sum of per-arm terms over
    # 32*sd_pooled^4). Reduces to (1+r^2)/(4N)'s n-1 analogue when SD_pre = SD_post.
    fm_exp <- mean_pre_sd_exp^4 + mean_sd_exp^4 +
              2 * r_pre_post_exp^2 * mean_pre_sd_exp^2 * mean_sd_exp^2
    fm_nexp <- mean_pre_sd_nexp^4 + mean_sd_nexp^4 +
               2 * r_pre_post_nexp^2 * mean_pre_sd_nexp^2 * mean_sd_nexp^2
    g2_coef <- (fm_exp / (n_exp - 1) + fm_nexp / (n_nexp - 1)) / (32 * sd_pooled^4)
    var_g <- T1 + g^2 * g2_coef
    var_d <- var_g / J^2
  }

  # CI on the pooled standardizer df (N - 2 = m) for every branch. dav's Cousineau
  # effective df nu = 2m/(1 + r2_avg) feeds the bias correction J alone, as in the
  # single-group kernel, whose dav CI likewise sits on n - 1 rather than on its
  # Cousineau df. The variance is built on m, so the interval is too. For bonett, dz
  # and drm nu equals m, so this matters only on the dav path.
  d_ci_lo <- d - sqrt(var_d) * qt(.975, m)
  d_ci_up <- d + sqrt(var_d) * qt(.975, m)
  g_ci_lo <- g - sqrt(var_g) * qt(.975, m)
  g_ci_up <- g + sqrt(var_g) * qt(.975, m)

  res <- cbind(
    d, var_d, d_ci_lo, d_ci_up,
    g, var_g, g_ci_lo, g_ci_up
  )

  return(res)
}

################# SINGLE GROUP PRE POST to SMD ##############
#' Calculate within-group standardized mean difference for a single group
#'
#' @param mean_pre mean at baseline (pre-test)
#' @param mean_post mean at follow-up (post-test)
#' @param mean_pre_sd standard deviation at baseline
#' @param mean_post_sd standard deviation at follow-up
#' @param n sample size
#' @param r_pre_post pre-post correlation
#' @param pre_post_to_smd method to use: "bonett" or "cooper"
#'
#' @return matrix with columns: d, var_d, d_ci_lo, d_ci_up, g, var_g, g_ci_lo, g_ci_up
#'
#' @noRd
.single_group_pre_post_to_smd <- function(mean_pre, mean_post,
                                           mean_pre_sd, mean_post_sd,
                                           n, r_pre_post,
                                           pre_post_to_smd) {
  # "cooper" alias kept for direct internal calls
  if (pre_post_to_smd == "cooper") {
    pre_post_to_smd <- "morris_drm"
  }

  r_pre_post <- .guard_r_pre_post(r_pre_post)

  # A negative SD is invalid input. NA it out so the ES is NA rather than quietly
  # corrupted: a negative raw SD flips the sign of the -2*r*sd_pre*sd_post cross
  # term in sd_diff, yielding a wrong-but-finite standardizer that the .guard_standardizer
  # on the assembled denominator cannot catch. Zero is allowed: the mean-change wrappers
  # legitimately pass mean_pre_sd = 0 (the pre slot is zeroed by construction).
  mean_pre_sd  <- ifelse(is.finite(mean_pre_sd)  & mean_pre_sd  >= 0, mean_pre_sd,  NA_real_)
  mean_post_sd <- ifelse(is.finite(mean_post_sd) & mean_post_sd >= 0, mean_post_sd, NA_real_)

  # An arm needs n >= 2 for a variance/CI to exist. NA the standardizing SDs below 2 so
  # d, var and CI are all NA. morris_drm's var_d = 2*(1-r)/n + d^2/(2n) carries no J, so
  # the J(n-1) = NA path that nulls the other three branches leaves it finite at
  # n = 1 (and Inf at n = 0); nulling the SDs here closes that leak and matches the
  # pooled kernel's explicit n >= 2 guard.
  bad_n <- !(is.finite(n) & n >= 2)
  mean_pre_sd  <- ifelse(bad_n, NA_real_, mean_pre_sd)
  mean_post_sd <- ifelse(bad_n, NA_real_, mean_post_sd)

  if (pre_post_to_smd == "bonett") {
    # Bonett method: standardize by baseline SD
    # Matches metafor SMCRH (heteroscedastic-robust variance formula)
    J <- .d_j(n - 1)

    var_change <- mean_pre_sd^2 + mean_post_sd^2 -
      .pre_post_cross(r_pre_post, mean_pre_sd, mean_post_sd)

    sd_std <- .guard_standardizer(mean_pre_sd)
    d <- (mean_post - mean_pre) / sd_std
    g <- d * J

    # metafor SMCRH heteroscedastic variance formula (Bonett 2008)
    # var = sd_change^2 / (sd1i^2 * (n-1)) + g^2 / (2 * (n-1))
    var_g <- var_change / (sd_std^2 * (n - 1)) + g^2 / (2 * (n - 1))
    var_d <- var_g / (J^2)

    d_ci_lo <- d - sqrt(var_d) * qt(.975, n - 1)
    d_ci_up <- d + sqrt(var_d) * qt(.975, n - 1)
    g_ci_lo <- g - sqrt(var_g) * qt(.975, n - 1)
    g_ci_up <- g + sqrt(var_g) * qt(.975, n - 1)

    res <- cbind(
      d = d, var_d = var_d, d_ci_lo = d_ci_lo, d_ci_up = d_ci_up,
      g = g, var_g = var_g, g_ci_lo = g_ci_lo, g_ci_up = g_ci_up
    )

    return(res)
  } else if (pre_post_to_smd == "morris_drm") {
    # Morris d_rm (alias "cooper"): d_rm = d_z * sqrt(2(1-r)) (Caldwell & Vigotsky 2020),
    # that is, the change SD rescaled onto the raw-score metric. It equals metafor's SMCR
    # (baseline/raw-SD standardizer) only under homoscedasticity (SD_pre = SD_post); when
    # the two SDs differ they diverge, so d_rm is defined here by the Caldwell rescaling
    # rather than by SMCR. Its variance below is 2(1-r) * Var(d_z), the exact scaling of
    # the SMCC form.
    # Note that the mean_change wrappers pass sd_change through mean_post_sd
    # (mean_pre = 0 and mean_pre_sd = 0, so sd_diff reduces to sd_change).
    J <- .d_j(n - 1)

    sd_diff <- .guard_standardizer(sqrt(mean_pre_sd^2 + mean_post_sd^2 -
                    .pre_post_cross(r_pre_post, mean_pre_sd, mean_post_sd)))

    d <- (mean_post - mean_pre) / sd_diff * sqrt(2 * (1 - r_pre_post))
    g <- d * J

    # Viechtbauer corrected formula
    var_d <- 2 * (1 - r_pre_post) / n + d^2 / (2 * n)
    var_g <- J^2 * var_d

    d_ci_lo <- d - sqrt(var_d) * qt(.975, n - 1)
    d_ci_up <- d + sqrt(var_d) * qt(.975, n - 1)
    g_ci_lo <- g - sqrt(var_g) * qt(.975, n - 1)
    g_ci_up <- g + sqrt(var_g) * qt(.975, n - 1)

    res <- cbind(
      d = d, var_d = var_d, d_ci_lo = d_ci_lo, d_ci_up = d_ci_up,
      g = g, var_g = var_g, g_ci_lo = g_ci_lo, g_ci_up = g_ci_up
    )

    return(res)
  } else if (pre_post_to_smd == "morris_dz") {
    # Morris & DeShon d_z: standardize by change score SD
    # Matches metafor SMCC (change score standardization)
    # Most conservative when r is high
    J <- .d_j(n - 1)

    sd_diff <- .guard_standardizer(sqrt(mean_pre_sd^2 + mean_post_sd^2 -
                    .pre_post_cross(r_pre_post, mean_pre_sd, mean_post_sd)))

    d <- (mean_post - mean_pre) / sd_diff
    g <- d * J

    # metafor SMCC variance formula (uses corrected g, not d)
    # var = 1/n + g^2/(2n)
    var_g <- 1 / n + g^2 / (2 * n)
    var_d <- var_g / (J^2)

    d_ci_lo <- d - sqrt(var_d) * qt(.975, n - 1)
    d_ci_up <- d + sqrt(var_d) * qt(.975, n - 1)
    g_ci_lo <- g - sqrt(var_g) * qt(.975, n - 1)
    g_ci_up <- g + sqrt(var_g) * qt(.975, n - 1)

    res <- cbind(
      d = d, var_d = var_d, d_ci_lo = d_ci_lo, d_ci_up = d_ci_up,
      g = g, var_g = var_g, g_ci_lo = g_ci_lo, g_ci_up = g_ci_up
    )

    return(res)
  } else if (pre_post_to_smd == "morris_dav") {
    # Morris (2008) d_av (d_ppc3): standardize by the quadratic mean of the pre and
    # post SDs. Its sampling variance is metafor SMCRPH, equivalently Bonett (2008)
    # eq. 10: a heteroscedasticity-robust change-SD leading term plus a fourth-moment
    # g^2 term. It is not the homoscedastic SMCRP form 2(1-r)/n + g^2(1+r^2)/(4n),
    # which assumes SD_pre = SD_post and understates var(g) when they differ, with CI
    # coverage falling to about 0.93 as SD_pre/SD_post moves away from 1. Morris (2008,
    # p.384) advised against d_av because its variance was unknown to him; Bonett (2008)
    # supplies it and this branch uses it.
    # nu = 2(n-1)/(1+r^2) is used for the Hedges bias correction J alone (Cousineau
    # 2020 eq. 2; metafor SMCRP/SMCRPH); the variance itself is on df = n-1.
    mi <- 2 * (n - 1) / (1 + r_pre_post^2)
    J <- .d_j(mi)

    sd_av <- .guard_standardizer(sqrt((mean_pre_sd^2 + mean_post_sd^2) / 2))

    d <- (mean_post - mean_pre) / sd_av
    g <- d * J

    # metafor SMCRPH "LS" is Bonett (2008) eq. 10, verified as a formula identity with
    # zero difference against metafor::escalc(measure = "SMCRPH"). It does not reproduce
    # the var = 0.0148 printed for Bonett's worked Example 2 (n = 60): this branch gives
    # 0.0146674131, and redoing Bonett's own arithmetic for that example gives
    # 0.0147449109, so 0.0148 appears to be a rounding slip in the paper rather than a
    # package error. Recomputing eq. 10 by hand with the bias-corrected g used here
    # returns 0.0146674131, identical to 12 digits.
    sd_diff2 <- mean_pre_sd^2 + mean_post_sd^2 -
      .pre_post_cross(r_pre_post, mean_pre_sd, mean_post_sd)
    fm <- mean_pre_sd^4 + mean_post_sd^4 + 2 * r_pre_post^2 * mean_pre_sd^2 * mean_post_sd^2
    var_g <- sd_diff2 / (sd_av^2 * (n - 1)) + g^2 * fm / (8 * sd_av^4 * (n - 1))
    var_d <- var_g / (J^2)

    d_ci_lo <- d - sqrt(var_d) * qt(.975, n - 1)
    d_ci_up <- d + sqrt(var_d) * qt(.975, n - 1)
    g_ci_lo <- g - sqrt(var_g) * qt(.975, n - 1)
    g_ci_up <- g + sqrt(var_g) * qt(.975, n - 1)

    res <- cbind(
      d = d, var_d = var_d, d_ci_lo = d_ci_lo, d_ci_up = d_ci_up,
      g = g, var_g = var_g, g_ci_lo = g_ci_lo, g_ci_up = g_ci_up
    )

    return(res)
  }
}

################# R/Z to SMD ###################
# Vectorised front end for .cor_to_smd().
#
# .cor_to_smd() is scalar, and calling it per row through mapply() is slow. Its
# "viechtbauer" branch, the default for cor_to_smd, builds a fresh
# metafor::conv.delta() call per row even though conv.delta is itself vectorised,
# which costs 18.1 s per 10,000 rows against 0.34 s for a single vectorised call.
#
# Memoisation (.mapply_memo) does not help here as it does on the tetrachoric path,
# because correlations are continuous and rows are effectively all distinct. Real
# vectorisation is the only option. Rows are grouped by their cor_to_smd value, which
# is a per-row column; each group is computed in one shot, and the groups are then
# reassembled in the original order.
#
# All three branches agree with the per-row path to within 1-2 ULP. The "viechtbauer"
# d is bit-identical on every row, but its SE differs by 1 ULP on about 0.02% of rows
# (44 out of 200,000, max |rel| 2.2e-16), because conv.delta's numerical derivative is
# evaluated at slightly different precision in a vector call. The closed-form "cooper"
# and "mathur" branches differ by 1-2 ULP as well (max |rel| 4.3e-16, about 0.02% of
# rows), because R evaluates length-1 and length-n arithmetic on different code paths.
#
# All arguments must be per-row vectors of the same length. Unlike mapply, this does
# not recycle a scalar, so every call site passes rep(x, length.out = n).
.cor_to_smd_vec <- function(r, r_se, unit_increase_iv, sd_iv, unit_type,
                            n_sample, cor_to_smd) {
  n <- length(r)
  out <- matrix(NA_real_, nrow = n, ncol = 2)
  meth <- as.character(cor_to_smd)

  for (m in unique(meth[!is.na(meth)])) {
    i <- which(meth == m)
    out[i, ] <- if (m == "viechtbauer") {
      # transf.rtod output is the g; d is backed out as g / J. See the rationale
      # (and the simulation that settles it) in .cor_to_smd() below.
      rg <- metafor::conv.delta(yi = r[i], vi = r_se[i]^2,
                                transf = metafor::transf.rtod,
                                var.names = c("g", "g_var"))
      J <- .d_j(n_sample[i] - 2)
      cbind(rg$g / J, sqrt(rg$g_var / J^2))
    } else if (m == "cooper") {
      cbind(2 * r[i] / sqrt(1 - r[i]^2),
            sqrt(4 * r_se[i]^2 / ((1 - r[i]^2)^3)))
    } else if (m == "mathur") {
      increase <- ifelse(unit_type[i] == "sd",
                         unit_increase_iv[i] * sd_iv[i], unit_increase_iv[i])
      d <- r[i] * increase / (sd_iv[i] * sqrt(1 - r[i]^2))
      cbind(d, abs(d) * sqrt(1 / (r[i]^2 * (n_sample[i] - 3)) +
                              1 / (2 * (n_sample[i] - 1))))
    } else {
      matrix(NA_real_, nrow = length(i), ncol = 2)
    }
  }
  out
}

.cor_to_smd <- function(r, r_se,
                        unit_increase_iv, sd_iv, unit_type,
                        n_sample, cor_to_smd) {
  if (cor_to_smd == "mathur") {
    increase <- ifelse(unit_type == "sd",
                       unit_increase_iv * sd_iv,
                       unit_increase_iv)

    d <- r * increase / (sd_iv * sqrt((1 - r^2)))
    d_se <- abs(d) * sqrt(1 / (r^2 * (n_sample - 3)) + 1 / (2 * (n_sample - 1)))
    res <- cbind(d, d_se)
    return(res)
  } else if (cor_to_smd == "viechtbauer") {
    # transf.rtod(r_hat) lands in the g slot rather than the d slot. Since
    # transf.rtod is a population map with no df term, its output looks as though
    # it ought to be an uncorrected d, but the choice of slot is a question about
    # estimator bias rather than about the algebra of the transform.
    #
    # Checked by simulation (bivariate normal, rho = 0.5, median split, true
    # delta = 0.870126, nrep = 2e5), bias of the reported g against delta:
    #
    #   n     Hedges g from raw data    g = transf.rtod(r_hat)   g = J*transf.rtod
    #   10          -0.004                    +0.033                  -0.055
    #   25          -0.000                    +0.012                  -0.017
    #  200          -0.000                    +0.002                  -0.002
    #
    # Putting transf.rtod in g is closer to unbiased at every n, because this
    # route's small-sample bias is not the Hedges bias. E[r_hat] is biased
    # downwards, which partly cancels the upward bias that the pooled-SD
    # denominator induces in a directly computed d, so applying J on top
    # over-corrects. Neither convention is exactly unbiased; the residual is the
    # uncancelled remainder, and it stays below 0.02 for n >= 25.
    res_g <- metafor::conv.delta(
      yi = r, vi = r_se^2, transf = metafor::transf.rtod, var.names = c("g", "g_var")
    )
    J <- .d_j(n_sample - 2)
    res_g$d <- res_g$g / J
    res_g$d_se <- sqrt(res_g$g_var / J^2)
    res <- cbind(res_g$d, res_g$d_se)

    return(res)
  } else if (cor_to_smd == "cooper") {
    d <- 2 * r / sqrt(1 - r^2)
    d_se <- sqrt(4 * r_se^2 / ((1 - r^2)^3))
    res <- cbind(d, d_se)

    return(res)
  }
}



.validate_pre_post_to_smd <- function(pre_post_to_smd, allowed_methods, context, func_name) {
  pre_post_to_smd_normalized <- ifelse(pre_post_to_smd == "cooper", "morris_drm", pre_post_to_smd)

  invalid <- !pre_post_to_smd_normalized %in% allowed_methods
  if (any(invalid)) {
    stop(paste0(
      "Invalid 'pre_post_to_smd' argument in ", func_name, "().\n\n",
      "Provided: '", paste(unique(pre_post_to_smd[invalid]), collapse = "', '"), "'\n",
      "Allowed methods: '", paste(allowed_methods, collapse = "', '"), "'\n\n",
      .get_method_restriction_rationale(context)
    ), call. = FALSE)
  }

  return(pre_post_to_smd_normalized)
}

.get_method_restriction_rationale <- function(context) {
  if (context == "mean_change") {
    paste0(
      "For mean change data, two standardization methods are available:\n\n",
      "  - 'morris_drm' (alias: 'cooper'): Raw score standardizer [DEFAULT]\n",
      "      d_rm = (mean_change / sd_change) * sqrt(2*(1-r))\n",
      "      REQUIRES pre-post correlation (r). Converts d_z to d_rm.\n\n",
      "  - 'morris_dz': Change score standardizer\n",
      "      d_z = mean_change / sd_change\n",
      "      INDEPENDENT of r. Direct standardization by change SD.\n\n",
      "Other methods not applicable:\n",
      "  - 'bonett': requires baseline SD (not available with mean change data)\n",
      "  - 'morris_dav': requires separate pre/post SDs (not available)\n\n",
      "See Morris & DeShon (2002) for guidance on choosing between d_rm and d_z."
    )
  } else if (context == "paired_t") {
    paste0(
      "For paired t-test data, two standardization methods are available:\n\n",
      "  - 'morris_drm' (alias: 'cooper'): Raw score standardizer [DEFAULT]\n",
      "      d_rm = t * sqrt(2*(1-r)/n)\n",
      "      REQUIRES pre-post correlation (r). Most common in meta-analysis.\n\n",
      "  - 'morris_dz': Change score standardizer\n",
      "      d_z = t / sqrt(n)\n",
      "      INDEPENDENT of r. Use when correlation is unknown or when\n",
      "      synthesizing with other d_z estimates.\n\n",
      "Note: d_rm and d_z are on different scales. See Morris & DeShon (2002)\n",
      "for guidance on choosing between them."
    )
  } else if (context == "pre_post_means") {
    paste0(
      "For pre-post means data, four standardization methods are available:\n\n",
      "  - 'bonett': baseline-SD standardizer (Morris 2008 d_ppc2). Sensitive to a\n",
      "      baseline SD that is restricted relative to the endpoint SD (e.g. by\n",
      "      eligibility cut-offs), which inflates the effect size.\n",
      "  - 'morris_drm' (alias: 'cooper'): change SD rescaled by sqrt(2(1-r)) onto\n",
      "      the raw-score metric. Depends on r_pre_post; assumes SD_pre = SD_post.\n",
      "  - 'morris_dz': change-SD standardizer. Independent of r, but NOT on the\n",
      "      same metric as an endpoint SMD -- do not pool the two (Cochrane 10.5.2).\n",
      "  - 'morris_dav': average-SD standardizer (Morris 2008 d_ppc3). Morris\n",
      "      recommends AGAINST it (2008, p.384): its sampling variance was unknown\n",
      "      and it is downward-biased when the post-treatment SD inflates.\n\n",
      "No single method is best in all cases; the choice changes the estimand.\n",
      "See Morris (2008) and Morris & DeShon (2002) for detailed comparisons."
    )
  } else {
    ""
  }
}


# tryCatch({
#   validate_positive(n_cases_exp, n_cases_nexp, n_controls_exp, n_controls_nexp,
#                     error_message = "The number of cases/controls in the exposed/non-exposed groups should be >0.")
# }, error = function(e) {
#   stop("Validation failed:", conditionMessage(e), "\n")
# })
#
#
# tryCatch({
#   validate_ci_symmetry(value, ci_lo, ci_up, func = "example_function",
#                        max_asymmetry_percent = 5)
# }, error = function(e) {
#   stop("Validation failed:", conditionMessage(e), "\n")
# })
#


# **A.** First, Cooper et al. (2019) - \code{table_2x2_to_cor = "cooper"} -
# proposes to convert the
# 2x2 table into a OR (formula above), to convert this OR into a SMD
# (see formula in \code{\link{es_from_or_se}()}), and to convert this
# SMD into a correlation coefficient (see formula in \code{\link{es_from_cohen_d}()},
# with the option \code{"smd_to_cor = 'lipsey_cooper'"}).
#
# **B.** Second, a correlation coefficient (more precisely - a phi coefficient)
# can be obtained from the contingency table using the formula given in
# Lipsey and Wilson (2001) - \code{table_2x2_to_cor = "lipsey"}.
# The formulas used to estimate the r and z are:
# \deqn{r = \frac{(n\_cases\_exp*n\_controls\_nexp - n\_controls\_exp*n\_cases\_nexp)}{\sqrt{(n\_exp) * (n\_nexp) * (n\_cases) * (n\_controls\_exp+n\_cases\_nexp)}}}
# \deqn{z = atanh(r)}
# \deqn{z\_se = logor\_se^2 * \frac{z^2}{\log(or)^2}}
# \deqn{z\_ci\_lo = z - qnorm(.975)*z\_se}
# \deqn{z\_ci\_up = z + qnorm(.975)*z\_se}
# \deqn{r\_ci\_lo = tanh(z\_ci\_lo)}
# \deqn{r\_ci\_up = tanh(z\_ci\_up)}
# \deqn{effective\_n = \frac{1}{z\_se^2 + 3}}
# \deqn{r\_se = \frac{(1 - r^2)^2}{effective\_n - 1}}

