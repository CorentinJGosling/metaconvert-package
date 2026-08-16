# # !!! for alll internal formulas, check whether nn_miss does not prevent the code to be ran
# # important to check the input of the method before calling these functions

# Row-wise mapply with memoisation on the argument tuple.
#
# The tetrachoric route (.contingency_to_cor -> .tet_r -> metafor::escalc with
# measure = "RTET") has no closed form and is solved row by row by numerical ML over a
# bivariate normal CDF (optim + mvtnorm::pmvnorm). Profiling es_from_2x2() and
# es_from_or_se() puts ~83% of the total time inside pmvnorm and ~93% inside .tet_r, so
# that single call dominates both. (.or_to_cor itself is closed-form arithmetic; it is
# routed through here only because it is scalar-valued and shares the call shape.)
#
# Because each result is a deterministic function of that row's arguments alone, rows
# with identical arguments are solved once and the solution reused. Counts repeat
# heavily in real datasets and overwhelmingly in simulation grids (~250 distinct tables
# per 1,500 draws at n = 50/arm): measured 44.75 s -> 1.31 s (34x) on 2,000 rows drawn
# from 60 distinct tables. With all-distinct rows it costs one paste() and two match()
# calls and is not slower (53.29 s -> 44.30 s on the same benchmark).
#
# Output is identical to t(mapply(FUN, ...)) for every input reachable from the exported
# functions. One theoretical caveat: paste() renders doubles at 15 significant digits, so
# two doubles agreeing to ~15 digits but not identical() share a key and hence a result.
# The induced relative error is ~1e-15 -- far below the ML solver's own tolerance -- and
# the memoised arguments are 2x2 counts, method names and reverse flags on one path and
# or / logor_se / margins on the other. Integer counts cannot collide.
# Which of estimraw::estim_raw()'s candidate 2x2 reconstructions are actually usable.
#
# When the reconstruction quadratic has a negative discriminant -- routine whenever the
# outcome is not rare -- estim_raw() still returns its candidate list, but with NaN cell
# counts. Distinguishing "no real solution" from "several real solutions" is what lets
# the two callers below emit an accurate message instead of one that over-counts the
# candidates and recommends a remedy that cannot work.
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
# mapply() only simplifies to a numeric matrix when EVERY call returns the same
# type. The conversion helpers (.or_to_rr, .rr_to_or, ...) mostly return a 1x4
# cbind() matrix, but .or_to_rr()'s 'dipietrantonj' branch returns a data.frame,
# on the successful reconstruction as much as on the no-real-solution one, and a
# single such row is enough to turn the whole result into a LIST-matrix. The
# column then extracts as a list, and assigning it into a numeric vector
# silently recycles values across rows: with or_to_rr supplied as a per-row
# column, one 'dipietrantonj' row made every other row inherit the FIRST row's
# standard error, leaving each SE inconsistent with its own confidence interval.
# Zero-length elements are mapped to NA rather than dropped, so the result can
# never be short (and hence can never be recycled).
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

    # Non-identifiability gate (see .rotation_tied). Only when the exact solve did
    # NOT fire: a solved table is identified whatever the margins look like.
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

    res <- cbind(
      logrr = calc_meta_cases$logrr,
      logrr_se = calc_meta_cases$logrr_se,
      logrr_ci_lo = calc_meta_cases$logrr_ci_lo,
      logrr_ci_up = calc_meta_cases$logrr_ci_up
    )

    return(res)
  } else if (or_to_rr == "metaumbrella_exp") {
    # n_cases / n_controls / baseline_risk are forwarded so the reconstruction can be
    # solved exactly rather than searched. All three are already in this function's
    # signature; before this they were received and discarded.
    contingency_meta_exp <- .estimate_n_from_or_and_n_exp(
      or = or, var = logor_se^2, n_exp = n_exp, n_nexp = n_nexp,
      n_cases = n_cases, n_controls = n_controls, baseline_risk = baseline_risk
    )

    # Non-identifiability gate (see .rotation_tied). Only when the exact solve did
    # NOT fire: a solved table is identified whatever the margins look like.
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

    res <- cbind(
      logrr = calc_meta_exp$logrr,
      logrr_se = calc_meta_exp$logrr_se,
      logrr_ci_lo = calc_meta_exp$logrr_ci_lo,
      logrr_ci_up = calc_meta_exp$logrr_ci_up
    )

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
      # of candidate solutions. It can also return NO usable candidate: for a
      # common outcome the reconstruction quadratic has a negative discriminant,
      # so the cell counts come back NaN. In that case which.min() below sees only
      # non-finite values and returns integer(0), and estim[[integer(0)]] raises
      # "attempt to select less than one element" -- aborting the whole
      # convert_df() run over a single unreconstructable row. Return NA instead:
      # one bad cell must never stop a run.
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
#' Solve a 2x2 table exactly from an odds ratio and ALL FOUR margins
#'
#' When both margin pairs are known the table is over-determined by one degree of
#' freedom, so the odds ratio pins it down exactly: no enumeration, no search, and --
#' critically -- \strong{no use of the reported variance}.
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
#' WHY THIS MATTERS. The two enumerating helpers below search for the candidate table
#' whose reconstructed variance best matches the reported one. Two failure modes follow.
#' (1) The 180-degree rotation \eqn{(a,b,c,d) \to (d,c,b,a)} preserves the odds ratio AND
#' \eqn{1/a+1/b+1/c+1/d} \emph{exactly}, and is admissible with the same arm sizes
#' precisely when \code{n_exp == n_nexp} (resp. \code{n_cases == n_controls}), so the
#' search faces an exact tie broken only by enumeration order -- measured 25-32% correct
#' at rare event rates. (2) On the \code{es_from_or()} route the \code{var} being matched
#' is itself imputed by \code{\link{.se_from_or}} and runs ~1.4x wide, so the search lands
#' on the wrong table \emph{before} any tie arises; with a realistically rounded OR the
#' current rule's MAE |logRR| is 0.10-0.19 even at unequal arms. Solving ignores
#' \code{var} entirely and removes both: MAE |logRR| 0.0002, branch hit 1.000 at every
#' event rate from 0.03 to 0.97.
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
  # report a case margin covering arms that are not in n_exp + n_nexp (flagged
  # [UNUSUAL], not [INVALID], by V15). Such a margin does not describe THIS table, so
  # the solve must decline rather than return a table built from mismatched inputs.
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

  a <- round(a)
  b <- n_exp - a
  cc <- n_cases - a
  d <- n_controls - b

  # Zero-cell guard. Rounding can land the solve on a table with an empty cell (2.14%
  # of solved tables over or in [0.1, 10]). The enumeration below has a purpose-built
  # +0.5 branch for those, so hand them back to it rather than emitting a cell of 0
  # that would make var(logOR) infinite.
  if (!is.finite(a) || min(a, b, cc, d) < 1) return(NULL)

  data.frame(n_cases_exp = a, n_cases_nexp = cc,
             n_controls_exp = b, n_controls_nexp = d)
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
#' admissible with the SAME reported inputs precisely when those two margins are equal --
#' \code{n_exp == n_nexp} for the \code{_exp} parameterisation, \code{n_cases ==
#' n_controls} for \code{_cases}. In that configuration the enumerating search faces an
#' exact tie broken only by enumeration order, so its answer is arbitrary.
#'
#' Measured on the \code{_exp} parameterisation, 782 usable draws at
#' \code{n_exp == n_nexp == 50} with no second margin: the returned table is the true
#' one 37.0% of the time and its rotation 66.6% (4.4% of true tables are their own
#' rotation, so those two overlap; 0.8% land on neither). Genuinely wrong: 62.3%.
#'
#' The cliff is sharp rather than gradual -- one participant of imbalance makes the
#' rotation inadmissible -- which is why the gate tests exact equality and nothing
#' looser. Hit rate by \code{|n_exp - n_nexp|}: 0.384 (0) -> 0.995 (1) -> 0.987 (2)
#' -> 0.995 (5) -> 0.982 (10). The \code{_cases} mirror behaves the same way at
#' \code{|n_cases - n_controls|}: 0.745 (0) -> 0.995 (1) -> 0.995 (2) -> 0.993 (5)
#' -> 0.993 (10). The two tied-cell rates are not comparable to each other (they come
#' from different draw distributions, and the published 0.352 for \code{_cases} from a
#' third); what reproduces across all of them is the size of the drop at zero.
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
#' and which columns would break the tie. No hit-rate percentage is quoted: the two
#' parameterisations were measured at different rates and both figures are properties
#' of a draw distribution, not of the method. What is true of every measurement, and
#' is what the message says, is that the tie is EXACT -- so the search has nothing to
#' decide on.
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
  # Rung 2 (baseline_risk) is deliberately NOT applied here: unlike the _n_exp mirror
  # it is not fully identifying on this parameterisation (measured 0.9965, because
  # c/(c+d) == b/(a+b) admits b + c == n_cases as a second solution), whereas n_exp is
  # exact (1.000 at every exposure prevalence 0.1-0.9).
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

  # ---- Otherwise: fall through to the enumeration below, UNCHANGED. No prior is
  # applied. Where neither rung fires the row is genuinely non-identified, and every
  # candidate rule tested was a bet on outcome coding: "assume events are the minority"
  # is 3.65x worse than the status quo on common outcomes, wrong on 74% of those rows,
  # and produced +44% pooled-RR bias end-to-end on a common-outcome review. Recoding an
  # outcome from "response" to "non-response" flips its answer on identical data, so
  # there is no principled default. Output here is bit-identical to pre-cascade.
  #
  # What DID change (item 1.8): the caller now refuses to publish a risk ratio built
  # on this path when the two supplied margins are equal, because there the search is
  # not merely imprecise but exactly tied. See .rotation_tied.
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
    # which is only defined while RR * BR < 1. It is applied to THREE values -- the
    # point estimate and both interval bounds -- and the upper bound is the largest of
    # them, so it leaves the domain first. The failure was therefore asymmetric and
    # silent: a row could return a perfectly plausible log OR beside a NaN standard
    # error and a half-open interval, with every log() wrapped in suppressWarnings().
    # A finite estimate with no usable variance is worse than no estimate at all,
    # because it looks poolable and is not.
    #
    # The caller (es_from_stand_RR.R:153) only routes rows here when rr, baseline_risk
    # and both CI bounds are non-NA, so a non-finite result is always a domain
    # violation and never a missing input. Return the whole quartet as NA and say why.
    #
    # NB the mirror direction is safe and is deliberately left alone: .or_to_rr()'s
    # grant branch computes or / (1 - BR + BR*or), whose denominator is positive for
    # every or > 0 and BR in (0, 1), so it has no domain boundary to cross.
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
    # The previous message advertised 'grant_2x2' and 'grant_CI'. Neither is accepted
    # anywhere -- the branch above is plain 'grant' -- so it named two values that could
    # not be selected while omitting two that could ('grant', 'dipietrantonj').
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
    # var(log RR) = 1/a - 1/n1 + 1/c - 1/n2: the two ARM-TOTAL terms are SUBTRACTED
    # (delta method for log of a binomial proportion; identical to es_from_2x2()'s
    # se_rr and the Cochrane Handbook). This is NOT the log-OR pattern 1/a+1/b+1/c+1/d
    # (all added). Matching the candidate table to logrr_se^2 with the arm-total terms
    # ADDED selected the wrong table along the RR-constrained family and biased the
    # reconstructed OR (round-trip recovery of a known OR failed). Mapping of the sweep
    # variables: a = n_cases - n_cases_nexp_sim, c = n_cases_nexp_sim,
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

  # NB there is deliberately no "lipsey" (phi) branch. A commented-out one lived here
  # until it was removed: an if-branch containing nothing but comments and, because it
  # had no return(), silently falling through to NULL had anyone re-enabled it.
  #
  # It is not coming back, and the reason is not that the formula was wrong (the version
  # kept here had already been corrected -- its denominator originally used the diagonal
  # (n_controls_exp + n_cases_nexp) where phi needs the controls MARGIN
  # (n_controls_exp + n_controls_nexp)). It is that phi is not a poolable estimand:
  # its attainable range is bounded by the margins, so studies of the same association
  # with different event rates report different phi values, and pooling them manufactures
  # heterogeneity that is pure margin artefact (I^2 rising from 14% to 88% as the primary
  # studies get LARGER, because tau^2 is pinned by the margins while the within-study
  # variance falls as 1/n). And the choice could not be offered responsibly even if it
  # were wanted: a 2x2 table with fixed n has three free parameters and the
  # dichotomised-bivariate-normal family also has three, so the latent-normal model is
  # saturated and no goodness-of-fit test can tell a user which estimand applies to their
  # data. See ?convert_df. Genuinely dichotomous variables should use a binary measure
  # (logor / rr / rd), not a correlation.
  if (table_2x2_to_cor == "tetrachoric") {
    res <- .tet_r(as.numeric(n_cases_exp),
                  as.numeric(n_controls_exp),
                  as.numeric(n_cases_nexp),
                  as.numeric(n_controls_nexp))
    res[res == "calculation failure"] <- NA

    # res[] <- lapply(res, function(x) as.numeric(as.character(x)))
    # On reverse: negate AND swap each CI (new_lo = -old_up, new_up = -old_lo). Negating
    # the bounds in place left lo > up -- a wrong-signed, inverted interval that did not
    # bracket the negated point estimate. This is the same defect that was fixed on the
    # OR path (see es_from_stand_OR.R, the .or_to_cor result block). The z interval is a
    # symmetric Wald interval (z +- z*sqrt(vz)) and the r interval is its tanh
    # back-transform; both are therefore odd-symmetric about 0, since tanh(-x) = -tanh(x).
    # So the reflected interval is exactly the interval recomputed around the negated
    # estimate, in both cases. Old bounds are saved first, since res is overwritten.
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
  # over a bivariate normal CDF, which needs 'mvtnorm'. mvtnorm is a Suggests of metafor,
  # NOT an Imports, so installing metafor does not bring it in and a perfectly ordinary
  # installation can lack it. Without this guard escalc() raised an error, the tryCatch
  # below swallowed it, and every 2x2-derived correlation and Fisher's z came back NA
  # with nothing said -- indistinguishable from data that genuinely cannot support them.
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
      # The r-scale interval is the BACK-TRANSFORMED z interval, not a symmetric Wald
      # interval on r. A symmetric interval r +- z*sqrt(vr) is unbounded and routinely
      # escapes the parameter space: measured over a realistic grid it left [-1, 1] on
      # 45.8% of tables (and flag B1b -- the [INFO] disclosure that exists only because
      # of this -- fired on 20.6%, rising to 67% at n = 50 with |r| >= 0.6). Since vz is
      # already the delta-method transform of vr (line above), tanh() of the z bounds is
      # the same interval expressed on a scale where it cannot overshoot: escape falls to
      # 0.0% with coverage of the true correlation unchanged (96.1% -> 96.2%).
      r_lo <- tanh(z_lo)
      r_up <- tanh(z_up)

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
      c <- (1 - abs(n_exp/n_sample - n_cases/n_sample) / 5 - (1 / 2 - small_margin_prop)^2) / 2
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
  # closed form for the biserial correlation, which is a function of n, p and r ONLY.
  # That closed form is not an independent estimator: to leading order it IS the delta
  # propagation of the CRUDE d variance (they agree to ~0.3% at n = 400, and exactly up
  # to (n-1)/(n-2) at r = 0). So using it verbatim silently throws away whatever
  # precision the d actually has -- the Cooper eq. 12.26 (1 - R^2) shrink on an ANCOVA
  # row, the 2(1 - r_pre_post) factor on a pre-post row, the control-group df on a
  # Glass row, or a user-reported standard error. Rescaling by vd / vd_crude restores
  # it while leaving the crude case BIT-IDENTICAL (the ratio is exactly 1, since
  # vd_crude below is the same expression .es_from_d() uses), so agreement with
  # metafor's measure = "RBIS" is preserved for the rows it applies to.
  vd_crude <- (n_exp + n_nexp) / (n_exp * n_nexp) + d^2 / (2 * (n_exp + n_nexp))
  prec_ratio <- ifelse(!is.na(vd) & is.finite(vd) & !is.na(vd_crude) & vd_crude > 0,
                       vd / vd_crude, 1)

  if (smd_to_cor == "viechtbauer") {
    # h encodes the finite-sample point-biserial identity r_pb = t / sqrt(t^2 + df),
    # which holds for the MARGINAL two-group design at df = n_exp + n_nexp - 2. The d
    # arriving here is always on the marginal (unadjusted) SD scale -- that is the
    # whole point of Cooper's eq. 12.23/12.24 convention, which keeps ANCOVA studies
    # poolable with unadjusted ones. Subtracting n_cov_ancova used to shrink h and
    # inflate r into a quantity that is NEITHER the marginal r_pb NOR the partial one
    # (the partial would additionally require h * (1 - cov_outcome_r^2)), and made the
    # reported correlation depend on how many covariates the source study happened to
    # adjust for. n_cov_ancova now enters only the confidence-interval degrees of
    # freedom, where it belongs.
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
    # this function). Both are scaled by the SAME factor, so the variance-stabilising
    # relation between them -- z is built so that (dz/dr)^2 * vr = 1/(n-1) under the
    # standard design -- is preserved exactly, whatever the design.
    vr_viechtbauer <- vr_viechtbauer * prec_ratio
    vz_viechtbauer <- vz_viechtbauer * prec_ratio
    # ========= 95% CI ===== #


    z_lo_viechtbauer <- z_viechtbauer - qnorm(.975) * sqrt(vz_viechtbauer)
    z_up_viechtbauer <- z_viechtbauer + qnorm(.975) * sqrt(vz_viechtbauer)
    r_lo_viechtbauer <- (1 / a_viechtbauer) * ((exp(2 * z_lo_viechtbauer / a_viechtbauer) - 1) / (exp(2 * z_lo_viechtbauer / a_viechtbauer) + 1))
    r_up_viechtbauer <- (1 / a_viechtbauer) * ((exp(2 * z_up_viechtbauer / a_viechtbauer) - 1) / (exp(2 * z_up_viechtbauer / a_viechtbauer) + 1))

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
    # delta method for z = atanh(r), r = d / sqrt(d^2 + a) with a = 1/(p*(1-p)):
    # dz/dd = 1/sqrt(d^2 + a), so Var(z) = vd / (d^2 + a). The additive denominator
    # term is the SQUARED point estimate d^2, NOT the sampling variance vd. This makes
    # vz_lipsey the EXACT Fisher transform of vr_lipsey above (vz = vr/(1-r^2)^2).
    # NB: esc::convert_d2r() returns vd/(vd + a) here, which is NOT Fisher-consistent
    # with its own r variance; metaConvert intentionally uses the consistent value, so
    # the lipsey_cooper z-SE differs from esc for large |d| (they agree as d -> 0).
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

    # Negate AND swap each CI on reverse -- see the identical block in
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

    # Negate AND swap each CI on reverse -- see the identical block in
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
#' resulting SMD is NA rather than a silent Inf/NaN. Scoped deliberately to the
#' DENOMINATOR: the mean-change wrappers legitimately pass mean_pre_sd = 0 (the
#' pre slot is zeroed by construction), so guarding raw inputs would break them,
#' whereas a zero *standardizer* is always degenerate. Reachable when a study
#' reports SD = 0, or when r_pre_post = 1 with equal pre/post SDs makes the
#' change SD collapse to 0.
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
#' Each branch's variance follows metafor's LS *pattern* for the corresponding measure
#' -- an empirical leading term plus a g^2 term -- but note this is a pattern, not a
#' single closed-form rule: metafor puts n (or n1+n2), NOT the standardizer df, in the
#' g^2 denominator of its non-heteroscedastic LS forms, while evaluating J at the
#' standardizer df; the two df deliberately differ. So the g^2 term below is over 2N,
#' and J is J(nu) with nu the standardizer df. Var(d) = Var(g)/J^2 exactly, since
#' g = J*d with J a deterministic constant.
#'
#' The leading term is heteroscedasticity-robust: it is built from the EMPIRICAL pooled
#' change SD, not from Var(change) = 2*sigma^2*(1-r) (valid only when SD_pre = SD_post).
#' For bonett this reduces EXACTLY to Viechtbauer's published two-group pre/post
#' variance vi = 2(1-r)(1/nT + 1/nC) + g^2/(2N) when SD_pre = SD_post within each arm
#' (the df-weighted r_avg makes the reduction exact even when n1 != n2); under
#' heteroscedasticity it departs from that homoscedastic form, which understates the
#' variance (its coverage falls to ~0.86 at SD_pre/SD_post ~ 0.64). d_av takes its
#' fourth-moment g^2 COEFFICIENT from Bonett (2008) eq. 19 but NOT its leading term --
#' see the d_av entry below.
#'
#' Publication status of each pooled variance:
#'   - bonett: reduces to Viechtbauer's published vi (metafor-project Morris-2008 page)
#'             under homoscedasticity; the robust departure is metafor SMCRH / Bonett
#'             (2008) generalized across arms.
#'   - d_z:    exactly metafor::escalc(measure = "SMD", vtype = "LS") on the change
#'             scores (Hedges 1981) -- a published two-group variance.
#'   - d_rm:   d_rm = d_z * sqrt(2(1-r)) (Caldwell & Vigotsky 2020 eq. 13, a definition),
#'             so Var(d_rm) = 2(1-r)*Var(d_z) for known r.
#'   - d_av:   PARTLY Bonett (2008) eq. 19 (two-group mixed design, all-four-SD
#'             standardizer). Only the fourth-moment g^2 coefficient is eq. 19's, and it
#'             is reproduced exactly: eq. 19's first bracket
#'             [(s1^4+s2^4+2 r12^2 s1^2 s2^2)/df1 + (s3^4+s4^4+2 r34^2 s3^2 s4^2)/df2]
#'             /(32 s^4) is identical to g2_coef below (checked to 1e-14 over a grid of
#'             arm sizes). The LEADING term is NOT eq. 19's. Eq. 19 uses the df-based
#'             sum Sc1^2/(n1-1) + Sc2^2/(n2-1); this code uses metafor's n-based pooled
#'             form (sd_change_pooled^2/sd_pooled^2) * N/(n1*n2), i.e. the two-sample
#'             SMD leading term. Under homoscedasticity the ratio of the two is
#'             (1/(n1-1) + 1/(n2-1)) / (1/n1 + 1/n2), so eq. 19 is LARGER by ~2% at
#'             n = 50/50, ~3% at 30/30, ~10% at 12/11 and ~32% at 100/4 -- with
#'             heteroscedastic SDs the departure reaches ~38% at n = 100/10 and is
#'             unbounded as min(n1, n2) -> 2. The n-based form is the conventional and
#'             more conservative (smaller-variance) choice, but at strongly unequal arm
#'             sizes it is not eq. 19 and should not be described as such.
#' The pooled forms without a directly published two-group source (the robust
#' bonett departure, d_rm's algebra, and d_av's leading term) are Monte-Carlo calibrated
#' in tests_save/checked/test-pooled-variance-calibration.R, whose coverage grid is
#' n1, n2 in {10, 29, 90} -- imbalance up to 9:1, where the departure above is ~10%.
#' Ratios beyond that (and small min(n) with a large imbalance) are uncalibrated.
#'
#' CIs use qt(.975, m) with m = N - 2 the pooled standardizer df, for EVERY branch --
#' including d_av, whose Cousineau effective df nu = 2m/(1 + r2_avg) feeds ONLY the bias
#' correction J, not the interval (the single-group kernel follows the same split: its d_av
#' CI is on n - 1, not on its Cousineau df). qt is a small-sample-conservative choice (it
#' slightly over-covers at n ~ 10); metafor and Bonett (2008) use a normal (z) critical value.
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
  # what makes the robust bonett/dz variances reduce EXACTLY to Viechtbauer's published
  # vi = 2(1-r)(1/nT + 1/nC) + g^2/(2N) under homoscedasticity when the arms' SDs are
  # equal but n1 != n2 (an n-weighted mean leaves a residual there).
  r_avg <- ((n_exp - 1) * r_pre_post_exp + (n_nexp - 1) * r_pre_post_nexp) / m

  # Pooled change SD: the empirical scale of the numerator. Used by every branch,
  # so that no branch assumes SD_pre = SD_post.
  sd_change_exp <- sqrt(mean_pre_sd_exp^2 + mean_sd_exp^2 -
                        2 * r_pre_post_exp * mean_pre_sd_exp * mean_sd_exp)
  sd_change_nexp <- sqrt(mean_pre_sd_nexp^2 + mean_sd_nexp^2 -
                         2 * r_pre_post_nexp * mean_pre_sd_nexp * mean_sd_nexp)
  sd_change_pooled <- sqrt(((n_exp - 1) * sd_change_exp^2 +
                            (n_nexp - 1) * sd_change_nexp^2) / m)

  # Var(mean_diff) on the change scale, in units of sd_change_pooled^2
  T_change <- N / (n_exp * n_nexp)

  if (pre_post_to_smd == "bonett") {
    # Morris (2008) d_ppc2: numerator = difference in mean change, standardizer =
    # BASELINE SD pooled across arms (eq. 8-9). Variance is the SMCRH form: the
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
    # is metafor::escalc(measure = "SMD", vtype = "LS") -- no r, no 2(1-r) term.
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
    # deterministic rescaling, NOT a homoscedasticity assumption.
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
    # standardized by all four cell SDs) == metafor SMCRPH carried across two
    # independent arms: a robust change-SD leading term plus a per-arm fourth-moment
    # g^2 term. This is NOT the homoscedastic SMCRP coefficient (1+r^2)/(4N), which is
    # only its SD_pre = SD_post special case.
    sd_av_exp <- sqrt((mean_pre_sd_exp^2 + mean_sd_exp^2) / 2)
    sd_av_nexp <- sqrt((mean_pre_sd_nexp^2 + mean_sd_nexp^2) / 2)
    sd_pooled <- .guard_standardizer(sqrt(((n_exp - 1) * sd_av_exp^2 +
                                           (n_nexp - 1) * sd_av_nexp^2) / m))
    # nu is used ONLY for the bias correction J. The (1+r^2) effective df is
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

  # CI on the pooled standardizer df (N - 2 = m) for EVERY branch. dav's Cousineau
  # effective df nu = 2m/(1 + r2_avg) feeds ONLY the bias correction J (as in the
  # single-group kernel, whose dav CI is likewise on n - 1, not on its Cousineau df) --
  # the variance is built on m, so the interval is too. For bonett/dz/drm nu == m, so
  # this changes nothing there; it only homogenises the dav path.
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

  # A negative SD is invalid input. NA it out so the ES is NA rather than a silently
  # corrupted one: a negative raw SD flips the sign of the -2*r*sd_pre*sd_post cross
  # term in sd_diff, yielding a wrong-but-finite standardizer that the .guard_standardizer
  # on the assembled denominator cannot catch. Zero is allowed: the mean-change wrappers
  # legitimately pass mean_pre_sd = 0 (the pre slot is zeroed by construction).
  mean_pre_sd  <- ifelse(is.finite(mean_pre_sd)  & mean_pre_sd  >= 0, mean_pre_sd,  NA_real_)
  mean_post_sd <- ifelse(is.finite(mean_post_sd) & mean_post_sd >= 0, mean_post_sd, NA_real_)

  # An arm needs n >= 2 for a variance/CI to exist. NA the standardizing SDs below 2 so
  # d, var and CI are ALL NA. morris_drm's var_d = 2*(1-r)/n + d^2/(2n) carries no J, so
  # the J(n-1) = NA path that silently nulls the other three branches leaves it finite at
  # n = 1 (and Inf at n = 0); nulling the SDs here closes that leak and matches the
  # pooled kernel's explicit n >= 2 guard.
  bad_n <- !(is.finite(n) & n >= 2)
  mean_pre_sd  <- ifelse(bad_n, NA_real_, mean_pre_sd)
  mean_post_sd <- ifelse(bad_n, NA_real_, mean_post_sd)

  if (pre_post_to_smd == "bonett") {
    # Bonett method: standardize by baseline SD
    # Matches metafor SMCRH (heteroscedastic-robust variance formula)
    J <- .d_j(n - 1)

    var_change <- mean_pre_sd^2 + mean_post_sd^2 - 2 * r_pre_post * mean_pre_sd * mean_post_sd

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
    # i.e. the change SD rescaled onto the raw-score metric. This equals metafor's SMCR
    # (baseline/raw-SD standardizer) ONLY under homoscedasticity (SD_pre = SD_post); when
    # they differ the two diverge, so d_rm is defined by the Caldwell rescaling here, not
    # by SMCR. Its variance below is 2(1-r) * Var(d_z), the exact scaling of the SMCC form.
    # nb: the mean_change wrappers pass sd_change through mean_post_sd
    # (mean_pre = 0 and mean_pre_sd = 0, so sd_diff reduces to sd_change)
    J <- .d_j(n - 1)

    sd_diff <- .guard_standardizer(sqrt(mean_pre_sd^2 + mean_post_sd^2 -
                    (2 * r_pre_post * mean_pre_sd * mean_post_sd)))

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
                    2 * r_pre_post * mean_pre_sd * mean_post_sd))

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
    # post SDs. Its sampling variance is metafor SMCRPH == Bonett (2008) eq. 10: a
    # heteroscedasticity-robust change-SD leading term PLUS a fourth-moment g^2 term.
    # This is NOT the homoscedastic SMCRP form 2(1-r)/n + g^2(1+r^2)/(4n), which
    # assumes SD_pre = SD_post and understates var(g) when they differ (its CI
    # coverage falls to ~0.93 as SD_pre/SD_post moves away from 1). Morris (2008,
    # p.384) recommended AGAINST d_av precisely because its variance was unknown to
    # him; Bonett (2008) supplies it and this branch uses it.
    # nu = 2(n-1)/(1+r^2) is used ONLY for the Hedges bias correction J (Cousineau
    # 2020 eq. 2; metafor SMCRP/SMCRPH); the variance itself is on df = n-1.
    mi <- 2 * (n - 1) / (1 + r_pre_post^2)
    J <- .d_j(mi)

    sd_av <- .guard_standardizer(sqrt((mean_pre_sd^2 + mean_post_sd^2) / 2))

    d <- (mean_post - mean_pre) / sd_av
    g <- d * J

    # metafor SMCRPH "LS" == Bonett (2008) eq. 10. Verified as a FORMULA IDENTITY:
    # 0 difference against metafor::escalc(measure = "SMCRPH"). It does NOT reproduce
    # the var = 0.0148 printed for Bonett's worked Example 2 (n = 60): this branch
    # gives 0.0146674131, and re-doing Bonett's own arithmetic for that example gives
    # 0.0147449109, so 0.0148 is a rounding slip in the paper rather than a package
    # error. Recomputing eq. 10 by hand with the bias-corrected g used here returns
    # 0.0146674131, identical to 12 digits.
    sd_diff2 <- mean_pre_sd^2 + mean_post_sd^2 - 2 * r_pre_post * mean_pre_sd * mean_post_sd
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
# .cor_to_smd() is scalar and was reached through a per-row mapply() at four call
# sites. Its "viechtbauer" branch -- the DEFAULT for cor_to_smd -- builds a fresh
# metafor::conv.delta() call per row, although conv.delta is itself vectorised:
# 18.1 s per 10,000 rows versus 0.34 s for a single vectorised call.
#
# Memoisation (.mapply_memo) does not help here the way it does on the tetrachoric
# path: correlations are continuous, so rows are effectively all distinct. The fix
# has to be real vectorisation. Rows are grouped by their cor_to_smd value (it is a
# per-row column) and each group is computed in one shot, then reassembled in the
# original order.
#
# Agreement with the per-row path: all three branches agree to within 1-2 ULP. The
# "viechtbauer" d is bit-identical on every row, but its SE differs by 1 ULP on ~0.02%
# of rows (44 / 200,000, max |rel| 2.2e-16), because conv.delta's numerical derivative
# is evaluated at slightly different precision in a vector call. The closed-form
# "cooper" and "mathur" branches likewise differ by 1-2 ULP (max |rel| 4.3e-16, ~0.02%
# of rows) because R evaluates length-1 and length-n arithmetic on different code paths.
#
# NB all arguments must be per-row vectors of the same length; unlike mapply this
# does NOT recycle a scalar. Every call site passes rep(x, length.out = n).
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
    # transf.rtod(r_hat) lands in the g slot, NOT the d slot, and this is
    # deliberate. It is tempting to reason that transf.rtod is a population map
    # with no df term, so its output "must" be an uncorrected d -- but which slot
    # it belongs in is a question about ESTIMATOR BIAS, not about the algebra of
    # the transform.
    #
    # Checked by simulation (bivariate normal, rho = 0.5, median split, true
    # delta = 0.870126, nrep = 2e5), bias of the reported g against delta:
    #
    #   n     Hedges g from raw data    g = transf.rtod(r_hat)   g = J*transf.rtod
    #   10          -0.004                    +0.033                  -0.055
    #   25          -0.000                    +0.012                  -0.017
    #  200          -0.000                    +0.002                  -0.002
    #
    # Putting transf.rtod in g is closer to unbiased at every n. The reason is
    # that this route's small-sample bias is not the Hedges bias: E[r_hat] is
    # biased DOWN, which partly cancels the upward bias the pooled-SD denominator
    # induces in a directly computed d. Applying J on top over-corrects.
    # Neither convention is exactly unbiased -- the residual is the uncancelled
    # remainder, and it is small (<0.02 for n >= 25).
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

