# ---------------------------------------------------------------------------
# Sensitivity of effect sizes to the choice of conversion formula.
#
# metaConvert addresses two distinct sources of variability, which require
# differentiation:
#
#   1. The source statistics: the reported quantities from which an effect
#      size is derived (means and standard deviations, a t statistic, a 2x2
#      table, ...). For a given estimand, alternative derivations are
#      algebraically equivalent; discrepancies therefore indicate an extraction
#      or primary reporting error. That is the domain of the metaDETECT
#      [DISCORDANT] checks, exposed by convert_df(main_es = FALSE).
#
#   2. The conversion formula: the analytical transformation applied once the
#      source statistics have been read (or_to_rr, cor_to_smd, pre_post_to_smd,
#      ...). Alternative formulas encode different statistical assumptions and
#      are therefore not algebraically equivalent: discrepancies reflect the
#      methodological assumptions of the analyst rather than inconsistencies in
#      the underlying data. This source of structural uncertainty was previously
#      unaccounted for: each convert_df() execution applied a single
#      user-specified value per conversion parameter and returned a point
#      estimate, without indicating how far that estimate depended on the choice.
#
# es_formulas() addresses the second source. It re-evaluates the conversion under
# every alternative formula and returns one row per (comparison, parameter,
# formula). Results are labelled [FORMULA] rather than [DISCORDANT] so that they
# cannot be interpreted as an error signal.
# ---------------------------------------------------------------------------

# Conversion parameters that admit more than one formula, with their complete
# option sets. Parameters with a single tolerated value (e.g. table_2x2_to_cor)
# are omitted, since there is no alternative to evaluate. Options documented
# under two names are listed once only: pre_post_to_smd = "cooper" is rewritten
# to "morris_drm" before any computation (internal_multiple_formulas.R), so
# listing both would return the same estimate twice.
.formula_parameters <- function() {
  list(
    or_to_rr        = c("metaumbrella_cases", "metaumbrella_exp", "transpose",
                        "grant", "dipietrantonj"),
    rr_to_or        = c("metaumbrella", "transpose", "grant", "dipietrantonj"),
    or_to_cor       = c("pearson", "digby", "bonett", "lipsey_cooper"),
    cor_to_smd      = c("viechtbauer", "cooper", "mathur"),
    smd_to_cor      = c("viechtbauer", "lipsey_cooper"),
    pre_post_to_smd = c("bonett", "morris_dz", "morris_drm", "morris_dav"),
    smd_denom       = c("pooled", "glass", "glass_robust"),
    smd_var         = c("borenstein", "hedges_olkin"),
    prop_to_es      = c("raw", "logit", "freeman_tukey"),
    alpha_to_es     = c("bonett", "raw", "hakstian_whalen"),
    icc_to_es       = c("bonett", "raw"),
    omega_to_es     = c("bonett", "raw", "hakstian_whalen")
  )
}

# Parameters whose options alter the analysis scale rather than the formula
# applied on a fixed scale. alpha_to_es = "bonett" returns ln(1 - alpha) whereas
# "raw" returns alpha itself; prop_to_es = "logit" returns a log-odds and "raw" a
# proportion. Direct mathematical comparison between these scales is invalid, so
# the values are reported side by side but `deviation` and `spread` are returned
# as NA: the choice determines which quantity is meta-analysed, not how precisely
# a single quantity is computed.
.scale_changing_parameters <- function() c("alpha_to_es", "icc_to_es", "prop_to_es",
                                           "omega_to_es")

# Identifies the conversion parameters relevant to a specific target effect size.
# Iterating over irrelevant parameters (e.g. varying alpha_to_es for Hedges' g)
# would trigger redundant computations that reproduce identical values, so these
# are excluded a priori.
.parameters_for_measure <- function(measure) {
  smd_like <- c("d", "g", "md", "dw", "gw", "mdw")
  cor_like <- c("r", "z", "rp", "zp")
  ratio    <- c("logor", "logrr", "logirr", "loghr", "nnt", "rd")
  k <- character(0)
  if (measure %in% smd_like) {
    k <- c("pre_post_to_smd", "smd_denom", "smd_var", "cor_to_smd",
           "or_to_cor", "or_to_rr", "rr_to_or")
  } else if (measure %in% cor_like) {
    # smd_denom, or_to_rr and rr_to_or reach a correlation too: an SMD standardized
    # on the control arm gives a different point-biserial r, and a reported RR is
    # routed through the OR before or_to_cor sees it.
    k <- c("smd_to_cor", "or_to_cor", "cor_to_smd", "smd_var", "pre_post_to_smd",
           "smd_denom", "or_to_rr", "rr_to_or")
  } else if (measure %in% ratio) {
    # A correlation or a pre/post design reaches an odds ratio through the SMD, so
    # the parameters governing that leg move the ratio as well: on a single
    # pearson_r row cor_to_smd moves the OR by 53%, and on a paired-t row
    # pre_post_to_smd moves it by 8.5%. Over-inclusion is free -- a parameter is
    # reported only where the estimates actually differ -- whereas omitting one
    # prints an affirmative all-clear.
    k <- c("or_to_rr", "rr_to_or", "or_to_cor", "smd_to_cor",
           "cor_to_smd", "pre_post_to_smd", "smd_denom", "smd_var")
  } else if (measure == "prop") {
    k <- "prop_to_es"
  } else if (measure == "alpha") {
    k <- "alpha_to_es"
  } else if (measure == "icc") {
    k <- "icc_to_es"
  } else if (measure == "omega") {
    k <- "omega_to_es"
  }
  intersect(k, names(.formula_parameters()))
}

# Columns appended by .check_data() to hold the results of an estimation. The
# "raw_data" attribute is the checked frame, so it already carries them; feeding
# it back to convert_df() would append a second copy under the same names and
# the estimation would abort on the duplicated names. They are therefore removed
# before every re-evaluation. All of them are empty in "raw_data", so nothing
# supplied by the user is lost.
.formula_result_columns <- function() {
  base <- c(
    "all_info", "measure", "info_measure", "n_estimations", "es_selected", "info_used",
    "es", "se", "es_ci_lo", "es_ci_up",
    "min_info", "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
    "max_info", "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up",
    "diff_min_max", "overlap_min_max", "dispersion_es"
  )
  c("adjusted_input", base, paste0(base, "_crude"), paste0(base, "_adjusted"))
}

# Maps a value recorded by convert_df() onto the option set of a conversion
# parameter. Several options are documented under two names ("control" is
# Glass's delta, "viechtbauer" the Hedges & Olkin variance); unname() is required
# because the lookup returns a named element, which does not compare equal to
# the plain option name.
.resolve_formula_alias <- function(value, opts) {
  alias <- c(control = "glass", control_robust = "glass_robust",
             viechtbauer = "hedges_olkin", cooper = "morris_drm")
  value <- as.character(value)
  hit <- !is.na(value) & !value %in% opts & value %in% names(alias)
  value[hit] <- unname(alias[value[hit]])
  value[!is.na(value) & !value %in% opts] <- NA_character_
  value
}

# The formula applied by the original convert_df() call, one value per
# comparison. convert_df() reads a per-row column named after the parameter in
# preference to its scalar argument (main_convert_df.R), so the applied formula
# is a property of the row, not of the call.
.resolve_formula <- function(parameter, cargs, raw_data, opts) {
  scalar <- .resolve_formula_alias(as.character(cargs[[parameter]])[1], opts)
  out <- rep(scalar, nrow(raw_data))
  if (parameter %in% colnames(raw_data)) {
    col <- .resolve_formula_alias(raw_data[[parameter]], opts)
    out[!is.na(col)] <- col[!is.na(col)]
  }
  out
}

# Re-evaluates convert_df() with a single conversion parameter modified. A
# per-row column of the same name takes precedence over the scalar argument
# within convert_df(), so that column must be modified as well; otherwise the
# user-specified value would be disregarded.
.run_with_parameter <- function(raw_data, cargs, parameter, value, measure, exp,
                                es_selected, hierarchy) {
  dat <- raw_data[, setdiff(colnames(raw_data), .formula_result_columns()),
                  drop = FALSE]
  if (parameter %in% colnames(dat)) dat[[parameter]] <- value
  args <- cargs[intersect(names(cargs), names(formals(convert_df)))]
  args[[parameter]] <- value
  args$x <- dat
  args$measure <- if (exp && measure == "logor") "or"
                  else if (exp && measure == "logrr") "rr"
                  else if (exp && measure == "logirr") "irr"
                  else if (exp && measure == "loghr") "hr"
                  else measure
  args$verbose <- FALSE
  args$main_es <- TRUE
  args$split_adjusted <- FALSE
  args$es_selected <- es_selected
  args$hierarchy <- hierarchy
  out <- try(suppressWarnings(suppressMessages(
    summary(do.call(convert_df, args), flags = FALSE, guidance = FALSE,
            include_raw = FALSE))), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

#' Sensitivity of effect sizes to alternative conversion formulas
#'
#' @description
#' Recomputes the effect sizes of a \dQuote{metaConvert} object under all
#' alternative conversion formulas, and returns the results side by side, with
#' one row per comparison per formula.
#'
#' Whereas \code{\link{convert_df}(main_es = FALSE)} evaluates variability in the
#' \emph{source statistics} from which an effect size is derived,
#' \code{es_formulas()} evaluates the uncertainty of the
#' \emph{analytical conversion} applied to those statistics. For example, it
#' assesses the five methodological approaches available for converting an odds
#' ratio into a risk ratio (\code{or_to_rr}), or the four pre-post
#' standardisation methods (\code{pre_post_to_smd}). Because alternative formulas
#' encode different statistical assumptions rather than algebraic equivalencies,
#' discrepancies among them constitute methodological uncertainty rather than
#' data extraction errors. Results are therefore labelled \code{[FORMULA]} rather
#' than with the metaDETECT \code{[DISCORDANT]} flag.
#'
#' @param object an object of class \dQuote{metaConvert} produced by
#'   \code{\link{convert_df}}.
#' @param parameters a character vector specifying which conversion parameters to
#'   vary. The default (\code{NULL}) evaluates every conversion parameter that can
#'   influence the requested effect size measure. See Details for the available
#'   parameters.
#' @param digits an integer specifying the number of decimal places used for
#'   rounding. Default is 3.
#' @param verbose a logical value indicating whether a brief computation summary
#'   should be printed. Default is TRUE.
#'
#' @details
#' The conversion parameters that support several methodological formulas are:
#' \tabular{ll}{
#'  \code{or_to_rr} \tab metaumbrella_cases, metaumbrella_exp, transpose, grant, dipietrantonj\cr
#'  \code{rr_to_or} \tab metaumbrella, transpose, grant, dipietrantonj\cr
#'  \code{or_to_cor} \tab pearson, digby, bonett, lipsey_cooper\cr
#'  \code{cor_to_smd} \tab viechtbauer, cooper, mathur\cr
#'  \code{smd_to_cor} \tab viechtbauer, lipsey_cooper\cr
#'  \code{pre_post_to_smd} \tab bonett, morris_dz, morris_drm (alias cooper), morris_dav\cr
#'  \code{smd_denom} \tab pooled, glass (alias control), glass_robust (alias control_robust)\cr
#'  \code{smd_var} \tab borenstein, hedges_olkin (alias viechtbauer)\cr
#'  \code{prop_to_es} \tab raw, logit, freeman_tukey\cr
#'  \code{alpha_to_es} \tab bonett, raw, hakstian_whalen\cr
#'  \code{icc_to_es} \tab bonett, raw\cr
#'  \code{omega_to_es} \tab bonett, raw, hakstian_whalen\cr
#' }
#' Options documented under two names are evaluated once, under the name listed
#' above.
#'
#' Each alternative estimate is obtained by re-evaluating \code{\link{convert_df}}
#' with a single conversion parameter modified. When the original call supplied a
#' conversion parameter as a column of the dataset rather than as a single value,
#' the \code{deviation} of each comparison is measured against the formula that
#' comparison itself received.
#'
#' For each comparison, \code{\link{convert_df}} returns one effect size, derived
#' from one type of input data, the type selected by the hierarchy and reported
#' in the \code{info_used} column. A conversion parameter can therefore change the
#' result only if that particular type of input data is converted through it. A
#' comparison is reported only if the parameter reached it: at least one formula
#' must have produced an estimate, and the formulas must differ in the effect
#' size or in its standard error. For example,
#' \code{pre_post_to_smd} governs the conversion of pre-post data: if a comparison
#' also reports endpoint means and standard deviations, and if the hierarchy selects
#' \code{means_sd}, then no pre-post conversion is performed and
#' \code{pre_post_to_smd} is absent from the output. To examine its influence,
#' request the pre-post method explicitly through the \code{es_selected} and
#' \code{hierarchy} arguments of \code{\link{convert_df}}.
#'
#' Some options leave the effect size untouched and alter only its standard error,
#' and hence the weight the comparison receives in a meta-analysis:
#' \code{smd_var} selects the sampling-variance formula, and \code{glass} and
#' \code{glass_robust} share a standardizer but not a variance. These comparisons
#' are reported with a \code{deviation} of zero, and the range of standard errors
#' is given by the \code{se_min} and \code{se_max} columns of the
#' \code{"summary"} attribute.
#'
#' Parameters whose options alter the analysis scale rather than the formula
#' applied on a fixed scale (\code{alpha_to_es}, \code{omega_to_es},
#' \code{icc_to_es}, \code{prop_to_es}) are reported with \code{deviation} set to
#' \code{NA}, as are the \code{es_min}, \code{es_max} and \code{spread} columns of the
#' \code{"summary"} attribute. For these parameters, \code{"bonett"} returns a
#' transformed quantity, such as \eqn{\ln(1 - \alpha)}, whereas \code{"raw"}
#' returns the coefficient itself, so direct mathematical comparison between the
#' two is invalid.
#'
#' For the ratio measures reported on their natural scale (\code{or}, \code{rr},
#' \code{irr}, \code{hr}), \code{deviation} and the \code{spread} column of the
#' \code{"summary"} attribute are differences of ratios, not of their logarithms.
#' They are comparable across the formulas of one comparison, but not across
#' comparisons of different magnitude: a given ratio of two risk ratios yields a
#' larger difference the larger the risk ratios are.
#'
#' \strong{Warning regarding meta-analytic pooling.} The rows returned by this
#' function must never be pooled in a meta-analysis. They are structurally
#' dependent re-evaluations of identical primary study data under different
#' mathematical assumptions, not independent effect size estimates.
#'
#' @return
#' A dataframe of class \code{metaConvert_formulas} containing one row per
#' (comparison, parameter, formula) combination:
#' \tabular{ll}{
#'  \code{row_id} \tab row index in the original dataset.\cr
#'  \code{study_id} \tab study identifier, when available.\cr
#'  \code{parameter} \tab the conversion parameter being varied.\cr
#'  \code{formula} \tab the formula applied on this row.\cr
#'  \code{is_active} \tab a logical value indicating whether this formula was the one applied to this comparison in the original \code{convert_df()} call.\cr
#'  \code{info_used} \tab type of source statistics from which the estimate was derived.\cr
#'  \code{es}, \code{se}, \code{es_ci_lo}, \code{es_ci_up} \tab the effect size, its standard error and its 95% confidence interval under this formula.\cr
#'  \code{deviation} \tab the difference between \code{es} and the effect size obtained under the formula applied to this comparison in the original call.\cr
#'  \code{scale_change} \tab a logical value indicating whether the options of this parameter lie on different analysis scales.\cr
#' }
#' One row per (comparison, parameter) combination is attached as the
#' \code{"summary"} attribute, ordered by decreasing \code{spread}:
#' \tabular{ll}{
#'  \code{row_id}, \code{study_id}, \code{parameter} \tab as above.\cr
#'  \code{n_formulas} \tab number of formulas evaluated.\cr
#'  \code{n_estimated} \tab number of these formulas that produced an estimate. The remaining ones require input data the comparison does not report.\cr
#'  \code{active} \tab the formula applied to this comparison in the original call.\cr
#'  \code{active_estimated} \tab a logical value indicating whether that formula produced an estimate.\cr
#'  \code{es_min}, \code{es_max}, \code{spread} \tab the smallest and largest effect size, and their difference.\cr
#'  \code{se_min}, \code{se_max}, \code{se_spread} \tab the smallest and largest standard error, and their difference.\cr
#'  \code{scale_change} \tab as above.\cr
#' }
#'
#' @seealso
#' \code{\link{convert_df}} for the estimation of effect sizes\cr
#' \code{\link{summary.metaConvert}} for the main effect size table\cr
#'
#' @export
#'
#' @md
#'
#' @examples
#' \donttest{
#' ### sensitivity of the risk ratio to the odds ratio conversion formula
#' ### the case and control margins are required by two of the five formulas
#' dat <- data.frame(or = 2.5, or_ci_lo = 1.6, or_ci_up = 3.9,
#'                   n_exp = 100, n_nexp = 100,
#'                   n_cases = 60, n_controls = 140, baseline_risk = 0.2)
#' es_formulas(convert_df(dat, measure = "rr", verbose = FALSE))
#' }
es_formulas <- function(object, parameters = NULL, digits = 3, verbose = TRUE) {
  if (!inherits(object, "metaConvert")) {
    stop("'object' must be an object of class 'metaConvert' (see convert_df()).")
  }
  cargs <- attr(object, "conversion_args")
  if (is.null(cargs)) {
    stop("This 'metaConvert' object does not record its conversion parameters, ",
         "as it was created by an earlier version of the package. Please re-run ",
         "convert_df() before calling es_formulas().")
  }
  raw_data <- attr(object, "raw_data")
  measure  <- attr(object, "measure")
  exp      <- attr(object, "exp")
  es_selected <- attr(object, "es_selected")
  hierarchy   <- attr(object, "hierarchy")

  all_parameters <- .formula_parameters()
  user_specified <- !is.null(parameters)
  if (!user_specified) {
    parameters <- .parameters_for_measure(measure)
  } else {
    bad <- setdiff(parameters, names(all_parameters))
    if (length(bad) > 0) {
      stop("Unknown conversion parameter(s): ", paste(bad, collapse = ", "),
           ". Available parameters are: ",
           paste(names(all_parameters), collapse = ", "), ".")
    }
  }
  if (length(parameters) == 0) {
    if (verbose) {
      message("No conversion formula can influence the '", measure,
              "' measure, so no alternative was evaluated.")
    }
    return(.empty_formula_frame())
  }

  study_id <- if ("study_id" %in% colnames(raw_data)) {
    as.character(raw_data$study_id)
  } else rep(NA_character_, nrow(raw_data))

  blocks <- list()
  failed <- character(0)
  unestimable <- character(0)
  for (pm in parameters) {
    opts <- all_parameters[[pm]]
    # The formula applied by the original call, resolved separately for every
    # comparison: convert_df() reads a per-row column of the same name in
    # preference to the scalar argument, so a dataset may apply a different
    # formula to each row.
    active <- .resolve_formula(pm, cargs, raw_data, opts)
    runs <- list()
    for (op in opts) {
      r <- .run_with_parameter(raw_data, cargs, pm, op, measure, exp,
                               es_selected, hierarchy)
      if (!is.null(r)) runs[[op]] <- r
    }
    if (length(runs) < 2) {
      failed <- c(failed, pm)
      next
    }
    # per comparison, the estimate obtained under the formula that comparison
    # actually received: the reference against which `deviation` is measured
    ref_es <- rep(NA_real_, nrow(raw_data))
    for (op in names(runs)) {
      r <- runs[[op]]
      p <- match(r$row_id, raw_data$row_id)
      sel <- !is.na(p) & !is.na(active[p]) & active[p] == op
      if (!any(sel)) next
      ref_es[p[sel]] <- suppressWarnings(as.numeric(as.character(r$es)))[sel]
    }

    param_rows <- do.call(rbind, lapply(names(runs), function(op) {
      r <- runs[[op]]
      p <- match(r$row_id, raw_data$row_id)
      es_v <- suppressWarnings(as.numeric(as.character(r$es)))
      data.frame(
        row_id   = r$row_id,
        study_id = study_id[p],
        parameter = pm,
        formula  = op,
        is_active = !is.na(active[p]) & active[p] == op,
        info_used = as.character(r$info_used),
        es  = round(es_v, digits),
        se  = round(suppressWarnings(as.numeric(as.character(r$se))), digits),
        es_ci_lo = round(suppressWarnings(as.numeric(as.character(r$es_ci_lo))), digits),
        es_ci_up = round(suppressWarnings(as.numeric(as.character(r$es_ci_up))), digits),
        deviation = if (pm %in% .scale_changing_parameters()) NA_real_ else
          round(es_v - ref_es[p], digits),
        scale_change = pm %in% .scale_changing_parameters(),
        stringsAsFactors = FALSE
      )
    }))
    # A comparison is reported only if the parameter reached it: at least one
    # formula must have produced an estimate, and the formulas must not all
    # return the same value. Identical values prove that the parameter played no
    # part in the estimation, whatever the analysis scale of its options, so
    # scale-changing parameters are filtered on the same rule. The comparison is
    # also retained when only the standard error differs, as `smd_var` and
    # `glass` versus `glass_robust` leave the point estimate untouched and alter
    # the weight the comparison receives in a meta-analysis.
    tol <- 10^(-digits)
    rng <- function(v) {
      v <- v[is.finite(v)]
      if (length(v) == 0) return(NA_real_)
      diff(range(v))
    }
    by_row     <- as.character(param_rows$row_id)
    es_spread  <- vapply(split(param_rows$es, by_row), rng, numeric(1))
    se_spread  <- vapply(split(param_rows$se, by_row), rng, numeric(1))
    estimable  <- vapply(split(param_rows$es, by_row),
                         function(v) any(is.finite(v)), logical(1))
    changed <- (!is.na(es_spread) & es_spread > tol) |
               (!is.na(se_spread) & se_spread > tol)
    keep <- names(es_spread)[estimable & changed]
    if (length(keep) == 0) {
      # A parameter under which NO comparison produced an estimate has not been
      # examined either. That is the shape the reliability SE-source switches make:
      # under alpha_se_source = "reported", omega_se_source = "reported" or
      # icc_agreement_se = "drop", a row with no reported SE has es = NA on every
      # alternative, so the filter above drops the block for want of a value rather
      # than for want of a difference. Folding that into the "changed nothing"
      # all-clear below would answer a question that was never asked.
      if (!any(estimable)) unestimable <- c(unestimable, pm)
      next
    }
    blocks[[pm]] <- param_rows[by_row %in% keep, , drop = FALSE]
  }

  # A parameter whose re-evaluations aborted has not been examined; reporting it
  # as having changed nothing would be an unwarranted all-clear.
  if (length(failed) > 0) {
    warning("The alternative formulas of the following conversion parameter(s) ",
            "could not be evaluated, and their influence is therefore unknown: ",
            paste(failed, collapse = ", "), ".", call. = FALSE)
  }
  if (verbose && length(unestimable) > 0) {
    message("No comparison could be estimated under any alternative of: ",
            paste(unestimable, collapse = ", "), ". Their influence is unknown, not ",
            "nil: the analysis recorded in this object returns no effect size for ",
            "these comparisons, which is what alpha_se_source = 'reported', ",
            "omega_se_source = 'reported' and icc_agreement_se = 'drop' are for.")
  }
  if (length(blocks) == 0) {
    if (verbose && length(unestimable) < length(parameters)) {
      message("No conversion formula altered any effect size in this dataset.")
      message("Each comparison is estimated from one type of input data, ",
              "reported in the 'info_used' column of summary(). A conversion ",
              "parameter can change the result only if that type of input data ",
              "is converted through it. For example, 'pre_post_to_smd' governs ",
              "the conversion of pre-post data, so it has no effect on a ",
              "comparison estimated from endpoint means and standard deviations ",
              "(info_used = 'means_sd'). Request a different type of input data ",
              "through the 'es_selected' and 'hierarchy' arguments of convert_df().")
    }
    return(.empty_formula_frame())
  }
  dropped <- setdiff(parameters, names(blocks))
  if (verbose && length(dropped) > 0 && user_specified) {
    message("The following parameter(s) did not influence the estimation ",
            "methods selected in this dataset, and are therefore not reported: ",
            paste(dropped, collapse = ", "), ".")
  }

  out <- do.call(rbind, blocks)
  out <- out[order(out$row_id, out$parameter, !out$is_active, out$formula), ]
  rownames(out) <- NULL

  smry <- .formula_spread_summary(out, digits)
  attr(out, "summary") <- smry
  attr(out, "measure") <- measure
  class(out) <- c("metaConvert_formulas", "data.frame")

  if (verbose) {
    message("\n-- Sensitivity to the choice of conversion formula --")
    message("Conversion parameters evaluated: ",
            paste(names(blocks), collapse = ", "))
    message("Comparisons affected: ", length(unique(smry$row_id)),
            " (", nrow(smry), " comparison-by-parameter combination",
            if (nrow(smry) > 1) "s" else "", ")")
    message("These rows are re-estimations of the same studies under different ",
            "assumptions and must not be pooled.")
  }
  out
}

.empty_formula_frame <- function() {
  out <- data.frame(row_id = numeric(0), study_id = character(0),
                    parameter = character(0), formula = character(0),
                    is_active = logical(0), info_used = character(0),
                    es = numeric(0), se = numeric(0),
                    es_ci_lo = numeric(0), es_ci_up = numeric(0),
                    deviation = numeric(0), scale_change = logical(0),
                    stringsAsFactors = FALSE)
  attr(out, "summary") <- out[0, c("row_id", "parameter")]
  class(out) <- c("metaConvert_formulas", "data.frame")
  out
}

# One row per (comparison, conversion parameter): the range of estimates spanned
# by the alternative formulas.
.formula_spread_summary <- function(out, digits = 3) {
  key <- paste(out$row_id, out$parameter, sep = "\r")
  parts <- split(seq_len(nrow(out)), key)
  res <- do.call(rbind, lapply(parts, function(idx) {
    v <- out$es[idx]
    s <- out$se[idx]
    act <- unique(out$formula[idx][out$is_active[idx]])
    sc <- isTRUE(out$scale_change[idx][1])
    finite_v <- v[is.finite(v)]
    finite_s <- s[is.finite(s)]
    data.frame(
      row_id  = out$row_id[idx][1],
      study_id = out$study_id[idx][1],
      parameter = out$parameter[idx][1],
      n_formulas = length(idx),
      # formulas that produced an estimate: es_min, es_max and spread are
      # computed over these alone, so n_formulas alone overstates the range
      n_estimated = length(finite_v),
      active  = if (length(act) == 1) act else NA_character_,
      active_estimated = length(act) == 1 &&
        any(is.finite(v[out$is_active[idx]])),
      # NA for scale-changing parameters: a minimum, maximum or range computed
      # across different analysis scales would not be interpretable.
      es_min  = if (sc || !length(finite_v)) NA_real_ else round(min(finite_v), digits),
      es_max  = if (sc || !length(finite_v)) NA_real_ else round(max(finite_v), digits),
      spread  = if (sc || !length(finite_v)) NA_real_ else round(diff(range(finite_v)), digits),
      se_min  = if (sc || !length(finite_s)) NA_real_ else round(min(finite_s), digits),
      se_max  = if (sc || !length(finite_s)) NA_real_ else round(max(finite_s), digits),
      se_spread = if (sc || !length(finite_s)) NA_real_ else round(diff(range(finite_s)), digits),
      scale_change = sc,
      stringsAsFactors = FALSE
    )
  }))
  res <- res[order(-ifelse(is.na(res$spread), -Inf, res$spread)), ]
  rownames(res) <- NULL
  res
}

#' Print the sensitivity of effect sizes to alternative conversion formulas
#'
#' @param x an object of class \dQuote{metaConvert_formulas}
#' @param ... other arguments that can be passed to the function
#'
#' @return \code{x}, invisibly. Called for its printed output.
#'
#' @exportS3Method
#'
#' @md
print.metaConvert_formulas <- function(x, ...) {
  if (nrow(x) == 0) {
    cat("No conversion formula altered any effect size.\n")
    return(invisible(x))
  }
  cat("-- Effect size under each alternative conversion formula --\n")
  cat("These rows are re-estimations of the same studies under different\n")
  cat("statistical assumptions. They are not independent estimates and must\n")
  cat("not be pooled in a meta-analysis.\n\n")
  print(utils::head(as.data.frame(x), 40))
  if (nrow(x) > 40) cat("... ", nrow(x) - 40, " further row(s)\n", sep = "")
  smry <- attr(x, "summary")
  if (!is.null(smry) && nrow(smry) > 0) {
    cat("\nComparisons most sensitive to the choice of conversion formula:\n")
    print(utils::head(smry, 10))
  }
  invisible(x)
}

# Compact per-row "[FORMULA] ..." messages, indexed by row_id, used by
# summary(formulas = TRUE).
.formula_tokens <- function(fx, digits = 3) {
  if (is.null(fx) || nrow(fx) == 0) return(NULL)
  smry <- attr(fx, "summary")
  if (is.null(smry) || nrow(smry) == 0) return(NULL)
  parts <- split(seq_len(nrow(smry)), smry$row_id)
  vapply(parts, function(idx) {
    # the formula this comparison actually received; when it produced no
    # estimate the reported range comes entirely from the alternatives
    applied <- vapply(idx, function(i) {
      if (is.na(smry$active[i])) return("")
      paste0(" (applied: '", smry$active[i],
             if (isTRUE(smry$active_estimated[i])) "'" else "', which yields no estimate here",
             ")")
    }, character(1))
    # how many of the formulas produced an estimate
    over <- vapply(idx, function(i) {
      if (smry$n_estimated[i] < smry$n_formulas[i]) {
        paste0(smry$n_estimated[i], " of ", smry$n_formulas[i], " formulas")
      } else {
        paste0(smry$n_formulas[i], " formulas")
      }
    }, character(1))
    es_change <- rep(FALSE, length(idx))
    bits <- vapply(seq_along(idx), function(k) {
      i <- idx[k]
      if (isTRUE(smry$scale_change[i]) || is.na(smry$spread[i])) {
        paste0(smry$parameter[i], " offers ", smry$n_formulas[i],
               " options lying on different analysis scales", applied[k])
      } else if (smry$spread[i] > 0) {
        es_change[k] <<- TRUE
        paste0(smry$parameter[i], " ranges from ", smry$es_min[i], " to ",
               smry$es_max[i], " across ", over[k], applied[k])
      } else {
        # the point estimate is common to every formula and only its precision
        # differs, as with smd_var or glass versus glass_robust
        paste0(smry$parameter[i], " leaves the estimate unchanged but its ",
               "standard error ranges from ", smry$se_min[i], " to ",
               smry$se_max[i], " across ", over[k], applied[k])
      }
    }, character(1))
    lead <- if (any(es_change) || any(smry$scale_change[idx])) {
      "The estimate depends on the choice of conversion formula: "
    } else {
      "The precision of the estimate depends on the choice of conversion formula: "
    }
    # contains no "; ": .flag_es_quality() splits merged messages on that string
    paste0("[FORMULA] ", lead, paste(bits, collapse = " / "),
           " - a methodological choice of the analyst, not a data inconsistency")
  }, character(1))
}
