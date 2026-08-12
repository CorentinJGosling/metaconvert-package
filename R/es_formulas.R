# ---------------------------------------------------------------------------
# Sensitivity of effect sizes to the choice of conversion formula.
#
# metaConvert addresses two distinct sources of variability, which require
# differentiation:
#
#   1. The SOURCE STATISTICS -- the reported quantities from which an effect
#      size is derived (means and standard deviations, a t statistic, a 2x2
#      table, ...). For a given estimand, alternative derivations are
#      algebraically equivalent; discrepancies therefore indicate an extraction
#      or primary reporting error. That is the domain of the metaDETECT
#      [DISCORDANT] checks, exposed by convert_df(main_es = FALSE).
#
#   2. The CONVERSION FORMULA -- the analytical transformation applied once the
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
# are omitted, since there is no alternative to evaluate.
.formula_parameters <- function() {
  list(
    or_to_rr        = c("metaumbrella_cases", "metaumbrella_exp", "transpose",
                        "grant", "dipietrantonj"),
    rr_to_or        = c("metaumbrella", "transpose", "grant", "dipietrantonj"),
    or_to_cor       = c("pearson", "digby", "bonett", "lipsey_cooper"),
    cor_to_smd      = c("viechtbauer", "cooper", "mathur"),
    smd_to_cor      = c("viechtbauer", "lipsey_cooper"),
    pre_post_to_smd = c("bonett", "morris_dz", "morris_drm", "morris_dav", "cooper"),
    smd_denom       = c("pooled", "glass", "glass_robust"),
    smd_var         = c("borenstein", "hedges_olkin"),
    prop_to_es      = c("raw", "logit", "freeman_tukey"),
    alpha_to_es     = c("bonett", "raw"),
    icc_to_es       = c("bonett", "raw")
  )
}

# Parameters whose options alter the ANALYSIS SCALE rather than the formula
# applied on a fixed scale. alpha_to_es = "bonett" returns ln(1 - alpha) whereas
# "raw" returns alpha itself; prop_to_es = "logit" returns a log-odds and "raw" a
# proportion. Direct mathematical comparison between these scales is invalid, so
# the values are reported side by side but `deviation` and `spread` are returned
# as NA: the choice determines which quantity is meta-analysed, not how precisely
# a single quantity is computed.
.scale_changing_parameters <- function() c("alpha_to_es", "icc_to_es", "prop_to_es")

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
    k <- c("smd_to_cor", "or_to_cor", "cor_to_smd", "smd_var", "pre_post_to_smd")
  } else if (measure %in% ratio) {
    k <- c("or_to_rr", "rr_to_or", "or_to_cor", "smd_to_cor")
  } else if (measure == "prop") {
    k <- "prop_to_es"
  } else if (measure == "alpha") {
    k <- "alpha_to_es"
  } else if (measure == "icc") {
    k <- "icc_to_es"
  }
  intersect(k, names(.formula_parameters()))
}

# Re-evaluates convert_df() with a single conversion parameter modified. A
# per-row column of the same name takes precedence over the scalar argument
# within convert_df(), so that column must be modified as well; otherwise the
# user-specified value would be disregarded.
.run_with_parameter <- function(raw_data, cargs, parameter, value, measure, exp,
                                es_selected, hierarchy) {
  dat <- raw_data
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
#' ratio into a risk ratio (\code{or_to_rr}), or the five pre-post
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
#'  \code{pre_post_to_smd} \tab bonett, cooper, morris_dz, morris_drm, morris_dav\cr
#'  \code{smd_denom} \tab pooled, glass, glass_robust\cr
#'  \code{smd_var} \tab borenstein, hedges_olkin\cr
#'  \code{prop_to_es} \tab raw, logit, freeman_tukey\cr
#'  \code{alpha_to_es} \tab bonett, raw\cr
#'  \code{icc_to_es} \tab bonett, raw\cr
#' }
#'
#' Each alternative estimate is obtained by re-evaluating \code{\link{convert_df}}
#' with a single conversion parameter modified.
#'
#' For each comparison, \code{\link{convert_df}} returns one effect size, derived
#' from one type of input data -- the type selected by the hierarchy and reported
#' in the \code{info_used} column. A conversion parameter can therefore change the
#' result only if that particular type of input data is converted through it, and
#' parameters leaving every estimate unchanged are not reported. For example,
#' \code{pre_post_to_smd} governs the conversion of pre-post data: if a comparison
#' also reports endpoint means and standard deviations, and if the hierarchy selects
#' \code{means_sd}, then no pre-post conversion is performed and
#' \code{pre_post_to_smd} is absent from the output. To examine its influence,
#' request the pre-post method explicitly through the \code{es_selected} and
#' \code{hierarchy} arguments of \code{\link{convert_df}}.
#'
#' Parameters whose options alter the analysis scale rather than the formula
#' applied on a fixed scale (\code{alpha_to_es}, \code{icc_to_es},
#' \code{prop_to_es}) are reported with \code{deviation} and \code{spread} set to
#' \code{NA}. For these parameters, \code{"bonett"} returns a transformed
#' quantity, such as \eqn{\ln(1 - \alpha)}, whereas \code{"raw"} returns the
#' coefficient itself, so direct mathematical comparison between the two is
#' invalid.
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
#'  \code{is_active} \tab a logical value indicating whether this formula was applied in the original \code{convert_df()} call.\cr
#'  \code{info_used} \tab type of source statistics from which the estimate was derived.\cr
#'  \code{es}, \code{se}, \code{es_ci_lo}, \code{es_ci_up} \tab the effect size, its standard error and its 95% confidence interval under this formula.\cr
#'  \code{deviation} \tab the difference between \code{es} and the effect size obtained under the formula applied in the original call.\cr
#'  \code{scale_change} \tab a logical value indicating whether the options of this parameter lie on different analysis scales.\cr
#' }
#' The range of estimates spanned by each (comparison, parameter) combination is
#' attached as the \code{"summary"} attribute.
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
#' dat <- data.frame(or = 2.5, or_ci_lo = 1.6, or_ci_up = 3.9,
#'                   n_exp = 100, n_nexp = 100, baseline_risk = 0.2)
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
  for (pm in parameters) {
    opts <- all_parameters[[pm]]
    active <- as.character(cargs[[pm]])[1]
    # map the recorded value onto the option set (aliases: control -> glass)
    if (!active %in% opts) {
      alias <- c(control = "glass", control_robust = "glass_robust",
                 viechtbauer = "hedges_olkin")
      if (!is.na(alias[active]) && alias[active] %in% opts) active <- alias[active]
    }
    runs <- list()
    for (op in opts) {
      r <- .run_with_parameter(raw_data, cargs, pm, op, measure, exp,
                               es_selected, hierarchy)
      if (!is.null(r)) runs[[op]] <- r
    }
    if (length(runs) < 2) next
    ref <- if (active %in% names(runs)) runs[[active]] else runs[[1]]
    ref_es <- suppressWarnings(as.numeric(as.character(ref$es)))

    param_rows <- do.call(rbind, lapply(names(runs), function(op) {
      r <- runs[[op]]
      es_v <- suppressWarnings(as.numeric(as.character(r$es)))
      data.frame(
        row_id   = r$row_id,
        study_id = study_id[match(r$row_id, raw_data$row_id)],
        parameter = pm,
        formula  = op,
        is_active = identical(op, active),
        info_used = as.character(r$info_used),
        es  = round(es_v, digits),
        se  = round(suppressWarnings(as.numeric(as.character(r$se))), digits),
        es_ci_lo = round(suppressWarnings(as.numeric(as.character(r$es_ci_lo))), digits),
        es_ci_up = round(suppressWarnings(as.numeric(as.character(r$es_ci_up))), digits),
        deviation = if (pm %in% .scale_changing_parameters()) NA_real_ else
          round(es_v - ref_es[match(r$row_id, ref$row_id)], digits),
        scale_change = pm %in% .scale_changing_parameters(),
        stringsAsFactors = FALSE
      )
    }))
    # Parameters that leave every estimate unchanged are excluded: an identical
    # column across all formulas carries no information. Scale-changing
    # parameters are always retained, since their options are not directly
    # comparable and a numerical test of equivalence does not apply to them.
    # A comparison that no formula could estimate yields an all-NA block, for
    # which range() would return c(Inf, -Inf) with a warning; that case is
    # therefore handled separately.
    spread <- tapply(param_rows$es, param_rows$row_id, function(v) {
      v <- v[is.finite(v)]
      if (length(v) == 0) return(0)
      diff(range(v))
    })
    spread[!is.finite(spread)] <- 0
    if (!pm %in% .scale_changing_parameters() &&
        all(spread <= 10^(-digits), na.rm = TRUE)) next
    keep <- if (pm %in% .scale_changing_parameters()) names(spread) else
      names(spread)[spread > 10^(-digits)]
    blocks[[pm]] <- param_rows[as.character(param_rows$row_id) %in% keep, , drop = FALSE]
  }

  if (length(blocks) == 0) {
    if (verbose) {
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
    act <- out$formula[idx][out$is_active[idx]]
    sc <- isTRUE(out$scale_change[idx][1])
    finite_v <- v[is.finite(v)]
    data.frame(
      row_id  = out$row_id[idx][1],
      study_id = out$study_id[idx][1],
      parameter = out$parameter[idx][1],
      n_formulas = length(idx),
      active  = if (length(act) == 1) act else NA_character_,
      # NA for scale-changing parameters: a minimum, maximum or range computed
      # across different analysis scales would not be interpretable.
      es_min  = if (sc || !length(finite_v)) NA_real_ else round(min(finite_v), digits),
      es_max  = if (sc || !length(finite_v)) NA_real_ else round(max(finite_v), digits),
      spread  = if (sc || !length(finite_v)) NA_real_ else round(diff(range(finite_v)), digits),
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
    bits <- vapply(idx, function(i) {
      if (isTRUE(smry$scale_change[i]) || is.na(smry$spread[i])) {
        paste0(smry$parameter[i], " offers ", smry$n_formulas[i],
               " options lying on different analysis scales (applied: '",
               smry$active[i], "')")
      } else {
        paste0(smry$parameter[i], " ranges from ", smry$es_min[i], " to ",
               smry$es_max[i], " across ", smry$n_formulas[i],
               " formulas (applied: '", smry$active[i], "')")
      }
    }, character(1))
    # contains no "; ": .flag_es_quality() splits merged messages on that string
    paste0("[FORMULA] The estimate depends on the choice of conversion formula: ",
           paste(bits, collapse = " / "),
           " - a methodological choice of the analyst, not a data inconsistency")
  }, character(1))
}
