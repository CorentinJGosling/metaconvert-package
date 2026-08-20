#' Convert McDonald's omega into an effect size measure
#'
#' @param omega McDonald's omega reliability coefficient (natural scale, usually between 0 and 1)
#' @param omega_se standard error of omega, \strong{on the natural (coefficient) scale}
#' @param omega_ci_lo lower bound of the 95% confidence interval of omega (natural scale)
#' @param omega_ci_up upper bound of the 95% confidence interval of omega (natural scale)
#' @param n_sample the total number of participants that completed the scale. Not used to
#'   compute the standard error (see \strong{Details}); carried through for reporting and
#'   moderator analyses.
#' @param n_items number of items in the scale. Carried through for reporting; not used in
#'   the standard error.
#' @param omega_type which omega is reported: \code{"total"} (default), \code{"hierarchical"},
#'   \code{"asymptotic"} or \code{"subscale"}. These are \strong{different estimands} and must
#'   not be pooled together (see \strong{Details}).
#' @param omega_estimator how omega was estimated in the primary study: one of
#'   \code{"cfa_bifactor"}, \code{"cfa_1factor"}, \code{"efa_schmid_leiman"},
#'   \code{"first_pc"}, \code{"first_pf"} or \code{"unspecified"} (default).
#'   Case-insensitive, with the common synonyms accepted. This is a
#'   \strong{provenance} field, not a computational one -- it changes nothing in the
#'   arithmetic, but a pool that mixes estimators is flagged, because the estimator
#'   moves omega far more than most moderators do (see \strong{Details}).
#' @param omega_to_es method used to compute the effect size from omega. One of
#'   \code{"bonett"} (default), \code{"raw"} or \code{"hakstian_whalen"}, matching
#'   \code{\link{es_from_cronbach_alpha}}.
#'
#' @details
#' \strong{Why this function needs a standard error or a confidence interval, and
#' \code{es_from_cronbach_alpha()} does not.} Coefficient alpha has a closed-form
#' sampling variance in \eqn{(n, k)} alone -- Bonett's \eqn{2k/((k-1)(n-2))} and the
#' Hakstian-Whalen variance both descend from the Feldt (1965) / Kristof (1963) result
#' that \eqn{(1 - r)/(1 - \alpha)} is distributed as \eqn{F}, which holds under
#' essential tau-equivalence, compound symmetry and multivariate normality. Omega exists
#' precisely to \emph{drop} tau-equivalence.
#'
#' It would be tidy to say that borrowing the alpha variance is therefore a category
#' error, but that is \strong{not what the numbers show}, and the honest reasons are
#' more specific. Under a correct unidimensional congeneric \emph{normal} model the
#' Bonett variance on the \eqn{\ln(1-\omega)} scale is within \eqn{[0.899, 1.003]} of
#' the full normal-theory asymptotic SE of \eqn{\omega_{total}} (1283 random loading
#' patterns, \eqn{k \ge 4}, \eqn{\omega \ge 0.60}, \eqn{n = 200}), and coincides with
#' it exactly under tau-equivalence. The reasons not to use it are that it is
#' anti-conservative exactly where omega is actually reached for:
#'
#' \itemize{
#'   \item \strong{7.6% and 25.6% narrower} than the two published bootstrap CIs in
#'     Flora (2020) -- inverse-variance weight inflated 1.17x and 1.81x;
#'   \item \strong{13% narrower} under a single correlated residual
#'     (\eqn{r = .30}, \eqn{k = 6}, \eqn{n = 200}), where \eqn{\hat\omega} is also
#'     biased \eqn{+0.028};
#'   \item \strong{14-24% narrower} for a bifactor \eqn{\omega_h}, whose sampling
#'     variance is governed by the number of \emph{group factors} \eqn{m}, not by
#'     \eqn{k} (the normal-theory identity is
#'     \eqn{SE(\ln(1-\omega_h)) = \sqrt{2m/((m-1)n)}});
#'   \item it \strong{diverges at \eqn{k = 3}} (a just-identified one-factor model);
#'   \item and \strong{no published source endorses it} -- every source consulted
#'     declines to give an omega variance at all.
#' }
#'
#' There is also no route through the loadings that adds information: conditioning on a
#' correctly specified model makes Raykov's (2002) delta method evaluable from
#' \eqn{(\lambda, n)}, but the result lands inside that same
#' \eqn{[0.899, 1.003]} band, so a printed loadings table buys essentially nothing over
#' \eqn{(k, n)} -- while the dominant source of variability \emph{across} studies is
#' which estimator was used (EFA / Schmid-Leiman / bifactor CFA / first principal
#' component), which no summary statistic identifies. The recommended interval remains a
#' bootstrap (BCa) one computed from raw data (Kelley & Pornprasertmanit, 2016).
#'
#' The standard error is therefore taken from what the primary study reported, in this
#' order:
#'
#' \enumerate{
#'   \item \code{omega_se}, interpreted on the \strong{natural} scale and mapped to the
#'     analysis scale by the delta method:
#'     \deqn{SE_T = SE_\omega / (1 - \omega)}  (bonett)
#'     \deqn{SE_T = SE_\omega / (3 (1 - \omega)^{2/3})}  (hakstian_whalen)
#'   \item otherwise the reported confidence interval, whose \strong{bounds} are transformed
#'     and differenced. Transforming the bounds rather than the point estimate is what makes
#'     this correct for the asymmetric bootstrap intervals omega usually carries, since a
#'     BCa interval is transformation-respecting.
#'   \item otherwise the effect size is returned with \code{omega_se = NA}. The row is then
#'     visible and countable in the pipeline (it appears in the extraction sheet and in the
#'     reporting-rate accounting) but is dropped by \code{summary()} and cannot be pooled --
#'     which is the honest outcome for a bare omega with no uncertainty attached.
#' }
#'
#' \strong{Lossiness.} Step 2 collapses a possibly asymmetric interval into a single
#' standard error, from which a symmetric Wald interval is rebuilt (the convention used by
#' every other \code{*_ci} route in this package). The point estimate and the interval
#' width are preserved; the asymmetry is not.
#'
#' \strong{omega_type is an estimand, not a label.} \code{"total"} is the proportion of
#' total score variance attributable to all common factors; \code{"hierarchical"} is the
#' proportion attributable to the \emph{general} factor only, and is systematically smaller.
#' They answer different questions and must not share a pool. \code{convert_df()} flags a
#' pool that mixes them.
#'
#' \strong{Do not pool omega with alpha either.} Under a congeneric model
#' \eqn{\omega_{total} \ge \alpha}, with equality only under exact tau-equivalence, so the
#' two are different quantities and their average is neither. There is a second, subtler
#' reason specific to this package: alpha's standard errors here are analytic functions of
#' \eqn{(n, k)} while omega's come from reported (often bootstrap) intervals, so the
#' \emph{weighting scheme} would differ systematically by metric, making metric perfectly
#' confounded with precision. Run two pools and report them side by side; if a study
#' reports both, estimate the gap as a within-study paired contrast rather than by
#' differencing two separately pooled means (the two coefficients are computed on the same
#' data and are near-perfectly dependent).
#'
#' \strong{Which estimator produced the number matters more than most moderators.}
#' Omega is not one computation. Zinbarg et al. (2006) applied several estimators to
#' identical data and found the first principal component overestimated
#' \eqn{\omega_h} by about \strong{+0.40} on average, against about +0.02 for a
#' hierarchical (bifactor) confirmatory factor analysis; Revelle & Zinbarg (2009)
#' record EFA returning \eqn{\omega_h = .04} where CFA returns exactly 0.0 on the
#' same data. A bias of that size is larger than any moderator effect a
#' reliability-generalization review is likely to report, so a pool that mixes
#' estimators can manufacture a "finding" that is purely an artefact of which
#' software each author ran. Record it in \code{omega_estimator}; a pool mixing two
#' \emph{known} estimators is flagged (rows left \code{"unspecified"} are not
#' counted as mixing, or the flag would fire on nearly every real dataset).
#'
#' In practice the estimator is usually codeable from the software a paper reports:
#' \code{psych::omega()} performs an EFA with a Schmid-Leiman transformation
#' (\code{"efa_schmid_leiman"}), \code{semTools::compRelSEM()} and
#' \code{MBESS::ci.reliability} are CFA-based (\code{"cfa_1factor"} or
#' \code{"cfa_bifactor"} according to the fitted model).
#'
#' \strong{Two extraction traps worth knowing before you copy a number.}
#'
#' \enumerate{
#'   \item \strong{Do not recompute omega from a printed loadings table unless the
#'     residual variances are printed too.} \code{lavaan}'s \code{std.lv = TRUE}
#'     fixes the \emph{factor} variance to 1, not the item variances, so the
#'     "Estimate" column is \strong{unstandardised}. Assuming
#'     \eqn{\theta_j = 1 - \lambda_j^2} from such a table gives 0.7725 for the
#'     worked example in Flora (2020) whose published omega is 0.5999 -- an error of
#'     +0.17 that nothing downstream can detect.
#'   \item \strong{On a bifactor model, the output row labelled "omega" is not
#'     \eqn{\omega_h}.} In \code{semTools} it is the total (0.97 in Flora's
#'     example) where \eqn{\omega_h} is 0.91. Copy the row you mean, and record
#'     which one in \code{omega_type}.
#' }
#'
#' \strong{Heywood cases.} An estimated \eqn{\omega \ge 1} arises from an improper solution
#' (a negative estimated error variance) and is not simply a transcription error. It has no
#' Bonett or Hakstian-Whalen transform, so those scales return \code{NA}; \code{"raw"}
#' returns the value so it stays inspectable.
#'
#' @export es_from_omega
#'
#' @references
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for population
#' reliability coefficients: Evaluation of methods, recommendations, and software for
#' composite measures. Psychological Methods, 21(1), 69-92.
#'
#' McDonald, R. P. (1999). Test theory: A unified treatment. Lawrence Erlbaum.
#'
#' Raykov, T. (2002). Analytic estimation of standard error and confidence interval for
#' scale reliability. Multivariate Behavioral Research, 37(1), 89-103.
#'
#' Revelle, W., & Zinbarg, R. E. (2009). Coefficients alpha, beta, omega, and the glb:
#' comments on Sijtsma. Psychometrika, 74(1), 145-154.
#'
#' Zinbarg, R. E., Yovel, I., Revelle, W., & McDonald, R. P. (2006). Estimating
#' generalizability to a latent variable common to all of a scale's indicators: A
#' comparison of estimators for omega_h. Applied Psychological Measurement, 30(2), 121-144.
#'
#' @md
#'
#' @return
#' This function estimates the standard error of McDonald's omega.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab omega\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab N/A\cr
#'  \tab \cr
#'  \code{required input data} \tab omega + (omega_se OR omega_ci_lo & omega_ci_up)\cr
#'  \tab \cr
#' }
#'
#' @examples
#' # from a reported standard error (natural scale)
#' es_from_omega(omega = 0.86, omega_se = 0.021, n_sample = 300, n_items = 12)
#'
#' # from a reported (possibly asymmetric) confidence interval
#' es_from_omega(omega = 0.86, omega_ci_lo = 0.81, omega_ci_up = 0.90,
#'               n_sample = 300, n_items = 12)
#'
#' # a bare omega: kept and visible, but with no SE it cannot be pooled
#' es_from_omega(omega = 0.86)
es_from_omega <- function(omega, omega_se, omega_ci_lo, omega_ci_up,
                          n_sample, n_items, omega_type = "total",
                          omega_estimator = "unspecified",
                          omega_to_es = "bonett") {

  len <- length(omega)
  if (missing(omega_se)) omega_se <- rep(NA_real_, len)
  if (missing(omega_ci_lo)) omega_ci_lo <- rep(NA_real_, len)
  if (missing(omega_ci_up)) omega_ci_up <- rep(NA_real_, len)
  if (missing(n_sample)) n_sample <- rep(NA_real_, len)
  if (missing(n_items)) n_items <- rep(NA_real_, len)
  if (missing(omega_type)) omega_type <- rep("total", len)
  if (missing(omega_estimator)) omega_estimator <- rep("unspecified", len)

  # Same recycling contract as es_from_cronbach_alpha(): a length-1 argument beside a
  # vector of coefficients is the natural way to write an extraction sheet, and without
  # this the row-index subsetting below would silently NA every row after the first.
  rec <- function(v, nm) {
    if (length(v) == 1) v <- rep(v, len)
    if (length(v) != len)
      stop(paste0("The length of the '", nm, "' argument is incorrectly specified."))
    v
  }
  omega_se <- rec(omega_se, "omega_se")
  omega_ci_lo <- rec(omega_ci_lo, "omega_ci_lo")
  omega_ci_up <- rec(omega_ci_up, "omega_ci_up")
  n_sample <- rec(n_sample, "n_sample")
  n_items <- rec(n_items, "n_items")
  omega_type <- rec(omega_type, "omega_type")
  omega_type[is.na(omega_type)] <- "total"
  omega_estimator <- rec(omega_estimator, "omega_estimator")

  # omega_type arrives from a DATA COLUMN, so it is user transcription rather than
  # an API argument, and must degrade per row: the package's own rule is that one
  # bad cell cannot abort a convert_df() run (cf. es_from_prop_single_group). A
  # capitalised "Total" -- the natural way a human extractor writes it -- used to
  # throw and take the whole run with it. Case is folded and the common synonyms
  # mapped; anything still unrecognised falls back to the default with a warning.
  omega_type <- .normalise_omega_type(omega_type)
  omega_estimator <- .normalise_omega_estimator(omega_estimator)
  if (!omega_to_es %in% c("bonett", "raw", "hakstian_whalen")) {
    stop(paste0("'", omega_to_es, "' not in tolerated values for the 'omega_to_es' argument. ",
                "Possible inputs are: 'bonett', 'raw', 'hakstian_whalen'"))
  }

  # A non-positive dispersion cannot be a standard error; a transposed interval is the
  # same interval (R/internal_guards.R).
  omega_se <- .positive_or_na(omega_se)
  ci_lo <- .ci_lower(omega_ci_lo, omega_ci_up)
  ci_up <- .ci_upper(omega_ci_lo, omega_ci_up)

  # forward transform and its derivative-based SE map
  fwd <- function(w) {
    if (identical(omega_to_es, "bonett")) return(log(1 - w))
    if (identical(omega_to_es, "hakstian_whalen")) return(1 - (1 - w)^(1 / 3))
    w
  }
  se_map <- function(w, se) {
    if (identical(omega_to_es, "bonett")) return(se / (1 - w))
    if (identical(omega_to_es, "hakstian_whalen")) return(se / (3 * (1 - w)^(2 / 3)))
    se
  }

  # omega >= 1 (a Heywood case) has no bonett/hakstian_whalen transform; keep it on the
  # raw scale so it stays inspectable rather than vanishing.
  transformable <- if (identical(omega_to_es, "raw")) {
    !is.na(omega)
  } else {
    !is.na(omega) & omega < 1
  }

  es <- rep(NA_real_, len)
  se <- rep(NA_real_, len)
  es[transformable] <- fwd(omega[transformable])

  # 1. reported SE (natural scale) -> delta method
  from_se <- which(transformable & !is.na(omega_se))
  if (length(from_se)) se[from_se] <- se_map(omega[from_se], omega_se[from_se])

  # 2. otherwise the reported CI, transformed at the BOUNDS so an asymmetric bootstrap
  #    interval maps correctly (BCa is transformation-respecting)
  # A bound of exactly 1.00 has no log(1 - .) or (1 - .)^(1/3) transform, so it must
  # be excluded on those scales -- but on the raw scale fwd() is the identity and a
  # reported "omega = .88 [.79, 1.00]" (a routine bootstrap print-out at high
  # reliability) is perfectly usable. The guard is therefore scale-conditional;
  # before, such a row silently lost its SE and became unpoolable for no reason.
  ci_usable <- if (identical(omega_to_es, "raw")) rep(TRUE, len) else (ci_lo < 1 & ci_up < 1)
  from_ci <- which(transformable & is.na(se) &
                     !is.na(ci_lo) & !is.na(ci_up) & ci_usable)
  if (length(from_ci)) {
    t_lo <- fwd(ci_lo[from_ci])
    t_up <- fwd(ci_up[from_ci])
    se[from_ci] <- abs(t_up - t_lo) / (2 * qnorm(0.975))
  }
  se <- .positive_or_na(se)

  omega_ci_lo_out <- es - qnorm(0.975) * se
  omega_ci_up_out <- es + qnorm(0.975) * se

  data.frame(
    omega = es,
    omega_se = se,
    omega_ci_lo = omega_ci_lo_out,
    omega_ci_up = omega_ci_up_out,
    n_sample = n_sample,
    n_items = n_items,
    omega_type = omega_type,
    omega_estimator = omega_estimator,
    info_used = "omega"
  )
}

# Fold case and map the common synonyms for omega_type. Unrecognised values fall
# back to the default with a warning rather than aborting the run (see the note in
# es_from_omega()).
.normalise_omega_type <- function(x, warn = TRUE) {
  raw <- as.character(x)
  key <- tolower(trimws(raw))
  key <- gsub("[^a-z]", "", key)
  map <- c(
    total = "total", omegat = "total", omegatotal = "total", tot = "total",
    omegatot = "total",
    hierarchical = "hierarchical", omegah = "hierarchical", hierarchique = "hierarchical",
    hier = "hierarchical", general = "hierarchical", omegahierarchical = "hierarchical",
    generalfactor = "hierarchical", omegageneral = "hierarchical",
    asymptotic = "asymptotic", omegalim = "asymptotic", asymp = "asymptotic",
    subscale = "subscale", subscales = "subscale", sub = "subscale",
    # the normaliser must be IDEMPOTENT: convert_df() normalises the stored column
    # and es_from_omega() normalises again, so every output must be a valid input
    unspecified = "unspecified", unknown = "unspecified", notreported = "unspecified",
    nr = "unspecified", none = "unspecified"
  )
  out <- unname(map[key])
  bad <- which(is.na(out) & !is.na(raw) & nzchar(key))
  if (length(bad) && isTRUE(warn)) {
    warning(paste0(
      "Unrecognised omega_type value(s): '", paste(unique(raw[bad]), collapse = "', '"),
      "'. Treated as 'unspecified'. Recognised values are 'total', 'hierarchical', ",
      "'asymptotic', 'subscale' (case-insensitive; 'omega_t'/'omega_h'/'omega.lim' ",
      "are also accepted)."))
  }
  # Fall back to "unspecified", NOT to "total". Mapping an unrecognised value onto a
  # real estimand is worse than leaving it unknown: it silently relabels a
  # hierarchical omega as a total one and DISARMS V38, the flag that exists for
  # exactly that error. (.normalise_omega_estimator() already fell back this way;
  # this one did not.) V38 ignores "unspecified" the same way V39 does.
  out[is.na(out)] <- "unspecified"
  out
}

# Fold case and map the common synonyms for omega_estimator. Same per-row contract as
# .normalise_omega_type(): this arrives from a data column, so an unrecognised value
# warns and falls back rather than aborting the run.
#
# The fallback is "unspecified" rather than any real estimator, deliberately: guessing
# a specific one would let the mixing flag fire on a guess.
.normalise_omega_estimator <- function(x, warn = TRUE) {
  raw <- as.character(x)
  key <- gsub("[^a-z]", "", tolower(trimws(raw)))
  map <- c(
    cfabifactor = "cfa_bifactor", bifactor = "cfa_bifactor", bifactorcfa = "cfa_bifactor",
    hicf = "cfa_bifactor", hierarchicalcfa = "cfa_bifactor", omegahcfa = "cfa_bifactor",
    cfafactor = "cfa_1factor", cfaonefactor = "cfa_1factor", cfa = "cfa_1factor",
    onefactor = "cfa_1factor", unidimensional = "cfa_1factor", congeneric = "cfa_1factor",
    factorcfa = "cfa_1factor", sem = "cfa_1factor",
    efaschmidleiman = "efa_schmid_leiman", schmidleiman = "efa_schmid_leiman",
    sl = "efa_schmid_leiman", efa = "efa_schmid_leiman", efasl = "efa_schmid_leiman",
    hocf = "efa_schmid_leiman", higherorder = "efa_schmid_leiman", psychomega = "efa_schmid_leiman",
    # psych::omega() defaults to an EFA with a Schmid-Leiman transformation, so the
    # package name alone determines the estimator. semTools / MBESS / lavaan do NOT:
    # they are CFA-based, but whether the fitted model is unidimensional or bifactor
    # is the analyst's, so those names are deliberately left to fall through to
    # "unspecified" rather than guessing a level the mixing flag would then act on.
    psych = "efa_schmid_leiman",
    firstpc = "first_pc", pc = "first_pc", principalcomponent = "first_pc",
    firstprincipalcomponent = "first_pc", pca = "first_pc",
    firstpf = "first_pf", pf = "first_pf", principalfactor = "first_pf",
    firstprincipalfactor = "first_pf", paf = "first_pf",
    unspecified = "unspecified", unknown = "unspecified", notreported = "unspecified",
    nr = "unspecified", na = "unspecified", none = "unspecified"
  )
  out <- unname(map[key])
  bad <- which(is.na(out) & !is.na(raw) & nzchar(key))
  if (length(bad) && isTRUE(warn)) {
    warning(paste0(
      "Unrecognised omega_estimator value(s): '", paste(unique(raw[bad]), collapse = "', '"),
      "'. Treated as 'unspecified'. Recognised values are 'cfa_bifactor', ",
      "'cfa_1factor', 'efa_schmid_leiman', 'first_pc', 'first_pf', 'unspecified' ",
      "(case-insensitive; 'bifactor', 'schmid_leiman', 'psych_omega', 'pca' and ",
      "similar synonyms are also accepted)."))
  }
  out[is.na(out)] <- "unspecified"
  out
}
