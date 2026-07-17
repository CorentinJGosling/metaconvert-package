# known es types
.known_es_types <- function() {
  c("d", "g", "md", "or", "logor", "rr", "logrr", "r", "z", "rd",
    "irr", "logirr", "nnt", "dw", "gw", "mdw", "hr", "loghr")
}

# output slots of the entered measure (dw mirrors d, etc.)
.user_slot_prefixes <- function(orig_type) {
  switch(orig_type,
    "d"   = c("d", "dw"),
    "dw"  = c("dw", "d"),
    "g"   = c("g", "gw"),
    "gw"  = c("gw", "g"),
    "md"  = c("md", "mdw"),
    "mdw" = c("mdw", "md"),
    "or"  = , "logor"  = "logor",
    "rr"  = , "logrr"  = "logrr",
    "irr" = , "logirr" = "logirr",
    "hr"  = , "loghr"  = "loghr",
    "r"   = "r",
    "z"   = "z",
    "rd"  = "rd",
    "nnt" = "nnt",
    NULL
  )
}

# CI to restore on the entered measure after dispatch (converters rebuild
# t-based CIs; the entered slot must stay wald-z coherent with the recovered SE)
.user_preserved_ci <- function(orig_measure_input, user_ci_lo, user_ci_up,
                                es_val, se_val) {
  orig <- tolower(as.character(orig_measure_input))
  user_gave_ci <- !is.na(user_ci_lo) & !is.na(user_ci_up)
  natural_ratio <- orig %in% c("or", "rr", "irr", "hr")
  # ratio originals without a user ci: converters already build log-scale CIs
  additive <- orig %in% c("d", "g", "md", "dw", "gw", "mdw", "r", "z", "rd")

  log_lo <- suppressWarnings(log(user_ci_lo))
  log_up <- suppressWarnings(log(user_ci_up))
  log_lo[!is.finite(log_lo)] <- NA_real_
  log_up[!is.finite(log_up)] <- NA_real_

  z_lo <- ifelse(additive & !is.na(es_val) & !is.na(se_val),
                 es_val - qnorm(.975) * se_val, NA_real_)
  z_up <- ifelse(additive & !is.na(es_val) & !is.na(se_val),
                 es_val + qnorm(.975) * se_val, NA_real_)

  lo <- ifelse(user_gave_ci, ifelse(natural_ratio, log_lo, user_ci_lo), z_lo)
  up <- ifelse(user_gave_ci, ifelse(natural_ratio, log_up, user_ci_up), z_up)
  list(lo = lo, up = up, user_gave_ci = user_gave_ci)
}

.dispatch_user_conversion <- function(original_type, es_val, se_val,
                                       n_exp, n_nexp, n_sample,
                                       n_cases, n_controls,
                                       baseline_risk, small_margin_prop,
                                       or_to_rr, or_to_cor, smd_to_cor,
                                       cor_to_smd, rr_to_or) {
  len <- length(es_val)
  no_reverse <- rep(FALSE, len)

  switch(original_type,
    "d" = {
      res_d <- .es_from_d(d = es_val, d_se = se_val,
                 n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample,
                 smd_to_cor = rep(smd_to_cor, len), reverse = no_reverse)
      # dw = d for user input
      res_d$dw       <- res_d$d
      res_d$dw_se    <- res_d$d_se
      res_d$dw_ci_lo <- res_d$d_ci_lo
      res_d$dw_ci_up <- res_d$d_ci_up
      res_d
    },
    "g" = {
      # g -> d, keeping the user se
      df <- n_exp + n_nexp - 2
      J <- .d_j(df)
      # no sample sizes: J = 1
      J <- ifelse(is.na(J), 1, J)
      d <- es_val / J
      d_se <- se_val / J
      res_g <- .es_from_d(d = d, d_se = d_se,
                 n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample,
                 smd_to_cor = rep(smd_to_cor, len), reverse = no_reverse)
      # keep user g/se when conversion returned NA
      res_g$g <- ifelse(is.na(res_g$g) & !is.na(es_val), es_val, res_g$g)
      res_g$g_se <- ifelse(is.na(res_g$g_se) & !is.na(se_val), se_val, res_g$g_se)
      res_g$g_ci_lo <- ifelse(is.na(res_g$g_ci_lo) & !is.na(es_val) & !is.na(se_val),
                               es_val - qnorm(.975) * se_val, res_g$g_ci_lo)
      res_g$g_ci_up <- ifelse(is.na(res_g$g_ci_up) & !is.na(es_val) & !is.na(se_val),
                               es_val + qnorm(.975) * se_val, res_g$g_ci_up)
      res_g$gw       <- res_g$g
      res_g$gw_se    <- res_g$g_se
      res_g$gw_ci_lo <- res_g$g_ci_lo
      res_g$gw_ci_up <- res_g$g_ci_up
      res_g
    },
    "logor" = {
      es_from_or_se(or = exp(es_val), logor_se = se_val,
                    baseline_risk = baseline_risk,
                    small_margin_prop = small_margin_prop,
                    n_exp = n_exp, n_nexp = n_nexp,
                    n_cases = n_cases, n_controls = n_controls,
                    n_sample = n_sample,
                    or_to_rr = or_to_rr, or_to_cor = or_to_cor,
                    reverse_or = no_reverse)
    },
    "or" = {
      # se on the log scale
      es_from_or_se(or = es_val, logor_se = se_val,
                    baseline_risk = baseline_risk,
                    small_margin_prop = small_margin_prop,
                    n_exp = n_exp, n_nexp = n_nexp,
                    n_cases = n_cases, n_controls = n_controls,
                    n_sample = n_sample,
                    or_to_rr = or_to_rr, or_to_cor = or_to_cor,
                    reverse_or = no_reverse)
    },
    "logrr" = {
      # A risk ratio is a ratio-family measure: convert to OR / NNT / RD only,
      # NOT to an SMD or correlation. Reaching a standardized family would route
      # through the OR and an assumed baseline risk (not identified from the RR
      # alone, anti-conservative SE); this matches the risk-ratio pipeline
      # measure (es_from_rr_se) and the risk-difference user path. To obtain a
      # D/G/R/Z from a risk ratio, enter it as an odds ratio instead.
      es_from_rr_se(rr = exp(es_val), logrr_se = se_val,
                    baseline_risk = baseline_risk,
                    n_exp = n_exp, n_nexp = n_nexp,
                    n_cases = n_cases, n_controls = n_controls,
                    rr_to_or = rr_to_or,
                    reverse_rr = no_reverse)
    },
    "rr" = {
      # See the logrr case: a risk ratio converts to OR / NNT / RD only.
      es_from_rr_se(rr = es_val, logrr_se = se_val,
                    baseline_risk = baseline_risk,
                    n_exp = n_exp, n_nexp = n_nexp,
                    n_cases = n_cases, n_controls = n_controls,
                    rr_to_or = rr_to_or,
                    reverse_rr = no_reverse)
    },
    "r" = {
      r <- es_val
      r_se_user <- se_val
      n_s <- ifelse(!is.na(n_sample), n_sample, n_exp + n_nexp)
      n_e <- ifelse(!is.na(n_exp), n_exp, n_s / 2)
      n_ne <- ifelse(!is.na(n_nexp), n_nexp, n_s / 2)

      nn <- which(!is.na(r) & !is.na(r_se_user) & !is.na(n_s))
      d_val <- d_se_val <- rep(NA_real_, len)
      if (length(nn) > 0) {
        d_res <- t(mapply(.cor_to_smd,
          r = r[nn], r_se = r_se_user[nn], n_sample = n_s[nn],
          sd_iv = NA, unit_increase_iv = NA, unit_type = NA,
          cor_to_smd = cor_to_smd))
        d_val[nn] <- unlist(d_res[, 1])
        d_se_val[nn] <- unlist(d_res[, 2])
      }

      es <- .es_from_d(d = d_val, d_se = d_se_val,
                       n_exp = n_e, n_nexp = n_ne, n_sample = n_s,
                       reverse = no_reverse)

      es$r <- r
      es$r_se <- r_se_user
      es$r_ci_lo <- r - qt(.975, n_s - 2) * r_se_user
      es$r_ci_up <- r + qt(.975, n_s - 2) * r_se_user
      z_val <- atanh(r)
      z_se_val <- r_se_user / (1 - r^2)
      es$z <- z_val
      es$z_se <- z_se_val
      es$z_ci_lo <- z_val - qnorm(.975) * z_se_val
      es$z_ci_up <- z_val + qnorm(.975) * z_se_val
      es
    },
    "z" = {
      z_user <- es_val
      z_se_user <- se_val
      r <- tanh(z_user)
      r_se_user <- z_se_user * (1 - r^2)  # inverse delta method
      n_s <- ifelse(!is.na(n_sample), n_sample, n_exp + n_nexp)
      n_e <- ifelse(!is.na(n_exp), n_exp, n_s / 2)
      n_ne <- ifelse(!is.na(n_nexp), n_nexp, n_s / 2)

      nn <- which(!is.na(r) & !is.na(r_se_user) & !is.na(n_s))
      d_val <- d_se_val <- rep(NA_real_, len)
      if (length(nn) > 0) {
        d_res <- t(mapply(.cor_to_smd,
          r = r[nn], r_se = r_se_user[nn], n_sample = n_s[nn],
          sd_iv = NA, unit_increase_iv = NA, unit_type = NA,
          cor_to_smd = cor_to_smd))
        d_val[nn] <- unlist(d_res[, 1])
        d_se_val[nn] <- unlist(d_res[, 2])
      }

      es <- .es_from_d(d = d_val, d_se = d_se_val,
                       n_exp = n_e, n_nexp = n_ne, n_sample = n_s,
                       reverse = no_reverse)

      es$r <- r
      es$r_se <- r_se_user
      es$r_ci_lo <- r - qt(.975, n_s - 2) * r_se_user
      es$r_ci_up <- r + qt(.975, n_s - 2) * r_se_user
      es$z <- z_user
      es$z_se <- z_se_user
      es$z_ci_lo <- z_user - qnorm(.975) * z_se_user
      es$z_ci_up <- z_user + qnorm(.975) * z_se_user
      es
    },
    "rd" = {
      es_from_rd_se(rd = es_val, rd_se = se_val,
                    baseline_risk = baseline_risk,
                    n_exp = n_exp, n_nexp = n_nexp,
                    n_cases = n_cases, n_controls = n_controls,
                    n_sample = n_sample,
                    reverse_rd = no_reverse)
    },
    "md" = {
      can_convert <- !is.na(n_exp) & !is.na(n_nexp) & n_exp > 0 & n_nexp > 0
      if (any(can_convert)) {
        res_md <- es_from_md_se(md = es_val, md_se = se_val,
                      n_exp = n_exp, n_nexp = n_nexp,
                      smd_to_cor = smd_to_cor, reverse_md = no_reverse)
        no_convert <- which(!can_convert)
        if (length(no_convert) > 0) {
          res_md$md[no_convert] <- es_val[no_convert]
          res_md$md_se[no_convert] <- se_val[no_convert]
          res_md$md_ci_lo[no_convert] <- es_val[no_convert] - qnorm(.975) * se_val[no_convert]
          res_md$md_ci_up[no_convert] <- es_val[no_convert] + qnorm(.975) * se_val[no_convert]
        }
      } else {
        res_md <- data.frame(
          md = es_val,
          md_se = se_val,
          md_ci_lo = es_val - qnorm(.975) * se_val,
          md_ci_up = es_val + qnorm(.975) * se_val
        )
      }
      res_md$mdw       <- res_md$md
      res_md$mdw_se    <- res_md$md_se
      res_md$mdw_ci_lo <- res_md$md_ci_lo
      res_md$mdw_ci_up <- res_md$md_ci_up
      res_md
    },
    "logirr" = {
      data.frame(
        logirr = es_val,
        logirr_se = se_val,
        logirr_ci_lo = es_val - qnorm(.975) * se_val,
        logirr_ci_up = es_val + qnorm(.975) * se_val
      )
    },
    "loghr" = {
      data.frame(
        loghr = es_val,
        loghr_se = se_val,
        loghr_ci_lo = es_val - qnorm(.975) * se_val,
        loghr_ci_up = es_val + qnorm(.975) * se_val
      )
    },
    "nnt" = {
      # nnt -> rd
      rd <- ifelse(es_val == 0, NA, 1 / es_val)
      rd_se <- ifelse(es_val == 0 | se_val == 0, NA, se_val / es_val^2)
      res <- es_from_rd_se(rd = rd, rd_se = rd_se,
                           baseline_risk = baseline_risk,
                           n_exp = n_exp, n_nexp = n_nexp,
                           n_cases = n_cases, n_controls = n_controls,
                           n_sample = n_sample,
                           reverse_rd = no_reverse)
      # keep the RD-inverted CI (Altman 1998)
      res$nnt <- es_val
      res$nnt_se <- se_val
      res
    },
    "mdw" = {
      data.frame(
        mdw       = es_val,
        mdw_se    = se_val,
        mdw_ci_lo = es_val - qnorm(.975) * se_val,
        mdw_ci_up = es_val + qnorm(.975) * se_val,
        md        = es_val,
        md_se     = se_val,
        md_ci_lo  = es_val - qnorm(.975) * se_val,
        md_ci_up  = es_val + qnorm(.975) * se_val
      )
    },
    "dw" = {
      data.frame(
        dw       = es_val,
        dw_se    = se_val,
        dw_ci_lo = es_val - qnorm(.975) * se_val,
        dw_ci_up = es_val + qnorm(.975) * se_val,
        d        = es_val,
        d_se     = se_val,
        d_ci_lo  = es_val - qnorm(.975) * se_val,
        d_ci_up  = es_val + qnorm(.975) * se_val
      )
    },
    "gw" = {
      data.frame(
        gw       = es_val,
        gw_se    = se_val,
        gw_ci_lo = es_val - qnorm(.975) * se_val,
        gw_ci_up = es_val + qnorm(.975) * se_val,
        g        = es_val,
        g_se     = se_val,
        g_ci_lo  = es_val - qnorm(.975) * se_val,
        g_ci_up  = es_val + qnorm(.975) * se_val
      )
    },
    NULL
  )
}


#' Directly input a crude effect size value + variance, with optional conversion
#'
#' @param user_es_original_measure_crude the type of effect size entered (see details).
#' @param user_es_crude effect size value
#' @param user_se_crude standard error of the effect size (on the log scale for the or/rr/irr/hr measures)
#' @param user_ci_lo_crude lower bound of the 95% CI around the effect size value
#' @param user_ci_up_crude upper bound of the 95% CI around the effect size value
#' @param user_es_target_measure_crude the effect size measure of the output (automatically set to the \code{measure} argument when called by the \code{\link{convert_df}} function)
#' @param n_exp number of participants in the exposed group
#' @param n_nexp number of participants in the non-exposed group
#' @param n_sample total number of participants in the sample
#' @param n_cases number of cases/events across exposed/non-exposed groups
#' @param n_controls number of controls/no-event across exposed/non-exposed groups
#' @param baseline_risk proportion of cases in the non-exposed group
#' @param small_margin_prop smallest margin proportion of cases/events in the underlying 2x2 table
#' @param or_to_rr formula used to convert an odds ratio value into a risk ratio (see \code{\link{es_from_or_se}}).
#' @param or_to_cor formula used to convert an odds ratio value into a correlation coefficient (see \code{\link{es_from_or_se}}).
#' @param smd_to_cor formula used to convert a SMD value into a coefficient correlation (see \code{\link{es_from_cohen_d}}).
#' @param cor_to_smd formula used to convert a correlation coefficient value into a SMD (see \code{\link{es_from_pearson_r}}).
#' @param rr_to_or formula used to convert a risk ratio value into an odds ratio (see \code{\link{es_from_rr_se}}).
#' @param measure deprecated alias for \code{user_es_target_measure_crude}, kept for backward compatibility with metaConvert <= 1.0.3.
#' @param user_es_measure_crude deprecated alias for \code{user_es_original_measure_crude}, kept for backward compatibility with metaConvert <= 1.0.3.
#'
#' @details
#' This function is a generic function allowing to include any crude effect size measure value + variance.
#' Importantly, when the \code{user_es_original_measure_crude} is one of the known measures
#' (d, g, md, mdw, dw, gw, or, logor, rr, logrr, irr, logirr, hr, r, z, rd, nnt),
#' conversions towards the other effect size measures are performed.
#' Otherwise, no conversion is performed (the effect size value + variance you enter is
#' the value + variance exported by this function) and a warning is issued.
#' The sample sizes and baseline risk are used only to perform the conversions.
#'
#' For the or/rr/irr/hr measures, the standard error you enter must be on the log scale.
#' If you indicate the 95% CI bounds instead of the standard error, the log-transformation
#' is applied automatically.
#'
#' @md
#'
#' @export es_from_user_crude
#'
#' @return
#' This function allows to directly input any of the available effect size measures
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab Any of the available measures\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab Depends on the entered measure (see details)\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 23. User's input (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @examples
#' dat = data.frame(user_es_original_measure_crude = "or", user_es_crude = 2.5,
#'                  user_ci_lo_crude = 1.2, user_ci_up_crude = 5.2,
#'                  n_exp = 120, n_nexp = 44)
#' summary(convert_df(dat, measure = "g"))
es_from_user_crude <- function(user_es_original_measure_crude,
                                user_es_crude, user_se_crude,
                                user_ci_lo_crude, user_ci_up_crude,
                                user_es_target_measure_crude = "g",
                                n_exp, n_nexp, n_sample,
                                n_cases, n_controls, baseline_risk,
                                small_margin_prop,
                                or_to_rr = "metaumbrella_cases",
                                or_to_cor = "pearson",
                                smd_to_cor = "viechtbauer",
                                cor_to_smd = "viechtbauer",
                                rr_to_or = "metaumbrella",
                                measure, user_es_measure_crude) {

  # Backward compatibility with the pre-2.0 argument names (metaConvert <= 1.0.3,
  # still used by the published metaumbrella): 'user_es_measure_crude' was the
  # entered (original) measure and 'measure' was the requested (target) measure.
  if (!missing(user_es_measure_crude) && missing(user_es_original_measure_crude))
    user_es_original_measure_crude <- user_es_measure_crude
  # The alias must not silently override an explicitly supplied new-style
  # argument, and its use should be visible to the caller.
  if (!missing(measure) && missing(user_es_target_measure_crude)) {
    warning("The 'measure' argument of es_from_user_crude() is deprecated; use 'user_es_target_measure_crude'.",
            call. = FALSE)
    user_es_target_measure_crude <- measure
  }

  len <- length(user_es_original_measure_crude)

  if (missing(user_es_crude)) user_es_crude <- rep(NA_real_, len)
  if (missing(user_se_crude)) user_se_crude <- rep(NA_real_, len)
  if (missing(user_ci_lo_crude)) user_ci_lo_crude <- rep(NA_real_, len)
  if (missing(user_ci_up_crude)) user_ci_up_crude <- rep(NA_real_, len)
  if (missing(n_exp)) n_exp <- rep(NA_real_, len)
  if (missing(n_nexp)) n_nexp <- rep(NA_real_, len)
  if (missing(n_sample)) n_sample <- rep(NA_real_, len)
  if (missing(n_cases)) n_cases <- rep(NA_real_, len)
  if (missing(n_controls)) n_controls <- rep(NA_real_, len)
  if (missing(baseline_risk)) baseline_risk <- rep(NA_real_, len)
  if (missing(small_margin_prop)) small_margin_prop <- rep(NA_real_, len)

  # es/se from ci bounds
  es_val <- ifelse(is.na(user_es_crude) & !is.na(user_es_original_measure_crude) &
                    !is.na(user_ci_lo_crude) & !is.na(user_ci_up_crude),
    (user_ci_up_crude + user_ci_lo_crude) / 2,
    user_es_crude
  )

  se_val <- ifelse(is.na(user_se_crude) & !is.na(user_es_original_measure_crude) &
                    !is.na(user_ci_lo_crude) & !is.na(user_ci_up_crude),
    (user_ci_up_crude - user_ci_lo_crude) / (2 * qnorm(.975)),
    user_se_crude
  )

  ci_lo <- ifelse(is.na(user_ci_lo_crude) & !is.na(user_es_original_measure_crude) &
                   !is.na(user_es_crude) & !is.na(user_se_crude),
    user_es_crude - user_se_crude * qnorm(.975),
    user_ci_lo_crude
  )

  ci_up <- ifelse(is.na(user_ci_up_crude) & !is.na(user_es_original_measure_crude) &
                   !is.na(user_es_crude) & !is.na(user_se_crude),
    user_es_crude + user_se_crude * qnorm(.975),
    user_ci_up_crude
  )

  target <- unique(user_es_target_measure_crude)
  if (length(target) != 1) target <- target[1]

  orig_type <- tolower(user_es_original_measure_crude)
  known_types <- .known_es_types()

  # or/rr/irr/hr: se from ci on the log scale
  ratio_from_ci <- which(orig_type %in% c("or", "rr", "irr", "hr") &
                          is.na(user_se_crude) &
                          !is.na(user_ci_lo_crude) & !is.na(user_ci_up_crude) &
                          user_ci_lo_crude > 0 & user_ci_up_crude > 0)
  if (length(ratio_from_ci) > 0) {
    log_lo <- log(user_ci_lo_crude[ratio_from_ci])
    log_up <- log(user_ci_up_crude[ratio_from_ci])
    se_val[ratio_from_ci] <- (log_up - log_lo) / (2 * qnorm(.975))
    es_val[ratio_from_ci] <- ifelse(
      is.na(user_es_crude[ratio_from_ci]),
      (log_lo + log_up) / 2,
      log(user_es_crude[ratio_from_ci])
    )
    ci_lo[ratio_from_ci] <- log_lo
    ci_up[ratio_from_ci] <- log_up
    ot <- orig_type[ratio_from_ci]
    orig_type[ratio_from_ci] <- ifelse(ot == "or",  "logor",
                                ifelse(ot == "rr",  "logrr",
                                ifelse(ot == "irr", "logirr",
                                ifelse(ot == "hr",  "loghr", ot))))
  }

  # natural irr/hr -> log scale
  irr_with_se <- which(orig_type == "irr" & !is.na(se_val))
  if (length(irr_with_se) > 0) {
    es_val[irr_with_se] <- log(es_val[irr_with_se])
    orig_type[irr_with_se] <- "logirr"
  }

  hr_with_se <- which(orig_type == "hr" & !is.na(se_val))
  if (length(hr_with_se) > 0) {
    es_val[hr_with_se] <- log(es_val[hr_with_se])
    orig_type[hr_with_se] <- "loghr"
  }

  ratio_se_check <- which(orig_type %in% c("or", "rr") &
                           !is.na(user_se_crude) & !is.na(user_es_crude) &
                           user_es_crude > 0 & user_se_crude > user_es_crude)
  if (length(ratio_se_check) > 0) {
    warning(paste0(
      "For rows with user_es_original_measure_crude = 'or'/'rr', ",
      "user_se_crude must be on the LOG scale (SE of logOR/logRR), ",
      "not the natural scale. Detected SE > ES in ",
      length(ratio_se_check), " row(s), which suggests the SE may be ",
      "on the natural scale. If so, provide CI bounds instead and the ",
      "log-scale SE will be derived automatically."
    ))
  }

  has_known_type <- !is.na(orig_type) & orig_type %in% known_types
  has_unknown_type <- !is.na(orig_type) & !(orig_type %in% known_types)

  if (any(has_unknown_type)) {
    unknown_vals <- unique(orig_type[has_unknown_type])
    warning(paste0(
      "Unknown user_es_original_measure_crude values: ",
      paste(unknown_vals, collapse = ", "),
      ". No conversion applied for these rows. Known types: ",
      paste(known_types, collapse = ", "), "."
    ))
  }

  res <- data.frame(info_used = rep("user_input_crude", len))

  res[[target]] <- rep(NA_real_, len)
  res[[paste0(target, "_se")]] <- rep(NA_real_, len)
  res[[paste0(target, "_ci_lo")]] <- rep(NA_real_, len)
  res[[paste0(target, "_ci_up")]] <- rep(NA_real_, len)

  passthrough_idx <- which(!has_known_type & !is.na(es_val))
  if (length(passthrough_idx) > 0) {
    res[[target]][passthrough_idx] <- es_val[passthrough_idx]
    res[[paste0(target, "_se")]][passthrough_idx] <- se_val[passthrough_idx]
    res[[paste0(target, "_ci_lo")]][passthrough_idx] <- ci_lo[passthrough_idx]
    res[[paste0(target, "_ci_up")]][passthrough_idx] <- ci_up[passthrough_idx]
  }

  for (es_type in unique(orig_type[has_known_type])) {
    idx <- which(has_known_type & orig_type == es_type)
    if (length(idx) == 0) next

    converted <- .dispatch_user_conversion(
      original_type = es_type,
      es_val = es_val[idx], se_val = se_val[idx],
      n_exp = n_exp[idx], n_nexp = n_nexp[idx],
      n_sample = n_sample[idx],
      n_cases = n_cases[idx], n_controls = n_controls[idx],
      baseline_risk = baseline_risk[idx],
      small_margin_prop = small_margin_prop[idx],
      or_to_rr = or_to_rr, or_to_cor = or_to_cor,
      smd_to_cor = smd_to_cor, cor_to_smd = cor_to_smd,
      rr_to_or = rr_to_or
    )

    if (!is.null(converted)) {
      for (col in names(converted)) {
        if (col == "info_used") next
        if (!col %in% names(res)) res[[col]] <- rep(NA_real_, len)
        res[[col]][idx] <- converted[[col]]
      }
    }
  }

  # restore user ci on the entered measure
  pres <- .user_preserved_ci(user_es_original_measure_crude,
                              user_ci_lo_crude, user_ci_up_crude,
                              es_val, se_val)
  for (es_type in unique(orig_type[has_known_type])) {
    prefixes <- .user_slot_prefixes(es_type)
    if (is.null(prefixes)) next
    idx <- which(has_known_type & orig_type == es_type &
                 !is.na(pres$lo) & !is.na(pres$up))
    if (es_type == "nnt") idx <- idx[pres$user_gave_ci[idx]]
    if (length(idx) == 0) next
    for (prefix in prefixes) {
      lo_col <- paste0(prefix, "_ci_lo")
      up_col <- paste0(prefix, "_ci_up")
      if (!lo_col %in% names(res) || !up_col %in% names(res)) next
      res[[lo_col]][idx] <- pres$lo[idx]
      res[[up_col]][idx] <- pres$up[idx]
    }
  }

  attr(res, "measure") <- user_es_original_measure_crude
  return(res)
}


#' Directly input an adjusted effect size value + variance, with optional conversion
#'
#' @param user_es_original_measure_adj the type of effect size entered (see details).
#' @param user_es_adj adjusted effect size value
#' @param user_se_adj adjusted standard error of the effect size (on the log scale for the or/rr/irr/hr measures)
#' @param user_ci_lo_adj adjusted lower bound of the 95% CI around the effect size value
#' @param user_ci_up_adj adjusted upper bound of the 95% CI around the effect size value
#' @param user_es_target_measure_adj the effect size measure of the output (automatically set to the \code{measure} argument when called by the \code{\link{convert_df}} function)
#' @param n_exp number of participants in the exposed group
#' @param n_nexp number of participants in the non-exposed group
#' @param n_sample total number of participants in the sample
#' @param n_cases number of cases/events across exposed/non-exposed groups
#' @param n_controls number of controls/no-event across exposed/non-exposed groups
#' @param baseline_risk proportion of cases in the non-exposed group
#' @param small_margin_prop smallest margin proportion of cases/events in the underlying 2x2 table
#' @param or_to_rr formula used to convert an odds ratio value into a risk ratio (see \code{\link{es_from_or_se}}).
#' @param or_to_cor formula used to convert an odds ratio value into a correlation coefficient (see \code{\link{es_from_or_se}}).
#' @param smd_to_cor formula used to convert a SMD value into a coefficient correlation (see \code{\link{es_from_cohen_d}}).
#' @param cor_to_smd formula used to convert a correlation coefficient value into a SMD (see \code{\link{es_from_pearson_r}}).
#' @param rr_to_or formula used to convert a risk ratio value into an odds ratio (see \code{\link{es_from_rr_se}}).
#' @param measure deprecated alias for \code{user_es_target_measure_adj}, kept for backward compatibility with metaConvert <= 1.0.3.
#' @param user_es_measure_adj deprecated alias for \code{user_es_original_measure_adj}, kept for backward compatibility with metaConvert <= 1.0.3.
#'
#' @details
#' This function is a generic function allowing to include any adjusted effect size measure value + variance.
#' Importantly, when the \code{user_es_original_measure_adj} is one of the known measures
#' (d, g, md, mdw, dw, gw, or, logor, rr, logrr, irr, logirr, hr, r, z, rd, nnt),
#' conversions towards the other effect size measures are performed.
#' Otherwise, no conversion is performed (the effect size value + variance you enter is
#' the value + variance exported by this function) and a warning is issued.
#' The sample sizes and baseline risk are used only to perform the conversions.
#'
#' For the or/rr/irr/hr measures, the standard error you enter must be on the log scale.
#' If you indicate the 95% CI bounds instead of the standard error, the log-transformation
#' is applied automatically.
#'
#' @md
#'
#' @export es_from_user_adj
#'
#' @return
#' This function allows to directly input any of the available effect size measures
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab Any of the available measures\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab Depends on the entered measure (see details)\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 24. User's input (adjusted)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @examples
#' dat = data.frame(user_es_original_measure_adj = "or", user_es_adj = 2.5,
#'                  user_ci_lo_adj = 1.2, user_ci_up_adj = 5.2,
#'                  n_exp = 120, n_nexp = 44)
#' summary(convert_df(dat, measure = "g"))
es_from_user_adj <- function(user_es_original_measure_adj,
                              user_es_adj, user_se_adj,
                              user_ci_lo_adj, user_ci_up_adj,
                              user_es_target_measure_adj = "g",
                              n_exp, n_nexp, n_sample,
                              n_cases, n_controls, baseline_risk,
                              small_margin_prop,
                              or_to_rr = "metaumbrella_cases",
                              or_to_cor = "pearson",
                              smd_to_cor = "viechtbauer",
                              cor_to_smd = "viechtbauer",
                              rr_to_or = "metaumbrella",
                              measure, user_es_measure_adj) {

  # Backward compatibility with the pre-2.0 argument names (metaConvert <= 1.0.3):
  # 'user_es_measure_adj' was the entered (original) measure and 'measure' was the
  # requested (target) measure.
  if (!missing(user_es_measure_adj) && missing(user_es_original_measure_adj))
    user_es_original_measure_adj <- user_es_measure_adj
  # The alias must not silently override an explicitly supplied new-style
  # argument, and its use should be visible to the caller.
  if (!missing(measure) && missing(user_es_target_measure_adj)) {
    warning("The 'measure' argument of es_from_user_adj() is deprecated; use 'user_es_target_measure_adj'.",
            call. = FALSE)
    user_es_target_measure_adj <- measure
  }

  len <- length(user_es_original_measure_adj)

  if (missing(user_es_adj)) user_es_adj <- rep(NA_real_, len)
  if (missing(user_se_adj)) user_se_adj <- rep(NA_real_, len)
  if (missing(user_ci_lo_adj)) user_ci_lo_adj <- rep(NA_real_, len)
  if (missing(user_ci_up_adj)) user_ci_up_adj <- rep(NA_real_, len)
  if (missing(n_exp)) n_exp <- rep(NA_real_, len)
  if (missing(n_nexp)) n_nexp <- rep(NA_real_, len)
  if (missing(n_sample)) n_sample <- rep(NA_real_, len)
  if (missing(n_cases)) n_cases <- rep(NA_real_, len)
  if (missing(n_controls)) n_controls <- rep(NA_real_, len)
  if (missing(baseline_risk)) baseline_risk <- rep(NA_real_, len)
  if (missing(small_margin_prop)) small_margin_prop <- rep(NA_real_, len)

  es_val <- ifelse(is.na(user_es_adj) & !is.na(user_es_original_measure_adj) &
                    !is.na(user_ci_lo_adj) & !is.na(user_ci_up_adj),
    (user_ci_up_adj + user_ci_lo_adj) / 2,
    user_es_adj
  )

  se_val <- ifelse(is.na(user_se_adj) & !is.na(user_es_original_measure_adj) &
                    !is.na(user_ci_lo_adj) & !is.na(user_ci_up_adj),
    (user_ci_up_adj - user_ci_lo_adj) / (2 * qnorm(.975)),
    user_se_adj
  )

  ci_lo <- ifelse(is.na(user_ci_lo_adj) & !is.na(user_es_original_measure_adj) &
                   !is.na(user_es_adj) & !is.na(user_se_adj),
    user_es_adj - user_se_adj * qnorm(.975),
    user_ci_lo_adj
  )

  ci_up <- ifelse(is.na(user_ci_up_adj) & !is.na(user_es_original_measure_adj) &
                   !is.na(user_es_adj) & !is.na(user_se_adj),
    user_es_adj + user_se_adj * qnorm(.975),
    user_ci_up_adj
  )

  target <- unique(user_es_target_measure_adj)
  if (length(target) != 1) target <- target[1]

  orig_type <- tolower(user_es_original_measure_adj)
  known_types <- .known_es_types()

  ratio_from_ci <- which(orig_type %in% c("or", "rr", "irr", "hr") &
                          is.na(user_se_adj) &
                          !is.na(user_ci_lo_adj) & !is.na(user_ci_up_adj) &
                          user_ci_lo_adj > 0 & user_ci_up_adj > 0)
  if (length(ratio_from_ci) > 0) {
    log_lo <- log(user_ci_lo_adj[ratio_from_ci])
    log_up <- log(user_ci_up_adj[ratio_from_ci])
    se_val[ratio_from_ci] <- (log_up - log_lo) / (2 * qnorm(.975))
    es_val[ratio_from_ci] <- ifelse(
      is.na(user_es_adj[ratio_from_ci]),
      (log_lo + log_up) / 2,
      log(user_es_adj[ratio_from_ci])
    )
    ci_lo[ratio_from_ci] <- log_lo
    ci_up[ratio_from_ci] <- log_up
    ot <- orig_type[ratio_from_ci]
    orig_type[ratio_from_ci] <- ifelse(ot == "or",  "logor",
                                ifelse(ot == "rr",  "logrr",
                                ifelse(ot == "irr", "logirr",
                                ifelse(ot == "hr",  "loghr", ot))))
  }

  irr_with_se <- which(orig_type == "irr" & !is.na(se_val))
  if (length(irr_with_se) > 0) {
    es_val[irr_with_se] <- log(es_val[irr_with_se])
    orig_type[irr_with_se] <- "logirr"
  }

  hr_with_se <- which(orig_type == "hr" & !is.na(se_val))
  if (length(hr_with_se) > 0) {
    es_val[hr_with_se] <- log(es_val[hr_with_se])
    orig_type[hr_with_se] <- "loghr"
  }

  ratio_se_check <- which(orig_type %in% c("or", "rr") &
                           !is.na(user_se_adj) & !is.na(user_es_adj) &
                           user_es_adj > 0 & user_se_adj > user_es_adj)
  if (length(ratio_se_check) > 0) {
    warning(paste0(
      "For rows with user_es_original_measure_adj = 'or'/'rr', ",
      "user_se_adj must be on the LOG scale (SE of logOR/logRR), ",
      "not the natural scale. Detected SE > ES in ",
      length(ratio_se_check), " row(s), which suggests the SE may be ",
      "on the natural scale. If so, provide CI bounds instead and the ",
      "log-scale SE will be derived automatically."
    ))
  }

  has_known_type <- !is.na(orig_type) & orig_type %in% known_types
  has_unknown_type <- !is.na(orig_type) & !(orig_type %in% known_types)

  if (any(has_unknown_type)) {
    unknown_vals <- unique(orig_type[has_unknown_type])
    warning(paste0(
      "Unknown user_es_original_measure_adj values: ",
      paste(unknown_vals, collapse = ", "),
      ". No conversion applied for these rows. Known types: ",
      paste(known_types, collapse = ", "), "."
    ))
  }

  res <- data.frame(info_used = rep("user_input_adj", len))

  res[[target]] <- rep(NA_real_, len)
  res[[paste0(target, "_se")]] <- rep(NA_real_, len)
  res[[paste0(target, "_ci_lo")]] <- rep(NA_real_, len)
  res[[paste0(target, "_ci_up")]] <- rep(NA_real_, len)

  passthrough_idx <- which(!has_known_type & !is.na(es_val))
  if (length(passthrough_idx) > 0) {
    res[[target]][passthrough_idx] <- es_val[passthrough_idx]
    res[[paste0(target, "_se")]][passthrough_idx] <- se_val[passthrough_idx]
    res[[paste0(target, "_ci_lo")]][passthrough_idx] <- ci_lo[passthrough_idx]
    res[[paste0(target, "_ci_up")]][passthrough_idx] <- ci_up[passthrough_idx]
  }

  for (es_type in unique(orig_type[has_known_type])) {
    idx <- which(has_known_type & orig_type == es_type)
    if (length(idx) == 0) next

    converted <- .dispatch_user_conversion(
      original_type = es_type,
      es_val = es_val[idx], se_val = se_val[idx],
      n_exp = n_exp[idx], n_nexp = n_nexp[idx],
      n_sample = n_sample[idx],
      n_cases = n_cases[idx], n_controls = n_controls[idx],
      baseline_risk = baseline_risk[idx],
      small_margin_prop = small_margin_prop[idx],
      or_to_rr = or_to_rr, or_to_cor = or_to_cor,
      smd_to_cor = smd_to_cor, cor_to_smd = cor_to_smd,
      rr_to_or = rr_to_or
    )

    if (!is.null(converted)) {
      for (col in names(converted)) {
        if (col == "info_used") next
        if (!col %in% names(res)) res[[col]] <- rep(NA_real_, len)
        res[[col]][idx] <- converted[[col]]
      }
    }
  }

  # restore user ci on the entered measure
  pres <- .user_preserved_ci(user_es_original_measure_adj,
                              user_ci_lo_adj, user_ci_up_adj,
                              es_val, se_val)
  for (es_type in unique(orig_type[has_known_type])) {
    prefixes <- .user_slot_prefixes(es_type)
    if (is.null(prefixes)) next
    idx <- which(has_known_type & orig_type == es_type &
                 !is.na(pres$lo) & !is.na(pres$up))
    if (es_type == "nnt") idx <- idx[pres$user_gave_ci[idx]]
    if (length(idx) == 0) next
    for (prefix in prefixes) {
      lo_col <- paste0(prefix, "_ci_lo")
      up_col <- paste0(prefix, "_ci_up")
      if (!lo_col %in% names(res) || !up_col %in% names(res)) next
      res[[lo_col]][idx] <- pres$lo[idx]
      res[[up_col]][idx] <- pres$up[idx]
    }
  }

  attr(res, "measure") <- user_es_original_measure_adj
  return(res)
}
