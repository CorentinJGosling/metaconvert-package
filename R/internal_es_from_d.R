.normalize_smd_var <- function(smd_var) {
  # Map the user-facing names for the SMD sampling-variance convention onto the
  # internal canonical labels "LS" and "LS2" (metafor's vtype names). The
  # user-facing names follow the package convention of author names (as for
  # smd_to_cor = "viechtbauer", pre_post_to_smd = "bonett", etc.):
  #   "borenstein" (default, current metaConvert behaviour) = Borenstein et al.
  #     (2009, Introduction to Meta-Analysis, ch. 4): v_d = 1/n1 + 1/n2 + d^2/(2N)
  #     then v_g = cm^2 * v_d  (bit-exact with metafor vtype = "LS2").
  #   "hedges_olkin" (opt-in) = Hedges & Olkin (1985) / Viechtbauer (2007),
  #     metafor's own default: v_g = 1/n1 + 1/n2 + g^2/(2N), v_d = v_g / cm^2. This
  #     unifies the endpoint/between-group family with the pre/post family.
  # The metafor codes "LS2"/"LS" are kept as accepted synonyms.
  lut <- c(
    "borenstein"  = "LS2",
    "hedges_olkin" = "LS", "hedges-olkin" = "LS", "viechtbauer" = "LS"
  )
  key <- tolower(as.character(smd_var))
  # NA (e.g. a blank cell in a per-row smd_var column on a direct call) falls
  # back to the package default, mirroring .normalize_smd_denom's NA -> "pooled".
  key[is.na(smd_var)] <- "borenstein"
  out <- unname(lut[key])
  if (any(is.na(out))) {
    stop(paste0(
      "'", paste(unique(smd_var[is.na(out)]), collapse = "', '"),
      "' not in tolerated values for the 'smd_var' argument. Possible inputs are: ",
      "'borenstein' (the default) or 'hedges_olkin' (alias 'viechtbauer')."
    ), call. = FALSE)
  }
  out
}

.normalize_smd_denom <- function(smd_denom) {
  # Map the user-facing names for the SMD standardiser onto the internal canonical
  # labels "pooled" / "control" / "control_robust". Glass's delta (the control-SD
  # standardiser) follows the package's author/name convention as "glass"; the
  # descriptive "control" / "control_robust" spellings are kept as accepted synonyms.
  # Returns NA for unrecognised values so the caller can raise a targeted error.
  lut <- c(
    "pooled" = "pooled",
    "glass" = "control", "control" = "control",
    "glass_robust" = "control_robust", "control_robust" = "control_robust"
  )
  unname(lut[tolower(as.character(smd_denom))])
}

.d_j <- function(x) {
  # Hedges' small-sample correction J(df) is only defined for df > 1.
  # When df <= 1, the correction (and hence Hedges' g) is returned as NA.
  bad <- which(!is.na(x) & x <= 1)
  if (length(bad) > 0) {
    warning(sprintf(
      "Hedges' J correction is undefined for df <= 1; %d row(s) had df <= 1 and will produce NA for Hedges' g. This usually reflects extremely small samples (n_exp + n_nexp - 2 <= 1).",
      length(bad)
    ), call. = FALSE)
  }
  j <- ifelse(x <= 1, NA, 1) * exp(lgamma(x / 2) - 0.5 * log(x / 2) - lgamma((x - 1) / 2))
  return(j)
}

.glass_from_means <- function(mean_exp, mean_sd_exp, mean_nexp, mean_sd_nexp,
                              n_exp, n_nexp, smd_to_cor = "viechtbauer",
                              smd_var = "borenstein", reverse, robust = FALSE) {
  # Glass's delta: the endpoint mean difference standardised by the control
  # (non-experimental) endpoint SD instead of the pooled SD. Bit-exact with
  # metafor measure = "SMD1" (homoscedastic variance, with the vtype matched to
  # smd_var: "borenstein" = LS2, "hedges_olkin" = LS) and "SMD1H" (robust,
  # heteroscedasticity-consistent variance). glass_robust/SMD1H has a single
  # variance form, so smd_var is ignored for it. The Hedges small-sample
  # correction and every variance term use the control-group degrees of freedom
  # (n_nexp - 1) rather than the pooled df, which is what distinguishes SMD1
  # from SMD.
  if (missing(reverse)) reverse <- rep(FALSE, length(mean_exp))
  reverse[is.na(reverse)] <- FALSE
  smd_var <- .normalize_smd_var(smd_var)
  if (length(smd_var) == 1) smd_var <- rep(smd_var, length(mean_exp))

  di   <- (mean_exp - mean_nexp) / mean_sd_nexp           # uncorrected Glass delta
  df_g <- n_nexp - 1
  cm_g <- .d_j(df_g)                                      # cm(n_nexp - 1)
  g    <- cm_g * di

  # g-scale sampling variance
  vg_hom <- ifelse(smd_var == "LS",
    1 / n_exp + 1 / n_nexp + g^2 / (2 * n_nexp),                          # SMD1, vtype "LS"
    cm_g^2 * (1 / n_exp + 1 / n_nexp + di^2 / (2 * n_nexp))               # SMD1, vtype "LS2"
  )
  vg_rob <- (mean_sd_exp^2 / mean_sd_nexp^2) / (n_exp - 1) +
    1 / (n_nexp - 1) + g^2 / (2 * (n_nexp - 1))                           # SMD1H
  vg   <- if (robust) vg_rob else vg_hom
  g_se <- sqrt(vg)
  d_se <- g_se / cm_g                                     # d-scale SE (g_se = cm * d_se)

  # r / z / logOR / d scaffolding from the uncorrected Glass delta + its SE
  es <- .es_from_d(
    d = di, d_se = d_se, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse
  )

  # override g with the control-df correction (.es_from_d used the pooled df)
  g_signed <- ifelse(reverse, -g, g)
  es$g      <- g_signed
  es$g_se   <- g_se
  es$g_ci_lo <- g_signed - qt(.975, df_g) * g_se
  es$g_ci_up <- g_signed + qt(.975, df_g) * g_se
  # keep the d CI on the same (control) df as g for internal consistency
  es$d_ci_lo <- es$d - qt(.975, df_g) * es$d_se
  es$d_ci_up <- es$d + qt(.975, df_g) * es$d_se
  es
}

.es_from_d <- function(d, d_se, n_exp, n_nexp, n_sample, smd_to_cor = "viechtbauer",
                       adjusted, n_cov_ancova, cov_outcome_r, reverse,
                       smd_var = "borenstein") {
  if (missing(d_se)) d_se <- rep(NA_real_, length(d))
  if (missing(n_exp)) n_exp <- rep(NA_real_, length(d))
  if (missing(n_nexp)) n_nexp <- rep(NA_real_, length(d))
  if (missing(n_sample)) n_sample <- rep(NA_real_, length(d))
  if (missing(adjusted)) adjusted <- rep(FALSE, length(d))
  # An omitted n_cov_ancova means "no covariates beyond the grouping", i.e. q = 0: it
  # only costs degrees of freedom, so the crude default is a safe reading and keeps the
  # returned row internally coherent.
  if (missing(n_cov_ancova)) n_cov_ancova <- rep(0, length(d))
  # cov_outcome_r is different: it scales the Cooper eq. 12.26 (1 - R^2) variance
  # shrink, and there is no safe default value for it. It is read only on the adjusted
  # branch, so the 0.5 placeholder below is never used on the crude branch and stays.
  # On the adjusted branch it must not be invented, because R's missing() propagates
  # through an unforced promise: es_from_cohen_d_adj() is the one adjusted entry point
  # where cov_outcome_r never enters the point estimate and so is never forced, which
  # would let it fall through to R = 0.5 and return a confidently wrong SE (weight
  # inflated about 29%) while its 13 sibling routes hard-error on the same omission.
  # NA propagates to d_se and from there to every variance-bearing output, blanking the
  # row coherently, which matches what convert_df() returns for an absent or NA
  # column.
  .is_adj <- rep_len(adjusted, length(d)) %in% TRUE
  if (missing(cov_outcome_r)) {
    cov_outcome_r <- rep(0.5, length(d))
    cov_outcome_r[.is_adj] <- NA_real_
  }
  if (missing(reverse)) reverse <- rep(FALSE, length(d))
  reverse[is.na(reverse)] <- FALSE
  if (length(reverse) == 1) reverse = c(rep(reverse, length(d)))
  if (length(reverse) != length(d)) stop("The length of the 'reverse' argument is incorrectly specified.")
  if (length(adjusted) == 1) adjusted = c(rep(adjusted, length(d)))
  if (length(adjusted) != length(d)) stop("The length of the 'adjusted' argument is incorrectly specified.")


  if (!all(smd_to_cor %in% c("viechtbauer", "lipsey_cooper"))) {
    stop(paste0(
      "'",
      unique(smd_to_cor[!smd_to_cor %in% c("viechtbauer", "lipsey_cooper")]),
      "' not in tolerated values for the 'smd_to_cor' argument.",
      " Possible inputs are: 'viechtbauer', 'lipsey_cooper'"
    ))
  }

  smd_var <- .normalize_smd_var(smd_var)
  if (length(smd_var) == 1) smd_var <- rep(smd_var, length(d))
  if (length(smd_var) != length(d)) stop("The length of the 'smd_var' argument is incorrectly specified.")
  use_LS <- smd_var == "LS"

  # ========= FLIP THE EFFECT SIZE ========== #
  d <- ifelse(reverse, -d, d)
  # ========= homogeneize sample sizes ========== #
  n_sample <- ifelse(!is.na(n_sample), n_sample, n_exp + n_nexp)
  n_exp <- ifelse(!is.na(n_exp), n_exp, n_sample / 2)
  n_nexp <- ifelse(!is.na(n_nexp), n_nexp, n_sample / 2)

  df <- ifelse(adjusted,
    n_exp + n_nexp - 2 - n_cov_ancova,
    n_exp + n_nexp - 2
  )

  # ========= d_se ========= #
  # Large-sample default variance (only used for rows where d_se was not supplied).
  # The leading term is 1/n1 + 1/n2 (adjusted for the covariate multiple correlation
  # in the ANCOVA case). g_se = d_se * J holds in both conventions (see below), so the
  # only thing smd_var changes is the *default* d_se:
  #   LS2 (default): v_d = leading + d^2/(2N)                  -> v_g = cm^2 * v_d
  #   LS  (opt-in) : v_g = leading + g^2/(2N), v_d = v_g/cm^2  (g = d * J, metafor default)
  leading <- ifelse(adjusted,
    (n_exp + n_nexp) / (n_exp * n_nexp) * (1 - cov_outcome_r^2),
    (n_exp + n_nexp) / (n_exp * n_nexp)
  )
  if (any(use_LS)) {
    J_var <- suppressWarnings(.d_j(df))
    d_se_ls <- sqrt((leading + (d * J_var)^2 / (2 * (n_exp + n_nexp))) / J_var^2)
  } else {
    d_se_ls <- rep(NA_real_, length(d))
  }
  # Recorded before the default fills the gaps: the r/z conversion below has to tell a
  # standard error the caller supplied (which carries real design information) from one
  # this function invented (which carries none). vd_ls2 is that invented variance on the
  # LS2 convention, kept separately because the conversion needs it whichever convention
  # smd_var selected -- see the note on dat_r.
  d_se_supplied <- !is.na(d_se)
  vd_ls2 <- leading + d^2 / (2 * (n_exp + n_nexp))
  d_se <- ifelse(
    d_se_supplied,
    d_se,
    ifelse(use_LS,
      d_se_ls,
      sqrt(vd_ls2)
    )
  )

  d_ci_lo <- d - d_se * qt(.975, df)
  d_ci_up <- d + d_se * qt(.975, df)

  # ========= OR ========= #
  logor <- d * pi / sqrt(3)
  logor_se <- sqrt(d_se^2 * pi^2 / 3)
  logor_ci_lo <- logor - logor_se * qnorm(.975)
  logor_ci_up <- logor + logor_se * qnorm(.975)

  #######################################################
  # track non missing information for speedy conversion #
  #######################################################
  nn_miss <- which(!is.na(d) & !is.na(d_se) & !is.na(n_exp) & !is.na(n_nexp))

  g <- g_se <- g_ci_lo <- g_ci_up <-
    r <- r_se <- r_ci_lo <- r_ci_up <-
    z <- z_se <- z_ci_lo <- z_ci_up <- rep(NA, length(d))

  # ========= g ========= #
  J <- .d_j(df[nn_miss])
  g[nn_miss] <- d[nn_miss] * J
  g_se[nn_miss] <- sqrt(d_se[nn_miss]^2 * (J^2))
  g_ci_lo[nn_miss] <- g[nn_miss] - g_se[nn_miss] * qt(.975, df[nn_miss])
  g_ci_up[nn_miss] <- g[nn_miss] + g_se[nn_miss] * qt(.975, df[nn_miss])


  # ========= r/Z ========= #
  # .smd_to_cor() rescales Soper's closed-form biserial variance by vd / vd_crude, and
  # its vd_crude is the LS2 expression, so the ratio is exactly 1 on an ordinary crude
  # row only when vd is on that same convention. A supplied d_se carries precision the
  # closed form cannot know about (the ANCOVA (1 - R^2) shrink, the pre-post
  # 2(1 - r_pre_post) factor, the control-group df of a Glass row, a user-reported
  # standard error) and must pass through untouched. A defaulted one carries none, so it
  # is handed over on the LS2 convention: otherwise smd_var = "hedges_olkin" would leak
  # the pure estimator-convention factor 1/J^2 into r_se and z_se (11.6% of the variance
  # at N = 16, 5.5% at N = 30, 0.74% at N = 200), taking r_se off metafor's exact RBIS
  # variance and z_se^2 off the stabilised 1/(N - 1). Only what the conversion sees
  # changes; the returned d_se/g_se stay on whichever convention smd_var selected.
  dat_r <- data.frame(
    d = d, vd = ifelse(d_se_supplied, d_se^2, vd_ls2),
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, n_cov_ancova = n_cov_ancova
  )

  if (length(nn_miss) != 0) {
    cor <- t(mapply(.smd_to_cor,
      d = dat_r$d[nn_miss],
      vd = dat_r$vd[nn_miss],
      n_exp = dat_r$n_exp[nn_miss],
      n_nexp = dat_r$n_nexp[nn_miss],
      smd_to_cor = dat_r$smd_to_cor[nn_miss],
      n_cov_ancova = dat_r$n_cov_ancova[nn_miss]
    ))

    r[nn_miss] <- cor[, 1]
    r_se[nn_miss] <- sqrt(cor[, 2])
    r_ci_lo[nn_miss] <- cor[, 3]
    r_ci_up[nn_miss] <- cor[, 4]
    z[nn_miss] <- cor[, 5]
    z_se[nn_miss] <- sqrt(cor[, 6])
    z_ci_lo[nn_miss] <- cor[, 7]
    z_ci_up[nn_miss] <- cor[, 8]
  }

  res <- data.frame(
    d, d_se, d_ci_lo, d_ci_up,
    g, g_se, g_ci_lo, g_ci_up,
    r, r_se, r_ci_lo, r_ci_up,
    z, z_se, z_ci_lo, z_ci_up,
    logor, logor_se, logor_ci_lo, logor_ci_up
  )

  return(res)
}
