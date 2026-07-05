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

.es_from_d <- function(d, d_se, n_exp, n_nexp, n_sample, smd_to_cor = "viechtbauer",
                       adjusted, n_cov_ancova, cov_outcome_r, reverse) {
  if (missing(d_se)) d_se <- rep(NA_real_, length(d))
  if (missing(n_exp)) n_exp <- rep(NA_real_, length(d))
  if (missing(n_nexp)) n_nexp <- rep(NA_real_, length(d))
  if (missing(n_sample)) n_sample <- rep(NA_real_, length(d))
  if (missing(adjusted)) adjusted <- rep(FALSE, length(d))
  if (missing(n_cov_ancova)) n_cov_ancova <- rep(0, length(d))
  if (missing(cov_outcome_r)) cov_outcome_r <- rep(0.5, length(d))
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
  d_se <- ifelse(
    !is.na(d_se),
    d_se,
    ifelse(adjusted,
      sqrt(((n_exp + n_nexp) / (n_exp * n_nexp) * (1 - cov_outcome_r^2)) + d^2 / (2 * (n_exp + n_nexp))),
      sqrt((n_exp + n_nexp) / (n_exp * n_nexp) + d^2 / (2 * (n_exp + n_nexp)))
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
  dat_r <- data.frame(
    d = d, vd = d_se^2, n_exp = n_exp, n_nexp = n_nexp,
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
