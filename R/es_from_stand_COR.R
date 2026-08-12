#' Convert a Pearson's correlation coefficient to several effect size measures
#'
#' @param pearson_r a Pearson's correlation coefficient value
#' @param n_sample the total number of participants
#' @param sd_iv the standard deviation of the independent variable
#' @param unit_increase_iv a value of the independent variable that will be used to estimate the Cohen's d (see details).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Must be either "sd" or "value"
#' @param n_exp number of the experimental/exposed group
#' @param n_nexp number of the non-experimental/non-exposed group
#' @param cor_to_smd formula used to convert a \code{pearson_r} or \code{fisher_z} value into a SMD.
#' @param reverse_pearson_r a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function estimates the variance of a Pearson's correlation coefficient, and computes the Fisher's r-to-z
#' transformation.
#' Cohen's d (D), Hedges' g (G) are converted from the Pearson's r, and odds ratio (OR)
#' are converted from the Cohen's d.
#'
#' 1. **The formula used to estimate the standard error of the Pearson's correlation coefficient** and 95% CI
#' are (Formula 12.27 in Cooper):
#' \deqn{R\_se = \sqrt{\frac{(1 - pearson\_r^2)^2}{n\_sample - 1}}}
#' \deqn{R\_lo = pearson\_r - qt(.975, n\_sample - 2) * R\_se}
#' \deqn{R\_up = pearson\_r + qt(.975, n\_sample - 2) * R\_se}
#'
#' 2. **The formula used to estimate the Fisher's z** are (Formula 12.28 & 12.29 in Cooper):
#' \deqn{Z = atanh(r)}
#' \deqn{Z\_se = \sqrt{\frac{1}{n\_sample - 3}}}
#' \deqn{Z\_ci\_lo = Z - qnorm(.975) * Z\_se}
#' \deqn{Z\_ci\_up = Z + qnorm(.975) * Z\_se}
#'
#' 3. Several approaches can be used to convert a correlation coefficient to a SMD.
#'
#' **A.** Mathur proposes to use this formula (Formula 1.2 in Mathur, \code{cor_to_smd = "mathur"}):
#' \deqn{increase = ifelse(unit_type == "sd", unit\_increase\_iv * sd\_iv, unit\_increase\_iv)}
#' \deqn{d = \frac{r * increase}{sd_iv * \sqrt{1 - r^2}}}
#' \deqn{d\_se = abs(d) * \sqrt{\frac{1}{r^2 * (n\_sample - 3)} + \frac{1}{2*(n\_sample - 1))}}}
#' The resulting Cohen's d is the average increase in the dependent variable associated with an increase of x units in the independent variable (with x = \code{unit_increase_iv}).
#'
#' **B.** Viechtbauer proposes to use the delta method to derive a Cohen's d from a correlation coefficient (Viechtbauer, 2023, \code{cor_to_smd = "viechtbauer"})
#'
#' **C.** Cooper converts a point-biserial correlation to Cohen's d with the fixed
#' per-2-SD transformation (Formula 12.38 & 12.39 in Cooper, \code{cor_to_smd = "cooper"}).
#' Unlike the \code{"mathur"} option, this branch does NOT use \code{unit_increase_iv},
#' \code{unit_type} or \code{sd_iv}; it always returns the standardized mean difference
#' between two groups lying one predictor standard deviation apart:
#' \deqn{d = \frac{2 * r}{\sqrt{1 - r^2}}}
#' \deqn{d\_se = \sqrt{\frac{4 * R\_se^2}{(1 - r^2)^3}}}
#' where \eqn{R\_se} is the standard error of the correlation. It therefore coincides
#' with the \code{cor_to_smd = "mathur"} result only in the special case
#' \code{unit_type = "sd"} and \code{unit_increase_iv = 2}.
#'
#' To know how the Cohen's d value is converted to other effect measures (G/OR), see details of the \code{\link{es_from_cohen_d}} function.
#'
#' @export es_from_pearson_r
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Mathur, M. B., & VanderWeele, T. J. (2020). A Simple, Interpretable Conversion from Pearson's Correlation to Cohen's for d Continuous Exposures. Epidemiology (Cambridge, Mass.), 31(2), e16-e18. https://doi.org/10.1097/EDE.0000000000001105
#'
#' Viechtbauer W (2010). "Conducting meta-analyses in R with the metafor package." Journal of Statistical Software, 36(3), 1-48. doi:10.18637/jss.v036.i03.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab R + Z\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + OR\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 4. Pearson's r or Fisher's z'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @examples
#' es_from_pearson_r(
#'   pearson_r = .51, sd_iv = 0.24, n_sample = 214,
#'   unit_increase_iv = 1, unit_type = "sd"
#' )
es_from_pearson_r <- function(pearson_r, sd_iv, n_sample,
                              n_exp, n_nexp, cor_to_smd = "viechtbauer",
                              unit_increase_iv, unit_type = "raw_scale", reverse_pearson_r) {
  if (missing(reverse_pearson_r)) {
    reverse_pearson_r <- rep(FALSE, length(pearson_r))
  }
  reverse_pearson_r[is.na(reverse_pearson_r)] <- FALSE
  if (length(reverse_pearson_r) == 1) reverse_pearson_r = c(rep(reverse_pearson_r, length(pearson_r)))
  if (length(reverse_pearson_r) != length(pearson_r)) stop("The length of the 'reverse_pearson_r' argument is incorrectly specified.")

  if (missing(n_exp)) {
    n_exp <- rep(NA, length(pearson_r))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA, length(pearson_r))
  }
  if (missing(n_sample)) {
    n_sample <- rep(NA, length(pearson_r))
  }
  if (missing(sd_iv)) {
    sd_iv <- rep(NA, length(pearson_r))
  }
  if (missing(unit_increase_iv)) {
    unit_increase_iv <- rep(NA, length(pearson_r))
  }
  if (missing(unit_type)) {
    unit_type <- rep(NA, length(pearson_r))
  }

  if (!all(cor_to_smd %in% c("cooper", "mathur", "viechtbauer"))) {
    stop(paste0("'",
                unique(cor_to_smd[!cor_to_smd %in% c("cooper", "mathur", "viechtbauer")]),
                "' not in tolerated values for the 'cor_to_smd' argument.
                Possible inputs are: 'cooper', 'mathur', 'viechtbauer'"))
  }

  # pearson_r = dat$pearson_r; sd_iv = unit_increase_iv = rep(NA, length(pearson_r))#dat$sd_iv;
  # n_sample = dat$n_sample;
  # n_exp = dat$n_exp; n_nexp = dat$n_nexp; cor_to_smd = "viechtbauer";
  # unit_type = "raw_scale"; reverse_pearson_r = rep(FALSE, length(pearson_r))

  n_sample <- ifelse(is.na(n_sample), n_exp + n_nexp, n_sample)
  n_exp <- ifelse(is.na(n_exp), n_sample / 2, n_exp)
  n_nexp <- ifelse(is.na(n_nexp), n_sample / 2, n_nexp)

  r <- ifelse(reverse_pearson_r, -pearson_r, pearson_r)
  r_se <- sqrt((1 - r^2)^2 / (n_sample - 1))

  dat_cor <- data.frame(
    r = r, r_se = r_se,
    sd_iv = sd_iv, n_sample = n_sample,
    unit_increase_iv = unit_increase_iv,
    unit_type = unit_type,
    cor_to_smd = cor_to_smd
  )

  nn_miss <- with(dat_cor, which(
    (cor_to_smd == "mathur" & !is.na(r) & !is.na(sd_iv) &
      !is.na(n_sample) & !is.na(unit_increase_iv) & !is.na(unit_type)) |
      (cor_to_smd == "cooper" & !is.na(r) & !is.na(r_se)) |
      (cor_to_smd == "viechtbauer" & !is.na(r) & !is.na(r_se) & !is.na(n_sample))
  ))

  es <- data.frame(
    d = rep(NA, nrow(dat_cor)),
    d_se = rep(NA, nrow(dat_cor))
  )

  if (length(nn_miss) != 0) {
    res_d <- .cor_to_smd_vec(
      r = dat_cor$r[nn_miss],
      r_se = dat_cor$r_se[nn_miss],
      n_sample = dat_cor$n_sample[nn_miss],
      sd_iv = dat_cor$sd_iv[nn_miss],
      unit_increase_iv = dat_cor$unit_increase_iv[nn_miss],
      unit_type = dat_cor$unit_type[nn_miss],
      cor_to_smd = dat_cor$cor_to_smd[nn_miss]
    )

    es$d[nn_miss] <- unlist(res_d[, 1])
    es$d_se[nn_miss] <- unlist(res_d[, 2])
  }

  es <- .es_from_d(d = es$d, d_se = es$d_se, n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample)

  es$r <- r
  es$r_se <- r_se
  es$r_ci_lo <- r - qt(.975, n_sample - 2) * r_se
  es$r_ci_up <- r + qt(.975, n_sample - 2) * r_se

  es$z <- atanh(r)
  es$z_se <- sqrt(1 / (n_sample - 3))
  es$z_ci_lo <- es$z - qnorm(.975) * es$z_se
  es$z_ci_up <- es$z + qnorm(.975) * es$z_se

  es$info_used <- "pearson_r"
  return(es)
}

#' Convert a Fisher's z (r-to-z transformation) to several effect size measures
#'
#' @param fisher_z a Fisher's r-to-z transformed correlation coefficient
#' @param n_sample the total number of participants
#' @param sd_iv the standard deviation of the independent variable
#' @param unit_increase_iv a value of the independent variable that will be used to estimate the Cohen's d (see details).
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Must be either "sd" or "value"
#' @param reverse_fisher_z a logical value indicating whether the direction of the generated effect sizes should be flipped.
#' @param n_exp number of the experimental/exposed group
#' @param n_nexp number of the non-experimental/non-exposed group
#' @param cor_to_smd formula used to convert a \code{pearson_r} or \code{fisher_z} value into a SMD.
#'
#' @details
#' This function converts estimates the standard error of the Fisher's z and performs the z-to-r Fisher's transformation.
#'
#' Last, it converts this r value into a Cohen's d and OR (see details in \code{\link{es_from_pearson_r}()}).
#'
#'
#' @export es_from_fisher_z
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Mathur, M. B., & VanderWeele, T. J. (2020). A Simple, Interpretable Conversion from Pearson's Correlation to Cohen's for d Continuous Exposures. Epidemiology (Cambridge, Mass.), 31(2), e16-e18. https://doi.org/10.1097/EDE.0000000000001105
#'
#' Viechtbauer W (2010). "Conducting meta-analyses in R with the metafor package." Journal of Statistical Software, 36(3), 1-48. doi:10.18637/jss.v036.i03.
#'
#' @md
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab R + Z\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + OR\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 4. Pearson's r or Fisher's z'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @examples
#' es_from_fisher_z(
#'   fisher_z = .21, n_sample = 44,
#' )
es_from_fisher_z <- function(fisher_z, n_sample, unit_type = "raw_scale",
                             n_exp, n_nexp, cor_to_smd = "viechtbauer",
                             sd_iv, unit_increase_iv, reverse_fisher_z) {
  if (missing(reverse_fisher_z)) {
    reverse_fisher_z <- rep(FALSE, length(fisher_z))
  }
  reverse_fisher_z[is.na(reverse_fisher_z)] <- FALSE
  if (missing(n_exp)) {
    n_exp <- rep(NA, length(fisher_z))
  }
  if (missing(n_nexp)) {
    n_nexp <- rep(NA, length(fisher_z))
  }
  if (missing(n_sample)) {
    n_sample <- rep(NA, length(fisher_z))
  }
  if (missing(sd_iv)) {
    sd_iv <- rep(NA, length(fisher_z))
  }
  if (missing(unit_increase_iv)) {
    unit_increase_iv <- rep(NA, length(fisher_z))
  }
  if (missing(unit_type)) {
    unit_type <- rep(NA, length(fisher_z))
  }
  r <- tanh(fisher_z)

  es <- es_from_pearson_r(
    pearson_r = r, sd_iv = sd_iv, n_sample = n_sample,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    unit_increase_iv = unit_increase_iv, unit_type = unit_type,
    reverse_pearson_r = reverse_fisher_z
  )

  es$info_used <- "fisher_z"
  return(es)
}
