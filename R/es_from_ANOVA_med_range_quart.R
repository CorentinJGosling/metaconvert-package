# Bounds-safe lookup into a Wan (2014) / Hozo (2005) constant table. At index 0 R
# silently DROPS the element (and does negative indexing at index < 0), which shortens
# the result vector so the enclosing ifelse() recycles a constant onto the WRONG (valid)
# row -- silently corrupting a neighbouring study's SD/effect size. Mapping every
# out-of-range index to NA keeps the result the same length as the input, so each row
# keeps its own constant and an out-of-range row is NA (n <= 1 for the IQR/quartile
# tables, n <= 0 for the range table). This complements the Tier-1 V16 check in
# internal_flags.R, which flags n == 1 rows as invalid (within-group variance undefined).
.wan_table_lookup <- function(tbl, idx) {
  idx <- as.integer(idx)
  idx[!is.na(idx) & (idx < 1L | idx > length(tbl))] <- NA_integer_
  tbl[idx]
}

#' Convert median and interquartile range of two independent groups into several effect size measures
#'
#' @param q1_exp first quartile of the experimental/exposed group.
#' @param med_exp median value of the experimental/exposed group.
#' @param q3_exp third quartile of the experimental/exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param q1_nexp first quartile of the non-experimental/non-exposed group.
#' @param med_nexp median value of the non-experimental/non-exposed group.
#' @param q3_nexp third quartile of the non-experimental/non-exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
#' @param reverse_med a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first converts a Cohen's d (D), Hedges' g (G) and mean difference (MD)
#' from the medians and interquartile ranges of two independent groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **This function recreates means+SD of the two groups**  (Wan et al., 2014):
#' \deqn{mean\_exp = \frac{q1\_exp + med\_exp + q3\_exp}{3}}
#' \deqn{mean\_nexp = \frac{q1\_nexp + med\_nexp + q3\_nexp}{3}}
#' \deqn{mean\_sd\_exp = \frac{q3\_exp - q1\_exp}{2*qnorm(\frac{0.75*n\_exp - 0.125}{n\_exp+0.25})}}
#' \deqn{mean\_sd\_nexp = \frac{q3\_nexp - q1\_nexp}{2*qnorm(\frac{0.75*n\_nexp - 0.125}{n\_nexp+0.25})}}
#'
#' The denominator shown above is a large-sample approximation. When \eqn{Q = (n - 1)/4}
#' is no greater than 50 - that is, for group sample sizes up to 201 - the tabulated
#' constant of Wan et al. (2014) is used in its place.
#'
#' **From these means+SD, the function computes MD, D and G** using formulas
#' described in \code{\link{es_from_means_sd}()}.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @references
#' Wan, X., Wang, W., Liu, J. et al. Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. BMC Med Res Methodol 14, 135 (2014). https://doi.org/10.1186/1471-2288-14-135
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' @export es_from_med_quarts
#'
#' @md
#'
#' @examples
#' es_from_med_quarts(
#'   q1_exp = 1335, med_exp = 1400,
#'   q3_exp = 1765, n_exp = 40,
#'   q1_nexp = 1481, med_nexp = 1625,
#'   q3_nexp = 1800, n_nexp = 40
#' )
es_from_med_quarts <- function(q1_exp, med_exp, q3_exp, n_exp,
                               q1_nexp, med_nexp, q3_nexp, n_nexp,
                               smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_med) {
  if (missing(reverse_med)) reverse_med <- rep(FALSE, length(q1_exp))
  reverse_med[is.na(reverse_med)] <- FALSE


  e <- c(
    0.990, 1.144, 1.206, 1.239, 1.260,
    1.274, 1.284, 1.292, 1.298, 1.303,
    1.307, 1.311, 1.313, 1.316, 1.318,
    1.320, 1.322, 1.323, 1.324, 1.326,
    1.327, 1.328, 1.329, 1.330, 1.330,
    1.331, 1.332, 1.332, 1.333, 1.333,
    1.334, 1.334, 1.335, 1.335, 1.336,
    1.336, 1.336, 1.337, 1.337, 1.337,
    1.338, 1.338, 1.338, 1.338, 1.339,
    1.339, 1.339, 1.339, 1.339, 1.340
  )

  mean_exp <- (q1_exp + med_exp + q3_exp) / 3
  Q_exp <- (n_exp - 1) / 4
  cor_exp <- ifelse(Q_exp <= 50,
    .wan_table_lookup(e, ceiling(Q_exp)),
    (2 * qnorm((0.75 * n_exp - 0.125) / (n_exp + 0.25)))
  )
  sd_exp <- (q3_exp - q1_exp) / cor_exp

  mean_nexp <- (q1_nexp + med_nexp + q3_nexp) / 3
  Q_nexp <- (n_nexp - 1) / 4
  cor_nexp <- ifelse(Q_nexp <= 50,
    .wan_table_lookup(e, ceiling(Q_nexp)),
    (2 * qnorm((0.75 * n_nexp - 0.125) / (n_nexp + 0.25)))
  )
  sd_nexp <- (q3_nexp - q1_nexp) / cor_nexp

  es <- es_from_means_sd(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_sd_exp = sd_exp, mean_sd_nexp = sd_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_means = reverse_med
  )

  es$info_used <- "med_quarts"

  return(es)
}

#' Convert median, range and interquartile range of two independent groups into several effect size measures
#'
#' @param min_exp minimum value of the experimental/exposed group.
#' @param q1_exp first quartile of the experimental/exposed group.
#' @param med_exp median value of the experimental/exposed group.
#' @param q3_exp third quartile of the experimental/exposed group.
#' @param max_exp maximum value of the experimental/exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param min_nexp minimum value of the non-experimental/non-exposed group.
#' @param q1_nexp first quartile of the non-experimental/non-exposed group.
#' @param med_nexp median value of the non-experimental/non-exposed group.
#' @param q3_nexp third quartile of the non-experimental/non-exposed group.
#' @param max_nexp maximum value of the non-experimental/non-exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
#' @param reverse_med a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first converts a Cohen's d (D), Hedges' g (G) and mean difference (MD)
#' from the medians, ranges, and interquartile ranges of two independent groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **This function recreates means+SD of the two groups**  (Wan et al., 2014):
#' \deqn{mean\_exp = \frac{min\_exp + 2*q1\_exp + 2*med\_exp + 2*q3\_exp + max\_exp}{8}}
#' \deqn{mean\_nexp = \frac{min\_nexp + 2*q1\_nexp + 2*med\_nexp + 2*q3\_nexp + max\_nexp}{8}}
#' \deqn{mean\_sd\_exp = \frac{max\_exp - min\_exp}{4*qnorm(\frac{n\_exp-0.375}{n\_exp+0.25})} + \frac{q3\_exp-q1\_exp}{4*qnorm(\frac{0.75*n\_exp-0.125}{n\_exp+0.25})}}
#' \deqn{mean\_sd\_nexp = \frac{max\_nexp - min\_nexp}{4*qnorm(\frac{n\_nexp-0.375}{n\_nexp+0.25})} + \frac{q3\_nexp-q1\_nexp}{4*qnorm(\frac{0.75*n\_nexp-0.125}{n\_nexp+0.25})}}
#'
#' The two normal-quantile factors shown above are large-sample approximations, and each
#' is replaced by a tabulated constant of Wan et al. (2014) at small samples: the range
#' factor for group sample sizes no greater than 50, and the interquartile factor when
#' \eqn{Q = (n - 1)/4} is no greater than 50 - that is, for group sample sizes up to 201.
#'
#' **From these means+SD, the function computes MD, D and G** using formulas
#' described in \code{\link{es_from_means_sd}()}.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab \cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab MD + D + G\cr
#'  \tab OR + R + Z \cr
#'  \code{required input data} \tab See 'Section 12. Median, range and/or interquartile range'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Wan, X., Wang, W., Liu, J. et al. Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. BMC Med Res Methodol 14, 135 (2014). https://doi.org/10.1186/1471-2288-14-135
#'
#' @export es_from_med_min_max_quarts
#'
#' @md
#'
#' @examples
#' es_from_med_min_max_quarts(
#'   min_exp = 1102, q1_exp = 1335,
#'   med_exp = 1400, q3_exp = 1765,
#'   max_exp = 1899, n_exp = 40,
#'   min_nexp = 1181, q1_nexp = 1481,
#'   med_nexp = 1625, q3_nexp = 1800,
#'   max_nexp = 1910, n_nexp = 40
#' )
es_from_med_min_max_quarts <- function(q1_exp, med_exp, q3_exp,
                                       min_exp, max_exp, n_exp,
                                       q1_nexp, med_nexp, q3_nexp,
                                       min_nexp, max_nexp, n_nexp,
                                       smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_med) {
  if (missing(reverse_med)) reverse_med <- rep(FALSE, length(min_exp))
  reverse_med[is.na(reverse_med)] <- FALSE

  list1 <- c(
    0.000, 1.128, 1.693, 2.059, 2.326, 2.534, 2.704,
    2.847, 2.970, 3.078, 3.173, 3.259, 3.336, 3.407,
    3.472, 3.532, 3.588, 3.640, 3.689, 3.735, 3.778,
    3.819, 3.858, 3.895, 3.931, 3.964, 3.997, 4.027,
    4.057, 4.086, 4.113, 4.139, 4.165, 4.189, 4.213,
    4.236, 4.259, 4.280, 4.301, 4.322, 4.341, 4.361,
    4.379, 4.398, 4.415, 4.433, 4.450, 4.466, 4.482,
    4.498
  )

  list2 <- c(
    0.990, 1.144, 1.206, 1.239, 1.260,
    1.274, 1.284, 1.292, 1.298, 1.303,
    1.307, 1.311, 1.313, 1.316, 1.318,
    1.320, 1.322, 1.323, 1.324, 1.326,
    1.327, 1.328, 1.329, 1.330, 1.330,
    1.331, 1.332, 1.332, 1.333, 1.333,
    1.334, 1.334, 1.335, 1.335, 1.336,
    1.336, 1.336, 1.337, 1.337, 1.337,
    1.338, 1.338, 1.338, 1.338, 1.339,
    1.339, 1.339, 1.339, 1.339, 1.340
  )


  mean_exp <- (min_exp + 2 * q1_exp + 2 * med_exp + 2 * q3_exp + max_exp) / 8
  cor_m_exp <- ifelse(n_exp <= 50,
    .wan_table_lookup(list1, n_exp),
    2 * qnorm((n_exp - 0.375) / (n_exp + 0.25))
  )
  Q_exp <- (n_exp - 1) / 4
  cor_q_exp <- ifelse(Q_exp <= 50,
    .wan_table_lookup(list2, ceiling(Q_exp)),
    (2 * qnorm((0.75 * n_exp - 0.125) / (n_exp + 0.25)))
  )

  sd_exp <- ((max_exp - min_exp) / cor_m_exp + (q3_exp - q1_exp) / cor_q_exp) / 2

  mean_nexp <- (min_nexp + 2 * q1_nexp + 2 * med_nexp + 2 * q3_nexp + max_nexp) / 8
  cor_m_nexp <- ifelse(n_nexp <= 50,
    .wan_table_lookup(list1, n_nexp),
    2 * qnorm((n_nexp - 0.375) / (n_nexp + 0.25))
  )
  Q_nexp <- (n_nexp - 1) / 4
  cor_q_nexp <- ifelse(Q_nexp <= 50,
    .wan_table_lookup(list2, ceiling(Q_nexp)),
    (2 * qnorm((0.75 * n_nexp - 0.125) / (n_nexp + 0.25)))
  )

  sd_nexp <- ((max_nexp - min_nexp) / cor_m_nexp + (q3_nexp - q1_nexp) / cor_q_nexp) / 2

  es <- es_from_means_sd(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_sd_exp = sd_exp, mean_sd_nexp = sd_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_means = reverse_med
  )

  es$info_used <- "med_min_max_quarts"

  return(es)
}

#' Convert median and range of two independent groups into several effect size measures
#'
#' @param min_exp minimum value of the experimental/exposed group.
#' @param med_exp median value of the experimental/exposed group.
#' @param max_exp maximum value of the experimental/exposed group.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param min_nexp minimum value of the non-experimental/non-exposed group.
#' @param med_nexp median value of the non-experimental/non-exposed group.
#' @param max_nexp maximum value of the non-experimental/non-exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
#' @param reverse_med a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first converts a Cohen's d (D), Hedges' g (G)
#' and mean difference (MD) from the medians and ranges of two independent groups.
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **This function recreates means+SD of the two groups**  (Wan et al., 2014):
#' \deqn{mean\_exp = \frac{min\_exp + 2*med\_exp + max\_exp}{4} + \frac{min\_exp - 2*med\_exp + max\_exp}{4*n\_exp}}
#' \deqn{mean\_nexp = \frac{min\_nexp + 2*med\_nexp + max\_nexp}{4} + \frac{min\_nexp - 2*med\_nexp + max\_nexp}{4*n\_nexp}}
#' \deqn{mean\_sd\_exp = \frac{max\_exp - min\_exp}{2*qnorm((n\_exp-0.375) / (n\_exp+0.25))}}
#' \deqn{mean\_sd\_nexp = \frac{max\_nexp - min\_nexp}{2*qnorm((n\_nexp-0.375) / (n\_nexp+0.25))}}
#'
#' The denominator shown above is a large-sample approximation. For group sample sizes
#' no greater than 50, the tabulated constant of Wan et al. (2014) is used in its place.
#'
#' **From these means+SD, the function computes MD, D and G** using formulas
#' described in \code{\link{es_from_means_sd}()}.
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' The authors of the Cochrane Handbook stated
#' "As a general rule, we recommend that ranges should not be used
#' to estimate SDs." (see section 6.5.2.6).
#' It is thus a good practice to explore the consequences of
#' the use of this conversion in sensitivity analyses.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab \cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab MD + D + G\cr
#'  \tab OR + R + Z \cr
#'  \code{required input data} \tab See 'Section 12. Median, range and/or interquartile range'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Wan, X., Wang, W., Liu, J. et al. Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. BMC Med Res Methodol 14, 135 (2014). https://doi.org/10.1186/1471-2288-14-135
#'
#' Higgins JPT, Li T, Deeks JJ (editors). Chapter 6: Choosing effect size measures and computing estimates of effect. In: Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Cochrane, 2022. Available from www.training.cochrane.org/handbook.
#'
#' @export es_from_med_min_max
#'
#' @md
#'
#' @examples
#' es_from_med_min_max(
#'   min_exp = 1335, med_exp = 1400,
#'   max_nexp = 1765, n_exp = 40,
#'   min_nexp = 1481, med_nexp = 1625,
#'   max_exp = 1800, n_nexp = 40
#' )
es_from_med_min_max <- function(min_exp, med_exp, max_exp, n_exp,
                                min_nexp, med_nexp, max_nexp, n_nexp,
                                smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_med) {
  if (missing(reverse_med)) reverse_med <- rep(FALSE, length(min_exp))
  reverse_med[is.na(reverse_med)] <- FALSE

  list1 <- c(
    0.000, 1.128, 1.693, 2.059, 2.326, 2.534, 2.704,
    2.847, 2.970, 3.078, 3.173, 3.259, 3.336, 3.407,
    3.472, 3.532, 3.588, 3.640, 3.689, 3.735, 3.778,
    3.819, 3.858, 3.895, 3.931, 3.964, 3.997, 4.027,
    4.057, 4.086, 4.113, 4.139, 4.165, 4.189, 4.213,
    4.236, 4.259, 4.280, 4.301, 4.322, 4.341, 4.361,
    4.379, 4.398, 4.415, 4.433, 4.450, 4.466, 4.482,
    4.498
  )

  mean_exp <- (min_exp + 2 * med_exp + max_exp) / 4 + (min_exp - 2 * med_exp + max_exp) / (4 * n_exp)
  cor_exp <- ifelse(n_exp <= 50, .wan_table_lookup(list1, n_exp), 2 * qnorm((n_exp - 0.375) / (n_exp + 0.25)))
  sd_exp <- (max_exp - min_exp) / cor_exp

  mean_nexp <- (min_nexp + 2 * med_nexp + max_nexp) / 4 + (min_nexp - 2 * med_nexp + max_nexp) / (4 * n_nexp)
  cor_nexp <- ifelse(n_nexp <= 50, .wan_table_lookup(list1, n_nexp), 2 * qnorm((n_nexp - 0.375) / (n_nexp + 0.25)))
  sd_nexp <- (max_nexp - min_nexp) / cor_nexp

  es <- es_from_means_sd(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_sd_exp = sd_exp, mean_sd_nexp = sd_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_means = reverse_med
  )

  es$info_used <- "med_min_max"

  return(es)
}
