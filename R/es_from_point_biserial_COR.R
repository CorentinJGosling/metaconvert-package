#' Convert a point-biserial correlation coefficient into several effect size measures
#'
#' @param pt_bis_r value of a point-biserial correlation coefficient
#' @param n_exp total number of participants in the exposed group
#' @param n_nexp total number of participants in the non exposed group
#' @param smd_to_cor formula used to convert the \code{pt_bis_r} value into a coefficient correlation.
#' @param reverse_pt_bis_r a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function uses a point-biserial correlation coefficient to estimate a
#' Cohen's d (D) and Hedges' g (G).
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **The formula used to obtain the Cohen's d are (Viechtbauer, 2021)**:
#' \deqn{m = n\_exp + n\_nexp - 2}
#' \deqn{h = \frac{m}{n\_exp} + \frac{m}{n\_nexp}}
#' \deqn{d = \frac{pt\_bis\_r * \sqrt{h}}{\sqrt{1 - pt\_bis\_r^2}}}
#'
#' **To estimate other effect size measures**,
#' calculations of the \code{\link{es_from_cohen_d}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 11. ANOVA statistics, Student's t-test, or point-bis correlation'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Viechtbauer (2021). Accessed at: https://stats.stackexchange.com/questions/526789/convert-correlation-r-to-cohens-d-unequal-groups-of-known-size
#'
#' @export es_from_pt_bis_r
#'
#' @md
#'
#' @examples
#' es_from_pt_bis_r(pt_bis_r = 0.2, n_exp = 121, n_nexp = 121)
es_from_pt_bis_r <- function(pt_bis_r, n_exp, n_nexp, smd_to_cor = "viechtbauer", reverse_pt_bis_r) {
  if (missing(reverse_pt_bis_r)) reverse_pt_bis_r <- rep(FALSE, length(pt_bis_r))
  reverse_pt_bis_r[is.na(reverse_pt_bis_r)] <- FALSE
  if (length(reverse_pt_bis_r) == 1) reverse_pt_bis_r = c(rep(reverse_pt_bis_r, length(pt_bis_r)))
  if (length(reverse_pt_bis_r) != length(pt_bis_r)) stop("The length of the 'reverse_pt_bis_r' argument is incorrectly specified.")

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       error_message = paste0("The number of people exposed/non-exposed, ",
                                              "should be >0."),
                       func = "es_from_pt_bis_r")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  df <- n_exp + n_nexp - 2
  h <- df / n_exp + df / n_nexp
  p <- n_exp / (n_exp + n_nexp)
  q <- n_nexp / (n_exp + n_nexp)

  d <- pt_bis_r * sqrt(h) / sqrt(1 - pt_bis_r^2)
  d <- ifelse(reverse_pt_bis_r, -d, d)
  es <- .es_from_d(d = d, n_exp = n_exp, n_nexp = n_nexp)

  es$info_used <- "pt_bis_r"

  return(es)
}

#' Convert a p-value of a point-biserial correlation coefficient into several effect size measures
#'
#' @param pt_bis_r_pval p-value of a point-biserial correlation coefficient
#' @param n_exp total number of participants in the exposed group
#' @param n_nexp total number of participants in the non exposed group
#' @param smd_to_cor formula used to convert the \code{pt_bis_r_pval} value into a coefficient correlation.
#' @param reverse_pt_bis_r_pval a logical value indicating whether the direction of the generated effect sizes should be flipped.
#'
#' @details
#' This function converts the p-value of a point biserial correlation into a Student's t-value.
#'
#' **The formula used to obtain this Student's t-value is**:
#' \deqn{t = pt(\frac{pt\_bis\_r\_pval}{2}, df = n\_exp + n\_nexp - 2)}
#'
#' Calculations of the \code{\link{es_from_student_t}} function are then applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 11. ANOVA statistics, Student's t-test, or point-bis correlation'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' @export es_from_pt_bis_r_pval
#'
#' @md
#'
#' @examples
#' es_from_pt_bis_r_pval(pt_bis_r_pval = 0.2, n_exp = 121, n_nexp = 121)
es_from_pt_bis_r_pval <- function(pt_bis_r_pval, n_exp, n_nexp,
                                  smd_to_cor = "viechtbauer", reverse_pt_bis_r_pval) {
  if (missing(reverse_pt_bis_r_pval)) reverse_pt_bis_r_pval <- rep(FALSE, length(n_exp))
  reverse_pt_bis_r_pval[is.na(reverse_pt_bis_r_pval)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       error_message = paste0("The number of people exposed/non-exposed, and p-values ",
                                              "should be >0."),
                       func = "es_from_pt_bis_r_pval")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  t <- qt(p = pt_bis_r_pval / 2, df = n_exp + n_nexp - 2, lower.tail = FALSE)

  es <- es_from_student_t(
    student_t = t, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_student_t = reverse_pt_bis_r_pval
  )

  es$info_used <- "pt_bis_r_pval"

  return(es)
}
