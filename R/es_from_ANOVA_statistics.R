#' Convert a Student's t-test value to several effect size measures
#'
#' @param student_t Student's t-test value.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{student_t} value into a coefficient correlation (see details).
#' @param reverse_student_t a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the Student's t-test value into a Cohen's d (D) and Hedges' g (G),
#' Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.
#'
#' **To estimate a Cohen's d** the formula used is (table 12.1 in Cooper):
#' \deqn{cohen\_d = student\_t * \sqrt{\frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp}}}
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
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 11. ANOVA statistics, Student's t-test, or point-bis correlation'\cr
#'  \tab https://metaconvert.org/html/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_student_t
#'
#' @md
#'
#' @examples
#' es_from_student_t(student_t = 2.1, n_exp = 20, n_nexp = 22)
es_from_student_t <- function(student_t, n_exp, n_nexp,
                              smd_to_cor = "viechtbauer", reverse_student_t) {
  if (missing(reverse_student_t)) reverse_student_t <- rep(FALSE, length(student_t))
  reverse_student_t[is.na(reverse_student_t)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp,
                       error_message = paste0("The number of people exposed/non-exposed, ",
                                              "should be >0."),
                       func = "es_from_student_t")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })


  d <- student_t * sqrt(1 / n_exp + 1 / n_nexp)

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse = reverse_student_t
  )

  es$info_used <- "student_t"
  return(es)
}

#' Convert a Student's t-test p-value to several effect size measures
#'
#' @param student_t_pval p-value (two-tailed) from a Student's t-test. If your p-value is one-tailed, simply multiply it by two.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{student_t_pval} value into a coefficient correlation (see details).
#' @param reverse_student_t_pval a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the Student's t-test p-value into a t-value,
#' and then relies on the calculations of the \code{\link{es_from_student_t}()} function.
#'
#' **To convert the p-value into a t-value,** the following formula is used (table 12.1 in Cooper):
#' \deqn{student\_t = qt(\frac{student\_t\_pval}{2}, df = n\_exp + n\_nexp - 2)}
#' Then, calculations of the \code{\link{es_from_student_t}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 11. ANOVA statistics, Student's t-test, or point-bis correlation'\cr
#'  \tab https://metaconvert.org/html/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_student_t_pval
#'
#' @md
#'
#' @examples
#' es_from_student_t_pval(student_t_pval = 0.24, n_exp = 20, n_nexp = 22)
es_from_student_t_pval <- function(student_t_pval, n_exp, n_nexp,
                                   smd_to_cor = "viechtbauer", reverse_student_t_pval) {
  if (missing(reverse_student_t_pval)) reverse_student_t_pval <- rep(FALSE, length(student_t_pval))
  reverse_student_t_pval[is.na(reverse_student_t_pval)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp, student_t_pval,
                       error_message = paste0("The number of people exposed/non-exposed, ",
                                              " and the p-value ",
                                              "should be >0."),
                       func = "es_from_student_t_pval")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  t <- qt(p = student_t_pval / 2, df = n_exp + n_nexp - 2, lower.tail = FALSE)

  es <- es_from_student_t(
    student_t = t, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_student_t = reverse_student_t_pval
  )

  es$info_used <- "student_t_pval"

  return(es)
}

#' Convert a one-way independent ANOVA F-value to several effect size measures
#'
#' @param anova_f ANOVA F-value (one-way, binary predictor).
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{anova_f} value into a coefficient correlation (see details).
#' @param reverse_anova_f a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the F-value (one-way, binary predictor) into a t-value,
#' and then relies on the calculations of the \code{\link{es_from_student_t}()} function.
#'
#' **To convert the F-value into a t-value,** the following formula is used (table 12.1 in Cooper):
#' \deqn{student\_t = \sqrt{anova\_f}}
#' Then, calculations of the \code{\link{es_from_student_t}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 11. ANOVA statistics, Student's t-test, or point-bis correlation'\cr
#'  \tab https://metaconvert.org/html/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_anova_f
#'
#' @md
#'
#' @examples
#' es_from_anova_f(anova_f = 2.01, n_exp = 20, n_nexp = 22)
es_from_anova_f <- function(anova_f, n_exp, n_nexp, smd_to_cor = "viechtbauer", reverse_anova_f) {
  if (missing(reverse_anova_f)) reverse_anova_f <- rep(FALSE, length(anova_f))
  reverse_anova_f[is.na(reverse_anova_f)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp, anova_f,
                       error_message = paste0("The number of people exposed/non-exposed, ",
                                              " and the F-ANOVA ",
                                              "should be >0."),
                       func = "es_from_anova_f")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  t <- sqrt(anova_f)

  es <- es_from_student_t(
    student_t = t, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_student_t = reverse_anova_f
  )

  es$info_used <- "anova_f"

  return(es)
}

#' Convert a p-value from a one-way independent ANOVA to several effect size measures
#'
#' @param anova_f_pval p-value (two-tailed) from an ANOVA (binary predictor). If your p-value is one-tailed, simply multiply it by two.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the \code{anova_f_pval} value into a coefficient correlation (see details).
#' @param reverse_anova_f_pval a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function converts the p-value from the F-value of an ANOVA (one-way, binary predictor) into a t-value,
#' and then relies on the calculations of the \code{\link{es_from_student_t}()} function.
#'
#' **To convert the p-value into a t-value,** the following formula is used (table 12.1 in Cooper):
#' \deqn{student\_t = qt(\frac{anova\_f\_pval}{2}, df = n\_exp + n\_nexp - 2)}
#' Then, calculations of the \code{\link{es_from_student_t}()} are applied.
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab D + G\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab OR + R + Z\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 11. ANOVA statistics, Student's t-test, or point-bis correlation'\cr
#'  \tab https://metaconvert.org/html/input.html\cr
#'  \tab \cr
#' }
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' @export es_from_anova_pval
#'
#' @md
#'
#' @examples
#' es_from_anova_pval(anova_f_pval = 0.0012, n_exp = 20, n_nexp = 22)
es_from_anova_pval <- function(anova_f_pval, n_exp, n_nexp, smd_to_cor = "viechtbauer",
                               reverse_anova_f_pval) {
  if (missing(reverse_anova_f_pval)) reverse_anova_f_pval <- rep(FALSE, length(anova_f_pval))
  reverse_anova_f_pval[is.na(reverse_anova_f_pval)] <- FALSE

  tryCatch({
    .validate_positive(n_exp, n_nexp, anova_f_pval,
                       error_message = paste0("The number of people exposed/non-exposed, ",
                                              " and the p-value ",
                                              "should be >0."),
                       func = "es_from_anova_pval")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })


  es <- es_from_student_t_pval(
    student_t_pval = anova_f_pval, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_student_t_pval = reverse_anova_f_pval
  )

  es$info_used <- "anova_f_pval"

  return(es)
}
