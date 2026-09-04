#' Convert a Student's t-test value to several effect size measures
#'
#' @param student_t Student's t-test value.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
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
#' **Important - Student's t, not Welch's t.** This formula is the exact inverse of the
#' pooled-variance (Student) t-test, so \code{student_t} must be the *equal-variance* t.
#' It should not be used with a Welch (unequal-variance) t-test, which is the default of
#' R's \code{t.test()}. A Welch t uses \eqn{\sqrt{s_1^2/n_1 + s_2^2/n_2}} rather than
#' \eqn{s_{pooled}\sqrt{1/n_1 + 1/n_2}}, so when group sizes and variances both differ,
#' plugging a Welch t into this formula yields a biased Cohen's d (the bias can exceed
#' 50% and may flip sign depending on which arm carries the larger variance; it vanishes
#' only when \eqn{n_1 = n_2}). A Welch t cannot be converted to Cohen's d from
#' \code{(t, n_exp, n_nexp)} alone: recovering the pooled SD requires the two arm SDs, in
#' which case \code{\link{es_from_means_sd}} should be used directly.
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
#'  \tab https://metaconvert.org/input.html\cr
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
                              smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_student_t) {
  if (missing(reverse_student_t)) reverse_student_t <- rep(FALSE, length(student_t))
  reverse_student_t[is.na(reverse_student_t)] <- FALSE


  d <- student_t * sqrt(1 / n_exp + 1 / n_nexp)

  es <- .es_from_d(
    d = d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse = reverse_student_t
  )

  es$info_used <- "student_t"
  return(es)
}

#' Convert a Student's t-test p-value to several effect size measures
#'
#' @param student_t_pval p-value (two-tailed) from a Student's t-test. If your p-value is one-tailed, simply multiply it by two.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
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
#' Note that a two-sided p-value carries no direction, so the recovered t (and hence the
#' generated effect sizes) are always non-negative. Use \code{reverse_student_t_pval} to
#' encode the correct sign for effects that favour the non-experimental group. The same
#' equal-variance (Student, not Welch) assumption as \code{\link{es_from_student_t}} applies.
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
#'  \tab https://metaconvert.org/input.html\cr
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
                                   smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_student_t_pval) {
  if (missing(reverse_student_t_pval)) reverse_student_t_pval <- rep(FALSE, length(student_t_pval))
  reverse_student_t_pval[is.na(reverse_student_t_pval)] <- FALSE

  # p <= 0 inverts to t = Inf (d = Inf, se = Inf); p = 1 is the ordinary t = 0 and is kept.
  student_t_pval <- .pval_or_na(student_t_pval, se_from_ratio = FALSE)
  t <- qt(p = student_t_pval / 2, df = n_exp + n_nexp - 2, lower.tail = FALSE)

  es <- es_from_student_t(
    student_t = t, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_student_t = reverse_student_t_pval
  )

  es$info_used <- "student_t_pval"

  return(es)
}

#' Convert a one-way independent ANOVA F-value to several effect size measures
#'
#' @param anova_f ANOVA F-value (one-way, binary predictor).
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
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
#' **Important - single numerator degree of freedom only.** The identity
#' \eqn{\sqrt{F} = |t|} holds only when the F-test has a single numerator degree of freedom,
#' i.e. a one-way ANOVA comparing exactly two groups (a binary predictor). Supplying an
#' omnibus F from a factor with three or more levels (numerator df > 1) produces a
#' meaningless effect size and is not detected by the function. In addition, \eqn{\sqrt{F}}
#' discards the sign of the effect, so the generated effect sizes are always non-negative;
#' use \code{reverse_anova_f} to encode direction.
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
#'  \tab https://metaconvert.org/input.html\cr
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
es_from_anova_f <- function(anova_f, n_exp, n_nexp, smd_to_cor = "viechtbauer", smd_var = "borenstein", reverse_anova_f) {
  if (missing(reverse_anova_f)) reverse_anova_f <- rep(FALSE, length(anova_f))
  reverse_anova_f[is.na(reverse_anova_f)] <- FALSE

  t <- sqrt(anova_f)

  es <- es_from_student_t(
    student_t = t, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_student_t = reverse_anova_f
  )

  es$info_used <- "anova_f"

  return(es)
}

#' Convert a p-value from a one-way independent ANOVA to several effect size measures
#'
#' @param anova_f_pval p-value (two-tailed) from an ANOVA (binary predictor). If your p-value is one-tailed, simply multiply it by two.
#' @param n_exp number of participants in the experimental/exposed group.
#' @param n_nexp number of participants in the non-experimental/non-exposed group.
#' @param smd_to_cor formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).
#' @param smd_var name of the sampling-variance formula for the standardized mean difference: "borenstein" (default) or "hedges_olkin" (alias "viechtbauer"). The two differ by a squared small-sample-correction factor (J^2); "hedges_olkin" is a few percent larger at small samples. This choice affects the standard error only: the effect size itself is identical under both formulas. The standardizer is selected with \code{smd_denom}, which does change the effect size.
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
#' As for \code{\link{es_from_anova_f}}, this conversion is valid only for an F-test with a
#' single numerator degree of freedom (a two-group comparison): the p-value of a
#' multi-level (numerator df > 1) omnibus F is inverted here as if it were a two-sided
#' two-group t p-value, which is incorrect. The two-sided p-value also carries no
#' direction, so the generated effect sizes are always non-negative; use
#' \code{reverse_anova_f_pval} to encode the correct sign.
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
#'  \tab https://metaconvert.org/input.html\cr
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
                               smd_var = "borenstein", reverse_anova_f_pval) {
  if (missing(reverse_anova_f_pval)) reverse_anova_f_pval <- rep(FALSE, length(anova_f_pval))
  reverse_anova_f_pval[is.na(reverse_anova_f_pval)] <- FALSE


  es <- es_from_student_t_pval(
    student_t_pval = anova_f_pval, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, smd_var = smd_var, reverse_student_t_pval = reverse_anova_f_pval
  )

  es$info_used <- "anova_f_pval"

  return(es)
}
