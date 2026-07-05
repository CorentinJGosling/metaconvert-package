#' Meta-analytic dataset inspired from Haza and colleagues (2024)
#'
#' Dataset of a meta-analysis exploring the specificity of social functioning of children with ADHD (compared to healthy controls) in case-control studies.
#' This dataset contains:
#' 1. several information coming from the same participants (due to the completion of multiple outcomes).
#' 1. several information coming from the same study (due to the presence of multiple subgroups).
#' 1. overlapping information for the same comparison
#' 1. several information types from which a standardized mean difference can be estimated/converted
#'
#' @source Haza B, Gosling CJ, Conty L & Pinabiaux C (2024). Social Functioning in Children and Adolescents with ADHD: A Meta-analysis. Journal of Child Psychology and Psychiatry and Allied Disciplines.
"df.haza"

#' Short version of the df.haza dataset
#'
#' This dataset is a shoter version of the \code{\link{df.haza}} dataset.
#'
#' @source Haza B, Gosling CJ, Conty L & Pinabiaux C (2024). Social Functioning in Children and Adolescents with ADHD: A Meta-analysis. Journal of Child Psychology and Psychiatry and Allied Disciplines.
"df.short"


#' Fictitious dataset 1
#'
#' First fictitious dataset aiming to understand how the \code{\link{compare_df}} function works.
#' Slightly different from \code{\link{df.compare1}}
#'
"df.compare1"

#' Fictitious dataset 2
#'
#' First fictitious dataset aiming to understand how the \code{\link{compare_df}} function works.
#' Slightly different from \code{\link{df.compare2}}
#'
"df.compare2"

#' Simulated dataset for a COSMIN-based systematic review of PROMs
#'
#' Simulated dataset of a systematic review evaluating measurement properties
#' of a fictional patient-reported outcome measure (the Mental Health Wellbeing
#' Scale, MHWS) following the COSMIN framework. Designed to demonstrate the
#' psychometric functions of the metaConvert package.
#'
#' The dataset contains studies reporting:
#' \itemize{
#'   \item Internal consistency (Cronbach's alpha)
#'   \item Test-retest reliability (ICC)
#'   \item Criterion and construct validity (Pearson and Spearman correlations)
#'   \item Responsiveness (change-score correlations)
#'   \item Floor/ceiling effects (proportions)
#' }
#'
#' Additional columns support standalone psychometric utility functions
#' (SEM, SDC, disattenuation, change-score reliability).
#'
"df.psychom"
