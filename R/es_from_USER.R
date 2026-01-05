#' Directly input a value + variance of an effect size measure
#'
#' @param measure the effect size measure used in calculations (must be one of the 11 effect size measures available in metaConvert)
#' @param user_es_measure_crude the name of the effect size measure that will appear when this function is called by the \link{convert_df} function (can be any character string)
#' @param user_es_crude effect size value
#' @param user_se_crude standard error of the effect size
#' @param user_ci_lo_crude lower bound of the 95% CI around the effect size value
#' @param user_ci_up_crude upper bound of the 95% CI around the effect size value
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#'
#' @details
#' This function is a generic function allowing to include any crude effect size measure value + variance.
#' Importantly, with this function, no conversions are performed (i.e., the effect size value + variance
#' you enter is the value + variance exported by this function).
#'
#' @md
#'
#' @export es_from_user_crude
#'
#' @return
#' This function allows to directly input any of the 11 effect size measures
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab Any of the 11 available measures\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab No conversion performed\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 23. User's input (crude)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @examples
#' dat = data.frame(measure = "OR", user_es_measure_crude = "mortality rate ratio",
#'                  user_es_crude = -0.04, user_se_crude = 0.2)
#' summary(convert_df(dat, measure="logor"))
es_from_user_crude <- function(measure, user_es_measure_crude, user_es_crude, user_se_crude,
                                user_ci_lo_crude, user_ci_up_crude, max_asymmetry = 10) {
  if (missing(user_es_crude)) user_es_crude <- rep(NA_real_, length(user_es_measure_crude))
  if (missing(user_se_crude)) user_se_crude <- rep(NA_real_, length(user_es_measure_crude))
  if (missing(user_ci_lo_crude)) user_ci_lo_crude <- rep(NA_real_, length(user_es_measure_crude))
  if (missing(user_ci_up_crude)) user_ci_up_crude <- rep(NA_real_, length(user_es_measure_crude))

  tryCatch({
    .validate_positive(user_se_crude,
                       error_message = paste0("The standard error ",
                                              "should be >0."),
                       func = "es_from_user_crude")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  tryCatch({
    .validate_ci_symmetry(user_es_crude, user_ci_lo_crude, user_ci_up_crude,
                          func = "es_from_user_crude",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })


  es <- ifelse(is.na(user_es_crude) & !is.na(user_es_measure_crude) & !is.na(user_ci_lo_crude) & !is.na(user_ci_up_crude),
    (user_ci_up_crude + user_ci_lo_crude) / 2,
    user_es_crude
  )

  se <- ifelse(is.na(user_se_crude) & !is.na(user_es_measure_crude) & !is.na(user_ci_lo_crude) & !is.na(user_ci_up_crude),
    (user_ci_up_crude - user_ci_lo_crude) / (2 * qnorm(.975)),
    user_se_crude
  )

  ci_lo <- ifelse(is.na(user_ci_lo_crude) & !is.na(user_es_measure_crude) & !is.na(user_es_crude) & !is.na(user_se_crude),
    user_es_crude - user_se_crude * qnorm(.975),
    user_ci_lo_crude
  )

  ci_up <- ifelse(is.na(user_ci_up_crude) & !is.na(user_es_measure_crude) & !is.na(user_es_crude) & !is.na(user_se_crude),
    user_es_crude + user_se_crude * qnorm(.975),
    user_ci_up_crude
  )

  res <- data.frame(info_used = rep(NA, length(user_es_measure_crude)))
  res$info_used <- "user_input_crude"
  res[, paste(unique(measure))] <- es
  res[, paste(unique(measure), "_se", sep = "")] <- se
  res[, paste(unique(measure), "_ci_lo", sep = "")] <- ci_lo
  res[, paste(unique(measure), "_ci_up", sep = "")] <- ci_up

  attr(res, "measure") <- user_es_measure_crude
  return(res)
}

#' Directly input an adjusted value + variance of an effect size measure
#'
#' @param measure the effect size measure used in calculations (must be one of the 11 effect size measures available in metaConvert)
#' @param user_es_measure_adj the name of the effect size measure that will appear when this function is called by the \link{convert_df} function (can be any character string)
#' @param user_es_adj adjusted effect size value
#' @param user_se_adj adjusted standard error of the effect size
#' @param user_ci_lo_adj adjusted lower bound of the 95% CI around the effect size value
#' @param user_ci_up_adj adjusted upper bound of the 95% CI around the effect size value
#' @param max_asymmetry A percentage indicating the tolerance before detecting asymmetry in the 95% CI bounds.
#'
#' @details
#' This function is a generic function allowing to include any adjusted effect size measure value + variance.
#' Importantly, with this function, no conversions are performed (i.e., the effect size value + variance
#' you enter is the value + variance exported by this function).
#'
#' @md
#'
#' @export es_from_user_adj
#'
#' @return
#' This function allows to directly input any of the 11 effect size measures
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab Any of the 11 available measures\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab No conversion performed\cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 24. User's input (adjusted)'\cr
#'  \tab https://metaconvert.org/input.html\cr
#'  \tab \cr
#' }
#'
#' @examples
#' dat = data.frame(measure = "OR", user_es_measure_adj = "adjusted OR",
#'                  user_es_adj = -0.04, user_se_adj = 0.2)
#' summary(convert_df(dat, measure="logor"))
es_from_user_adj <- function(measure, user_es_measure_adj, user_es_adj, user_se_adj,
                              user_ci_lo_adj, user_ci_up_adj, max_asymmetry = 10) {
  if (missing(user_es_adj)) user_es_adj <- rep(NA_real_, length(user_es_measure_adj))
  if (missing(user_se_adj)) user_se_adj <- rep(NA_real_, length(user_es_measure_adj))
  if (missing(user_ci_lo_adj)) user_ci_lo_adj <- rep(NA_real_, length(user_es_measure_adj))
  if (missing(user_ci_up_adj)) user_ci_up_adj <- rep(NA_real_, length(user_es_measure_adj))

  tryCatch({
    .validate_positive(user_se_adj,
                       error_message = paste0("The standard error ",
                                              "should be >0."),
                       func = "es_from_user_adj")
  }, error = function(e) {
    stop("Data entry error: ", conditionMessage(e), "\n")
  })

  tryCatch({
    .validate_ci_symmetry(user_es_adj, user_ci_lo_adj, user_ci_up_adj,
                          func = "es_from_user_adj",
                          max_asymmetry_percent = max_asymmetry)
  }, error = function(e) {
    stop("Validation failed: ", conditionMessage(e), "\n")
  })

  es <- ifelse(is.na(user_es_adj) & !is.na(user_es_measure_adj) & !is.na(user_ci_lo_adj) & !is.na(user_ci_up_adj),
    (user_ci_up_adj + user_ci_lo_adj) / 2,
    user_es_adj
  )

  se <- ifelse(is.na(user_se_adj) & !is.na(user_es_measure_adj) & !is.na(user_ci_lo_adj) & !is.na(user_ci_up_adj),
    (user_ci_up_adj - user_ci_lo_adj) / (2 * qnorm(.975)),
    user_se_adj
  )

  ci_lo <- ifelse(is.na(user_ci_lo_adj) & !is.na(user_es_measure_adj) & !is.na(user_es_adj) & !is.na(user_se_adj),
    user_es_adj - user_se_adj * qnorm(.975),
    user_ci_lo_adj
  )

  ci_up <- ifelse(is.na(user_ci_up_adj) & !is.na(user_es_measure_adj) & !is.na(user_es_adj) & !is.na(user_se_adj),
    user_es_adj + user_se_adj * qnorm(.975),
    user_ci_up_adj
  )

  res <- data.frame(info_used = rep(NA, length(user_es_measure_adj)))

  res$info_used <- "user_input_adj"

  res[, paste(unique(measure))] <- es
  res[, paste(unique(measure), "_se", sep = "")] <- se
  res[, paste(unique(measure), "_ci_lo", sep = "")] <- ci_lo
  res[, paste(unique(measure), "_ci_up", sep = "")] <- ci_up

  attr(res, "measure") <- user_es_measure_adj
  return(res)
}
