# # !!! for alll internal formulas, check whether nn_miss does not prevent the code to be ran
# # important to check the input of the method before calling these functions
#
# "Re: Specific question 4. I've thought of a new strategy. Instead of pivoting on the p-value, we can pivot on the "effective n". In other words, we calculate the n that would result in vz, and then we calculated vr using its formula with this n:
# effective_n = 1 / vz + 3
# vr = (1 - r^2) / (effective_n - 2)"
# 'Re: Specific question 1. That is strange, the formula of the variance of r is sqrt((1-r_raw_data^2) / (n-2)), see it e.g., in wikipedia (https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Standard_error). Indeed, you can check that this formula coincides with the one returned by "cor.test":'
#
# 'Re: Specific questions 2 and 3. I don't know, maybe we could use qnorm but trigger a warning?
# '


#
# One last question about this. In some instances, authors can obtain R without the sample size
# (e.g., when converting OR => R using the pearson/digby approaches).
# What should we do to estimate the 95% CI in these cases?
# We apply the qt() to obtain the 95% CI if users give the n_sample,
# but we derive the 95% from the Z (tanh(z_ci_lo) ; tanh(z_ci_up))
# if we have no information on sample? (or we always use tanh /
# we force users to indicate a sample size?)
#
# ## Specific question 3
# This question 2 led me to the exact same question for SMD. When converting OR + SE => SMD, users can obtain SMD + SE without having access to a sample size. In these cases, we derive the 95% CI using qnorm()? (otherwise, using qt, we present a SMD+SE but without 95% CI)
#

################### OR to RR ###################
.or_to_rr <- function(or, logor_se, or_ci_lo, or_ci_up,
                      n_cases, n_controls, n_exp, n_nexp, baseline_risk, or_to_rr) {
   if (or_to_rr == "grant") {
    logrr_grant <- suppressWarnings(log(or / (1 - baseline_risk + (baseline_risk * or))))
    logrr_ci_lo_grant_CI <- suppressWarnings(log(or_ci_lo / (1 - baseline_risk + (baseline_risk * or_ci_lo))))
    logrr_ci_up_grant_CI <- suppressWarnings(log(or_ci_up / (1 - baseline_risk + (baseline_risk * or_ci_up))))
    logrr_se_CI <- (logrr_ci_up_grant_CI - logrr_ci_lo_grant_CI) / (2 * qnorm(.975))
    res <- cbind(
      logrr = logrr_grant,
      logrr_se = logrr_se_CI,
      logrr_ci_lo = logrr_ci_lo_grant_CI,
      logrr_ci_up = logrr_ci_up_grant_CI
    )

    return(res)
  } else if (or_to_rr == "metaumbrella_cases") {
    contingency_meta_cases <- .estimate_n_from_or_and_n_cases(
      or = or, var = logor_se^2,
      n_cases = n_cases, n_controls = n_controls
    )

    calc_meta_cases <- es_from_2x2(
      n_cases_exp = contingency_meta_cases$n_cases_exp,
      n_controls_exp = contingency_meta_cases$n_controls_exp,
      n_cases_nexp = contingency_meta_cases$n_cases_nexp,
      n_controls_nexp = contingency_meta_cases$n_controls_nexp
    )

    res <- cbind(
      logrr = calc_meta_cases$logrr,
      logrr_se = calc_meta_cases$logrr_se,
      logrr_ci_lo = calc_meta_cases$logrr_ci_lo,
      logrr_ci_up = calc_meta_cases$logrr_ci_up
    )

    return(res)
  } else if (or_to_rr == "metaumbrella_exp") {
    contingency_meta_exp <- .estimate_n_from_or_and_n_exp(or = or, var = logor_se^2, n_exp = n_exp, n_nexp = n_nexp)
    calc_meta_exp <- es_from_2x2(
      n_cases_exp = contingency_meta_exp$n_cases_exp,
      n_controls_exp = contingency_meta_exp$n_controls_exp,
      n_cases_nexp = contingency_meta_exp$n_cases_nexp,
      n_controls_nexp = contingency_meta_exp$n_controls_nexp
    )

    res <- cbind(
      logrr = calc_meta_exp$logrr,
      logrr_se = calc_meta_exp$logrr_se,
      logrr_ci_lo = calc_meta_exp$logrr_ci_lo,
      logrr_ci_up = calc_meta_exp$logrr_ci_up
    )

    return(res)
  } else if (or_to_rr == "transpose") {
    res <- cbind(
      logrr = log(or),
      logrr_se = logor_se,
      logrr_ci_lo = log(or_ci_lo),
      logrr_ci_up = log(or_ci_up)
    )

    return(res)
  } else if (or_to_rr == "dipietrantonj") {
    n_dec <- max(
      nchar(gsub("^.+[.]", "", or)),
      nchar(gsub("^.+[.]", "", or_ci_lo)),
      nchar(gsub("^.+[.]", "", or_ci_up))
    )

    estim <- estimraw::estim_raw(
      es = or, lb = or_ci_lo, ub = or_ci_up,
      m1 = n_exp, m2 = n_nexp, dec = n_dec, measure = "or"
    )

    if (length(estim) != 4) {
      numb = 1
      if (!is.na(baseline_risk)) {
        numb = which.min(
          c(
            abs(estim[[1]]$c[1] / (estim[[1]]$c[1] + estim[[1]]$d[1]) -
                  baseline_risk),
            abs(estim[[2]]$c[1] / (estim[[2]]$c[1] + estim[[2]]$d[1]) -
                  baseline_risk)
            )
        )
      }
      calc_dipie <- es_from_2x2(
        n_cases_exp = estim[[numb]]$a[1],
        n_controls_exp = estim[[numb]]$b[1],
        n_cases_nexp = estim[[numb]]$c[1],
        n_controls_nexp = estim[[numb]]$d[1]
      )

    } else {
      calc_dipie <- es_from_2x2(
        n_cases_exp = estim$a[1],
        n_controls_exp = estim$b[1],
        n_cases_nexp = estim$c[1],
        n_controls_nexp = estim$d[1]
      )
    }
    res <- data.frame(
      logrr = as.numeric(calc_dipie$logrr),
      logrr_se = as.numeric(calc_dipie$logrr_se),
      logrr_ci_lo = as.numeric(calc_dipie$logrr_ci_lo),
      logrr_ci_up = as.numeric(calc_dipie$logrr_ci_up)
    )
    res$logrr = res$logrr
    res$logrr_se = res$logrr_se
    res$logrr_ci_lo = res$logrr_ci_lo
    res$logrr_ci_up = res$logrr_ci_up

    return(res)
  }
}
################### OR to 2x2 ##################
# internal function
.estimate_n_from_or_and_n_cases <- function(or, var, n_cases, n_controls) {
  res <- data.frame(n_cases_exp = NA, n_cases_nexp = NA, n_controls_exp = NA, n_controls_nexp = NA)

  if (!is.na(or) & !is.na(var) & !is.na(n_cases) & !is.na(n_controls)) {
    # Create all possibilites of n
    n_cases_nexp_sim1 <- 0:n_cases
    n_controls_nexp_sim1 <- round(n_controls * (1 - (n_cases - n_cases_nexp_sim1) / (n_cases + (or - 1) * n_cases_nexp_sim1)))
    n_cases_exp_sim1 <- n_cases - n_cases_nexp_sim1
    n_controls_exp_sim1 <- n_controls - n_controls_nexp_sim1
    # sim1: possiblities with strictly positive n
    idx_non_zero <- which(
      n_cases_nexp_sim1 > 0 &
        n_controls_nexp_sim1 > 0 &
        n_cases_exp_sim1 > 0 &
        n_controls_exp_sim1 > 0
    )
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # sim2: possiblities with positive n and at least one zero (Add 0.5 to the possiblities with any 0)
    n_cases_nexp_sim2 <- 0:n_cases
    n_controls_nexp_sim2 <- round((n_controls + 0.5) - (n_controls + 1) * (n_cases - n_cases_nexp_sim2 + 0.5) / ((n_cases_nexp_sim2 + 0.5) * or + n_cases - n_cases_nexp_sim2 + 0.5))
    n_cases_exp_sim2 <- n_cases - n_cases_nexp_sim2
    n_controls_exp_sim2 <- n_controls - n_controls_nexp_sim2
    # (n_cases_exp_sim2 + 0.5) / (n_cases_nexp_sim2 + 0.5) / (n_controls_exp_sim2 + 0.5) * (n_controls_nexp_sim2 + 0.5)
    # select the ones with some 0 but non-negative
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    # join both previous vectors
    n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)
    n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)
    n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)

    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0

    var_sim <- ifelse(some_zero,
      1 / ((n_cases + 1) - (n_cases_nexp_sim + 0.5)) + 1 / ((n_controls + 1) - (n_controls_nexp_sim + 0.5)) + 1 / (n_cases_nexp_sim + 0.5) + 1 / (n_controls_nexp_sim + 0.5),
      1 / (n_cases - n_cases_nexp_sim) + 1 / (n_controls - n_controls_nexp_sim) + 1 / n_cases_nexp_sim + 1 / n_controls_nexp_sim
    )

    # var_sim2 = 1 / ((n_cases+1) - (n_cases_nexp_sim+0.5)) + 1 / ((n_controls+1) - (n_controls_nexp_sim+0.5)) + 1 / (n_cases_nexp_sim+0.5) + 1 / (n_controls_nexp_sim+0.5)

    best <- order((var_sim - var)^2)[1]

    res$n_cases_nexp <- n_cases_nexp_sim[best]
    res$n_controls_nexp <- n_controls_nexp_sim[best]
    res$n_cases_exp <- n_cases - res$n_cases_nexp
    res$n_controls_exp <- n_controls - res$n_controls_nexp
  }

  return(res)
}

#' Estimate the n, using the variance, the number of exposed and non-exposed subjects
#'
#' @param or OR
#' @param var variance
#' @param n_exp number of exposed participants
#' @param n_nexp number of non exposed participants
#'
#' @noRd
.estimate_n_from_or_and_n_exp <- function(or, var, n_exp, n_nexp) {
  res <- data.frame(n_cases_exp = NA, n_cases_nexp = NA, n_controls_exp = NA, n_controls_nexp = NA)

  if (!is.na(or) & !is.na(var) & !is.na(n_exp) & !is.na(n_nexp)) {
    # first: uncorrected values with 0
    n_controls_exp_sim1 <- 0:n_exp
    n_controls_nexp_sim1 <- round(n_nexp / (1 + (n_exp - n_controls_exp_sim1) / (or * n_controls_exp_sim1)))
    n_cases_exp_sim1 <- n_exp - n_controls_exp_sim1
    n_cases_nexp_sim1 <- n_nexp - n_controls_nexp_sim1
    # we take the ones without 0 and non-negative
    idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # correcting by 0.5
    n_controls_exp_sim2 <- 0:n_exp
    n_controls_nexp_sim2 <- round((n_nexp + 0.5) - ((n_nexp + 1) * (n_exp - n_controls_exp_sim2 + 0.5)) / ((n_controls_exp_sim2 + 0.5) * or + n_exp - n_controls_exp_sim2 + 0.5))
    n_cases_exp_sim2 <- n_exp - n_controls_exp_sim2
    n_cases_nexp_sim2 <- n_nexp - n_controls_nexp_sim2

    # SELECT THE ONES THAT HAS SOME 0 BUT NO NEGATIVE ONES
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)
    n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)
    n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)


    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
    var_sim <- ifelse(some_zero,
      1 / ((n_exp + 1) - (n_controls_exp_sim + 0.5) + 1 / (n_controls_exp_sim + 0.5) + 1 / ((n_nexp + 1) - (n_controls_nexp_sim + 0.5)) + 1 / (n_controls_nexp_sim + 0.5)),
      1 / (n_exp - n_controls_exp_sim) + 1 / n_controls_exp_sim + 1 / (n_nexp - n_controls_nexp_sim) + 1 / n_controls_nexp_sim
    )

    # var_sim = 1 / (n_exp - n_controls_exp_sim) + 1 / n_controls_exp_sim + 1 / (n_nexp - n_controls_nexp_sim) + 1 / n_controls_nexp_sim
    best <- order((var_sim - var)^2)[1]
    res$n_controls_exp <- n_controls_exp_sim[best]
    res$n_controls_nexp <- n_controls_nexp_sim[best]
    res$n_cases_exp <- n_exp - res$n_controls_exp
    res$n_cases_nexp <- n_nexp - res$n_controls_nexp
  }
  return(res)
}

################### SE of OR ###################
.se_from_or <- function(x) {
  or <- as.numeric(x[1])
  n_cases <- as.numeric(x[2])
  n_controls <- as.numeric(x[3])

  res <- data.frame(value = NA, var = NA, se = NA)

  if (!is.na(or) & !is.na(n_cases) & !is.na(n_controls)) {
    cases_exp <- 1:(n_cases - 1)
    cases_nexp <- n_cases - cases_exp
    controls_exp <- round(n_controls / (1 + cases_nexp * or / cases_exp))
    controls_exp[which(controls_exp < 1 | controls_exp > n_controls - 1)] <- NA
    controls_nexp <- n_controls - controls_exp
    v_or_mean <- mean(1 / cases_exp + 1 / cases_nexp + 1 / controls_exp + 1 / controls_nexp, na.rm = TRUE)

    res$value <- or
    res$var <- v_or_mean
    res$se <- sqrt(v_or_mean)
  }

  return(res)
}

################### RR to OR ###################
.rr_to_or <- function(rr, logrr_se, rr_ci_lo, rr_ci_up,
                      n_cases, n_controls, n_exp, n_nexp,
                      baseline_risk, rr_to_or) {
  if (rr_to_or == "grant") {
    logor_grant <- suppressWarnings(log(rr * (1 - baseline_risk) / (1 - rr * baseline_risk)))
    logor_ci_lo_grant <- suppressWarnings(log(rr_ci_lo * (1 - baseline_risk) / (1 - rr_ci_lo * baseline_risk)))
    logor_ci_up_grant <- suppressWarnings(log(rr_ci_up * (1 - baseline_risk) / (1 - rr_ci_up * baseline_risk)))
    logor_se <- (logor_ci_up_grant - logor_ci_lo_grant) / (2 * qnorm(.975))
    res <- cbind(
      logor = logor_grant,
      logor_se = logor_se,
      logor_ci_lo = logor_ci_lo_grant,
      logor_ci_up = logor_ci_up_grant
    )
    return(res)
  } else if (rr_to_or == "metaumbrella") {
    raw_res <- .metaumbrella_rr_se_to_or(rr = rr, logrr_se = logrr_se, n_cases = n_cases, n_controls = n_controls)
    res <- cbind(
      logor = log(raw_res$value),
      logor_se = raw_res$se,
      logor_ci_lo = log(raw_res$value) - qnorm(.975) * raw_res$se,
      logor_ci_up = log(raw_res$value) + qnorm(.975) * raw_res$se
    )
    return(res)
  } else if (rr_to_or == "transpose") {
    res <- cbind(
      logor = log(rr),
      logor_se = logrr_se,
      logor_ci_lo = log(rr_ci_lo),
      logor_ci_up = log(rr_ci_up)
    )
    return(res)
  } else if (rr_to_or == "dipietrantonj") {
    n_dec <- max(
      nchar(gsub("^.+[.]", "", rr)),
      nchar(gsub("^.+[.]", "", rr_ci_lo)),
      nchar(gsub("^.+[.]", "", rr_ci_up))
    )

    estim <- estimraw::estim_raw(
      es = rr, lb = rr_ci_lo, ub = rr_ci_up,
      m1 = n_exp, m2 = n_nexp, dec = n_dec, measure = "rr"
    )
    if (length(estim) != 4) {
      numb = 1
      if (!is.na(baseline_risk)) {
        numb = which.min(
          c(
            abs(estim[[1]]$c[1] / (estim[[1]]$c[1] + estim[[1]]$d[1]) - baseline_risk),
            abs(estim[[2]]$c[1] / (estim[[2]]$c[1] + estim[[2]]$d[1]) - baseline_risk)
          )
        )
      }
      calc_dipie <- es_from_2x2(
        n_cases_exp = estim[[numb]]$a[1],
        n_controls_exp = estim[[numb]]$b[1],
        n_cases_nexp = estim[[numb]]$c[1],
        n_controls_nexp = estim[[numb]]$d[1]
      )
    } else {
      calc_dipie <- es_from_2x2(
        n_cases_exp = estim$a[1],
        n_controls_exp = estim$b[1],
        n_cases_nexp = estim$c[1],
        n_controls_nexp = estim$d[1]
      )
    }
    res <- cbind(
      logor = calc_dipie$logor,
      logor_se = calc_dipie$logor_se,
      logor_ci_lo = calc_dipie$logor_ci_lo,
      logor_ci_up = calc_dipie$logor_ci_up
    )

    return(res)
  } else {
    stop(paste0("'", rr_to_or, "' not in tolerated values for the 'rr_to_or' argument. Possible inputs are: 'metaumbrella', 'transpose', 'grant_2x2', 'grant_CI'"))
  }
}

.metaumbrella_rr_se_to_or <- function(rr, logrr, logrr_se, n_cases, n_controls) {
  es <- data.frame(value = NA, se = NA)

  if (!is.na(rr) & !is.na(logrr_se) & !is.na(n_cases) & !is.na(n_controls)) {
    # uncorrected
    n_cases_nexp_sim1 <- 0:n_cases
    n_controls_nexp_sim1 <- round(n_cases_nexp_sim1 * ((rr * (n_cases + n_controls)) / (n_cases + (rr - 1) * n_cases_nexp_sim1) - 1))
    n_cases_exp_sim1 <- n_cases - n_cases_nexp_sim1
    n_controls_exp_sim1 <- n_controls - n_controls_nexp_sim1

    # we take only positives (no-zero)
    idx_non_zero <- which(n_cases_nexp_sim1 > 0 & n_controls_nexp_sim1 > 0 & n_cases_exp_sim1 > 0 & n_controls_exp_sim1 > 0) # Posem ">" i no "!=" per treure els negatius!
    n_cases_nexp_sim1 <- n_cases_nexp_sim1[idx_non_zero]
    n_controls_nexp_sim1 <- n_controls_nexp_sim1[idx_non_zero]
    n_cases_exp_sim1 <- n_cases_exp_sim1[idx_non_zero]
    n_controls_exp_sim1 <- n_controls_exp_sim1[idx_non_zero]

    # corregint 0.5
    n_cases_nexp_sim2 <- 0:n_cases
    n_controls_nexp_sim2 <- ((n_cases + n_controls - n_cases_nexp_sim2 + 1) - (n_cases + n_controls + 2) * (n_cases - n_cases_nexp_sim2 + 0.5) / ((n_cases_nexp_sim2 + 0.5) * rr + n_cases - n_cases_nexp_sim2 + 0.5))
    n_cases_exp_sim2 <- n_cases - n_cases_nexp_sim2
    n_controls_exp_sim2 <- n_controls - n_controls_nexp_sim2

    # we take the ones with some 0 but non negative
    idx_some_zero <- which(
      (n_cases_nexp_sim2 == 0 | n_controls_nexp_sim2 == 0 | n_cases_exp_sim2 == 0 | n_controls_exp_sim2 == 0) &
        (n_cases_nexp_sim2 >= 0 & n_controls_nexp_sim2 >= 0 & n_cases_exp_sim2 >= 0 & n_controls_exp_sim2 >= 0)
    )
    n_cases_nexp_sim2 <- n_cases_nexp_sim2[idx_some_zero]
    n_controls_nexp_sim2 <- n_controls_nexp_sim2[idx_some_zero]
    n_cases_exp_sim2 <- n_cases_exp_sim2[idx_some_zero]
    n_controls_exp_sim2 <- n_controls_exp_sim2[idx_some_zero]

    #
    n_controls_exp_sim <- append(n_controls_exp_sim1, n_controls_exp_sim2)
    n_controls_nexp_sim <- append(n_controls_nexp_sim1, n_controls_nexp_sim2)
    n_cases_nexp_sim <- append(n_cases_nexp_sim1, n_cases_nexp_sim2)
    n_cases_exp_sim <- append(n_cases_exp_sim1, n_cases_exp_sim2)

    some_zero <- n_cases_exp_sim == 0 | n_controls_exp_sim == 0 | n_cases_nexp_sim == 0 | n_controls_nexp_sim == 0
    var_sim <- ifelse(some_zero,
      1 / ((n_cases + 1) - (n_cases_nexp_sim + 0.5)) + 1 / ((n_cases + 1) + (n_controls + 1) - ((n_cases_nexp_sim + 0.5) + (n_controls_nexp_sim + 0.5))) +
        1 / (n_cases_nexp_sim + 0.5) + 1 / ((n_cases_nexp_sim + 0.5) + (n_controls_nexp_sim + 0.5)),
      1 / (n_cases - n_cases_nexp_sim) + 1 / (n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim)) +
        1 / n_cases_nexp_sim + 1 / (n_cases_nexp_sim + n_controls_nexp_sim)
    )

    var_sim <- 1 / (n_cases - n_cases_nexp_sim) + 1 / (n_cases + n_controls - (n_cases_nexp_sim + n_controls_nexp_sim)) +
      1 / n_cases_nexp_sim + 1 / (n_cases_nexp_sim + n_controls_nexp_sim)

    best <- order((var_sim - logrr_se^2)^2)[1]
    n_cases_nexp <- n_cases_nexp_sim[best]
    n_controls_nexp <- n_controls_nexp_sim[best]
    n_cases_exp <- n_cases - n_cases_nexp
    n_controls_exp <- n_controls - n_controls_nexp

    # es$n_cases_nexp = n_cases_nexp
    # es$n_controls_nexp = n_controls_nexp
    # es$n_cases_exp = n_cases_exp
    # es$n_controls_exp = n_controls_exp

    cont_table <- es_from_2x2(
      n_cases_exp = n_cases_exp, n_controls_exp = n_controls_exp,
      n_cases_nexp = n_cases_nexp, n_controls_nexp = n_controls_nexp
    )

    es$value <- exp(cont_table$logor)
    es$se <- cont_table$logor_se
  }
  return(es)
}



################# 2x2 to R/Z ###################
.contigency_to_cor <- function(n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp,
                               table_2x2_to_cor, reverse_2x2) {

  if (table_2x2_to_cor == "lipsey") {
    # TO DO
    # log_or <- log((n_cases_exp * n_controls_nexp) / (n_cases_nexp * n_controls_exp))
    # v.log_or <- 1 / n_cases_exp + 1 / n_cases_nexp + 1 / n_controls_exp + 1 / n_controls_nexp
    #
    # r <- (n_cases_exp * n_controls_nexp - n_controls_exp * n_cases_nexp) /
    #   sqrt((n_cases_exp + n_controls_exp) * (n_cases_nexp + n_controls_nexp) *
    #     (n_cases_exp + n_cases_nexp) * (n_controls_exp + n_cases_nexp))
    # r_lipsey <- ifelse(reverse_2x2, -r, r)
    # z_lipsey <- atanh(r_lipsey)
    # vz_lipsey <- v.log_or * (z_lipsey^2) / (log_or^2)
    #
    # z_lo_lipsey <- z_lipsey - qnorm(.975) * sqrt(vz_lipsey)
    # z_up_lipsey <- z_lipsey + qnorm(.975) * sqrt(vz_lipsey)
    # r_lo_lipsey <- tanh(z_lo_lipsey)
    # r_up_lipsey <- tanh(z_up_lipsey)
    #
    # effective_n = 1/vz_lipsey + 3
    # vr_lipsey = (1 - r_lipsey^2)^2 / (effective_n - 1)
    #
    # res <- cbind(
    #   r_lipsey, vr_lipsey, r_lo_lipsey, r_up_lipsey,
    #   z_lipsey, vz_lipsey, z_lo_lipsey, z_up_lipsey
    # )
    # return(res)
  } else if (table_2x2_to_cor == "tetrachoric") {
    res <- .tet_r(as.numeric(n_cases_exp),
                  as.numeric(n_controls_exp),
                  as.numeric(n_cases_nexp),
                  as.numeric(n_controls_nexp))
    res[res == "calculation failure"] <- NA

    # res[] <- lapply(res, function(x) as.numeric(as.character(x)))
    res[1] <- ifelse(reverse_2x2, -res[1], res[1])
    res[3] <- ifelse(reverse_2x2, -res[3], res[3])
    res[4] <- ifelse(reverse_2x2, -res[4], res[4])
    res[5] <- ifelse(reverse_2x2, -res[5], res[5])
    res[7] <- ifelse(reverse_2x2, -res[7], res[7])
    res[8] <- ifelse(reverse_2x2, -res[8], res[8])

    return(res)
  } else {
    stop(paste0("'", table_2x2_to_cor, "' not in tolerated values for the 'table_2x2_to_cor' argument. Possible inputs are: 'tetrachoric', 'cooper_delta', 'cooper_std' or 'lipsey'"))
  }
}

.tet_r <- function(n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp) {
  tryCatch(
    expr = {
      tet <- metafor::escalc(
        ai = n_cases_exp, bi = n_controls_exp,
        ci = n_cases_nexp, di = n_controls_nexp,
        measure = "RTET"
      )
      n_sample = n_cases_exp + n_controls_exp + n_cases_nexp + n_controls_nexp
      r <- tet$yi
      vr <- tet$vi
      z <- atanh(r)
      vz <- vr / ((1 - r^2)^2)

      z_lo <- z - qnorm(.975) * sqrt(vz)
      z_up <- z + qnorm(.975) * sqrt(vz)
      r_lo <- r - qnorm(.975) * sqrt(vr)
      r_up <- r + qnorm(.975) * sqrt(vr)

      dat <- cbind(r, vr, r_lo, r_up, z, vz, z_lo, z_up)
      return(dat)
    },
    error = function(e) {
      dat <- cbind(NA, NA, NA, NA, NA, NA, NA, NA)
      return(dat)
    },
    warning = function(w) {
      dat <- cbind(NA, NA, NA, NA, NA, NA, NA, NA)
      return(dat)
    },
    finally = {
    }
  )
}

################# OR to R/Z ###################
.or_to_cor <- function(or, logor_se,
                       n_cases,
                       n_exp,
                       small_margin_prop,
                       n_sample,
                       or_to_cor) {
  if (or_to_cor %in% c("bonett", "pearson")) {
    if (or_to_cor == "bonett") {
      c <- (1 - abs(n_exp/n_sample - n_cases/n_sample) / 5 - (1 / 2 - small_margin_prop)^2) / 2
    } else {
      c <- 1 / 2
    }

    r <- cos(pi / (1 + or^c))
    r_se <- logor_se * (pi * c * or^c) * sin(pi / (1 + or^c)) / (1 + or^c)^2

    or_ci_lo <- exp(log(or) - qnorm(.975) * logor_se)
    or_ci_up <- exp(log(or) + qnorm(.975) * logor_se)
    r_lo <- cos(pi / (1 + or_ci_lo^c))
    r_up <- cos(pi / (1 + or_ci_up^c))

    z <- atanh(r)
    z_se <- sqrt(r_se^2 / ((1 - r^2)^2)) # delta method
    z_lo <- atanh(r_lo)
    z_up <- atanh(r_up)

    res <- cbind(r, r_se, r_lo, r_up, z, z_se, z_lo, z_up)
    return(res)
  } else if (or_to_cor == "digby") {
    c <- 3 / 4

    r <- (or^c - 1) / (or^c + 1)
    r_se <- sqrt((c^2 / 4) * (1 - r^2)^2 * logor_se^2)

    z <- atanh(r)
    z_se <- sqrt(r_se^2 / ((1 - r^2)^2)) # delta method
    z_lo <- z - qnorm(.975) * sqrt(c^2 / 4 * logor_se^2)
    z_up <- z + qnorm(.975) * sqrt(c^2 / 4 * logor_se^2)

    r_lo <- tanh(z_lo)
    r_up <- tanh(z_up)

    res <- cbind(r, r_se, r_lo, r_up, z, z_se, z_lo, z_up)
    return(res)
  }
}

################# SMD to R/Z ###################
.smd_to_cor <- function(d, vd, n_exp, n_nexp, smd_to_cor, n_cov_ancova) {
  if (smd_to_cor == "viechtbauer") {
    df <- n_exp + n_nexp - 2 - n_cov_ancova
    h <- df / n_exp + df / n_nexp
    p <- n_exp / (n_exp + n_nexp)
    q <- n_nexp / (n_exp + n_nexp)
    r_pb <- d / sqrt(d^2 + h)

    f <- dnorm(qnorm(p, lower.tail = FALSE))
    r_viechtbauer <- sqrt(p * q) / f * r_pb
    r_trunc = ifelse(r_viechtbauer > 1, 1, ifelse(r_viechtbauer < -1, -1, r_viechtbauer))
    vr_viechtbauer <- 1 / (n_exp + n_nexp - 1) *
      (p * q / f^2 - (3 / 2 + (1 - p * qnorm(p, lower.tail = FALSE) / f) *
                              (1 + q * qnorm(p, lower.tail = FALSE) / f)) *
         r_trunc^2 + r_trunc^4)


    # ========= z ========= #
    fzp <- dnorm(qnorm(p))
    a_viechtbauer <- sqrt(fzp) / (p * (1 - p))^(1 / 4)
    z_viechtbauer <- (a_viechtbauer / 2) * log((1 + a_viechtbauer * r_trunc) /
                                               (1 - a_viechtbauer * r_trunc))
    vz_viechtbauer <- 1 / (n_exp + n_nexp - 1)
    # ========= 95% CI ===== #


    z_lo_viechtbauer <- z_viechtbauer - qnorm(.975) * sqrt(vz_viechtbauer)
    z_up_viechtbauer <- z_viechtbauer + qnorm(.975) * sqrt(vz_viechtbauer)
    r_lo_viechtbauer <- (1 / a_viechtbauer) * ((exp(2 * z_lo_viechtbauer / a_viechtbauer) - 1) / (exp(2 * z_lo_viechtbauer / a_viechtbauer) + 1))
    r_up_viechtbauer <- (1 / a_viechtbauer) * ((exp(2 * z_up_viechtbauer / a_viechtbauer) - 1) / (exp(2 * z_up_viechtbauer / a_viechtbauer) + 1))

    res <- cbind(
      r_viechtbauer, vr_viechtbauer, r_lo_viechtbauer, r_up_viechtbauer,
      z_viechtbauer, vz_viechtbauer, z_lo_viechtbauer, z_up_viechtbauer
    )

    return(res)
  } else if (smd_to_cor == "lipsey_cooper") {
    a <- ((n_exp + n_nexp)^2) / (n_exp * n_nexp)
    p <- n_exp / (n_exp + n_nexp)
    r_lipsey <- d / sqrt(d^2 + 1 / (p * (1 - p)))
    vr_lipsey <- a^2 * vd / ((d^2 + a)^3)
    z_lipsey <- atanh(r_lipsey)
    vz_lipsey <- vd / (vd + 1 / (p * (1 - p)))
    r_lo_lipsey <- r_lipsey - qt(.975, df = n_exp + n_nexp - 2) * sqrt(vr_lipsey)
    r_up_lipsey <- r_lipsey + qt(.975, df = n_exp + n_nexp - 2) * sqrt(vr_lipsey)
    z_lo_lipsey <- z_lipsey - qnorm(.975) * sqrt(vz_lipsey)
    z_up_lipsey <- z_lipsey + qnorm(.975) * sqrt(vz_lipsey)

    res <- cbind(
      r_lipsey, vr_lipsey, r_lo_lipsey, r_up_lipsey,
      z_lipsey, vz_lipsey, z_lo_lipsey, z_up_lipsey
    )

    return(res)
  }
}

################# PHI to R/Z ###################
.phi_to_cor <- function(phi, n_sample,
                        n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp,
                        phi_to_cor, reverse_phi) {
  if (phi_to_cor == "lipsey") {
    r <- ifelse(reverse_phi, -phi, phi)
    z <- atanh(r)
    vz <- z^2 / (r^2 * n_sample)
    z_lo <- z - qnorm(.975) * sqrt(vz)
    z_up <- z + qnorm(.975) * sqrt(vz)
    effective_n = 1 / vz + 3
    vr = (1 - r^2)^2 / (effective_n - 1)
    # t <- qt(pnorm(z / sqrt(vz)), n_sample - 2)
    # vr <- (r / t)^2
    # r_lo = r - qnorm(.975)*sqrt(vr)
    # r_up = r + qnorm(.975)*sqrt(vr)
    r_lo <- tanh(z_lo)
    r_up <- tanh(z_up)
    res <- cbind(r, vr, r_lo, r_up, z, vz, z_lo, z_up)
    return(res)
  } else if (phi_to_cor == "tetrachoric") {
    # if (!requireNamespace("mvtnorm", quietly = TRUE) & !requireNamespace("metafor", quietly = TRUE)) {
    #   stop("Please install the 'mvtnorm' and 'metafor' packages to compute tetrachoric correlation.")
    # } else if (!requireNamespace("metafor", quietly = TRUE)) {
    #   stop("Please install the 'metafor' package to compute tetrachoric correlation.")
    # } else if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    #   stop("Please install the 'mvtnorm' package to compute tetrachoric correlation.")
    # }

    res <- .tet_r(
      as.numeric(n_cases_exp), as.numeric(n_controls_exp),
      as.numeric(n_cases_nexp), as.numeric(n_controls_nexp)
    )
    res[res == "calculation failure"] <- NA

    res[1] <- ifelse(reverse_phi, -res[1], res[1])
    res[3] <- ifelse(reverse_phi, -res[3], res[3])
    res[4] <- ifelse(reverse_phi, -res[4], res[4])
    res[5] <- ifelse(reverse_phi, -res[5], res[5])
    res[7] <- ifelse(reverse_phi, -res[7], res[7])
    res[8] <- ifelse(reverse_phi, -res[8], res[8])

    return(res)
  } else {
    stop(paste0("'", phi_to_cor, "' not in tolerated values for the 'phi_to_cor' argument. Possible inputs are: 'tetrachoric', or 'lipsey'"))
  }
}

################# CHI to R/Z ###################
.chi_to_cor <- function(chisq, n_sample,
                        n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp,
                        chisq_to_cor, reverse_chisq) {
  if (chisq_to_cor == "lipsey") {
    r <- sqrt(chisq / n_sample)
    r <- ifelse(reverse_chisq, -r, r)
    z <- atanh(r)
    vz <- z^2 / (chisq)
    z_lo <- z - qnorm(.975) * sqrt(vz)
    z_up <- z + qnorm(.975) * sqrt(vz)
    # t <- qt(pnorm(z / sqrt(vz)), n_sample - 2)
    # vr <- (r / t)^2
    effective_n = 1 / vz + 3
    vr = (1 - r^2)^2 / (effective_n - 1)
    r_lo <- tanh(z_lo)
    r_up <- tanh(z_up)
    # r_lo = r - qt(.975, n_sample - 2)*sqrt(vr)
    # r_up = r + qt(.975, n_sample - 2)*sqrt(vr)
    # r_lo <- tanh(z_lo)
    # r_up <- tanh(z_up)
    res <- cbind(r, vr, r_lo, r_up, z, vz, z_lo, z_up)
    return(res)
  } else if (chisq_to_cor == "tetrachoric") {
    # if (!requireNamespace("mvtnorm", quietly = TRUE) & !requireNamespace("metafor", quietly = TRUE)) {
    #   stop("Please install the 'mvtnorm' and 'metafor' packages to compute tetrachoric correlation.")
    # } else if (!requireNamespace("metafor", quietly = TRUE)) {
    #   stop("Please install the 'metafor' package to compute tetrachoric correlation.")
    # } else if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    #   stop("Please install the 'mvtnorm' package to compute tetrachoric correlation.")
    # }

    res <- .tet_r(
      as.numeric(n_cases_exp), as.numeric(n_controls_exp),
      as.numeric(n_cases_nexp), as.numeric(n_controls_nexp)
    )
    res[res == "calculation failure"] <- NA

    res[1] <- ifelse(reverse_chisq, -res[1], res[1])
    res[3] <- ifelse(reverse_chisq, -res[3], res[3])
    res[4] <- ifelse(reverse_chisq, -res[4], res[4])
    res[5] <- ifelse(reverse_chisq, -res[5], res[5])
    res[7] <- ifelse(reverse_chisq, -res[7], res[7])
    res[8] <- ifelse(reverse_chisq, -res[8], res[8])

    return(res)
  } else {
    stop(paste0("'", chisq_to_cor, "' not in tolerated values for the 'chisq_to_cor' argument. Possible inputs are: 'tetrachoric', or 'lipsey'"))
  }
}


################# PRE POST to SMD ##############
.pre_post_to_smd <- function(mean_pre_exp, mean_pre_sd_exp,
                             mean_exp, mean_sd_exp,
                             mean_pre_nexp, mean_pre_sd_nexp,
                             mean_nexp, mean_sd_nexp,
                             n_exp, n_nexp,
                             r_pre_post_exp, r_pre_post_nexp,
                             pre_post_to_smd) {
  if (pre_post_to_smd == "bonett") {
    J_exp <- .d_j(n_exp - 1)
    var_exp <- mean_pre_sd_exp^2 + mean_sd_exp^2 -
      2 * r_pre_post_exp * mean_pre_sd_exp * mean_sd_exp
    d_exp <- (mean_pre_exp - mean_exp) / mean_pre_sd_exp
    g_exp <- d_exp * J_exp
    var_g_exp <- var_exp / (mean_pre_sd_exp^2 * (n_exp - 1)) + g_exp^2 / (2 * (n_exp - 1))
    var_d_exp <-  var_g_exp / (J_exp^2)

    J_nexp <- .d_j(n_nexp - 1)
    var_nexp <- mean_pre_sd_nexp^2 + mean_sd_nexp^2 - 2 * r_pre_post_nexp * mean_pre_sd_nexp * mean_sd_nexp
    d_nexp <- (mean_pre_nexp - mean_nexp) / mean_pre_sd_nexp
    g_nexp <- d_nexp * J_nexp
    var_g_nexp <- var_nexp / (mean_pre_sd_nexp^2 * (n_nexp - 1)) + g_nexp^2 / (2 * (n_nexp - 1))
    var_d_nexp <-  var_g_nexp / (J_nexp^2)

    d_bonett <- d_exp - d_nexp
    g_bonett <- g_exp - g_nexp
    vd_bonett <- var_d_exp + var_d_nexp
    vg_bonett <- var_g_exp + var_g_nexp
    d_lo_bonett <- d_bonett - sqrt(vd_bonett) * qt(.975, n_exp + n_nexp - 2)
    d_up_bonett <- d_bonett + sqrt(vd_bonett) * qt(.975, n_exp + n_nexp - 2)
    g_lo_bonett <- g_bonett - sqrt(vg_bonett) * qt(.975, n_exp + n_nexp - 2)
    g_up_bonett <- g_bonett + sqrt(vg_bonett) * qt(.975, n_exp + n_nexp - 2)

    res <- cbind(
      d_bonett, vd_bonett, d_lo_bonett, d_up_bonett,
      g_bonett, vg_bonett, g_lo_bonett, g_up_bonett
    )

    return(res)
  } else if (pre_post_to_smd == "cooper") {
    J_exp <- .d_j(n_exp - 1)
    J_nexp <- .d_j(n_nexp - 1)

    sd_diff_exp <- sqrt(mean_pre_sd_exp^2 + mean_sd_exp^2 -
      (2 * r_pre_post_exp * mean_pre_sd_exp * mean_sd_exp))

    sd_diff_nexp <- sqrt(mean_pre_sd_nexp^2 + mean_sd_nexp^2 -
      (2 * r_pre_post_nexp * mean_pre_sd_nexp * mean_sd_nexp))

    d_exp <- (mean_pre_exp - mean_exp) / sd_diff_exp * sqrt(2 * (1 - r_pre_post_exp))
    g_exp <- d_exp * J_exp

    d_nexp <- (mean_pre_nexp - mean_nexp) / sd_diff_nexp * sqrt(2 * (1 - r_pre_post_nexp))
    g_nexp <- d_nexp * J_nexp

    # COR  vi[i] <- (2*(1-ri[i])/ni[i] + di[i]^2 / (2*ni[i]))
    # WRONG vi[i] <- 2*(1-ri[i]) * (1/ni[i] + di[i]^2 / (2*ni[i]))
    # corrected formula by W Viechtbauer
    d_var_exp <- 2 * (1 - r_pre_post_exp)/n_exp + d_exp^2 / (2*n_exp)
    d_var_nexp <- 2 * (1 - r_pre_post_nexp)/n_nexp + d_nexp^2 / (2*n_nexp)

    var_g_exp <- J_exp^2 * d_var_exp
    var_g_nexp <- J_nexp^2 * d_var_nexp

    d_cooper <- d_exp - d_nexp
    g_cooper <- g_exp - g_nexp
    vd_cooper <- d_var_exp + d_var_nexp
    vg_cooper <- var_g_exp + var_g_nexp
    d_lo_cooper <- d_cooper - sqrt(vd_cooper) * qt(.975, n_exp + n_nexp - 2)
    d_up_cooper <- d_cooper + sqrt(vd_cooper) * qt(.975, n_exp + n_nexp - 2)
    g_lo_cooper <- g_cooper - sqrt(vg_cooper) * qt(.975, n_exp + n_nexp - 2)
    g_up_cooper <- g_cooper + sqrt(vg_cooper) * qt(.975, n_exp + n_nexp - 2)

    res <- cbind(
      d_cooper, vd_cooper, d_lo_cooper, d_up_cooper,
      g_cooper, vg_cooper, g_lo_cooper, g_up_cooper
    )
    return(res)
  }
}

################# R/Z to SMD ###################
.cor_to_smd <- function(r, r_se,
                        unit_increase_iv, sd_iv, unit_type,
                        n_sample, cor_to_smd) {
  if (cor_to_smd == "mathur") {
    increase <- ifelse(unit_type == "sd",
                       unit_increase_iv * sd_iv,
                       unit_increase_iv)

    d <- r * increase / (sd_iv * sqrt((1 - r^2)))
    d_se <- abs(d) * sqrt(1 / (r^2 * (n_sample - 3)) + 1 / (2 * (n_sample - 1)))
    res <- cbind(d, d_se)
    return(res)
  } else if (cor_to_smd == "viechtbauer") {
    res_g <- metafor::conv.delta(
      yi = r, vi = r_se^2, transf = metafor::transf.rtod, var.names = c("g", "g_var")
    )
    J <- .d_j(n_sample - 2)
    res_g$d <- res_g$g / J
    res_g$d_se <- sqrt(res_g$g_var / J^2)
    res <- cbind(res_g$d, res_g$d_se)

    return(res)
  } else if (cor_to_smd == "cooper") {
    d <- 2 * r / sqrt(1 - r^2)
    d_se <- sqrt(4 * sqrt(r_se) / ((1 - r^2)^3))
    res <- cbind(d, d_se)

    return(res)
  }
}
