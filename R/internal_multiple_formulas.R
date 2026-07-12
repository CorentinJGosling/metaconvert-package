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
      1 / ((n_exp + 1) - (n_controls_exp_sim + 0.5)) + 1 / (n_controls_exp_sim + 0.5) + 1 / ((n_nexp + 1) - (n_controls_nexp_sim + 0.5)) + 1 / (n_controls_nexp_sim + 0.5),
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
.contingency_to_cor <- function(n_cases_exp, n_controls_exp, n_cases_nexp, n_controls_nexp,
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


################# PRE POST helpers ##############
#' Guard a standardizing SD before it is used as a denominator
#'
#' Returns NA for any standardizer that is zero, negative or non-finite, so the
#' resulting SMD is NA rather than a silent Inf/NaN. Scoped deliberately to the
#' DENOMINATOR: the mean-change wrappers legitimately pass mean_pre_sd = 0 (the
#' pre slot is zeroed by construction), so guarding raw inputs would break them,
#' whereas a zero *standardizer* is always degenerate. Reachable when a study
#' reports SD = 0, or when r_pre_post = 1 with equal pre/post SDs makes the
#' change SD collapse to 0.
#'
#' @noRd
.guard_standardizer <- function(sd_value) {
  ifelse(is.finite(sd_value) & sd_value > 0, sd_value, NA_real_)
}

#' Guard a pre-post correlation
#'
#' |r| >= 1 is not a valid correlation and drives sd_change to 0 (or a negative
#' radicand) in every pre-post formula. Strictly at the bounds so that r = 0.99
#' still computes.
#'
#' @noRd
.guard_r_pre_post <- function(r) {
  ifelse(is.finite(r) & abs(r) < 1, r, NA_real_)
}

################# PRE POST to SMD ##############
.pre_post_to_smd <- function(mean_pre_exp, mean_pre_sd_exp,
                             mean_exp, mean_sd_exp,
                             mean_pre_nexp, mean_pre_sd_nexp,
                             mean_nexp, mean_sd_nexp,
                             n_exp, n_nexp,
                             r_pre_post_exp, r_pre_post_nexp,
                             pre_post_to_smd,
                             pool_sd = FALSE) {
  # pool_sd = TRUE: standardizing SD pooled across groups (Morris 2008, eq. 9/13)
  if (pool_sd) {
    return(.pooled_pre_post_to_smd(
      mean_pre_exp = mean_pre_exp, mean_pre_sd_exp = mean_pre_sd_exp,
      mean_exp = mean_exp, mean_sd_exp = mean_sd_exp,
      mean_pre_nexp = mean_pre_nexp, mean_pre_sd_nexp = mean_pre_sd_nexp,
      mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp,
      n_exp = n_exp, n_nexp = n_nexp,
      r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,
      pre_post_to_smd = pre_post_to_smd
    ))
  }

  res_exp <- .single_group_pre_post_to_smd(
    mean_pre = mean_pre_exp, mean_post = mean_exp,
    mean_pre_sd = mean_pre_sd_exp, mean_post_sd = mean_sd_exp,
    n = n_exp, r_pre_post = r_pre_post_exp,
    pre_post_to_smd = pre_post_to_smd
  )

  res_nexp <- .single_group_pre_post_to_smd(
    mean_pre = mean_pre_nexp, mean_post = mean_nexp,
    mean_pre_sd = mean_pre_sd_nexp, mean_post_sd = mean_sd_nexp,
    n = n_nexp, r_pre_post = r_pre_post_nexp,
    pre_post_to_smd = pre_post_to_smd
  )

  d_final <- res_exp[,"d"] - res_nexp[,"d"]
  g_final <- res_exp[,"g"] - res_nexp[,"g"]
  vd_final <- res_exp[,"var_d"] + res_nexp[,"var_d"]
  vg_final <- res_exp[,"var_g"] + res_nexp[,"var_g"]

  d_ci_lo <- d_final - sqrt(vd_final) * qt(.975, n_exp + n_nexp - 2)
  d_ci_up <- d_final + sqrt(vd_final) * qt(.975, n_exp + n_nexp - 2)
  g_ci_lo <- g_final - sqrt(vg_final) * qt(.975, n_exp + n_nexp - 2)
  g_ci_up <- g_final + sqrt(vg_final) * qt(.975, n_exp + n_nexp - 2)

  res <- cbind(
    d_final, vd_final, d_ci_lo, d_ci_up,
    g_final, vg_final, g_ci_lo, g_ci_up
  )

  return(res)
}

################# POOLED TWO-GROUP PRE POST to SMD ##############
#' Between-group pre/post SMD, standardizing SD pooled across arms
#'
#' All four branches follow the single rule that metafor's own escalc() encodes for
#' every pre/post measure (SMD, SMCC, SMCR, SMCRH, SMCRP):
#'
#'     Var(g) = Var(numerator) / SD_standardizer^2  +  g^2 / (2 * nu_std)
#'     J      = J(nu_std)
#'
#' with NO leading J^2 factor, and Var(d) obtained by substituting d for g (NOT by
#' dividing Var(g) by J^2 -- that is only an identity under the LS2 convention).
#' nu_std is the degrees of freedom of the STANDARDIZER, which is what the g^2 term
#' is divided by: N-2 for an SD pooled over both arms, and 2(N-2)/(1+r^2) for the
#' quadratic-mean-of-two-correlated-SDs standardizer used by d_av (Cousineau, 2020).
#'
#' Var(numerator) is taken from the EMPIRICAL pooled change SD, not from the identity
#' Var(change) = 2*sigma^2*(1-r), which is valid only when SD_pre = SD_post. This is
#' the heteroscedasticity-robust choice (metafor SMCRH/SMCRPH, attributed to Bonett
#' 2008) and it is already what this package's SINGLE-GROUP bonett branch does. It
#' matters: with SD_pre/SD_post = 0.64 (the median in the authors' own PETRA data) the
#' homoscedastic form understates Var(g_bonett) by ~43% and its CI coverage falls to
#' 0.86; the robust form holds coverage at 0.95 in every regime tested.
#'
#' References:
#'   Hedges (1981); Viechtbauer (2007) JEBS 32(1):39-60, eq. 31 -- two-sample SMD "LS"
#'   Bonett (2008) Psych Methods 13(2):99-109 -- heteroscedasticity-robust variance
#'   Cousineau (2020) TQMP 16(4):418-421, eq. 2 -- the (1 + r^2) effective df for d_av
#'   Morris (2008) ORM 11(2):364-386 -- d_ppc2 / d_ppc3 point estimates (eq. 8-14)
#'   Viechtbauer (2007, p.57) recommends the large-sample ("LS") forms over both the
#'   plug-in-exact and the unbiased forms; metafor defaults to LS.
#'
#' NOTE: no published sampling variance exists for the POOLED TWO-GROUP form of d_z,
#' d_rm or d_av. These are delta-method generalizations of the single-group results,
#' obtained by carrying the standardizer's degrees of freedom across the two arms.
#' They are Monte-Carlo calibrated in tests/testthat/test-pooled-variance-calibration.R.
#'
#' @noRd
.pooled_pre_post_to_smd <- function(mean_pre_exp, mean_pre_sd_exp,
                                     mean_exp, mean_sd_exp,
                                     mean_pre_nexp, mean_pre_sd_nexp,
                                     mean_nexp, mean_sd_nexp,
                                     n_exp, n_nexp,
                                     r_pre_post_exp, r_pre_post_nexp,
                                     pre_post_to_smd) {
  if (pre_post_to_smd == "cooper") {
    pre_post_to_smd <- "morris_drm"
  }

  r_pre_post_exp <- .guard_r_pre_post(r_pre_post_exp)
  r_pre_post_nexp <- .guard_r_pre_post(r_pre_post_nexp)

  N <- n_exp + n_nexp
  m <- N - 2 # pooled degrees of freedom
  # m <= 0 breaks J; the variance additionally needs m > 0
  m <- ifelse(is.finite(m) & m > 0, m, NA_real_)

  change_exp <- mean_exp - mean_pre_exp
  change_nexp <- mean_nexp - mean_pre_nexp
  mean_diff <- change_exp - change_nexp

  r_avg <- (n_exp * r_pre_post_exp + n_nexp * r_pre_post_nexp) / N

  # Pooled change SD: the empirical scale of the numerator. Used by every branch,
  # so that no branch assumes SD_pre = SD_post.
  sd_change_exp <- sqrt(mean_pre_sd_exp^2 + mean_sd_exp^2 -
                        2 * r_pre_post_exp * mean_pre_sd_exp * mean_sd_exp)
  sd_change_nexp <- sqrt(mean_pre_sd_nexp^2 + mean_sd_nexp^2 -
                         2 * r_pre_post_nexp * mean_pre_sd_nexp * mean_sd_nexp)
  sd_change_pooled <- sqrt(((n_exp - 1) * sd_change_exp^2 +
                            (n_nexp - 1) * sd_change_nexp^2) / m)

  # Var(mean_diff) on the change scale, in units of sd_change_pooled^2
  T_change <- N / (n_exp * n_nexp)

  if (pre_post_to_smd == "bonett") {
    # Morris (2008) d_ppc2: numerator = difference in mean change, standardizer =
    # BASELINE SD pooled across arms (eq. 8-9). Variance is the SMCRH form: the
    # numerator's empirical variance expressed in standardizer units.
    sd_pooled <- .guard_standardizer(sqrt(((n_exp - 1) * mean_pre_sd_exp^2 +
                                           (n_nexp - 1) * mean_pre_sd_nexp^2) / m))
    nu <- m
    J <- .d_j(nu)
    d <- mean_diff / sd_pooled
    g <- d * J

    T1 <- (sd_change_pooled^2 / sd_pooled^2) * T_change
    var_g <- T1 + g^2 / (2 * N)
    var_d <- T1 + d^2 / (2 * N)

  } else if (pre_post_to_smd == "morris_dz") {
    # Change-score metric. Once the change SD is pooled across arms this is exactly
    # an independent-groups Hedges g computed on the change scores, so its variance
    # is metafor::escalc(measure = "SMD", vtype = "LS") -- no r, no 2(1-r) term.
    sd_pooled <- .guard_standardizer(sd_change_pooled)
    nu <- m
    J <- .d_j(nu)
    d <- mean_diff / sd_pooled
    g <- d * J

    var_g <- T_change + g^2 / (2 * N)
    var_d <- T_change + d^2 / (2 * N)

  } else if (pre_post_to_smd == "morris_drm") {
    # Raw-score metric: d_rm = d_z * sqrt(2(1-r)) (Caldwell & Vigotsky 2020 eq. 13).
    # r is a known constant, so Var(d_rm) = 2(1-r) * Var(d_z): the 2(1-r) here is a
    # deterministic rescaling, NOT a homoscedasticity assumption.
    sd_pooled <- .guard_standardizer(sd_change_pooled)
    nu <- m
    J <- .d_j(nu)
    k <- sqrt(2 * (1 - r_avg))
    d <- (mean_diff / sd_pooled) * k
    g <- d * J

    T1 <- 2 * (1 - r_avg) * T_change
    var_g <- T1 + g^2 / (2 * N)
    var_d <- T1 + d^2 / (2 * N)

  } else if (pre_post_to_smd == "morris_dav") {
    # Morris (2008) d_ppc3: standardizer = quadratic mean of the pre and post SDs,
    # pooled across arms (eq. 12-13). Its effective df is nu = 2m/(1+r^2)
    # (Cousineau 2020 eq. 2; metafor SMCRP uses 2(n-1)/(1+r^2) per arm, and df add
    # across independent arms), so J is evaluated at nu and the g^2 term is divided
    # by 2*nu -- which is exactly the (1 + r^2)/(4N) coefficient below.
    sd_av_exp <- sqrt((mean_pre_sd_exp^2 + mean_sd_exp^2) / 2)
    sd_av_nexp <- sqrt((mean_pre_sd_nexp^2 + mean_sd_nexp^2) / 2)
    sd_pooled <- .guard_standardizer(sqrt(((n_exp - 1) * sd_av_exp^2 +
                                           (n_nexp - 1) * sd_av_nexp^2) / m))
    nu <- 2 * m / (1 + r_avg^2)
    J <- .d_j(nu)
    d <- mean_diff / sd_pooled
    g <- d * J

    T1 <- (sd_change_pooled^2 / sd_pooled^2) * T_change
    var_g <- T1 + g^2 * (1 + r_avg^2) / (4 * N)
    var_d <- T1 + d^2 * (1 + r_avg^2) / (4 * N)
  }

  d_ci_lo <- d - sqrt(var_d) * qt(.975, nu)
  d_ci_up <- d + sqrt(var_d) * qt(.975, nu)
  g_ci_lo <- g - sqrt(var_g) * qt(.975, nu)
  g_ci_up <- g + sqrt(var_g) * qt(.975, nu)

  res <- cbind(
    d, var_d, d_ci_lo, d_ci_up,
    g, var_g, g_ci_lo, g_ci_up
  )

  return(res)
}

################# SINGLE GROUP PRE POST to SMD ##############
#' Calculate within-group standardized mean difference for a single group
#'
#' @param mean_pre mean at baseline (pre-test)
#' @param mean_post mean at follow-up (post-test)
#' @param mean_pre_sd standard deviation at baseline
#' @param mean_post_sd standard deviation at follow-up
#' @param n sample size
#' @param r_pre_post pre-post correlation
#' @param pre_post_to_smd method to use: "bonett" or "cooper"
#'
#' @return matrix with columns: d, var_d, d_ci_lo, d_ci_up, g, var_g, g_ci_lo, g_ci_up
#'
#' @noRd
.single_group_pre_post_to_smd <- function(mean_pre, mean_post,
                                           mean_pre_sd, mean_post_sd,
                                           n, r_pre_post,
                                           pre_post_to_smd) {
  # "cooper" alias kept for direct internal calls
  if (pre_post_to_smd == "cooper") {
    pre_post_to_smd <- "morris_drm"
  }

  r_pre_post <- .guard_r_pre_post(r_pre_post)

  if (pre_post_to_smd == "bonett") {
    # Bonett method: standardize by baseline SD
    # Matches metafor SMCRH (heteroscedastic-robust variance formula)
    J <- .d_j(n - 1)

    var_change <- mean_pre_sd^2 + mean_post_sd^2 - 2 * r_pre_post * mean_pre_sd * mean_post_sd

    sd_std <- .guard_standardizer(mean_pre_sd)
    d <- (mean_post - mean_pre) / sd_std
    g <- d * J

    # metafor SMCRH heteroscedastic variance formula (Bonett 2008)
    # var = sd_change^2 / (sd1i^2 * (n-1)) + g^2 / (2 * (n-1))
    var_g <- var_change / (sd_std^2 * (n - 1)) + g^2 / (2 * (n - 1))
    var_d <- var_g / (J^2)

    d_ci_lo <- d - sqrt(var_d) * qt(.975, n - 1)
    d_ci_up <- d + sqrt(var_d) * qt(.975, n - 1)
    g_ci_lo <- g - sqrt(var_g) * qt(.975, n - 1)
    g_ci_up <- g + sqrt(var_g) * qt(.975, n - 1)

    res <- cbind(
      d = d, var_d = var_d, d_ci_lo = d_ci_lo, d_ci_up = d_ci_up,
      g = g, var_g = var_g, g_ci_lo = g_ci_lo, g_ci_up = g_ci_up
    )

    return(res)
  } else if (pre_post_to_smd == "morris_drm") {
    # Morris d_rm (alias "cooper"): standardize by change SD, raw-score correction
    # nb: the mean_change wrappers pass sd_change through mean_post_sd
    # (mean_pre = 0 and mean_pre_sd = 0, so sd_diff reduces to sd_change)
    J <- .d_j(n - 1)

    sd_diff <- .guard_standardizer(sqrt(mean_pre_sd^2 + mean_post_sd^2 -
                    (2 * r_pre_post * mean_pre_sd * mean_post_sd)))

    d <- (mean_post - mean_pre) / sd_diff * sqrt(2 * (1 - r_pre_post))
    g <- d * J

    # Viechtbauer corrected formula
    var_d <- 2 * (1 - r_pre_post) / n + d^2 / (2 * n)
    var_g <- J^2 * var_d

    d_ci_lo <- d - sqrt(var_d) * qt(.975, n - 1)
    d_ci_up <- d + sqrt(var_d) * qt(.975, n - 1)
    g_ci_lo <- g - sqrt(var_g) * qt(.975, n - 1)
    g_ci_up <- g + sqrt(var_g) * qt(.975, n - 1)

    res <- cbind(
      d = d, var_d = var_d, d_ci_lo = d_ci_lo, d_ci_up = d_ci_up,
      g = g, var_g = var_g, g_ci_lo = g_ci_lo, g_ci_up = g_ci_up
    )

    return(res)
  } else if (pre_post_to_smd == "morris_dz") {
    # Morris & DeShon d_z: standardize by change score SD
    # Matches metafor SMCC (change score standardization)
    # Most conservative when r is high
    J <- .d_j(n - 1)

    sd_diff <- .guard_standardizer(sqrt(mean_pre_sd^2 + mean_post_sd^2 -
                    2 * r_pre_post * mean_pre_sd * mean_post_sd))

    d <- (mean_post - mean_pre) / sd_diff
    g <- d * J

    # metafor SMCC variance formula (uses corrected g, not d)
    # var = 1/n + g^2/(2n)
    var_g <- 1 / n + g^2 / (2 * n)
    var_d <- var_g / (J^2)

    d_ci_lo <- d - sqrt(var_d) * qt(.975, n - 1)
    d_ci_up <- d + sqrt(var_d) * qt(.975, n - 1)
    g_ci_lo <- g - sqrt(var_g) * qt(.975, n - 1)
    g_ci_up <- g + sqrt(var_g) * qt(.975, n - 1)

    res <- cbind(
      d = d, var_d = var_d, d_ci_lo = d_ci_lo, d_ci_up = d_ci_up,
      g = g, var_g = var_g, g_ci_lo = g_ci_lo, g_ci_up = g_ci_up
    )

    return(res)
  } else if (pre_post_to_smd == "morris_dav") {
    # Morris & DeShon d_av: standardize by average SD
    # Recommended by Morris (2008) as general-purpose measure
    # Robust to variance heterogeneity
    # modified df mi = 2*(n-1)/(1+r^2), matches metafor::escalc(measure = "SMCRP")
    mi <- 2 * (n - 1) / (1 + r_pre_post^2)
    J <- .d_j(mi)

    sd_av <- .guard_standardizer(sqrt((mean_pre_sd^2 + mean_post_sd^2) / 2))

    d <- (mean_post - mean_pre) / sd_av
    g <- d * J

    # Variance formula from metafor SMCRP (uses corrected g in formula)
    # vi = 2*(1-r)/n + g^2*(1+r^2)/(4*n)
    var_g <- 2 * (1 - r_pre_post) / n + g^2 * (1 + r_pre_post^2) / (4 * n)
    var_d <- var_g / (J^2)

    d_ci_lo <- d - sqrt(var_d) * qt(.975, n - 1)
    d_ci_up <- d + sqrt(var_d) * qt(.975, n - 1)
    g_ci_lo <- g - sqrt(var_g) * qt(.975, n - 1)
    g_ci_up <- g + sqrt(var_g) * qt(.975, n - 1)

    res <- cbind(
      d = d, var_d = var_d, d_ci_lo = d_ci_lo, d_ci_up = d_ci_up,
      g = g, var_g = var_g, g_ci_lo = g_ci_lo, g_ci_up = g_ci_up
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
    d_se <- sqrt(4 * r_se^2 / ((1 - r^2)^3))
    res <- cbind(d, d_se)

    return(res)
  }
}



.validate_pre_post_to_smd <- function(pre_post_to_smd, allowed_methods, context, func_name) {
  pre_post_to_smd_normalized <- ifelse(pre_post_to_smd == "cooper", "morris_drm", pre_post_to_smd)

  invalid <- !pre_post_to_smd_normalized %in% allowed_methods
  if (any(invalid)) {
    stop(paste0(
      "Invalid 'pre_post_to_smd' argument in ", func_name, "().\n\n",
      "Provided: '", paste(unique(pre_post_to_smd[invalid]), collapse = "', '"), "'\n",
      "Allowed methods: '", paste(allowed_methods, collapse = "', '"), "'\n\n",
      .get_method_restriction_rationale(context)
    ), call. = FALSE)
  }

  return(pre_post_to_smd_normalized)
}

.get_method_restriction_rationale <- function(context) {
  if (context == "mean_change") {
    paste0(
      "For mean change data, two standardization methods are available:\n\n",
      "  - 'morris_drm' (alias: 'cooper'): Raw score standardizer [DEFAULT]\n",
      "      d_rm = (mean_change / sd_change) * sqrt(2*(1-r))\n",
      "      REQUIRES pre-post correlation (r). Converts d_z to d_rm.\n\n",
      "  - 'morris_dz': Change score standardizer\n",
      "      d_z = mean_change / sd_change\n",
      "      INDEPENDENT of r. Direct standardization by change SD.\n\n",
      "Other methods not applicable:\n",
      "  - 'bonett': requires baseline SD (not available with mean change data)\n",
      "  - 'morris_dav': requires separate pre/post SDs (not available)\n\n",
      "See Morris & DeShon (2002) for guidance on choosing between d_rm and d_z."
    )
  } else if (context == "paired_t") {
    paste0(
      "For paired t-test data, two standardization methods are available:\n\n",
      "  - 'morris_drm' (alias: 'cooper'): Raw score standardizer [DEFAULT]\n",
      "      d_rm = t * sqrt(2*(1-r)/n)\n",
      "      REQUIRES pre-post correlation (r). Most common in meta-analysis.\n\n",
      "  - 'morris_dz': Change score standardizer\n",
      "      d_z = t / sqrt(n)\n",
      "      INDEPENDENT of r. Use when correlation is unknown or when\n",
      "      synthesizing with other d_z estimates.\n\n",
      "Note: d_rm and d_z are on different scales. See Morris & DeShon (2002)\n",
      "for guidance on choosing between them."
    )
  } else if (context == "pre_post_means") {
    paste0(
      "For pre-post means data, four standardization methods are available:\n\n",
      "  - 'bonett': baseline-SD standardizer (Morris 2008 d_ppc2). Sensitive to a\n",
      "      baseline SD that is restricted relative to the endpoint SD (e.g. by\n",
      "      eligibility cut-offs), which inflates the effect size.\n",
      "  - 'morris_drm' (alias: 'cooper'): change SD rescaled by sqrt(2(1-r)) onto\n",
      "      the raw-score metric. Depends on r_pre_post; assumes SD_pre = SD_post.\n",
      "  - 'morris_dz': change-SD standardizer. Independent of r, but NOT on the\n",
      "      same metric as an endpoint SMD -- do not pool the two (Cochrane 10.5.2).\n",
      "  - 'morris_dav': average-SD standardizer (Morris 2008 d_ppc3). Morris\n",
      "      recommends AGAINST it (2008, p.384): its sampling variance was unknown\n",
      "      and it is downward-biased when the post-treatment SD inflates.\n\n",
      "No single method is best in all cases; the choice changes the estimand.\n",
      "See Morris (2008) and Morris & DeShon (2002) for detailed comparisons."
    )
  } else {
    ""
  }
}


# tryCatch({
#   validate_positive(n_cases_exp, n_cases_nexp, n_controls_exp, n_controls_nexp,
#                     error_message = "The number of cases/controls in the exposed/non-exposed groups should be >0.")
# }, error = function(e) {
#   stop("Validation failed:", conditionMessage(e), "\n")
# })
#
#
# tryCatch({
#   validate_ci_symmetry(value, ci_lo, ci_up, func = "example_function",
#                        max_asymmetry_percent = 5)
# }, error = function(e) {
#   stop("Validation failed:", conditionMessage(e), "\n")
# })
#


# **A.** First, Cooper et al. (2019) - \code{table_2x2_to_cor = "cooper"} -
# proposes to convert the
# 2x2 table into a OR (formula above), to convert this OR into a SMD
# (see formula in \code{\link{es_from_or_se}()}), and to convert this
# SMD into a correlation coefficient (see formula in \code{\link{es_from_cohen_d}()},
# with the option \code{"smd_to_cor = 'lipsey_cooper'"}).
#
# **B.** Second, a correlation coefficient (more precisely - a phi coefficient)
# can be obtained from the contingency table using the formula given in
# Lipsey and Wilson (2001) - \code{table_2x2_to_cor = "lipsey"}.
# The formulas used to estimate the r and z are:
# \deqn{r = \frac{(n\_cases\_exp*n\_controls\_nexp - n\_controls\_exp*n\_cases\_nexp)}{\sqrt{(n\_exp) * (n\_nexp) * (n\_cases) * (n\_controls\_exp+n\_cases\_nexp)}}}
# \deqn{z = atanh(r)}
# \deqn{z\_se = logor\_se^2 * \frac{z^2}{\log(or)^2}}
# \deqn{z\_ci\_lo = z - qnorm(.975)*z\_se}
# \deqn{z\_ci\_up = z + qnorm(.975)*z\_se}
# \deqn{r\_ci\_lo = tanh(z\_ci\_lo)}
# \deqn{r\_ci\_up = tanh(z\_ci\_up)}
# \deqn{effective\_n = \frac{1}{z\_se^2 + 3}}
# \deqn{r\_se = \frac{(1 - r^2)^2}{effective\_n - 1}}

