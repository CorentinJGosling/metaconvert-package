.quick_val_ci <- function(value, ci_lo, ci_up) {
  # First check if any inputs are length 0
  if (length(value) == 0 || length(ci_lo) == 0 || length(ci_up) == 0) {
    return(FALSE)
  }

  lengths <- c(length(value), length(ci_lo), length(ci_up))
  if (length(unique(lengths)) > 1) {
    return(FALSE)
  }

  # Proceed with the original logic
  mapply(function(v, l, u) {
    if (is.na(v) || is.na(l) || is.na(u)) return(FALSE)
    return((v - l) > 0 && (u - v) > 0)
  }, value, ci_lo, ci_up)
}

.extract_info <- function(x) {
  row <- which((!is.na(x$d) & !is.na(x$d_se) & .quick_val_ci(x$d, x$d_ci_lo, x$d_ci_up)) |
    (!is.na(x$g) & !is.na(x$g_se) & .quick_val_ci(x$g, x$g_ci_lo, x$g_ci_up)) |
    (!is.na(x$md) & !is.na(x$md_se) & .quick_val_ci(x$md, x$md_ci_lo, x$md_ci_up)) |
    (!is.na(x$r) & !is.na(x$r_se) & .quick_val_ci(x$r, x$r_ci_lo, x$r_ci_up)) |
    (!is.na(x$z) & !is.na(x$z_se) & .quick_val_ci(x$z, x$z_ci_lo, x$z_ci_up)) |
    (!is.na(x$logor) & !is.na(x$logor_se) & .quick_val_ci(x$logor, x$logor_ci_lo, x$logor_ci_up)) |
    (!is.na(x$logrr) & !is.na(x$logrr_se) & .quick_val_ci(x$logrr, x$logrr_ci_lo, x$logrr_ci_up)) |
    (!is.na(x$logirr) & !is.na(x$logirr_se) & .quick_val_ci(x$logirr, x$logirr_ci_lo, x$logirr_ci_up)) |
    (!is.na(x$logcvr) & !is.na(x$logcvr_se) & .quick_val_ci(x$logcvr, x$logcvr_ci_lo, x$logcvr_ci_up)) |
    (!is.na(x$logvr) & !is.na(x$logvr_se) & .quick_val_ci(x$logvr, x$logvr_ci_lo, x$logvr_ci_up)))
  info <- rep(NA, nrow(x))
  if (length(row) > 0) {
    info[row] <- x$info_used[row]
  }
  return(info)
}

.add_columns <- function(x, y) {
  cols <- y[which(!y %in% colnames(x))]

  if (length(cols) > 0) {
    x[, cols] <- NA
  }
  return(x)
}

.extract_es <- function(x, measure, exp) {
  if (measure == "d") {
    if ("d" %in% colnames(x)) {
      x$value <- x$d
      x$se_value <- x$d_se
      x$value_ci_lo <- x$d_ci_lo
      x$value_ci_up <- x$d_ci_up
    } else {
      x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
    }
  } else if (measure == "g") {
    if ("g" %in% colnames(x)) {
      x$value <- x$g
      x$se_value <- x$g_se
      x$value_ci_lo <- x$g_ci_lo
      x$value_ci_up <- x$g_ci_up
    } else {
      x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
    }
  } else if (measure == "z") {
    if ("z" %in% colnames(x)) {
      x$value <- x$z
      x$se_value <- x$z_se
      x$value_ci_lo <- x$z_ci_lo
      x$value_ci_up <- x$z_ci_up
    } else {
      x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
    }
  } else if (measure == "r") {
    if ("r" %in% colnames(x)) {
      x$value <- x$r
      x$se_value <- x$r_se
      x$value_ci_lo <- x$r_ci_lo
      x$value_ci_up <- x$r_ci_up
    } else {
      x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
    }
  } else if (measure == "logor") {
    if (exp) {
      if ("logor" %in% colnames(x)) {
        x$value <- exp(as.numeric(as.character(x$logor)))
        x$se_value <- x$logor_se
        x$value_ci_lo <- exp(as.numeric(as.character(x$logor_ci_lo)))
        x$value_ci_up <- exp(as.numeric(as.character(x$logor_ci_up)))
      } else {
        x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
      }
    } else {
      if ("logor" %in% colnames(x)) {
        x$value <- x$logor
        x$se_value <- x$logor_se
        x$value_ci_lo <- x$logor_ci_lo
        x$value_ci_up <- x$logor_ci_up
      } else {
        x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
      }
    }
  } else if (measure == "logrr") {
    if (exp) {
      if ("logrr" %in% colnames(x)) {
        x$value <- exp(as.numeric(as.character(x$logrr)))
        x$se_value <- x$logrr_se
        x$value_ci_lo <- exp(as.numeric(as.character(x$logrr_ci_lo)))
        x$value_ci_up <- exp(as.numeric(as.character(x$logrr_ci_up)))
      } else {
        x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
      }
    } else {
      if ("logrr" %in% colnames(x)) {
        x$value <- x$logrr
        x$se_value <- x$logrr_se
        x$value_ci_lo <- x$logrr_ci_lo
        x$value_ci_up <- x$logrr_ci_up
      } else {
        x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
      }
    }
  } else if (measure == "logirr") {
    if (exp) {
      if ("logirr" %in% colnames(x)) {
        x$value <- exp(as.numeric(as.character(x$logirr)))
        x$se_value <- x$logirr_se
        x$value_ci_lo <- exp(as.numeric(as.character(x$logirr_ci_lo)))
        x$value_ci_up <- exp(as.numeric(as.character(x$logirr_ci_up)))
      } else {
        x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
      }
    } else {
      if ("logirr" %in% colnames(x)) {
        x$value <- x$logirr
        x$se_value <- x$logirr_se
        x$value_ci_lo <- x$logirr_ci_lo
        x$value_ci_up <- x$logirr_ci_up
      } else {
        x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
      }
    }
  } else if (measure == "logcvr") {
    if ("logcvr" %in% colnames(x)) {
      x$value <- x$logcvr
      x$se_value <- x$logcvr_se
      x$value_ci_lo <- x$logcvr_ci_lo
      x$value_ci_up <- x$logcvr_ci_up
    } else {
      x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
    }
  } else if (measure == "logvr") {
    if ("logvr" %in% colnames(x)) {
      x$value <- x$logvr
      x$se_value <- x$logvr_se
      x$value_ci_lo <- x$logvr_ci_lo
      x$value_ci_up <- x$logvr_ci_up
    } else {
      x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
    }
  } else if (measure == "md") {
    if ("md" %in% colnames(x)) {
      x$value <- x$md
      x$se_value <- x$md_se
      x$value_ci_lo <- x$md_ci_lo
      x$value_ci_up <- x$md_ci_up
    } else {
      x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
    }
  } else if (measure == "nnt") {
    if ("nnt" %in% colnames(x)) {
      x$value <- x$nnt
      x$se_value <- rep(NA, nrow(x))
      x$value_ci_lo <- rep(NA, nrow(x))
      x$value_ci_up <- rep(NA, nrow(x))
    } else {
      x$value <- x$se_value <- x$value_ci_lo <- x$value_ci_up <- rep(NA_real_, nrow(x))
    }
  }
  return(x)
}


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
.generate_df <- function(x, main_es = TRUE, list_df_es_enh, list_df, ordering,
                         measure, suffix = "",
                         es_selected,
                         exp, digits) {
  list_keep <- do.call(cbind, lapply(list_df, function(x) x$info_used))
  cols_keep <- which(as.character(list_keep[1, ]) %in% ordering)
  list_df_restrict <- list_df[cols_keep]

  # build datasets with all ES / SE / information ----------------------------
  df_value_min_max <- do.call(cbind, lapply(list_df_restrict, function(x) x$value))
  df_se_min_max <- do.call(cbind, lapply(list_df_restrict, function(x) x$se_value))
  df_ci_lo_min_max <- do.call(cbind, lapply(list_df_restrict, function(x) x$value_ci_lo))
  df_ci_up_min_max <- do.call(cbind, lapply(list_df_restrict, function(x) x$value_ci_up))
  df_info_min_max <- do.call(cbind, lapply(list_df_restrict, function(x) x$info_used))


  # select the requested ES to highlight BASED ON HIERARCHY -------------------------------------------------
  # we start by the first ES requested, and if NA, try the second one, etc...
  for (i in 1:length(ordering)) {
    # i = 1
    # print(i);
    # print(ordering[i]);
    df_temp <- list_df_restrict[[which(df_info_min_max[1, ] == ordering[i])]]
    row <- which(
      # info non-missing in the selected results
      !is.na(df_temp$value) & (!is.na(df_temp$se) | measure == "nnt") &
      # info missing in user datasets
        ((((is.na(x[, paste0("es", suffix)]) | is.na(x[, paste0("se", suffix)]))) & measure != "nnt") |
        (is.na(x[, paste0("es", suffix)]) & measure == "nnt"))
      )

    # print(length(row))
    if (length(row) >= 1) {
      x[row, paste0("es", suffix)] <- df_temp$value[row]
      x[row, paste0("se", suffix)] <- df_temp$se_value[row]
      x[row, paste0("es_ci_lo", suffix)] <- df_temp$value_ci_lo[row]
      x[row, paste0("es_ci_up", suffix)] <- df_temp$value_ci_up[row]
      x[row, paste0("info_used", suffix)] <- df_temp$info_used[row]
      x[row, paste0("measure", suffix)] <- measure
      if (grepl("user_input", unique(df_temp$info_used), fixed = TRUE)) {
        x[row, paste0("measure", suffix)] <- attr(df_temp, "measure")[row]
      }
    }
  }

  # rename the measure in case of exp selection
  if (exp & measure %in% c("logrr", "logor", "logirr")) {
    x[, paste0("measure", suffix)][x[, paste0("measure", suffix)] == "logrr" & !is.na(x[, paste0("measure", suffix)])] <- "rr"
    x[, paste0("measure", suffix)][x[, paste0("measure", suffix)] == "logor" & !is.na(x[, paste0("measure", suffix)])] <- "or"
    x[, paste0("measure", suffix)][x[, paste0("measure", suffix)] == "logirr" & !is.na(x[, paste0("measure", suffix)])] <- "irr"
  }

  # From there, the info_used+ES+SE+CIs have been inserted
  # this is for info_measure
  # message("OK")
  for (col in 1:length(ordering)) {
    # for (col in 1:38) {
    for (row in 1:nrow(df_value_min_max)) {
      # print(col)
      # for a given information, if a calculation (ES or SE) failed, we delete the others (CI are necessarily NA)----------------------------
      if ((!is.na(df_value_min_max[row, col]) & is.na(df_se_min_max[row, col]) & measure != "nnt") |
        (is.na(df_value_min_max[row, col]) & !is.na(df_se_min_max[row, col]))) {
        df_value_min_max[row, col] <- df_se_min_max[row, col] <- NA
      }
      # we identify all information from which we can derive an effect size ----------------------------
      if (!is.na(df_value_min_max[row, col]) & (!is.na(df_se_min_max[row, col]) | measure == "nnt")) {
        x[row, paste0("info_measure", suffix)] <- paste(x[row, paste0("info_measure", suffix)],
                                                        df_info_min_max[row, col], sep = " + ")
      }
    }
  }

  x[, paste0("info_measure", suffix)] <- sub("NA", "", x[, paste0("info_measure", suffix)], fixed = TRUE)
  x[, paste0("info_measure", suffix)] <- ifelse(grepl(" + ", substring(x[, paste0("info_measure", suffix)], 1, 3), fixed = TRUE),
    sub(" \\+ ", "", x[, paste0("info_measure", suffix)]),
    x[, paste0("info_measure", suffix)]
  )

  # This is for "all_info"
  info_list_all <- do.call(cbind, lapply(list_df_es_enh[cols_keep], .extract_info))
  all_info <- apply(info_list_all, 1, function(x) paste(na.omit(x), collapse = " + "))
  all_info[all_info == ""] <- "-- None available --"
  no_avlble <- c("", "-- None available --", NA)

  x[, paste0("all_info", suffix)] <- ifelse(
    x[, paste0("all_info", suffix)] %in% no_avlble & all_info %in% no_avlble,
   "-- None available --",
   ifelse(x[, paste0("all_info", suffix)] %in% no_avlble & !all_info %in% no_avlble,
          all_info,
          ifelse(!x[, paste0("all_info", suffix)] %in% no_avlble & all_info %in% no_avlble,
                 x[, paste0("all_info", suffix)],
                 ifelse(!x[, paste0("all_info", suffix)] %in% no_avlble & !all_info %in% no_avlble,
                        paste0(x[, paste0("all_info", suffix)], " + ", all_info),
                        "Error - contact corentin.gosling@parisnanterre.fr for more information"
                 )
          )
     )
    )

  # This is for "n_estimation"
  x[, paste0("n_estimations", suffix)] <- ncol(df_value_min_max) - rowSums(is.na(df_value_min_max))

  # This is for min_** and max_***
  col_min <- apply(df_value_min_max, 1, function(x) which.min(x)[1])
  col_max <- apply(df_value_min_max, 1, function(x) which.max(x)[1])

  cols_equal = which(col_min == col_max)
  col_min[cols_equal] <- apply(data.frame(df_ci_lo_min_max)[cols_equal, ], 1, function(x) which.min(x)[1])
  col_max[cols_equal] <- apply(data.frame(df_ci_lo_min_max)[cols_equal, ], 1, function(x) which.max(x)[1])

  for (i in 1:length(col_min)) {
    x[i, paste0("min_info", suffix)] <- df_info_min_max[i, col_min[i]]
    x[i, paste0("min_es_value", suffix)] <- as.numeric(as.character(df_value_min_max[i, col_min[i]]))
    x[i, paste0("min_es_se", suffix)] <- as.numeric(as.character(df_se_min_max[i, col_min[i]]))
    x[i, paste0("min_es_ci_lo", suffix)] <- as.numeric(as.character(df_ci_lo_min_max[i, col_min[i]]))
    x[i, paste0("min_es_ci_up", suffix)] <- as.numeric(as.character(df_ci_up_min_max[i, col_min[i]]))

    x[i, paste0("max_info", suffix)] <- df_info_min_max[i, col_max[i]]
    x[i, paste0("max_es_value", suffix)] <- as.numeric(as.character(df_value_min_max[i, col_max[i]]))
    x[i, paste0("max_es_se", suffix)] <- as.numeric(as.character(df_se_min_max[i, col_max[i]]))
    x[i, paste0("max_es_ci_lo", suffix)] <- as.numeric(as.character(df_ci_lo_min_max[i, col_max[i]]))
    x[i, paste0("max_es_ci_up", suffix)] <- as.numeric(as.character(df_ci_up_min_max[i, col_max[i]]))
  }

  # This is for "diff_min_max"
  x[, paste0("diff_min_max", suffix)] <- ifelse(x[, paste0("n_estimations", suffix)] > 1,
                                                x[, paste0("max_es_value", suffix)] - x[, paste0("min_es_value", suffix)],
                                                "< 2 types of input data available"
  )

  # This is for "overlap_min_max"
  x[, paste0("overlap_min_max", suffix)] <- ifelse(
    x[, paste0("n_estimations", suffix)] > 1,
    ifelse(x[, paste0("max_es_ci_lo", suffix)] > x[, paste0("min_es_ci_up", suffix)],
           0,
           abs(x[, paste0("max_es_ci_lo", suffix)] - x[, paste0("min_es_ci_up", suffix)]) /
             abs(x[, paste0("max_es_ci_up", suffix)] - x[, paste0("min_es_ci_lo", suffix)])
    ),
    "< 2 types of input data available"
  )
  # This is for dispersion_es
  dat_long = rep(NA, nrow(df_value_min_max))
  i = 0
  for (dat in c("df_value_min_max", "df_se_min_max",
                "df_ci_lo_min_max", "df_ci_up_min_max",
                "df_info_min_max")) {
    # print(dat)
    i = i+1
    df_tot = data.frame(mget(dat))
    df_tot$row_id = x$row_id
    dat_long_transit = reshape(
      df_tot, direction = "long",
      varying = names(df_tot)[which(names(df_tot) != "row_id")],
      times = names(df_tot)[which(names(df_tot) != "row_id")],
      v.names = "es",
      timevar = "row_id")
    colnames(dat_long_transit) <- c(c("info_measure_es_used", "info_measure_se",
                                      "info_measure_ci_lo", "info_measure_ci_up", "info_measure_used")[i],
                                    c("es", "se", "es_ci_lo", "es_ci_up", "info_used")[i],
                                    c("row_id", "row_id2", "row_id3", "row_id4", "row_id5")[i])
    dat_long = cbind(dat_long, dat_long_transit)
  }
  dat_long = dat_long[(!is.na(dat_long$es) & !is.na(dat_long$se) &
                         !is.na(dat_long$es_ci_lo) & !is.na(dat_long$es_ci_up) &
                         rep(measure, nrow(dat_long)) != "nnt") |
                        (!is.na(dat_long$es) & measure == "nnt"), ]
  dat_long = dat_long[order(dat_long$row_id),]

  res_dispersion = data.frame(tapply(dat_long$es, dat_long$row_id, sd))
  dispersion = data.frame(dispersion_es = res_dispersion[,1],
                          row_id = rownames(res_dispersion))

  dat_long = merge(x = dat_long,
                   y = dispersion)

  if (main_es == TRUE & nrow(dat_long) != 0) {
    dispersion$blank = NA
    dispersion$row_id = as.numeric(as.character(dispersion$row_id))
    x_transit = merge(x = dispersion[,c("row_id", "blank")],
                      y = x)
    if (all(dispersion$row_id == x_transit$row_id)) {
      x_transit[, paste0("dispersion_es", suffix)] <- dispersion$dispersion_es
    } else {
      stop("an error occured when estimating the 'dispersion_es' variable")
    }

    x_transit <- x_transit[, -which(names(x_transit) == "blank")]

    x_empty = x[!x$row_id %in% x_transit$row_id, ]

    x = rbind(x_transit, x_empty)

    x = x[order(x$row_id),]

  } else if (nrow(dat_long) != 0) {
    x_transit = merge(x = dat_long[, c("dat_long", "row_id")], y = x)

    x_transit[, paste0("info_used", suffix)] <- dat_long$info_used
    x_transit[, paste0("es", suffix)] <- dat_long$es
    x_transit[, paste0("se", suffix)] <- dat_long$se
    x_transit[, paste0("es_ci_up", suffix)] <- dat_long$es_ci_up
    x_transit[, paste0("es_ci_lo", suffix)] <- dat_long$es_ci_lo
    x_transit[, paste0("dispersion_es", suffix)] <- dat_long$dispersion_es

    x_transit <- x_transit[, -which(names(x_transit) == "dat_long")]

    x_empty = x[!x$row_id %in% x_transit$row_id, ]

    x = rbind(x_transit, x_empty)

    x = x[order(x$row_id),]

  }

  # INSERT FINAL ES if users prefer using minimum/maximum rather than hierarchy
  if (es_selected == "minimum") {
    x[, paste0("es_selected", suffix)] <- "minimum"
    x[, paste0("es", suffix)] <- x[, paste0("min_es_value", suffix)]
    x[, paste0("se", suffix)] <- x[, paste0("min_es_se", suffix)]
    x[, paste0("es_ci_lo", suffix)] <- x[, paste0("min_es_ci_lo", suffix)]
    x[, paste0("es_ci_up", suffix)] <- x[, paste0("min_es_ci_up", suffix)]
    x[, paste0("info_used", suffix)] <- x[, paste0("min_info", suffix)]
    x[, paste0("measure", suffix)] <- x[, paste0("measure", suffix)]
  } else if (es_selected == "maximum") {
    x[, paste0("es_selected", suffix)] <- "maximum"
    x[, paste0("es", suffix)] <- x[, paste0("max_es_value", suffix)]
    x[, paste0("se", suffix)] <- x[, paste0("max_es_se", suffix)]
    x[, paste0("es_ci_lo", suffix)] <- x[, paste0("max_es_ci_lo", suffix)]
    x[, paste0("es_ci_up", suffix)] <- x[, paste0("max_es_ci_up", suffix)]
    x[, paste0("info_used", suffix)] <- x[, paste0("max_info", suffix)]
    x[, paste0("measure", suffix)] <- x[, paste0("measure", suffix)]
  } else {
    x[, paste0("es_selected", suffix)] <- "hierarchy"
  }

  # this is for securizing certain columns.
  # set NA when n_estimations = 0
  row_miss <- which(x[, paste0("n_estimations", suffix)] == 0)
  x[row_miss, paste0(
    c(
      "overlap_min_max", "diff_min_max",
      "min_info", "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
      "max_info", "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up"
    ),
    suffix
  )] <- NA
  # set "< 2 types of input data available" when n_estimations = 1
  row_1 <- which(x[, paste0("n_estimations", suffix)] == 1)
  x[row_1, paste0(
    c("dispersion_es",
      "min_info", "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
      "max_info", "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up"
    ),
    suffix
  )] <- "< 2 types of input data available"

  # this is for rounding
  for (cols in paste0(
    c(
      "es", "se", "es_ci_lo", "es_ci_up",
      "overlap_min_max", "diff_min_max", "dispersion_es",
      "min_es_value", "min_es_se", "min_es_ci_lo", "min_es_ci_up",
      "max_es_value", "max_es_se", "max_es_ci_lo", "max_es_ci_up"
    ),
    suffix
  )) {
    for (rows in which(x[, cols] != "< 2 types of input data available")) {
      x[rows, cols] <- # as.numeric(as.character(
        # sprintf(paste0("%.", digits, "f"),
        round(as.numeric(as.character(
          x[rows, cols]
        )), digits) # )))
    }
  }

  return(x)
}
