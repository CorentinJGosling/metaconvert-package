#' Aggregate a dataframe containing dependent effect sizes
#'
#' @param x a dataframe that should be aggregated (must contain effect size values and standard errors).
#' @param dependence The type of dependence in your dataframe (can be "outcomes", "times", or "subgroups"). See details.
#' @param agg_fact A character string identifying the column name that contains the clustering units (all rows with the same \code{agg_fact} value will be aggregated together).
#' @param es A character string identifying the column name containing the effect size values. Default is "es".
#' @param se A character string identifying the column name containing the standard errors of the effect size. Default is "se".
#' @param cor_unit The correlation between effect sizes coming from the same clustering unit. Used directly when \code{dependence = "outcomes"}. When \code{dependence = "times"} the within-cluster correlation is taken from a per-row \code{cor_unit} column if present (which must be constant within each cluster), and this argument is only used as a fallback when that column is absent; \code{dependence = "times"} additionally requires a \code{time_agg} column giving the time-point of each effect size. Ignored when \code{dependence = "subgroups"}.
#' @param col_mean a vector of character strings identifying the column names for which the dependent values are summarized by taking their mean.
#' @param col_weighted_mean a vector of character strings identifying the column names for which the dependent values are summarized by taking their weighted mean.
#' @param weights The weights that will be used to estimated the weighted means.
#' @param col_sum a vector of character strings identifying the column names for which the dependent values are summarized by taking their sum.
#' @param col_min a vector of character strings identifying the column names for which the dependent values are summarized by taking their minimum.
#' @param col_max a vector of character strings identifying the column names for which the dependent values are summarized by taking their maximum.
#' @param col_fact a vector of character strings identifying the column names that are factors (different values will be separated by a "/" character).
#' @param na.rm a logical vector indicating whether missing values should be ignored in the calculations for the \code{col_mean}, \code{col_weighted_mean}, \code{col_sum}, \code{col_min} and \code{col_max} arguments.
#'
#' @details
#' 1. In the \code{dependence} argument, you should indicate "outcomes" if the dependence within the same clustering unit (e.g., study) is due to the presence of multiple effect sizes produced from the same participants at the same time-point (e.g., multiple outcome measures)
#' 2. In the \code{dependence} argument, you should indicate "times" if the dependence within the same clustering unit (e.g., study) is due to the presence of multiple effect sizes produced from the same participants at the different time-points (e.g., an RCT with several follow-up waves). This option requires a \code{time_agg} column (the time-point of each effect size) and a \code{cor_unit} column giving the within-cluster correlation (constant within each cluster; falls back to the \code{cor_unit} argument if the column is absent).
#' 3. In the \code{dependence} argument, you should indicate "subgroups" if the dependence within the same clustering unit (e.g., study) is due to the presence of multiple effect sizes produced by independent subgroups (e.g., one effect size for boys, and one for girls).
#'
#' **Which columns survive.** `aggregate_df()` returns the aggregating factor, the
#' aggregated effect size and its standard error, and *nothing else*. Every other
#' column of the input is dropped, including the moderators a meta-regression will
#' need. This is deliberate - there is no defined value for a column that varies
#' within a cluster, so the function refuses to guess one - but it is silent, so it is
#' worth stating: a 40-column `summary()` frame comes back with four columns.
#'
#' Two ways to keep a moderator, and the choice depends on the column, not on taste:
#'
#' - name it in `col_mean` / `col_weighted_mean` / `col_sum` / `col_min` / `col_max`
#'   (numeric, summarised across the cluster) or `col_fact` (categorical, and it must
#'   be constant within the cluster);
#' - or merge it back afterwards on the aggregating factor, which is the simpler route
#'   for study-level variables that cannot vary within a study anyway (publication
#'   year, country, risk-of-bias rating).
#'
#' If you are working with ratio measures, make sure that the information on
#' the effect size estimates (i.e., the column passed to the es argument of the function)
#' is presented on the log scale.
#'
#' @return
#' A dataframe containing the aggregating factor, the aggregated effect size and its
#' standard error, plus any columns named in the \code{col_*} arguments -- and
#' \strong{no other column of the input}. See the note in the details: everything not
#' listed below is dropped, so a moderator must either be named in a \code{col_*}
#' argument or merged back afterwards on the aggregating factor.
#' \tabular{ll}{
#'  \code{row_index} \tab the row number in the original dataset.\cr
#'  \tab \cr
#'  \code{es} \tab the aggregated effect size value.\cr
#'  \tab \cr
#'  \code{se} \tab the standard error of the aggregated effect size.\cr
#'  \tab \cr
#'  \code{...} \tab any columns indicated in the \code{col_*} arguments.\cr
#'  \tab \cr
#' }
#'
#' @export aggregate_df
#'
#' @md
#'
#' @examples
#' res <- summary(convert_df(df.haza, measure = "d"))
#' aggregate_df(res, dependence = "outcomes", cor_unit = 0.8,
#'              agg_fact = "study_id", es = "es_crude", se = "se_crude",
#'              col_fact = c("outcome", "type_publication"))
#'
aggregate_df <- function(x, dependence = "outcomes", cor_unit = 0.8,
                         agg_fact, es = "es", se = "se",
                         col_mean = NA, col_weighted_mean = NA, weights = NA,
                         col_sum = NA, col_min = NA, col_max = NA,
                         col_fact = NA, na.rm = TRUE) {
  # x=x_mult_outcome; dependence = "outcomes"; cor_unit = 0.8;
  # agg_fact = "ID"; es = "value"; se = "se";
  # col_mean = NA; col_weighted_mean = NA; weights = NA;
  # col_sum = NA; col_min = NA; col_max = NA;
  # col_fact = "row_index"; na.rm = TRUE
  if (!dependence %in% c("outcomes", "times", "subgroups")) {
    stop(paste0("'", dependence, "' not in tolerated measures. Possible inputs are: 'outcomes', 'times', or 'subgroups'."))
  }

  if ((missing(es) & !"es" %in% colnames(x)) | (missing(se) & !"se" %in% colnames(x))) {
    stop(paste0("Please indicate the column names storing information on the effect size values and standard errors using the 'es' & 'se' arguments"))
  }
  if ((missing(es) & "es" %in% colnames(x))) {
    es <- "es"
  }
  if ((missing(se) & "se" %in% colnames(x))) {
    se <- "se"
  }

  if (length(es) > 1 | length(se) > 1) {
    stop(paste0("The '", es, "' and the '", se, "' arguments should contain only 1 element"))
  }
  if (!es %in% colnames(x) | !se %in% colnames(x)) {
    stop(paste0("The '", es, "' or the '", se, "' variable has not been found in your dataset"))
  }

  if (dependence == "outcomes") {
    res <- .agg.outcomes(
      x = x, agg_fact = agg_fact, es = es, se = se,
      cor_unit = cor_unit,
      col_mean = col_mean,
      col_weighted_mean = col_weighted_mean, weights = weights,
      col_sum = col_sum, col_min = col_min, col_max = col_max,
      col_fact = col_fact, na.rm = na.rm
    )
  } else if (dependence == "times") {
    res <- .agg.times(
      x = x, agg_fact = agg_fact, es = es, se = se,
      cor_unit = cor_unit,
      col_mean = col_mean,
      col_weighted_mean = col_weighted_mean, weights = weights,
      col_sum = col_sum, col_min = col_min, col_max = col_max,
      col_fact = col_fact, na.rm = na.rm
    )
  } else {
    res <- .agg.subgroups(
      x = x, agg_fact = agg_fact, es = es, se = se,
      col_mean = col_mean,
      col_weighted_mean = col_weighted_mean, weights = weights,
      col_sum = col_sum, col_min = col_min, col_max = col_max,
      col_fact = col_fact, na.rm = na.rm
    )
  }

  rownames(res) <- 1:nrow(res)
  return(res)
}

.agg.subgroups <- function(x, agg_fact, es, se,
                           col_mean = NA, col_weighted_mean = NA, weights = NA,
                           col_sum = NA, col_min = NA, col_max = NA,
                           col_fact = NA, na.rm = TRUE) {
  cols <- c(col_mean, col_weighted_mean, col_sum, col_min, col_max, col_fact)

  dup_col <- ifelse(length(na.omit(cols)) > 0, any(duplicated(na.omit(cols))), FALSE)

  if (dup_col) {
    stop("Each column of your dataset can be aggregated only once.")
  }

  if (any(!is.na(col_weighted_mean) & is.na(weights))) {
    stop("You should specify the 'weights' argument when requesting to calculate a weighted mean for some columns.")
  }

  df_mean <- df_w_mean <- df_sum <- df_min <- df_max <- df_fact <- NA

  dup_row <- duplicated(x[, agg_fact]) | duplicated(x[, agg_fact], fromLast = TRUE)

  x_tot <- cbind(x, es = x[, es], se = x[, se], row_index = 1:nrow(x))

  x_unique <- x_tot[!dup_row, ]

  x <- x_tot[dup_row, ]

  if (nrow(x) == length(unique(x[, agg_fact]))) {
    warning("the number of rows and unique 'agg_fact' values are equal")
    return(x_unique)
  }

  # Split on a character key: split() on a factor keeps an (empty) group for
  # every unused level, which would slip past the guards and inject phantom
  # all-NA rows into the output.
  x_split <- split(x, as.character(x[, agg_fact]))

  agg_bor <- do.call(rbind, lapply(x_split, .unique_es_subgroups, na.rm = na.rm))

  agg_bor$agg <- row.names(agg_bor)
  # agg_bor$agg <- unique(x[, agg_fact])

  if (any(!is.na(col_mean))) {
    df_mean <- aggregate.data.frame(x[, col_mean],
                                    by = list(agg = x[, agg_fact]),
                                    FUN = .mean_na, na.rm = na.rm)
    colnames(df_mean) <- c("agg", col_mean)
  }
  if (any(!is.na(col_weighted_mean))) {
    x_pond <- x[, col_weighted_mean, drop = FALSE] * x[, weights]

    sum_x_pond <- aggregate.data.frame(x_pond, by = list(agg = x[, agg_fact]),
                                       FUN = .sum_na, na.rm = na.rm)
    sum_N <- aggregate.data.frame(x[, weights], by = list(agg = x[, agg_fact]),
                                  FUN = .sum_na, na.rm = na.rm)

    # Carry the aggregate-produced (sorted) key alongside its own values -
    # unique(x[, agg_fact]) is in APPEARANCE order and would mis-pair labels
    # with values whenever the input is not pre-sorted. Indexing by name (not
    # $x) also supports multi-column col_weighted_mean.
    df_w_mean <- data.frame(agg = sum_x_pond$agg,
                            sum_x_pond[, col_weighted_mean, drop = FALSE] / sum_N$x,
                            check.names = FALSE)
    colnames(df_w_mean) <- c("agg", col_weighted_mean)
  }

  if (any(!is.na(col_sum))) {
    df_sum <- aggregate.data.frame(x[, col_sum], by = list(agg = x[, agg_fact]),
                                   FUN = .sum_na, na.rm = na.rm)
    colnames(df_sum) <- c("agg", col_sum)
  }

  if (any(!is.na(col_min))) {
    df_min <- aggregate.data.frame(x[, col_min], by = list(agg = x[, agg_fact]),
                                   FUN = .min_na, na.rm = na.rm)
    colnames(df_min) <- c("agg", col_min)
  }

  if (any(!is.na(col_max))) {
    df_max <- aggregate.data.frame(x[, col_max], by = list(agg = x[, agg_fact]),
                                   FUN = .max_na, na.rm = na.rm)
    colnames(df_max) <- c("agg", col_max)
  }

  if (any(!is.na(col_fact))) {
    df_fact <- aggregate.data.frame(x[, col_fact], by = list(agg = x[, agg_fact]),
                                    FUN = .agg.fact, na.rm = na.rm)
    colnames(df_fact) <- c("agg", col_fact)
  }

  # Merge only the frames that were actually built. The unused col_* placeholders
  # default to scalar NA; merging against them used to append an artifact column
  # named 'y' (later stripped with subset(-c(y))) - a landmine that collided with
  # and destroyed any real user column named 'y'.
  parts <- Filter(is.data.frame, list(agg_bor, df_mean, df_w_mean,
                                      df_sum, df_min, df_max, df_fact))
  df_add <- Reduce(function(a, b) merge(a, b), parts)

  first <- c(
    which(colnames(df_add) == "row_index"),
    which(colnames(df_add) == "agg")
  )

  df_add_clean <- cbind(df_add[, first], df_add[, -first])

  colnames(df_add_clean)[colnames(df_add_clean) == "agg"] <- agg_fact

  df_return <- rbind(x_unique[, colnames(df_add_clean)], df_add_clean)

  res <- df_return[order(df_return$row_index), ]
  rownames(res) <- 1:nrow(res)

  return(res)
}
.sum_na <- function(x, na.rm) {
  if (sum(is.na(x)) == length(x) & na.rm == TRUE) {
    NA
  } else {
    sum(x, na.rm = na.rm)
  }
}
.mean_na <- function(x, na.rm) {
  if (sum(is.na(x)) == length(x) & na.rm == TRUE) {
    NA
  } else {
    mean(x, na.rm = na.rm)
  }
}
.min_na <- function(x, na.rm) {
  if (sum(is.na(x)) == length(x) & na.rm == TRUE) {
    NA
  } else {
    min(x, na.rm = na.rm)
  }
}
.max_na <- function(x, na.rm) {
  if (sum(is.na(x)) == length(x) & na.rm == TRUE) {
    NA
  } else {
    max(x, na.rm = na.rm)
  }
}

.agg.outcomes <- function(x, agg_fact, es, se,
                          col_mean, col_weighted_mean, weights,
                          col_sum, col_min, col_max,
                          col_fact, cor_unit, na.rm) {
  cols <- c(col_mean, col_weighted_mean, col_sum, col_min, col_max, col_fact)

  dup_col <- ifelse(length(na.omit(cols)) > 0, any(duplicated(na.omit(cols))), FALSE)

  if (dup_col) {
    stop("Each column of your dataset can be aggregated only once.")
  }

  if (any(!is.na(col_weighted_mean) & is.na(weights))) {
    stop("You should specify the 'weights' argument when requesting to calculate a weighted mean for some columns.")
  }

  df_mean <- df_w_mean <- df_sum <- df_min <- df_max <- df_fact <- NA


  dup_row <- duplicated(x[, agg_fact]) | duplicated(x[, agg_fact], fromLast = TRUE)

  x_tot <- cbind(x,
    es = x[, es], se = x[, se],
    row_index = 1:nrow(x)
  )

  x_unique <- x_tot[!dup_row, ]

  x <- x_tot[dup_row, ]

  if (nrow(x) == length(unique(x[, agg_fact]))) {
    warning("the number of rows and unique 'agg_fact' values are equal")
    return(x_unique)
  }

  # Split on a character key: split() on a factor keeps an (empty) group for
  # every unused level, which would slip past the guards and inject phantom
  # all-NA rows into the output.
  x_split <- split(x, as.character(x[, agg_fact]))
  # measure = "SMD"
  agg_bor <- do.call(rbind, lapply(x_split, .unique_es_outcomes, cor_unit = cor_unit, na.rm = na.rm))

  agg_bor$agg <- row.names(agg_bor)

  if (any(!is.na(col_mean))) {
    df_mean <- aggregate.data.frame(x[, col_mean], by = list(agg = x[, agg_fact]), FUN = .mean_na, na.rm = na.rm)
    colnames(df_mean) <- c("agg", col_mean)
  }

  if (any(!is.na(col_weighted_mean))) {
    x_pond <- x[, col_weighted_mean, drop = FALSE] * x[, weights]
    sum_x_pond <- aggregate.data.frame(x_pond, by = list(agg = x[, agg_fact]), FUN = .sum_na, na.rm = na.rm)
    sum_N <- aggregate.data.frame(x[, weights], by = list(agg = x[, agg_fact]), FUN = .sum_na, na.rm = na.rm)
    df_w_mean <- data.frame(agg = sum_x_pond$agg,
                            sum_x_pond[, col_weighted_mean, drop = FALSE] / sum_N$x,
                            check.names = FALSE)
    colnames(df_w_mean) <- c("agg", col_weighted_mean)
  }

  if (any(!is.na(col_sum))) {
    df_sum <- aggregate.data.frame(x[, col_sum], by = list(agg = x[, agg_fact]),
                                   FUN = .sum_na, na.rm = na.rm)
    colnames(df_sum) <- c("agg", col_sum)
  }

  if (any(!is.na(col_min))) {
    df_min <- aggregate.data.frame(x[, col_min],
                                   by = list(agg = x[, agg_fact]),
                                   FUN = .min_na, na.rm = na.rm)
    colnames(df_min) <- c("agg", col_min)
  }

  if (any(!is.na(col_max))) {
    df_max <- aggregate.data.frame(x[, col_max],
                                   by = list(agg = x[, agg_fact]),
                                   FUN = .max_na, na.rm = na.rm)
    colnames(df_max) <- c("agg", col_max)
  }

  if (any(!is.na(col_fact))) {
    df_fact <- aggregate.data.frame(x[, col_fact],
                                    by = list(agg = x[, agg_fact]),
                                    FUN = .agg.fact, na.rm = na.rm)
    colnames(df_fact) <- c("agg", col_fact)
  }

  # Merge only the frames that were actually built (see .agg.subgroups for the
  # rationale - avoids the 'y' artifact-column collision).
  parts <- Filter(is.data.frame, list(agg_bor, df_mean, df_w_mean,
                                      df_sum, df_min, df_max, df_fact))
  df_add <- Reduce(function(a, b) merge(a, b), parts)

  first <- c(
    which(colnames(df_add) == "row_index"),
    which(colnames(df_add) == "agg")
  )

  df_add_clean <- cbind(df_add[, first], df_add[, -first])

  colnames(df_add_clean)[colnames(df_add_clean) == "agg"] <- agg_fact

  df_return <- rbind(x_unique[, colnames(df_add_clean)], df_add_clean)

  res <- df_return[order(df_return$row_index), ]
  return(res)
}
.agg.times <- function(x, agg_fact, es, se,
                          col_mean, col_weighted_mean, weights,
                          col_sum, col_min, col_max,
                          col_fact, cor_unit, na.rm) {
  cols <- c(col_mean, col_weighted_mean, col_sum, col_min, col_max, col_fact)

  dup_col <- ifelse(length(na.omit(cols)) > 0, any(duplicated(na.omit(cols))), FALSE)

  if (dup_col) {
    stop("Each column of your dataset can be aggregated only once.")
  }

  if (any(!is.na(col_weighted_mean) & is.na(weights))) {
    stop("You should specify the 'weights' argument when requesting to calculate a weighted mean for some columns.")
  }

  df_mean <- df_w_mean <- df_sum <- df_min <- df_max <- df_fact <- NA


  dup_row <- duplicated(x[, agg_fact]) | duplicated(x[, agg_fact], fromLast = TRUE)

  x_tot <- cbind(x,
                 es = x[, es], se = x[, se],
                 row_index = 1:nrow(x)
  )

  x_unique <- x_tot[!dup_row, ]

  x <- x_tot[dup_row, ]

  if (nrow(x) == length(unique(x[, agg_fact]))) {
    warning("the number of rows and unique 'agg_fact' values are equal")
    return(x_unique)
  }

  # dependence = "times" needs a per-row time index and a within-cluster
  # correlation. Require the time column, and fall back to the cor_unit ARGUMENT
  # when the data carries no cor_unit column (a cor_unit column takes precedence
  # when present).
  if (!"time_agg" %in% colnames(x)) {
    stop("dependence = 'times' requires a 'time_agg' column giving the time-point (wave) of each effect size.")
  }
  if (!"cor_unit" %in% colnames(x)) {
    x$cor_unit <- cor_unit
  }

  # Split on a character key: split() on a factor keeps an (empty) group for
  # every unused level, which would slip past the guards and inject phantom
  # all-NA rows into the output.
  x_split <- split(x, as.character(x[, agg_fact]))
  agg_bor <- do.call(rbind, lapply(x_split, .unique_es_times, na.rm = na.rm))

  agg_bor$agg <- row.names(agg_bor)

  if (any(!is.na(col_mean))) {
    df_mean <- aggregate.data.frame(x[, col_mean], by = list(agg = x[, agg_fact]), FUN = .mean_na, na.rm = na.rm)
    colnames(df_mean) <- c("agg", col_mean)
  }

  if (any(!is.na(col_weighted_mean))) {
    x_pond <- x[, col_weighted_mean, drop = FALSE] * x[, weights]
    sum_x_pond <- aggregate.data.frame(x_pond, by = list(agg = x[, agg_fact]), FUN = .sum_na, na.rm = na.rm)
    sum_N <- aggregate.data.frame(x[, weights], by = list(agg = x[, agg_fact]), FUN = .sum_na, na.rm = na.rm)
    df_w_mean <- data.frame(agg = sum_x_pond$agg,
                            sum_x_pond[, col_weighted_mean, drop = FALSE] / sum_N$x,
                            check.names = FALSE)
    colnames(df_w_mean) <- c("agg", col_weighted_mean)
  }

  if (any(!is.na(col_sum))) {
    df_sum <- aggregate.data.frame(x[, col_sum], by = list(agg = x[, agg_fact]),
                                   FUN = .sum_na, na.rm = na.rm)
    colnames(df_sum) <- c("agg", col_sum)
  }

  if (any(!is.na(col_min))) {
    df_min <- aggregate.data.frame(x[, col_min],
                                   by = list(agg = x[, agg_fact]),
                                   FUN = .min_na, na.rm = na.rm)
    colnames(df_min) <- c("agg", col_min)
  }

  if (any(!is.na(col_max))) {
    df_max <- aggregate.data.frame(x[, col_max],
                                   by = list(agg = x[, agg_fact]),
                                   FUN = .max_na, na.rm = na.rm)
    colnames(df_max) <- c("agg", col_max)
  }

  if (any(!is.na(col_fact))) {
    df_fact <- aggregate.data.frame(x[, col_fact],
                                    by = list(agg = x[, agg_fact]),
                                    FUN = .agg.fact, na.rm = na.rm)
    colnames(df_fact) <- c("agg", col_fact)
  }

  # Merge only the frames that were actually built (see .agg.subgroups for the
  # rationale - avoids the 'y' artifact-column collision).
  parts <- Filter(is.data.frame, list(agg_bor, df_mean, df_w_mean,
                                      df_sum, df_min, df_max, df_fact))
  df_add <- Reduce(function(a, b) merge(a, b), parts)

  first <- c(
    which(colnames(df_add) == "row_index"),
    which(colnames(df_add) == "agg")
  )

  df_add_clean <- cbind(df_add[, first], df_add[, -first])

  colnames(df_add_clean)[colnames(df_add_clean) == "agg"] <- agg_fact

  df_return <- rbind(x_unique[, colnames(df_add_clean)], df_add_clean)

  res <- df_return[order(df_return$row_index), ]
  return(res)
}

.unique_es_subgroups <- function(x, na.rm = TRUE) {
  row_index <- x$row_index[1]

  # Drop effect sizes with a missing es/se before combining, matching the
  # behaviour of metafor::aggregate.escalc() (which discards NA rows). Without
  # this, a single NA nullifies the whole cluster and can leave an incoherent
  # (es = NA, se = <valid>) pair.
  if (na.rm) {
    x <- x[!is.na(x$es) & !is.na(x$se), , drop = FALSE]
  }
  # A non-positive or non-finite SE cannot enter inverse-variance weighting: 1/se^2
  # becomes Inf/NaN, which turns the combined estimate into Inf/Inf = NaN paired with
  # se = 0 (an incoherent result). Drop such rows -- se = 0 is itself pathological and
  # is flagged upstream by the Tier-2 'SE is zero' check -- and warn so it is visible.
  bad_se <- !is.finite(x$se) | x$se <= 0
  if (any(bad_se)) {
    warning(sprintf(
      "aggregate_df: dropped %d effect size(s) with a non-positive/non-finite SE from a subgroup cluster before inverse-variance combination.",
      sum(bad_se)))
    x <- x[!bad_se, , drop = FALSE]
  }
  if (nrow(x) == 0) {
    return(data.frame(es = NA_real_, se = NA_real_, row_index = row_index))
  }

  weights <- 1 / (x$se^2)
  mean_es <- sum(weights * x$es) / sum(weights)
  se <- sqrt(1 / sum(weights))

  data.frame(es = mean_es, se = se, row_index = row_index)
}
.unique_es_outcomes <- function(x, cor_unit, na.rm = TRUE) {
  row_index <- x$row_index[1]

  if (na.rm) {
    x <- x[!is.na(x$es) & !is.na(x$se), , drop = FALSE]
  }
  if (nrow(x) == 0) {
    return(data.frame(es = NA_real_, se = NA_real_, row_index = row_index))
  }

  var_es <- x$se^2
  prod_se <- x$se %*% t(x$se)
  prod_se_r <- prod_se * cor_unit
  prod_se_r[lower.tri(prod_se_r)] <- 0
  diag(prod_se_r) <- 0

  mean_es <- mean(x$es)
  var <- (1 / length(x$es))^2 * (sum(var_es) + 2 * sum(prod_se_r))
  se <- sqrt(var)

  data.frame(es = mean_es, se = se, row_index = row_index)
}


.unique_es_times <- function(x, na.rm = TRUE) {
  row_index <- x[1, "row_index"]

  # NB: length(unique(.)) - a *count* of distinct values - not unique(.) > 1,
  # which tests the value and errors ("condition has length > 1") when a cluster
  # legitimately carries more than one cor_unit value.
  if (length(unique(x[, "cor_unit"])) > 1) {
    stop("cor_unit values should be constant within a cluster unit when 'dependence=times' is called")
  }
  cor_unit <- unique(x[, "cor_unit"])

  if (na.rm) {
    x <- x[!is.na(x[, "es"]) & !is.na(x[, "se"]), , drop = FALSE]
  }
  if (nrow(x) == 0) {
    return(data.frame(es = NA_real_, se = NA_real_, row_index = row_index))
  }

  time <- x[, "time_agg"]
  R <- outer(time, time, function(a, b) cor_unit^abs(a - b))

  # nrow= guards the length-1 case: diag(scalar) would build an identity matrix.
  S <- diag(sqrt(x[, "se"]^2), nrow = nrow(x))
  V <- S %*% R %*% S

  W <- try(solve(V), silent = TRUE)
  if (inherits(W, "try-error")) {
    stop("Cannot take inverse of the covariance matrix 'V'.")
  }

  sum_W <- sum(W)
  mean_es <- sum(W %*% x[, "es"]) / sum_W
  se <- sqrt(1 / sum_W)

  data.frame(es = mean_es, se = se, row_index = row_index)
}

.agg.fact <- function(x, na.rm = na.rm) {
  if (na.rm) {
    ifelse(length(unique(na.omit(x))) == 0,
      "N/A",
      ifelse(length(unique(na.omit(x))) == 1,
        unique(na.omit(x)),
        paste(na.omit(sort(unique(x))), collapse = " / ")
      )
    )
  } else {
    ifelse(length(unique(x)) == 0,
      "N/A",
      ifelse(length(unique(x)) == 1,
        unique(x),
        paste(sort(unique(x)), collapse = " / ")
      )
    )
  }
}
