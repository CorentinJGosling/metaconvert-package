#' Flag the differences between two dataframes.
#'
#' @param df_extractor_1 a first dataset. Differences with the second dataset will be flagged in green.
#' @param df_extractor_2 a second dataset. Differences with the first dataset will be flagged in red.
#' @param tolerance the cut-off value used to flag differences between two numeric values
#' @param tolerance_type must be either 'difference' or 'ratio'
#' @param output type of object returned by the function (see 'Value' section). Must be either 'wide', 'long', 'html', 'html2' or 'xlsx'.
#' @param file_name the name of the generated file (only used when \code{output="xlsx"})
#' @param ordering_columns column names that should be used to re-order the two datasets before running the comparisons
#'
#' @details
#' This function aims to facilitate the comparison of two datasets created by blind data extractors during a systematic review.
#' It is a wrapper of several functions from the 'compareDF' package.
#'
#' @return
#' This function returns a dataframe composed of the rows that include a
#' difference (all identical rows are removed).
#' Several outputs can be requested :
#' 1. setting \code{output="xlsx"} returns an excel file. A message indicates the location of the generated file on your computer.
#' 1. setting \code{output="html"} returns an html file
#' 1. setting \code{output="html2"} returns an html file (only useful when the "html" command did not make the html pane appear in R studio).
#' 1. setting \code{output="wide"} a wide dataframe
#' 1. setting \code{output="long"} a long dataframe
#'
#' @md
#'
#' @export compare_df
#'
#' @references
#' Alex Joseph (2022). compareDF: Do a Git Style Diff of the Rows Between Two Dataframes with Similar Structure. R package version 2.3.3. https://CRAN.R-project.org/package=compareDF
#'
#' @examples
#' df.compare1 = df.compare1[order(df.compare1$author), ]
#' df.compare2 = df.compare2[order(df.compare2$year), ]
#'
#' compare_df(
#'   df_extractor_1 = df.compare1,
#'   df_extractor_2 = df.compare2,
#'   ordering_columns = c("author", "year")
#' )

compare_df <- function(df_extractor_1, df_extractor_2,
                       ordering_columns = NULL,
                       tolerance = 0,
                       tolerance_type = "ratio",
                       output = "html",
                       file_name = "comparison.xlsx") {
  # An unrecognised value used to fall through every branch and fail later with
  # "object 'res' not found", which says nothing about the argument that caused it.
  output <- match.arg(output, c("wide", "long", "html", "html2", "xlsx"))

  df_extractor_1 = data.frame(df_extractor_1)
  df_extractor_2 = data.frame(df_extractor_2)


  cols_available = unique(append(names(df_extractor_1), names(df_extractor_2)))
  cols_retained = intersect(names(df_extractor_1), names(df_extractor_2))
  cols_unique = cols_available[!cols_available %in% cols_retained]
  if (length(cols_unique) > 0) {
    warning("\nThe columns ", paste(cols_unique, collapse = " / "),
            " were removed from the comparison since they were included only in one of the two datasets.\n")
  }

  if (!is.null(ordering_columns)) {
    if (!all(ordering_columns %in% cols_retained)) {
      stop("The ordering columns should be present in the two datasets")
    } else if (length(ordering_columns) > 1) {

      create_unique_id <- function(df, columns) {
        do.call(paste, c(df[, columns, drop = FALSE], sep = "_"))
      }

      df_extractor_1$ID_metaConvert <- create_unique_id(df_extractor_1, ordering_columns)
      df_extractor_2$ID_metaConvert <- create_unique_id(df_extractor_2, ordering_columns)

      ordering_columns = "ID_metaConvert"
      cols_retained = append(cols_retained, ordering_columns)

    }
  }



  if (!is.null(ordering_columns)) {
    # The alignment below matches each dataset's rows to a shared, DE-DUPLICATED list of
    # key values and assigns them one at a time. A key that repeats therefore selects
    # several rows for a single destination: base R assigns the first, warns
    # "replacement element 1 has 2 rows to replace 1 row", and the other rows are
    # silently dropped from the comparison. The warning is opaque and the data loss is
    # invisible, so the ambiguity is refused here instead. ordering_columns exists to
    # line the two sheets up row for row, which a non-unique key cannot do.
    for (.d in list(df_extractor_1, df_extractor_2)) {
      .k <- if (length(ordering_columns) > 1) {
        do.call(paste, c(.d[, ordering_columns, drop = FALSE], sep = "_"))
      } else .d[, ordering_columns]
      .dup <- unique(.k[duplicated(.k)])
      if (length(.dup)) {
        stop("'ordering_columns' must identify each row uniquely, but ",
             paste0("'", utils::head(.dup, 5), "'", collapse = ", "),
             if (length(.dup) > 5) paste0(" (and ", length(.dup) - 5L, " more)") else "",
             " occur more than once. Add a column that distinguishes them, or drop ",
             "ordering_columns to compare the sheets in their given row order.",
             call. = FALSE)
      }
    }
    df_extractor_1 = suppressWarnings(df_extractor_1[order(df_extractor_1[, ordering_columns]), cols_retained])
    df_extractor_2 = suppressWarnings(df_extractor_2[order(df_extractor_2[, ordering_columns]), cols_retained])

    combined_values <- unique(c(df_extractor_1[, ordering_columns], df_extractor_2[, ordering_columns]))
    combined_values <- combined_values[order(combined_values)]

    insert_blank_rows <- function(df, combined_values, ordering_columns) {
      new_df <- data.frame(matrix(ncol = ncol(df), nrow = length(combined_values)))
      colnames(new_df) <- colnames(df)

      for (i in seq_along(combined_values)) {
        match_row <- which(df[, ordering_columns] == combined_values[i])
        if (length(match_row) > 0) {
          new_df[i, ] <- df[match_row, ]
        }
      }

      return(new_df)
    }

    # Apply the function to both dataframes
    df_extractor_1 <- insert_blank_rows(df_extractor_1, combined_values, ordering_columns)
    df_extractor_2 <- insert_blank_rows(df_extractor_2, combined_values, ordering_columns)
    df_extractor_1$rowname <- 1:nrow(df_extractor_1)
    df_extractor_2$rowname <- 1:nrow(df_extractor_2)

    } else {
    df_extractor_1 = df_extractor_1[, cols_retained]
    df_extractor_2 = df_extractor_2[, cols_retained]
    df_extractor_1$rowname <- 1:nrow(df_extractor_1)
    df_extractor_2$rowname <- 1:nrow(df_extractor_2)
    }


  comp <- suppressMessages(compareDF::compare_df(
    df_new = df_extractor_1,
    df_old = df_extractor_2,
    group_col = "rowname",
    exclude = NULL,
    tolerance = tolerance,
    tolerance_type = tolerance_type,
    stop_on_error = TRUE,
    change_markers = c("df_extractor_1", "df_extractor_2", "="),
    round_output_to = 3
  ))

  if (output == "long") {
    res <- comp$comparison_df[order(as.numeric(as.character(comp$comparison_df$rowname))), ]
  } else if (output == "wide") {
    res <- compareDF::create_wide_output(comp)
    res <- res[order(as.numeric(as.character(res$rowname))), ]
  } else if (output == "html") {
    res <- suppressMessages(compareDF::create_output_table(comp, output_type = "html", limit = 5000))
  } else if (output == "html2") {
    # view_html() drives the RStudio viewer pane; outside RStudio getOption("viewer")
    # is NULL and calling it raises "attempt to apply non-function".
    if (is.null(getOption("viewer"))) {
      stop("output = 'html2' needs the RStudio viewer pane; use output = 'html' instead.",
           call. = FALSE)
    }
    res <- suppressMessages(compareDF::view_html(comp))
  } else if (output == "xlsx") {
    if (!grepl(".xlsx", file_name)) {
      file_name <- paste0(file_name, ".xlsx")
    }
    message("File generated at '", paste0(getwd(), "/", file_name, "'"))
    res <- compareDF::create_output_table(comp,
      output_type = "xlsx",
      file_name = file_name, limit = 5000
    )
  }
  return(res)
}
