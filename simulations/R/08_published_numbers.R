## =============================================================================
## NUMBERS PUBLISHED IN README.md, TIED BACK TO THE CSV THEY CAME FROM
##
## WHY THIS EXISTS. data/aggregated/PROVENANCE.csv answers "was this CSV produced by
## the code now in the tree?". It says nothing about the numbers COPIED OUT of that CSV
## into prose -- and that is where the rot actually happened:
##
##   * README's study-09 table published bonett at 0.0185 / 0.2571 / 1.051 while the
##     shipped aggregate said 0.0197 / 0.2868 / 1.053;
##   * it published 2x2_tetrachoric's coverage as 0.928 (worst 0.849) where the
##     aggregate said 0.961 (worst 0.939) -- and THAT figure carried a claim reused
##     further down, "the best method on bias is the WORST on coverage", which was
##     false on the data as shipped.
##
## Both had been wrong since the item 3.5 regeneration. Nothing noticed, because
## nothing connected the paragraph to the file.
##
## THE DESIGN, and why it is not a list of expected numbers. A hand-written table of
## "the value should be 0.0195" would rot exactly as the README did -- it is the same
## copied number, moved. So the EXPECTED SIDE IS DERIVED: each entry below is a recipe
## that recomputes the table from the aggregate, and the test compares it against the
## markdown table parsed out of README.md. Neither side is a transcription.
##
## Adding a table is one function here plus one `<!-- pinned: id -->` marker above the
## table in README.md. HTML comments do not render, so the document is unaffected.
##
## KNOWN LIMIT, stated rather than hidden: tables computed from the per-replication
## data/raw/*.rds (e.g. the RMSE-vs-EmpSE illustration in the "own target" section,
## which needs cor(est, target)) cannot be pinned this way, because data/raw/ is
## gitignored and regenerable rather than shipped. Only aggregate-derived tables are
## reachable.
## =============================================================================

#' Recompute one published README table from its aggregate
#'
#' @param id marker id, as written in the README's `<!-- pinned: id -->` comment
#' @return a data.frame whose columns are the table's columns, in published order and
#'   at published precision, and whose first column is the row key
published_table <- function(id) {
  f <- switch(id,
    "study09a-result"  = .pub_09a_result,
    "study09a-ranking" = .pub_09a_ranking,
    stop("no recipe for pinned table '", id, "'"))
  f()
}


#' Every id this module can recompute
published_table_ids <- function() c("study09a-result", "study09a-ranking")


## The (r)-scale rows of 09a scored on each method's own estimand -- "computational
## error, not estimand mismatch", which is what the README section says it reports.
.pub_09a_slice <- function() {
  d <- utils::read.csv(dir_agg("09a_or_to_cor_CONT_nrep1000.csv"), stringsAsFactors = FALSE)
  d <- d[d$target == "own" & grepl("[(]r[)]", d$method), , drop = FALSE]
  d$m <- sub(" [(]r[)]", "", d$method)
  d
}


.pub_09a_result <- function() {
  d <- .pub_09a_slice()
  o <- do.call(rbind, lapply(split(d, d$m), function(x) data.frame(
    method   = x$m[1],
    mean_abs = round(mean(abs(x$bias), na.rm = TRUE), 4),
    worst    = round(max(abs(x$bias), na.rm = TRUE), 4),
    coverage = round(mean(x$coverage, na.rm = TRUE), 3),
    se_ratio = round(mean(x$se_ratio, na.rm = TRUE), 3),
    stringsAsFactors = FALSE)))
  o[order(o$mean_abs), ]
}


## The ranking illustration: bias rank against coverage rank, the latter on
## |coverage - 0.95| because a coverage is judged by distance from nominal, not by size.
## Only the three rows the README prints.
.pub_09a_ranking <- function() {
  d <- .pub_09a_slice()
  o <- do.call(rbind, lapply(split(d, d$m), function(x) data.frame(
    method   = x$m[1],
    mean_abs = mean(abs(x$bias), na.rm = TRUE),
    coverage = mean(x$coverage, na.rm = TRUE),
    worst_cov = min(x$coverage, na.rm = TRUE),
    stringsAsFactors = FALSE)))
  o$rank_bias <- rank(o$mean_abs)
  o$rank_cov  <- rank(abs(o$coverage - 0.95))
  shown <- c("2x2_tetrachoric", "bonett", "lipsey_cooper")
  o <- o[match(shown, o$method), ]
  data.frame(
    method    = o$method,
    rank_bias = .pub_ordinal(o$rank_bias),
    coverage  = ifelse(o$method == "2x2_tetrachoric",
                       sprintf("%.3f (worst %.3f)", o$coverage, o$worst_cov),
                       sprintf("%.3f", o$coverage)),
    rank_cov  = .pub_ordinal(o$rank_cov),
    stringsAsFactors = FALSE)
}


.pub_ordinal <- function(k) {
  k <- as.integer(round(k))
  suffix <- ifelse(k %% 100 %in% 11:13, "th",
                   ifelse(k %% 10 == 1, "st",
                          ifelse(k %% 10 == 2, "nd",
                                 ifelse(k %% 10 == 3, "rd", "th"))))
  paste0(k, suffix)
}


#' Parse the markdown table that follows a `<!-- pinned: id -->` marker in README.md
#'
#' @param id marker id
#' @param path README to read
#' @return a character data.frame of the table's cells, header names kept verbatim
read_pinned_table <- function(id, path = sim_path("README.md")) {
  l <- readLines(path, warn = FALSE)
  i <- grep(paste0("<!--\\s*pinned:\\s*", id, "\\s*-->"), l)
  if (length(i) != 1L)
    stop("marker '", id, "' found ", length(i), " times in ", basename(path))
  ## the table starts at the next line beginning with "|" and runs while lines do
  j <- i + 1L
  while (j <= length(l) && !startsWith(trimws(l[j]), "|")) j <- j + 1L
  k <- j
  while (k <= length(l) && startsWith(trimws(l[k]), "|")) k <- k + 1L
  rows <- l[j:(k - 1L)]
  rows <- rows[!grepl("^\\s*\\|[\\s:|-]*\\|\\s*$", rows, perl = TRUE)]  # drop the ---- rule
  ## Split on UNESCAPED pipes only. These headers contain `mean \|bias\|`, where the
  ## pipes are markdown-escaped and belong to the cell rather than separating cells.
  ## Splitting on them produced rows of unequal length and a cryptic recycling error
  ## ("number of columns of result is not a multiple of vector length") rather than
  ## anything pointing at the cause.
  split_row <- function(r) {
    r <- sub("^\\s*\\|", "", r)
    r <- sub("(?<!\\\\)\\|\\s*$", "", r, perl = TRUE)
    cells <- strsplit(r, "(?<!\\\\)\\|", perl = TRUE)[[1]]
    trimws(gsub("\\\\\\|", "|", cells))
  }
  m <- do.call(rbind, lapply(rows, split_row))
  out <- as.data.frame(m[-1, , drop = FALSE], stringsAsFactors = FALSE)
  names(out) <- m[1, ]
  out
}


#' Strip the markdown a published cell may carry, leaving the value
#'
#' Bold (`**0.0195**`), code (`` `bonett` ``) and the minus sign U+2212 all appear in
#' these tables; none of them is part of the number or the name.
.pub_clean <- function(x) {
  x <- gsub("\\*\\*", "", x)
  x <- gsub("`", "", x)
  x <- gsub("−", "-", x)
  trimws(x)
}
