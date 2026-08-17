# =============================================================================
# Roadmap items 2.3 and 2.4 -- an argument's contract must match its behaviour.
#
# 2.3  THE SAME TABLE CONVERTED DIFFERENTLY DEPENDING ON WHICH DOOR YOU CAME IN.
#      convert_df() and four of five OR entry points defaulted to
#      or_to_cor = "bonett"; es_from_or_se(), es_from_user_crude() and
#      es_from_user_adj() shipped "pearson". A formals scan over every exported
#      function found 16 divergences in total, of which 13 are legitimate (below).
#
# 2.4  unit_type is read in exactly ONE place, as `unit_type == "sd"`, so every
#      other value -- a typo included -- silently selected raw units. The shipped
#      default "raw_scale" was not even among the documented values: six @param
#      blocks said "sd" or "raw_scale", two said "sd" or "value".
#
# WHY THIS FILE IS A SCAN AND NOT A LIST. Nineteen arguments are shared between
# convert_df() and the es_from_*() family. Checking them by hand is how they
# drifted. The scan below fails on ANY new divergence, so the next one is caught
# when it is introduced rather than by a later audit.
#
# THE 13 JUSTIFIED DIVERGENCES ARE DERIVED, NOT HARDCODED. The mean-change and
# paired routes default to pre_post_to_smd = "cooper" against convert_df()'s
# "bonett" -- correctly: bonett needs a baseline SD that change-score and paired
# test-statistic data do not carry, and those functions ALREADY refuse it through
# .validate_pre_post_to_smd(allowed_methods = ...). So a divergence is legitimate
# exactly when the function's own allowed_methods excludes convert_df()'s default.
# That rule is self-maintaining: a hand-written allowlist would let genuine drift
# hide inside it, and a blanket "all defaults must match" would have demanded 13
# changes that break every direct call to those routes.
# =============================================================================

NS <- asNamespace("metaConvert")

# An argument with no default is the empty symbol. ANY function that receives it
# as an argument forces it and raises "argument is missing", so it has to be
# detected structurally, by comparing sub-lists.
.arg_defaults <- function(f) {
  fm <- formals(f)
  if (!length(fm)) return(list())
  keep <- vapply(seq_along(fm), function(i)
    !identical(unname(as.list(fm)[i]), unname(alist(q = ))), logical(1))
  lapply(as.list(fm)[keep], function(v) tryCatch(eval(v), error = function(e) NULL))
}

.exported_functions <- function() {
  sort(Filter(function(n) is.function(get(n, envir = NS)),
              getNamespaceExports("metaConvert")))
}

# Walk a function body for `.validate_pre_post_to_smd(allowed_methods = c(...))`
# and return that set.
#
# ONE LEVEL OF DELEGATION IS FOLLOWED. es_from_paired_f(), es_from_paired_t_pval()
# and es_from_paired_f_pval() never call the validator themselves -- they convert
# their statistic and hand off to es_from_paired_t(), which does. Measured: all
# three still reject pre_post_to_smd = "bonett", but the error names the DELEGATE,
# so a structural check that stopped at the function itself would call three
# legitimate defaults "unjustified" and demand they be broken.
#
# NULL when nothing is found even after delegation, or when allowed_methods is not
# a literal -- both count as "not justified", the safe direction for a guard.
.allowed_pre_post_local <- function(f) {
  found <- NULL
  walk <- function(e) {
    if (is.call(e)) {
      fn <- e[[1]]
      if (is.name(fn) && identical(as.character(fn), ".validate_pre_post_to_smd")) {
        al <- as.list(e)
        am <- al$allowed_methods
        if (is.null(am) && length(al) >= 3) am <- al[[3]]
        v <- tryCatch(eval(am), error = function(err) NULL)
        if (is.character(v)) found <<- unique(c(found, v))
      }
      for (i in seq_along(e)) {
        el <- tryCatch(e[[i]], error = function(err) NULL)
        if (!is.null(el)) walk(el)
      }
    }
  }
  tryCatch(walk(body(f)), error = function(err) NULL)
  found
}

.called_package_functions <- function(f) {
  out <- character(0)
  walk <- function(e) {
    if (is.call(e)) {
      fn <- e[[1]]
      if (is.name(fn)) {
        nm <- as.character(fn)
        if (exists(nm, envir = NS, inherits = FALSE) &&
            is.function(get(nm, envir = NS))) out <<- c(out, nm)
      }
      for (i in seq_along(e)) {
        el <- tryCatch(e[[i]], error = function(err) NULL)
        if (!is.null(el)) walk(el)
      }
    }
  }
  tryCatch(walk(body(f)), error = function(err) NULL)
  unique(out)
}

## Delegation is followed TRANSITIVELY, not one level: es_from_paired_f_pval()
## hands off to es_from_paired_t_pval(), which hands off to es_from_paired_t(),
## which is where the validator finally lives. A one-level search called that
## function's default drift. Depth is bounded and visited nodes are tracked, so a
## cycle in the call graph cannot hang the suite.
.allowed_pre_post <- function(f, depth = 4L, .seen = character(0)) {
  local <- .allowed_pre_post_local(f)
  if (!is.null(local)) return(local)
  if (depth <= 0L) return(NULL)
  found <- NULL
  for (callee in setdiff(.called_package_functions(f), .seen)) {
    a <- .allowed_pre_post(get(callee, envir = NS), depth - 1L, c(.seen, callee))
    if (!is.null(a)) found <- unique(c(found, a))
  }
  found
}

.convert_df_defaults <- .arg_defaults(convert_df)
.shared_method_args <- {
  a <- names(.convert_df_defaults)[vapply(.convert_df_defaults, function(v)
    is.character(v) && length(v) >= 1 && nzchar(v[1]), logical(1))]
  setdiff(a, c("measure", "format"))
}

.scan_divergences <- function() {
  rows <- list()
  for (fn in .exported_functions()) {
    fd <- .arg_defaults(get(fn, envir = NS))
    for (a in intersect(names(fd), .shared_method_args)) {
      got <- fd[[a]]
      if (is.null(got) || !length(got)) next
      want <- .convert_df_defaults[[a]][1]
      if (!identical(got[1], want))
        rows[[length(rows) + 1L]] <- data.frame(fn = fn, arg = a,
                                                its_default = got[1], convert_df = want,
                                                stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) do.call(rbind, rows) else
    data.frame(fn = character(0), arg = character(0),
               its_default = character(0), convert_df = character(0))
}

# --- 2.3: the scan ----------------------------------------------------------

test_that("the shared method arguments are actually shared", {
  # If this shrinks, the scan below silently stops covering something.
  expect_gte(length(.shared_method_args), 15L)
  for (a in c("or_to_cor", "or_to_rr", "rr_to_or", "smd_to_cor", "cor_to_smd",
              "pre_post_to_smd", "smd_var", "smd_denom", "unit_type",
              "prop_to_es", "alpha_to_es", "icc_to_es", "table_2x2_to_cor"))
    expect_true(a %in% .shared_method_args, info = a)
})

test_that("every default divergence from convert_df() is a justified one", {
  d <- .scan_divergences()
  unjustified <- d[!vapply(seq_len(nrow(d)), function(i) {
    if (d$arg[i] != "pre_post_to_smd") return(FALSE)
    allowed <- .allowed_pre_post(get(d$fn[i], envir = NS))
    !is.null(allowed) && !(d$convert_df[i] %in% allowed)
  }, logical(1)), , drop = FALSE]

  expect_equal(nrow(unjustified), 0L,
               info = paste0("undocumented default drift: ",
                             paste(sprintf("%s(%s = '%s') vs convert_df '%s'",
                                           unjustified$fn, unjustified$arg,
                                           unjustified$its_default,
                                           unjustified$convert_df),
                                   collapse = "; ")))
})

test_that("the justified divergences are exactly the pre_post_to_smd family", {
  d <- .scan_divergences()
  expect_true(all(d$arg == "pre_post_to_smd"))
  expect_true(all(d$its_default == "cooper"))
  # 13 of them at the time of writing: mean-change and paired routes
  expect_gte(nrow(d), 13L)
  expect_true(all(grepl("mean_change|paired", d$fn)))
})

test_that("each justified divergence really does refuse convert_df()'s default", {
  # The rule that licenses the allowlist, asserted rather than trusted -- so the
  # allowlist cannot become a place for genuine drift to hide.
  d <- .scan_divergences()
  for (i in seq_len(nrow(d))) {
    allowed <- .allowed_pre_post(get(d$fn[i], envir = NS))
    expect_false(is.null(allowed),
                 info = paste(d$fn[i], "must declare allowed_methods"))
    expect_false("bonett" %in% allowed,
                 info = paste(d$fn[i], "excludes bonett, so 'cooper' is its only sane default"))
    # and its own default must be inside its own allowed set ("cooper" normalises
    # to "morris_drm")
    norm <- ifelse(d$its_default[i] == "cooper", "morris_drm", d$its_default[i])
    expect_true(norm %in% allowed, info = d$fn[i])
  }
})

test_that("a route that excludes bonett errors when asked for it", {
  # End-to-end confirmation that allowed_methods is not decorative.
  expect_error(
    es_from_mean_change_sd(mean_change_exp = 2, mean_change_sd_exp = 3,
                           mean_change_nexp = 1, mean_change_sd_nexp = 3,
                           n_exp = 30, n_nexp = 30,
                           r_pre_post_exp = .8, r_pre_post_nexp = .8,
                           pre_post_to_smd = "bonett"),
    "bonett")
})

test_that("the three DELEGATING routes refuse bonett too, via es_from_paired_t()", {
  # es_from_paired_f(), es_from_paired_t_pval() and es_from_paired_f_pval() declare
  # no allowed_methods of their own; they convert their statistic and hand off.
  # Structural inspection alone would mark their "cooper" default unjustified, so
  # the delegation is pinned behaviourally as well.
  base <- list(n_exp = 30, n_nexp = 30, r_pre_post_exp = .8, r_pre_post_nexp = .8)
  expect_error(do.call(es_from_paired_f,
                       c(list(paired_f_exp = 9, paired_f_nexp = 4,
                              pre_post_to_smd = "bonett"), base)), "bonett")
  expect_error(do.call(es_from_paired_t_pval,
                       c(list(paired_t_pval_exp = .01, paired_t_pval_nexp = .04,
                              pre_post_to_smd = "bonett"), base)), "bonett")
  expect_error(do.call(es_from_paired_f_pval,
                       c(list(paired_f_pval_exp = .01, paired_f_pval_nexp = .04,
                              pre_post_to_smd = "bonett"), base)), "bonett")
  # ... and their own default works
  expect_false(is.na(do.call(es_from_paired_f,
                             c(list(paired_f_exp = 9, paired_f_nexp = 4), base))$g))
})

test_that("or_to_cor no longer depends on which entry point you use", {
  # The 2.3 fix. These three shipped "pearson" while convert_df() and the other
  # OR routes shipped "bonett".
  for (fn in c("es_from_or_se", "es_from_user_crude", "es_from_user_adj")) {
    got <- .arg_defaults(get(fn, envir = NS))$or_to_cor
    expect_equal(got, .convert_df_defaults$or_to_cor,
                 info = paste(fn, "must match convert_df()"))
  }
  # and no exported function is left on the old default
  d <- .scan_divergences()
  expect_equal(sum(d$arg == "or_to_cor"), 0L)
})

test_that("aligning or_to_cor makes the entry points agree with convert_df()", {
  # Behaviour, not just formals: same row, both doors, same answer.
  row <- data.frame(or = 2.5, logor_se = 0.3, n_sample = 200, n_cases = 60,
                    n_controls = 140, n_exp = 100, n_nexp = 100)
  direct <- suppressMessages(
    es_from_or_se(or = row$or, logor_se = row$logor_se, n_sample = row$n_sample,
                  n_cases = row$n_cases, n_controls = row$n_controls,
                  n_exp = row$n_exp, n_nexp = row$n_nexp))
  via_df <- suppressMessages(summary(suppressMessages(convert_df(row, measure = "r"))))
  # summary() rounds es_crude for display (0.32 against the direct 0.320469), so the
  # comparison is made at the precision summary() actually publishes.
  expect_equal(round(direct$r, 2), via_df$es_crude, tolerance = 1e-8)
})

test_that("the documented consequence: the minimal input now yields NA, both ways", {
  # bonett cannot be computed from (or, logor_se) alone, and the lipsey_cooper
  # fallback needs n_sample, so the minimal row has no correlation by the default
  # method. That is a REAL change for direct callers -- es_from_or_se() used to
  # return 0.3463 here via pearson -- and it is exactly what convert_df() has
  # always done with the same row. Pinned so the trade-off cannot be reversed by
  # accident.
  direct <- suppressMessages(es_from_or_se(or = 2.5, logor_se = 0.3))
  expect_true(is.na(direct$r))
  via_df <- suppressMessages(summary(suppressMessages(
    convert_df(data.frame(or = 2.5, logor_se = 0.3), measure = "r"))))
  expect_true(is.na(via_df$es_crude))
  # the capability is still one explicit argument away
  explicit <- suppressMessages(es_from_or_se(or = 2.5, logor_se = 0.3,
                                             or_to_cor = "pearson"))
  expect_false(is.na(explicit$r))
  expect_equal(round(explicit$r, 4), 0.3463)
})

# --- 2.4: unit_type ---------------------------------------------------------

test_that("unit_type's tolerated set is documented and includes the shipped default", {
  vals <- metaConvert:::.unit_type_values()
  expect_true("sd" %in% vals)
  expect_true("raw_scale" %in% vals)          # the shipped default, previously undocumented
  expect_true("value" %in% vals)              # the name two @param blocks promised
  expect_true("raw_data" %in% vals)           # see below
  expect_equal(.convert_df_defaults$unit_type, "raw_scale")
  expect_true(.convert_df_defaults$unit_type %in% vals)
})

test_that("the tolerated set covers every spelling the repository actually uses", {
  # A first version of this guard was built from the documentation and listed three
  # names. tests_save/checked/test-ES-COR.R promptly failed: it passes "raw_data",
  # a fourth spelling that appears in no @param block and worked only because
  # everything that is not "sd" means raw units. Rebuilt from a census of the repo.
  # This test keeps the set evidence-based: if a new spelling enters the sources,
  # it fails here rather than in the archived suite.
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = FALSE)
  files <- unlist(lapply(c("R", "tests/testthat", "tests_save/checked", "vignettes"),
                         function(d) list.files(file.path(root, d),
                                                pattern = "[.](R|Rmd)$",
                                                full.names = TRUE, recursive = TRUE)))
  skip_if(length(files) == 0, "not a source checkout")
  used <- unlist(lapply(files, function(f) {
    txt <- readLines(f, warn = FALSE)
    m <- regmatches(txt, gregexpr('unit_type *= *"[^"]*"', txt))
    sub('^unit_type *= *"', "", sub('"$', "", unlist(m)))
  }))
  used <- setdiff(unique(used), "banana")     # this file's own negative fixture
  unknown <- setdiff(used, metaConvert:::.unit_type_values())
  expect_equal(length(unknown), 0L,
               info = paste("unit_type spellings in the sources that the guard would",
                            "reject:", paste(unknown, collapse = ", ")))
})

test_that("every exported taker of unit_type agrees on the default", {
  takers <- Filter(function(fn) "unit_type" %in% names(formals(get(fn, envir = NS))),
                   .exported_functions())
  expect_gte(length(takers), 8L)
  for (fn in takers) {
    got <- .arg_defaults(get(fn, envir = NS))$unit_type
    if (is.null(got)) next
    expect_equal(got, "raw_scale", info = fn)
  }
})

test_that("an invalid unit_type is rejected instead of silently meaning raw units", {
  d <- data.frame(n_exp = 20, n_nexp = 20, mean_exp = 1, mean_nexp = 0,
                  mean_sd_exp = 1, mean_sd_nexp = 1)
  expect_error(convert_df(d, measure = "g", unit_type = "banana"), "unit_type")
  expect_error(es_from_pearson_r(pearson_r = 0.3, n_sample = 100,
                                 unit_type = "banana"), "unit_type")
  expect_error(es_from_fisher_z(fisher_z = 0.3, n_sample = 100,
                                unit_type = "banana"), "unit_type")
  expect_error(es_from_linreg_t(linreg_t = 2, n_sample = 100,
                                unit_type = "banana"), "unit_type")
  expect_error(es_from_spearman_rho(spearman_r = 0.3, n_sample = 100,
                                    unit_type = "banana"), "unit_type")
})

test_that("all three tolerated values, and NA, are accepted", {
  for (v in c("sd", "raw_scale", "value"))
    expect_silent(metaConvert:::.validate_unit_type(v))
  expect_silent(metaConvert:::.validate_unit_type(NA))
  expect_silent(metaConvert:::.validate_unit_type(c("sd", NA, "raw_scale")))
  expect_silent(metaConvert:::.validate_unit_type(character(0)))
  expect_silent(metaConvert:::.validate_unit_type(NULL))
})

test_that("'value' and 'raw_scale' are synonyms, and neither is 'sd'", {
  # The whole argument is one comparison against "sd", so this is the behaviour
  # the documentation now promises.
  base <- list(pearson_r = 0.4, n_sample = 120, sd_iv = 2, unit_increase_iv = 1,
               cor_to_smd = "mathur")
  a <- do.call(es_from_pearson_r, c(base, list(unit_type = "raw_scale")))
  b <- do.call(es_from_pearson_r, c(base, list(unit_type = "value")))
  s <- do.call(es_from_pearson_r, c(base, list(unit_type = "sd")))
  expect_equal(a$d, b$d)
  expect_false(isTRUE(all.equal(a$d, s$d)))
})

test_that("a bad unit_type COLUMN warns but does not abort the run", {
  # One bad cell must never stop a convert_df() run -- the convention the
  # dipietrantonj and proportion routes already follow.
  d <- data.frame(n_exp = c(20, 20), n_nexp = c(20, 20), mean_exp = c(1, 1),
                  mean_nexp = c(0, 0), mean_sd_exp = c(1, 1), mean_sd_nexp = c(1, 1),
                  unit_type = c("sd", "banana"))
  expect_warning(res <- convert_df(d, measure = "g"), "unit_type")
  expect_equal(nrow(res[[1]]), 2L)
})
