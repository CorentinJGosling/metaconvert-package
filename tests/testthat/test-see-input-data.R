# see_input_data(): every documented `measure` value must work, and an unsupported one
# must say so.
#
# Both switch() blocks had a "Z" branch but no "z" branch, and no default. A value with no
# branch left `dat` NULL, and the following order(dat$list_input_data) died with the opaque
# "argument 1 is not a vector". "z" is listed in the function's own default `choices`
# vector and in its \usage block, so copying a value straight out of the signature crashed.

test_that("every documented measure value returns a data frame", {
  documented <- eval(formals(see_input_data)$measure)
  for (m in documented) {
    out <- see_input_data(measure = m, extension = "data.frame", verbose = FALSE)
    expect_s3_class(out, "data.frame")
    expect_gt(nrow(out), 0)
  }
})


test_that("lowercase 'z' works and agrees with 'Z' and with 'r'", {
  z_lower <- see_input_data(measure = "z", extension = "data.frame", verbose = FALSE)
  z_upper <- see_input_data(measure = "Z", extension = "data.frame", verbose = FALSE)
  r_out   <- see_input_data(measure = "r", extension = "data.frame", verbose = FALSE)

  expect_identical(z_lower, z_upper)
  # r and z are reported together in the natural/converted tags, so they select the same rows
  expect_identical(z_lower, r_out)
})


test_that("both type_of_measure settings work for every documented measure", {
  documented <- eval(formals(see_input_data)$measure)
  for (m in documented) {
    for (tm in c("natural", "natural+converted")) {
      out <- see_input_data(measure = m, type_of_measure = tm,
                            extension = "data.frame", verbose = FALSE)
      expect_s3_class(out, "data.frame")
    }
  }
})


test_that("a measure this table does not cover is rejected by name", {
  # These are real metaConvert measures, but see_input_data()'s table is built from the
  # measure = "d" hierarchy and carries no RD/ALPHA/ICC/PROP/HR tags, so it cannot report
  # them. The error must name the value and point at data_extraction_sheet().
  for (m in c("rd", "alpha", "icc", "prop", "hr")) {
    expect_error(see_input_data(measure = m, extension = "data.frame", verbose = FALSE),
                 "unsupported 'measure' value", fixed = TRUE)
    expect_error(see_input_data(measure = m, extension = "data.frame", verbose = FALSE),
                 "data_extraction_sheet()", fixed = TRUE)
  }
})


test_that("an unrecognised measure fails with the same clear error, not an internal one", {
  msg <- tryCatch(
    see_input_data(measure = "not_a_measure", extension = "data.frame", verbose = FALSE),
    error = function(e) conditionMessage(e))

  expect_true(grepl("unsupported 'measure' value", msg, fixed = TRUE))
  # and NOT the opaque internal failure this used to surface
  expect_false(grepl("argument 1 is not a vector", msg, fixed = TRUE))
  expect_false(grepl("n'est pas un vecteur", msg, fixed = TRUE))
})


test_that("every function name advertised in the sheet actually exists", {
  # The sheet's `corresponding_R_function` column is the user's map from "what I
  # extracted" to "what to call". Seven of the names were wrong -- renamed functions
  # (es_from_means_plot -> es_from_plot_means, es_from_variability_means_* ->
  # es_variab_from_means_*), a name that never existed (es_from_anova_f_pval; the
  # function is es_from_anova_pval), and one that dropped its _crude suffix
  # (es_from_ancova_means_sd_pooled). Copying any of them out of the sheet gave
  # "could not find function".
  d <- see_input_data(extension = "data.frame", verbose = FALSE)
  txt <- unlist(lapply(d, function(col) {
    if (is.character(col) || is.factor(col)) as.character(col) else NULL
  }))
  advertised <- unique(sub("[(][)]$", "",
    unlist(regmatches(txt, gregexpr("[A-Za-z_.][A-Za-z0-9_.]*[(][)]", txt)))))
  expect_gt(length(advertised), 50)
  exported <- getNamespaceExports("metaConvert")
  expect_equal(sort(setdiff(advertised, exported)), character(0))
})
