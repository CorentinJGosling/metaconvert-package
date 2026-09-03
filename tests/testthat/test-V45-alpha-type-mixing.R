# V45 -- a pool mixing raw (covariance-matrix) and standardised (correlation-matrix)
# alpha. The alpha member of the V38 / V43 family.
#
# Raw and standardised alpha are different coefficients, equal only when the item
# variances are. psych::alpha() prints both (raw_alpha, std.alpha), so an extractor has
# the pair in front of them with nothing in the sheet recording which was taken -- and
# alpha_type enters neither the estimate nor the (n, k) standard error, so no numeric
# check anywhere can reveal the mix.

pool <- function(types) {
  data.frame(study_id = paste0("S", seq_along(types)),
             cronbach_alpha = c(0.81, 0.84, 0.88, 0.79)[seq_along(types)],
             n_sample = c(200, 260, 180, 240)[seq_along(types)],
             n_items = 10, alpha_type = types, stringsAsFactors = FALSE)
}
flags_of <- function(d) {
  s <- as.data.frame(suppressMessages(suppressWarnings(
    summary(convert_df(d, measure = "alpha", verbose = FALSE), flags = TRUE))))
  if ("flags_crude" %in% names(s)) s$flags_crude else s$flags
}

test_that("V45 fires only when two KNOWN alpha types coexist", {
  f <- flags_of(pool(c("raw", "std.alpha", "raw", "raw")))
  expect_true(all(grepl("Pool mixes alpha types", f)))
  expect_match(f[1], "covariance-matrix", fixed = TRUE)
  expect_match(f[2], "correlation-matrix", fixed = TRUE)
  expect_match(f[1], "[INFO]", fixed = TRUE)
})

test_that("a homogeneous pool is silent, whichever type it is", {
  for (ty in c("raw", "std.alpha", "covariance", "correlation")) {
    expect_false(any(grepl("Pool mixes alpha types", flags_of(pool(rep(ty, 4))))),
                 info = ty)
  }
})

test_that("unspecified is not a level, so blanks plus one known type stay silent", {
  # V39's rule, not V43's: there is no package default alpha_type, the coefficient
  # comes from the study, and most primary reports never say which they computed.
  # Counting blanks as raw would fire on nearly every real dataset, on a guess.
  expect_false(any(grepl("Pool mixes alpha types", flags_of(pool(c(NA, NA, NA, "raw"))))))
  expect_false(any(grepl("Pool mixes alpha types",
                         flags_of(pool(c("unspecified", "not reported", NA, "std.alpha"))))))
  # ...but two known types still fire even with blanks present
  expect_true(any(grepl("Pool mixes alpha types", flags_of(pool(c(NA, "raw", "std.alpha", NA))))))
})

test_that("the column is absent-safe and does not change any estimate", {
  d <- pool(rep("raw", 4)); d$alpha_type <- NULL
  expect_false(any(grepl("Pool mixes alpha types", flags_of(d))))
  a <- as.data.frame(suppressMessages(suppressWarnings(
        summary(convert_df(d, measure = "alpha", verbose = FALSE)))))
  b <- as.data.frame(suppressMessages(suppressWarnings(
        summary(convert_df(pool(c("raw", "std.alpha", "raw", "correlation")),
                           measure = "alpha", verbose = FALSE)))))
  expect_equal(a$es_crude, b$es_crude, tolerance = 1e-12)
  expect_equal(a$se_crude, b$se_crude, tolerance = 1e-12)
})

test_that("the normaliser maps what extractors actually see, and never aborts", {
  f <- metaConvert:::.normalise_alpha_type
  expect_equal(f(c("raw", "raw_alpha", "Cronbach", "cov")), rep("covariance", 4))
  expect_equal(f(c("std", "std.alpha", "Standardised", "correlation")), rep("correlation", 4))
  expect_equal(suppressWarnings(f("nonsense")), "unspecified")
  expect_warning(f("nonsense"), "Unrecognised alpha_type")
})
