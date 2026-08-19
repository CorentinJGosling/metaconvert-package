# Regression tests for the Wan (2014) SD-constant table lookup in the median/quartile
# converters. A row with n <= 1 makes the table index 0 (Q = (n-1)/4 <= 0, or n itself
# for the range table); R DROPS an index-0 element, shortening the constant vector so the
# enclosing ifelse() recycles a constant onto the WRONG (valid) neighbouring row --
# silently corrupting that study's SD/effect size with no error. The fix
# (.wan_table_lookup) maps out-of-range indices to NA so each row keeps its own constant.
#
# The decisive assertion is INVARIANCE: computing a batch that also contains an n <= 1
# row must leave every OTHER row's effect size exactly equal to its standalone value.

test_that("es_from_med_quarts: an n=1 row does not corrupt neighbouring rows", {
  standalone <- es_from_med_quarts(
    q1_exp = c(10, 10), med_exp = c(20, 20), q3_exp = c(40, 40), n_exp = c(9, 13),
    q1_nexp = c(12, 12), med_nexp = c(22, 22), q3_nexp = c(41, 41), n_nexp = c(9, 13)
  )
  batch <- es_from_med_quarts(
    q1_exp = c(10, 10, 10), med_exp = c(20, 20, 20), q3_exp = c(40, 40, 40), n_exp = c(1, 9, 13),
    q1_nexp = c(12, 12, 12), med_nexp = c(22, 22, 22), q3_nexp = c(41, 41, 41), n_nexp = c(1, 9, 13)
  )
  # neighbours (rows 2,3 = n 9,13) must be identical to their standalone values
  expect_equal(batch$d[2:3], standalone$d[1:2])
  expect_equal(batch$g[2:3], standalone$g[1:2])
  # the n=1 row itself is undefined -> NA (not a corrupted finite number)
  expect_true(is.na(batch$d[1]))
})

test_that("es_from_med_quarts: corruption also blocked via the control-arm (n_nexp=1) path", {
  standalone <- es_from_med_quarts(
    q1_exp = c(10, 10), med_exp = c(20, 20), q3_exp = c(40, 40), n_exp = c(9, 13),
    q1_nexp = c(12, 12), med_nexp = c(22, 22), q3_nexp = c(41, 41), n_nexp = c(9, 13)
  )
  batch <- es_from_med_quarts(
    q1_exp = c(10, 10, 10), med_exp = c(20, 20, 20), q3_exp = c(40, 40, 40), n_exp = c(9, 9, 13),
    q1_nexp = c(12, 12, 12), med_nexp = c(22, 22, 22), q3_nexp = c(41, 41, 41), n_nexp = c(1, 9, 13)
  )
  expect_equal(batch$d[2:3], standalone$d[1:2])
  expect_true(is.na(batch$d[1]))
})

test_that("es_from_med_min_max_quarts: an n=1 row does not corrupt neighbouring rows", {
  standalone <- es_from_med_min_max_quarts(
    min_exp = c(5, 5), q1_exp = c(10, 10), med_exp = c(20, 20), q3_exp = c(40, 40), max_exp = c(60, 60), n_exp = c(9, 13),
    min_nexp = c(6, 6), q1_nexp = c(12, 12), med_nexp = c(22, 22), q3_nexp = c(41, 41), max_nexp = c(62, 62), n_nexp = c(9, 13)
  )
  batch <- es_from_med_min_max_quarts(
    min_exp = c(5, 5, 5), q1_exp = c(10, 10, 10), med_exp = c(20, 20, 20), q3_exp = c(40, 40, 40), max_exp = c(60, 60, 60), n_exp = c(1, 9, 13),
    min_nexp = c(6, 6, 6), q1_nexp = c(12, 12, 12), med_nexp = c(22, 22, 22), q3_nexp = c(41, 41, 41), max_nexp = c(62, 62, 62), n_nexp = c(1, 9, 13)
  )
  expect_equal(batch$d[2:3], standalone$d[1:2])
  expect_true(is.na(batch$d[1]))
})

test_that("es_from_med_min_max: an n=0 row does not corrupt neighbouring rows (range table)", {
  # The range-table lookup list1[n] drops at n=0 (index 0) and does negative indexing
  # at n<0; guard it too. n=1 (list1[1]=0 -> Inf) is left to its own row (V16 territory).
  standalone <- es_from_med_min_max(
    min_exp = c(5, 5), med_exp = c(20, 20), max_exp = c(60, 60), n_exp = c(9, 13),
    min_nexp = c(6, 6), med_nexp = c(22, 22), max_nexp = c(62, 62), n_nexp = c(9, 13)
  )
  batch <- es_from_med_min_max(
    min_exp = c(5, 5, 5), med_exp = c(20, 20, 20), max_exp = c(60, 60, 60), n_exp = c(0, 9, 13),
    min_nexp = c(6, 6, 6), med_nexp = c(22, 22, 22), max_nexp = c(62, 62, 62), n_nexp = c(0, 9, 13)
  )
  expect_equal(batch$d[2:3], standalone$d[1:2])
})

test_that("convert_df: n=1 median row is flagged (V16) and neighbours stay correct", {
  # Interaction with the Tier-1 flag layer (internal_flags.R): the n=1 row must be
  # flagged 'within-group variance is undefined' AND yield NA, while the pipeline still
  # returns the correct effect sizes for the valid neighbouring studies.
  df <- data.frame(
    q1_exp = c(10, 10, 10), med_exp = c(20, 20, 20), q3_exp = c(40, 40, 40), n_exp = c(1, 9, 13),
    q1_nexp = c(12, 12, 12), med_nexp = c(22, 22, 22), q3_nexp = c(41, 41, 41), n_nexp = c(1, 9, 13)
  )
  s <- suppressWarnings(suppressMessages(
    summary(convert_df(df, measure = "d"), flags = TRUE)
  ))
  es_col <- grep("^es_crude$|^es$", colnames(s), value = TRUE)[1]
  fl_col <- grep("flag", colnames(s), value = TRUE)[1]

  # neighbours computed, n=1 row NA
  expect_true(is.na(s[[es_col]][1]))
  expect_false(any(is.na(s[[es_col]][2:3])))

  # neighbours match a direct standalone computation (no corruption through the pipeline)
  direct <- es_from_med_quarts(
    q1_exp = c(10, 10), med_exp = c(20, 20), q3_exp = c(40, 40), n_exp = c(9, 13),
    q1_nexp = c(12, 12), med_nexp = c(22, 22), q3_nexp = c(41, 41), n_nexp = c(9, 13)
  )
  # summary() used to round es_crude IN the returned data.frame, so this had to
  # compare against round(direct$d, 3) at a 1e-3 tolerance. The analysis columns
  # are now returned at full precision, so the pipeline and the standalone route
  # agree to machine precision -- a strictly stronger assertion.
  expect_equal(unname(s[[es_col]][2:3]), direct$d[1:2], tolerance = 1e-10)

  # V16 flag present on the n=1 row
  expect_match(s[[fl_col]][1], "n_exp = 1|within-group variance is undefined")
})
