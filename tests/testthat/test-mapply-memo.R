## .mapply_memo() memoises the row-wise correlation routes (.contingency_to_cor,
## .or_to_cor) whose tetrachoric step is solved by numerical ML and costs ~20 ms
## per row. Deduplication must be an OPTIMISATION ONLY: identical output, in the
## original row order, for every input shape.

test_that(".mapply_memo reproduces t(mapply()) exactly", {
  f <- function(x, y) c(x + y, x * y, x - y)

  # all rows distinct -> takes the passthrough branch
  x <- c(1, 2, 3, 4); y <- c(5, 6, 7, 8)
  expect_identical(.mapply_memo(f, x = x, y = y), t(mapply(f, x = x, y = y)))

  # heavy duplication, duplicates NOT adjacent -> order must be preserved
  x <- c(1, 2, 1, 3, 2, 1, 3); y <- c(9, 8, 9, 7, 8, 9, 7)
  expect_identical(.mapply_memo(f, x = x, y = y), t(mapply(f, x = x, y = y)))

  # a single row
  expect_identical(.mapply_memo(f, x = 2, y = 3), t(mapply(f, x = 2, y = 3)))

  # mixed types in the key (character argument, as table_2x2_to_cor is)
  g <- function(a, lab) c(a * 2, nchar(lab))
  a <- c(1, 1, 2, 2, 1); lab <- c("aa", "aa", "b", "ccc", "aa")
  expect_identical(.mapply_memo(g, a = a, lab = lab), t(mapply(g, a = a, lab = lab)))

  # NA in the arguments must not collapse distinct rows together
  h <- function(x, y) c(x, y)
  x <- c(1, NA, 1, NA); y <- c(NA, 1, NA, 2)
  expect_identical(.mapply_memo(h, x = x, y = y), t(mapply(h, x = x, y = y)))
})

test_that("es_from_2x2 is unaffected by deduplication of repeated tables", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  # the same table repeated must give the same row as the table on its own
  one <- es_from_2x2(n_cases_exp = 15, n_controls_exp = 25,
                     n_cases_nexp = 8, n_controls_nexp = 32)
  many <- es_from_2x2(n_cases_exp    = c(15, 4, 15, 20, 15),
                      n_controls_exp = c(25, 36, 25, 20, 25),
                      n_cases_nexp   = c(8, 9, 8, 11, 8),
                      n_controls_nexp = c(32, 31, 32, 29, 32))

  expect_equal(many$r[c(1, 3, 5)], rep(one$r, 3), tolerance = 0)
  expect_equal(many$z[c(1, 3, 5)], rep(one$z, 3), tolerance = 0)
  expect_equal(many$r_se[c(1, 3, 5)], rep(one$r_se, 3), tolerance = 0)

  # the interleaved distinct rows must equal their own standalone computation
  alone2 <- es_from_2x2(n_cases_exp = 4, n_controls_exp = 36,
                        n_cases_nexp = 9, n_controls_nexp = 31)
  expect_equal(many$r[2], alone2$r, tolerance = 0)
})

test_that("reverse_2x2 is part of the memoisation key", {
  skip_if_not_installed("mvtnorm")  # metafor RTET (tetrachoric) requires it
  # identical counts, opposite reverse flags: must NOT collapse to one answer
  res <- es_from_2x2(n_cases_exp     = c(15, 15),
                     n_controls_exp  = c(25, 25),
                     n_cases_nexp    = c(8, 8),
                     n_controls_nexp = c(32, 32),
                     reverse_2x2     = c(FALSE, TRUE))
  expect_equal(res$r[1], -res$r[2], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(res$r[1], res$r[2])))
})
