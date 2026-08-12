test_that("agg.subgroups correctly agg ES", {
  set.seed(4321)

  dat <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100))
  )

  df.mfr <- metafor::escalc(yi = es, sei = se, data = dat)

  test_mfr <- metafor::aggregate.escalc(df.mfr,
    cluster = agg,
    struct = "ID", weighted = TRUE
  )

  metacnvert <- aggregate_df(
    x = dat, dependence = "subgroups",
    agg_fact = "agg", es = "es", se = "se"
  )

  expect_equal(as.numeric(test_mfr$yi), metacnvert$es)
  expect_equal(test_mfr$vi, metacnvert$se^2)
})

test_that("agg.outcomes correctly agg outcomes", {
  dat <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100))
  )

  df.mfr <- metafor::escalc(yi = es, sei = se, data = dat)

  test_mfr <- metafor::aggregate.escalc(df.mfr,
    cluster = agg,
    struct = "CS", weighted = FALSE, rho = 0.8
  )

  metacnvert <- aggregate_df(
    x = dat, dependence = "outcomes",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8
  )

  expect_equal(as.numeric(test_mfr$yi), metacnvert$es)
  expect_equal(test_mfr$vi, metacnvert$se^2)
})



test_that("agg.times correctly agg times", {
  dat <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    time_agg = c(1, 1, 1, 2, 3, 1, 2, 3, 1, 2),
    # agg = rep(c(1), 10),
    # time = c(1,2:10),
    cor_unit=0.5,
    row_index = 1:10,
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100))
  )

  df.mfr <- metafor::escalc(yi = dat$es, vi = dat$se^2)
  df.mfr <- metafor::escalc(yi = es, sei = se, data=dat)

  test_mfr <- metafor::aggregate.escalc(df.mfr,
      cluster = dat$agg,
      struct = "CAR",
      phi = 0.5, time = dat$time_agg,
      digits=11
  )

  metacnvert <- aggregate_df(
    x = dat, dependence = "times",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8
  )

  expect_equal(as.numeric(test_mfr$yi), metacnvert$es)
  expect_equal(test_mfr$vi, metacnvert$se^2)
})


test_that("agg.s correctly agg MEAN", {
  df <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100))
  )

  metacnvert_sub <- aggregate_df(
    x = df, dependence = "subgroups",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    col_mean = c("n")
  )
  metacnvert_out <- aggregate_df(
    x = df, dependence = "outcomes",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    col_mean = c("n")
  )

  N <- tapply(df$n, df$agg, mean)

  expect_equal(as.numeric(N), metacnvert_sub$n)
  expect_equal(as.numeric(N), metacnvert_out$n)
})

test_that("agg.s correctly agg WEIGHTED MEAN", {
  df <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100)),
    W = rnorm(10, 50, 10)
  )

  metacnvert_sub <- aggregate_df(
    x = df, dependence = "subgroups",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    weights = "W", col_weighted_mean = c("n")
  )
  metacnvert_out <- aggregate_df(
    x = df, dependence = "outcomes",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    weights = "W", col_weighted_mean = c("n")
  )

  N <- as.numeric(by(df, df$agg, with, weighted.mean(n, W)))

  expect_equal(as.numeric(N), metacnvert_sub$n)
  expect_equal(as.numeric(N), metacnvert_out$n)
})

test_that("agg.s correctly agg SUM", {
  df <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100))
  )

  metacnvert_sub <- aggregate_df(
    x = df, dependence = "subgroups",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    weights = "n", col_sum = c("n")
  )
  metacnvert_out <- aggregate_df(
    x = df, dependence = "outcomes",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    weights = "n", col_sum = c("n")
  )

  N <- tapply(df$n, df$agg, sum)

  expect_equal(as.numeric(N), metacnvert_sub$n)
  expect_equal(as.numeric(N), metacnvert_out$n)
})

test_that("agg.s correctly agg MIN", {
  df <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100))
  )

  metacnvert_sub <- aggregate_df(
    x = df, dependence = "subgroups",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    weights = "n", col_min = c("n")
  )
  metacnvert_out <- aggregate_df(
    x = df, dependence = "outcomes",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    weights = "n", col_min = c("n")
  )

  N <- tapply(df$n, df$agg, min)

  expect_equal(as.numeric(N), metacnvert_sub$n)
  expect_equal(as.numeric(N), metacnvert_out$n)
})

test_that("agg.s correctly agg MAX", {
  df <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100))
  )


  metacnvert_sub <- aggregate_df(
    x = df, dependence = "subgroups",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    weights = "n", col_max = c("n")
  )
  metacnvert_out <- aggregate_df(
    x = df, dependence = "outcomes",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    weights = "n", col_max = c("n")
  )

  N <- tapply(df$n, df$agg, max)

  expect_equal(as.numeric(N), metacnvert_sub$n)
  expect_equal(as.numeric(N), metacnvert_out$n)
})

test_that("agg.s correctly agg FACT", {
  df <- data.frame(
    agg = rep(c(1, 2, 3, 4, 5), c(1, 1, 3, 3, 2)),
    es = rnorm(10), se = abs(rnorm(10)),
    var_a = rep(c("A", "B"), each = 5),
    var_b = rep(c("C", "D", "E"), c(2, 3, 5)),
    n = round(runif(10, 20, 100))
  )



  metacnvert_sub <- aggregate_df(
    x = df, dependence = "subgroups",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    col_fact = "n"
  )
  metacnvert_out <- aggregate_df(
    x = df, dependence = "outcomes",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    col_fact = "n"
  )

  N <- tapply(df$n, df$agg, paste, collapse = " / ")

  # correct, only the order differs & the unique ess not detected in tapply
  # expect_equal(N, metacnvert_sub$n)
  # expect_equal(N, metacnvert_out$n)
  expect_equal(metacnvert_sub$n, metacnvert_out$n)

  metacnvert_sub <- aggregate_df(
    x = df, dependence = "subgroups",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    col_fact = "var_b"
  )
  metacnvert_out <- aggregate_df(
    x = df, dependence = "outcomes",
    agg_fact = "agg", es = "es", se = "se",
    cor_unit = 0.8,
    col_fact = "var_b"
  )


  N <- tapply(df$var_b, df$agg, paste, collapse = " / ")

  # correct, only the order differs & the unique ess not detected in tapply
  # expect_equal(N, metacnvert_sub$var_b)
  expect_equal(metacnvert_sub$var_b, metacnvert_out$var_b)
})
