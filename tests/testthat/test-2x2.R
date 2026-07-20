
# 2x2 to OR -----
test_that("OR from 2x2", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp

  mfr1 <- metafor::escalc(
    ai = n_cases_exp,
    bi = n_controls_exp,
    ci = n_cases_nexp,
    di = n_controls_nexp, data = dat, measure = "OR", digits = 11
  )
  mfr2 <- metafor::escalc(
    ai = n_cases_exp,
    bi = n_cases_nexp,
    ci = n_controls_exp,
    di = n_controls_nexp, data = dat, measure = "OR", digits = 11
  )

  es.mcv_or <- summary(convert_df(dat, verbose = FALSE,
                                  es_selected = "hierarchy", hierarchy = "2x2", measure = "logor"),
                       digits = 11)

  ## OR
  expect_equal(unique(es.mcv_or$info_used_crude), "2x2")
  expect_equal(es.mcv_or$es_crude, as.numeric(mfr2$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude^2, as.numeric(mfr2$vi), tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, as.numeric(mfr1$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude^2, as.numeric(mfr1$vi), tolerance = 1e-10)
})

test_that("OR from 2x2-sum", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE,
                                  es_selected = "hierarchy", hierarchy = "2x2",
                                  measure = "logor"), digits = 11)
  es.mcv_or2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "logor"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_or$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_or2$info_used_crude), "2x2_sum")
  expect_equal(es.mcv_or$es_crude, es.mcv_or2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or2$se_crude, tolerance = 1e-10)
})

test_that("OR from prop", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "logor"), digits = 11)
  es.mcv_or2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "logor"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_or$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_or2$info_used_crude), "2x2_prop")
  expect_equal(es.mcv_or$es_crude, es.mcv_or2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or2$se_crude, tolerance = 1e-10)
})

# 2x2 to RR -----
test_that("RR from 2x2", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp

  mfr1 <- metafor::escalc(
    ai = n_cases_exp,
    bi = n_controls_exp,
    ci = n_cases_nexp,
    di = n_controls_nexp, data = dat, measure = "RR", digits = 11
  )

  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "logrr"), digits = 11)

  ## OR
  expect_equal(unique(es.mcv_or$info_used_crude), "2x2")
  expect_equal(es.mcv_or$es_crude, as.numeric(mfr1$yi), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, sqrt(as.numeric(mfr1$vi)), tolerance = 1e-10)
})

test_that("RR from 2x2-sum", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_rr <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "logrr"), digits = 11)
  es.mcv_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "logrr"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_rr$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_rr2$info_used_crude), "2x2_sum")
  expect_equal(es.mcv_rr$es_crude, es.mcv_rr2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr$se_crude, es.mcv_rr2$se_crude, tolerance = 1e-10)
})

test_that("RR from prop", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_rr <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "logrr"), digits = 11)
  es.mcv_rr2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "logrr"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_rr$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_rr2$info_used_crude), "2x2_prop")
  expect_equal(es.mcv_rr$es_crude, es.mcv_rr2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_rr$se_crude, es.mcv_rr2$se_crude, tolerance = 1e-10)
})

# 2x2 to SMD -----
test_that("SMD from 2x2", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE,
                                 es_selected = "hierarchy", hierarchy = "2x2",
                                 measure = "d"), digits = 11)

  comp_res_d = NULL
  for (i in 1:nrow(dat)) {
    comp_res_d_i <- data.frame(esc::esc_2x2(
      grp1yes = dat$n_cases_exp[i],
      grp1no = dat$n_controls_exp[i],
      grp2yes = dat$n_cases_nexp[i],
      grp2no = dat$n_controls_nexp[i],
      es.type = "d", digits = 11
    ))
    comp_res_d = rbind(comp_res_d, comp_res_d_i)
  }
  ## SMD
  expect_equal(unique(es.mcv_d$info_used_crude), "2x2")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 1e-10)
})

test_that("SMD from 2x2-sum", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "d"), digits = 11)
  es.mcv_d2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "d"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_d$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_d2$info_used_crude), "2x2_sum")
  expect_equal(es.mcv_d$es_crude, es.mcv_d2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d2$se_crude, tolerance = 1e-10)
})

test_that("SMD from prop", {
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "d"), digits = 11)
  es.mcv_d2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "d"), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_d$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_d2$info_used_crude), "2x2_prop")
  expect_equal(es.mcv_d$es_crude, es.mcv_d2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d2$se_crude, tolerance = 1e-10)
})

# 2x2 to COR -----
# test_that("2x2 to R (COOPER)", {
#   dat <- metaumbrella::df.OR
#   dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
#   dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
#
#   es.mcv_or <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2",
#     measure = "logor"
#   ), digits = 11)
#
#   es.mcv_r <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2",
#     measure = "r",
#     table_2x2_to_cor = "cooper"
#   ), digits = 11)
#
#   comp_r <- compute.es::lores(
#     lor = es.mcv_or$es_crude,
#     var.lor = es.mcv_or$se_crude^2,
#     verbose = FALSE,
#     n.1 = es.mcv_or$n_exp, n.2 = es.mcv_or$n_nexp,
#     dig = 12
#   )
#   esc_r = NULL
#   for (i in 1:nrow(dat)) {
#     esc_r_i <- data.frame(esc::esc_2x2(
#       grp1yes = dat$n_cases_exp[i],
#       grp1no = dat$n_controls_exp[i],
#       grp2yes = dat$n_cases_nexp[i],
#       grp2no = dat$n_controls_nexp[i],
#       es.type = "r", digits = 11
#     ))
#     esc_r = rbind(esc_r, esc_r_i)
#   }
#
#   expect_equal(unique(es.mcv_r$info_used_crude), "2x2")
#   expect_equal(es.mcv_r$es_crude, comp_r$r, tolerance = 1e-10)
#   expect_equal(es.mcv_r$se_crude^2, comp_r$var.r, tolerance = 1e-10)
#   expect_equal(es.mcv_r$es_crude, esc_r$es, tolerance = 1e-10)
#   # expect_equal(es.mcv_r$se_crude, esc_r$se, tolerance = 1e-10)
# })

# test_that("2x2 to Z (COOPER)", {
#   dat <- metaumbrella::df.OR
#   dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
#   dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
#   es.mcv_or <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2",
#     measure = "logor"
#   ), digits = 11)
#
#   es.mcv_r <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2",
#     measure = "z",
#     table_2x2_to_cor = "cooper"
#   ), digits = 11)
#
#   comp_r <- compute.es::lores(
#     lor = es.mcv_or$es_crude,
#     var.lor = es.mcv_or$se_crude^2,
#     verbose = FALSE,
#     n.1 = es.mcv_or$n_exp, n.2 = es.mcv_or$n_nexp,
#     dig = 12
#   )
#   esc_r = NULL
#   for (i in 1:nrow(dat)) {
#     esc_r_i <- data.frame(esc::esc_2x2(
#       grp1yes = dat$n_cases_exp[i],
#       grp1no = dat$n_controls_exp[i],
#       grp2yes = dat$n_cases_nexp[i],
#       grp2no = dat$n_controls_nexp[i],
#       es.type = "r", digits = 11
#     ))
#     esc_r = rbind(esc_r, esc_r_i)
#   }
#   expect_equal(unique(es.mcv_r$info_used_crude), "2x2")
#   expect_equal(es.mcv_r$es_crude, comp_r$fisher.z, tolerance = 1e-10)
#   # expect_equal(es.mcv_r$se_crude^2, comp_r$var.z, tolerance = 1e-2)
#   expect_equal(es.mcv_r$es_crude, esc_r$fishers.z, tolerance = 1e-10)
#   expect_equal(es.mcv_r$se_crude, (esc_r$ci.hi.z - esc_r$ci.lo.z)/(2*qnorm(.975)), tolerance = 1e-10)
# })

test_that("2x2 to R (TETRACHORIC)", {
  skip_if_not_installed("mvtnorm")
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
    table_2x2_to_cor = "tetrachoric"
  ), digits = 11)
  comp_res_r <- summary(metafor::escalc(
    ai = dat$n_cases_exp,
    bi = dat$n_controls_exp,
    ci = dat$n_cases_nexp,
    di = dat$n_controls_nexp,
    measure = "RTET", digits = 11
  ))

  expect_equal(unique(es.mcv_r$info_used_crude), "2x2")
  expect_equal(es.mcv_r$es_crude, as.numeric(comp_res_r$yi), tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, sqrt(comp_res_r$vi), tolerance = 1e-10)
  expect_equal(es.mcv_r$es_ci_lo_crude, comp_res_r$ci.lb, tolerance = 1e-10)
  expect_equal(es.mcv_r$es_ci_up_crude, comp_res_r$ci.ub, tolerance = 1e-10)
})

test_that("2x2 to Z (TETRACHORIC)", {
  skip_if_not_installed("mvtnorm")
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
    table_2x2_to_cor = "tetrachoric"
  ), digits = 11)
  comp_res_r <- summary(metafor::escalc(
    ai = dat$n_cases_exp,
    bi = dat$n_controls_exp,
    ci = dat$n_cases_nexp,
    di = dat$n_controls_nexp,
    measure = "ZTET", digits = 11)
  )

  expect_equal(unique(es.mcv_r$info_used_crude), "2x2")
  expect_equal(es.mcv_r$es_crude, as.numeric(comp_res_r$yi), tolerance = 1e-10)
  expect_equal(es.mcv_r$se_crude, sqrt(comp_res_r$vi), tolerance = 1e-10)
  expect_equal(es.mcv_r$es_ci_lo_crude, comp_res_r$ci.lb, tolerance = 1e-10)
  expect_equal(es.mcv_r$es_ci_up_crude, comp_res_r$ci.ub, tolerance = 1e-10)
})

# test_that("2x2 to R (lipsey)", {
#   dat <- metaumbrella::df.OR
#   dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
#   dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
#   es.mcv_r <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
#     table_2x2_to_cor = "lipsey"
#   ), digits = 11)
#   es.mcv_rt <- summary(convert_df(dat,
#                                  verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
#                                  table_2x2_to_cor = "tetrachoric"
#   ), digits = 11)
#
#   expect_equal(unique(es.mcv_r$info_used_crude), "2x2")
#   expect_equal(es.mcv_r$es_crude, es.mcv_rt$es_crude, tolerance = 3e-1)
#   expect_equal(es.mcv_r$se_crude, es.mcv_rt$se_crude, tolerance = 2e-1)
#   expect_true(mean(abs(es.mcv_r$es_crude-es.mcv_rt$es_crude)) < 0.10)
#   expect_true(mean(abs(es.mcv_r$se_crude-es.mcv_rt$se_crude)) < 0.10)
#   # != formulas (esc : 2x2 => d => r)
# })
# test_that("2x2 to Z (Lipsey)", {
#   dat <- metaumbrella::df.OR
#   dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
#   dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
#   es.mcv_r <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
#     table_2x2_to_cor = "lipsey"
#   ), digits = 11)
#   es.mcv_rt <- summary(convert_df(dat,
#                                   verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2",
#                                   measure = "z",
#                                   table_2x2_to_cor = "tetrachoric"
#   ), digits = 11)
#
#   expect_equal(unique(es.mcv_r$info_used_crude), "2x2")
#   expect_equal(es.mcv_r$es_crude, es.mcv_rt$es_crude, tolerance = 3e-1)
#   expect_equal(es.mcv_r$se_crude, es.mcv_rt$se_crude, tolerance = 2e-1)
#   expect_true(mean(abs(es.mcv_r$es_crude-es.mcv_rt$es_crude)) < 0.10)
#   expect_true(mean(abs(es.mcv_r$se_crude-es.mcv_rt$se_crude)) < 0.10)
#   # != formulas (esc : 2x2 => d => r)
# })



# INTERNAL for COR------
# test_that("cooper R - 2x2 & sum & prop", {
#   dat <- metaumbrella::df.OR
#   dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
#   dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
#
#   es.mcv_or1 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
#     table_2x2_to_cor = "cooper"
#   ), digits = 11)
#   es.mcv_or2 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
#     table_2x2_to_cor = "cooper"
#   ), digits = 11)
#   es.mcv_or3 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
#     table_2x2_to_cor = "cooper"
#   ), digits = 11)
#
#   ## test ES
#   expect_equal(unique(es.mcv_or1$info_used_crude), "2x2")
#   expect_equal(unique(es.mcv_or2$info_used_crude), "2x2_sum")
#   expect_equal(unique(es.mcv_or3$info_used_crude), "2x2_prop")
#   expect_equal(es.mcv_or1$es_crude, es.mcv_or2$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$se_crude, es.mcv_or2$se_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$es_crude, es.mcv_or3$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$se_crude, es.mcv_or3$se_crude, tolerance = 1e-10)
# })
# test_that("cooper Z - 2x2 & 2x2 & prop", {
#   dat <- metaumbrella::df.OR
#   dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
#   dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
#
#   es.mcv_or1 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
#     table_2x2_to_cor = "cooper"
#   ), digits = 11)
#   es.mcv_or2 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
#     table_2x2_to_cor = "cooper"
#   ), digits = 11)
#   es.mcv_or3 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
#     table_2x2_to_cor = "cooper"
#   ), digits = 11)
#
#   ## test ES
#   expect_equal(unique(es.mcv_or1$info_used_crude), "2x2")
#   expect_equal(unique(es.mcv_or2$info_used_crude), "2x2_sum")
#   expect_equal(unique(es.mcv_or3$info_used_crude), "2x2_prop")
#   expect_equal(es.mcv_or1$es_crude, es.mcv_or2$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$se_crude, es.mcv_or2$se_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$es_crude, es.mcv_or3$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$se_crude, es.mcv_or3$se_crude, tolerance = 1e-10)
# })

test_that("tetra R - 2x2 & 2x2 & prop", {
  skip_if_not_installed("mvtnorm")
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp

  es.mcv_or1 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
    table_2x2_to_cor = "tetrachoric"
  ), digits = 11)
  es.mcv_or2 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
    table_2x2_to_cor = "tetrachoric"
  ), digits = 11)
  es.mcv_or3 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
    table_2x2_to_cor = "tetrachoric"
  ), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_or1$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_or2$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_or3$info_used_crude), "2x2_prop")
  expect_equal(es.mcv_or1$es_crude, es.mcv_or2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or1$se_crude, es.mcv_or2$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or1$es_crude, es.mcv_or3$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or1$se_crude, es.mcv_or3$se_crude, tolerance = 1e-10)
})
test_that("tetra Z - 2x2 & 2x2 & prop", {
  skip_if_not_installed("mvtnorm")
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp

  es.mcv_or1 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
    table_2x2_to_cor = "tetrachoric"
  ), digits = 11)
  es.mcv_or2 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
    table_2x2_to_cor = "tetrachoric"
  ), digits = 11)
  es.mcv_or3 <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
    table_2x2_to_cor = "tetrachoric"
  ), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_or1$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_or2$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_or3$info_used_crude), "2x2_prop")
  expect_equal(es.mcv_or1$es_crude, es.mcv_or2$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or1$se_crude, es.mcv_or2$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or1$es_crude, es.mcv_or3$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or1$se_crude, es.mcv_or3$se_crude, tolerance = 1e-10)
})

# test_that("lipsey R - 2x2 & 2x2 & prop", {
#   dat <- metaumbrella::df.OR
#   dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
#   dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
#
#   es.mcv_or1 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
#     table_2x2_to_cor = "lipsey"
#   ), digits = 11)
#   es.mcv_or2 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
#     table_2x2_to_cor = "lipsey"
#   ), digits = 11)
#   es.mcv_or3 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
#     table_2x2_to_cor = "lipsey"
#   ), digits = 11)
#
#   ## test ES
#   expect_equal(unique(es.mcv_or1$info_used_crude), "2x2")
#   expect_equal(unique(es.mcv_or2$info_used_crude), "2x2_sum")
#   expect_equal(unique(es.mcv_or3$info_used_crude), "2x2_prop")
#   expect_equal(es.mcv_or1$es_crude, es.mcv_or2$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$se_crude, es.mcv_or2$se_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$es_crude, es.mcv_or3$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$se_crude, es.mcv_or3$se_crude, tolerance = 1e-10)
# })
# test_that("lipsey Z - 2x2 & 2x2 & prop", {
#   dat <- metaumbrella::df.OR
#   dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
#   dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp
#
#   es.mcv_or1 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
#     table_2x2_to_cor = "lipsey"
#   ), digits = 11)
#   es.mcv_or2 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
#     table_2x2_to_cor = "lipsey"
#   ), digits = 11)
#   es.mcv_or3 <- summary(convert_df(dat,
#     verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
#     table_2x2_to_cor = "lipsey"
#   ), digits = 11)
#
#   ## test ES
#   expect_equal(unique(es.mcv_or1$info_used_crude), "2x2")
#   expect_equal(unique(es.mcv_or2$info_used_crude), "2x2_sum")
#   expect_equal(unique(es.mcv_or3$info_used_crude), "2x2_prop")
#   expect_equal(es.mcv_or1$es_crude, es.mcv_or2$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$se_crude, es.mcv_or2$se_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$es_crude, es.mcv_or3$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_or1$se_crude, es.mcv_or3$se_crude, tolerance = 1e-10)
# })

# REVERSE -----
test_that("2x2 - Reverse", {
  skip_if_not_installed("mvtnorm")
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp

  dat$reverse_2x2 <- FALSE
  es.mcv_2x2_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "d"), digits = 11)
  es.mcv_2x2_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "g"), digits = 11)
  es.mcv_2x2_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "logor"), digits = 11)
  es.mcv_2x2_rr <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "logrr"), digits = 11)
  # es.mcv_2x2_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
  #                                     table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
  #                                     table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_r3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
                                      table_2x2_to_cor = "tetrachoric"), digits = 11)
  # es.mcv_2x2_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
  #                                    table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
  #                                    table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_z3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
                                     table_2x2_to_cor = "tetrachoric"), digits = 11)


  dat$reverse_2x2 <- TRUE
  es.mcv_2x2_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "d"), digits = 11)
  es.mcv_2x2_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "g"), digits = 11)
  es.mcv_2x2_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "logor"), digits = 11)
  es.mcv_2x2_rr_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "logrr"), digits = 11)
  # es.mcv_2x2_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
  #                                     table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
  #                                     table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_r3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "r",
                                      table_2x2_to_cor = "tetrachoric"), digits = 11)
  # es.mcv_2x2_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
  #                                     table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
  #                                     table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_z3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2", measure = "z",
                                         table_2x2_to_cor = "tetrachoric"), digits = 11)

  expect_equal(unique(es.mcv_2x2_d$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_g$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_or$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_rr$info_used_crude), "2x2")
  # expect_equal(unique(es.mcv_2x2_r1$info_used_crude), "2x2")
  # expect_equal(unique(es.mcv_2x2_r2$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_r3$info_used_crude), "2x2")
  # expect_equal(unique(es.mcv_2x2_z1$info_used_crude), "2x2")
  # expect_equal(unique(es.mcv_2x2_z2$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_z3$info_used_crude), "2x2")

  expect_equal(unique(es.mcv_2x2_d_rv$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_g_rv$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_or_rv$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_rr_rv$info_used_crude), "2x2")
  # expect_equal(unique(es.mcv_2x2_r1_rv$info_used_crude), "2x2")
  # expect_equal(unique(es.mcv_2x2_r2_rv$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_r3_rv$info_used_crude), "2x2")
  # expect_equal(unique(es.mcv_2x2_z1_rv$info_used_crude), "2x2")
  # expect_equal(unique(es.mcv_2x2_z2_rv$info_used_crude), "2x2")
  expect_equal(unique(es.mcv_2x2_z3_rv$info_used_crude), "2x2")

  expect_equal(es.mcv_2x2_d$es_crude, -es.mcv_2x2_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_d$se_crude, es.mcv_2x2_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_g$es_crude, -es.mcv_2x2_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_g$se_crude, es.mcv_2x2_g_rv$se_crude, tolerance = 1e-10)
#
#   expect_equal(es.mcv_2x2_r1$es_crude, -es.mcv_2x2_r1_rv$es_crude, tolerance = 1e-10)
#   expect_equal(es.mcv_2x2_r1$se_crude, es.mcv_2x2_r1_rv$se_crude, tolerance = 1e-10)

  # expect_equal(es.mcv_2x2_z1$es_crude, -es.mcv_2x2_z1_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_z1$se_crude, es.mcv_2x2_z1_rv$se_crude, tolerance = 1e-10)
  #
  # expect_equal(es.mcv_2x2_r2$es_crude, -es.mcv_2x2_r2_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_r2$se_crude, es.mcv_2x2_r2_rv$se_crude, tolerance = 1e-10)

  # expect_equal(es.mcv_2x2_z2$es_crude, -es.mcv_2x2_z2_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_z2$se_crude, es.mcv_2x2_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_r3$es_crude, -es.mcv_2x2_r3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_r3$se_crude, es.mcv_2x2_r3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_z3$es_crude, -es.mcv_2x2_z3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_z3$se_crude, es.mcv_2x2_z3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_or$es_crude, -es.mcv_2x2_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_or$se_crude, es.mcv_2x2_or_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_rr$es_crude, -es.mcv_2x2_rr_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_rr$se_crude, es.mcv_2x2_rr_rv$se_crude, tolerance = 1e-10)
})

test_that("2x2 sum - Reverse", {
  skip_if_not_installed("mvtnorm")
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp

  dat$reverse_2x2 <- FALSE
  es.mcv_2x2_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "d"), digits = 11)
  es.mcv_2x2_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "g"), digits = 11)
  es.mcv_2x2_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "logor"), digits = 11)
  es.mcv_2x2_rr <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "logrr"), digits = 11)
  # es.mcv_2x2_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
  #                                     table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
  #                                     table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_r3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
                                      table_2x2_to_cor = "tetrachoric"), digits = 11)
  # es.mcv_2x2_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
  #                                     table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
  #                                     table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_z3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
                                      table_2x2_to_cor = "tetrachoric"), digits = 11)

  dat$reverse_2x2 <- TRUE
  es.mcv_2x2_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "d"), digits = 11)
  es.mcv_2x2_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "g"), digits = 11)
  es.mcv_2x2_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "logor"), digits = 11)
  es.mcv_2x2_rr_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "logrr"), digits = 11)
  # es.mcv_2x2_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
  #                                        table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
  #                                        table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_r3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "r",
                                         table_2x2_to_cor = "tetrachoric"), digits = 11)
  # es.mcv_2x2_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
  #                                        table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
  #                                        table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_z3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_sum", measure = "z",
                                         table_2x2_to_cor = "tetrachoric"), digits = 11)

  expect_equal(unique(es.mcv_2x2_d$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_g$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_or$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_rr$info_used_crude), "2x2_sum")
  # expect_equal(unique(es.mcv_2x2_r1$info_used_crude), "2x2_sum")
  # expect_equal(unique(es.mcv_2x2_r2$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_r3$info_used_crude), "2x2_sum")
  # expect_equal(unique(es.mcv_2x2_z1$info_used_crude), "2x2_sum")
  # expect_equal(unique(es.mcv_2x2_z2$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_z3$info_used_crude), "2x2_sum")

  expect_equal(unique(es.mcv_2x2_d_rv$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_g_rv$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_or_rv$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_rr_rv$info_used_crude), "2x2_sum")
  # expect_equal(unique(es.mcv_2x2_r1_rv$info_used_crude), "2x2_sum")
  # expect_equal(unique(es.mcv_2x2_r2_rv$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_r3_rv$info_used_crude), "2x2_sum")
  # expect_equal(unique(es.mcv_2x2_z1_rv$info_used_crude), "2x2_sum")
  # expect_equal(unique(es.mcv_2x2_z2_rv$info_used_crude), "2x2_sum")
  expect_equal(unique(es.mcv_2x2_z3_rv$info_used_crude), "2x2_sum")

  expect_equal(es.mcv_2x2_d$es_crude, -es.mcv_2x2_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_d$se_crude, es.mcv_2x2_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_g$es_crude, -es.mcv_2x2_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_g$se_crude, es.mcv_2x2_g_rv$se_crude, tolerance = 1e-10)

  # expect_equal(es.mcv_2x2_r1$es_crude, -es.mcv_2x2_r1_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_r1$se_crude, es.mcv_2x2_r1_rv$se_crude, tolerance = 1e-10)
  #
  # expect_equal(es.mcv_2x2_z1$es_crude, -es.mcv_2x2_z1_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_z1$se_crude, es.mcv_2x2_z1_rv$se_crude, tolerance = 1e-10)
  #
  # expect_equal(es.mcv_2x2_r2$es_crude, -es.mcv_2x2_r2_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_r2$se_crude, es.mcv_2x2_r2_rv$se_crude, tolerance = 1e-10)
  #
  # expect_equal(es.mcv_2x2_z2$es_crude, -es.mcv_2x2_z2_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_z2$se_crude, es.mcv_2x2_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_r3$es_crude, -es.mcv_2x2_r3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_r3$se_crude, es.mcv_2x2_r3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_z3$es_crude, -es.mcv_2x2_z3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_z3$se_crude, es.mcv_2x2_z3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_or$es_crude, -es.mcv_2x2_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_or$se_crude, es.mcv_2x2_or_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_rr$es_crude, -es.mcv_2x2_rr_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_rr$se_crude, es.mcv_2x2_rr_rv$se_crude, tolerance = 1e-10)
})

test_that("prop - Reverse", {
  skip_if_not_installed("mvtnorm")
  dat <- metaumbrella::df.OR
  dat$prop_cases_exp <- dat$n_cases_exp / dat$n_exp
  dat$prop_cases_nexp <- dat$n_cases_nexp / dat$n_nexp

  dat$reverse_prop <- FALSE
  es.mcv_prop_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "d"), digits = 11)
  es.mcv_prop_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "g"), digits = 11)
  es.mcv_prop_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "logor"), digits = 11)
  es.mcv_prop_rr <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "logrr"), digits = 11)
  # es.mcv_2x2_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
  #                                     table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
  #                                     table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_r3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
                                      table_2x2_to_cor = "tetrachoric"), digits = 11)
  # es.mcv_2x2_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
  #                                     table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
  #                                     table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_z3 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
                                      table_2x2_to_cor = "tetrachoric"), digits = 11)

  dat$reverse_prop <- TRUE
  es.mcv_prop_d_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "d"), digits = 11)
  es.mcv_prop_g_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "g"), digits = 11)
  es.mcv_prop_or_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "logor"), digits = 11)
  es.mcv_prop_rr_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "logrr"), digits = 11)
  # es.mcv_2x2_r1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
  #                                        table_2x2_to_cor = "lipsey"), digits = 11)
  # es.mcv_2x2_r2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
  #                                        table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_r3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "r",
                                         table_2x2_to_cor = "tetrachoric"), digits = 11)
  es.mcv_2x2_z1_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
                                         table_2x2_to_cor = "lipsey"), digits = 11)
  es.mcv_2x2_z2_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
                                         table_2x2_to_cor = "cooper"), digits = 11)
  es.mcv_2x2_z3_rv <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "2x2_prop", measure = "z",
                                         table_2x2_to_cor = "tetrachoric"), digits = 11)

  expect_equal(unique(es.mcv_prop_d$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_prop_g$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_prop_or$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_prop_rr$info_used_crude), "2x2_prop")
  # expect_equal(unique(es.mcv_2x2_r1$info_used_crude), "2x2_prop")
  # expect_equal(unique(es.mcv_2x2_r2$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_2x2_r3$info_used_crude), "2x2_prop")
  # expect_equal(unique(es.mcv_2x2_z1$info_used_crude), "2x2_prop")
  # expect_equal(unique(es.mcv_2x2_z2$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_2x2_z3$info_used_crude), "2x2_prop")

  expect_equal(unique(es.mcv_prop_d_rv$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_prop_g_rv$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_prop_or_rv$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_prop_rr_rv$info_used_crude), "2x2_prop")
  # expect_equal(unique(es.mcv_2x2_r1_rv$info_used_crude), "2x2_prop")
  # expect_equal(unique(es.mcv_2x2_r2_rv$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_2x2_r3_rv$info_used_crude), "2x2_prop")
  # expect_equal(unique(es.mcv_2x2_z1_rv$info_used_crude), "2x2_prop")
  # expect_equal(unique(es.mcv_2x2_z2_rv$info_used_crude), "2x2_prop")
  expect_equal(unique(es.mcv_2x2_z3_rv$info_used_crude), "2x2_prop")

  expect_equal(es.mcv_prop_d$es_crude, -es.mcv_prop_d_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_prop_d$se_crude, es.mcv_prop_d_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_prop_g$es_crude, -es.mcv_prop_g_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_prop_g$se_crude, es.mcv_prop_g_rv$se_crude, tolerance = 1e-10)

  # expect_equal(es.mcv_2x2_r1$es_crude, -es.mcv_2x2_r1_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_r1$se_crude, es.mcv_2x2_r1_rv$se_crude, tolerance = 1e-10)
  #
  # expect_equal(es.mcv_2x2_z1$es_crude, -es.mcv_2x2_z1_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_z1$se_crude, es.mcv_2x2_z1_rv$se_crude, tolerance = 1e-10)
  #
  # expect_equal(es.mcv_2x2_r2$es_crude, -es.mcv_2x2_r2_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_r2$se_crude, es.mcv_2x2_r2_rv$se_crude, tolerance = 1e-10)
  #
  # expect_equal(es.mcv_2x2_z2$es_crude, -es.mcv_2x2_z2_rv$es_crude, tolerance = 1e-10)
  # expect_equal(es.mcv_2x2_z2$se_crude, es.mcv_2x2_z2_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_r3$es_crude, -es.mcv_2x2_r3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_r3$se_crude, es.mcv_2x2_r3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_2x2_z3$es_crude, -es.mcv_2x2_z3_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_2x2_z3$se_crude, es.mcv_2x2_z3_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_prop_or$es_crude, -es.mcv_prop_or_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_prop_or$se_crude, es.mcv_prop_or_rv$se_crude, tolerance = 1e-10)

  expect_equal(es.mcv_prop_rr$es_crude, -es.mcv_prop_rr_rv$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_prop_rr$se_crude, es.mcv_prop_rr_rv$se_crude, tolerance = 1e-10)
})
