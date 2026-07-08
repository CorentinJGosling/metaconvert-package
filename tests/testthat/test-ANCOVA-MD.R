## ANC MEANS v T (raw data) ---------------------------------

test_that("ancova means and ancova t converge", {
  x = rep(c("A", "B"), each=25)
  c = rnorm(50)
  y = rnorm(50) + c
  r = cor.test(~c+y)$estimate
  res = lm(y~0+x+c)
  dat=data.frame(study=1)
  dat$ancova_mean_exp = summary(res)$coefficients[1,1]
  dat$ancova_mean_se_exp = summary(res)$coefficients[1,2]
  dat$ancova_mean_nexp = summary(res)$coefficients[2,1]
  dat$ancova_mean_se_nexp = summary(res)$coefficients[2,2]
  dat$n_exp=25; dat$n_nexp=25
  dat$cov_outcome_r=r
  dat$n_cov_ancova = 1

  res2 = lm(y~x+c)
  dat$ancova_t = summary(res2)$coefficients[2,3]

  es.mcv_d_t <- summary(convert_df(dat,
                                   verbose = FALSE,
                                   es_selected = "hierarchy",
                                   hierarchy = "ancova_t",
                                   measure = "d",
                                   pre_post_to_smd = "cooper"
  ), digits = 11)
  es.mcv_d <- summary(convert_df(dat,
                                 verbose = FALSE,
                                 es_selected = "hierarchy", hierarchy = "ancova_means_se",
                                 measure = "d",
                                 pre_post_to_smd = "cooper"
  ), digits = 11)
  ## d
  expect_equal(unique(es.mcv_d$info_used_adjusted), "ancova_means_se")
  expect_equal(unique(es.mcv_d_t$info_used_adjusted), "ancova_t")
  expect_equal(abs(es.mcv_d$es_adjusted), abs(es.mcv_d_t$es_adjusted),
               tolerance = 1e-1)
  expect_equal(es.mcv_d$se_adjusted, es.mcv_d_t$se_adjusted, tolerance = 1e-1)
})

## ANC MD v T (raw data) ---------------------------------

test_that("ancova MD+SE and ancova t converge", {
  x = rep(c("A", "B"), each=25)
  c = rnorm(50)
  y = rnorm(50) + c
  r = cor.test(~c+y)$estimate
  res = lm(y~x+c)

  dat=data.frame(study=1)
  dat$ancova_md = summary(res)$coefficients[2,1]
  dat$ancova_md_se = summary(res)$coefficients[2,2]
  dat$n_exp=25; dat$n_nexp=25
  dat$cov_outcome_r=r
  dat$n_cov_ancova = 1
  dat$ancova_t = summary(res)$coefficients[2,3]

  es.mcv_d_t <- summary(convert_df(dat,
                                   verbose = FALSE,
                                   es_selected = "hierarchy", hierarchy = "ancova_t",
                                   measure = "d",
                                   pre_post_to_smd = "cooper"
  ), digits = 11)
  es.mcv_d <- summary(convert_df(dat,
                                 verbose = FALSE,
                                 es_selected = "hierarchy", hierarchy = "ancova_md_se",
                                 measure = "d",
                                 pre_post_to_smd = "cooper"
  ), digits = 11)
  ## d
  expect_equal(unique(es.mcv_d$info_used_adjusted), "ancova_md_se")
  expect_equal(unique(es.mcv_d_t$info_used_adjusted), "ancova_t")
  expect_equal(abs(es.mcv_d$es_adjusted), abs(es.mcv_d_t$es_adjusted),
               tolerance = 1e-10)
  expect_equal(es.mcv_d$se_adjusted, es.mcv_d_t$se_adjusted, tolerance = 1e-10)
})

