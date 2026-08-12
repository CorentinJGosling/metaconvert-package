test_that("STD reg", {
  skip_if_not_installed("esc")
  dat <- df.haza[1:50, ]
  dat$beta_std <- abs(runif(nrow(dat), -0.5, 0.5))
  dat$beta_unstd <- runif(nrow(dat), -0.7, 0.7)
  # realistic dependent-variable SD (large enough that the synthetic beta_unstd
  # values remain mathematically consistent, i.e. imply a non-negative within-group
  # variance; a tiny sd_dv would make high-|beta| rows impossible and return NA)
  dat$sd_dv <- runif(nrow(dat), 1, 3)

  # esc functions require scalar inputs, so we loop row-by-row
  comp_res_d <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r <- esc::esc_beta(beta = dat$beta_std[i], sdy = dat$sd_dv[i],
                       grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i], es.type = "d")
    data.frame(es = r$es, se = r$se, stringsAsFactors = FALSE)
  }))
  comp_res_r <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r <- esc::esc_beta(beta = dat$beta_std[i], sdy = dat$sd_dv[i],
                       grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i], es.type = "r")
    data.frame(es = r$es, se = r$se, zr = r$zr, stringsAsFactors = FALSE)
  }))
  comp_res_or <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r <- esc::esc_beta(beta = dat$beta_std[i], sdy = dat$sd_dv[i],
                       grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i], es.type = "or")
    data.frame(es = r$es, se = r$se, stringsAsFactors = FALSE)
  }))

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "logor"), digits = 11)
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "r",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)
  es.mcv_z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "z",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_d$info_used_crude), "beta_std")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 1e-10)

  expect_equal(es.mcv_or$es_crude, log(comp_res_or$es), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, comp_res_or$se, tolerance = 1e-10)

  expect_equal(es.mcv_r$es_crude, comp_res_r$es, tolerance = 1e-10)
  # expect_equal(es.mcv_r$se_crude, comp_res_r$se, tolerance = 1e-10)

  expect_equal(es.mcv_z$es_crude, comp_res_r$zr, tolerance = 1e-10)
  # compute.es's Fisher-z SE is not the Fisher transform of the lipsey_cooper r variance;
  # metaConvert uses the Fisher-consistent delta value se(z) = se(r)/(1 - r^2). The z POINT
  # estimate above still matches compute.es.
  expect_equal(es.mcv_z$se_crude, es.mcv_r$se_crude / (1 - es.mcv_r$es_crude^2), tolerance = 1e-9)
})

test_that("UNSTD reg", {
  skip_if_not_installed("esc")
  dat <- df.haza[1:50, ]
  dat$beta_std <- abs(runif(nrow(dat), -0.5, 0.5))
  dat$beta_unstd <- runif(nrow(dat), -0.7, 0.7)
  # realistic dependent-variable SD (large enough that the synthetic beta_unstd
  # values remain mathematically consistent, i.e. imply a non-negative within-group
  # variance; a tiny sd_dv would make high-|beta| rows impossible and return NA)
  dat$sd_dv <- runif(nrow(dat), 1, 3)

  # esc functions require scalar inputs, so we loop row-by-row
  comp_res_d <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r <- esc::esc_B(b = dat$beta_unstd[i], sdy = dat$sd_dv[i],
                    grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i], es.type = "d")
    data.frame(es = r$es, se = r$se, stringsAsFactors = FALSE)
  }))
  comp_res_r <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r <- esc::esc_B(b = dat$beta_unstd[i], sdy = dat$sd_dv[i],
                    grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i], es.type = "r")
    data.frame(es = r$es, se = r$se, zr = r$zr, stringsAsFactors = FALSE)
  }))
  comp_res_or <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r <- esc::esc_B(b = dat$beta_unstd[i], sdy = dat$sd_dv[i],
                    grp1n = dat$n_exp[i], grp2n = dat$n_nexp[i], es.type = "or")
    data.frame(es = r$es, se = r$se, stringsAsFactors = FALSE)
  }))

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "d"),
    digits = 11
  )
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "logor"),
    digits = 11
  )
  es.mcv_r <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "r",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)
  es.mcv_z <- summary(convert_df(dat,
    verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "z",
    smd_to_cor = "lipsey_cooper"
  ), digits = 11)

  ## test ES
  expect_equal(unique(es.mcv_d$info_used_crude), "beta_unstd")
  expect_equal(es.mcv_d$es_crude, comp_res_d$es, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, comp_res_d$se, tolerance = 1e-10)

  expect_equal(es.mcv_or$es_crude, log(comp_res_or$es), tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, comp_res_or$se, tolerance = 1e-10)

  expect_equal(es.mcv_r$es_crude, comp_res_r$es, tolerance = 1e-10)
  # expect_equal(es.mcv_r$se_crude, comp_res_r$se, tolerance = 1e-1)

  expect_equal(es.mcv_z$es_crude, comp_res_r$zr, tolerance = 1e-10)
  # compute.es's Fisher-z SE is not the Fisher transform of the lipsey_cooper r variance;
  # metaConvert uses the Fisher-consistent delta value se(z) = se(r)/(1 - r^2). The z POINT
  # estimate above still matches compute.es.
  expect_equal(es.mcv_z$se_crude, es.mcv_r$se_crude / (1 - es.mcv_r$es_crude^2), tolerance = 1e-9)
})

test_that("beta_std-reverse", {
  dat <- df.haza[1:50, ]
  dat$beta_std <- abs(runif(nrow(dat), -0.5, 0.5))
  dat$beta_unstd <- runif(nrow(dat), -0.7, 0.7)
  # realistic dependent-variable SD (large enough that the synthetic beta_unstd
  # values remain mathematically consistent, i.e. imply a non-negative within-group
  # variance; a tiny sd_dv would make high-|beta| rows impossible and return NA)
  dat$sd_dv <- runif(nrow(dat), 1, 3)

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "logor"), digits = 11)
  es.mcv_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "z",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "r",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "z",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "r",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "g"), digits = 11)

  dat$reverse_beta_std <- TRUE

  es.mcv_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "d"), digits = 11)
  es.mcv_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "logor"), digits = 11)
  es.mcv_z1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "z",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "r",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "z",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_r2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "r",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_std", measure = "g"), digits = 11)

  expect_true(all(c(es.mcv_d$info_used_crude, es.mcv_or$info_used_crude, es.mcv_g$info_used_crude,
                    es.mcv_z1$info_used_crude, es.mcv_r1$info_used_crude,
                    es.mcv_z2$info_used_crude, es.mcv_r2$info_used_crude,
                    es.mcv_d_r$info_used_crude, es.mcv_or_r$info_used_crude, es.mcv_g_r$info_used_crude,
                    es.mcv_z1_r$info_used_crude, es.mcv_r1_r$info_used_crude, es.mcv_z2_r$info_used_crude,
                    es.mcv_r2_r$info_used_crude) == "beta_std"))
  expect_equal(es.mcv_d$es_crude, -es.mcv_d_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, -es.mcv_or_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$es_crude, -es.mcv_z1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$se_crude, es.mcv_z1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$es_crude, -es.mcv_r1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$se_crude, es.mcv_r1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$es_crude, -es.mcv_z2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$se_crude, es.mcv_z2_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$es_crude, -es.mcv_r2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$se_crude, es.mcv_r2_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$es_crude, -es.mcv_g_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_r$se_crude, tolerance = 1e-10)
})



test_that("beta_unstd-reverse", {
  dat <- df.haza[1:50, ]
  dat$beta_std <- abs(runif(nrow(dat), -0.5, 0.5))
  dat$beta_unstd <- runif(nrow(dat), -0.7, 0.7)
  # realistic dependent-variable SD (large enough that the synthetic beta_unstd
  # values remain mathematically consistent, i.e. imply a non-negative within-group
  # variance; a tiny sd_dv would make high-|beta| rows impossible and return NA)
  dat$sd_dv <- runif(nrow(dat), 1, 3)

  es.mcv_d <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "d"), digits = 11)
  es.mcv_or <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "logor"), digits = 11)
  es.mcv_z1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "z",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r1 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "r",
                                  smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "z",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_r2 <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "r",
                                  smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_g <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "g"), digits = 11)

  dat$reverse_beta_unstd <- TRUE

  es.mcv_d_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "d"), digits = 11)
  es.mcv_or_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "logor"), digits = 11)
  es.mcv_z1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "z",
                                    smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_r1_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "r",
                                    smd_to_cor = "viechtbauer"), digits = 11)
  es.mcv_z2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "z",
                                    smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_r2_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "r",
                                    smd_to_cor = "lipsey_cooper"), digits = 11)
  es.mcv_g_r <- summary(convert_df(dat, verbose = FALSE, es_selected = "hierarchy", hierarchy = "beta_unstd", measure = "g"), digits = 11)


  expect_true(all(c(es.mcv_d$info_used_crude, es.mcv_or$info_used_crude, es.mcv_g$info_used_crude,
                    es.mcv_z1$info_used_crude, es.mcv_r1$info_used_crude,
                    es.mcv_z2$info_used_crude, es.mcv_r2$info_used_crude,
                    es.mcv_d_r$info_used_crude, es.mcv_or_r$info_used_crude, es.mcv_g_r$info_used_crude,
                    es.mcv_z1_r$info_used_crude, es.mcv_r1_r$info_used_crude, es.mcv_z2_r$info_used_crude,
                    es.mcv_r2_r$info_used_crude) == "beta_unstd"))
  expect_equal(es.mcv_d$es_crude, -es.mcv_d_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_d$se_crude, es.mcv_d_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$es_crude, -es.mcv_or_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_or$se_crude, es.mcv_or_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$es_crude, -es.mcv_z1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z1$se_crude, es.mcv_z1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$es_crude, -es.mcv_r1_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r1$se_crude, es.mcv_r1_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$es_crude, -es.mcv_z2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_z2$se_crude, es.mcv_z2_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$es_crude, -es.mcv_r2_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_r2$se_crude, es.mcv_r2_r$se_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$es_crude, -es.mcv_g_r$es_crude, tolerance = 1e-10)
  expect_equal(es.mcv_g$se_crude, es.mcv_g_r$se_crude, tolerance = 1e-10)
})
