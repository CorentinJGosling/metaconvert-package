# Per-arm direction flags on the paired-statistic routes.
#
# A paired t or F carries no sign, so |d_arm| is all that is recoverable and the
# contrast d = d_exp - d_nexp silently assumes both arms moved the SAME way. The
# per-arm flags exist to state otherwise, and the whole-contrast flag cannot express
# it: when treatment improves and control deteriorates the true contrast is
# d_exp + d_nexp, which is not +-(d_exp - d_nexp) for any flag setting.
#
# These routes shipped with no test at all, and the six columns were unreachable from
# convert_df() until they were registered in .check_data() and data_extraction_sheet().

test_that("per-arm reversal generalises the whole-contrast flag", {
  n <- 40; r <- 0.5
  g <- function(...) es_from_paired_f(paired_f_exp = 9, paired_f_nexp = 4,
                                      n_exp = n, n_nexp = n,
                                      r_pre_post_exp = r, r_pre_post_nexp = r, ...)$d
  de <- sqrt(9) * sqrt(2 * (1 - r) / n)   # |d| of the exposed arm
  dn <- sqrt(4) * sqrt(2 * (1 - r) / n)   # |d| of the control arm

  # the four sign combinations, each derived by hand
  expect_equal(g(),                                de - dn, tolerance = 1e-12)
  expect_equal(g(reverse_paired_f = TRUE),      -(de - dn), tolerance = 1e-12)
  expect_equal(g(reverse_paired_f_exp = TRUE),  -de - dn,   tolerance = 1e-12)
  expect_equal(g(reverse_paired_f_nexp = TRUE),  de + dn,   tolerance = 1e-12)

  # THE INVARIANT: flipping both arms is the same statement as flipping the contrast.
  # If this ever fails, one of the two mechanisms has picked up a sign of its own.
  expect_equal(g(reverse_paired_f_exp = TRUE, reverse_paired_f_nexp = TRUE),
               g(reverse_paired_f = TRUE), tolerance = 1e-12)

  # ...and the arms-diverged value is reachable by NO setting of the whole-contrast
  # flag, which is the reason the per-arm flags exist.
  expect_false(isTRUE(all.equal(g(reverse_paired_f_nexp = TRUE), g())))
  expect_false(isTRUE(all.equal(g(reverse_paired_f_nexp = TRUE),
                                g(reverse_paired_f = TRUE))))

  # the SE is a magnitude and must not move with any of them
  se <- function(...) es_from_paired_f(paired_f_exp = 9, paired_f_nexp = 4,
                                       n_exp = n, n_nexp = n,
                                       r_pre_post_exp = r, r_pre_post_nexp = r, ...)$d_se
  expect_equal(se(reverse_paired_f_exp = TRUE),  se(), tolerance = 1e-12)
  expect_equal(se(reverse_paired_f_nexp = TRUE), se(), tolerance = 1e-12)
  expect_equal(se(reverse_paired_f = TRUE),      se(), tolerance = 1e-12)
})

test_that("the per-arm flags are reachable from convert_df(), not only from the route", {
  n <- 40; r <- 0.5
  base <- data.frame(study_id = "s1", n_exp = n, n_nexp = n,
                     r_pre_post_exp = r, r_pre_post_nexp = r,
                     paired_f_exp = 9, paired_f_nexp = 4)
  pull <- function(d) summary(convert_df(d, measure = "d", verbose = FALSE,
                                         es_selected = "hierarchy",
                                         hierarchy = "paired_f"))$es_crude

  with_flag <- base; with_flag$reverse_paired_f_nexp <- TRUE
  expect_false(isTRUE(all.equal(pull(base), pull(with_flag))))
  expect_equal(pull(with_flag),
               es_from_paired_f(paired_f_exp = 9, paired_f_nexp = 4,
                                n_exp = n, n_nexp = n,
                                r_pre_post_exp = r, r_pre_post_nexp = r,
                                reverse_paired_f_nexp = TRUE)$d,
               tolerance = 1e-10)

  # all six columns must survive .check_data() as logicals and be in the sheet
  six <- c("reverse_paired_t_pval_exp", "reverse_paired_t_pval_nexp",
           "reverse_paired_f_exp", "reverse_paired_f_nexp",
           "reverse_paired_f_pval_exp", "reverse_paired_f_pval_nexp")
  expect_true(all(six %in% names(data_extraction_sheet("all"))))
  d <- base
  for (k in six) d[[k]] <- TRUE
  cd <- metaConvert:::.check_data(d)
  expect_true(all(six %in% names(cd)))
  expect_true(all(vapply(six, function(k) is.logical(cd[[k]]), logical(1))))
})
