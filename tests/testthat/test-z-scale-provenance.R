# =============================================================================
# Roadmap item 2.5 -- `measure = "z"` holds two different transforms.
#
# The z column is supposed to hold ONE quantity, Fisher's z, so that it can be
# pooled. It does not:
#
#   pearson_r, fisher_z, or_se, 2x2 ......... z = atanh(r)          (Fisher)
#   means_sd, cohen_d, student_t ............ a variance-stabilising transform
#
# and the second group is the DEFAULT (smd_to_cor = "viechtbauer"). Measured on
# identical data at a point-biserial rho = 0.75, the SMD route reports z = 1.0925
# while atanh() of its own r is 1.7468. A review holding both SMD studies and
# correlation studies -- the ordinary case for a z meta-analysis -- mixes the two
# with no user choice involved.
#
# NOT THE SAME THING as the estimand difference between viechtbauer (which
# estimates the biserial correlation) and lipsey_cooper (the point-biserial).
# That is a documented choice about WHICH correlation to estimate, and study 01 of
# the simulation programme measures it. This file is about the z column not being
# atanh() of the r column on every route.
#
# WHY A FLAG AND NOT JUST DOCUMENTATION, AND WHY ONLY z. The package already
# implements exactly this check for the SMD family -- E6 fires when change-SD rows
# (morris_dz) and raw-score-SD rows share a pool, and E7 when paired t/F rows
# coexist with pooled-standardizer rows. ANCOVA cannot mix silently at all: its
# estimates go to the *adjusted* columns while endpoint means go to the crude ones.
# So z was the one place the pattern was missing, and E8 restores consistency
# rather than singling z out.
# =============================================================================

# A two-group dataset whose point-biserial correlation is `rho`.
.smd_row <- function(rho = 0.5, n = 400) {
  d <- 2 * rho / sqrt(1 - rho^2)
  data.frame(n_exp = n / 2, n_nexp = n / 2, mean_exp = d, mean_nexp = 0,
             mean_sd_exp = 1, mean_sd_nexp = 1)
}

.tokens <- function(res, col = "flags_crude") {
  fl <- res[[col]]
  fl <- fl[!is.na(fl) & nzchar(fl)]
  if (!length(fl)) return(character(0))
  unique(unlist(strsplit(paste(fl, collapse = "; "), "; ")))
}
.e8 <- function(res, col = "flags_crude")
  grep("Mixed z transforms", .tokens(res, col), value = TRUE)

# --- the fact itself --------------------------------------------------------

test_that("the SMD family's z is NOT atanh of its own r, while the correlation family's is", {
  s <- .smd_row(0.5)
  smd <- es_from_means_sd(mean_exp = s$mean_exp, mean_nexp = s$mean_nexp,
                          mean_sd_exp = s$mean_sd_exp, mean_sd_nexp = s$mean_sd_nexp,
                          n_exp = s$n_exp, n_nexp = s$n_nexp)
  expect_false(isTRUE(all.equal(smd$z, atanh(smd$r), tolerance = 1e-8)))

  cor_row <- es_from_pearson_r(pearson_r = 0.5, n_sample = 400)
  expect_equal(cor_row$z, atanh(cor_row$r))

  fz <- es_from_fisher_z(fisher_z = atanh(0.5), n_sample = 400)
  expect_equal(fz$z, atanh(fz$r))

  tt <- es_from_2x2(n_cases_exp = 60, n_controls_exp = 140,
                    n_cases_nexp = 40, n_controls_nexp = 160)
  expect_equal(tt$z, atanh(tt$r))
})

test_that("the two smd_to_cor routes return different z for identical data", {
  for (rho in c(0.25, 0.5, 0.75)) {
    s <- .smd_row(rho)
    v <- es_from_means_sd(mean_exp = s$mean_exp, mean_nexp = s$mean_nexp,
                          mean_sd_exp = s$mean_sd_exp, mean_sd_nexp = s$mean_sd_nexp,
                          n_exp = s$n_exp, n_nexp = s$n_nexp,
                          smd_to_cor = "viechtbauer")
    l <- es_from_means_sd(mean_exp = s$mean_exp, mean_nexp = s$mean_nexp,
                          mean_sd_exp = s$mean_sd_exp, mean_sd_nexp = s$mean_sd_nexp,
                          n_exp = s$n_exp, n_nexp = s$n_nexp,
                          smd_to_cor = "lipsey_cooper")
    expect_false(isTRUE(all.equal(v$z, l$z, tolerance = 1e-6)),
                 info = paste("rho =", rho))
    # lipsey_cooper's z IS Fisher's z of its own r; viechtbauer's is not
    expect_equal(l$z, atanh(l$r))
    expect_false(isTRUE(all.equal(v$z, atanh(v$r), tolerance = 1e-8)))
  }
  # the gap grows with rho -- 0.0021 / 0.0170 / 0.1195 at 0.25 / 0.50 / 0.75
  gap <- function(rho) {
    s <- .smd_row(rho)
    v <- es_from_means_sd(mean_exp = s$mean_exp, mean_nexp = s$mean_nexp,
                          mean_sd_exp = s$mean_sd_exp, mean_sd_nexp = s$mean_sd_nexp,
                          n_exp = s$n_exp, n_nexp = s$n_nexp, smd_to_cor = "viechtbauer")
    l <- es_from_means_sd(mean_exp = s$mean_exp, mean_nexp = s$mean_nexp,
                          mean_sd_exp = s$mean_sd_exp, mean_sd_nexp = s$mean_sd_nexp,
                          n_exp = s$n_exp, n_nexp = s$n_nexp, smd_to_cor = "lipsey_cooper")
    abs(v$z - l$z)
  }
  expect_lt(gap(0.25), gap(0.50))
  expect_lt(gap(0.50), gap(0.75))
})

# --- the provenance marker --------------------------------------------------

test_that(".z_transform_by_route() classifies from the output, not a route list", {
  d <- cbind(.smd_row(0.5), pearson_r = NA_real_, n_sample = NA_real_)
  d <- rbind(d, data.frame(n_exp = NA, n_nexp = NA, mean_exp = NA, mean_nexp = NA,
                           mean_sd_exp = NA, mean_sd_nexp = NA,
                           pearson_r = 0.5, n_sample = 400))
  o <- suppressMessages(convert_df(d, measure = "z"))
  zt <- attr(o, "z_transform")

  expect_true(is.character(zt))
  expect_true(length(zt) >= 2L)
  expect_equal(unname(zt["means_sd"]), "vst")
  expect_equal(unname(zt["pearson_r"]), "fisher")
  # nothing is classified "mixed" on well-behaved data
  expect_false(any(zt == "mixed", na.rm = TRUE))
  # names are info_used values -- the key summary() reports -- not frame names
  expect_false(any(grepl("^es_", names(zt))))
})

test_that("the marker follows smd_to_cor rather than being fixed per route", {
  d <- .smd_row(0.5)
  v <- attr(suppressMessages(convert_df(d, measure = "z", smd_to_cor = "viechtbauer")),
            "z_transform")
  l <- attr(suppressMessages(convert_df(d, measure = "z", smd_to_cor = "lipsey_cooper")),
            "z_transform")
  expect_equal(unname(v["means_sd"]), "vst")
  expect_equal(unname(l["means_sd"]), "fisher")
})

# --- E8, the cross-row flag -------------------------------------------------

.mixed_df <- function() data.frame(
  n_exp = c(200, NA), n_nexp = c(200, NA),
  mean_exp = c(1.15, NA), mean_nexp = c(0, NA),
  mean_sd_exp = c(1, NA), mean_sd_nexp = c(1, NA),
  pearson_r = c(NA, 0.5), n_sample = c(NA, 400))

test_that("E8 fires on a pool mixing the two transforms, and names both sides", {
  res <- suppressMessages(summary(suppressMessages(
    convert_df(.mixed_df(), measure = "z")), flags = TRUE))
  hits <- .e8(res)
  expect_length(hits, 2L)                       # one message per side
  expect_true(any(grepl("variance-stabilising transform, not atanh", hits)))
  expect_true(any(grepl("Fisher's z, atanh\\(r\\)", hits)))
  expect_true(all(grepl("^\\[INFO\\]", hits)))  # a pool property, not a row defect
  expect_true(any(grepl("from means_sd", hits, fixed = TRUE)))
  expect_true(any(grepl("from pearson_r", hits, fixed = TRUE)))
  # and it points at the fix
  expect_true(all(grepl("lipsey_cooper", hits, fixed = TRUE)))
})

test_that("E8 is silent when the pool is uniform", {
  only_cor <- data.frame(pearson_r = c(0.5, 0.4), n_sample = c(400, 300))
  r1 <- suppressMessages(summary(suppressMessages(
    convert_df(only_cor, measure = "z")), flags = TRUE))
  expect_length(.e8(r1), 0L)

  only_smd <- data.frame(n_exp = c(200, 150), n_nexp = c(200, 150),
                         mean_exp = c(1.15, 0.9), mean_nexp = c(0, 0),
                         mean_sd_exp = c(1, 1), mean_sd_nexp = c(1, 1))
  r2 <- suppressMessages(summary(suppressMessages(
    convert_df(only_smd, measure = "z")), flags = TRUE))
  expect_length(.e8(r2), 0L)
})

test_that("E8 goes quiet once the pool is put on one transform", {
  # The remedy the message recommends must actually work.
  res <- suppressMessages(summary(suppressMessages(
    convert_df(.mixed_df(), measure = "z", smd_to_cor = "lipsey_cooper")), flags = TRUE))
  expect_length(.e8(res), 0L)
  zt <- attr(suppressMessages(
    convert_df(.mixed_df(), measure = "z", smd_to_cor = "lipsey_cooper")), "z_transform")
  expect_true(all(zt[c("means_sd", "pearson_r")] == "fisher"))
})

test_that("E8 respects the cross-row toggle", {
  res <- suppressMessages(summary(suppressMessages(
    convert_df(.mixed_df(), measure = "z")), flags = TRUE,
    flag_options = list(enable_cross_row = FALSE)))
  expect_length(.e8(res), 0L)
})

test_that("E8 does not fire for measures other than z", {
  # The r column is a correlation on every route; only z mixes transforms.
  res <- suppressMessages(summary(suppressMessages(
    convert_df(.mixed_df(), measure = "r")), flags = TRUE))
  expect_length(.e8(res), 0L)
})

test_that("E8 survives the route view (main_es = FALSE)", {
  # Cross-row checks decide on the route that would be SELECTED; under the route
  # view a comparison occupies several rows and must not flag itself.
  res <- suppressMessages(summary(suppressMessages(
    convert_df(.mixed_df(), measure = "z")), flags = TRUE))
  expect_length(.e8(res), 2L)
  only_smd <- data.frame(n_exp = 200, n_nexp = 200, mean_exp = 1.15,
                         mean_nexp = 0, mean_sd_exp = 1, mean_sd_nexp = 1)
  res2 <- suppressMessages(summary(suppressMessages(
    convert_df(only_smd, measure = "z")), flags = TRUE))
  expect_length(.e8(res2), 0L)
})

# --- the reason this is E8 and not a new idea --------------------------------

test_that("the SMD families the same argument applies to are already flagged (E6)", {
  # Guards the consistency claim in this file's header: if E6 ever stopped firing,
  # z would indeed be singled out and this file's premise would be wrong.
  d <- data.frame(
    n_exp = c(50, 50), n_nexp = c(50, 50),
    mean_change_exp = c(2, NA), mean_change_nexp = c(0, NA),
    mean_change_sd_exp = c(3, NA), mean_change_sd_nexp = c(3, NA),
    r_pre_post_exp = c(0.8, NA), r_pre_post_nexp = c(0.8, NA),
    mean_exp = c(NA, 12), mean_nexp = c(NA, 10),
    mean_sd_exp = c(NA, 3), mean_sd_nexp = c(NA, 3))
  res <- suppressMessages(summary(suppressMessages(
    convert_df(d, measure = "g", pre_post_to_smd = "morris_dz")), flags = TRUE))
  expect_gte(length(grep("Mixed SMD standardizers", .tokens(res))), 2L)
})

test_that("ANCOVA cannot mix with endpoint means: they are in different columns", {
  # The other family named when this item was scoped. ANCOVA estimates are
  # ADJUSTED and never share a column with crude endpoint estimates, so there is
  # no silent pooling for a flag to catch.
  d <- data.frame(
    n_exp = c(50, 50), n_nexp = c(50, 50),
    ancova_mean_exp = c(12, NA), ancova_mean_nexp = c(10, NA),
    ancova_mean_sd_exp = c(3, NA), ancova_mean_sd_nexp = c(3, NA),
    cov_outcome_r = c(0.5, NA), n_cov_ancova = c(1, NA),
    mean_exp = c(NA, 12), mean_nexp = c(NA, 10),
    mean_sd_exp = c(NA, 3), mean_sd_nexp = c(NA, 3))
  res <- suppressMessages(summary(suppressMessages(convert_df(d, measure = "g"))))
  expect_true(is.na(res$es_crude[1]))       # the ANCOVA row is not crude ...
  expect_false(is.na(res$es_adjusted[1]))   # ... it is adjusted
  expect_false(is.na(res$es_crude[2]))      # the endpoint row is crude
  expect_true(is.na(res$es_adjusted[2]))
})
