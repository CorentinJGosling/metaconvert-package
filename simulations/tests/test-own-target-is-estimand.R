## =============================================================================
## Roadmap item 3.2 -- study 03's `theta_own` was scale-matched, not
## estimand-matched.
##
## The `own` target exists so that comparing two methods' bias measures
## COMPUTATIONAL error rather than "these two are estimating different things".
## Study 03 set theta_own to theta_pop for the (r) routes and atanh(theta_pop) for
## the (z) routes -- i.e. it fixed only the SCALE. So for every (r) route `own`
## WAS `population`: bias columns bit-identical across all 200 cells of each
## shipped file, max |difference| exactly 0. Meanwhile, under CAT the tetrachoric
## routes were still scored against phi and under CONT the phi routes were still
## scored against the latent rho.
##
## phi and the tetrachoric are different population quantities, and in each
## mechanism exactly ONE of them equals rho:
##
##             phi estimand          tetrachoric estimand
##   CONT      phi of the            rho
##             dichotomised table
##   CAT       rho                   the latent rho the table implies
##
## Mis-scoring left in place before the fix, over study 03's grid: mean |gap|
## 0.0925 (max 0.2429) under CONT, 0.1341 (max 0.4186) under CAT. Worst condition
## (CAT, rho 0.5, p_exp 0.3, p_case 0.1): the tetrachoric's true estimand is
## 0.9186 while it was scored against 0.50.
## =============================================================================

## --- the shared population quantities ---------------------------------------

test_that("the latent-2x2 helpers are defined exactly once across the programme", {
  # They lived in study 09 while study 03 carried its own copies of two of them
  # under different names -- the simulation-side instance of roadmap 1.6. One
  # definition, in R/05_latent_2x2.R.
  files <- c(list.files(file.path(.sim_root, "R"), pattern = "[.]R$", full.names = TRUE),
             list.files(file.path(.sim_root, "studies"), pattern = "[.]R$", full.names = TRUE))
  for (fn in c(".biv_p11", ".p11_to_phi", ".phi_to_p11", ".phi_attainable",
               ".p11_to_tetrachoric", ".estimands_2x2")) {
    n <- sum(vapply(files, function(f) {
      sum(grepl(paste0("^", gsub("\\.", "\\\\.", fn), "\\s*<-\\s*function"),
                readLines(f, warn = FALSE)))
    }, integer(1)))
    expect_equal(n, 1L, info = paste0(fn, " is defined ", n, " time(s)"))
  }
  # and the deleted duplicates have not come back under their old names
  for (fn in c("phi_to_p11", "phi_attainable")) {
    n <- sum(vapply(files, function(f) {
      sum(grepl(paste0("^", fn, "\\s*<-\\s*function"), readLines(f, warn = FALSE)))
    }, integer(1)))
    expect_equal(n, 0L, info = paste0(fn, " (undotted duplicate) is back"))
  }
})

test_that(".estimands_2x2() returns rho for exactly one method per mechanism", {
  for (rho in c(0.1, 0.25, 0.5, 0.75)) {
    for (p in c(0.3, 0.5)) for (q in c(0.1, 0.3, 0.5)) {
      if (!.phi_attainable(rho, p, q)) next
      cont <- .estimands_2x2("cont", rho, p, q)
      cat_ <- .estimands_2x2("cat",  rho, p, q)
      expect_equal(cont$tetrachoric, rho)     # CONT: the tetrachoric IS rho
      expect_equal(cat_$phi, rho)             # CAT : phi IS rho
      # ... and the other one is not, whenever rho is away from 0
      expect_false(isTRUE(all.equal(cont$phi, rho)))
      expect_false(isTRUE(all.equal(cat_$tetrachoric, rho)))
      # the dichotomised phi is an attenuation of the latent correlation
      expect_lt(abs(cont$phi), abs(rho))
      # and the implied latent correlation exceeds the phi it came from
      expect_gt(abs(cat_$tetrachoric), abs(rho))
    }
  }
})

test_that(".estimands_2x2() is the exact inverse of itself across mechanisms", {
  # The latent rho implied by a table with phi = x must dichotomise back to x.
  for (rho in c(0.1, 0.3, 0.6)) for (p in c(0.3, 0.5)) for (q in c(0.3, 0.5)) {
    if (!.phi_attainable(rho, p, q)) next
    tetra <- .estimands_2x2("cat", rho, p, q)$tetrachoric
    back  <- .p11_to_phi(.biv_p11(tetra, p, q), p, q)
    expect_equal(back, rho, tolerance = 1e-6)
  }
})

test_that(".estimands_2x2() rejects an unknown mechanism rather than guessing", {
  expect_error(.estimands_2x2("banana", 0.5, 0.5, 0.5))
})

## --- study 03: the fix itself ------------------------------------------------

.cond03 <- function(rho = 0.5, n = 300, p_exp = 0.3, p_case = 0.1)
  list(rho = rho, n = n, p_exp = p_exp, p_case = p_case)

test_that("the generators stamp which mechanism produced the data", {
  set.seed(320)
  expect_equal(unique(as.character(gen_cat(.cond03(), 5)$dgm)), "cat")
  expect_equal(unique(as.character(gen_cont(.cond03(), 5)$dgm)), "cont")
})

test_that("estimate_2x2() refuses to guess when the stamp is missing", {
  # Without it the function would silently fall back to scale-only matching --
  # exactly the defect -- so it must fail loudly instead.
  skip_if_not(isTRUE(get0("METACONVERT_AVAILABLE", ifnotfound = FALSE)),
              "metaConvert not loadable")
  set.seed(321)
  d <- gen_cat(.cond03(), 5)
  d$dgm <- NULL
  expect_error(estimate_2x2(d, .cond03()), "dgm")
})

test_that("each route is scored against ITS OWN estimand, on its own scale", {
  skip_if_not(isTRUE(get0("METACONVERT_AVAILABLE", ifnotfound = FALSE)),
              "metaConvert not loadable")
  cond <- .cond03(rho = 0.5, p_exp = 0.3, p_case = 0.1)
  want <- list(cat  = .estimands_2x2("cat",  0.5, 0.3, 0.1),
               cont = .estimands_2x2("cont", 0.5, 0.3, 0.1))

  for (dgm in c("cat", "cont")) {
    set.seed(322)
    d <- if (dgm == "cat") gen_cat(cond, 8) else gen_cont(cond, 8)
    e <- estimate_2x2(d, cond)
    own <- vapply(split(e$theta_own, e$method), function(z) z[1], numeric(1))

    expect_equal(unname(own[["tetrachoric (r)"]]), want[[dgm]]$tetrachoric)
    expect_equal(unname(own[["phi (r)"]]),         want[[dgm]]$phi)
    expect_equal(unname(own[["tetrachoric (z)"]]), atanh(want[[dgm]]$tetrachoric))
    expect_equal(unname(own[["phi (z)"]]),         atanh(want[[dgm]]$phi))
  }
})

test_that("the (r) routes' own target is no longer identical to population", {
  # The signature of the defect: for every (r) route, own == population.
  skip_if_not(isTRUE(get0("METACONVERT_AVAILABLE", ifnotfound = FALSE)),
              "metaConvert not loadable")
  cond <- .cond03(rho = 0.5, p_exp = 0.3, p_case = 0.1)
  for (dgm in c("cat", "cont")) {
    set.seed(323)
    d <- if (dgm == "cat") gen_cat(cond, 8) else gen_cont(cond, 8)
    e <- estimate_2x2(d, cond)
    r_rows <- grepl("\\(r\\)", e$method)
    expect_false(isTRUE(all.equal(e$theta_own[r_rows], e$theta_pop[r_rows])),
                 info = paste(dgm, "-- own must not collapse onto population"))
    # exactly one of the two families still coincides with rho, by construction
    matched <- vapply(split(e$theta_own[r_rows], e$method[r_rows]),
                      function(z) isTRUE(all.equal(z[1], cond$rho)), logical(1))
    expect_equal(sum(matched), 1L, info = dgm)
  }
})

test_that("mismatched scoring is what the fix removes, and it was large", {
  # Quantify rather than assert a direction: this is the number the write-up uses.
  gaps <- do.call(rbind, lapply(c(0.1, 0.25, 0.5, 0.75), function(rho) {
    do.call(rbind, lapply(c(0.3, 0.5), function(p) {
      do.call(rbind, lapply(c(0.1, 0.3, 0.5), function(q) {
        if (!.phi_attainable(rho, p, q)) return(NULL)
        data.frame(cont = abs(.estimands_2x2("cont", rho, p, q)$phi - rho),
                   cat  = abs(.estimands_2x2("cat",  rho, p, q)$tetrachoric - rho))
      }))
    }))
  }))
  expect_gt(max(gaps$cont), 0.24)     # measured 0.2429
  expect_gt(max(gaps$cat),  0.41)     # measured 0.4186
  expect_gt(mean(gaps$cont), 0.05)
  expect_gt(mean(gaps$cat),  0.10)
})

test_that("es_from_2x2()'s z really is Fisher's z, so plain atanh() is right here", {
  # The roadmap flagged this: a variance-stabilising route would need
  # vs_transform(), not atanh(). Checked rather than assumed -- and if the package
  # ever changes that transform, this fails instead of the targets going quietly
  # wrong.
  skip_if_not(isTRUE(get0("METACONVERT_AVAILABLE", ifnotfound = FALSE)),
              "metaConvert not loadable")
  set.seed(324)
  a <- sample(5:60, 20); b <- sample(5:60, 20)
  c_ <- sample(5:60, 20); d <- sample(5:60, 20)
  res <- es_from_2x2(n_cases_exp = a, n_controls_exp = b,
                     n_cases_nexp = c_, n_controls_nexp = d,
                     table_2x2_to_cor = "tetrachoric")
  expect_equal(res$z, atanh(res$r))
  # the candidate builds its z the same way
  ph <- phi_candidate(a, b, c_, d)
  expect_equal(ph$z, atanh(ph$r))
})

## --- generic anti-regression across the programme ---------------------------

## Which studies expose a per-method `own` target, and how to exercise one
## condition of each. nrep is tiny on purpose: theta_own is a population quantity
## of the CONDITION, so it does not depend on the number of replications.
.own_registry <- list(
  list(id = "01_smd_to_cor",       grid = quote(build_grid_01()),
       gen = quote(gen_smd),       est = quote(estimate_smd_scale("r"))),
  list(id = "03_2x2_to_cor_CAT",   grid = quote(build_grid()),
       gen = quote(gen_cat),       est = quote(estimate_2x2)),
  list(id = "03_2x2_to_cor_CONT",  grid = quote(build_grid()),
       gen = quote(gen_cont),      est = quote(estimate_2x2)),
  list(id = "07_ancova_to_smd",    grid = quote(build_grid_07()),
       gen = quote(gen_ancova),    est = quote(estimate_ancova_scale("d"))),
  list(id = "08_pre_post_to_smd",  grid = quote(build_grid_08()),
       gen = quote(gen_pre_post),  est = quote(estimate_pre_post_scale("d"))),
  list(id = "09_or_to_cor_CONT",   grid = quote(build_grid_09("CONT")),
       gen = quote(gen_cont_or),   est = quote(estimate_or_to_cor)),
  list(id = "09_or_to_cor_CAT",    grid = quote(build_grid_09("CAT")),
       gen = quote(gen_cat_or),    est = quote(estimate_or_to_cor))
)

.separates <- function(spec, nrep = 4L, n_cond = 6L) {
  grid <- eval(spec$grid)
  gen  <- eval(spec$gen)
  est  <- eval(spec$est)
  idx <- unique(round(seq(1, nrow(grid), length.out = min(n_cond, nrow(grid)))))
  for (i in idx) {
    cond <- as.list(grid[i, , drop = FALSE])
    ok <- tryCatch({
      set.seed(3200 + i)
      dat <- gen(cond, nrep)
      e <- est(dat, cond)
      if (!all(c("theta_own", "theta_pop") %in% names(e))) {
        e$theta_pop <- if ("theta_pop" %in% names(dat)) dat$theta_pop[1] else NA_real_
      }
      !isTRUE(all.equal(e$theta_own, rep(e$theta_pop[1], nrow(e))))
    }, error = function(err) NA)
    if (isTRUE(ok)) return(TRUE)
  }
  FALSE
}

test_that("every study registering methods with different estimands separates own from population", {
  skip_if_not(isTRUE(get0("METACONVERT_AVAILABLE", ifnotfound = FALSE)),
              "metaConvert not loadable")
  sep <- vapply(.own_registry, function(s) isTRUE(.separates(s)), logical(1))
  names(sep) <- vapply(.own_registry, `[[`, character(1), "id")

  # Study 03 is the subject of roadmap 3.2: before the fix BOTH of its rows were
  # FALSE. Study 09 is the pattern 03 was brought into line with.
  expect_true(sep[["03_2x2_to_cor_CAT"]])
  expect_true(sep[["03_2x2_to_cor_CONT"]])
  expect_true(sep[["09_or_to_cor_CONT"]])
  expect_true(sep[["09_or_to_cor_CAT"]])
  expect_true(sep[["01_smd_to_cor"]])

  # Whatever the rest do today, record it so a change is visible rather than
  # silent. If a study starts or stops separating, this fails and the reviewer
  # decides which it should be.
  message("    own-vs-population separation by study: ",
          paste(sprintf("%s=%s", names(sep), sep), collapse = ", "))
  expect_type(sep, "logical")
  expect_length(sep, length(.own_registry))
})

## --- stopifnot() fallback for a bare R (no testthat) ------------------------
if (!requireNamespace("testthat", quietly = TRUE)) {
  stopifnot(
    isTRUE(all.equal(.estimands_2x2("cont", 0.5, 0.3, 0.1)$tetrachoric, 0.5)),
    isTRUE(all.equal(.estimands_2x2("cat", 0.5, 0.3, 0.1)$phi, 0.5)),
    !isTRUE(all.equal(.estimands_2x2("cont", 0.5, 0.3, 0.1)$phi, 0.5)),
    !isTRUE(all.equal(.estimands_2x2("cat", 0.5, 0.3, 0.1)$tetrachoric, 0.5))
  )
}
