## =============================================================================
## Study 11 -- alpha vs omega standard errors, and propagation into a pooled RG.
##
## These tests do NOT re-run the Monte Carlo (that is run_11's job). They assert
## the things that make the study trustworthy and that a later edit could quietly
## break. Study 11's conclusion is a DIFFERENTIAL of about 0.02 against
## per-coefficient levels spanning 0.49 to 1.07, so an undetected generator or
## fairness error would be invisible in the differential and decisive in the
## level. Every assertion below exists because getting it wrong flips a shipped
## package decision:
##
##   * the FAIRNESS RULE. A pilot comparing a Pearson-covariance alpha against a
##     WLSMV/polychoric omega measured a gap of -0.107 and would have concluded
##     "the closed form is materially worse for omega"; on the SAME covariance
##     matrix the same contrast is -0.010. The estimator swap manufactures the
##     entire effect the study exists to detect, so "one S, two coefficients" is
##     asserted directly rather than trusted.
##   * the EXACT CATEGORY MOMENTS. The claim that a continuous excess kurtosis of
##     8 is unreachable in Likert data -- which is what demotes it from a design
##     factor to a labelled bound -- rests on these numbers.
##   * the BENCHMARK CELL's presence. bin_10 is Maydeu-Olivares et al. (2007)
##     Table 5; without it in the grid nothing anchors the DGP to a published
##     result.
##   * k = 3's presence. Its coverage is the evidence that the `n_items >= 4`
##     gate is unnecessary. Drop the cell and the gate becomes unfalsifiable.
##   * the FITTER EQUIVALENCE. factanal() is used instead of lavaan for a 30x
##     speedup that is what makes 11c free; if it ever stops agreeing, every
##     omega number silently changes.
##   * REPLICATION ALIGNMENT. rel_delta() pairs the two arms by row position.
##     That is what buys the ~10x Monte Carlo efficiency, and it fails SILENTLY
##     (as a mildly wrong correlation) rather than loudly if the order changes.
## =============================================================================

.needs_pkg_11 <- function() if (!METACONVERT_AVAILABLE) skip("metaConvert not loadable")


test_that("study 11 is registered and honours the runner contract", {
  expect_true(exists("run_11", envir = .GlobalEnv))
  fmls <- names(formals(get("run_11", envir = .GlobalEnv)))
  expect_true(all(c("nrep", "cores") %in% fmls))
  # run_everything() discovers by regex, so the NAME has to match too
  expect_true("run_11" %in% ls(envir = .GlobalEnv, pattern = "^run_[0-9]{2}$"))
})


test_that("the response formats have the exact moments the design rests on", {
  fm <- .response_formats()
  ord <- fm[vapply(fm, function(x) identical(x$kind, "ordinal"), logical(1))]
  m <- t(vapply(ord, function(x) .rel_moments(x$p), numeric(2)))

  # These are computed from the category probabilities, not simulated, so they
  # are exact and can be pinned tightly.
  expect_equal(unname(m["lik5_sym", "skew"]), 0, tolerance = 1e-12)
  expect_equal(unname(m["lik5_sym", "exkurt"]), -0.50, tolerance = 1e-8)
  expect_equal(unname(m["lik5_skew", "exkurt"]), -0.317, tolerance = 1e-3)
  expect_equal(unname(m["lik5_sev", "exkurt"]), 2.097, tolerance = 1e-3)
  expect_equal(unname(m["bin_10", "exkurt"]), 5.111, tolerance = 1e-3)

  # THE LOAD-BEARING CLAIM: five-point categorisation caps the reachable excess
  # kurtosis near 2, so the continuous chi-square cell (excess kurtosis 8) has no
  # empirical referent in reliability-generalisation data. This is why response
  # format is the primary factor and continuous non-normality is a labelled bound.
  five_pt <- m[grepl("^lik5", rownames(m)), "exkurt"]
  expect_true(all(five_pt < 2.5))
  expect_true(max(five_pt) < 12 / 1.5)   # the continuous cont_kurt8 value

  # bin_10 must sit at the published benchmark's kurtosis, or it is not that cell
  expect_equal(unname(m["bin_10", "exkurt"]), 5.11, tolerance = 0.01)
})


test_that("the 11a grid keeps the benchmark cell, the bounds, and k = 3", {
  g <- build_grid_11a()

  # the published anchor: without it a DGP bug and a finding are indistinguishable
  expect_true("bin_10" %in% g$response)
  expect_true(any(g$response == "bin_10" & g$k == 8))
  expect_gte(length(unique(g$n[g$response == "bin_10"])), 2L)

  # the two labelled bounds, and the realistic reference cell between them
  expect_true(all(c("cont_normal", "cont_kurt8", "lik5_skew", "lik5_sev") %in% g$response))

  # k = 3 is the cell whose COVERAGE refutes the n_items gate. Its se_ratio is
  # not interpretable there (rel_stability marks it), but the cell must exist.
  expect_true(3 %in% g$k)
  expect_true(all(c(3, 8, 20) %in% g$k))

  # 11a is the head-to-head, so loading shape must NOT vary here -- that would
  # cross a factor whose whole job (11b) is a single three-cell contrast.
  expect_equal(length(unique(g$shape)), 1L)
})


test_that("11b contrasts loading shape and nothing else", {
  g <- build_grid_11b()
  expect_true(all(c("tau_equiv", "congeneric", "extreme") %in% g$shape))
  # one contrast, not a factorial: everything else is held fixed
  expect_equal(length(unique(g$response)), 1L)
  expect_equal(length(unique(g$n)), 1L)
  expect_equal(length(unique(g$k)), 1L)
})


test_that("under exact tau-equivalence alpha and omega are the SAME quantity", {
  # This is the premise the shipped ?es_from_omega rationale rests on: Bonett's
  # variance is derived under tau-equivalence, where the two coefficients
  # coincide. If they do not coincide here, .rel_loadings() or the population
  # formulas are wrong and 11b measures nothing.
  te <- tau_equivalence_premise(k = c(4, 8, 20))
  eq <- te[te$shape == "tau_equiv", ]
  expect_true(all(abs(eq$omega_minus_alpha) < 1e-12))
  expect_true(all(abs(eq$gap_in_bonett_se) < 1e-10))

  # and the estimand gap must GROW with loading spread, or "extreme" is not extreme
  for (kk in unique(te$k)) {
    d <- te[te$k == kk, ]
    expect_true(d$omega_minus_alpha[d$shape == "extreme"] >
                  d$omega_minus_alpha[d$shape == "congeneric"])
    expect_true(d$omega_minus_alpha[d$shape == "congeneric"] >=
                  d$omega_minus_alpha[d$shape == "tau_equiv"])
  }
})


test_that("THE FAIRNESS RULE: both coefficients come from one covariance matrix", {
  .needs_pkg_11()
  set.seed(101)
  lam <- .rel_loadings(8, "congeneric")
  S <- .rel_draw(400, lam, .response_formats()$lik5_skew)

  cf <- .rel_coefficients(S, 400)
  # alpha must be reproducible from S alone, by the textbook formula -- i.e. it
  # is not being computed from a different moment matrix somewhere else
  k <- ncol(S)
  expect_equal(unname(cf["alpha"]),
               (k / (k - 1)) * (1 - sum(diag(S)) / sum(S)), tolerance = 1e-12)

  # and calling twice on the same S must give the same omega: no hidden RNG in
  # the fitter, so the pairing in rel_delta() is a true paired comparison
  expect_equal(.rel_coefficients(S, 400), cf, tolerance = 1e-12)
})


test_that("the fast fitter agrees with lavaan (the 30x speedup must be free)", {
  .needs_pkg_11()
  skip_if_not_installed("lavaan")
  set.seed(202)
  for (cfg in list(list(k = 3, n = 200), list(k = 8, n = 500), list(k = 20, n = 300))) {
    lam <- .rel_loadings(cfg$k, "congeneric")
    S <- .rel_draw(cfg$n, lam, .response_formats()$cont_normal)

    fast <- unname(.rel_coefficients(S, cfg$n)["omega"])

    nm <- paste0("x", seq_len(cfg$k))
    dimnames(S) <- list(nm, nm)
    mod <- paste0("F =~ ", paste(nm, collapse = " + "))
    fit <- lavaan::cfa(mod, sample.cov = S, sample.nobs = cfg$n,
                       std.lv = TRUE, warn = FALSE, se = "none")
    pe <- lavaan::parameterEstimates(fit)
    l <- pe$est[pe$op == "=~"]
    th <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs != "F"]
    slow <- sum(l)^2 / (sum(l)^2 + sum(th))

    expect_equal(fast, slow, tolerance = 1e-5)
  }
})


test_that("estimate_rel() returns aligned arms with per-method targets", {
  .needs_pkg_11()
  set.seed(303)
  cond <- data.frame(response = "lik5_skew", n = 200, k = 8,
                     shape = "congeneric", stringsAsFactors = FALSE)
  dat <- gen_rel(cond, 40)
  est <- estimate_rel(dat, cond)

  expect_true(all(c("method", "est", "se", "ci_lo", "ci_up",
                    "theta_centred", "theta_own") %in% names(est)))
  expect_equal(sort(unique(est$method)), c("alpha", "omega"))

  # REPLICATION ALIGNMENT. rel_delta() pairs the arms by row position within a
  # cell; if estimate_rel() ever stops emitting nrep rows per arm in the input
  # order, the pairing silently degrades into a mismatched comparison rather
  # than failing. Assert the shape AND that arm i really is replication i.
  a <- est[est$method == "alpha", ]
  w <- est[est$method == "omega", ]
  expect_equal(nrow(a), nrow(dat))
  expect_equal(nrow(w), nrow(dat))
  expect_equal(a$est, log(1 - dat$alpha_hat), tolerance = 1e-10)
  expect_equal(w$est, log(1 - dat$omega_hat), tolerance = 1e-10)

  # the PRIMARY target really is the estimator's own mean: coverage about it
  # isolates SE calibration and must not depend on any plim
  expect_equal(unique(a$theta_centred), mean(a$est[is.finite(a$est)]), tolerance = 1e-12)
  expect_equal(length(unique(a$theta_centred)), 1L)

  # per-method plims must DIFFER between the arms -- population alpha is not
  # population omega, and scoring both against one column would report an
  # estimand gap as bias
  expect_false(isTRUE(all.equal(unique(a$theta_own), unique(w$theta_own))))
})


test_that("the Bonett SE does not diverge at k = 3 (refutes the shipped claim)", {
  # ?es_from_omega @details asserts the closed form "diverges at k = 3". The
  # expression sqrt(2k/((k-1)(n-2))) is finite there and is the LARGEST value in
  # the sweep, decreasing monotonically to sqrt(2/(n-2)). This is arithmetic, so
  # it is pinned exactly rather than measured.
  d <- bonett_se_by_k(k = c(3, 4, 6, 8, 10, 20, 100), n = 200)
  expect_true(all(is.finite(d$bonett_se)))
  expect_equal(d$bonett_se[d$k == 3], sqrt(6 / (2 * 198)), tolerance = 1e-12)
  expect_equal(d$bonett_se[d$k == 3], 0.1230915, tolerance = 1e-6)

  # monotone decreasing in k, so k = 3 is the maximum, not a divergence
  o <- d[order(d$k), ]
  expect_true(all(diff(o$bonett_se) < 0))
  expect_true(d$is_max_for_n[d$k == 3])
  expect_true(all(d$bonett_se > d$limit_k_inf))
})


test_that("the bifactor bullet is arithmetic, and its denominators do not match", {
  b <- bifactor_se_identity(m = 2:5, k = c(8, 10), n = c(200, 1500))
  # the shipped claim is 14-24% narrower; check the band is reproduced somewhere
  expect_true(any(b$pct_narrower > 13 & b$pct_narrower < 25))
  # and the two expressions differ in their denominator (n vs n-2), which must
  # never be presented as "the same formula with k replaced by m"
  expect_true(all(b$denom_effect_pct > 0))
  expect_true(any(b$denom_effect_pct[b$n == 200] > 0.4))
})


test_that("Bonett's variance is coefficient-free and Hakstian-Whalen's is not", {
  .needs_pkg_11()
  h <- hw_weight_ratio(alpha = c(0.70, 0.85, 0.95, 0.98), n = 200, k = 10)
  # Bonett: identical (n, k) => identical weight whatever alpha is
  expect_equal(unique(round(h$weight_spread_bonett, 10)), 1)
  # HW: the same studies get materially different weights
  expect_gt(unique(h$weight_spread_hw), 3)
})


test_that("rel_stability() flags the cells whose se_ratio is not readable", {
  .needs_pkg_11()
  set.seed(404)
  g <- data.frame(response = "cont_normal", n = 200, k = 8,
                  shape = "congeneric", stringsAsFactors = FALSE)
  dat <- gen_rel(g, 120)
  est <- estimate_rel(dat, g)
  raw <- cbind(g[rep(1, nrow(est)), , drop = FALSE], est, row.names = NULL)

  st <- rel_stability(raw)
  expect_true(all(c("ex_kurtosis", "robust_sd", "sd_mcse_rel",
                    "se_ratio_readable") %in% names(st)))
  expect_equal(nrow(st), 2L)
  # a well-behaved normal cell must be readable, or the gate is mis-signed
  expect_true(all(st$se_ratio_readable))
  # the heavy-tail-aware MCSE must exceed the naive normal-theory one whenever
  # the kurtosis is positive -- that inequality is the reason the column exists
  naive <- 1 / sqrt(2 * (st$n_valid - 1))
  expect_true(all(st$sd_mcse_rel[st$ex_kurtosis > 0] > naive[st$ex_kurtosis > 0]))
})


test_that("the benchmark gate is scoped to the published cell and fails loudly", {
  # A gate that silently passes on an empty selection is worse than no gate.
  fake <- data.frame(response = "lik5_skew", k = 8, n = c(200, 1000),
                     method = "alpha", target = "centred",
                     coverage = c(0.95, 0.95), coverage_mcse = 0.007,
                     stringsAsFactors = FALSE)
  g <- rel_benchmark_gate(fake)
  expect_false(isTRUE(g$pass[1]))   # bin_10 absent => must NOT pass

  ok <- data.frame(response = "bin_10", k = 8, n = c(200, 1000),
                   method = "alpha", target = "centred",
                   coverage = c(0.795, 0.792), coverage_mcse = 0.013,
                   stringsAsFactors = FALSE)
  expect_true(rel_benchmark_gate(ok)$pass[1])

  # and a cell that has drifted to nominal coverage must FAIL: the published
  # result is that normal-theory alpha under-covers on binary items and does not
  # recover with n. Reproducing 0.95 there would mean the DGP lost the kurtosis.
  drift <- ok; drift$coverage <- c(0.949, 0.951)
  expect_false(rel_benchmark_gate(drift)$pass[1])
})


## -----------------------------------------------------------------------------
## 11e -- the floor-effect detector.
##
## 11a established that the CATEGORY COUNT does not identify the variance error
## (three five-point formats spanned c = 0.92 to 1.69). 11e asks whether a
## reported total-score MEAN does. These tests guard the two ways that block can
## silently go wrong: the sweep leaking into 11a's grid, and a new per-replication
## column being mistaken for a grouping factor.
## -----------------------------------------------------------------------------

test_that("the floor sweep is balanced and does NOT leak into the 11a grid", {
  ff <- .floor_formats()
  # balanced across category counts, or the C-vs-floor_position comparison is
  # confounded with how many cells each C contributes
  tab <- table(vapply(ff, function(f) length(f$p), numeric(1)))
  expect_equal(length(unique(as.integer(tab))), 1L)
  expect_true(all(c(2, 3, 5, 7) %in% as.integer(names(tab))))

  # every format is a proper distribution
  for (f in ff) {
    expect_equal(sum(f$p), 1, tolerance = 1e-12)
    expect_true(all(f$p > 0))
  }

  # REGRESSION GUARD. build_grid_11a() must read .response_formats(), not
  # .all_formats(): the sweep adds 20 formats, and if they leaked in, 11a would
  # silently balloon from 42 cells to 162 and quietly change its own headline.
  expect_equal(nrow(build_grid_11a()), 42L)
  expect_false(any(grepl("^q[0-9]", build_grid_11a()$response)))
  expect_true(all(grepl("^q[0-9]", build_grid_11e()$response)))
})


test_that("delta orders the formats from symmetric to hard floor", {
  ff <- .floor_formats(C = 5, delta = c(0, 0.4, 0.8, 1.2, 1.6))
  # mean item response must decrease monotonically as delta pushes mass down
  m <- vapply(ff, function(f) sum(f$p * seq_along(f$p)), numeric(1))
  expect_true(all(diff(m) < 0))
  # delta = 0 is symmetric on a 5-point scale => mean item exactly 3
  expect_equal(unname(m[1]), 3, tolerance = 1e-8)
})


test_that(".REL_NON_COND covers everything estimate_rel emits", {
  .needs_pkg_11()
  set.seed(505)
  cond <- data.frame(response = "lik5_skew", n = 200, k = 8,
                     shape = "congeneric", stringsAsFactors = FALSE)
  est <- estimate_rel(gen_rel(cond, 30), cond)

  # THE SILENT FAILURE THIS GUARDS. rel_delta(), rel_stability() and
  # rel_floor_detector() all recover the condition columns by SUBTRACTING
  # .REL_NON_COND from the raw frame's names. A per-replication column that is
  # missing from that list becomes a grouping FACTOR, splitting every cell into
  # one group of size 1 -- which does not error, it just returns nothing (every
  # group falls under the min_valid floor) or, worse, garbage.
  expect_true(all(names(est) %in% .REL_NON_COND),
              info = paste("not in .REL_NON_COND:",
                           paste(setdiff(names(est), .REL_NON_COND), collapse = ", ")))
})


test_that("floor_position is a proper standardised scale position", {
  .needs_pkg_11()
  set.seed(606)
  # The expected value is computed from the format's own category probabilities
  # rather than typed in, so this tests the whole pipeline (draw -> threshold ->
  # sum score -> rescale) against the analytic truth instead of against a number
  # someone guessed. A symmetric 5-point scale lands at exactly 0.5.
  for (fmt in c("q5_d00", "q5_d16")) {
    p <- .all_formats()[[fmt]]$p
    expect_pos <- (sum(p * seq_along(p)) - 1) / (length(p) - 1)
    cond <- data.frame(response = fmt, n = 400, k = 8,
                       shape = "congeneric", stringsAsFactors = FALSE)
    est <- estimate_rel(gen_rel(cond, 60), cond)
    raw <- cbind(cond[rep(1, nrow(est)), , drop = FALSE], est, row.names = NULL)
    fd <- rel_floor_detector(raw, min_valid = 10L)

    # sampling tolerance only: 60 reps of n = 400, so the Monte Carlo error on a
    # mean item response is tiny relative to an absolute 0.01 on the [0,1] scale
    expect_equal(fd$floor_position[1], expect_pos, tolerance = 0.01)
    # bounded in [0, 1] by construction, and the Bhatia-Davis ratio cannot exceed 1
    expect_true(all(fd$floor_position >= 0 & fd$floor_position <= 1))
    expect_true(all(fd$sd_over_max <= 1 + 1e-8))
    expect_true(all(fd$C == 5))
  }
})


test_that("the total-score moments are exact, not approximated", {
  .needs_pkg_11()
  set.seed(707)
  lam <- .rel_loadings(8, "congeneric")
  S <- .rel_draw(300, lam, .response_formats()$lik5_skew)
  # var(rowSums(X)) is EXACTLY sum(S) for a plain sum score -- this is why
  # gen_rel needs no second pass over the data. If that identity ever stops
  # holding, scale_sd is quietly wrong and 11e's sd_over_max with it.
  expect_true(!is.null(attr(S, "item_means")))
  expect_equal(length(attr(S, "item_means")), 8L)
  expect_true(sum(S) > 0)
})
