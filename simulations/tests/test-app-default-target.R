## =============================================================================
## Roadmap item 4.3 -- the app opened study 01b on the one benchmark that neither
## of its routes estimates.
##
## THE DEFECT. The app picks the "Scored against" default with a generic rule:
## population-scale targets first, then the other shared targets, then `own`. Study
## 01b records no target named "population" (its targets are `fisherz_biserial`,
## `fisherz_pointbiserial`, `own`), so the rule fell through to the alphabetically
## first shared target -- `fisherz_biserial` -- while the sidebar told the reader to
## prefer shared benchmarks over `own`.
##
## `fisherz_biserial` is atanh() of the biserial correlation. On the z scale
## viechtbauer returns a VARIANCE-STABILISING transform, not atanh() of anything, and
## lipsey_cooper returns atanh() of the POINT-biserial. At rho = 0.75, p_exp = 0.5 the
## three quantities are 0.973, 0.724 and 0.691 -- so the opening view scored both
## routes against a number neither is trying to produce, and both looked broken
## (coverage 0.736 / 0.773 against 0.944 / 0.934 on `fisherz_pointbiserial`).
##
## WHY THE FIX IS A DECLARED MAP AND NOT A RULE. The same shape -- a shared target
## every method misses by a margin that is large and flat in n -- is the intended
## FINDING in studies 07, 08 and 09b. No statistic in the aggregates separates "wrong
## scale" from "real estimand gap"; both are flat in n. So the judgement is declared
## per study in STUDY_RANK_ON, and these tests pin it against the shipped files so a
## declaration naming a target that does not exist cannot sit there doing nothing.
##
## STUDIES 03 AND 09 WERE MEASURED AND DELIBERATELY LEFT ALONE. Their `population`
## target is attained exactly by whichever route is correct for the mechanism, and
## where nothing attains it (09b: phi is the estimand and metaConvert ships no phi
## route) that gap is the study's result, not a defect. The tests below assert that,
## so a later "fix for consistency" cannot quietly delete a finding.
## =============================================================================

.app_file <- file.path(.sim_root, "app", "app.R")

## The app is a Shiny script, not a package: source it for the pure helpers only.
## Everything it does at load time is constant definitions plus a scan of
## data/aggregated, so this is safe -- but shiny::shinyApp() at the end must not run,
## hence the parse-and-evaluate-top-level-definitions approach rather than source().
.app_env <- local({
  e <- new.env(parent = globalenv())
  exprs <- parse(.app_file)
  for (ex in exprs) {
    ## take assignments only; skip library() calls and the shinyApp() launch
    if (is.call(ex) && length(ex) >= 3 &&
        as.character(ex[[1]])[1] %in% c("<-", "=") && is.name(ex[[2]])) {
      tryCatch(eval(ex, envir = e), error = function(err) NULL)
    }
  }
  e
})

.get <- function(nm) get(nm, envir = .app_env)

.agg <- function(stem) utils::read.csv(
  file.path(.sim_root, "data", "aggregated", paste0(stem, "_nrep1000.csv")),
  stringsAsFactors = FALSE)

test_that("the app's helpers were picked up from app.R", {
  for (nm in c("STUDY_RANK_ON", "STUDY_SHARED_SCALE", "target_order",
               "default_target", "unattained_targets", ".scale_mismatch"))
    expect_true(exists(nm, envir = .app_env), info = nm)
})

test_that("01b no longer opens on a target neither route estimates", {
  tg <- unique(.agg("01b_smd_to_cor_z")$target)
  d <- .get("default_target")
  expect_equal(d(tg, "01b_smd_to_cor_z"), "fisherz_pointbiserial")
  # and the old behaviour is what it replaced, so the regression is visible
  expect_equal(d(tg, NULL), "fisherz_biserial")
})

test_that("fisherz_biserial is genuinely nobody's estimand, and gets labelled", {
  a <- .agg("01b_smd_to_cor_z")
  u <- .get("unattained_targets")(a)
  expect_true("fisherz_biserial" %in% u)
  expect_false("fisherz_pointbiserial" %in% u)

  # the numbers behind the item, so a regeneration that changed them would show up
  cov_at <- function(t) mean(a$coverage[a$target == t], na.rm = TRUE)
  expect_lt(cov_at("fisherz_biserial"), 0.80)
  expect_gt(cov_at("fisherz_pointbiserial"), 0.93)
  expect_gt(cov_at("own"), 0.93)
})

test_that("every STUDY_RANK_ON entry names a target the study actually records", {
  # This is the guard the phantom column list of roadmap 1.2 did not have: a
  # declaration that names something absent would silently do nothing.
  m <- .get("STUDY_RANK_ON")
  for (stem in names(m)) {
    tg <- unique(.agg(stem)$target)
    expect_true(m[[stem]] %in% tg,
                info = paste0(stem, ": '", m[[stem]], "' not among ",
                              paste(tg, collapse = ", ")))
  }
})

test_that("no study opens on an unattained target while an attained one exists", {
  # The property the item is really about, asserted across every shipped aggregate
  # rather than for 01b alone: opening on a target nothing estimates is only tolerated
  # when there is no fixed-quantity alternative. That exception is real -- in 08 and
  # 09b `population` is what the review wants and no route delivers it, which is those
  # studies' finding -- and their only other shared target is the same-sample one,
  # against which coverage cannot be read at all.
  d <- .get("default_target"); u <- .get("unattained_targets")
  stems <- sub("_nrep1000[.]csv$", "",
               list.files(file.path(.sim_root, "data", "aggregated"),
                          pattern = "_nrep1000[.]csv$"))
  checked <- 0L
  for (stem in stems) {
    a <- .agg(stem)
    if (!"own" %in% a$target) next          # nothing to be wrong about
    tg <- unique(a$target)
    sel <- d(tg, stem); un <- u(a)
    fixed_shared <- setdiff(tg, c(grep("^own", tg, value = TRUE),
                                  grep("sample", tg, value = TRUE)))
    checked <- checked + 1L
    if (sel %in% un)
      expect_length(setdiff(fixed_shared, un), 0)
    else
      expect_false(sel %in% un, info = stem)
  }
  expect_gt(checked, 5L)
})

test_that("a same-sample target is never labelled unattained", {
  # It cannot be attained by construction -- `own` is a population quantity in every
  # study -- so the test would have fired on all 12 studies and said nothing. A method
  # DOES estimate that statistic; it is random rather than fixed, which the app warns
  # about separately. This is a regression guard: before it, 01a labelled both of its
  # sample targets and the label became noise.
  u <- .get("unattained_targets")
  stems <- sub("_nrep1000[.]csv$", "",
               list.files(file.path(.sim_root, "data", "aggregated"),
                          pattern = "_nrep1000[.]csv$"))
  for (stem in stems) {
    a <- .agg(stem)
    if (!"own" %in% a$target) next
    expect_length(grep("sample", u(a)), 0)
  }
  # 01a is the case that exposed it: nothing there should be labelled at all
  expect_length(u(.agg("01a_smd_to_cor_r")), 0)
})

test_that("08 and 09b keep population as the default, labelled for what it is", {
  # Their gap to it is the finding (the standardiser choice; no OR-to-correlation
  # route recovers phi), so the label informs rather than demotes.
  d <- .get("default_target"); u <- .get("unattained_targets")
  for (stem in c("08a_pre_post_to_smd_d", "08b_pre_post_to_smd_g", "09b_or_to_cor_CAT")) {
    a <- .agg(stem)
    expect_equal(d(unique(a$target), stem), "population", info = stem)
    expect_true("population" %in% u(a), info = stem)
  }
})

test_that("studies 03 and 09 keep their defaults -- the gap there is the finding", {
  d <- .get("default_target")
  for (stem in c("03a_2x2_to_cor_CAT", "03b_2x2_to_cor_CONT",
                 "09a_or_to_cor_CONT", "09b_or_to_cor_CAT")) {
    a <- .agg(stem)
    expect_equal(d(unique(a$target), stem), "population", info = stem)
  }
  # 03a: phi (r) IS the population target under a categorical mechanism, and the
  # tetrachoric's distance from it is the study's headline. 09b: nothing attains phi,
  # because the package ships no phi route -- which is what that study reports.
  a03 <- .agg("03a_2x2_to_cor_CAT")
  expect_false("population" %in% .get("unattained_targets")(a03))
  a09 <- .agg("09b_or_to_cor_CAT")
  expect_true("population" %in% .get("unattained_targets")(a09))
})

test_that("the scale-mismatch warning fires only where the scales really differ", {
  f <- .get(".scale_mismatch")
  # 09 keeps its shared targets on the r scale; the app opens on (r)
  expect_false(f("09b_or_to_cor_CAT", "r", "population"))
  expect_true(f("09b_or_to_cor_CAT", "z", "population"))
  expect_true(f("03a_2x2_to_cor_CAT", "z", "sample"))
  # `own` is scale-matched by construction, so it never fires
  expect_false(f("09b_or_to_cor_CAT", "z", "own"))
  # a study with no scale switch is never affected
  expect_false(f("01b_smd_to_cor_z", "z", "fisherz_biserial"))
  expect_false(f("07a_ancova_to_smd_d", "z", "population"))
  # missing inputs must not error
  expect_false(f(NULL, "z", "population"))
  expect_false(f("09b_or_to_cor_CAT", NULL, "population"))
})

test_that("STUDY_SHARED_SCALE names only studies that really carry a scale switch", {
  # The switch exists when the method names carry an (r)/(z) suffix; if a study were
  # listed here without one, the warning could never fire and the entry would rot.
  for (stem in names(.get("STUDY_SHARED_SCALE"))) {
    m <- unique(.agg(stem)$method)
    sfx <- sub(".*\\s\\(([^()]+)\\)$", "\\1", m[grepl("\\s\\([^()]+\\)$", m)])
    expect_gt(length(unique(sfx)), 1L, label = paste0(stem, " scale suffixes"))
    expect_true(unname(.get("STUDY_SHARED_SCALE")[[stem]]) %in% sfx, info = stem)
  }
})

test_that("target_order still puts shared benchmarks first and own last", {
  o <- .get("target_order")
  expect_equal(o(c("own", "sample", "population")), c("population", "sample", "own"))
  expect_equal(o(c("own_true_r", "own", "sample", "population")),
               c("population", "sample", "own", "own_true_r"))
  expect_equal(o(c("own", "fisherz_pointbiserial", "fisherz_biserial")),
               c("fisherz_biserial", "fisherz_pointbiserial", "own"))
})

## --- stopifnot() fallback for a bare R (no testthat) ------------------------
if (!requireNamespace("testthat", quietly = TRUE)) {
  .tg <- unique(.agg("01b_smd_to_cor_z")$target)
  stopifnot(
    .get("default_target")(.tg, "01b_smd_to_cor_z") == "fisherz_pointbiserial",
    "fisherz_biserial" %in% .get("unattained_targets")(.agg("01b_smd_to_cor_z")),
    !.get(".scale_mismatch")("09b_or_to_cor_CAT", "r", "population"),
    .get(".scale_mismatch")("09b_or_to_cor_CAT", "z", "population")
  )
}
