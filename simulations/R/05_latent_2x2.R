## =============================================================================
## Population quantities for a 2x2 table and its latent-normal model.
##
## WHY THESE LIVE IN R/ AND NOT IN A STUDY FILE. Studies 03 and 09 both need to
## state, for a given (rho, p_exp, p_case), what each method is actually
## ESTIMATING -- phi and the tetrachoric are different population quantities, and
## which one a method targets depends on the data-generating mechanism, not on the
## method alone. Until 2026-08 study 09 defined these helpers privately and study
## 03 carried its own copies of two of them under different names
## (phi_to_p11 / phi_attainable vs .phi_to_p11 / .phi_attainable), verified
## identical on 500 random inputs. Two definitions of one function is the defect
## roadmap item 1.6 removed on the package side; this is the simulation-side
## instance. One definition, one place, both studies.
##
## They are pure functions of the CONDITION, not of any replication, so they are
## cheap: each is evaluated once per condition, never per replication.
## =============================================================================

## Joint upper-quadrant probability of a standard bivariate normal with
## correlation rho, above the (1 - p) and (1 - q) quantiles. This is the cell
## probability p11 of the dichotomised table.
.biv_p11 <- function(rho, p, q) {
  cut_x <- stats::qnorm(1 - p)
  cut_y <- stats::qnorm(1 - q)
  if (abs(rho) < 1e-12) return(p * q)
  mvtnorm::pmvnorm(
    lower = c(cut_x, cut_y), upper = c(Inf, Inf),
    mean  = c(0, 0), corr = rbind(c(1, rho), c(rho, 1))
  )[[1]]
}

## phi implied by a cell probability p11 with margins p, q
.p11_to_phi <- function(p11, p, q) (p11 - p * q) / sqrt(p * q * (1 - p) * (1 - q))

## inverse: the p11 that produces a target phi
.phi_to_p11 <- function(phi, p, q) phi * sqrt(p * q * (1 - p) * (1 - q)) + p * q

.phi_attainable <- function(phi, p, q) {
  p11 <- .phi_to_p11(phi, p, q)
  p11 >= 0 && p11 <= min(p, q) && (p11 - p - q + 1) >= 0 &&
    (p - p11) >= 0 && (q - p11) >= 0
}

## The tetrachoric correlation of a table whose cell probability is p11: solve
## .biv_p11(rho) = p11 for rho. Used for the CAT mechanism, where the data are
## genuinely categorical but the tetrachoric methods still have a well-defined
## target (the latent correlation a normal-latent model would infer).
.p11_to_tetrachoric <- function(p11, p, q) {
  f <- function(rho) .biv_p11(rho, p, q) - p11
  lo <- f(-0.999); hi <- f(0.999)
  if (!is.finite(lo) || !is.finite(hi) || lo * hi > 0) return(NA_real_)
  stats::uniroot(f, c(-0.999, 0.999), tol = 1e-9)$root
}

## ---- the two estimands of a 2x2 table, per mechanism -------------------------

#' Population estimand of the phi and tetrachoric families, given the mechanism
#'
#' The point of the `own` target: score each method against the quantity IT is
#' estimating, so that what is left over is computational error rather than the
#' two methods measuring different things.
#'
#' @param dgm  "cat"  genuinely 4-cell multinomial with population phi = rho, or
#'             "cont" bivariate normal with correlation rho, then dichotomised
#' @param rho  the correlation the data were generated from (phi under "cat",
#'             the latent correlation under "cont")
#' @param p,q  the exposure and case marginal proportions
#' @return list(phi =, tetrachoric =), both on the r scale
#'
#' Under "cont" the tetrachoric estimand IS rho and phi's is the phi of the
#' dichotomised table; under "cat" phi's estimand IS rho and the tetrachoric's is
#' the latent correlation the table implies. In each mechanism exactly one of the
#' two equals rho, which is why a single shared target column cannot serve both.
.estimands_2x2 <- function(dgm, rho, p, q) {
  dgm <- match.arg(dgm, c("cat", "cont"))
  if (dgm == "cont") {
    list(phi = .p11_to_phi(.biv_p11(rho, p, q), p, q),
         tetrachoric = rho)
  } else {
    list(phi = rho,
         tetrachoric = .p11_to_tetrachoric(.phi_to_p11(rho, p, q), p, q))
  }
}
