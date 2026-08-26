## =============================================================================
## metaConvert -- toutes les verifications, depuis une console R.
##
##   setwd("c:/Users/coren/Documents/sideprojet/metaConvert")
##   source("checks.R")
##
##   check_all()                       # les trois suites            (~30 min)
##   check_archived()                  # juste tests_save/checked    (~25 min)
##   check_archived("test-flags.R")    # un seul fichier de tests_save
##   check_archived(verbose = TRUE)    # voir QUELS tests echouent
##   check_deps()                      # hygiene des dependances     (instantane)
##   check_cran()                      # R CMD check --as-cran       (~10 min)
##   check_nosuggests()                # la flavour qui a regresse deux fois
##   check_win("devel")                # win-builder  (envoi externe, reponse par mail)
##   check_revdeps()                   # depuis une session R PROPRE
##   check_before_cran()               # la sequence complete, du moins cher au plus cher
##
## LE PIEGE N°1 : sans NOT_CRAN="true", 12 fichiers / ~1585 assertions sont
## silencieusement sautes et la suite passe au vert en testant beaucoup moins.
## Chaque fonction ci-dessous le pose elle-meme.
## =============================================================================

PKG_ROOT <- "c:/Users/coren/Documents/sideprojet/metaConvert"

## References au 2026-08-26 (branche dev/2.1.0). Mets-les a jour quand elles bougent
## POUR UNE RAISON QUE TU PEUX EXPLIQUER -- un total qui BAISSE sans echec signifie que
## des assertions ont cesse de s'executer, pas que tout va bien.
BASELINE <- c(main = 3859L, archived = 7202L, sims = 792L)


.report <- function(label, r, key = NULL) {
  d <- as.data.frame(r)
  n <- c(pass = sum(d$passed), fail = sum(d$failed),
         err = sum(d$error), skip = sum(d$skipped))
  msg <- sprintf("%-9s PASS %d  FAIL %d  ERROR %d  SKIP %d",
                 label, n[["pass"]], n[["fail"]], n[["err"]], n[["skip"]])
  if (!is.null(key) && key %in% names(BASELINE))
    msg <- paste0(msg, sprintf("   (reference %d, ecart %+d)",
                               BASELINE[[key]], n[["pass"]] - BASELINE[[key]]))
  message(msg)
  ## REGARDE LA COLONNE error, PAS SEULEMENT failed : un bloc test_that qui ERREUR
  ## rapporte failed = 0 tout en abandonnant chaque assertion qui suit.
  bad <- d[d$failed > 0 | d$error > 0, c("file", "test"), drop = FALSE]
  if (nrow(bad)) { message("  -- blocs en echec --"); print(bad, row.names = FALSE) }
  invisible(n)
}


## ---------------------------------------------------------------- les trois suites

#' Suite principale : tests/testthat/  (~4 min)
check_main <- function() {
  withr::with_dir(PKG_ROOT, withr::with_envvar(c(NOT_CRAN = "true"),
    .report("MAIN", devtools::test(reporter = "silent"), "main")))
}


#' Suite archivee : tests_save/checked/  (~25 min)
#'
#' TROIS PARTICULARITES, et chacune la casse si on l'oublie :
#'   * devtools::test() NE LA LANCE PAS -- elle est hors de tests/testthat/, donc
#'     invisible pour lui. C'est une suite a part entiere, et c'est elle qui a
#'     rattrape une regression de 19 assertions que la principale laissait passer.
#'   * pkgload::load_all() est obligatoire : ces tests appellent des fonctions
#'     INTERNES (.d_j, les kernels pre/post...).
#'   * package = "metaConvert" sinon les helpers ne se resolvent pas.
#'
#' @param files un ou plusieurs noms de fichiers pour n'en lancer qu'une partie
#' @param verbose TRUE pour voir quels tests echouent plutot que le seul total
check_archived <- function(files = NULL, verbose = FALSE) {
  withr::with_dir(PKG_ROOT, withr::with_envvar(c(NOT_CRAN = "true"), {
    pkgload::load_all(PKG_ROOT, quiet = TRUE)
    rep <- if (verbose) "summary" else "silent"
    if (is.null(files)) {
      .report("ARCHIVED", testthat::test_dir("tests_save/checked", reporter = rep,
                                             package = "metaConvert"), "archived")
    } else {
      res <- lapply(files, function(f) {
        p <- if (file.exists(f)) f else file.path("tests_save/checked", f)
        as.data.frame(testthat::test_file(p, reporter = rep, package = "metaConvert"))
      })
      .report("ARCHIVED", do.call(rbind, res))
    }
  }))
}


#' Suite simulations : une TROISIEME suite que ni l'une ni l'autre n'execute (~1 min)
check_sims <- function() {
  ## run_tests.R deduit sa racine du repertoire courant, d'ou le with_dir
  withr::with_dir(file.path(PKG_ROOT, "simulations"),
    system2(file.path(R.home("bin"), "Rscript.exe"), "tests/run_tests.R"))
}


#' Les trois suites, du moins cher au plus cher (pour echouer tot)
check_all <- function() { check_sims(); check_main(); check_archived(); invisible(NULL) }


## ---------------------------------------------------------------- R CMD check

#' R CMD check --as-cran
check_cran <- function(build_vignettes = TRUE) {
  tar <- pkgbuild::build(PKG_ROOT, dest_path = tempdir(), vignettes = build_vignettes)
  rcmdcheck::rcmdcheck(tar, args = "--as-cran", error_on = "never")
}


#' La flavour noSuggests, reproduite localement et SANS RIEN DESINSTALLER
#'
#' _R_CHECK_DEPENDS_ONLY_ restreint le check aux Depends/Imports/LinkingTo. Ce n'est
#' pas optionnel ici : la route tetrachorique a besoin de mvtnorm, un *Suggests*, et
#' sans lui elle renvoie NA avec une simple notice -- donc un test qui affirme une
#' valeur tetrachorique echoue EN SILENCE, et R CMD check ERREUR sur la flavour
#' noSuggests de CRAN. CLAUDE.md note que ca a regresse DEUX fois.
check_nosuggests <- function() {
  tar <- pkgbuild::build(PKG_ROOT, dest_path = tempdir(), vignettes = FALSE)
  rcmdcheck::rcmdcheck(tar, args = "--no-manual",
                       libpath = .libPaths(),
                       env = c(`_R_CHECK_DEPENDS_ONLY_` = "true"),
                       error_on = "never")
}


## ---------------------------------------------------------------- plateformes CRAN

#' win-builder : R-devel / release / oldrelease sur Windows
#'
#' CE SONT DES ENVOIS VERS UN SERVICE EXTERNE. La tarball part sur
#' win-builder.r-project.org et le resultat arrive PAR MAIL a l'adresse Maintainer de
#' DESCRIPTION, sous ~30 min. Rien ne s'affiche dans la console : si le mail n'arrive
#' pas, regarde les indesirables.
#'
#' R-devel est celui qui compte pour CRAN -- c'est la version contre laquelle ils
#' testeront. release/oldrelease verifient que tu ne casses pas les utilisateurs qui
#' n'ont pas mis a jour.
#'
#' @param which "devel" (defaut), "release", "oldrelease" ou "all"
check_win <- function(which = c("devel", "release", "oldrelease", "all")) {
  which <- match.arg(which)
  if (which %in% c("devel", "all"))      devtools::check_win_devel(PKG_ROOT)
  if (which %in% c("release", "all"))    devtools::check_win_release(PKG_ROOT)
  if (which %in% c("oldrelease", "all")) devtools::check_win_oldrelease(PKG_ROOT)
  message("Envoye. Resultat par mail sous ~30 min a l'adresse Maintainer de DESCRIPTION.")
}


#' macOS builder -- l'autre plateforme que CRAN teste (envoi externe lui aussi)
check_mac <- function() devtools::check_mac_release(PKG_ROOT)


## ---------------------------------------------------------------- dependances

#' Hygiene des dependances : declare vs reellement utilise
#'
#' R CMD check attrape "utilise mais non declare". Il n'attrape PAS l'inverse : un
#' Imports declare que plus rien n'appelle reste invisible, et c'est une installation
#' imposee aux utilisateurs pour rien. Il ne dit pas non plus si un Suggests manque
#' SUR CETTE MACHINE -- or un Suggests absent fait sauter des tests en silence.
check_deps <- function() {
  d <- read.dcf(file.path(PKG_ROOT, "DESCRIPTION"))
  field <- function(f) {
    if (!f %in% colnames(d)) return(character(0))
    x <- strsplit(gsub("[[:space:]]+", "", d[1, f]), ",")[[1]]
    setdiff(sub("[(].*", "", x), c("R", ""))
  }
  imports  <- field("Imports")
  suggests <- field("Suggests")

  ## ce que le code appelle explicitement, sous la forme pkg::fun
  src <- unlist(lapply(list.files(file.path(PKG_ROOT, "R"), pattern = "[.]R$",
                                  full.names = TRUE), readLines, warn = FALSE))
  src <- src[!grepl("^[[:space:]]*#", src)]
  used <- unique(sub("::.*", "", unlist(regmatches(
    src, gregexpr("[A-Za-z][A-Za-z0-9._]*::", src)))))
  ## ... plus ce que NAMESPACE importe en bloc
  ns <- readLines(file.path(PKG_ROOT, "NAMESPACE"), warn = FALSE)
  imp <- grep("^import", ns, value = TRUE)
  used <- unique(c(used, sub("^import[a-zA-Z]*\\(([^,)]+).*", "\\1", imp)))

  unused <- setdiff(imports, used)
  miss   <- suggests[!vapply(suggests, requireNamespace, logical(1), quietly = TRUE)]

  cat("Imports declares  :", paste(imports, collapse = ", "), "\n")
  cat("  jamais vus en pkg::fun ni en import() :",
      if (length(unused)) paste(unused, collapse = ", ") else "aucun", "\n")
  cat("  (verifie a la main avant d'en retirer un : il peut etre utilise via une\n")
  cat("   methode S3 ou un appel dynamique que ce grep ne voit pas)\n\n")
  cat("Suggests declares :", length(suggests), "\n")
  cat("  NON installes sur cette machine :",
      if (length(miss)) paste(miss, collapse = ", ") else "aucun", "\n")
  if (length(miss))
    cat("  -> des tests vont sauter en silence. Installe-les avant de conclure.\n")
  invisible(list(unused_imports = unused, missing_suggests = miss))
}


#' Dependances inverses (metaumbrella appelle les es_from_* en direct)
#'
#' A LANCER DEPUIS UNE SESSION R QUI N'A PAS CHARGE metaConvert : pkgload::load_all()
#' compte comme un chargement et fait echouer l'installation. Redemarre R d'abord.
check_revdeps <- function(tarball = NULL) {
  if ("metaConvert" %in% loadedNamespaces())
    stop("metaConvert est charge dans cette session -- redemarre R avant de lancer ceci")
  withr::with_dir(PKG_ROOT, {
    source("revdep/check-revdeps.R")
    get("check_revdeps", envir = globalenv())(tarball)
  })
}


## ---------------------------------------------------------------- la sequence

#' Tout ce qui se fait en local avant une soumission CRAN, du moins cher au plus cher
check_before_cran <- function() {
  check_deps()
  check_sims(); check_main(); check_archived()
  print(check_cran())
  print(check_nosuggests())
  message("\nRestent, et chacun depuis un etat propre :")
  message("  check_win('devel')   # envoi externe, reponse par mail sous ~30 min")
  message("  check_revdeps()      # REDEMARRE R d'abord")
}
