.onAttach <- function(libname, pkgname) {
  installed_version <- utils::packageVersion("metaConvert")

  tryCatch({
    cran_version <- NULL
    available_pkgs <- utils::available.packages(repos = "https://cloud.r-project.org")
    if ("metaConvert" %in% rownames(available_pkgs)) {
      cran_version <- package_version(available_pkgs["metaConvert", "Version"])
    }
    if (!is.null(cran_version) && installed_version < cran_version) {
      packageStartupMessage(
        sprintf(
          "\nmetaConvert Update Available!\nCurrently using version %s, but version %s is available on CRAN.\n",
          installed_version,
          cran_version
        ),
        "To update, run: install.packages('metaConvert')\n",
        "See NEWS for changes: https://cran.r-project.org/web/packages/metaConvert/news/news.html\n"
      )
    } else if (!is.null(cran_version) && installed_version == cran_version) {
      packageStartupMessage(
        paste0(
          "metaConvert v", installed_version, ":\n",
          "This version is still considered as a beta.\n",
          "\n-> Email us at cgosling@parisnanterre.fr if you notice any bug. <-\n",
          "If you use the package in the publication, you can cite the asssociated reference",
          "Gosling CJ, Cortese S, Solmi M, et al. metaConvert: an automatic suite for estimation of 11 different effect size measures and flexible conversion across them. Research Synthesis Methods. Published online 2025:1-12. doi:10.1017/rsm.2025.11."
        )
      )
    }
  }, error = function(e) {
    return(NULL)
  })

  invisible()
}


# test_version_check <- function() {
#   installed_version <- utils::packageVersion("metaConvert")
#
#   cran_version <- package_version("999.0.0")
#
#   if (!is.null(cran_version) && installed_version < cran_version) {
#     packageStartupMessage(
#       sprintf(
#         "\nmetaConvert Update Available!\nCurrently using version %s, but version %s is available on CRAN.\n",
#         installed_version,
#         cran_version
#       ),
#       "To update, run: install.packages('metaConvert')\n",
#       "See NEWS for changes: https://cran.r-project.org/web/packages/metaConvert/news/news.html\n"
#     )
#   } else if (!is.null(cran_version) && installed_version == cran_version) {
#     packageStartupMessage(
#       paste0(
#         "metaConvert v", installed_version, ":\n",
#         "\nSince the companion paper of the package is under review, ",
#         "this version is still considered as a beta.\n",
#         "\n-> Email us at cgosling@parisnanterre.fr if you notice any bug. <-\n",
#         "\nTo help you start with this package, you can retrieve a tutorial online ",
#         "(https://metaconvert.org/tutorial.html)",
#         "\n\nYou can also quickly generate data extraction sheets using ",
#         "(https://metaconvert.org/input.html)\n"
#       )
#     )
#   }
# }
#
# test_version_check()

