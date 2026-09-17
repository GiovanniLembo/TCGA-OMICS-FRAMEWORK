#!/usr/bin/env Rscript
## install_dependencies.R ------------------------------------------------------
## Installs everything the framework needs. Run once:
##   Rscript install/install_dependencies.R
##   Rscript install/install_dependencies.R --optional   # + enrichment extras
## -----------------------------------------------------------------------------

optional <- "--optional" %in% commandArgs(trailingOnly = TRUE)

cran <- c("yaml", "optparse", "ggplot2", "ggrepel", "scales", "matrixStats", "R.utils")

bioc <- c("TCGAbiolinks", "SummarizedExperiment", "S4Vectors", "GenomicRanges",
          "DESeq2", "apeglm", "ashr", "limma", "maftools")

bioc_optional <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "sesame", "sva")

install_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    message("already installed: ", paste(pkgs, collapse = ", "))
    return(invisible(TRUE))
  }
  message("installing: ", paste(missing, collapse = ", "))
  installer(missing)
}

repos <- getOption("repos")
if (is.null(repos[["CRAN"]]) || repos[["CRAN"]] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

install_missing(cran, function(p) install.packages(p))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install_missing(bioc, function(p) BiocManager::install(p, ask = FALSE, update = FALSE))

if (optional) {
  install_missing(bioc_optional, function(p) BiocManager::install(p, ask = FALSE, update = FALSE))
}

## report --------------------------------------------------------------------
check <- c(cran, bioc, if (optional) bioc_optional)
status <- vapply(check, requireNamespace, logical(1), quietly = TRUE)
message("\n--- dependency status ---")
for (i in seq_along(check)) {
  message(sprintf("%-24s %s", check[i], if (status[i]) "ok" else "MISSING"))
}
if (any(!status)) {
  message("\nSome packages are missing. On Linux the usual cause is missing system\n",
          "libraries; install them and re-run:\n",
          "  sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev\n")
  quit(status = 1)
}
message("\nAll dependencies available.")
