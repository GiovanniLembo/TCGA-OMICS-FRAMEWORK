## load_framework.R -----------------------------------------------------------
## Sources every module. All CLI scripts start with:
##   source("R/load_framework.R")
## -----------------------------------------------------------------------------

.framework_root <- function() {
  ## works whether the script is run from the repo root or from scripts/
  for (cand in c(".", "..", "../..")) {
    if (file.exists(file.path(cand, "R", "utils.R"))) return(normalizePath(cand))
  }
  stop("cannot locate the repository root - run from the repo directory", call. = FALSE)
}

local({
  root <- .framework_root()
  files <- c("utils.R", "config.R", "gdc_query.R", "download.R",
             "cohort.R", "explore.R", "de_rnaseq.R", "mutations.R",
             "methylation.R", "cnv.R", "mirna.R", "report.R")
  for (f in files) source(file.path(root, "R", f))
  assign("FRAMEWORK_ROOT", root, envir = globalenv())
})

suppressPackageStartupMessages({
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
})

message("TCGA omics framework loaded (root: ", FRAMEWORK_ROOT, ")")
