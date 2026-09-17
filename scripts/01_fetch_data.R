#!/usr/bin/env Rscript
## 01_fetch_data.R -------------------------------------------------------------
## Download + prepare + cache every layer requested by a config.
## Safe to re-run: existing cache entries are reused unless --force is given.
##
##   Rscript scripts/01_fetch_data.R --config config/brca_tumor_vs_normal.yml
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
source("R/load_framework.R")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character", help = "path to the YAML config"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "ignore the cache and re-download"),
  make_option("--chunk", type = "integer", default = 20,
              help = "files per GDC API chunk; lower it on unstable networks [default %default]")
)))
if (is.null(opt$config)) stop("--config is required", call. = FALSE)

cfg  <- load_config(opt$config)
data <- fetch_all(cfg, force = opt$force)

log_step("cache summary for ", cfg$project)
for (nm in names(data)) {
  obj <- data[[nm]]
  if (is.null(obj)) { log_warn(nm, ": not available"); next }
  dims <- if (inherits(obj, "SummarizedExperiment")) sprintf("%d features x %d samples", nrow(obj), ncol(obj))
          else sprintf("%d rows", nrow(obj))
  log_ok(sprintf("%-12s %s", nm, dims))
}
log_ok("cache directory: ", cfg$cache_dir)
