#!/usr/bin/env Rscript
## 02_run_analysis.R -----------------------------------------------------------
## The full pipeline: cohort -> exploration -> DESeq2 -> mutations -> methylation
## -> integration -> HTML report. Driven entirely by the config.
##
##   Rscript scripts/02_run_analysis.R --config config/brca_tumor_vs_normal.yml
##   Rscript scripts/02_run_analysis.R --config config/luad_stage.yml --only rnaseq
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
source("R/load_framework.R")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character", help = "path to the YAML config"),
  make_option("--only", type = "character", default = NULL,
              help = "restrict to some steps: explore,rnaseq,mutation,methylation"),
  make_option("--force-download", action = "store_true", default = FALSE,
              help = "bypass the cache")
)))
if (is.null(opt$config)) stop("--config is required", call. = FALSE)

cfg   <- load_config(opt$config)

steps <- if (is.null(opt$only)) {
  c("explore", "rnaseq", "mutation", "methylation", "cnv", "mirna")
} else {
  trimws(strsplit(opt$only, ",")[[1]])
}

log_step("run: ", cfg$contrast_name)
data <- fetch_all(cfg, force = isTRUE(opt$`force-download`))

de_table <- NULL; dmp_table <- NULL; cohort_summary <- NULL; vsd <- NULL

## --- RNA-seq cohort is the reference cohort for everything else --------------
rna_cohort <- NULL
if (!is.null(data$rnaseq)) {
  rna_cohort <- build_cohort(data$rnaseq, cfg, data$clinical, data$subtypes)
  cohort_summary <- rna_cohort$summary
  write_table(rna_cohort$summary, file.path(cfg$tables_dir, "cohort_summary.tsv"))
  print(rna_cohort$summary, row.names = FALSE)
}

if ("explore" %in% steps) {
  avail <- try(scan_availability(cfg$project,
                                 names(cfg$omics)[vapply(cfg$omics, isTRUE, logical(1))]),
               silent = TRUE)
  if (inherits(avail, "try-error")) avail <- NULL
  explore_project(cfg, avail = avail,
                  se = if (!is.null(rna_cohort)) rna_cohort$se else NULL)
}

if ("rnaseq" %in% steps && !is.null(rna_cohort)) {
  de <- run_deseq2(rna_cohort, cfg)
  de_table <- de$table
  vsd <- de$vsd          # reused by the copy-number dosage analysis
}

if ("mutation" %in% steps && !is.null(data$mutation)) {
  maf <- data$mutation
  if (!is.null(rna_cohort)) maf <- cohort_maf(maf, rna_cohort$se)
  analyse_mutations(maf, cfg, compare = isTRUE(cfg$mutation$compare_groups))
}

if ("methylation" %in% steps && !is.null(data$methylation)) {
  meth_cohort <- try(build_cohort(data$methylation, cfg, data$clinical, data$subtypes), silent = TRUE)
  if (inherits(meth_cohort, "try-error")) {
    log_warn("methylation cohort could not be built: ", as.character(meth_cohort))
  } else {
    dmp_table <- run_dmp(meth_cohort, cfg)
    integrate_meth_expression(dmp_table, de_table, cfg)
  }
}

if ("cnv" %in% steps && !is.null(data$cnv)) {
  cnv_cohort <- try(build_cohort(data$cnv, cfg, data$clinical, data$subtypes), silent = TRUE)
  if (inherits(cnv_cohort, "try-error")) {
    log_warn("CNV cohort could not be built: ", as.character(cnv_cohort))
  } else {
    analyse_cnv(cnv_cohort, cfg, vsd = vsd)
  }
}

if ("mirna" %in% steps && !is.null(data$mirna)) {
  mirna_cohort <- try(build_cohort(data$mirna, cfg, data$clinical, data$subtypes), silent = TRUE)
  if (inherits(mirna_cohort, "try-error")) {
    log_warn("miRNA cohort could not be built: ", as.character(mirna_cohort))
  } else {
    run_mirna_de(mirna_cohort, cfg)
  }
}

snapshot_run(cfg)
build_report(cfg, cohort_summary = cohort_summary, de_table = de_table, dmp_table = dmp_table)
log_ok("finished: ", cfg$results_dir)
