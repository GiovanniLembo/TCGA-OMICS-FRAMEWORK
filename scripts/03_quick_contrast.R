#!/usr/bin/env Rscript
## 03_quick_contrast.R ---------------------------------------------------------
## Escape hatch for exploratory work: run a contrast straight from the command
## line, no YAML needed. The resolved config is still written next to the
## results, so an ad-hoc run stays reproducible.
##
##   # tumour vs normal in three projects
##   Rscript scripts/03_quick_contrast.R --projects TCGA-BRCA,TCGA-LUAD,TCGA-KIRC
##
##   # PAM50 subtype contrast
##   Rscript scripts/03_quick_contrast.R --projects TCGA-BRCA \
##       --group-by subtype_BRCA_Subtype_PAM50 --treatment Basal --reference LumA
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
source("R/load_framework.R")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--projects", type = "character", help = "comma-separated project ids"),
  make_option("--group-by", type = "character", default = "tissue_class",
              help = "grouping column [default %default]"),
  make_option("--treatment", type = "character", default = "tumor", help = "numerator level"),
  make_option("--reference", type = "character", default = "normal", help = "denominator level"),
  make_option("--paired", action = "store_true", default = FALSE,
              help = "block on patient (only meaningful for tumour/normal pairs)"),
  make_option("--omics", type = "character", default = "rnaseq,mutation",
              help = "layers to include [default %default]"),
  make_option("--cache-dir", type = "character", default = "cache"),
  make_option("--results-dir", type = "character", default = "results")
)))
if (is.null(opt$projects)) stop("--projects is required", call. = FALSE)

projects <- trimws(strsplit(opt$projects, ",")[[1]])
layers   <- trimws(strsplit(opt$omics, ",")[[1]])

tmp_config <- function(project) {
  cfg <- list(
    project     = project,
    cache_dir   = opt$`cache-dir`,
    results_dir = opt$`results-dir`,
    omics       = stats::setNames(as.list(names(OMIC_SPECS) %in% layers), names(OMIC_SPECS)),
    cohort      = list(group_by = opt$`group-by`, treatment = opt$treatment,
                       reference = opt$reference, paired = opt$paired)
  )
  path <- tempfile(fileext = ".yml")
  yaml::write_yaml(cfg, path)
  path
}

results <- list()
for (p in projects) {
  log_step("project ", p)
  out <- try({
    cfg <- load_config(tmp_config(p))
    data <- fetch_all(cfg)
    stopifnot(!is.null(data$rnaseq))
    coh <- build_cohort(data$rnaseq, cfg, data$clinical, data$subtypes)
    de  <- run_deseq2(coh, cfg)
    if (!is.null(data$mutation)) {
      analyse_mutations(cohort_maf(data$mutation, coh$se), cfg)
    }
    snapshot_run(cfg)
    build_report(cfg, cohort_summary = coh$summary, de_table = de$table)
    data.frame(project = p,
               n_treatment = coh$summary$n_samples[coh$summary$role == "treatment"],
               n_reference = coh$summary$n_samples[coh$summary$role == "reference"],
               n_significant = sum(de$table$significant),
               results = cfg$results_dir, stringsAsFactors = FALSE)
  }, silent = TRUE)

  if (inherits(out, "try-error")) {
    log_warn("project ", p, " failed: ", as.character(out))
    results[[p]] <- data.frame(project = p, n_treatment = NA, n_reference = NA,
                               n_significant = NA, results = "FAILED", stringsAsFactors = FALSE)
  } else {
    results[[p]] <- out
  }
}

summary_df <- do.call(rbind, results)
print(summary_df, row.names = FALSE)
write_table(summary_df, file.path(opt$`results-dir`, "batch_summary.tsv"))
