#!/usr/bin/env Rscript
## 00_explore_projects.R -------------------------------------------------------
## Discover what a project offers before committing to a download.
##
##   Rscript scripts/00_explore_projects.R --list
##   Rscript scripts/00_explore_projects.R --list --all-programs
##   Rscript scripts/00_explore_projects.R --project TCGA-BRCA
##
## For --project, two things are produced:
##   1. omics_availability.tsv + figures 01/02 - how many samples per omic
##      layer and sample type (tumour, normal, ...), from the GDC file
##      listing. No assay data is downloaded for this.
##   2. available_contrasts.tsv + figures 03_* - every clinical and molecular
##      subtype field usable as cohort$group_by, with the exact level
##      spelling and sample/patient counts on each side. This needs the
##      clinical + subtype tables, which are a light, fast download (a few
##      hundred KB) - still no assay data. Skip it with --no-clinical.
##
## These are a PREVIEW: counts are not de-duplicated or cohort-filtered the
## way build_cohort() will do it once you actually run an analysis, so they
## can differ by a handful of samples from the final cohort.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
source("R/load_framework.R")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--list", action = "store_true", default = FALSE,
              help = "list projects with their case count"),
  make_option("--all-programs", action = "store_true", default = FALSE,
              help = "with --list: include every GDC program (TARGET, CPTAC, MMRF, ...), not just TCGA-* [default: TCGA only]"),
  make_option("--project", type = "character", default = NULL,
              help = "project id to scan, e.g. TCGA-BRCA"),
  make_option("--omics", type = "character", default = "rnaseq,mutation,methylation",
              help = "comma-separated layers to scan [default %default]"),
  make_option("--out", type = "character", default = "results/scan",
              help = "output directory [default %default]"),
  make_option("--cache-dir", type = "character", default = "cache",
              help = "where to cache the clinical/subtype lookups; shared with the full pipeline's cache [default %default]"),
  make_option("--no-clinical", action = "store_true", default = FALSE,
              help = "skip clinical/subtype profiling - availability figures only (old behaviour)")
)))

## --- --list -------------------------------------------------------------------
if (opt$list) {
  p <- list_tcga_projects(all_programs = opt$`all-programs`)
  print(p[, c("project", "name")], right = FALSE)
  ensure_dir(opt$out)
  out_file <- if (opt$`all-programs`) "gdc_projects.tsv" else "tcga_projects.tsv"
  write_table(p, file.path(opt$out, out_file))
  if (opt$`all-programs`) {
    log_msg("NOTE: only project ids starting with 'TCGA-' are supported end-to-end by ",
            "this framework - tumour/normal and sample-type calls rely on the TCGA ",
            "barcode format (tcga_patient()/tcga_sample_type() in R/utils.R). Other ",
            "GDC programs (", paste(setdiff(unique(p$program), "TCGA"), collapse = ", "),
            ") are listed for reference; using one would need matching barcode/submitter-id logic.")
  }
  quit(status = 0)
}

if (is.null(opt$project)) stop("provide --project or --list", call. = FALSE)

## --- omics availability (unchanged from before) --------------------------------
omics <- trimws(strsplit(opt$omics, ",")[[1]])
avail <- scan_availability(opt$project, omics)
if (is.null(avail)) quit(status = 1)

cfg <- list(project = opt$project,
            cache_dir   = ensure_dir(file.path(opt$`cache-dir`, opt$project)),
            results_dir = ensure_dir(file.path(opt$out, opt$project)),
            figures_dir = ensure_dir(file.path(opt$out, opt$project, "figures")),
            tables_dir  = ensure_dir(file.path(opt$out, opt$project, "tables")))

## --- clinical / subtype categories (new, additive; skip with --no-clinical) ----
ann <- NULL
if (!opt$`no-clinical`) {
  clinical <- try(fetch_clinical(cfg), silent = TRUE)
  if (inherits(clinical, "try-error")) { log_warn("clinical fetch failed: ", as.character(clinical)); clinical <- NULL }
  subtypes <- try(fetch_subtypes(cfg), silent = TRUE)
  if (inherits(subtypes, "try-error")) { log_warn("subtype fetch failed: ", as.character(subtypes)); subtypes <- NULL }
  ann <- scan_annotation(avail, clinical, subtypes)
}

explore_project(cfg, avail = avail, sample_annotation = ann)

if (!is.null(ann) && "tissue_class" %in% names(ann) &&
    length(unique(stats::na.omit(ann$tissue_class))) >= 2) {
  log_step("checking known batch fields against tissue_class (tumor vs normal)")
  for (bf in c("tss", "plate", "center")) {
    if (!bf %in% names(ann)) next
    confounded <- covariate_is_confounded(ann[[bf]], ann$tissue_class, bf)
    if (!confounded) log_ok(bf, ": not confounded with tissue_class - safe to add as a covariate")
  }
}

if (!is.null(ann)) {
  cats <- summarise_categories(ann)
  if (!is.null(cats)) {
    write_table(cats, file.path(cfg$tables_dir, "available_contrasts.tsv"))

    ## compact per-variable overview for the console: which columns exist,
    ## from where, and how many levels each has
    ov <- do.call(rbind, lapply(split(cats, cats$variable), function(d) {
      data.frame(variable = d$variable[1], source = d$source[1], n_levels = nrow(d),
                 top_levels = paste(sprintf("%s(n=%d)", utils::head(d$level, 3),
                                            utils::head(d$n_samples, 3)), collapse = ", "),
                 stringsAsFactors = FALSE)
    }))
    ov <- ov[order(ov$source, ov$variable), ]
    rownames(ov) <- NULL
    log_step("grouping variables available for ", opt$project)
    print(ov, row.names = FALSE)
    log_ok("full level-by-level counts: ", file.path(cfg$tables_dir, "available_contrasts.tsv"))
  } else {
    log_msg("no usable categorical grouping variables found (clinical table may be sparse for this project)")
  }
}

print(avail, row.names = FALSE)
log_ok("done - see ", cfg$results_dir)
