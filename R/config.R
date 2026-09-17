## config.R -------------------------------------------------------------------
## Every analysis in this framework is described by a single YAML file, so a
## result can always be traced back to the exact parameters that produced it.
## -----------------------------------------------------------------------------

DEFAULT_CONFIG <- list(
  project      = NULL,          # e.g. "TCGA-BRCA"
  cache_dir    = "cache",
  results_dir  = "results",
  omics        = list(
    rnaseq      = TRUE,
    mutation    = TRUE,
    methylation = FALSE,        # 450k arrays are heavy; opt-in
    cnv         = FALSE,        # gene-level copy number (ASCAT3)
    mirna       = FALSE         # miRNA expression
  ),
  cohort = list(
    group_by      = "tissue_class",   # colData column, or "tissue_class"/"sample_type"
    reference     = "normal",         # denominator of the contrast
    treatment     = "tumor",          # numerator of the contrast
    keep_levels   = NULL,             # optional explicit whitelist of levels
    filters       = NULL,             # named list: colData column -> allowed values
    paired        = FALSE,            # block on patient (tumour/normal pairs)
    min_group_n   = 3,
    drop_duplicate_aliquots = TRUE
  ),
  deseq2 = list(
    assay             = "unstranded", # raw integer counts from STAR
    min_count         = 10,
    min_samples       = NULL,         # default: size of the smallest group
    protein_coding_only = TRUE,
    covariates        = NULL,         # e.g. list("tss", "clin_gender") - see docs/config_reference.md
    auto_sva          = FALSE,        # opt-in: estimate hidden batch effects with SVA (needs the 'sva' package)
    n_sv              = NULL,         # surrogate variable count; NULL = auto-estimated
    alpha             = 0.05,
    lfc_threshold     = 1,
    shrink            = TRUE,
    top_n_heatmap     = 50,
    run_enrichment    = FALSE
  ),
  mutation = list(
    top_genes   = 25,
    min_mut     = 5,              # minimum mutated samples for group comparison
    compare_groups = TRUE
  ),
  cnv = list(
    assay            = NULL,    # auto-detected ("copy_number", GISTIC scores, ...)
    ploidy_correct   = TRUE,    # call gains/losses relative to sample ploidy
    gain_ratio       = 1.4,     # CN/ploidy >= this  -> gain
    amp_ratio        = 2.0,     # CN/ploidy >= this  -> amplification
    loss_ratio       = 0.6,     # CN/ploidy <= this  -> loss
    min_freq         = 0.05,    # only test genes altered in >=5% of a group
    top_genes        = 25,
    correlate_expression = TRUE # dosage effect vs the RNA-seq layer
  ),
  mirna = list(
    min_count     = 10,
    min_samples   = NULL,
    covariates    = NULL,
    alpha         = 0.05,
    lfc_threshold = 1,
    shrink        = TRUE
  ),
  methylation = list(
    platform        = "Illumina Human Methylation 450",
    drop_sex_chr    = TRUE,
    max_na_fraction = 0.1,
    covariates      = NULL,           # e.g. list("tss") - dropped automatically if constant/confounded
    p_adjust        = 0.05,
    delta_beta      = 0.2,
    top_n_heatmap   = 50
  ),
  seed = 1234
)

#' Recursively merge a user config over the defaults.
merge_config <- function(default, user) {
  for (nm in names(user)) {
    if (is.list(default[[nm]]) && is.list(user[[nm]]) && !is.null(names(user[[nm]]))) {
      default[[nm]] <- merge_config(default[[nm]], user[[nm]])
    } else {
      default[[nm]] <- user[[nm]]
    }
  }
  default
}

#' Read and validate an analysis config.
#'
#' @param path path to a YAML file
#' @return a validated config list with absolute output paths
load_config <- function(path) {
  need_pkg("yaml")
  if (!file.exists(path)) log_die("config file not found: ", path)
  user <- yaml::read_yaml(path)
  cfg  <- merge_config(DEFAULT_CONFIG, user)

  if (is.null(cfg$project) || !nzchar(cfg$project)) {
    log_die("config must define 'project' (e.g. TCGA-BRCA)")
  }
  if (!grepl("^TCGA-", cfg$project)) {
    log_warn("project '", cfg$project, "' does not look like a TCGA project id")
  }
  if (is.null(cfg$cohort$treatment) || is.null(cfg$cohort$reference)) {
    log_die("config must define cohort$treatment and cohort$reference")
  }
  if (identical(cfg$cohort$treatment, cfg$cohort$reference)) {
    log_die("cohort$treatment and cohort$reference must differ")
  }

  cfg$contrast_name <- sprintf("%s_%s_vs_%s", slug(cfg$project),
                               slug(cfg$cohort$treatment), slug(cfg$cohort$reference))
  cfg$cache_dir   <- ensure_dir(file.path(cfg$cache_dir, cfg$project))
  cfg$results_dir <- ensure_dir(file.path(cfg$results_dir, cfg$contrast_name))
  cfg$figures_dir <- ensure_dir(file.path(cfg$results_dir, "figures"))
  cfg$tables_dir  <- ensure_dir(file.path(cfg$results_dir, "tables"))
  cfg$config_path <- normalizePath(path)

  set.seed(cfg$seed %||% 1234)
  cfg
}

#' Write the resolved config plus the session info next to the results.
#' This is what makes a run reproducible six months later.
snapshot_run <- function(cfg) {
  need_pkg("yaml")
  dump <- cfg[setdiff(names(cfg), c("config_path"))]
  yaml::write_yaml(dump, file.path(cfg$results_dir, "run_config.resolved.yml"))
  si <- utils::capture.output(utils::sessionInfo())
  writeLines(si, file.path(cfg$results_dir, "sessionInfo.txt"))
  log_ok("run snapshot written to ", cfg$results_dir)
  invisible(TRUE)
}
