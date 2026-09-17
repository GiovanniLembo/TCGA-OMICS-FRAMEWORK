## download.R -----------------------------------------------------------------
## Download once, reuse forever. Every layer is prepared into a Bioconductor
## object and cached as .rds under <cache_dir>/<PROJECT>/.
## -----------------------------------------------------------------------------

cache_path <- function(cfg, name) file.path(cfg$cache_dir, name)

#' Download + prepare one omic layer, with caching.
#'
#' @param cfg config list from load_config()
#' @param omic one of names(OMIC_SPECS)
#' @param force re-download even if a cache entry exists
#' @param files_per_chunk chunk size for the GDC API (lower it on flaky networks)
#' @return SummarizedExperiment (rnaseq / methylation / cnv) or data.frame (mutation)
fetch_omic <- function(cfg, omic, force = FALSE, files_per_chunk = 20) {
  need_pkg("TCGAbiolinks")
  rds <- cache_path(cfg, paste0(omic, ".rds"))
  if (file.exists(rds) && !force) {
    log_msg("cache hit: ", basename(rds))
    return(readRDS(rds))
  }

  log_step("fetching ", omic, " for ", cfg$project)
  q <- build_query(cfg$project, omic)
  if (is.null(q)) return(NULL)

  gdc_dir <- ensure_dir(file.path(cfg$cache_dir, "GDCdata"))
  ok <- try(TCGAbiolinks::GDCdownload(q, method = "api",
                                      files.per.chunk = files_per_chunk,
                                      directory = gdc_dir), silent = TRUE)
  if (inherits(ok, "try-error")) {
    log_warn("download failed for ", omic, ": ", as.character(ok))
    log_warn("tip: lower files_per_chunk, or re-run - GDCdownload resumes from disk")
    return(NULL)
  }

  obj <- try(TCGAbiolinks::GDCprepare(q, directory = gdc_dir, summarizedExperiment = TRUE),
             silent = TRUE)
  if (inherits(obj, "try-error")) {
    log_warn("GDCprepare failed for ", omic, ": ", as.character(obj))
    return(NULL)
  }

  obj <- annotate_object(obj, omic)
  saveRDS(obj, rds)
  log_ok("cached ", omic, " -> ", rds)
  obj
}

#' Add barcode-derived columns that we rely on downstream.
annotate_object <- function(obj, omic) {
  if (omic == "mutation") {
    obj <- as.data.frame(obj)
    bc  <- obj$Tumor_Sample_Barcode
    obj$patient      <- tcga_patient(bc)
    obj$sample_type  <- tcga_sample_type(bc)
    obj$tissue_class <- as.character(tcga_tissue_class(bc))
    obj$tss          <- tcga_tss(bc)
    obj$plate        <- tcga_plate(bc)
    obj$center       <- tcga_center(bc)
    obj$portion      <- tcga_portion(bc)
    return(obj)
  }
  if (omic == "mirna" && !inherits(obj, "SummarizedExperiment")) {
    ## GDC ships miRNA as a wide data.frame - reshape so the cohort machinery applies
    obj <- mirna_to_se(obj)
  }
  if (inherits(obj, "SummarizedExperiment")) {
    need_pkg("SummarizedExperiment")
    cd <- SummarizedExperiment::colData(obj)
    bc <- colnames(obj)
    cd$barcode_full  <- bc
    cd$patient       <- tcga_patient(bc)
    cd$sample_id     <- tcga_sample(bc)
    cd$sample_type   <- cd$sample_type %||% tcga_sample_type(bc)
    cd$tissue_class  <- as.character(tcga_tissue_class(bc))
    ## batch proxies - populated whenever the column name carries the full
    ## aliquot barcode; NA (and dropped as constant) otherwise, so this is
    ## harmless even for layers where it's not available.
    cd$tss           <- tcga_tss(bc)
    cd$plate         <- tcga_plate(bc)
    cd$center        <- tcga_center(bc)
    cd$portion       <- tcga_portion(bc)
    SummarizedExperiment::colData(obj) <- cd
  }
  obj
}

#' Patient-level clinical table (indexed / harmonised clinical data).
fetch_clinical <- function(cfg, force = FALSE) {
  need_pkg("TCGAbiolinks")
  rds <- cache_path(cfg, "clinical.rds")
  if (file.exists(rds) && !force) return(readRDS(rds))

  log_step("fetching clinical metadata for ", cfg$project)
  cl <- try(TCGAbiolinks::GDCquery_clinic(project = cfg$project, type = "clinical"),
            silent = TRUE)
  if (inherits(cl, "try-error")) {
    log_warn("clinical download failed: ", as.character(cl))
    return(NULL)
  }
  cl$patient <- cl$submitter_id
  saveRDS(cl, rds)
  log_ok("cached clinical -> ", rds)
  cl
}

#' Published molecular subtypes curated by TCGAbiolinks (PAM50, CMS, ...).
#' Returns NULL when the disease has no curated subtype table.
fetch_subtypes <- function(cfg, force = FALSE) {
  need_pkg("TCGAbiolinks")
  rds <- cache_path(cfg, "subtypes.rds")
  if (file.exists(rds) && !force) return(readRDS(rds))

  disease <- tolower(sub("^TCGA-", "", cfg$project))
  st <- try(TCGAbiolinks::TCGAquery_subtype(tumor = disease), silent = TRUE)
  if (inherits(st, "try-error") || is.null(st)) {
    log_msg("no curated subtype table for ", cfg$project)
    return(NULL)
  }
  st <- as.data.frame(st)
  if ("patient" %in% names(st)) st$patient <- as.character(st$patient)
  saveRDS(st, rds)
  log_ok("cached subtypes -> ", rds)
  st
}

#' Fetch every layer requested in the config.
#' @return named list; missing layers are NULL rather than an error
fetch_all <- function(cfg, force = FALSE) {
  wanted <- names(cfg$omics)[vapply(cfg$omics, isTRUE, logical(1))]
  out <- list()
  for (om in wanted) out[[om]] <- fetch_omic(cfg, om, force = force)
  out$clinical <- fetch_clinical(cfg, force = force)
  out$subtypes <- fetch_subtypes(cfg, force = force)
  out
}
