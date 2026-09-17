## gdc_query.R ----------------------------------------------------------------
## Thin, opinionated wrappers around TCGAbiolinks::GDCquery.
##
## The GDC harmonised database changes its workflow labels from time to time.
## Centralising the query definitions here means a single edit keeps the whole
## framework working instead of hunting through a dozen scripts.
## -----------------------------------------------------------------------------

#' Supported omics layers and their GDC coordinates (harmonised / hg38).
OMIC_SPECS <- list(
  rnaseq = list(
    label         = "Gene expression (STAR raw counts)",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  ),
  mutation = list(
    label         = "Somatic mutations (MAF, masked)",
    data.category = "Simple Nucleotide Variation",
    data.type     = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  ),
  methylation = list(
    label         = "DNA methylation (beta values)",
    data.category = "DNA Methylation",
    data.type     = "Methylation Beta Value",
    platform      = "Illumina Human Methylation 450"
  ),
  mirna = list(
    label         = "miRNA expression",
    data.category = "Transcriptome Profiling",
    data.type     = "miRNA Expression Quantification",
    workflow.type = "BCGSC miRNA Profiling"
  ),
  cnv = list(
    label         = "Copy number (gene level)",
    data.category = "Copy Number Variation",
    data.type     = "Gene Level Copy Number",
    workflow.type = "ASCAT3"
  )
)

#' List all TCGA projects available on the GDC.
#'
#' @return data.frame with project id, name, primary site and case count
#' List projects available on the GDC.
#'
#' @param all_programs FALSE (default, unchanged behaviour): only the 33 TCGA
#'   disease-type projects, which is what this framework's barcode parsing
#'   (tcga_patient(), tcga_sample_type(), tcga_tissue_class() in R/utils.R)
#'   assumes. TRUE: every GDC program (TARGET, CPTAC, MMRF, HCMI, CGCI, ...),
#'   for reference - those use different submitter-id schemes and are NOT
#'   supported end-to-end by this framework without extending the barcode
#'   helpers in R/utils.R.
#' @return data.frame with project id, program, name, primary site, disease
#'   type and case count
list_tcga_projects <- function(all_programs = FALSE) {
  need_pkg("TCGAbiolinks")
  p <- TCGAbiolinks::getGDCprojects()
  if (!all_programs) p <- p[grepl("^TCGA-", p$project_id), ]
  out <- data.frame(
    project      = p$project_id,
    program      = sub("-.*$", "", p$project_id),
    name         = p$name,
    primary_site = vapply(p$primary_site, function(x) paste(unlist(x), collapse = "; "), character(1)),
    disease_type = vapply(p$disease_type, function(x) paste(unlist(x), collapse = "; "), character(1)),
    barcode_compatible = grepl("^TCGA-", p$project_id),
    stringsAsFactors = FALSE
  )
  out[order(out$program, out$project), ]
}

#' Build a GDCquery for one omic layer of one project.
#'
#' @param project GDC project id, e.g. "TCGA-LUAD"
#' @param omic one of names(OMIC_SPECS)
#' @param sample_types optional character vector of GDC sample type labels
#' @param barcodes optional character vector restricting the query
build_query <- function(project, omic, sample_types = NULL, barcodes = NULL) {
  need_pkg("TCGAbiolinks")
  spec <- OMIC_SPECS[[omic]]
  if (is.null(spec)) log_die("unknown omic '", omic, "'. Available: ",
                             paste(names(OMIC_SPECS), collapse = ", "))

  args <- list(project = project,
               data.category = spec$data.category,
               data.type     = spec$data.type,
               access        = "open")
  if (!is.null(spec$workflow.type)) args$workflow.type <- spec$workflow.type
  if (!is.null(spec$platform))      args$platform      <- spec$platform
  if (!is.null(sample_types))       args$sample.type   <- sample_types
  if (!is.null(barcodes))           args$barcode       <- barcodes

  q <- try(do.call(TCGAbiolinks::GDCquery, args), silent = TRUE)
  if (inherits(q, "try-error")) {
    log_warn("no ", omic, " data for ", project, " (", as.character(q), ")")
    return(NULL)
  }
  q
}

#' Inventory of what exists on the GDC for a project, without downloading it.
#'
#' Answers "how many samples do I actually have, and of which type" before you
#' commit to a multi-gigabyte download.
#'
#' @return data.frame: omic, sample_type, n_samples, n_patients
scan_availability <- function(project, omics = names(OMIC_SPECS)) {
  need_pkg("TCGAbiolinks")
  rows <- list()
  for (om in omics) {
    log_msg("scanning ", project, " / ", om)
    q <- build_query(project, om)
    if (is.null(q)) next
    res <- try(TCGAbiolinks::getResults(q), silent = TRUE)
    if (inherits(res, "try-error") || is.null(res) || nrow(res) == 0) next

    bc <- res$cases %||% res$sample.submitter_id
    if (is.null(bc)) next
    bc <- as.character(bc)

    df <- data.frame(
      omic        = om,
      barcode     = bc,
      patient     = tcga_patient(bc),
      sample_type = tcga_sample_type(bc),
      tissue      = as.character(tcga_tissue_class(bc)),
      tss         = tcga_tss(bc),
      plate       = tcga_plate(bc),
      center      = tcga_center(bc),
      stringsAsFactors = FALSE
    )
    rows[[om]] <- df
  }
  if (length(rows) == 0) {
    log_warn("no data found for ", project)
    return(NULL)
  }
  long <- do.call(rbind, rows)
  key <- paste(long$omic, long$sample_type, long$tissue, sep = "\r")
  agg <- do.call(rbind, lapply(split(long, key), function(d) {
    data.frame(omic = d$omic[1], sample_type = d$sample_type[1], tissue = d$tissue[1],
               n_samples  = length(unique(d$barcode)),
               n_patients = length(unique(d$patient)),
               stringsAsFactors = FALSE)
  }))
  rownames(agg) <- NULL
  agg <- agg[order(agg$omic, -agg$n_samples), ]
  rownames(agg) <- NULL
  ## keep the sample-level table attached: the overlap plot needs it.
  ## (set AFTER subsetting - data.frame subsetting drops custom attributes)
  attr(agg, "long") <- long
  agg
}
