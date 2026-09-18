## utils.R --------------------------------------------------------------------
## Small helpers shared across the framework: logging, TCGA barcode parsing,
## aliquot de-duplication and safe IO.
## -----------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) return(y)
  if (length(x) == 0L) return(y)
  if (length(x) == 1L && is.na(x)) return(y)
  x
}

.ts <- function() format(Sys.time(), "%H:%M:%S")

log_msg  <- function(...) message(sprintf("[%s] %s", .ts(), paste0(...)))
log_step <- function(...) message(sprintf("\n[%s] ==> %s", .ts(), paste0(...)))
log_ok   <- function(...) message(sprintf("[%s]  ok  %s", .ts(), paste0(...)))
log_warn <- function(...) warning(sprintf("[%s] WARN %s", .ts(), paste0(...)), call. = FALSE, immediate. = TRUE)
log_die  <- function(...) stop(sprintf("[%s] FAIL %s", .ts(), paste0(...)), call. = FALSE)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(normalizePath(path, mustWork = FALSE))
}

need_pkg <- function(pkg, why = "") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_die(sprintf("package '%s' is required%s. Run: Rscript install/install_dependencies.R",
                    pkg, if (nzchar(why)) paste0(" ", why) else ""))
  }
  invisible(TRUE)
}

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

## --- TCGA barcode handling ---------------------------------------------------
## Barcode layout: TCGA-02-0001-01C-01D-0182-01
##                 [ 1 ][2 ][ 3 ][4 ][5 ][6 ][7]
##   1-12 : patient (submitter_id)
##   1-15 : sample  (patient + sample type code)
##   14-15: sample type code (01 = primary tumour, 11 = solid tissue normal, ...)
##   16   : vial ; 17-18 portion ; 19 analyte ; 20-23 plate ; 24-25 centre

SAMPLE_TYPE_CODES <- data.frame(
  code  = c("01", "02", "03", "05", "06", "07", "10", "11", "12", "13", "14", "20"),
  label = c("Primary Tumor", "Recurrent Tumor", "Primary Blood Derived Cancer - Peripheral Blood",
            "Additional - New Primary", "Metastatic", "Additional Metastatic",
            "Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal",
            "EBV Immortalized Normal", "Bone Marrow Normal", "Control Analyte"),
  class = c("tumor", "tumor", "tumor", "tumor", "tumor", "tumor",
            "normal", "normal", "normal", "normal", "normal", "control"),
  stringsAsFactors = FALSE
)

tcga_patient <- function(barcode) substr(as.character(barcode), 1, 12)
tcga_sample  <- function(barcode) substr(as.character(barcode), 1, 15)
tcga_code    <- function(barcode) substr(as.character(barcode), 14, 15)

tcga_sample_type <- function(barcode) {
  idx <- match(tcga_code(barcode), SAMPLE_TYPE_CODES$code)
  out <- SAMPLE_TYPE_CODES$label[idx]
  out[is.na(out)] <- "Unknown"
  out
}

#' Coarse tumour / normal class derived from the barcode itself.
#' Never trust a single clinical column for this - the barcode is authoritative.
tcga_tissue_class <- function(barcode) {
  idx <- match(tcga_code(barcode), SAMPLE_TYPE_CODES$code)
  out <- SAMPLE_TYPE_CODES$class[idx]
  out[is.na(out)] <- "unknown"
  factor(out, levels = c("tumor", "normal", "control", "unknown"))
}

#' Drop replicate aliquots so that each (patient, sample type) appears once.
#'
#' GDC ships several aliquots for some samples (different vial / portion /
#' plate). Keeping them silently inflates group sizes and breaks the
#' independence assumption of every test downstream. We keep the aliquot with
#' the lexicographically smallest barcode tail, which corresponds to the
#' earliest vial/portion - the convention used by most TCGA papers.
#'
#' @param barcodes character vector of full aliquot barcodes
#' @return logical vector, TRUE for the aliquots to keep
dedup_aliquots <- function(barcodes) {
  barcodes <- as.character(barcodes)
  key  <- tcga_sample(barcodes)
  tail <- substr(barcodes, 16, nchar(barcodes))
  ord  <- order(key, tail)
  keep <- rep(FALSE, length(barcodes))
  keep[ord[!duplicated(key[ord])]] <- TRUE
  n_drop <- sum(!keep)
  if (n_drop > 0) log_msg(sprintf("de-duplicated %d replicate aliquot(s)", n_drop))
  keep
}

## --- IO ----------------------------------------------------------------------

write_table <- function(x, path, ...) {
  ensure_dir(dirname(path))
  df <- flatten_list_columns(as.data.frame(x))
  utils::write.table(df, file = path, sep = "\t",
                     quote = FALSE, row.names = FALSE, ...)
  log_ok("wrote ", path)
  invisible(path)
}

#' Coerce any list-type column of a data.frame to plain character.
#'
#' TCGAbiolinks' GDCprepare() sometimes returns colData with list columns -
#' a case can have more than one recorded value for a field (e.g. several
#' treatments), so it's stored as one list per sample. Harmless for
#' Bioconductor objects, but base write.table() cannot serialise it
#' ("unimplemented type 'list' in 'EncodeElement'"), and as.character() on a
#' whole list column produces deparsed noise like 'c("a", "b")' rather than
#' a readable value. This is the one place that gets called wherever GDC
#' metadata is turned into a plain data.frame, so nothing downstream has to
#' know or care that the column used to be a list.
flatten_list_columns <- function(df) {
  is_list_col <- vapply(df, is.list, logical(1))
  if (!any(is_list_col)) return(df)
  for (cc in names(df)[is_list_col]) {
    df[[cc]] <- vapply(df[[cc]], function(v) {
      if (is.null(v) || length(v) == 0 || all(is.na(v))) return(NA_character_)
      paste(as.character(v), collapse = "; ")
    }, character(1))
  }
  df
}

save_plot <- function(plot, path, width = 8, height = 6, dpi = 300) {
  ensure_dir(dirname(path))
  ok <- try({
    ggplot2::ggsave(filename = path, plot = plot, width = width,
                    height = height, dpi = dpi, limitsize = FALSE)
  }, silent = TRUE)
  if (inherits(ok, "try-error")) {
    log_warn("could not save plot ", path, ": ", as.character(ok))
    return(invisible(NULL))
  }
  log_ok("wrote ", path)
  invisible(path)
}

#' Cache any expensive object as .rds, recomputing only when asked.
cached <- function(path, expr, refresh = FALSE) {
  if (!refresh && file.exists(path)) {
    log_msg("cache hit: ", basename(path))
    return(readRDS(path))
  }
  obj <- expr   # promise is evaluated here, only when the cache misses
  ensure_dir(dirname(path))
  saveRDS(obj, path)
  log_ok("cached: ", path)
  obj
}

#' Sanitise a string so it can be used as a file name / factor level.
slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

#' TCGA barcode positions beyond the sample-type code, decoded when the full
#' aliquot barcode is available (e.g. "TCGA-02-0001-01C-01D-0182-01"):
#'   [1-4]=TCGA [6-7]=TSS [9-12]=participant [14-15]=sample [16]=vial
#'   [18-19]=portion [20]=analyte [22-25]=plate [27-28]=center
#'
#' These are the standard proxies for TCGA batch effects: samples sharing a
#' Tissue Source Site (TSS) were collected at the same institution, samples
#' sharing a plate were processed together, and both are well-documented
#' sources of non-biological variation in TCGA data (see e.g. the original
#' MAQC/TCGA batch-effect reports). NA is returned where the barcode is too
#' short to contain that field, rather than an error - not every layer's
#' barcode always carries the full 28 characters.

.tcga_field <- function(barcode, start, stop) {
  barcode <- as.character(barcode)
  out <- rep(NA_character_, length(barcode))
  ok  <- nchar(barcode) >= stop
  out[ok] <- substr(barcode[ok], start, stop)
  out
}

tcga_tss    <- function(barcode) .tcga_field(barcode, 6, 7)     # Tissue Source Site
tcga_vial   <- function(barcode) .tcga_field(barcode, 16, 16)
tcga_portion<- function(barcode) .tcga_field(barcode, 18, 19)
tcga_analyte<- function(barcode) .tcga_field(barcode, 20, 20)
tcga_plate  <- function(barcode) .tcga_field(barcode, 22, 25)
tcga_center <- function(barcode) .tcga_field(barcode, 27, 28)

#' TRUE if every level of `covariate` co-occurs with every level of `group`
#' - a necessary (not sufficient) condition for an identifiable design.
#' When it fails, the covariate is either fully or partially aliased with
#' the contrast: DESeq2/limma will either drop it, report it as NA, or the
#' whole design matrix can become rank-deficient. Warns with the offending
#' levels named, so the fix is obvious, and returns TRUE so the caller can
#' drop the covariate rather than let the fit fail opaquely later.
covariate_is_confounded <- function(covariate, group, name = "covariate") {
  tab <- table(as.character(covariate), as.character(group))
  if (any(tab == 0)) {
    culprits <- rownames(tab)[apply(tab == 0, 1, any)]
    log_warn(sprintf(
      "'%s' is confounded with the contrast - level(s) {%s} of it appear in only one group. ",
      name, paste(culprits, collapse = ", ")),
      "Dropped from the design. If this batch genuinely only exists on one side, ",
      "it cannot be separated from the biological effect; restrict the cohort with ",
      "cohort$filters instead of trying to adjust for it.")
    return(TRUE)
  }
  FALSE
}

#' Filter a list of candidate covariates down to the ones that are usable:
#' present in cd, have >=2 levels, and are not confounded with `group_col`.
#' Shared by the RNA-seq, miRNA and methylation design-building code so the
#' three don't drift apart.
clean_covariates <- function(covars, cd, group_col = "cohort_group") {
  covars <- unique(covars[covars %in% names(cd)])
  Filter(function(cv) {
    lv <- droplevels(factor(as.character(cd[[cv]])))
    if (nlevels(lv) < 2) { log_warn("covariate '", cv, "' is constant - dropped"); return(FALSE) }
    if (covariate_is_confounded(lv, cd[[group_col]], cv)) return(FALSE)
    TRUE
  }, covars)
}
#' column, prefixing every added column so its source is always visible
#' (clin_, subtype_, ...). Shared by the full per-sample annotation
#' (cohort.R) and the lightweight scan preview (explore.R) so the two never
#' drift apart.
#'
#' @param df data.frame with a `patient` column
#' @param ... named patient-level tables, e.g. join_patient_tables(df, clin_ = clinical)
join_patient_tables <- function(df, ...) {
  tabs <- list(...)
  for (nm in names(tabs)) {
    tab <- tabs[[nm]]
    if (is.null(tab) || !"patient" %in% names(tab)) next
    tab <- tab[!duplicated(tab$patient), , drop = FALSE]
    keep <- setdiff(names(tab), c(names(df), "patient"))
    if (length(keep) == 0) next
    add <- tab[match(df$patient, tab$patient), keep, drop = FALSE]
    names(add) <- paste0(nm, names(add))
    df <- cbind(df, add)
  }
  df
}

theme_tcga <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey30")
    )
}
