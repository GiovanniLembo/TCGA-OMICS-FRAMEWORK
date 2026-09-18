## cohort.R -------------------------------------------------------------------
## Turns a config block into a validated two-group comparison.
##
## Design rule: the cohort is defined ONCE and the same definition is applied to
## every omic layer, so the RNA-seq, mutation and methylation results always
## describe the same patients.
## -----------------------------------------------------------------------------

#' Attach clinical + subtype annotation to the colData of a SummarizedExperiment.
#'
#' Clinical data is patient-level, assays are sample-level, so we join on the
#' 12-character patient barcode.
annotate_colData <- function(se, clinical = NULL, subtypes = NULL) {
  need_pkg("SummarizedExperiment")
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cd$patient <- cd$patient %||% tcga_patient(rownames(cd))
  cd <- join_patient_tables(cd, clin_ = clinical, subtype_ = subtypes)
  cd <- flatten_list_columns(cd)
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd, row.names = rownames(cd))
  se
}

#' Resolve the grouping variable into a plain character vector.
#'
#' "tissue_class" and "sample_type" are always available because they are
#' derived from the barcode; anything else must be a colData column.
resolve_group <- function(cd, group_by) {
  if (!group_by %in% names(cd)) {
    log_die(sprintf("group_by '%s' is not a column of the sample annotation.\n  Available (first 40): %s",
                    group_by, paste(utils::head(names(cd), 40), collapse = ", ")))
  }
  as.character(cd[[group_by]])
}

#' Apply config$cohort$filters, e.g. list(clin_gender = "female").
apply_filters <- function(cd, filters) {
  if (is.null(filters) || length(filters) == 0) return(rep(TRUE, nrow(cd)))
  keep <- rep(TRUE, nrow(cd))
  for (col in names(filters)) {
    if (!col %in% names(cd)) {
      log_warn("filter column '", col, "' not found - ignored")
      next
    }
    allowed <- unlist(filters[[col]])
    hit <- as.character(cd[[col]]) %in% as.character(allowed)
    log_msg(sprintf("filter %s in {%s}: %d/%d samples kept",
                    col, paste(allowed, collapse = ", "), sum(hit & keep), nrow(cd)))
    keep <- keep & hit
  }
  keep
}

#' Build the analysis cohort from a SummarizedExperiment.
#'
#' @return list(se = subsetted SE with cohort_group factor, summary = data.frame)
build_cohort <- function(se, cfg, clinical = NULL, subtypes = NULL) {
  need_pkg("SummarizedExperiment")
  cc <- cfg$cohort
  se <- annotate_colData(se, clinical, subtypes)
  cd <- as.data.frame(SummarizedExperiment::colData(se))

  log_step("building cohort: ", cc$treatment, " vs ", cc$reference,
           " (group_by = ", cc$group_by, ")")

  ## 1. drop replicate aliquots ------------------------------------------------
  keep <- rep(TRUE, nrow(cd))
  if (isTRUE(cc$drop_duplicate_aliquots)) {
    keep <- keep & dedup_aliquots(rownames(cd))
  }

  ## 2. arbitrary clinical filters ---------------------------------------------
  keep <- keep & apply_filters(cd, cc$filters)

  ## 3. group assignment -------------------------------------------------------
  grp <- resolve_group(cd, cc$group_by)
  grp[is.na(grp) | grp %in% c("", "NA", "not reported", "Not Reported")] <- NA

  wanted <- c(cc$reference, cc$treatment)
  in_grp <- grp %in% wanted
  if (!any(in_grp & keep)) {
    log_die(sprintf("no samples match the requested levels.\n  Observed levels of '%s': %s",
                    cc$group_by, paste(sort(unique(stats::na.omit(grp))), collapse = " | ")))
  }
  keep <- keep & in_grp

  ## 4. optional pairing -------------------------------------------------------
  if (isTRUE(cc$paired)) {
    pat  <- cd$patient
    tab  <- table(pat[keep], grp[keep])
    both <- rownames(tab)[rowSums(tab[, wanted, drop = FALSE] > 0) == 2]
    log_msg(sprintf("paired design: %d patient(s) have both groups", length(both)))
    if (length(both) < 2) log_die("paired analysis needs at least 2 complete pairs")
    keep <- keep & pat %in% both
  }

  se  <- se[, keep]
  cd  <- as.data.frame(SummarizedExperiment::colData(se))
  grp <- resolve_group(cd, cc$group_by)

  ## reference first: DESeq2/limma use the first level as the denominator
  cd$cohort_group <- factor(grp, levels = wanted)
  cd$patient      <- factor(cd$patient)
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd, row.names = rownames(cd))

  ## 5. sanity checks ----------------------------------------------------------
  n <- table(cd$cohort_group)
  log_msg("cohort sizes: ", paste(sprintf("%s=%d", names(n), as.integer(n)), collapse = ", "))
  if (any(n < cc$min_group_n)) {
    log_die(sprintf("group smaller than min_group_n (%d). Differential analysis on n<3 is not interpretable.",
                    cc$min_group_n))
  }

  summary_df <- data.frame(
    group    = names(n),
    role     = ifelse(names(n) == cc$reference, "reference", "treatment"),
    n_samples = as.integer(n),
    n_patients = as.integer(tapply(cd$patient, cd$cohort_group,
                                   function(x) length(unique(as.character(x))))[names(n)]),
    stringsAsFactors = FALSE
  )
  list(se = se, summary = summary_df, paired = isTRUE(cc$paired))
}

#' Subset a MAF data.frame to the patients of an existing cohort, carrying the
#' group labels across so mutations can be compared with the same contrast.
cohort_maf <- function(maf_df, cohort_se) {
  need_pkg("SummarizedExperiment")
  cd <- as.data.frame(SummarizedExperiment::colData(cohort_se))
  map <- stats::setNames(as.character(cd$cohort_group), as.character(cd$patient))
  maf_df$cohort_group <- unname(map[maf_df$patient])
  out <- maf_df[!is.na(maf_df$cohort_group), , drop = FALSE]
  log_msg(sprintf("MAF restricted to cohort: %d variants, %d patients",
                  nrow(out), length(unique(out$patient))))
  out
}
