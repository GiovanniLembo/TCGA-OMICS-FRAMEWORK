## explore.R ------------------------------------------------------------------
## "What do I actually have?" - the step that should always come before any
## differential analysis. Produces the sample-inventory figures and tables.
## -----------------------------------------------------------------------------

#' Barplot of sample counts per sample type for each omic layer.
plot_sample_types <- function(avail, project = "") {
  need_pkg("ggplot2")
  df <- avail
  ## avail has one row per (omic, sample_type, tissue), so sample_type repeats
  ## across rows - order the LEVELS by total count, not the column itself.
  totals <- tapply(df$n_samples, df$sample_type, sum)
  lvls <- names(sort(totals))
  df$sample_type <- factor(df$sample_type, levels = lvls)
  ggplot2::ggplot(df, ggplot2::aes(x = sample_type, y = n_samples, fill = tissue)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = n_samples), hjust = -0.15, size = 3) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ omic, scales = "free_x") +
    ggplot2::scale_fill_manual(values = c(tumor = "#C0392B", normal = "#2E86C1",
                                          control = "#7F8C8D", unknown = "#BDC3C7")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(title = paste("Sample availability", project),
                  subtitle = "Sample types per omic layer (GDC harmonised, open access)",
                  x = NULL, y = "Number of samples", fill = "Tissue") +
    theme_tcga()
}

#' How many patients are covered by which combination of omic layers.
#' A compact stand-in for an UpSet plot with no extra dependency.
plot_omics_overlap <- function(avail, project = "") {
  need_pkg("ggplot2")
  long <- attr(avail, "long")
  if (is.null(long)) return(NULL)

  tum <- long[long$tissue == "tumor", c("patient", "omic")]
  tum <- unique(tum)
  omics <- sort(unique(tum$omic))
  mat <- table(tum$patient, tum$omic) > 0
  combo <- apply(mat, 1, function(r) paste(omics[r], collapse = " + "))
  df <- as.data.frame(table(combo), stringsAsFactors = FALSE)
  names(df) <- c("combination", "n_patients")
  df <- df[order(-df$n_patients), ]
  df$combination <- factor(df$combination, levels = rev(df$combination))

  ggplot2::ggplot(df, ggplot2::aes(x = combination, y = n_patients)) +
    ggplot2::geom_col(fill = "#34495E") +
    ggplot2::geom_text(ggplot2::aes(label = n_patients), hjust = -0.15, size = 3) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(title = paste("Multi-omic coverage", project),
                  subtitle = "Tumour patients per available combination of layers",
                  x = NULL, y = "Number of patients") +
    theme_tcga()
}

#' Distribution of a clinical variable, split by tissue class when present.
plot_clinical_variable <- function(cd, column, fill_by = "tissue_class") {
  need_pkg("ggplot2")
  if (!column %in% names(cd)) return(NULL)
  v <- cd[[column]]
  v_num <- suppressWarnings(as.numeric(as.character(v)))
  numeric_like <- !all(is.na(v_num)) && length(unique(stats::na.omit(v_num))) > 10

  fill_ok <- fill_by %in% names(cd)
  base <- if (numeric_like) {
    d <- data.frame(value = v_num, fill = if (fill_ok) as.character(cd[[fill_by]]) else "all")
    d <- d[!is.na(d$value), , drop = FALSE]
    ggplot2::ggplot(d, ggplot2::aes(x = value, fill = fill)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.8, position = "identity") +
      ggplot2::labs(x = column, y = "Samples")
  } else {
    d <- data.frame(value = as.character(v), fill = if (fill_ok) as.character(cd[[fill_by]]) else "all")
    d$value[is.na(d$value) | d$value == ""] <- "not reported"
    keep <- names(sort(table(d$value), decreasing = TRUE))[seq_len(min(20, length(unique(d$value))))]
    d <- d[d$value %in% keep, , drop = FALSE]
    d$value <- factor(d$value, levels = rev(keep))
    ggplot2::ggplot(d, ggplot2::aes(y = value, fill = fill)) +
      ggplot2::geom_bar() +
      ggplot2::labs(y = column, x = "Samples")
  }
  base +
    ggplot2::labs(title = paste("Distribution of", column), fill = fill_by) +
    theme_tcga()
}

#' Build a lightweight, patient-annotated sample table straight from an
#' availability scan - no assay data is downloaded, only the clinical and
#' subtype lookups (a few hundred KB, seconds to fetch).
#'
#' This is a PREVIEW, not the analysis cohort: it includes every sample the
#' GDC lists for the requested omics, without replicate-aliquot
#' de-duplication or cohort filtering. Counts can differ slightly - usually
#' by a handful of duplicate aliquots - from what build_cohort() produces
#' once you actually fetch and analyse the data.
#'
#' @param avail output of scan_availability()
#' @param clinical output of fetch_clinical(), or NULL
#' @param subtypes output of fetch_subtypes(), or NULL
scan_annotation <- function(avail, clinical = NULL, subtypes = NULL) {
  long <- attr(avail, "long")
  if (is.null(long)) return(NULL)
  long$sample_id <- tcga_sample(long$barcode)
  long <- long[!duplicated(long$sample_id), ]
  ann <- data.frame(
    barcode      = long$barcode,
    sample_id    = long$sample_id,
    patient      = long$patient,
    sample_type  = long$sample_type,
    tissue_class = long$tissue,
    tss          = long$tss,
    plate        = long$plate,
    center       = long$center,
    stringsAsFactors = FALSE
  )
  ann <- join_patient_tables(ann, clin_ = clinical, subtype_ = subtypes)
  flatten_list_columns(ann)
}

#' One row per (grouping variable, level): how many samples and patients
#' carry it. This is the table you read to fill in cohort$group_by,
#' cohort$treatment and cohort$reference - column name, exact level spelling,
#' and the count on each side of the contrast, all in one place.
#'
#' Columns are kept only if they are plausible grouping variables: at least 2
#' levels, not almost-unique (drops patient/case identifier columns), and not
#' continuous-looking (drops things like age or days-to-event, which need
#' bucketing before they are useful as a `group_by`).
#'
#' @param sample_annotation data.frame with a `patient` column, as returned
#'   by scan_annotation() or as.data.frame(colData(cohort_se))
#' @param max_levels drop columns with more levels than this (free-text /
#'   near-identifier columns that slipped past the near-unique filter)
#' @return data.frame: variable, source, level, n_samples, n_patients
summarise_categories <- function(sample_annotation, max_levels = 30) {
  cd <- sample_annotation
  batch_fields <- c("tss", "plate", "center", "portion")
  candidates <- c("tissue_class", "sample_type", batch_fields,
                  grep("^(clin_|subtype_|paper_)", names(cd), value = TRUE))
  candidates <- unique(candidates[candidates %in% names(cd)])

  rows <- lapply(candidates, function(v) {
    is_batch <- v %in% batch_fields
    x <- as.character(cd[[v]])
    x[x %in% c("", "NA", "not reported", "Not Reported", "[Not Available]",
              "[Not Evaluated]", "[Unknown]", "[Not Applicable]")] <- NA
    u <- unique(stats::na.omit(x))
    ## batch identifiers (tss/plate/center) legitimately have many levels and
    ## are alphanumeric codes, not free text or a continuous measurement, so
    ## none of the "looks like noise" filters below apply to them - the whole
    ## point is to see every one, to check it isn't confounded with the
    ## contrast (see covariate_is_confounded() in R/utils.R for the same
    ## check applied automatically once you actually run an analysis)
    if (length(u) < 2) return(NULL)                                  # constant
    if (!is_batch) {
      if (length(u) > max_levels) return(NULL)                       # free text
      if (length(u) > 0.9 * nrow(cd)) return(NULL)                   # identifier column
      xn <- suppressWarnings(as.numeric(u))
      if (!any(is.na(xn)) && length(u) > 10) return(NULL)             # continuous
    }

    source <- if (is_batch) "batch"
              else if (grepl("^clin_", v)) "clinical"
              else if (grepl("^(subtype_|paper_)", v)) "subtype"
              else "barcode"

    do.call(rbind, lapply(u, function(lv) {
      idx <- which(x == lv)
      data.frame(variable = v, source = source, level = lv,
                 n_samples = length(idx),
                 n_patients = length(unique(cd$patient[idx])),
                 stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(NULL)
  out <- out[order(out$source, out$variable, -out$n_samples), ]
  rownames(out) <- NULL
  out
}


#'
#' @param cfg config
#' @param avail output of scan_availability()
#' @param se optional SummarizedExperiment (annotated) for clinical plots -
#'   used by the full pipeline, where the exact analysis cohort is available
#' @param sample_annotation optional plain data.frame with a `patient` column
#'   and clinical/subtype columns - used by the lightweight scan preview
#'   (scan_annotation()). Takes precedence over `se` when both are given.
#' @param clinical_vars columns to profile; auto-detected when NULL
explore_project <- function(cfg, avail = NULL, se = NULL, sample_annotation = NULL,
                            clinical_vars = NULL) {
  log_step("exploring ", cfg$project)
  fig <- cfg$figures_dir; tab <- cfg$tables_dir

  if (!is.null(avail)) {
    write_table(avail[, c("omic", "sample_type", "tissue", "n_samples", "n_patients")],
                file.path(tab, "omics_availability.tsv"))
    save_plot(plot_sample_types(avail, cfg$project),
              file.path(fig, "01_sample_availability.png"), width = 10, height = 6)
    p <- plot_omics_overlap(avail, cfg$project)
    if (!is.null(p)) save_plot(p, file.path(fig, "02_multiomic_overlap.png"), width = 9, height = 5)
  }

  cd <- NULL
  if (!is.null(sample_annotation)) {
    cd <- sample_annotation
  } else if (!is.null(se)) {
    need_pkg("SummarizedExperiment")
    cd <- as.data.frame(SummarizedExperiment::colData(se))
  }
  if (!is.null(cd)) cd <- flatten_list_columns(cd)

  if (!is.null(cd)) {
    write_table(cd, file.path(tab, "sample_annotation.tsv"))

    if (is.null(clinical_vars)) {
      candidates <- c("sample_type", "tissue_class",
                      grep("^(clin_|subtype_|paper_)", names(cd), value = TRUE))
      ## keep variables that are informative: >1 level, not almost-unique ids
      informative <- vapply(candidates, function(c) {
        u <- length(unique(stats::na.omit(as.character(cd[[c]]))))
        u > 1 && u <= max(25, nrow(cd) * 0.5)
      }, logical(1))
      clinical_vars <- utils::head(candidates[informative], 12)
    }
    log_msg("profiling clinical variables: ", paste(clinical_vars, collapse = ", "))
    for (i in seq_along(clinical_vars)) {
      p <- plot_clinical_variable(cd, clinical_vars[i])
      if (!is.null(p)) {
        save_plot(p, file.path(fig, sprintf("03_clinical_%02d_%s.png", i, slug(clinical_vars[i]))),
                  width = 8, height = 5)
      }
    }
    ## cross-tab of every profiled variable against tissue class
    xt <- do.call(rbind, lapply(clinical_vars, function(c) {
      if (!c %in% names(cd)) return(NULL)
      t <- as.data.frame(table(variable = c, level = as.character(cd[[c]]),
                               tissue = as.character(cd$tissue_class)))
      t[t$Freq > 0, ]
    }))
    if (!is.null(xt)) write_table(xt, file.path(tab, "clinical_crosstab.tsv"))
  }
  log_ok("exploration written to ", cfg$results_dir)
  invisible(TRUE)
}
