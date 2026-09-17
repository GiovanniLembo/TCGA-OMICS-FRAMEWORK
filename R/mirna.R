## mirna.R --------------------------------------------------------------------
## miRNA expression.
##
## GDCprepare returns miRNA data as a wide data.frame, not a
## SummarizedExperiment: one row per miRNA, and three columns per sample
## (read_count_*, reads_per_million_*, cross-mapped_*). We reshape it into an SE
## so that the exact same cohort machinery applies to this layer too.
##
## read_count is a raw integer count, so DESeq2 is appropriate. Note that miRNA
## libraries have far fewer features than mRNA (~1900 vs ~20000), which makes
## the multiple-testing burden much lighter and the dispersion estimates
## noisier - hence the gentler default filter.
## -----------------------------------------------------------------------------

#' Reshape the GDC miRNA data.frame into a SummarizedExperiment.
mirna_to_se <- function(df) {
  need_pkg("SummarizedExperiment")
  df <- as.data.frame(df)
  id_col <- intersect(c("miRNA_ID", "miRNA_id"), names(df))
  if (length(id_col) == 0) log_die("unexpected miRNA table: no miRNA_ID column")

  count_cols <- grep("^read_count_", names(df), value = TRUE)
  if (length(count_cols) == 0) log_die("unexpected miRNA table: no read_count_ columns")

  counts <- as.matrix(df[, count_cols, drop = FALSE])
  storage.mode(counts) <- "integer"
  rownames(counts) <- df[[id_col[1]]]
  colnames(counts) <- sub("^read_count_", "", count_cols)

  rpm_cols <- grep("^reads_per_million", names(df), value = TRUE)
  assays <- list(read_count = counts)
  if (length(rpm_cols) == length(count_cols)) {
    rpm <- as.matrix(df[, rpm_cols, drop = FALSE])
    dimnames(rpm) <- dimnames(counts)
    assays$rpm <- rpm
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = assays,
    colData = S4Vectors::DataFrame(row.names = colnames(counts))
  )
  log_ok(sprintf("miRNA matrix: %d miRNAs x %d samples", nrow(se), ncol(se)))
  se
}

#' Differential miRNA expression with DESeq2.
#'
#' Reuses the RNA-seq plotting functions, so the figures are directly
#' comparable with the mRNA ones.
run_mirna_de <- function(cohort, cfg) {
  need_pkg("DESeq2"); need_pkg("SummarizedExperiment")
  mm  <- cfg$mirna
  fig <- cfg$figures_dir; tab <- cfg$tables_dir
  log_step("differential miRNA expression: ", cfg$cohort$treatment, " vs ", cfg$cohort$reference)

  se <- cohort$se
  counts <- SummarizedExperiment::assay(se, "read_count")
  storage.mode(counts) <- "integer"
  cd <- as.data.frame(SummarizedExperiment::colData(se))

  min_n <- mm$min_samples %||% min(table(cd$cohort_group))
  keep <- rowSums(counts >= (mm$min_count %||% 10)) >= min_n
  log_msg(sprintf("filter: %d/%d miRNAs retained", sum(keep), length(keep)))
  counts <- counts[keep, , drop = FALSE]
  if (nrow(counts) < 10) log_die("too few miRNAs survive filtering - check the counts")

  covars <- unlist(mm$covariates %||% character(0))
  covars <- covars[covars %in% names(cd)]
  for (cv in covars) cd[[cv]] <- factor(as.character(cd[[cv]]))
  if (isTRUE(cohort$paired)) covars <- unique(c("patient", covars))
  covars <- clean_covariates(covars, cd)
  design_str <- paste("~", paste(c(covars, "cohort_group"), collapse = " + "))
  log_msg("design: ", design_str)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = cd[, unique(c("cohort_group", covars, "patient")), drop = FALSE],
    design    = stats::as.formula(design_str)
  )
  dds$cohort_group <- stats::relevel(droplevels(dds$cohort_group), ref = cfg$cohort$reference)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  coef_name <- utils::tail(DESeq2::resultsNames(dds), 1)
  res <- DESeq2::results(dds, name = coef_name, alpha = mm$alpha %||% 0.05)
  if (isTRUE(mm$shrink %||% TRUE)) {
    sh <- try(DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE), silent = TRUE)
    if (inherits(sh, "try-error")) sh <- try(DESeq2::lfcShrink(dds, coef = coef_name, type = "ashr", quiet = TRUE), silent = TRUE)
    if (!inherits(sh, "try-error")) { res$log2FoldChange <- sh$log2FoldChange; res$lfcSE <- sh$lfcSE }
  }

  ## same column names as the mRNA table so the shared plots just work
  out <- data.frame(
    gene_id = rownames(res), gene_id_versioned = rownames(res),
    gene_name = rownames(res), gene_type = "miRNA",
    baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
    pvalue = res$pvalue, padj = res$padj, stringsAsFactors = FALSE
  )
  out$significant <- !is.na(out$padj) & out$padj < (mm$alpha %||% 0.05) &
    abs(out$log2FoldChange) >= (mm$lfc_threshold %||% 1)
  out$direction <- ifelse(!out$significant, "ns", ifelse(out$log2FoldChange > 0, "up", "down"))
  out <- out[order(out$padj, -abs(out$log2FoldChange)), ]

  write_table(out, file.path(tab, "mirna_results_all.tsv"))
  write_table(out[out$significant, ], file.path(tab, "mirna_results_significant.tsv"))
  log_ok(sprintf("%d significant miRNAs (%d up / %d down)", sum(out$significant),
                 sum(out$direction == "up"), sum(out$direction == "down")))

  cfg_mi <- cfg; cfg_mi$deseq2$alpha <- mm$alpha %||% 0.05
  cfg_mi$deseq2$lfc_threshold <- mm$lfc_threshold %||% 1
  save_plot(plot_volcano(out, cfg_mi), file.path(fig, "40_mirna_volcano.png"), width = 8, height = 7)
  save_plot(plot_ma(out, cfg_mi), file.path(fig, "41_mirna_ma.png"), width = 8, height = 6)

  vsd <- try(DESeq2::vst(dds, blind = FALSE, nsub = min(1000, nrow(dds))), silent = TRUE)
  if (inherits(vsd, "try-error")) vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  save_plot(plot_pca(vsd, cfg_mi, ntop = min(500, nrow(dds))),
            file.path(fig, "42_mirna_pca.png"), width = 7, height = 6)
  p <- plot_top_boxplots(dds, out, cfg_mi)
  if (!is.null(p)) save_plot(p, file.path(fig, "43_mirna_top_boxplots.png"), width = 10, height = 8)

  invisible(list(dds = dds, table = out, vsd = vsd))
}
