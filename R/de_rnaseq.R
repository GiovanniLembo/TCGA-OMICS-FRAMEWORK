## de_rnaseq.R ----------------------------------------------------------------
## Differential expression with DESeq2 on GDC STAR raw counts.
##
## Notes that matter for TCGA specifically:
##  * use the "unstranded" assay - it holds raw integer counts. TPM/FPKM assays
##    are normalised and must never be fed to DESeq2.
##  * Ensembl ids carry a version suffix (ENSG00000141510.16) - strip it before
##    any join with external annotation.
##  * TCGA cohorts are large and heterogeneous; a patient block (paired design)
##    or at least a covariate for plate/gender is usually worth the extra df.
## -----------------------------------------------------------------------------

strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)

#' Estimate hidden batch effects with SVA and add them as covariates SV1..SVn.
#'
#' Complements, rather than replaces, explicit covariates like tss/plate/
#' center: pass those as `covars` so SVA looks for structure they do NOT
#' already explain, instead of rediscovering the same thing. Use this when
#' you suspect batch effects but don't have (or don't trust) an explicit
#' label for them - e.g. samples processed over several years with no
#' recorded plate id.
#'
#' @param counts raw integer count matrix (post gene-filtering)
#' @param cd colData; must contain cohort_group and every column in `covars`
#' @param covars known covariates already in the design
#' @param n_sv number of surrogate variables; NULL auto-estimates with
#'   sva::num.sv() (Buja-Eyuboglu permutation method)
#' @return list(cd, sv_names), or NULL if sva is unavailable, fails, or
#'   estimates zero surrogate variables
estimate_surrogate_variables <- function(counts, cd, covars = character(0), n_sv = NULL) {
  if (!has_pkg("sva")) {
    log_warn("auto_sva requested but the 'sva' package is not installed - skipped. ",
            "Install with: BiocManager::install('sva')")
    return(NULL)
  }
  need_pkg("DESeq2")
  log_step("SVA: estimating hidden surrogate variables")

  sf <- DESeq2::estimateSizeFactorsForMatrix(counts)
  normcounts <- sweep(counts, 2, sf, "/")
  rv <- rowSums((normcounts - rowMeans(normcounts))^2) / (ncol(normcounts) - 1)
  normcounts <- normcounts[rv > 0 & !is.na(rv), , drop = FALSE]

  mod_str  <- paste("~", paste(c(covars, "cohort_group"), collapse = " + "))
  mod0_str <- if (length(covars) > 0) paste("~", paste(covars, collapse = " + ")) else "~ 1"
  mod  <- stats::model.matrix(stats::as.formula(mod_str), data = cd)
  mod0 <- stats::model.matrix(stats::as.formula(mod0_str), data = cd)

  if (is.null(n_sv)) {
    n_sv <- tryCatch(sva::num.sv(normcounts, mod, method = "be"), error = function(e) NA_integer_)
    if (is.na(n_sv)) { log_warn("sva::num.sv() failed to estimate a count - skipped"); return(NULL) }
  }
  if (n_sv < 1) { log_msg("SVA: 0 surrogate variables estimated - nothing added"); return(NULL) }
  n_sv <- min(n_sv, ncol(mod) - 1)

  sv <- try(sva::svaseq(as.matrix(normcounts), mod, mod0, n.sv = n_sv)$sv, silent = TRUE)
  if (inherits(sv, "try-error") || is.null(sv)) { log_warn("sva::svaseq() failed: ", as.character(sv)); return(NULL) }

  sv <- as.data.frame(sv)
  names(sv) <- paste0("SV", seq_len(ncol(sv)))
  cd2 <- cbind(cd, sv)
  log_ok(sprintf("SVA: added %d surrogate variable(s) to the design", ncol(sv)))
  list(cd = cd2, sv_names = names(sv))
}

#' Build a DESeqDataSet from a cohort SummarizedExperiment.
make_dds <- function(cohort, cfg) {
  need_pkg("DESeq2"); need_pkg("SummarizedExperiment")
  se <- cohort$se
  dd <- cfg$deseq2

  assays_available <- SummarizedExperiment::assayNames(se)
  if (!dd$assay %in% assays_available) {
    log_die(sprintf("assay '%s' not found. Available: %s", dd$assay,
                    paste(assays_available, collapse = ", ")))
  }
  counts <- SummarizedExperiment::assay(se, dd$assay)
  storage.mode(counts) <- "integer"
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  rd <- as.data.frame(SummarizedExperiment::rowData(se))

  ## gene filtering ------------------------------------------------------------
  if (isTRUE(dd$protein_coding_only) && "gene_type" %in% names(rd)) {
    keep <- rd$gene_type == "protein_coding"
    log_msg(sprintf("keeping %d protein-coding genes of %d", sum(keep), length(keep)))
    counts <- counts[keep, , drop = FALSE]; rd <- rd[keep, , drop = FALSE]
  }
  min_n <- dd$min_samples %||% min(table(cd$cohort_group))
  expressed <- rowSums(counts >= dd$min_count) >= min_n
  log_msg(sprintf("independent pre-filter: %d/%d genes with >=%d counts in >=%d samples",
                  sum(expressed), nrow(counts), dd$min_count, min_n))
  counts <- counts[expressed, , drop = FALSE]; rd <- rd[expressed, , drop = FALSE]

  ## design --------------------------------------------------------------------
  covars <- unlist(dd$covariates %||% character(0))
  covars <- covars[covars %in% names(cd)]
  for (cv in covars) cd[[cv]] <- factor(as.character(cd[[cv]]))
  if (isTRUE(cohort$paired)) covars <- unique(c("patient", covars))

  ## drop covariates that are constant or confounded with the contrast -
  ## a confounded covariate isn't just weaker, it can make the whole design
  ## matrix rank-deficient or silently absorb the group effect
  covars <- clean_covariates(covars, cd)

  if (isTRUE(dd$auto_sva)) {
    sva_out <- estimate_surrogate_variables(counts, cd, covars = setdiff(covars, "patient"),
                                            n_sv = dd$n_sv)
    if (!is.null(sva_out)) { cd <- sva_out$cd; covars <- c(covars, sva_out$sv_names) }
  }

  design_str <- paste("~", paste(c(covars, "cohort_group"), collapse = " + "))
  log_msg("design: ", design_str)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = cd[, unique(c("cohort_group", covars, "patient", "sample_type")), drop = FALSE],
    design    = stats::as.formula(design_str)
  )
  SummarizedExperiment::rowData(dds) <- S4Vectors::DataFrame(rd)
  dds$cohort_group <- stats::relevel(droplevels(dds$cohort_group), ref = cfg$cohort$reference)
  dds
}

#' Run the DESeq2 pipeline and write results + figures.
#'
#' @return list(dds, res, vsd, results_table)
run_deseq2 <- function(cohort, cfg) {
  need_pkg("DESeq2")
  dd  <- cfg$deseq2
  fig <- cfg$figures_dir; tab <- cfg$tables_dir

  log_step("DESeq2: ", cfg$cohort$treatment, " vs ", cfg$cohort$reference)
  dds <- make_dds(cohort, cfg)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  coef_name <- paste0("cohort_group_", cfg$cohort$treatment, "_vs_", cfg$cohort$reference)
  if (!coef_name %in% DESeq2::resultsNames(dds)) {
    ## non-syntactic level names get mangled by model.matrix; fall back on position
    coef_name <- utils::tail(DESeq2::resultsNames(dds), 1)
    log_warn("using coefficient '", coef_name, "'")
  }
  res <- DESeq2::results(dds, name = coef_name, alpha = dd$alpha)

  if (isTRUE(dd$shrink)) {
    shrunk <- try(DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE), silent = TRUE)
    if (inherits(shrunk, "try-error")) {
      shrunk <- try(DESeq2::lfcShrink(dds, coef = coef_name, type = "ashr", quiet = TRUE), silent = TRUE)
    }
    if (!inherits(shrunk, "try-error")) {
      log_ok("applied LFC shrinkage")
      res$log2FoldChange <- shrunk$log2FoldChange
      res$lfcSE <- shrunk$lfcSE
    } else {
      log_warn("shrinkage unavailable - reporting MLE fold changes")
    }
  }

  ## results table -------------------------------------------------------------
  rd <- as.data.frame(SummarizedExperiment::rowData(dds))
  out <- data.frame(
    gene_id     = strip_ensembl_version(rownames(res)),
    gene_id_versioned = rownames(res),
    gene_name   = rd$gene_name %||% NA_character_,
    gene_type   = rd$gene_type %||% NA_character_,
    baseMean    = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE       = res$lfcSE,
    pvalue      = res$pvalue,
    padj        = res$padj,
    stringsAsFactors = FALSE
  )
  out$significant <- !is.na(out$padj) & out$padj < dd$alpha &
    abs(out$log2FoldChange) >= dd$lfc_threshold
  out$direction <- ifelse(!out$significant, "ns",
                          ifelse(out$log2FoldChange > 0, "up", "down"))
  out <- out[order(out$padj, -abs(out$log2FoldChange)), ]

  write_table(out, file.path(tab, "deseq2_results_all.tsv"))
  write_table(out[out$significant, ], file.path(tab, "deseq2_results_significant.tsv"))
  log_ok(sprintf("%d significant genes (padj < %.3g, |LFC| >= %.2g): %d up / %d down",
                 sum(out$significant), dd$alpha, dd$lfc_threshold,
                 sum(out$direction == "up"), sum(out$direction == "down")))

  ## QC + result figures -------------------------------------------------------
  vsd <- try(DESeq2::vst(dds, blind = FALSE), silent = TRUE)
  if (inherits(vsd, "try-error")) vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)

  save_plot(plot_pca(vsd, cfg), file.path(fig, "10_pca.png"), width = 7, height = 6)
  save_plot(plot_sample_distance(vsd), file.path(fig, "11_sample_distance.png"), width = 8, height = 7)
  save_plot(plot_volcano(out, cfg), file.path(fig, "12_volcano.png"), width = 8, height = 7)
  save_plot(plot_ma(out, cfg), file.path(fig, "13_ma_plot.png"), width = 8, height = 6)
  save_plot(plot_pvalue_hist(out), file.path(fig, "14_pvalue_histogram.png"), width = 7, height = 5)
  p <- plot_top_heatmap(vsd, out, cfg); if (!is.null(p)) save_plot(p, file.path(fig, "15_top_genes_heatmap.png"), width = 9, height = 10)
  p <- plot_top_boxplots(dds, out, cfg); if (!is.null(p)) save_plot(p, file.path(fig, "16_top_genes_boxplots.png"), width = 10, height = 8)

  if (isTRUE(dd$run_enrichment)) run_enrichment(out, cfg)

  invisible(list(dds = dds, res = res, vsd = vsd, table = out))
}

## --- figures -----------------------------------------------------------------

plot_pca <- function(vsd, cfg, ntop = 2000) {
  need_pkg("ggplot2")
  d <- DESeq2::plotPCA(vsd, intgroup = "cohort_group", ntop = ntop, returnData = TRUE)
  pv <- round(100 * attr(d, "percentVar"))
  ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = cohort_group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::stat_ellipse(level = 0.9, linewidth = 0.4, na.rm = TRUE) +
    ggplot2::labs(title = "PCA on variance-stabilised counts",
                  subtitle = sprintf("%s | top %d variable genes", cfg$project, ntop),
                  x = sprintf("PC1 (%d%%)", pv[1]), y = sprintf("PC2 (%d%%)", pv[2]),
                  colour = NULL) +
    theme_tcga()
}

plot_sample_distance <- function(vsd, max_n = 80) {
  need_pkg("ggplot2")
  m <- SummarizedExperiment::assay(vsd)
  if (ncol(m) > max_n) m <- m[, sample(ncol(m), max_n)]
  d <- as.matrix(stats::dist(t(m)))
  grp <- as.character(SummarizedExperiment::colData(vsd)[colnames(d), "cohort_group"])
  ord <- order(grp)
  d <- d[ord, ord]
  df <- expand.grid(x = colnames(d), y = rownames(d), stringsAsFactors = FALSE)
  df$dist <- as.vector(d)
  df$x <- factor(df$x, levels = colnames(d)); df$y <- factor(df$y, levels = rownames(d))
  ggplot2::ggplot(df, ggplot2::aes(x, y, fill = dist)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(direction = -1) +
    ggplot2::labs(title = "Sample-to-sample Euclidean distance",
                  subtitle = "samples ordered by group", x = NULL, y = NULL, fill = "distance") +
    theme_tcga() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
}

plot_volcano <- function(tab, cfg, label_n = 20) {
  need_pkg("ggplot2")
  d <- tab[!is.na(tab$padj), ]
  d$neglog10p <- -log10(pmax(d$padj, .Machine$double.xmin))
  top <- utils::head(d[d$significant, ], label_n)
  p <- ggplot2::ggplot(d, ggplot2::aes(log2FoldChange, neglog10p, colour = direction)) +
    ggplot2::geom_point(size = 1, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * cfg$deseq2$lfc_threshold, linetype = 2, colour = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(cfg$deseq2$alpha), linetype = 2, colour = "grey40") +
    ggplot2::scale_colour_manual(values = c(up = "#C0392B", down = "#2E86C1", ns = "grey75")) +
    ggplot2::labs(title = "Volcano plot",
                  subtitle = sprintf("%s: %s vs %s", cfg$project, cfg$cohort$treatment, cfg$cohort$reference),
                  x = "log2 fold change", y = "-log10 adjusted p", colour = NULL) +
    theme_tcga()
  if (has_pkg("ggrepel") && nrow(top) > 0) {
    p <- p + ggrepel::geom_text_repel(data = top,
                                      ggplot2::aes(label = gene_name), size = 3,
                                      max.overlaps = 30, show.legend = FALSE)
  }
  p
}

plot_ma <- function(tab, cfg) {
  need_pkg("ggplot2")
  d <- tab[!is.na(tab$padj) & tab$baseMean > 0, ]
  ggplot2::ggplot(d, ggplot2::aes(baseMean, log2FoldChange, colour = direction)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.55) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(up = "#C0392B", down = "#2E86C1", ns = "grey75")) +
    ggplot2::labs(title = "MA plot", x = "mean of normalised counts",
                  y = "log2 fold change", colour = NULL) +
    theme_tcga()
}

plot_pvalue_hist <- function(tab) {
  need_pkg("ggplot2")
  d <- tab[!is.na(tab$pvalue), ]
  ggplot2::ggplot(d, ggplot2::aes(pvalue)) +
    ggplot2::geom_histogram(binwidth = 0.02, fill = "#34495E") +
    ggplot2::labs(title = "Raw p-value distribution",
                  subtitle = "a flat tail with a spike near 0 is what you want; anything else means the model is misspecified",
                  x = "p-value", y = "genes") +
    theme_tcga()
}

plot_top_heatmap <- function(vsd, tab, cfg) {
  need_pkg("ggplot2")
  sig <- tab[tab$significant, ]
  if (nrow(sig) < 2) { log_warn("not enough significant genes for a heatmap"); return(NULL) }
  top <- utils::head(sig, cfg$deseq2$top_n_heatmap)
  m <- SummarizedExperiment::assay(vsd)[top$gene_id_versioned, , drop = FALSE]
  m <- t(scale(t(m)))                       # z-score per gene
  m[is.na(m)] <- 0
  grp <- as.character(SummarizedExperiment::colData(vsd)$cohort_group)
  ord <- order(grp)
  m <- m[, ord, drop = FALSE]; grp <- grp[ord]
  hc <- stats::hclust(stats::dist(m))
  labels <- ifelse(is.na(top$gene_name) | top$gene_name == "", top$gene_id, top$gene_name)

  df <- expand.grid(sample = colnames(m), gene = rownames(m), stringsAsFactors = FALSE)
  df$z <- m[cbind(match(df$gene, rownames(m)), match(df$sample, colnames(m)))]
  df$gene <- factor(df$gene, levels = rownames(m)[hc$order],
                    labels = labels[match(rownames(m)[hc$order], top$gene_id_versioned)])
  df$sample <- factor(df$sample, levels = colnames(m))
  df$group <- grp[match(df$sample, colnames(m))]

  ggplot2::ggplot(df, ggplot2::aes(sample, gene, fill = z)) +
    ggplot2::geom_raster() +
    ggplot2::facet_grid(~ group, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(low = "#2E86C1", mid = "white", high = "#C0392B", limits = c(-3, 3),
                                  oob = scales::squish) +
    ggplot2::labs(title = sprintf("Top %d differentially expressed genes", nrow(top)),
                  subtitle = "row z-scores of variance-stabilised counts",
                  x = NULL, y = NULL, fill = "z") +
    theme_tcga(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

plot_top_boxplots <- function(dds, tab, cfg, n = 9) {
  need_pkg("ggplot2")
  sig <- utils::head(tab[tab$significant, ], n)
  if (nrow(sig) == 0) return(NULL)
  nc <- DESeq2::counts(dds, normalized = TRUE)[sig$gene_id_versioned, , drop = FALSE]
  grp <- as.character(SummarizedExperiment::colData(dds)$cohort_group)
  df <- do.call(rbind, lapply(seq_len(nrow(sig)), function(i) {
    data.frame(gene = sig$gene_name[i] %||% sig$gene_id[i],
               value = log2(nc[i, ] + 1), group = grp, stringsAsFactors = FALSE)
  }))
  ggplot2::ggplot(df, ggplot2::aes(group, value, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
    ggplot2::facet_wrap(~ gene, scales = "free_y") +
    ggplot2::labs(title = "Top differentially expressed genes",
                  x = NULL, y = "log2(normalised count + 1)") +
    theme_tcga() + ggplot2::theme(legend.position = "none")
}

#' Optional over-representation analysis (GO BP) on the significant genes.
#' Silently skipped when clusterProfiler / org.Hs.eg.db are not installed.
run_enrichment <- function(tab, cfg) {
  if (!has_pkg("clusterProfiler") || !has_pkg("org.Hs.eg.db")) {
    log_warn("clusterProfiler / org.Hs.eg.db missing - enrichment skipped")
    return(invisible(NULL))
  }
  log_step("GO over-representation analysis")
  universe <- unique(stats::na.omit(tab$gene_id))
  for (dir in c("up", "down")) {
    genes <- unique(tab$gene_id[tab$direction == dir])
    if (length(genes) < 10) next
    ego <- try(clusterProfiler::enrichGO(gene = genes, universe = universe,
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                         keyType = "ENSEMBL", ont = "BP",
                                         pAdjustMethod = "BH", qvalueCutoff = 0.05,
                                         readable = TRUE), silent = TRUE)
    if (inherits(ego, "try-error") || is.null(ego) || nrow(as.data.frame(ego)) == 0) next
    write_table(as.data.frame(ego), file.path(cfg$tables_dir, sprintf("go_bp_%s.tsv", dir)))
    p <- try(clusterProfiler::dotplot(ego, showCategory = 20) +
               ggplot2::ggtitle(sprintf("GO:BP enriched in %s-regulated genes", dir)), silent = TRUE)
    if (!inherits(p, "try-error")) {
      save_plot(p, file.path(cfg$figures_dir, sprintf("17_go_bp_%s.png", dir)), width = 9, height = 8)
    }
  }
  invisible(NULL)
}
