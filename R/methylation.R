## methylation.R --------------------------------------------------------------
## Differentially methylated probes (DMPs) from Illumina 450k/EPIC beta values.
##
## Statistics are run on M-values (logit of beta) because beta values are
## heteroscedastic and bounded, which violates the linear-model assumptions.
## Effect sizes are reported as delta-beta because that is the interpretable
## scale.
## -----------------------------------------------------------------------------

beta_to_m <- function(beta, eps = 1e-3) {
  beta[beta < eps] <- eps
  beta[beta > 1 - eps] <- 1 - eps
  log2(beta / (1 - beta))
}

#' Filter probes: missingness, sex chromosomes, zero variance.
filter_probes <- function(se, cfg) {
  need_pkg("SummarizedExperiment")
  mm <- cfg$methylation
  b <- SummarizedExperiment::assay(se)
  n0 <- nrow(b)

  na_frac <- rowMeans(is.na(b))
  keep <- na_frac <= mm$max_na_fraction
  log_msg(sprintf("probes with <=%.0f%% NA: %d/%d", 100 * mm$max_na_fraction, sum(keep), n0))

  if (isTRUE(mm$drop_sex_chr)) {
    chr <- NULL
    rr <- try(SummarizedExperiment::rowRanges(se), silent = TRUE)
    if (!inherits(rr, "try-error") && length(rr) == nrow(b)) {
      chr <- as.character(GenomicRanges::seqnames(rr))
    } else if ("chr" %in% names(SummarizedExperiment::rowData(se))) {
      chr <- as.character(SummarizedExperiment::rowData(se)$chr)
    }
    if (!is.null(chr)) {
      sexy <- chr %in% c("chrX", "chrY", "X", "Y")
      log_msg(sprintf("dropping %d sex-chromosome probes", sum(sexy & keep)))
      keep <- keep & !sexy
    }
  }
  se <- se[keep, ]
  b <- SummarizedExperiment::assay(se)
  ## impute the few remaining NAs with the probe mean so limma keeps the probe
  idx <- which(is.na(b), arr.ind = TRUE)
  if (nrow(idx) > 0) b[idx] <- rowMeans(b, na.rm = TRUE)[idx[, 1]]
  SummarizedExperiment::assay(se) <- b
  vr <- matrixStats_rowVars(b)
  se <- se[vr > 0 & !is.na(vr), ]
  log_ok(sprintf("%d probes retained", nrow(se)))
  se
}

matrixStats_rowVars <- function(x) {
  if (has_pkg("matrixStats")) return(matrixStats::rowVars(x))
  rowSums((x - rowMeans(x))^2) / (ncol(x) - 1)
}

#' Differential methylation between the two cohort groups with limma.
#'
#' @return data.frame of probe-level statistics
run_dmp <- function(cohort, cfg) {
  need_pkg("limma"); need_pkg("SummarizedExperiment")
  mm <- cfg$methylation
  fig <- cfg$figures_dir; tab <- cfg$tables_dir

  log_step("differential methylation: ", cfg$cohort$treatment, " vs ", cfg$cohort$reference)
  se <- filter_probes(cohort$se, cfg)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  beta <- SummarizedExperiment::assay(se)
  mval <- beta_to_m(beta)

  grp <- droplevels(cd$cohort_group)
  cd$grp <- grp   # convenience alias used in the design formula below

  covars <- unlist(mm$covariates %||% character(0))
  covars <- covars[covars %in% names(cd)]
  for (cv in covars) cd[[cv]] <- factor(as.character(cd[[cv]]))
  if (isTRUE(cohort$paired)) covars <- unique(c("patient", covars))
  covars <- clean_covariates(covars, cd)   # drops constant / confounded covariates

  design_str <- paste("~", paste(c(covars, "grp"), collapse = " + "))
  log_msg("design: ", design_str)
  design <- stats::model.matrix(stats::as.formula(design_str), data = cd)
  coef_i <- ncol(design)   # the group term is always last
  fit <- limma::eBayes(limma::lmFit(mval, design))
  tt <- limma::topTable(fit, coef = coef_i, number = Inf, sort.by = "P")

  ## delta-beta: mean(treatment) - mean(reference), the reportable effect size
  is_trt <- grp == cfg$cohort$treatment
  db <- rowMeans(beta[, is_trt, drop = FALSE]) - rowMeans(beta[, !is_trt, drop = FALSE])

  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  out <- data.frame(
    probe      = rownames(tt),
    gene       = rd[rownames(tt), intersect(c("gene_HGNC", "Gene_Symbol", "gene"), names(rd))[1]] %||% NA,
    delta_beta = db[rownames(tt)],
    logFC_M    = tt$logFC,
    pvalue     = tt$P.Value,
    padj       = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
  out$significant <- out$padj < mm$p_adjust & abs(out$delta_beta) >= mm$delta_beta
  out$direction <- ifelse(!out$significant, "ns",
                          ifelse(out$delta_beta > 0, "hyper", "hypo"))
  out <- out[order(out$padj), ]

  write_table(out, file.path(tab, "dmp_results_all.tsv"))
  write_table(out[out$significant, ], file.path(tab, "dmp_results_significant.tsv"))
  log_ok(sprintf("%d significant DMPs (%d hyper / %d hypo)", sum(out$significant),
                 sum(out$direction == "hyper"), sum(out$direction == "hypo")))

  save_plot(plot_dmp_volcano(out, cfg), file.path(fig, "30_dmp_volcano.png"), width = 8, height = 6)
  p <- plot_dmp_heatmap(beta, out, grp, cfg)
  if (!is.null(p)) save_plot(p, file.path(fig, "31_dmp_heatmap.png"), width = 9, height = 9)

  invisible(out)
}

plot_dmp_volcano <- function(tab, cfg) {
  need_pkg("ggplot2")
  d <- tab[!is.na(tab$padj), ]
  ggplot2::ggplot(d, ggplot2::aes(delta_beta, -log10(pmax(padj, .Machine$double.xmin)),
                                  colour = direction)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * cfg$methylation$delta_beta, linetype = 2, colour = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(cfg$methylation$p_adjust), linetype = 2, colour = "grey40") +
    ggplot2::scale_colour_manual(values = c(hyper = "#C0392B", hypo = "#2E86C1", ns = "grey75")) +
    ggplot2::labs(title = "Differentially methylated probes",
                  subtitle = sprintf("%s: %s vs %s", cfg$project, cfg$cohort$treatment, cfg$cohort$reference),
                  x = expression(Delta * beta), y = "-log10 adjusted p", colour = NULL) +
    theme_tcga()
}

plot_dmp_heatmap <- function(beta, tab, grp, cfg) {
  need_pkg("ggplot2")
  sig <- tab[tab$significant, ]
  if (nrow(sig) < 2) { log_warn("not enough significant DMPs for a heatmap"); return(NULL) }
  top <- utils::head(sig, cfg$methylation$top_n_heatmap)
  m <- beta[top$probe, , drop = FALSE]
  ord <- order(grp)
  m <- m[, ord, drop = FALSE]; g <- as.character(grp)[ord]
  hc <- stats::hclust(stats::dist(m))
  df <- expand.grid(sample = colnames(m), probe = rownames(m), stringsAsFactors = FALSE)
  df$beta <- m[cbind(match(df$probe, rownames(m)), match(df$sample, colnames(m)))]
  df$probe <- factor(df$probe, levels = rownames(m)[hc$order])
  df$sample <- factor(df$sample, levels = colnames(m))
  df$group <- g[match(df$sample, colnames(m))]

  ggplot2::ggplot(df, ggplot2::aes(sample, probe, fill = beta)) +
    ggplot2::geom_raster() +
    ggplot2::facet_grid(~ group, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradientn(colours = c("#2E86C1", "#F7F7F7", "#C0392B"), limits = c(0, 1)) +
    ggplot2::labs(title = sprintf("Top %d differentially methylated probes", nrow(top)),
                  x = NULL, y = NULL, fill = expression(beta)) +
    theme_tcga(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

#' Integrate methylation with expression: genes that are both hypermethylated
#' and down-regulated (or the reverse) are the classic epigenetic candidates.
integrate_meth_expression <- function(dmp, de, cfg) {
  if (is.null(dmp) || is.null(de)) return(invisible(NULL))
  if (!"gene" %in% names(dmp) || all(is.na(dmp$gene))) {
    log_warn("no gene annotation on the methylation probes - integration skipped")
    return(invisible(NULL))
  }
  d <- dmp[dmp$significant & !is.na(dmp$gene), ]
  d$gene <- sub(";.*$", "", as.character(d$gene))   # keep the first mapped gene
  e <- de[de$significant & !is.na(de$gene_name), ]
  mg <- merge(d, e, by.x = "gene", by.y = "gene_name", suffixes = c("_meth", "_expr"))
  if (nrow(mg) == 0) { log_msg("no overlap between DMPs and DE genes"); return(invisible(NULL)) }

  mg$relationship <- ifelse(mg$delta_beta > 0 & mg$log2FoldChange < 0, "hyper_down",
                     ifelse(mg$delta_beta < 0 & mg$log2FoldChange > 0, "hypo_up", "concordant_other"))
  write_table(mg, file.path(cfg$tables_dir, "methylation_expression_integration.tsv"))
  log_ok(sprintf("%d gene-probe pairs, %d classic epigenetic-silencing candidates",
                 nrow(mg), sum(mg$relationship == "hyper_down")))

  need_pkg("ggplot2")
  p <- ggplot2::ggplot(mg, ggplot2::aes(delta_beta, log2FoldChange, colour = relationship)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
    ggplot2::labs(title = "Methylation vs expression",
                  subtitle = "each point is a significant probe mapped to a significant gene",
                  x = expression(Delta * beta), y = "log2 fold change (expression)", colour = NULL) +
    theme_tcga()
  save_plot(p, file.path(cfg$figures_dir, "32_meth_vs_expression.png"), width = 8, height = 6)
  invisible(mg)
}
