## cnv.R ----------------------------------------------------------------------
## Gene-level copy number.
##
## The GDC ships two flavours of gene-level CNV and they need different
## handling, so the module detects which one it was given:
##
##   * ASCAT3 "Gene Level Copy Number"  -> absolute integer copy number per gene
##     (0, 1, 2, 3, ...). Must be interpreted RELATIVE TO SAMPLE PLOIDY: copy
##     number 3 in a near-tetraploid genome is a relative loss, not a gain.
##   * "Gene Level Copy Number Scores"  -> GISTIC-style discrete calls already
##     in {-2,-1,0,1,2}. Use as-is.
##
## Everything downstream works on a discretised matrix with values
##   -2 deep deletion, -1 loss, 0 neutral, +1 gain, +2 amplification.
## -----------------------------------------------------------------------------

#' Pick the assay holding copy number, tolerating GDC naming changes.
pick_cnv_assay <- function(se, requested = NULL) {
  need_pkg("SummarizedExperiment")
  av <- SummarizedExperiment::assayNames(se)
  if (!is.null(requested)) {
    if (!requested %in% av) log_die("assay '", requested, "' not found. Available: ", paste(av, collapse = ", "))
    return(requested)
  }
  for (cand in c("copy_number", "CNV", "score", "segment_mean")) {
    if (cand %in% av) return(cand)
  }
  log_warn("no known CNV assay name; using the first one: ", av[1])
  av[1]
}

#' Is this matrix already GISTIC-style discrete calls?
looks_like_gistic <- function(m) {
  v <- stats::na.omit(as.vector(m[seq_len(min(2000, nrow(m))), , drop = FALSE]))
  length(v) > 0 && all(v %in% -2:2)
}

#' Per-sample ploidy, estimated as the median copy number across all genes.
#' Crude but robust, and it is what makes a "gain" call meaningful in an
#' aneuploid genome.
estimate_ploidy <- function(m) {
  p <- apply(m, 2, stats::median, na.rm = TRUE)
  p[is.na(p) | p <= 0] <- 2
  p
}

#' Discretise absolute copy number into -2/-1/0/1/2 relative to sample ploidy.
#'
#' @param m gene x sample matrix of absolute copy number
#' @param gain_ratio  CN/ploidy at or above which a gene is a gain
#' @param amp_ratio   CN/ploidy at or above which it is an amplification
#' @param loss_ratio  CN/ploidy at or below which it is a loss
discretise_cnv <- function(m, gain_ratio = 1.4, amp_ratio = 2.0,
                           loss_ratio = 0.6, use_ploidy = TRUE) {
  ploidy <- if (use_ploidy) estimate_ploidy(m) else rep(2, ncol(m))
  ratio  <- sweep(m, 2, ploidy, "/")
  out <- matrix(0L, nrow = nrow(m), ncol = ncol(m),
                dimnames = dimnames(m))
  out[ratio >= gain_ratio] <- 1L
  out[ratio >= amp_ratio]  <- 2L
  out[ratio <= loss_ratio] <- -1L
  out[!is.na(m) & m == 0]  <- -2L     # homozygous deletion is absolute
  out[is.na(m)] <- NA_integer_
  attr(out, "ploidy") <- ploidy
  out
}

#' Alteration frequency per gene, per group.
#' @return data.frame: gene, group, n, freq_gain, freq_loss
cnv_frequencies <- function(calls, groups) {
  groups <- as.character(groups)
  do.call(rbind, lapply(unique(groups), function(g) {
    sub <- calls[, groups == g, drop = FALSE]
    n <- rowSums(!is.na(sub))
    data.frame(
      gene      = rownames(calls),
      group     = g,
      n         = n,
      freq_gain = rowSums(sub >= 1, na.rm = TRUE) / pmax(n, 1),
      freq_loss = rowSums(sub <= -1, na.rm = TRUE) / pmax(n, 1),
      stringsAsFactors = FALSE
    )
  }))
}

#' Fisher test per gene: is this gene altered more often in one group?
#' Run separately for gains and losses, BH-corrected across genes.
cnv_group_test <- function(calls, groups, cfg) {
  lv <- c(cfg$cohort$reference, cfg$cohort$treatment)
  groups <- as.character(groups)
  keep <- groups %in% lv
  calls <- calls[, keep, drop = FALSE]; groups <- groups[keep]
  is_trt <- groups == cfg$cohort$treatment

  min_freq <- cfg$cnv$min_freq %||% 0.05
  res <- do.call(rbind, lapply(c("gain", "loss"), function(dir) {
    alt <- if (dir == "gain") calls >= 1 else calls <= -1
    n_trt <- rowSums(alt[, is_trt, drop = FALSE], na.rm = TRUE)
    n_ref <- rowSums(alt[, !is_trt, drop = FALSE], na.rm = TRUE)
    tot_trt <- sum(is_trt); tot_ref <- sum(!is_trt)
    ## test only genes altered often enough somewhere - saves ~50k useless tests
    testable <- (n_trt / tot_trt >= min_freq) | (n_ref / tot_ref >= min_freq)
    if (!any(testable)) return(NULL)
    idx <- which(testable)
    pv <- vapply(idx, function(i) {
      tb <- matrix(c(n_trt[i], tot_trt - n_trt[i], n_ref[i], tot_ref - n_ref[i]), nrow = 2)
      stats::fisher.test(tb)$p.value
    }, numeric(1))
    data.frame(
      gene = rownames(calls)[idx], alteration = dir,
      n_treatment = n_trt[idx], freq_treatment = n_trt[idx] / tot_trt,
      n_reference = n_ref[idx], freq_reference = n_ref[idx] / tot_ref,
      pvalue = pv, stringsAsFactors = FALSE
    )
  }))
  if (is.null(res)) return(NULL)
  res$padj <- stats::p.adjust(res$pvalue, method = "BH")
  res$delta_freq <- res$freq_treatment - res$freq_reference
  res[order(res$padj, -abs(res$delta_freq)), ]
}

#' Genome-wide alteration frequency: gains above the axis, losses below.
#' The familiar "GISTIC-like" view of which chromosome arms are unstable.
plot_cnv_landscape <- function(freq, rowranges, cfg) {
  need_pkg("ggplot2")
  if (is.null(rowranges)) return(NULL)
  pos <- data.frame(
    gene = names(rowranges),
    chr  = as.character(GenomicRanges::seqnames(rowranges)),
    start = GenomicRanges::start(rowranges),
    stringsAsFactors = FALSE
  )
  d <- merge(freq, pos, by = "gene")
  d <- d[d$chr %in% paste0("chr", c(1:22, "X")), ]
  if (nrow(d) == 0) return(NULL)
  d$chr <- factor(d$chr, levels = paste0("chr", c(1:22, "X")))
  d <- d[order(d$chr, d$start), ]

  long <- rbind(
    data.frame(d[, c("gene", "group", "chr", "start")], freq = d$freq_gain, type = "gain"),
    data.frame(d[, c("gene", "group", "chr", "start")], freq = -d$freq_loss, type = "loss")
  )
  ggplot2::ggplot(long, ggplot2::aes(start, freq, colour = type)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::facet_grid(group ~ chr, scales = "free_x", space = "free_x", switch = "x") +
    ggplot2::scale_colour_manual(values = c(gain = "#C0392B", loss = "#2E86C1")) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(abs(round(100 * x)), "%")) +
    ggplot2::labs(title = "Copy-number landscape",
                  subtitle = sprintf("%s | gains above, losses below", cfg$project),
                  x = NULL, y = "Altered samples", colour = NULL) +
    theme_tcga(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.spacing.x = grid::unit(0.05, "lines"),
                   strip.text.x = ggplot2::element_text(size = 6, angle = 90))
}

plot_cnv_top_genes <- function(test, cfg, n = 25) {
  need_pkg("ggplot2")
  if (is.null(test) || nrow(test) == 0) return(NULL)
  d <- utils::head(test, n)
  d$label <- paste0(d$gene, " (", d$alteration, ")")
  d$label <- factor(d$label, levels = rev(d$label))
  long <- rbind(
    data.frame(label = d$label, freq = d$freq_treatment, group = cfg$cohort$treatment),
    data.frame(label = d$label, freq = d$freq_reference, group = cfg$cohort$reference)
  )
  ggplot2::ggplot(long, ggplot2::aes(freq, label, fill = group)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(title = "Most differentially altered genes",
                  subtitle = "Fisher test per gene, BH-corrected",
                  x = "Altered samples", y = NULL, fill = NULL) +
    theme_tcga()
}

#' Does copy number actually drive expression?
#'
#' Spearman correlation per gene between the copy-number call and the
#' variance-stabilised expression, over the samples that have both layers.
#' A strong positive tail is the expected dosage effect; genes that are
#' amplified but NOT over-expressed are often the interesting ones.
correlate_cnv_expression <- function(calls, vsd, cfg, min_samples = 20) {
  need_pkg("SummarizedExperiment")
  if (is.null(vsd)) return(invisible(NULL))
  expr <- SummarizedExperiment::assay(vsd)
  rownames(expr) <- strip_ensembl_version(rownames(expr))
  cn <- calls
  rownames(cn) <- strip_ensembl_version(rownames(cn))

  ## match on patient: CNV and RNA-seq aliquots differ
  p_cn <- tcga_patient(colnames(cn)); p_ex <- tcga_patient(colnames(expr))
  common_p <- intersect(p_cn, p_ex)
  if (length(common_p) < min_samples) {
    log_warn(sprintf("only %d patients have both CNV and expression - correlation skipped",
                     length(common_p)))
    return(invisible(NULL))
  }
  cn <- cn[, match(common_p, p_cn), drop = FALSE]
  expr <- expr[, match(common_p, p_ex), drop = FALSE]

  genes <- intersect(rownames(cn), rownames(expr))
  genes <- genes[!duplicated(genes)]
  log_msg(sprintf("CNV-expression correlation over %d genes x %d patients",
                  length(genes), length(common_p)))

  rho <- vapply(genes, function(g) {
    x <- as.numeric(cn[g, ]); y <- as.numeric(expr[g, ])
    if (length(unique(stats::na.omit(x))) < 2) return(NA_real_)
    suppressWarnings(stats::cor(x, y, method = "spearman", use = "complete.obs"))
  }, numeric(1))

  out <- data.frame(gene = genes, spearman_rho = rho, stringsAsFactors = FALSE)
  out <- out[!is.na(out$spearman_rho), ]
  out <- out[order(-out$spearman_rho), ]
  write_table(out, file.path(cfg$tables_dir, "cnv_expression_correlation.tsv"))

  need_pkg("ggplot2")
  p <- ggplot2::ggplot(out, ggplot2::aes(spearman_rho)) +
    ggplot2::geom_histogram(bins = 60, fill = "#34495E") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40", linetype = 2) +
    ggplot2::labs(title = "Copy-number dosage effect on expression",
                  subtitle = sprintf("Spearman rho per gene, %d patients with both layers",
                                     length(common_p)),
                  x = "rho (copy number vs expression)", y = "genes") +
    theme_tcga()
  save_plot(p, file.path(cfg$figures_dir, "52_cnv_expression_correlation.png"), width = 7, height = 5)
  invisible(out)
}

#' Full copy-number analysis for a cohort.
#'
#' @param cohort output of build_cohort() on the CNV SummarizedExperiment
#' @param vsd optional vst object from the RNA-seq run, for the dosage analysis
analyse_cnv <- function(cohort, cfg, vsd = NULL) {
  need_pkg("SummarizedExperiment")
  cc <- cfg$cnv
  fig <- cfg$figures_dir; tab <- cfg$tables_dir
  log_step("copy-number analysis")

  se <- cohort$se
  assay_name <- pick_cnv_assay(se, cc$assay)
  m <- SummarizedExperiment::assay(se, assay_name)
  storage.mode(m) <- "double"

  ## prefer gene symbols as row names when the object carries them
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  sym_col <- intersect(c("gene_name", "Gene_Symbol", "symbol"), names(rd))
  if (length(sym_col) > 0) {
    sym <- as.character(rd[[sym_col[1]]])
    ok <- !is.na(sym) & sym != "" & !duplicated(sym)
    m <- m[ok, , drop = FALSE]; rownames(m) <- sym[ok]
  }

  if (looks_like_gistic(m)) {
    log_msg("input looks like GISTIC-style scores - using them directly")
    calls <- m; storage.mode(calls) <- "integer"
  } else {
    log_msg("input looks like absolute copy number - discretising relative to sample ploidy")
    calls <- discretise_cnv(m,
                            gain_ratio = cc$gain_ratio %||% 1.4,
                            amp_ratio  = cc$amp_ratio  %||% 2.0,
                            loss_ratio = cc$loss_ratio %||% 0.6,
                            use_ploidy = isTRUE(cc$ploidy_correct %||% TRUE))
    pl <- attr(calls, "ploidy")
    write_table(data.frame(sample = names(pl), estimated_ploidy = as.numeric(pl)),
                file.path(tab, "cnv_estimated_ploidy.tsv"))
  }

  groups <- as.character(SummarizedExperiment::colData(se)$cohort_group)

  freq <- cnv_frequencies(calls, groups)
  write_table(freq, file.path(tab, "cnv_alteration_frequency.tsv"))

  test <- cnv_group_test(calls, groups, cfg)
  if (!is.null(test)) {
    write_table(test, file.path(tab, "cnv_group_comparison.tsv"))
    sig <- sum(test$padj < 0.05, na.rm = TRUE)
    log_ok(sprintf("%d gene(s) differentially altered at FDR < 0.05", sig))
    p <- plot_cnv_top_genes(test, cfg, n = cc$top_genes %||% 25)
    if (!is.null(p)) save_plot(p, file.path(fig, "51_cnv_top_genes.png"), width = 9, height = 8)
  }

  rr <- try(SummarizedExperiment::rowRanges(se), silent = TRUE)
  if (!inherits(rr, "try-error") && length(rr) > 0) {
    if (length(sym_col) > 0) names(rr) <- as.character(rd[[sym_col[1]]])
    p <- plot_cnv_landscape(freq, rr, cfg)
    if (!is.null(p)) save_plot(p, file.path(fig, "50_cnv_landscape.png"), width = 14, height = 6)
  }

  if (isTRUE(cc$correlate_expression) && !is.null(vsd)) {
    correlate_cnv_expression(calls, vsd, cfg)
  }
  invisible(list(calls = calls, frequency = freq, test = test))
}
