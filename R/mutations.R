## mutations.R ----------------------------------------------------------------
## Somatic mutation layer, built on maftools.
##
## Important: TCGA somatic MAFs are tumour-only calls (matched normal used as a
## filter, not reported as a sample). A "tumour vs normal" mutation comparison
## is therefore meaningless here - mutation contrasts must be between two groups
## of TUMOURS (subtype A vs B, stage I vs IV, mutant vs wild-type, ...).
## -----------------------------------------------------------------------------

#' Convert a (possibly cohort-restricted) MAF data.frame into a maftools object.
make_maf <- function(maf_df, clinical = NULL) {
  need_pkg("maftools")
  if (is.null(maf_df) || nrow(maf_df) == 0) { log_warn("empty MAF"); return(NULL) }

  clin <- NULL
  if (!is.null(clinical)) {
    clin <- clinical
    clin$Tumor_Sample_Barcode <- clin$Tumor_Sample_Barcode %||% clin$patient
    clin <- clin[!duplicated(clin$Tumor_Sample_Barcode), , drop = FALSE]
  }
  ## maftools keys on Tumor_Sample_Barcode; align it to the patient id so the
  ## clinical join and the cohort_group labels line up.
  maf_df$Tumor_Sample_Barcode <- maf_df$patient %||% maf_df$Tumor_Sample_Barcode

  if (!is.null(maf_df$cohort_group)) {
    grp <- unique(maf_df[, c("Tumor_Sample_Barcode", "cohort_group")])
    clin <- if (is.null(clin)) grp else merge(clin, grp, by = "Tumor_Sample_Barcode", all.y = TRUE)
  }
  m <- try(maftools::read.maf(maf = maf_df, clinicalData = clin, verbose = FALSE), silent = TRUE)
  if (inherits(m, "try-error")) { log_warn("read.maf failed: ", as.character(m)); return(NULL) }
  m
}

#' Tumour mutational burden per sample (non-synonymous variants per Mb).
#' The 38 Mb denominator is the usual approximation of the exome target size.
compute_tmb <- function(maf_df, exome_size_mb = 38) {
  nonsyn <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
              "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
              "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")
  d <- maf_df[maf_df$Variant_Classification %in% nonsyn, , drop = FALSE]
  tb <- as.data.frame(table(patient = d$patient), stringsAsFactors = FALSE)
  names(tb)[2] <- "n_nonsyn"
  tb$tmb <- tb$n_nonsyn / exome_size_mb
  if (!is.null(maf_df$cohort_group)) {
    map <- unique(maf_df[, c("patient", "cohort_group")])
    tb <- merge(tb, map, by = "patient", all.x = TRUE)
  }
  tb[order(-tb$tmb), ]
}

plot_tmb <- function(tmb, cfg) {
  need_pkg("ggplot2")
  if (!"cohort_group" %in% names(tmb)) return(NULL)
  d <- tmb[!is.na(tmb$cohort_group), ]
  pv <- try(stats::wilcox.test(tmb ~ cohort_group, data = d)$p.value, silent = TRUE)
  sub <- if (inherits(pv, "try-error")) NULL else sprintf("Wilcoxon rank-sum p = %.3g", pv)
  ggplot2::ggplot(d, ggplot2::aes(cohort_group, tmb, fill = cohort_group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Tumour mutational burden", subtitle = sub,
                  x = NULL, y = "non-synonymous mutations / Mb") +
    theme_tcga() + ggplot2::theme(legend.position = "none")
}

#' Full mutation analysis for a cohort.
#'
#' @param maf_df MAF data.frame, ideally already passed through cohort_maf()
#' @param cfg config
#' @param compare when TRUE and cohort_group is present, run mafCompare
analyse_mutations <- function(maf_df, cfg, compare = TRUE) {
  need_pkg("maftools")
  fig <- cfg$figures_dir; tab <- cfg$tables_dir
  log_step("mutation analysis")

  maf <- make_maf(maf_df)
  if (is.null(maf)) return(invisible(NULL))

  gs <- maftools::getGeneSummary(maf)
  write_table(gs, file.path(tab, "mutation_gene_summary.tsv"))
  write_table(maftools::getSampleSummary(maf), file.path(tab, "mutation_sample_summary.tsv"))

  tmb <- compute_tmb(maf_df)
  write_table(tmb, file.path(tab, "tumor_mutational_burden.tsv"))
  p <- plot_tmb(tmb, cfg)
  if (!is.null(p)) save_plot(p, file.path(fig, "21_tmb_by_group.png"), width = 6, height = 5)

  ## base figures (maftools draws with base graphics -> wrap in png device)
  png_plot <- function(path, expr, width = 1800, height = 1400, res = 160) {
    ensure_dir(dirname(path))
    grDevices::png(path, width = width, height = height, res = res)
    ok <- try(force(expr), silent = TRUE)
    grDevices::dev.off()
    if (inherits(ok, "try-error")) log_warn("plot failed: ", basename(path), " - ", as.character(ok))
    else log_ok("wrote ", path)
  }

  clin_features <- if ("cohort_group" %in% names(maf@clinical.data)) "cohort_group" else NULL
  png_plot(file.path(fig, "20_oncoplot.png"),
           maftools::oncoplot(maf = maf, top = cfg$mutation$top_genes,
                              clinicalFeatures = clin_features,
                              sortByAnnotation = !is.null(clin_features)))
  png_plot(file.path(fig, "22_maf_summary.png"), maftools::plotmafSummary(maf = maf))

  ## differential mutation between the two cohort groups ----------------------
  if (compare && "cohort_group" %in% names(maf_df)) {
    lv <- c(cfg$cohort$reference, cfg$cohort$treatment)
    lv <- lv[lv %in% unique(maf_df$cohort_group)]
    if (length(lv) == 2) {
      m1 <- make_maf(maf_df[maf_df$cohort_group == lv[2], , drop = FALSE])
      m2 <- make_maf(maf_df[maf_df$cohort_group == lv[1], , drop = FALSE])
      if (!is.null(m1) && !is.null(m2)) {
        cmp <- try(maftools::mafCompare(m1 = m1, m2 = m2, m1Name = lv[2], m2Name = lv[1],
                                        minMut = cfg$mutation$min_mut), silent = TRUE)
        if (!inherits(cmp, "try-error")) {
          write_table(cmp$results, file.path(tab, "mutation_group_comparison.tsv"))
          png_plot(file.path(fig, "23_forest_plot.png"),
                   maftools::forestPlot(mafCompareRes = cmp, pVal = 0.05))
          png_plot(file.path(fig, "24_cooncoplot.png"),
                   maftools::coOncoplot(m1 = m1, m2 = m2, m1Name = lv[2], m2Name = lv[1],
                                        genes = utils::head(cmp$results$Hugo_Symbol, 20)))
          log_ok(sprintf("%d gene(s) differentially mutated at p < 0.05",
                         sum(cmp$results$pval < 0.05, na.rm = TRUE)))
        } else {
          log_warn("mafCompare failed: ", as.character(cmp))
        }
      }
    } else {
      log_warn("mutation comparison skipped: both groups must contain tumour samples")
    }
  }
  invisible(maf)
}
