# Config reference

Every option, its default, and what it actually changes. Anything you omit falls back to `DEFAULT_CONFIG` in `R/config.R`.

## Top level

| Key | Default | Meaning |
|---|---|---|
| `project` | *(required)* | GDC project id, e.g. `TCGA-BRCA` |
| `cache_dir` | `cache` | downloads land in `<cache_dir>/<project>/` |
| `results_dir` | `results` | outputs land in `<results_dir>/<contrast_name>/` |
| `seed` | `1234` | set before anything stochastic |

`contrast_name` is derived automatically as `<project>_<treatment>_vs_<reference>`, so two different contrasts on the same project never overwrite each other.

## `omics`

Which layers to download. All default to `false` except `rnaseq` and `mutation`.

| Key | Data downloaded |
|---|---|
| `rnaseq` | STAR raw counts (`Gene Expression Quantification`) |
| `mutation` | masked somatic MAF (ensemble caller) |
| `methylation` | Illumina 450k beta values — large, opt in deliberately |
| `mirna` | miRNA quantification — downloads, no analysis module yet |
| `cnv` | gene-level copy number (ASCAT3) — downloads, no analysis module yet |

## `cohort`

This block is the one that matters. It defines the comparison.

| Key | Default | Meaning |
|---|---|---|
| `group_by` | `tissue_class` | any column of the sample annotation |
| `treatment` | `tumor` | **numerator** — positive log2FC means higher here |
| `reference` | `normal` | **denominator** — the baseline |
| `paired` | `false` | keep only patients with both groups, block on patient |
| `filters` | none | `column: [allowed, values]`, applied before grouping |
| `min_group_n` | `3` | hard stop below this; n<3 is not interpretable |
| `drop_duplicate_aliquots` | `true` | one aliquot per sample; leave this on |

**Direction matters.** `treatment: tumor, reference: normal` gives positive fold changes for genes up in tumour. Swap them and every sign flips.

**Always filter to tumours for tumour-vs-tumour contrasts.** Without `filters: {tissue_class: [tumor]}`, a subtype contrast can silently pull in adjacent normals that happen to carry a subtype label inherited from their patient.

```yaml
cohort:
  group_by: subtype_BRCA_Subtype_PAM50
  treatment: Basal
  reference: LumA
  min_group_n: 10
  filters:
    tissue_class: [tumor]
    clin_gender: [female]
```

## `deseq2`

| Key | Default | Meaning |
|---|---|---|
| `assay` | `unstranded` | must be raw integer counts |
| `min_count` | `10` | pre-filter threshold |
| `min_samples` | smallest group size | how many samples must clear `min_count` |
| `protein_coding_only` | `true` | drops lncRNA, pseudogenes, etc. |
| `covariates` | none | extra terms added *before* the group term |
| `alpha` | `0.05` | adjusted-p cutoff (BH) |
| `lfc_threshold` | `1` | effect-size cutoff for calling significance |
| `shrink` | `true` | apeglm, falling back to ashr |
| `top_n_heatmap` | `50` | genes in the heatmap |
| `run_enrichment` | `false` | GO:BP over-representation, needs clusterProfiler |

Notes:

- **Pre-filtering is independent filtering**, done before testing and without looking at group labels, so it does not inflate the false discovery rate. It mainly speeds things up and improves the multiple-testing correction.
- **`lfc_threshold` is applied post hoc** here (filtering the results table), not passed to `results(lfcThreshold=)`. That is the common convention but it is *not* a formal test against the threshold. If you need that guarantee, call `DESeq2::results(dds, lfcThreshold = 1)` directly.
- **Covariates cost degrees of freedom.** Add them because you have a reason (a known confounder), not defensively. A covariate confounded with the contrast makes the model unidentifiable; `make_dds()` drops both constant *and* confounded covariates automatically and warns which and why (see [Batch correction](#batch-correction) below).

## Batch correction

**Nothing is corrected automatically.** `covariates` (RNA-seq, miRNA, methylation) and `auto_sva` (RNA-seq only) are both opt-in. This is a deliberate choice, not a missing feature: a batch correction applied blindly can silently absorb real biological signal when the batch happens to align with the contrast, which is common in TCGA (a given collection site often contributed disproportionately more tumours or more normals).

**What's available to correct for:**

| Source | Where it comes from | Cost |
|---|---|---|
| `tss`, `plate`, `center`, `portion` | parsed directly from the barcode (`R/utils.R`) — always present, no extra download | free |
| any `clin_*` column | the clinical table | free (clinical download is small) |
| hidden / unknown structure | `deseq2$auto_sva: true` — estimated from the count matrix with `sva::svaseq()` | needs the `sva` package; adds runtime |

`tss` (Tissue Source Site) and `plate` are the standard TCGA batch proxies — samples sharing a TSS were collected at the same institution, samples sharing a plate were processed together. Both are documented sources of non-biological variation in TCGA specifically. They're computed for every layer automatically; add one with `covariates: [tss]`.

**Confounding is checked automatically**, for every covariate, every time: `covariate_is_confounded()` (`R/utils.R`) builds a covariate × group contingency table and looks for a zero cell — a batch level that only exists on one side of the contrast. If it finds one, that covariate is dropped and a warning names the offending level, because a fully confounded covariate cannot be adjusted for: DESeq2/limma will alias it out, report it as `NA`, or the whole design can become rank-deficient. The fix in that case is `cohort$filters`, not a covariate — restrict the cohort so the batch appears on both sides, or accept that batch and biology cannot be separated for these samples.

You can see this *before* downloading anything: `scripts/00_explore_projects.R` runs the same confounding check against `tissue_class` and prints the result for `tss`/`plate`/`center`.

**`auto_sva`** (RNA-seq only) is for when you don't have — or don't trust — an explicit batch label. It estimates hidden structure directly from the expression matrix with `sva::svaseq()` and adds it as `SV1..SVn`, ahead of the group term. Pass known covariates alongside it (`covariates: [tss]`, `auto_sva: true`) so SVA looks for what they *don't* already explain rather than rediscovering the same thing. `n_sv` fixes the count; leave it `NULL` to auto-estimate. This needs the `sva` package (`install/install_dependencies.R --optional`, or it's in `environment.yml` under the optional block) and is silently skipped with a warning if it isn't installed.

**Mutation and copy-number group comparisons do not support covariates.** `mafCompare()` and the CNV Fisher test (`R/cnv.R`) are simple 2×2 contingency tests with no covariate slot. If you suspect a batch effect there, the only lever is `cohort$filters` — restrict the cohort so the confound doesn't exist in the first place. Extending these to a covariate-adjusted logistic model (`glm(altered ~ group + batch, family = binomial)`) is a reasonable future addition, not currently implemented.

See also: [FAQ — batch effects](faq.md#should-i-worry-about-batch-effects).

## `mutation`

| Key | Default | Meaning |
|---|---|---|
| `top_genes` | `25` | genes shown in the oncoplot |
| `min_mut` | `5` | minimum mutated samples for `mafCompare` |
| `compare_groups` | `true` | run the Fisher test between groups |

Set `compare_groups: false` for tumour-vs-normal contrasts — TCGA MAFs contain tumour-only calls, so the normal arm has no mutations to compare.

## `methylation`

| Key | Default | Meaning |
|---|---|---|
| `platform` | `Illumina Human Methylation 450` | array to query |
| `drop_sex_chr` | `true` | removes chrX/chrY probes |
| `max_na_fraction` | `0.1` | drop probes missing in >10% of samples |
| `p_adjust` | `0.05` | adjusted-p cutoff |
| `delta_beta` | `0.2` | minimum absolute Δβ to call a DMP |
| `top_n_heatmap` | `50` | probes in the heatmap |

`delta_beta: 0.2` is a conventional, fairly strict threshold for tumour vs normal. For subtler contrasts (stage, subtype) `0.1` is more realistic — see `config/luad_stage_i_vs_iv.yml`.

## `cnv`

Gene-level copy number. Enable with `omics: {cnv: true}`.

| Key | Default | Meaning |
|---|---|---|
| `assay` | auto | assay holding copy number; auto-detects `copy_number` / GISTIC scores |
| `ploidy_correct` | `true` | call gains/losses relative to each sample's ploidy |
| `gain_ratio` | `1.4` | CN / ploidy at or above this is a gain |
| `amp_ratio` | `2.0` | CN / ploidy at or above this is an amplification |
| `loss_ratio` | `0.6` | CN / ploidy at or below this is a loss |
| `min_freq` | `0.05` | only test genes altered in ≥5% of one group |
| `top_genes` | `25` | genes in the comparison barplot |
| `correlate_expression` | `true` | Spearman CNV vs expression, needs the RNA-seq step |

`ploidy_correct` matters more than it looks. ASCAT3 gives absolute copy number, so a fixed "CN ≥ 3 = gain" rule mislabels every near-tetraploid tumour — in a genome with ploidy 4, three copies is a relative loss. With correction on, calls are made on the CN/ploidy ratio and only homozygous deletion (CN = 0) stays absolute.

`min_freq` exists because testing all ~20,000 genes when most are altered in two samples wastes the FDR budget on noise. Raise it to `0.1` for a stricter, cleaner table.

The dosage correlation matches CNV and RNA-seq **by patient**, not by barcode, because the two layers are profiled from different aliquots of the same tumour.

## `mirna`

miRNA expression, tested with DESeq2 exactly like mRNA. Enable with `omics: {mirna: true}`.

| Key | Default | Meaning |
|---|---|---|
| `min_count` | `10` | pre-filter threshold on raw read counts |
| `min_samples` | smallest group size | how many samples must clear it |
| `covariates` | none | extra design terms |
| `alpha` | `0.05` | adjusted-p cutoff |
| `lfc_threshold` | `1` | effect-size cutoff |
| `shrink` | `true` | apeglm, falling back to ashr |

There are only ~1,900 miRNAs versus ~20,000 genes, so the multiple-testing burden is far lighter — but dispersion estimates are noisier and a handful of highly abundant miRNAs dominate the library. Treat a miRNA hit as a lead, not a result.
