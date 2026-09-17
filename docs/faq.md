# FAQ and troubleshooting

## Errors you will probably hit

**`group_by 'X' is not a column of the sample annotation`**
The error prints the first 40 available columns. Run `--only explore` first and read `tables/sample_annotation.tsv` for the full list. Clinical columns are prefixed `clin_`, subtypes `subtype_`.

**`no samples match the requested levels`**
The message lists every level actually observed for your grouping column. Levels are matched literally, so `Stage I` ≠ `stage i` ≠ `STAGE I`. Copy the spelling from the error output.

**`group smaller than min_group_n`**
Genuinely too few samples. Either widen the grouping (merge `Stage I` and `Stage II` into an "early" category by pre-processing, or pick a different variable), or accept that this contrast cannot be run. Lowering `min_group_n` does not create statistical power.

**`GDCdownload` fails or times out**
The GDC API is occasionally flaky. Lower the chunk size (`--chunk 5`) and re-run — downloads resume from disk, so nothing is lost. Corporate proxies and VPNs block it fairly often.

**`GDCprepare` runs out of memory**
Almost always methylation. Options: run methylation separately from the RNA-seq steps (`--only methylation`), give R more memory, or use the Docker image with a raised container memory limit.

**`error in evaluating the argument 'x' in selecting a method for function 'assay'`**
Usually a stale cache entry from an interrupted download. Delete the offending `cache/<PROJECT>/<layer>.rds` and re-fetch.

**Empty or missing subtype data**
`TCGAquery_subtype()` only covers the diseases with published marker papers. Missing subtypes are not an error — the framework logs it and continues, you just cannot use `subtype_*` as a grouping variable for that project.

## Methodological questions

**Should I use `paired: true` for tumour vs normal?**
If you have enough pairs, usually yes. Blocking on patient removes inter-individual variability and typically finds more genes despite the smaller n. Check the pair count first — the cohort builder reports it. Below ~10 pairs the loss of samples usually outweighs the gain.

**Why are my adjacent normals not really normal?**
They are not. Tissue adjacent to a tumour carries field-effect changes and is transcriptionally distinct from healthy tissue from a healthy donor. This is a well-documented limitation of every TCGA tumour-vs-normal analysis, not a bug in the pipeline. If it matters for your question, GTEx is the comparison you want — with the substantial caveat that GTEx and TCGA were processed by different pipelines, so a naive merge produces mostly batch effect.

**Can I compare across tumour types?**
Not directly through the config, which is scoped to one project. `scripts/03_quick_contrast.R` runs the *same* contrast across several projects and writes a `batch_summary.tsv`, which is usually what people actually want. A true pan-cancer joint model needs the project as a covariate and a merged object — a reasonable extension, not currently implemented.

**Why is my p-value histogram not flat?**
A hump near 1 usually means overly conservative dispersion estimates or a mis-specified design. A slope toward 1 often means a confounder. Look at the PCA coloured by candidate confounders (sex, plate, stage) before adding terms blindly.

**Why so many significant genes?**
With n in the hundreds, tiny differences become significant. That is why `lfc_threshold` exists. For tumour vs normal in a large cohort, several thousand significant genes is expected and largely real — rank by effect size, not by p-value, when choosing what to follow up.

## Practical

**Disk usage.** One tumour type with RNA-seq + mutations is roughly 2–5 GB; adding 450k methylation pushes it to 15–20 GB. `cache/` is gitignored — never commit it.

**Can I share results?** Yes. TCGA open-access data (counts, masked MAFs, methylation, de-identified clinical) can be redistributed freely. Controlled-access data — raw BAMs, unmasked germline variants — is not touched by this framework at all, which is why no dbGaP token is needed.

**Can I use my own counts instead of TCGA?** Yes, with a little work: build a `SummarizedExperiment` with an integer assay named `unstranded` and a `colData` containing your grouping column, then call `build_cohort()` and `run_deseq2()` directly. The DESeq2 module has no TCGA-specific assumptions beyond the assay name.

## Should I worry about batch effects?

**Usually yes, for TCGA specifically.** Samples were collected across dozens of institutions (Tissue Source Sites) over years and processed on many sequencing plates — neither is biological, and both are well-documented sources of structure in TCGA data that can masquerade as a real effect, or mask one.

**Nothing is corrected by default.** See [Batch correction](config_reference.md#batch-correction) in the config reference for the full picture. In short:

1. `tss`, `plate` and `center` are parsed from every barcode automatically — no extra download. Check `results/scan/<project>/tables/available_contrasts.tsv` (from `00_explore_projects.R`) to see their distribution, and the console output of that same script for whether they're confounded with tumour/normal before you commit to a contrast.
2. Add one as a covariate: `deseq2: {covariates: [tss]}` (also supported for `mirna` and `methylation`).
3. A confounded covariate is dropped automatically with a warning — this is not something you need to check for by hand.
4. For hidden batch structure with no clean label, `deseq2: {auto_sva: true}` estimates it directly from the data.

**Why isn't this on by default?** Because a batch correction applied without checking can remove real signal — if one collection site happens to have sent mostly tumours and another mostly normals, "correcting for site" and "removing the tumour/normal effect" become close to the same operation. The confounding check exists specifically to catch that case and refuse to proceed with it silently.

**Mutation and copy-number contrasts** (`mafCompare`, the CNV Fisher test) have no covariate mechanism at all — they're simple contingency-table tests. If batch matters there, filter the cohort (`cohort$filters`) rather than trying to adjust for it after the fact.
