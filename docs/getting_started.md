# Getting started

A complete first run, from an empty clone to an HTML report. Budget about 30–60 minutes, most of it download time.

## 0. Prerequisites

- Either conda/mamba, or R ≥ 4.3 with Bioconductor ≥ 3.18
- ~20 GB free disk for one tumour type with all three layers
- ~8 GB RAM for RNA-seq + mutations; ~16 GB if you enable 450k methylation
- A working internet connection to `api.gdc.cancer.gov` (no account or token needed — this framework only touches open-access data)

### Option A — conda (recommended)

```bash
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate tcga-omics
```

This brings R, every Bioconductor package **and** the system libraries that normally break an `install.packages()` run. Use `mamba` if you have it; the Bioconductor dependency graph is large and plain conda's solver is slow on it.

Verify:

```bash
Rscript -e 'library(TCGAbiolinks); library(DESeq2); library(maftools); cat("ok\n")'
```

### Option B — an R you already have

```bash
Rscript install/install_dependencies.R          # add --optional for enrichment
```

If a package fails on Linux, it is almost always a missing system library:

```bash
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev \
                     libfontconfig1-dev libharfbuzz-dev libfribidi-dev
```

### Option C — Docker

```bash
docker build -t tcga-omics .
```

## 1. Decide what to analyse

```bash
Rscript scripts/00_explore_projects.R --list
```

Then inventory the tumour type you care about:

```bash
Rscript scripts/00_explore_projects.R --project TCGA-BRCA
```

Open `results/scan/TCGA-BRCA/figures/01_sample_availability.png`. The question you are answering is: **do I have enough normals?** Several TCGA projects have fewer than ten adjacent normals, which makes a tumour-vs-normal contrast underpowered no matter how many tumours there are. BRCA, LUAD, LUSC, KIRC, THCA, PRAD, LIHC, COAD, STAD and HNSC all have reasonable normal counts; OV, LAML, GBM, SKCM and UCS have very few or none.

If the normals are thin, pivot to a tumour-vs-tumour contrast (subtype, stage, mutation status) instead — see `config/brca_pam50_basal_vs_luma.yml`.

## 2. Download once

```bash
cp config/brca_tumor_vs_normal.yml config/my_analysis.yml
$EDITOR config/my_analysis.yml
Rscript scripts/01_fetch_data.R --config config/my_analysis.yml
```

This writes `cache/TCGA-BRCA/rnaseq.rds`, `mutation.rds`, `clinical.rds` and, if the project has one, `subtypes.rds`. It is safe to interrupt and re-run — `GDCdownload` resumes from what is already on disk, and prepared objects are only rebuilt when missing.

If the download keeps timing out, lower the chunk size:

```bash
Rscript scripts/01_fetch_data.R --config config/my_analysis.yml --chunk 5
```

## 3. Find your grouping variable

This is the point of `00_explore_projects.R` — it fetches the clinical and molecular-subtype tables (a few hundred KB, not the actual omics data) and writes `available_contrasts.tsv`: one row per (variable, level), with sample and patient counts on each side. This is what you actually read to fill in `cohort.group_by` / `treatment` / `reference`:

```bash
Rscript scripts/00_explore_projects.R --project TCGA-BRCA
cat results/scan/TCGA-BRCA/tables/available_contrasts.tsv
```

```
variable                     source     level        n_samples   n_patients
tissue_class                 barcode    tumor        1123        1099
tissue_class                 barcode    normal       113         113
clin_ajcc_pathologic_stage   clinical   Stage II     626         620
clin_ajcc_pathologic_stage   clinical   Stage I      183         181
clin_ajcc_pathologic_stage   clinical   Stage III    236         233
clin_ajcc_pathologic_stage   clinical   Stage IV     20          20
subtype_BRCA_Subtype_PAM50   subtype    LumA         566         566
subtype_BRCA_Subtype_PAM50   subtype    Basal        190         190
subtype_BRCA_Subtype_PAM50   subtype    LumB         203         203
subtype_BRCA_Subtype_PAM50   subtype    Her2         82          82
subtype_BRCA_Subtype_PAM50   subtype    Normal       40          40
```

The console also prints a compact one-row-per-variable overview so you can see at a glance which columns exist and where they come from (`barcode` = derived from the TCGA barcode, always available; `clinical` = `GDCquery_clinic()`; `subtype` = the disease's curated molecular-subtype table, when one exists). Copy the exact `variable` name and `level` spelling into your YAML — matching is literal, so `Stage I` ≠ `stage i`.

This table is a **preview**: it comes straight from the GDC file listing and the clinical download, before any replicate-aliquot de-duplication or cohort filtering, so counts can differ by a handful of samples from what `build_cohort()` reports once you actually run the analysis. Skip it with `--no-clinical` if you only want the old availability-only scan.

For the exact, post-filtering counts of your final cohort, `sample_annotation.tsv` from a full run has every column too:

| Prefix | Source | Example |
|---|---|---|
| *(none)* | derived from the barcode | `tissue_class`, `sample_type` |
| `clin_` | `GDCquery_clinic()` | `clin_gender`, `clin_ajcc_pathologic_stage` |
| `subtype_` | `TCGAquery_subtype()` | `subtype_BRCA_Subtype_PAM50` |
| `paper_` | shipped inside the GDC object | `paper_BRCA_Subtype_PAM50` |

Copy the column name and the **exact** spelling of the two levels you want into your config. Level names are matched literally: `Stage I` will not match `stage i`. If you get them wrong, the cohort builder fails loudly and prints every level it actually observed — which is the fastest way to find the right spelling.

## 4. Run everything

```bash
Rscript scripts/02_run_analysis.R --config config/my_analysis.yml
```

Open `results/<contrast>/report.html`.

## 5. Read the output in the right order

1. **`cohort_summary.tsv`** — are the group sizes what you expected? If not, a filter is wrong.
2. **`10_pca.png`** — do the groups separate at all? If tumour and normal overlap completely, something is mislabelled.
3. **`14_pvalue_histogram.png`** — flat with a spike at zero is healthy. A hump anywhere else means the model is misspecified: add a covariate, or check for a confounded batch.
4. **`12_volcano.png`** and `deseq2_results_significant.tsv` — only now look at the genes.

Skipping to step 4 is how people publish batch effects.

## Typical runtimes

| Step | TCGA-BRCA, ~1200 samples |
|---|---|
| Availability scan | under a minute |
| RNA-seq download + prepare | 15–30 min |
| Methylation 450k download + prepare | 45–90 min, heavy RAM |
| DESeq2 (tumour vs normal) | 5–15 min |
| Mutation analysis | 2–5 min |
| limma DMPs | 5–10 min |
| Copy number download + analysis | 10–20 min |
| miRNA download + DESeq2 | 5–10 min |

Every step after the first download reads from the cache and is dramatically faster.
