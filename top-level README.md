# GM_2026_Obesity.project

Multi-omics profiling of visceral adipose tissue (VAT) adipocytes, comparing dysglycaemic (DG) and normoglycaemic (NG) obese women matched by age and BMI via propensity-score matching. The aim is to identify molecular signatures distinguishing DG from NG using integrated transcriptomics (RNA-seq), DNA methylation (RRBS), and lipidomics, with MOFA2 multi-omics factor analysis and random forest classification as integrative layers.

## Repository structure

```
RNAseq/           Differential expression, pathway enrichment, snoRNA diagnostics
RRBS/             DNA methylation (locus- and gene-level), MspI/TSS reference files
  reference/        MspI fragment + TSS annotation files, and the script to regenerate MspI fragments
  wholegenome_sensitivity/   Genome-wide sensitivity branch (no TSS filter), HPC-run
Lipidomics/       Lipid abundance differential testing, clinical correlations
MOFA/             Multi-omics factor analysis, integration, combined biomarker panel
Integration/      Analyses spanning two or more omics layers (methylation x expression)
utils/            Shared plotting colours/theme used across all pipelines
```

Each folder's scripts are numbered and meant to be run **in order**: every script saves the R objects the next one needs into that folder's `11.06.run/rds/`, so you don't need to keep everything in one long session. Each folder has its own `README.md` with the run order and any folder-specific caveats — read those before running anything, especially RRBS and MOFA, which have several flagged issues (see below).

## Recommended run order across folders

```
RNAseq  →  RRBS  →  Lipidomics  →  MOFA  →  Integration
```

RNA-seq and RRBS can technically run independently of each other, but MOFA needs the outputs of all three single-omics pipelines, and Integration needs RNA-seq + RRBS (or RNA-seq + Lipidomics, depending on the script).

## Data availability

Raw sequencing data, clinical metadata, and other patient-derived files are **not included** in this repository (see `.gitignore`) — both for size and because they're not appropriate to version-control publicly. To rerun any pipeline from scratch you'll need, alongside the code:

- Raw counts, gene annotation, and sample metadata for RNA-seq
- Bismark coverage files (or the pre-built locus matrices — see `RRBS/reference/README.md` for the one unresolved dependency here) for RRBS
- The lipidomics Excel workbook (`Lipidomicsdata.xlsx`)
- `clinical_metadata.xlsx`, used by RNA-seq, MOFA, and (in a cut-down form) Lipidomics

Reference files that are small and non-patient-derived (e.g. `RRBS/reference/hg38_TSS.bed`) are committed directly.

## Canonical analysis thresholds

These are the thresholds the analysis should converge on — several scripts in this repo are flagged where the source notebooks didn't consistently apply them (see each folder's README for specifics):

| Layer | Threshold |
|---|---|
| RNA-seq DEGs | FDR < 0.05 & \|log2FC\| > log2(1.5) |
| RRBS, primary | FDR < 0.25 & \|Δβ\| > 0.05 |
| RRBS, exploratory | p < 0.001 & \|Δβ\| > 0.01 |
| Lipidomics | FDR < 0.1, per-class BH correction, no fold-change filter |
| MOFA CV AUCs (confirmed) | all 15 factors = 0.780, Factors 1+2+3 = 0.746 |

## Known open items

- **RRBS coverage matrices**: `RRBS/01_data_loading.R` depends on a pre-built `RRBS_genome_results.RData` whose construction script (reading raw Bismark `.cov.gz` files) isn't in this repo — see `RRBS/reference/README.md`.
- **RRBS threshold inconsistencies**: locus-level vs gene-level primary/exploratory delta-beta cutoffs don't fully match the canonical thresholds above in the source code — see `RRBS/README.md`.
- **Lipidomics significance flag**: an earlier, superseded significance definition (nominal p<0.05 + fold-change filter) is kept alongside the canonical one (FDR<0.1, no FC filter) for comparison only — see `Lipidomics/README.md`.
- **MOFA AUC values**: most of the MOFA pipeline's hardcoded comparison numbers (0.750/0.769) don't match the confirmed canonical values (0.780/0.746) — see `MOFA/README.md` for the full explanation, and recommend regenerating these numbers programmatically before final write-up.
- **Cohort count (DG=31 vs 32)**: a discrepancy in the DG sample count appears across sections of the dissertation and is visible in the code too (`MOFA` expects DG=31 from the clinical sheet; the hardcoded `dg_subjects` list used in `Lipidomics` and `MOFA/11_combined_omics_panel.R` has 32 entries) — needs reconciling before submission.

## Software

R packages (confirmed versions): edgeR 4.8.2, limma 3.66.0, Rtsne 0.17, clusterProfiler 4.18.4, org.Hs.eg.db 3.22.0, ggplot2 4.0.3, randomForest, pROC, MOFA2, DiagrammeR, data.table, patchwork, enrichplot.

Plotting conventions: DG = `#D94F3D`, NG = `#4C72B0` (see `utils/plotting_helpers.R`) — though note the colour-inconsistency flag in `Lipidomics/README.md`, where some figures use a different palette.
