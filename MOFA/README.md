# MOFA2 and random-forest multi-omics analysis

This repository contains cleaned, methods-aligned scripts for integrating
transcriptomics, promoter-focused RRBS methylation, and lipidomics data with
MOFA2, followed by random-forest classification of DG versus NG participants.

## Files

- `01_MOFA2_multiomics.R`
- `02_random_forest_classification.R`
- `README.md`

## MOFA2 workflow

`01_MOFA2_multiomics.R`:

1. loads the processed transcriptomic, methylation, and lipidomic matrices;
2. removes `DG`/`NG` prefixes from sample names;
3. retains the three-way sample intersection;
4. selects the 5,000 most variable transcriptomic and methylation features;
5. removes features with more than 50% missing values;
6. retains the processed lipidomics features;
7. regresses sample-wise total signal from RNA and methylation features;
8. scales each omics view to unit variance;
9. trains a 15-factor MOFA2 model with slow convergence and a fixed seed;
10. calculates factor-wise variance explained;
11. identifies the dominant factor for each omics layer;
12. tests DG versus NG factor-score differences with two-sample Welch tests;
13. tests factor-clinical associations using Spearman correlation and BH
    correction;
14. exports top feature weights and feature annotations;
15. performs GO Biological Process and KEGG enrichment for top Factor 1 genes.

The script warns when the run does not reproduce the reported 64 common
samples, 5,000 RNA features, 5,000 methylation features, or 314 lipid species.

## Random-forest workflow

`02_random_forest_classification.R` uses:

- 500 trees;
- the default `mtry`;
- stratified 5-fold cross-validation;
- 10 repeats;
- 50 total train/test splits;
- pooled held-out probabilities for ROC and AUC calculation.

Eight models are evaluated:

1. all 15 MOFA factors;
2. the three layer-dominant factors together;
3. transcriptomic-dominant factor alone;
4. lipidomic-dominant factor alone;
5. methylation-dominant factor alone;
6. transcriptomic features;
7. lipidomic features;
8. methylation features.

The same folds are generated deterministically for every feature set, enabling
direct performance comparisons.

## Required inputs

The default filenames are:

```text
RRBS_genome_results_promoter.RData
RNAseq_DG_vs_NG_results.RData
lipidomics_DG_vs_NG_results.RData
clinical_metadata.xlsx
hg38_refGene.txt.gz
```

Expected objects:

- RRBS RData: `mvalue_mat`
- RNA-seq RData: filtered and TMM-normalised edgeR object `y`
- lipidomics RData: `log2_matrix`

Update the paths at the top of `01_MOFA2_multiomics.R` when needed.

## Required R packages

```r
BiocManager::install(c(
  "MOFA2",
  "edgeR",
  "GenomicRanges",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler"
))

install.packages(c(
  "data.table",
  "readxl",
  "dplyr",
  "randomForest",
  "pROC"
))
```

MOFA2 uses basilisk-managed Python dependencies through `run_mofa()`.

## Run

```bash
Rscript 01_MOFA2_multiomics.R
Rscript 02_random_forest_classification.R
```

## Main MOFA outputs

- `MOFA_model.hdf5`
- `MOFA_model.rds`
- `variance_explained.tsv`
- `layer_dominant_factors.tsv`
- `factor_scores.tsv`
- `factor_group_tests.tsv`
- `factor_clinical_correlations.tsv`
- `top20_weights_per_factor_view.tsv`
- `Factor1_GO_BP_FDR_0.05.tsv`
- `Factor1_KEGG_FDR_0.05.tsv`
- `MOFA_RF_inputs.RData`

## Main random-forest outputs

- `cross_validated_predictions.tsv`
- `pooled_ROC_coordinates.tsv`
- `model_performance_AUC.tsv`
- `random_forest_results.RData`


