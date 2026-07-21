# RNA-seq DG versus NG analysis

This repository contains a methods-aligned RNA-seq analysis for comparing
dysglycaemic (`DG`) and normoglycaemic (`NG`) participants using edgeR.

## Workflow

The script:

1. loads raw gene counts, gene annotations, and the sequencing-ID mapping;
2. renames samples to subject identifiers and excludes samples that failed
   upstream quality control;
3. creates an edgeR `DGEList`;
4. retains genes with CPM >= 1 in at least 32 samples;
5. applies TMM normalisation with `calcNormFactors()`;
6. tests treatment-associated transcriptomic variation among DG participants
   using PERMANOVA on log2-CPM expression;
7. fits a no-intercept edgeR model (`~0 + group`);
8. tests the `DG - NG` contrast using `glmFit()` and `glmLRT()`;
9. applies Benjamini-Hochberg correction;
10. performs GO Biological Process and KEGG over-representation analyses using
    FDR-significant DEGs and the complete filtered gene set as the universe;
11. tests snoRNA over-representation using Fisher's exact test.

## Files

- `RNAseq_DG_vs_NG.R`: complete analysis script.
- `README.md`: repository documentation.

## Required input files

Place these files in the working directory or update the paths at the top of
the script:

- `raw_count.tsv`
- `human_gene_names.txt`
- `new_metatable copy.txt`
- `clinical_metadata.xlsx`

The sample metadata file is expected to contain `Index`, `Index.1`,
`New.Grouping`, and `Name`. The clinical workbook is expected to use the
two-header-row structure present in the original analysis.

## Required R packages

```r
BiocManager::install(c(
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db"
))

install.packages(c(
  "readxl",
  "dplyr",
  "stringr",
  "vegan",
  "data.table"
))
```

The methods report edgeR 4.8.2 and clusterProfiler 4.18.4. Save the generated
`sessionInfo.txt` with the analysis outputs to document the versions actually
used.

## Run

```bash
Rscript RNAseq_DG_vs_NG.R
```

Results are written to `RNAseq_results/`.

## Main outputs

- `filtered_gene_universe.tsv`
- `DG_treatment_PERMANOVA.txt`
- `DG_vs_NG_all_genes.tsv`
- `DG_vs_NG_significant_DEGs_FDR_0.05.tsv`
- `GO_BP_enrichment_all.tsv`
- `GO_BP_enrichment_FDR_0.05.tsv`
- `KEGG_enrichment_all.tsv`
- `KEGG_enrichment_FDR_0.05.tsv`
- `snoRNA_overrepresentation_Fisher.txt`
- `RNAseq_DG_vs_NG_results.RData`
- `sessionInfo.txt`

## Deliberately excluded from the cleaned script

The original notebook also contained PCA, t-SNE, volcano plots, manually
selected pathway tables, gene-specific boxplots, direction-stratified
enrichment, random-forest classification, and repeated plotting/export chunks.
Those analyses were removed because they are not part of the supplied methods
section.

The DEG definition follows the supplied methods exactly: FDR < 0.05, without
an additional fold-change threshold.
