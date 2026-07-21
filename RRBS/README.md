# RRBS differential methylation analysis

This repository contains the minimal analysis code corresponding to the RRBS methods:

- `01_primary_promoter.R`: MspI fragment filter (40–220 bp), strand-aware promoter filter (−1500/+500 bp), CpG-level testing, promoter/gene aggregation, and the ZNF423 complete-case sensitivity analysis.
- `02_secondary_whole_genome.R`: whole-genome CpG analysis without MspI or promoter filtering, nearest-gene annotation, and gene-level aggregation using CpGs within ±1500 bp of a TSS.

## Input

Both scripts start from the same coverage-filtered `.RData` file containing:

- `meth_50`: methylated-read counts
- `cov_50`: total coverage

Rows must be named as `chromosome:position` (for example, `20:123456`) and columns must be sample identifiers beginning with `DG` or `NG`.

The shared matrix is expected to contain CpGs with at least 50× coverage in at least 50% of samples.

## Required R packages

```r
install.packages("data.table")
BiocManager::install(c("edgeR", "GenomicRanges"))
```

## Configuration

Edit the paths at the top of each script:

```r
RDATA_IN <- "/path/to/RRBS_genome_results_wholegenome.RData"
MSPI_BED <- "/path/to/MspI_fragments_hg38.bed"
REFGENE  <- "/path/to/hg38_TSS.bed"
```

Then run:

```bash
Rscript 01_primary_promoter.R
Rscript 02_secondary_whole_genome.R
```

## Statistical workflow

Missing methylated-count and coverage values are replaced by the rounded locus-wise mean. Beta values are calculated as methylated counts divided by total coverage and clipped to `[0.001, 0.999]`; M-values are `log2(beta / (1 - beta))`.

Methylated counts are analysed in edgeR using total coverage as the sample library size. No TMM normalization is applied. The model uses `~0 + group`, and the `DG - NG` contrast is tested with `glmFit()` and `glmLRT()` after `filterByExpr()` and dispersion estimation.

CpG candidates are reported at:

- FDR < 0.25 and |delta-beta| > 0.05
- exploratory p < 0.001 and |delta-beta| > 0.05

Aggregated gene/promoter candidates are reported at exploratory p < 0.001 and |delta-beta| > 0.05.

## Deliberately excluded

The scripts omit PCA/t-SNE, volcano plots, boxplots, RNA-seq cross-referencing, random-forest modelling, hard-coded candidate corrections, and dissertation-specific text because these are not part of the stated differential-methylation methods.

