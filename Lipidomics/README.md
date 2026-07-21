
# Lipidomics DG versus NG analysis

This repository contains a methods-aligned analysis of DNA-normalised
lipidomics measurements from dysglycaemic (`DG`) and normoglycaemic (`NG`)
participants.

## Workflow

The script:

1. loads the DNA-normalised lipidomics worksheet;
2. restricts the data to subjects present in the predefined DG and NG lists;
3. removes lipid species with more than five missing values;
4. derives lipid classes from lipid names;
5. converts each lipid to its percentage of the total abundance of its class
   within each sample;
6. imputes remaining missing or non-positive values with half the minimum
   positive value in the percent-normalised dataset;
7. applies a log2 transformation;
8. assesses normality separately in DG and NG using the Shapiro-Wilk test;
9. uses Welch's t-test when both groups pass normality at p > 0.05;
10. otherwise uses the Mann-Whitney U test;
11. applies Benjamini-Hochberg correction independently within each lipid
    class.

No global BH correction, nominal-p-value significance rule, fold-change
threshold, FDR 0.1 threshold, PCA, t-SNE, volcano plot, heatmap, or HbA1c
correlation is included because these are not part of the supplied methods.

## Files

- `Lipidomics_DG_vs_NG.R`: complete analysis script.
- `README.md`: repository documentation.

## Required input

Place this workbook in the working directory, or change `INPUT_FILE` at the
top of the script:

```text
Lipidomicsdata.xlsx
```

The default data worksheet is:

```text
Normalized female data with NA
```

The first column must contain unique lipid names. All other relevant columns
must be named with subject identifiers.

## Required R packages

```r
install.packages(c(
  "readxl",
  "data.table"
))
```

## Run

```bash
Rscript Lipidomics_DG_vs_NG.R
```

Results are written to `lipidomics_results/`.

## Main outputs

- `lipid_missingness_summary.tsv`
- `lipid_percent_normalised_imputed.tsv`
- `lipid_log2_final_matrix.tsv`
- `normality_tests.tsv`
- `DG_vs_NG_all_lipids_per_class_BH.tsv`
- `DG_vs_NG_significant_lipids_FDR_0.05.tsv`
- `lipid_class_summary.tsv`
- `lipidomics_DG_vs_NG_results.RData`
- `sessionInfo.txt`

## Reproducibility checks

The supplied methods report 65 matched samples (`DG = 32`, `NG = 33`) and 315
retained lipid species from 666 measured species. The script prints the
observed counts and warns when the retained-lipid count differs from 315.

The significant-results export uses within-class FDR < 0.05 as the conventional
reporting cutoff. The complete results file retains every raw p-value and
within-class BH-adjusted p-value, so another reporting threshold can be applied
without rerunning the statistical tests.
