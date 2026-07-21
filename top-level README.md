# GM_2026_Obesity.project

Multi-omics profiling of visceral adipose tissue (VAT) adipocytes, comparing dysglycaemic (DG) and normoglycaemic (NG) obese women matched by age and BMI via propensity-score matching. The aim is to identify molecular signatures distinguishing DG from NG using integrated transcriptomics (RNA-seq), DNA methylation (RRBS), and lipidomics, with MOFA2 multi-omics factor analysis and random forest classification as integrative layers.

## Repository structure

```
.
├── README.md
├── scripts/
│   ├── RNAseq_DG_vs_NG.R
│   ├── RRBS_TSS.R
│   ├── RRBS_WholeGenome.R
│   ├── Lipidomics_DG_vs_NG.R
│   ├── MOFA2_multiomics.R
│   └── RandomForest_classification.R
│
├── data/
│   ├── reference/
│   └── results/
│       ├── RNAseq/
│       ├── RRBS/
│       ├── Lipidomics/
│       └── MOFA/
│
└── figures/
```

---

## Overview of analyses

### RNA sequencing

Differential gene expression analysis performed using **edgeR**.

Main steps include:

- quality-controlled sample selection
- low-expression filtering
- TMM normalisation
- dispersion estimation
- negative binomial GLM
- likelihood ratio testing
- Benjamini–Hochberg correction
- GO and KEGG pathway enrichment
- snoRNA over-representation analysis

---

### RRBS methylation

Two complementary methylation analyses were performed.

**1. Promoter (TSS) analysis**

- MspI fragment filtering
- promoter annotation (-1500/+500 bp)
- M-value transformation
- edgeR differential methylation analysis
- gene-level aggregation

**2. Whole-genome analysis**

- genome-wide CpG analysis
- nearest-gene annotation
- gene-level aggregation
- differential methylation testing

---

### Lipidomics

DNA-normalised lipid abundances were analysed following:

- missing-value filtering
- lipid-class assignment
- within-class percent normalisation
- half-minimum imputation
- log2 transformation
- Shapiro-Wilk normality testing
- Welch's t-test or Mann-Whitney U test
- Benjamini-Hochberg correction within lipid classes

---

### Multi-omics integration

Multi-omics integration was performed using **MOFA2**.

The workflow includes:

- harmonisation of transcriptomic, methylation and lipidomic datasets
- feature filtering
- sequencing-depth correction
- training of a 15-factor MOFA2 model
- variance explained analysis
- factor association with clinical variables
- functional enrichment of dominant factors

---

### Machine learning

Random Forest classifiers were constructed to compare predictive performance of

- transcriptomics
- methylation
- lipidomics
- MOFA latent factors
- combined multi-omics features

Performance was assessed using

- stratified 5-fold cross-validation
- 10 repeats
- ROC analysis
- Area Under the Curve (AUC)

---

## Software

Analyses were performed in **R** using packages including:

- edgeR
- MOFA2
- clusterProfiler
- randomForest
- pROC
- vegan
- GenomicRanges
- org.Hs.eg.db
- data.table
- dplyr

Package versions are recorded within each analysis script using `sessionInfo()`.

---

## Data

Raw sequencing and metabolomics datasets are **not included** because they contain participant-derived data.

The repository includes:

- analysis scripts
- processed statistical result tables
- supplementary outputs required to reproduce figures and downstream analyses

Reference genome annotation files (e.g. hg38) should be obtained from their original sources.

---

## Reproducibility

Each analysis can be run independently using the scripts in the `scripts/` directory.

Typical workflow:

1. RNA-seq analysis
2. RRBS promoter analysis
3. RRBS whole-genome analysis
4. Lipidomics analysis
5. MOFA2 integration
6. Random Forest classification

Outputs from earlier analyses are used as inputs for the downstream multi-omics analyses.

