# ============================================================================
# 01_data_loading.R
# TSS-filtered RRBS pipeline (this is the PRIMARY analysis used as MOFA
# methylation input — see wholegenome_sensitivity/ for the genome-wide
# sensitivity branch, which skips the TSS filter).
#
# Loads the pre-built locus-level coverage/methylation matrices, applies the
# MspI fragment-size filter, restricts to TSS windows (-1500/+500bp) with
# nearest-gene annotation, imputes missing values, and builds beta/M-value
# matrices.
#
# Inputs:  RRBS_genome_results.RData, MspI_fragments_hg38.bed, hg38_TSS.bed
# Outputs: 11.06.run/rds/{meth_50,cov_50,meth_clean,cov_clean,beta_full,
#          mvalue_full,locus_anno,group,dg_idx,ng_idx}.rds
# ============================================================================

library(data.table)
library(dplyr)

OUT_DIR <- "11.06.run"
dir.create(file.path(OUT_DIR, "rds"), recursive = TRUE, showWarnings = FALSE)

# ── 1. Load pre-built locus matrices ────────────────────────────────────────
load("RRBS_genome_results.RData")   # provides meth_50, cov_50
cat("meth_50:", dim(meth_50), "\n")
cat("cov_50: ", dim(cov_50),  "\n")

# ── 2. MspI fragment size filter (40-220bp) ─────────────────────────────────
# MspI_fragments_hg38.bed itself isn't committed to the repo (~112MB, over
# GitHub's 100MB limit) — regenerate it with reference/build_mspi_fragments.py
# before running this, or point this path at wherever you keep it locally.
# FLAG (not fixed, unverified): the actual file has a header row
# ("chr\tfrag_start\tfrag_end\tfrag_id"). fread()'s default header="auto"
# should detect and skip this automatically even with col.names supplied —
# but I don't have an R runtime here to confirm that, so worth checking
# nrow(msp) against (file line count - 1) the first time this runs.
msp <- fread("reference/MspI_fragments_hg38.bed", col.names = c("chr", "start", "end", "frag_id"))
msp[, chr       := sub("chr", "", chr)]
msp[, frag_size := end - start]
msp_filt <- msp[frag_size >= 40 & frag_size <= 220]

locus_dt <- data.table(
  locus = rownames(meth_50),
  chr   = sub(":.*", "", rownames(meth_50)),
  start = as.integer(sub(".*:", "", rownames(meth_50)))
)
locus_dt[, end := start]
setkey(locus_dt,  chr, start, end)
setkey(msp_filt,  chr, start, end)

keep_msp <- unique(foverlaps(locus_dt, msp_filt, nomatch = 0L)$locus)
meth_50  <- meth_50[rownames(meth_50) %in% keep_msp, ]
cov_50   <- cov_50[rownames(cov_50)   %in% keep_msp, ]
cat("After MspI filter:", nrow(meth_50), "loci\n")

# ── 3. TSS proximity filter (-1500/+500bp) + nearest-gene annotation ───────
tss <- fread("hg38_TSS.bed",
             col.names = c("bin","name","chr","strand","txStart","txEnd",
                           "cdsStart","cdsEnd","exonCount","exonStarts",
                           "exonEnds","score","gene","cdsStartStat",
                           "cdsEndStat","exonFrames"))
tss[, tss_pos := ifelse(strand == "+", txStart, txEnd)]
tss[, `:=`(
  win_start = ifelse(strand == "+", tss_pos - 1500, tss_pos - 500),
  win_end   = ifelse(strand == "+", tss_pos + 500,  tss_pos + 1500)
)]
tss[, chr := sub("chr", "", chr)]
tss_windows <- tss[, .(chr, win_start, win_end, gene)]
setkey(tss_windows, chr, win_start, win_end)

locus_dt2 <- data.table(
  locus = rownames(meth_50),
  chr   = sub(":.*", "", rownames(meth_50)),
  start = as.integer(sub(".*:", "", rownames(meth_50)))
)
locus_dt2[, end := start]
setkey(locus_dt2, chr, start, end)

hits_tss   <- foverlaps(locus_dt2, tss_windows, nomatch = 0L)
keep_tss   <- unique(hits_tss$locus)
locus_anno <- hits_tss[, .(gene = gene[1]), by = locus]

meth_50 <- meth_50[rownames(meth_50) %in% keep_tss, ]
cov_50  <- cov_50[rownames(cov_50)   %in% keep_tss, ]

cat("After TSS filter:", nrow(meth_50), "loci\n")
cat("Unique genes:    ", length(unique(locus_anno$gene)), "\n")

# ── 4. Impute NAs, build beta / M-value matrices ────────────────────────────
# CRITICAL: impute_rowmean() must be applied BEFORE any gene-level
# aggregation, not after — applying it after aggregation caused all-NA
# gene/sample combinations to return Cov = 0 and Inf delta-beta values for
# ~1,174 genes in an earlier whole-genome run. This ordering (impute here at
# locus level, aggregate to genes later in 03_gene_level_DE.R) is correct.
impute_rowmean <- function(mat) {
  mat       <- as.matrix(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  for (i in which(rowSums(is.na(mat)) > 0)) {
    mat[i, is.na(mat[i, ])] <- round(row_means[i])
  }
  mat
}

meth_clean <- impute_rowmean(meth_50)
cov_clean  <- impute_rowmean(cov_50)
cat("NA in meth_clean:", sum(is.na(meth_clean)), "\n")
cat("Loci available:  ", nrow(meth_clean), "\n")

beta_full <- meth_clean / cov_clean
beta_full[beta_full <= 0] <- 0.001
beta_full[beta_full >= 1] <- 0.999
mvalue_full <- log2(beta_full / (1 - beta_full))
cat("Non-finite in mvalue_full:", sum(!is.finite(mvalue_full)), "\n")

# ── 5. Sample groups ─────────────────────────────────────────────────────────
group  <- factor(ifelse(grepl("^DG", colnames(meth_clean)), "DG", "NG"))
dg_idx <- which(group == "DG")
ng_idx <- which(group == "NG")
cat("Group breakdown:\n")
print(table(group))

# ── 6. Save for downstream scripts ──────────────────────────────────────────
saveRDS(meth_50,      file.path(OUT_DIR, "rds", "meth_50.rds"))
saveRDS(cov_50,       file.path(OUT_DIR, "rds", "cov_50.rds"))
saveRDS(meth_clean,   file.path(OUT_DIR, "rds", "meth_clean.rds"))
saveRDS(cov_clean,    file.path(OUT_DIR, "rds", "cov_clean.rds"))
saveRDS(beta_full,    file.path(OUT_DIR, "rds", "beta_full.rds"))
saveRDS(mvalue_full,  file.path(OUT_DIR, "rds", "mvalue_full.rds"))
saveRDS(locus_anno,   file.path(OUT_DIR, "rds", "locus_anno.rds"))
saveRDS(group,        file.path(OUT_DIR, "rds", "group.rds"))
saveRDS(dg_idx,       file.path(OUT_DIR, "rds", "dg_idx.rds"))
saveRDS(ng_idx,       file.path(OUT_DIR, "rds", "ng_idx.rds"))
