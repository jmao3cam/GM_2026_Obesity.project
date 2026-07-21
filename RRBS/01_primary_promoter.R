#!/usr/bin/env Rscript

# RRBS primary analysis: MspI-filtered, promoter-focused
# DG versus NG

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

# ---- Configuration ---------------------------------------------------------
RDATA_IN <- "/path/to/RRBS_genome_results_wholegenome.RData"
MSPI_BED <- "/path/to/MspI_fragments_hg38.bed"
REFGENE  <- "/path/to/hg38_TSS.bed"
OUT_DIR  <- "results/primary_promoter"

MSPI_MIN <- 40L
MSPI_MAX <- 220L
TSS_UPSTREAM <- 1500L
TSS_DOWNSTREAM <- 500L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ---------------------------------------------------------------
impute_row_mean <- function(x) {
  x <- as.matrix(x)
  means <- rowMeans(x, na.rm = TRUE)
  if (any(!is.finite(means))) stop("At least one locus is missing in every sample.")
  missing <- which(is.na(x), arr.ind = TRUE)
  if (nrow(missing)) x[missing] <- round(means[missing[, 1]])
  x
}

clip_beta <- function(x) pmin(pmax(x, 0.001), 0.999)

run_edger <- function(meth, coverage, group) {
  y <- DGEList(
    counts = meth,
    lib.size = colSums(coverage),
    group = group
  )
  keep <- filterByExpr(y, group = group)
  y <- y[keep, , keep.lib.sizes = FALSE]

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  # No calcNormFactors(): total coverage is used as the sample library size.
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  test <- glmLRT(fit, contrast = c(1, -1)) # DG - NG

  list(table = as.data.table(topTags(test, n = Inf)$table, keep.rownames = "feature"),
       keep = keep)
}

aggregate_rows <- function(mat, mapping) {
  dt <- as.data.table(mat, keep.rownames = "locus")
  dt <- merge(mapping, dt, by = "locus", allow.cartesian = TRUE)
  sample_cols <- setdiff(names(dt), c("locus", "gene"))
  out <- dt[, lapply(.SD, sum), by = gene, .SDcols = sample_cols]
  rn <- out$gene
  out[, gene := NULL]
  ans <- as.matrix(out)
  rownames(ans) <- rn
  ans
}

# ---- Load shared coverage-filtered matrices --------------------------------
loaded <- load(RDATA_IN)
required <- c("meth_50", "cov_50")
if (!all(required %in% loaded)) {
  stop("RData must contain: ", paste(required, collapse = ", "))
}
stopifnot(identical(dim(meth_50), dim(cov_50)))
stopifnot(identical(rownames(meth_50), rownames(cov_50)))

coverage_pass <- rowSums(cov_50 >= 50, na.rm = TRUE) >= ceiling(0.5 * ncol(cov_50))
if (!all(coverage_pass)) stop("Input is not the expected coverage-filtered matrix.")

# ---- MspI fragment filter: 40-220 bp ---------------------------------------
msp <- fread(MSPI_BED, header = FALSE)
if (ncol(msp) < 3L) stop("MspI BED must contain at least three columns.")
setnames(msp, names(msp)[1:3], c("chr", "frag_start", "frag_end"))
msp[, chr := sub("^chr", "", chr)]
msp[, frag_size := frag_end - frag_start]
msp <- msp[frag_size >= MSPI_MIN & frag_size <= MSPI_MAX,
           .(chr, frag_start, frag_end)]

loci <- data.table(
  locus = rownames(meth_50),
  chr = sub(":.*", "", rownames(meth_50)),
  start = as.integer(sub(".*:", "", rownames(meth_50)))
)
loci[, end := start]

setkey(msp, chr, frag_start, frag_end)
setkey(loci, chr, start, end)
msp_hits <- foverlaps(
  loci, msp,
  by.x = c("chr", "start", "end"),
  by.y = c("chr", "frag_start", "frag_end"),
  nomatch = 0L
)
keep_msp <- unique(msp_hits$locus)

meth_primary <- meth_50[keep_msp, , drop = FALSE]
cov_primary <- cov_50[keep_msp, , drop = FALSE]

# ---- Strand-aware promoter filter: -1500/+500 bp ---------------------------
refgene <- fread(REFGENE, header = FALSE)
if (ncol(refgene) < 13L) stop("RefSeq table must contain at least 13 columns.")
setnames(
  refgene,
  names(refgene)[1:13],
  c("bin", "transcript", "chr", "strand", "txStart", "txEnd",
    "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
    "score", "gene")
)
refgene[, chr := sub("^chr", "", chr)]
refgene[, tss := fifelse(strand == "+", txStart, txEnd)]
refgene[, promoter_start := fifelse(
  strand == "+", tss - TSS_UPSTREAM, tss - TSS_DOWNSTREAM
)]
refgene[, promoter_end := fifelse(
  strand == "+", tss + TSS_DOWNSTREAM, tss + TSS_UPSTREAM
)]
promoters <- unique(
  refgene[!is.na(gene) & gene != "",
          .(chr, promoter_start, promoter_end, gene)]
)

promoter_loci <- data.table(
  locus = rownames(meth_primary),
  chr = sub(":.*", "", rownames(meth_primary)),
  start = as.integer(sub(".*:", "", rownames(meth_primary)))
)
promoter_loci[, end := start]

setkey(promoters, chr, promoter_start, promoter_end)
setkey(promoter_loci, chr, start, end)
promoter_hits <- foverlaps(
  promoter_loci, promoters,
  by.x = c("chr", "start", "end"),
  by.y = c("chr", "promoter_start", "promoter_end"),
  nomatch = 0L
)

# A CpG may overlap multiple transcripts/genes. Keep unique locus-gene pairs.
locus_gene <- unique(promoter_hits[, .(locus, gene)])
keep_promoter <- unique(locus_gene$locus)
meth_primary <- meth_primary[keep_promoter, , drop = FALSE]
cov_primary <- cov_primary[keep_promoter, , drop = FALSE]

# ---- Imputation and beta/M-values ------------------------------------------
meth_clean <- impute_row_mean(meth_primary)
cov_clean <- impute_row_mean(cov_primary)
if (any(cov_clean <= 0)) stop("Non-positive coverage remains after imputation.")

beta <- clip_beta(meth_clean / cov_clean)
mvalue <- log2(beta / (1 - beta))

group <- factor(
  ifelse(grepl("^DG", colnames(meth_clean)), "DG", "NG"),
  levels = c("DG", "NG")
)
if (anyNA(group) || any(table(group) == 0)) stop("Could not define DG and NG groups.")
dg <- group == "DG"
ng <- group == "NG"

# ---- CpG-level differential methylation ------------------------------------
cpg_fit <- run_edger(meth_clean, cov_clean, group)
cpg_res <- cpg_fit$table
setnames(cpg_res, "feature", "locus")
cpg_res[, beta_DG := rowMeans(beta[locus, dg, drop = FALSE])]
cpg_res[, beta_NG := rowMeans(beta[locus, ng, drop = FALSE])]
cpg_res[, delta_beta := beta_DG - beta_NG]

# Add all overlapping promoter gene symbols for interpretation.
cpg_gene <- locus_gene[, .(gene = paste(sort(unique(gene)), collapse = ";")), by = locus]
cpg_res <- merge(cpg_res, cpg_gene, by = "locus", all.x = TRUE)

cpg_primary <- cpg_res[FDR < 0.25 & abs(delta_beta) > 0.05]
cpg_exploratory <- cpg_res[PValue < 0.001 & abs(delta_beta) > 0.05]

# ---- Promoter/gene-level aggregation ---------------------------------------
gene_meth <- aggregate_rows(meth_clean, locus_gene)
gene_cov <- aggregate_rows(cov_clean, locus_gene)
stopifnot(identical(rownames(gene_meth), rownames(gene_cov)))

gene_beta <- clip_beta(gene_meth / gene_cov)
gene_fit <- run_edger(gene_meth, gene_cov, group)
gene_res <- gene_fit$table
setnames(gene_res, "feature", "gene")
gene_res[, beta_DG := rowMeans(gene_beta[gene, dg, drop = FALSE])]
gene_res[, beta_NG := rowMeans(gene_beta[gene, ng, drop = FALSE])]
gene_res[, delta_beta := beta_DG - beta_NG]
gene_exploratory <- gene_res[PValue < 0.001 & abs(delta_beta) > 0.05]

# ---- ZNF423 complete-case sensitivity analysis -----------------------------
znf_loci <- locus_gene[gene == "ZNF423", unique(locus)]
znf_sensitivity <- NULL
if (length(znf_loci)) {
  raw_m <- meth_primary[znf_loci, , drop = FALSE]
  raw_c <- cov_primary[znf_loci, , drop = FALSE]
  complete <- colSums(is.na(raw_m) | is.na(raw_c)) == 0
  beta_cc <- clip_beta(raw_m[, complete, drop = FALSE] /
                       raw_c[, complete, drop = FALSE])
  beta_cc_gene <- colMeans(beta_cc)
  group_cc <- droplevels(group[complete])
  if (length(unique(group_cc)) == 2L) {
    wt <- wilcox.test(beta_cc_gene ~ group_cc)
    znf_sensitivity <- data.table(
      n = sum(complete),
      n_DG = sum(group_cc == "DG"),
      n_NG = sum(group_cc == "NG"),
      delta_beta = mean(beta_cc_gene[group_cc == "DG"]) -
                   mean(beta_cc_gene[group_cc == "NG"]),
      wilcoxon_p = wt$p.value
    )
  }
}

# ---- Save ------------------------------------------------------------------
fwrite(cpg_res, file.path(OUT_DIR, "primary_cpg_all.tsv"), sep = "\t")
fwrite(cpg_primary, file.path(OUT_DIR, "primary_cpg_FDR025_db005.tsv"), sep = "\t")
fwrite(cpg_exploratory, file.path(OUT_DIR, "primary_cpg_p0001_db005.tsv"), sep = "\t")
fwrite(gene_res, file.path(OUT_DIR, "primary_gene_all.tsv"), sep = "\t")
fwrite(gene_exploratory, file.path(OUT_DIR, "primary_gene_p0001_db005.tsv"), sep = "\t")
if (!is.null(znf_sensitivity)) {
  fwrite(znf_sensitivity, file.path(OUT_DIR, "ZNF423_complete_case.tsv"), sep = "\t")
}

save(
  beta, mvalue, meth_clean, cov_clean, locus_gene,
  cpg_res, cpg_primary, cpg_exploratory,
  gene_res, gene_exploratory, znf_sensitivity,
  file = file.path(OUT_DIR, "primary_promoter_results.RData")
)

cat("Coverage-filtered CpGs:", nrow(meth_50), "\n")
cat("After MspI filter:", length(keep_msp), "\n")
cat("After promoter filter:", nrow(meth_clean), "\n")
cat("CpGs tested:", nrow(cpg_res), "\n")
cat("Genes tested:", nrow(gene_res), "\n")
cat("Finished:", format(Sys.time()), "\n")
