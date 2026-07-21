#!/usr/bin/env Rscript

# RRBS secondary analysis: whole-genome
# DG versus NG; no MspI filter and no promoter filter at CpG level

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(GenomicRanges)
})

# ---- Configuration ---------------------------------------------------------
RDATA_IN <- "/path/to/RRBS_genome_results_wholegenome.RData"
REFGENE  <- "/path/to/hg38_TSS.bed"
OUT_DIR  <- "results/secondary_whole_genome"
GENE_WINDOW <- 1500L

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

# No MspI or TSS filter is applied here.
meth_clean <- impute_row_mean(meth_50)
cov_clean <- impute_row_mean(cov_50)
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

# ---- Nearest-gene annotation ------------------------------------------------
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
gene_tss <- unique(refgene[!is.na(gene) & gene != "", .(chr, tss, gene)])

locus_pos <- data.table(
  locus = cpg_res$locus,
  chr = sub(":.*", "", cpg_res$locus),
  pos = as.integer(sub(".*:", "", cpg_res$locus))
)
locus_gr <- GRanges(locus_pos$chr, IRanges(locus_pos$pos, locus_pos$pos))
tss_gr <- GRanges(gene_tss$chr, IRanges(gene_tss$tss, gene_tss$tss))
idx <- nearest(locus_gr, tss_gr, ignore.strand = TRUE)

locus_pos[, gene := gene_tss$gene[idx]]
locus_pos[, dist_tss := pos - gene_tss$tss[idx]]
cpg_res <- merge(cpg_res, locus_pos[, .(locus, gene, dist_tss)],
                 by = "locus", all.x = TRUE)

cpg_primary <- cpg_res[FDR < 0.25 & abs(delta_beta) > 0.05]
cpg_exploratory <- cpg_res[PValue < 0.001 & abs(delta_beta) > 0.05]

# ---- Gene-level aggregation within +/-1500 bp of a TSS ---------------------
all_loci <- data.table(
  locus = rownames(meth_clean),
  chr = sub(":.*", "", rownames(meth_clean)),
  pos = as.integer(sub(".*:", "", rownames(meth_clean)))
)
all_loci[, `:=`(start = pos, end = pos)]

tss_windows <- gene_tss[, .(
  chr,
  win_start = pmax(0L, tss - GENE_WINDOW),
  win_end = tss + GENE_WINDOW,
  gene
)]
setkey(all_loci, chr, start, end)
setkey(tss_windows, chr, win_start, win_end)
gene_hits <- foverlaps(
  all_loci, tss_windows,
  by.x = c("chr", "start", "end"),
  by.y = c("chr", "win_start", "win_end"),
  nomatch = 0L
)
locus_gene <- unique(gene_hits[, .(locus, gene)])

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

# ---- Save ------------------------------------------------------------------
fwrite(cpg_res, file.path(OUT_DIR, "whole_genome_cpg_all.tsv"), sep = "\t")
fwrite(cpg_primary, file.path(OUT_DIR, "whole_genome_cpg_FDR025_db005.tsv"), sep = "\t")
fwrite(cpg_exploratory, file.path(OUT_DIR, "whole_genome_cpg_p0001_db005.tsv"), sep = "\t")
fwrite(locus_gene, file.path(OUT_DIR, "whole_genome_TSS1500_locus_gene_map.tsv"), sep = "\t")
fwrite(gene_res, file.path(OUT_DIR, "whole_genome_gene_all.tsv"), sep = "\t")
fwrite(gene_exploratory, file.path(OUT_DIR, "whole_genome_gene_p0001_db005.tsv"), sep = "\t")

save(
  beta, mvalue, meth_clean, cov_clean,
  cpg_res, cpg_primary, cpg_exploratory,
  locus_gene, gene_res, gene_exploratory,
  file = file.path(OUT_DIR, "secondary_whole_genome_results.RData")
)

cat("Coverage-filtered CpGs carried forward:", nrow(meth_clean), "\n")
cat("CpGs tested:", nrow(cpg_res), "\n")
cat("CpGs within +/-1500 bp for gene aggregation:", nrow(locus_gene), "\n")
cat("Genes tested:", nrow(gene_res), "\n")
cat("Finished:", format(Sys.time()), "\n")
