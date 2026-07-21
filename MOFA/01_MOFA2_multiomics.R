#!/usr/bin/env Rscript
# ============================================================================
# MOFA2 multi-omics integration: DG versus NG
# Methods-aligned GitHub version
# ============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(edgeR)
  library(data.table)
  library(readxl)
  library(dplyr)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

# ---- User-configurable paths -------------------------------------------------
RRBS_RDATA <- "RRBS_genome_results_promoter.RData"
RNA_RDATA <- "RNAseq_DG_vs_NG_results.RData"
LIPID_RDATA <- "lipidomics_DG_vs_NG_results.RData"
CLINICAL_FILE <- "clinical_metadata.xlsx"
REFGENE_FILE <- "hg38_refGene.txt.gz"
OUT_DIR <- "MOFA_results"

N_FACTORS <- 15L
TOP_VARIABLE <- 5000L
MAX_MISSING_FRACTION <- 0.50
TOP_FACTOR1_FEATURES <- 200L
SEED <- 42L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers -----------------------------------------------------------------
strip_group_prefix <- function(x) {
  x <- sub("^(DG|NG)[-_]", "", x)
  sub("^X", "", x)
}

top_variable_features <- function(mat, n = TOP_VARIABLE) {
  vars <- apply(mat, 1, var, na.rm = TRUE)
  vars[!is.finite(vars)] <- -Inf
  ord <- order(vars, decreasing = TRUE)
  mat[ord[seq_len(min(n, nrow(mat)))], , drop = FALSE]
}

filter_missing <- function(mat, max_fraction = MAX_MISSING_FRACTION) {
  mat[rowMeans(is.na(mat)) <= max_fraction, , drop = FALSE]
}

regress_total_signal <- function(mat) {
  total_signal <- colSums(mat, na.rm = TRUE)

  corrected <- t(apply(mat, 1, function(feature) {
    ok <- is.finite(feature) & is.finite(total_signal)
    out <- rep(NA_real_, length(feature))

    if (sum(ok) >= 3L && length(unique(total_signal[ok])) >= 2L) {
      out[ok] <- residuals(lm(feature[ok] ~ total_signal[ok]))
    } else {
      out[ok] <- feature[ok] - mean(feature[ok], na.rm = TRUE)
    }
    out
  }))

  dimnames(corrected) <- dimnames(mat)
  list(matrix = corrected, total_signal = total_signal)
}

read_clinical_metadata <- function(path) {
  raw <- read_excel(path, col_names = FALSE)
  top <- as.character(raw[1, ])
  sub <- as.character(raw[2, ])

  for (i in seq_along(top)) {
    if ((is.na(top[i]) || top[i] == "NA") && i > 1L) top[i] <- top[i - 1L]
  }

  nm <- ifelse(is.na(sub) | sub == "NA" | sub == top,
               top, paste(top, sub, sep = "_"))
  nm <- make.names(nm, unique = TRUE)
  dat <- read_excel(path, skip = 2, col_names = nm)

  candidate <- list(
    subject = c("Subject", "subject"),
    group = c("Group", "group"),
    hba1c = c("HbA1c", "hba1c_pct"),
    glucose = c("Glucose", "glucose"),
    homa_ir = c("HOMA.Insulin.Resistance", "homa_ir"),
    bmi = c("BMI", "bmi"),
    age = c("Age", "age"),
    triglycerides = c("Triglycerides", "triglycerides"),
    hdl = c("HDL", "hdl"),
    ldl = c("LDL", "ldl"),
    total_cholesterol = c("Total.Cholesterol", "total_cholesterol", "total_chol")
  )

  pull_first <- function(keys) {
    hit <- intersect(keys, names(dat))
    if (length(hit) == 0L) return(rep(NA, nrow(dat)))
    dat[[hit[1L]]]
  }

  data.frame(
    subject = trimws(as.character(as.integer(pull_first(candidate$subject)))),
    group = as.character(pull_first(candidate$group)),
    hba1c = as.numeric(pull_first(candidate$hba1c)),
    glucose = as.numeric(pull_first(candidate$glucose)),
    homa_ir = as.numeric(pull_first(candidate$homa_ir)),
    bmi = as.numeric(pull_first(candidate$bmi)),
    age = as.numeric(pull_first(candidate$age)),
    triglycerides = as.numeric(pull_first(candidate$triglycerides)),
    hdl = as.numeric(pull_first(candidate$hdl)),
    ldl = as.numeric(pull_first(candidate$ldl)),
    total_cholesterol = as.numeric(pull_first(candidate$total_cholesterol)),
    stringsAsFactors = FALSE
  )
}

annotate_methylation <- function(loci, refgene_path) {
  refgene <- fread(cmd = paste("zcat", shQuote(refgene_path)), header = FALSE)
  setnames(
    refgene,
    c("bin", "name", "chrom", "strand", "txStart", "txEnd",
      "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
      "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  )
  refgene[, tss := ifelse(strand == "+", txStart, txEnd)]
  gene_tss <- unique(refgene[, .(chrom, tss, gene = name2)])

  locus_chr <- ifelse(grepl("^chr", loci),
                      sub(":.*", "", loci),
                      paste0("chr", sub(":.*", "", loci)))
  locus_pos <- as.integer(sub(".*:", "", loci))

  locus_gr <- GRanges(locus_chr, IRanges(locus_pos, locus_pos))
  gene_gr <- GRanges(gene_tss$chrom, IRanges(gene_tss$tss, gene_tss$tss))
  idx <- nearest(locus_gr, gene_gr, ignore.strand = TRUE)

  data.frame(
    feature = loci,
    gene_symbol = gene_tss$gene[idx],
    distance_to_tss = locus_pos - gene_tss$tss[idx],
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, name) {
  fwrite(as.data.table(x), file.path(OUT_DIR, name), sep = "\t", na = "NA")
}

# ---- 1. Load the three processed omics views --------------------------------
cat("Loading processed omics matrices...\n")

rrbs <- new.env(parent = emptyenv())
load(RRBS_RDATA, envir = rrbs)
if (!exists("mvalue_mat", envir = rrbs, inherits = FALSE)) {
  stop("RRBS file must contain 'mvalue_mat'.")
}
methylation <- as.matrix(rrbs$mvalue_mat)

rna <- new.env(parent = emptyenv())
load(RNA_RDATA, envir = rna)
if (!exists("y", envir = rna, inherits = FALSE)) {
  stop("RNA file must contain the filtered, TMM-normalised edgeR object 'y'.")
}
transcriptomics <- cpm(rna$y, log = TRUE, prior.count = 1)

lipid <- new.env(parent = emptyenv())
load(LIPID_RDATA, envir = lipid)
if (!exists("log2_matrix", envir = lipid, inherits = FALSE)) {
  stop("Lipidomics file must contain 'log2_matrix'.")
}
lipidomics <- as.matrix(lipid$log2_matrix)

storage.mode(methylation) <- "numeric"
storage.mode(transcriptomics) <- "numeric"
storage.mode(lipidomics) <- "numeric"

# ---- 2. Harmonise sample identifiers ----------------------------------------
colnames(methylation) <- strip_group_prefix(colnames(methylation))
colnames(transcriptomics) <- strip_group_prefix(colnames(transcriptomics))
colnames(lipidomics) <- strip_group_prefix(colnames(lipidomics))

common_samples <- Reduce(
  intersect,
  list(colnames(transcriptomics), colnames(methylation), colnames(lipidomics))
)
common_samples <- sort(common_samples)

if (length(common_samples) == 0L) stop("No samples were shared by all three views.")
if (length(common_samples) != 64L) {
  warning("Methods report 64 complete multi-omics samples; observed ",
          length(common_samples), ".")
}

transcriptomics <- transcriptomics[, common_samples, drop = FALSE]
methylation <- methylation[, common_samples, drop = FALSE]
lipidomics <- lipidomics[, common_samples, drop = FALSE]

cat("Common samples:", length(common_samples), "\n")

# ---- 3. Feature filtering ----------------------------------------------------
transcriptomics <- filter_missing(top_variable_features(transcriptomics))
methylation <- filter_missing(top_variable_features(methylation))
lipidomics <- filter_missing(lipidomics)

cat("Transcriptomic features:", nrow(transcriptomics), "\n")
cat("Methylation features:", nrow(methylation), "\n")
cat("Lipidomic features:", nrow(lipidomics), "\n")

if (nrow(transcriptomics) != 5000L) warning("Expected 5,000 RNA features.")
if (nrow(methylation) != 5000L) warning("Expected 5,000 methylation features.")
if (nrow(lipidomics) != 314L) warning("Methods report 314 lipid features for MOFA.")

# ---- 4. Regress sequencing depth / total signal -----------------------------
rna_depth <- regress_total_signal(transcriptomics)
meth_depth <- regress_total_signal(methylation)

transcriptomics_corrected <- rna_depth$matrix
methylation_corrected <- meth_depth$matrix

write_tsv(
  data.frame(sample = common_samples,
             RNA_total_signal = rna_depth$total_signal,
             methylation_total_signal = meth_depth$total_signal),
  "depth_covariates.tsv"
)

# ---- 5. Train MOFA2 ----------------------------------------------------------
cat("Training MOFA2 model...\n")

mofa <- create_mofa(list(
  transcriptomics = transcriptomics_corrected,
  lipidomics = lipidomics,
  methylation = methylation_corrected
))

data_options <- get_default_data_options(mofa)
data_options$scale_views <- TRUE

model_options <- get_default_model_options(mofa)
model_options$num_factors <- N_FACTORS

training_options <- get_default_training_options(mofa)
training_options$convergence_mode <- "slow"
training_options$seed <- SEED

mofa <- prepare_mofa(
  mofa,
  data_options = data_options,
  model_options = model_options,
  training_options = training_options
)

model <- run_mofa(
  mofa,
  outfile = file.path(OUT_DIR, "MOFA_model.hdf5"),
  use_basilisk = TRUE
)
saveRDS(model, file.path(OUT_DIR, "MOFA_model.rds"))

# ---- 6. Variance explained and dominant factors -----------------------------
variance <- get_variance_explained(model)$r2_per_factor[[1]]
variance_df <- data.frame(
  factor = rownames(variance),
  variance,
  total_variance = rowSums(variance),
  check.names = FALSE
)
write_tsv(variance_df, "variance_explained.tsv")

dominant_factor <- sapply(colnames(variance), function(view) {
  as.integer(sub("Factor", "", rownames(variance)[which.max(variance[, view])]))
})
dominant_df <- data.frame(
  view = names(dominant_factor),
  factor = as.integer(dominant_factor)
)
write_tsv(dominant_df, "layer_dominant_factors.tsv")

# ---- 7. Attach group and clinical metadata ----------------------------------
clinical <- read_clinical_metadata(CLINICAL_FILE)
metadata <- clinical[match(common_samples, clinical$subject), , drop = FALSE]

if (anyNA(metadata$group)) stop("Missing DG/NG group labels after metadata matching.")
metadata$sample <- common_samples
metadata$group <- factor(metadata$group, levels = c("DG", "NG"))

factor_scores <- as.data.frame(get_factors(model)[[1]])
factor_scores$sample <- rownames(factor_scores)
factor_scores <- factor_scores[match(common_samples, factor_scores$sample), , drop = FALSE]
factor_scores$group <- metadata$group

write_tsv(factor_scores, "factor_scores.tsv")
write_tsv(metadata, "harmonised_clinical_metadata.tsv")

# ---- 8. Group differences in factor scores ----------------------------------
factor_names <- grep("^Factor", names(factor_scores), value = TRUE)
group_results <- lapply(factor_names, function(f) {
  tt <- t.test(factor_scores[[f]] ~ factor_scores$group)
  data.frame(
    factor = f,
    mean_DG = mean(factor_scores[[f]][factor_scores$group == "DG"]),
    mean_NG = mean(factor_scores[[f]][factor_scores$group == "NG"]),
    difference_DG_minus_NG =
      mean(factor_scores[[f]][factor_scores$group == "DG"]) -
      mean(factor_scores[[f]][factor_scores$group == "NG"]),
    p_value = tt$p.value
  )
})
group_results <- rbindlist(group_results)
group_results[, FDR := p.adjust(p_value, method = "BH")]
write_tsv(group_results, "factor_group_tests.tsv")

# ---- 9. Factor-clinical Spearman correlations -------------------------------
clinical_vars <- c(
  "hba1c", "glucose", "homa_ir", "bmi", "age",
  "triglycerides", "hdl", "ldl", "total_cholesterol"
)

cor_results <- rbindlist(lapply(factor_names, function(f) {
  rbindlist(lapply(clinical_vars, function(v) {
    ok <- complete.cases(factor_scores[[f]], metadata[[v]])
    if (sum(ok) < 3L) {
      return(data.frame(factor = f, variable = v,
                        rho = NA_real_, p_value = NA_real_))
    }
    ct <- suppressWarnings(cor.test(
      factor_scores[[f]][ok],
      metadata[[v]][ok],
      method = "spearman",
      exact = FALSE
    ))
    data.frame(
      factor = f,
      variable = v,
      rho = unname(ct$estimate),
      p_value = ct$p.value
    )
  }))
}))
cor_results[, FDR := p.adjust(p_value, method = "BH")]
write_tsv(cor_results, "factor_clinical_correlations.tsv")

# ---- 10. Extract top-weighted features --------------------------------------
weights <- get_weights(model, as.data.frame = TRUE)
weights[, abs_weight := abs(value)]
setorder(weights, factor, view, -abs_weight)

top_weights <- weights[, head(.SD, 20L), by = .(factor, view)]
write_tsv(top_weights, "top20_weights_per_factor_view.tsv")

meth_annot <- annotate_methylation(
  unique(weights[view == "methylation", feature]),
  REFGENE_FILE
)
write_tsv(meth_annot, "methylation_feature_annotation.tsv")

rna_ids <- unique(weights[view == "transcriptomics", feature])
rna_symbols <- mapIds(
  org.Hs.eg.db,
  keys = sub("\\..*$", "", rna_ids),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
rna_annot <- data.frame(
  feature = rna_ids,
  gene_symbol = unname(rna_symbols),
  stringsAsFactors = FALSE
)
write_tsv(rna_annot, "transcriptomic_feature_annotation.tsv")

# ---- 11. GO and KEGG enrichment for top Factor 1 features -------------------
factor1 <- "Factor1"
f1_rna <- weights[
  factor == factor1 & view == "transcriptomics"
][order(-abs_weight)]
f1_meth <- weights[
  factor == factor1 & view == "methylation"
][order(-abs_weight)]

f1_rna <- head(f1_rna, TOP_FACTOR1_FEATURES)
f1_meth <- head(f1_meth, TOP_FACTOR1_FEATURES)

rna_genes <- na.omit(rna_annot$gene_symbol[
  match(f1_rna$feature, rna_annot$feature)
])
meth_genes <- na.omit(meth_annot$gene_symbol[
  match(f1_meth$feature, meth_annot$feature)
])
factor1_genes <- unique(c(rna_genes, meth_genes))

background_symbols <- unique(c(
  na.omit(rna_annot$gene_symbol),
  na.omit(meth_annot$gene_symbol)
))

gene_map <- bitr(
  factor1_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
background_map <- bitr(
  background_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

if (nrow(gene_map) > 0L && nrow(background_map) > 0L) {
  go <- enrichGO(
    gene = unique(gene_map$ENTREZID),
    universe = unique(background_map$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  kegg <- enrichKEGG(
    gene = unique(gene_map$ENTREZID),
    universe = unique(background_map$ENTREZID),
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )

  go_df <- as.data.frame(go)
  kegg_df <- as.data.frame(kegg)
  write_tsv(go_df, "Factor1_GO_BP_all.tsv")
  write_tsv(kegg_df, "Factor1_KEGG_all.tsv")
  write_tsv(go_df[go_df$p.adjust < 0.05, , drop = FALSE],
            "Factor1_GO_BP_FDR_0.05.tsv")
  write_tsv(kegg_df[kegg_df$p.adjust < 0.05, , drop = FALSE],
            "Factor1_KEGG_FDR_0.05.tsv")
}

# ---- 12. Save matrices needed by the RF script -------------------------------
save(
  transcriptomics_corrected,
  methylation_corrected,
  lipidomics,
  common_samples,
  metadata,
  factor_scores,
  dominant_factor,
  file = file.path(OUT_DIR, "MOFA_RF_inputs.RData")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

cat("MOFA2 analysis complete. Results saved to:", OUT_DIR, "\n")
