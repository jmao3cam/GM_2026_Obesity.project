#!/usr/bin/env Rscript
# ============================================================================
# Lipidomics analysis: DG versus NG
# Methods-aligned GitHub version
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})

# ---- User-configurable settings ---------------------------------------------
INPUT_FILE <- "Lipidomicsdata.xlsx"
DATA_SHEET <- "Normalized female data with NA"
OUT_DIR <- "lipidomics_results"

MAX_MISSING <- 5
NORMALITY_ALPHA <- 0.05
SIGNIFICANCE_FDR <- 0.05

# Clinical subject lists used to restrict the lipidomics dataset.
DG_SUBJECTS <- as.character(c(
  69, 109, 112, 115, 128, 129, 132, 133, 134, 143, 156,
  164, 170, 181, 184, 194, 195, 198, 214, 222, 226, 232,
  238, 242, 248, 259, 262, 265, 268, 273, 276, 288
))

NG_SUBJECTS <- as.character(c(
  82, 94, 100, 102, 106, 108, 114, 120, 127, 135, 139,
  147, 152, 155, 158, 168, 169, 174, 178, 187, 190,
  201, 202, 205, 215, 223, 228, 231, 237, 245, 252, 277, 292
))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helper functions --------------------------------------------------------
parse_lipid_class <- function(lipid_names) {
  # Reproduces the class parsing used in the source analysis:
  # retain the portion before "_(" or the first whitespace.
  sub("_\\(.*|\\s.*", "", lipid_names)
}

safe_shapiro <- function(x) {
  x <- x[is.finite(x)]

  # shapiro.test() requires 3-5000 observations and non-constant values.
  if (length(x) < 3L || length(x) > 5000L || length(unique(x)) < 2L) {
    return(NA_real_)
  }

  tryCatch(
    shapiro.test(x)$p.value,
    error = function(e) NA_real_
  )
}

safe_group_test <- function(x_dg, x_ng, use_t_test) {
  x_dg <- x_dg[is.finite(x_dg)]
  x_ng <- x_ng[is.finite(x_ng)]

  if (length(x_dg) < 2L || length(x_ng) < 2L) {
    return(list(p.value = NA_real_, method = "Insufficient observations"))
  }

  if (use_t_test) {
    out <- tryCatch(
      t.test(x_dg, x_ng, var.equal = FALSE),
      error = function(e) NULL
    )
    method <- "Welch t-test"
  } else {
    out <- tryCatch(
      wilcox.test(x_dg, x_ng, exact = FALSE),
      error = function(e) NULL
    )
    method <- "Mann-Whitney U test"
  }

  list(
    p.value = if (is.null(out)) NA_real_ else out$p.value,
    method = method
  )
}

write_tsv <- function(x, filename, row_names = FALSE) {
  out <- as.data.table(x, keep.rownames = row_names)
  fwrite(out, file.path(OUT_DIR, filename), sep = "\t", na = "NA")
}

# ---- 1. Load DNA-normalised lipidomics measurements --------------------------
cat("Loading lipidomics data...\n")

raw <- read_excel(INPUT_FILE, sheet = DATA_SHEET)

if (ncol(raw) < 2L) {
  stop("The selected worksheet must contain a lipid-name column and sample columns.")
}

lipid_names <- as.character(raw[[1]])
if (anyNA(lipid_names) || any(lipid_names == "")) {
  stop("Missing or empty lipid names were found in the first column.")
}
if (anyDuplicated(lipid_names)) {
  stop("Lipid names must be unique. Duplicates detected: ",
       paste(unique(lipid_names[duplicated(lipid_names)]), collapse = ", "))
}

sample_columns <- colnames(raw)[-1]
requested_subjects <- c(DG_SUBJECTS, NG_SUBJECTS)
keep_samples <- sample_columns[sample_columns %in% requested_subjects]

missing_dg <- setdiff(DG_SUBJECTS, keep_samples)
missing_ng <- setdiff(NG_SUBJECTS, keep_samples)

if (length(missing_dg) > 0L) {
  warning("DG subjects absent from the lipidomics sheet: ",
          paste(missing_dg, collapse = ", "))
}
if (length(missing_ng) > 0L) {
  warning("NG subjects absent from the lipidomics sheet: ",
          paste(missing_ng, collapse = ", "))
}
if (length(keep_samples) == 0L) {
  stop("No lipidomics samples matched the DG/NG subject lists.")
}

lipid_mat <- as.matrix(
  data.frame(lapply(raw[, keep_samples, drop = FALSE], as.numeric),
             check.names = FALSE)
)
rownames(lipid_mat) <- lipid_names
colnames(lipid_mat) <- keep_samples
storage.mode(lipid_mat) <- "numeric"

groups <- factor(
  ifelse(keep_samples %in% DG_SUBJECTS, "DG", "NG"),
  levels = c("DG", "NG")
)
names(groups) <- keep_samples

cat("Matched samples:", ncol(lipid_mat), "\n")
cat("DG:", sum(groups == "DG"), "| NG:", sum(groups == "NG"), "\n")
cat("Raw lipid species:", nrow(lipid_mat), "\n")

# ---- 2. Exclude lipids with more than five missing values --------------------
cat("\nApplying missing-value filter...\n")

na_per_lipid <- rowSums(is.na(lipid_mat))
missingness_summary <- data.frame(
  lipid = rownames(lipid_mat),
  missing_values = na_per_lipid,
  retained = na_per_lipid <= MAX_MISSING
)
write_tsv(missingness_summary, "lipid_missingness_summary.tsv")

keep_lipids <- na_per_lipid <= MAX_MISSING
lipid_filtered <- lipid_mat[keep_lipids, , drop = FALSE]

if (nrow(lipid_filtered) == 0L) {
  stop("No lipids passed the missing-value filter.")
}

cat("Lipids retained:", nrow(lipid_filtered), "\n")
if (nrow(lipid_filtered) != 315L) {
  warning(
    "The methods report 315 retained lipids, but this run retained ",
    nrow(lipid_filtered),
    ". Check the input worksheet, subject matching, and missing-value coding."
  )
}

# ---- 3. Assign lipid classes and apply within-class percent normalisation ----
cat("\nApplying within-class percent normalisation...\n")

lipid_classes <- parse_lipid_class(rownames(lipid_filtered))
if (anyNA(lipid_classes) || any(lipid_classes == "")) {
  stop("One or more lipid classes could not be parsed from lipid names.")
}

lipid_percent <- matrix(
  NA_real_,
  nrow = nrow(lipid_filtered),
  ncol = ncol(lipid_filtered),
  dimnames = dimnames(lipid_filtered)
)

for (class_name in unique(lipid_classes)) {
  idx <- which(lipid_classes == class_name)
  class_values <- lipid_filtered[idx, , drop = FALSE]

  # The denominator is the observed total abundance for the class in a sample.
  # Missing component values are ignored here and imputed after normalisation.
  class_total <- colSums(class_values, na.rm = TRUE)
  class_total[class_total <= 0 | !is.finite(class_total)] <- NA_real_

  lipid_percent[idx, ] <- sweep(class_values, 2, class_total, "/") * 100
}

# ---- 4. Impute residual missing values and log2 transform --------------------
cat("\nImputing residual missing values...\n")

positive_values <- lipid_percent[
  is.finite(lipid_percent) & lipid_percent > 0
]

if (length(positive_values) == 0L) {
  stop("No positive values were available after percent normalisation.")
}

half_minimum <- min(positive_values) / 2

# The methods specify imputation of remaining NA values. Non-positive values
# cannot be log2-transformed and are therefore treated as missing at this step.
to_impute <- is.na(lipid_percent) |
             !is.finite(lipid_percent) |
             lipid_percent <= 0
lipid_percent[to_impute] <- half_minimum

log2_matrix <- log2(lipid_percent)

if (any(!is.finite(log2_matrix))) {
  stop("The final log2 matrix contains non-finite values.")
}

write_tsv(
  data.frame(
    lipid = rownames(lipid_percent),
    class = lipid_classes,
    lipid_percent,
    check.names = FALSE
  ),
  "lipid_percent_normalised_imputed.tsv"
)

write_tsv(
  data.frame(
    lipid = rownames(log2_matrix),
    class = lipid_classes,
    log2_matrix,
    check.names = FALSE
  ),
  "lipid_log2_final_matrix.tsv"
)

cat("Half-minimum imputation value:", format(half_minimum, scientific = TRUE), "\n")
cat("Final matrix:", nrow(log2_matrix), "lipids x",
    ncol(log2_matrix), "samples\n")

# ---- 5. Group-wise Shapiro-Wilk normality testing ----------------------------
cat("\nAssessing normality within DG and NG...\n")

dg_samples <- names(groups)[groups == "DG"]
ng_samples <- names(groups)[groups == "NG"]

normality_results <- data.frame(
  lipid = rownames(log2_matrix),
  class = lipid_classes,
  shapiro_p_DG = apply(
    log2_matrix[, dg_samples, drop = FALSE],
    1,
    safe_shapiro
  ),
  shapiro_p_NG = apply(
    log2_matrix[, ng_samples, drop = FALSE],
    1,
    safe_shapiro
  ),
  stringsAsFactors = FALSE
)

normality_results$normal_in_both <- with(
  normality_results,
  !is.na(shapiro_p_DG) &
  !is.na(shapiro_p_NG) &
  shapiro_p_DG > NORMALITY_ALPHA &
  shapiro_p_NG > NORMALITY_ALPHA
)

write_tsv(normality_results, "normality_tests.tsv")

cat("Welch t-test lipids:", sum(normality_results$normal_in_both), "\n")
cat("Mann-Whitney lipids:", sum(!normality_results$normal_in_both), "\n")

# ---- 6. Differential testing: DG versus NG ----------------------------------
cat("\nTesting DG versus NG differences...\n")

results <- lapply(seq_len(nrow(log2_matrix)), function(i) {
  x_dg <- as.numeric(log2_matrix[i, dg_samples])
  x_ng <- as.numeric(log2_matrix[i, ng_samples])
  normal_both <- normality_results$normal_in_both[i]

  test <- safe_group_test(
    x_dg = x_dg,
    x_ng = x_ng,
    use_t_test = isTRUE(normal_both)
  )

  data.frame(
    lipid = rownames(log2_matrix)[i],
    class = lipid_classes[i],
    mean_log2_DG = mean(x_dg),
    mean_log2_NG = mean(x_ng),
    mean_difference_DG_minus_NG = mean(x_dg) - mean(x_ng),
    median_log2_DG = median(x_dg),
    median_log2_NG = median(x_ng),
    shapiro_p_DG = normality_results$shapiro_p_DG[i],
    shapiro_p_NG = normality_results$shapiro_p_NG[i],
    normal_in_both = normal_both,
    test = test$method,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
})

results <- rbindlist(results)

# Benjamini-Hochberg correction is performed independently within each class.
results[, FDR_within_class := p.adjust(p_value, method = "BH"), by = class]
results[, significant_FDR_0.05 :=
          !is.na(FDR_within_class) &
          FDR_within_class < SIGNIFICANCE_FDR]
results[, direction := fifelse(
  mean_difference_DG_minus_NG > 0,
  "Higher in DG",
  fifelse(mean_difference_DG_minus_NG < 0, "Lower in DG", "No difference")
)]

setorder(results, class, FDR_within_class, p_value)

significant_results <- results[significant_FDR_0.05 == TRUE]

write_tsv(results, "DG_vs_NG_all_lipids_per_class_BH.tsv")
write_tsv(
  significant_results,
  "DG_vs_NG_significant_lipids_FDR_0.05.tsv"
)

class_summary <- results[, .(
  lipids_tested = .N,
  significant_FDR_0.05 = sum(significant_FDR_0.05, na.rm = TRUE)
), by = class][order(class)]

write_tsv(class_summary, "lipid_class_summary.tsv")

cat("Lipids tested:", nrow(results), "\n")
cat("Significant lipids (within-class FDR < 0.05):",
    nrow(significant_results), "\n")

# ---- 7. Save reproducible analysis objects ----------------------------------
sample_metadata <- data.frame(
  sample = names(groups),
  group = as.character(groups),
  stringsAsFactors = FALSE
)

save(
  lipid_filtered,
  lipid_classes,
  lipid_percent,
  half_minimum,
  log2_matrix,
  sample_metadata,
  normality_results,
  results,
  significant_results,
  file = file.path(OUT_DIR, "lipidomics_DG_vs_NG_results.RData")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

cat("\nAnalysis complete. Results saved to:", OUT_DIR, "\n")
