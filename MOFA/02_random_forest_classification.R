#!/usr/bin/env Rscript
# ============================================================================
# Random-forest classification for single- and multi-omics feature sets
# 5-fold cross-validation repeated 10 times
# ============================================================================

suppressPackageStartupMessages({
  library(randomForest)
  library(pROC)
  library(data.table)
})

INPUT_RDATA <- "MOFA_results/MOFA_RF_inputs.RData"
OUT_DIR <- "random_forest_results"

N_TREES <- 500L
N_FOLDS <- 5L
N_REPEATS <- 10L
SEED <- 42L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
load(INPUT_RDATA)

# ---- Helpers -----------------------------------------------------------------
make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  folds <- integer(length(y))

  for (level in levels(y)) {
    idx <- which(y == level)
    idx <- sample(idx, length(idx))
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

safe_rf_predict <- function(train_x, train_y, test_x) {
  train_x <- as.data.frame(train_x, check.names = FALSE)
  test_x <- as.data.frame(test_x, check.names = FALSE)

  # Remove zero-variance predictors using training data only.
  keep <- vapply(train_x, function(x) {
    x <- x[is.finite(x)]
    length(x) > 1L && length(unique(x)) > 1L
  }, logical(1))

  train_x <- train_x[, keep, drop = FALSE]
  test_x <- test_x[, keep, drop = FALSE]

  if (ncol(train_x) == 0L) stop("No usable predictors in this training split.")

  fit <- randomForest(
    x = train_x,
    y = train_y,
    ntree = N_TREES
  )

  predict(fit, newdata = test_x, type = "prob")[, "DG"]
}

evaluate_feature_set <- function(name, x, y) {
  stopifnot(identical(rownames(x), names(y)))
  predictions <- vector("list", N_REPEATS * N_FOLDS)
  cursor <- 1L

  for (repeat_id in seq_len(N_REPEATS)) {
    fold_id <- make_stratified_folds(
      y,
      k = N_FOLDS,
      seed = SEED + repeat_id
    )

    for (fold in seq_len(N_FOLDS)) {
      test_idx <- which(fold_id == fold)
      train_idx <- setdiff(seq_along(y), test_idx)

      probability <- safe_rf_predict(
        train_x = x[train_idx, , drop = FALSE],
        train_y = y[train_idx],
        test_x = x[test_idx, , drop = FALSE]
      )

      predictions[[cursor]] <- data.frame(
        feature_set = name,
        repeat = repeat_id,
        fold = fold,
        sample = names(y)[test_idx],
        truth = y[test_idx],
        probability_DG = probability,
        stringsAsFactors = FALSE
      )
      cursor <- cursor + 1L
    }
  }

  pred <- rbindlist(predictions)
  roc_obj <- roc(
    response = pred$truth,
    predictor = pred$probability_DG,
    levels = c("NG", "DG"),
    direction = "<",
    quiet = TRUE
  )

  list(
    predictions = pred,
    auc = as.numeric(auc(roc_obj)),
    roc = data.frame(
      feature_set = name,
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      threshold = roc_obj$thresholds
    )
  )
}

# ---- Build aligned feature matrices -----------------------------------------
samples <- common_samples
group <- factor(metadata$group[match(samples, metadata$subject)],
                levels = c("NG", "DG"))
names(group) <- samples

factor_names <- grep("^Factor", names(factor_scores), value = TRUE)
factor_matrix <- as.matrix(
  factor_scores[match(samples, factor_scores$sample), factor_names, drop = FALSE]
)
rownames(factor_matrix) <- samples

transcriptomic_x <- t(transcriptomics_corrected[, samples, drop = FALSE])
methylation_x <- t(methylation_corrected[, samples, drop = FALSE])
lipidomic_x <- t(lipidomics[, samples, drop = FALSE])

dominant <- unique(as.integer(dominant_factor))
if (length(dominant) < 3L) {
  warning("Fewer than three unique layer-dominant factors were identified.")
}
dominant <- dominant[seq_len(min(3L, length(dominant)))]

factor_sets <- list(
  MOFA_all_15 = factor_matrix,
  MOFA_three_dominant = factor_matrix[, paste0("Factor", dominant), drop = FALSE]
)

for (i in seq_along(dominant)) {
  factor_sets[[paste0("MOFA_dominant_", i)]] <-
    factor_matrix[, paste0("Factor", dominant[i]), drop = FALSE]
}

# Five factor configurations are expected:
# all 15 factors, the three dominant factors together, and each dominant factor.
if (length(factor_sets) != 5L) {
  warning("Expected five MOFA factor configurations; constructed ",
          length(factor_sets), ".")
}

feature_sets <- c(
  factor_sets,
  list(
    transcriptomics = transcriptomic_x,
    lipidomics = lipidomic_x,
    methylation = methylation_x
  )
)

# ---- Repeated stratified cross-validation -----------------------------------
cat("Evaluating", length(feature_sets), "feature sets...\n")

evaluations <- lapply(names(feature_sets), function(name) {
  cat("  ", name, "\n")
  evaluate_feature_set(name, feature_sets[[name]], group)
})

predictions <- rbindlist(lapply(evaluations, `[[`, "predictions"))
roc_points <- rbindlist(lapply(evaluations, `[[`, "roc"))
performance <- data.frame(
  feature_set = names(feature_sets),
  AUC = vapply(evaluations, `[[`, numeric(1), "auc"),
  predictors = vapply(feature_sets, ncol, integer(1)),
  stringsAsFactors = FALSE
)
performance <- performance[order(performance$AUC, decreasing = TRUE), ]

fwrite(predictions,
       file.path(OUT_DIR, "cross_validated_predictions.tsv"),
       sep = "\t")
fwrite(roc_points,
       file.path(OUT_DIR, "pooled_ROC_coordinates.tsv"),
       sep = "\t")
fwrite(performance,
       file.path(OUT_DIR, "model_performance_AUC.tsv"),
       sep = "\t")

save(
  feature_sets,
  predictions,
  roc_points,
  performance,
  file = file.path(OUT_DIR, "random_forest_results.RData")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

print(performance)
cat("Random-forest analysis complete. Results saved to:", OUT_DIR, "\n")
