###############################################################
#' 06_SVM_LOOCV.R
#'
#' @description
#' Leave-one-out cross-validated (LOOCV) SVM (radial kernel, via caret)
#' classifier for Case vs Control, using rlog-transformed expression of a
#' candidate biomarker gene set, with ROC/AUC, confusion matrix, and
#' prediction distribution plots. Companion script to 05_LOOCV_ElasticNet.R,
#' using the same gene set / CV design but an SVM instead of elastic net.
#'
#' @details
#' - Cohort: 24 Cases / 24 Controls
#' - 19 differentially expressed genes from the stringent cutoff, or the
#'   original 19 candidate biomarkers plus related/borderline genes when
#'   `RELAX_FLAG` is TRUE (see below).
#' - `RELAX_FLAG = FALSE`: use only the original 19 biomarkers passing the
#'   stringent cutoff (no relaxed/borderline genes added).
#' - Hyperparameters are fixed (not tuned) at `sigma = 1/ncol(X)`, `C = 1`,
#'   since LOOCV leaves no room for a safe inner tuning loop on top of the
#'   outer leave-one-out loop.
#'
#' @inputs
#' - ./data/calcDEG.Rdata
#'     Must contain `rld` (DESeqTransform, rlog-transformed data) and
#'     `dds_meta` (fitted DESeqDataSet with Group in colData). Produced by
#'     04_CalcDEG.R.
#' - ./results/04_CalcDEG/deg_0.05_1/04_GeneExpr_Shrunk_Filtered.tsv
#'     Filtered, shrunken DESeq2 results (the candidate gene list). NOTE:
#'     this path is hardcoded and does not match the timestamped
#'     `out_dir` naming pattern produced by the current 04_CalcDEG.R
#'     (`deg_<pAdj_cutoff>_<log2FC_cutoff>_<date>`) — update this path if
#'     the folder name from your 04_CalcDEG.R run differs.
#'
#' @outputs (under out_dir, see RELAX_FLAG below)
#' - SVM_CARET_LOOCV_AUC_CM.txt      Confusion matrix (0.5 threshold) + AUC
#' - SVM_CARET_LOOCV_ROC.pdf         ROC curve
#' - SVM_CARET_LOOCV_Prob.pdf        Predicted probability jitter plot
#' - SVM_CARET_LOOCV_Prob_Box.pdf    Predicted probability box plot
#' - SVM_CARET_LOOCV_Predictions.csv Per-sample truth + predicted probability
#' - ./data/SVM_LOOCV.Rdata          Saved workspace (all objects)
#'
#' @param RELAX_FLAG
#' If TRUE, uses the relaxed gene set (19 core + 5 borderline candidates) and
#' writes to a separate "_RELAXED" output directory instead of overwriting
#' the strict-set results.
###############################################################

# ----------------------------
# Load required packages
# ----------------------------
# packages <- c("caret", "pROC", "ggplot2", "dplyr")
# for (p in packages) {
#   if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
# }
library(DESeq2)
library(caret)
library(pROC)
library(knitr)
library(tidyverse)

# ------------------------------------------------------------------------------
# 0. Config
# ------------------------------------------------------------------------------

RELAX_FLAG <- FALSE  # TRUE = include the 5 borderline candidate genes

out_dir <- if (RELAX_FLAG) {
  "./results/06_SVM_CARET_LOOCV_RELAXED"
} else {
  "./results/06_SVM_CARET_LOOCV"
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load upstream DESeq2 results
# ------------------------------------------------------------------------------

env <- new.env()
load("./data/calcDEG.Rdata", envir = env)

if (!"rld" %in% ls(env)) {
  stop("Object 'rld' not found in RData file")
}

rld <- env[["rld"]]
dds_meta <- env[["dds_meta"]]

# ------------------------------------------------------------------------------
# 2. Input data: candidate gene set + expression matrix
# ------------------------------------------------------------------------------

sel_genes <- read_tsv("./results/04_CalcDEG/deg_0.05_1/04_GeneExpr_Shrunk_Filtered.tsv") %>%
  select(Gene)

if (RELAX_FLAG) {
  # Borderline candidates: pAdj < 0.05; |FC| > 1.86 (did not meet the
  # stricter log2FC_Cutoff used in 04_CalcDEG.R)
  borderline_genes <- c("SGPP1", "SKP1", "ZNF706", "TRIM39", "FCHO1")
  sel_genes <- bind_rows(sel_genes, tibble(Gene = borderline_genes)) %>%
    distinct() %>%
    arrange(Gene)
}

# Subset rlog matrix to candidate genes, transpose to samples x genes, and
# z-score each gene across samples for the SVM (caret also re-centers /
# re-scales per training fold via preProcess below)
rlog_mat <- assay(rld)
X <- t(rlog_mat[rownames(rlog_mat) %in% sel_genes$Gene, ])
X <- scale(X)

# Class labels from colData
y <- colData(dds_meta)$Group

# Combine into a single modeling data frame (samples x [genes, condition])
data <- cbind(X, condition = tibble(condition = y))

# ------------------------------------------------------------------------------
# 3. LOOCV setup
# ------------------------------------------------------------------------------

n <- nrow(data)

loocv_probs <- numeric(n)
loocv_truth <- data$condition

# ------------------------------------------------------------------------------
# 4. caret training control
# ------------------------------------------------------------------------------
# method = "none": caret fits the model as-is with the fixed hyperparameters
# below, with no inner resampling. This is intentional under LOOCV — an inner
# CV loop per outer fold would be expensive and, with this few samples per
# fold, unstable.

ctrl <- trainControl(
  method = "none",
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "none"
)

# Fixed hyperparameters (safe under LOOCV, not tuned per fold)
svm_grid <- expand.grid(
  sigma = 1 / ncol(data[, -ncol(data)]),
  C = 1
)

# ------------------------------------------------------------------------------
# 5. Leave-one-out CV loop
# ------------------------------------------------------------------------------
# For each sample: fit a radial-kernel SVM on the remaining n-1 samples at
# the fixed hyperparameters above, then predict P(Case) for the held-out
# sample.

set.seed(14532)

for (i in 1:n) {

  train_data <- data[-i, ]
  test_data  <- data[i, , drop = FALSE]

  svm_fit <- train(
    condition ~ .,
    data = train_data,
    method = "svmRadial",
    trControl = ctrl,
    tuneGrid = svm_grid,
    preProcess = c("center", "scale"),
    metric = "ROC"
  )

  loocv_probs[i] <- predict(
    svm_fit,
    newdata = test_data,
    type = "prob"
  )[, "Case"]
}

# ------------------------------------------------------------------------------
# 6. ROC curve and AUC
# ------------------------------------------------------------------------------

roc_obj <- roc(
  response = loocv_truth,
  predictor = loocv_probs,
  levels = c("Control", "Case"),
  direction = "<",
  smooth = TRUE
)

auc_value <- auc(roc_obj)
cat("LOOCV AUC =", round(auc_value, 3), "\n")

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

# ------------------------------------------------------------------------------
# 7. ROC plot
# ------------------------------------------------------------------------------

roc_plot <- ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(color = "blue", linewidth = 1.3) +
  geom_abline(linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 14) +
  labs(
    title = sprintf("SVM (caret) LOOCV ROC (AUC = %.3f)", auc_value),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )
roc_plot

# ------------------------------------------------------------------------------
# 8. Confusion matrix at 0.5 threshold
# ------------------------------------------------------------------------------

pred_class <- factor(
  ifelse(loocv_probs >= 0.5, "Case", "Control"),
  levels = c("Control", "Case")
)

conf_mat <- confusionMatrix(pred_class, loocv_truth)

sink(paste0(out_dir, "/SVM_CARET_LOOCV_AUC_CM.txt"))
print(conf_mat)
cat("\nSVM Caret LOOCV AUC =", round(auc_value, 3), "\n")
sink()

# ------------------------------------------------------------------------------
# 9. Prediction distribution plots
# ------------------------------------------------------------------------------

pred_df <- data.frame(
  Sample = rownames(data),
  Truth = loocv_truth,
  Probability = loocv_probs
)

# --- 9a. Jitter plot -----------------------------------------------------------
prob_plot <- ggplot(pred_df, aes(Truth, Probability, color = Truth)) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "LOOCV Predicted Probabilities (SVM \u2013 caret)",
    y = "Predicted Probability (Case)",
    x = "True Class"
  )
prob_plot

# --- 9b. Box plot with jittered points and median labels -----------------------
box_plot <- ggplot(pred_df, aes(Truth, Probability, color = Truth)) +
  geom_boxplot(
    width = 0.4,
    outlier.shape = NA,
    alpha = 0.25
  ) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "LOOCV Predicted Probabilities (SVM \u2013 caret)",
    y = "Predicted Probability (Case)",
    x = "True Class"
  ) +
  # Median value labels
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y))),
    vjust = 0.0,
    hjust = -1.8,
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  guides(color = "none")
box_plot

# ------------------------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------------------------

ggsave(paste0(out_dir, "/SVM_CARET_LOOCV_ROC.pdf"), roc_plot, width = 6, height = 5)
ggsave(paste0(out_dir, "/SVM_CARET_LOOCV_Prob.pdf"), prob_plot, width = 6, height = 5)
ggsave(paste0(out_dir, "/SVM_CARET_LOOCV_Prob_Box.pdf"), box_plot, width = 6, height = 5)

write.csv(
  pred_df,
  paste0(out_dir, "/SVM_CARET_LOOCV_Predictions.csv"),
  row.names = FALSE
)

#############################
## Save environment
#############################
# Persist the full workspace so downstream scripts can pick up where this
# one left off.
save(list = ls(), file = "./data/SVM_LOOCV.Rdata")

cat("SVM (caret) LOOCV analysis complete.\n")
