###############################################################
#' 05_LOOCV_ElasticNet.R
#'
#' @description
#' Leave-one-out cross-validated (LOOCV) elastic net logistic regression to
#' classify Case vs Control from rlog-transformed expression of a candidate
#' biomarker gene set, with ROC/AUC, confusion matrix, and prediction
#' distribution plots.
#'
#' @details
#' - Cohort: 24 Cases / 24 Controls
#' - Candidate biomarkers: 19 genes (see input file below), optionally
#'   extended with 5 borderline candidates when `RELAX_FLAG` is TRUE:
#'   SGPP1, SKP1, ZNF706, TRIM39, FCHO1 (padj < 0.05; |FC| > 1.86)
#' - alpha = 0.9 (elastic net mixing parameter, chosen via prior cv.glmnet
#'   optimization; see 05c_ElasticNet_cv_optimization or similar)
#' - At each LOOCV fold, `cv.glmnet` re-selects lambda via internal CV on the
#'   training fold, then predicts the held-out sample at lambda.min.
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
#' - ElasticNet_LOOCV_Coeff.tsv           Per-fold elastic net coefficients
#' - ElasticNet_LOOCV_AUC_CM.txt          Confusion matrix (0.5 threshold) + AUC
#' - ElasticNet_LOOCV_ROC.pdf             ROC curve
#' - ElasticNet_LOOCV_ProbBox.pdf         Predicted probability box plot
#' - ElasticNet_LOOCV_Prob.pdf            Predicted probability jitter plot
#' - ElasticNet_LOOCV_Predictions_GLM.csv Per-sample truth + predicted probability
#' - ./data/ElasticNet_LOOCV.Rdata        Saved workspace (all objects)
#'
#' @param RELAX_FLAG
#' If TRUE, uses the relaxed gene set (19 core + 5 borderline candidates) and
#' writes to a separate "_RELAXED" output directory instead of overwriting
#' the strict-set results.
###############################################################

# ----------------------------
# Load required packages
# ----------------------------
# packages <- c("ggplot2", "pROC", "dplyr", "tidyr", "caret", "glmnet")
# for (p in packages) {
#   if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
# }
library(DESeq2)
library(pROC)
library(caret)
library(knitr)
library(tidyverse)
library(glmnet)

# ------------------------------------------------------------------------------
# 0. Config
# ------------------------------------------------------------------------------

RELAX_FLAG <- FALSE  # TRUE = include the 5 borderline candidate genes

out_dir <- if (RELAX_FLAG) {
  "./results/10_ElasticNet_LOOCV_RELAXED"
} else {
  "./results/10_ElasticNet_LOOCV"
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
# z-score each gene across samples for the elastic net model
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

# Storage for predictions and per-fold coefficients
loocv_probs <- numeric(n)
loocv_truth <- data$condition
coef_tb <- tibble(Gene = c("Intercept", colnames(X)))

# ------------------------------------------------------------------------------
# 4. Leave-one-out cross-validation
# ------------------------------------------------------------------------------
# For each sample: fit an elastic net logistic model (alpha = 0.9) on the
# remaining n-1 samples, with lambda chosen by cv.glmnet's internal CV, then
# predict the probability of "Case" for the held-out sample at lambda.min.

set.seed(15432)

for (i in 1:n) {

  # Split data: leave sample i out for testing
  train_data <- data[-i, ]
  test_data  <- data[i, , drop = FALSE]

  x_train <- as.matrix(train_data[, colnames(train_data) != "condition"])
  y_train <- ifelse(train_data$condition == "Case", 1, 0)

  # NOTE: lambda is not fixed across folds — cv.glmnet re-runs internal CV on
  # each training fold and lambda.min is used for prediction below. (An
  # earlier version explored a fixed lambda = 0.000148516, but a shared
  # lambda across folds can cause convergence issues in some GLM fits.)
  cv_fit <- cv.glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = 0.9,
    type.measure = "deviance"
  )

  x_test <- as.matrix(test_data[, colnames(test_data) != "condition"])

  # Record this fold's coefficients (including intercept) for later inspection
  coef_tb <- coef_tb %>%
    mutate(!!paste0("Coeff_", i) := as.vector(as.matrix(coef(cv_fit))))

  # Predict probability of "Case" for the held-out sample
  loocv_probs[i] <- as.numeric(
    predict(
      cv_fit,
      newx = x_test,
      s = "lambda.min",
      type = "response"
    )
  )
}

write_tsv(coef_tb, paste0(out_dir, "/ElasticNet_LOOCV_Coeff.tsv"))

# ------------------------------------------------------------------------------
# 5. ROC curve and AUC
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

# ROC data frame for plotting
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

# ------------------------------------------------------------------------------
# 6. ROC plot (publication-ready)
# ------------------------------------------------------------------------------

roc_plot <- ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(color = "#3B4CC0", linewidth = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 14) +
  labs(
    title = sprintf("LOOCV ROC Curve (Elastic Net, AUC = %.3f)", auc_value),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )
roc_plot

# ------------------------------------------------------------------------------
# 7. Confusion matrix at 0.5 threshold
# ------------------------------------------------------------------------------

pred_class <- factor(
  ifelse(loocv_probs >= 0.5, "Case", "Control"),
  levels = c("Control", "Case")
)

conf_mat <- confusionMatrix(pred_class, loocv_truth)

sink(paste0(out_dir, "/ElasticNet_LOOCV_AUC_CM.txt"))
print(conf_mat)
cat("\nElastic Net LOOCV AUC =", round(auc_value, 3), "\n")
sink()

# ------------------------------------------------------------------------------
# 8. Prediction distribution plots
# ------------------------------------------------------------------------------
# Shared data frame of per-sample truth + predicted probability, used for
# both the jitter plot and the box plot below.

pred_df <- data.frame(
  Sample = rownames(data),
  Truth = loocv_truth,
  Probability = loocv_probs
)

# --- 8a. Jitter plot -----------------------------------------------------------
prob_plot <- ggplot(pred_df, aes(x = Truth, y = Probability, color = Truth)) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "LOOCV Predicted Probabilities (ElasticNet Regression)",
    y = "Predicted Probability (Case)",
    x = "True Class"
  )
prob_plot

# --- 8b. Box plot with jittered points and median labels -----------------------
box_plot <- ggplot(pred_df, aes(x = Truth, y = Probability, color = Truth)) +
  geom_boxplot(
    width = 0.4,
    outlier.shape = NA,
    alpha = 0.2
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8
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
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(
    title = "LOOCV Predicted Probabilities (ElasticNet Regression)",
    y = "Predicted Probability (Case)",
    x = "True Class"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = "none")

box_plot

# ------------------------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------------------------

ggsave(paste0(out_dir, "/ElasticNet_LOOCV_ROC.pdf"), roc_plot, width = 6, height = 5)
ggsave(paste0(out_dir, "/ElasticNet_LOOCV_ProbBox.pdf"), box_plot, width = 6, height = 5)
ggsave(paste0(out_dir, "/ElasticNet_LOOCV_Prob.pdf"), prob_plot, width = 12, height = 10)

write.csv(
  data.frame(
    Sample = rownames(data),
    Truth = loocv_truth,
    Predicted_Probability = loocv_probs
  ),
  paste0(out_dir, "/ElasticNet_LOOCV_Predictions_GLM.csv"),
  row.names = FALSE
)

#############################
## Save environment
#############################
# Persist the full workspace so downstream scripts can pick up where this
# one left off.
save(list = ls(), file = "./data/ElasticNet_LOOCV.Rdata")

cat("ElasticNet Regression LOOCV analysis complete.\n")
