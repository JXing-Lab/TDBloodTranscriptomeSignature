###
#' 04_CalcDEG.R — Differential Expression Analysis with DESeq2
#'
#' @description
#' Uses DESeq2 to compute differentially expressed genes (DEGs) for Case vs
#' Control and generates volcano plots, heatmaps, PCA, dispersion, Wald vs
#' Cook's distance diagnostics, and violin/box plots for selected genes.
#' Results and intermediate tables are written to TSV files.
#'
#' @details
#' - Thresholds:
#'   - `pAdj_cutoff` (currently 0.05) and `log2FC_Cutoff` (currently log2(2),
#'     i.e. |FC| >= 2) define "significant" throughout this script.
#'     NOTE: the original exploratory-discovery threshold discussed for this
#'     project was padj <= 0.10; this script currently runs at padj <= 0.05.
#'     Update `pAdj_cutoff` below if the looser exploratory threshold is
#'     intended.
#' - Design: `~ Group + Sex + Age + Group*Sex + Group*Age + Sex*Age + Group*Sex*Age`
#' - Notes: as per the 03_PreprocInput.r header, the 3 smallest samples
#'   (cases 217, 215; control 218) are excluded upstream in preprocessing,
#'   not in this script.
#'
#' @author
#' Last Updated: 2026-04-16
#'
#' @seealso
#' DESeq2 vignette, "Methods changes since the 2014 DESeq2 paper":
#' https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#methods-changes-since-the-2014-deseq2-paper
#'
#' @examples
#' # Run the full analysis (paths are relative to the project):
#' # source("04_CalcDEG.R")
#' # The script runs end-to-end and saves outputs under
#' # ./results/04_CalcDEG/deg_<pAdj_cutoff>_<log2FC_cutoff>_<date>/
#'
#' @keywords DESeq2, RNA-seq, differential expression, volcano plot, PCA, heatmap
###

#############################
## Libraries
#############################
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(stringr)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(reshape2)

#load(file = "./data/calcDEG.Rdata")

#############################
## Global Parameters
#############################
pAdj_cutoff   <- 0.05     # Significance threshold on adjusted p-value
log2FC_Cutoff <- log2(2)  # => abs(FC) >= 2

cat(paste0("pAdj Cutoff:   ",  pAdj_cutoff, "\n"))
cat(paste0("log2FC Cutoff: ",  log2FC_Cutoff, "\n"))

inp_dir <- "./input/"

# Output directory is timestamped and tagged with the cutoffs used, so runs
# with different thresholds don't overwrite each other.
out_dir <- paste0(
  "./results/04_CalcDEG/deg_", pAdj_cutoff, "_",
  log2FC_Cutoff, "_",
  format(Sys.time(), "%Y%m%d")
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# External sample-ID lookup table used to generate human-readable PCA plot
# labels in additional_plots(). NOTE: this is a machine-specific absolute
# path (Windows) — update it or move the file under ./input/ before running
# on a different machine.
sample_id_map_file <- "D:/IQB_PhD/RNASeq/paper/sampleIDmap.tsv"

#############################
## Functions: Plotting DESeq Results
#############################

#' Plot DESeq2 results (Case vs Control), shrink LFC, and export filtered tables
#'
#' Generates: MA plots (raw + shrunken), top-gene count plot, volcano plot
#' (ggplot2, with gene labels), and EnhancedVolcano plot. Writes raw and
#' apeglm-shrunken results tables (full and filtered) to TSV.
#'
#' Relies on the following objects from the calling environment:
#' `out_dir`, `pAdj_cutoff`, `log2FC_Cutoff`.
#'
#' @param desout DESeqDataSet after DESeq() has been run.
#' @param dsname Character label for dataset used in plot titles and file names.
#' @return Invisibly returns NULL. Writes several TSVs and produces plots.
plot_deseq_results <- function(desout, dsname) {
  # List coefficients present
  DESeq2::resultsNames(desout)

  # Case vs Control results (no shrinkage)
  res <- DESeq2::results(desout, name = "Group_Case_vs_Control")

  cat("DESeq results w/o shrinkage\n")
  summary(res)
  sum(res$padj < pAdj_cutoff, na.rm = TRUE)

  # MA plot (raw)
  pdf(file.path(out_dir, "04_MA_Plot_Raw.pdf"), width = 6, height = 4)
  plotMA(res, ylim = c(-2, 2))
  dev.off()

  # Write raw results (all genes)
  as.data.frame(res) %>%
    rownames_to_column(var = "gene") %>%
    write_tsv(file = file.path(out_dir, "04_DESeqResults_Raw_All.tsv"))

  # Filtered (padj & log2FC) and export
  resf <- res %>%
    as.data.frame() %>%
    dplyr::filter(padj <= pAdj_cutoff) %>%
    dplyr::filter(abs(log2FoldChange) >= log2FC_Cutoff)

  write_tsv(
    as.data.frame(resf) %>% tibble::rownames_to_column(var = "Gene"),
    file = file.path(out_dir, "04_DESeqResultsFilt_Raw.tsv")
  )

  # Shrink LFC using apeglm (see citation below)
  sres <- lfcShrink(
    desout,
    coef = "Group_Case_vs_Control",
    res = res,
    type = "apeglm"
  )
  summary(sres)
  sum(sres$padj < pAdj_cutoff, na.rm = TRUE)

  # MA plot (shrunken)
  pdf(file.path(out_dir, "04_MA_Plot_Shrunken.pdf"), width = 6, height = 4)
  plotMA(sres, ylim = c(-2, 2))
  dev.off()

  # Export shrunken full results
  write_tsv(
    as.data.frame(sres) %>% tibble::rownames_to_column(var = "Gene"),
    file = file.path(out_dir, "04_GeneExpr_Shrunk_All.tsv")
  )

  # Filtered shrunken results
  sresf <- sres %>%
    as.data.frame() %>%
    dplyr::filter(padj <= pAdj_cutoff) %>%
    dplyr::filter(abs(log2FoldChange) >= log2FC_Cutoff)

  write_tsv(
    as.data.frame(sresf) %>% tibble::rownames_to_column(var = "Gene"),
    file = file.path(out_dir, "04_GeneExpr_Shrunk_Filtered.tsv")
  )

  # Plot counts for the most significant gene (unshrunken alternative, unused)
  # plotCounts(desout, gene = which.min(res$padj), intgroup = "Group")

  # Customized count plot (log scale) for the top gene by shrunken padj
  d <- plotCounts(
    desout,
    gene = which.min(sres$padj),
    intgroup = "Group",
    returnData = TRUE
  )
  # Create short sample labels from rownames (last 8 characters)
  d$sample_short <- str_sub(rownames(d), -8)

  pdf(file.path(out_dir, "04_CountPlot_CaseVsCtl.pdf"), width = 6, height = 4)
  print(
    ggplot(d, aes(x = Group, y = count)) +
      geom_point(
        position = position_jitter(width = 0.1, height = 0),
        size = 3
      ) +
      geom_text_repel(
        aes(label = sample_short),
        size = 3,
        show.legend = FALSE
      ) +
      scale_y_log10(breaks = c(25, 100, 400)) +
      labs(
        title = paste0(dsname, " Counts (Top Gene)"),
        x = "Group",
        y = "Top Gene Count (log scale)"
      ) +
      theme_minimal()
  )
  dev.off()

  # Gene descriptions (interactive inspection only, not saved to file)
  mcols(sres)$description

  # Volcano with gene labels (subset to |log2FC| > 0.8 to reduce clutter)
  df <- tibble::rownames_to_column(as.data.frame(sres), var = "gene") %>%
    dplyr::rename(log2FC = log2FoldChange, p_value = padj)

  df$neg_log10_p_value <- -log10(df$p_value)

  df$significant <- ifelse(
    df$p_value < pAdj_cutoff & abs(df$log2FC) > log2FC_Cutoff,
    "Significant", "Not Significant"
  )

  pdf(file.path(out_dir, "04_VolcanoWithGenes_CaseVsCtl.pdf"), width = 6, height = 4)
  print(
    ggplot(
      df %>% dplyr::filter(abs(log2FC) > 0.8),
      aes(x = log2FC, y = neg_log10_p_value, label = gene, color = significant)
    ) +
      geom_point() +
      scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
      geom_hline(yintercept = -log10(pAdj_cutoff), linetype = "dashed", color = "red") +
      geom_vline(
        xintercept = c(-log2FC_Cutoff, log2FC_Cutoff),
        linetype = "dashed", color = "blue"
      ) +
      ggrepel::geom_text_repel() +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12)
      ) +
      labs(
        title = paste0(dsname, " Volcano Plot"),
        subtitle = paste0(
          "Significant if: log2FC cutoff >= ", log2FC_Cutoff, "; ",
          "padj <= ", pAdj_cutoff
        ),
        x = "log2 Fold Change", y = "-log10(p-value)"
      )
  )
  dev.off()

  ## EnhancedVolcano plot, with dynamic axis limits sized to the data
  pdf(file.path(out_dir, "04_EnhancedVolcanoPlot_CaseVsCtl.pdf"), width = 11, height = 6)

  # X-axis: log2FoldChange, padded to 1.1x the observed range
  x_rng <- range(sres$log2FoldChange, na.rm = TRUE)
  xlim_dyn <- c(1.1 * min(x_rng), 1.1 * max(x_rng))

  # Y-axis: -log10(padj), the value EnhancedVolcano uses internally
  y_rng <- range(-log10(sres$padj), na.rm = TRUE)
  ylim_dyn <- c(0, 1.1 * max(y_rng))

  plot(EnhancedVolcano(
    sres,
    lab            = rownames(res),
    x              = "log2FoldChange",
    y              = "padj",
    xlim           = xlim_dyn,
    ylim           = ylim_dyn,
    pCutoff        = pAdj_cutoff,
    FCcutoff       = log2FC_Cutoff,
    labSize        = 4,
    drawConnectors = TRUE,
    max.overlaps   = 60,
    title          = ""
  ))
  dev.off()

  invisible(NULL)
}

#############################
## Functions: Additional Plots
#############################

#' Additional visualization: heatmap (top 25 DE genes) and PCA (by Sex, by Age)
#'
#' Relies on `out_dir` and `sample_id_map_file` from the calling environment.
#'
#' @param desout DESeqDataSet after DESeq() has been run.
#' @param dsname Character dataset label used in plot titles.
#' @return Invisibly returns NULL. Produces plots.
additional_plots <- function(desout, dsname) {
  # Case vs Control results and rlog transform
  res <- DESeq2::results(desout, name = "Group_Case_vs_Control")
  rld <- rlog(desout)

  ## Heatmap: top 25 DE genes by padj
  # topGenes <- head(order(res$padj), 50)  # earlier version used top 50
  topGenes <- head(order(res$padj), 25)
  mat <- assay(rld)[topGenes, ]
  mat <- t(scale(t(mat)))  # Z-score normalization per gene

  pdf(file.path(out_dir, "04_HeatMap_CaseVsCtl_T25.pdf"), width = 16, height = 16)
  print(pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    annotation_col = as.data.frame(colData(desout)[, "Group", drop = FALSE])
  ))
  dev.off()

  ## PCA
  pcaData <- plotPCA(rld, intgroup = "Group", returnData = TRUE)
  pcaData$Sex <- colData(rld)$Sex[match(rownames(pcaData), colnames(rld))]
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData$Age <- colData(rld)$Age[match(rownames(pcaData), colnames(rld))]
  pcaData$label <- str_sub(rownames(pcaData), -8)

  # Map internal sample IDs to publication-friendly labels via external lookup
  y <- as.data.frame(rownames(pcaData)) %>%
    rename(SubjID_Display = `rownames(pcaData)`)
  z <- read_tsv(sample_id_map_file) %>%
    inner_join(as.data.frame(y)) %>%
    select(-SubjID_Display)
  pcaData$label <- z[[1]]

  # Labelled PCA plot (square, presentation-sized), colored by Group/Sex
  pdf(file.path(out_dir, "04_PCA_Dim_1_2_Plot_CaseVsCtl_Sex.pdf"), width = 16, height = 16)
  print(
    ggplot(pcaData, aes(PC1, PC2, color = Group, shape = Sex)) +
      geom_point(size = 3) +
      geom_text_repel(
        aes(label = label),
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      ggtitle(paste0(dsname, "TD Case vs Controls - PCA Dim1 vs Dim2 Plot - Sex")) +
      theme_bw()
  )
  dev.off()

  # Same PCA (Group/Sex), portrait layout with larger text/points for
  # publication figures. NOTE: originally written to the same filename as
  # the plot above (silently overwriting it) — given a distinct filename
  # here so both versions are kept.
  pdf(file.path(out_dir, "04_PCA_Dim_1_2_Plot_CaseVsCtl_Sex_Portrait.pdf"), width = 8, height = 11)
  print(
    ggplot(pcaData, aes(PC1, PC2, color = Group, shape = Sex)) +
      geom_point(size = 5) +
      geom_text_repel(
        aes(label = label),
        size = 8,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      ggtitle("") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold")
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 8)),
        shape = guide_legend(override.aes = list(size = 8))
      )
  )
  dev.off()

  pcaData <- pcaData %>%
    mutate(Age_bin = ifelse(Age < 18, "<18", ">=18"))

  # Labelled PCA plot, colored by Group/Age bin
  pdf(file.path(out_dir, "04_PCA_Dim_1_2_Plot_CaseVsCtl_Age.pdf"), width = 16, height = 16)
  print(
    ggplot(pcaData, aes(PC1, PC2, color = Group, shape = Age_bin)) +
      geom_point(size = 3) +
      geom_text_repel(
        aes(label = label),
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      ggtitle(paste0(dsname, "TD Case vs Controls - PCA Dim1 vs Dim2 Plot - Age Cutoff 18")) +
      theme_bw()
  )
  dev.off()

  invisible(NULL)
}

#############################
## Data Loading & Preprocessing
#############################

# Load consolidated gene expression counts (produced by 03_PreprocInput.r)
ecnt <- read_tsv(file.path(inp_dir, "gene_expr.tsv"))

# Load Case vs Control metadata (produced by 03_PreprocInput.r)
meta <- read_tsv(file.path(inp_dir, "CaseCtlMeta.tsv")) %>%
  mutate(
    Group = as.factor(Group),
    Sex = as.factor(Sex)
  ) %>%
  arrange(sample)
meta$Group <- relevel(meta$Group, ref = "Control")

# Convert long-format counts to a gene x sample count matrix.
# expected_count is scaled by 100 and truncated to an integer to preserve
# the original script's behavior (DESeq2 requires integer counts).
cnts <- ecnt %>%
  dplyr::select(gene_id, expected_count, sample) %>%
  dplyr::mutate(expected_count = as.integer(expected_count * 100)) %>%
  dplyr::group_by(sample) %>%
  tidyr::pivot_wider(names_from = sample, values_from = expected_count) %>%
  tibble::column_to_rownames("gene_id")

# Create DESeq dataset with full factorial design (Group, Sex, Age, and all
# their interactions)
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(cnts),
  colData = meta,
  design = ~ Group + Sex + Age +
    Group:Sex + Group:Age + Sex:Age + Group:Sex:Age
)

# Prefilter low-count genes: keep genes with >= 10 counts in at least 3 samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

# Fit the full model
dds_meta <- DESeq(dds)

# Report genes whose beta coefficients did not converge
# (e.g. "22 rows did not converge in beta" — see mcols(object)$betaConv;
# consider a larger `maxit` argument to nbinomWaldTest if this list is large)
rownames(dds_meta[(mcols(dds_meta)$betaConv == FALSE)])

# Case vs Control results and rlog transform (rld is saved to the workspace
# for reuse / re-plotting without recomputing)
res <- DESeq2::results(dds_meta, name = "Group_Case_vs_Control")
rld <- rlog(dds_meta)

# *** Note: only the Case vs Control coefficient is plotted here, though the
#     full interaction model above is used for fitting. ***
plot_deseq_results(dds_meta, "")

# apeglm citation (used for LFC shrinkage above):
# Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions
# for sequence count data: removing the noise and preserving large
# differences. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

#############################
## Reference output (from a prior run, preserved for comparison)
#############################
# out of 28304 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 16, 0.057%
# LFC < 0 (down)     : 23, 0.081%
# outliers [1]       : 4135, 15%
# low counts [2]     : 1701, 6%
# (mean count < 17)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# Warning message:
#   ggrepel: 54 unlabeled data points (too many overlaps). Consider increasing max.overlaps

additional_plots(dds_meta, "")

# Notes on transformation / typical warnings seen with this workflow:
# - rlog() may take a few minutes with 30+ samples; vst() is faster.
# - Heatmap/PCA use rlog-transformed data.
# - 1-2: ggrepel unlabeled points due to overlaps; consider raising max.overlaps.
# - 3: "one or more p-values is 0"; EnhancedVolcano substitutes the lowest
#      non-zero value for plotting.
# - 4-5: ggplot2 deprecation warnings (size -> linewidth) from upstream packages.
# - 6: pheatmap argument passthrough warning.

#############################
## Basic Summaries & Exports
#############################
# Ref: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

# Size factors per sample
size_factors <- as.data.frame(sizeFactors(dds_meta))
write_tsv(size_factors, file.path(out_dir, "size_factors.tsv"))

# Total raw counts per sample
raw_count_tot <- as.data.frame(colSums(counts(dds_meta)))
write_tsv(raw_count_tot, file.path(out_dir, "raw_count_tot.tsv"))

# Total normalized counts per sample
norm_count <- as.data.frame(colSums(counts(dds_meta, normalized = TRUE)))
write_tsv(norm_count, file.path(out_dir, "norm_count_tot.tsv"))

# Dispersion estimates
pdf(file.path(out_dir, "04_DispersionEstimates_Plot.pdf"), width = 6, height = 4)
plotDispEsts(dds_meta)
dev.off()

# Examine available coefficients
resultsNames(dds_meta)

# Case vs Control results from the full model
res_group <- results(dds_meta, name = "Group_Case_vs_Control")

# Significant genes table (padj <= pAdj_cutoff and |log2FC| >= log2FC_Cutoff)
sig_genes <- as.data.frame(
  res_group[!is.na(res_group$padj) &
    res_group$padj < pAdj_cutoff &
    !is.na(res_group$log2FoldChange) &
    abs(res_group$log2FoldChange) >= log2FC_Cutoff, ]
)

write_tsv(
  sig_genes %>% tibble::rownames_to_column(var = "Gene"),
  file.path(out_dir, "all_sig_genes.tsv")
)

#############################
## Violin & Box Plots for Selected Genes
#############################

# Genes to plot: the significant set from above
selected_genes <- sig_genes %>% row.names()
# Alternative source (from a previously filtered, shrunken results file):
# selected_genes <- read_tsv(file.path(out_dir,"04_GeneExpr_Shrunk_Filtered.tsv")) %>%
#   select(Gene) %>% as.character()

# Normalized counts
norm_counts <- counts(dds_meta, normalized = TRUE)

# Subset to selected genes
subset_counts <- norm_counts[selected_genes, ]

# Reshape to long format for ggplot2
df_long <- melt(subset_counts)
colnames(df_long) <- c("Gene", "Sample", "NormalizedCount")

# Attach Group from colData (preserving original script behavior of joining
# against `dds`, not `dds_meta`)
df_long$Group <- colData(dds)$Group[match(df_long$Sample, rownames(colData(dds)))]

# Violin plots, one facet per gene
pdf(file.path(out_dir, "04_SigGenes_ViolinPlots.pdf"), width = 12, height = 12)
print(
  ggplot(df_long, aes(x = Group, y = NormalizedCount, fill = Group)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    facet_wrap(~ Gene, scales = "free_y") +
    labs(
      title = "Violin Plots of Normalized Counts for Selected Genes",
      x = "Condition", y = "Normalized Counts"
    ) +
    theme_minimal()
)
dev.off()

# Box plots, one facet per gene
pdf(file.path(out_dir, "04_SigGenes_BoxPlots.pdf"), width = 12, height = 12)
print(
  ggplot(df_long, aes(x = Group, y = NormalizedCount, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    facet_wrap(~ Gene, scales = "free_y") +
    labs(
      title = "Box Plots of Normalized Counts for Selected Genes",
      x = "Condition", y = "Normalized Counts"
    ) +
    theme_minimal()
)
dev.off()

#############################
## Quick Lookup of Specific Genes
#############################
res_tbl <- as_tibble(res_group, rownames = "gene")
res_tbl %>%
  dplyr::filter(gene %in% sig_genes$Gene)

#############################
## Impact Summary per Coefficient
#############################
# For every model coefficient except the intercept, count how many genes are
# significant and report the mean |log2FC|, as a rough sense of which terms
# (Group, Sex, Age, interactions) drive the most variation.
coef_names <- resultsNames(dds_meta)
impact_list <- lapply(coef_names[-1], function(n) {
  res <- results(dds_meta, name = n)
  data.frame(
    Coefficient      = n,
    SignificantGenes = sum(res$padj < pAdj_cutoff, na.rm = TRUE),
    MeanLFC          = mean(abs(res$log2FoldChange), na.rm = TRUE)
  )
})
impact_summary <- do.call(rbind, impact_list)
impact_summary

#############################
## Wald Statistic vs Cook's Distance (Diagnostics)
#############################
# Genes with a high Wald statistic and low Cook's distance represent robust
# differential expression, while high Cook's distance indicates a signal
# driven by influential outlier samples. DESeq2 flags genes with extreme
# Cook's distance: their p-values are set to NA (small n) or down-weighted
# (larger n).

W <- res_group$stat
maxCooks <- apply(assays(dds_meta)[["cooks"]], 1, max)
idx <- !is.na(W)

pdf(file.path(out_dir, "04_WaldStatsVsCooksDist_Plot.pdf"), width = 6, height = 4)

# Base plot: rank of Wald statistic vs max Cook's distance, with an F-based
# reference cutoff line (p=3 parameters estimated per gene; see DESeq2 docs)
plot(
  rank(W[idx]), maxCooks[idx],
  xlab = "rank of Wald statistic",
  ylab = "maximum Cook's distance per gene",
  ylim = c(0, 5), cex = 0.4, col = rgb(0, 0, 0, 0.3)
)
m <- ncol(dds)
p <- 3
abline(h = qf(0.99, p, m - p))

## Plot variant: log-scaled y-axis, descending rank
cooks <- assays(dds_meta)[["cooks"]]
max_cooks <- apply(cooks, 1, max, na.rm = TRUE)
wald_stat <- res_group$stat
rank_wald <- rank(-wald_stat, ties.method = "first")  # descending

plot_data <- data.frame(
  Rank     = rank_wald,
  MaxCooks = max_cooks
)

plot(
  plot_data$Rank, plot_data$MaxCooks,
  log  = "y", pch = 20, col = "blue",
  xlab = "Rank of Wald Statistic",
  ylab = "Maximum Cook's Distance per Gene",
  main = "Cook's Distance vs Wald Statistic Rank"
)
abline(h = 1, col = "red", lty = 2)

## Cook's distance diagnostic, per sample, to identify influential samples

# Extract Cook's distance matrix (genes x samples)
cd <- assays(dds_meta)[["cooks"]]

# Rule-of-thumb cutoff: 4 / (number of samples)
cutoff <- 4 / ncol(cd)

cd_df <- as.data.frame(cd) %>%
  mutate(gene = rownames(cd)) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cooks") %>%
  group_by(sample) %>%
  summarise(frac_outliers = mean(cooks > cutoff, na.rm = TRUE))

print(
  ggplot(cd_df, aes(x = sample, y = frac_outliers)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
    theme_bw() +
    labs(
      x = "Sample",
      y = "Fraction of genes with high Cook's distance",
      title = "Sample-level Cook's Distance Impact"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
)
dev.off()

#############################
## Save Environment
#############################
# Persist the full workspace so downstream scripts can pick up where this
# one left off.
save(list = ls(), file = "./data/calcDEG.Rdata")

####
# Alternative approach: contribution of individual factors via LRT
####

# ---- Sex contribution: drop Sex + all Sex-involving interactions ----
dds_sex_red <- DESeq(dds, test = "LRT", reduced = ~ Group + Age + Group:Age)
res_sex_LRT <- results(dds_sex_red)
sex_frac <- mean(res_sex_LRT$padj < 0.05, na.rm = TRUE)  # fraction significant

# ---- Age contribution: drop Age + all Age-involving interactions ----
dds_age_red <- DESeq(dds, test = "LRT", reduced = ~ Group + Sex + Group:Sex)
res_age_LRT <- results(dds_age_red)
age_frac <- mean(res_age_LRT$padj < 0.05, na.rm = TRUE)  # fraction significant

sex_frac
age_frac
