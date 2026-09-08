#' ==============================================================================
#' 03_PreprocInput.r
#'
#' Purpose:
#'   Preprocess input data and consolidate samples for differential expression
#'   analysis: build sample metadata (group/sex/age), summarize sample
#'   characteristics, and merge per-sample RSEM gene/isoform count files into
#'   single consolidated expression tables.
#'
#' Last updated: 4/16/2026
#'
#' Notes:
#'   1. 24 cases and 24 controls have sufficient reads and are used here.
#'   2. Library prep and processing excluded globins and ribosomal RNAs.
#'   3. Ribodetect was used to further filter ribosomal RNA reads.
#'   4. Globin % was confirmed to be < 5% and globin genes were excluded from
#'      analysis.
#'
#' Inputs:
#'   - ./input/CaseCtlMeta.tsv
#'       Case/control list merged with sex and age attributes.
#'   - ./data/genecnts/*.txt
#'       Per-sample RSEM gene-level count files.
#'   - ./data/isocnts/*.txt
#'       Per-sample RSEM isoform-level count files.
#'   - ./data/preproc.Rdata
#'       Previously saved workspace (if re-running from a cached state).
#'
#' Outputs:
#'   - ./results/03_PreprocInput/01_AgeByCaseCtl_Histogram.pdf
#'       Stacked histogram of age by case/control group.
#'   - ./results/03_PreprocInput/01_SexGroup.html
#'       Summary table of sex breakdown by group.
#'   - ./input/gene_expr.tsv
#'       Consolidated gene-level counts across all samples.
#'   - ./input/isoform_expr.tsv
#'       Consolidated isoform-level counts across all samples.
#'   - ./data/preproc.Rdata
#'       Updated workspace (all objects).
#'
#' Dependencies:
#'   - xlsx        (reading the sample attributes Excel file)
#'   - ggplot2     (age histogram)
#'   - gt          (sex/group summary table)
#'   - tidyverse   (dplyr, readr, stringr, tidyr, etc.)
#' ==============================================================================

library("xlsx")
library(ggplot2)
library(gt)
library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

load(file = "./data/preproc.Rdata")

inp_dir <- "./input/"
out_dir <- "./results/03_PreprocInput/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Functions
# ------------------------------------------------------------------------------

#' Read and combine all per-sample RSEM count files in a directory.
#'
#' Each ".txt" file in `datadir` is assumed to be a single sample's RSEM
#' output. The sample ID is derived from the file path: after splitting on
#' "/", the 4th path element is taken and truncated to its first 18
#' characters (this assumes a consistent directory depth / naming scheme,
#' e.g. ".../genecnts/600-01-TS0234-5002.genes.results.txt").
#'
#' @param datadir Path to a directory of per-sample ".txt" count files.
#' @return A data frame combining all samples' counts, with a `sample` column.
get_rnaseq_data <- function(datadir) {
  ecnt <- data.frame()

  for (fn in list.files(path = datadir, pattern = ".txt$", full.names = TRUE)) {
    df <- read_tsv(fn) %>%
      mutate(sample = substr(str_split(fn, "/")[[1]][4], 1, 18))
    if (nrow(ecnt) == 0) {
      ecnt <- df
    } else {
      ecnt <- rbind(ecnt, df)
    }
  }
  return(ecnt)
}

# ------------------------------------------------------------------------------
# 3. Build sample metadata (Group + Sex + Age)
# ------------------------------------------------------------------------------

meta <- read_tsv(inp_dir, "CaseCtlMeta.tsv") %>%
  mutate(Sex = as.factor(Sex),
         Age = as.numeric(Age),
         Group = as.factor(Group))

# ------------------------------------------------------------------------------
# 4. Sample characteristics: stats / plots
# ------------------------------------------------------------------------------

# --- 4a. Stacked histogram of age by case/control -----------------------------
data <- meta %>%
  select(Group, Age)

pdf(file.path(out_dir, "01_AgeByCaseCtl_Histogram.pdf"), width = 6, height = 4)
ggplot(data, aes(x = Age, fill = Group)) +
  geom_histogram(binwidth = 5, position = "stack", color = "black") +
  labs(
    title = "Stacked Histogram of Age by Case/Control",
    x = "Age",
    y = "Count"
  ) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  theme_minimal()
dev.off()

# --- 4b. Summary table of sex breakdown by group -------------------------------
summary_table <- meta %>%
  select(Group, Sex) %>%
  group_by(Group, Sex) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Sex, values_from = Count, values_fill = 0) %>%
  mutate(Total = Male + Female)

gt(summary_table) %>%
  tab_header(title = "Breakdown of Sex by Group") %>%
  fmt_number(columns = c("Male", "Female", "Total"), decimals = 0) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels(everything())
  ) %>%
  gtsave(file.path(out_dir, "01_SexGroup.html"))

# ------------------------------------------------------------------------------
# 5. Consolidate RNA-seq gene & isoform counts across samples
# ------------------------------------------------------------------------------
# Combine per-sample RSEM output into single wide tables (one row per
# gene/isoform per sample), for use in downstream differential expression
# analysis. Group/Sex/Age are dropped here since these files hold only counts.

RNASeqGeneCnts <- "./data/genecnts/"
RNASeqIsoCnts <- "./data/isocnts/"

ecnt <- get_rnaseq_data(RNASeqGeneCnts) %>%
  inner_join(meta) %>%
  select(-c(Group, Sex, Age))
write_tsv(ecnt, "./input/gene_expr.tsv")

icnt <- get_rnaseq_data(RNASeqIsoCnts) %>%
  inner_join(meta) %>%
  select(-c(Group, Sex, Age))
write_tsv(icnt, "./input/isoform_expr.tsv")

# ------------------------------------------------------------------------------
# 6. Save workspace
# ------------------------------------------------------------------------------
# Persist the full workspace so downstream scripts can pick up where this
# one left off.

save(list = ls(), file = "./data/preproc.Rdata")
