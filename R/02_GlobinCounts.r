#' ==============================================================================
#' 02_GlobinCounts.r
#'
#' Purpose:
#'   Check the percentage of globin gene expression (HB* genes) in RSEM
#'   quantification output, as a QC metric for globin depletion efficiency.
#'   Expected value after successful depletion is ~0.5% +/- 0.6%.
#'   Reference: https://doi.org/10.1186/s12864-020-07304-4
#'
#' Inputs:
#'   - ./data/Globins.Rdata
#'       Must contain a data frame `df` with (at minimum) the columns:
#'         sample   - sample identifier (character/factor)
#'         gene_id  - gene identifier (character), globin genes start with "HB"
#'         TPM      - transcripts per million (numeric) from RSEM
#'   - Alternatively, uncomment the block below to read `df` directly from
#'     ./input/globin_counts.txt (raw RSEM TPM output) instead of the .Rdata
#'     file.
#'
#' Outputs:
#'   - ./results/02_Globins/globinsPct.txt   Per-sample globin % (TSV)
#'   - ./results/02_Globins/GlobinPct.pdf    Scatter plot of globin % by sample
#'   - ./data/Globins.Rdata                  Updated workspace (all objects)
#'
#' Dependencies:
#'   - tidyverse (includes ggplot2, dplyr, readr, etc.)
#' ==============================================================================

library(tidyverse)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

inp_dir <- "./input/"
out_dir <- "./results/02_Globins/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------------------------

# Primary source: previously saved workspace containing `df`.
load(file = "./data/Globins.Rdata")

# Alternative source: read directly from raw RSEM TPM output instead of the
# cached .Rdata. Uncomment if `df` needs to be regenerated from source
# (e.g. after re-running RSEM on unfiltered RNA-seq data).
# df <- read_tsv(paste0(inp_dir, "globin_counts.txt"))

stopifnot(
  "df is missing required columns" =
    all(c("sample", "gene_id", "TPM") %in% names(df))
)

# ------------------------------------------------------------------------------
# 3. Compute globin percentage per sample
# ------------------------------------------------------------------------------
# Globin genes are identified by gene_id starting with "HB" (e.g. HBA1, HBB).
# TPM values sum to 1,000,000 per sample, so dividing the summed globin TPM
# by 10,000 converts it to a percentage of total expression.

globinPct <- df %>%
  filter(startsWith(gene_id, "HB")) %>%
  group_by(sample) %>%
  summarize(globpct = sum(TPM)) %>%
  mutate(globpct = globpct / 10000) %>%
  arrange(-globpct)

write_tsv(globinPct, paste0(out_dir, "globinsPct.txt"))

# ------------------------------------------------------------------------------
# 4. Plot globin percentage per sample
# ------------------------------------------------------------------------------
# Samples are ordered along the x-axis from highest to lowest globin %, so
# poorly-depleted samples stand out on the left.

glob_plot <- ggplot(
  globinPct,
  aes(x = reorder(sample, -globpct), y = globpct)
) +
  geom_point(color = "blue", size = 2) +
  labs(title = "", x = "Sample", y = "Globins %") +
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

glob_plot

ggsave(
  plot = glob_plot,
  filename = paste0(out_dir, "GlobinPct.pdf"),
  width = 8,
  height = 4
)

# ------------------------------------------------------------------------------
# 5. Save workspace
# ------------------------------------------------------------------------------
# Persist the full workspace so downstream scripts can pick up where this
# one left off.

save(list = ls(), file = "./data/Globins.Rdata")
