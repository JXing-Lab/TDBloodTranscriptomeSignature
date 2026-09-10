#' ==============================================================================
#' 01_RiboStats.r
#'
#' Purpose:
#'   Assess ribosomal RNA (rRNA) contamination in RNA-seq samples by plotting
#'   the percentage of non-rRNA reads per sample.
#'
#' Inputs:
#'   - ./data/RiboStats.Rdata
#'       Must contain a data frame `df` with (at minimum) the columns:
#'         Sample        - sample identifier (character/factor)
#'         Percent_rRNA  - percentage of reads mapping to rRNA (numeric, 0-100)
#'   - Alternatively, read `df` directly from ./input/RiboCounts.txt instead of .Rdata file.
#'    RiboCounts.txt is a tsv file created using the log file from runribodetector.sh
#'    containing the following fields: Sample	rRNA_reads	non_rRNA	rRNA%
#'
#' Outputs:
#'   - ./results/01_RiboStats/RiboStats.pdf   Scatter plot of non-rRNA % by sample
#'   - ./results/01_RiboStats/RiboStats.tsv   Tab-separated copy of `df`
#'   - ./data/RiboStats.Rdata                 Updated workspace (all objects)
#'
#' Dependencies:
#'   - tidyverse (dplyr, ggplot2, readr, etc.)
#'   - xlsx      (only needed if reading directly from the Excel source file)
#' ==============================================================================

library(tidyverse)
#library(xlsx)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

inp_dir <- "./input/"
out_dir <- "./results/01_RiboStats/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------------------------

# Primary source: previously saved workspace containing `df`.
# load(file = "./data/RiboStats.Rdata")

# Alternative source: read directly from the raw Excel stats file instead of
# relying on the cached .Rdata. Uncomment if `df` needs to be regenerated
# from source (e.g. after a re-run of the pipeline).
df <- read_tsv(paste0(inp_dir, "RiboCounts.txt"))

stopifnot(
  "df is missing required columns" =
    all(c("Sample", "Percent_rRNA") %in% names(df))
)

# ------------------------------------------------------------------------------
# 3. Plot non-rRNA percentage per sample
# ------------------------------------------------------------------------------
# Samples are ordered along the x-axis by their rRNA percentage so that the
# gradient from cleanest to most contaminated sample is easy to read.

ribo_plot <- ggplot(
  df,
  aes(x = reorder(Sample, Percent_rRNA), y = 100 - Percent_rRNA)
) +
  geom_point(color = "blue", size = 2) +
  labs(
    title = "",
    x = "Sample",
    y = "non_rRNA reads percentage"
  ) +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ribo_plot

ggsave(
  plot = ribo_plot,
  filename = paste0(out_dir, "RiboStats.pdf"),
  width = 8,
  height = 4
)

# ------------------------------------------------------------------------------
# 4. Export table + save workspace
# ------------------------------------------------------------------------------

write_tsv(df, paste0(out_dir, "RiboStats.tsv"))

# Persist the full workspace so downstream scripts can pick up where this
# one left off.
save(list = ls(), file = "./data/RiboStats.Rdata")
