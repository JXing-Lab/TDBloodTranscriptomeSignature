# RNA-seq Case/Control Differential Expression & Biomarker Pipeline

RNA-seq analysis pipeline for a Case vs Control study (globin- and rRNA-depleted
whole blood/tissue RNA-seq, STAR/RSEM quantification): QC, differential
expression (DESeq2), and biomarker classification (ElasticNet and SVM, both
with leave-one-out cross-validation).

## Pipeline overview

Scripts are numbered and **must be run in order** — each one reads outputs
produced by an earlier script.

| # | Script | Purpose |
|---|--------|---------|
| 01 | [`01_RiboStats.r`](01_RiboStats.r) | QC: plot residual ribosomal RNA (rRNA) % per sample |
| 02 | [`02_GlobinCounts.r`](02_GlobinCounts.r) | QC: compute and plot globin gene expression % per sample (depletion efficiency) |
| 03 | [`03_PreprocInput.R`](03_PreprocInput.R) | Build sample metadata (Group/Sex/Age) and consolidate per-sample RSEM gene/isoform counts into single expression tables |
| 04 | [`04_CalcDEG.R`](04_CalcDEG.R) | Differential expression analysis with DESeq2 (Case vs Control, adjusting for Sex/Age); volcano, heatmap, PCA, diagnostics |
| 05 | [`05_LOOCV_ElasticNet.R`](05_LOOCV_ElasticNet.R) | LOOCV elastic net classifier on candidate biomarker genes; ROC/AUC |
| 06 | [`06_SVM_LOOCV.R`](06_SVM_LOOCV.R) | LOOCV SVM (radial kernel, caret) classifier on the same gene set; ROC/AUC |

Scripts 01 and 02 are independent QC checks and can be run any time after raw
counts exist; 03 must run before 04; 04 must run before 05 and 06. 05 and 06
are independent of each other (both consume 04's output).

```
        01_RiboStats.r   02_GlobinCounts.r (QC — independent)
                        |
                        v
              03_PreprocInput.R
                        |
                        v
                04_CalcDEG.R
                    /       \
                   v         v
   05_LOOCV_ElasticNet.R   06_SVM_LOOCV.R
```

## Requirements

R (developed against a recent 4.x release) with the following packages:

```r
install.packages(c(
  "tidyverse", "xlsx", "gt", "ggrepel", "pheatmap", "reshape2",
  "caret", "pROC", "glmnet", "knitr"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
```

| Script | Key packages |
|---|---|
| 01 | tidyverse, xlsx |
| 02 | tidyverse, ggplot2 |
| 03 | xlsx, ggplot2, gt, tidyverse |
| 04 | DESeq2, ggplot2, ggrepel, stringr, tidyverse, pheatmap, EnhancedVolcano, reshape2 |
| 05 | DESeq2, pROC, caret, knitr, tidyverse, glmnet |
| 06 | DESeq2, caret, pROC, knitr, tidyverse |

## Directory layout

Scripts assume they are run from a project root with this structure (created
automatically where noted):

```
.
├── input/                        # Excel/TSV inputs you supply
│   ├── CaseCtlMeta.tsv
│                                 # optional, only if re-reading raw stats
├── data/                         # per-sample RSEM output + cached .Rdata workspaces
│   ├── genecnts/*.txt
│   ├── isocnts/*.txt
│   ├── RiboStats.Rdata
│   ├── Globins.Rdata
│   ├── preproc.Rdata
│   ├── calcDEG.Rdata
│   ├── ElasticNet_LOOCV.Rdata
│   └── SVM_LOOCV.Rdata
└── results/                      # created automatically by each script
    ├── 01_RiboStats/
    ├── 02_Globins/
    ├── 03_PreprocInput/
    ├── 04_CalcDEG/deg_<pAdj>_<log2FC>_<date>/
    ├── 05_ElasticNet_LOOCV/ (or _RELAXED)
    └── 06_SVM_CARET_LOOCV/ (or _RELAXED)
```

Each script loads a `.Rdata` workspace at the start (if one exists) and saves
its own updated workspace at the end, so intermediate state persists between
runs.

## Usage

Run in order from the project root:

```r
source("./R/01_RiboStats.r")
source("./R/02_GlobinCounts.r")
source("./R/03_PreprocInput.R")
source("./R/04_CalcDEG.R")
source("./R/05_LOOCV_ElasticNet.R")
source("./R/06_SVM_LOOCV.R")
```
