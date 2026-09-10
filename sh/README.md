# RNA-seq Processing Pipeline

Scripts for processing paired-end bulk RNA-seq samples on the Rutgers Amarel
SLURM cluster: rRNA filtering, adapter trimming, alignment/quantification
(two independent methods — STAR+RSEM and Salmon), and a globin-transcript QC
check for blood-derived samples.

## Pipeline overview

```
raw FASTQ (*_R1.fq.gz / *_R2.fq.gz)
        │
        ▼
 [1] runribodetector.sh   ── remove rRNA reads (ribodetector)
        │  (submitted at scale via runarray_ribodet.sh)
        ▼
 *.nonrrna.1/2.fq.gz
        │
        ├──────────────────────────────┐
        ▼                              ▼
 [2] runrspl.sh                  [2'] salmoncnt.sh
     Trim Galore → STAR → RSEM        Salmon (pseudo-alignment)
     (submitted at scale via          quant.sf per sample
      runarray_rspl.sh)
        │
        ▼
 counts.genes.results (per sample)
        │
        ▼
 [3] globinCount.sh  ── QC: check globin / housekeeping gene counts
                         (uses gene list in `globinGenes`)
```

STAR+RSEM and Salmon are two **independent** quantification routes run on
the same rRNA-filtered input — not sequential steps. Use whichever
(or both, for cross-checking) your downstream analysis needs.

## Files in this folder

| File | Type | Purpose |
|---|---|---|
| `runribodetector.sh` | per-sample script | Removes rRNA reads from one sample's paired FASTQ files using `ribodetector_cpu`. |
| `runarray_ribodet.sh` | SLURM array wrapper | Distributes `runribodetector.sh`-style commands across many samples (54, throttled to 18 concurrent) as one SLURM array job. |
| `runrspl.sh` | per-sample script | Trims adapters (Trim Galore), aligns to the T2T genome (STAR), and quantifies expression (RSEM) for one sample. |
| `runarray_rspl.sh` | SLURM array wrapper | Distributes `runrspl.sh`-style commands across many samples (54, throttled to 12 concurrent) as one SLURM array job. |
| `salmoncnt.sh` | per-sample script | Alternative quantification of one sample against the hg38 transcriptome using Salmon (no STAR alignment needed). |
| `globinCount.sh` | QC utility | Recursively searches for RSEM `counts.genes.results` files and extracts counts for the genes listed in `globinGenes`. |
| `globinGenes` | data file | Plain list of gene symbols (one per line) used by `globinCount.sh`: hemoglobin/globin genes plus `ACTB`/`GAPDH` as housekeeping references. **Not a script — do not add comments to it**, since every line is read as a literal `grep` match pattern. |

## The SLURM array-of-scripts pattern

`runarray_ribodet.sh` and `runarray_rspl.sh` don't call the per-sample
scripts directly. Instead, each expects 54 pre-generated command files named
`x1`, `x2`, ... `x54` in `./work` — one per sample/array index —
where each `xN` file contains the actual command to run for that sample
(e.g. a line invoking `runribodetector.sh <SAMPLE>` or `runrspl.sh <SAMPLE>`).
SLURM's `--array=1-54` then runs each `xN` file as its own task, each in its
own scratch working directory, with concurrency capped by the `%N` suffix
(`%18`, `%12`) to avoid overloading the cluster/filesystem.

**To use these wrappers:**
1. Generate `./work/x1` ... `x54`, one file per sample,
   each containing the command for that sample.
2. Submit with `sbatch runarray_ribodet.sh` (or `runarray_rspl.sh`).

## Typical usage

```bash
# 1. Filter rRNA reads (per sample, or at scale via the array wrapper)
./runribodetector.sh SAMPLE01

# 2a. Trim + align + quantify with STAR/RSEM (per sample, or at scale)
./runrspl.sh SAMPLE01

# 2b. Or/also quantify with Salmon
./salmoncnt.sh SAMPLE01

# 3. QC: check globin/housekeeping gene counts across all samples
cd /path/to/rsem/output/root
./globinCount.sh > globin_counts.tsv
```

## Requirements

- **Cluster:** Rutgers Amarel (SLURM).
- **Conda environment:** `ribodetector` (must be activated manually before
  running `runribodetector.sh` — it is not activated inside the script).
- **Environment modules:** `TrimGalor`, `STAR/2.7.5a`, `samtools/1.3.1`,
  `RSEM/1.3.3-yc759`, `salmon` (loaded automatically within the relevant
  scripts via `module load`).
- **Reference data:**
  - T2T genome + STAR index + RSEM reference: `./DataSets/rsem/`
    (GTF: `GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf`)
  - Salmon hg38 index + GENCODE v38 GTF: `./DataSets/salmon/`
  
## Notes

- **ribodetector logs:**
  - rRNA and non-rRNA read counts by sample have to be extracted from ribodetector log  files into a tab separated file that is required as input for R script that checks 	rRNA filtering per sample

