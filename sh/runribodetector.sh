#!/usr/bin/bash
###############################################################################
# runribodetector.sh
#
# Purpose:
#   Remove ribosomal RNA (rRNA) reads from a paired-end RNA-seq sample using
#   ribodetector (CPU version). This is typically the first QC/filtering step
#   applied to raw FASTQ files before adapter trimming and alignment.
#
# Usage:
#   ./runribodetector.sh <SAMPLE>
#
#   <SAMPLE>  Sample name/prefix. Input FASTQ files are expected at:
#             ${DATADIR}/<SAMPLE>_R1.fq.gz
#             ${DATADIR}/<SAMPLE>_R2.fq.gz
#
# Requirements:
#   - The `ribodetector` conda environment must be activated before running
#     this script, e.g.:
#         conda activate ribodetector
#
# Outputs (written to ${OUTDIR}):
#   <SAMPLE>.nonrrna.1.fq.gz / <SAMPLE>.nonrrna.2.fq.gz  - reads WITHOUT rRNA
#                                                           (used downstream)
#   <SAMPLE>.rrna.1.fq.gz    / <SAMPLE>.rrna.2.fq.gz      - reads flagged as rRNA
#
# Log:
#   ${LOGDIR}/<SAMPLE>.rdet.log
###############################################################################

set -euo pipefail

# --- Path configuration -----------------------------------------------------
DATADIR=/scratch/ks1437/tsdata     # location of raw input FASTQ files
OUTDIR=/scratch/ks1437/tsout       # location for filtered FASTQ output
LOGDIR=/scratch/ks1437/tslog       # location for ribodetector log files

# --- Input argument ----------------------------------------------------------
SMPL="$1"   # sample name, used to build input/output file paths

# --- Run ribodetector ---------------------------------------------------------
# -t                 : number of threads
# -l                 : read length (bp)
# -i                 : input paired-end FASTQ files (R1 R2)
# -e both            : evaluate both reads of the pair for rRNA content
# -o                 : output FASTQ files with rRNA reads removed (R1 R2)
# --chunk_size       : number of reads processed per batch
# -r                 : output FASTQ files containing only the detected rRNA reads
# --log              : path to write the run log
ribodetector_cpu -t 20 \
    -l 100 \
    -i "${DATADIR}/${SMPL}_R1.fq.gz" "${DATADIR}/${SMPL}_R2.fq.gz" \
    -e both \
    -o "${OUTDIR}/${SMPL}.nonrrna.1.fq.gz" "${OUTDIR}/${SMPL}.nonrrna.2.fq.gz" \
    --chunk_size 3072 \
    -r "${OUTDIR}/${SMPL}.rrna.1.fq.gz" "${OUTDIR}/${SMPL}.rrna.2.fq.gz" \
    --log "${LOGDIR}/${SMPL}.rdet.log"
