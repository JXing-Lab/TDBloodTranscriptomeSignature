#!/usr/bin/bash
###############################################################################
# salmoncnt.sh
#
# Purpose:
#   Alternative, pseudo-alignment-based quantification of a single RNA-seq
#   sample using Salmon against the hg38 transcriptome, run directly on the
#   ribosomal-RNA-filtered FASTQ files. This is an independent path to
#   RSEM/STAR (see runrspl.sh) - Salmon does not require a prior STAR
#   alignment step, since it maps reads to transcripts itself.
#
# Usage:
#   ./salmoncnt.sh <SAMPLE>
#
#   <SAMPLE>  Sample name/prefix. Input FASTQ files are expected at:
#             ${DATADIR}/<SAMPLE>.nonrrna.1.fq.gz
#             ${DATADIR}/<SAMPLE>.nonrrna.2.fq.gz
#             (i.e. the "non-rRNA" output of runribodetector.sh)
#
# Output:
#   ${WORKDIR}/<SAMPLE>/            - Salmon quantification output
#                                      (quant.sf, lib_format_counts.json, etc.)
#   ${WORKDIR}/<SAMPLE>/logs/       - SLURM/log output directory
#
# NOTE: module load can be sensitive to `set -e`/`set -u` on some module
# systems, so strict mode is intentionally not enabled in this script.
###############################################################################

module load salmon

# --- Path configuration -----------------------------------------------------
WORKDIR='/scratch/ks1437/rnaseqout'
DATADIR='/scratch/ks1437/RNASeq'

INDEXDIR='/projectsp/f_jx76_1/DataSets/salmon/hg38_idx'          # pre-built Salmon index
GTFFILE='/projectsp/f_jx76_1/DataSets/salmon/gencode.v38.annotation.gtf'  # hg38 gene annotation (for --geneMap)

# --- Input argument ----------------------------------------------------------
SAMPLE="$1"
# SAMPLE=600_01_TS0108_5003   # example sample name (left here for reference)

R1="${DATADIR}/${SAMPLE}.nonrrna.1.fq.gz"
R2="${DATADIR}/${SAMPLE}.nonrrna.2.fq.gz"

MAPDIR="${WORKDIR}/${SAMPLE}"
LOGDIR="${WORKDIR}/${SAMPLE}/logs"

# --- Echo run parameters (useful for debugging in the log) ------------------
echo "datadir= ""$DATADIR"
echo "workdir= ""$WORKDIR"
echo "logdir= ""$LOGDIR"
echo "gtf= ""$GTFFILE"
echo "sample= ""$SAMPLE"
echo "R1= ""$R1"
echo "R2= ""$R2"

# --- Create output directories -----------------------------------------------
mkdir -p "${MAPDIR}"
mkdir -p "${LOGDIR}"   # SLURM reports will go here

# --- Run Salmon quantification -----------------------------------------------
# -i                    : path to pre-built Salmon index
# --libType A           : auto-detect library type (strandedness)
# -1 / -2               : paired-end FASTQ input (R1 / R2)
# --threads             : number of threads
# --validateMappings    : use Salmon's more accurate selective-alignment mode
# --geneMap             : GTF used to summarize transcript-level to gene-level counts
# --gcBias              : correct for GC-content bias
# --seqBias             : correct for sequence-specific bias
# --minAssignedFrags    : minimum number of assigned fragments required (else warn/fail)
# --output              : output directory for this sample
salmon quant \
    -i "$INDEXDIR" \
    --libType A \
    -1 "${R1}" \
    -2 "${R2}" \
    --threads 64 \
    --validateMappings \
    --geneMap "$GTFFILE" \
    --gcBias \
    --seqBias \
    --minAssignedFrags 1 \
    --output "$MAPDIR"
