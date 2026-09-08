#!/usr/bin/bash
###############################################################################
# runrspl.sh  (RNA-Seq PipeLine: trim -> align -> quantify)
#
# Purpose:
#   Run the "downstream" portion of the RNA-seq processing pipeline for a
#   single sample, on top of ribosomal-RNA-filtered reads produced by
#   runribodetector.sh:
#     1. Trim adapters / low-quality bases with Trim Galore
#     2. Align trimmed reads to the T2T reference genome with STAR
#     3. Quantify gene/transcript expression with RSEM
#
# Usage:
#   ./runrspl.sh <SAMPLE>
#
#   <SAMPLE>  Sample name/prefix. Input FASTQ files are expected at:
#             ${DATADIR}/<SAMPLE>.nonrrna.1.fq.gz
#             ${DATADIR}/<SAMPLE>.nonrrna.2.fq.gz
#             (i.e. the "non-rRNA" output of runribodetector.sh)
#
# Notes:
#   - This script performs `module purge` / `module load` calls, which can
#     misbehave under `set -e` / `set -u` on some cluster module systems, so
#     strict-mode is intentionally NOT enabled here (unlike the other scripts
#     in this folder). Validate paths/arguments carefully before submitting.
#   - Each of the three stages below was originally a separate script
#     (trimgalor.sh, staraln.sh, rsemcount.sh); they have been kept together
#     here but the section banners preserve that boundary for reference.
###############################################################################

export SAMPLE="$1"

# --- Path configuration -----------------------------------------------------
export DATADIR="/scratch/ks1437/RNASeq/"
export GENOMEDIR="/projectsp/f_jx76_1/DataSets/rsem/"
export RSEMREF="${GENOMEDIR}/T2T"
export WORKDIR="/home/ks1437/jgd/TSP3/RNASeq/out/"
export OUTDIR="${WORKDIR}/${SAMPLE}/"

###############################################################################
# Stage 1: trimgalor.sh - Trim adapters and filter low-quality reads
###############################################################################
module purge
module load TrimGalor

# -q               : Phred quality threshold for trimming (Q20)
# --phred33        : input quality encoding
# --fastqc         : run FastQC on the trimmed output
# --stringency     : minimum overlap required to trim an adapter (lenient: 1bp)
# -e               : max allowed error rate for adapter matching
# --gzip           : gzip-compress trimmed output
# --basename       : prefix for output file names
# --cores          : number of threads
# --paired         : paired-end mode, followed by R1 and R2 input files
trim_galore \
    --output_dir "${OUTDIR}" \
    -q 20 \
    --phred33 \
    --fastqc \
    --stringency 1 \
    -e 0.1 \
    --gzip \
    --basename "${SAMPLE}" \
    --cores 16 \
    --paired \
    "${DATADIR}/${SAMPLE}.nonrrna.1.fq.gz" "${DATADIR}/${SAMPLE}.nonrrna.2.fq.gz"

###############################################################################
# Stage 2: staraln.sh - Align trimmed reads to the T2T genome with STAR
###############################################################################
module purge
module load STAR/2.7.5a
module load samtools/1.3.1

# Trim Galore output naming convention: <basename>_val_1.fq.gz / _val_2.fq.gz
export fq1="${OUTDIR}/${SAMPLE}_val_1.fq.gz"
export fq2="${OUTDIR}/${SAMPLE}_val_2.fq.gz"   # NOTE: added missing "/" before
                                                #       ${SAMPLE} (original had
                                                #       "${OUTDIR}${SAMPLE}...")
export GTFFILE="${GENOMEDIR}/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf"
export PFX="${OUTDIR}/${SAMPLE}_"

# --runMode alignReads          : standard alignment run
# --genomeDir                   : path to pre-built STAR genome index
# --outFileNamePrefix           : prefix for all STAR output files
# --readFilesIn                 : R1 R2 trimmed FASTQ input
# --readFilesCommand zcat       : decompress gzipped FASTQ on the fly
# --outSAMtype BAM SortedByCoordinate : output coordinate-sorted BAM
# --outSAMunmapped Within       : keep unmapped reads in the output BAM
# --outReadsUnmapped Fastx      : also write unmapped reads as FASTQ
# --seedSearchStartLmax 30      : seed search parameter (tuned for read length)
# --quantMode TranscriptomeSAM GeneCounts : produce transcriptome BAM (for RSEM)
#                                            and per-gene STAR counts
# --sjdbGTFfile                 : gene annotation used for splice junctions
STAR \
    --runThreadN 16 \
    --runMode alignReads \
    --genomeDir "${GENOMEDIR}" \
    --outFileNamePrefix "${PFX}" \
    --readFilesIn "$fq1" "$fq2" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outReadsUnmapped Fastx \
    --seedSearchStartLmax 30 \
    --quantMode TranscriptomeSAM GeneCounts \
    --sjdbGTFfile "$GTFFILE"

###############################################################################
# Stage 3: rsemcount.sh - Quantify gene expression for the sample with RSEM
###############################################################################
export RSEMREF="${GENOMEDIR}/T2T"
export COUNTDIR="${OUTDIR}/counts"
export LOGDIR="${COUNTDIR}/logs"

module purge
module load RSEM/1.3.3-yc759

# --alignments                  : input is already an aligned BAM
#                                  (STAR's Aligned.toTranscriptome.out.bam)
# --paired-end                  : paired-end mode
# --append-names                : append gene symbols to gene IDs in output
# --calc-pme                    : calculate posterior mean estimates
# --calc-ci                     : calculate 95% credibility intervals
rsem-calculate-expression \
    -p 16 \
    --alignments \
    --paired-end "${OUTDIR}/${SAMPLE}_Aligned.toTranscriptome.out.bam" \
    --append-names \
    --calc-pme \
    --calc-ci \
    "${RSEMREF}" "${COUNTDIR}"
