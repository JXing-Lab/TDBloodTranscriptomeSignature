#!/usr/bin/bash
#SBATCH --partition=p_ch977_1,nonpre,p_am2820_1,p_ak1833_1,genetics_1,main
#SBATCH --job-name=rseqpl
#SBATCH --array=1-54%12            # 54 tasks total, max 12 running concurrently
#SBATCH --requeue
#SBATCH --cpus-per-task=40
#SBATCH --mem=120GB
#SBATCH --time=15:30:00
#SBATCH -o %j.out                  # stdout: <SLURM_JOB_ID>.out
#SBATCH -e %j.err                  # stderr: <SLURM_JOB_ID>.err
#SBATCH --export=ALL               # forward current environment to compute nodes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ks1437@rutgers.edu
###############################################################################
# runarray_rspl.sh
#
# Purpose:
#   SLURM array-job wrapper that fans out the trim/align/quantify pipeline
#   (runrspl.sh: Trim Galore -> STAR -> RSEM) across many samples in
#   parallel on the Rutgers Amarel cluster.
#
#   Like runarray_ribodet.sh, this script does NOT call runrspl.sh directly.
#   It expects a set of pre-generated per-task shell scripts named
#   x1, x2, ... x54 in /scratch/ks1437/work (one per SLURM array index),
#   each of which invokes `runrspl.sh <SAMPLE>` (or equivalent) for a single
#   sample. This "array of scripts" pattern lets many independent, heavier
#   pipeline runs be queued as one SLURM array job while limiting how many
#   run concurrently (here, 12 at a time via the %12 throttle).
#
#   Array index -> input script mapping:
#     SLURM_ARRAY_TASK_ID = 1  -->  /scratch/ks1437/work/x1
#     SLURM_ARRAY_TASK_ID = 2  -->  /scratch/ks1437/work/x2
#     ...                            ...
#     SLURM_ARRAY_TASK_ID = 54 -->  /scratch/ks1437/work/x54
#
# Usage:
#   1. Generate /scratch/ks1437/work/x1 .. x54, one per sample (run AFTER
#      ribodetector filtering has completed for all samples).
#   2. Submit with: sbatch runarray_rspl.sh
#
# Output:
#   Each task gets its own working directory:
#     /scratch/ks1437/work/<JOBID>/<ARRAY_TASK_ID>/
#   containing <ARRAY_TASK_ID>.output (stdout of the per-task script) plus
#   the standard SLURM %j.out / %j.err files.
#
# NOTE: module purge/use calls can be sensitive to `set -e`/`set -u` on some
# module systems, so strict mode is intentionally not enabled in this script.
###############################################################################

module purge
module use /projects/community/modulefiles

# --- Job/environment info (useful for debugging in the .out log) -----------
echo "Start Date/Time              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo ""
echo "SLURM JOB ID                   = $SLURM_JOB_ID"
echo "SLURM ARRAY TASK ID            = $SLURM_ARRAY_TASK_ID"

# --- Set up a per-task scratch working directory ----------------------------
cd /scratch/ks1437/work
mkdir -p "/scratch/ks1437/work/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}"
cd "/scratch/ks1437/work/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}"

# --- Run this array task's per-sample command file --------------------------
# Reads and executes /scratch/ks1437/work/x<ARRAY_TASK_ID> as a bash script,
# redirecting its stdout to <ARRAY_TASK_ID>.output in the per-task directory.
srun /usr/bin/bash < "/scratch/ks1437/work/x${SLURM_ARRAY_TASK_ID}" > "${SLURM_ARRAY_TASK_ID}.output"

echo "End Date/Time              = $(date)"
