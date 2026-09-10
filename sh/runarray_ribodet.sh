#!/usr/bin/bash
#SBATCH --partition=main
#SBATCH --job-name=ribodet
#SBATCH --array=1-54%18            # 54 tasks total, max 18 running concurrently
#SBATCH --requeue
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB
#SBATCH --time=2:55:00
#SBATCH -o %j.out                  # stdout: <SLURM_JOB_ID>.out
#SBATCH -e %j.err                  # stderr: <SLURM_JOB_ID>.err
#SBATCH --export=ALL               # forward current environment to compute nodes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
###############################################################################
# runarray_ribodet.sh
#
# Purpose:
#   SLURM array-job wrapper that fans out ribodetector runs across many
#   samples in parallel on the Rutgers Amarel cluster.
#
#   This script does NOT call ribodetector directly. Instead, it expects a
#   set of pre-generated per-task shell scripts named x1, x2, ... x54 in
#   ./work (one per SLURM array index), each of which contains
#   the actual `runribodetector.sh <SAMPLE>` command (or equivalent) for a
#   single sample. This "array of scripts" pattern lets many independent
#   commands be queued as one SLURM array job.
#
#   Array index -> input script mapping:
#     SLURM_ARRAY_TASK_ID = 1  -->  ./work/x1
#     SLURM_ARRAY_TASK_ID = 2  -->  ./work/x2
#     ...                            ...
#     SLURM_ARRAY_TASK_ID = 54 -->  ./work/x54
#
# Usage:
#   1. Generate ./work/x1 .. x54, one per sample.
#   2. Submit with: sbatch runarray_ribodet.sh
#
# Output:
#   Each task gets its own working directory:
#     ./work/<JOBID>/<ARRAY_TASK_ID>/
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
cd ./work
mkdir -p "./work/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}"
cd "./work/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}"

# --- Run this array task's per-sample command file --------------------------
# Reads and executes ./work/x<ARRAY_TASK_ID> as a bash script,
# redirecting its stdout to <ARRAY_TASK_ID>.output in the per-task directory.
srun /usr/bin/bash < "./work/x${SLURM_ARRAY_TASK_ID}" > "${SLURM_ARRAY_TASK_ID}.output"

echo "End Date/Time              = $(date)"
