#!/usr/bin/bash
###############################################################################
# globinCount.sh
#
# Purpose:
#   QC utility that scans all RSEM output (recursively, from the current
#   directory) for counts belonging to a fixed list of genes of interest -
#   primarily globin genes, used to check for globin transcript contamination
#   in blood-derived RNA-seq samples, plus two housekeeping genes (ACTB,
#   GAPDH) as a normalization/sanity reference. See the `globinGenes` file
#   in this folder for the exact gene list.
#
# Usage:
#   Run from a directory tree containing one or more RSEM
#   "counts.genes.results" files (e.g. the top-level output directory from
#   runrspl.sh), with no arguments:
#
#       ./globinCount.sh
#
#   (Optionally redirect to a file: ./globinCount.sh > globin_counts.tsv)
#
# Requirements:
#   - The gene list file must exist at: ./globinGenes
#     (one gene symbol per line; see globinGenes in this folder). Update this
#     hardcoded path if the gene list is moved or this script is run by a
#     different user.
#
# Output:
#   For every counts.genes.results file found under the current directory,
#   prints any lines whose gene symbol/ID exactly matches a whole word in
#   the gene list, prefixed with the source file path
#   (grep -H behavior), e.g.:
#
#       ./sampleA/counts/counts.genes.results:HBB&ENSG00000244734.4  123  45.6 ...
###############################################################################

set -euo pipefail

# -H            : always print the matching file name (useful with `find -exec`)
# -F            : treat patterns in the gene list as fixed strings, not regex
# -w            : match whole words only (avoid partial matches, e.g. "HBB"
#                 inside a longer identifier)
# -f <file>     : read match patterns (gene symbols) from the given file
find . -name counts.genes.results \
    -exec grep -H -Fwf ./globinGenes {} \;
