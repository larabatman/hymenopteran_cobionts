#!/usr/bin/env bash
# Wrapper for merge_nuwt_regions.R
# Loads bedtools and R, then runs the merge.
#
# Usage:
# bash merge_nuwt_regions.sh <kept_nuwt_tsv> <fragments_tsv> <out_prefix>

set -euo pipefail

SCRIPT_DIR="/data/users/lland/cobionts/scripts/nuwt/supergroups"

module load BEDTools/2.30.0-GCC-10.3.0
module load R/4.2.1-foss-2021a

Rscript "${SCRIPT_DIR}/merge_nuwt_regions.R" "$@"