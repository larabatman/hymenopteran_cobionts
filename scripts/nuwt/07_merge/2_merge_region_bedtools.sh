#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=merge_nuwt_regions
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/merge_nuwt_regions_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/merge_nuwt_regions_%j.err

# Loads bedtools and launches 2_merge_nuwt_regions.R
# Usage:
# sbash 2_merge_nuwt_regions_bedtools.sh

set -euo pipefail
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
cd "$PROJECT_ROOT"
NUWT_DIR="cobionts/nuwtscan"
SCRIPT_DIR="/scripts/nuwtt/07_merge/"
# Fragments to keep given by 4_og_mito_ank.R
KEPT_NUWT="${NUWT_DIR}/placement/aggregate/hym_final_supergroup_nuwt.tsv"
# Fragment coordinates
FRAGMENTS="${NUWT_DIR}/nuwt_fragments.tsv"
OUT_PREFIX="${NUWT_DIR}/placement/aggregate/hym"

module load BEDTools/2.30.0-GCC-10.3.0
module load R/4.2.1-foss-2021a

Rscript "${SCRIPT_DIR}/2_merge_nuwt_regions.R" "$KEPT_NUWT" "$FRAGMENTS" "$OUT_PREFIX"