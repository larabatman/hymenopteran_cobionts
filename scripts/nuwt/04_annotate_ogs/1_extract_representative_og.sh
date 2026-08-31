#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=extract_representative_og
#SBATCH --output=logs/extract_representative_og_%j.out
#SBATCH --error=logs/extract_representative_og_%j.err

# Launches python file to select a representative protein for each OG
# Usage:
# sbatch 1_extract_representative_og.sh

set -euo pipefail

# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
cd "$PROJECT_ROOT"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
HMM_DIR="${NUWT_DIR}/hmm_database"
# Deposit with the alignments
ALN_DIR="${HMM_DIR}/NUWTs_Release_v3.0/Data/4_WolbachiaOrthoAlignments"
# Output
OUT_FASTA="${HMM_DIR}/og_representatives.fna"
# Script directory
SCRIPT_DIR="${PROJECT_ROOT}/scripts/nuwt/04_annotation"
EXTRACT_PY="${SCRIPT_DIR}/extract_og_representatives.py"

module load Python/3.9.5-GCCcore-10.3.0

# For each OG, extract one representative sequence
python3 "${EXTRACT_PY}" "${ALN_DIR}" "${OUT_FASTA}"