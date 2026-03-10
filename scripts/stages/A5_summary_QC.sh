#!/usr/bin/env bash
#SBATCH --job-name=summarize_QC
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/summarize_QC_%j.out
#SBATCH --error=logs/summarize_QC_%j.err

# Stage A: Assembly
# A5: Summarize assembly QC
# Output: species-level and global QC summary tables
# This script parses stage A assembly outputs and prepares a standardized QC table that is then consumed by the R summarization script.

set -euo pipefail

# 1) Setup
SPECIES="$1"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

# 2) Launch R summary
module load R/4.2.1-foss-2021a

echo "[INFO] Launching QC summarization"

Rscript scripts/stages/A5_summarize_QC.R "$SPECIES" "$WORKDIR"

echo "[INFO] QC summary completed for $SPECIES"