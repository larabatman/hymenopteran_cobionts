#!/usr/bin/env bash
#SBATCH --job-name=identify_dominant_cov_mode
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/identify_dominant_cov_mode_%j.out
#SBATCH --error=logs/identify_dominant_cov_mode_%j.err

# Stage B: Orthogonal Evidence Layers
# Coverage-based partition for host backbone definition
# B2: Define host backbone with coverage as primary discriminator 
# Inputs: gc_cov.tsv
# Output: coverage_classification.tsv, host_backbone.tsv and coverage_backbone_summary.tsv
# Model host-like coverage weighted by contig length through MAD 

# 1) Setup
set -euo pipefail

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

# Inputs
SPECIES="$1"
INPUT="results/${SPECIES}_stages/gc_cov/gc_cov.tsv"
OUTDIR="results/${SPECIES}_stages/host_backbone"
mkdir -p "$OUTDIR"

# Logs
echo "[INFO] Species: $SPECIES"
echo "[INFO] Input file: $INPUT"
echo "[INFO] Output dir: $OUTDIR"

# Load R
module load R/4.2.1-foss-2021a

# 2) Launch MAD COV analysis
Rscript scripts/stages/B1b_identify_dominant_cov_mode.R "${SPECIES}" "${INPUT}" "${OUTDIR}"