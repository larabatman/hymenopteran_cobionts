#!/usr/bin/env bash
#SBATCH --job-name=identify_dominant_cov_mode
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/identify_dominant_cov_mode_%j.out
#SBATCH --error=logs/identify_dominant_cov_mode_%j.err

# Wrapper for coverage modelling in R script B2_identify_dominant_cov_mode.R
# Uses gc_cov.tsv and outputs the classification and tables
# Usage:
# sbatch species

set -euo pipefail

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

SPECIES="$1"
# gc_cov.tsv with contig, length, GC and mean depth for each contig
INPUT="results/${SPECIES}_stages/gc_cov/gc_cov.tsv"
OUTDIR="results/${SPECIES}_stages/host_backbone"
mkdir -p "$OUTDIR"

# Load R
module load R/4.2.1-foss-2021a

# Lauch Rscript
Rscript scripts/exploratory_phase/COVERAGE/B2_identify_dominant_cov_mode.R "${SPECIES}" "${INPUT}" "${OUTDIR}"