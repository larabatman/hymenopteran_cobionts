#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=summarize_cobiont_analysis
#SBATCH --output=logs/summarize_cobiont_analysis_%j.out
#SBATCH --error=logs/summarize_cobiont_analysis_%j.err

# This script is a wrapper to launch the R script 03d_summarize_cobiont_analysis.R

set -euo pipefail

# Load R
module load R/4.2.1-foss-2021a

# Project root
WORKDIR="/data/projects/p2025-0083_a_cross-species_pipeline_for_mining_cobionts_in_hymenopteran_genomes"
cd "$WORKDIR"

# Run QC script
Rscript scripts/03d_integrate_GC_cov_Kraken.R