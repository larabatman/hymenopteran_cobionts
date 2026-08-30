#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=busco_analysis
#SBATCH --output=logs/busco_analysis_%j.out
#SBATCH --error=logs/busco_analysis_%j.err


# Wrapper that launches the three Rscript to parse BUSCO full_table.tsv, make density blobplots and cross-tabulations with coverage classes
# Usage: 
# sbatch C2_run_busco_analysis.sh species

set -euo pipefail

SPECIES="$1"

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
SCRIPTS="scripts/exploratory_phase/02_busco"
# Input directory
BUSCO_HYM="results/${SPECIES}_stages/busco/${SPECIES}_busco_hymenoptera/run_hymenoptera_odb10/full_table.tsv"
BUSCO_ARTH="results/${SPECIES}_stages/busco/${SPECIES}_busco_arthropoda/run_arthropoda_odb10/full_table.tsv"
BUSCO_BACT="results/${SPECIES}_stages/busco/${SPECIES}_busco_bacteria/run_bacteria_odb10/full_table.tsv"
COV="results/${SPECIES}_stages/host_backbone/coverage_classification.tsv"
# Result directory
OUTDIR="results/${SPECIES}_stages/busco_analysis"
mkdir -p "$OUTDIR"

module load R/4.2.1-foss-2021a

Rscript "$SCRIPTS/C2a_busco_analysis.R" "$SPECIES" "$BUSCO_HYM" "$BUSCO_ARTH" "$BUSCO_BACT" "$COV" "$OUTDIR"

# Results from that parsing:
BUSCO_SUMMARY="${OUTDIR}/busco_per_contig_summary.tsv"
# Blobplots
Rscript "$SCRIPTS/C2b_busco_plots.R" "$SPECIES" "$BUSCO_SUMMARY" "$OUTDIR"
# Cross tables
Rscript "$SCRIPTS/C2c_busco_anchors.R" "$SPECIES" "$BUSCO_SUMMARY" "$OUTDIR"