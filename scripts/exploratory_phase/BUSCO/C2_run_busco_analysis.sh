#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=busco_analysis
#SBATCH --output=logs/busco_analysis_%j.out
#SBATCH --error=logs/busco_analysis_%j.err

# =============================================================================
# Stage C2: Unified BUSCO contig analysis
# Contig-level collapse + 6 blobplots + all summary tables
#
# Requires: C1_busco.sh completed + B1b coverage classification
# Usage: sbatch C2_run_busco_analysis.sh <species>
# =============================================================================

set -euo pipefail

SPECIES="$1"

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

# ── Inputs ──
BUSCO_HYM="results/${SPECIES}_stages/busco/${SPECIES}_busco_hymenoptera/run_hymenoptera_odb10/full_table.tsv"
BUSCO_ARTH="results/${SPECIES}_stages/busco/${SPECIES}_busco_arthropoda/run_arthropoda_odb10/full_table.tsv"
BUSCO_BACT="results/${SPECIES}_stages/busco/${SPECIES}_busco_bacteria/run_bacteria_odb10/full_table.tsv"
COV="results/${SPECIES}_stages/host_backbone/coverage_classification.tsv"

OUTDIR="results/${SPECIES}_stages/busco_analysis"
mkdir -p "$OUTDIR"

# ── Checks ──
for f in "$BUSCO_HYM" "$BUSCO_ARTH" "$BUSCO_BACT" "$COV"; do
    [[ -s "$f" ]] || { echo "[ERROR] Missing: $f" >&2; exit 1; }
done

echo "=============================="
echo "C2 Unified BUSCO analysis"
echo "Species: $SPECIES"
echo "Output:  $OUTDIR"
echo "=============================="

module load R/4.2.1-foss-2021a

Rscript scripts/stages/C2_busco_analysis.R \
    "$SPECIES" \
    "$BUSCO_HYM" \
    "$BUSCO_ARTH" \
    "$BUSCO_BACT" \
    "$COV" \
    "$OUTDIR"

echo "[OK] Done. Results in: $OUTDIR"