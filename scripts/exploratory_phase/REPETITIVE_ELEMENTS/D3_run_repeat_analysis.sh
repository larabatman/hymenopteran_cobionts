#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=repeat_analysis
#SBATCH --output=logs/repeat_analysis_%j.out
#SBATCH --error=logs/repeat_analysis_%j.err

# =============================================================================
# Stage D3: Unified repeat element analysis
# Parses EDTA GFF3 + RepeatMasker + TRF, produces 3 blobplots and summaries.
#
# Requires: D1_edta.sh + D2_rm_trf.sh completed + B1b coverage classification
# Usage: sbatch D3_run_repeat_analysis.sh <species> <asm_mode>
# =============================================================================

set -euo pipefail

SPECIES="$1"
ASM_MODE="$2"

WORKDIR=/data/projects/p2025-0083_mining_cobionts
cd "$WORKDIR"

GENOME="assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
PREFIX="$(basename "$GENOME")"

# ── Inputs ──
EDTA_GFF3="results/${SPECIES}_stages/edta/${PREFIX}.mod.EDTA.TEanno.gff3"
RM_OUT="results/${SPECIES}_stages/repeat_masking/repeatmasker/${PREFIX}.out"
TRF_DAT="results/${SPECIES}_stages/repeat_masking/trf/${PREFIX}.2.7.7.80.10.50.500.dat"
COV="results/${SPECIES}_stages/host_backbone/coverage_classification.tsv"

OUTDIR="results/${SPECIES}_stages/repeat_analysis"
mkdir -p "$OUTDIR"

# ── Checks ──
for f in "$GENOME" "$EDTA_GFF3" "$RM_OUT" "$TRF_DAT" "$COV"; do
    [[ -s "$f" ]] || { echo "[ERROR] Missing: $f" >&2; exit 1; }
done

echo "=============================="
echo "D3 Unified repeat analysis"
echo "Species: $SPECIES | Mode: $ASM_MODE"
echo "Output:  $OUTDIR"
echo "=============================="

module load R/4.2.1-foss-2021a

Rscript scripts/stages/exploratory_phase/REPETITIVE_ELEMENTS/D3_repeat_analysis.R \
    "$SPECIES" \
    "$EDTA_GFF3" \
    "$RM_OUT" \
    "$TRF_DAT" \
    "$GENOME" \
    "$COV" \
    "$OUTDIR"

echo "[OK] Done. Results in: $OUTDIR"