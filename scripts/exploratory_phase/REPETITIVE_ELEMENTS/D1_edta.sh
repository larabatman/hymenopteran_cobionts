#!/usr/bin/env bash
#SBATCH --job-name=EDTA
#SBATCH --partition=pibu_el8
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --output=logs/edta_%j.out
#SBATCH --error=logs/edta_%j.err

# =============================================================================
# Stage D1: EDTA — de novo TE annotation
# Runs the EDTA pipeline for whole-genome transposable element annotation.
# Produces a classified TE library and genome-wide GFF3 annotation.
#
# Usage: sbatch D1_edta.sh <species> <asm_mode>
# =============================================================================

set -euo pipefail

SPECIES="$1"
ASM_MODE="$2"

WORKDIR=/data/projects/p2025-0083_mining_cobionts
GENOME=$WORKDIR/assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta
IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif

RUN_DIR="$WORKDIR/results/${SPECIES}_stages/edta"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

[ -s "$GENOME" ] || { echo "ERROR: genome not found: $GENOME" >&2; exit 1; }
[ -s "$IMG" ]    || { echo "ERROR: EDTA image not found: $IMG" >&2; exit 1; }

echo "[INFO] Species: $SPECIES | Mode: $ASM_MODE"
echo "[INFO] Running EDTA..."

apptainer exec --bind "$WORKDIR","/data/courses/assembly-annotation-course" \
    "$IMG" EDTA.pl \
    --genome "$GENOME" \
    --species others \
    --step all \
    --sensitive 1 \
    --anno 1 \
    --threads "${SLURM_CPUS_PER_TASK}"

echo "[OK] EDTA completed."