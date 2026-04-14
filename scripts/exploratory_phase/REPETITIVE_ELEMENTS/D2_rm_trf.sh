#!/usr/bin/env bash
#SBATCH --job-name=RM_TRF
#SBATCH --partition=pibu_el8
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=logs/rm_trf_%j.out
#SBATCH --error=logs/rm_trf_%j.err

# =============================================================================
# Stage D2: RepeatMasker + TRF
# RepeatMasker uses the EDTA-derived TE library for classified repeat masking.
# TRF detects tandem repeats and microsatellites.
# Together they capture simple repeats and low-complexity regions not covered
# by the TE library alone.
#
# Requires: D1_edta.sh completed (needs .TElib.fa)
# Usage: sbatch D2_rm_trf.sh <species> <asm_mode>
# =============================================================================

set -euo pipefail

SPECIES="$1"
ASM_MODE="$2"

WORKDIR=/data/projects/p2025-0083_mining_cobionts
GENOME=$WORKDIR/assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta
PREFIX="$(basename "$GENOME")"

EDTA_DIR=$WORKDIR/results/${SPECIES}_stages/edta
TELIB="${EDTA_DIR}/${PREFIX}.mod.EDTA.TElib.fa"

OUTDIR=$WORKDIR/results/${SPECIES}_stages/repeat_masking
RM_DIR="${OUTDIR}/repeatmasker"
TRF_DIR="${OUTDIR}/trf"
mkdir -p "$RM_DIR" "$TRF_DIR"

[ -s "$GENOME" ] || { echo "ERROR: assembly not found: $GENOME" >&2; exit 1; }
[ -s "$TELIB" ]  || { echo "ERROR: EDTA TE library not found: $TELIB" >&2; exit 1; }

module load RepeatMasker/4.1.5-foss-2021a
module load TRF/4.09.1-GCC-10.3.0

# ── RepeatMasker ──
echo "[INFO] Running RepeatMasker on $PREFIX"
cp "$GENOME" "$RM_DIR/"
cd "$RM_DIR"

RepeatMasker \
    -lib "$TELIB" \
    -pa "${SLURM_CPUS_PER_TASK}" \
    -gff \
    -xsmall \
    -no_is \
    -dir "$RM_DIR" \
    "${RM_DIR}/${PREFIX}"

echo "[OK] RepeatMasker done."

# ── TRF ──
echo "[INFO] Running TRF on $PREFIX"
cd "$TRF_DIR"

trf "$GENOME" 2 7 7 80 10 50 500 -d -h || true

echo "[OK] TRF done."
echo "[INFO] RepeatMasker output: $RM_DIR"
echo "[INFO] TRF output:          $TRF_DIR"