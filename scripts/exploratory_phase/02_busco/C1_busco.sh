#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=busco
#SBATCH --output=logs/busco_%j.out
#SBATCH --error=logs/busco_%j.err

# =============================================================================
# Stage C1: BUSCO marker gene detection
# Run BUSCO against three lineage databases on the hifiasm assembly
#   - hymenoptera_odb10 : strong host anchors
#   - arthropoda_odb10  : supportive host anchors
#   - bacteria_odb10    : cobiont / anti-host signal
#
# Usage: sbatch C1_busco.sh <species> <asm_mode>

set -euo pipefail

# Setup
# Arguments
SPECIES="$1"
ASM_MODE="$2"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

ASM="assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
BUSCO_OUT="results/${SPECIES}_stages/busco"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$BUSCO_OUT" logs
# Sanity check assembly existence
[[ -s "$ASM" ]] || { echo "[ERROR] Assembly not found: $ASM" >&2; exit 1; }

# Modules
module purge
module load BUSCO/5.4.2-foss-2021a

echo "[INFO] Species: $SPECIES | Mode: $ASM_MODE"
echo "[INFO] Assembly: $ASM"

# Hymenoptera
# genome mode
echo "[INFO] Running BUSCO — hymenoptera_odb10"
busco \
    --in "$ASM" \
    --out "${SPECIES}_busco_hymenoptera" \
    --out_path "$BUSCO_OUT" \
    --mode genome \
    --lineage hymenoptera_odb10 \
    --cpu "$THREADS" \
    -f

# Arthropoda
echo "[INFO] Running BUSCO — arthropoda_odb10"
busco \
    --in "$ASM" \
    --out "${SPECIES}_busco_arthropoda" \
    --out_path "$BUSCO_OUT" \
    --mode genome \
    --lineage arthropoda_odb10 \
    --cpu "$THREADS" \
    -f

# Bacteria
echo "[INFO] Running BUSCO — bacteria_odb10"
busco \
    --in "$ASM" \
    --out "${SPECIES}_busco_bacteria" \
    --out_path "$BUSCO_OUT" \
    --mode genome \
    --lineage bacteria_odb10 \
    --cpu "$THREADS" \
    -f

echo "[OK] BUSCO complete for $SPECIES"