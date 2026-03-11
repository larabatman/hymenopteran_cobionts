#!/usr/bin/env bash
#SBATCH --job-name=hifiasm_assembly
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --output=logs/hifiasm_assembly_%j.out
#SBATCH --error=logs/hifiasm_assembly_%j.err

# Stage A: Assembly
# A2: hifiasm assembly
# Assemble raw HiFi reads, exploiting Hi-C contact information for phasing when available 
# Supports:
#   - HiFi-only assembly (ASM_MODE=bp)
#   - HiFi + Hi-C assembly (ASM_MODE=hic)

# 1) Setup
set -euo pipefail

# Variables
SPECIES="$1"
ASM_MODE="$2"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

OUTDIR="assemblies/hifiasm_clean/${SPECIES}"
PREFIX="${OUTDIR}/asm"

mkdir -p "$OUTDIR" logs

HIFI=(reads/pacbio_hifi/${SPECIES}/*.fastq.gz)

THREADS=${SLURM_CPUS_PER_TASK}

# Modules
module load hifiasm/0.16.1-GCCcore-10.3.0

# Logs and checks
echo "[INFO] Species: $SPECIES"
echo "[INFO] Assembly mode: $ASM_MODE"

for f in "${HIFI[@]}"; do
    [[ -s "$f" ]] || { echo "ERROR: missing HiFi file $f"; exit 1; }
done

# 2a) Run hifiasm, HiFi only mode
if [[ "$ASM_MODE" == "bp" ]]; then

    echo "[INFO] Running HiFi-only assembly"

    hifiasm \
      -o "$PREFIX.bp" \
      -t "$THREADS" \
      "${HIFI[@]}"

    GFA="${PREFIX}.bp.p_ctg.gfa"
    FASTA="${PREFIX}.bp.p_ctg.fasta"

fi

# 2b) Run hifiasm, HiFi + Hi-C mode
if [[ "$ASM_MODE" == "hic" ]]; then

    R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
    R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"

    [[ -s "$R1" ]] || { echo "ERROR: missing Hi-C R1"; exit 1; }
    [[ -s "$R2" ]] || { echo "ERROR: missing Hi-C R2"; exit 1; }

    echo "[INFO] Running HiFi + Hi-C assembly"

    hifiasm \
      -o "$PREFIX.hic" \
      -t "$THREADS" \
      --h1 "$R1" \
      --h2 "$R2" \
      "${HIFI[@]}"

    GFA="${PREFIX}.hic.p_ctg.gfa"
    FASTA="${PREFIX}.hic.p_ctg.fasta"

fi

# Convert GFA to FASTA
echo "[INFO] Converting GFA to FASTA"

[[ -s "$GFA" ]] || { echo "ERROR: GFA missing: $GFA"; exit 1; }
# awk command to select sequences S lines from GFA
awk '/^S/{print ">"$2"\n"$3}' "$GFA" > "$FASTA"

echo "[INFO] Assembly finished"
echo "[INFO] FASTA: $FASTA"