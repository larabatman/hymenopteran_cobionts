#!/usr/bin/env bash
#SBATCH --job-name=hifiasm_assembly
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=164G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hifiasm_assembly_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hifiasm_assembly_%j.err

# =============================================================================
# Stage A2: hifiasm assembly
# Assemble filtered HiFi reads, exploiting Hi-C contact information for phasing when available
# Supports:
#   - HiFi-only assembly (ASM_MODE=bp)
#   - HiFi + Hi-C assembly (ASM_MODE=hic)
# It converts the primary GFA outputed by hifiasm to FASTA and produces basic assemblies statistics
# Usage: sbatch A2_hifiasm.sh <species> <asm_mode>

set -euo pipefail

# Setup
# Arguments
SPECIES="$1"
ASM_MODE="$2"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

OUTDIR="assemblies/hifiasm/${SPECIES}"
PREFIX="${OUTDIR}/asm"

mkdir -p "$OUTDIR" logs

HIFI=(results/${SPECIES}_stages/read_qc/hifi_qc/hifi_filtered/*.fastq.gz)

# Hi-C cleaned reads from A1
HIC_R1="reads/hic_clean/${SPECIES}/hic_R1.clean.fastq.gz"
HIC_R2="reads/hic_clean/${SPECIES}/hic_R2.clean.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Load modules
module purge
module load hifiasm/0.16.1-GCCcore-10.3.0
module load SeqKit/2.6.1

# Logs and checks on inputs
echo "[INFO] Species: $SPECIES"
echo "[INFO] Assembly mode: $ASM_MODE"
echo "[INFO] Threads: $THREADS"

for f in "${HIFI[@]}"; do
    [[ -s "$f" ]] || { echo "[ERROR] Missing or empty HiFi file: $f"; exit 1; }
done

# HiFi-only assembly
if [[ "$ASM_MODE" == "bp" ]]; then
    echo "[INFO] Running HiFi-only assembly"
    hifiasm \
        -o "${PREFIX}" \
        -t "$THREADS" \
        "${HIFI[@]}"

    GFA="${PREFIX}.bp.p_ctg.gfa"
    FASTA="${PREFIX}.bp.p_ctg.fasta"
fi

# HiFi + Hi-C assembly
if [[ "$ASM_MODE" == "hic" ]]; then
    [[ -s "$HIC_R1" ]] || { echo "[ERROR] Missing or empty Hi-C R1: $HIC_R1"; exit 1; }
    [[ -s "$HIC_R2" ]] || { echo "[ERROR] Missing or empty Hi-C R2: $HIC_R2"; exit 1; }

    echo "[INFO] Running HiFi + Hi-C assembly"
    hifiasm \
        -o "${PREFIX}" \
        -t "$THREADS" \
        --h1 "$HIC_R1" \
        --h2 "$HIC_R2" \
        "${HIFI[@]}"

    GFA="${PREFIX}.hic.p_ctg.gfa"
    FASTA="${PREFIX}.hic.p_ctg.fasta"
fi

# Convert GFA to FASTA
# awk command grabs lines contianing the sequences (S) in the GFA file
echo "[INFO] Converting GFA to FASTA"
[[ -s "$GFA" ]] || { echo "[ERROR] GFA missing or empty: $GFA"; exit 1; }
awk '/^S/{print ">"$2"\n"$3}' "$GFA" > "$FASTA"

# Assembly stats
echo "[INFO] Assembly stats"
seqkit stats -T "$FASTA" > "${OUTDIR}/assembly_basic_stats.tsv"
cat "${OUTDIR}/assembly_stats.tsv"

echo "[INFO] Assembly finished: $FASTA"