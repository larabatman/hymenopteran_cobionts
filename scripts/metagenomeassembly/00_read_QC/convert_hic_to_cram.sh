#!/usr/bin/env bash
#SBATCH --job-name=hic_to_cram
#SBATCH --partition=pibu_el8
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hic_to_cram_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hic_to_cram_%j.err

# A2b: Convert clean Hi-C paired FASTQs to unaligned CRAM
#
# Input:
# reads/hic_clean/<species>/hic_R1.clean.fastq.gz   (from A1b)
# reads/hic_clean/<species>/hic_R2.clean.fastq.gz
# Output: reads/processed/<species>/hic.cram  
#
# Usage:  sbatch A2b_hic_to_cram.sh <species>

set -euo pipefail

SPECIES="$1"

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

mkdir -p logs

R1="reads/hic_clean/${SPECIES}/hic_R1.clean.fastq.gz"
R2="reads/hic_clean/${SPECIES}/hic_R2.clean.fastq.gz"

[[ -f "$R1" ]] || { echo "[ERROR] Missing R1: $R1"; exit 1; }
[[ -f "$R2" ]] || { echo "[ERROR] Missing R2: $R2"; exit 1; }

PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$PROC_DIR"
OUT="${PROC_DIR}/hic.cram"

echo "[INFO] ${SPECIES}"
echo " R1  = ${R1}"
echo " R2  = ${R2}"
echo " OUT = ${OUT}"

# Modules
module purge
module load SAMtools

# Convert
# Define the RG read group line to stamp onto every read
# ID set to $SPECIES
# SM sample name set to $SPECIES
# LB library set to hic
# PL platform set to ILLUMINA
# PU platform unit set to $SPECIES
RG=$(printf '@RG\tID:%s\tSM:%s\tLB:hic\tPL:ILLUMINA\tPU:%s' \
     "$SPECIES" "$SPECIES" "$SPECIES")

# samtools import converts FASTQ files to ulaigned SAM/BAM/CRAM format
# -1 & -2 specify the first and second read files for paired-end data
# -o sets the output file, set to .cram extention
# -@ controls the number of additional threads to use
# -r to stamp the reads
# unmapped CRAM: no reference genome
samtools import \
    -1 "$R1" \
    -2 "$R2" \
    -o "$OUT" \
    -@ $((SLURM_CPUS_PER_TASK - 1)) \
    -r "$RG"

# Check
READ_COUNT=$(samtools view -c "$OUT")
echo "[INFO] Read count: ${READ_COUNT}"
[[ "$READ_COUNT" -gt 0 ]] || { echo "[ERROR] CRAM has no reads"; exit 1; }

echo ""
echo "[OK] Hi-C CRAM: ${OUT}"
ls -lh "$OUT"