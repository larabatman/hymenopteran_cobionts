#!/usr/bin/env bash
#SBATCH --job-name=hic_to_cram
#SBATCH --partition=pibu_el8
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hic_to_cram_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hic_to_cram_%j.err

# Convert clean HiC FASTQs to unaligne CRAM
# Usage:
# sbatch A2b_hic_to_cram.sh species

set -euo pipefail

SPECIES="$1"

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

mkdir -p logs

R1="reads/hic_clean/${SPECIES}/hic_R1.clean.fastq.gz"
R2="reads/hic_clean/${SPECIES}/hic_R2.clean.fastq.gz"

PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$PROC_DIR"
OUT="${PROC_DIR}/hic.cram"

# Modules
module purge
module load SAMtools

# Convert. first define the read groups line which needs to be stamp onto every read
# ID to $SPECIES, like sample name SN and platform unit PU. Then library Lb to hic,  and platform PL to ILLUMINA
RG=$(printf '@RG\tID:%s\tSM:%s\tLB:hic\tPL:ILLUMINA\tPU:%s' "$SPECIES" "$SPECIES" "$SPECIES")

# samtools import: converts FASTQ files to unaligned CRAM format
# -o with .cram extention and the -r with the RG stamp above
samtools import \
    -1 "$R1" \
    -2 "$R2" \
    -o "$OUT" \
    -@ $((SLURM_CPUS_PER_TASK - 1)) \
    -r "$RG"