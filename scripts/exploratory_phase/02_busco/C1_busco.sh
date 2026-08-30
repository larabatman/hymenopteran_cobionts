#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=busco
#SBATCH --output=logs/busco_%j.out
#SBATCH --error=logs/busco_%j.err

# Run BUSCO for conserved gene detection in three lineages: hymenoptera_odb10, arthropofa_odb10, bacteria_odb10
# Usage:
# sbatch C1_busco.sh species ASM_MODE==bp or ASM_MODE==hic

set -euo pipefail

# Arguments
SPECIES="$1"
ASM_MODE="$2"

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
# Assembly files, fasta
ASM="assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
# Output directory
BUSCO_OUT="results/${SPECIES}_stages/busco"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$BUSCO_OUT" logs

# Modules
module purge
module load BUSCO/5.4.2-foss-2021a

# BUSCO runs: in genome mode
# hymenoptera_odb10
busco \
    --in "$ASM" \
    --out "${SPECIES}_busco_hymenoptera" \
    --out_path "$BUSCO_OUT" \
    --mode genome \
    --lineage hymenoptera_odb10 \
    --cpu "$THREADS" \
    -f
# arthropoda_odb10
busco \
    --in "$ASM" \
    --out "${SPECIES}_busco_arthropoda" \
    --out_path "$BUSCO_OUT" \
    --mode genome \
    --lineage arthropoda_odb10 \
    --cpu "$THREADS" \
    -f
# bacteria_odb10
busco \
    --in "$ASM" \
    --out "${SPECIES}_busco_bacteria" \
    --out_path "$BUSCO_OUT" \
    --mode genome \
    --lineage bacteria_odb10 \
    --cpu "$THREADS" \
    -f