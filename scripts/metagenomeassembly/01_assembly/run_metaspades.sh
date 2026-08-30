#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=768G
#SBATCH --job-name=metaspades
#SBATCH --output=logs/metaspades_%x_%j.out
#SBATCH --error=logs/metaspades_%x_%j.err

# Run metaSPAdes metagenome assembler on Illumina reads 
# Usage:
# sbatch run_metaspades.sh species

set -euo pipefail

SPECIES="${1}"

# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"

# Cleaned Illumina paired end reads directory
FILTERED_R1="${PROJECT_ROOT}/reads/processed/${SPECIES}/illumina_filtered_R1.fastq.gz"
FILTERED_R2="${PROJECT_ROOT}/reads/processed/${SPECIES}/illumina_filtered_R2.fastq.gz"
# Output directory for metaspades and selecting scaffolds.fasta as the assembly 
ASM_DIR="${PROJECT_ROOT}/assemblies/${SPECIES}/metaspades"
ASM_OUT="${ASM_DIR}/scaffolds.fasta"
ASM_GZ="${ASM_DIR}/${SPECIES}_metaspades.fasta.gz"
# Use scratch
export TMPDIR="/scratch/${USER}/metaspades_${SLURM_JOB_ID}"
mkdir -p logs "${ASM_DIR}" "${TMPDIR}"
# Modules
module purge
module load SPAdes/3.15.3-GCC-10.3.0

# metaSPAdes: --meta for metagenome assembly mode 
# -m set manually to include a small memory buffer and conversion to GB
MEM_GB=$(( SLURM_MEM_PER_NODE / 1024 - 8 ))
spades.py \
    --meta \
    -1 "${FILTERED_R1}" \
    -2 "${FILTERED_R2}" \
    -o "${ASM_DIR}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -m "${MEM_GB}"

# zip output
gzip -c "${ASM_OUT}" > "${ASM_GZ}"