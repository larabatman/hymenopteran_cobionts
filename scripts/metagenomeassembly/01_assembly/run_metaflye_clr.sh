#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=700G
#SBATCH --job-name=metaflye
#SBATCH --output=logs/metaflye_%x_%j.out
#SBATCH --error=logs/metaflye_%x_%j.err

# Run metaFlye metagenome assembler on CLR cleaned reads
# Usage:
# sbatch run_metaflye.sh species 


set -euo pipefail

SPECIES="${1}"
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
# Cleaned FASTA files
CLR_FASTA="${PROJECT_ROOT}/reads/processed/${SPECIES}/clr_filtered.fasta.gz"
# Output directory: for metaflye intermediate files and the final assembly assembly.fasta
ASM_DIR="${PROJECT_ROOT}/assemblies/${SPECIES}/metaflye"
ASM_OUT="${ASM_DIR}/assembly.fasta"
ASM_GZ="${ASM_DIR}/${SPECIES}_metaflye.fasta.gz"
# Use scratch 
export TMPDIR="/scratch/${USER}/metaflye_${SLURM_JOB_ID}"

# Modules
module purge
module load Flye/2.9-GCC-10.3.0

mkdir -p logs "${ASM_DIR}"  "${TMPDIR}"

# Run metaFlye: --pacbio-raw for CLR reads that handles noisy reads
# --meta for metagenome 
# --min-overlap between reads set at 2000 bp
flye \
    --pacbio-raw "${CLR_FASTA}" \
    --meta \
    --out-dir "${ASM_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --min-overlap 2000

# zip and rename output:
# rename with species name 
gzip -c "${ASM_OUT}" > "${ASM_GZ}"