#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=300:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --job-name=myloasm
#SBATCH --output=logs/myloasm_%x_%j.out
#SBATCH --error=logs/myloasm_%x_%j.err

# Run myloasm on HiFi reads then converts FASTQ to FASTA.
# For reads that could not be assembled with metaMDBG
# Usage:
# sbatch run_myloasm species

set -euo pipefail

SPECIES="${1}"

# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
# Downloaded binary file to run myloasm
MYLOASM="${PROJECT_ROOT}/.conda_envs/myloasm_bin/myloasm-0.5.1-x86_64-avx2"
# Input reads directory
HIFI_FASTQ="${PROJECT_ROOT}/reads/processed/${SPECIES}/hifi_filtered.fastq.gz"
# Output directory
HIFI_FASTA="${PROJECT_ROOT}/reads/processed/${SPECIES}/hifi_filtered.fasta.gz"
ASM_DIR="${PROJECT_ROOT}/assemblies/${SPECIES}/myloasm"
ASM_OUT="${ASM_DIR}/assembly_primary.fa"
ASM_GZ="${ASM_DIR}/${SPECIES}_myloasm.fasta.gz"

mkdir -p logs "${ASM_DIR}"

# Modules
module purge
module load SeqKit/2.6.1

# myloasm: --hifi for HiFi mode
"${MYLOASM}" \
    "${HIFI_FASTQ}" \
    --hifi \
    -o "${ASM_DIR}" \
    -t "${SLURM_CPUS_PER_TASK}"

# Compress assembly
gzip -c "${ASM_OUT}" > "${ASM_GZ}"

# Convert fastq to fasta: seqkit fq2fa
# takes the + separator line and quality scores away 
zcat "${HIFI_FASTQ}" | seqkit fq2fa -j "${SLURM_CPUS_PER_TASK}" | gzip -c > "${HIFI_FASTA}"