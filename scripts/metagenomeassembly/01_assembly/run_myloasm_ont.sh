#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --job-name=myloasm_ont
#SBATCH --output=logs/myloasm_ont_%x_%j.out
#SBATCH --error=logs/myloasm_ont_%x_%j.err

# Runs myloasm on ONT reads 
# Usage
# sbatch run_myloasm_ont.sh species

set -euo pipefail

# Arguments
SPECIES="${1}"

# Paths
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
# Downloaded binary
MYLOASM="${PROJECT_ROOT}/.conda_envs/myloasm_bin/myloasm-0.5.1-x86_64-avx2"
# Read directory
READ_DIR="${PROJECT_ROOT}/reads/processed/${SPECIES}"
ONT_FASTQ="${READ_DIR}/ont_filtered.fastq.gz"
# Output directory
ASM_DIR="${PROJECT_ROOT}/assemblies/${SPECIES}/myloasm_ont"
ASM_GZ="${ASM_DIR}/${SPECIES}_myloasm_ont.fasta.gz"
ONT_FASTA="${READ_DIR}/ont_filtered.fasta.gz"

mkdir -p "${ASM_DIR}"

# modules
module purge
module load SeqKit/2.6.1

# Run myloasm
"${MYLOASM}" \
    "${ONT_FASTQ}" \
    -o "${ASM_DIR}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    --clean-dir

# Compress
gzip -c "${ASM_DIR}/assembly_primary.fa" > "${ASM_GZ}"
# Convert fastq reads to fasta with fq2fa
zcat "${ONT_FASTQ}" | seqkit fq2fa -j "${SLURM_CPUS_PER_TASK}" | gzip -c > "${ONT_FASTA}"
