#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=wol_prokka
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/wol_prokka_%A_%a.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/wol_prokka_%A_%a.err

# Run Prokka on one Wolbachia reference genome
# Submitted as a job array to run on every reference genome
# This regenerates the .ffn CDS nucleotide files for the references, which we need to find the original nucleotide sequences from which each OG came from
# Usage: 
# Each job reads the nth line of the genome list file
# sbatch --array=1-${N}%20 2a_wolb_prokka_cds_array.sh

set -euo pipefail
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
# File with the path of every genome for which we need to find CDS
GENOME_LIST="${NUWT_DIR}/wolbachia_genome_list.txt"
# Output directory
PROKKA_DIR="${NUWT_DIR}/wolbachia_prokka"

mkdir -p "${PROKKA_DIR}"

# Pick the genome for this job from the line number in the genome list
GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${GENOME_LIST}")

# Retrieve its accession from the filename
ACCESSION=$(basename "${GENOME}" .fna)
OUT_DIR="${PROKKA_DIR}/${ACCESSION}"

module load prokka/1.14.5-gompi-2021a

# Run Prokka: --locustag accession to keep track of the strain genome source
prokka \
    --quiet \
    --cpus 4 \
    --kingdom Bacteria \
    --locustag "${ACCESSION}" \
    --prefix "${ACCESSION}" \
    --outdir "${OUT_DIR}" \
    --force \
    "${GENOME}"