#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=build_cds_index
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/build_cds_index_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/build_cds_index_%j.err

# Build CDS index: concatenate every Prokka .ffn CDS nucleotide file from the reference Wolbachia genomes into one fasta and index with samtools faidx
# Usage:
# sbatch 2b_build_cds_idex.sh

set -euo pipefail
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
PROKKA_DIR="${NUWT_DIR}/wolbachia_prokka"
CDS_DIR="${NUWT_DIR}/wolbachia_cds"
# Output directory
CDS_FASTA="${CDS_DIR}/all_wolbachia_cds.ffn"

mkdir -p "${CDS_DIR}" "${NUWT_DIR}/logs"

module load SAMtools/1.13-GCC-10.3.0

# Concatenate the .ffn files without the descriptions, just the gene ID
# /^>/: header line, print only field 1 (">GCA_000008025_00359") and drop the description
# else: sequence line, keep as is
# Input: one Prokka .ffn:
# >GC_000008025_00001 Chromosomal replication initiator ptoein DnaA
# ATGCTCTACGATCGACTAGCACTAGC
# >Ace_00001 Glutamate--tRNA ligase
# ATGAACGTAGATACAGTA
# Output: all_wolbachia_cds.ffn
# >Ace_00001
# ATGAACGTAGATACAGTA
# >GC_000008025_00001
# ATGCTCTACGATCGACTAGCACTAGC
N_FFN=0
while IFS= read -r FFN; do
    awk '/^>/ { print $1; next } { print }' "${FFN}" >> "${CDS_FASTA}"
    N_FFN=$(( N_FFN + 1 ))
done < <(find "${PROKKA_DIR}" -name "*.ffn" | sort)

# Run samtools faidx for indexation of the CDS pool file
# Ace_00001 60  11  60  61
# GC_000008025_00001   120 83  60  61
samtools faidx "${CDS_FASTA}"