#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=build_cds_index
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/build_cds_index_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/build_cds_index_%j.err

# Setup for the placement step
#
# Concatenates every Prokka .ffn CDS nucleotide file from all 1'444 Wolbachia reference genomes into a single FASTA, then indexes it with samtools faidx.
# This pooled file is the CDS pool the placement searches. prep_reference_cds.sh subsets it to the dereplicated reference genomes + OGs, and place_one_og.sh runs each OG's nucleotide HMM profile against that subset with nhmmer to find the OG's reference orthologues
# Also his becomes a lookup random-access searched by gene ID instead of opening and scanning a per-genome file. 
# Output:
# nuwt_scan/wolbachia_cds/all_wolbachia_cds.ffn with the concatenated CDS
# nuwt_scan/wolbachia_cds/all_wolbachia_cds.ffn.fai with the samtools index
#
# Usage:
# sbatch scripts/nuwt/supergroups/build_cds_index.sh

set -euo pipefail

PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
PROKKA_DIR="${NUWT_DIR}/wolbachia_prokka"
CDS_DIR="${NUWT_DIR}/wolbachia_cds"
CDS_FASTA="${CDS_DIR}/all_wolbachia_cds.ffn"

mkdir -p "${CDS_DIR}" "${NUWT_DIR}/logs"

module load SAMtools/1.13-GCC-10.3.0

# Concatenate all .ffn files without their headers:
# Find every .ffn under the Prokka output directory.
# awk keeps only the first token of each header line (the gene ID), drops the description. Non-header lines aka sequences pass through unchanged.
# /^>/: header line, print only field 1 (">GCA_000008025_00359") and drop the description
# else: sequence line, keep as is
# Input: one Prokka .ffn:
# >GC_000008025_00001 Chromosomal replication initiator ptoein DnaA
# ATGCTCTACGATCGACTAGCACTAGC
# >Ace_00001 Glutamate--tRNA ligase
# ATGAACGTAGATACAGTA
# Output: all_wolbachia_cds.ffn
# >Ace_00001 Glutamate--tRNA ligase
# ATGAACGTAGATACAGTA
# >GC_000008025_00001
# ATGCTCTACGATCGACTAGCACTAGC
N_FFN=0
while IFS= read -r FFN; do
    awk '/^>/ { print $1; next } { print }' "${FFN}" >> "${CDS_FASTA}"
    N_FFN=$(( N_FFN + 1 ))
done < <(find "${PROKKA_DIR}" -name "*.ffn" | sort)

# SAMtools index
# Builds the .fai index enabling random-access retrieval of any sequence by its ID (the first header token).
# The index: all_wolbachia_cds.ffn.fai one tab-separated line per sequence
# Ace_00001 60  11  60  61
# GC_000008025_00001   120 83  60  61
# Columns: name, sequence length, byte offset where the sequence starts, bases per line, bytes per line
samtools faidx "${CDS_FASTA}"