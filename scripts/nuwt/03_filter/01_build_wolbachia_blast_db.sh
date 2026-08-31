#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=build_wolbachia_db
#SBATCH --output=logs/build_wolbachia_db_%j.out
#SBATCH --error=logs/build_wolbachia_db_%j.err

# Concatenate the downloaded Wolbachia reference genomes into a single fasta 
# Build blast nucleotide database
# Usage:
# sbatch 01_build_wolbachia_db.sh

set -euo pipefail
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
# Directory with the downloaded Wolbachia references
REF_DIR="${PROJECT_ROOT}/nuwt_scan/wolbachia_references"
# Outpu directory
COMBINED="${REF_DIR}/Wolbachia_refs_all.fna"
DB_PREFIX="${REF_DIR}/Wolbachia_refs_db"

module load BLAST+/2.15.0-gompi-2021a

# Concatenate every Wolbachia genome to a single fasta file
# Each genome keeps its header to keep track of the blast hit
cat "${REF_DIR}"/*.fna > "${COMBINED}"

# Build the blast database
# -dbtype nucl for nucleotide blastn
# -parse_seqids to index the sequences id
makeblastdb -in "${COMBINED}" -dbtype nucl -out "${DB_PREFIX}" -parse_seqids