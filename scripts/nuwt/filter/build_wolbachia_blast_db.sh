#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=build_wolbachia_db
#SBATCH --output=logs/build_wolbachia_db_%j.out
#SBATCH --error=logs/build_wolbachia_db_%j.err

# Concatenates the downloaded Wolbachia reference genomes into a single FASTA and builds a BLAST nucleotide database from it.
#
# The database is used by nuwt_filter_vancaester.sh. A host scaffold covered over >=80% of its length at >=99% identity is flagged as living-Wolbachia  contamination mis-assembled into the host, not a true NUWT.
#
# Files
# in: one nucleotide FASTA per genome, GCA_*.fna
# out: 
# Wolbachia_refs_all.fna as concatenated
# Wolbachia_refs_db.* with BLAST binary index: .nhr .nin .nsq ...
# Usage:
# sbatch build_wolbachia_db.sh

set -euo pipefail

PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
REF_DIR="${PROJECT_ROOT}/nuwt_scan/wolbachia_references"
COMBINED="${REF_DIR}/Wolbachia_refs_all.fna"
DB_PREFIX="${REF_DIR}/Wolbachia_refs_db"

module load BLAST+/2.15.0-gompi-2021a

# Concatenate everything
# Each genome keeps its own headers, so a BLAST hit can be traced back to its source genome. Sorted order makes the combined file reproducible.
cat "${REF_DIR}"/*.fna > "${COMBINED}"

# Build BLAST DB
# -dbtype nucl: nucleotide database for blastn
# -parse_seqids: index by sequence id, so a matched reference contig can be retrieved by ID 
makeblastdb \
    -in     "${COMBINED}" \
    -dbtype nucl \
    -out    "${DB_PREFIX}" \
    -parse_seqids
