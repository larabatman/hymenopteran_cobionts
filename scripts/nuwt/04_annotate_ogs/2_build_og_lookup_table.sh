#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --job-name=build_og_lookup
#SBATCH --output=logs/build_og_lookup_%j.out
#SBATCH --error=logs/build_og_lookup_%j.err

# Annotate the orthogroups by blasting the protein representatives that have been selected against UniProt/SwissProt
# Usage:
# sbatch 2_build_og_lookup_table.sh

set -euo pipefail

# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
HMM_DIR="${NUWT_DIR}/hmm_database"
# The chosen proteins as representatives of their OG
REP_FASTA="${HMM_DIR}/og_representatives.fna"
# Output directory
ANNOT_DIR="${HMM_DIR}/og_annotation"
BLAST_SP="${ANNOT_DIR}/og_vs_swissprot.blast6"
LOOKUP_TSV="${HMM_DIR}/og_lookup_table.tsv"
# Rscript to classify the annotation table 
R_SCRIPT="${PROJECT_ROOT}/scripts/nuwt/database/build_og_lookup_table.R"
# SwissProt database
SWISSPROT_DB="${ANNOT_DIR}/swissprot_db"


mkdir -p "${ANNOT_DIR}"

module load BLAST+/2.15.0-gompi-2021a
module load R/4.2.1-foss-2021a

# Download SwissProt database: UniProt subset
wget -q -P "${ANNOT_DIR}" "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
# Decompress the database
gunzip "${ANNOT_DIR}/uniprot_sprot.fasta.gz"

# Run blastp: build the blast database first
# -dbtype prot to build proteins 
makeblastdb -in "${ANNOT_DIR}/uniprot_sprot.fasta" -dbtype prot -out "${SWISSPROT_DB}" -parse_seqids

rm "${ANNOT_DIR}/uniprot_sprot.fasta"

# Run blastp: -evalue 1e-5, -max_target_seqs 5 to return top 5 hits
# Output format columns -outfmt
# qseqid is the OG identifier, as OG0001392
# sseqid is the SwissProt accession
# pident is the identity 
# length is the alignment length
# evalue is the Evalue of the match
# bitscore is the bit score
# stitle is the SwissProt entry title
blastp \
    -query "${REP_FASTA}" \
    -db "${SWISSPROT_DB}" \
    -out "${BLAST_SP}" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
    -evalue 1e-5 \
    -num_threads 16 \
    -max_target_seqs 5

# Launch Rscript to classify the protein products
Rscript "${R_SCRIPT}" "${REP_FASTA}" "${BLAST_SP}" "${LOOKUP_TSV}"