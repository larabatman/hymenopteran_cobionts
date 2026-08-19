#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --job-name=build_og_lookup
#SBATCH --output=logs/build_og_lookup_%j.out
#SBATCH --error=logs/build_og_lookup_%j.err

# Step 2 of 3 in the OG annotation database build.
# Run after download_og_alignments.sh and before build_og_gene_annotations.py.
#
# Annotates each Wolbachia orthogroup by BLASTing its protein representative (from og_representatives.fna) against SwissProt, then classifies each OG into one of four categories based on the hit description and source organism. 
# The result is og_lookup_table.tsv, which is joined to every NUWT hit during the annotation step.
# SwissProt entries are manually reviewed and have curated, human-readable product names.
#
# blastp: og_representatives.fna contains amino acid sequences (protein alignments from the Zenodo archive). blastp is the correct flavour for
#
# Requires:
# nuwt_scan/hmm_database/og_representatives.fna
# scripts/nuwt/database/build_og_lookup_table.R
#
# Outputs:
# nuwt_scan/hmm_database/og_annotation/og_vs_swissprot.blast6
# nuwt_scan/hmm_database/og_lookup_table.tsv
#
# Usage:
# sbatch scripts/nuwt/database/build_og_lookup_table.sh

set -euo pipefail

PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
HMM_DIR="${NUWT_DIR}/hmm_database"
ANNOT_DIR="${HMM_DIR}/og_annotation"
REP_FASTA="${HMM_DIR}/og_representatives.fna"
R_SCRIPT="${PROJECT_ROOT}/scripts/nuwt/database/build_og_lookup_table.R"
SWISSPROT_DB="${ANNOT_DIR}/swissprot_db"
BLAST_SP="${ANNOT_DIR}/og_vs_swissprot.blast6"
LOOKUP_TSV="${HMM_DIR}/og_lookup_table.tsv"

mkdir -p "${ANNOT_DIR}"

# Module load 
module load BLAST+/2.15.0-gompi-2021a
module load R/4.2.1-foss-2021a

# Step 1: SwissProt BLAST database
# makeblastbd: mandatory first step of any BLAST+ workflow
# blastp cannot search a FASTA file: it needs the sequences indexed into BLAST's binary format
# Download the reviewed UniProt subset, decompress and index it. Keep only the index files.
wget -q -P "${ANNOT_DIR}" \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

gunzip "${ANNOT_DIR}/uniprot_sprot.fasta.gz"

# makeblastdb options:
# -dbtype prot : build a protein database (required for blastp queries)
# -parse_seqids: index sequences by their ID, enabling later retrieval
#                  by accession if needed for debugging
makeblastdb \
    -in "${ANNOT_DIR}/uniprot_sprot.fasta" \
    -dbtype prot \
    -out "${SWISSPROT_DB}" \
    -parse_seqids

rm "${ANNOT_DIR}/uniprot_sprot.fasta"

# blastp
# blastp: protein query vs protein database
# -evalue 1e-5: Permissive threshold for annotation
# -max_target_seqs 5: Return up to 5 hits per query to have some alternatives
# -num_threads 16: blastp parallelises well across queries.
#
# Output format columns -outfmt "6 ...":
# qseqid: OG identifier (e.g. OG0001392)
# sseqid: SwissProt accession
# pident: percent identity 
# length: alignment length
# evalue: E-value of the match
# bitscore: bit score, database-size independent
# stitle: full SwissProt entry title including gene name and organism
# sscinames: scientific name of the source organism (used to distinguish bacterial Wolbachia IS elements from eukaryotic TEs)
blastp \
    -query "${REP_FASTA}" \
    -db "${SWISSPROT_DB}" \
    -out "${BLAST_SP}" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore stitle sscinames" \
    -evalue 1e-5 \
    -num_threads 16 \
    -max_target_seqs 5

# The R script reads og_representatives.fna to get the complete list of OG IDs, reads the blastp results, and classifies each OG into one of four categories: eukaryotic_TE, bacterial_IS, wolbachia_gene, no_swissprot_hit.
Rscript "${R_SCRIPT}" \
    "${REP_FASTA}" \
    "${BLAST_SP}" \
    "${LOOKUP_TSV}"