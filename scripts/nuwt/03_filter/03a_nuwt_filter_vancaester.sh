#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=nuwt_vanc
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/cobionts/nuwt_scan/logs/nuwt_vanc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/cobionts/nuwt_scan/logs/nuwt_vanc_%j.err

# Filter each dfam.tbl against the reference database 
# blastn the host assembly against the full Wolbachia reference blast database
# Scaffold or contigs that match a reference at more than 99% over 80% of the contig length is treated as living Wolbachia contamination and not a nuwt
# Usage:
# sbatch 03a_nuwt_filter_vancaester.sh species assembly.fna


set -euo pipefail

SPECIES="${1}"
ASSEMBLY="${2}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
cd "$PROJECT_ROOT"
# Path to Rscript
R_FILTER="scripts/nuwt/03_filter/03c_filter_dfam.R"
# Directory with dfam tables
RESULTS_DIR="cobionts/nuwt_scan/results/${SPECIES}"
RAW_DFAM="${RESULTS_DIR}/${SPECIES}_nuwt_hits_dfam.tbl"
# Reference Wolbachia database with the files that blastn needs 
WOLBACHIA_DB="nuwt_scan/wolbachia_references/Wolbachia_refs_db"
# Output directories
FILTER_DIR="${RESULTS_DIR}/filter"
BLAST_REFS="${FILTER_DIR}/assembly_vs_wolbachia_refs.blast6"
FILTERED_DFAM="${FILTER_DIR}/${SPECIES}_nuwt_hits_dfam_filtered.tbl"
FLAGS_FILE="${FILTER_DIR}/${SPECIES}_contamination_flags.txt"

mkdir -p "${FILTER_DIR}"

module load BLAST+/2.15.0-gompi-2021a
module load SAMtools/1.13-GCC-10.3.0
module load R/4.2.1-foss-2021a

# Index the assembly
samtools faidx "${ASSEMBLY}"

# blastn host assembly against the Wolbachia reference database
# -per_identity set at 90 and filter afterwards
blastn \
    -query "${ASSEMBLY}" \
    -db "${WOLBACHIA_DB}" \
    -out "${BLAST_REFS}" \
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
    -perc_identity 90 \
    -max_target_seqs 5 \
    -num_threads "${THREADS}"

# Launch the filtering script
Rscript "${R_FILTER}" "${RAW_DFAM}" "${FILTER_DIR}" "${FILTERED_DFAM}" "${FLAGS_FILE}"