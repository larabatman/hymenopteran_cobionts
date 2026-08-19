#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=nuwt_vanc
#SBATCH --output=/data/users/lland/cobionts/nuwt_scan/logs/%x_%j.out
#SBATCH --error=/data/users/lland/cobionts/nuwt_scan/logs/%x_%j.err

# Vancaester reference-database contamination filter.
#
# This script BLASTn the host assembly against the full Wolbachia reference database (Vancaester et al. 2025). 
# Any host scaffold matching a reference at >=99% identity over >=80% of its length is flagged as living-Wolbachia contamination (mis-assembled into the host) rather than a true NUWT, and is dropped from the dfam table.
#
# Conversions:
# host .fna samtools faidx to .fai
# host .fna + ref db through blastn to .blast6
# raw _dfam.tbl to filter_dfam.R  -> _dfam_filtered.tbl
#
# Usage
# sbatch nuwt_filter_vancaester.sh <Species> <Accession> <assembly.fna>

set -euo pipefail

SPECIES="${1}"
ACCESSION="${2}"
ASSEMBLY="${3}"

WOLBACHIA_DB="/data/projects/p2025-0083_mining_cobionts/nuwt_scan/wolbachia_references/Wolbachia_refs_db"

R_MODULE="R/4.2.1-foss-2021a"

SCRIPTS_DIR="/data/users/lland/cobionts/scripts/nuwt/filter"


COBIONTS_ROOT="/data/users/lland/cobionts"
RESULTS_DIR="${COBIONTS_ROOT}/nuwt_scan/results/${SPECIES}"
FILTER_DIR="${RESULTS_DIR}/filter"
RAW_DFAM="${RESULTS_DIR}/${SPECIES}_nuwt_hits_dfam.tbl"
BLAST_REFS="${FILTER_DIR}/assembly_vs_wolbachia_refs.blast6"
FILTERED_DFAM="${FILTER_DIR}/${SPECIES}_nuwt_hits_dfam_filtered.tbl"
FLAGS_FILE="${FILTER_DIR}/${SPECIES}_contamination_flags.txt"
R_FILTER="${SCRIPTS_DIR}/filter_dfam.R"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${FILTER_DIR}"

# Modules
module load BLAST+/2.15.0-gompi-2021a
module load SAMtools/1.13-GCC-10.3.0
module load "${R_MODULE}"

# Index assembly
samtools faidx "${ASSEMBLY}"

# BLASTn host assembly vs Wolbachia reference database
# -perc_identity 90: broader capture, the strict 99% bar is applied in R so that the same threshold logic is shared with the MAG filter.
blastn \
    -query "${ASSEMBLY}" \
    -db "${WOLBACHIA_DB}" \
    -out "${BLAST_REFS}" \
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
    -perc_identity 90 \
    -max_target_seqs 5 \
    -num_threads "${THREADS}"

# Filtered dfam
Rscript "${R_FILTER}" "${RAW_DFAM}" "${FILTER_DIR}" "${FILTERED_DFAM}" "${FLAGS_FILE}"