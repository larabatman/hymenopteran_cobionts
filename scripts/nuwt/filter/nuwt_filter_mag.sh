#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=nuwt_mag
#SBATCH --output=/data/users/lland/cobionts/nuwt_scan/logs/%x_%j.out
#SBATCH --error=/data/users/lland/cobionts/nuwt_scan/logs/%x_%j.err

# MAG contamination filter.
#
# Screens the host assembly against the Wolbachia MAG recovered from this specific individual by the metagenomic pipeline. The MAG is the resident strain, so it catches mis-assembled living Wolbachia that the reference database may miss.
#
# The MAG is located from the curated Wolbachia MAG sheet, which lists one row per recovered Wolbachia genome with its assembler and refinement tool. Those two fields plus the species name determine the bin filename:
# results/<Species>/full_binners/bins/fasta/<tool>/<Species>_<assembler>_<tool>_*.fa.gz
# The sheet lookup and the FASTA assembly are done by prepare_mag_fasta.R.
#
# Usage:
# sbatch nuwt_filter_mag.sh <Species> <Accession> <assembly.fna>

set -euo pipefail

SPECIES="${1}"
ACCESSION="${2}"
ASSEMBLY="${3}"


SHEET="/data/projects/p2025-0083_mining_cobionts/sample_sheets/Wolbachia_MAGs_1607.tsv"
COBIONTS_ROOT="/data/users/lland/cobionts"
SCRIPTS_DIR="${COBIONTS_ROOT}/scripts/nuwt/filter"
R_MODULE="R/4.2.1-foss-2021a"


RESULTS_DIR="${COBIONTS_ROOT}/nuwt_scan/results/${SPECIES}"
FILTER_DIR="${RESULTS_DIR}/filter"
RAW_DFAM="${RESULTS_DIR}/${SPECIES}_nuwt_hits_dfam.tbl"
BLAST_MAG="${FILTER_DIR}/assembly_vs_MAG.blast6"
MAG_FASTA="${FILTER_DIR}/${SPECIES}_wolbachia_mags.fa"
MAG_DB="${FILTER_DIR}/${SPECIES}_wolbachia_mag_db"
FILTERED_DFAM="${FILTER_DIR}/${SPECIES}_nuwt_hits_dfam_filtered.tbl"
FLAGS_FILE="${FILTER_DIR}/${SPECIES}_contamination_flags.txt"
R_PREPARE="${SCRIPTS_DIR}/prepare_mag_fasta.R"
R_FILTER="${SCRIPTS_DIR}/filter_dfam.R"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "[$(date)] MAG filter: ${SPECIES} (${ACCESSION})"

module load BLAST+/2.15.0-gompi-2021a
module load "${R_MODULE}"

# Build the combined MAG FASTA
Rscript "${R_PREPARE}" "${SPECIES}" "${SHEET}" "${COBIONTS_ROOT}" "${MAG_FASTA}"

# BLATn against the MAGs
# -perc_identity 90 is a broad capture
makeblastdb -in "${MAG_FASTA}" -dbtype nucl -out "${MAG_DB}"

blastn \
    -query           "${ASSEMBLY}" \
    -db              "${MAG_DB}" \
    -out             "${BLAST_MAG}" \
    -outfmt          "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
    -perc_identity   90 \
    -max_target_seqs 5 \
    -num_threads     "${THREADS}"

# Re-make the filtered dfam from BOTH reference and MAG evidence
Rscript "${R_FILTER}" "${RAW_DFAM}" "${FILTER_DIR}" "${FILTERED_DFAM}" "${FLAGS_FILE}"