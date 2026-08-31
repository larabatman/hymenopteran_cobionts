#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=nuwt_mag
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/cobionts/nuwt_scan/logs/nuwt_mag_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/cobionts/nuwt_scan/logs/nuwt_mag_%j.err

# blastn the host assembly against the Wolbachia MAGs that was retrieved from it
# Usage
# sbatch 03b_nuwt_filter_mag.sh species assembly

set -euo pipefail

SPECIES="${1}"
ASSEMBLY="${2}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
cd "$PROJECT_ROOT"
# Wolbachia MAG table
WOLBACHIA_MAG="/sample_sheet/Wolbachia_MAGs.tsv"
# Scripts directory
SCRIPTS_DIR="/scripts/nuwt/03_filter"
R_PREPARE="${SCRIPTS_DIR}/02_prepare_mag_fasta.R"
R_FILTER="${SCRIPTS_DIR}/03c_filter_dfam.R"
# Output directory
RESULTS_DIR="cobionts/nuwt_scan/results/${SPECIES}"
RAW_DFAM="${RESULTS_DIR}/${SPECIES}_nuwt_hits_dfam.tbl"
# Directory for the filtered files
FILTER_DIR="${RESULTS_DIR}/filter"
BLAST_MAG="${FILTER_DIR}/assembly_vs_MAG.blast6"
MAG_FASTA="${FILTER_DIR}/${SPECIES}_wolbachia_mags.fa"
MAG_DB="${FILTER_DIR}/${SPECIES}_wolbachia_mag_db"
FILTERED_DFAM="${FILTER_DIR}/${SPECIES}_nuwt_hits_dfam_filtered.tbl"
FLAGS_FILE="${FILTER_DIR}/${SPECIES}_contamination_flags.txt"

module load BLAST+/2.15.0-gompi-2021a
module load R/4.2.1-foss-2021a

# Launch the script to build the combined mag fasta file, containing all the mags for one species
Rscript "${R_PREPARE}" "${SPECIES}" "${SHEET}" "${COBIONTS_ROOT}" "${MAG_FASTA}"

# blastn the genomes that had mags against their own mags
# Make a blast database of the mags fasta file
makeblastdb -in "${MAG_FASTA}" -dbtype nucl -out "${MAG_DB}"

blastn \
    -query "${ASSEMBLY}" \
    -db "${MAG_DB}" \
    -out "${BLAST_MAG}" \
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
    -perc_identity 90 \
    -max_target_seqs 5 \
    -num_threads "${THREADS}"

# Launch the script to filter the resulting blast6 and dfam tables
Rscript "${R_FILTER}" "${RAW_DFAM}" "${FILTER_DIR}" "${FILTERED_DFAM}" "${FLAGS_FILE}"