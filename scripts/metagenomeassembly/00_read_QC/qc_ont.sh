#!/usr/bin/env bash
#SBATCH --job-name=ont_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/ont_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/ont_qc_%j.err

# Length filtering long read ONT, qulaity filtering is left to the assembler to handle 
# Minimal length: 1000 kbp, maximal length 2xN50 
# Usage
# sbatch qc_ont.sh species
set -euo pipefail

SPECIES="$1"
ONT_MIN_LEN=1000
# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
# Output directory
QC_DIR="results/${SPECIES}_stages/read_qc/ont_qc"
PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$QC_DIR" "$PROC_DIR" logs
ONT_FILTERED_FQ="${PROC_DIR}/ont_filtered.fastq.gz"
# Raw reads directory
ONT_RAW="${WORKDIR}/reads/ont/${SPECIES}/ont_raw.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Modules
module purge
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a

# Raw read length
zcat "${ONT_RAW}" | awk 'NR%4==2 {print length($0)}' > "${QC_DIR}/ont_raw_read_lengths.txt"

# Seqkit stats to extract N50
zcat "${ONT_RAW}" | seqkit stats -T -a > "${QC_DIR}/ont_raw_stats.tsv"

RAW_N50=$(awk -F '\t' 'NR==2 {print $13}' "${QC_DIR}/ont_raw_stats.tsv")
ONT_MAX_LEN=$(( RAW_N50 * 2 ))

# Filter with seqkit for intervals below -m and above -M
zcat "${ONT_RAW}" | seqkit seq -m "${ONT_MIN_LEN}" -M "${ONT_MAX_LEN}" -g --threads "${THREADS}" | gzip -c > "${ONT_FILTERED_FQ}"

# Cleaned read length
zcat "${ONT_FILTERED_FQ}" | awk 'NR%4==2 {print length($0)}' > "${QC_DIR}/ont_filtered_read_lengths.txt"

# Cleaned statistics
seqkit stats -T -a "${ONT_FILTERED_FQ}" > "${QC_DIR}/ont_filtered_stats.tsv"

# Diagnostics plot
Rscript scripts/metagenomeassembly/helpers/reads_diagnostics.R "$SPECIES" "${QC_DIR}/ont_raw_read_lengths.txt" "${QC_DIR}/ont_filtered_read_lengths.txt" "$QC_DIR"