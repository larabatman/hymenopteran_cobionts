#!/usr/bin/env bash
#SBATCH --job-name=clr_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/clr_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/clr_qc_%j.err

# PacBio long read CLR QC: filter by length only, the quality is left to the assembler to handle
# Minimum: 1000 bp, maximum 2xN50 
# Usage:
# sbatch qc_clr.sh species

set -euo pipefail

SPECIES="$1"
CLR_MIN_LEN=1000

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
# Results directory 
QC_DIR="results/${SPECIES}_stages/read_qc/clr_qc"
PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$QC_DIR" "$PROC_DIR" logs
CLR_FILTERED_FQ="${PROC_DIR}/clr_filtered.fastq.gz"
# Read directory
CLR_RAW="${WORKDIR}/reads/pacbio_clr/${SPECIES}/clr_raw.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
# Modules
module purge
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a

# Raw read lengths and stats 
zcat "${CLR_RAW}" | awk 'NR%4==2 {print length($0)}' > "${QC_DIR}/clr_raw_read_lengths.txt"

# Raw statistics: -t -a for column 13 that is the N50
zcat "${CLR_RAW}" | seqkit stats -T -a > "${QC_DIR}/clr_raw_stats.tsv"

RAW_N50=$(awk -F '\t' 'NR==2 {print $13}' "${QC_DIR}/clr_raw_stats.tsv")
CLR_MAX_LEN=$(( RAW_N50 * 2 ))

# Filter length with seqkit seq
# Cut reads in length window -m minimum and -M maximum
# -g drops zeros
zcat "${CLR_RAW}" | seqkit seq -m "${CLR_MIN_LEN}" -M "${CLR_MAX_LEN}" -g --threads "${THREADS}" | gzip -c > "${CLR_FILTERED_FQ}"

# Read lengths for R histograms
zcat "${CLR_FILTERED_FQ}" | awk 'NR%4==2 {print length($0)}' > "${QC_DIR}/clr_filtered_read_lengths.txt"

# Compute filtered statistics
seqkit stats -T -a "${CLR_FILTERED_FQ}" > "${QC_DIR}/clr_filtered_stats.tsv"

# Diagnostics plot
Rscript scripts/metagenomeassembly/helpers/reads_diagnostics.R "$SPECIES" "${QC_DIR}/clr_raw_read_lengths.txt" "${QC_DIR}/clr_filtered_read_lengths.txt" "$QC_DIR"