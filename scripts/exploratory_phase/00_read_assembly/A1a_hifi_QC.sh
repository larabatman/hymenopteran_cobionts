#!/usr/bin/env bash
#SBATCH --job-name=hifi_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hifi_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hifi_qc_%j.err

# This script filters the raw hifi reads with fastplong:
# minimum 3000 bp long, maximum 2xN50 
# FastQC report
# Read and length extraction for R histograms in A1_reads_diagnostics
# Usage:
# sbatch A1a_hifi_qc.sh species  

set -euo pipefail

# Arguments
SPECIES="$1"

# Min and max lengths
HIFI_MIN_LEN=3000
HIFI_MIN_QUAL=20

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

# Read directory
HIFI=reads/pacbio_hifi/${SPECIES}/*.fastq.gz
HIFI_FILTERED="${HIFI_CLEAN_DIR}/hifi.filtered.fastq.gz"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Outdir: clean rean directory and plots in species directory
HIFI_QC_DIR="results/${SPECIES}_stages/read_qc/hifi_qc"
HIFI_CLEAN_DIR="${HIFI_QC_DIR}/hifi_filtered"
mkdir -p "$HIFI_QC_DIR" "$HIFI_CLEAN_DIR" logs

# Modules
module purge
module load FastQC/0.11.9-Java-11
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a
module load Anaconda3/2022.05

FASTPLONG_ENV="${WORKDIR}/.conda_envs/fastplong"

# Raw statistics: grab the N50 for the maximum using seqkit
# -T tap separated output for awk parsing
# -a includes all stats and column 13 is the N50
# Output to hifi_raw_stats.tsv
seqkit stats -T -a "$HIFI" > "${HIFI_QC_DIR}/hifi_raw_stats.tsv"

# Select max length: 2xN50 
# From hifi_raw_stats.tsv:
# -F '\t' specifies the tab separator of seqkit -T
# NR==2 skips the header written by seqkit
# columns 13 $13 is the N50
RAW_N50=$(awk -F '\t' 'NR==2 {print $13}' "${HIFI_QC_DIR}/hifi_raw_stats.tsv")
HIFI_MAX_LEN=$(( RAW_N50 * 2 ))

# Extract read length for histogram:
# In FASTQ, there are 4 lines: 
# Line 1: header that starts with @ 
# Line 2: sequence which is what we want awk to grab and read length with length($0)
# Line 3: separator that starts with +
# Line 4: quality string with a character for each base
# NR==2 matches the second line in one file, but we need to match the second line in eery redor: NR%4==2 matches 2, 6, 10, 14 keeping the line if it is the second of every group of 4
zcat "$HIFI" | awk 'NR%4==2 {print length($0)}' > "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt"

# Clean raw reads: fastplong for length and quality
# --mean_qual: Phred score threshold, Q20 keeps bases tha are 99% accurate
conda run -p "$FASTPLONG_ENV" fastplong \
    --in "$HIFI" \
    --out "$HIFI_FILTERED" \
    --length_required "$HIFI_MIN_LEN" \
    --length_limit "$HIFI_MAX_LEN" \
    --mean_qual "$HIFI_MIN_QUAL" \
    --thread "$THREADS" \
    --html "${HIFI_QC_DIR}/fastplong_report.html" \
    --json "${HIFI_QC_DIR}/fastplong_report.json"

# Extract filter stats with seqkit
seqkit stats -T "$HIFI_FILTERED" > "${HIFI_QC_DIR}/hifi_filtered_stats.tsv"

# Write FastQC report 
fastqc -t "$THREADS" -o "$HIFI_QC_DIR" "$HIFI_FILTERED"

# Extract clean read length for histograms
zcat "$HIFI_FILTERED" | awk 'NR%4==2 {print length($0)}' > "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt"

# Diagnostics plot
Rscript scripts/exploratory_phase/00_read_assembly/A1_reads_diagnostics.R "$SPECIES" "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt" "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt" "$HIFI_QC_DIR"