#!/usr/bin/env bash
#SBATCH --job-name=illumina_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/illumina_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/illumina_qc_%j.err

# Illumina QC using fastp
# Usage:
# sbatch illumina_qc.sh species

set -euo pipefail

SPECIES="$1"

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
# Output directory
QC_DIR="results/${SPECIES}_stages/read_qc/illumina_qc"
PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$QC_DIR" "$PROC_DIR" logs
FILTERED_R1="${PROC_DIR}/illumina_filtered_R1.fastq.gz"
FILTERED_R2="${PROC_DIR}/illumina_filtered_R2.fastq.gz"
# Input directory 
RAW_R1="${WORKDIR}/reads/illumina/${SPECIES}/illumina_raw_R1.fastq.gz"
RAW_R2="${WORKDIR}/reads/illumina/${SPECIES}/illumina_raw_R2.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Modules
module purge
module load fastp/0.23.4-GCC-10.3.0
module load SeqKit/2.6.1

# Compute raw statistics
seqkit stats -T -a "${RAW_R1}" > "${QC_DIR}/illumina_raw_stats_R1.tsv"
seqkit stats -T -a "${RAW_R2}" > "${QC_DIR}/illumina_raw_stats_R2.tsv"

# Filtering: fastp 
# --detect_adapter_for_pe autodetects adapter sequences in paired end data 
# -q 20 for Q20 Pred scores,  -u 40 and -l 50 
# --cut_right: sliding window from 3' end with cut_mean_quality 20
# --trim_poly_g nabled 
fastp \
    -i "${RAW_R1}" \
    -I "${RAW_R2}" \
    -o "${FILTERED_R1}" \
    -O "${FILTERED_R2}" \
    --detect_adapter_for_pe \
    -q 20 \
    -u 40 \
    -l 50 \
    --cut_right \
    --cut_mean_quality 20 \
    --trim_poly_g \
    -w "${THREADS}" \
    --html "${QC_DIR}/fastp_report.html" \
    --json "${QC_DIR}/fastp_report.json"

# Cleaned statistics
seqkit stats -T -a "${FILTERED_R1}" > "${QC_DIR}/illumina_filtered_stats_R1.tsv"
seqkit stats -T -a "${FILTERED_R2}" > "${QC_DIR}/illumina_filtered_stats_R2.tsv"