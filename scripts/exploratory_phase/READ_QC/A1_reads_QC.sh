#!/usr/bin/env bash
#SBATCH --job-name=reads_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/reads_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/reads_qc_%j.err

# Stage A: Assembly
# A1: General read QC
# HiFi: fastplong filtering (min 3kb, max 2x N50 if not provided) + FastQC + read length extraction for R curation script
# Hi-C: fastp cleaning (adapter trim, Q20, min 50bp) + FastQC pre/post
# Length distribution curation delegated to A1_length_curation.R
#
# Usage: sbatch A1_reads_qc.sh <species> <asm_mode> [hifi_max_len]
#   asm_mode     : 'hic' or 'bp' (bp skips Hi-C steps)
#   hifi_max_len : max HiFi read length in bp (optional)
#                  if not provided, estimated automatically as 2x raw read N50
#                  override manually for species flagged WARN in curation TSV

set -euo pipefail

# ------------------------------------------------------------
# Args
# ------------------------------------------------------------
SPECIES="$1"
ASM_MODE="$2"
HIFI_MIN_LEN=3000

echo "[INFO] species=$SPECIES"
echo "[INFO] mode=$ASM_MODE"

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

QC_DIR="results/${SPECIES}_stages/read_qc"
HIFI_QC_DIR="${QC_DIR}/hifi_qc"
HIC_QC_DIR="${QC_DIR}/hic_qc"

HIFI_CLEAN_DIR="${HIFI_QC_DIR}/hifi_filtered"
HIC_CLEAN_DIR="reads/hic_clean/${SPECIES}"

mkdir -p "$QC_DIR" "$HIFI_QC_DIR" "$HIC_QC_DIR" "$HIFI_CLEAN_DIR" "$HIC_CLEAN_DIR" logs

# Single HiFi file (as promised)
HIFI=$(ls reads/pacbio_hifi/${SPECIES}/*.fastq.gz)

HIFI_FILTERED="${HIFI_CLEAN_DIR}/hifi.filtered.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ------------------------------------------------------------
# Modules
# ------------------------------------------------------------
module purge
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-10.3.0
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a
module load Anaconda3/2022.05

source "$(conda info --base)/etc/profile.d/conda.sh"
FASTPLONG_ENV="${WORKDIR}/.conda_envs/fastplong"

# ------------------------------------------------------------
# RAW stats
# ------------------------------------------------------------
echo "[INFO] RAW stats"
seqkit stats -T "$HIFI" > "${HIFI_QC_DIR}/hifi_raw_stats.tsv"

# ------------------------------------------------------------
# N50 → max length (optional override)
# ------------------------------------------------------------
if [[ -n "${3:-}" ]]; then
    HIFI_MAX_LEN="$3"
    LENGTH_ARGS=(--length_limit "$HIFI_MAX_LEN")
else
    RAW_N50=$(seqkit stats -T "$HIFI" | awk -F '\t' 'NR==2 {print $8}')
    HIFI_MAX_LEN=$(( RAW_N50 * 2 ))
    LENGTH_ARGS=(--length_limit "$HIFI_MAX_LEN")
fi

echo "[INFO] length filter: min=${HIFI_MIN_LEN} max=${HIFI_MAX_LEN}"

echo -e "species\thifi_min_len\thifi_max_len" > "${HIFI_QC_DIR}/hifi_length_cutoffs.tsv"
echo -e "${SPECIES}\t${HIFI_MIN_LEN}\t${HIFI_MAX_LEN}" >> "${HIFI_QC_DIR}/hifi_length_cutoffs.tsv"

# ------------------------------------------------------------
# RAW length distribution
# ------------------------------------------------------------
echo "[INFO] extracting RAW lengths"
zcat "$HIFI" | awk 'NR%4==2 {print length($0)}' \
> "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt"

# ------------------------------------------------------------
# fastplong
# ------------------------------------------------------------
echo "[INFO] running fastplong"

conda run -p "$FASTPLONG_ENV" fastplong \
    --in "$HIFI" \
    --out "$HIFI_FILTERED" \
    --length_required "$HIFI_MIN_LEN" \
    "${LENGTH_ARGS[@]}" \
    --thread "$THREADS" \
    --html "${HIFI_QC_DIR}/fastplong_report.html" \
    --json "${HIFI_QC_DIR}/fastplong_report.json"

# ------------------------------------------------------------
# FILTERED stats
# ------------------------------------------------------------
echo "[INFO] filtered stats"

seqkit stats -T "$HIFI_FILTERED" > "${HIFI_QC_DIR}/hifi_filtered_stats.tsv"

fastqc -t "$THREADS" -o "$HIFI_QC_DIR" "$HIFI_FILTERED"

# ------------------------------------------------------------
# FILTERED length distribution
# ------------------------------------------------------------
echo "[INFO] extracting FILTERED lengths"

zcat "$HIFI_FILTERED" | awk 'NR%4==2 {print length($0)}' \
> "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt"

# ------------------------------------------------------------
# Diagnostics plot
# ------------------------------------------------------------
echo "[INFO] running diagnostics"

Rscript scripts/stages/A1_reads_diagnostics.R \
    "$SPECIES" \
    "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt" \
    "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt" \
    "$HIFI_QC_DIR"

# ------------------------------------------------------------
# Skip Hi-C if needed
# ------------------------------------------------------------
if [[ "$ASM_MODE" == "bp" ]]; then
    echo "[INFO] skipping Hi-C"
    exit 0
fi

# ------------------------------------------------------------
# Hi-C (unchanged)
# ------------------------------------------------------------
HIC_R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
HIC_R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"

fastqc -t "$THREADS" -o "$HIC_QC_DIR" "$HIC_R1" "$HIC_R2"

HIC_CLEAN_R1="${HIC_CLEAN_DIR}/hic_R1.clean.fastq.gz"
HIC_CLEAN_R2="${HIC_CLEAN_DIR}/hic_R2.clean.fastq.gz"

fastp \
    --in1 "$HIC_R1" \
    --in2 "$HIC_R2" \
    --out1 "$HIC_CLEAN_R1" \
    --out2 "$HIC_CLEAN_R2" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 50 \
    --correction \
    --thread "$THREADS" \
    --html "${HIC_QC_DIR}/fastp_hic_report.html" \
    --json "${HIC_QC_DIR}/fastp_hic_report.json"

fastqc -t "$THREADS" -o "$HIC_QC_DIR" "$HIC_CLEAN_R1" "$HIC_CLEAN_R2"

seqkit stats -T "$HIC_R1" "$HIC_R2" > "${HIC_QC_DIR}/hic_raw_stats.tsv"
seqkit stats -T "$HIC_CLEAN_R1" "$HIC_CLEAN_R2" > "${HIC_QC_DIR}/hic_clean_stats.tsv"

