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

# ------------------------------------------------------------
# 1) Setup
# ------------------------------------------------------------
set -euo pipefail
# Ensure proper usage
usage() {
    echo "Usage: $0 <species> <asm_mode> [hifi_max_len]"
    exit 1
}

[[ $# -lt 2 ]] && usage
# Arguments
SPECIES="$1"
ASM_MODE="$2"
# Hardcoded HIFI_MIN_LEN
HIFI_MIN_LEN=3000

# Log arguments
echo "[INFO] job started on $(hostname)"
echo "[INFO] species=$SPECIES"
echo "[INFO] assembly mode=$ASM_MODE"
# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

QC_DIR="results/${SPECIES}_stages/read_qc"
HIFI_QC_DIR="${QC_DIR}/hifi_qc"
HIC_QC_DIR="${QC_DIR}/hic_qc"

HIFI_CLEAN_DIR="${HIFI_QC_DIR}/hifi_filtered"
HIC_CLEAN_DIR="reads/hic_clean/${SPECIES}"

mkdir -p "$QC_DIR" "$HIFI_QC_DIR" "$HIC_QC_DIR" "$HIFI_CLEAN_DIR" "$HIC_CLEAN_DIR" logs

# Check on HiFi
HIFI=(reads/pacbio_hifi/${SPECIES}/*.fastq.gz)
HIC_R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
HIC_R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"

HIFI_FILTERED="${HIFI_CLEAN_DIR}/hifi.filtered.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "[INFO] Using $THREADS threads"
# Load modules
module purge
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-10.3.0
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a
# fastplong is used through Conda env (see ENV_create_fastplong.sh)
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"

FASTPLONG_ENV="${WORKDIR}/.conda_envs/fastplong"

# ------------------------------------------------------------
# 2) RAW HiFi stats
# ------------------------------------------------------------

# This we could use as well to get the N50 information for later stages
echo "[INFO] Computing HiFi raw stats"
seqkit stats -T "${HIFI[@]}" > "${HIFI_QC_DIR}/hifi_raw_stats.tsv"

# ------------------------------------------------------------
# 3) Determine max length cutoff
# ------------------------------------------------------------

# If we don't have given HIFI_MAX_LEN argument, default to 2XN50
# seqkit stats -N 50 to find N50 and set to 2XN50
if [[ -n "${3:-}" ]]; then
    HIFI_MAX_LEN="$3"
else
    RAW_N50=$(seqkit stats -N 50 -T "${HIFI[@]}" \
        | tail -n +2 \
        | awk '{print $NF}')
    HIFI_MAX_LEN=$(( RAW_N50 * 2 ))
fi

# Log length filtering
echo "[INFO] HiFi filter: min=${HIFI_MIN_LEN} max=${HIFI_MAX_LEN}"

echo -e "species\thifi_min_len\thifi_max_len" > "${HIFI_QC_DIR}/hifi_length_cutoffs.tsv"
echo -e "${SPECIES}\t${HIFI_MIN_LEN}\t${HIFI_MAX_LEN}" >> "${HIFI_QC_DIR}/hifi_length_cutoffs.tsv"

# ------------------------------------------------------------
# 4) RAW read length distribution
# ------------------------------------------------------------

# From raw reads, select line 2 that contains the sequence and extract length
echo "[INFO] Extracting RAW read lengths"

zcat "${HIFI[@]}" \
| awk 'NR%4==2 {print length($0)}' \
> "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt"

# ------------------------------------------------------------
# 5) Filter HiFi reads
# ------------------------------------------------------------

# By default, fastplong trims adapters

echo "[INFO] Running fastplong"

conda run -p "$FASTPLONG_ENV" fastplong \
    --in "${HIFI[@]}" \
    --out "$HIFI_FILTERED" \
    --length_required "$HIFI_MIN_LEN" \
    --length_limit "$HIFI_MAX_LEN" \
    --thread "$THREADS" \
    --html "${HIFI_QC_DIR}/fastplong_report.html" \
    --json "${HIFI_QC_DIR}/fastplong_report.json"

# ------------------------------------------------------------
# 6) Filtered stats
# ------------------------------------------------------------

echo "[INFO] Computing filtered stats"

seqkit stats -T "$HIFI_FILTERED" > "${HIFI_QC_DIR}/hifi_filtered_stats.tsv"

fastqc \
    -t "$THREADS" \
    -o "$HIFI_QC_DIR" \
    "$HIFI_FILTERED"

# ------------------------------------------------------------
# 7) Filtered read length distribution
# ------------------------------------------------------------

# From filtered reads, select line 2 that contains the sequence and extract length
echo "[INFO] Extracting FILTERED read lengths"

zcat "$HIFI_FILTERED" \
| awk 'NR%4==2 {print length($0)}' \
> "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt"

# ------------------------------------------------------------
# 8) Run diagnostics R script
# ------------------------------------------------------------

# This creates before vs after length filtering plots
echo "[INFO] Running length diagnostics"

Rscript scripts/stages/A1_reads_diagnostics.R \
    "$SPECIES" \
    "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt" \
    "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt" \
    "$HIFI_QC_DIR"

# ------------------------------------------------------------
# 9) Skip Hi-C if bp mode
# ------------------------------------------------------------

if [[ "$ASM_MODE" == "bp" ]]; then
    echo "[INFO] ASM_MODE=bp: skipping Hi-C QC"
    exit 0
fi

# ------------------------------------------------------------
# 10) Hi-C QC
# ------------------------------------------------------------

echo "[INFO] Hi-C raw FastQC"

fastqc -t "$THREADS" -o "$HIC_QC_DIR" "$HIC_R1" "$HIC_R2"

echo "[INFO] Cleaning Hi-C reads"

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

echo "[INFO] Hi-C post-clean FastQC"

fastqc -t "$THREADS" -o "$HIC_QC_DIR" "$HIC_CLEAN_R1" "$HIC_CLEAN_R2"

seqkit stats -T "$HIC_R1" "$HIC_R2" > "${HIC_QC_DIR}/hic_raw_stats.tsv"
seqkit stats -T "$HIC_CLEAN_R1" "$HIC_CLEAN_R2" > "${HIC_QC_DIR}/hic_clean_stats.tsv"

echo "[INFO] A1 Read QC completed successfully"