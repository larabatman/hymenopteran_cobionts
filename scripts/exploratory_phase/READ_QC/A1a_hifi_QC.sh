#!/usr/bin/env bash
#SBATCH --job-name=hifi_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hifi_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hifi_qc_%j.err

# =============================================================================
# Stage A1a: HiFi read QC
# fastplong filtering (min 3kb, max 2x N50 or manual override) + FastQC +
# read length extraction + R diagnostics plot
#
# Usage: sbatch A1a_hifi_qc.sh <species> [hifi_max_len]
#   hifi_max_len : max HiFi read length in bp (optional)
#                  if not provided, estimated automatically as 2x raw read N50
# =============================================================================

set -euo pipefail

SPECIES="$1"
HIFI_MIN_LEN=3000

echo "[INFO] species=$SPECIES"

# ── Paths ──
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

HIFI_QC_DIR="results/${SPECIES}_stages/read_qc/hifi_qc"
HIFI_CLEAN_DIR="${HIFI_QC_DIR}/hifi_filtered"
mkdir -p "$HIFI_QC_DIR" "$HIFI_CLEAN_DIR" logs

HIFI=$(ls reads/pacbio_hifi/${SPECIES}/*.fastq.gz)
HIFI_FILTERED="${HIFI_CLEAN_DIR}/hifi.filtered.fastq.gz"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ── Modules ──
module purge
module load FastQC/0.11.9-Java-11
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a
module load Anaconda3/2022.05

source "$(conda info --base)/etc/profile.d/conda.sh"
FASTPLONG_ENV="${WORKDIR}/.conda_envs/fastplong"

# ── RAW stats ──
echo "[INFO] RAW stats"
seqkit stats -T "$HIFI" > "${HIFI_QC_DIR}/hifi_raw_stats.tsv"

# ── N50 → max length ──
if [[ -n "${2:-}" ]]; then
    HIFI_MAX_LEN="$2"
else
    RAW_N50=$(seqkit stats -T "$HIFI" | awk -F '\t' 'NR==2 {print $8}')
    HIFI_MAX_LEN=$(( RAW_N50 * 2 ))
fi

echo "[INFO] length filter: min=${HIFI_MIN_LEN} max=${HIFI_MAX_LEN}"

echo -e "species\thifi_min_len\thifi_max_len" > "${HIFI_QC_DIR}/hifi_length_cutoffs.tsv"
echo -e "${SPECIES}\t${HIFI_MIN_LEN}\t${HIFI_MAX_LEN}" >> "${HIFI_QC_DIR}/hifi_length_cutoffs.tsv"

# ── RAW length distribution ──
echo "[INFO] extracting RAW lengths"
zcat "$HIFI" | awk 'NR%4==2 {print length($0)}' \
    > "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt"

# ── fastplong ──
echo "[INFO] running fastplong"
conda run -p "$FASTPLONG_ENV" fastplong \
    --in "$HIFI" \
    --out "$HIFI_FILTERED" \
    --length_required "$HIFI_MIN_LEN" \
    --length_limit "$HIFI_MAX_LEN" \
    --thread "$THREADS" \
    --html "${HIFI_QC_DIR}/fastplong_report.html" \
    --json "${HIFI_QC_DIR}/fastplong_report.json"

# ── FILTERED stats ──
echo "[INFO] filtered stats"
seqkit stats -T "$HIFI_FILTERED" > "${HIFI_QC_DIR}/hifi_filtered_stats.tsv"
fastqc -t "$THREADS" -o "$HIFI_QC_DIR" "$HIFI_FILTERED"

# ── FILTERED length distribution ──
echo "[INFO] extracting FILTERED lengths"
zcat "$HIFI_FILTERED" | awk 'NR%4==2 {print length($0)}' \
    > "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt"

# ── Diagnostics plot ──
echo "[INFO] running diagnostics"
Rscript scripts/exploratory_phase/READ_QC/A1_reads_diagnostics.R \
    "$SPECIES" \
    "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt" \
    "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt" \
    "$HIFI_QC_DIR"

echo "[OK] HiFi QC complete for $SPECIES"