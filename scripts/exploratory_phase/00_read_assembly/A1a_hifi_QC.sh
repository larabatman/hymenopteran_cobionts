#!/usr/bin/env bash
#SBATCH --job-name=hifi_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hifi_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hifi_qc_%j.err

# Stage A1a: HiFi read QC
# fastplong filtering (min 3kb, max 2x N50 or manual override) + FastQC + read length extraction + R diagnostics plot
#
# Usage: sbatch A1a_hifi_qc.sh <species> [hifi_max_len]
# hifi_max_len : max HiFi read length in bp (optional). If not provided, estimated automatically as 2x raw read N50

set -euo pipefail

SPECIES="$1"
HIFI_MIN_LEN=3000
HIFI_MIN_QUAL=20

echo "[INFO] species=$SPECIES"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

HIFI_QC_DIR="results/${SPECIES}_stages/read_qc/hifi_qc"
HIFI_CLEAN_DIR="${HIFI_QC_DIR}/hifi_filtered"
mkdir -p "$HIFI_QC_DIR" "$HIFI_CLEAN_DIR" logs

HIFI=$(ls reads/pacbio_hifi/${SPECIES}/*.fastq.gz)
HIFI_FILTERED="${HIFI_CLEAN_DIR}/hifi.filtered.fastq.gz"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Modules
module purge
module load FastQC/0.11.9-Java-11
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a
module load Anaconda3/2022.05

source "$(conda info --base)/etc/profile.d/conda.sh"
FASTPLONG_ENV="${WORKDIR}/.conda_envs/fastplong"

# RAW stats
# -T: tab-separated output for consistent awk parsing downstream.
# -a: include ALL statistics, column 13 is N50, used for the max cutoff.
echo "[INFO] RAW stats"
seqkit stats -T -a "$HIFI" > "${HIFI_QC_DIR}/hifi_raw_stats.tsv"

# Max length from N50
if [[ -n "${2:-}" ]]; then
    HIFI_MAX_LEN="$2"
else
    # Read N50 from the stats file already written above, avoids a second pass over the raw data.
    # -F '\t': tab separator to match seqkit stats -T output.
    # NR==2: skip the header row.
    # $13: column 13 in seqkit stats -T -a output is N50.
    RAW_N50=$(seqkit stats -T "$HIFI" | awk -F '\t' 'NR==2 {print $13}')
    HIFI_MAX_LEN=$(( RAW_N50 * 2 ))
fi

echo "[INFO] length filter: min=${HIFI_MIN_LEN} max=${HIFI_MAX_LEN}"

# Record all cutoffs for the R diagnostics script.
echo -e "species\thifi_min_len\thifi_max_len" > "${HIFI_QC_DIR}/hifi_length_cutoffs.tsv"
echo -e "${SPECIES}\t${HIFI_MIN_LEN}\t${HIFI_MAX_LEN}" >> "${HIFI_QC_DIR}/hifi_length_cutoffs.tsv"

# Raw read lengths (for diagnostics plot)
echo "[INFO] extracting RAW lengths"
# awk 'NR%4==2': FASTQ is 4-line records; line 2 of each record is the sequence. length($0) = read length in bases.
zcat "$HIFI" | awk 'NR%4==2 {print length($0)}' \
    > "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt"

# fastplong filtering
echo "[INFO] running fastplong"
# fastplong applies both length and quality filters simultaneously.
# --mean_qual: per-read mean Phred score threshold. Q20 = 99% base accuracy, appropriate for HiFi
# --out: single output FASTQ going directly to reads/processed/
conda run -p "$FASTPLONG_ENV" fastplong \
    --in "$HIFI" \
    --out "$HIFI_FILTERED" \
    --length_required "$HIFI_MIN_LEN" \
    --length_limit "$HIFI_MAX_LEN" \
    --mean_qual "$HIFI_MIN_QUAL" \
    --thread "$THREADS" \
    --html "${HIFI_QC_DIR}/fastplong_report.html" \
    --json "${HIFI_QC_DIR}/fastplong_report.json"

# Filtered stats + FastQC 
echo "[INFO] filtered stats"
# seqkit stats on the filtered FASTQ, same flags as raw for consistent columns.
seqkit stats -T "$HIFI_FILTERED" > "${HIFI_QC_DIR}/hifi_filtered_stats.tsv"
# FastQC: per-base and per-read quality overview of the filtered reads.
# -t: parallel threads for FastQC's internal processing.
# -o: output directory for the HTML report and zipped data.
fastqc -t "$THREADS" -o "$HIFI_QC_DIR" "$HIFI_FILTERED"

# Filtered read lengths
echo "[INFO] extracting FILTERED lengths"
zcat "$HIFI_FILTERED" | awk 'NR%4==2 {print length($0)}' \
    > "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt"

# Diagnostics plot
echo "[INFO] running diagnostics"
Rscript scripts/exploratory_phase/00_read_assembly/A1_reads_diagnostics.R \
    "$SPECIES" \
    "${HIFI_QC_DIR}/hifi_raw_read_lengths.txt" \
    "${HIFI_QC_DIR}/hifi_filtered_read_lengths.txt" \
    "$HIFI_QC_DIR"

echo "[OK] HiFi QC complete for $SPECIES"