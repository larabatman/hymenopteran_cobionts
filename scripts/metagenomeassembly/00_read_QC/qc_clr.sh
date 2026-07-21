#!/usr/bin/env bash
#SBATCH --job-name=clr_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/clr_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/clr_qc_%j.err

# PacBio CLR subread QC
# Filters by length only (no quality filter, CLR Phred scores are unreliable).
# min: 1 kbp, reads below this cannot form useful overlaps
# max: 2× N50, removes outliers likely to be chimeric or contaminated
#
# Usage: sbatch Aq_clr.sh <species> [max_length_override]

set -euo pipefail

SPECIES="$1"
CLR_MIN_LEN=1000

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

QC_DIR="results/${SPECIES}_stages/read_qc/clr_qc"
PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$QC_DIR" "$PROC_DIR" logs

CLR_RAW="${WORKDIR}/reads/pacbio_clr/${SPECIES}/clr_raw.fastq.gz"
[[ -f "${CLR_RAW}" ]] || { echo "[ERROR] Missing CLR FASTQ: ${CLR_RAW}"; exit 1; }

CLR_FILTERED_FQ="${PROC_DIR}/clr_filtered.fastq.gz"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Modules
module purge
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a

# 1: Raw read lengths and stats
echo "[INFO] species=${SPECIES}"
echo "[INFO] PASS 1 — raw stats and read lengths"

# Select sequences, discard Phred scores
# FASTQ records are 4 lines: header, sequence, "+", quality
# NR%4==2 selects only sequence lines; length($0) gives the read length in bases.
zcat "${CLR_RAW}" \
    | awk 'NR%4==2 {print length($0)}' \
    > "${QC_DIR}/clr_raw_read_lengths.txt"

# Get raw statistics 
# -T: tab-separated output  
# -a: include all stats (N50 is column 13)
zcat "${CLR_RAW}" \
    | seqkit stats -T -a \
    > "${QC_DIR}/clr_raw_stats.tsv"

# Get length filter to apply: either given argument, or default 2XN50
if [[ -n "${2:-}" ]]; then
    CLR_MAX_LEN="$2"
else
    RAW_N50=$(awk -F '\t' 'NR==2 {print $13}' "${QC_DIR}/clr_raw_stats.tsv")
    CLR_MAX_LEN=$(( RAW_N50 * 2 ))
fi

echo "[INFO] length filter: min=${CLR_MIN_LEN}  max=${CLR_MAX_LEN}"
echo "[INFO] NO quality filter (CLR — Phred scores unreliable)"

# Get a record of the length cutoffs
# Using printf for tab separation explicit
printf "species\tclr_min_len\tclr_max_len\tclr_min_qual\n" \
    > "${QC_DIR}/clr_length_cutoffs.tsv"
printf "%s\t%d\t%d\tNA\n" "${SPECIES}" "${CLR_MIN_LEN}" "${CLR_MAX_LEN}" \
    >> "${QC_DIR}/clr_length_cutoffs.tsv"

# 2: filter, write FASTQ, extract filtered lengths
echo "[INFO] PASS 2 — filter and write FASTQ"

# Filter length with seqkit seq by cutting the reads directly
# -m and-M: length window  
# -g: drop zero-length records
zcat "${CLR_RAW}" \
    | seqkit seq \
        -m "${CLR_MIN_LEN}" \
        -M "${CLR_MAX_LEN}" \
        -g \
        --threads "${THREADS}" \
    | gzip -c > "${CLR_FILTERED_FQ}"

# Sequential step: use gzip before we read back for lengths, to avoid file corruption
zcat "${CLR_FILTERED_FQ}" \
    | awk 'NR%4==2 {print length($0)}' \
    > "${QC_DIR}/clr_filtered_read_lengths.txt"

# Filtered stats
echo "[INFO] Computing filtered stats"

FILT_COUNT=$(wc -l < "${QC_DIR}/clr_filtered_read_lengths.txt")

[[ "${FILT_COUNT}" -eq 0 ]] && {
    echo "[WARN] No reads survived filtering for ${SPECIES}. Check length cutoffs."
    exit 1
}

# Filtered statistics
seqkit stats -T -a "${CLR_FILTERED_FQ}" > "${QC_DIR}/clr_filtered_stats.tsv"

# Diagnostics plot
echo "[INFO] Running diagnostics"
Rscript scripts/metagenomeassembly/helpers/reads_diagnostics.R \
    "$SPECIES" \
    "${QC_DIR}/clr_raw_read_lengths.txt" \
    "${QC_DIR}/clr_filtered_read_lengths.txt" \
    "$QC_DIR"

# Summary
RAW_READS=$(wc -l < "${QC_DIR}/clr_raw_read_lengths.txt")
echo ""
echo "[OK] CLR QC complete for ${SPECIES}"
echo " Raw reads: ${RAW_READS}"
echo " Filtered reads: ${FILT_COUNT}"
echo " Length cutoffs: ${QC_DIR}/clr_length_cutoffs.tsv"
echo " Output FASTQ: ${CLR_FILTERED_FQ}"
echo " QC reports: ${QC_DIR}/"
echo " Next step: sbatch scripts/metagenomeassembly/01_assembly/run_metaflye.sh ${SPECIES}"