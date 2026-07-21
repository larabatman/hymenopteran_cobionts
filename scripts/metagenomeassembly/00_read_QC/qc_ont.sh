#!/usr/bin/env bash
#SBATCH --job-name=ont_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/ont_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/ont_qc_%j.err

# ONT read QC: length filter only, FASTQ output preserved for myloasm.
#
# No quality filter since ONT chemistry is not confirmed per species; R10 a Q15 cutoff would be safe, while R9 median Q ~12-15; a Q15 cutoff would discard most reads.
# myloasm handles this internally at runtime
#
# Usage: sbatch qc_ont.sh <species> [max_length_override]

set -euo pipefail

SPECIES="$1"
ONT_MIN_LEN=1000

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

QC_DIR="results/${SPECIES}_stages/read_qc/ont_qc"
PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$QC_DIR" "$PROC_DIR" logs

ONT_RAW="${WORKDIR}/reads/ont/${SPECIES}/ont_raw.fastq.gz"
[[ -f "${ONT_RAW}" ]] || { echo "[ERROR] Missing ONT FASTQ: ${ONT_RAW}"; exit 1; }

ONT_FILTERED_FQ="${PROC_DIR}/ont_filtered.fastq.gz"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

module purge
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a

# 1: raw read lengths and stats

echo "[INFO] species=${SPECIES}"
echo "[INFO] PASS 1 — raw stats and read lengths"

# FASTQ records are 4 lines: header, sequence, "+", quality.
# NR%4==2 selects only sequence lines; length($0) gives the read length in bases.
zcat "${ONT_RAW}" \
    | awk 'NR%4==2 {print length($0)}' \
    > "${QC_DIR}/ont_raw_read_lengths.txt"

# -T: tab-separated output  
# -a: include all stats (N50 is column 13)
zcat "${ONT_RAW}" \
    | seqkit stats -T -a \
    > "${QC_DIR}/ont_raw_stats.tsv"

if [[ -n "${2:-}" ]]; then
    ONT_MAX_LEN="$2"
else
    RAW_N50=$(awk -F '\t' 'NR==2 {print $13}' "${QC_DIR}/ont_raw_stats.tsv")
    ONT_MAX_LEN=$(( RAW_N50 * 2 ))
fi

echo "[INFO] length filter: min=${ONT_MIN_LEN}  max=${ONT_MAX_LEN}"
echo "[INFO] NO quality filter (ONT chemistry unconfirmed — see script header)"

printf "species\tont_min_len\tont_max_len\tont_min_qual\n" \
    > "${QC_DIR}/ont_length_cutoffs.tsv"
printf "%s\t%d\t%d\tNA\n" "${SPECIES}" "${ONT_MIN_LEN}" "${ONT_MAX_LEN}" \
    >> "${QC_DIR}/ont_length_cutoffs.tsv"

# 2: filter, write FASTQ, extract filtered lengths
echo "[INFO] PASS 2 — filter and write FASTQ"

# -m/-M: length window  
# -g: drop zero-length records
zcat "${ONT_RAW}" \
    | seqkit seq \
        -m "${ONT_MIN_LEN}" \
        -M "${ONT_MAX_LEN}" \
        -g \
        --threads "${THREADS}" \
    | gzip -c > "${ONT_FILTERED_FQ}"

# Sequential: gzip has fully flushed before we read back for lengths.
zcat "${ONT_FILTERED_FQ}" \
    | awk 'NR%4==2 {print length($0)}' \
    > "${QC_DIR}/ont_filtered_read_lengths.txt"

# Filtered stats
echo "[INFO] Computing filtered stats"

FILT_COUNT=$(wc -l < "${QC_DIR}/ont_filtered_read_lengths.txt")

[[ "${FILT_COUNT}" -eq 0 ]] && {
    echo "[WARN] No reads survived filtering for ${SPECIES}. Check length cutoffs."
    exit 1
}

seqkit stats -T -a "${ONT_FILTERED_FQ}" > "${QC_DIR}/ont_filtered_stats.tsv"


# Diagnostics plot
echo "[INFO] Running diagnostics"
Rscript scripts/metagenomeassembly/helpers/reads_diagnostics.R \
    "$SPECIES" \
    "${QC_DIR}/ont_raw_read_lengths.txt" \
    "${QC_DIR}/ont_filtered_read_lengths.txt" \
    "$QC_DIR"


# Summary
RAW_READS=$(wc -l < "${QC_DIR}/ont_raw_read_lengths.txt")
echo ""
echo "[OK] ONT QC complete for ${SPECIES}"
echo " Raw reads:      ${RAW_READS}"
echo " Filtered reads: ${FILT_COUNT}"
echo " Length cutoffs: ${QC_DIR}/ont_length_cutoffs.tsv"
echo " Output FASTQ:   ${ONT_FILTERED_FQ}"
echo " QC reports:     ${QC_DIR}/"
echo " Next step: sbatch scripts/metagenomeassembly/ont/run_myloasm_ont.sh ${SPECIES}"