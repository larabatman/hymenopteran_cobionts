#!/usr/bin/env bash
#SBATCH --job-name=hifi_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hifi_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hifi_qc_%j.err

# qc_hifi.sh
# HiFi read QC: fastplong filtering, FastQC, seqkit stats, R diagnostics
#
# OUTPUT FORMAT: FASTQ (not FASTA)
# Keep FASTQ so run_myloasm.sh can use quality scores for SNPmer calling and polishing 
# FASTA conversion for the pipeline is done at the end of run_myloasm.sh, matching the CLR and ONT workflows.
#
# Input: reads/pacbio_hifi/<species>/*.fastq.gz
# Output: 
# reads/processed/<species>/hifi_filtered.fastq.gz (pipeline-ready)
# results/<species>_stages/read_qc/hifi_qc/ (QC reports)
#
# Usage:  qc_hifi.sh <species> [max_length_override]


set -euo pipefail

SPECIES="$1"

# Min length: 3 kbp for HiFi
# Min quality: Q20 mean per read
# per-read mean Phred scores are accurate.
HIFI_MIN_LEN=3000
HIFI_MIN_QUAL=20

echo "[INFO] species=${SPECIES}"

# Paths 
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

QC_DIR="results/${SPECIES}_stages/read_qc/hifi_qc"
PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$QC_DIR" "$PROC_DIR" logs

# Find all raw HiFi FASTQs for this species
HIFI=reads/pacbio_hifi/${SPECIES}/*.fastq.gz

# FASTQ output: quality scores preserved for myloasm
HIFI_FILTERED_FQ="${PROC_DIR}/hifi_filtered.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Modules
module purge
module load FastQC/0.11.9-Java-11
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a
module load Anaconda3/2022.05

source "$(conda info --base)/etc/profile.d/conda.sh"
FASTPLONG_ENV="${WORKDIR}/.conda_envs/fastplong"

# Raw stats
echo "[INFO] Raw stats"

# -T: tab-separated output for consistent awk parsing downstream.
# -a: include ALL statistics, column 13 is N50, used for the max cutoff.
seqkit stats -T -a $HIFI > "${QC_DIR}/hifi_raw_stats.tsv"

# Max length from N50
if [[ -n "${2:-}" ]]; then
    HIFI_MAX_LEN="$2"
else
    # Read N50 from the stats file already written above, avoids a second pass over the raw data.
    # -F '\t': tab separator to match seqkit stats -T output.
    # NR==2: skip the header row.
    # $13: column 13 in seqkit stats -T -a output is N50.
    RAW_N50=$(awk -F '\t' 'NR==2 {print $13}' "${QC_DIR}/hifi_raw_stats.tsv")
    HIFI_MAX_LEN=$(( RAW_N50 * 2 ))
fi

echo "[INFO] length filter : min=${HIFI_MIN_LEN}  max=${HIFI_MAX_LEN}"
echo "[INFO] quality filter: avg Q >= ${HIFI_MIN_QUAL}"

# Record all cutoffs for the R diagnostics script.
echo -e "species\thifi_min_len\thifi_max_len\thifi_min_qual" \
    > "${QC_DIR}/hifi_length_cutoffs.tsv"
echo -e "${SPECIES}\t${HIFI_MIN_LEN}\t${HIFI_MAX_LEN}\t${HIFI_MIN_QUAL}" \
    >> "${QC_DIR}/hifi_length_cutoffs.tsv"

# Raw read lengths (for diagnostics plot)
echo "[INFO] Extracting raw lengths"

# awk 'NR%4==2': FASTQ is 4-line records; line 2 of each record is the sequence. length($0) = read length in bases.
zcat $HIFI | awk 'NR%4==2 {print length($0)}' \
    > "${QC_DIR}/hifi_raw_read_lengths.txt"

# fastplong filtering
echo "[INFO] Running fastplong"

# fastplong applies both length and quality filters simultaneously.
# --mean_qual: per-read mean Phred score threshold. Q20 = 99% base accuracy, appropriate for HiFi
# --out: single output FASTQ going directly to reads/processed/
conda run -p "$FASTPLONG_ENV" fastplong \
    --in $HIFI \
    --out "$HIFI_FILTERED_FQ" \
    --length_required "$HIFI_MIN_LEN" \
    --length_limit "$HIFI_MAX_LEN" \
    --mean_qual "$HIFI_MIN_QUAL" \
    --thread "$THREADS" \
    --html "${QC_DIR}/fastplong_report.html" \
    --json "${QC_DIR}/fastplong_report.json"

# Filtered stats + FastQC 
echo "[INFO] Filtered stats and FastQC"

# seqkit stats on the filtered FASTQ, same flags as raw for consistent columns.
seqkit stats -T -a "$HIFI_FILTERED_FQ" > "${QC_DIR}/hifi_filtered_stats.tsv"

# FastQC: per-base and per-read quality overview of the filtered reads.
# -t: parallel threads for FastQC's internal processing.
# -o: output directory for the HTML report and zipped data.
fastqc -t "$THREADS" -o "$QC_DIR" "$HIFI_FILTERED_FQ"

# Filtered read lengths
echo "[INFO] Extracting filtered lengths"
zcat "$HIFI_FILTERED_FQ" | awk 'NR%4==2 {print length($0)}' \
    > "${QC_DIR}/hifi_filtered_read_lengths.txt"

# Diagnostics plot
echo "[INFO] Running diagnostics"

Rscript scripts/metagenomeassembly/helpers/reads_diagnostics.R \
    "$SPECIES" \
    "${QC_DIR}/hifi_raw_read_lengths.txt" \
    "${QC_DIR}/hifi_filtered_read_lengths.txt" \
    "$QC_DIR"

echo ""
echo "[OK] HiFi QC complete for ${SPECIES}"
echo " Output FASTQ:  ${HIFI_FILTERED_FQ}"
echo " QC reports: ${QC_DIR}/"
echo " Next step: sbatch scripts/metagenomeassembly/hifi/run_myloasm.sh ${SPECIES}"