#!/usr/bin/env bash
#SBATCH --job-name=illumina_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/illumina_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/illumina_qc_%j.err

# qc_illumina.sh
# Illumina paired-end WGS read QC
# fastp filtering:
#   - Adapter trimming: auto-detection (--detect_adapter_for_pe)
#   - Quality: mean Q20 per read
#   - Length: minimum 50bp after trimming (reads shorter than this produce no useful overlaps in metaSPAdes assembly graphs)
#   - PolyG trimming: enabled
#   - No maximum length filter (Illumina reads are uniform)
#
# Input:  reads/illumina/<species>/illumina_raw_R1/R2.fastq.gz
# Output: reads/processed/<species>/illumina_filtered_R1/R2.fastq.gz
#         results/<species>_stages/read_qc/illumina_qc/
#
# Usage:  sbatch illumina_qc.sh <species> 

set -euo pipefail

SPECIES="$1"

echo "[INFO] species=${SPECIES}"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

QC_DIR="results/${SPECIES}_stages/read_qc/illumina_qc"
PROC_DIR="reads/processed/${SPECIES}"
mkdir -p "$QC_DIR" "$PROC_DIR" logs

# Input: merged raw Illumina FASTQ pair from download_illumina_sra.sh.
RAW_R1="${WORKDIR}/reads/illumina/${SPECIES}/illumina_raw_R1.fastq.gz"
RAW_R2="${WORKDIR}/reads/illumina/${SPECIES}/illumina_raw_R2.fastq.gz"

[[ -f "${RAW_R1}" ]] || { echo "[ERROR] Missing R1: ${RAW_R1}. Run download_illumina_sra.sh first."; exit 1; }
[[ -f "${RAW_R2}" ]] || { echo "[ERROR] Missing R2: ${RAW_R2}. Run download_illumina_sra.sh first."; exit 1; }

# Final pipeline-ready outputs: filtered R1 and R2.
# Both files are needed by metaSPAdes for paired-end assembly.
FILTERED_R1="${PROC_DIR}/illumina_filtered_R1.fastq.gz"
FILTERED_R2="${PROC_DIR}/illumina_filtered_R2.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Modules
module purge
module load fastp/0.23.4-GCC-10.3.0
module load SeqKit/2.6.1
module load R/4.2.1-foss-2021a

# Raw stats
echo "[INFO] RAW stats"
# -T: tab-separated; -a: all statistics including N50 and Q20/Q30 rates.
# We run stats on R1 only for speed (R2 stats are similar and the read count is what matters for the diagnostic report)
seqkit stats -T -a "${RAW_R1}" > "${QC_DIR}/illumina_raw_stats_R1.tsv"
seqkit stats -T -a "${RAW_R2}" > "${QC_DIR}/illumina_raw_stats_R2.tsv"

# Raw read lengths
echo "[INFO] Extracting RAW lengths"
# For Illumina reads NR%4==2 extracts sequence lines from FASTQ.
# We use R1 only for the length distribution, all reads are the same length (202bp for HiSeq 4000) so R1 is representative.
zcat "${RAW_R1}" | awk 'NR%4==2 {print length($0)}' \
    > "${QC_DIR}/illumina_raw_read_lengths.txt"

# fastp filtering 
echo "[INFO] Running fastp QC..."


# -i/-I: input R1 and R2 (paired-end)
# -o/-O: output filtered R1 and R2
# --detect_adapter_for_pe: auto-detect adapter sequences for paired-end data. fastp analyses read overlap to identify adapters without needing the adapter sequence to be specified
# -q 20: per-base quality threshold. Bases below Q20 are considered unqualified. Q20 = 99% base accuracy
# -u 40: allow up to 40% of bases per read to be below Q20 before discarding the read. fastp default; keeps reads with a few low-quality bases rather than discarding aggressively.
# -l 50: discard reads shorter than 50bp after trimming. metaSPAdes needs reads long enough to form assembly graph overlaps; 50bp is the consensus minimum for metagenome assembly.
# --cut_right: sliding window trimming from the 3' end. Removes the common Illumina quality drop-off at the end of reads. We use cut_right rather than cut_tail (stricter) to preserve more read length while still removing bad 3' bases.
# --cut_mean_quality 20: quality threshold for the sliding window. A window is cut when mean quality drops below Q20, consistent with -q 20.
# --trim_poly_g: remove polyG tails.
# -w: number of worker threads for parallel processing.
# --html / --json: quality reports consistent with fastplong for HiFi data.
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
    --json "${QC_DIR}/fastp_report.json" \
    ${TENX_FLAGS}

# Filtered stats
echo "[INFO] Filtered stats"
seqkit stats -T -a "${FILTERED_R1}" > "${QC_DIR}/illumina_filtered_stats_R1.tsv"
seqkit stats -T -a "${FILTERED_R2}" > "${QC_DIR}/illumina_filtered_stats_R2.tsv"

# Filtered read lengths
echo "[INFO] Extracting FILTERED lengths"
zcat "${FILTERED_R1}" | awk 'NR%4==2 {print length($0)}' \
    > "${QC_DIR}/illumina_filtered_read_lengths.txt"

# Record QC cutoffs for consistency with other data types
echo -e "species\tmin_qual\tmin_len\tadapter_detection\tpoly_g_trim" \
    > "${QC_DIR}/illumina_qc_cutoffs.tsv"
echo -e "${SPECIES}\t20\t50\tauto\ttrue" \
    >> "${QC_DIR}/illumina_qc_cutoffs.tsv"

# Summary
RAW_READS=$(wc -l < "${QC_DIR}/illumina_raw_read_lengths.txt")
FILT_READS=$(wc -l < "${QC_DIR}/illumina_filtered_read_lengths.txt")
echo ""
echo "[OK] Illumina QC complete for ${SPECIES}"
echo " Raw read pairs: ${RAW_READS}"
echo " Filtered read pairs: ${FILT_READS}"
echo " QC reports: ${QC_DIR}/"
echo " Next step: sbatch run_metaspades.sh ${SPECIES}"