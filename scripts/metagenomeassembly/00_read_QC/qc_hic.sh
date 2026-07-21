#!/usr/bin/env bash
#SBATCH --job-name=hic_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hic_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hic_qc_%j.err

# Hi-C read QC
#
# fastp cleaning (adapter trim, Q20, min 50 bp), FastQC pre/post, seqkit stats
#
# Input: reads/hic/<species>/*_1.fastq.gz, *_2.fastq.gz (from download_reads.sh)
# Output:
# reads/hic_clean/<species>/hic_R1.clean.fastq.gz (for CRAM conversion)
# reads/hic_clean/<species>/hic_R2.clean.fastq.gz
# results/<species>_stages/read_qc/hic_qc/ (QC reports)
#
# Usage: sbatch A1b_hic_qc.sh <species>

set -euo pipefail

SPECIES="$1"

echo "[INFO] species=${SPECIES}"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

QC_DIR="results/${SPECIES}_stages/read_qc/hic_qc"
HIC_CLEAN_DIR="reads/hic_clean/${SPECIES}"
mkdir -p "$QC_DIR" "$HIC_CLEAN_DIR" logs

# Find Hi-C paired reads (ENA names: *_1.fastq.gz / *_2.fastq.gz)
HIC_R1=$(ls reads/hic/${SPECIES}/*_1.fastq.gz 2>/dev/null | head -1) \
    || { echo "[ERROR] No Hi-C R1 in reads/hic/${SPECIES}/"; exit 1; }
HIC_R2=$(ls reads/hic/${SPECIES}/*_2.fastq.gz 2>/dev/null | head -1) \
    || { echo "[ERROR] No Hi-C R2 in reads/hic/${SPECIES}/"; exit 1; }

echo "[INFO] R1: ${HIC_R1}"
echo "[INFO] R2: ${HIC_R2}"

HIC_CLEAN_R1="${HIC_CLEAN_DIR}/hic_R1.clean.fastq.gz"
HIC_CLEAN_R2="${HIC_CLEAN_DIR}/hic_R2.clean.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Modules
module purge
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-10.3.0
module load SeqKit/2.6.1

# Raw stats + FastQC
echo "[INFO] RAW Hi-C stats"
# Summary table
seqkit stats -T "$HIC_R1" "$HIC_R2" > "${QC_DIR}/hic_raw_stats.tsv"
fastqc -t "$THREADS" -o "$QC_DIR" "$HIC_R1" "$HIC_R2"

# fastp cleaning
echo "[INFO] running fastp"
# in2, in2, out1, out2 for clean input and output
# --detect_adapter_for_pe to autodetect the adapters for paired-end reads and remove them during trimming
# --qualitfied_quality_phred 20: defines base quality at Q20 rather than default Q15
# --unqualitfied_percent_limit 40: a read with more than 40% of its bases below Q20 are removed
# --length_required 50: discards reads shorter than 50 after trimming

fastp \
    --in1 "$HIC_R1" \
    --in2 "$HIC_R2" \
    --out1 "$HIC_CLEAN_R1" \
    --out2 "$HIC_CLEAN_R2" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 50 \
    --thread "$THREADS" \
    --html "${QC_DIR}/fastp_hic_report.html" \
    --json "${QC_DIR}/fastp_hic_report.json"

# Clean stats + FastQC
echo "[INFO] CLEAN Hi-C stats"
seqkit stats -T "$HIC_CLEAN_R1" "$HIC_CLEAN_R2" > "${QC_DIR}/hic_clean_stats.tsv"
fastqc -t "$THREADS" -o "$QC_DIR" "$HIC_CLEAN_R1" "$HIC_CLEAN_R2"

echo ""
echo "[OK] Hi-C QC complete for ${SPECIES}"
echo " Clean FASTQ: ${HIC_CLEAN_R1}, ${HIC_CLEAN_R2}"
echo " QC reports: ${QC_DIR}/"