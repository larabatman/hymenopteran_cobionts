#!/usr/bin/env bash
#SBATCH --job-name=hic_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hic_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hic_qc_%j.err

# =============================================================================
# Stage A1b: Hi-C read QC
# fastp cleaning (adapter trim, Q20, min 50bp) + FastQC pre/post + seqkit stats
#
# Usage: sbatch A1b_hic_qc.sh <species>

set -euo pipefail

SPECIES="$1"

echo "[INFO] species=$SPECIES"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

HIC_QC_DIR="results/${SPECIES}_stages/read_qc/hic_qc"
HIC_CLEAN_DIR="reads/hic_clean/${SPECIES}"
mkdir -p "$HIC_QC_DIR" "$HIC_CLEAN_DIR" logs

HIC_R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
HIC_R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"

[[ -s "$HIC_R1" ]] || { echo "[ERROR] Missing Hi-C R1: $HIC_R1" >&2; exit 1; }
[[ -s "$HIC_R2" ]] || { echo "[ERROR] Missing Hi-C R2: $HIC_R2" >&2; exit 1; }

HIC_CLEAN_R1="${HIC_CLEAN_DIR}/hic_R1.clean.fastq.gz"
HIC_CLEAN_R2="${HIC_CLEAN_DIR}/hic_R2.clean.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Modules
module purge
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-10.3.0
module load SeqKit/2.6.1

# RAW stats + FastQC
echo "[INFO] RAW Hi-C stats"
seqkit stats -T "$HIC_R1" "$HIC_R2" > "${HIC_QC_DIR}/hic_raw_stats.tsv"
fastqc -t "$THREADS" -o "$HIC_QC_DIR" "$HIC_R1" "$HIC_R2"

# fastp cleaning
echo "[INFO] running fastp"
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

# Clean stats + FastQC
echo "[INFO] CLEAN Hi-C stats"
seqkit stats -T "$HIC_CLEAN_R1" "$HIC_CLEAN_R2" > "${HIC_QC_DIR}/hic_clean_stats.tsv"
fastqc -t "$THREADS" -o "$HIC_QC_DIR" "$HIC_CLEAN_R1" "$HIC_CLEAN_R2"

echo "[OK] Hi-C QC complete for $SPECIES"