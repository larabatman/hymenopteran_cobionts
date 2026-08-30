#!/usr/bin/env bash
#SBATCH --job-name=hic_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hic_qc_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hic_qc_%j.err

# Clean HiC raw reads: fastp cleaning for adapter trimming, Q20 and minimum 50bp
# FastQC reports and seqkit stats
# Usage: 
# sbatch A1b_hic_qc.sh species

set -euo pipefail

SPECIES="$1"

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
# Clean directory
HIC_QC_DIR="results/${SPECIES}_stages/read_qc/hic_qc"
HIC_CLEAN_DIR="reads/hic_clean/${SPECIES}"
mkdir -p "$HIC_QC_DIR" "$HIC_CLEAN_DIR" logs
# Raw reads directory
HIC_R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
HIC_R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"
# Clean reads directory
HIC_CLEAN_R1="${HIC_CLEAN_DIR}/hic_R1.clean.fastq.gz"
HIC_CLEAN_R2="${HIC_CLEAN_DIR}/hic_R2.clean.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Modules
module purge
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-10.3.0
module load SeqKit/2.6.1

# Extract raw statistics
# Launch FastQC report 
seqkit stats -T "$HIC_R1" "$HIC_R2" > "${HIC_QC_DIR}/hic_raw_stats.tsv"
fastqc -t "$THREADS" -o "$HIC_QC_DIR" "$HIC_R1" "$HIC_R2"

# Clean reads with fastp: 
# --detect_adapter_for_pe autodetects adapters in pair end reads
# set minimal length at 50 bp, quality score at Q20, and minimal qualification over 40% of the read
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
    --html "${HIC_QC_DIR}/fastp_hic_report.html" \
    --json "${HIC_QC_DIR}/fastp_hic_report.json"

# Extract clean statistics and FastQC report
seqkit stats -T "$HIC_CLEAN_R1" "$HIC_CLEAN_R2" > "${HIC_QC_DIR}/hic_clean_stats.tsv"
fastqc -t "$THREADS" -o "$HIC_QC_DIR" "$HIC_CLEAN_R1" "$HIC_CLEAN_R2"