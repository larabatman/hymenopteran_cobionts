#!/usr/bin/bash
#SBATCH --job-name=reads_qc
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/reads_qc_%j.out
#SBATCH --error=logs/reads_qc_%j.err

# Stage A: Assembly 
# A1: General read QC
# Overall HiFi and Hi-C statistics
# Output: fastqc reports

set -euo pipefail

# 1) Setup
SPECIES="$1"
ASM_MODE="$2"

echo "[INFO] job started on $(hostname)"
echo "[INFO] species=$SPECIES"
echo "[INFO] assembly mode=$ASM_MODE"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

QC_DIR="results/${SPECIES}_stages/read_qc"
HIC_CLEAN="reads/hic_clean/${SPECIES}"

mkdir -p "$QC_DIR" "$HIC_CLEAN" logs

HIFI_FILES=(reads/pacbio_hifi/${SPECIES}/*.fastq.gz)

HIC_R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
HIC_R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"

THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "[INFO] Using $THREADS threads"

# Load modules
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-10.3.0
module load SeqKit/2.6.1

# 2) HiFi QC
echo "[INFO] Running HiFi QC"

seqkit stats "${HIFI_FILES[@]}" > "$QC_DIR/hifi_stats.tsv"

fastqc \
  -t "$THREADS" \
  -o "$QC_DIR" \
  "${HIFI_FILES[@]}"

seqkit fx2tab -l "${HIFI_FILES[@]}" \
| awk '{print $2}' \
| sort -n \
> "$QC_DIR/hifi_read_lengths.txt"

# End of script for ASM_MODE=="bp"
if [[ "$ASM_MODE" == "bp" ]]; then
    echo "[INFO] ASM_MODE=bp → skipping Hi-C QC"
    echo "[INFO] HiFi QC completed successfully"
    exit 0
fi

# 3) Hi-C QC
echo "[INFO] Running Hi-C raw QC"

fastqc \
  -t "$THREADS" \
  -o "$QC_DIR" \
  "$HIC_R1" "$HIC_R2"

# Optional trim, but dangerous
# fastp can desynchronize pairs → disabled for now

# fastp \
#   -i "$HIC_R1" \
#   -I "$HIC_R2" \
#   -o "$HIC_CLEAN/hic_R1.clean.fastq.gz" \
#   -O "$HIC_CLEAN/hic_R2.clean.fastq.gz" \
#   --detect_adapter_for_pe \
#   --trim_poly_g \
#   --thread "$THREADS" \
#   --html "$QC_DIR/fastp_report.html" \
#   --json "$QC_DIR/fastp_report.json"

#echo "[INFO] Running Hi-C QC after trimming"

#fastqc \
#  -t "$THREADS" \
#  -o "$QC_DIR" \
#  "$HIC_CLEAN/hic_R1.clean.fastq.gz" \
#  "$HIC_CLEAN/hic_R2.clean.fastq.gz"

# Done 
echo "[INFO] Read QC completed successfully"