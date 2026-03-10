#!/usr/bin/env bash
#SBATCH --job-name=read_mapping
#SBATCH --partition=pibu_el8
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/read_mapping_%j.out
#SBATCH --error=logs/read_mapping_%j.err

# Stage A: Assembly 
# A4: Read mapping validation
# Output: assembly statistics, read mapping metrics and coverage estimates

# 1) Setup
set -euo pipefail

# Variables
SPECIES="$1"
ASM_MODE="$2"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
ASM="assemblies/hifiasm_clean/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
HIFI=reads/pacbio_hifi/${SPECIES}/*.fastq.gz
OUTDIR="results/${SPECIES}_stages/assembly_qc"
mkdir -p "$OUTDIR"

# Load modules
module purge
module load SeqKit/2.6.1
module load minimap2/2.20-GCCcore-10.3.0
module load SAMtools/1.13-GCC-10.3.0

# 2) Assembly stats
seqkit stats -T "$ASM" > "$OUTDIR/assembly_basic_stats.tsv"

# 3) Mapping
minimap2 -t 12 -ax map-hifi "$ASM" $HIFI \
  | samtools sort -@ 4 -o "$OUTDIR/reads.bam"

samtools index "$OUTDIR/reads.bam"

# Quick statistics
samtools flagstat "$OUTDIR/reads.bam" > "$OUTDIR/mapping_flagstat.txt"
samtools stats "$OUTDIR/reads.bam" > "$OUTDIR/mapping_stats.txt"

# 4) Per-contig coverage computation
samtools coverage "$OUTDIR/reads.bam" > "$OUTDIR/coverage_per_contig.tsv"

# Start position is on column 2 and end position on column 3 is multiplied by per-contig coverage on column 7 to obtain the total coverag, which is then divided by the total length to obtain the mean 
# -v takes the shell variable $SPECIES inside awk as sp
# Set output field as separator for species mean_depth
# Compute interval between start and end with inclusive coordinates + 1
# Add each contig length to the running total with total_len += len
# Build length-weighted coverage sum: if length = 1000 and depth = 20, contribution = 20000
# Finally, weighted mean cov = sum(length*depth)/sum(length)
awk -v sp="$SPECIES" '
BEGIN{
  OFS="\t"
}
NR>1{
  len = $3 - $2 + 1
  total_len += len
  total_cov += len * $7
}
END{
  print "species", "mean_depth"
  if (total_len > 0) {
    print sp, total_cov / total_len
  } else {
    print sp, "NA"
  }
}' "$OUTDIR/coverage_per_contig.tsv" > "$OUTDIR/coverage_summary.tsv"

echo "[INFO] Assembly QC complete for $SPECIES"