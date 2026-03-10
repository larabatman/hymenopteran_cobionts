#!/usr/bin/env bash
#SBATCH --job-name=map_hic
#SBATCH --partition=pibu_el8
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --chdir=/data/projects/p2025-0083_mining_cobionts
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/map_hic_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/map_hic_%j.err

# Stage C: Orthogonal evidence layer
# C2: Hi-C alignement, mapping Hi-C reads to assembly for contact evidence
# Align Hi-C reads to assembly, clean alignements and remove PCR duplicates to produce final BAM
# Inputs: HiC R1 and R2, assembly

# 1) Setup 
set -euo pipefail


SPECIES="$1"

# Inputs
ASM="assemblies/hifiasm_clean/${SPECIES}/asm.hic.p_ctg.fasta"
R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"

# Outputs
OUTDIR="results/${SPECIES}_stages/hic_map"
mkdir -p "$OUTDIR"

# Load modules
module purge
module load BWA-MEM2/2.2.1-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0

# 2) Index assembly 
bwa-mem2 index "$ASM"

THREADS="${SLURM_CPUS_PER_TASK}"

# 3) Mapping Hi-C reads to assembly
# Filter for primary reads, removing secondary and supplementary reads and sorting by name
# -F 256 flags secondary elignment, 2048 flags supplementary alignment
# Keep primry alignements only 
# sort -n sort reads by read name such that mates are adjacent
bwa-mem2 mem -t "$THREADS" "$ASM" "$R1" "$R2" \
  | samtools view -@ "$THREADS" -b -F 2304 \
  | samtools sort -@ "$THREADS" -n -o "$OUTDIR/hic.namesort.bam" -

# fixmate: fills in mate coordinates and insert size fields
# adds metadata on mate position insert size pair flags for duplicate marking
samtools fixmate -@ "$THREADS" -m "$OUTDIR/hic.namesort.bam" "$OUTDIR/hic.fixmate.bam"

# Position sort for markdup to remove duplicates 
# samtools sort by genomic position, required for duplicate dtection as well
samtools sort -@ "$THREADS" -o "$OUTDIR/hic.possort.bam" "$OUTDIR/hic.fixmate.bam"
samtools index -@ "$THREADS" "$OUTDIR/hic.possort.bam"

# Remove duplicates 
# samtools markdup -r removes PCR duplicates to avoid artificailly inflating contact counts
samtools markdup -@ "$THREADS" -r "$OUTDIR/hic.possort.bam" "$OUTDIR/hic.filtered.bam"
samtools index -@ "$THREADS" "$OUTDIR/hic.filtered.bam"

# Quick statistics
samtools flagstat "$OUTDIR/hic.filtered.bam" > "$OUTDIR/hic.filtered.flagstat.txt"

echo "[INFO] C2 done for $SPECIES"