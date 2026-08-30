#!/usr/bin/env bash
#SBATCH --job-name=read_mapping
#SBATCH --partition=pibu_el8
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/read_mapping_%j.out
#SBATCH --error=logs/read_mapping_%j.err

# This script maps the reads back to the hifiasm assembly with minimap 2
# Usage
# sbatch A3_read_mapping.sh species ASM_MODE==bp or ASM_MODE==hic dependong on HiC availability

set -euo pipefail
# Arguments:
SPECIES="$1"
ASM_MODE="$2"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Working directory:
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
# Assembly directory
ASM="assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
# Reads directory
HIFI="results/archive_exploratory_phase/${SPECIES}_stages/read_qc/hifi_qc/hifi_filtered/hifi.filtered.fastq.gz"
OUTDIR="results/${SPECIES}_stages/assembly_qc"
mkdir -p "$OUTDIR"
# Modules
module purge
module load minimap2/2.20-GCCcore-10.3.0 SAMtools/1.13-GCC-10.3.0

# minimap2: 
# -a for SAM output, text alignment for samtools
# -x map-hifi: present for HiFi reads 
# Piped directily into samtools to sort that orders the alignments by cooredinate
minimap2 -t "$THREADS" -ax map-hifi "$ASM" "$HIFI" | samtools sort -@ 4 -o "$OUTDIR/reads.bam"
# Index: .bai file to go directly to given contigs on a sorted BAM file 
samtools index "$OUTDIR/reads.bam"
# samtools statistics: 
samtools flagstat "$OUTDIR/reads.bam" > "$OUTDIR/mapping_flagstat.txt"
samtools stats "$OUTDIR/reads.bam" > "$OUTDIR/mapping_stats.txt"
# samtools coverage
samtools coverage "$OUTDIR/reads.bam" > "$OUTDIR/coverage_per_contig.tsv"