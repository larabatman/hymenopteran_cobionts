#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=gc_cov
#SBATCH --output=logs/gc_cov_%j.out
#SBATCH --error=logs/gc_cov_%j.err

# Create a single table with contig, lentgth, GC and mean coverage for each contig
# For each contig, compute its length and GC content from the assembly fasta and mearge into a single table gc_cov.tsv
# Usage:
# sbatch B1_make_gc_cov.sh species ASM_MODE==bp or ASM_MODE==hic

set -euo pipefail 

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

# Arguments: 
SPECIES="$1"
ASM_MODE="$2"
# Primary contigs assembly fasta
ASM="${WORKDIR}/assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
# BAM file
BAM="${WORKDIR}/results/${SPECIES}_stages/assembly_qc/reads.bam"
# Coverage
COV="results/${SPECIES}_stages/assembly_qc/coverage_per_contig.tsv"
# Outputs 
OUTDIR="${WORKDIR}/results/${SPECIES}_stages/gc_cov"
mkdir -p "$OUTDIR"

# Modules
module load SeqKit/2.6.1
module load SAMtools/1.13-GCC-10.3.0

# GC and length for each contig: seqkit fx2tab
# -n for contig name
# -l for length
# -g for GC content
# Piped to awk: "\t" tab field separator, print the header contig len gc and the three output columns for each contig
# gc_len.tsv: contig  len gc
seqkit fx2tab -n -l -g "$ASM" | awk 'BEGIN{OFS="\t"; print "contig","len","gc"} {print $1,$2,$3}'> "$OUTDIR/gc_len.tsv"

# From samtools coverage: extract contig name and mean depth
# awk: print header tab separated contig  mean_cov and print the two columns
awk 'BEGIN{OFS="\t"; print "contig","mean_cov"} NR==1{next} {print $1,$7}' ${COV} > "$OUTDIR/cov.tsv"

# Merge gc_len.tsv with cov.tsv by contig name, column 1 in both cases
# First read gc_len.tsv: NR==NFR is true only for the first file 
# For each contig name: store name[contig] = len tab gc
# Read cov.tsv
# Print contig  len gc  mean_cov
awk 'BEGIN{OFS="\t"} NR==FNR{name[$1]=$2 OFS $3; next} FNR==1{next} ($1 in name){print $1, name[$1], $2}' "$OUTDIR/gc_len.tsv" "$OUTDIR/cov.tsv" \
  | awk 'BEGIN{OFS="\t"; print "contig","len","gc","mean_cov"} {print}' > "$OUTDIR/gc_cov.tsv"