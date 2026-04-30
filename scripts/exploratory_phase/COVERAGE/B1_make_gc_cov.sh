#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=gc_cov
#SBATCH --output=logs/gc_cov_%j.out
#SBATCH --error=logs/gc_cov_%j.err

# =============================================================================
# Stage B1: Compute per-contig mean coverage
# Coverage-based partition for host backbone definition
# Define host backbone with coverage as primary discriminator 
# Compute mean_cov per contig

# This script computes per-contig length and GC content from the assembly fasta, retrieves mean coverage and merges them into a single table: gc_cov.tsv

set -euo pipefail 

# Setup 
# Project root
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

# Inputs
# Sample identifier 
SPECIES="$1"
# Assembly mode: either Hi-C information (hic) or not (bp)
ASM_MODE="$2"
# Path to primary contigs assembly FASTA, reference for extracting contig GC/ length stats and mapping reads to compute coverage
ASM="${WORKDIR}/assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
# Path to BAM produced in A4
BAM="${WORKDIR}/results/${SPECIES}_stages/assembly_qc/reads.bam"

# Outputs 
OUTDIR="${WORKDIR}/results/${SPECIES}_stages/gc_cov"
mkdir -p "$OUTDIR"

echo "[INFO] Assembly: $ASM"
echo "[INFO] BAM reads: $BAM"
echo "[INFO] Output dir: $OUTDIR"

# Load modules
# SeqKit is used to extract contig states (name, length, GC) from FASTA
module load SeqKit/2.6.1
# samtools: used to sort/ index BAM and compute coverage per contig
module load SAMtools/1.13-GCC-10.3.0

# Compute GC and length per contig
# SeqKit fx2tab outputs tabular stats per FASTA record:
# -n: contig name
# -l: length
# -g: GC content
# The awk command sets output field separator to TAB, prints a header once and prints 3 columns per contig. 
# The output file gc_len.tsv looks like: contig len gc
seqkit fx2tab -n -l -g "$ASM" \
  | awk 'BEGIN{OFS="\t"; print "contig","len","gc"} {print $1,$2,$3}' \
  > "$OUTDIR/gc_len.tsv"

# Compute per contig mean coverage
# Per-contig coverage with samtools
# samtools coverage summarizes coverage per reference sequence (per contig here). It output a header line and then one line per contig
# The column 7 is the mean depth, and the output file is the raw samtools coverage output coverage.tsv
samtools coverage "$BAM" > "$OUTDIR/coverage.tsv"

# Extract only the contig name and mean depth:
# The awk logic prints the header "contig mean_cov", skips the first line which is the coverage.tsv header and prints contigs and mean depth
# The output file cov.tsv contains: contig  mean_cov
awk 'BEGIN{OFS="\t"; print "contig","mean_cov"}
     NR==1{next}
     {print $1,$7}' "$OUTDIR/coverage.tsv" > "$OUTDIR/cov.tsv"

# Merge GC and length with coverage into gc_cov.tsv
# Join gc_len and cov.tsv by contig name which is column 1 in both. 
# The first awk reads gc_len.tsv (NR==FNR is true only for the first file). It then stors, for each contig name, a[contig] = len<TAB>gc , then reads cov.tsv by skipping its header and printin conitg, (len, gc), mean_cov
# The second awk adds the final header row. 
# The final output file gc_civ.tsv has columns: contig  len gc  mean_cov
awk 'BEGIN{OFS="\t"}
     NR==FNR{a[$1]=$2 OFS $3; next}
     FNR==1{next}
     ($1 in a){print $1, a[$1], $2}' \
  "$OUTDIR/gc_len.tsv" "$OUTDIR/cov.tsv" \
| awk 'BEGIN{OFS="\t"; print "contig","len","gc","mean_cov"} {print}' \
> "$OUTDIR/gc_cov.tsv"

# Done
echo "[INFO] Wrote: $OUTDIR/gc_cov.tsv"
