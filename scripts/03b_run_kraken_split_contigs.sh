#!/usr/bin/env bash
#SBATCH --job-name=kraken_windowed
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/kraken_windowed_%j.out
#SBATCH --error=logs/kraken_windowed_%j.err

# This script performs window-based taxonomic classification of an assembled genome using Kraken2. 
# The aim is to detect host, cobiont and mixed signals at sub-contig resolution. 


set -euo pipefail

# Sample identifier: 
SAMPLE_ID="Bombus_terrestris"

# Input
# GFA file produced by hifiasm with contig sequences in "S" (segment) lines
GFA="assemblies/hifiasm/${SAMPLE_ID}/asm.hic.p_ctg.gfa"
# FASTA version of the primary contigs generated from the GFA
FASTA="assemblies/hifiasm/${SAMPLE_ID}/asm.hic.p_ctg.fasta"
# Path to a reduced Kraken 2 database (MinKraken) used for signal detection
KRAKEN_DB="/data/users/lfalquet/SOFTS/minikraken_8GB_20200312"

# Output
# Output directory for all Kraken window results
OUTDIR="results/${SAMPLE_ID}_kraken_windows"
#FASTA file containing the sliding windows extracted from contigs, anming encodes the window size and step size. Windows overlap by 2 kbp. 
WINFA="$OUTDIR/windows_W10000_S8000.fa"
mkdir -p "$OUTDIR"

# Setup
# Load modules
# SeqKit is used to generate the sliding windows along contigs
module load SeqKit/2.6.1
# Kraken2 performs k-mer based taxonomic classification
module load Kraken2/2.1.2-gompi-2021a

# Convert GFA to FASTA when needed 
# If FASTA file is not already existing, extract contig sequences from the GFA file. 
# In GFA froamt, lines starting with "S" are segments/ contigs. Column 2 contain the contig name and column 3, the contig sequence. 
if [[ ! -f "$FASTA" ]]; then
    echo "[INFO] Converting GFA to FASTA"
    awk '/^S/{print ">"$2"\n"$3}' "$GFA" > "$FASTA"
fi

# Windowing contigs
# Split each contig into overlapping windows:
# - Window size: 10'000bp
# - Step size: 8'000 bp
# Overlap of 2'000 bp betweeen consecutive windows. Each window becomes an independent FASTA entry. 
seqkit sliding  -W 10000 -s 8000 "$FASTA" > "$WINFA"

# Kraken classification
# Run kraken on all windows 
# --confidence: minimum threshold for taxonomic assignement, low value hosen to keep weak, but potentially informative signals
# --report: produces a human-readable summary report with counts per taxon 
kraken2 \
  --threads ${SLURM_CPUS_PER_TASK} \
  --db "$KRAKEN_DB" \
  --confidence 0.01 \
  --report "$OUTDIR/kraken_windows.report" \
  --output "$OUTDIR/kraken_windows.out" \
  "$WINFA"

# Post-processing: window to taxid table 
# Extarct key fields from Kraken output: classification status, window ID and assigned taxonomic ID
# Output format: window_id  taxid status
awk '{print $2"\t"$3"\t"$1}' \
  "$OUTDIR/kraken_windows.out" \
  > "$OUTDIR/window_taxid_status.tsv"

# Done 
echo "[INFO] Done. Outputs in: $OUTDIR"