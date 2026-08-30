#!/usr/bin/env bash
#SBATCH --job-name=RM_TRF
#SBATCH --partition=pibu_el8
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=logs/rm_trf_%j.out
#SBATCH --error=logs/rm_trf_%j.err

# Run RepeatMasker and TRF on the TE library outputted by EDTA
# TRF is able to detect randem repeats and microsatellites using pattern detection without any library
# Usage:
# sbatch D2_rm_trf,sh species ASM_MODE==bp or ASM_MODE==hic

set -euo pipefail

SPECIES="$1"
ASM_MODE="$2"

# Working directory
WORKDIR=/data/projects/p2025-0083_mining_cobionts
# Assembly fasta
GENOME=$WORKDIR/assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta
# basename: filename from full path
PREFIX="$(basename "$GENOME")"
# EDTA run directory
EDTA_DIR=$WORKDIR/results/${SPECIES}_stages/edta
# EDTA TE library
TELIB="${EDTA_DIR}/${PREFIX}.mod.EDTA.TElib.fa"
# Output directories for RepeatMasker and TRF
OUTDIR=$WORKDIR/results/${SPECIES}_stages/repeat_masking
RM_DIR="${OUTDIR}/repeatmasker"
TRF_DIR="${OUTDIR}/trf"
mkdir -p "$RM_DIR" "$TRF_DIR"

# Modules
module load RepeatMasker/4.1.5-foss-2021a
module load TRF/4.09.1-GCC-10.3.0

# RepeatMasker: 
# Copy assembly to output 
cp "$GENOME" "$RM_DIR/"
cd "$RM_DIR"

# -lib for custom library to search against the EDTA TElib
# -gff to produce a GFF output file
# -xsmall preserves the masked sequences
# -no_is skips bacterial insertion sequence search
# Outputs a .out tab table with one row for each repeat hit
# Columns: score    divergence  deletion    insertion   query_name  query_start query_end   query_left  strand  repeat_name repeat_class    
RepeatMasker \
    -lib "$TELIB" \
    -pa "${SLURM_CPUS_PER_TASK}" \
    -gff \
    -xsmall \
    -no_is \
    -dir "$RM_DIR" \
    "${RM_DIR}/${PREFIX}"

# Tandem Repeat Finder 
cd "$TRF_DIR"
# positional arguments: $GENOME with input FASTA, match weight 2, 7 mismatch penalty, 7 indel penalty, 80 match probability percentage, 10 indel probability percentage, 50 min alignment score, 500 max period size of a repeat unit
# -d for .dat output with one line for each repeat
trf "$GENOME" 2 7 7 80 10 50 500 -d -h || true