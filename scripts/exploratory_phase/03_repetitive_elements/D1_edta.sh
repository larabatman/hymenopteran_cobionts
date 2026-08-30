#!/usr/bin/env bash
#SBATCH --job-name=EDTA
#SBATCH --partition=pibu_el8
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --output=logs/edta_%j.out
#SBATCH --error=logs/edta_%j.err

# Run EDTA transposable elements annotation
# TE library and GFF3 annotation:
# .TElib.fa for the TE library
# .TEanno.gff3: TE annotation
# Usage:
# sbatch D1_edta.sh species ASM_MODE==bp or ASM_MODE==hic

set -euo pipefail

SPECIES="$1"
ASM_MODE="$2"
# Working directory
WORKDIR=/data/projects/p2025-0083_mining_cobionts
# Assembly fasta
GENOME=$WORKDIR/assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta
# Apptainer image
IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif
# Output directory
RUN_DIR="$WORKDIR/results/${SPECIES}_stages/edta"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

# EDTA with apptainer exec --bind for the image and EDTA.pk the main EDTA Perl script
# --genome input assembly fasta
# --species others 
# --step all to run the full pipeline
# --sensitive 1 for RepeatModeler
# --anno 1 produces the GFF3 annotation after the TE library is built
apptainer exec --bind "$WORKDIR","/data/courses/assembly-annotation-course" \
    "$IMG" EDTA.pl \
    --genome "$GENOME" \
    --species others \
    --step all \
    --sensitive 1 \
    --anno 1 \
    --threads "${SLURM_CPUS_PER_TASK}"