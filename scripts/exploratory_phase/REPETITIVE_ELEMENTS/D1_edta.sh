#!/usr/bin/env bash
#SBATCH --job-name=EDTA
#SBATCH --partition=pibu_el8
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --output=logs/edta_%j.out
#SBATCH --error=logs/edta_%j.err

# =============================================================================
# Stage D1: EDTA de novo TE annotation
# Runs the EDTA pipeline for whole-genome transposable element annotation.
# Produces a classified TE library and genome-wide GFF3 annotation.
#
# EDTA performs structural TE discovery, filtering, library building and genome annotation. 
# EDTA uses stuctural signatures to finde TE:
# - LTR retotransposons via LTR_FINDER + LTR_harvest
# - TIR DNA transposons via TIR-Learner
# - Helitrons via HelitronScanner
# - LINEs via its own module 
# A non-redundant TE library is built, and the genome is annotated. 

# Key outputs:
# - *.mod.EDTA.TElib.fa contains the species-specific TE library that is used by D2
# - *.mod.EDTA.TEanno.gff3 contains the genome-wide TE annotation used by D3

# Usage: sbatch D1_edta.sh <species> <asm_mode>


# Safety flags
set -euo pipefail

# Arguments
# Species name 
SPECIES="$1"
# Assembly mode (hic or bp)
ASM_MODE="$2"

# Paths
WORKDIR=/data/projects/p2025-0083_mining_cobionts
# Assembly FASTA
GENOME=$WORKDIR/assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta
# Apptainer conainer image for EDTA and dependencies
IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif
# Output directory under the results tree
RUN_DIR="$WORKDIR/results/${SPECIES}_stages/edta"

mkdir -p "$RUN_DIR"
# EDTA writes output to the current working directory 
cd "$RUN_DIR"

# Sanity checks
# [ -s FILE ] tests that the file exists AND is non-empty
# If either input is missing, print an error to stderr (<&2) and exit with code 1
[ -s "$GENOME" ] || { echo "ERROR: genome not found: $GENOME" >&2; exit 1; }
[ -s "$IMG" ]    || { echo "ERROR: EDTA image not found: $IMG" >&2; exit 1; }

echo "[INFO] Species: $SPECIES | Mode: $ASM_MODE"
echo "[INFO] Running EDTA..."

# Run EDTA
# apptainer exec: run a command inside a container 
# --bind: make host directories visible inside the container, comma-separated list of host paths to mount
# "$IMG": the container image file 
# EDTA.pl: the main EDTA Perl script, installed inside the container 
# EDTA arguments: 
# --genome: input assembly FASTA
# --species: "others" for generic models use
# --step all: run the full pipeline, aka structure discovery + filtering + annotation
# --sensitive 1: run RepeatModeler for an extra round of TE discovery that catches fragments missed by the structural methods
# --anno 1: after TE library, run RepeatMasker internally to produce the genome-wide GFF3 annotation
# --threads: number of CPU threads read from SLURM allocation 
apptainer exec --bind "$WORKDIR","/data/courses/assembly-annotation-course" \
    "$IMG" EDTA.pl \
    --genome "$GENOME" \
    --species others \
    --step all \
    --sensitive 1 \
    --anno 1 \
    --threads "${SLURM_CPUS_PER_TASK}"

echo "[OK] EDTA completed."