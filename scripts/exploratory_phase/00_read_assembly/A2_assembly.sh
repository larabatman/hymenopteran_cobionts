#!/usr/bin/env bash
#SBATCH --job-name=hifiasm_assembly
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=164G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/hifiasm_assembly_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/hifiasm_assembly_%j.err

# This script launches hifiasm assembly on hifi clean reads, with two options: either HiC data present and helps assembly, either HiC data not available 
# ASM_MODE=bp for no HiC, and ASM_MODE=hic when HiC is available 
# Also converts the primary assembly output GFA to FASTA for downstream analysis
# Usage:
# sbatch A2_hifiasm.sh species ASM_MODE

set -euo pipefail

# Arguments
SPECIES="$1"
ASM_MODE="$2"

# Working directory
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"
# Out directory
OUTDIR="assemblies/hifiasm/${SPECIES}"
# Prefic for hifiasm
PREFIX="${OUTDIR}/asm"

mkdir -p "$OUTDIR" logs
# HiFi clean reads
HIFI="results/${SPECIES}_stages/read_qc/hifi_qc/hifi_filtered/*.fastq.gz"
# HiC clean reads
HIC_R1="reads/hic_clean/${SPECIES}/hic_R1.clean.fastq.gz"
HIC_R2="reads/hic_clean/${SPECIES}/hic_R2.clean.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
# Modules
module purge
module load hifiasm/0.16.1-GCCcore-10.3.0
module load SeqKit/2.6.1

# If HiFi only assembly ASM_MODE==bp
# define the GFA and FASTA files
if [[ "$ASM_MODE" == "bp" ]]; then
    hifiasm \
        -o "${PREFIX}" \
        -t "$THREADS" \
        "${HIFI[@]}"

    GFA="${PREFIX}.bp.p_ctg.gfa"
    FASTA="${PREFIX}.bp.p_ctg.fasta"
fi

# If HiFi and HiC reads: ASM_MODE==hic
if [[ "$ASM_MODE" == "hic" ]]; then
    hifiasm \
        -o "${PREFIX}" \
        -t "$THREADS" \
        --h1 "$HIC_R1" \
        --h2 "$HIC_R2" \
        "${HIFI[@]}"

    GFA="${PREFIX}.hic.p_ctg.gfa"
    FASTA="${PREFIX}.hic.p_ctg.fasta"
fi

# Converting GFA to FASTA: awk to grab lines that contain the sequences in the GFA file
# The lines that start with S contain the fasta sequence: 
# /^S/ describes the pattern 
# Field 2 of an S line is the segment name: prefix with > to make a FASTA header
# And field 3 is the sequence
# Add "\n" for new line between them
awk '/^S/{print ">"$2"\n"$3}' "$GFA" > "$FASTA"

# Assembly statistics with seqkit
seqkit stats -T "$FASTA" > "${OUTDIR}/assembly_basic_stats.tsv"
