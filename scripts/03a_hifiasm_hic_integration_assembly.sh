#!/usr/bin/env bash
#SBATCH --job-name=hifiasm_hic_assembly_Lasi
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --output=logs/hifiasm_hic_assembly_Lasi_%j.out
#SBATCH --error=logs/hifiasm_hic_assembly_Lasi_%j.err

# This script assembles a genome using PacBio HiFi reads and leverages Hi-C paired-ends reads to produce a HI-C-interated assembly using hifiasm. 

set -euo pipefail

# Sample identifier
SPECIES="Lasioglossum_pauxillum"

# Input 
# HiFi read as fastq.gz, the primary sequence data used for assembly graph construction
HIFI="reads/pacbio_hifi/${SPECIES}/ERR9081702.fastq.gz"
# Hi-C reads 1 and 2
# For technical replicates, concatenate HiC via concatenate_hic_replicate.sh
HIC_R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
HIC_R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"

# Outputs
OUTDIR="assemblies/hifiasm/${SPECIES}"
# Prefix used ny hifiasm to name output files. 
PREFIX="${OUTDIR}/asm"
mkdir -p logs "$OUTDIR"

# Setup
# Load modules
module load hifiasm/0.16.1-GCCcore-10.3.0

# Sanity checks: files exist and are not empty with -s, otherwise print error and exit
for f in "$HIFI" "$HIC_R1" "$HIC_R2"; do
  [[ -s "$f" ]] || { echo "ERROR: missing file $f"; exit 1; }
done

# Run hifiasm with Hi-C support
# --h1 and --h2: add long-range linkage information for mate 1 and mate 2 reads
hifiasm \
  -o "$PREFIX" \
  -t "$SLURM_CPUS_PER_TASK" \
  --h1 "$HIC_R1" \
  --h2 "$HIC_R2" \
  "$HIFI"
