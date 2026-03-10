#!/usr/bin/env bash
#SBATCH --job-name=hic_score
#SBATCH --partition=pibu_el8
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=logs/hic_score_%j.out
#SBATCH --error=logs/hic_score_%j.err

# Stage C: Orthogonal evidence layer
# C2: Contact tendency evaluation
# Inputs: contig_links.tsv host_anchor.list assembly.chrom.sizes
# Outputs: hic_host_link_enrichment.tsv

# 1) Setup
SPECIES="$1"
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

OUTDIR="${WORKDIR}/results/${SPECIES}_stages/hic_validation"

LINKS="${OUTDIR}/contig_links.tsv"
HOST_LIST="${OUTDIR}/host_anchor.list"
CHROMSIZES="${OUTDIR}/assembly.chrom.sizes"

# Load modules
module load R/4.2.1-foss-2021a

# 2) Launch R script 
echo "[INFO] Running R script"
Rscript scripts/stages/C2_Hi-C_scoring.R "$SPECIES" "$LINKS" "$HOST_LIST" "$OUTDIR" "$CHROMSIZES"