#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=hifi_busco
#SBATCH --output=logs/hifi_busco_%j.out
#SBATCH --error=logs/hifi_busco_%j.err

# Stage C: Orthogonal evidence layer
# C1: BUSCO validation
# Inputs: assembly
# Outputs: BUSCO reports and parsed busco_anchor_contigs.tsv, busco_anchor_coverage_distribution.tsv
# Run BUSCO on assembly: BUSCO-positive contigs should fall in host backbone


set -euo pipefail

# 1) Set up
SPECIES="$1"
ASM_MODE="$2"

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
ASM="${WORKDIR}/assemblies/hifiasm_clean/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"

BUSCO_OUT="results/${SPECIES}_stages/busco_anchor"
VALIDATION_OUT="results/${SPECIES}_stages/busco_validation"
mkdir -p "${BUSCO_OUT}" "${VALIDATION_OUT}"

# Load modules
module load BUSCO/5.4.2-foss-2021a
module load R/4.2.1-foss-2021a

# 2) Run BUSCO with hymenoptera_odb10 reference genes
busco \
    --in "${ASM}" \
    --out "${SPECIES}_busco" \
    --out_path "${BUSCO_OUT}" \
    --mode genome \
    --lineage hymenoptera_odb10 \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    -f

# 3) Parse full_table.tsv for contig analysis against coverage backbone 
BUSCO_TABLE="${BUSCO_OUT}/${SPECIES}_busco/run_hymenoptera_odb10/full_table.tsv"
COVERAGE_TABLE="results/${SPECIES}_stages/host_backbone/coverage_classification.tsv"

Rscript scripts/stages/C1_BUSCO_anchors.R "${SPECIES}" "${BUSCO_TABLE}" "${COVERAGE_TABLE}" "${VALIDATION_OUT}"


