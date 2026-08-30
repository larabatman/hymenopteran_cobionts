#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=repeat_analysis
#SBATCH --output=logs/repeat_analysis_%j.out
#SBATCH --error=logs/repeat_analysis_%j.err

# Repeat element analysis: parse EDTA, RepeatMasker and TRF into an interval table and merge overlapping intervals, then make blobplots
# Usage: 
# sbatch D3_run_repeat_analysis.sh species ASM_MODE==bp or ASM_MODE==hic

set -euo pipefail

SPECIES="$1"
ASM_MODE="$2"

# Working directory
WORKDIR=/data/projects/p2025-0083_mining_cobionts
cd "$WORKDIR"
SCRIPTS="scripts/exploratory_phase/03_repeats"
prefix="asm.${ASM_MODE}.p_ctg.fasta"
# EDTA GFF3 input directory
EDTA_GFF3="results/${SPECIES}_stages/edta/${PREFIX}.mod.EDTA.TEanno.gff3"
# RepeatMasker .out input directory
RM_OUT="results/${SPECIES}_stages/repeat_masking/repeatmasker/${PREFIX}.out"
# TRF .dat input directory
TRF_DAT="results/${SPECIES}_stages/repeat_masking/trf/${PREFIX}.2.7.7.80.10.50.500.dat"
# Coverage table input directory
COV="results/${SPECIES}_stages/host_backbone/coverage_classification.tsv"
# Output directory
OUTDIR="results/${SPECIES}_stages/repeat_analysis"
mkdir -p "$OUTDIR"

# Module
module load R/4.2.1-foss-2021a

Rscript "$SCRIPTS/D3a_parse_repeats.R" "$SPECIES" "$EDTA_GFF3" "$RM_OUT" "$TRF_DAT" "$OUTDIR"

# Results from that script:
INTERVALS="${OUTDIR}/repeat_intervals.tsv.gz"

Rscript "$SCRIPTS/D3b_merge_intervals.R" "$SPECIES" "$INTERVALS" "$COV" "$OUTDIR"
# Results from that script:
REPEAT_TABLE="${OUTDIR}/contig_repeat_coverage.tsv"

Rscript "$SCRIPTS/D3c_repeat_plots.R" "$SPECIES" "$REPEAT_TABLE" "$OUTDIR"