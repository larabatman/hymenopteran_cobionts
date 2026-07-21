#!/usr/bin/env bash
#SBATCH --job-name=RM_TRF
#SBATCH --partition=pibu_el8
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=logs/rm_trf_%j.out
#SBATCH --error=logs/rm_trf_%j.err

# =============================================================================
# Stage D2: RepeatMasker + TRF
# RepeatMasker uses the EDTA-derived TE library for classified repeat masking.
# EDTA's internal RM run only annotates TEs, and running RM separately allows 
# for detection of simple repeats, low-complexity and satellites DNA which are 
# non-TE repetitive elements that EDTA ignores!
#
# TRF detects tandem repeats and microsatellites. It adds another layer by 
# finding tandem repeats using pure pattern detection, without any library. 
# Together they capture simple repeats and low-complexity regions not covered
# by the TE library alone. 
#
# Requires: D1_edta.sh completed (it needs .TElib.fa)
# Usage: sbatch D2_rm_trf.sh <species> <asm_mode>

set -euo pipefail

# Same arguments as in D1
SPECIES="$1"
ASM_MODE="$2"

# Paths 
WORKDIR=/data/projects/p2025-0083_mining_cobionts
GENOME=$WORKDIR/assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta
# basename extracts the filename from the full path as EDTA names its outputs based on prefix
PREFIX="$(basename "$GENOME")"

EDTA_DIR=$WORKDIR/results/${SPECIES}_stages/edta
# The TE library built by EDTA in D1 is what makes RM species-specific rather than using a generic database
TELIB="${EDTA_DIR}/${PREFIX}.mod.EDTA.TElib.fa"

# Output directories for separate subdirs for RM and TRF
OUTDIR=$WORKDIR/results/${SPECIES}_stages/repeat_masking
RM_DIR="${OUTDIR}/repeatmasker"
TRF_DIR="${OUTDIR}/trf"
mkdir -p "$RM_DIR" "$TRF_DIR"

# Checks
[ -s "$GENOME" ] || { echo "ERROR: assembly not found: $GENOME" >&2; exit 1; }
[ -s "$TELIB" ]  || { echo "ERROR: EDTA TE library not found: $TELIB" >&2; exit 1; }

# Load modules
module load RepeatMasker/4.1.5-foss-2021a
module load TRF/4.09.1-GCC-10.3.0

# Run RepeatMasker
echo "[INFO] Running RepeatMasker on $PREFIX"
# RM modifies files in-place, so the assembly is copied to the output directory to avoid touching the original
cp "$GENOME" "$RM_DIR/"
cd "$RM_DIR"
# RM arguments: 
# -lib: custom repeat library to search against, aka our EDTA TElib. Without this, RM would use Dfam/ RepBase which are generic
# -pa: number of parallele processes
# -gff: produce GFF output file 
# -xsmall: soft-mask repeats instead of hard-mask, which preserves the sequence for downstream tools
# -no_is: skip bacterial insertion sequence search
# -dir: where to write output files 
# Key outputs: ${PREFIX}.out which is a tab-limited table with one row per repeat hit:
# columns: score, divergence, deletion, insertion, query_name, query_start, query_end, query_left, strand, repeat_name, repeat_class, ....
# RM also detects non-TE repeaty by default (Simple_repeat, Low_complexity, Satellite). These appear in the class column on field 11 of the .out file. 
RepeatMasker \
    -lib "$TELIB" \
    -pa "${SLURM_CPUS_PER_TASK}" \
    -gff \
    -xsmall \
    -no_is \
    -dir "$RM_DIR" \
    "${RM_DIR}/${PREFIX}"

echo "[OK] RepeatMasker done."

# Run TRF
echo "[INFO] Running TRF on $PREFIX"
cd "$TRF_DIR"

# TRF arguments are positional, no flags
# Using the recommended default patameters from the TRF documentation 
# $GNENOME: input FASTA
# 2: match weight
# 7: msmatch penalty
# 7: indel penalty
# 80: math probability %
# 10: indel probability %
# 50: minimum alignment score to report
# 500: maximum period size, as in repeat unit length up to 500 bp
# Flags: 
# -d: produces a .dat output file that is machine-readable, one line per repeat
# -h: suppresses HTML ouput
# TRF is single-threaded and processes one contig at a time, sequentially. For each contig, it outputs start, end, period, copies, consensus, ...

# The output file is named automatically: ${input_filename}.2.7.7.80.10.50.500.dat which are the parameters embeddded in the filename
# || true: do not fail if the TRF exits witht non-zero, which happens on edge cases when output is fine!
trf "$GENOME" 2 7 7 80 10 50 500 -d -h || true

echo "[OK] TRF done."
echo "[INFO] RepeatMasker output: $RM_DIR"
echo "[INFO] TRF output:          $TRF_DIR"