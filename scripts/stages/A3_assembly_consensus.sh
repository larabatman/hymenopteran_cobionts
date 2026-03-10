#!/usr/bin/env bash
#SBATCH --job-name=assembly_merqury
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=logs/assembly_merqury_%j.out
#SBATCH --error=logs/assembly_merqury_%j.err

# Stage A: Assembly
# A3: Merqury consensus
# Input: hifiasm primary assembly
# This script aims to further produce useful assembly metrics, such as QV and completeness

# 1) Setup 
set -euo pipefail

SPECIES="$1"
ASM_MODE="$2"

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

ASM="${WORKDIR}/assemblies/hifiasm_clean/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"

OUTDIR="${WORKDIR}/results/${SPECIES}_stages/assembly_qc/merqury"
mkdir -p "$OUTDIR"

CONTAINER="/containers/apptainer/merqury_1.3.sif"
export MERQURY="/usr/local/share/merqury"

RAW_DB="$OUTDIR/hifi_raw.meryl"
MERYL_DB="$OUTDIR/hifi.meryl"
HIFI_READ="reads/pacbio_hifi/${SPECIES}"/*.fastq.gz

echo "[INFO] Species: $SPECIES"
echo "[INFO] Assembly: $ASM"
echo "[INFO] Output dir: $OUTDIR"

# 2) Build raw k-mer database
if [[ ! -d "$RAW_DB" ]]; then
    echo "[INFO] Building raw meryl database"

    apptainer exec --bind /data "$CONTAINER" \
        meryl count k=21 \
        memory=90 \
        threads=$SLURM_CPUS_PER_TASK \
        output "$RAW_DB" \
        "$HIFI_READ"
fi

# 3) Remove singleton kmers from sequencing errors
if [[ ! -d "$MERYL_DB" ]]; then
    echo "[INFO] Filtering singleton kmers"

    apptainer exec --bind /data "$CONTAINER" \
        meryl greater-than 1 "$RAW_DB" \
        output "$MERYL_DB"
fi

# 4) Run Merqury
RUN_DIR="$OUTDIR/${SPECIES}"
mkdir -p "$RUN_DIR"

echo "[INFO] Running Merqury"

apptainer exec --bind /data "$CONTAINER" bash -lc '
export MERQURY=/usr/local/share/merqury
cd "'"$RUN_DIR"'"
"$MERQURY/merqury.sh" "'"$MERYL_DB"'" "'"$ASM"'" "'"$SPECIES"'"
'

# Done
echo "[INFO] Merqury finished for $SPECIES"

# 5) Parse Merqury output
# Paths
QC_DIR="${WORKDIR}/results/${SPECIES}_stages/assembly_qc"
MERQURY_DIR="${OUTDIR}/${SPECIES}"
QV_FILE="${MERQURY_DIR}/${SPECIES}.qv"
COMP_FILE="${MERQURY_DIR}/${SPECIES}.completeness.stats"
SUMMARY_FILE="${QC_DIR}/merqury_summary.tsv"

# Parse Merqury output 
echo "[INFO] Parsing Merqury output"

if [[ -s "$QV_FILE" && -s "$COMP_FILE" ]]; then

    # QV is column 4 in the .qv file
    QV=$(awk 'NR==1 {print $4}' "$QV_FILE")

    # completeness percentage is column 5 where column 2 == "all"
    COMP=$(awk '$2=="all"{print $5}' "$COMP_FILE")
    # write in summary file
    echo -e "species\tQV\tcompleteness_percent" > "$SUMMARY_FILE"
    echo -e "${SPECIES}\t${QV}\t${COMP}" >> "$SUMMARY_FILE"

    echo "[INFO] Merqury metrics parsed"

else
    echo "[WARN] Merqury files missing, writing NA"

    echo -e "species\tQV\tcompleteness_percent" > "$SUMMARY_FILE"
    echo -e "${SPECIES}\tNA\tNA" >> "$SUMMARY_FILE"
fi