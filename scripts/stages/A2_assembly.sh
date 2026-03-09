#!/usr/bin/env bash
#SBATCH --job-name=hifiasm_assembly
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --output=logs/hifiasm_assembly_%j.out
#SBATCH --error=logs/hifiasm_assembly_%j.err

# Stage A: Assembly
# A2: hifiasm assembly, integrating Hi-C information for phasing when available
# Outputs: hifiasm GFA and FASTA assemblies 
set -euo pipefail

# 1) Setup
SPECIES="$1"
ASM_MODE="$2"

# Paths
OUTDIR="assemblies/hifiasm_clean/${SPECIES}"
PREFIX="${OUTDIR}/asm"

mkdir -p "$OUTDIR" logs
# Input
HIFI=(reads/pacbio_hifi/${SPECIES}/*.fastq.gz)

HIC_R1="reads/hic/${SPECIES}/hic_R1.fastq.gz"
HIC_R2="reads/hic/${SPECIES}/hic_R2.fastq.gz"

THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "[INFO] Species: $SPECIES"
echo "[INFO] Assembly mode: $ASM_MODE"
echo "[INFO] HiFi reads: ${HIFI}"
echo "[INFO] Output dir: $OUTDIR"

# Load modules
module load hifiasm/0.16.1-GCCcore-10.3.0

# Checks
[[ -s "$HIFI" ]] || { echo "ERROR: missing HiFi file $f"; exit 1; }

if [[ "$ASM_MODE" == "hic" ]]; then
    [[ -s "$HIC_R1" ]] || { echo "ERROR: missing Hi-C R1 $HIC_R1"; exit 1; }
    [[ -s "$HIC_R2" ]] || { echo "ERROR: missing Hi-C R2 $HIC_R2"; exit 1; }
fi

# 2) Run hifiasm with or without Hi-C raw reads 
echo "[INFO] Running hifiasm"

if [[ "$ASM_MODE" == "hic" ]]; then

    hifiasm \
        -o "$PREFIX" \
        -t "$THREADS" \
        --h1 "$HIC_R1" \
        --h2 "$HIC_R2" \
        "$HIFI"

    GFA="${OUTDIR}/asm.hic.p_ctg.gfa"
    FASTA="${OUTDIR}/asm.hic.p_ctg.fasta"

else

    hifiasm \
        -o "$PREFIX" \
        -t "$THREADS" \
        "$HIFI"

    GFA="${OUTDIR}/asm.bp.p_ctg.gfa"
    FASTA="${OUTDIR}/asm.bp.p_ctg.fasta"

fi

echo "[INFO] Assembly finished"

# 3) Convert GFA to FASTA

[[ -s "$GFA" ]] || { echo "ERROR: GFA not found: $GFA"; exit 1; }

# Grab GFA sequence segments starting with S
if [[ ! -s "$FASTA" ]]; then
    echo "[INFO] Converting GFA to FASTA"
    awk '/^S/{print ">"$2"\n"$3}' "$GFA" > "$FASTA"
fi

echo "[INFO] A2 Assembly complete for $SPECIES"