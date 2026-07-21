#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=700G
#SBATCH --job-name=metaflye
#SBATCH --output=logs/metaflye_%x_%j.out
#SBATCH --error=logs/metaflye_%x_%j.err

# run_metaflye.sh
# Assemble CLR metagenome reads with metaFlye. This produces an assembly to feed into the sanger-tol metagenomeassembly pipeline (run_metagenomeassembly_metaflye.sh).
#
# This script uses --pacbio-raw, which sets Flye's error model for noisy (~85%) CLR subreads. Only point this at species confirmed to be CLR.
#
# Read argument(optional 2nd arg): by default this assembles clr_filtered.fasta.gz. For very large read sets that OOM metaFlye (~240-270 Gbp CLR), first run subsample_clr.sh and pass the subsampled file as the 2nd argument:
# sbatch run_metaflye.sh Leptopilina_drosophilae clr_filtered.subsampled.fasta.gz
# The argument is a filename inside reads/processed/<SAMPLE_ID>/.
#
# Usage: sbatch run_metaflye.sh <SAMPLE_ID> [READS_FILENAME]
# sbatch run_metaflye.sh Muscidifurax_raptorellus
# sbatch run_metaflye.sh Leptopilina_drosophilae clr_filtered.subsampled.fasta.gz


set -euo pipefail


# Arguments

SAMPLE_ID="${1:?ERROR: SAMPLE_ID required as first argument. Usage: sbatch run_metaflye.sh <SAMPLE_ID> [READS_FILENAME]}"
# Optional 2nd arg: reads filename inside reads/processed/<SAMPLE_ID>/.
# Defaults to the full filtered CLR set.
READS_FILENAME="${2:-clr_filtered.fasta.gz}"

# Paths

PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"

# Input: CLR reads (full set by default, or a subsampled file if passed as $2).
CLR_FASTA="${PROJECT_ROOT}/reads/processed/${SAMPLE_ID}/${READS_FILENAME}"

# Output directory for metaFlye. metaFlye writes many intermediate files here; we keep the directory intact after assembly so flye.log is available for troubleshooting and --resume works if the job is interrupted.
ASM_DIR="${PROJECT_ROOT}/assemblies/${SAMPLE_ID}/metaflye"

# metaFlye writes its final assembly to assembly.fasta in the output directory.
ASM_OUT="${ASM_DIR}/assembly.fasta"
# Renamed and gzipped output matching the naming convention used by myloasm and metaMDBG, ready to pass to the pipeline YAML.
ASM_GZ="${ASM_DIR}/${SAMPLE_ID}_metaflye.fasta.gz"

# Modules
module purge
module load Flye/2.9-GCC-10.3.0

# Sanity checks

[[ -f "${CLR_FASTA}" ]] || { echo "[ERROR] Missing reads file: ${CLR_FASTA}. Run A1b_clr_qc.sh (and subsample_clr.sh if subsampling) first."; exit 1; }

echo "[INFO] Sample: ${SAMPLE_ID}"
echo "[INFO] Reads: ${CLR_FASTA}"
echo "[INFO] Output dir: ${ASM_DIR}"
echo "[INFO] Flye: $(flye --version)"

mkdir -p logs
mkdir -p "${ASM_DIR}"


# TMPDIR on fast scratch
export TMPDIR="/scratch/${USER}/metaflye_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"

# Run metaFlye
echo "[INFO] Starting metaFlye assembly at $(date)"

# --pacbio-raw: input is PacBio CLR subreads (<20% error rate). This sets Flye's internal error model appropriately for noisy reads.
# --meta: enable metagenome mode. Without this, Flye assumes uniform coverage (single genome) and will discard low-coverage organisms.
# In metagenome mode, Flye uses a k-mer selection strategy that handles the wide range of coverages across organisms.
# --out-dir: all Flye output goes here, including assembly.fasta, assembly_graph.gfa, assembly_info.txt, and flye.log.
# --threads: use all SLURM-allocated CPUs.
# --min-overlap:  minimum overlap length between reads to consider them connected in the assembly graph. 2000bp is the recommended value for CLR metagenomes from the metaFlye paper (Kolmogorov et al. 2020).
# The default auto-set value can be too high for short CLR reads; 2000 balances sensitivity and specificity for this data type.

# --resume: only pass this flag when a previous Flye run exists in the output directory; Flye aborts with "Can't find save file" if --resume is given on a fresh run. Flye 2.9 writes params.json plus flye.log when a run begins, so we detect either.
if [[ -f "${ASM_DIR}/params.json" || -f "${ASM_DIR}/flye.log" ]]; then
    echo "[INFO] Found existing Flye run - resuming from last checkpoint."
    RESUME_FLAG="--resume"
else
    echo "[INFO] Starting fresh Flye run."
    RESUME_FLAG=""
fi

flye \
    --pacbio-raw "${CLR_FASTA}" \
    --meta \
    --out-dir "${ASM_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --min-overlap 2000 \
    ${RESUME_FLAG}


# Compress and rename output
[[ -f "${ASM_OUT}" ]] || { echo "[ERROR] metaFlye did not produce assembly.fasta - check ${ASM_DIR}/flye.log"; exit 1; }

echo "[INFO] Compressing assembly..."
# -k: keep the original uncompressed assembly.fasta so --resume still works if we need to rerun, and so flye.log references remain valid.
# -f: force-overwrite a stale assembly.fasta.gz from a previous run; without -f, gzip exits non-zero on an existing target and set -e kills the job.
gzip -kf "${ASM_OUT}"
# Rename to include sample ID, matching the convention of other assemblers.
# -f so a re-run cleanly overwrites a previous renamed copy.
mv -f "${ASM_DIR}/assembly.fasta.gz" "${ASM_GZ}"

# Summary
echo "[INFO] Assembly complete at $(date)"
echo "[INFO] Output: ${ASM_GZ}"

# Quick assembly statistics from assembly_info.txt which Flye writes automatically.
# This file has per-contig coverage, length, and circularity - more informative
# than just counting FASTA headers.
if [[ -f "${ASM_DIR}/assembly_info.txt" ]]; then
    N_CONTIGS=$(tail -n +2 "${ASM_DIR}/assembly_info.txt" | wc -l)
    TOTAL_BP=$(tail -n +2 "${ASM_DIR}/assembly_info.txt" | \
        awk '{sum += $2} END {print sum}')
    echo "[INFO] Contigs: ${N_CONTIGS}"
    echo "[INFO] Total bp: ${TOTAL_BP}"
fi

echo "[INFO] Done. Next step: sbatch run_metagenomeassembly_metaflye.sh ${SAMPLE_ID}"