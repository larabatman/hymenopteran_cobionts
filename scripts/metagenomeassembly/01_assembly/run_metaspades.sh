#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=768G
#SBATCH --job-name=metaspades
#SBATCH --output=logs/metaspades_%x_%j.out
#SBATCH --error=logs/metaspades_%x_%j.err

# run_metaspades.sh
# Assemble Illumina paired-end WGS reads with metaSPAdes for a single species. Output is fed into the sanger-tol pipeline via run_metagenomeassembly_metaspades.sh.
#
# Usage:
#   sbatch run_metaspades.sh <SAMPLE_ID>

set -euo pipefail

# Arguments
SAMPLE_ID="${1:?ERROR: SAMPLE_ID required as first argument. Usage: sbatch run_metaspades.sh <SAMPLE_ID>}"

# Paths
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"

# Input: filtered paired-end reads from A1c_illumina_qc.sh.
FILTERED_R1="${PROJECT_ROOT}/reads/processed/${SAMPLE_ID}/illumina_filtered_R1.fastq.gz"
FILTERED_R2="${PROJECT_ROOT}/reads/processed/${SAMPLE_ID}/illumina_filtered_R2.fastq.gz"

# metaSPAdes output directory. SPAdes writes many intermediate files here; we keep them for --restart-from if the job is killed mid-assembly.
ASM_DIR="${PROJECT_ROOT}/assemblies/${SAMPLE_ID}/metaspades"

# metaSPAdes writes contigs.fasta and scaffolds.fasta into the output dir.
# We use scaffolds.fasta as the final assembly: scaffolds are contigs joined by paired-end read information, giving longer sequences for binning.
# (for metagenomes, scaffolding can occasionally chimerically join contigs across organisms; acceptable here as these MAGs are used for presence/identity, not phylogeny)
ASM_GZ="${ASM_DIR}/${SAMPLE_ID}_metaspades.fasta.gz"

# Modules
module purge
# SPAdes 3.15.3: includes metaSPAdes mode via the --meta flag.
module load SPAdes/3.15.3-GCC-10.3.0

# Sanity checks
[[ -f "${FILTERED_R1}" ]] || { echo "[ERROR] Missing R1: ${FILTERED_R1}. Run A1c_illumina_qc.sh first."; exit 1; }
[[ -f "${FILTERED_R2}" ]] || { echo "[ERROR] Missing R2: ${FILTERED_R2}. Run A1c_illumina_qc.sh first."; exit 1; }

echo "[INFO] Sample:  ${SAMPLE_ID}"
echo "[INFO] R1: ${FILTERED_R1}"
echo "[INFO] R2: ${FILTERED_R2}"
echo "[INFO] Out dir: ${ASM_DIR}"

mkdir -p logs
mkdir -p "${ASM_DIR}"

# TMPDIR
export TMPDIR="/scratch/${USER}/metaspades_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"

# Run metaSPAdes
echo "[INFO] Starting metaSPAdes assembly at $(date)"

# --meta: enable metagenome assembly mode (metaSPAdes). Without this SPAdes assumes a single isolate genome with uniform coverage and will collapse or discard low-coverage organisms.
# -1 / -2: paired-end R1 and R2 input files. metaSPAdes currently supports only a single paired-end library
# -o: output directory. SPAdes creates it if absent.
# -t: thread count from SLURM allocation.
# -m: memory limit in GB. We pass the SLURM allocation minus a small buffer. SPAdes will fail if it exceeds this rather than being killed by the OOM killer mid-assembly. SLURM_MEM_PER_NODE is in MB, so we convert to GB.
MEM_GB=$(( SLURM_MEM_PER_NODE / 1024 - 8 ))

# --restart-from last: if a previous run was interrupted, SPAdes can resume from the last completed stage rather than restarting from scratch. We check if a params.txt file exists (written by SPAdes on first run) to decide whether to resume or start fresh.
if [[ -f "${ASM_DIR}/params.txt" ]]; then
    echo "[INFO] Found existing SPAdes run: resuming from last checkpoint."
    # --restart-from cannot be combined with input file arguments; SPAdes reads the original input paths (and k-mer set) from params.txt.
    spades.py \
        --restart-from last \
        -o "${ASM_DIR}" \
        -t "${SLURM_CPUS_PER_TASK}" \
        -m "${MEM_GB}"
else
    echo "[INFO] Starting fresh SPAdes run."
    # --only-assembler: skip SPAdes' BayesHammer read error-correction stage
    spades.py \
        --meta \
        -1 "${FILTERED_R1}" \
        -2 "${FILTERED_R2}" \
        -o "${ASM_DIR}" \
        -t "${SLURM_CPUS_PER_TASK}" \
        -m "${MEM_GB}"
fi

# Compress and rename output
[[ -f "${ASM_OUT}" ]] || { echo "[ERROR] metaSPAdes did not produce scaffolds.fasta — check ${ASM_DIR}/spades.log"; exit 1; }

echo "[INFO] Compressing assembly..."
# -k: keep original scaffolds.fasta so --restart-from still works if needed.
gzip -k "${ASM_OUT}"
mv "${ASM_DIR}/scaffolds.fasta.gz" "${ASM_GZ}"

# Summary
echo "[INFO] Assembly complete at $(date)"
echo "[INFO] Output: ${ASM_GZ}"

N_CONTIGS=$(zcat "${ASM_GZ}" | grep -c "^>")
TOTAL_BP=$(zcat "${ASM_GZ}" | grep -v "^>" | tr -d '\n' | wc -c)
echo "[INFO] Scaffolds: ${N_CONTIGS}"
echo "[INFO] Total bp: ${TOTAL_BP}"
echo "[INFO] Done. Next step: sbatch run_metagenomeassembly_metaspades.sh ${SAMPLE_ID}"