#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=nuwt_extract
#SBATCH --output=logs/nuwt_extract_%j.out
#SBATCH --error=logs/nuwt_extract_%j.err

# Launch the NUWT fragment extraction from all the host genome assemblies and writes one fasta for each orthogroup
# Usage:
# sbatch nuw_extrac_seq.sh

set -euo pipefail

# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
cd "$PROJECT_ROOT"
NUWT_DIR="cobionts/nuwt_scan"
RESULTS_DIR="${NUWT_DIR}/results"
ASSEMBLY_DIR="${NUWT_DIR}/host_assemblies"
# Script directory
SCRIPT="${PROJECT_ROOT}/scripts/nuwt/supergroups/Fetchregionseq_ins.py"
# Output directory
OG_SEQ_DIR="${NUWT_DIR}/og_sequences"
FRAG_TABLE="${NUWT_DIR}/nuwt_fragments.tsv"
MANIFEST="${NUWT_DIR}/og_hit_manifest.txt"

mkdir -p "${OG_SEQ_DIR}"

module load Biopython/1.79-foss-2021a

# Write a header for the table
printf "name\tspecies\tscaffold\tog\tstrand\tstart\tend\tlength\n" > "${FRAG_TABLE}"

n_done=0
for DFAM in "${RESULTS_DIR}"/*/filter/*_nuwt_hits_dfam_filtered.tbl; do
    SPECIES=$(basename "$(dirname "$(dirname "${DFAM}")")")
    ASSEMBLY=$(ls "${ASSEMBLY_DIR}/${SPECIES}_"*.fna 2>/dev/null | grep -v '\.fai$' | head -1)

    python3 "${SCRIPT}" -i "${DFAM}" -f "${ASSEMBLY}" -d "${OG_SEQ_DIR}" -t "${FRAG_TABLE}"
    n_done=$(( n_done + 1 ))
done

# Build a manifest with one line for each OG and the sequence count
for FA in "${OG_SEQ_DIR}"/OG*.fa; do
    OG=$(basename "${FA}" .fa)
    printf "%s\t%s\n" "${OG}" "$(grep -c '^>' "${FA}")"
done | sort > "${MANIFEST}"