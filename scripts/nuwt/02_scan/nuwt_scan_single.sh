#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --job-name=nuwt_scan  # overridden per-species at submission with --job-name
#SBATCH --output=logs/nuwt_scan_%j.out
#SBATCH --error=logs/nuwt_scan_%j.err

# Run nhmmscan on the insect gneomes against the 2'647 Wolbachia HMM profiles
# Usage:
# sbatch nuwt_scan_single species assembly_path
#
# dfam.tbl has 15 columns: 
# target name is the profile that matched, with the OG ID
# acc is the accession registered inside the HMM file
# query name i s the host sequence as a scaffold or contig from its assembly
# bits is the bitscore
# e-value as the expected number of hits scoring at least that well by chance
# hmm-st start position on the profile
# hmm-en end position on the profile
# strand + or -
# ali-st is the alignment start on the contig
# ali-en is the alignment end on the contig
# env-st is the enveloppe start on the contig
# env-en is the envelope end in the contig
# modlen is the length of the HMM profile
# description of the target 

set -euo pipefail

SPECIES="${1}"
ASSEMBLY="${2}"
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
# Profile library
HMM_DB="${NUWT_DIR}/hmm_database/Wolbachia_all.hmm"
# Output directory
OUT_DIR="${NUWT_DIR}/results/${SPECIES}"
mkdir -p "${OUT_DIR}"

# Module
module load HMMER/3.3.2-gompi-2021a

# Run nhmmscan: -E 1e-30 for reporting only values below or at that threshold
# Search nucleotide sequences in genome assemlby against the library of HMM profiles
# --tblout for a compact table of each hit
# --dfamtblout for the _dfam.tbl 
nhmmscan \
    --tblout "${OUT_DIR}/${SPECIES}_nuwt_hits.tbl" \
    --dfamtblout "${OUT_DIR}/${SPECIES}_nuwt_hits_dfam.tbl" \
    -E 1e-30 \
    --cpu 16 \
    "${HMM_DB}" \
    "${ASSEMBLY}" \
    > "${OUT_DIR}/${SPECIES}_nuwt_hits.out"