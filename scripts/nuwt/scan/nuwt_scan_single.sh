#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --job-name=nuwt_scan  # overridden per-species at submission with --job-name
#SBATCH --output=logs/nuwt_scan_%j.out
#SBATCH --error=logs/nuwt_scan_%j.err

# This script verifies that the assembly and pressed HMM database exist, run nhmmscan with the 2'647 Wolbachia nucleotide profiles against one host genome and writes .tbl, _dfam.tbl and .out
# Usage
# sbatch --job-name=scan_<Species> nuwt_scan_single.sh <Species> <Accession> <assembly.fna>
#
# Arguments:
# $1: Species name
# $2: Accession
# $3: Assembly path, full path to the decompressed .fna file
# The _dfam.tbl output 15 columns:
# target name: profile that matched, OG ID
# acc: accession registered inside the HMM file
# query name: the host sequence, a scaffold from its assembly
# bits: bit score, og-ods of the seuqence under the profile vs under random background model. Higher the better
# e-value: expected number of hits scoring at least this well by chance in a database this big
# bias: a correction already substracted from the score accounting for biased composition that would otherwise inflate it
# hmm-st: start position on the profile
# hmm-en: end position on the profile
# strand: + or -, the strand of the scaffold that the profile matched
# ali-st: alignment start on the scaffold
# ali-en: alignment end on the scaffold
# env-st: envelope start on the scaffold, a wider estimate of where homology plausibly extends
# env-en: envelope end
# modlen: total length of the HMM profile in model positions
# description of target: description of the HMM file 

set -euo pipefail


SPECIES="${1}"
ACCESSION="${2}"
ASSEMBLY="${3}"


PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
HMM_DB="${NUWT_DIR}/hmm_database/Wolbachia_all.hmm"
OUT_DIR="${NUWT_DIR}/results/${SPECIES}"
# hmmpress builds the binary index files (.h3f .h3i .h3m .h3p) that nhmmscan
# requires to search efficiently. Only needs to happen once per HMM database.
[[ -f "${HMM_DB}.h3i" ]] || { echo "ERROR: HMM index missing — run hmmpress on ${HMM_DB}"; exit 1; }

mkdir -p "${OUT_DIR}"

# Module
module load HMMER/3.3.2-gompi-2021a

# nhmmscan
# Searches nucleotide sequences in the genome assembly against a database of nucleotide HMM profiles (one profile per Wolbachia orthogroup)
# --tblout: compact per-hit table (one line per hit)
# --dfamtblout: Dfam-format table, same hits but with genomic coordinates
# -E 1e-30: only report hits with E-value at or below this threshold, matches Vancaester et al.
# --cpu 16: parallelise across all allocated cores; nhmmscan scales well here
# HMM database first, then the query nucleotide sequences
nhmmscan \
    --tblout "${OUT_DIR}/${SPECIES}_nuwt_hits.tbl" \
    --dfamtblout "${DFAM_OUT}" \
    -E 1e-30 \
    --cpu 16 \
    "${HMM_DB}" \
    "${ASSEMBLY}" \
    > "${OUT_DIR}/${SPECIES}_nuwt_hits.out"
