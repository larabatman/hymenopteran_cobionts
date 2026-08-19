#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=nuwt_extract
#SBATCH --output=/data/users/lland/cobionts/nuwt_scan/logs/nuwt_extract_%j.out
#SBATCH --error=/data/users/lland/cobionts/nuwt_scan/logs/nuwt_extract_%j.err

# Extracts NUWT nucleotide sequences from all host genome assemblies and writes one FASTA file per Wolbachia orthogroup.
#
# Outputs:
# /data/users/lland/cobionts/nuwt_scan/og_sequences/OG*.fa: One FASTA per OG, sequences from all species with hits to that OG
# /data/users/lland/cobionts/nuwt_scan/og_hit_manifest.txt: Tab-separated list of OG_ID and sequence count
#
# Usage:
# sbatch scripts/nuwt/supergroups/nuwt_extract_seqs.sh

set -euo pipefail

COBIONTS_ROOT="/data/users/lland/cobionts"
SCRIPT="${COBIONTS_ROOT}/scripts/nuwt/supergroups/extract_nuwt_seqs.py"

# We need samtools faidx to extract sequences by coordinate from indexed assemblies.
module load SAMtools/1.13-GCC-10.3.0
module load Python/3.9.5-GCCcore-10.3.0

python3 "${SCRIPT}"
