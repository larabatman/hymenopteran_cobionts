#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=download_og_alignments
#SBATCH --output=logs/download_og_alignments_%j.out
#SBATCH --error=logs/download_og_alignments_%j.err

# This script downloads the Wolbachia orthogroup protein alignement files from the Vancaester et al. Zenodo archive and extracts one representative protein sequence per orthogroup into og_representatives.fna
# That file will be the quey input for the SiwssProt blastp step
# First script to run before build_og_lookup_table.sh and build_og_gene_annotations.py.
#
# Minimum sequence length threshold is 10 aa 
#
#  On the sequence type:
# Despite the .fna extension on the output file, the sequences in 4_WolbachiaOrthoAlignments are protein (amino acid), not nucleotide.
# og_representatives.fna therefore contains protein sequences, but it should really be .faa
#
# Outputs:
# nuwt_scan/hmm_database/og_representatives.fna One protein sequence per orthogroup, header is >OG_ID
# nuwt_scan/hmm_database/4_WolbachiaOrthoAlignments/ Raw alignment files
#
# Usage:
# sbatch scripts/nuwt/database/download_og_alignments.sh

set -euo pipefail

# Setup
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
HMM_DIR="${NUWT_DIR}/hmm_database"
ALN_DIR="${HMM_DIR}/4_WolbachiaOrthoAlignments"
ZENODO_ZIP="${HMM_DIR}/NUWTs_Release_v3.0_alignments.zip"
OUT_FASTA="${HMM_DIR}/og_representatives.fna"

SCRIPT_DIR="/data/users/lland/cobionts/scripts/nuwt/database"
EXTRACT_PY="${SCRIPT_DIR}/extract_og_representatives.py"
 

mkdir -p "${ALN_DIR}" "${NUWT_DIR}/logs"
mkdir -p "${ALN_DIR}"


# Download the Zenodo archive 
# We use Zenodo record 16833941
#
# wget options:
# -O: write to a specific output path rather than using the remote filename
wget --no-verbose \
    -O "${ZENODO_ZIP}" \
    "https://zenodo.org/records/16833941/files/NUWTs_Release_v3.0.zip?download=1"


# Extract folder 4 which contain the protein alignments from zip
# unzip options:
# -o: overwrite existing files without prompting
# Note: the zip was created on macOS and contains __MACOSX/ resource fork
# entries alongside the real files (e.g. ._OG0000001.aln.fa.gz). These are
# harmless metadata files but would clutter the alignment directory. The
# find command below excludes them with "! -name '._*'".
unzip -o "${ZENODO_ZIP}" \
    "NUWTs_Release_v3.0/Data/4_WolbachiaOrthoAlignments/*" \
    -d "${HMM_DIR}/"

# The zip extracts to NUWTs_Release_v3.0/Data/4_WolbachiaOrthoAlignments/ inside HMM_DIR. Move the contents to ALN_DIR for consistency with the
mv "${HMM_DIR}/NUWTs_Release_v3.0/Data/4_WolbachiaOrthoAlignments/"* "${ALN_DIR}/"

# Remove zip to save disk space.
rm "${ZENODO_ZIP}"

# Extact one representative sequence per orthogroup
# The Python helper globs ALN_DIR itself It

python3 "${EXTRACT_PY}" "${ALN_DIR}" "${OUT_FASTA}"