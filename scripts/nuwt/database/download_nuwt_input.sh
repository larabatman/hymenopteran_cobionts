#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=download_nuwt_inputs
#SBATCH --output=logs/download_nuwt_inputs_%j.out
#SBATCH --error=logs/download_nuwt_inputs_%j.err

# Downloads the Wolbachia HMM library from the Vancaester et al. 2025 Zenodo archive and builds the concatenated HMM database used by nhmmscan.
#
# Also extracts the protein alignment files (directories 1, 2, 4) needed to build the OG annotation table
#
# Outputs:
# nuwt_scan/hmm_database/Wolbachia_all.hmm with the concatenated HMMER
# nuwt_scan/hmm_database/Wolbachia_all.hmm.h3{f,i,m,p} with all the hmmpress index files
# nuwt_scan/hmm_database/NUWTs_Release_v3.0/Data/1_WolbachiaAnnotation/
# nuwt_scan/hmm_database/NUWTs_Release_v3.0/Data/2_WolbachiaOrthoFinder/
# nuwt_scan/hmm_database/NUWTs_Release_v3.0/Data/4_WolbachiaOrthoAlignments/
# nuwt_scan/hmm_database/NUWTs_Release_v3.0/Data/5_WolbachiaHMM/

# Usage:
#  sbatch download_nuwt_inputs.sh

set -euo pipefail

PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
HMM_DIR="${NUWT_DIR}/hmm_database"
ZENODO_ZIP="${HMM_DIR}/NUWTs_Release_v3.0.zip"

mkdir -p "${HMM_DIR}" "${NUWT_DIR}/logs"

module load HMMER/3.3.2-gompi-2021a

# Download Zenodo archive
# The Zenodo record 15032066 (Vancaester et al. 2025) contains the full NUWT release including HMMs, protein alignments, Prokka annotations, and OrthoFinder results. We download the full zip and extract selectively.
# ?download=1 forces Zenodo to serve raw file bytes rather than a preview page

wget -O "${ZENODO_ZIP}" "https://zenodo.org/records/15032066/files/NUWTs_Release_v3.0.zip?download=1"

# Extract  directories
# The zip internal structure is: NUWTs_Release_v3.0/Data/<directory>/
# We extract four directories:
# 1_WolbachiaAnnotation: Prokka FAA files per Wolbachia genome
# 2_WolbachiaOrthoFinder: OrthoFinder Orthogroups.tsv.gz
# 4_WolbachiaOrthoAlignments: Protein alignments per OG for blastp
# 5_WolbachiaHMM: Nucleotide HMM files, one per OG

unzip -o "${ZENODO_ZIP}" \
    "NUWTs_Release_v3.0/Data/1_WolbachiaAnnotation/*" \
    "NUWTs_Release_v3.0/Data/2_WolbachiaOrthoFinder/*" \
    "NUWTs_Release_v3.0/Data/4_WolbachiaOrthoAlignments/*" \
    "NUWTs_Release_v3.0/Data/5_WolbachiaHMM/*" \
    -d "${HMM_DIR}/"

rm "${ZENODO_ZIP}"


# Concatenate HMM files into a single database 
# nhmmscan requires a single concatenated .hmm file as its database, but aach OG has its own .hmm file
# They need to be concatenated, in sorted order for reproducibility

HMM_DB="${HMM_DIR}/Wolbachia_all.hmm"

# find dir -name "*.hmm" walks the directory and prints the path of every file whose name matches
# sort to put them in alphabetical order
find "${HMM_DIR}/NUWTs_Release_v3.0/Data/5_WolbachiaHMM/" \
    -name "*.hmm" | \
    sort | \
    xargs cat > "${HMM_DB}"

# Compress the HMM database
# hmmpress builds four binary index files: .h3f .h3i .h3m .h3p that nhmmscan requires for efficient searching

hmmpress "${HMM_DB}"