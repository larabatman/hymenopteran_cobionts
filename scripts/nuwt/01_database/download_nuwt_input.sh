#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=download_nuwt_inputs
#SBATCH --output=logs/download_nuwt_inputs_%j.out
#SBATCH --error=logs/download_nuwt_inputs_%j.err

# Download the Wolbachia deposite from Vancaester et al. Zenodo archives
# Builds the concatenated HMM library for nhmmscan
# Usage
# sbatch download_nuwt_input.sh

set -euo pipefail
# Working directory 
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
# Output directory
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
HMM_DIR="${NUWT_DIR}/hmm_database"
ZENODO_ZIP="${HMM_DIR}/NUWTs_Release_v3.0.zip"
mkdir -p "${HMM_DIR}" "${NUWT_DIR}/logs"

module load HMMER/3.3.2-gompi-2021a

# Zenodo archive record 15032066
wget -O "${ZENODO_ZIP}" "https://zenodo.org/records/15032066/files/NUWTs_Release_v3.0.zip?download=1"

# Extract the directories of interest: 
# 1_WolbachiaAnnotation contains the Prokka .faa files for each Wolbachia genome
# 2_WolbachiaOrthoFinder contains the Orthogroups.tsv.gz files
# 4_WolbachiaOrthoAlignments contains the protein alignments for each OG
# 5_WolbachiaHMM contains the nucleotide HMM files for each OG
unzip -o "${ZENODO_ZIP}" "NUWTs_Release_v3.0/Data/1_WolbachiaAnnotation/*" "NUWTs_Release_v3.0/Data/2_WolbachiaOrthoFinder/*" "NUWTs_Release_v3.0/Data/4_WolbachiaOrthoAlignments/*" "NUWTs_Release_v3.0/Data/5_WolbachiaHMM/*" -d "${HMM_DIR}/"
rm "${ZENODO_ZIP}"

# Concatenate the HMM files into a library
# Output directory: 
HMM_DB="${HMM_DIR}/Wolbachia_all.hmm"
# Find all the .hmm file paths inside the deposit, sort them and cat the files into Wolbachia_all.hmm, one after the other
find "${HMM_DIR}/NUWTs_Release_v3.0/Data/5_WolbachiaHMM/" \
    -name "*.hmm" | sort | xargs cat > "${HMM_DB}"
# Compress the profiles into a searchable library .h files
hmmpress "${HMM_DB}"