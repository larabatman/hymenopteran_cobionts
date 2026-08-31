#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --job-name=align_markers
#SBATCH --output=logs/align_markers_%j.out
#SBATCH --error=logs/align_markers_%j.err

# Insert each MAG marker protein into its marker reference alignment, such that the final alignment contains all the marker proteins from the alignment and our MAGs
# ref_aln/OG000046.faa
# >Ace_01009
#-MWSVDRRFFAWJDWKWDCFI
# >GCA_00008025_00412
#-MWSVDRRFFAWJDWKWDCFI
# ...
# Is augmented to aug/OG000046.faa
# # >Ace_01009
#-MWSVDRRFFAWJDWKWDCFI
# >GCA_00008025_00412
#-MWSVDRRFFAWJDWKWDCFI
# >Acromyrmex_octospinosus_metamdbg_magscot_1
#-MWSVDRRFFAWJDWKWDCFI
#... 
# Usage: 
# sbatch 05_mafft_add.sh

set -euo pipefail
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
OUT="${PROJECT_ROOT}/placement_out"
# Input file: the cleaned list of marker proteins
OGS="${PROJECT_ROOT}/scripts/metagenomeassembly/files/tree_ogs.clean.txt"

mkdir -p "${OUT}"/{add,aug}

module purge
module load MAFFT/7.490-GCC-10.3.0-with-extensions
module load SeqKit/2.6.1

# Loop: go through clean protein marker list one at a time
# Pull one line from OGS
while read -r og; do
  ref_aln="${OUT}/ref_aln/${og}.faa" # decompressed reference alignments for that og
  add="${OUT}/add/${og}.faa" # diretory with the MAG proteins 

# Collect the protein marker corresponding to that of in every MAG .faa file
# Like this: 
# orthologs/Acromyrmex_octospinosus_...faa
# >Acromyrmex_octospinosus_...|OG0000046
# >Acromyrmex_octospinosus_...|OG0000061
# and 
# orthologs/Aleiodes_alternator_...faa
# >Aleiodes_alternator_...|OG0000046 
# >Aleiodes_alternator_...|OG0000061
# But to build the alignment of OG0000046, we need every OG0000046 protein representative from each MAG
# so: for OG0000046, collect the records marker |OG0000046 and write them into add/OG0000046.faa
# seqkit reads all the per-MAG files in one call
# -n matches the full header
# -r treats the pattern as a regex and \$ anchors it so OG0000046 cannot also match OG00000461
# The sed then strips |OG, leaving the MAG id alone as the sequence name which becomes the tips name
seqkit grep -nrp "\|${og}\$" "${OUT}"/orthologs/*.faa 2>/dev/null | sed -E "s/^>(.*)\|${og}\$/>\1/" > "${add}" || true

# Run MAFFT: safeguard for empty files if no MAG carried the marker
# --add for the sequence to insert and align, --keeplength to keep the same number of columns
# --add: sequence to insert, aligning these against an existing alignment mode
if [ -s "${add}" ]; then
mafft --add "${add}" --keeplength --thread "${SLURM_CPUS_PER_TASK}" \
"${ref_aln}" > "${OUT}/aug/${og}.faa" 2>/dev/null
else
cp "${ref_aln}" "${OUT}/aug/${og}.faa"
fi
done < "${OGS}"