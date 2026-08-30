#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=build_hmm_library
#SBATCH --output=logs/build_hmm_library_%j.out
#SBATCH --error=logs/build_hmm_library_%j.err

# Prepares files for phylogenomic placement of the Wolbachia into their supergroups
# prot2strain.tsv: contains the map of the marker proteins to their strains. When doing the supermatrix concatenation, we need to know which strains each alignments belongs to
# all_ogs.hmm: build one HMM profile for each selected protein marker, and compress it into a searchable library which will be used by hmmscan

# Usage: 
# sbatch 01_build_hmm_library.sh

set -euo pipefail

module purge
module load HMMER/3.3.2-gompi-2021a
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
HMMDB="${PROJECT_ROOT}/nuwt_scan/hmm_database"

# OrthoFinder clustered all proteins across 1451 genomes into orthogroups, then aligned each group
# the Orthogroup.tsv tells which proteins are in an orthogroup, and the .fna give the actual proteins aligned
# so basically one is a list of IDs and the other contains the alignments for each OG and that's what we need for the HMM and matrix concatenation
# The deposit contains the alignments for each OG, with one file for each OG, each in an aligned amino acid fasta file
# OG000046.aln.fa.gz contains every protein across the set of genomes that belonged to OG0000046, already aligned to each other:
# >Ace_00025
# MKALIV--DBIADNKRCKE
# >GCA_00008025_00412
# MKALIVTEDBIADNKRCKE
# >GCA_00002285_00412
# MKALIVTEDBIADNK---E
# The headers are protein ids, not the strain names: prot2strain.tsv will help preserve the mapping from protein ids to strains

# Directory with the protein alignements for each orthogroup, from the deposite
REF_ALN_DIR="${HMMDB}/4_WolbachiaOrthoAlignments"
# Directory with the protein annotations. These are gzipped strain.gaa.gz where the basename is the strain name, matching the column names in Orthogroups.tsv meaning GCA_ accessions for Wolbachia, and Ace, Ama, Ech, Eru for outgroups
REF_ANNOT_DIR="${HMMDB}/NUWTs_Release_v3.0/Data/1_WolbachiaAnnotation"
# Marker proteins: one line for one OG
TREE_OGS="${PROJECT_ROOT}/scripts/metagenomeassembly/files/tree_ogs.txt"
# Output directory: one for each hmm profile and one for the ref_aln decompressed
OUT="${PROJECT_ROOT}/placement_out"
mkdir -p "${OUT}"/{hmm,ref_aln}

# Build a protein to strain map using the headers 
# Alignment rows are protein ids Ace_00025, GCA_000008025_00001. Every id in a strain .faa.gz belongs to that strain, so the map is just header to basename.
# Loop through the .faa.gz files in the directory, then strip its suffix which leaves the strain name GCA_00008025
# then zcat decompresses and the grep keeps the FASTA header lines and sed rewrites each into two tab separated fields 
# Then we have one line for each protein
# Ace_01009 Ace
# GCA_00008025_00412    GCA_00008025
for faa in "${REF_ANNOT_DIR}"/*.faa.gz; do
  strain=$(basename "${faa}" .faa.gz)
  zcat "${faa}" | grep '^>' | sed -E "s/^>([^[:space:]]+).*/\1\t${strain}/"
done > "${OUT}/prot2strain.tsv"

# Build one HMM profile for each OG
# We need to turn each of the marker alignments into a statistical model of that protein family and stack them into one file
# An HMM profile is build column by column from the alignment: for each position, it records which aa occur and how often, and how often that position is gapped
# We need those because then, for each MAG, when we predict its proteins, we need to map back which of those proteins are actually the ortholog marker OG000046
# zcat "${REF_ALN_DIR}/${og}.aln.fa.gz" > "${ref_aln}" hmmbuild cannot read gzip, we need decompression and we keep them for the alignment anyways
# hmmbuild --amino forces protein alphabet 
# -n "{og}": names the profile that we need to track hits
while read -r og; do
  ref_aln="${OUT}/ref_aln/${og}.faa"
  zcat "${REF_ALN_DIR}/${og}.aln.fa.gz" > "${ref_aln}"
  hmmbuild --amino -n "${og}" /dev/stdout "${ref_aln}"
done < "${TREE_OGS}" > "${OUT}/all_ogs.hmm"

# Index the library for hmmscan, producig the .h3 files
hmmpress "${OUT}/all_ogs.hmm"