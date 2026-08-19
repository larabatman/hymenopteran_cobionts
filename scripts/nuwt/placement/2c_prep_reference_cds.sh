#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=nuwt_place_prep
#SBATCH --output=/data/users/lland/cobionts/nuwt_scan/placement/logs/prep_%j.out
#SBATCH --error=/data/users/lland/cobionts/nuwt_scan/placement/logs/prep_%j.err


# Builds the inputs the placement array reuses:
# ref_cds.ffn: the reference CDS we put in every tree: the dRep (wolbachia_reference.txt) plus 7 outgroups, minus GCA_023661085
# hit_ogs.txt: the OGs that actually have NUWTs and we want to place, from the manifest)
# manual_skip_ogs.txt — OGs the assignment script should not emit calls for
#      (empty by default; matches the authors' -m argument).
# the .ssi indexes hmmfetch and esl-sfetch need.
# Build the reference CDS subset and build the indexes the array will need
# .ffn are CDS nucleotide FASTA
# .fai are samtools index
# .ssi are Easel/HMMER sequence index
# .hmm(.ssi) is profile index.
#
# Usage  
# sbatch prep_reference_cds.sh

set -euo pipefail

PROJ_NUWT="/data/projects/p2025-0083_mining_cobionts/nuwt_scan"
USER_NUWT="/data/users/lland/cobionts/nuwt_scan"

CDS_FULL="${PROJ_NUWT}/wolbachia_cds/all_wolbachia_cds.ffn" # all 1444 genomes' CDS
REF_CDS="${PROJ_NUWT}/wolbachia_cds/ref_cds.ffn" # dRep + outgroups subset built here
DREP="${PROJ_NUWT}/wolbachia_database/wolbachia_reference.txt" # 156 dRep accessions
HMM="${PROJ_NUWT}/hmm_database/Wolbachia_all.hmm" # nucleotide HMM library
MANIFEST="${USER_NUWT}/og_hit_manifest.txt" # og_id nb_hits

PLACE_DIR="${USER_NUWT}/placement"
HIT_OGS="${PLACE_DIR}/hit_ogs.txt"
MANUAL="${PLACE_DIR}/manual_skip_ogs.txt"
SUBSET="${PLACE_DIR}/reference_subset.txt"

OUTGROUPS="Ace Ama Ech Eru mAli mBlo mPat" # non-Wolbachia outgroups for rooting
DROP="GCA_023661085"

mkdir -p "${PLACE_DIR}/logs"

module load HMMER/3.3.2-gompi-2021a

# Build reference accession: the list of genomes that go in every tree
# {...; ...; } groups commands such that their combined output can be piped as one streamand sort -u sees everything
# grep -v -x "{DROP} "${DREP}" prints every line of wolbachia_reference.txt with the dRep accession minus that one we don't want
# -v inverts the match to see lines that don't match: one is dropped 
# for loop: OUTGROUPS is one string split into 7 items 
# Sorts and remove duplicates
{
    grep -v -x "${DROP}" "${DREP}"
    for o in ${OUTGROUPS}; do echo "$o"; done
} | sort -u > "${SUBSET}"

# CDS selection: the ones that belong to these genomes and extract them 
# CDS IDs are <ACCESSION>_NNNNN in Prokka .ffn: strip the trailing _number to get the accession and keep the line if that accession is in the subset.
# NR is the total record number across all input and FNR is the record number within the current file
# NR==FNR is true while reading file 1 only
# Build a set: keep holds GCA_000008025, Ace, Ama, ...
# Then read all_wolbachia_cds.ffn.fai: $1 of a .fai line is the sequence name, so id: GCA_000008025_00359
# acc=id copies it as sub() would modify it
# sub() substitutes the first match of underscore followed by one or more digits andchored at the end of the string
# GCA_000008025_00359 becomes GCA_000008025
# if (acc in keep) tests the membership: that's why block 1 only needed the create the keys, hat we then compare to acc
# print id gives the full CDS ID, not the stripped accession 
awk 'NR==FNR { keep[$1]; next }
     { id=$1; acc=id; sub(/_[0-9]+$/,"",acc); if (acc in keep) print id }' \
    "${SUBSET}" "${CDS_FULL}.fai" > "${PLACE_DIR}/ref_cds.ids"

# Use Easel sequence-handling library 
# esl-sfetch: sequence fetch, does the same as samtools faidy aka retrieving sequence from a FASTA by name, using an index so it can look up rather than scan
# index format .ssi instead of .fai
# esl-sfetch needs an .ssi index on the source FASTA to pull sequences by ID.
# Builds all_wolbachia_cds.ffn.ssi a binary map from every sequence name to its byte position in the file needed for the next line
esl-sfetch --index "${CDS_FULL}"
# This is the extraction: 
# -f: the last argument is a file containing names, one per line rather than a single name 
# -0 write to a file
# This reads all the CDS gene names from ref_cds.ids, looks into their ssi and writes it out
# the result is ref_cds.ffn
esl-sfetch -o "${REF_CDS}" -f "${CDS_FULL}" "${PLACE_DIR}/ref_cds.ids"

# index the subset for fast per-OG fetching in the array
# same operation as the first index but on the file we just created, for downstream steps
esl-sfetch --index "${REF_CDS}"

# Index the HMM library as well
hmmfetch --index "${HMM}"