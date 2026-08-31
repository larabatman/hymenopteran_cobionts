#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=nuwt_place_prep
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/prep_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/prep_%j.err

# Dereplicates the pool to the reference accession list, select the CDS from those genomes and find the sequences corresponding to each OG profile
# Usage:
# 2c_prep_reference_cds.sh

set -euo pipefail
# Working directory
PROJ_ROOT="/data/projects/p2025-0083_mining_cobionts/nuwt_scan"
cd "$PROJ_ROOT"
# Database and references
DB_DIR="nuwt_scan"
# Al 1444 genome CDS pool
CDS_FULL="${DB_DIR}/wolbachia_cds/all_wolbachia_cds.ffn"
# Dereplicated reference list
DREP="${DB_DIR}/wolbachia_database/wolbachia_reference.txt"
# HMM profiles
HMM="${DB_DIR}/hmm_database/Wolbachia_all.hmm"
# Results
NUWT_DIR="cobionts/nuwt_scan"
# List of OG to place
MANIFEST="${NUWT_DIR}/og_hit_manifest.txt" 
# Output directory
REF_CDS="${PROJ_NUWT}/wolbachia_cds/ref_cds.ffn"
PLACE_DIR="${USER_NUWT}/placement"
HIT_OGS="${PLACE_DIR}/hit_ogs.txt"
MANUAL="${PLACE_DIR}/manual_skip_ogs.txt"
SUBSET="${PLACE_DIR}/reference_subset.txt"
# Outgroups and dropped accessions
OUTGROUPS="Ace Ama Ech Eru mAli mBlo mPat"
DROP="GCA_023661085"

mkdir -p "${PLACE_DIR}/logs"

module load HMMER/3.3.2-gompi-2021a

# Build the list of genomes whose CDS will go into the trees as references
# grep -v -x "{DROP} "${DREP}" prints every line of wolbachia_reference.txt with the dRep accession minus that one we don't want
# Sorts and remove duplicates
{
    grep -v -x "${DROP}" "${DREP}"
    for o in ${OUTGROUPS}; do echo "$o"; done
} | sort -u > "${SUBSET}"

# Select the CDS that belong to the subset 
# CDS IDs are GCA001232_00001 in Prokka .ffn: strip the trailing _00001 to get the accession and keep the line if that accession is in the subset.
# Build a set: keep holds GCA_000008025, Ace, Ama, ...
# Then read all_wolbachia_cds.ffn.fai: $1 of a .fai line is the sequence name, so id: GCA_000008025_00359
# sub() substitutes the first match of underscore followed by one or more digits anchored at the end of the string and GCA_000008025_00359 becomes GCA_000008025
# if (acc in keep) tests membership
awk 'NR==FNR { keep[$1]; next }
     { id=$1; acc=id; sub(/_[0-9]+$/,"",acc); if (acc in keep) print id }' \
    "${SUBSET}" "${CDS_FULL}.fai" > "${PLACE_DIR}/ref_cds.ids"

# Index the full CDS into a library and search within it: Easel
# esl-sfetch: sequence fetch, does the same as samtools faidx aka retrieving sequence from a fasta by name using an index, format .ssi instead of .fai
esl-sfetch --index "${CDS_FULL}"
# Extract the CDS from the reference genomes
esl-sfetch -o "${REF_CDS}" -f "${CDS_FULL}" "${PLACE_DIR}/ref_cds.ids"
# Index the subset
esl-sfetch --index "${REF_CDS}"
# Index the HMM library as well
hmmfetch --index "${HMM}"