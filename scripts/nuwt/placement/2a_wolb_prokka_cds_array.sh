#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=wol_prokka
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/wol_prokka_%A_%a.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/wol_prokka_%A_%a.err

# Job array: run Prokka on one Wolbachia reference genome per array task.
#
# The aim is to regenerate the .ffn CDS nucleotide files for all Wolbachia reference genomes. The Zenodo archive only deposited .faa (protein) files, but the tress are built from nucleotide sequence and a nucleotide reference is needed.
# These .ffn become the CDS pool that the plcamenet sep searches: for each OG, place_one_og.sh runs that OG's nucleotide HMM profile against the pool of CDS with hmmer and keeps the bes-scoring CDS per genome. 
# Those hits are the OG's reference orthologues for that tree
# 
# Array indexing: each task reads the Nth line of the genome list file, where N = SLURM_ARRAY_TASK_ID.

set -euo pipefail

PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR="${PROJECT_ROOT}/nuwt_scan"
PROKKA_DIR="${NUWT_DIR}/wolbachia_prokka"
GENOME_LIST="${NUWT_DIR}/wolbachia_genome_list.txt"

# Pick the genome of the task at hand 
# sed -n "${N}p" prints only the Nth line of the file: this is how we map the array task ID to a specific genome
GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${GENOME_LIST}")

# Derive the accession from the filename: GCA_000008025.fna to GCA_000008025
# basename strips the directory, %.fna strips the extension
ACCESSION=$(basename "${GENOME}" .fna)

mkdir -p "{PROKKA_DIR}"

OUT_DIR="${PROKKA_DIR}/${ACCESSION}"

# Skip if already done, for resubmission of failed jobs
if [[ -s "${OUT_DIR}/${ACCESSION}.ffn" ]]; then exit 0; fi

module load prokka/1.14.5-gompi-2021a

# Prokka
# --locustag ${ACCESSION}: names genes <ACCESSION>_NNNNN. That trailing number is stripped downstream to recover which genome a CDS came from, so it's important
# --prefix ${ACCESSION}: names all output files <ACCESSION>.*
# --outdir: one directory per genome as Prokka requires a dedicated output dir
# --cpus 4: matches the SBATCH allocation
# --kingdom Bacteria: Wolbachia are bacteria
# --quiet: suppress per-gene progress chatter
prokka \
    --quiet \
    --cpus 4 \
    --kingdom Bacteria \
    --locustag "${ACCESSION}" \
    --prefix "${ACCESSION}" \
    --outdir "${OUT_DIR}" \
    --force \
    "${GENOME}"