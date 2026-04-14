#!/bin/bash

#SBATCH --job-name=query_ncbi_sequencing
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/query_ncbi_sequencing_%j.out
#SBATCH --error=logs/query_ncbi_sequencing_%j.err

# This script is a wrapper that launches 01a_map_gca_to_runs.py
# -e stops on error, -u fail on unset variable and pipefail propagate failures
set -euo pipefail 

mkdir -p logs
mkdir -p data

# Load Python module
module load Python/3.9.5-GCCcore-10.3.0

# Run script
python3 -u scripts/01a_map_gca_to_runs.py < species/Hymenopteran_genomes.csv > data/data_inventory.tsv