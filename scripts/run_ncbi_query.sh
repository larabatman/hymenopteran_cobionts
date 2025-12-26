#!/bin/bash

#SBATCH --job-name=ncbi_meta
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/ncbi_meta_%j.out
#SBATCH --error=logs/ncbi_meta_%j.err

# This script is a wrapper that launches query_ncbi_assembly.py
set -euo pipefail 

# Load Python module
module load Python/3.9.5-GCCcore-10.3.0

# Run script
python3 -u scripts/ncbi_GCA_query.py > data/raw_data_inventory.tsv