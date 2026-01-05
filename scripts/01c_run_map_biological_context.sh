#!/bin/bash

#SBATCH --job-name=ncbi_biosample_meta
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/ncbi_biosample_meta%j.out
#SBATCH --error=logs/ncbi_biosample_meta%j.err

# This script is a wrapper that launches 01c_map_biological_context.py
set -euo pipefail 

mkdir -p logs
mkdir -p data/annotation

# Load Python module
module load Python/3.9.5-GCCcore-10.3.0

# Input and output paths
INPUT="data/data_inventory.tsv"
OUTPUT="data/annotation/data_inventory_biosample_annotation.tsv"

# Run script
python3 -u scripts/01c_map_biological_context.py "$INPUT" > "$OUTPUT"