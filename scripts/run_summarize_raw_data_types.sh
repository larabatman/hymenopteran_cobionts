#!/bin/bash

#SBATCH --job-name=summarize_ncbi_query
#SBATCH --partition=pibu_el8
#SBATCH --time=0:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=logs/summarize_ncbi_query_%j.out
#SBATCH --error=logs/summarize_ncbi_query_%j.err

# This script is a wrapper that launches summarize_raw_data_types.R
set -euo pipefail 

mkdir -p logs

# Load R module
module load R/4.2.1-foss-2021a

# Run script
Rscript scripts/summarize_raw_data_types.R