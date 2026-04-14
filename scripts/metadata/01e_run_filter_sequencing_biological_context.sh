#!/bin/bash

#SBATCH --job-name=biosample_metadata_filtering
#SBATCH --partition=pibu_el8
#SBATCH --time=0:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=logs/biosample_metadata_filtering_%j.out
#SBATCH --error=logs/biosample_metadata_filtering_%j.err

# This script is a wrapper that launches 01e_filter_sequencing_biological_context.R
set -euo pipefail 

mkdir -p logs

# Load R module
module load R/4.2.1-foss-2021a

# Run script
Rscript scripts/01e_filter_sequencing_biological_context.R