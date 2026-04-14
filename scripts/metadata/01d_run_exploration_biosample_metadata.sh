#!/bin/bash

#SBATCH --job-name=explore_biosample_metadata
#SBATCH --partition=pibu_el8
#SBATCH --time=0:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=logs/explore_biosample_metadata_%j.out
#SBATCH --error=logs/explore_biosample_metadata_%j.err

# This script is a wrapper that launches 01d_explore_biosample_metadata.R
set -euo pipefail 

mkdir -p logs

# Load R module
module load R/4.2.1-foss-2021a

# Run script
Rscript scripts/01d_explore_biosample_metadata.R