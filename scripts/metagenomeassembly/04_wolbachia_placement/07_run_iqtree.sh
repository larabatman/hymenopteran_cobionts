#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=28-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=run_iqtree
#SBATCH --output=logs/run_iqtree_%j.out
#SBATCH --error=logs/run_iqtree_%j.err

# Infer maximum-likelihood tree with IQTREE2
# GTR20+G4 and 1000 ultrafast bootstraps
# Usage:
# sbatch 07_run_iqtree.sh

set -euo pipefail
# Working directory
OUT="placement_out"
# Supermatrix file
SUPERMATRIX="${OUT}/super.all.fasta"
# Output prefix
PREFIX="${OUT}/wolb_mag_tree"
# Outgroups
OUTGROUPS="Ace,Ama,Ech,Eru"

module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /data/projects/p2025-0083_mining_cobionts/.conda_envs/iqtree

# Run IQTREE: -s for input alignment supermatrix, -m model, -B bootstraps, -o outgrops
iqtree -s "${SUPERMATRIX}" -m GTR20+G4 -B 1000 -T "${SLURM_CPUS_PER_TASK}" -o "${OUTGROUPS}" --prefix "${PREFIX}"