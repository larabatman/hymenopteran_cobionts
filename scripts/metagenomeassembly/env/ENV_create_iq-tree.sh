#!/usr/bin/env bash
#SBATCH --job-name=iqtree_env
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/iqtree_env_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/iqtree_env_%j.err
set -euo pipefail

PROJECT="/data/projects/p2025-0083_mining_cobionts"
ENV_DIR="${PROJECT}/.conda_envs/iqtree"
ENV_YML="${PROJECT}/envs/iqtree_env.yml"

module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"

mkdir -p "$(dirname "$ENV_DIR")"
mkdir -p "$(dirname "$ENV_YML")"
mkdir -p "${PROJECT}/logs"

cat > "$ENV_YML" << 'YAML'
name: iqtree
channels:
  - conda-forge
  - bioconda
channel_priority: strict
dependencies:
  - iqtree=2.*
YAML

conda env create -p "$ENV_DIR" -f "$ENV_YML"