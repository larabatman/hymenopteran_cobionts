#!/usr/bin/env bash
#SBATCH --job-name=fastplong_env
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/logs/fastplong_env_%j.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/logs/fastplong_env_%j.err

set -euo pipefail

PROJECT="/data/projects/p2025-0083_mining_cobionts"
ENV_DIR="${PROJECT}/.conda_envs/fastplong"
ENV_YML="${PROJECT}/envs/fastplong_env.yml"

module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"

mkdir -p "$(dirname "$ENV_DIR")"
mkdir -p "$(dirname "$ENV_YML")"

cat > "$ENV_YML" << 'YAML'
name: fastplong
channels:
  - conda-forge
  - bioconda
channel_priority: strict
dependencies:
  - fastplong
YAML

conda env create -p "$ENV_DIR" -f "$ENV_YML"