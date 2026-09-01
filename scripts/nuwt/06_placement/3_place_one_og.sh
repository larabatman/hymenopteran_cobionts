#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --job-name=nuwt_place
#SBATCH --output=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/place_%A_%a.out
#SBATCH --error=/data/projects/p2025-0083_mining_cobionts/nuwt_scan/logs/place_%A_%a.err

# Job array: one orthogroup for each job
# For one OG with hits, builds a reference and NUWT gene tree to assign each NUWT to its Wolbachia supergroup of origin
# hmmfetch: take the OG nucleotide profile from the HMM library 
# nhmmer: scan that profile against the dereplicated reference CDS and keep the best scoring CDS of ech genome
# esl-sfetch: fetch the best CDS in each genome, concatenate it with the hit fragment of the insect genome, and make headers 
# hmmalign: align the reference and fragment combined set to the OG profile
# IQ-TREE: infer gene tree
# reroot_rename_tree_INS.py: makes the supergroup call
# Usage: 
# sbatch --array=1-${N}%50 3_place_one_og.sh

set -uo pipefail

# Working directory
PROJ_NUWT="/data/projects/p2025-0083_mining_cobionts"
RENAME_PY="${PROJ_NUWT}/scripts/nuwt/06_placement/reroot_rename_tree_INS.py"
DB_DIR="nuwt_scan"
NUWT_DIR="cobionts/nuwt_scan"
# Database directory
HMM="${DB_DIR}/hmm_database/Wolbachia_all.hmm"
REF_CDS="${DB_DIR}/wolbachia_cds/ref_cds.ffn"
CONV="${DB_DIR}/wolbachia_database/accession_to_supergroup.tsv"
# Results directory
PLACE_DIR="${NUWT_DIR}/placement"
OG_SEQ_DIR="${NUWT_DIR}/og_sequences"
SKIP_OG="${NUWT_DIR}/manual_skip_ogs.txt"
ASSIGN_DIR="${NUWT_DIR}/assignments"
OG_LIST="${NUWT_DIR}/hit_ogs.txt"
# Env for IQTREE and ETE3
IQTREE_ENV="${PROJ_NUWT}/.conda_envs/iqtree"
ETE_PY="${PROJ_NUWT}/.conda_envs/ete3/bin/python"

mkdir -p "${ASSIGN_DIR}" "${PLACE_DIR}/logs"

# Modules
module load HMMER/3.3.2-gompi-2021a
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${IQTREE_ENV}"

# Threshold for CDS and OG search
EVALUE="1e-5"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Select the OG for this task
OG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${OG_LIST}")
# Define a working directory for each OG
WORK="${PLACE_DIR}/og/${OG}"
mkdir -p "${WORK}"

# Find the corresponding OG fasta
NUWT="${OG_SEQ_DIR}/${OG}.fa"

# hmmfetch
# Find the corresponding OG profile
hmmfetch -o "${WORK}/${OG}.hmm" "${HMM}" "${OG}"

# nhmmer
# Find the corresponding OG sequence in each dereplicated reference genome 
# The profile of the OG is used to search against the pool of reference CDS
nhmmer --tblout "${WORK}/${OG}.tbl" -E "${EVALUE}" --cpu "${THREADS}" "${WORK}/${OG}.hmm" "${REF_CDS}" > /dev/null 2>&1

# Keep the best scoring CDS for each reference genome
# id=$1 is the target name as GCA_000008025_00359. 
# acc=id to copy that id, then strip the trailing name left by prokka keeping only GCA_000008025
# Find maximum: if (s > best[acc]) { best[acc]=s; keep[acc]=id } best holds the highest scroe seen so far for that genome and keep holds the corresponding CDS ID
# END { for (a in keep) print keep[a] } after all lines are read, print the winning CDS ID for each genome for iterates the keys and keep[a] has the value
# Result: one lie for each dereplicated reference genome, with the best scoring nucleotide for that OG profile
grep -v '^#' "${WORK}/${OG}.tbl" \
  | awk '{ id=$1; acc=id; sub(/_[0-9]+$/,"",acc); s=$14+0;
           if (s > best[acc]) { best[acc]=s; keep[acc]=id } }
         END { for (a in keep) print keep[a] }' \
  > "${WORK}/ref.ids"

# est-fetch
# Take the reference CDS and put them together
esl-sfetch -o "${WORK}/ref.fa" -f "${REF_CDS}" "${WORK}/ref.ids"
# Concatenate the reference CDS with the insect fragment hits
cat "${WORK}/ref.fa" "${NUWT}" > "${WORK}/combined.raw.fa"
# Clean the headers
# The problem is that the nuwt headers contain a column >INS_Acromyrmex_lobicornis_CM107051.1_OG0000000:25114388-25114868
# IQ-TREE writes in Newick syntax where the columns mean branch length so need to get rid of it
# gsub thus converts : ( ) ; , [ ] to _
awk '/^>/ { h=$1; gsub(/[:();,\[\]]/,"_",h); print h; next } { print }' "${WORK}/combined.raw.fa" > "${WORK}/combined.fa"

# hmmalign
# Align the reference sequences and insect fragments to their OG profile
# This produces an alignmeent where all sequences have identical length
# >GCA_00008025_00359
# ATCGATCGATATAGATCTGATCGA
# >GCA_000022285_00417
# ATCGATCGATATAGATCTGATCGA
# >INS_Abia_candens_OZ125706.1_OG0001339_15060639-15061970
# A-CGA---ATATAGA----TCGA
# --trim discards the columns on top of the profile 
hmmalign --dna --trim -o "${WORK}/${OG}.sto" "${WORK}/${OG}.hmm" "${WORK}/combined.fa"
# Convert the Stockholm alignement to aligned FASTA for IQ-TREE
esl-reformat afa "${WORK}/${OG}.sto" > "${WORK}/${OG}.afa"

# IQ-TREE
# Build the gene tree for this OG using the alignment and -m MFP ModelFinder Plus
iqtree2 -s "${WORK}/${OG}.afa" -m MFP -T "${THREADS}" --prefix "${WORK}/${OG}" > "${WORK}/${OG}.iqtree.log" 2>&1

# reroot_rename_ins.py
# Assign supergroup using adapted python script
# -t for the Newick tree, -a alignment, -og OG ID -c accession to supergroup table, -m the skipped lise
# ete3 tree library needed 
"${ETE_PY}" "${RENAME_PY}" \
    -t "${WORK}/${OG}.treefile" \
    -a "${WORK}/${OG}.afa" \
    -og "${OG}" \
    -c "${CONV}" \
    -m "${SKIP_OG}" \
    > "${ASSIGN_DIR}/${OG}.assign.tsv" \
    2> "${WORK}/${OG}.assign.log"