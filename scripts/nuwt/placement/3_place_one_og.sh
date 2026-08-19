#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --job-name=nuwt_place
#SBATCH --output=/data/users/lland/cobionts/nuwt_scan/placement/logs/place_%A_%a.out
#SBATCH --error=/data/users/lland/cobionts/nuwt_scan/placement/logs/place_%A_%a.err

# job array: one orthogroup per task
# For one OG with NUWT hits, build the reference+NUWT nucleotide tree and assign each NUWT to a Wolbachia supergroup:
# 1. hmmfetch the OG's nucleotide profile from the library.
# 2. nhmmer that profile against the dRep reference CDS in ref_cds.ffn, keeping the best-scoring CDS per genome. These represent the OG's reference orthologues aka CDS 
# 3. esl-sfetch those CDS, concatenate with this OG's NUWTs, make proper headers
# 4. hmmalign the combined set to the OG profile (--trim keeps match columns), convert Stockholm to an aligned FASTA with esl-reformat.
# 5. IQ-TREE with model search (-m MFP)
# 6. reroot_rename_tree_INS.py for per-NUWT supergroup call.
#
#
# File types: 
# .hmm profile
# .ffn CDS nucleotide FASTA
# .sto Stockholm alignment
# .afa aligned FASTA
# .ckp.gz IQ-TREE checkpoint
# .treefile Newick
# .assign.tsv result table with OG, INS, supergroup, nuwt_name
#
# Submission
# sbatch --array=1-$(wc -l < .../placement/hit_ogs.txt)%50 place_one_og.sh to send 50 jobs at a time

set -uo pipefail

PROJ_NUWT="/data/projects/p2025-0083_mining_cobionts/nuwt_scan"
USER_NUWT="/data/users/lland/cobionts/nuwt_scan"

HMM="${PROJ_NUWT}/hmm_database/Wolbachia_all.hmm"
REF_CDS="${PROJ_NUWT}/wolbachia_cds/ref_cds.ffn"
CONV="${PROJ_NUWT}/wolbachia_database/accession_to_supergroup.tsv"

PLACE_DIR="${USER_NUWT}/placement"
OG_SEQ_DIR="${USER_NUWT}/og_sequences"
MANUAL="${PLACE_DIR}/manual_skip_ogs.txt"
ASSIGN_DIR="${PLACE_DIR}/assignments"
RENAME_PY="/data/users/lland/cobionts/scripts/nuwt/supergroups/reroot_rename_tree_INS.py"

# OG list: hit_ogs.txt by default, overridable from the environment so a retry subset can be run with the same script.
OG_LIST="${OG_LIST:-${PLACE_DIR}/hit_ogs.txt}"

IQTREE_ENV="/data/projects/p2025-0083_mining_cobionts/.conda_envs/iqtree"
ETE_PY="/data/projects/p2025-0083_mining_cobionts/.conda_envs/ete3/bin/python"

EVALUE="1e-5" # CDS-to-OG membership cutoff.
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${ASSIGN_DIR}" "${PLACE_DIR}/logs"

# Select the OG for this task
OG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${OG_LIST}")

WORK="${PLACE_DIR}/og/${OG}"
DONE="${WORK}/.done"
mkdir -p "${WORK}"
# Find its hit fasta
NUWT="${OG_SEQ_DIR}/${OG}.fa"

# Modules
module load HMMER/3.3.2-gompi-2021a
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${IQTREE_ENV}"
IQT=$(command -v iqtree2)

# 1. Find this OG profile
# wolbachia_all.hmm holds all 2647 HMM profiles concatenated and this pulls it out by name using the .ssi index built by prep_reference_cds.sh
hmmfetch -o "${WORK}/${OG}.hmm" "${HMM}" "${OG}"

# 2. Find this OG's gene in each reference genome
# ref_cds.ffn holfs all the CDS from the dRep reference genomes, but we don't know which OG each CDS belongs to.
# Here we use the profile as a definition of the OG, so searching it against the pool of reference must recover its original CDS
# For every genome in ref_cds.ffn that hit the profile, keep its best-scoring CDS as that genome's copy og this OG's gene
# nhmmscan: many profiles, one sequence set: take this genome and scan against the howle HMM library
# nhmmer: one profile, many sequences: take this one profile and serach it through this pile of sequences
nhmmer --tblout "${WORK}/${OG}.tbl" -E "${EVALUE}" --cpu "${THREADS}" \
       "${WORK}/${OG}.hmm" "${REF_CDS}" > /dev/null 2>&1

# Clean the output: reduce the nhmmer hit table to one CDS per gene, the best scoring one 
# first: grep -v to get rid of the # comments and header lines: only data rows remain
# id=$1 is the target name as GCA_000008025_00359. The sequence in --tblout is the target name as the profile is the query, which is not the same in nhmmscan were column 1 was the OG
# acc=id to copy that id, then strip the trailing name left by prokka keeping only GCA_000008025
# then column 14 is the bitscore: +0 to treat it as a number rather than a string, so > compares numerically 
# running maximum patter: if (s > best[acc]) { best[acc]=s; keep[acc]=id } best holds the hgihest scroe seen so far, for that genome and keep hold the CDS ID that holds it
# both are updated together if a hit happens
# END { for (a in keep) print keep[a] } after all lines are read, print the winning CDS ID for each genome for iterates the keys and keep[a] has the value
grep -v '^#' "${WORK}/${OG}.tbl" \
  | awk '{ id=$1; acc=id; sub(/_[0-9]+$/,"",acc); s=$14+0;
           if (s > best[acc]) { best[acc]=s; keep[acc]=id } }
         END { for (a in keep) print keep[a] }' \
  > "${WORK}/ref.ids"
# This results in max one line per genome aka 160 for a well-conserved OG and fewer if some lacked the gene 

# 3. Fetch the reference sequences, put them together with the NUWT sequences, and clean the headers
# The list of IDs form step 2 are turned into sequences by searching in the .ssi pre-built
# The output: ref.fa, the reference tips for this OG's tree
esl-sfetch -o "${WORK}/ref.fa" -f "${REF_CDS}" "${WORK}/ref.ids"
# Concatenate the reference CDS with its insect fragment hits
# Everyhing meets here!
cat "${WORK}/ref.fa" "${NUWT}" > "${WORK}/combined.raw.fa"
# h=$1 keeps only the first whitespace delimited piece, as NUWT headers could carry text
# Need to apply Newick syntax characters: the tree format uses parenthese for clades, commas between siblings, and colons before branch length and semicolon to terminate
# So a name must not carry any of these or that would corrupt the tree
# gsub thus converts : ( ) ; , [ ] to _
# NUWT headers contain: >INS_Abia_candens_OZ125706.1_OG0001339:15060639-15061970 so the colon must go
awk '/^>/ { h=$1; gsub(/[:();,\[\]]/,"_",h); print h; next } { print }' \
    "${WORK}/combined.raw.fa" > "${WORK}/combined.fa"

# 4. Align the profile
# hmmalign aligns every sequence in combined.fa, references and NUWTs to the OG's profile rather than to each other
# We chose that rather than direct MAFFT because the NUWT might be short and diverged and fragmented: they'd change the column structure
# This produces an alignmeent where all sequences have identical length, for instance:
# >GCA_00008025_00359
# ATCGATCGATATAGATCTGATCGA
# >GCA_000022285_00417
# ATCGATCGATATAGATCTGATCGA
# >INS_Abia_candens_OZ125706.1_OG0001339_15060639-15061970
# A-CGA---ATATAGA----TCGA
# --dna: nucleotide sequences
# --trim: discards the columns outside the profile's match states aka the insertions each sequence has relative to the model
hmmalign --dna --trim -o "${WORK}/${OG}.sto" "${WORK}/${OG}.hmm" "${WORK}/combined.fa"
# Convert the Stockholm alignement to aligned FASTA for IQ-TREE
esl-reformat afa "${WORK}/${OG}.sto" > "${WORK}/${OG}.afa"

# 5. Build the tree
# -s "${WORK}/${OG}.afa" contains the input alignemnt for sequence alignement reading FASTA
# -m MFP: ModelFinder Plus which means that IQ-TREE tests a range of them, scores them, and uses the winner 
# -prefix names every output file: IQ-TREE writes lots of outputs: .treefile which is the Newick of interst, .iqtree for a report of the model, -log, -skp.gz checopoint, .model.gz, .mdlist
"${IQT}" -s "${WORK}/${OG}.afa" -m MFP -T "${THREADS}" \
         --prefix "${WORK}/${OG}" > "${WORK}/${OG}.iqtree.log" 2>&1

# 6. Assign supergroup using Vancaester et al. script 
# -t: Newicj tree
# -a: alignment to classify each sequence as outgroup, reference or NUWT
# -og: the OG IS so it can be put on each row
# -c the accession to supergroup table
# -m the manually skipped list
# specifically: the scripts needs ete3 tree library to run 
"${ETE_PY}" "${RENAME_PY}" \
    -t "${WORK}/${OG}.treefile" \
    -a "${WORK}/${OG}.afa" \
    -og "${OG}" \
    -c "${CONV}" \
    -m "${MANUAL}" \
    > "${ASSIGN_DIR}/${OG}.assign.tsv" \
    2> "${WORK}/${OG}.assign.log" \
    || { echo "[$(date)] assignment failed for ${OG} (see ${WORK}/${OG}.assign.log)"; exit 1; }

N_CALL=$(wc -l < "${ASSIGN_DIR}/${OG}.assign.tsv")
echo "[$(date)] ${OG}: ${N_CALL} NUWT supergroup calls -> ${ASSIGN_DIR}/${OG}.assign.tsv"
touch "${DONE}"