#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=scan_mags
#SBATCH --output=logs/scan_mags_%A_%a.out
#SBATCH --error=logs/scan_mags_%A_%a.err

# Array job: one for each MAG
# Predict proteins with Prodigal
# Scan them against the marker protein HMM library
# Keep the best scoring one for each marker protein
# Each MAG gives:
# .tbl with the full hmmscan table 
# .best: each line has one OG with OG tab chosen protein id
# .faa: the actual chosen protein files, with headers as mag_id|OG to preserve the OG alignment later
# Usage: submitting with an array range that matches manifest_wolbachia.tsv
# sbatch --array=1-${N}%20 03_scan_mags.sh to have 20 jobs queued at once and not launch the whole thing

set -euo pipefail
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
OUT="${PROJECT_ROOT}/placement_out"
# Input file: the curated wolbachia manifest
# The first two columns of manifest_wolbachia.tsv hold the mag_id and then its absolute path 
MANIFEST="${PROJECT_ROOT}/scripts/metagenomeassembly/files/manifest_wolbachia.tsv"
# Genetic code for Prodigal set to bacterial
GENETIC_CODE=11

mkdir -p "${OUT}"/{proteomes,scan,orthologs}

module purge
module load prodigal/2.6.3-GCCcore-10.3.0
module load HMMER/3.3.2-gompi-2021a
module load SeqKit/2.6.1

# Find the MAG for this job: the SLURM_ARRAY_TASK_ID becomes the row number of the manifest
# +1 to skip the header, p command prints the line
ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")
# Select mag_id
MAG_ID=$(echo "${ROW}" | cut -f1)
# Select absolute path
MAG_PATH=$(echo "${ROW}" | cut -f2)

FA="${OUT}/proteomes/${MAG_ID}.input.fa" 
# Amino acid proteins
FAA="${OUT}/proteomes/${MAG_ID}.faa"
# Prodigal: predict proteins, needs unzipped FASTA so decompress the given MAG into .input.fa
gzip -dc "${MAG_PATH}" > "${FA}"

# Run Prodigal: -i input decompressed mag fasta, -a writes the predicted proteins as aa, -p for single procedure mode and -q fo suppress the progression
prodigal -i "${FA}" -a "${FAA}" -p single -g "${GENETIC_CODE}" -q
# Delete decompressed files
rm -f "${FA}"

# Prodigal output example: 
# Acromyrmex_octospinosus_metamdbg_magscot_1.faa 
# >ctg1111_1 # 3 # 779 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.345
# FVGRVMKDEQWWKILSAIEKEKDLNRDNVIEKIKEELKAEDKNWYKKWEKAGFGVNYLFE

# Run hmmscan
# Find the marker proteins within the proteomes that we predicted
# Every preditced protein gets scored against the marker protein HMM library profile
# Many profiles for one protein: need hmmscan
# "${OUT}/all_ogs.hmm" "${faa}" is database first, query second 
# --noali to suppress the alignment blocks for each hit
# --tblout a tabular summary for each sequence, one line for each protein profile hit in which columns 1, 3 and 6 are the profile name, protein id and bitscore
hmmscan --noali --cpu "${SLURM_CPUS_PER_TASK}" --tblout "${OUT}/scan/${MAG_ID}.tbl" "${OUT}/all_ogs.hmm" "${FAA}" >/dev/null

# Keep the best-scoring protein marker: 
# tblout columns: 1 = profile name with the OG from hmmbuild -n, 3 = protein id, 6 = bitscore
# Sort by OG then descending bitscore, and keep the first line seen for one OG: that is the best-scoring protein for that marker.
# grep -v to drop the # headers
# sort -k1,1 -k6,6gr: -k1,1 groups by field 1 te OG name and 1,1 means we stop at field 1 instead of the end of line
# then -k6,6 orders within each group by field 6, the bitscore and g to handle numeric scientific notations and r reverses so highest is first
# This means that all OG000046 hits are adjacent in the end, with the best-scoring at the top of the block
# This is filtered by awk '!seen[$1]++ {print $1"\t"$3}' Line 1 for OG000046: seen ["OG000046"] is unset and treated as 0, so !0 is true and prints the ++ makes it 1. In line 2, seen["OG000046"] is 1 and !1 is false so skipped
grep -v '^#' "${OUT}/scan/${MAG_ID}.tbl" | sort -k1,1 -k6,6gr \
  | awk '!seen[$1]++ {print $1"\t"$3}' > "${OUT}/scan/${MAG_ID}.best"

# Retrieve the chosen proteins with SeqKit: .best only contains names as OG000046 tab PROT_00412 but the supermatrix needs the amino acids so we need to fetch them
# So here: og=OG0000046, prot=PROT_00412 and mag_id=Ametastegia_equiseti_bin3
# .best goes in one line at a time
# in the loop, the .best looks like: OG000046    ctg34_12 which is split by read -r og prot into two variables on tab, and iteration one sets og=OG000046 prot=ctg34_12
# so in one iteration: seqkit grep -p ctg34_12 "${faa}" finds the record in the MAG proteome and prints it
# which gets piped in the sed and rewrites line 1 only: >Aleiodes_alternator_metamdbg_dastool_2|OG0000046
# then it repeats once for each line in .best, appending each rewritten record to the same output
# this results in one FASTA holding this MAG's chosen protein for every marker, each labelled with which MAG it came from and which marker it represents 
# this is necessary because the ids are ctg1111_1 and they can repeat between MAGs so we must have unique IDs
while IFS=$'\t' read -r og prot; do
  seqkit grep -p "${prot}" "${FAA}" | sed -E "1 s/^>.*/>${MAG_ID}|${og}/"
done < "${OUT}/scan/${MAG_ID}.best" > "${OUT}/orthologs/${MAG_ID}.faa"