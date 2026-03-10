#!/usr/bin/env bash
#SBATCH --job-name=hic_pairs
#SBATCH --partition=pibu_el8
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=logs/hic_pairs_%j.out
#SBATCH --error=logs/hic_pairs_%j.err

# Stage C: Orthogonal evidence layer
# C2: Hi-C pairing information, convert alignements to contacts using pairtools
# Inputs: hic.filtered.bam and host_backbone.tsv produced in earlier steps, that we know are host
# Produces an edge list contig_links.tsv, computing the total contacts per contig
set -euo pipefail
# 1) Setup
SPECIES="$1"
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

BAM="results/${SPECIES}_stages/hic_map/hic.filtered.bam"
ANCHORS_TSV="results/${SPECIES}_stages/host_backbone/host_backbone.tsv"

OUTDIR="results/${SPECIES}_stages/hic_validation"
mkdir -p "$OUTDIR"

PAIRS="${OUTDIR}/hic.pairs.gz"
LINKS="${OUTDIR}/contig_links.tsv"
HOST_LIST="${OUTDIR}/host_anchor.list"
CHROMSIZES="${OUTDIR}/assembly.chrom.sizes"


[[ -s "$BAM" ]] || { echo "ERROR: BAM missing"; exit 1; }
[[ -s "$ANCHORS_TSV" ]] || { echo "ERROR: anchors missing"; exit 1; }


# Load Conda
module purge
module load Anaconda3/2022.05
module load R/4.2.1-foss-2021a

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${WORKDIR}/.conda_envs/hicpairs"

unset PYTHONPATH PYTHONHOME LD_LIBRARY_PATH

echo "[INFO] Using pairtools at:"
which pairtools
pairtools --version

# Scratch tmp
TMP_BASE="/scratch/${USER}/hicpairs_${SPECIES}_${SLURM_JOB_ID}"
mkdir -p "$TMP_BASE"
export TMPDIR="$TMP_BASE"
trap 'rm -rf "$TMP_BASE"' EXIT

echo "[INFO] TMPDIR=$TMPDIR"

# 2) Buil host list
# This extracts host_contig_ids for later
awk -F'\t' 'NR==1{next} NF{print $1}' "$ANCHORS_TSV" | sort -u > "$HOST_LIST"
[[ -s "$HOST_LIST" ]] || { echo "ERROR: host list empty"; exit 1; }

# 3) Extract contig size
# Build chrom.sizes from BAM header for pairtools
# @SQ SN:contig LN:length
samtools view -H "$BAM" \
  | awk -F'\t' '$1=="@SQ"{
      sn=""; ln="";
      for(i=1;i<=NF;i++){
        if($i ~ /^SN:/) sn=substr($i,4)
        if($i ~ /^LN:/) ln=substr($i,4)
      }
      if(sn!="" && ln!="") print sn"\t"ln
    }' > "$CHROMSIZES"

[[ -s "$CHROMSIZES" ]] || { echo "ERROR: chrom sizes empty"; exit 1; }

# 4) Name sort BAM for pairtools
NAME_SORTED="${OUTDIR}/hic.name_sorted.bam"
samtools sort -n -@ 4 -o "$NAME_SORTED" "$BAM"
BAM="$NAME_SORTED"

# 5) Parse Hi-C pairs using pairtools
# Construction of the actual Hi-C contacts
PARSED="${OUTDIR}/hic.parsed.pairs"
# Filters:
# --min-mapq 30 keeps MAPQ read > 30
# -- walks-policy 5unique for complex Hi-C molecules
pairtools parse \
  --chroms-path "$CHROMSIZES" \
  --assembly "$SPECIES" \
  --nproc-in 2 \
  --nproc-out 1 \
  --min-mapq 30 \
  --walks-policy 5unique \
  --drop-sam \
  --drop-seq \
  "$BAM" \
> "$PARSED"

# Sanity checks
echo "[INFO] Parsed pairs size:"
ls -lh "$PARSED"

head -n 30 "$PARSED" > "${OUTDIR}/hic.parsed.preview.txt"

grep -q '^#' "$PARSED" || { echo "ERROR: no header in parsed pairs"; exit 1; }
BODY_LINES=$(grep -vc '^#' "$PARSED" || true)
[[ "$BODY_LINES" -ge 1 ]] || { echo "ERROR: parsed pairs empty"; exit 1; }

# 7) Sorting for downstream analysis
pairtools sort \
  --nproc 1 \
  --tmpdir "$TMPDIR" \
  "$PARSED" \
| bgzip -c > "$PAIRS"

echo "[INFO] Finished sorting"
ls -lh "$PAIRS"

# 8) Build contig contact table from pairs
# Keep inter-contig contacts, count them with awk, creating an edge list
# contig_links.tsv columns: contig1 contig2 n_links
echo "[INFO] Building contig link table: $LINKS"

pairtools select '(chrom1!=chrom2)' "$PAIRS" \
  | awk 'BEGIN{OFS="\t"}
         $1 ~ /^#/ {next}
         {
           c1=$2; c2=$4;
           if (c1=="" || c2=="") next;
           # make undirected edges canonical
           if (c1>c2) {t=c1; c1=c2; c2=t}
           key=c1"\t"c2;
           cnt[key]++
         }
         END{for (k in cnt) print k, cnt[k]}' \
  | sort -k3,3nr > "$LINKS"

[[ -s "$LINKS" ]] || { echo "ERROR: contig_links.tsv empty"; exit 1; }
echo "[INFO] Link table done:"
head -n 5 "$LINKS"
