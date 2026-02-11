#!/usr/bin/env bash
#SBATCH --job-name=candidate_kraken_intersect
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/candidate_kraken_%j.out
#SBATCH --error=logs/candidate_kraken_%j.err

# This script generates GC and coverage-based non-host candidate contigs, applies a minimum length cutoff and intersects candidates with Kraken contig classes using R. 
# Usage: sbatch scripts/03h_run_non-host_candidates.sh <SPECIES>

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <SPECIES_NAME>"
  exit 1
fi

SPECIES="$1"

# Input
GC_COV="results/${SPECIES}_gc_cov/gc_cov.tsv"
# Output
QC_DIR="results/${SPECIES}_qc_summary"
LONG_TMP="$QC_DIR/${SPECIES}_long_gc_cov.txt"

mkdir -p "$QC_DIR"

# Define host baseline as contigs ≥ 1 Mbp
# Extract GC and coverage from long contigs (≥ 1 Mbp) and sort by GC 
awk -F'\t' 'NR>1 && $2>=1000000 {print $3, $4}' "$GC_COV" \
| sort -k1,1n \
> "$LONG_TMP"

# Compute GC thresholds on sorted file
read GC_LOW GC_HIGH < <(
awk '{
  gc[NR]=$1
}
END{
  n=NR
  p10=int(n*0.1); if(p10<1)p10=1
  p90=int(n*0.9); if(p90<1)p90=1
  print gc[p10], gc[p90]
}' "$LONG_TMP"
)

# Compute coverage thresholds, internally sorted
read COV_LOW COV_HIGH < <(
awk '{
  cov[NR]=$2
}
END{
  n=NR
  asort(cov)
  p10=int(n*0.1); if(p10<1)p10=1
  p90=int(n*0.9); if(p90<1)p90=1
  print cov[p10], cov[p90]
}' "$LONG_TMP"
)

# Done
echo "[INFO] n_long_contigs: $(wc -l < "$LONG_TMP")"
echo "[INFO] GC host baseline: $GC_LOW – $GC_HIGH"
echo "[INFO] COV host baseline: $COV_LOW – $COV_HIGH"

# Outlier selection: candidates outside GC or Cov thresholds
awk -F'\t' \
'NR==1 || ($3 < '"$GC_LOW"' || $3 > '"$GC_HIGH"' || $4 < '"$COV_LOW"' || $4 > '"$COV_HIGH"')' \
"$GC_COV" \
> "$QC_DIR/nonhost_candidates.tsv"

# Length filter 50000 bp
awk -F'\t' 'NR==1 || $2>=50000' \
"$QC_DIR/nonhost_candidates.tsv" \
> "$QC_DIR/nonhost_50kb.tsv"

# Intersect with assembly-wide Kraken, GC and cov analysis

module load R/4.2.1-foss-2021a
Rscript scripts/03e_non_host_candidates_Kraken.R "$SPECIES"

echo "[INFO] Candidate Kraken intersection finished for $SPECIES"