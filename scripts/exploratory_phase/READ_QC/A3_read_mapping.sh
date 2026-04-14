#!/usr/bin/env bash
#SBATCH --job-name=read_mapping
#SBATCH --partition=pibu_el8
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/read_mapping_%j.out
#SBATCH --error=logs/read_mapping_%j.err

# A3: Map filtered HiFi reads back to assembly for mapping rate and coverage
# Usage: sbatch A3_read_mapping.sh <species> <asm_mode>
set -euo pipefail

SPECIES="$1"
ASM_MODE="$2"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

ASM="assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
HIFI="results/archive_exploratory_phase/${SPECIES}_stages/read_qc/hifi_qc/hifi_filtered/hifi.filtered.fastq.gz"
OUTDIR="results/${SPECIES}_stages/assembly_qc"
mkdir -p "$OUTDIR"

module purge
module load minimap2/2.20-GCCcore-10.3.0 SAMtools/1.13-GCC-10.3.0

# ── Map, sort, index ──
minimap2 -t "$THREADS" -ax map-hifi "$ASM" "$HIFI" \
    | samtools sort -@ 4 -o "$OUTDIR/reads.bam"
samtools index "$OUTDIR/reads.bam"

# ── Mapping stats ──
samtools flagstat "$OUTDIR/reads.bam" > "$OUTDIR/mapping_flagstat.txt"
samtools stats   "$OUTDIR/reads.bam" > "$OUTDIR/mapping_stats.txt"

# ── Coverage ──
samtools coverage "$OUTDIR/reads.bam" > "$OUTDIR/coverage_per_contig.tsv"

# Length-weighted mean depth across all contigs
awk -v sp="$SPECIES" 'BEGIN{OFS="\t"}
NR>1 { len=$3-$2+1; tl+=len; tc+=len*$7 }
END  { print "species","mean_depth"; print sp, (tl>0 ? tc/tl : "NA") }
' "$OUTDIR/coverage_per_contig.tsv" > "$OUTDIR/coverage_summary.tsv"

echo "[INFO] Read mapping QC complete for $SPECIES"