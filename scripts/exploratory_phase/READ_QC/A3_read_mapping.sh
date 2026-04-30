#!/usr/bin/env bash
#SBATCH --job-name=read_mapping
#SBATCH --partition=pibu_el8
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/read_mapping_%j.out
#SBATCH --error=logs/read_mapping_%j.err

# =============================================================================
# Stage A3: Read mapping
# Map filtered HiFi reads back to assembly for mapping rate and coverage
# Usage: sbatch A3_read_mapping.sh <species> <asm_mode>

set -euo pipefail

# Setup
# Arguments
SPECIES="$1"
ASM_MODE="$2"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Paths
WORKDIR="/data/projects/p2025-0083_mining_cobionts"
cd "$WORKDIR"

ASM="assemblies/hifiasm/${SPECIES}/asm.${ASM_MODE}.p_ctg.fasta"
HIFI="results/archive_exploratory_phase/${SPECIES}_stages/read_qc/hifi_qc/hifi_filtered/hifi.filtered.fastq.gz"
OUTDIR="results/${SPECIES}_stages/assembly_qc"
mkdir -p "$OUTDIR"

module purge
module load minimap2/2.20-GCCcore-10.3.0 SAMtools/1.13-GCC-10.3.0

# Map, sort, index
minimap2 -t "$THREADS" -ax map-hifi "$ASM" "$HIFI" \
    | samtools sort -@ 4 -o "$OUTDIR/reads.bam"
samtools index "$OUTDIR/reads.bam"

# Mapping stats
samtools flagstat "$OUTDIR/reads.bam" > "$OUTDIR/mapping_flagstat.txt"
samtools stats   "$OUTDIR/reads.bam" > "$OUTDIR/mapping_stats.txt"

# Coverage 
samtools coverage "$OUTDIR/reads.bam" > "$OUTDIR/coverage_per_contig.tsv"

# Length-weighted mean depth across all contigs
# samtools produces TSV with col $7=meandepth, the mean read depth for that contig
# -v sp="$SPECIES": -v passes a shell variable into awk as an internal awk variable
# BEGIN block runs once before any line is read
# OFS="\t" sets Output Filed Separator to a tab character
# Block 2: runs once per line, only when NR>1, skipping the header row in line 1
# For each contig: len = $3 -$2 + 1 computes the length of the contig in bp, as $2 is the startpos and $3 the endpos. The +1 is needed to have inclusive coordinates on both ends
# tl += len accumulates the total length of all contigs seen so far, t = total length. Running sum across all rows
# tc += len * $7 accumulates the length-weighted depth as tc = total coverage. Each contig computes its depth as $7 = meandepth, multiplied by its longth
# Longer contigs contibute proportionally more to the final mean
# END block runs once after all lines have been processed
# prints species and mean_depth as header row of the output TSV, the comma being replaced by OFS set in BEGIN
# print result row: sp = species passed in via -v, tl>0 ? tc/tl : "NA" is a guard if the total length is greater than 0, compute the weighted mean depth, otherwise prin "NA" to avoid error

awk -v sp="$SPECIES" 'BEGIN{OFS="\t"}
NR>1 { len=$3-$2+1; tl+=len; tc+=len*$7 }
END  { print "species","mean_depth"; print sp, (tl>0 ? tc/tl : "NA") }
' "$OUTDIR/coverage_per_contig.tsv" > "$OUTDIR/coverage_summary.tsv"

echo "[INFO] Read mapping QC complete for $SPECIES"