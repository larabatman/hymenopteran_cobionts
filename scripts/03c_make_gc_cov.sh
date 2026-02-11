#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=nomada_gc_cov
#SBATCH --output=logs/nomada_gc_cov_%j.out
#SBATCH --error=logs/nomada_gc_cov_%j.err

# This script computes per-contig length and GC content from the assembly fasta, and mean read depth (coverage) by mapping reads back to the assembly. It then merges them into a single table: gc_cov.tsv

set -euo pipefail 

# Project root
WORKDIR="/data/projects/p2025-0083_a_cross-species_pipeline_for_mining_cobionts_in_hymenopteran_genomes"
cd "$WORKDIR"

# Inputs
# Sample identifier 
SPECIES="Bombus_terrestris"
# Path to primary contigs assembly FASTA, reference for extracting contig GC/ length stats and mapping reads to compute coverage
ASM="${WORKDIR}/assemblies/hifiasm/QC/${SPECIES}/${SPECIES}.p_ctg.fa"
# Path to HiFi FASTQ.gz reads 
HIFI="${WORKDIR}/reads/pacbio_hifi/${SPECIES}/ERR6558189.fastq.gz"

# Outputs 
OUTDIR="${WORKDIR}/results/${SPECIES}_gc_cov"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Logs
echo "[INFO] Assembly: $ASM"
echo "[INFO] HiFi reads: $HIFI"
echo "[INFO] Output dir: $OUTDIR"

# Setup
# Load modules
# SeqKit is used to extract contig states (name, length, GC) from FASTA
module load SeqKit/2.6.1
# minimap2 is used to align HiFi reads back to the assembly with the map-hifi preset
module load minimap2/2.20-GCCcore-10.3.0
# samtools: used to sort/ index BAM and compute coverage per contig
module load SAMtools/1.13-GCC-10.3.0


# GC and length per contig
# SeqKit fx2tab outputs tabular stats per FASTA record:
# -n: contig name
# -l: length
# -g: GC content
# The awk command sets output field separator to TAB, prints a header once and prints 3 columns per contig. 
# The output file gc_len.tsv looks like: contig len gc
seqkit fx2tab -n -l -g "$ASM" \
  | awk 'BEGIN{OFS="\t"; print "contig","len","gc"} {print $1,$2,$3}' \
  > gc_len.tsv

# Map HiFi reads to assembly and produce sorted BAM
# minimap2 alignement with -a for output in the SAM format and -x map-hifi for preset tuned for PacBion HiFi reads.
# The inputs are the aseembly and the HiFi reads, piped into samtools sort to produce a coordinate-sorted BAM. 
# The final output is reads.bam
minimap2 -t "$SLURM_CPUS_PER_TASK" -ax map-hifi "$ASM" "$HIFI" \
  | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o reads.bam
# Build BAM index reads.bam.bai for fast per-contig queries
samtools index reads.bam

# Per-contig coverage with samtools
# samtools coverage summarizes coverage per reference sequence (per contig here). It output a header line and then one line per contig
# The column 7 is the mean depth, and the output file is the raw samtools coverage output coverage.tsv
samtools coverage reads.bam > coverage.tsv

# Now, extract inly the contig name and mean depth:
# The awk logic prints the header "contig mean_cov", skips the first ine which is the coverage.tsv header and prints contigs and mean depth
# The output file cov.tsv contains: contig  mean_cov
awk 'BEGIN{OFS="\t"; print "contig","mean_cov"}
     NR==1{next}
     {print $1,$7}' coverage.tsv > cov.tsv

# Merge GC and length with coverage into gc_cov.tsv
# Join gc_len and cov.tsv by contig name which is column 1 in both. 
# The first awk reads gc_len.tsv (NR==FNR is true only for the first file). It then stors, for each contig name, a[contig] = len<TAB>gc , then reads cov.tsv by skipping its header and printin conitg, (len, gc), mean_cov
# The second awk adds the final header row. 
# The final output file gc_civ.tsv has columns: contig  len gc  mean_cov
awk 'BEGIN{OFS="\t"}
     NR==FNR{a[$1]=$2 OFS $3; next}
     FNR==1{next}
     ($1 in a){print $1, a[$1], $2}' \
  gc_len.tsv cov.tsv \
| awk 'BEGIN{OFS="\t"; print "contig","len","gc","mean_cov"} {print}' \
> gc_cov.tsv

# Done
echo "[INFO] Wrote: $OUTDIR/gc_cov.tsv"
