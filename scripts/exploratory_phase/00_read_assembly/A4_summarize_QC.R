#!/usr/bin/env Rscript

# R script to collect assembly statistics, launched by A4_summarize_Qc.sh
# Usage:
# Rscript A4_summarize_QC.R species workdir

library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
WORKDIR <- args[2]
# Working directories
QC_DIR <- file.path(WORKDIR, "results", paste0(SPECIES, "_stages"), "assembly_qc")
ASM_DIR <- file.path(WORKDIR, "assemblies", "hifiasm", SPECIES)

# Assembly statistics from A2 seqkit
basic_stats <- read_tsv(file.path(ASM_DIR, "assembly_basic_stats.tsv"),show_col_types = FALSE)

# Mapping rate from samtools flagstat
flagstat <- read_lines(file.path(QC_DIR, "mapping_flagstat.txt"))
mapped_line <- flagstat[str_detect(flagstat, " mapped \\(")][1]
mapped_pct <- as.numeric(str_extract(mapped_line, "[0-9.]+(?=%)"))

# Assembly QC table:
qc_row <- data.frame(
  species = SPECIES,
  assembly_size = basic_stats$sum_len,
  n_contigs = basic_stats$num_seqs,
  avg_len = basic_stats$avg_len,
  mapped_pct = mapped_pct
)

write_tsv( qc_row, file.path(QC_DIR, "assembly_qc_summary.tsv"))