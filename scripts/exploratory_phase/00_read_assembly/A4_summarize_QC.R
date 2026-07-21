#!/usr/bin/env Rscript

# This R script merges assembly statistics per species and globally

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

# Arguments
SPECIES <- args[1]
WORKDIR <- args[2]

QC_DIR <- file.path(WORKDIR, "results", paste0(SPECIES, "_stages"), "assembly_qc")

# 1) Basic assembly stats (seqkit stats)
basic_stats <- read_tsv(
  file.path(QC_DIR, "assembly_basic_stats.tsv"),
  show_col_types = FALSE
)

assembly_size <- basic_stats$sum_len
n_contigs <- basic_stats$num_seqs
avg_len <- basic_stats$avg_len

# 2) Mapping rate (samtools flagstat)
flagstat <- read_lines(file.path(QC_DIR, "mapping_flagstat.txt"))
mapped_line <- flagstat[str_detect(flagstat, " mapped \\(")][1]
mapped_pct <- str_extract(mapped_line, "[0-9.]+(?=%)") %>% as.numeric()

# 3) Per-contig coverage (samtools coverage)
depth_tbl  <- read_tsv(
  file.path(QC_DIR, "coverage_summary.tsv"),
  show_col_types = FALSE
)
mean_depth <- depth_tbl$mean_depth

# 4) QC thresholds based on mapping rate
# Mapping rate is important since coverage is a strong signal for the pipeline
status <- case_when(
  mapped_pct < 95 ~ "FAIL_low_mapping",
  TRUE ~ "PASS"
)

# 5) Build final table
qc_row <- tibble(
  species       = SPECIES,
  assembly_size = assembly_size,
  n_contigs     = n_contigs,
  avg_len       = avg_len,
  mapped_pct    = mapped_pct,
  mean_depth    = mean_depth,
  status        = status
)

# 6) Write species-specific table
write_tsv(
  qc_row,
  file.path(QC_DIR, "assembly_qc_status.tsv")
)

# 7) Update global summary
GLOBAL_FILE <- file.path(WORKDIR, "results", "global_qc_summary.tsv")

if (file.exists(GLOBAL_FILE)) {
  global <- read_tsv(GLOBAL_FILE, show_col_types = FALSE) %>%
    filter(species != SPECIES) %>%
    bind_rows(qc_row)
} else {
  global <- qc_row
}

write_tsv(global, GLOBAL_FILE)

# Done logs
message("QC summary updated")
message("Status: ", status)