#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

# Variables
SPECIES <- args[1]
WORKDIR <- args[2]

QC_DIR <- file.path(WORKDIR, "results", paste0(SPECIES, "_stages"), "assembly_qc")

# 1) Retrieve stage A4 assembly statistics from seqkit stats output
basic_stats <- read_tsv(
  file.path(QC_DIR, "assembly_basic_stats.tsv"),
  show_col_types = FALSE
)

assembly_size <- basic_stats$sum_len
n_contigs <- basic_stats$num_seqs
avg_len <- basic_stats$avg_len

# 2) Retrieve stage A4 assembly statistics from samtools flagstat
flagstat <- read_lines(file.path(QC_DIR, "mapping_flagstat.txt"))
mapped_line <- flagstat[str_detect(flagstat, " mapped \\(")][1]

mapped_pct <- str_extract(mapped_line, "[0-9.]+(?=%)") %>% as.numeric()

# 3) Retrive stage A4 per-contig coverage statistics from samtools coverage
depth_tbl <- read_tsv(
  file.path(QC_DIR, "coverage_summary.tsv"),
  show_col_types = FALSE
)

mean_depth <- depth_tbl$mean_depth

# 4) Retrieve Merqury parsed output from stage A3 assembly consensus

merqury <- read_tsv(
  file.path(QC_DIR, "merqury_summary.tsv"),
  show_col_types = FALSE
)

QV <- merqury$QV
completeness <- merqury$completeness_percent

# 5) QC thresholds
# Flag low mapping, low QV/ completeness
status <- case_when(
  mapped_pct < 95 ~ "FAIL_low_mapping",
  QV < 30 ~ "WARN_low_QV",
  completeness < 90 ~ "WARN_low_completeness",
  TRUE ~ "PASS"
)

# 6) Build final table
qc_row <- tibble(
  species = SPECIES,
  assembly_size = assembly_size,
  n_contigs = n_contigs,
  avg_len = avg_len,
  mapped_pct = mapped_pct,
  mean_depth = mean_depth,
  QV = QV,
  completeness = completeness,
  status = status
)

# Write species specific table
write_tsv(
  qc_row,
  file.path(QC_DIR, "assembly_qc_status.tsv")
)

# Write global table by binding rows on already existing global table
GLOBAL_FILE <- file.path(WORKDIR, "results", "global_qc_summary.tsv")

if (file.exists(GLOBAL_FILE)) {
  global <- read_tsv(GLOBAL_FILE, show_col_types = FALSE)
  global <- global %>%
    filter(species != SPECIES) %>%
    bind_rows(qc_row)
} else {
  global <- qc_row
}

write_tsv(global, GLOBAL_FILE)

# Done
message("QC summary updated")
message("Status: ", status)