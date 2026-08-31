#!/usr/bin/env Rscript

# Join the lookup table of protein annotations to the hits in the filtered dfam table
# Usage
# Rscript 3_annotate_nuwt_hits.R

library(dplyr)
library(readr)
library(tidyr)

PROJECT_ROOT <- "/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR <- file.path(PROJECT_ROOT, "nuwt_scan")
RESULTS_DIR <- file.path(NUWT_DIR, "results")
SWISSPROT_TSV <- file.path(NUWT_DIR, "hmm_database", "og_lookup_table.tsv")
SUMMARY_OUT <- file.path(RESULTS_DIR, "annotation_summary.tsv")

# Read annotation table 
swissprot <- read_tsv(SWISSPROT_TSV, show_col_types = FALSE) %>%
  select(og_id, og_class, sp_description)

# Define dfam columns
dfam_col_names <- c( "og_id", "acc", "scaffold", "score", "evalue", "bias", "hmm_from", "hmm_to", "strand", "ali_from", "ali_to", "env_from", "env_to", "modlen", "description")

# Function to parse one dfam table 
parse_dfam <- function(filepath) {
  lines <- readLines(filepath)
  # Select the lines that are not headers
  data_lines <- lines[!startsWith(lines, "#")]
  if (length(data_lines) == 0) return(NULL)
  df <- read.table(
    text = data_lines,
    header = FALSE,
    quote = "",
    comment = "",
    col.names = dfam_col_names,
    colClasses = c(og_id = "character", scaffold = "character", score = "numeric", evalue = "numeric", bias = "numeric", strand = "character", ali_from = "numeric", ali_to = "numeric", env_from = "numeric", env_to = "numeric"))
  if (nrow(df) == 0) return(NULL)
  # Normalize env coordinates
  df %>%
    mutate(hit_start = pmin(env_from, env_to), hit_end = pmax(env_from, env_to), hit_len = hit_end - hit_start + 1)
}

# For each species, find the contamination filtered dfam tables 
dfam_files <- Sys.glob(file.path(RESULTS_DIR, "*", "filter", "*_nuwt_hits_dfam_filtered.tbl"))

# Main loop
# Empty list and counters
summary_rows <- list()
# For one file
for (dfam_file in dfam_files) {
  # Extract species name from path
  species <- basename(dirname(dirname(dfam_file)))
  filter_dir <- dirname(dfam_file)
  out_file <- file.path(filter_dir, paste0(species, "_nuwt_hits_annotated.tsv"))
  # Parse its dfam table
  hits <- parse_dfam(dfam_file)
  if(is.null(hits)){
    next
  }
  # Join annotation column based on og id, NA for OGs that have no swissprot table
  annotated <- hits %>%
    left_join(swissprot, by = "og_id") %>%
    mutate(
      og_class = replace_na(og_class, "no_swissprot_hit"),
      sp_description = replace_na(sp_description, "no_hit"),
      species = species) %>%
    select(species, og_id, og_class, sp_description, scaffold, strand, hit_start, hit_end, hit_len, score, evalue, bias, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, modlen)
  write_tsv(annotated, out_file)
  summary_rows[[species]] <- count(annotated, species, og_class)
  }
bind_rows(summary_rows) %>% 
  write_tsv(SUMMARY_OUT)