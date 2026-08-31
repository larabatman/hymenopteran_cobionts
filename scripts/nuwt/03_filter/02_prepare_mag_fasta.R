#!/usr/bin/env Rscript

# Find a species name and produce a fasta containing all the Wolbachia MAGs for that species
# Writes unique headers with prefix mag0_, mag1_, ...
# Usage:
# Rscript 02_prepare_mag_fasta.R species wolbachia_mags.tsv project_directory out_fa binner


library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

SPECIES <- args[1]
WOL_TSV <- args[2]
PROJECT_DIR <- args[3]
OUT_FA <- args[4]
RUN_BINNER <- ifelse(is.na(args[5]), "full_binners", args[5])

# Select a species row
mags <- read_tsv(WOL_TSV, show_col_types = FALSE) %>%
  rename(sheet_species = `Host species`, assembler = Assembler, tool = `Refinement Tool`) %>%
  mutate(sheet_species = gsub(" ", "_", sheet_species)) %>%
  # Cut to rows for one species, containing as many rows as there are MAGs for that species
  filter(sheet_species == SPECIES) %>%
  select(assembler, tool)

# Build the path directory to the bin
bin_dir <- file.path(PROJECT_DIR, "results", SPECIES, RUN_BINNER, "bins", "fasta")

# Find the MAG paths for each mag_files in the species
mag_files <- mags %>%
# Path pattern for each row
  mutate(pattern = file.path(bin_dir, tool, sprintf("%s_%s_%s_*.fa.gz", species, assembler, tool))) %>%
  pull(pattern) %>%
  # Extend to paths with Sy.glob
  lapply(Sys.glob) %>%
  # One character vector: mag_files are the path to each MAGs in the species
  unlist()

# Write the combined fasta
# seq_along(): 1, 2, 3 and the function is applied once for each file, unlist to concatenate everything to a character vector
lines <- unlist(lapply(seq_along(mag_files), function(i) {
  # Prefix with i -1L so it starts at 0, >mag0_, >mag_1 as i iterates
  prefix <- sprintf(">mag%d_", i - 1L)
  # Decompress each fasta inte a character vector, rewrite the ones starting with > to the prefix
  sub("^>", prefix, readLines(gzfile(mag_files[i])))
}))
writeLines(lines, OUT_FA)