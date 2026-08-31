#!/usr/bin/env Rscript

# Build the OG list for the placement
# Some OGs were removed as they were thought to be hitting transposable elements in insect genomes
# Usage:
# Rscript 1_make_placement_list.R

library(readr)
library(dplyr)

NUWT_USER <-  "/data/projects/p2025-0083_mining_cobionts/cobionts/nuwt_scan"
NUWT_PROJ <- "/data/projects/p2025-0083_mining_cobionts/nuwt_scan"

LOOKUP <- file.path(NUWT_PROJ, "hmm_database", "og_lookup_table.tsv")
MANIFEST <- file.path(NUWT_USER, "og_hit_manifest.txt")
PLACE <- file.path(NUWT_USER, "placement")

dir.create(PLACE, showWarnings = FALSE, recursive = TRUE)

# Manually discarded OGs
MANUAL_TE <- c("OG0001384", "OG0001260", "OG0002151")
# Read annotation table
lookup   <- read_tsv(LOOKUP, show_col_types = FALSE)
# Read manifest, two columns: og_id and n_hits
manifest <- read_tsv(MANIFEST, col_names = c("og_id", "n_hits"), col_types = "ci")
# Make a discarded list
euk_te <- lookup %>% 
  filter(og_class == "eukaryotic_TE") %>% 
  pull(og_id)
exclude_all <- union(euk_te, MANUAL_TE)
# Make placement list: everything that was not discarded
hit_ogs <- manifest %>% 
  filter(!og_id %in% exclude_all) %>% 
  pull(og_id) %>% 
  sort()
# Write output
writeLines(hit_ogs, file.path(PLACE, "hit_ogs.txt"))