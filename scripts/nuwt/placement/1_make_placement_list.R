#!/usr/bin/env Rscript
#
# Builds the OG lists that feed the placement step, by removing orthogroups judged to be transposable elements
# Exclusion came from "eukaryotic_TE" as well as manual inspection
# Here we filter the OGs themselves rather than the hits because for placing each OG into its reference tree against its own reference set is easier
# bacterial_IS OGs are kept excluded
#
# Files:
# in: og_lookup_table.tsv with the annotated OGs
# in: og_hit_manifest.txt the table with the raw hit counts: og_id n_hits with no header
# out: plain text, one OG id per line
#
# Usage:
# Rscript filter_te_ogs.R

library(readr)
library(dplyr)

NUWT_USER <- "/data/users/lland/cobionts/nuwt_scan"
NUWT_PROJ <- "/data/projects/p2025-0083_mining_cobionts/nuwt_scan"

LOOKUP   <- file.path(NUWT_PROJ, "hmm_database", "og_lookup_table.tsv")
MANIFEST <- file.path(NUWT_USER, "og_hit_manifest.txt")
PLACE    <- file.path(NUWT_USER, "placement")

dir.create(PLACE, showWarnings = FALSE, recursive = TRUE)

# Manually discarded OGs
MANUAL_TE <- c("OG0001384", "OG0001260", "OG0002151")

lookup   <- read_tsv(LOOKUP, show_col_types = FALSE)
# The manifest has no headers and just two columns: og_id and n_hits
manifest <- read_tsv(MANIFEST, col_names = c("og_id", "n_hits"), col_types = "ci")

# Discard list
euk_te <- lookup %>% 
  filter(og_class == "eukaryotic_TE") %>% 
  pull(og_id)
exclude_all <- union(euk_te, MANUAL_TE)

# Placement task list: all that were not discarded
hit_ogs <- manifest %>% 
  filter(!og_id %in% exclude_all) %>% 
  pull(og_id) %>% 
  sort()
writeLines(hit_ogs, file.path(PLACE, "hit_ogs.txt"))