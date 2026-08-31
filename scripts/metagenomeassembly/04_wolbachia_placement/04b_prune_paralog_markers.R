#!/usr/bin/env Rscript

# Drop protein markers when the hits are ambiguous or duplicated in many MAGs
# Write a new cleaned_tree_ogs.txt file with the kept marker proteins
# Usage:
# Rscript 04b_prune_paralog_markers.R

library(readr)
library(dplyr)
# Table with the og assignments and their evaluation 
one_og <- "placement_out/qc/per_og_qc.tsv"
# Text with the complete set single-marker proteins 
tree_ogs <- "tree_ogs.txt"
# Output file with the cleaned og list
clean_ogs <- "tree_ogs.clean.txt"
# Output file with the dropped marker proteins
dropped_tsv <- "placement_out/qc/dropped_markers.tsv"
# The fraction of MAGs for which the marker protein is ambiguous or duplicated, marker protein is dropped if above this
FLAG_FRAC <- 0.30

qc <- read_tsv(one_og, show_col_types = FALSE) %>%
  mutate(FLAG_FRAC = (dup + ambiguous) / n_mags, drop = frac_flagged >= FLAG_FRAC )

dropped <- qc %>% filter(drop) %>%
  transmute(og, n_mags, dup, ambiguous, frac_flagged = round(frac_flagged, 3), reason = sprintf("ambiguous/dup in >=%.0f%% of MAGs", 100 * FLAG_FRAC))
write_tsv(dropped, dropped_tsv)

# Read complete marker protein list, strip whitespace and drop blanks
all_ogs  <- trimws(readLines(tree_ogs))
all_ogs <- all_ogs[all_ogs != ""]
# setdiff to choose the kept ogs: everything in all_ogs that is not in dropped
keep_ogs <- setdiff(all_ogs, dropped$og)
# Write the clean_tree_ogs.txt
writeLines(keep_ogs, clean_ogs)