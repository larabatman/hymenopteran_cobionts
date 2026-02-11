#!/usr/bin/env Rscript

# This script aims at intersecting GC and coverage-based non-host candidate contigs (≥50 kb) with Kraken-derived contig-level taxonomy.

# Outputs are base-pair-weighted summary of Kraken contig classes among candidates and of putative non-host genomes as well as bar plots for visualization

library(dplyr)
library(readr)

# Argument parsing: species name is the first argument
args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]

# Paths
qc_dir <- file.path("results", paste0(SPECIES, "_qc_summary"))
# Inputs
# Contig-level Kraken + GC + coverage table (with all contigs)
contig_file <- file.path(qc_dir, "generic_contig_tax_gc_cov.tsv")
# Non-host candidates filtered to ≥50 kb (awk-produced, no header)
candidate_file <- file.path(qc_dir, "nonhost_50kb.tsv")
# Output prefix
out_prefix <- file.path(qc_dir, "candidate_kraken")

# Sanity checks
if (!file.exists(contig_file)) stop("Missing contig file: ", contig_file)
if (!file.exists(candidate_file)) stop("Missing contig file: ", candidate_file)

# Read files
# Contig-level annotation table
contigs <- read_tsv(contig_file, show_col_types = FALSE)
# Candidate contigs (no header, columns are positional)
candidates <- read_tsv(
  candidate_file,
  show_col_types = FALSE
) %>%
  select(contig)

if (nrow(candidates) == 0){
  stop("No candidates contigs found for ", SPECIES) 
}

# Intersection
# Keep only contigs that are both: non-host candidates (GC/coverage + length) and annotated by Kraken
cand_annot <- contigs %>%
  semi_join(candidates, by = "contig")

# Kraken contig class
cand_class_summary <- cand_annot %>%
  group_by(label) %>%
  summarise(
    n_contigs = n(),
    bp = sum(len, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(percent = 100 * bp / sum(bp))

write_tsv(
  cand_class_summary,
  paste0(out_prefix, "_class_bp_summary.tsv")
)

# Putative non-host contigs
# Dominant-taxon contigs grouped by Kraken taxid represent putative genomes (not guaranteed complete).
cand_genomes <- cand_annot %>%
  filter(label == "dominant_taxon") %>%
  group_by(top_taxid) %>%
  summarise(
    n_contigs = n(),
    bp = sum(len, na.rm = TRUE),
    mean_gc = mean(gc, na.rm = TRUE),
    mean_cov = mean(mean_cov, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(percent_bp = 100 * bp / sum(bp))

write_tsv(
  cand_genomes,
  paste0(out_prefix, "_putative_genomes.tsv")
)

# Done
message("Intersection analysis finished for ", SPECIES)
message("Outputs:")
message(" - ", paste0(out_prefix, "_class_bp_summary.tsv"))
message(" - ", paste0(out_prefix, "_putative_genomes.tsv"))
