#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript script.R SPECIES LINKS HOST_LIST OUTDIR CHROMSIZES")
}

SPECIES <- args[1]
LINKS <- args[2]
HOST_LIST <- args[3]
OUTDIR <- args[4]
CHROMSIZES <- args[5]

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# 1) Read inputs
# contig_links.tsv: ptg000405l      ptg000423l      1
links <- read_tsv(
  LINKS,
  col_names = c("contig1", "contig2", "n_links"),
  show_col_types = FALSE
)

host_contigs <- read_lines(HOST_LIST) # Known host genome pieces
# assembly.chrom.sizes: ptg000523l      19750
chrom_sizes <- read_tsv(
  CHROMSIZES,
  col_names = c("contig", "length"),
  show_col_types = FALSE
)


# 2) Expand symmetric link table
# contig1 contig2 n_links
#  A  B 10
# Becomes: 
# contig  n_links
# A 10
# B 10
# Transform columns into rowa
links_long <- links %>%
  pivot_longer(
    cols = c(contig1, contig2),
    names_to = "side",
    values_to = "contig"
  ) %>%
  select(contig, n_links)

# 3) Total links per contig
# Compute total contacts per contig, like contigA total_links = 200
# Group rows belonging to the same contig and add all contacts 
total_links <- links_long %>%
  group_by(contig) %>%
  summarise(
    total_links = sum(n_links),
    .groups = "drop"
  )

# 4) Host-link counts
# Compute host contacts, like contigA <-> host contigs = 160
# Count how many of above contacts ivnolve host contigs
# Check if the evaluated contig is in the host_anchor list
links_host <- links %>%
  mutate(
    contig1_is_host = contig1 %in% host_contigs,
    contig2_is_host = contig2 %in% host_contigs
  ) %>%
  filter(contig1_is_host | contig2_is_host) %>% # Keep rows where at least one contig is host
  pivot_longer(
    cols = c(contig1, contig2),
    names_to = "side",
    values_to = "contig"
  ) %>%
  select(contig, n_links)
# Sum host contacts
host_links <- links_host %>%
  group_by(contig) %>%
  summarise(
    host_links = sum(n_links),
    .groups = "drop"
  )

# 5) Merge base metrics
# COmbine tables to add host contacts and contig length
# contig  total_links host_links  length
hic_scores <- total_links %>%
  left_join(host_links, by = "contig") %>%
  mutate(
    host_links = replace_na(host_links, 0)
  ) %>%
  left_join(chrom_sizes, by = "contig")

# 6) Compute expected host fraction (length-based)
# Genome length:
total_genome_length <- sum(chrom_sizes$length, na.rm = TRUE)

# Host genome length
host_total_length <- chrom_sizes %>%
  filter(contig %in% host_contigs) %>%
  summarise(sum(length)) %>%
  pull()
# Expected random contact is the fraction of genome that represents the host
expected_host_fraction <- host_total_length / total_genome_length

# 7) Compute fractions and enrichment
# The enrichment is then the host_links_fraction/expected_host_fraction
# Threshold for sufficient contacts: at least 20, otherwise dat ais unreliable
hic_scores <- hic_scores %>%
  mutate(
    host_links_fraction = ifelse(
      total_links > 0,
      host_links / total_links,
      0
    ),
    host_enrichment = ifelse(
      total_links > 0,
      host_links_fraction / expected_host_fraction,
      NA_real_
    ),
    sufficient_contacts = total_links >= 20
  )

# 8) Classification 
# Categories
# <20 contacts: low support
# enrichment ≥ 2: host_enriched
# enrichment ≤ 0.5: host_depleted
# otherwise: neutral
hic_scores <- hic_scores %>%
  mutate(
    hic_structural_class = case_when(
      !sufficient_contacts ~ "low_support",
      host_enrichment >= 2 ~ "host_enriched",
      host_enrichment <= 0.5 ~ "host_depleted",
      TRUE ~ "neutral"
    )
  ) %>%
  arrange(desc(host_enrichment)) # Sort results ro have strongest host signal at top

# 9) Write output
out_file <- file.path(OUTDIR, "hic_host_link_enrichment.tsv")

write_tsv(hic_scores, out_file)

message("Hi-C structural enrichment scoring complete.")
message("Species: ", SPECIES)
message("Total contigs scored: ", nrow(hic_scores))
message("Expected host fraction: ", round(expected_host_fraction, 4))
message("Output: ", out_file)