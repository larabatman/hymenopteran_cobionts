#!/usr/bin/env Rscript

# Filter each species dfam table according to their blastn results
# Discard entire scaffolds or contigs that matched a Wolbachia reference genome with high identity over its length
# If 80% or more of the scaffold or contig matches a Wolbachia genome at more than 99% identity, discard the hits from the dfam table 
# Dfam tables with candidate nuwts and blast6 tables share the host scaffold name: qseqid in blast6 and col3 in dfam.tbl
# Usage: 
# Rscript 03c_filter_dfam.R raw_dfam.tbl filter_dir output_filtered.tbl discarded_hits.txt

library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

RAW_DFAM <- args[1]
FILTER_DIR <- args[2]
OUT_FILTER <- args[3] 
OUT_FLAGS <- args[4]

# Set contamination rule
# Percent identity to apply to alignment
PIDENT_MIN <- 99
# Percentage of the host scaffold that is covered by the hits
COV_MIN <- 80 

# Header rows for the BLAST format 6 output
blast6_cols <- c("qseqid", "sseqid", "pident", "length", "qlen","qstart", "qend", "sstart", "send", "evalue", "bitscore")
# Column simplification: c for character, d for double and i for integer
blast6_types <- "ccdiiiiiidd"

# Function to read a blast6 file and make it a table
# Each alignment receives one row, which carries the scaffold name with length and identity aligned
# Interval of alignment: qs-qe normalized to minus strands
# Also handles missing files, and the reference only or reference and mag cases
read_blast6 <- function(path, source_label) {
  # If there is no file, or if it is empty, then create an empty table with the right columns and types
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(tibble(qseqid = character(), pident = double(), qlen = integer(), qs = integer(), qe = integer(), source = character()))
  }
  # Pass one blast6 file 
  # Normalize minus strands: qstart = pmin(qstart, qend) and qe = pmax(qstart, qend)
  read_tsv(path, col_names = blast6_cols, col_types = blast6_types) %>%
    mutate(qs = pmin(qstart, qend), qe = pmax(qstart, qend), source = source_label) %>%
    select(qseqid, pident, qlen, qs, qe, source)
}

# Functio to merge overlapping intervals
# Return the total base pairs covered
# Sort by start and walk through conting to extend or close the interval
merge_and_sum <- function(starts, ends){
    if(length(starts) == 0) return(0)
    ord <- order(starts)
    s <- starts[ord]
    e <- ends[ord]
    current_start <- s[1]
    current_end <- e[1]
    total <- 0
    if(length(s) > 1){
        for(i in 2:length(s)) {
            if( s[i] <= current_end){ # overla with current interval: extend
                current_end <- max(current_end, e[i])
            } else { # no overlap: close the current interval and start new one
                total <- total + (current_end - current_start + 1)
                current_start <- s[i]
                current_end <- e[i]
            }
        }
    }
    total + (current_end - current_start + 1) # for the last interval
}

# Read blast6 files from mag or ref screens
refs <- read_blast6(file.path(FILTER_DIR, "assembly_vs_wolbachia_refs.blast6"), "Wolbachia_reference_db")
mag <- read_blast6(file.path(FILTER_DIR, "assembly_vs_MAG.blast6"), "Species_Wolbachia_MAG")

# Bind the datasets: apply the identity filter, merge and divide by scaffold length
one_scaffold <- bind_rows(refs, mag) %>%
# Filter rows below 99% out
  filter(pident >= PIDENT_MIN) %>% 
  # arrange by scaffold/contig
  arrange(qseqid) %>%
  group_by(qseqid) %>%
  # one for each contig
  summarise(scaffold_len = first(qlen), 
            n_alignments = n(),
            max_pident = max(pident),
            covered_bp = merge_and_sum(qs, qe),
            sources = paste(sort(unique(source)), collapse = ","),
            .groups = "drop") %>%
  mutate(pct_scaffold_covered = covered_bp / scaffold_len * 100)

flagged <- one_scaffold %>% filter(pct_scaffold_covered >= COV_MIN)
flagged_scaffolds <- flagged$qseqid

# Keep the discarded record
flagged %>%
  transmute(scaffold = qseqid, scaffold_len, covered_bp, pct_scaffold_covered = round(pct_scaffold_covered, 2), max_pident = round(max_pident, 2), n_alignments, filter = "CONTAMINATION", source = sources) %>%
  write_tsv(OUT_FLAGS)

# Read dfam table to apply decisions on the flagged contigs
# Plan text with readLines to pass the whole file as character vector
lines <- readLines(RAW_DFAM)
# Subset the comments and #
is_comment <- grepl("^#", lines) 
# Select col3 of every line:
columns <- strsplit(trimws(lines), "[[:space:]]+")
scaffold <- vapply(columns, function(f) if (length(f) >=3) f[3] else NA_character_, "")
# Keep comment lines and not the flagged scaffolds
keep <- is_comment | !(scaffold %in% flagged_scaffolds)
# Regenerate the filtered table
writeLines(lines[keep], OUT_FILTER)