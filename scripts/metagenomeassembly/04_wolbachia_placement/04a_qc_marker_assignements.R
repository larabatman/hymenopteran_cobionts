#!/usr/bin/env Rscript

# Quality Control of the marker proteins that were found with hmmscan
# For each MAG, we have their .tbl output of the score each of their proteins got from the profile library
# Control that we have a clear best protein and not an ambiguous case:
# For every MAG and OG pair, look at the best, the second-best protein, and how far we were from it
# Classify based on the bitscore ratio of second/best
# OK
# AMBIGUOUS: second best scores between 0.50 and 0.90
# DUP: second-best scores >= 0.90 of the best
# WEAK: best hit is above the hit threshold
# ABSENT: no hit at all, the gene missing and that column full of gaps
# Usage:
# Rscript 04a_qc_marker_assignments.R

library(dplyr)
library(readr)
# Directory with the hmmscan .tbl of the MAGs
scan_dir <- "placement_out/scan"
# File with the single-copy marker proteins
tree_ogs <- "tree_ogs.txt"
# Output directory
out_dir <- "placement_out/qc"
# Threshold for similarity
THRESHOLD <- 1e-5
# Ratio for the second best bitscore
DUP_RATIO <- 0.90
# Ratio for the second-best bitscore ambiguity
AMB_RATIO <- 0.50

# Collect the .tbl files as a character vector
table_files <- list.files(scan_dir, pattern = "\\.tbl$", full.names = TRUE)
# The expected ogs in the .tbl are the ones from the single-copy protein marker list
expected_ogs <- trimws(readLines(tree_ogs))
expected_ogs <- expected_ogs[expected_ogs != ""]

# Function to read one hmmscan .tbl and parse og, protein, evalue and bits
#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
#OG0000239            -          ctg1111_1            -            1.4e-05   19.9   6.8    0.0025   12.5   0.0   3.1   2   1   1   3   3   3   2 -
#OG0000730            -          ctg1111_1            -              0.039    8.3  10.3    0.0024   12.3   3.2   2.3   2   1   0   2   2   2   0 -
# .tbl is delimited with whitespaces, only the description field has spaces
# col1: target og
# col3: query protein
# col5: E-value
# col6: bitscore
read_hmmscan <- function(f) {
  # Take one file, pass it as one string for each line
  lines <- readLines(f)
  # Drop the headers and blank rows
  lines <- lines[!grepl("^#", lines) & trimws(lines) != ""]
  if (length(lines) == 0)
    return(tibble(og = character(), protein = character(), evalue = numeric(), bits = numeric()))
  # Split each line by whitespace
  p <- strsplit(trimws(lines), "\\s+")
  # Take elements 1, 3, 5 and 6 of each line 
  tibble(og = vapply(p, `[`, "", 1),
         protein = vapply(p, `[`, "", 3),
         evalue = as.numeric(vapply(p, `[`, "", 5)),
         bits = as.numeric(vapply(p, `[`, "", 6)))
}

# Function to rank hits within each og ans set the flags: 
qc_one_mag <- function(f) {
  mag  <- sub("\\.tbl$", "", basename(f))
  hits <- read_hmmscan(f)
  # Number the hits within each OG, ranked so the best two can be collected
  ranked <- hits %>%
    group_by(og) %>% # within each og
    arrange(desc(bits), .by_group = TRUE) %>% # sort by bitscore, do not ignore grouping
    mutate(rank = row_number()) %>%
    ungroup()
  # The best matching protein is the first in rank
  best <- ranked %>% filter(rank == 1) %>%
    select(og, best_protein = protein, best_bits = bits, best_evalue = evalue)
  # The second best matching protein is the second in rank
  second <- ranked %>% filter(rank == 2) %>%
    select(og, second_protein = protein, second_bits = bits)
  # Join tables with expected ogs, the best and the sceond best
  tibble(og = expected_ogs) %>%
    left_join(best, by = "og") %>%
    left_join(second, by = "og") %>%
    mutate(
      mag = mag,
      score_ratio = second_bits / best_bits,
      flag = case_when(
        is.na(best_bits) ~ "ABSENT", # if there is no best, no matching protein
        best_evalue > THRESHOLD ~ "WEAK", # if the best evalue is above the threshold, no similarity
        !is.na(score_ratio) & score_ratio >= DUP_RATIO ~ "DUP", # if the second best to best raio is higher than DUP_RATIO, the two matches are too close
        !is.na(score_ratio) & score_ratio >= AMB_RATIO ~ "AMBIGUOUS", # if the second best to best ratio is between AMB_RATIO and DUP_RATIO, it is not clear which match is the best
        TRUE ~ "OK" )) %>%
    select(mag, og, flag, best_protein, best_bits, best_evalue, second_protein, second_bits, score_ratio)
}
# Apply the qc_one_mag function for each .tbl file
# Returns a list of MAG tibbles where each row is a marker protein with bind_rows to sstack them together
assignment <- bind_rows(lapply(table_files, qc_one_mag))
write_tsv(assignment, file.path(out_dir, "assignment_qc.tsv"))

# Summary for each MAG: 
one_mag <- assignment %>%
  group_by(mag) %>%
  summarise(
    expected = n(),
    ok = sum(flag == "OK"),
    ambiguous = sum(flag == "AMBIGUOUS"),
    dup = sum(flag == "DUP"),
    weak = sum(flag == "WEAK"),
    absent = sum(flag == "ABSENT"),
    .groups = "drop")
write_tsv(one_mag, file.path(out_dir, "mag_qc.tsv"))

# Summary for each OG 
one_og <- assignment %>%
  group_by(og) %>%
  summarise(
    n_mags = n(),
    weak = sum(flag == "WEAK"),
    dup = sum(flag == "DUP"),
    ambiguous = sum(flag == "AMBIGUOUS"),
    ok = sum(flag == "OK"),
    .groups = "drop")
write_tsv(one_og, file.path(out_dir, "og_qc.tsv"))