#!/usr/bin/env Rscript
# filter_dfam.R
#
# This script filters each species dfam table according to their blastn results. 
# It is implemented such that it discards entire scaffolds that matched a Wolbachia reference genome with high identity over its length
# The aime is to judge each scaffold:taking the whole host scaffold, ask if it is an actual piece of living Wolbachia that got mis-assembled into the host assembly, using the evidence of how much of its length aligned to Wolbachia at ≥ 99% identity.
# If 80% or more of the scaffold is near-identical to a Wolbachia genome, it is a Wolbachia contig sitting and must be discarded from the dfam
#
# Usage
# Rscript filter_dfam.R <raw_dfam.tbl> <filter_dir> <out_filtered.tbl> <out_flags.txt>

library(readr)
library(dplyr)

# Contamination rule (Vancaester et al. 2025)
PIDENT_MIN <- 99   # percent identity, applied per alignment
COV_MIN    <- 80   # percent of the host scaffold covered, UNION of alignments

# Arguments
args <- commandArgs(trailingOnly = TRUE)

raw_dfam <- args[1] #dfam to filter
filter_dir <- args[2] #directory
out_filt <- args[3] 
out_flags <- args[4]

# Header rows for the BLAST format 6 output
blast6_cols  <- c("qseqid", "sseqid", "pident", "length", "qlen","qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast6_types <- "ccdiiiiiidd"# readr compact column-type specification: c character, d double, i integer

# Read one blast6 file and make it a table
# Gives one row per alignement, each carrying the host scaffold's name with its length and identity aligned
# The interval it covers is thus qs-qe
# qstart/qend are normalised so a minus-strand alignment is still an increasing interval.
# Turns one blast6 file into a uniform table with scaffold name, its full length, the identity of each alignement and the interval that alignment covers
# This allows to make a missing file behave like no evidence from a source and handles the reference-only and reference + MAG cases
read_blast6 <- function(path, source_label) {
  if (!file.exists(path) || file.info(path)$size == 0) { # if there is no file, or if it is empty, then return an empty tibble with the right columns and typess
    return(tibble(qseqid = character(), pident = double(), qlen = integer(), # this is needed for bind_rows
                  qs = integer(), qe = integer(), source = character()))
  }
  # Read the blast6 file
  # qs = pmin(qstart, qend) and qe = pmax(qstart, qend) is a normalization for the minus strand alignments
  read_tsv(path, col_names = blast6_cols, col_types = blast6_types) %>%
    mutate(qs = pmin(qstart, qend), qe = pmax(qstart, qend), source = source_label) %>%
    select(qseqid, pident, qlen, qs, qe, source)
}

# Union length of a set of intervals
# Given a set of intervals, how many distinct bases do they cover: count overlaping alignments once
# Walks the sorted intervals keeping a running block, closing it whenever the next interval starts beyond its end. This is what bedtools merge does.
union_len <- function(s, e) {
  tot <- 0; rs <- s[1]; re <- e[1]
  if (length(s) > 1) for (i in 2:length(s)) {
    if (s[i] > re + 1) { tot <- tot + (re - rs + 1); rs <- s[i]; re <- e[i] }
    else                 re <- max(re, e[i])
  }
  tot + (re - rs + 1)
}

# Read the two files: gather evidence from both screens, the general reference database and its individual own Wolbachia MAG
refs <- read_blast6(file.path(filter_dir, "assembly_vs_wolbachia_refs.blast6"), "Wolbachia_reference_db")
mag  <- read_blast6(file.path(filter_dir, "assembly_vs_MAG.blast6"), "Species_Wolbachia_MAG")

# Union of evidence from both sources. Identity is applied per alignment first so only high-identity blocks contribute to the covered fraction.
# Bind rows: the two tables become one, and MAG evidence for the same scaffold sit side by side 
# Convert many alignment rows into one row per scaffold: what frction of this scaffold is near-identical to Wolbachia?
# Identitiy filter
# Merge
# divide by scaffold length
per_scaffold <- bind_rows(refs, mag) %>%
  filter(pident >= PIDENT_MIN) %>% # filtering the rows belo 99% identity
  arrange(qseqid, qs, qe) %>% # sort by scaffold, then by start position
  group_by(qseqid) %>% # every calculation is one per scaffold
  summarise(scaffold_len = first(qlen), 
            n_alignments = n(),
            max_pident = max(pident),
            covered_bp = union_len(qs, qe),
            sources = paste(sort(unique(source)), collapse = ","),
            .groups = "drop") %>%
  mutate(pct_scaffold_covered = covered_bp / scaffold_len * 100)

flagged <- per_scaffold %>% filter(pct_scaffold_covered >= COV_MIN)
flagged_scaffolds <- flagged$qseqid

# Write a discarded record
flagged %>%
  transmute(scaffold = qseqid, scaffold_len, covered_bp, pct_scaffold_covered = round(pct_scaffold_covered, 2), max_pident = round(max_pident, 2), n_alignments, filter = "CONTAMINATION", source = sources) %>%
  write_tsv(out_flags)

# Filter dfam table to apply the decisions
# The dfam table is handled as plain text as to not change the output for downstream analysis
lines <- readLines(raw_dfam) # the whole file as character vector, one element per line including the # header blocks
is_comment <- grepl("^#", lines) # true when a line is part f the # header block 

# Return the scaffold name in field 3 from one dfam line or NA if the line is too short to have one
# trimws() strips leading whitespace and then "[[:space:]]+" splits on any whitespace, as the dfam format is column-aligned with padding rather than tab-separated
scaffold_of <- function(line) { 
  f <- strsplit(trimws(line), "[[:space:]]+")[[1]] 
  if (length(f) >= 3) f[3] else NA_character_ # seelct field 3 which is the query name, the host scaffold 
}

# One true/false per line of the file to know if it must be kept or not: 
# vapply walks the line indices with seq_along rather than the lines themselves. Every iteration returns exacly one logical value
keep <- vapply(seq_along(lines), function(i) {
  if (is_comment[i]) return(TRUE) # headers are always kept
  !(scaffold_of(lines[i]) %in% flagged_scaffolds) # a data line is kept only if its scaffold was not flagged as contination
}, logical(1))

writeLines(lines[keep], out_filt) # subset tby the logical vector and write 