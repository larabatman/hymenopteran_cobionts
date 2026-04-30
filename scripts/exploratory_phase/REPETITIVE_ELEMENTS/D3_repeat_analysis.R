#!/usr/bin/env Rscript

# =============================================================================
# Stage D3: Unified repeat element analysis
#
# This script consolidates all repeat data (EDTA, RepeatMasker, TRF) into a single per-contig analysis. 
# It: 
# - Reads contig lengths from the assembly FASTA
# - Parses EDTA GFF3 for TE coordinates
# - Parses RM .out for TE and non-TE repeat coordinates
# - Parses TRF .dat for tandem repeat coordinates
# - Merges overlapping intervals to avoid double-counting
# - Computes per-contig fractions for TEs, non-TEs repeats and totals
# - Joins with coverage data and produces blobplots and summary tables
#
# Usage:
#   Rscript D3_repeat_analysis.R <SPECIES> <EDTA_GFF3> <RM_OUT> <TRF_DAT> \
#                                 <ASSEMBLY> <COVERAGE_TABLE> <OUTDIR> [--threshold 0.5]

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# Parse CLI arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: D3_repeat_analysis.R SPECIES EDTA_GFF3 RM_OUT TRF_DAT ASSEMBLY COVERAGE_TABLE OUTDIR [threshold]")
}

SPECIES <- args[1] # species name for plot titles and filenames
EDTA_GFF3 <- args[2] # path to EDTA's TEanno.gff3
RM_OUT_FILE <- args[3] # path to RM's .out file
TRF_DAT_FILE <- args[4] # path to TRF's .dat file
ASSEMBLY <- args[5] # path to assembly FASTA for contig lengths
COVERAGE_TABLE <- args[6] # path to coverage_classification.tsv from B2
OUTDIR <- args[7] # output dir
TE_THRESHOLD <- if (length(args) >= 8) as.numeric(args[8]) else 0.5 # optional 8th argument, threshold for flagging contigs as repeptitive and default set to 0.5 ad 50% of the contig must be repeat-annotated for flagging

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Helper functions: merge overlapping intervals and return total bp

# merge_and_sum: core interval-merging step
# Repeat annotations from different tools overlap: if contig X has a TE from position 100-500 in EDTA, and randep repeat from 400-500 in TRF, adding them naively will give 400 + 200 = 600 bp. 
# But positions 400-500 are covered by both, and the real covered region is 100 - 600 = 500 bp! 
# To avoid that: sort intervals by start position, and walk through them. 
# If the next interval overlaps the current one, extend the current interval, otherwise close the current interval by adding its length to total and start a new one. 

# Classic merge overlapping intervals problem 
merge_and_sum <- function(starts, ends) {
  if (length(starts) == 0) return(0) # no intervals: 0 bp covered

  ord <- order(starts) # sort by start position
  s <- starts[ord]; e <- ends[ord] # reorder both vectors 

  ms <- s[1]; me <- e[1]; total <- 0 # initialize: first interval is the current one
  
  if (length(s) > 1) {
    for (i in 2:length(s)) {
      if (s[i] <= me) { # this interval overlaps or is adjacent to the current one, thus extend
        me <- max(me, e[i])
      } else { # no overlap: close current interval and add its length
        total <- total + (me - ms + 1)
        ms <- s[i]; me <- e[i] # start new current interval
      }
    }
  }
  total + (me - ms + 1) # add the las interval
}

# compute_bp: apply merge_and_sum per contig 
# Takes a dataframe of intervals with columns contig, start, end and a vector of contig names. 
# It returns a named vector of masked bp per contig. 
compute_bp <- function(records, contigs) {
  by_ctg <- split(records, records$contig) # split() groups the df by contig name to a list of sub-dataframes
  # For each contig in the assembly, look up its intervals and merge them 
  sapply(contigs, function(ctg) {
    if (ctg %in% names(by_ctg)) {
      d <- by_ctg[[ctg]]
      merge_and_sum(d$start, d$end)
    } else { 0 } # not intervals for that contig: 0 bp covered
  })
}

# 1) Contig lengths from FASTA

# Contig lengths are used as the denominators for all fraction calculations. 
# A FASTA file alternates between header lines and sequence lines. Read all lines, find the headers, and sum the character count of sequence lines between conscutive headers. 

message("[INFO] Reading assembly: ", ASSEMBLY)
lines <- readLines(ASSEMBLY) # read entire file into memory
hdr_idx <- grep("^>", lines) # indices of header lines
ctg_names <- sub("^>(\\S+).*", "\\1", lines[hdr_idx]) # extract contig name, the first word after ">"
# For each header, the sequence runs from header_line + 1 to next_header -1
# For the last contig, it runs to the end of the file 
ends <- c(hdr_idx[-1] - 1, length(lines)) 
# Sum character counts of all sequence lines for each contig
ctg_lens <- mapply(function(s, e) sum(nchar(lines[(s + 1):e])), hdr_idx, ends)
# Put in df
contig_len <- data.frame(contig = ctg_names, length = ctg_lens, stringsAsFactors = FALSE)
# Free the large character vector from memory!
rm(lines); gc(verbose = FALSE)
message("[INFO] Loaded ", nrow(contig_len), " contigs")

# 2) Parse EDTA GFF3 for TE annotations

# EDTA's GFF3 contains two kinds of annotations: 
# - structurally intact LTR elements found de novo:
# repeat_ region (parent)
#   Gypsy_LTR_retrotransposon (actual TE)
#     long_terminal_repeat (LTR arm, subfeature)
#     long_termina_repeat (LTR arm, subfeature)
#     target_site_duplication (TSD, sub-feature)
#     target_site_duplication (TSD, sub-feature)
# - homogoly-based annotations from RepeatMasker --anno 1:
# Gypsy_LTR_retrotransposon (no parent)
# Jockey_LINE_retrotransposon (no parent)
# repeat_fragment (no parent)
# Filtering by keeping repeat_region thus discarded all of the homology-based annotations that actually count for the vast majority of annotations. 
# Keep everything except sub-features that cause double counting: long_terminal_repeat, target_site_duplication, repeat_region
# All remaining features are individual TE annotations, and overlaps are resolved by merge_and_sum. 

message("[INFO] Reading EDTA GFF3: ", EDTA_GFF3)
gff <- read.delim(EDTA_GFF3, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE,
                  col.names = c("seqid", "source", "type", "start", "end",
                                "score", "strand", "phase", "attributes"))
# Print feature type distirbution for file verification and logging
feat_counts <- sort(table(gff$type), decreasing = TRUE)
message("[INFO] Feature types in GFF3:")
print(feat_counts)

# Filter to exclude sub-features instead of keeping only repeat_region
# repeat_region: wrapper whose coordinates overlap the child TE
# long_terminal_repeat: LTR arms already within the parent TE's span
# target_site_duplication: TSDs already within the parent TE's span 
exclude_types <- c("long_terminal_repeat", "target_site_duplication", "repeat_region")
n_before <- nrow(gff)
gff <- gff[!gff$type %in% exclude_types, ]
message("[INFO] Excluded sub-features (", paste(exclude_types, collapse = ", "), ")")
message("[INFO] Using ", nrow(gff), " TE features (from ", n_before, " total)")

# Extract columns needed for interval operations
edta_intervals <- data.frame(
  contig = gff$seqid, start = gff$start, end = gff$end,
  stringsAsFactors = FALSE
)

# 3) Parse RepeatMasker .out
# the .out file has 3 header lines, then one line per repeat hit. Each line is whitespace-delimited with key fields:
# - field 5: query sequence name, contig
# field 6: dtart position in query
# field 6: end position in query
# field 11: repeat class/ family (DNA/hAT, Simple_repeat, LTR/Gypsy, ...) 

message("[INFO] Reading RepeatMasker output: ", RM_OUT_FILE)
rm_lines <- readLines(RM_OUT_FILE)
# Skip the 3 header lines with column labels and blank lines
rm_lines <- rm_lines[-(1:3)]
rm_lines <- rm_lines[nchar(trimws(rm_lines)) > 0]
# initialize empty df in case the parsing is empty
rm_records <- data.frame(contig = character(), start = integer(),
                         end = integer(), class = character(),
                         stringsAsFactors = FALSE)

if (length(rm_lines) > 0) {
  # Parse each line by splitting on whitespace
  parsed <- lapply(rm_lines, function(line) {
    f <- strsplit(trimws(line), "\\s+")[[1]] # split on any whitespace
    if (length(f) >= 11) {
      data.frame(contig = f[5], start = as.integer(f[6]), end = as.integer(f[7]),
                 class = f[11], stringsAsFactors = FALSE)
    }
    # Lineas wiht fewer than 11 fields are silently skipped (malformed rows)
  })
  # Combine all parsed rows into one df, dropping NULLs
  rm_records <- do.call(rbind, parsed[!sapply(parsed, is.null)])
}

# Classify RM hits into broad categories based on the class string. This allows for separation between TE hits and non-TE repeat hits. 
# "DNA/hAT-Ac" contains non of the non-TE eywords and goes to "te"
# "Simple_repeat" metches and goes to "simple_repeat"
# "Low_complexity" matches and goes to "low_complexity"
# "Satellite/centr" matches and goes to "satellite"
# "rRNA" matches and goes to "rna"
rm_records$category <- ifelse(grepl("Simple_repeat", rm_records$class), "simple_repeat",
                       ifelse(grepl("Low_complexity", rm_records$class), "low_complexity",
                       ifelse(grepl("Satellite", rm_records$class), "satellite",
                       ifelse(grepl("rRNA|tRNA|snRNA|scRNA", rm_records$class), "rna",
                              "te"))))

message("[INFO] RepeatMasker hits: ", nrow(rm_records))

# 4) Parse TRF .dat
# TRF's .dat format alternates between header and data lines:
# Sequence: contig_name tells which contig
# ... (parameter lines)
# 100 250 15 10.0 15 85 5 ... is a data line with start = 100, end = 250, period = 15 etc
# 300 500 20 9.5 20 90 3 is another randem repeat, on the same contig
# Sequence: next_contig marks the start of a new contig 
# Track which contig we are on with current_ctg and extract start/ end from lines that begin with integers
message("[INFO] Reading TRF output: ", TRF_DAT_FILE)
trf_lines <- readLines(TRF_DAT_FILE)

trf_records <- data.frame(contig = character(), start = integer(),
                          end = integer(), stringsAsFactors = FALSE)
current_ctg <- ""

for (line in trf_lines) {
  if (grepl("^Sequence:", line)) {
    # Extract contig name from "Seqnece: contig_name" line
    current_ctg <- sub("^Sequence:\\s+(\\S+).*", "\\1", line)
  } else {
    f <- strsplit(trimws(line), "\\s+")[[1]]
    # Data lines start with two integers, the start and end positions
    if (length(f) >= 2 && !is.na(suppressWarnings(as.integer(f[1])))) {
      s <- as.integer(f[1]); e <- as.integer(f[2])
      if (!is.na(s) && !is.na(e) && nchar(current_ctg) > 0) {
        trf_records <- rbind(trf_records,
          data.frame(contig = current_ctg, start = s, end = e,
                     stringsAsFactors = FALSE))
      }
    }
  }
}
message("[INFO] TRF hits: ", nrow(trf_records))

# 5) Compute per-contig repeat fractions by category
# The three data sources get merged into one per-contig table, computing the coverage at three levels:
# - EDTA TEs only from the GFF3
# - Non-TE repeats only (simple + low-complexity + satellite + TRF)
# - Total: everything merged, overlaps resolved
message("[INFO] Computing per-contig repeat fractions...")

# Split RM hits by category
rm_te      <- rm_records[rm_records$category == "te", ]
rm_simple  <- rm_records[rm_records$category == "simple_repeat", ]
rm_lowc    <- rm_records[rm_records$category == "low_complexity", ]
rm_sat     <- rm_records[rm_records$category == "satellite", ]

# Start building the result table from contig lengths
result <- contig_len

# EDTA TE bp: how many bases per contig are annotated as TEs by EDTA
result$edta_te_bp <- compute_bp(edta_intervals, result$contig)

# RM category-specific bp
result$rm_te_bp             <- compute_bp(rm_te, result$contig)
result$simple_repeat_bp     <- compute_bp(rm_simple, result$contig)
result$low_complexity_bp    <- compute_bp(rm_lowc, result$contig)
result$satellite_bp         <- compute_bp(rm_sat, result$contig)
result$tandem_repeat_bp     <- compute_bp(trf_records, result$contig)

# Non-TE repeat bp: merge simple + low_complexity + satellite + TRF intervals
# Avoid double-counting where a TRF randem repeats overlaps with a Simple_repeat from RM annotation on the same locus, for instance
non_te_intervals <- rbind(
  rm_simple[, c("contig", "start", "end")],
  rm_lowc[, c("contig", "start", "end")],
  rm_sat[, c("contig", "start", "end")],
  trf_records[, c("contig", "start", "end")]
)
result$non_te_repeat_bp <- compute_bp(non_te_intervals, result$contig)

# Total repeat bp: merge ALL intervals (EDTA + RM + TRF) to avoid double-counting
# This tells what fraction of each contig is covered by any kind of repetitive element. Merging ensures that if EDTA, RM nad TRF all annnotate the same 1000 bp region, it is counted once and not thrice!
message("[INFO] Computing total repeat coverage (merged EDTA + RM + TRF)...")
all_intervals <- rbind(
  edta_intervals[, c("contig", "start", "end")],
  rm_records[, c("contig", "start", "end")],
  trf_records[, c("contig", "start", "end")]
)
result$total_repeat_bp <- compute_bp(all_intervals, result$contig)

# Compute the fractions: the proportion of contig length covered by repeats
result$edta_te_fraction <- result$edta_te_bp / result$length
result$non_te_repeat_fraction <- result$non_te_repeat_bp / result$length
result$total_repeat_fraction <- result$total_repeat_bp / result$length

# Flagging contigs where repeats content exceeds the threshold (default 0.5, 50%)
result$flagged_repetitive <- ifelse(result$total_repeat_fraction >= TE_THRESHOLD, "yes", "no")

result <- result[order(-result$length), ]

# 6) Merge with coverage classification
# The coverage classification table from B2 has per-contig: mean_cov, gc, len, coverage_clas (host_like/outlier/ambiguous)
# Merging: lets corss-tabulate repeat content with coverage behaviour
message("[INFO] Loading coverage classification...")
cov <- read_tsv(COVERAGE_TABLE, show_col_types = FALSE)

# Verify required columns exist
required_cols <- c("contig", "mean_cov", "gc", "len", "coverage_class")
missing <- setdiff(required_cols, colnames(cov))
if (length(missing) > 0) {
  stop("Missing columns in coverage file: ", paste(missing, collapse = ", "))
}

# Left-join: keep all contigs from coverage table and add repeat columns
# coalesce(..., 0) replaces NA, contigs with no repeat annotation, with 0
df <- cov %>%
  left_join(
    as_tibble(result) %>%
      select(contig, edta_te_fraction, non_te_repeat_fraction, total_repeat_fraction,
             flagged_repetitive, edta_te_bp, non_te_repeat_bp, total_repeat_bp),
    by = "contig"
  ) %>%
  mutate(
    edta_te_fraction = coalesce(edta_te_fraction, 0),
    non_te_repeat_fraction = coalesce(non_te_repeat_fraction, 0),
    total_repeat_fraction = coalesce(total_repeat_fraction, 0),
    flagged_repetitive = coalesce(flagged_repetitive, "no")
  )

# 7) Blobplots
# Three blobplots: each coloured by a different repeat metric. All share same structure: 
# GC on x, coverage on y, point size = contig length and gradient shows repeat fraction from 0 in grey to high, coloured. 

# Filter out problematic contigs with 0 coverage or length and extreme GC 
df_plot <- df %>% filter(mean_cov > 0, len > 0, gc >= 13, gc <= 80)

# Reustable blobplot function
make_repeat_blob <- function(data, frac_col, label, color_high, outfile) {
  # Sort by the fraction column so hig.repeat contigs are drawn last, on top, which prevents them from being burried under grey mass
  data <- data %>% arrange(.data[[frac_col]])

  p <- ggplot(data, aes(x = gc, y = mean_cov, size = log10(len),
                         colour = .data[[frac_col]])) +
    geom_point(alpha = 0.7) + # semi-transparent points
    scale_y_log10() + # log-scale for coverage axis
    scale_size(range = c(1, 8)) + # map conti length to point size
    scale_colour_gradient(
      low    = "grey80", # 0% repeat: light grey
      high   = color_high, # high repeat: strong colour
      name   = label,
      limits = c(0, 1) # fix sxale from 0 to 1 (0% to 100%)
    ) +
    labs(
      title = paste0(label, ": ", SPECIES),
      x     = "GC (%)",
      y     = "Mean coverage (log10)",
      size  = "log10(contig length)"
    ) +
    theme_bw()

  ggsave(outfile, p, width = 8, height = 5, dpi = 300)
  message("[INFO] Written: ", outfile)
}

# Plot 1: TE fraction (EDTA)
# Shows where TE concentrate, red = high TE
make_repeat_blob(
  df_plot, "edta_te_fraction", "TE fraction (EDTA)", "firebrick",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_te.png"))
)

# Plot 2: non-TE repeats (simple + low_complexity + satellite + tandem)
# Shows where non-TE repetitive elements concentrate, orange = high non-TE repeat
make_repeat_blob(
  df_plot, "non_te_repeat_fraction", "Non-TE repeat fraction", "darkorange",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_non_te_repeats.png"))
)

# Plot 3: total repeat content (all merged)
# Combined picture, purple = hish total repeat content
make_repeat_blob(
  df_plot, "total_repeat_fraction", "Total repeat fraction", "purple4",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_total_repeats.png"))
)

# 8) Summary tables

message("[INFO] Writing summary tables...")

# Table 1: full per-contig repeat table (all categories, all bp counts)
# Master table, one row per contig, all metrics
write.table(result, file = file.path(OUTDIR, "contig_repeat_coverage.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Table 2: Repeat × coverage class cross-tabulation
# To know if coverage outlier contigs are predominantly repetitive. If yes, supports the argument that coverage heterogeneity are repeats rather than cobionts
repeat_vs_cov <- df %>%
  mutate(repeat_status = ifelse(flagged_repetitive == "yes", "high_repeat", "low_repeat")) %>%
  count(repeat_status, coverage_class) %>% # count contigs per combination
  group_by(repeat_status) %>%
  mutate(percent = n / sum(n) * 100) %>% # percentage within each repeat status
  ungroup()

write_tsv(repeat_vs_cov, file.path(OUTDIR, "repeat_coverage_distribution.tsv"))

# Table 3: Genome-wide summary, one row per species
# Quick overview numbers for the report 
genome_summary <- tibble(
  species = SPECIES,
  n_contigs = nrow(result),
  genome_te_pct = 100 * sum(result$edta_te_bp) / sum(result$length),
  genome_simple_repeat_pct = 100 * sum(result$simple_repeat_bp) / sum(result$length),
  genome_tandem_repeat_pct = 100 * sum(result$tandem_repeat_bp) / sum(result$length),
  genome_total_repeat_pct = 100 * sum(result$total_repeat_bp) / sum(result$length),
  n_flagged_repetitive = sum(result$flagged_repetitive == "yes"),
  flagged_bp_pct = 100 * sum(result$length[result$flagged_repetitive == "yes"]) / sum(result$length)
)

write_tsv(genome_summary, file.path(OUTDIR, "repeat_genome_summary.tsv"))
# Print to stdout for the log
cat("\n=== Genome repeat summary ===\n")
print(genome_summary)
cat("\n=== Repeat status × coverage class ===\n")
print(repeat_vs_cov)

message("[OK] Repeat analysis complete. Output: ", OUTDIR)