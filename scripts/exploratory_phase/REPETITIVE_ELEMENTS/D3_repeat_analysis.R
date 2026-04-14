#!/usr/bin/env Rscript

# =============================================================================
# Stage D3: Unified repeat element analysis
#
# Consolidates EDTA GFF3 parsing, RepeatMasker + TRF parsing, blobplot
# generation and summary tables into a single script.
#
# Produces:
#   3 blobplots: TE fraction, non-TE repeat fraction, total repeat fraction
#   Per-contig repeat tables with category breakdown
#   Summary statistics linking repeat content to coverage heterogeneity
#
# Usage:
#   Rscript D3_repeat_analysis.R <SPECIES> <EDTA_GFF3> <RM_OUT> <TRF_DAT> \
#                                 <ASSEMBLY> <COVERAGE_TABLE> <OUTDIR> [--threshold 0.5]
#
# Inputs:
#   EDTA_GFF3       : .mod.EDTA.TEanno.gff3 from EDTA
#   RM_OUT          : RepeatMasker .out file
#   TRF_DAT         : TRF .dat file
#   ASSEMBLY        : Assembly FASTA (for contig lengths)
#   COVERAGE_TABLE  : coverage_classification.tsv from B1b
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: D3_repeat_analysis.R SPECIES EDTA_GFF3 RM_OUT TRF_DAT ASSEMBLY COVERAGE_TABLE OUTDIR [threshold]")
}

SPECIES        <- args[1]
EDTA_GFF3      <- args[2]
RM_OUT_FILE    <- args[3]
TRF_DAT_FILE   <- args[4]
ASSEMBLY       <- args[5]
COVERAGE_TABLE <- args[6]
OUTDIR         <- args[7]
TE_THRESHOLD   <- if (length(args) >= 8) as.numeric(args[8]) else 0.5

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Helper: merge overlapping intervals and return total bp
# =============================================================================

merge_and_sum <- function(starts, ends) {
  if (length(starts) == 0) return(0)
  ord <- order(starts)
  s <- starts[ord]; e <- ends[ord]
  ms <- s[1]; me <- e[1]; total <- 0
  if (length(s) > 1) {
    for (i in 2:length(s)) {
      if (s[i] <= me) {
        me <- max(me, e[i])
      } else {
        total <- total + (me - ms + 1)
        ms <- s[i]; me <- e[i]
      }
    }
  }
  total + (me - ms + 1)
}

compute_bp <- function(records, contigs) {
  by_ctg <- split(records, records$contig)
  sapply(contigs, function(ctg) {
    if (ctg %in% names(by_ctg)) {
      d <- by_ctg[[ctg]]
      merge_and_sum(d$start, d$end)
    } else { 0 }
  })
}

# =============================================================================
# 1) Contig lengths from FASTA
# =============================================================================

message("[INFO] Reading assembly: ", ASSEMBLY)
lines    <- readLines(ASSEMBLY)
hdr_idx  <- grep("^>", lines)
ctg_names <- sub("^>(\\S+).*", "\\1", lines[hdr_idx])
ends <- c(hdr_idx[-1] - 1, length(lines))
ctg_lens <- mapply(function(s, e) sum(nchar(lines[(s + 1):e])), hdr_idx, ends)
contig_len <- data.frame(contig = ctg_names, length = ctg_lens, stringsAsFactors = FALSE)
rm(lines); gc(verbose = FALSE)
message("[INFO] Loaded ", nrow(contig_len), " contigs")

# =============================================================================
# 2) Parse EDTA GFF3 for TE annotations
# =============================================================================

message("[INFO] Reading EDTA GFF3: ", EDTA_GFF3)
gff <- read.delim(EDTA_GFF3, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE,
                  col.names = c("seqid", "source", "type", "start", "end",
                                "score", "strand", "phase", "attributes"))

feat_counts <- sort(table(gff$type), decreasing = TRUE)
message("[INFO] Feature types in GFF3:")
print(feat_counts)

if ("repeat_region" %in% names(feat_counts)) {
  gff <- gff[gff$type == "repeat_region", ]
  message("[INFO] Using 'repeat_region' features (n = ", nrow(gff), ")")
} else {
  message("[WARN] 'repeat_region' not found; using all features")
}

edta_intervals <- data.frame(
  contig = gff$seqid, start = gff$start, end = gff$end,
  stringsAsFactors = FALSE
)

# =============================================================================
# 3) Parse RepeatMasker .out
# =============================================================================

message("[INFO] Reading RepeatMasker output: ", RM_OUT_FILE)
rm_lines <- readLines(RM_OUT_FILE)
rm_lines <- rm_lines[-(1:3)]
rm_lines <- rm_lines[nchar(trimws(rm_lines)) > 0]

rm_records <- data.frame(contig = character(), start = integer(),
                         end = integer(), class = character(),
                         stringsAsFactors = FALSE)

if (length(rm_lines) > 0) {
  parsed <- lapply(rm_lines, function(line) {
    f <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(f) >= 11) {
      data.frame(contig = f[5], start = as.integer(f[6]), end = as.integer(f[7]),
                 class = f[11], stringsAsFactors = FALSE)
    }
  })
  rm_records <- do.call(rbind, parsed[!sapply(parsed, is.null)])
}

# Classify RM hits into broad categories
rm_records$category <- ifelse(grepl("Simple_repeat", rm_records$class), "simple_repeat",
                       ifelse(grepl("Low_complexity", rm_records$class), "low_complexity",
                       ifelse(grepl("Satellite", rm_records$class), "satellite",
                       ifelse(grepl("rRNA|tRNA|snRNA|scRNA", rm_records$class), "rna",
                              "te"))))

message("[INFO] RepeatMasker hits: ", nrow(rm_records))

# =============================================================================
# 4) Parse TRF .dat
# =============================================================================

message("[INFO] Reading TRF output: ", TRF_DAT_FILE)
trf_lines <- readLines(TRF_DAT_FILE)

trf_records <- data.frame(contig = character(), start = integer(),
                          end = integer(), stringsAsFactors = FALSE)
current_ctg <- ""
for (line in trf_lines) {
  if (grepl("^Sequence:", line)) {
    current_ctg <- sub("^Sequence:\\s+(\\S+).*", "\\1", line)
  } else {
    f <- strsplit(trimws(line), "\\s+")[[1]]
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

# =============================================================================
# 5) Compute per-contig repeat fractions by category
# =============================================================================

message("[INFO] Computing per-contig repeat fractions...")

rm_te      <- rm_records[rm_records$category == "te", ]
rm_simple  <- rm_records[rm_records$category == "simple_repeat", ]
rm_lowc    <- rm_records[rm_records$category == "low_complexity", ]
rm_sat     <- rm_records[rm_records$category == "satellite", ]

result <- contig_len

# EDTA TE bp (from GFF3 directly)
result$edta_te_bp <- compute_bp(edta_intervals, result$contig)

# RM category bp
result$rm_te_bp             <- compute_bp(rm_te, result$contig)
result$simple_repeat_bp     <- compute_bp(rm_simple, result$contig)
result$low_complexity_bp    <- compute_bp(rm_lowc, result$contig)
result$satellite_bp         <- compute_bp(rm_sat, result$contig)
result$tandem_repeat_bp     <- compute_bp(trf_records, result$contig)

# Non-TE repeat bp: merge simple + low_complexity + satellite + TRF intervals
non_te_intervals <- rbind(
  rm_simple[, c("contig", "start", "end")],
  rm_lowc[, c("contig", "start", "end")],
  rm_sat[, c("contig", "start", "end")],
  trf_records[, c("contig", "start", "end")]
)
result$non_te_repeat_bp <- compute_bp(non_te_intervals, result$contig)

# Total repeat bp: merge ALL intervals (EDTA + RM + TRF) to avoid double-counting
message("[INFO] Computing total repeat coverage (merged EDTA + RM + TRF)...")
all_intervals <- rbind(
  edta_intervals[, c("contig", "start", "end")],
  rm_records[, c("contig", "start", "end")],
  trf_records[, c("contig", "start", "end")]
)
result$total_repeat_bp <- compute_bp(all_intervals, result$contig)

# Fractions
result$edta_te_fraction       <- result$edta_te_bp / result$length
result$non_te_repeat_fraction <- result$non_te_repeat_bp / result$length
result$total_repeat_fraction  <- result$total_repeat_bp / result$length

# Flagging
result$flagged_repetitive <- ifelse(result$total_repeat_fraction >= TE_THRESHOLD, "yes", "no")

result <- result[order(-result$length), ]

# =============================================================================
# 6) Merge with coverage classification
# =============================================================================

message("[INFO] Loading coverage classification...")
cov <- read_tsv(COVERAGE_TABLE, show_col_types = FALSE)

required_cols <- c("contig", "mean_cov", "gc", "len", "coverage_class")
missing <- setdiff(required_cols, colnames(cov))
if (length(missing) > 0) {
  stop("Missing columns in coverage file: ", paste(missing, collapse = ", "))
}

df <- cov %>%
  left_join(
    as_tibble(result) %>%
      select(contig, edta_te_fraction, non_te_repeat_fraction, total_repeat_fraction,
             flagged_repetitive, edta_te_bp, non_te_repeat_bp, total_repeat_bp),
    by = "contig"
  ) %>%
  mutate(
    edta_te_fraction       = coalesce(edta_te_fraction, 0),
    non_te_repeat_fraction = coalesce(non_te_repeat_fraction, 0),
    total_repeat_fraction  = coalesce(total_repeat_fraction, 0),
    flagged_repetitive     = coalesce(flagged_repetitive, "no")
  )

# =============================================================================
# 7) Blobplots (3 total)
# =============================================================================

df_plot <- df %>% filter(mean_cov > 0, len > 0, gc >= 13, gc <= 80)

make_repeat_blob <- function(data, frac_col, label, color_high, outfile) {

  data <- data %>% arrange(.data[[frac_col]])

  p <- ggplot(data, aes(x = gc, y = mean_cov, size = log10(len),
                         colour = .data[[frac_col]])) +
    geom_point(alpha = 0.7) +
    scale_y_log10() +
    scale_size(range = c(1, 8)) +
    scale_colour_gradient(
      low    = "grey80",
      high   = color_high,
      name   = label,
      limits = c(0, 1)
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
make_repeat_blob(
  df_plot, "edta_te_fraction", "TE fraction (EDTA)", "firebrick",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_te.png"))
)

# Plot 2: non-TE repeats (simple + low_complexity + satellite + tandem)
make_repeat_blob(
  df_plot, "non_te_repeat_fraction", "Non-TE repeat fraction", "darkorange",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_non_te_repeats.png"))
)

# Plot 3: total repeat content (all merged)
make_repeat_blob(
  df_plot, "total_repeat_fraction", "Total repeat fraction", "purple4",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_total_repeats.png"))
)

# =============================================================================
# 8) Summary tables
# =============================================================================

message("[INFO] Writing summary tables...")

# Full per-contig repeat table (all categories)
write.table(result, file = file.path(OUTDIR, "contig_repeat_coverage.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Repeat × coverage class cross-tabulation
repeat_vs_cov <- df %>%
  mutate(repeat_status = ifelse(flagged_repetitive == "yes", "high_repeat", "low_repeat")) %>%
  count(repeat_status, coverage_class) %>%
  group_by(repeat_status) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

write_tsv(repeat_vs_cov, file.path(OUTDIR, "repeat_coverage_distribution.tsv"))

# Genome-wide summary
genome_summary <- tibble(
  species                = SPECIES,
  n_contigs              = nrow(result),
  genome_te_pct          = 100 * sum(result$edta_te_bp) / sum(result$length),
  genome_simple_repeat_pct = 100 * sum(result$simple_repeat_bp) / sum(result$length),
  genome_tandem_repeat_pct = 100 * sum(result$tandem_repeat_bp) / sum(result$length),
  genome_total_repeat_pct  = 100 * sum(result$total_repeat_bp) / sum(result$length),
  n_flagged_repetitive   = sum(result$flagged_repetitive == "yes"),
  flagged_bp_pct         = 100 * sum(result$length[result$flagged_repetitive == "yes"]) / sum(result$length)
)

write_tsv(genome_summary, file.path(OUTDIR, "repeat_genome_summary.tsv"))

cat("\n=== Genome repeat summary ===\n")
print(genome_summary)
cat("\n=== Repeat status × coverage class ===\n")
print(repeat_vs_cov)

message("[OK] Repeat analysis complete. Output: ", OUTDIR)