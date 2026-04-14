#!/usr/bin/env Rscript

# =============================================================================
# Stage C2: Unified BUSCO contig analysis
#
# Consolidates contig-level BUSCO collapse, anchor validation, blobplots
# and summary tables into a single script.
#
# Produces:
#   6 blobplots (3 coverage-class + BUSCO overlay, 3 gradient by BUSCO/Mb)
#   Summary tables: gene-level, contig-level, bacteria, mixed, anchors,
#                   anchor×coverage, bacteria×coverage distributions
#
# Usage:
#   Rscript C2_busco_analysis.R <SPECIES> <BUSCO_HYM> <BUSCO_ARTH> <BUSCO_BACT> \
#                                <COVERAGE_TABLE> <OUTDIR>
#
# Inputs:
#   BUSCO_HYM / ARTH / BACT : full_table.tsv from each BUSCO lineage run
#   COVERAGE_TABLE           : coverage_classification.tsv from B1b
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: C2_busco_analysis.R SPECIES BUSCO_HYM BUSCO_ARTH BUSCO_BACT COVERAGE_TABLE OUTDIR")
}

SPECIES        <- args[1]
BUSCO_HYM      <- args[2]
BUSCO_ARTH     <- args[3]
BUSCO_BACT     <- args[4]
COVERAGE_TABLE <- args[5]
OUTDIR         <- args[6]

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1) Read BUSCO full_table.tsv — robust parser (handles broken TSV)
# =============================================================================

read_busco <- function(path, lineage) {

  raw <- readLines(path)
  raw <- raw[!grepl("^#", raw)]

  if (length(raw) == 0) {
    message("[INFO] No BUSCO hits for ", lineage)
    return(tibble())
  }

  # Fix escaped tabs
  raw <- gsub("\\\\t", "\t", raw)
  split <- strsplit(raw, "\t")

  df <- tibble(
    busco_id    = sapply(split, `[`, 1),
    status      = sapply(split, `[`, 2),
    sequence    = sapply(split, `[`, 3),
    description = sapply(split, function(x) {
      if (length(x) >= 10) x[10] else NA_character_
    })
  )

  df %>%
    filter(!is.na(status), status != "Missing") %>%
    transmute(
      contig = str_trim(sequence),
      busco_id,
      status,
      lineage = !!lineage,
      description
    )
}

message("[INFO] Reading BUSCO tables...")
busco <- bind_rows(
  read_busco(BUSCO_HYM,  "hymenoptera"),
  read_busco(BUSCO_ARTH, "arthropoda"),
  read_busco(BUSCO_BACT, "bacteria")
)

cat("\n=== BUSCO gene counts ===\n")
print(busco %>% count(lineage))
cat("=========================\n\n")

# =============================================================================
# 2) Contig-level aggregation
# =============================================================================

busco_contig <- busco %>%
  group_by(contig) %>%
  summarise(
    # Total
    n_busco = n(),

    # Per-lineage counts
    n_hym  = sum(lineage == "hymenoptera"),
    n_arth = sum(lineage == "arthropoda"),
    n_bact = sum(lineage == "bacteria"),

    has_hym  = n_hym > 0,
    has_arth = n_arth > 0,
    has_bact = n_bact > 0,

    # Status counts (global)
    n_complete   = sum(status == "Complete"),
    n_fragmented = sum(status == "Fragmented"),
    n_duplicated = sum(status == "Duplicated"),

    # Status per lineage
    hym_complete   = sum(lineage == "hymenoptera" & status == "Complete"),
    hym_fragmented = sum(lineage == "hymenoptera" & status == "Fragmented"),
    hym_duplicated = sum(lineage == "hymenoptera" & status == "Duplicated"),

    arth_complete   = sum(lineage == "arthropoda" & status == "Complete"),
    arth_fragmented = sum(lineage == "arthropoda" & status == "Fragmented"),
    arth_duplicated = sum(lineage == "arthropoda" & status == "Duplicated"),

    bact_complete   = sum(lineage == "bacteria" & status == "Complete"),
    bact_fragmented = sum(lineage == "bacteria" & status == "Fragmented"),
    bact_duplicated = sum(lineage == "bacteria" & status == "Duplicated"),

    # Functional annotations
    busco_ids          = paste(unique(busco_id), collapse = "; "),
    busco_lineages     = paste(unique(lineage), collapse = "; "),
    busco_descriptions = paste(unique(na.omit(description)), collapse = " | "),
    hym_functions      = paste(unique(description[lineage == "hymenoptera"]), collapse = " | "),
    arth_functions     = paste(unique(description[lineage == "arthropoda"]), collapse = " | "),
    bact_functions     = paste(unique(description[lineage == "bacteria"]), collapse = " | "),

    .groups = "drop"
  )

# =============================================================================
# 3) Merge with coverage classification
# =============================================================================

message("[INFO] Loading coverage classification...")
cov <- read_tsv(COVERAGE_TABLE, show_col_types = FALSE)

required_cols <- c("contig", "mean_cov", "gc", "len", "coverage_class")
missing <- setdiff(required_cols, colnames(cov))
if (length(missing) > 0) {
  stop("Missing columns in coverage file: ", paste(missing, collapse = ", "))
}

df <- cov %>%
  left_join(busco_contig, by = "contig") %>%
  mutate(
    across(starts_with("n_"),   ~coalesce(.x, 0L)),
    across(starts_with("has_"), ~coalesce(.x, FALSE))
  )

# =============================================================================
# 4) Anchor tier classification + coverage cross-tabulations
# =============================================================================

anchors <- df %>%
  select(contig, has_hym, has_arth, has_bact) %>%
  mutate(
    anchor_type = case_when(
      has_hym              ~ "strong_hymenoptera",
      !has_hym & has_arth  ~ "supportive_arthropoda",
      TRUE                 ~ "none"
    ),
    bacteria_signal = has_bact
  )

anchor_eval <- df %>%
  left_join(anchors %>% select(contig, anchor_type, bacteria_signal), by = "contig") %>%
  mutate(
    anchor_type     = coalesce(anchor_type, "none"),
    bacteria_signal = coalesce(bacteria_signal, FALSE)
  )

# Cross-tabulation: anchor type × coverage class
distribution_host <- anchor_eval %>%
  count(anchor_type, coverage_class) %>%
  group_by(anchor_type) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

# Cross-tabulation: bacteria signal × coverage class
distribution_bact <- anchor_eval %>%
  count(bacteria_signal, coverage_class) %>%
  group_by(bacteria_signal) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

cat("\n=== Anchor × coverage distribution ===\n")
print(distribution_host)
cat("\n=== Bacteria × coverage distribution ===\n")
print(distribution_bact)
cat("\n")

# =============================================================================
# 5) Blobplots — coverage class background + BUSCO overlay (3 plots)
# =============================================================================

df_plot <- df %>%
  filter(mean_cov > 0, len > 0, gc >= 13, gc <= 80)

make_blobplot_overlay <- function(data, flag_col, title, outfile) {

  highlight <- data %>% filter(.data[[flag_col]])

  p <- ggplot(data, aes(gc, mean_cov)) +
    geom_point(
      aes(size = log10(pmax(len, 1)), fill = coverage_class),
      shape = 21, color = "black", alpha = 0.4
    ) +
    geom_point(
      data = highlight,
      aes(size = log10(pmax(len, 1))),
      shape = 21, fill = NA, color = "black", stroke = 1.5
    ) +
    scale_y_log10() +
    scale_fill_manual(
      values = c(
        host_like        = "steelblue",
        coverage_outlier = "firebrick",
        ambiguous        = "grey70"
      ),
      na.value = "grey85"
    ) +
    theme_bw() +
    labs(
      title = title,
      x     = "GC (%)",
      y     = "Mean coverage",
      fill  = "Coverage class",
      size  = "log10(contig length)"
    )

  ggsave(outfile, p, width = 7, height = 5)
  message("[INFO] Written: ", outfile)
}

make_blobplot_overlay(
  df_plot, "has_hym",
  paste("Hymenoptera BUSCO contigs:", SPECIES),
  file.path(OUTDIR, "blobplot_hym_coverage.png")
)

make_blobplot_overlay(
  df_plot, "has_arth",
  paste("Arthropoda BUSCO contigs:", SPECIES),
  file.path(OUTDIR, "blobplot_arth_coverage.png")
)

make_blobplot_overlay(
  df_plot, "has_bact",
  paste("Bacteria BUSCO contigs:", SPECIES),
  file.path(OUTDIR, "blobplot_bact_coverage.png")
)

# =============================================================================
# 6) Blobplots — BUSCO density gradient (3 plots)
# =============================================================================

df_grad <- df_plot %>%
  mutate(
    hym_per_mb  = n_hym  / (len / 1e6),
    arth_per_mb = n_arth / (len / 1e6),
    bact_per_mb = n_bact / (len / 1e6)
  )

make_gradient_blob <- function(data, density_col, lineage_label, outfile) {

  data <- data %>% arrange(.data[[density_col]])

  vals <- data[[density_col]]
  cap  <- quantile(vals[vals > 0], 0.99, na.rm = TRUE)
  if (is.na(cap) || cap == 0) cap <- max(vals, na.rm = TRUE)

  p <- ggplot(data, aes(x = gc, y = mean_cov, size = log10(len),
                         colour = pmin(.data[[density_col]], cap))) +
    geom_point(alpha = 0.7) +
    scale_y_log10() +
    scale_size(range = c(1, 8)) +
    scale_colour_gradient(
      low  = "grey80",
      high = "steelblue",
      name = "BUSCO / Mb"
    ) +
    labs(
      title = paste0(lineage_label, " BUSCO density: ", SPECIES),
      x     = "GC (%)",
      y     = "Mean coverage (log10)",
      size  = "log10(contig length)"
    ) +
    theme_bw()

  ggsave(outfile, p, width = 8, height = 5, dpi = 300)
  message("[INFO] Written: ", outfile)
}

make_gradient_blob(
  df_grad, "hym_per_mb", "Hymenoptera",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_hym_gradient.png"))
)

make_gradient_blob(
  df_grad, "arth_per_mb", "Arthropoda",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_arth_gradient.png"))
)

make_gradient_blob(
  df_grad, "bact_per_mb", "Bacteria",
  file.path(OUTDIR, paste0(SPECIES, "_blobplot_bact_gradient.png"))
)

# =============================================================================
# 7) Write all summary tables
# =============================================================================

message("[INFO] Writing summary tables...")

# Gene-level: every BUSCO hit with contig, lineage, status, description
write_tsv(busco,
  file.path(OUTDIR, "busco_gene_level.tsv"))

# Contig-level: merged with coverage, all counts and annotations
write_tsv(df,
  file.path(OUTDIR, "busco_per_contig_summary.tsv"))

# Bacteria-positive contigs
write_tsv(df %>% filter(has_bact),
  file.path(OUTDIR, "busco_bacteria_contigs.tsv"))

# Mixed contigs (both hymenoptera + bacteria hits)
write_tsv(df %>% filter(has_hym & has_bact),
  file.path(OUTDIR, "busco_mixed_hym_bact.tsv"))

# Anchor tiers
write_tsv(anchors,
  file.path(OUTDIR, "busco_anchor_tiers.tsv"))

# Anchor × coverage cross-tabulation
write_tsv(distribution_host,
  file.path(OUTDIR, "busco_anchor_coverage_distribution.tsv"))

# Bacteria × coverage cross-tabulation
write_tsv(distribution_bact,
  file.path(OUTDIR, "busco_bacteria_coverage_distribution.tsv"))

# Full evaluation table (all contigs with anchor type + bacteria signal)
write_tsv(anchor_eval,
  file.path(OUTDIR, "busco_anchor_full_eval.tsv"))

# Blobplot counts summary
blobplot_summary <- df_plot %>%
  summarise(
    n_contigs      = n(),
    hym            = sum(has_hym),
    arth           = sum(has_arth),
    bact           = sum(has_bact),
    mixed_hym_bact = sum(has_hym & has_bact)
  )

write_tsv(blobplot_summary,
  file.path(OUTDIR, "blobplot_summary.tsv"))

cat("\n=== Blobplot summary ===\n")
print(blobplot_summary)

# BUSCO density summary
cat("\n=== BUSCO density summary ===\n")
cat(sprintf("Contigs with hym BUSCO:  %d\n",  sum(df$n_hym > 0)))
cat(sprintf("Contigs with arth BUSCO: %d\n",  sum(df$n_arth > 0)))
cat(sprintf("Contigs with bact BUSCO: %d\n",  sum(df$n_bact > 0)))
cat(sprintf("Contigs with no BUSCO:   %d\n",
  sum(df$n_hym == 0 & df$n_arth == 0 & df$n_bact == 0)))
cat("=============================\n")

message("[OK] BUSCO analysis complete. Output: ", OUTDIR)