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
#   COVERAGE_TABLE           : coverage_classification.tsv from B2

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: C2_busco_analysis.R SPECIES BUSCO_HYM BUSCO_ARTH BUSCO_BACT COVERAGE_TABLE OUTDIR")
}

SPECIES <- args[1]
BUSCO_HYM <- args[2]
BUSCO_ARTH <- args[3]
BUSCO_BACT <- args[4]
COVERAGE_TABLE <- args[5]
OUTDIR <- args[6]

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# 1) Read BUSCO full_table.tsv — robust parser (handles broken TSV)

# Custom parser instead of read_tsv() because BUSCO full_table.tsv files have a variable number of comment lines starting with "#" at the top that read_tsv() cannot handle cleanly
# readLines() lets us strip those before any parsing happens.
read_busco <- function(path, lineage) {

  raw <- readLines(path)

  # Drop every line beginning with "#" (BUSCO header/comment lines)
  raw <- raw[!grepl("^#", raw)]

  if (length(raw) == 0) {
    message("[INFO] No BUSCO hits for ", lineage)
    return(tibble())
  }

  # Some BUSCO versions write the literal two-character string "\t" instead of a real tab character in certain fields. Replace those before splitting on tabs
  # otherwise column boundaries will be wrong for affected rows.
  raw <- gsub("\\\\t", "\t", raw)

  # Split each line on tab to get a list of character vectors (one per line)
  split <- strsplit(raw, "\t")

  # Pull out the four columns by position; sapply(..., `[`, N) extracts the Nth element from each vector in the list.
  # Column 10 is the gene description, not always present in every BUSCO version, so check length before accessing to avoid out-of-bounds NAs.
  df <- tibble(
    busco_id = sapply(split, `[`, 1),
    status = sapply(split, `[`, 2),
    sequence = sapply(split, `[`, 3),
    description = sapply(split, function(x) {
      if (length(x) >= 10) x[10] else NA_character_
    })
  )

  # Drop genes with no status and genes scored "Missing": those were not placed on any contig and carry no positional information.
  # str_trim() strips whitespace from contig names to prevent silent join failures later (a trailing space would make "ctg1 " != "ctg1").
  # transmute() keeps only the renamed/derived columns, dropping the originals
  # !!lineage forces the function argument to be evaluated as a value rather than treated as a column name inside tidy evaluation.
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

# Read all three lineages and stack into a single gene-level table.
# bind_rows() handles the case where any tibble is empty (no hits).
busco <- bind_rows(
  read_busco(BUSCO_HYM, "hymenoptera"),
  read_busco(BUSCO_ARTH, "arthropoda"),
  read_busco(BUSCO_BACT, "bacteria")
)

cat("\n=== BUSCO gene counts ===\n")
print(busco %>% count(lineage))
cat("=========================\n\n")

# 2) Contig-level aggregation

# The raw BUSCO table has one row per gene per contig. A contig carrying 50 BUSCO genes has 50 rows, so collapse to one row per contig to join with the coverage table (also one row per contig)
# Each contig gets a single classification in the anchor tier step below.
busco_contig <- busco %>%
  group_by(contig) %>%
  summarise(

    # Total BUSCO hits across all lineages on this contig
    n_busco = n(),

    # Per-lineage hit counts. A contig with n_hym = 30 carries 30 Hymenoptera-conserved genes, strong host signal. A contig with n_hym = 2 and n_bact = 1 is a chimera candidate worth inspecting.
    n_hym = sum(lineage == "hymenoptera"),
    n_arth = sum(lineage == "arthropoda"),
    n_bact = sum(lineage == "bacteria"),

    # Boolean presence flags derived directly from the counts above
    has_hym = n_hym > 0,
    has_arth = n_arth > 0,
    has_bact = n_bact > 0,

    # Global status counts across all lineages combined
    n_complete = sum(status == "Complete"),
    n_fragmented = sum(status == "Fragmented"),
    n_duplicated = sum(status == "Duplicated"),

    # Per-lineage status breakdown — useful for spotting fragmented host signal or duplicated bacterial hits separately from global counts
    hym_complete = sum(lineage == "hymenoptera" & status == "Complete"),
    hym_fragmented = sum(lineage == "hymenoptera" & status == "Fragmented"),
    hym_duplicated = sum(lineage == "hymenoptera" & status == "Duplicated"),

    arth_complete = sum(lineage == "arthropoda" & status == "Complete"),
    arth_fragmented = sum(lineage == "arthropoda" & status == "Fragmented"),
    arth_duplicated = sum(lineage == "arthropoda" & status == "Duplicated"),

    bact_complete = sum(lineage == "bacteria" & status == "Complete"),
    bact_fragmented = sum(lineage == "bacteria" & status == "Fragmented"),
    bact_duplicated = sum(lineage == "bacteria" & status == "Duplicated"),

    # Collapse all gene annotations for this contig into single strings.
    # unique() prevents the same name appearing twice if a gene was detected n multiple lineages. 
    # na.omit() drops missing descriptions before pasting so we don't get the literal string "NA" in the output.
    busco_ids = paste(unique(busco_id), collapse = "; "),
    busco_lineages = paste(unique(lineage), collapse = "; "),
    busco_descriptions = paste(unique(na.omit(description)), collapse = " | "),
    hym_functions = paste(unique(description[lineage == "hymenoptera"]), collapse = " | "),
    arth_functions = paste(unique(description[lineage == "arthropoda"]), collapse = " | "),
    bact_functions = paste(unique(description[lineage == "bacteria"]), collapse = " | "),

    .groups = "drop"
  )

# 3) Merge with coverage classification

message("[INFO] Loading coverage classification...")
cov <- read_tsv(COVERAGE_TABLE, show_col_types = FALSE)

# Fail early with a clear message if the coverage table is missing expected columns — easier to debug than a cryptic join error further downstream
required_cols <- c("contig", "mean_cov", "gc", "len", "coverage_class")
missing <- setdiff(required_cols, colnames(cov))
if (length(missing) > 0) {
  stop("Missing columns in coverage file: ", paste(missing, collapse = ", "))
}

# left_join keeps ALL contigs from the coverage table, including those with no BUSCO hits at all. Contigs absent from busco_contig will have NA in every BUSCO column after the join. 
# Then replace those NAs explicitly:
#   - count columns (n_*): NA → 0L  ("no hits found", not "data missing")
#   - boolean flags (has_*): NA → FALSE
# coalesce(x, replacement) returns x where non-NA, replacement where NA.
# across(starts_with("n_"), ...) applies the same coalesce to every column
# whose name starts with "n_", avoiding repetitive column-by-column code.
df <- cov %>%
  left_join(busco_contig, by = "contig") %>%
  mutate(
    across(starts_with("n_"), ~coalesce(.x, 0L)),
    across(starts_with("has_"), ~coalesce(.x, FALSE))
  )

# 4) Anchor tier classification + coverage cross-tabulations

# Anchor tiers implement a hierarchy of host-evidence strength based on which BUSCO lineages placed genes on each contig:
# strong_hymenoptera: carries Hymenoptera BUSCOs. Hymenoptera is a narrow lineage, so these genes are highly specific to the host. Near-certain host sequence.
# supportive_arthropoda: carries Arthropoda but NOT Hymenoptera BUSCOs. Consistent with host but less specific — could also be another arthropod contaminant.
# none: no eukaryotic BUSCO signal; uninformative for host classification from gene content alone.

# bacteria_signal is kept separate rather than folded into anchor_type because bacterial hits are not a lower tier of host evidence.
# They are a categorically different signal indicating potential contamination.
# The TRUE ~ "none" line is the catch-all (equivalent to an else clause).
anchors <- df %>%
  select(contig, has_hym, has_arth, has_bact) %>%
  mutate(
    anchor_type = case_when(
      has_hym ~ "strong_hymenoptera",
      !has_hym & has_arth  ~ "supportive_arthropoda",
      TRUE ~ "none"
    ),
    bacteria_signal = has_bact
  )

# Join anchor classifications back onto the full contig table.
# coalesce() again handles any contigs that were missing from anchors (edge case guard that should not occur, but prevents silent NAs).
anchor_eval <- df %>%
  left_join(anchors %>% select(contig, anchor_type, bacteria_signal), by = "contig") %>%
  mutate(
    anchor_type = coalesce(anchor_type, "none"),
    bacteria_signal = coalesce(bacteria_signal, FALSE)
  )

# 2D contingency tables: anchor tier x coverage class, and bacteria x coverage.
# These answer the key diagnostic question: do host-anchored contigs fall in the host_like coverage class from B2? If strong_hymenoptera contigs are mostly host_like, the coverage model is well-calibrated. 
# If many fall in coverage_outlier, the assembly may have structural issues worth investigating.
# group_by + mutate (not summarise) keeps all rows while adding a
# within-group percentage column alongside the raw count.
distribution_host <- anchor_eval %>%
  count(anchor_type, coverage_class) %>%
  group_by(anchor_type) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

distribution_bact <- anchor_eval %>%
  count(bacteria_signal, coverage_class) %>%
  group_by(bacteria_signal) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

cat("\n=== Anchor x coverage distribution ===\n")
print(distribution_host)
cat("\n=== Bacteria x coverage distribution ===\n")
print(distribution_bact)
cat("\n")

# 5) Blobplots — coverage class background + BUSCO overlay (3 plots)

# Pre-filter before plotting: remove zero-coverage and zero-length contigs (cannot be log-scaled) and restrict GC to 13–80% to exclude biologically implausible values that are almost certainly assembly artifacts
df_plot <- df %>%
  filter(mean_cov > 0, len > 0, gc >= 13, gc <= 80)

# Produces a blobplot (GC vs coverage) with two layers:
# Layer 1: all contigs, filled by coverage_class from B2
# Layer 2: BUSCO-positive contigs only, drawn as an unfilled ring on top

# The two-layer design lets you see the coverage class color underneath AND the BUSCO signal simultaneously
# useful for asking whether BUSCO-positive contigs land in the expected coverage class.
# flag_col is passed as a string and accessed via .data[[flag_col]], which is the tidy evaluation idiom for using a variable as a column name inside aes()
make_blobplot_overlay <- function(data, flag_col, title, outfile) {

  # Subset to only the BUSCO-positive contigs for the highlight layer
  highlight <- data %>% filter(.data[[flag_col]])

  p <- ggplot(data, aes(gc, mean_cov)) +

    # Layer 1: all contigs colored by coverage class.
    # shape = 21 is a filled circle with an independent border color, lets fill (interior) and color (border) be set separately, so points stay visible when they overlap. alpha = 0.4 adds transparency.
    geom_point(
      aes(size = log10(pmax(len, 1)), fill = coverage_class),
      shape = 21, color = "black", alpha = 0.4
    ) +

    # Layer 2: BUSCO-positive contigs drawn as unfilled rings (fill = NA) with a thick border (stroke = 1.5) on top of layer 1. 
    # This highlights them without covering the coverage class color underneath.
    # pmax(len, 1) guards against log10(0) for any zero-length edge cases.
    geom_point(
      data = highlight,
      aes(size = log10(pmax(len, 1))),
      shape = 21, fill = NA, color = "black", stroke = 1.5
    ) +

    # Log scale on Y: without this, a cobiont at 500x would compress the host cluster at 40x into an unreadable band at the bottom of the plot
    scale_y_log10() +

    scale_fill_manual(
      values = c(
        host_like = "steelblue",
        coverage_outlier = "firebrick",
        ambiguous = "grey70"
      ),
      na.value = "grey85"
    ) +
    theme_bw() +
    labs(
      title = title,
      x = "GC (%)",
      y = "Mean coverage",
      fill = "Coverage class",
      size = "log10(contig length)"
    )

  ggsave(outfile, p, width = 7, height = 5)
  message("[INFO] Written: ", outfile)
}

# Call the function once per lineage. Same plot logic, different flag column
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

# 6) Blobplots — BUSCO density gradient (3 plots)

# Normalize BUSCO hit counts by contig length to get hits per megabase.
# Raw counts are length-biased: a 10 Mbp contig accumulates more hits than a 100 kbp contig simply by being larger. Dividing by (len / 1e6) puts all contigs on a comparable scale regardless of size.
df_grad <- df_plot %>%
  mutate(
    hym_per_mb  = n_hym  / (len / 1e6),
    arth_per_mb = n_arth / (len / 1e6),
    bact_per_mb = n_bact / (len / 1e6)
  )

# Produces a gradient blobplot where point color encodes BUSCO density (hits per Mb) rather than coverage class. High-density contigs appear in deep steelblue; zero-density contigs appear in grey.

# If a lineage produced no hits at all, falls back to a plain grey blobplot with "(no hits)" in the title to avoid a meaningless all-grey gradient.
make_gradient_blob <- function(data, density_col, lineage_label, outfile) {

  has_hits <- any(data[[density_col]] > 0, na.rm = TRUE)

  if (!has_hits) {
    message("[INFO] No ", lineage_label, " BUSCO hits — producing baseline blobplot")

    p <- ggplot(data, aes(x = gc, y = mean_cov, size = log10(len))) +
      geom_point(colour = "grey80", alpha = 0.7) +
      scale_y_log10() +
      scale_size(range = c(1, 8)) +
      labs(
        title = paste0(lineage_label, " BUSCO density: ", SPECIES, "  (no hits)"),
        x     = "GC (%)",
        y     = "Mean coverage (log10)",
        size  = "log10(contig length)"
      ) +
      theme_bw()

  } else {

    # Sort ascending by density so high-density contigs are drawn last. They appear on top of low-density grey points rather than being buried
    data <- data %>% arrange(.data[[density_col]])

    # Cap the color scale at the 99th percentile of non-zero values. Some extremely high-density contigs would otherwise compressthe gradient so that all other contigs appear uniformly grey.
    # pmin(value, cap) clips anything above the cap to the cap itself, so the top 1% all render at maximum color intensity without distorting the scale for the remaining 99%.
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
  }

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

# 7) Write all summary tables

message("[INFO] Writing summary tables...")

# Gene-level: every BUSCO hit with contig, lineage, status, description
write_tsv(busco,
  file.path(OUTDIR, "busco_gene_level.tsv"))

# Contig-level: all contigs merged with coverage and BUSCO counts/annotations
write_tsv(df,
  file.path(OUTDIR, "busco_per_contig_summary.tsv"))

# Bacteria-positive contigs — input for cobiont analysis downstream
write_tsv(df %>% filter(has_bact),
  file.path(OUTDIR, "busco_bacteria_contigs.tsv"))

# Mixed contigs: both Hymenoptera AND bacteria BUSCO hits on the same contig.
# Either assembly artifacts (two organisms fused into one contig), horizontal gene transfer regions, or bacterial insertions within host sequence. Warrant manual inspection before filtering.
write_tsv(df %>% filter(has_hym & has_bact),
  file.path(OUTDIR, "busco_mixed_hym_bact.tsv"))

# Anchor tier assignments (strong / supportive / none + bacteria flag)
write_tsv(anchors,
  file.path(OUTDIR, "busco_anchor_tiers.tsv"))

# 2D contingency: anchor type x coverage class (with within-group percentages)
write_tsv(distribution_host,
  file.path(OUTDIR, "busco_anchor_coverage_distribution.tsv"))

# 2D contingency: bacteria signal x coverage class
write_tsv(distribution_bact,
  file.path(OUTDIR, "busco_bacteria_coverage_distribution.tsv"))

# Full evaluation table: all contigs with anchor type + bacteria signal merged in
write_tsv(anchor_eval,
  file.path(OUTDIR, "busco_anchor_full_eval.tsv"))

# High-level counts of how many contigs carry each lineage's BUSCO signal and how many carry both Hymenoptera and bacteria (chimera candidates)
blobplot_summary <- df_plot %>%
  summarise(
    n_contigs = n(),
    hym = sum(has_hym),
    arth = sum(has_arth),
    bact = sum(has_bact),
    mixed_hym_bact = sum(has_hym & has_bact)
  )

write_tsv(blobplot_summary,
  file.path(OUTDIR, "blobplot_summary.tsv"))

cat("\n=== Blobplot summary ===\n")
print(blobplot_summary)

cat("\n=== BUSCO density summary ===\n")
cat(sprintf("Contigs with hym BUSCO:  %d\n",  sum(df$n_hym > 0)))
cat(sprintf("Contigs with arth BUSCO: %d\n",  sum(df$n_arth > 0)))
cat(sprintf("Contigs with bact BUSCO: %d\n",  sum(df$n_bact > 0)))
cat(sprintf("Contigs with no BUSCO:   %d\n",
  sum(df$n_hym == 0 & df$n_arth == 0 & df$n_bact == 0)))
cat("=============================\n")

message("[OK] BUSCO analysis complete. Output: ", OUTDIR)