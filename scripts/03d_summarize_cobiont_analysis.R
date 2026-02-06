#!/usr/bin/env Rscript

# This script integrates window-level Kraken classification and contig level GC content and HiFi coverage to profuce a contig-level taxonomy/ GC/ coverage table, a bp-weighted summary by contig class and a GC vs coverage QC plot

# Load libraries
library(dplyr) # data manipulation
library(readr) # TSV reeading and writing
library(stringr) # string manipulation
library(ggplot2) # plottinh


# Input paths
# TSV produced by the windowed Kraken script, no header in file: window_id  taxid status (C or U)
window_file <- "results/Lasioglossum_pauxillum_kraken_windows/window_taxid_status.tsv"
# TSV produced by GC and coverage script, columns: contig len gc  mean_cov
gc_cov_file <- "results/Lasioglossum_pauxillum_gc_cov/gc_cov.tsv"
# Prefix for output files 
out_prefix <- "lasioglossum_generic"
# Output path
out_dir <- "results/Lasioglossum_pauxillum_qc_summary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters controlling contig classification 
# Minimum number of windows required for a contig to be considered long enough for reliable classification
MIN_WINDOWS <- 10
# Minimum fraction of classified windows that must agree on the same taxid for a contig to be labeled as dominant_taxon
DOMINANCE_THRESHOLD <- 0.8
# Fraction of windows with taxid == 0 required to label a contig as mostly_unclassified
UNCLASSIFIED_THRESHOLD <- 0.95

# Helper functions
# Normalize contig name function
normalize_contig <- function(x){
  x %>%
    str_remove(":.*$") %>% # remove window coordinate suffix added by seqkit sliding
    str_remove("_sliding$") %>% # remove trialing "_sliding" if present
    str_trim() # remove leading/ trailing whitespace
}
# Funcion to build full output paths
out_path <- function(fname) file.path(out_dir, fname)

# Read window-level Kraken output
# No header in file, expected structure: window_id  taxid C|U
win <- read_tsv(
  window_file,
  col_names = c("window", "taxid", "status"),
  show_col_types = FALSE
) %>%
  mutate(
    contig = normalize_contig(window), # extract original contig name from window ID
    taxid = as.character(taxid) # ensure taxid is treated as character for comparison
  )

# Count windows per contig and taxid
# For each contig-taxid pair, n = number of windos assigned to that taxid
counts <- win %>%
  count(contig, taxid, name = "n")

# Collapse to contig-level taxonomic summary 
contig_tax <- counts %>%
  group_by(contig) %>%
  summarise(
    # Total number of windows in this contig
    n_windows = sum(n),
    
    # Windows with taxid == 0 are unclassified
    n_unclassified = sum(n[taxid == "0"]),
    frac_unclassified = n_unclassified / n_windows,

    # Windows assigned to any non-zero windows
    n_classified = sum(n[taxid != "0"]),

    # Dominant taxid among classified windows 
    top_taxid = if (any(taxid != "0")) {
      taxid[taxid != "0"][which.max(n[taxid != "0"])]
    } else {
      "0"
    },
    
    # Number of windows assigned to that dominant taxid 
    top_n = if (any(taxid != "0")) {
      max(n[taxid != "0"])
    } else {
      0
    },
    # Fraction of classified windows supporting the dominant taxid
    top_frac_classified =
      if (n_classified > 0) top_n / n_classified else 0,

    .groups = "drop"
  )

# Assign taxonomy labels
contig_tax <- contig_tax %>%
  mutate(
    label = case_when(
      # Long contigs with a clear dominant taxon
      n_windows >= MIN_WINDOWS &
        top_frac_classified >= DOMINANCE_THRESHOLD &
        top_taxid != "0" ~ "dominant_taxon",
      # Long contigs mostly unclassified by Kraken 
      n_windows >= MIN_WINDOWS &
        frac_unclassified >= UNCLASSIFIED_THRESHOLD ~ "mostly_unclassified",
      # Short contigs or mixed signals
      TRUE ~ "mixed_or_short"
    )
  )

# Read contig GC and coverage table
gcov <- read_tsv(
  gc_cov_file,
  col_types = cols(
    contig   = col_character(),
    len      = col_double(),
    gc       = col_double(),
    mean_cov = col_double()
  )
) %>%
  mutate(contig = normalize_contig(contig)) # Normalize contig names to ensure matching with Kraken-derived contigs

# Merge taxonomy with GC and coverage
contigs <- contig_tax %>%
  left_join(gcov, by = "contig") # Keep all contigs classified by Kraken, even if GC/ coverage is missing (NA)

# Write contig table
# Full contig-level table
contig_tsv <- out_path(paste0(out_prefix, "_contig_tax_gc_cov.tsv"))
write_tsv(contigs, contig_tsv)
# Basepair-weighted summary by contig class
bp_tsv     <- out_path(paste0(out_prefix, "_bp_summary.tsv"))
# Compute how much of the assembly, in bp, falls into each contig class
bp_summary <- contigs %>%
  group_by(label) %>%
  summarise(
    bp = sum(len, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(percent = 100 * bp / sum(bp))

write_tsv(bp_summary, bp_tsv)

# GC vs coverage scatter plot
plot_png   <- out_path(paste0(out_prefix, "_gc_cov_bloblike.png"))
# GC and log-scaled coverage
p <- ggplot(contigs, aes(x = gc, y = mean_cov)) +
  geom_point(
    aes(size = len, colour = label),
    alpha = 0.6
  ) +
  scale_y_log10() +
  scale_size(range = c(1, 10)) +
  scale_colour_manual(
    values = c(
      dominant_taxon = "firebrick",
      mostly_unclassified = "steelblue",
      mixed_or_short = "grey70"
    )
  ) +
  labs(
    x = "GC (fraction)",
    y = "Mean HiFi coverage (log scale)",
    colour = "Contig class",
    size = "Contig length"
  ) +
  theme_bw()

ggsave(
  filename = plot_png,
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)

# Done
message("QC finished. Outputs:")
message(" - ", contig_tsv)
message(" - ", bp_tsv)
message(" - ", plot_png)
