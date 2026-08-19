#!/usr/bin/env Rscript

# Stacked bar chart of the ten species with the most NUWT hits, coloured by annotation class.
#
# Here, wolbachia_gene and bacterial_IS are combined into one "NUWT signal" bar as both are genuine Wolbachia sequence and could legitimately have transferred. 
# eukaryotic_TE and no_annotation are kept separate
#
# Inputs
# annotation_summary.tsv as tab-separated, header row, one row per species with hit counts per og_class_final category
#
# Usage:
# Rscript plot_nuwt_top10.R

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(forcats)

RESULTS_DIR <- "/data/users/lland/cobionts/nuwt_scan/results"
SUMMARY_TSV <- file.path(RESULTS_DIR, "annotation_summary.tsv")
PLOT_DIR <- file.path(RESULTS_DIR, "plots")
OUTPUT_PDF <- file.path(PLOT_DIR, "nuwt_hits_top10_raw.pdf")

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# Read and select top 10 species
summary_tbl <- read_tsv(SUMMARY_TSV, show_col_types = FALSE)

top10 <- summary_tbl %>%
  mutate(nuwt_signal = wolbachia_gene + bacterial_IS) %>%
  slice_max(order_by = total_hits, n = 10) %>%
  select(species, total_hits, nuwt_signal, eukaryotic_TE, no_annotation)

# For stacking: need pivot longer to expand the segments
# ggplot stacks one bar segment per row, so the three count columns become three rows per species.
top10_long <- top10 %>%
  pivot_longer(
    cols = c(nuwt_signal, eukaryotic_TE, no_annotation),
    names_to = "hit_group",
    values_to = "n_hits"
  ) %>%
  mutate(species= fct_reorder(species, -total_hits), # Order species by total hits descending: tallest bar on the left. The minus sign reverses fct_reorder's default ascending order.
         hit_group = factor(hit_group, levels = c("nuwt_signal", "eukaryotic_TE", "no_annotation")))

# Plot
GROUP_COLOURS <- c(
  nuwt_signal   = "#4E79A7",
  eukaryotic_TE = "#E15759",
  no_annotation = "#BAB0AC"
)

GROUP_LABELS <- c(
  nuwt_signal   = "NUWT signal (Wolbachia gene + IS)",
  eukaryotic_TE = "Eukaryotic TE (noise)",
  no_annotation = "No annotation"
)

p <- ggplot(top10_long, aes(x = species, y = n_hits, fill = hit_group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL) +
  labs( title = "Top 10 species by total NUWT hit count", x = NULL, y = "Number of hits") +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11))

pdf(OUTPUT_PDF, width = 8, height = 6)
print(p)
invisible(dev.off())