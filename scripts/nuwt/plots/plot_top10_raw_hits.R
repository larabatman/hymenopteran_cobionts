#!/usr/bin/env Rscript

# Stacked bar of the top ten species with most raw nuwt hits coloured by annotation class
# Usage: 
# Rscript plot_top10_raw_hits.R annotation.tsv out_prefix

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

ANNOT_TSV <- args[1]
OUT_PREFIX <- args[2]

# Read annotation table
annot_table <- read_tsv(ANNOT_TSV, show_col_types = FALSE)
# Treat wolbachia_gene and bacterial_IS as nuwt signal together
# Select the top 10 species with most hits
top10 <- annot_table %>%
  mutate(nuwt_signal = wolbachia_gene + bacterial_IS) %>%
  slice_max(order_by = total_hits, n = 10) %>%
  select(species, total_hits, nuwt_signal, eukaryotic_TE, no_annotation)

# Stack the bars with pivot_longer to have three count columns for each orthogroup type
top10_long <- top10 %>%
  pivot_longer(cols = c(nuwt_signal, eukaryotic_TE, no_annotation), names_to = "hit_group", values_to = "n_hits") %>%
  # Reorder the levels descending so the largest comes first
  mutate(species = factor(species, levels = top10$species[order(-top10$total_hits)]),
         hit_group = factor(hit_group, levels = c("nuwt_signal", "eukaryotic_TE", "no_annotation")))

# Palette colour and labels
group_colours <- c(nuwt_signal = "#4E79A7",eukaryotic_TE = "#E15759", no_annotation = "#BAB0AC")
group_labels <- c(nuwt_signal = "NUWT signal (Wolbachia gene + IS)", eukaryotic_TE = "Eukaryotic TE (noise)",no_annotation = "No annotation")

top10_plot <- ggplot(top10_long, aes(x = species, y = n_hits, fill = hit_group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = group_colours, labels = group_labels, name = NULL) +
  labs(title = "Top 10 species by total NUWT hit count", x = NULL, y = "Number of hits") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 11))

ggsave(paste0(OUT_PREFIX, ".pdf"), top10_plot, width = 8, height = 6, units = "in")