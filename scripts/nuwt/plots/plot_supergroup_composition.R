#!/usr/bin/env Rscript

# Stacked bar chart: for each Wolbachia supergroup, the number of NUWT fragment assignments, coloured by which kind of orthogroup the fragments came from.
#
# Input
# _annotated.tsv from 5_g_mito_ank.R: one row per assigned fragment, with species, supergroup, og, og_type. It reads the annotated table, not the one retained for merging regions, so the excluded orthogroups still appear as their own bar segment.
#
# Output
# <out.pdf> the distribution  figure
# <out>_counts.tsv the table
#
# USAGE
# Rscript plot_supergroup_composition.R <annotated_tsv> <out.pdf>

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
annot_f  <- args[1]
out_file <- args[2]

# Fragments the tree could not resolve to a single supergroup ("?") are dropped from the figure
DROP_UNRESOLVED <- TRUE

# Read annotated file
d <- read_tsv(annot_f, quote = "", show_col_types = FALSE)
# Count fragments per supergroup times orthogroup type
agg <- d %>%
  { if (DROP_UNRESOLVED) filter(., supergroup != "?") else . } %>%
  count(supergroup, og_type, name = "n")

# Supergroups ordered by total size, largest on the left
sg_levels <- agg %>% group_by(supergroup) %>%
  summarise(tot = sum(n), .groups = "drop") %>%
  arrange(desc(tot)) %>% pull(supergroup)

# Fixed stacking, legend order and colours
type_levels <- c("other Wolbachia gene", "no annotation", "bacterial IS", "mitochondrial (residual)", "mitochondrial (removed)", "ankyrin")

pal <- c("other Wolbachia gene" = "grey72",
         "no annotation" = "grey88",
         "bacterial IS" = "#3A8DDE",
         "mitochondrial (residual)" = "#9AD9B5",
         "mitochondrial (removed)" = "#2FB573",
         "ankyrin" = "#e65d38ff")

agg <- agg %>%
  mutate(supergroup = factor(supergroup, levels = sg_levels), og_type = factor(og_type,    levels = type_levels))

p <- ggplot(agg, aes(x = supergroup, y = n, fill = og_type)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = pal, name = "Orthogroup type", drop = FALSE) +
  labs(x = "Wolbachia supergroup", y = "NUWT fragment assignments") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), legend.position    = "right")

ggsave(out_file, p, width = 9, height = 5.5, units = "in")

# Keep a record of the numbers
wide <- agg %>% pivot_wider(names_from = og_type, values_from = n, values_fill = 0)
counts_file <- sub("\\.pdf$", "_counts.tsv", out_file)
write_tsv(wide, counts_file)