#!/usr/bin/env Rscript

# Stacked bars with supergroup composition: for each supergroup, the number of nuwt fragments that were assigned
# Usage:
# Rscript plot_supergroup_composition annotated_tsv out_prefix

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

ANNOTATION_TSV  <- args[1]
OUT_PREFIX <- args[2]

# Read annotated file
annotation <- read_tsv(ANNOTATION_TSV, quote = "", show_col_types = FALSE)
# Count the fragments for each supergroup by orthogroup type
annotation <- filter(annotation, supergroup != "?")
aggregation <- count(annotation, supergroup, og_type, name = "n")

# Order supergroups by total size with the largest first
sg_levels <- aggregation %>% 
  group_by(supergroup) %>%
  summarise(tot = sum(n), .groups = "drop") %>%
  arrange(desc(tot)) %>% 
  # Take the supergroups as character vectors
  pull(supergroup)

# Colours and legends
type_levels <- c("other Wolbachia gene", "no annotation", "bacterial IS", "mitochondrial (residual)", "mitochondrial (removed)", "ankyrin")
pal <- c("other Wolbachia gene" = "grey72", "no annotation" = "grey88", "bacterial IS" = "#3A8DDE", "mitochondrial (residual)" = "#9AD9B5", "mitochondrial (removed)" = "#2FB573", "ankyrin" = "#e65d38ff")
aggregation <- aggregation %>%
  mutate(supergroup = factor(supergroup, levels = sg_levels), og_type = factor(og_type, levels = type_levels))

# Plot 
supergroup_composition <- ggplot(aggregation, aes(x = supergroup, y = n, fill = og_type)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = pal, name = "Orthogroup type", drop = FALSE) +
  labs(x = "Wolbachia supergroup", y = "NUWT fragment assignments") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), legend.position = "right")

ggsave(paste0(OUT_PREFIX, ".pdf"), supergroup_composition, width = 9, height = 5.5, units = "in")

# Print the numbers
wide <- aggregation %>% 
  pivot_wider(names_from = og_type, values_from = n, values_fill = 0)
write_tsv(wide, paste0(OUT_PREFIX, "_counts.tsv"))