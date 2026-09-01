#!/usr/bin/env Rscript

# Figure with host taxonomy and a ring on top that shows the nuwt regions for each species stacked by supergroup
# Usage:
# module load R-bundle-Bioconductor/3.15-foss-2021a-R-4.2.1
# Rscript plot_cladogram_with_regions_supergroups.R

library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)

args <- commandArgs(trailingOnly = TRUE)

TREE_NWK <- args[1]
NUWT_REGIONS <- args[2]
OUT_PREFIX <- args[3]

# Figure parameters
fig_size <- 24
open_angle <- 8
branch_width <- 0.3
label_size <- 1.5
label_gap <- 0.02
ring_offset <- 0.16
ring_width <- 0.30

# Colour supergroups F and L 
carrier_supergroups <- c("L", "F")
carrier_alpha <- 0.30

# Read the species cladogram
cladogram <- read.tree(TREE_NWK)

# Read the region file
# species sg_vec  Freq
regions <- read.table(NUWT_REGIONS, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
regions_df <- data.frame(species = regions$species,supergroup = regions$sg_vec, value = as.numeric(regions$Freq), stringsAsFactors = FALSE)
# Drop NA and 0s
regions_df <- regions_df[!is.na(regions_df$value) & regions_df$value > 0, ]

# Tips that have no row: empty
n_empty <- length(setdiff(cladogram$tip.label, unique(regions_df$species)))

# Keep raw values for supergorups L and F annotation before factoring by supergroup
regions_df_raw <- regions_df

# Colour palette
# Palette
palette_sg <- c(A = "#0072B2", B = "#E69F00", C = "#009E73", D = "#D55E00", F = "#CC79A7", J = "#56B4E9", E = "#8C6D31", G = "#7A8DDB", I = "#B79F00", K = "#9AA0A6", L = "#C4A29E", M = "#6FB59A", P = "#B888C4", S = "#9DB07C", T = "#C98A6A", V = "#7FB0C4", W = "#B0A57C", "?" = "#BDBDBD", unknown = "#DDDDDD")
sg_levels <- c("A","B","C","D","E","F","G","I","J","K","L","M","P","S","T","V","W","?","unknown")
all_sg <- intersect(sg_levels, unique(regions_df$supergroup))
regions_df$supergroup <- factor(regions_df$supergroup, levels = all_sg)
cols <- palette_sg[all_sg]

# Define which species carry L and F
# Define which species carries L or F:
carrier_species <- list()
for (sg in carrier_supergroups) {
  species_L_F <- sort(unique(regions_df_raw$species[regions_df_raw$supergroup == sg]))
  # Convert tip.label into positions
  carrier_species[[sg]] <- which(cladogram$tip.label %in% species_L_F)
}

# Plot with four layers: the cladogram ggtree fan, the coloured wedges for supergroups L and F geom_highlight, the tip labels with species name geom_tiplab and the ring of stacked bars geom_fruit
plot_clad <- ggtree(cladogram, layout = "fan", open.angle = open_angle, size = branch_width)
# Convert label gap
tree_width <- max(plot_clad$data$x, na.rm = TRUE) 
label_offset <- tree_width * label_gap 
# Colour the wedges for supergroup L and F 
for (sg in carrier_supergroups) {
  fill_sg <- palette_sg[[sg]]
  for (n in carrier_species[[sg]]) {
    plot_clad <- plot_clad + geom_hilight(node = n, fill = fill_sg, alpha = carrier_alpha)
  }
}
# Add tip labels
plot_clad <- plot_clad +
  geom_tiplab(size = label_size, offset = label_offset, align = TRUE, linesize = 0.05) +
  # geom_fruit: add the ring with the regions 
  # orientation = "y" draws horizontal stacked bars for each species 
  geom_fruit(data = regions_df, geom = geom_col, mapping = aes(y = species, x = value, fill = supergroup), orientation = "y", pwidth = ring_width, offset = ring_offset, axis.params = list(axis = "x", text.size = 1.3, nbreak = 3))

plot_clad <- plot_clad +
  scale_fill_manual(values = cols, drop = FALSE, name = "supergroup") + theme(legend.position = "right")

ggsave(paste0(OUT_PREFIX, "_composite.pdf"), plot_clad, width = fig_size, height = fig_size, limitsize = FALSE)

plot_data <- regions_df_raw[, c("species", "supergroup", "value")]
write.table(plot_data, paste0(OUT_PREFIX, "_composite_plotdata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)