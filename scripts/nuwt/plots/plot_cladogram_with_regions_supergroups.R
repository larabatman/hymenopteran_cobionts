#!/usr/bin/env Rscript


# This script builds the main figure of the host taxonomy with one ring on top against the tips, showing the per-specie NUWT counts stacked by supergroup
# Supergorups L and F are also coloured
#
# Input files:
# 1. tree.nwk of the host taxonomy tree from build_host_tree.py.
# 2. _region_counts.tsv from the bedtools merge: long with species sg_vec Freq which is the count. It contains one row per species times supergroup
#
# Outputs
# <out_prefix>_composite.{pdf,png} 
# <out_prefix>_composite_plotdata.tsv the summary with species, supergroup, value
#
# Usage
# module load R-bundle-Bioconductor/3.15-foss-2021a-R-4.2.1
# Rscript plot_fig2_composite.R <tree_nwk> <values_tsv> <metadata_tsv> <out_prefix>

library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)


args <- commandArgs(trailingOnly = TRUE)

tree_nwk <- args[1]
per_nuwt_f <- args[2]
out_prefix <- args[3]
# Plotting parameters
fig_size <- 24  # the whole canvas
open_angle  <- 8 # gap in the fan,
branch_width <- 0.3
label_size <- 1.5 # tip labels that od not grow with the figure
label_gap <- 0.02 # gap from the tips to the labels, as a ratio
ring_offset <- 0.16 # gap from the tree to the ring, as a ratio
ring_width <- 0.30 # thickness of the ring, as a ratio of tree width

# Colour the supergroups of interest
# As epcies is coloured it at least one NUWT region was present and assigned to that supergroup
carrier_supergroups <- c("L", "F")
carrier_alpha <- 0.30

# Read the species tree
tr <- read.tree(tree_nwk)

# Read the values: regions values with species sg_vec Freq
inp <- read.table(per_nuwt_f, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
r1 <- data.frame(species = inp$species,supergroup = inp$sg_vec, value = as.numeric(inp$Freq), stringsAsFactors = FALSE)
r1 <- r1[!is.na(r1$value) & r1$value > 0, ] # drop NA and 0s

# Tips with no row at all render as an empty spoke
n_empty <- length(setdiff(tr$tip.label, unique(r1$species)))

# Keeps the raw values and supergroup as plain text for the plotdata file and for the carrier lookup, before supergroup becomes a factor
r1_raw <- r1

# Palette
palette_sg <- c(A = "#0072B2", B = "#E69F00", C = "#009E73", D = "#D55E00",
                F = "#CC79A7", J = "#56B4E9",
                E = "#8C6D31", G = "#7A8DDB", I = "#B79F00", K = "#9AA0A6",
                L = "#C4A29E", M = "#6FB59A", P = "#B888C4", S = "#9DB07C",
                T = "#C98A6A", V = "#7FB0C4", W = "#B0A57C",
                "?" = "#BDBDBD", unknown = "#DDDDDD")
sg_levels <- c("A","B","C","D","E","F","G","I","J","K","L","M","P","S","T","V","W","?","unknown")
all_sg <- intersect(sg_levels, unique(r1$supergroup))
r1$supergroup <- factor(r1$supergroup, levels = all_sg)
cols <- palette_sg[all_sg]

# Define which species carries L or F:
carrier_nodes <- list()
for (sg in carrier_supergroups) {
  spp <- sort(unique(r1_raw$species[r1_raw$supergroup == sg])) # subset the species column to rows where L or F are found
  carrier_nodes[[sg]] <- which(tr$tip.label %in% spp) # this converts the tip.label into positions
}
# Draw
# This plot contains four layers of drawinf:
# 0: the tree, which is a ggtree layout "fan" object that bends the tree into a circle with open.angle giving the degree gap to ensure the ends don't collide
p <- ggtree(tr, layout = "fan", open.angle = open_angle, size = branch_width)
# p$data holds the x and y position for every node 
tree_w <- max(p$data$x, na.rm = TRUE) 
label_offset <- tree_w * label_gap # label gap must be converted because geom_tiplab's offset is absolute, while geom_fruit is a ratio

# Layer 1: the coloured wedges for each supergroup of interest
# The wedges are added first so they sit under the tree.
# Draw the carrier wedges, one at a time each in its own  colour
for (sg in carrier_supergroups) {
  fill_sg <- palette_sg[[sg]]
  for (n in carrier_nodes[[sg]]) {
    p <- p + geom_hilight(node = n, fill = fill_sg, alpha = carrier_alpha)
  }
}
# Layer 2: the tip labels
p <- p +
  geom_tiplab(size = label_size, offset = label_offset, align = TRUE, linesize = 0.05) +
  # the ring: counts, outside the labels, with an axis because bar length
  # carries a real quantity
  # Layer three: the ring with ggExtra mechanism of attaching a plot to a tree which takes a separate df and joins it to the tree by matching y = species against the tip labels
  geom_fruit(data = r1, geom = geom_col,
             mapping = aes(y = species, x = value, fill = supergroup),
             orientation = "y", pwidth = ring_width, offset = ring_offset,
             axis.params = list(axis = "x", text.size = 1.3, nbreak = 3))
# Orientation = "y" is what draws horizontal stacked bars, one per species, segments coloured by supergroup and sized by valze
p <- p +
  scale_fill_manual(values = cols, drop = FALSE, name = "supergroup") + theme(legend.position = "right")

ggsave(paste0(out_prefix, "_composite.pdf"), p,
       width = fig_size, height = fig_size, limitsize = FALSE)

# Plotted data, for reproducibility
pd <- r1_raw[, c("species", "supergroup", "value")]
write.table(pd, paste0(out_prefix, "_composite_plotdata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)