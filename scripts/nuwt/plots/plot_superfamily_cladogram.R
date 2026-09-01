#!/usr/bin/env Rscript

# Build a cladogram with the host superfamilies
# Usage
# Rscript plot_superfamily_cladogram.R tree.nwk metadata.tsv out_prefix

library(ape)
library(ggplot2)
library(ggtree)

args <- commandArgs(trailingOnly = TRUE)

TREE_NWK <- args[1]
METADATA_TSV <- args[2]
OUT_PREFIX <- args[3]

tree <- read.tree(TREE_NWK)
meta <- read.table(METADATA_TSV, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)

# Species to superfamily
superfamily <- setNames(meta$superfamily, meta$species)[tree$tip.label]

# Keep one tip for each superfamily and drop the rest
keep <- sapply(split(names(superfamily), superfamily), `[`, 1)
tree_superfamily <- keep.tip(tree, keep)
tree_superfamily$tip.label <- superfamily[tree_superfamily$tip.label]

# Plot as a rectangular cladogram with ggtree
cladogram <- ggtree(tree_superfamily, branch.length = "none", size = 0.4) +
    geom_tiplab(size = 5, offset = 0.2) +
    theme_void()

cladogram <- cladogram + xlim(0, max(cladogram$data$x)*2.0)
ggsave(paste0(out_prefix, ".pdf"), cladogram)