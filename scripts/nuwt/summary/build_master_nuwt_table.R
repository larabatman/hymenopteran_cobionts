#!/usr/bin/env Rscript

# Builds the master table: one row per (species, supergroup, orthogroup), with the fragment count, total aligned bp, and the orthogroup'sSwissProt class and product name.
#
# This is the input to plot_supergroup_composition.R, and the source for the supergroup ranking, the ankyrin fractions, and F's composition.
#
# Input
# _annotated.tsv from annotate_nuwt_fragments.R with one row per assigned fragment, with species, superfroup, og, length, og_clas, product, og_type, excluded.
# Output
# tsv with species, supergroup, OG_id, og_class, product, n_fragments, total_bp
#
# USAGE
# Rscript build_master_nuwt_table.R <_annotated.tsv> <out_tsv>

args <- commandArgs(trailingOnly = TRUE)
per_nuwt <- args[1]
out_tsv  <- args[2]

read_tsv0 <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)

# Count fragments and sum bp per species, supergroup and orthogroup
d <- read_tsv0(per_nuwt)
key <- list(species = d$species, supergroup = d$supergroup, OG_id = d$og)
n_frag  <- aggregate(d$length, key, length)
total_bp <- aggregate(d$length, key, sum)
master <- data.frame(n_frag[, c("species", "supergroup", "OG_id")], n_fragments = n_frag$x, total_bp    = total_bp$x,  stringsAsFactors = FALSE)

# Carry the orthogroup annotation across: og_class, product and og_type are properties of the orthogroup, so they are the same on every fragment of an OG: one row per OG is enough, and the merge cannot duplicate rows.
og_annot <- unique(d[, c("og", "og_class", "product", "og_type")])
names(og_annot)[1] <- "OG_id"
master <- merge(master, og_annot, by = "OG_id")

master <- master[order(master$species, master$supergroup, -master$n_fragments), c("species", "supergroup", "OG_id", "og_class", "product", "og_type", "n_fragments", "total_bp")]

write.table(master, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)