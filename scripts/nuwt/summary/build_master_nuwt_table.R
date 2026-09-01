#!/usr/bin/env Rscript

# Build a master table with one row for each species, supergroup and orthogroup, with th fragment count, the total aligned bp and the swissprot class and product name
# Usage:
# Rscript build_master_nuwt_table.R annotated.tsv out_prfix

args <- commandArgs(trailingOnly = TRUE)

NUWT <- args[1]
OUT_PREFIX  <- args[2]

# Function to read files
read_tsv0 <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)

# Count the fragments and sum the base pairs for each species, supergroup and orthogroup
data <- read_tsv0(NUWT)
all_groups <- list(species = data$species, supergroup = data$supergroup, OG_id = data$og)
n_frag <- aggregate(list(n_fragments = data$length), all_groups, length)
total_bp <- aggregate(list(total_bp = data$length), all_groups, sum)
master <- merge(n_frag, total_bp, by = c("species", "supergroup", "OG_id"))
# Add the orthogroup annotation across the table
# og_class, product and og_type are the same for every fragment that has that OG
og_annot <- unique(data[, c("og", "og_class", "product", "og_type")])
names(og_annot)[1] <- "OG_id"
master <- merge(master, og_annot, by = "OG_id")
# Reorder rows and columns
master <- master[order(master$species, master$supergroup, -master$n_fragments), ]
master <- master[ , c("species", "supergroup", "OG_id", "og_class", "product", "og_type", "n_fragments", "total_bp")]
write.table(master, paste0(OUT_PREFIX, "_master.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)