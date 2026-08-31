#!/usr/bin/env Rscript

# For each fragment that got a supergroup assigned, adds another layer of annotation
# Usage:
# Rscript 4_og_mito_ank.R per_nuwt.tsv og_lookup_table.tsv removed_files_dir out_prefix
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

NUWTS <- args[1]
LOOKUP_TABLE <- args[2]
FILES_DIR <- args[3]
OUT_PREFIX <- args[4]

# Function to parse tsv tables
read_tsv0 <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
# Function to read the removed and residula mito ids
read_ids <- function(f) {
  x <- trimws(readLines(file.path(FILES_DIR, f), warn = FALSE))
  x[nzchar(x)]
}
removed_mito <- read_ids("removed_mito.txt")
residual_mito <- read_ids("residual_mito.txt")

# Read files 
data <- read_tsv0(NUWTS)
lookup <- read_tsv0(LOOKUP_TABLE)[, c("og_id", "og_class", "sp_description")]
names(lookup) <- c("og", "og_class", "product")
# Merge the nuwt fragments and the annotation
data <- merge(data, lookup, by = "og", all.x = TRUE)
data$og_class[is.na(data$og_class)] <- "no_swissprot_hit"
data$product[is.na(data$product)] <- "no_hit"

# Classify the orthogroups
UNINFROMATIVE <- "^no_hit$|hypothetical|uncharacteri[sz]ed|unkown function|^$"
has_product <- !grepl(UNINFROMATIVE, data$product, ignore.case = TRUE)
# The first match wins
data$og_type <- case_when(
  data$og %in% removed_mito ~ "mitochondrial (removed)",
  data$og %in% residual_mito ~ "mitochondrial (residual)",
  grepl("ankyrin", data$product, ignore.case = TRUE) ~ "ankyrin",
  data$og_class == "bacterial_IS" ~ "bacterial IS",
  !has_product ~ "no annotation",
  TRUE ~ "other Wolbachia gene")
# The removed mitochondrial OG
data$excluded <- data$og %in% removed_mito
# Before and after
cols <- c("name", "species", "supergroup", "og", "length", "og_class", "product", "og_type", "excluded")
data <- data[order(data$species, data$supergroup, data$og), cols]

write.table(data, paste0(OUT_PREFIX, "_annotated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data[!data$excluded, ], paste0(OUT_PREFIX, "_kept.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)