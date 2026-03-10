#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
BUSCO_TABLE <- args[2] 
COVERAGE_TABLE <- args[3]
VALIDATION_DIR <- args[4]

dir.create(VALIDATION_DIR, showWarnings = FALSE, recursive = TRUE)

#---------------------------------------------------------
# 1) Read BUSCO full_table.tsv
# full_table.tsv actually contains 3 lines with #, of which the headers. 
# The columns names are thus given again

busco <- read_tsv(
    BUSCO_TABLE, 
    comment = "#", 
    col_names = c(
        "busco_id", 
        "status", 
        "sequence",
        "gene_start",
        "gene_end",
        "strand",
        "score",
        "length",
        "orthodb_url",
        "description"), 
    show_col_types = FALSE
)

#---------------------------------------------------------
# 2) Define anchors: contigs that show complete conserved genes
# Extract unique contig names
anchors <- busco %>%
    filter(status == "Complete") %>%
    transmute(contig = str_trim(sequence)) %>%
    distinct()

#---------------------------------------------------------
# 3) Evaluate coverage of host anchored contigs 
# Read coverage summary file
coverage <- read_tsv(COVERAGE_TABLE, show_col_types = FALSE)

# Join files by contig, retrieving host contigs
anchor_eval <- coverage %>%
    inner_join(anchors, by = "contig")

# Assess coverage distribution of the anchored contigs
distribution <- anchor_eval %>%
    count(coverage_class) %>%
    mutate(percent = n/sum(n) * 100)

#---------------------------------------------------------
# 4) Write files
write_tsv(anchors, file.path(VALIDATION_DIR, "busco_anchor_contigs.tsv"))
write_tsv(distribution, file.path(VALIDATION_DIR, "busco_anchor_coverage_distribution.tsv"))
print(distribution)

message("BUSCO backbone validation complete.")