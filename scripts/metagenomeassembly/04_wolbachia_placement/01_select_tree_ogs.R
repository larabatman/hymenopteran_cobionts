#!/usr/bin/env Rscript

# Find marker protein set for the Wolbachia strain tree by applying the single-copy in at least 75% of strains criteria on the OrthoFinder table
# In Orthogroups.tsv: column 1 is the OG id and every other column is one strain
# Each cell is a comma space list of one strain total proteins in one OG
# The input format: Orthogroup.tsv.gz has 1452 columns
# Column 1: OG id OG0000000, OG0000001, ...
# Column 2 and after: one for each genome, 1451
# Each cell contains the genome proteins for that OG. They are comma space separated if multiple, and the cell is empty if it had no member protein for that OG

# Here we have all the 1451 clustered genomes and we need to restrict the columns to those named in the dereplicated set: ref_strains.txt
# The outgroups are not in ref_strains.txt and they are excluded from the scoring. They will be used for placement only. 

# So for each orthogroup and for each kept strain:
# present = cell is non-empty
# single-copy = cell has exactly one protein aka no comma
# Orthogroup  Ace GCA_001752665 GCA_02283975  GCA_030441635
# OG000046  Ace_00025 ..._00412 ..._00390 ..._00355
# OG000047    ..._00511 ..._00390, ..._00402  
# 0G000048  Ace_00031 ..._0733    ..._00688

# OG000046 has three kept strains all one protein: the ratio is 3/3 = 1, the OG is kept
# OG000047 appears twice in GCA_02283975 and none in GCA_030441635: the ratio is 1/3 = 0.33 and dropped
# 0G000048 is not in GCA_02283975 : the ratio is 2/3 = 0.66 and dropped
# The output is tree_ogs.txt which will be the gene set for place_mags.sh and it will contain one OG ID for each line.
# Usage: 
# Rscript select_tree_ogs.R  Orthogroups.tsv.gz  ref_strains.txt  tree_ogs.txt

library(readr)

args <- commandArgs(trailingOnly = TRUE)

ORTHO_TSV <- args[1]
REF_STRAINS <- args[2]
OUTFILE <- args[3]
# An OG is kept if it is present in single copy in at least 75% of the strains
min_fraction <- 0.75 

# Read Orthogroup.tsv.gz: plain columns, the cells are protein id lists
ortho <- read_tsv(ORTHO_TSV, show_col_types = FALSE, progress = FALSE, col_types = cols(.default = col_character()))
# The first column: c("OG000000", "OG000001", "OG000002", ...) and its length is the number of OG
og_ids <- ortho[[1]]

# Keep only the reference strains:
ref_strains <- trimws(readLines(REF_STRAINS)) # readLines reads the file into a character vector, for each line and strips newlines: c("GCA_001752665", "GCA_022836975", "GCA_030441635", ...)
ref_strains <- ref_strains[ref_strains != ""] # produce a logical vector the same length as ref_strains which is TRUE for every accession and FALSE for any element that is an empty string, to subset keeps only the TRUE positions
present_cols <- intersect(ref_strains, names(ortho)) # this is what the selection will use to pull the dereplicated strain columns out of the 1451
strain_cols <- ortho[present_cols] # only the dereplicated strains
n_strains <- ncol(strain_cols) # count them to get the denominator

# Score single copy OGs: a cell is single copy if it is not NA and not empty, and it it contains no comma
# Function: test columns one at a time, and return TRUE or FALSE for a given OG
is_single <- function(col) {
  x <- ifelse(is.na(col), "", trimws(col)) #read_tsv turns empty cells into NA such that NA and "" become the same thing before testing
  x != "" & !grepl(",", x, fixed = TRUE) # conditions: non-empty aka the strain had the OG and comma free aka it has exactly one protein
}
# Apply for each column: iterate ocer strain_cols. Call is_single once for each strain to get all OGs at once, then test all rows of that column reutrning a logical vector as long as the number of OGs
single_copy_matrix <- vapply(strain_cols, is_single, logical(nrow(ortho))) # vapply our is_single function across every kept strain column and assemble the results into a matrix. the logical(nrow(tab)) is a tempplate that declares what each call must return aka a logical vector of that length
# single_mat is OG x kept strains, TRUE when that strain has one protein in that OG exactly
frac_single <- rowSums(single_copy_matrix) / n_strains #rowSums(single_copy_matrix): counts for each OG how many kept strains are single-copy by coercing the logical to 0/1 
keep <- frac_single >= min_frac # keep the OGs that appear in at least 0.75 of them 
writeLines(og_ids[keep], out_file)