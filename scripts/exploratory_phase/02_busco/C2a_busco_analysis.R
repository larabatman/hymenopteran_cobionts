#!/usr/bin/env Rscript

# This Rscript parses the BUSCO table outputs full_table.tsv to summarize the BUSCO content of each contig
# Usage:
# Rscript C2a_busco_contigs.R species busco_hym busco_arth busco_bact gc_cov outdir

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

SPECIES <- args[1]
BUSCO_HYM <- args[2]
BUSCO_ARTH <- args[3]
BUSCO_BACT <- args[4]
COVERAGE_TABLE <- args[5]
OUTDIR <- args[6]

# Function to parse BUSCO full_table.tsv: problem was that there were different numbers of # comment lines and Missing rows don't have the same fileds as Complete or other genes
read_busco <- function(path, lineage) {
  raw <- readLines(path)
  raw <- raw[!grepl("^#", raw)] # drop the comment and header lines
  if (length(raw) == 0) {
    message("No BUSCO hits for ", lineage) # added check because one of the species had no bacterial marker gene set
    return(tibble())
  }
  raw <- gsub("\\\\t", "\t", raw) # replace the literal "\t"
  split <- strsplit(raw, "\t") # split on tab: list of character vectors
  # columns to keep: the busco id, the status as in complete, fragmented, duplicated, missing, and the description of the BUSCO gene
  # extract by position: sapply(..., `[`, N) takes the Nth element from each vector in the list
  df <- tibble(
    busco_id = sapply(split, `[`, 1),
    status = sapply(split, `[`, 2),
    sequence = sapply(split, `[`, 3),
    description = sapply(split, function(x) { # the description column is not there in Missing rows. It is the tenth column
      if (length(x) >= 10) x[10] else NA_character_
    })
  )
  # Drop missing rows
  df %>%
    filter(!is.na(status), status != "Missing") %>%
    transmute(contig = str_trim(sequence), busco_id, status, lineage = .env$lineage, description)
}
# Put the three tables together
busco <- bind_rows(read_busco(BUSCO_HYM, "hymenoptera"), read_busco(BUSCO_ARTH, "arthropoda"), read_busco(BUSCO_BACT, "bacteria"))

# Aggregate for each contig
# ful_table.tsv has onw roe for each gene in a contig, but we need to ollapse to 1 to join with the coverage table 
busco_contig <- busco %>%
  group_by(contig) %>% # summarize for each contig
  summarise(
    n_busco = n(), # BUSCO hits across all lineages for that contig
    n_hym = sum(lineage == "hymenoptera"), # number of hymenoptera BUSCO
    n_arth = sum(lineage == "arthropoda"), # number of arthropoda BUSCO
    n_bact = sum(lineage == "bacteria"), # number of bacteria BUSCO
    has_hym = n_hym > 0,
    has_arth = n_arth > 0,
    has_bact = n_bact > 0,
    hym_complete = sum(lineage == "hymenoptera" & status == "Complete"), # count different status for hymenoptera
    hym_fragmented = sum(lineage == "hymenoptera" & status == "Fragmented"),
    hym_duplicated = sum(lineage == "hymenoptera" & status == "Duplicated"),
    arth_complete = sum(lineage == "arthropoda" & status == "Complete"), # count different status for arthropoda
    arth_fragmented = sum(lineage == "arthropoda" & status == "Fragmented"),
    arth_duplicated = sum(lineage == "arthropoda" & status == "Duplicated"),
    bact_complete = sum(lineage == "bacteria" & status == "Complete"), # count different status for bacteria
    bact_fragmented = sum(lineage == "bacteria" & status == "Fragmented"),
    bact_duplicated = sum(lineage == "bacteria" & status == "Duplicated"),
    busco_ids = paste(unique(busco_id), collapse = "; "), # get the different busco gene ids 
    busco_lineages = paste(unique(lineage), collapse = "; "), # get the different lineages that each contig got a set for
    busco_descriptions = paste(unique(na.omit(description)), collapse = " | "), # get all the descriptions
    hym_functions = paste(unique(description[lineage == "hymenoptera"]), collapse = " | "), # breakdown each of these descriptions in each lineages
    arth_functions = paste(unique(description[lineage == "arthropoda"]), collapse = " | "),
    bact_functions = paste(unique(description[lineage == "bacteria"]), collapse = " | "),
    .groups = "drop"
  )

# Read coverage classification table 
cov <- read_tsv(COVERAGE_TABLE, show_col_types = FALSE)
# Merge coverage and busco tables: left_join so that we keep every contig from the coverage table including those that got no BUSCO hits
# contigs that are not in busco_contigs: NA in their BUSCO columns fater the join and this is problematic for the plots
# Replacing the NAs by: 0 for NA counts, FALSE for the boolean 
# coalesce(x, replacement): returns wither x when the value is not NA, or the replacement when the value is NA. across() to apply the same to each column
df <- cov %>%
  left_join(busco_contig, by = "contig") %>%
  mutate(
    across(c(n_busco, n_hym, n_arth, n_bact), function(x) coalesce(x, 0L)),
    across(c(has_hym, has_arth, has_bact), function(x) coalesce(x, FALSE))
  )

# Tables
write_tsv(busco, file.path(OUTDIR, "busco_gene_level.tsv"))
write_tsv(df, file.path(OUTDIR, "busco_per_contig_summary.tsv"))

bacterial <- df %>%
  filter(has_bact)

write_tsv(bacterial, file.path(OUTDIR, "busco_bacteria_contigs.tsv"))

both_hym_bact <- df %>%
  filter(has_hym & has_bact)

write_tsv(both_hym_bact, file.path(OUTDIR, "busco_mixed_hym_bact.tsv"))