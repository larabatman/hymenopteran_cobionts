#!/usr/bin/env Rscript

# Joins the OG swissprot table and Prokka gene annotations to every hit in every species filtered dfam table. 
# Produces per-species annotated hit tables
#
# The annotation layers:
# og_class: classification from the swissprot table:
# eukaryotic_TE
# bacterial_IS
# wolbachia_gene
# no_swissprot_hit
# prokka_name: consensus gene product name from Prokka annotations across all member proteins of this OG
# sp_description: top SwissProt hit description (from blastp) "no_hit" if no SwissProt match
#
# Usage:
#   Rscript annotate_nuwt_hits.R
#
# Outputs per species, in results/<Species>/filter/:
#   <Species>_nuwt_hits_annotated.tsv
#    Tab-separated, one row per hit, all dfam columns plus annotation columns


library(dplyr)
library(tidyr)
library(readr)
library(stringr)


PROJECT_ROOT <- "/data/projects/p2025-0083_mining_cobionts"
NUWT_DIR <- file.path(PROJECT_ROOT, "nuwt_scan")
RESULTS_DIR <- file.path(NUWT_DIR, "results")
SWISSPROT_TSV <- file.path(NUWT_DIR, "hmm_database", "og_lookup_table.tsv")
PROKKA_TSV <- file.path(NUWT_DIR, "hmm_database", "og_gene_annotations.tsv")
SUMMARY_OUT <- file.path(RESULTS_DIR, "annotation_summary.tsv")

# Read annotation tables
swissprot <- read_tsv(SWISSPROT_TSV, show_col_types = FALSE) %>%
  select(og_id, og_class, sp_description)

prokka <- read_tsv(PROKKA_TSV, show_col_types = FALSE) |>
  select(og_id, prokka_name = consensus_name, prokka_is_te = is_te_by_name)

# Join swissprot and Prokka into a single annotation table
# This is the reference we join against every species hit table
og_annot <- swissprot %>%
  left_join(prokka, by = "og_id")

# dfam columns
# dfam table is space-delimited with these columns (verified from actual output)
# We name all columns so the output TSV is self-documenting
dfam_col_names <- c( "og_id", "acc", "scaffold", "score", "evalue", "bias", "hmm_from", "hmm_to", "strand", "ali_from", "ali_to", "env_from", "env_to", "modlen", "description")

# Parse one dfam table. 
parse_dfam <- function(filepath) {
  lines      <- readLines(filepath)
  data_lines <- lines[!startsWith(lines, "#")]

  if (length(data_lines) == 0) return(NULL) # filter removed every hit

  df <- read.table(
    text = data_lines,
    header = FALSE,
    quote = "", #quote="" and comment="": the description field is free text and could contain quote or "#" characters that would otherwise corrupt the parse.
    comment = "",
    col.names = dfam_col_names,
    colClasses = c(og_id = "character", scaffold = "character",
                   score = "numeric", evalue = "numeric", bias = "numeric",
                   strand = "character", ali_from = "numeric", ali_to = "numeric",
                   env_from = "numeric", env_to = "numeric")
  )

  if (nrow(df) == 0) return(NULL)

  # Envelope coordinates, normalised: on a minus-strand hit env_from > env_to,
  # so pmin/pmax give a consistent start < end and a positive length.
  df %>%
    mutate(hit_start = pmin(env_from, env_to), hit_end = pmax(env_from, env_to), hit_len = hit_end - hit_start + 1)
}

# For each species, use the contamination-filtered dfam tables that contain all hits after removing assembly contamination, before any
dfam_files <- Sys.glob(file.path(RESULTS_DIR, "*", "filter", "*_nuwt_hits_dfam_filtered.tbl"))

summary_rows <- list()
n_done <- 0
n_skipped <- 0

for (dfam_file in dfam_files) {
  # Extract species name from path
  species <- basename(dirname(dirname(dfam_file)))
  filter_dir <- dirname(dfam_file)
  out_file <- file.path(filter_dir, paste0(species, "_nuwt_hits_annotated.tsv"))

  hits <- parse_dfam(dfam_file)
  # Join annotation columns to every hit
  # left_join: all hits retained, NAs for OGs not in swissprot table
  annotated <- hits %>%
    left_join(og_annot, by = "og_id") %>%
    mutate(
      # Fill NAs for OGs not in the swissprot table
      og_class = replace_na(og_class, "no_swissprot_hit"),
      sp_description = replace_na(sp_description, "no_hit"),
      prokka_name = replace_na(prokka_name, "hypothetical protein"),
      prokka_is_te = replace_na(prokka_is_te,    FALSE),

      # Combined classification that merges SwissProt and Prokka evidence:
      # If Prokka directly annotated this OG as a transposase (prokka_is_te) but SwissProt classified it as bacterial_IS, keep bacterial_IS
      # If Prokka says transposase but SwissProt said wolbachia_gene or  no_swissprot_hit, upgrade to bacterial_IS
      og_class_final = case_when(
        og_class == "eukaryotic_TE" ~ "eukaryotic_TE",
        og_class == "bacterial_IS" ~ "bacterial_IS",
        prokka_is_te & og_class == "wolbachia_gene" ~ "bacterial_IS",
        prokka_is_te & og_class == "no_swissprot_hit" ~ "bacterial_IS",
        og_class == "wolbachia_gene" ~ "wolbachia_gene",
        TRUE ~ "no_annotation"
      ),

      # Include a species column
      species = species
    ) %>%
    # Reorder columns: species and annotation first, then hit coordinates, then full dfam fields for reference
    select(species, og_id, og_class_final, prokka_name, sp_description, scaffold, strand, hit_start, hit_end, hit_len, score, evalue, bias, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, modlen)

  # Write annotated table
  write_tsv(annotated, out_file)

