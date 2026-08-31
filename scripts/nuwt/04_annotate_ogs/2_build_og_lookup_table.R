#!/usr/bin/env Rscript

# Annotation table for the representative proteins of each orthogroup
# Output format columns -outfmt
# qseqid is the OG identifier, as OG0001392: og_id
# sseqid is the SwissProt accession: sp_hit
# pident is the identity: sp_pident
# length is the alignment length: aln_length
# evalue is the Evalue of the match: sp_evalue
# bitscore is the bit score: sp_bitscore
# stitle is the SwissProt entry title: sp_description
# Usage:
# Rscript 2_build_og_lookup_table.R representative_fasta blastp_results out_tsv

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

REP_FASTA <- args[1]
BLASTP <- args[2]
OUT_TSV <- args[3]

# Retrieve all the OG ids from the representative fasta file, since many of them did not get a blast hit and thus are not in the results table
# Read the whole fasta as a character vector, keep header lines, strips > to retrieve the IDs
fasta_lines <- readLines(REP_FASTA)
all_ogs <- fasta_lines[startsWith(fasta_lines, ">")] %>%
  str_remove("^>") %>% str_trim()

# Read blastp results: qseqid sseqid pident length evalue bitscore stitle
sp_col_names <- c("og_id", "sp_hit", "sp_pident", "aln_length","sp_evalue", "sp_bitscore", "sp_description")
sp_raw <- read.table(BLASTP, sep = "\t", header = FALSE, col.names = sp_col_names, fill = TRUE, quote = "", comment   = "", stringsAsFactors = FALSE)

# Select the best hit for one OG on the bitscore 
# Sort, slice and select top 1
sp_hits <- sp_raw %>%
    mutate(sp_evalue = as.numeric(sp_evalue), sp_pident = as.numeric(sp_pident), sp_bitscore = as.numeric(sp_bitscore)) %>%
    group_by(og_id) %>%
    arrange(desc(sp_bitscore), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(og_id, sp_hit, sp_description, sp_evalue, sp_pident, sp_bitscore)

# Classification:
# eukaryotic_TE when a description matches TE keywords, and source organism is eukaryote
# bacterial_IS when a description matches TE keywords, and source organism is bacteria
# wolbachia gene when a description does not match a TE keyword
# no_swissprot_hit when no blastp hit 

# Build a TE regex:
te_keywords <- regex(paste("transpos", "reverse.transcriptase", "gag.pol", "pol.polyprotein", "retrovir", "\\bLTR\\b", "\\bSINE\\b", "\\bLINE\\b", "\\bTE\\b", "Ty[0-9]", "\\bcopia\\b", "\\bgypsy\\b", "mariner", "helitron", "piggyBac", "\\bIS[0-9]\\b", sep = "|"), ignore_case = TRUE)

# Build bacteria regex:
bacterial_organisms <- regex(paste("escherichia", "coxiella", "brucella", "sinorhizobium", "deinococcus", "staphylococcus", "burkholderia", "microchaete", "bacillus", "rhizobium", "mesorhizobium", "agrobacterium", "caulobacter", "rickettsia", "wolbachia", "alphaproteobac", "bacteriophage", "phage", sep = "|"), ignore_case = TRUE)

# Join all OG ids and classify
lookup <- tibble(og_id = all_ogs) %>%
  left_join(sp_hits, by = "og_id") %>% 
  mutate(
    # Define the source organism: take anything between OS and OX
    sp_organism = str_match(sp_description, "OS=(.+?)\\s+OX=")[, 2],
    # Match to TE regex
    sp_keyword_match = !is.na(sp_description) & str_detect(sp_description, te_keywords),
    # Match to bacterial regex
    is_bacterial_source = !is.na(sp_organism) & str_detect(sp_organism, bacterial_organisms),
    # Assign classes
    og_class = case_when(
      # No SwissProt hit
      is.na(sp_description) ~ "no_swissprot_hit",
      # Match to TE, source match bacterial
      sp_keyword_match & is_bacterial_source ~ "bacterial_IS", 
      # Match to TE, source did not match bacteria
      sp_keyword_match & !is_bacterial_source ~ "eukaryotic_TE", 
      # No match to TE
      !sp_keyword_match ~ "wolbachia_gene"),
    # NA are replaced by no_hit all round
    sp_description = replace_na(sp_description, "no_hit"),
    sp_organism = replace_na(sp_organism, "no_hit"),
    sp_hit = replace_na(sp_hit,"no_hit")) %>%
  select(og_id, sp_hit, sp_description, sp_organism, sp_evalue, sp_pident, sp_bitscore, og_class) %>%
  arrange(og_id)

# Write output
write_tsv(lookup, OUT_TSV)