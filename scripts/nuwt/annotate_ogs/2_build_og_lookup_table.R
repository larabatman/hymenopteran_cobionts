#!/usr/bin/env Rscript
# build_og_lookup_table.R
#
# Builds an orthogroup annotation lookup table from blastp results against SwissProt. 
#
#
# Usage:
# Rscript build_og_lookup_table.R <og_reps.fna> <blast_sp.blast6> <output.tsv>
#
# Output columns:
# og_id: Orthogroup identifier (OG0001392)
# sp_hit: Top SwissProt accession (no_hit if no match)
# sp_description: Full SwissProt entry title (gene name + description)
# sp_organism: Source organism of the SwissProt hit
# sp_evalue: blastp E-value
# sp_pident: blastp percent identity
# sp_bitscore: blastp bit score
# og_class: One of four categories:
#           eukaryotic_TE: match insect/yeast/algal TEs
#           bacterial_IS: bacterial insertion sequence 
#           wolbachia_gene: Wolbachia or bacterial gene 
#           no_swissprot_hit: No SwissProt hit at e-value 1e-5 

library(dplyr)
library(tidyr) # for replace_na()
library(readr)
library(stringr)

# Find arguments
args <- commandArgs(trailingOnly = TRUE)

REP_FASTA <- args[1]
BLAST_SP  <- args[2]
OUT_TSV   <- args[3]

# Read OG IDs from the representative FASTA
# Every OG appear in the output table, even if it got no BLAST hit (annotated as "no_hit").
# Everything downstream is joined to this, which guarantees an OG with no BLAST still appears in the output 
fasta_lines <- readLines(REP_FASTA) # read the whole FASTA as a character vector, keep header lines, strips > to retrieve the IDs
all_ogs <- fasta_lines[startsWith(fasta_lines, ">")] %>%
  str_remove("^>") %>% str_trim()

# Read blastp results
# outfmt "6 qseqid sseqid pident length evalue bitscore stitle sscinames"
# read.table with sep="\t" and fill=TRUE handles extra tabs by padding with NA
# quote="" disables quote interpretation as SwissProt descriptions can contain single and double quotes that would otherwise confuse the parser.
# similarly for comment = ""
# Renaming the columns
sp_col_names <- c("og_id", "sp_hit", "sp_pident", "aln_length","sp_evalue", "sp_bitscore", "sp_description", "sp_organism")

sp_raw <- read.table(BLAST_SP, sep = "\t", header = FALSE, col.names = sp_col_names, fill = TRUE, quote = "", comment   = "", stringsAsFactors = FALSE)

# Select the best hit per OG by bitscore: sort by bitscore and select the top 1 with slice
sp_hits <- sp_raw %>%
    mutate(sp_evalue = as.numeric(sp_evalue), sp_pident   = as.numeric(sp_pident), sp_bitscore = as.numeric(sp_bitscore)) %>%
    group_by(og_id) %>% # split, sort, take first 
    arrange(desc(sp_bitscore), .by_group = TRUE) %>% # arrange descending by bitscore within each OG, then take the top row
    slice(1) %>%
    ungroup() %>%
    select(og_id, sp_hit, sp_description, sp_organism, sp_evalue, sp_pident, sp_bitscore)

# Classification 
# eukaryotic_TE: description matches TE keywords AND source organism is a eukaryote (insect, yeast, plant, worm, alga etc.)
# bacterial_IS: description matches TE keywords AND source organism is a bacterium (E. coli, Coxiella, Brucella etc.)
# wolbachia_gene : description does NOT match TE keywords
# no_swissprot_hit: no blastp hit at e-value 1e-5

# TE regex: Keywords that indicate a TE-related protein in the SwissProt description.
# Covers LTR retrotransposons (gypsy, copia, Ty), non-LTR (LINE, SINE), DNA transposons (mariner, Tc1, piggyBac, helitron), and IS elements.
te_keywords <- regex(
  paste(
    "transpos", # transposon, transposase, retrotransposon
    "reverse.transcriptase",
    "gag.pol", # gag-pol, gag/pol polyprotein
    "pol.polyprotein",
    "retrovir", # retrovirus, retroviral
    "\\bLTR\\b", # word-boundary
    "\\bSINE\\b",
    "\\bLINE\\b",
    "\\bTE\\b",
    "Ty[0-9]", # Ty1, Ty3 yeast/insect elements
    "\\bcopia\\b",
    "\\bgypsy\\b",
    "mariner",
    "helitron",
    "piggyBac",
    "\\bIS[0-9]\\b", # IS1, IS3, IS4 etc.
    sep = "|"
  ),
  ignore_case = TRUE
)

# Bacterial origin
# Use the OS= field from the SwissProt stitle, which is already in sp_description (e.g. "OS=Escherichia coli").
bacterial_organisms <- regex(
  paste(
    "escherichia", "coxiella", "brucella", "sinorhizobium",
    "deinococcus", "staphylococcus", "burkholderia", "microchaete",
    "bacillus", "rhizobium", "mesorhizobium", "agrobacterium",
    "caulobacter", "rickettsia", "wolbachia", "alphaproteobac",
    "bacteriophage", "phage",
    sep = "|"
  ),
  ignore_case = TRUE
)

# Join and classify
lookup <- tibble(og_id = all_ogs) %>%
  left_join(sp_hits, by = "og_id") %>% # Left join: all OGs appear in output regardless of BLAST result
  mutate(
    sp_keyword_match = !is.na(sp_description) & str_detect(sp_description, te_keywords), # Does the description contain a TE keyword?
    is_bacterial_source = !is.na(sp_description) & str_detect(sp_description, bacterial_organisms), # Is the source organism a bacterium?
    og_class = case_when( # class assignment
      is.na(sp_description) ~ "no_swissprot_hit", # No SwissProt hit at all
      sp_keyword_match & is_bacterial_source ~ "bacterial_IS", # TE keyword matched and source is a bacterium
      sp_keyword_match & !is_bacterial_source ~ "eukaryotic_TE",  # TE keyword matched and source is not a bacterium
      !sp_keyword_match ~ "wolbachia_gene"  # No TE keyword
    ),
    # Replace NA fields with explicit strings
    sp_description = replace_na(sp_description, "no_hit"),
    sp_organism = replace_na(sp_organism, "no_hit"),
    sp_hit = replace_na(sp_hit,"no_hit")) %>%
  select(og_id, sp_hit, sp_description, sp_organism, sp_evalue, sp_pident, sp_bitscore, og_class) %>%
  arrange(og_id)

# Output
write_tsv(lookup, OUT_TSV)