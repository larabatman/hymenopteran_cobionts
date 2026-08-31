#!/usr/bin/env Rscript

# Concatenate the augmented alignments for each marker into a supermatrix with one row for each strain 
# Keeps reference strains in the dereplicated list and keeps outgroups
# Fills a missing marker with gaps to keep the rectangular shape of the matrix
# Usage:
# Rscript 06_concat_supermatrix.R prot2strain.tsv aug_dir tree_ogs.clean.txt out_prefix keep_strains.txt

library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

# prot2strain.tsv has two columns: protein_id tab strain
PROT2STRAIN <- args[1]
# Directory with the MAFFT alignments
AUG_DIR <- args[2] 
# List of single-copy marker proteins
CLEAN_OG_FILE <- args[3]
PREFIX <- args[4]
# List of dereplicated strains
DEREP_STRAINS <- args[5] 
# Outgroups:
OUTGROUPS <- c("Ace", "Ama", "Ech", "Eru")

# Read prot2strain.tsv
# ACE_01009 Ace
# GCA_000008025_00412 GCA_000008025
# GCA_000008025_00733 GCA_000008025 
# and so on 
map_prot_id <- read.table(PROT2STRAIN, sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("protein_id", "strain"))
prot2strain <- setNames(map_prot_id$strain, map_prot_id$protein_id) # builds a name character vector where the values are the strains and the names are protien ids: prot2strain["GCA_00008025_00412"] returns "GCA_00008025"
# A protein id found in this map becomes its strains, and an id that was not stays as is: MAGs pass through Acromyrmex_octospinosus_metamdbg_magscot_1 and no such id exists in the reference map
ref_strains_all <- unique(map_prot_id$strain) # collapse the protein rows down to the 1451 distinct strain names

# Set the reference strains from the dereplicated list and the outgroups
keep_strains <- trimws(readLines(DEREP_STRAINS))
keep_strains <- keep_strains[keep_strains != ""]
keep_strains <- union(keep_strains, OUTGROUPS) # adds the four outgroups to the keep list as they are not in wolbachia_reference.txt 

# Function to keep MAGs that are not reference strains 
keep_taxon <- function(taxon) (!(taxon %in% ref_strains_all)) | (taxon %in% keep_strains)
# Build file path for the og alignements that are in the clean list
ogs <- trimws(readLines(CLEAN_OG_FILE))
ogs <- ogs[ogs != ""]
aln_files <- file.path(AUG_DIR, paste0(ogs, ".faa"))
missing <- aln_files[!file.exists(aln_files)]

# Read each marker alignment and convert it from rows labelled by protein to rows labelled by strain: each protein marker needs to be mapped to its strain 
# for aug/OG0000046.faa: readAAStringSet(f) reads the alignment FASTA: names(aln) gives the headers and as.character the sequences as plain strings
# ids <- sub("\\s.*$", "", names(aln)) takes everything before the whitespace so ids: GCA_000008025_00412 and some description become GCA_000008025_00412
# then: taxa ifelse as translation: for each id, if it's in prot2strain, replace it with the strain name 
# References becomes their genome, and MAGs pass untouched which makes two proteins from the same strain in different markers land in the same supermatrix row
# aug/OG0000046.faa: the keep list contains GCA_001752665 but not GCA_900000001
# names (aln) as.character(aln)
# GCA_001752665_00412 MXOWJAEHDUEWEN--
# GCA_001752665_00733 MXOWJAEIDUEWEN--
# GCA_900000001_00120 --OWJAEHDUEWEN--
# Ace_01009 -XOWJAEHDUEWENA-
# Acromyrmex_octospinosus_..._magscot_1 MXOWJAEIDUEWEN--
# row 1 and row 2 are proteins from the same strain
# after ids and taxa:
# ids taxa
# GCA_001752665_00412 GCA_001752665
# GCA_001752665_00733 GCA_001752665
# GCA_900000001_00120 GCA_900000001
# Ace_01009 Ace
# Acromyrmex_octospinosus_..._magscot_1 Acromyrmex_octospinosus_..._magscot_1
# The first two carry an identical label and MAG passes through because it was never in prot2strain
# The the filters: row 2 is dropped duplicates from being second in the file and row 3 is dropped from not being in the keep list
one_og <- list()
all_taxa <- character(0)
for (f in aln_files) {
  aln <- readAAStringSet(f)
  # Strip description before the first space
  ids <- sub("\\s.*$", "", names(aln))
   #if ids match exactly names in prot2strain, keep the prote2strain ID which is the GCA accession, otherwise keep the ids from here which are the MAG IDs. Vector with strain names
  taxa <- ifelse(ids %in% names(prot2strain), prot2strain[ids], ids)
  # !duplicated keeps the first copy for each strain: reference strain that carry many proteins in a marker, the retain one is by file order
  keep <- keep_taxon(taxa) & !duplicated(taxa) 
  # subset the sequence, relabel by strain
  seqs <- as.character(aln)[keep] 
  # keep record of the width to know how many - gap characters if protein is missing 
  one_og[[f]] <- list(seqs = seqs, width = width(aln)[1])
  all_taxa  <- union(all_taxa, names(seqs))
}
# seqs after the subset and relabelling: 
# GCA_001752665   MWDRRFFAWFSVAA--
# Ace -WDRRFFAWFSVAA--
# Acromyrmex_octospinosus_..._magscot_1 MWDRRFFAWFSVAA--
# three rows relabelled per taxon stored at one_og[[f]]

# Build supermatrix concatenation and fill with gaps 
# if all_taxa = four names and there are three markers:
# all_taxa = c("GCA_001752665", "Ace", "MAG_A", "MAG_B") with three markers of width 18, 9 and 10, and MAG_B is missing marker 1, MAG_A is missing marker 2
# cols <- vector("list", length(aln_files)) makes an empty list of the right length: cols [[1]] NULL [[2]] NULL [[3]] NULL 
# f <- aln_files[i] and s <- one_og[[f]]$seqs where the filename is the lookup key, on i = 1:
# s
# GCA_0001752666  Ace MAG_A
# "MWSVDRR--" "-WSVDRR--" "MWSVDRR-A"
# MAG_B is not there because it had no protein for marker 1
# gap <- paste(rep("-", one_og[[f]]$width), collapse = "") rep("-", 18) gives 18 separate - characters, and collapse = "" joins them into one string which makes it as wide as marker 1 which is why we needed width above
# cols[[i]] <- ifelse(all_taxa %in% names(s), s[all_taxa], gap) ifelse that walks through all_taxa and picks: all_taxa is in s? TRUE -> sequence, FALSE -> ------
# After three iterations: 
# cols[[1]] cols[[2]] cols[[3]]
# "MWXIAAFNEKJF--"  "--HAIDMF"  "TTGACRRW-"
# "MWXIAAFNEKJF--"  "--HAIDMF"  "TTTGCRRW-"
# "MWXIAAFNEKJF--"  "--------"  "TTGACRQW-"
# "--------------"  "--HAIDMF"  "TTGFCRRW-"
# Join elements: 
# "MWXIAAFNEKJF----HAIDMFTTGACRRW-"
# "MWXIAAFNEKJF----HAIDMFTTTGCRRW-"
# "MWXIAAFNEKJF----------TTGACRQW-"
# "----------------HAIDMFTTGFCRRW-"
# Then: super <- AAStringSet(do.call(paste0, cols)) puts it into a Biostrings object and reattaches the labels since we still have all_taxa order in every columns
# writeXStringSet(super, paste0(PREFIX, "all.fasta")):
# > GCA_001752665
# MWXIAAFNEKJF----HAIDMFTTGACRRW-
# >Ace
# MWXIAAFNEKJF----HAIDMFTTTGCRRW-
# >MAG_A
# "MWXIAAFNEKJF----------TTGACRQW-"
# >MAG_B
# "----------------HAIDMFTTGFCRRW-"
# Supermatrix in the end is one long chain: pasting every marker protein end to end into a single sequence for each strain
cols <- vector("list", length(aln_files))
for (i in seq_along(aln_files)) {
  f <- aln_files[i]
  s <- one_og[[f]]$seqs
  gap <- paste(rep("-", one_og[[f]]$width), collapse = "")
  cols[[i]] <- ifelse(all_taxa %in% names(s), s[all_taxa], gap)
}
super <- AAStringSet(do.call(paste0, cols))
names(super) <- all_taxa

writeXStringSet(super, paste0(PREFIX, ".all.fasta"))