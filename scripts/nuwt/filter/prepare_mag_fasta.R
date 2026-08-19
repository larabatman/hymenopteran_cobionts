#!/usr/bin/env Rscript


# Looks up one species in the curated Wolbachia MAG sheet and writes its MAG ready for makeblastdb.
#
# The sheet lists one row per recovered Wolbachia genome. Species, assembler and refinement tool together determine the bin path:
# results/<Species>/full_binners/bins/fasta/<tool>/<Species>_<assembler>_<tool>_*.fa.gz
#
# If the species is absent from the sheet it has no Wolbachia MAG.
#
# Files:
# in: sheet .tsv, tab-separated with a header row.
# in: MAG .fa.gz, gzipped nucleotide FASTA, read directly via gzfile()
# out: .fa, plain nucleotide FASTA.
#
# Headers are prefixed mag0_, mag1_ ... because dastool and magscot bin from the same assembly and so produce identical contig IDs.
# Usage:
# Rscript prepare_mag_fasta.R <Species> <sheet.tsv> <cobionts_root> <out.fa>

library(readr)
library(dplyr)

species <- args[1]
sheet <- args[2]
root <- args[3]
out_fa <- args[4]

# Species row
# Columns are selected by name, so the sheet can be regenerated with a
# different column order without breaking this.
mags <- read_tsv(sheet, show_col_types = FALSE) %>%
  rename(sheet_species = `Host species`,
         assembler = Assembler,
         tool = `Refinement Tool`) %>%
  mutate(sheet_species = gsub(" ", "_", sheet_species))  %>%
  filter(sheet_species == species)  %>%
  select(assembler, tool)

# Bin files

bin_dir <- file.path(root, "results", species, "full_binners", "bins", "fasta")

mag_files <- mags %>%
  mutate(pattern = file.path(bin_dir, tool, sprintf("%s_%s_%s_*.fa.gz", species, assembler, tool))) %>%
  pull(pattern) %>%
  lapply(Sys.glob) %>%
  unlist()

# Combined FASTA
# gzfile() decompresses and sub() rewrites only the leading ">" of each header
lines <- unlist(lapply(seq_along(mag_files), function(i) {
  prefix <- sprintf(">mag%d_", i - 1L)
  sub("^>", prefix, readLines(gzfile(mag_files[i])))
}))

writeLines(lines, out_fa)