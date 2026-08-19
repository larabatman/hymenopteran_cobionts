#!/usr/bin/env Rscript

# Collects the each assigned supergroup call into one long table: one row per assigned NUWT. 
# That table is the input to annotate_nuwt_fragments.R, which attaches the OG annotation and marks the fragments excluded from the region analysis

# The "?" label means the tree could not resolve the NUWT to a single supergroup clade
#
# Input files:
# 1: <assignments_dir>/<OG>.assign.tsv with for columns as og tab INS tab supergroup tab nuwt_name and the nuwt_name = INS_<Species>_<scaffold>_<OG>_<start>-<stop>
# 2: <results_dir>/ Needed because both species names and scaffold names contain '_', so <Species> cannot be split out of nuwt_name
#
# Output files_
# <prefix>_per_nuwt.tsv with one row per assigned NUWT: name, species, supergroup, og, length
# Where the NUWT length is a string to integer extraction from its header in nuwt_name<start>-<stop>
#
# Usage
# Rscript aggregate_supergroups.R <assignments_dir> <results_dir> <out_prefix>
args <- commandArgs(trailingOnly = TRUE)

assign_dir <- args[1] # placement/assignments
results_dir <- args[2] # nuwt_scan/results where sub-dirs are species names
out_prefix <- args[3]

# Species set
# Longest first so the regex alternation below prefers the longest matching species name, rather than a short name matching as a prefix of a longer one.
species_all <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
species_all <- species_all[nchar(species_all) > 0]
species_all <- species_all[order(nchar(species_all), decreasing = TRUE)]

# Read assignments
files <- list.files(assign_dir, pattern = "\\.assign\\.tsv$", full.names = TRUE)
files <- files[file.size(files) > 0]

# Stack all the files together: appy the function to each element of files, collect the results in a list
# A list of small data frame, one per OG that was placed
# Each gets read.tabled() and they are all headerless 4 column TSVs with wuote and comment disabled , everything as text
dat <- do.call(rbind, lapply(files, function(f) # do.call(rbind) calls rbind with the list's elements as its arguments: rbind(df1, df2, ... df1871) which stacks them all in one call rather than looping
  read.table(f, sep = "\t", header = FALSE, quote = "", comment.char = "",
             col.names = c("og", "ins", "supergroup", "name"),
             colClasses = "character")))
# Result: one data frame, one row per assigned NUWT across the whole run

# Retrieve species via header
sp_pat <- paste0("^INS_(", paste(species_all, collapse = "|"), ")_.*$")
dat$species <- sub(sp_pat, "\\1", dat$name)

# Summarize length from header as well: from the coordinates in the name that ends with _<start>-<stop>
ok_fmt <- grepl("_[0-9]+-[0-9]+$", dat$name)
if (!all(ok_fmt))
  stop(sprintf("nuwt_name(s) do not end in _<start>-<stop>, e.g.:\n  %s",
               paste(head(dat$name[!ok_fmt], 3), collapse = "\n  ")))

frag_start <- as.numeric(sub(".*_([0-9]+)-([0-9]+)$", "\\1", dat$name))
frag_stop <- as.numeric(sub(".*_([0-9]+)-([0-9]+)$", "\\2", dat$name))
dat$length <- abs(frag_stop - frag_start) + 1

# Output
out_file <- paste0(out_prefix, "_per_nuwt.tsv")
write.table( dat[, c("name", "species", "supergroup", "og", "length")], out_file, sep = "\t", quote = FALSE, row.names = FALSE)