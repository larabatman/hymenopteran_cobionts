#!/usr/bin/env Rscript

# Collect the assigned supergroup calls into a table with one row for each assigned NUWT
# the "?" labels are from the reference genomes that had no supergroup information 
# Usage
# Rscript 1_aggregate_supergroups.R assignment_directoy results_directory out_prefix

args <- commandArgs(trailingOnly = TRUE)

# Directory with the assignments
ASSIGN_DIR <- args[1]
# Result folder with the species names 
RESULTS_DIR <- args[2]
OUT_PREFIX <- args[3]

# Build the list of species name from the sub folders in the result directories
species_all <- list.dirs(RESULTS_DIR, full.names = FALSE, recursive = FALSE)
# Drop empty strings
species_all <- species_all[nchar(species_all) > 0]

# Read assignment files
files <- list.files(ASSIGN_DIR, pattern = "\\.assign\\.tsv$", full.names = TRUE)
files <- files[file.size(files) > 0]

# Put all the assignment files together
# Apply rbind to each element of the files collecting the results in a list
# Get one data frame with one row for each assigned nuwt across all assignments
data <- do.call(rbind, lapply(files, function(f) read.table(f, sep = "\t", header = FALSE, quote = "", comment.char = "", col.names = c("og", "ins", "supergroup", "name"), colClasses = "character")))

# Parse the headers to get the species
sp_pat <- paste0("^INS_(", paste(species_all, collapse = "|"), ")_.*$")
data$species <- sub(sp_pat, "\\1", data$name)

# Get the length from the coordinates: name ends with _start-stop
frag_start <- as.numeric(sub(".*_([0-9]+)-([0-9]+)$", "\\1", data$name))
frag_stop <- as.numeric(sub(".*_([0-9]+)-([0-9]+)$", "\\2", data$name))
data$length <- abs(frag_stop - frag_start) + 1

# Output
out_file <- paste0(OUT_PREFIX, "_supergroup_nuwt.tsv")
write.table( data[, c("name", "species", "supergroup", "og", "length")], out_file, sep = "\t", quote = FALSE, row.names = FALSE)