#!/usr/bin/env Rscript

# Collapses NUWT fragments into regions
#
# Fragments within 2.5 kbp on the same scaffold are merged. 
# Originally, in Vancaester et al., the merging rule also needed no coding gene sits between them but we don't have the gene annotations 
#
# Input
# 1: the kept_per_nuwt.tsv that came from 5_og_mito_ank.R with : name, species, supergroup, length. One row per OG gene assigned fragment
# 2: nuwt_fragments.tsv from extract_nuwt_seqs.py: name, scaffold, start, end. This is joined on name for the coordinates.
#
# Output
# 1: _regions.tsv with one row per region
# 2: _region_counts.tsv a long table with species, sg_vec, Freq which is the input plot_fig2_composite.R expects
#
# USAGE 
# Rscript merge_nuwt_regions.R <kept_per_nuwt.tsv> <fragments_tsv> <out_prefix>

args <- commandArgs(trailingOnly = TRUE)


per_nuwt <- args[1]
frags_tsv <- args[2]
out_prefix <- args[3]
gap <- 2500L

# Read tables function
read_tsv0 <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "",comment.char = "", stringsAsFactors = FALSE)

# Join the supergroups to theri coordinates by their names
calls <- read_tsv0(per_nuwt)[, c("name", "species", "supergroup", "length")]
coord <- read_tsv0(frags_tsv)[, c("name", "scaffold", "start", "end")]
d <- merge(calls, coord, by = "name")

# Prepare BED input:
# BED is a table of intervals. bedtool merge akes that table and glues overlaping or nearby intervals together
# A BED file goes like this: chrom  tsart   end
# What our bed should look like:
# Abia_candens::OZ125706.1  15060638    15061970    A   1332
# Abia_candens::OZ125706.1  15062100    15062940    A   840
# Abia_candens::OZ125706.1  15071500    15072300    B   800
# chrom is species::scaffold so intervals from different species or scaffolds are not merged together
# What merge -d 2500 does is it walks the intervals in order, and if the next one starts within 2500 bp of where the current block ends, then it gets absorbed otherwise the block is closed and we start a new one
# -c 4,5 -o collapse,collapse: -c 4, 5 tells to carry columns 4 and 5 through otherwise they are thrown away, then -o collapse makes the as a comma-separated lis, so that we have the record of what was there as in: # Abia_candens::OZ125706.1  15060638    15061970    A,A   1332,840 for that scaffold that has been merged
# which also allows after to choose the dominant supergroup!
bed <- data.frame(chrom = paste0(d$species, "::", d$scaffold), start = d$start - 1L, end = d$end, sg = d$supergroup, len = d$length, stringsAsFactors = FALSE)
bed <- bed[order(bed$chrom, bed$start), ] # bedtools requires sorted input
bed_f <- tempfile(fileext = ".bed")
out_f <- tempfile(fileext = ".bed")
write.table(bed, bed_f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

if (system2("bedtools", c("merge", "-d", gap, "-i", bed_f, "-c", "4,5", "-o", "collapse,collapse"), stdout = out_f) != 0)
  stop("bedtools merge failed")

m <- read.table(out_f, sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE, col.names = c("chrom", "start0", "end", "sgs", "lens"))

# Annotate its region with whichever contributed the most aligne bp
# Flagged mixed when there is more than one supergroup present 
# This takes one region that has collapsed strings and returns a single supergroup label:
# sgs = "A,B,A"
# lens = "1332,840,1200"
# Had three fragments merged, so strsplit(sgs, ",")[[1]] splits on commas and gives c("A", "B", "A")
# as.numerics(strsplit(lens, ","))[[1]] is the same then converted to numbers for the sum
# then we split values by grouping vector and apply the sum to each group: tapply(lens, sgs, sum): A gets 2532 and B gets 840
dominant <- function(sgs, lens) {
  bp <- tapply(as.numeric(strsplit(lens, ",")[[1]]), strsplit(sgs, ",")[[1]], sum)
  names(bp)[which.max(bp)] # and we chose the max: which.max is the position, and names is its name 
}
# One row per merged region, seven columns
# species: we need to delete :: and everything after
# scaffold: the other way, we keep everything after ::
# start: Convert BED 0-based back to 1-based convention that we have been using before, L for integers
# end stays the same
# spans_bp is end-start0, the region width 
# n_fragments splits the collapsed fragments and counts the items: lengths give the length of each element of a list so it counts how many supergroups where there
# supergroup: the dominant one 
regions <- data.frame(
  species = sub("::.*$", "", m$chrom),
  scaffold = sub("^.*::", "", m$chrom),
  start = m$start0 + 1L,                   # back to 1-based
  end = m$end,
  span_bp = m$end - m$start0,
  n_fragments = lengths(strsplit(m$sgs, ",")),
  supergroup = mapply(dominant, m$sgs, m$lens, USE.NAMES = FALSE),
  stringsAsFactors = FALSE)
# add mixed column with the strsplit pass
regions$mixed <- vapply(strsplit(m$sgs, ","), function(x) length(unique(x)) > 1, logical(1))
# order for pretty
regions <- regions[order(regions$species, regions$scaffold, regions$start), ]
write.table(regions, paste0(out_prefix, "_regions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Long format: region counts per species per supergroup
# Turns the region table ino the count table for the figures
# cross-tabulates two vectors: for every combination of species and supergroups, how many regions have that pair
# Result: matrix with species across supergroups
# species   sg_vec  Freq
# Abia_candense A   14
# Nomada_fabriciana A   3
# where Freq is the count as in, Abia candens has 14 merged insertion regiosn with dominant supergorup A
# Each region is counted once under its dominant supergruop only
# Here we also make the choice to drop zeros as to not overfeed ggplot with 0s
counts <- as.data.frame(table(species = regions$species,
                              sg_vec  = regions$supergroup),
                        stringsAsFactors = FALSE)
counts <- counts[counts$Freq > 0, ]
write.table(counts, paste0(out_prefix, "_region_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)