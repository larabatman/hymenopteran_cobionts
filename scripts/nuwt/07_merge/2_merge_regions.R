#!/usr/bin/env Rscript

# Merge nuwt fragments into regions: fragments within 2.5 kbp on the same contig or scaffold are merged 
# Usage:
# Rscript 2_merge_regions.R kept_nuwts.tsv fragments.tsv out_prefix

args <- commandArgs(trailingOnly = TRUE)

KEPT_NUWTS <- args[1]
FRAGMENTS_TSV <- args[2]
OUT_PREFIX <- args[3]
# Merging gap 
GAP <- 2500L

# Function to read tables
read_tsv0 <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "",comment.char = "", stringsAsFactors = FALSE)

# Join supergroups to their coordinates
calls <- read_tsv0(KEPT_NUWTS)[, c("name", "species", "supergroup", "length")]
coord <- read_tsv0(FRAGMENTS_TSV)[, c("name", "scaffold", "start", "end")]
data <- merge(calls, coord, by = "name")

# Prepare input for bedtools
# Needs a table of intervals as chrom  start end that bedtools merge will take to merge intervals
# Abia_candens::OZ125706.1  15060638    15061970    A   1332
# Abia_candens::OZ125706.1  15062100    15062940    A   840
# Abia_candens::OZ125706.1  15071500    15072300    B   800
bed <- data.frame(chrom = paste0(data$species, "::", data$scaffold), start = data$start - 1L, end = data$end, superg = data$supergroup, len = data$length, stringsAsFactors = FALSE)
bed <- bed[order(bed$chrom, bed$start), ]
bed_f <- tempfile(fileext = ".bed")
out_f <- tempfile(fileext = ".bed")
write.table(bed, bed_f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# bedtools call 
# -c 4,5 carries columns 4 and 5 through 
# -o collapse makes a comma separated list
# Abia_candens::OZ125706.1  15060638    15061970    A,A   1332,840 for a merged scaffold
if (system2("bedtools", c("merge", "-d", GAP, "-i", bed_f, "-c", "4,5", "-o", "collapse,collapse"), stdout = out_f) != 0)
  stop("bedtools merge failed")

bed_merged <- read.table(out_f, sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE, col.names = c("chrom", "start0", "end", "sgs", "lens"))

# Function to annotate the regions
# The region gets the supergroup that contributed the most aligned bp, mixed when there is more than one supergroup present 
# sgs = "A,B,A" lens = "1332,840,1200" has three fragments merged, so strsplit(sgs, ",")[[1]] splits on commas and gives c("A", "B", "A")
# as.numerics(strsplit(lens, ","))[[1]] is the same then converted to numbers for the sum
# then we split values by grouping vector and apply the sum to each group: tapply(lens, sgs, sum): A gets 2532 and B gets 840
dominant <- function(sgs, lens) {
  bp <- tapply(as.numeric(strsplit(lens, ",")[[1]]), strsplit(sgs, ",")[[1]], sum)
  names(bp)[which.max(bp)] }

# BuilD regions data frame
# Onerow for each merged region with seven columns
# start is converted to 1 based 
regions <- data.frame(
  species = sub("::.*$", "", bed_merged$chrom),
  scaffold = sub("^.*::", "", bed_merged$chrom),
  start = bed_merged$start0 + 1L,
  end = bed_merged$end,
  span_bp = bed_merged$end - bed_merged$start0,
  n_fragments = lengths(strsplit(bed_merged$sgs, ",")),
  supergroup = mapply(dominant, bed_merged$sgs, bed_merged$lens, USE.NAMES = FALSE),
  stringsAsFactors = FALSE)
# Mixed when more than one supergroup
regions$mixed <- vapply(strsplit(bed_merged$sgs, ","), function(x) length(unique(x)) > 1, logical(1))
regions <- regions[order(regions$species, regions$scaffold, regions$start), ]
write.table(regions, paste0(OUT_PREFIX, "_regions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Build long format table of region counts for each specis and for each supergroup, a count table 
# Cross tabulation of every combination of species and supergroups
# species   sg_vec  Freq
# Abia_candense A   14
# Nomada_fabriciana A   3
counts <- as.data.frame(table(species = regions$species, sg_vec = regions$supergroup), stringsAsFactors = FALSE)
counts <- counts[counts$Freq > 0, ]
write.table(counts, paste0(OUT_PREFIX, "_region_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)