#!/usr/bin/env Rscript

# Merge overlapping intervals and compute the repeat content of each contig
# Takes the normalized intervals and the coverage classification to produce one row for each contig with base pair and fractions for each repeat category
# Rscript D3b_merge_intervals.R species intervals coverage_table outdir

library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
INTERVALS <- args[2]
COVERAGE_TABLE <- args[3]
OUTDIR <- args[4]

# Fraction of a contig annotated as repetitive:
TE_THRESHOLD <- 0.5

# Merge overlapping intervals
# Return the total base pairs covered
# Sort by start and wak trhough conting to extend or close the interval
merge_and_sum <- function(starts, ends){
    if(length(starts) == 0) return(0)
    ord <- order(starts)
    s <- starts[ord]
    e <- ends[ord]
    current_start <- s[1]
    current_end <- e[1]
    total <- 0
    if(length(s) > 1){
        for(i in 2:length(s)) {
            if( s[i] <= current_end){ # overla with current interval: extend
                current_end <- max(current_end, e[i])
            } else { # no overlap: close the current interval and start new one
                total <- total + (current_end - current_start + 1)
                current_start <- s[i]
                current_end <- e[i]
            }
        }
    }
    total + (current_end - current_start + 1) # for the last interval
}
# Apply to contigs and return a vector numeric aligned to contigs
base_pair_by_contig <- function(records, contigs) {
    if(nrow(records) == 0) return(rep(0, length(contigs)))
    by_contig <- split(records[, c("start", "end")], records$contig)
    vapply(contigs, function(contig){
        if (contig %in% names(by_contig)){
            distance <- by_contig[[contig]]
            merge_and_sum(distance$start, distance$end)
        } else 0
    }, numeric(1))
}

# Read intervals
intervals <- read_tsv(INTERVALS, show_col_types = FALSE)
cov <- read_tsv(COVERAGE_TABLE, show_col_types = FALSE)
contigs <- cov$contig

# Base per for each category
te <- intervals %>%
    filter(category == "te")
non_te <- intervals %>%
    filter(category != "te")

base_pair_contig <- data.frame(contig = contigs, stringsAsFactors = FALSE)
base_pair_contig$te_bp <- base_pair_by_contig(te, contigs)
base_pair_contig$non_te_bp <- base_pair_by_contig(non_te, contigs)
base_pair_contig$total_repeat_bp <- base_pair_by_contig(intervals, contigs)

# Join to coverage and fractions
repeat_cov <- cov %>%
    left_join(base_pair_contig, by = "contig")
base_pair_cols <- c("te_bp", "non_te_bp", "total_repeat_bp")
for (col in base_pair_cols) repeat_cov[[col]] <- pmin(repeat_cov[[col]], repeat_cov$len)

repeat_cov <- repeat_cov %>% 
    mutate(
        te_fraction = te_bp / len,
        non_te_fraction = non_te_bp / len,
        total_repeat_fraction = total_repeat_bp / len,
        flagged_repetitive = ifelse(total_repeat_fraction >= TE_THRESHOLD, "yes", "no")
    ) %>%
    arrange(desc(len))

write_tsv(repeat_cov, file.path(OUTDIR, "contig_repeat_coverage.tsv"))
