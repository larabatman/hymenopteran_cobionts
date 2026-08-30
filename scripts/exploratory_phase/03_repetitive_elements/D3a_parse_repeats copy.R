#!/usr/bin/env Rscript

# Parse three repeat annotation tools into one table with interals
# EDTA GFF3, RepeatMasker .out and TRF .dat describe intervals on contigs.
# Normalize to the same four columns for merging downstream
# Usage:
# Rscript D3a_parse_repeats.R species EDTA_GFF3 RM.out TRF.dat outdir

library(readr)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
EDTA_GFF3 <- args[2]
RM_OUT_FILE <- args[3]
TRF_DAT_FILE <- args[4]
OUTDIR <- args[5]

# EDTA GFF3 parsing: drop the subfeatures 
gff <- read.delim(EDTA_GFF3, header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE, col.names = c("seqid", "source", "tpye", "start", "end", "score", "strand", "phase", "attributes"))
exclude_types <- c("long_terminal_repeat", "target_site_duplication", "repeat_region")
gff <- gff[!gff$type %in% exclude_types, ]

edta <- data.frame(
    contig = gff$seqid,
    start = as.integer(gff$start),
    end = as.integer(gff$end),
    source = "edta",
    category = "te",
    stringsAsFactors = FALSE
)

# RepeatMasker .out: three header lines then one whitespace hits for each line 
# column 5 contains the query name, 6 the start, 7 the end and 11 repeat class
rm_lines <- readLines(RM_OUT_FILE)
rm_lines <- rm_lines[-(1:3)]
rm_lines <- rm_lines[nchar(trimws(rm_lines)) > 0]
# Split on withespace and keep tje complete lines to build a matrix
fields <- strsplit(trimws(rm_lines), "\\s+")
fields <- fields[lengths(fields) >= 11] # lengths gives the length of each element
mat <- matrix(unlist(lapply(fields, `[`, c(5, 6, 7, 11))), ncol = 4, byrow = TRUE)

# Classify hits: if it does not match a keyword that is TE, it is treated as TE
rm_class <- mat[, 4]
rm_category <- ifelse(grepl("Simple_repeat", rm_class), "simple_repeat",
                ifelse(grepl("Low_complexity", rm_class), "low_complexity",
                ifelse(grepl("Satellite", rm_class), "satellite",
                ifelse(grepl("rRNA|tRNA|snRNA|scRNA", rm_class), "rna",
                "te"))))
repeat_masker_table <- data.frame(
    contig = mat[, 1],
    start = as.integer(mat[, 2]),
    end = as.integer(mat[, 3]),
    source = "rm",
    category = rm_category,
    stringsAsFactors = FALSE
)

# TRF .dat parsing: has "Sequence: contig" markers with blocks of line that start with two integers for the start and end of the tandem repeat
# cumsum
trf_lines <- readLines(TRF_DAT_FILE)
is_seq <- grepl("^Sequence:", trf_lines)
seq_names <- sub("^Sequence:\\s+(\\S+).*", "\\1", trf_lines[is_seq])
block <- cumsum(is_seq) # 0 before the first marker, then 1 inside the first block of lines, 2 inside the second and so on
contig_of_line <- c(NA_character_, seq_names)[block + 1]
is_data <- grepl("^\\s*[0-9]+\\s+[0-9]+\\s", trf_lines) & !is.na(contig_of_line)
trf_columns <- strsplit(trimws(trf_lines[is_data]), "\\s+")
trf_table <- data.frame(
    contig = contig_of_line[is_data],
    start = as.integer(sapply(trf_columns, `[`, 1)),
    end = as.integer(sapply(trf_columns, `[`, 2)),
    source = "trf",
    category = "tandem_repeat",
    stringsAsFactors = FALSE
)

# Put everything together
intervals <- rbind(edta, repeat_masker_table, trf_table)
intervals <- intervals[!is.na(intervals$start) & !is.na(intervals$end) & intervals$end >= intervals$start, ]

write_tsv(intervals, file.path(OUTDIR, "repeat_intervals.tsv.gz"))