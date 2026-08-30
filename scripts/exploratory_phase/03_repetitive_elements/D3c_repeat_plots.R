#!/usr/bin/env Rscript

# Make three blobplots with GC, coverage and point size by contig length. The colours are added for the repeat fraction
# Usage:
# Rscript D3c_repeat_plots.R species repeat_table outdir

library(dplyr)
library(readr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
REPEAT_TABLE <- args[2]
OUTDIR <- args[3]

repeats <- read_tsv(REPEAT_TABLE, show_col_types = FALSE)

# Drop contigs with 0 coverage or length and GC below 13 or above 80%
repeat_plot <- repeats %>%
    filter(mean_cov > 0, len > 0, gc >= 13, gc <= 80)

make_repeat_blobplot <- function(data, repeat_fraction, label, colour_high, outfile){
    data <- data %>% arrange(.data[[repeat_fraction]])
    bobplot <- ggplot(data, aes(x = gc, y = mean_cov, size = log10(pmax(len, 1)), colour = .data[[repeat_fraction]])) +
        geom_point(alpha = 0.7) +
        scale_y_log10() +
        scale_size(range = c(1, 8)) +
        scale_colour_gradient( low = "grey80", high = colour_high, name = label, limits = c(0, 1)) +
        theme_bw() +
        labs(title = paste0(label,": ", SPECIES), x = "GC (%)", y = "Mean coverage (log10)", size = "log10(contig length)")
    ggsave(outfile, bobplot, width = 8, height = 5, dpi = 300)
}

make_repeat_blobplot(repeat_plot, "te_fraction", "TE fraction (EDTA)", "firebrick", file.path(OUTDIR, "blobplot_te.pdf"))
make_repeat_blobplot(repeat_plot, "non_te_fraction", "Non-TE repeat fraction", "darkorange", file.path(OUTDIR, "blobplot_non_te_repeats.pdf"))
make_repeat_blobplot(repeat_plot, "total_repeat_fraction", "Total repeat fraction", "purple4", file.path(OUTDIR, "blobplot_total_repeats.pdf"))