#!/usr/bin/env Rscript

# Builds three blobplots from the BUSCO contig table, coloured by BUSCO density per Mbp
# Usage:
# Rscript C2b_busco_plots.R species busco_per_contig.tsv outdir

library(dplyr)
library(readr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
BUSCO <- args[2]
OUTDIR <- args[3]

busco <- read_tsv(BUSCO, show_col_types = FALSE)

# Fill colors:
COV_COLOURS <- c(host_like = "steelblue", coverage_outlier = "firebrick", ambiguous = "grey70")
BP_PER_MBP <- 1e6

# Drop contigs with 0 coverage or length
# Drop contigs below 13 or above 80% GC
busco_plot <- busco %>%
    filter(mean_cov > 0, len > 0, gc >= 13, gc <=80)

# Density blobplots: divid by length to have comparable contigs:
df_gradient <- busco_plot %>%
    mutate( 
        hym_per_mbp = n_hym / (len / BP_PER_MBP),
        arth_per_mbp = n_arth / (len / BP_PER_MBP),
        bact_per_mbp = n_bact / (len / BP_PER_MBP)
    )
# Contigs without BUSCO are grey in background
make_gradient_blobplot <- function(data, density_colour, lineage_name, outfile) {
    zero_density <- data %>%
        filter(.data[[density_colour]] == 0 | is.na(.data[[density_colour]]))
    nonzero_density <- data %>% # sort ascending: the contigs with most density are drawn last to stay on top
        filter(.data[[density_colour]] > 0) %>%
        arrange(.data[[density_colour]])
    if(nrow(nonzero_density) == 0) {
        message("No ", lineage_name, " BUSCOs, skip ", outfile)
        return(invisible(NULL))
    }
    # Colour scale cap at the 99th percentile, so that the gradient stays
    cap <- quantile(nonzero_density[[density_colour]], 0.99, na.rm = TRUE)
    if (is.na(cap) || cap == 0) cap <- max(nonzero_density[[density_colour]], na.rm = TRUE)

    bloplot <- ggplot() +
        geom_point(data = zero_density, aes(gc, mean_cov, size = log10(pmax(len, 1))), colour = "grey85", alpha = 0.4) +
        geom_point(data = nonzero_density, aes(gc, mean_cov, size = log10(pmax(len, 1)), colour = pmin(.data[[density_colour]], cap)), alpha = 0.8) +
        scale_y_log10() +
        scale_size(range = c(1, 8)) +
        scale_colour_viridis_c(name = "BUSCO / Mbp", option = "viridis") +
        theme_bw() +
        labs(title = paste0(lineage_name, " BUSCO density: ", SPECIES), x = "GC (%)", y = "Mean coverage", size = "log10(contig length)")
    ggsave(outfile, blobplot, width = 8, height = 5)
}
# Make blobplots
make_gradient_blobplot(df_gradient, "hym_per_mbp", "Hymenoptera", file.path(OUTDIR, "blobplot_hym_gradient.pdf"))
make_gradient_blobplot(df_gradient, "arth_per_mbp", "Arthropoda", file.path(OUTDIR, "blobplot_arth_gradient.pdf"))
make_gradient_blobplot(df_gradient, "bact_per_mbp", "Bacteria", file.path(OUTDIR, "blobplot_bact_gradient.pdf"))

# Tables:
blobplot_summary <- busco_plot %>%
    summarise(
        species = SPECIES, 
        n_plotted = n(), 
        n_excluded = nrow(busco) - n(),
        hym = sum(has_hym),
        arth = sum(has_arth),
        bact = sum(has_bact),
        mixed_hym_bact = sum(has_hym & has_bact),
        no_busco = sum(n_hym == 0 & n_arth == 0 & n_bact == 0)
    )
write_tsv(blobplot_summary, file.path(OUTDIR, "blobplot_summary.tsv"))

