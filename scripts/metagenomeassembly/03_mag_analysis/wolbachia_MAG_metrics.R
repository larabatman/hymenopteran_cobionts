#!/usr/bin/env Rscript

# Builds the Wolbachi MAG figure 
# Takes the genomic properties of the 98 MAGs in the Wolbachia_MAGs.tsv file
# Draws the MIMAG completeness vs contamination score, the GC distribution, the genome vs GC plot and mean coverage depth vs fragmentation
# Usage:
# Rscript wolbachia_MAG_metrics.R wolbachia_MAGs outdir wolbachia_ref_metrics

library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

WOL_TSV <- args[1]
OUTDIR <- args[2]
REF_TSV <- args[3]

COL_POINT <- "#0072B2"

# Read Wolbachi MAG
wol_mag_table <- read.delim(WOL_TSV, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA",""))
cols <- c(species = "Host species",
          assembler = "Assembler",
          refiner = "Refinement Tool",
          n_contigs = "Number of cobiont contigs",
          length_bp = "Total cobiont genome length [bp]",
          gc = "GC [%]",
          depth = "Mean depth",
          completeness = "Completeness [%]",
          contamination = "Contamination [%]",
          r16s = "Total number of 16S rRNA",
          r23s = "Total number of 23S rRNA")

data_wol_mags <- wol_mag_table[cols]
names(data_wol_mags) <- names(cols)
clean_col <- c("n_contigs","length_bp","gc","depth","completeness","contamination","r16s","r23s")
data_wol_mags[clean_col] <- lapply(data_wol_mags[clean_col], as.numeric)

# Set a theme for the figure 
theme_set(theme_bw(base_size = 9) + 
  theme(panel.grid.minor = element_blank(), 
  panel.grid.major = element_line(linewidth = .25, colour = "grey92"), 
  panel.border = element_rect(linewidth = .4, colour = "grey30"),
  plot.title = element_text(size = 9, face = "bold", hjust = 0),
  legend.key.size = unit(3.2, "mm"), legend.title = element_text(size = 8),
  legend.text = element_text(size = 7.5)))

# Published complete Wolbachia reference genome
ref_wolb <- read.delim(REF_TSV, stringsAsFactors = FALSE)
# Mean GC content
ref_wolb_gc <- mean(ref_wolb$gc)

# MIMAG space: completeness vs contamination
plot_mimag <- ggplot(data_wol_mags, aes(completeness, contamination)) +
  annotate("rect", xmin = 90, xmax = 100, ymin = -.3, ymax = 5, fill = "grey88", alpha = .55) +
  geom_hline(yintercept = 5, linetype = "dashed", linewidth = .35) +
  geom_vline(xintercept = 90, linetype = "dashed", linewidth = .35) +
  geom_point(colour = COL_POINT, size = 1.6, alpha = .8) +
  scale_x_continuous(breaks = c(90, 95, 100)) +
  coord_cartesian(xlim = c(89, 100), ylim = c(0, max(data_wol_mags$contamination) + .4)) +
  labs(title = "A  Completeness vs contamination", x = "Completeness (%)", y = "Contamination (%)")

# GC distribution
plot_gc <- ggplot(data_wol_mags, aes(gc)) +
  geom_histogram(binwidth = .15, fill = COL_POINT, colour = "white", linewidth = .2) +
  geom_vline(xintercept = ref_wolb_gc, linetype = "dashed", linewidth = .35) +
  labs(title = "B  Distribution of GC content", x = "GC (%)", y = "MAG count")

# Genome vs GC
plot_genome_gc <- ggplot(data_wol_mags, aes(length_bp/1e6, gc)) +
  geom_point(aes(size = n_contigs), colour = COL_POINT, alpha = .8) +
  geom_point(data = ref_wolb, aes(length_bp/1e6, gc), inherit.aes = FALSE, shape = 23,size = 2.4, fill = "white", colour = "black", stroke = .6) +
  scale_size_continuous(range = c(1, 4.2), name = "Contig count") +
  labs(title = "C  Genome length vs GC", x = "Genome length (Mb)", y = "GC (%)")

# Mean coverag depth vs fragmentation
plot_mean_frag <- ggplot(data_wol_mags, aes(depth, n_contigs)) +
  geom_point(colour = COL_POINT, size = 1.6, alpha = .8) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "D  Mean coverage vs fragmentation", x = "Mean depth (log10)", y = "Contig count")

fig <- (plot_mimag | plot_gc) / (plot_genome_gc | plot_mean_frag)
ggsave(file.path(FIG_DIR, "Fig1_wolbachia_metrics.pdf"), fig, width = 180, height = 140,units = "mm")